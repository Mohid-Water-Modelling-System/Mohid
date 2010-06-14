!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : LevitusFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert LevitusFormat files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleLevitusFormat

    use ModuleTime
    use ModuleGlobalData
    use ModuleFunctions
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertLevitusFormat
    private ::      ReadGlobalOptions
    private ::      ReadLevitusAnnualFile
    private ::      ReadLevitusMonthlyFile
    private ::          AddField
    private ::          CheckName
    private ::      OutputFields
    private ::          Open_HDF5_OutPut_File
    private ::      KillLevitusFormat

    !Parameters----------------------------------------------------------------
    integer, parameter                          :: NLevitusAnnual  = 33
    integer, parameter                          :: NLevitusMonthly = 10

    real,    parameter                          :: MaxDepth = 12000

    character(LEN = StringLength), parameter    :: field_block_begin   = '<<beginfield>>'
    character(LEN = StringLength), parameter    :: field_block_end     = '<<endfield>>'
    character(LEN = StringLength), parameter    :: input_files_begin   = '<<<begin_input_files>>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<<end_input_files>>>'


    !Types---------------------------------------------------------------------
    

    private :: T_Field
    type       T_Field
        type(T_PropertyID)                      :: ID
        character(len=StringLength)             :: Units
        integer                                 :: GridLocation
        real, dimension(:,:,:  ),   pointer     :: AnnualValues3D
        real, dimension(:,:,:  ),   pointer     :: MonthlyValues3D
        character(len=StringLength), dimension(:), pointer :: InputFiles
        character(len=StringLength)             :: AnnualFile
        logical                                 :: FirstField
        type(T_Field),              pointer     :: Next, Prev
    end type  T_Field

    
    private :: T_LevitusFormat
    type       T_LevitusFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: GeometryFileName
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: BaseBulletin
        character(len=PathLength)               :: DatesFile
        integer                                 :: NumDates = FillValueInt
        type(T_Time), dimension(:), pointer     :: OutPutDates
        integer                                 :: FieldsNumber
        integer                                 :: ClientNumber
        logical                                 :: FirstInputFiles = .false., FirstMapping = .true.
        real                                    :: SpatialResolution
        real                                    :: FillValue
        real,    dimension(2)                   :: UpperRightCornerXY, LowerLeftCornerXY
        integer, dimension(2)                   :: UpperRightCornerJI, LowerLeftCornerJI
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:  ),       pointer     :: XX, CenterX
        real, dimension(:  ),       pointer     :: YY, CenterY 
        real, dimension(:,:,:),     pointer     :: SZZ
        real, dimension(NLevitusAnnual)         :: LevitusLevels
        real,    dimension(:,:,:),  pointer     :: Aux3D, Aux3DLevitus
        integer, dimension(:,:,:),  pointer     :: WaterPoints3D
        type(T_Size3D)                          :: LevitusSize, LevitusWSize, WorkSize, Size
        type(T_Field),              pointer     :: FirstField         
        type(T_Field),              pointer     :: LastField              
    end type  T_LevitusFormat

    type(T_LevitusFormat), pointer              :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertLevitusFormat(EnterDataID, ClientNumber, STAT)

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

        call WriteLevitusGeometry

        call AllocateVariables

        call ConstructGrid
  
        call Open_HDF5_OutPut_File

        call Construct_FieldList

        call KillLevitusFormat


        STAT = SUCCESS_


    end subroutine ConvertLevitusFormat

    !------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        character(LEN=StringLength)                 :: period
        integer                                     :: STAT_CALL
        integer                                     :: iflag, i
        integer                                     :: lines, columns, layers, WindowLines, WindowColumns

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                             &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR15'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR20'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR25'

        call GetData(Me%GeometryFileName,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GEOMETRY_FILENAME',                         &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR35'


        call GetData(period,                                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PERIODICITY',                                      &
                     default      = 'monthly',                                          &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR60'

        if      (trim(period)=='annual' ) then

            Me%NumDates = 1
            allocate(Me%OutPutDates(Me%NumDates))  

            call SetDate(Me%OutPutDates(1),                                             &
                        Year    = CyclicTime,                                           & 
                        Month   = CyclicTime,                                           & 
                        Day     = CyclicTime,                                           &
                        Hour    = CyclicTime,                                           &
                        Minute  = CyclicTime,                                           &
                        Second  = CyclicTime) 

        else if (trim(period)=='monthly') then

            Me%NumDates = 12
            allocate(Me%OutPutDates(Me%NumDates))  
            do i = 1, 12
                call SetDate(Me%OutPutDates(i),                                         &
                            Year    = CyclicTime,                                       &
                            Month   = real(i),                                          & 
                            Day     = 1.,                                               &
                            Hour    = 0.,                                               &
                            Minute  = 0.,                                               &
                            Second  = 0.) 
            enddo

        else 
            
            write(*,*) 'PERIODICITY : ', trim(period), ' this is not a valid option' 
            stop       'ReadGlobalOptions - ModuleLevitusFormat - ERR70'

        endif


        call GetData(Me%SpatialResolution,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SPATIAL_RESOLUTION',                               &
                     default      = 0.25,                                               &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR80'

        if (Me%SpatialResolution /= 0.25 .and. Me%SpatialResolution /= 1. .and. Me%SpatialResolution /= 5.) then
            stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR90'
        endif

        lines   = int(180./Me%SpatialResolution)
        columns = int(360./Me%SpatialResolution)
        layers  = NLevitusAnnual


        call GetData(Me%LowerLeftCornerXY,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'LOWER_LEFT_CORNER',                                &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR84'

        if (iflag == 0) then

            Me%LowerLeftCornerXY(1) =   0.
            Me%LowerLeftCornerXY(2) = -90.

        else

            if (Me%LowerLeftCornerXY(1) >  360.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR85'
            if (Me%LowerLeftCornerXY(1) < -360.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR85a'
            if (Me%LowerLeftCornerXY(2) >   90.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR85b'
            if (Me%LowerLeftCornerXY(2) <  -90.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR85c'

        endif

        call GetData(Me%UpperRightCornerXY,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'UPPER_RIGHT_CORNER',                               &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR86'

        if (iflag == 0) then

            Me%UpperRightCornerXY(1) =  360.
            Me%UpperRightCornerXY(2) =  90.

        else

            if (Me%UpperRightCornerXY(1) >  360.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR87'
            if (Me%UpperRightCornerXY(1) < -360.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR87a'
            if (Me%UpperRightCornerXY(2) >   90.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR88'
            if (Me%UpperRightCornerXY(2) <  -90.) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR88a'

            if (Me%UpperRightCornerXY(1) < Me%LowerLeftCornerXY(1)) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR89'
            if (Me%UpperRightCornerXY(2) < Me%LowerLeftCornerXY(2)) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR89a'

        endif


!        if (Me%LowerLeftCornerXY(1)  < 0. ) Me%LowerLeftCornerXY(1)  = 360. - Me%LowerLeftCornerXY(1)
!        if (Me%UpperRightCornerXY(1) < 0. ) Me%UpperRightCornerXY(1) = 360. - Me%UpperRightCornerXY(1)

        Me%LowerLeftCornerXY (:) = int(Me%LowerLeftCornerXY (:)/Me%SpatialResolution)*Me%SpatialResolution
        Me%UpperRightCornerXY(:) = int(Me%UpperRightCornerXY(:)/Me%SpatialResolution)*Me%SpatialResolution


        
        WindowLines   = (Me%UpperRightCornerXY(2)- Me%LowerLeftCornerXY(2))/Me%SpatialResolution
        WindowColumns = (Me%UpperRightCornerXY(1)- Me%LowerLeftCornerXY(1))/Me%SpatialResolution

        if (WindowLines   < 1) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR87'
        if (WindowColumns < 1) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR88'


        if (Me%LowerLeftCornerXY(1) >= 0.) then

            Me%LowerLeftCornerJI (1) =   Me%LowerLeftCornerXY (1) / Me%SpatialResolution + 1

        else

            Me%LowerLeftCornerJI (1) =   (Me%LowerLeftCornerXY(1) + 360.)/ Me%SpatialResolution + 1

        endif

        Me%UpperRightCornerJI(1) =   Me%LowerLeftCornerJI (1) + WindowColumns - 1

        if (Me%UpperRightCornerJI(1) > columns) then
                
            Me%UpperRightCornerJI(1) = Me%UpperRightCornerJI(1) - columns

        endif


        Me%LowerLeftCornerJI (2) =   (Me%LowerLeftCornerXY (2) + 90.)/ Me%SpatialResolution + 1

        Me%UpperRightCornerJI(2) =   Me%LowerLeftCornerJI (2) + WindowLines - 1

        
        Me%LevitusSize%ILB = 0; Me%LevitusSize%JLB = 0;  Me%LevitusSize%KLB = 0
        Me%LevitusSize%IUB = lines + 1; Me%LevitusSize%JUB = columns + 1;  Me%LevitusSize%KUB = layers + 1

        Me%LevitusWSize%ILB = 1; Me%LevitusWSize%JLB = 1;  Me%LevitusWSize%KLB = 1
        Me%LevitusWSize%IUB = lines; Me%LevitusWSize%JUB = columns;  Me%LevitusWSize%KUB = layers

        Me%WorkSize%ILB = 1; Me%WorkSize%JLB = 1;  Me%WorkSize%KLB = 1
        Me%WorkSize%IUB =Windowlines; Me%WorkSize%JUB = Windowcolumns;  Me%WorkSize%KUB = layers

        Me%Size%ILB = Me%WorkSize%ILB - 1; Me%Size%JLB = Me%WorkSize%JLB - 1;  Me%Size%KLB = Me%WorkSize%KLB - 1
        Me%Size%IUB = Me%WorkSize%IUB + 1; Me%Size%JUB = Me%WorkSize%JUB + 1;  Me%Size%KUB = Me%WorkSize%KUB + 1


        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleLevitusFormat - ERR90'

        !Levitus Layers
        Me%LevitusLevels(33) = 0.             
        Me%LevitusLevels(32) = 10.            
        Me%LevitusLevels(31) = 20.            
        Me%LevitusLevels(30) = 30.            
        Me%LevitusLevels(29) = 50.            
        Me%LevitusLevels(28) = 75.          
        Me%LevitusLevels(27) = 100.
        Me%LevitusLevels(26) = 125.
        Me%LevitusLevels(25) = 150.
        Me%LevitusLevels(24) = 200.
        Me%LevitusLevels(23) = 250.
        Me%LevitusLevels(22) = 300.
        Me%LevitusLevels(21) = 400.
        Me%LevitusLevels(20) = 500.
        Me%LevitusLevels(19) = 600.
        Me%LevitusLevels(18) = 700.
        Me%LevitusLevels(17) = 800.
        Me%LevitusLevels(16) = 900.
        Me%LevitusLevels(15) = 1000.
        Me%LevitusLevels(14) = 1100.
        Me%LevitusLevels(13) = 1200.
        Me%LevitusLevels(12) = 1300.
        Me%LevitusLevels(11) = 1400.
        Me%LevitusLevels(10) = 1500.
        Me%LevitusLevels( 9) = 1750.
        Me%LevitusLevels( 8) = 2000.
        Me%LevitusLevels( 7) = 2500.
        Me%LevitusLevels( 6) = 3000.
        Me%LevitusLevels( 5) = 3500.
        Me%LevitusLevels( 4) = 4000.
        Me%LevitusLevels( 3) = 4500.
        Me%LevitusLevels( 2) = 5000.
        Me%LevitusLevels( 1) = 5500.

!0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700,800,900,1000,
! 1100,1200,1300,1400,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500


    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer             :: i, j
        !----------------------------------------------------------------------

        allocate(Me%Bathymetry    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Me%CenterX       (Me%Size%JLB:Me%Size%JUB))
        allocate(Me%XX            (Me%Size%JLB:Me%Size%JUB))
        allocate(Me%CenterY       (Me%Size%ILB:Me%Size%IUB))
        allocate(Me%YY            (Me%Size%ILB:Me%Size%IUB))
        allocate(Me%WaterPoints3D (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))

        Me%WaterPoints3D(:,:,:) = 1

        allocate(Me%SZZ           (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
        allocate(Me%Aux3D         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
        allocate(Me%Aux3DLevitus  (Me%LevitusSize%ILB:Me%LevitusSize%IUB,   &
                                   Me%LevitusSize%JLB:Me%LevitusSize%JUB,   &
                                   Me%LevitusSize%KLB:Me%LevitusSize%KUB))

        Me%Bathymetry (:,:) = MaxDepth


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            Me%SZZ(i,j,Me%WorkSize%KLB:Me%WorkSize%KUB) = Me%LevitusLevels(1:NLevitusAnnual)

            Me%SZZ(i,j,Me%WorkSize%KLB-1              ) = MaxDepth

        enddo
        enddo

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

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
!                        stop 'Construct_FieldList - ModuleLevitusFormat - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                stop 'Construct_FieldList - ModuleLevitusFormat - ERR02'
            else cd1
                stop 'Construct_FieldList - ModuleLevitusFormat - ERR03'
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

        nullify(NewField%MonthlyValues3D)

        nullify(NewField%AnnualValues3D)

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
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: line, l, NFiles, STAT_CALL, iflag

        !----------------------------------------------------------------------

        call GetData(NewField%AnnualFile,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'ANNUAL_FILE',                                      &
                     ClientModule = 'LevitusFormat',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleLevitusFormat - ERR10'
        if (iflag     == 0)        stop 'Construct_InputFiles - ModuleLevitusFormat - ERR20'


        allocate (NewField%AnnualValues3D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))

        NewField%AnnualValues3D(:,:,:) = FillValueReal

        call ReadLevitusAnnualFile    (NewField)

        !call MappingAnnualWaterPoints (NewField)


ID:     if (Me%NumDates == 12) then

        call ExtractBlockFromBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,           &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockInBlockFound= BlockFound,                &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_)then


BF:         if (BlockFound) then

                NFiles =  LastLine - FirstLine - 1

                !This module is only prepared to read annual or monthly climatologic values
                if (NFiles /= 12) stop 'Construct_InputFiles - ModuleLevitusFormat - ERR30'

                
                !This error happens when a different of instants is defined in one of the properties
                if (Me%FirstInputFiles .and. Me%NumDates /= NFiles) stop 'Construct_InputFiles - ModuleLevitusFormat - ERR40'

                !Allocates auxiliar variables
                allocate (NewField%InputFiles(NFiles))

                allocate (NewField%MonthlyValues3D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))
            
                l = 1
                do line = FirstLine + 1, LastLine - 1

                    call GetData(NewField%InputFiles(l), EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleLevitusFormat - ERR50'

                    NewField%MonthlyValues3D(:,:,:) = NewField%AnnualValues3D(:,:,:)

                    call ReadLevitusMonthlyFile    (NewField, l)

                    call OutputFields              (NewField, l)

                    l = l + 1

                enddo

                call ConvertToMohidUnits    (NewField)

                Me%FirstInputFiles = .true.

            else BF

                stop 'Construct_InputFiles - ModuleLevitusFormat - ERR60'

            end if BF

        else   IS

            stop 'Construct_InputFiles - ModuleLevitusFormat - ERR70'

        end if IS

        else ID

            call OutputFields              (NewField, 1)

        endif ID

    end subroutine Construct_InputFiles

    !--------------------------------------------------------------------------

    subroutine ReadLevitusMonthlyFile (NewField, l)

        !Arguemtns-------------------------------------------------------------
        type(T_Field), pointer                      :: NewField
        integer                                     :: l

        !Local----------------------------------------------------------------
        real                                        :: Aux1, Aux2
        integer                                     :: i, j, k, STAT_CALL, Jupper, iw, jw
        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading Levitus output file...',trim(NewField%InputFiles(l))
        


        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevitusMonthlyFile - ModuleLevitusFormat - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = NewField%InputFiles(l),                                           &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevitusMonthlyFile - ModuleLevitusFormat - ERR20'
        

        Me%Aux3DLevitus(:,:,:) = FillValueReal
        
        do k = NLevitusAnnual, NLevitusMonthly, - 1
            read(Me%Unit,'(10(f8.0))') ((Me%Aux3DLevitus(i,j,k), &
                j=Me%LevitusWSize%JLB,Me%LevitusWSize%JUB),i=Me%LevitusWSize%ILB,Me%LevitusWSize%IUB)
        enddo

        if (Me%UpperRightCornerJI(1) > Me%LowerLeftCornerJI(1)) then
            Jupper = Me%UpperRightCornerJI(1)
        else
            Jupper = int(360 / Me%SpatialResolution) 
        endif

        do k = NLevitusMonthly+1, NLevitusAnnual
        do j = Me%LowerLeftCornerJI(1), Jupper
        do i = Me%LowerLeftCornerJI(2), Me%UpperRightCornerJI(2)

            Aux1 = 0.
            Aux2 = 0.

            jw = j - Me%LowerLeftCornerJI(1) + 1
            iw = i - Me%LowerLeftCornerJI(2) + 1

            if (Me%Aux3DLevitus(i, j, k  ) > Me%FillValue / 2.) Aux1 = 1.
            if (Me%Aux3DLevitus(i, j, k-1) > Me%FillValue / 2.) Aux2 = 1.
            if ((Aux1 + Aux2) > 0.) then
                NewField%MonthlyValues3D(iw, jw, k) = (Me%Aux3DLevitus(i, j, k) * Aux1 + &
                                                       Me%Aux3DLevitus(i, j, k-1) * Aux2) / (Aux1 + Aux2)
            else
                NewField%MonthlyValues3D(iw, jw, k) = FillValueReal
            endif

        enddo
        enddo

        if (jUpper /= Me%UpperRightCornerJI(1)) then
            do j = 1, Me%UpperRightCornerJI(1)
            do i = Me%LowerLeftCornerJI(2), Me%UpperRightCornerJI(2)

                Aux1 = 0.
                Aux2 = 0.

                jw = j + Jupper - Me%LowerLeftCornerJI(1) + 1
                iw = i - Me%LowerLeftCornerJI(2) + 1

                if (Me%Aux3DLevitus(i, j, k  ) > Me%FillValue / 2.) Aux1 = 1.
                if (Me%Aux3DLevitus(i, j, k-1) > Me%FillValue / 2.) Aux2 = 1.
                if ((Aux1 + Aux2) > 0.) then
                    NewField%MonthlyValues3D(iw, jw, k) = (Me%Aux3DLevitus(i, j, k  ) * Aux1 + &
                                                           Me%Aux3DLevitus(i, j, k-1) * Aux2) / (Aux1 + Aux2)
                else
                    NewField%MonthlyValues3D(iw, jw, k) = FillValueReal
                endif

            enddo
            enddo

        endif

        enddo



        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevitusMonthlyFile - ModuleLevitusFormat - ERR30'

    end subroutine ReadLevitusMonthlyFile
    

    !----------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine ReadLevitusAnnualFile (NewField)

        !Arguemtns-------------------------------------------------------------
        type(T_Field), pointer                      :: NewField

        !Local----------------------------------------------------------------
        real                                        :: Aux1, Aux2
        integer                                     :: i, j, k, STAT_CALL, jupper, iw, jw
        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading Levitus output file...',trim(NewField%AnnualFile)
        

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevitusAnnualFile - ModuleLevitusFormat - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = NewField%AnnualFile,                                              &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevitusAnnualFile - ModuleLevitusFormat - ERR20'
        

        Me%Aux3DLevitus(:,:,:) = FillValueReal
        
        do k = NLevitusAnnual, 1, -1
            read(Me%Unit,'(10(f8.0))') ((Me%Aux3DLevitus(i,j,k),              &
                j=Me%LevitusWSize%JLB,Me%LevitusWSize%JUB),i=Me%LevitusWSize%ILB,Me%LevitusWSize%IUB)
        enddo


        if (Me%UpperRightCornerJI(1) > Me%LowerLeftCornerJI(1)) then
            Jupper = Me%UpperRightCornerJI(1)
        else
            Jupper = int(360 / Me%SpatialResolution) 
        endif


        do k = 1, NLevitusAnnual
        do j = Me%LowerLeftCornerJI(1), Jupper
        do i = Me%LowerLeftCornerJI(2), Me%UpperRightCornerJI(2)

            Aux1 = 0.
            Aux2 = 0.

            jw = j - Me%LowerLeftCornerJI(1) + 1
            iw = i - Me%LowerLeftCornerJI(2) + 1

            if (Me%Aux3DLevitus(i, j, k  ) > Me%FillValue / 2.) Aux1 = 1.
            if (Me%Aux3DLevitus(i, j, k-1) > Me%FillValue / 2.) Aux2 = 1.
            if ((Aux1 + Aux2) > 0.) then
                NewField%AnnualValues3D(iw, jw, k) = (Me%Aux3DLevitus(i, j, k  ) * Aux1 + &
                                                      Me%Aux3DLevitus(i, j, k-1) * Aux2) / (Aux1 + Aux2)
            else
                NewField%AnnualValues3D(iw, jw, k) = FillValueReal
            endif
        enddo
        enddo

        if (Jupper /= Me%UpperRightCornerJI(1)) then

            do j = 1, Me%UpperRightCornerJI(1)
            do i = Me%LowerLeftCornerJI(2), Me%UpperRightCornerJI(2)
                Aux1 = 0.
                Aux2 = 0.

                jw = j + Jupper - Me%LowerLeftCornerJI(1) + 1
                iw = i - Me%LowerLeftCornerJI(2) + 1

                if (Me%Aux3DLevitus(i, j, k  ) > Me%FillValue / 2.) Aux1 = 1.
                if (Me%Aux3DLevitus(i, j, k-1) > Me%FillValue / 2.) Aux2 = 1.
                if ((Aux1 + Aux2) > 0.) then
                    NewField%AnnualValues3D(iw, jw, k) = (Me%Aux3DLevitus(i, j, k  ) * Aux1 + &
                                                          Me%Aux3DLevitus(i, j, k-1) * Aux2) / (Aux1 + Aux2)
                else
                    NewField%AnnualValues3D(iw, jw, k) = FillValueReal
                endif

            enddo
            enddo

        endif

        enddo


        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLevitusAnnualFile - ModuleLevitusFormat - ERR30'

    end subroutine ReadLevitusAnnualFile
    


    !--------------------------------------------------------------------------

    subroutine MappingAnnualWaterPoints(NewField)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer       :: NewField

        !Local-----------------------------------------------------------------
        integer                      :: i, j, k, WaterPoint
        !----------------------------------------------------------------------

        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (NewField%AnnualValues3D(i,j,k) < Me%FillValue / 2.) then

                WaterPoint = 0

            else

                WaterPoint = 1

            endif

            if (Me%FirstMapping) then

                Me%WaterPoints3D(i, j, k) = WaterPoint

            else

                if (Me%WaterPoints3D(i, j, k) /= WaterPoint) then

                    write (*,*) 'Mapping error ', trim(NewField%ID%Name),' - Annual'

                    stop 'MappingAnnualWaterPoints - ModuleLevitusFormat - ERR10'

                endif

            endif
        enddo
        enddo
        enddo

        if (Me%FirstMapping) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                Me%SZZ(i,j,Me%WorkSize%KLB:Me%WorkSize%KUB) = Me%LevitusLevels(1:NLevitusAnnual)

                Me%SZZ(i,j,Me%WorkSize%KLB-1              ) = MaxDepth

                if (Me%WaterPoints3D(i, j, Me%WorkSize%KUB) == 0) then 
                    Me%Bathymetry(i, j) = -99.
                else

                    do k = Me%WorkSize%KUB-1, Me%WorkSize%KLB, -1

                        if (Me%WaterPoints3D(i, j, k) == 0) then
                            Me%Bathymetry(i, j) = Me%LevitusLevels(k)
                            exit
                        endif
                    enddo

                endif

            enddo
            enddo

            call ConstructGrid

       
            call Open_HDF5_OutPut_File

            Me%FirstMapping = .false.
        endif


    end subroutine MappingAnnualWaterPoints

    !--------------------------------------------------------------------------


    
    !------------------------------------------------------------------------

    
    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        type (T_Size2D)                             :: WorkSize2D
        real                                        :: Xorig, Yorig, Latitude, Longitude
        logical                                     :: ContinuesCompute
        integer                                     :: STAT_CALL, SIMPLE_GEOG, i, j
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'

        ContinuesCompute = .FALSE.

        WorkSize2D%ILB = Me%WorkSize%ILB
        WorkSize2D%IUB = Me%WorkSize%IUB
        WorkSize2D%JLB = Me%WorkSize%JLB
        WorkSize2D%JUB = Me%WorkSize%JUB

        Me%XX(1) = 0.
        Me%YY(1) = 0.
        do j = Me%WorkSize%JLB + 1, Me%WorkSize%JUB + 1
            Me%XX(j) = Me%XX(j-1) + Me%SpatialResolution
        enddo

        do i = Me%WorkSize%ILB + 1, Me%WorkSize%IUB + 1
            Me%YY(i) = Me%YY(i-1) + Me%SpatialResolution
        enddo

        Xorig =  Me%LowerLeftCornerXY(1)
        Yorig =  Me%LowerLeftCornerXY(2)

        call GetCoordTypeList(SIMPLE_GEOG = SIMPLE_GEOG)

        call WriteGridData (FileName        = Me%GridFileName,              &
                            XX              = Me%XX,                        &
                            YY              = Me%YY,                        &
                            COMENT1         = 'Levitus Grid based on file :',  &
                            COMENT2         = trim(Me%FileName),            &
                            WorkSize        = WorkSize2D,                   & 
                            CoordType       = SIMPLE_GEOG,                  &
                            Xorig           = Xorig,                        &
                            Yorig           = Yorig,                        &
                            Zone            = 0,                            &
                            Grid_Angle      = 0.,                           &
                            Latitude        = Latitude,                     &
                            Longitude       = Longitude,                    &
                            FillValue       = -99.,                         &
                            Overwrite       = ON,                           &
                            GridData2D_Real = Me%Bathymetry,                &
                            STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleLevitusFormat - ERR01'


        call StartComputeTime(Me%ObjTime, Me%OutPutDates(1), Me%OutPutDates(1), DT = 0.0,    &
                                 VariableDT = .false., STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleLevitusFormat - ERR02a'

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleLevitusFormat - ERR02b'


    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    
    subroutine WriteLevitusGeometry
        
        !Local-----------------------------------------------------------------
        integer                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Writting Geometry...'

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteLevitusGeometry - ModuleLevitusFormat - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = Me%GeometryFilename,                                              &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'UNKNOWN',                                                        &
             Action = 'WRITE',                                                          &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteLevitusGeometry - ModuleLevitusFormat - ERR20'

        write(Me%Unit,*) '<begindomain>'
        write(Me%Unit,*) 'ID : 1'
        write(Me%Unit,*) 'TYPE : CARTESIAN'
        write(Me%Unit,*) 'DOMAINDEPTH : 0.'
        write(Me%Unit,*) 'LAYERS : 33'
        write(Me%Unit,'(A138)') 'LAYERTHICKNESS : 6500 500 500 500 500 500 500 500 250 250 100 100 100 100 100 100 100 100 100 100 100 100 50 50 50 25 25 25 25 20 10 10 10'
        write(Me%Unit,*) 'MININITIALLAYERTHICKNESS : 1'
        write(Me%Unit,*) '<enddomain>'
        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteLevitusGeometry - ModuleLevitusFormat - ERR30'


    end subroutine WriteLevitusGeometry

    
    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits(NewField)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        
        !Begin-----------------------------------------------------------------

        select case(NewField%ID%IDNumber)

            case(Temperature_)

                NewField%ID%Units     = 'ºC'

            case(Salinity_)

                NewField%ID%Units     = 'psu'


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
    subroutine OutputFields(Field, i)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                          :: Field
        integer                                         :: i
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL
        !Begin-----------------------------------------------------------------
        

        if (Field%FirstField) then

    !           Dados para escriver uma soa vez cada date:
            call ExtractDate   (Me%OutPutDates(i),                                      &
                                AuxTime(1), AuxTime(2), AuxTime(3),                     &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleLevitusFormat - ERR10'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleLevitusFormat - ERR20'

        endif

        !Writes SZZ
        !call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
        !                     Me%WorkSize%JUB, LevitusKLB, LevitusKUB, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR30'

        !call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",             &
        !                    "m", Array3D = Me%SZZ, OutputNumber = i,                &
        !                    STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR40'


        !call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",          &
        !                    "m", Array3D = Me%WaterPoints3D, OutputNumber = i,      &
        !                    STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR50'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

   
        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,                            &
                           Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleLevitusFormat - ERR60'

        if (Me%NumDates == 1) then
            Me%Aux3D(:,:,:) = Field%AnnualValues3D (:,:,:)
        else
            Me%Aux3D(:,:,:) = Field%MonthlyValues3D(:,:,:)
        endif

        call HDF5WriteData(Me%ObjHDF5,                                      &
                           "/Results/"//Field%ID%Name,                      &
                           Field%ID%Name,                                   &
                           Field%ID%Units,                                  &
                           Array3D      = Me%Aux3D,                         &
                           OutputNumber = i,                                &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleLevitusFormat - ERR70'



        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleLevitusFormat - ERR80'


    end subroutine OutputFields


    !----------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR02'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR04'            
   

        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = Me%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR08'

        !Writes SZZ - Constant
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                             Me%WorkSize%JUB, Me%WorkSize%KLB-1, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR30'

        call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",             &
                            "m", Array3D = Me%SZZ,                                  &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR40'


        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleLevitusFormat - ERR09'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    logical function CheckName(LevitusName, MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: LevitusName
        character(Len=StringLength)     :: MohidName
        
        !Begin-----------------------------------------------------------------


        select case(trim(LevitusName))

            case('u')

                MohidName = GetPropertyName(VelocityU_)
                CheckName = .true.

            case('v')

                MohidName = GetPropertyName(VelocityV_)
                CheckName = .true.

            case('temperature')

                MohidName = GetPropertyName(Temperature_)
                CheckName = .true.

            case('salinity')

                MohidName = GetPropertyName(Salinity_)
                CheckName = .true.

            case default
                
                CheckName = .false.

        end select


    end function CheckName
    
    
    !--------------------------------------------------------------------------

    
    subroutine KillLevitusFormat
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer      :: Field
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillLevitusFormat - ModuleLevitusFormat - ERR01d'
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillLevitusFormat - ModuleLevitusFormat - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillLevitusFormat - ModuleLevitusFormat - ERR03'


        deallocate(Me%Bathymetry)
        deallocate(Me%CenterX)
        deallocate(Me%XX)
        deallocate(Me%CenterY)
        deallocate(Me%YY)
        deallocate(Me%WaterPoints3D)
        deallocate(Me%SZZ)
        deallocate(Me%Aux3D)
        deallocate(Me%Aux3DLevitus)

        Field => Me%FirstField 

        do while (associated(Field)) 

            deallocate(Field%AnnualValues3D  )

            if (Me%NumDates == 12) then
                deallocate(Field%MonthlyValues3D )
                deallocate(Field%InputFiles)
            endif

            Field => Field%Next

        enddo

        deallocate(Me%FirstField)
        deallocate(Me%LastField )
        nullify   (Me%FirstField)
        nullify   (Me%LastField )

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillLevitusFormat

    !--------------------------------------------------------------------------
 
end module ModuleLevitusFormat
