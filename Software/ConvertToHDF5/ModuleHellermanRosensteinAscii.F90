!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : HellermanRosensteinAscii
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert HellermanRosensteinAscii files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleHellermanRosensteinAscii

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
    public  :: ConvertHellermanRosensteinAscii
    private ::      ReadGlobalOptions
    private ::      ReadHRMonthlyFile
    private ::          AddField
    private ::      OutputFields
    private ::          Open_HDF5_OutPut_File
    private ::      KillHellermanRosensteinAscii

    !Parameters----------------------------------------------------------------
    integer, parameter                          :: NHellermanRosensteinAnnual  = 33
    integer, parameter                          :: NHellermanRosensteinMonthly = 10

    real,    parameter                          :: MaxDepth = 0.

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
        real, dimension(:,:    ),   pointer     :: AnnualValues2D
        real, dimension(:,:    ),   pointer     :: MonthlyValues2D
        character(len=StringLength)             :: File
        logical                                 :: FirstField
        type(T_Field),              pointer     :: Next, Prev
    end type  T_Field

    
    private :: T_HellermanRosensteinAscii
    type       T_HellermanRosensteinAscii
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
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
        real, dimension(NHellermanRosensteinAnnual)         :: HellermanRosensteinLevels
        real,    dimension(:,:  ),  pointer     :: Aux2D, Aux2DHellermanRosenstein
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D
        type(T_Size2D)                          :: HellermanRosensteinSize, HellermanRosensteinWSize, WorkSize, Size
        type(T_Field),              pointer     :: FirstField         
        type(T_Field),              pointer     :: LastField              
    end type  T_HellermanRosensteinAscii

    type(T_HellermanRosensteinAscii), pointer              :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertHellermanRosensteinAscii(EnterDataID, ClientNumber, STAT)

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

        call AllocateVariables

        call ConstructGrid
  
        call Open_HDF5_OutPut_File

        call Construct_FieldList

        call KillHellermanRosensteinAscii


        STAT = SUCCESS_


    end subroutine ConvertHellermanRosensteinAscii

    !------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: AuxTime
        character(LEN=StringLength)                 :: period
        integer                                     :: STAT_CALL
        integer                                     :: iflag, i
        integer                                     :: lines, columns, layers, WindowLines, WindowColumns

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                             &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR15'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR20'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR25'


        call GetData(period,                                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PERIODICITY',                                      &
                     default      = 'monthly',                                          &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR60'

        if      (trim(period)=='annual' ) then

            write(*,*) 'There is only monthly periodicity' 
            stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR65'

            Me%NumDates = 1
            allocate(Me%OutPutDates(Me%NumDates))  

            call SetDate(AuxTime,                                                       &
                        Year    = CyclicTime,                                           & 
                        Month   = 1.,                                                   & 
                        Day     = 1.,                                                   &
                        Hour    = 0.,                                                   &
                        Minute  = 0.,                                                   &
                        Second  = 0.) 

            Me%OutPutDates(1) = AuxTime

        else if (trim(period)=='monthly') then

            Me%NumDates = 12
            allocate(Me%OutPutDates(Me%NumDates))  
            do i = 1, 12
                call SetDate(AuxTime,                                                   &
                            Year    = CyclicTime,                                       &
                            Month   = real(i),                                          & 
                            Day     = 1.,                                               &
                            Hour    = 0.,                                               &
                            Minute  = 0.,                                               &
                            Second  = 0.) 
                Me%OutPutDates(i) = AuxTime
            enddo

        endif


        call GetData(Me%SpatialResolution,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SPATIAL_RESOLUTION',                               &
                     default      = 2.,                                                 &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR80'

        if (Me%SpatialResolution /= 2.) then
            stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR90'
        endif

        lines   = int(180./Me%SpatialResolution)
        columns = int(360./Me%SpatialResolution)
        layers  = NHellermanRosensteinAnnual


        call GetData(Me%LowerLeftCornerXY,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'LOWER_LEFT_CORNER',                                &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR84'

        if (iflag == 0) then

            Me%LowerLeftCornerXY(1) =   0.
            Me%LowerLeftCornerXY(2) = -90.

        else

            if (Me%LowerLeftCornerXY(1) >  360.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR85'
            if (Me%LowerLeftCornerXY(1) < -360.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR85a'
            if (Me%LowerLeftCornerXY(2) >   90.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR85b'
            if (Me%LowerLeftCornerXY(2) <  -90.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR85c'

        endif

        call GetData(Me%UpperRightCornerXY,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'UPPER_RIGHT_CORNER',                               &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR86'

        if (iflag == 0) then

            Me%UpperRightCornerXY(1) =  360.
            Me%UpperRightCornerXY(2) =  90.

        else

            if (Me%UpperRightCornerXY(1) >  360.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR87'
            if (Me%UpperRightCornerXY(1) < -360.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR87a'
            if (Me%UpperRightCornerXY(2) >   90.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR88'
            if (Me%UpperRightCornerXY(2) <  -90.) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR88a'

            if (Me%UpperRightCornerXY(1) < Me%LowerLeftCornerXY(1)) &
                stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR89'
            if (Me%UpperRightCornerXY(2) < Me%LowerLeftCornerXY(2)) &
                stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR89a'

        endif


!        if (Me%LowerLeftCornerXY(1)  < 0. ) Me%LowerLeftCornerXY(1)  = 360. - Me%LowerLeftCornerXY(1)
!        if (Me%UpperRightCornerXY(1) < 0. ) Me%UpperRightCornerXY(1) = 360. - Me%UpperRightCornerXY(1)

        Me%LowerLeftCornerXY (:) = int(Me%LowerLeftCornerXY (:)/Me%SpatialResolution)*Me%SpatialResolution
        Me%UpperRightCornerXY(:) = int(Me%UpperRightCornerXY(:)/Me%SpatialResolution)*Me%SpatialResolution


        
        WindowLines   = (Me%UpperRightCornerXY(2)- Me%LowerLeftCornerXY(2))/Me%SpatialResolution
        WindowColumns = (Me%UpperRightCornerXY(1)- Me%LowerLeftCornerXY(1))/Me%SpatialResolution

        if (WindowLines   < 1) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR87'
        if (WindowColumns < 1) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR88'


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

        
        Me%HellermanRosensteinSize%ILB = 0; Me%HellermanRosensteinSize%JLB = 0;
        Me%HellermanRosensteinSize%IUB = lines + 1; Me%HellermanRosensteinSize%JUB = columns + 1

        Me%HellermanRosensteinWSize%ILB = 1; Me%HellermanRosensteinWSize%JLB = 1
        Me%HellermanRosensteinWSize%IUB = lines; Me%HellermanRosensteinWSize%JUB = columns

        Me%WorkSize%ILB = 1; Me%WorkSize%JLB = 1
        Me%WorkSize%IUB =Windowlines; Me%WorkSize%JUB = Windowcolumns

        Me%Size%ILB = Me%WorkSize%ILB - 1; Me%Size%JLB = Me%WorkSize%JLB - 1
        Me%Size%IUB = Me%WorkSize%IUB + 1; Me%Size%JUB = Me%WorkSize%JUB + 1


        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleHellermanRosensteinAscii - ERR90'



    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        allocate(Me%Bathymetry    (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(Me%CenterX       (Me%Size%JLB:Me%Size%JUB))
        allocate(Me%XX            (Me%Size%JLB:Me%Size%JUB))
        allocate(Me%CenterY       (Me%Size%ILB:Me%Size%IUB))
        allocate(Me%YY            (Me%Size%ILB:Me%Size%IUB))
        allocate(Me%WaterPoints2D (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

        Me%WaterPoints2D(:,:) = 1

        allocate(Me%Aux2D         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%Aux2DHellermanRosenstein  (Me%HellermanRosensteinSize%ILB:Me%HellermanRosensteinSize%IUB,&
                                               Me%HellermanRosensteinSize%JLB:Me%HellermanRosensteinSize%JUB))

        Me%Bathymetry (:,:) = MaxDepth

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
!                        stop 'Construct_FieldList - ModuleHellermanRosensteinAscii - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                stop 'Construct_FieldList - ModuleHellermanRosensteinAscii - ERR02'
            else cd1
                stop 'Construct_FieldList - ModuleHellermanRosensteinAscii - ERR03'
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

        nullify(NewField%MonthlyValues2D)

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
        integer                                     :: l, STAT_CALL, iflag

        !----------------------------------------------------------------------

        call GetData(NewField%File,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlockInBlock,                                   &
                     keyword      = 'FILE',                                             &
                     ClientModule = 'HellermanRosensteinAscii',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleHellermanRosensteinAscii - ERR10'
        if (iflag     == 0)        stop 'Construct_InputFiles - ModuleHellermanRosensteinAscii - ERR20'

        allocate (NewField%MonthlyValues2D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading HellermanRosenstein ascii file...',trim(NewField%File)
        

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleHellermanRosensteinAscii - ERR30'

        open(Unit   = Me%Unit,                                                          &
             File   = NewField%File,                                                    &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleHellermanRosensteinAscii - ERR40'

    
        do l = 1, 12

            if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleHellermanRosensteinAscii - ERR50'

            call ReadHRMonthlyFile      (NewField)

            call ConvertToMohidUnits    (NewField)

            call OutputFields           (NewField, l)

        enddo


        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Construct_InputFiles - ModuleHellermanRosensteinAscii - ERR60'


    end subroutine Construct_InputFiles


    !--------------------------------------------------------------------------

    subroutine ReadHRMonthlyFile (NewField)

        !Arguemtns-------------------------------------------------------------
        type(T_Field), pointer                      :: NewField

        !Local----------------------------------------------------------------
        real                                        :: Aux1, Aux2
        integer                                     :: i, j, jupper, iw, jw
        !Begin-----------------------------------------------------------------
        

        Me%Aux2DHellermanRosenstein(:,:) = FillValueReal
        !The data is in 10 * Pa
        read(Me%Unit,*) ((Me%Aux2DHellermanRosenstein(i,j),              &
            j=Me%HellermanRosensteinWSize%JLB,Me%HellermanRosensteinWSize%JUB), &
            i=Me%HellermanRosensteinWSize%ILB,Me%HellermanRosensteinWSize%IUB)

        !Transform the data is Pascals 
        Me%Aux2DHellermanRosenstein(:,:) = 0.1 * Me%Aux2DHellermanRosenstein(:,:)


        if (Me%UpperRightCornerJI(1) > Me%LowerLeftCornerJI(1)) then
            Jupper = Me%UpperRightCornerJI(1)
        else
            Jupper = int(360 / Me%SpatialResolution) 
        endif


        do j = Me%LowerLeftCornerJI(1), Jupper
        do i = Me%LowerLeftCornerJI(2), Me%UpperRightCornerJI(2)

            Aux1 = 0.
            Aux2 = 0.

            jw = j - Me%LowerLeftCornerJI(1) + 1
            iw = i - Me%LowerLeftCornerJI(2) + 1

            NewField%MonthlyValues2D(iw, jw) = Me%Aux2DHellermanRosenstein(i, j)
        enddo
        enddo

        if (Jupper /= Me%UpperRightCornerJI(1)) then

            do j = 1, Me%UpperRightCornerJI(1)
            do i = Me%LowerLeftCornerJI(2), Me%UpperRightCornerJI(2)

                jw = j + Jupper - Me%LowerLeftCornerJI(1) + 1
                iw = i - Me%LowerLeftCornerJI(2) + 1

                NewField%MonthlyValues2D(iw, jw) = Me%Aux2DHellermanRosenstein(i, j)

            enddo
            enddo

        endif

    end subroutine ReadHRMonthlyFile
    

    !----------------------------------------------------------------------------



    
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
                            COMENT1         = 'HellermanRosenstein Grid based on file :',  &
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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHellermanRosensteinAscii - ERR01'


        call StartComputeTime(Me%ObjTime, Me%OutPutDates(1), Me%OutPutDates(1), Me%OutPutDates(1), DT = 0.0,    &
                                 VariableDT = .false., STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHellermanRosensteinAscii - ERR02a'

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleHellermanRosensteinAscii - ERR02b'


    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------

    
    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits(NewField)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        
        !Begin-----------------------------------------------------------------

        select case(NewField%ID%IDNumber)

            case(WindStressX_)

                NewField%ID%Units     = 'Pa'

            case(WindStressY_)

                NewField%ID%Units     = 'Pa'


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
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHellermanRosensteinAscii - ERR10'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHellermanRosensteinAscii - ERR20'

        endif

   
        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,                            &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHellermanRosensteinAscii - ERR60'

        if (Me%NumDates == 1) then
            !Me%Aux2D(:,:) = Field%AnnualValues2D (:,:)
        else
            Me%Aux2D(:,:) = Field%MonthlyValues2D(:,:)
        endif

        call HDF5WriteData(Me%ObjHDF5,                                      &
                           "/Results/"//Field%ID%Name,                      &
                           Field%ID%Name,                                   &
                           Field%ID%Units,                                  &
                           Array2D      = Me%Aux2D,                         &
                           OutputNumber = i,                                &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHellermanRosensteinAscii - ERR70'



        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleHellermanRosensteinAscii - ERR80'


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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR02'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR04'            
   

        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",    &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR08'


        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHellermanRosensteinAscii - ERR09'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

   
    !--------------------------------------------------------------------------

    
    subroutine KillHellermanRosensteinAscii
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer      :: Field
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------


        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHellermanRosensteinAscii - ModuleHellermanRosensteinAscii - ERR01d'
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHellermanRosensteinAscii - ModuleHellermanRosensteinAscii - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillHellermanRosensteinAscii - ModuleHellermanRosensteinAscii - ERR03'


        deallocate(Me%Bathymetry)
        deallocate(Me%CenterX)
        deallocate(Me%XX)
        deallocate(Me%CenterY)
        deallocate(Me%YY)
        deallocate(Me%WaterPoints2D)
        deallocate(Me%Aux2D)
        deallocate(Me%Aux2DHellermanRosenstein)

        Field => Me%FirstField 

        do while (associated(Field)) 

            !deallocate(Field%AnnualValues2D  )

            if (Me%NumDates == 12) then
                deallocate(Field%MonthlyValues2D )
            endif

            Field => Field%Next

        enddo

        deallocate(Me%FirstField)
        deallocate(Me%LastField )
        nullify   (Me%FirstField)
        nullify   (Me%LastField )

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillHellermanRosensteinAscii

    !--------------------------------------------------------------------------
 
end module ModuleHellermanRosensteinAscii
