program HDF5Operator

    !Modules needed in the project: 
    !ModuleGlobalData, ModuleTime, ModuleEnterData, ModuleFunctions, ModuleHDF5

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData

    implicit none

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: Differences          = 1
    integer, parameter                  :: Sum                  = 2

    !Variables-----------------------------------------------------------------
    integer                             :: HDF5_READ, HDF5_READWRITE, HDF5_CREATE
    integer                             :: Operation

    !Options-------------------------------------------------------------------
    logical                             :: DataSetWithNumber    = .true.
    logical                             :: CopyOpenPoints       = .true.
    logical                             :: CopyTime             = .true.
    logical                             :: RelativeDifferences  = .true.
    logical                             :: IgnoreInstants       = .false.

    logical                             :: PropertiesInFacesU   = .false.
    logical                             :: PropertiesInFacesV   = .false.
    logical                             :: Properties2D         = .false. 
    
    !Input data file-----------------------------------------------------------
    character(len=PathLength)           :: DataFile             = "HDF5Operator.dat"
    
    !HDF5 groups/datasets paths
    character(len=StringLength)         :: OutputGroup
    character(len=StringLength)         :: OutputDataSet
    character(len=StringLength)         :: InstantsPath         = "/Time"
    character(len=StringLength)         :: MappingGroup
    character(len=StringLength)         :: MappingDataSet
    
    !Mapping-------------------------------------------------------------------
    integer, dimension(:,:,:), pointer  :: Mapping3D

    !Types---------------------------------------------------------------------
    type     T_HDF5File
        integer                         :: ObjHDF5      = 0
        character(len=PathLength)       :: FileName
        type(T_Size3D)                  :: Size
        integer                         :: nInstants
    end type T_HDF5File

    type     T_Property
        character(len=PathLength)       :: Name
        character(len=PathLength)       :: Path
        character(len=PathLength)       :: OutputName
        real                            :: Factor
        integer                         :: TypeZUV
        character(StringLength)         :: Units
        integer                         :: Rank
        type(T_Property), pointer       :: Next
    end type T_Property

    !Property------------------------------------------------------------------
    type(T_Property), pointer           :: FirstProperty, Property
    
    !HDF5 files----------------------------------------------------------------
    type(T_HDF5File)                    :: FirstHDF5File
    type(T_HDF5File)                    :: SecondHDF5File
    type(T_HDF5File)                    :: OutputHDF5File
     
    !--------------------------------------------------------------------------

    call StartUpMohid("HDF5Operator")

    call ReadOptions
        
    call OpenHDF5File(FirstHDF5File)

    call CheckPropertiesInFile(FirstHDF5File)

    if(Operation .ne. Sum)then 

        call OpenHDF5File(SecondHDF5File)

        call CheckFiles

    end if

    call AllocateVariables

    call OpenOutputFile

    call DoOperations

    call Finish


    contains
    
    !--------------------------------------------------------------------------
    
    subroutine ReadOptions
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ClientNumber
        integer                                     :: ObjEnterData = 0
        logical                                     :: BlockFound
        type(T_Property), pointer                   :: NewProperty
        character(len=StringLength)                 :: Char_TypeZUV
        
        !Begin-----------------------------------------------------------------

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR01'

        call GetData(Operation,                                                         &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OPERATION',                                        &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR02'
        
        if (iflag == 0)then
            write(*,*)'You must specify the input file path'
            stop 'ReadOptions - HDFOperator - ERR03'
        end if

        call GetData(FirstHDF5File%FileName,                                            &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'FIRST_HDF_FILE',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR04'

        if (iflag == 0)then
            write(*,*)'You must specify the input file path #1'
            stop 'ReadOptions - HDFOperator - ERR05'
        end if

        if(Operation .ne. Sum)then

            call GetData(SecondHDF5File%FileName,                                       &
                         ObjEnterData, iflag,                                           &
                         SearchType   = FromFile,                                       &
                         keyword      = 'SECOND_HDF_FILE',                              &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR06'

            if (iflag == 0)then
                write(*,*)'You must specify the input file path #2'
                stop 'ReadOptions - HDFOperator - ERR07'
            end if

        end if

        call GetData(OutputHDF5File%FileName,                                           &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUT_FILE',                                      &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR08'

        if (iflag == 0)then
            write(*,*)'You must specify the output file path'
            stop 'ReadOptions - HDFOperator - HDFOperator - ERR09'
        end if

        call GetData(CopyOpenPoints,                                                    &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'COPY_OPEN_POINTS',                                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR10'
        
        call GetData(CopyTime,                                                          &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'COPY_TIME',                                        &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR11'


        call GetData(InstantsPath,                                                      &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'INSTANTS_PATH',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR12'

        if(Operation .eq. Sum)then

            call GetData(OutputGroup,                                                   &
                         ObjEnterData, iflag,                                           &
                         SearchType   = FromFile,                                       &
                         keyword      = 'OUTPUT_GROUP',                                 &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR13'

            call GetData(OutputDataSet,                                                 &
                         ObjEnterData, iflag,                                           &
                         SearchType   = FromFile,                                       &
                         keyword      = 'OUTPUT_DATASET',                               &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR14'
        
        end if

        call GetData(MappingGroup,                                                      &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MAPPING_GROUP',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR15'

        call GetData(MappingDataSet,                                                    &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MAPPING_DATASET',                                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR16'

        call GetData(DataSetWithNumber,                                                 &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'READ_DATASET_WITH_NUMBER',                         &
                     default      = .true.,                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR17'

        call GetData(RelativeDifferences,                                               &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'COMPUTE_RELATIVE_DIFFERENCES',                     &
                     default      = .true.,                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR18'

        call GetData(IgnoreInstants,                                                    &
                     ObjEnterData, iflag,                                               &
                     SearchType   = FromFile,                                           &
                     keyword      = 'IGNORE_INSTANTS',                                  &
                     default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR19'

do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                     &
                                        '<beginproperty>', '<endproperty>', BlockFound, &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR20'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then
if2 :           if (BlockFound) then

                    nullify(NewProperty)

                    call AddProperty(NewProperty)

                    call GetData(NewProperty%Name,                                      &
                                 ObjEnterData, iflag,                                   &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'NAME',                                 &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR21'
                    if (iflag == 0)then
                        write(*,*)'Must specify name of the property'
                        stop 'ReadOptions - HDFOperator - ERR15'
                    end if

                    call GetData(NewProperty%Path,                                      &
                                 ObjEnterData, iflag,                                   &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'PATH',                                 &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR22'
                    if (iflag == 0)then
                        write(*,*)'Must specify path of the property'
                        stop 'ReadOptions - HDFOperator - ERR21'
                    end if

                    call GetData(NewProperty%OutputName,                                &
                                 ObjEnterData, iflag,                                   &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'OUTPUT_NAME',                          &
                                 Default      = NewProperty%Name,                       &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR23'

                    call GetData(Char_TypeZUV, ObjEnterData, iflag,                     &
                                 keyword        = 'TYPE_ZUV',                           &
                                 SearchType     = FromBlock,                            &
                                 default        = "Z",                                  &
                                 STAT           = STAT_CALL)            
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'ReadOptions - HDFOperator - ERR23b'

                    NewProperty%TypeZUV  = TranslateTypeZUV(Char_TypeZUV)

                    if (NewProperty%TypeZUV == TypeU_ .and. .not. PropertiesInFacesU)   &
                        PropertiesInFacesU = .true.

                    if (NewProperty%TypeZUV == TypeV_ .and. .not. PropertiesInFacesV)   &
                        PropertiesInFacesV = .true.

                    if(Operation .eq. Sum)then

                        call GetData(NewProperty%Factor,                                &
                                     ObjEnterData, iflag,                               &
                                     SearchType   = FromBlock,                          &
                                     keyword      = 'FACTOR',                           &
                                     Default      = 1.,                                 &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR24'

                    end if

                else

                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'ReadOptions - HDFOperator - ERR25'
                        
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'ReadOptions - HDFOperator - ERR26'
                    
            end if if1
        end do do1

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - HDFOperator - ERR27'

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine OpenOutputFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, instant
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        real,    dimension(:    ), pointer          :: Array1D
        real,    dimension(:,:  ), pointer          :: Array2D
        integer, dimension(:,:,:), pointer          :: Array3D

        !Begin-----------------------------------------------------------------

        ILB = FirstHDF5File%Size%ILB
        IUB = FirstHDF5File%Size%IUB
        JLB = FirstHDF5File%Size%JLB
        JUB = FirstHDF5File%Size%JUB
        KLB = FirstHDF5File%Size%KLB
        KUB = FirstHDF5File%Size%KUB

        call GetHDF5FileAccess  (HDF5_READ      = HDF5_READ,                            &
                                 HDF5_READWRITE = HDF5_READWRITE,                       &
                                 HDF5_CREATE    = HDF5_CREATE)

        call ConstructHDF5 (OutputHDF5File%ObjHDF5, OutputHDF5File%FileName,            &
                            HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR01'
        
        allocate(Array1D  (1:6                         ))
        allocate(Array2D  (ILB:IUB+1, JLB:JUB+1        )) 
        allocate(Array3D  (ILB:IUB,   JLB:JUB,  KLB:KUB))
        

        !Bathymetry
        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR02'

        call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Grid",                              &
                          "Bathymetry",                                                 &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR03'
        
        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR04'

        call HDF5WriteData (OutputHDF5File%ObjHDF5,"/Grid",                             &
                            "Bathymetry",  "m",                                         &
                            Array2D      = Array2D,                                     &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR05'


        !ConnectionX
        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR06'

        call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Grid",                              &
                          "ConnectionX",                                                &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR07'
        
        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR08'

        call HDF5WriteData(OutputHDF5File%ObjHDF5, "/Grid",                             &
                          "ConnectionX", "m",                                           &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR09'

        !ConnectionX
        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR10'

        call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Grid",                              &
                          "ConnectionY",                                                &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR11'
        
        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR12'

        call HDF5WriteData(OutputHDF5File%ObjHDF5, "/Grid",                             &
                          "ConnectionY", "m",                                           &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR13'


        !Latitude
        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR14'

        call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Grid",                              &
                          "Latitude",                                                   &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR15'

        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR16'

        call HDF5WriteData(OutputHDF5File%ObjHDF5, "/Grid",                             &
                          "Latitude",  "º",                                             &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR17'

        !Longitude
        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR18'

        call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Grid",                              &
                          "Longitude",                                                  &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR19'

        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB+1, JLB, JUB+1)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR20'

        call HDF5WriteData(OutputHDF5File%ObjHDF5, "/Grid",                             &
                          "Longitude",  "º",                                            &
                          Array2D      = Array2D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR21'


        !WaterPoints
        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR22'
        
        call HDF5ReadData (FirstHDF5File%ObjHDF5, trim(MappingGroup),                   &
                          trim(MappingDataSet),                                         &
                          Array3D      = Array3D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR23'
        
        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR24'

        call HDF5WriteData(OutputHDF5File%ObjHDF5, trim(MappingGroup),                  &
                          trim(MappingDataSet), "-",                                    &
                          Array3D      = Array3D,                                       &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR25'

        do instant = 1, FirstHDF5File%nInstants

            if(CopyOpenPoints)then
            
            !OpenPoints
                call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR26'

                call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Grid/OpenPoints",           &
                                  "OpenPoints",                                         &
                                  Array3D      = Array3D,                               &
                                  OutputNumber = instant,                               &
                                  STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR05'

                call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR27'

                call HDF5WriteData(OutputHDF5File%ObjHDF5, "/Grid/OpenPoints",          &
                                  "OpenPoints", "-",                                    &
                                  Array3D      = Array3D,                               &
                                  OutputNumber = instant,                               &
                                  STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR28'
            
            end if

            if(CopyTime)then
            
                !Time
                call HDF5SetLimits(FirstHDF5File%ObjHDF5, 1, 6)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR29'

                call HDF5ReadData (FirstHDF5File%ObjHDF5, "/Time",                      &
                                  "Time",                                               &
                                  Array1D      = Array1D,                               &
                                  OutputNumber = instant,                               &
                                  STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR30'

                call HDF5SetLimits(OutputHDF5File%ObjHDF5, 1, 6)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR32'

                call HDF5WriteData(OutputHDF5File%ObjHDF5, "/Time",                     &
                                  "Time", "-",                                          &
                                  Array1D      = Array1D,                               &
                                  OutputNumber = instant,                               &
                                  STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR33'
            
            end if

        end do

        call HDF5FlushMemory(OutputHDF5File%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenOutputFile - HDF5Operator - ERR34'

        deallocate(Array1D, Array2D, Array3D)

    end subroutine OpenOutputFile

    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        type(T_Size3D)                  :: Size

        !Begin-----------------------------------------------------------------
        
        Size = FirstHDF5File%Size

        allocate(Mapping3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB, Size%KLB:Size%KUB))

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------
    
    subroutine CheckFiles

        !Local-----------------------------------------------------------------
        logical                                     :: Exist = ON
        type(T_Property), pointer                   :: Property
        character(len=PathLength)                   :: PropertyName
        character(StringLength)                     :: PropertyUnits
        integer                                     :: PropertyRank
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        if(FirstHDF5File%Size%ILB .ne. SecondHDF5File%Size%ILB)then
            write(*,*)
            write(*,*)"ILB in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR01'
        end if

        if(FirstHDF5File%Size%JLB .ne. SecondHDF5File%Size%JLB)then
            write(*,*)
            write(*,*)"JLB in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR02'
        end if
        
        if(FirstHDF5File%Size%KLB .ne. SecondHDF5File%Size%KLB)then
            write(*,*)
            write(*,*)"KLB in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR03'
        end if
        
        if(FirstHDF5File%Size%IUB .ne. SecondHDF5File%Size%IUB)then
            write(*,*)
            write(*,*)"IUB in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR04'
        end if
        
        if(FirstHDF5File%Size%JUB .ne. SecondHDF5File%Size%JUB)then
            write(*,*)
            write(*,*)"JUB in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR05'
        end if
        
        if(FirstHDF5File%Size%KUB .ne. SecondHDF5File%Size%KUB)then
            write(*,*)
            write(*,*)"KUB in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR06'
        end if

        if(FirstHDF5File%nInstants .ne. SecondHDF5File%nInstants)then
            write(*,*)
            write(*,*)"Number of instants in file 1 and file 2 are different"
            stop 'CheckFiles - HDF5Operator - ERR07'
        end if

        !Check if requested properties exist in second HDF5 and check rank and units
        Property => FirstProperty

        do while(associated(Property))

            call GetHDF5GroupExist (SecondHDF5File%ObjHDF5, Property%Path, Exist)

            if (.not. Exist) then

                write(*,*)
                write(*,*)'Expected property variable which does not exist in group:'
                write(*,*) trim(Property%Path)
                write(*,*)'HDF5 file:'
                write(*,*) trim(SecondHDF5File%FileName)
                stop 'CheckFiles - HDF5Operator - ERR08'
            
            else
                !Get rank of variable
                call GetHDF5GroupID(SecondHDF5File%ObjHDF5, Property%Path,              &
                                    1, PropertyName, PropertyUnits,                     &
                                    PropertyRank, STAT = STAT_CALL)             
                if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'CheckFiles - HDF5Operator - ERR09'

                if (PropertyRank .NE. Property%Rank) then

                    write(*,*)
                    write(*,*)'Property with different rank in the two HDF5 files:'
                    write(*,*) trim(Property%Path)
                    stop 'CheckFiles - HDF5Operator - ERR11'
                endif

                if (PropertyUnits .NE. Property%Units) then
                    write(*,*)'Property with different units in the two HDF5 files:'
                    write(*,*) trim(Property%Path)
                    stop 'CheckFiles - HDF5Operator - ERR12'
                end if
            endif

            Property => Property%Next
        end do

    end subroutine CheckFiles
    
    !--------------------------------------------------------------------------

    subroutine AddProperty (ObjProperty)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                   :: ObjProperty

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                   :: NewProperty
        type (T_Property), pointer                   :: PreviousProperty
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewProperty)
        nullify  (NewProperty%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstProperty)) then
            FirstProperty         => NewProperty
            ObjProperty           => NewProperty
        else
            PreviousProperty      => FirstProperty
            ObjProperty           => FirstProperty%Next
            do while (associated(ObjProperty))
                PreviousProperty  => ObjProperty
                ObjProperty       => ObjProperty%Next
            enddo
            ObjProperty           => NewProperty
            PreviousProperty%Next => NewProperty
        endif

    end subroutine AddProperty

    !--------------------------------------------------------------------------

    subroutine OpenHDF5File(HDF5File)
        
        !Arguments-------------------------------------------------------------
        type(T_HDF5File)                            :: HDF5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetHDF5FileAccess  (HDF5_READ      = HDF5_READ,                            &
                                 HDF5_READWRITE = HDF5_READWRITE,                       &
                                 HDF5_CREATE    = HDF5_CREATE)

        call ConstructHDF5 (HDF5File%ObjHDF5, HDF5File%FileName, HDF5_READ,             &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenHDF5File - HDF5Operator - ERR01'
        
        if (.not. IgnoreInstants) then
            call GetHDF5GroupNumberOfItems(HDF5File%ObjHDF5, trim(InstantsPath),        &
                                           HDF5File%nInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'OpenHDF5File - HDF5Operator - ERR02'
        else
            HDF5File%nInstants = 1
        endif

        call ObtainDimensions(HDF5File)

    end subroutine OpenHDF5File

    !--------------------------------------------------------------------------

    subroutine ObtainDimensions(HDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File)                            :: HDF5File


        !Local-----------------------------------------------------------------
        integer  , dimension(7)                     :: Dimensions
        integer                                     :: nItems, item
        character(len=StringLength)                 :: AuxName
        logical                                     :: Exist = ON
        integer                                     :: STAT_CALL, Rank
        
        !Begin-----------------------------------------------------------------

        !Check if provided name for water points data set exist
        call GetHDF5DataSetExist(HDF5File%ObjHDF5,                                      &
            trim(MappingGroup)//trim(MappingDataSet), exist, STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ObtainDimensions - HDF5Operator - ERR01'

        if (.not. Exist) then

            write(*,*)
            write(*,*)'Expected mapping variable which does not exist in group:'
            write(*,*) trim(MappingGroup)//trim(MappingDataSet)
            stop 'ObtainDimensions - HDF5Operator - ERR02'

        else

            !find position for data set
            call GetHDF5GroupNumberOfItems(HDF5File%ObjHDF5, "/Grid", nItems, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ObtainDimensions - HDF5Operator - ERR03'

            do item = 1, nItems

                call GetHDF5GroupID(HDF5File%ObjHDF5,                                   &
                                    FatherGroupName = "/Grid",                          &
                                    GroupPosition   = item,                             &
                                    GroupName       = AuxName,                          & 
                                    Rank            = Rank,                             &
                                    Dimensions      = Dimensions,                       &
                                    STAT            = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ObtainDimensions - HDF5Operator - ERR04'

                if (trim(AuxName) == trim(MappingDataSet)) then
                    exit
                endif

            end do

            select case (Rank)

                case(2)

                    HDF5File%Size%ILB = 1
                    HDF5File%Size%IUB = Dimensions(1)
                    HDF5File%Size%JLB = 1
                    HDF5File%Size%JUB = Dimensions(2)
                    HDF5File%Size%KLB = 1
                    HDF5File%Size%KUB = Dimensions(2)

                case(3)
                    
                    HDF5File%Size%ILB = 1
                    HDF5File%Size%IUB = Dimensions(1)
                    HDF5File%Size%JLB = 1
                    HDF5File%Size%JUB = Dimensions(2)
                    HDF5File%Size%KLB = 1
                    HDF5File%Size%KUB = Dimensions(3)

            case default 

                write(*,*)'Mapping variable must be 2D or 3D.'
                stop 'ObtainDimensions - HDF5Operator - ERR05'

            end select

        end if

    end subroutine ObtainDimensions

    !--------------------------------------------------------------------------

    subroutine CheckPropertiesInFile(HDF5File)

        !Arguments-------------------------------------------------------------
        type(T_HDF5File)                            :: HDF5File


        !Local-----------------------------------------------------------------
        logical                                     :: Exist = ON
        integer                                     :: STAT_CALL
        character(len=PathLength)                   :: PropertyName
        
        !Begin-----------------------------------------------------------------

        !Check if requested properties exist in HDF5 and get rank and units
        Property => FirstProperty

        do while(associated(Property))

            call GetHDF5GroupExist (HDF5File%ObjHDF5, Property%Path, Exist)

            if (.not. Exist) then

                write(*,*)
                write(*,*)'Expected property variable which does not exist in group:'
                write(*,*) trim(Property%Path)
                write(*,*)'HDF5 file:'
                write(*,*) trim(HDF5File%FileName)
                stop 'CheckPropertiesInFile - HDF5Operator - ERR01'
            
            else
                !Get rank of variable
                call GetHDF5GroupID(HDF5File%ObjHDF5, Property%Path,                    &
                                    1, PropertyName, Property%Units,                    &
                                    Property%Rank, STAT = STAT_CALL)             
                if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'CheckPropertiesInFile - HDF5Operator - ERR02'

                if (Property%Rank == 2 .and. .not. Properties2D)                        &
                    Properties2D = .true.

            endif

            Property => Property%Next
        end do

    end subroutine CheckPropertiesInFile

    !--------------------------------------------------------------------------

    subroutine DoOperations

        !Local-----------------------------------------------------------------
        
        !Begin-----------------------------------------------------------------

        select case(Operation)

            case(Differences)

                call ComputeDifferences

            case(Sum)

                call ComputeSum

            case default

                stop 'An unknown operation was specified - HDF5Operator - ERR01'

        end select

    end subroutine DoOperations

    !--------------------------------------------------------------------------

    subroutine ComputeDifferences
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer          :: RelativeDiffField
        real,    dimension(:,:,:), pointer          :: AbsoluteDiffField
        real,    dimension(:,:,:), pointer          :: FirstPropertyField
        real,    dimension(:,:,:), pointer          :: SecondPropertyField
        real,    dimension(:,:,:), pointer          :: RelativeDiffFieldFacesU
        real,    dimension(:,:,:), pointer          :: AbsoluteDiffFieldFacesU
        real,    dimension(:,:,:), pointer          :: FirstPropertyFieldFacesU
        real,    dimension(:,:,:), pointer          :: SecondPropertyFieldFacesU
        real,    dimension(:,:,:), pointer          :: RelativeDiffFieldFacesV
        real,    dimension(:,:,:), pointer          :: AbsoluteDiffFieldFacesV
        real,    dimension(:,:,:), pointer          :: FirstPropertyFieldFacesV
        real,    dimension(:,:,:), pointer          :: SecondPropertyFieldFacesV
        real,    dimension(:,:)  , pointer          :: RelativeDiffField2D
        real,    dimension(:,:)  , pointer          :: AbsoluteDiffField2D
        real,    dimension(:,:)  , pointer          :: FirstPropertyField2D
        real,    dimension(:,:)  , pointer          :: SecondPropertyField2D
        integer                                     :: STAT_CALL, instant, i, j, k
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------
        
        ILB = FirstHDF5File%Size%ILB
        IUB = FirstHDF5File%Size%IUB
        JLB = FirstHDF5File%Size%JLB
        JUB = FirstHDF5File%Size%JUB
        KLB = FirstHDF5File%Size%KLB
        KUB = FirstHDF5File%Size%KUB

        if (RelativeDifferences) then
            allocate(RelativeDiffField  (ILB:IUB,JLB:JUB,KLB:KUB))
            RelativeDiffField   = 0.

            if (PropertiesInFacesU) then
                allocate(RelativeDiffFieldFacesU  (ILB:IUB,JLB:JUB+1,KLB:KUB))
                RelativeDiffFieldFacesU   = 0.
            endif

            if (PropertiesInFacesV) then
                allocate(RelativeDiffFieldFacesV  (ILB:IUB+1,JLB:JUB,KLB:KUB))
                RelativeDiffFieldFacesV   = 0.
            endif

            if (Properties2D) then
                allocate(RelativeDiffField2D  (ILB:IUB,JLB:JUB))
                RelativeDiffField2D   = 0.
            endif

        endif

        allocate(AbsoluteDiffField  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(FirstPropertyField (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(SecondPropertyField(ILB:IUB,JLB:JUB,KLB:KUB))

        AbsoluteDiffField   = 0.
        FirstPropertyField  = null_real
        SecondPropertyField = null_real

        if (PropertiesInFacesU) then
            allocate(AbsoluteDiffFieldFacesU  (ILB:IUB,JLB:JUB+1,KLB:KUB))
            allocate(FirstPropertyFieldFacesU (ILB:IUB,JLB:JUB+1,KLB:KUB))
            allocate(SecondPropertyFieldFacesU (ILB:IUB,JLB:JUB+1,KLB:KUB))

            AbsoluteDiffFieldFacesU   = 0.
            FirstPropertyFieldFacesU  = null_real
            SecondPropertyFieldFacesU = null_real
        endif

        if (PropertiesInFacesV) then
            allocate(AbsoluteDiffFieldFacesV  (ILB:IUB+1,JLB:JUB,KLB:KUB))
            allocate(FirstPropertyFieldFacesV (ILB:IUB+1,JLB:JUB,KLB:KUB))
            allocate(SecondPropertyFieldFacesV (ILB:IUB+1,JLB:JUB,KLB:KUB))

            AbsoluteDiffFieldFacesV   = 0.
            FirstPropertyFieldFacesV  = null_real
            SecondPropertyFieldFacesV = null_real
        endif

        if (Properties2D) then
            allocate(AbsoluteDiffField2D  (ILB:IUB,JLB:JUB))
            allocate(FirstPropertyField2D (ILB:IUB,JLB:JUB))
            allocate(SecondPropertyField2D (ILB:IUB,JLB:JUB))

            AbsoluteDiffField2D   = 0.
            FirstPropertyField2D  = null_real
            SecondPropertyField2D = null_real
        endif

        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputeDifferences - HDF5Operator - ERR01'

        !call HDF5SetLimits(SecondHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
        !if (STAT_CALL .NE. SUCCESS_) stop 'ComputeDifferences - HDF5Operator - ERR02'

        !call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
        !if (STAT_CALL .NE. SUCCESS_) stop 'ComputeDifferences - HDF5Operator - ERR01'

        call HDF5ReadData (FirstHDF5File%ObjHDF5, trim(MappingGroup),                   &
                          trim(MappingDataSet),                                         &
                          Array3D      = Mapping3D,                                     &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputeDifferences - HDF5Operator - ERR02'

        do instant = 1, FirstHDF5File%nInstants

            Property => FirstProperty

            do while(associated(Property))

                if (Property%TypeZUV == TypeZ_) then

                    if (Property%Rank == 2) then

                        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR03'

                        call HDF5SetLimits(SecondHDF5File%ObjHDF5, ILB, IUB, JLB, JUB)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR04'

                        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR05'

                        if(DataSetWithNumber)then

                            call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path), &
                                              trim(Property%Name),                      &
                                              Array2D      = FirstPropertyField2D,      &
                                              OutputNumber = instant,                   &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR06'
                    
                            call HDF5ReadData(SecondHDF5File%ObjHDF5, trim(Property%Path), &
                                              trim(Property%Name),                      &
                                              Array2D      = SecondPropertyField2D,     &
                                              OutputNumber = instant,                   &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR07'

                        else

                            call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path), &
                                              trim(Property%Name),                      &
                                              Array2D      = FirstPropertyField2D,      &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR08'
                    
                            call HDF5ReadData(SecondHDF5File%ObjHDF5, trim(Property%Path), &
                                              trim(Property%Name),                      &
                                              Array2D      = SecondPropertyField2D,     &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR09'

                        end if

                        !(assumes 2D properties at surface)
                        do j = JLB, JUB
                        do i = ILB, IUB

                            if(Mapping3D(i,j,KUB) == WaterPoint)then

                                !RelativeDiffField(i,j,k) = (SecondPropertyField(i,j,k) - FirstPropertyField(i,j,k)) / &
                                !                            FirstPropertyField(i,j,k)

                                AbsoluteDiffField2D(i,j) = FirstPropertyField2D(i,j) -  &
                                    SecondPropertyField2D(i,j)

                            end if

                        enddo
                        enddo

                        if (RelativeDifferences) then

                            do j = JLB, JUB
                            do i = ILB, IUB

                                if(Mapping3D(i,j,KUB) == WaterPoint)then

                                    RelativeDiffField2D(i,j) =                          &
                                        (SecondPropertyField2D(i,j) -                   &
                                        FirstPropertyField2D(i,j)) /                    &
                                        FirstPropertyField2D(i,j)
                                end if

                            enddo
                            enddo
                        endif

                        if(DataSetWithNumber)then

                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                           "/AbsoluteDifference/"//trim(Property%Name), &
                                             trim(Property%OutputName), Property%Units, &
                                               Array2D = AbsoluteDiffField2D,           &
                                               OutputNumber = instant,                  &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR10'
                
                            if (RelativeDifferences) then
                                call HDF5WriteData(OutputHDF5File%ObjHDF5,              &
                                           "/RelativeDifference/"//trim(Property%Name), &
                                                   trim(Property%OutputName), "-",      & 
                                                   Array2D = RelativeDiffField2D,       &
                                                   OutputNumber = instant,              &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                                    stop 'ComputeDifferences - HDF5Operator - ERR11'
                            endif

                        else

                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                           "/AbsoluteDifference/"//trim(Property%Name), &
                                               trim(Property%OutputName), Property%Units, &
                                               Array2D = AbsoluteDiffField2D,           &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR12'

                            if (RelativeDifferences) then
                                call HDF5WriteData(OutputHDF5File%ObjHDF5,              &
                                           "/RelativeDifference/"//trim(Property%Name), &
                                                   trim(Property%OutputName), "-",      &
                                                   Array2D = RelativeDiffField2D,       &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                                    stop 'ComputeDifferences - HDF5Operator - ERR13'
                            endif

                        end if

                    else

                        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB,   &
                                           KLB, KUB)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR14'

                        call HDF5SetLimits(SecondHDF5File%ObjHDF5, ILB, IUB, JLB, JUB,  &
                                           KLB, KUB)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR15'

                        call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB,  &
                                           KLB, KUB)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR16'

                        if(DataSetWithNumber)then

                            call HDF5ReadData(FirstHDF5File%ObjHDF5,                    &
                                              trim(Property%Path),                      &
                                              trim(Property%Name),                      &
                                              Array3D      = FirstPropertyField,        &
                                              OutputNumber = instant,                   &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR17'
                    
                            call HDF5ReadData(SecondHDF5File%ObjHDF5,                   &
                                              trim(Property%Path),                      &
                                              trim(Property%Name),                      &
                                              Array3D      = SecondPropertyField,       &
                                              OutputNumber = instant,                   &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR18'

                        else

                            call HDF5ReadData(FirstHDF5File%ObjHDF5,                    &
                                              trim(Property%Path),                      &
                                              trim(Property%Name),                      &
                                              Array3D      = FirstPropertyField,        &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR19'
                    
                            call HDF5ReadData(SecondHDF5File%ObjHDF5,                   &
                                              trim(Property%Path),                      &
                                              trim(Property%Name),                      &
                                              Array3D      = SecondPropertyField,       &
                                              STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR20'

                        end if

                        do j = JLB, JUB
                        do i = ILB, IUB
                        do k = JLB, KUB

                            if(Mapping3D(i,j,k) == WaterPoint)then

                                !RelativeDiffField(i,j,k) = (SecondPropertyField(i,j,k) - FirstPropertyField(i,j,k)) / &
                                !                            FirstPropertyField(i,j,k)

                                AbsoluteDiffField(i,j,k) = FirstPropertyField(i,j,k) -  &
                                    SecondPropertyField(i,j,k)

                            end if

                        enddo
                        enddo
                        enddo

                        if (RelativeDifferences) then

                            do j = JLB, JUB
                            do i = ILB, IUB
                            do k = JLB, KUB

                                if(Mapping3D(i,j,k) == WaterPoint)then

                                    RelativeDiffField(i,j,k) =                          &
                                        (SecondPropertyField(i,j,k) -                   &
                                        FirstPropertyField(i,j,k)) /                    &
                                        FirstPropertyField(i,j,k)
                                end if

                            enddo
                            enddo
                            enddo
                        endif

                        if(DataSetWithNumber)then

                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                               "/AbsoluteDifference/"//trim(Property%Name), &
                                               trim(Property%OutputName), Property%Units, &
                                               Array3D = AbsoluteDiffField,             &
                                               OutputNumber = instant,                  &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR21'
                
                            if (RelativeDifferences) then
                                call HDF5WriteData(OutputHDF5File%ObjHDF5,              &
                                               "/RelativeDifference/"//trim(Property%Name), &
                                                   trim(Property%OutputName), "-",      & 
                                                   Array3D = RelativeDiffField,         &
                                                   OutputNumber = instant,              &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                                    stop 'ComputeDifferences - HDF5Operator - ERR22'
                            endif

                        else

                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                               "/AbsoluteDifference/"//trim(Property%Name), &
                                               trim(Property%OutputName), Property%Units, &
                                               Array3D = AbsoluteDiffField,             &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR23'

                            if (RelativeDifferences) then
                                call HDF5WriteData(OutputHDF5File%ObjHDF5,              &
                                               "/RelativeDifference/"//trim(Property%Name), &
                                                   trim(Property%OutputName), "-",      &
                                                   Array3D = RelativeDiffField,         &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_)                            &
                                    stop 'ComputeDifferences - HDF5Operator - ERR24'
                            endif

                        end if
                    endif

                elseif (Property%TypeZUV == TypeU_) then

                    call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB + 1,   &
                                       KLB, KUB)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ComputeDifferences - HDF5Operator - ERR25'

                    call HDF5SetLimits(SecondHDF5File%ObjHDF5, ILB, IUB, JLB, JUB + 1,  &
                                       KLB, KUB)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ComputeDifferences - HDF5Operator - ERR26'

                    call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB + 1,  &
                                       KLB, KUB)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ComputeDifferences - HDF5Operator - ERR27'

                    if(DataSetWithNumber)then

                        call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path),   &
                                          trim(Property%Name),                          &
                                          Array3D      = FirstPropertyFieldFacesU,      &
                                          OutputNumber = instant,                       &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR28'
                    
                        call HDF5ReadData(SecondHDF5File%ObjHDF5, trim(Property%Path),  &
                                          trim(Property%Name),                          &
                                          Array3D      = SecondPropertyFieldFacesU,     &
                                          OutputNumber = instant,                       &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR29'

                    else

                        call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path),   &
                                          trim(Property%Name),                          &
                                          Array3D      = FirstPropertyFieldFacesU,      &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR30'
                    
                        call HDF5ReadData(SecondHDF5File%ObjHDF5, trim(Property%Path),  &
                                          trim(Property%Name),                          &
                                          Array3D      = SecondPropertyFieldFacesU,     &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR31'

                    end if

                    do i = ILB, IUB
                    do k = JLB, KUB

                        if (Mapping3D(i,JLB,k) == WaterPoint) then

                            AbsoluteDiffFieldFacesU(i,JLB,k) =                          &
                                FirstPropertyFieldFacesU(i,JLB,k) -                     &
                                SecondPropertyFieldFacesU(i,JLB,k)

                        end if

                        if (Mapping3D(i,JUB,k) == WaterPoint) then

                            AbsoluteDiffFieldFacesU(i,JUB + 1,k) =                      &
                                FirstPropertyFieldFacesU(i,JUB + 1,k) -                 &
                                SecondPropertyFieldFacesU(i,JUB + 1,k)

                        end if

                    enddo
                    enddo

                    do j = JLB + 1, JUB
                    do i = ILB, IUB
                    do k = JLB, KUB

                        if((Mapping3D(i,j-1,k) == WaterPoint) .and.                     &
                           (Mapping3D(i,j,k) == WaterPoint))then

                            AbsoluteDiffFieldFacesU(i,j,k) =                            &
                                FirstPropertyFieldFacesU(i,j,k) -                       &
                                SecondPropertyFieldFacesU(i,j,k)

                        end if

                    enddo
                    enddo
                    enddo

                    if (RelativeDifferences) then

                        do i = ILB, IUB
                        do k = JLB, KUB

                            if (Mapping3D(i,JLB,k) == WaterPoint) then

                                AbsoluteDiffFieldFacesU(i,JLB,k) =                      &
                                    FirstPropertyFieldFacesU(i,JLB,k) -                 &
                                    SecondPropertyFieldFacesU(i,JLB,k)

                            end if

                            if (Mapping3D(i,JUB,k) == WaterPoint) then

                                AbsoluteDiffFieldFacesU(i,JUB + 1,k) =                  &
                                    FirstPropertyFieldFacesU(i,JUB + 1,k) -             &
                                    SecondPropertyFieldFacesU(i,JUB + 1,k)
                            end if

                        enddo
                        enddo


                        do j = JLB + 1, JUB
                        do i = ILB, IUB
                        do k = JLB, KUB

                            if((Mapping3D(i,j-1,k) == WaterPoint) .and.                 &
                           (Mapping3D(i,j,k) == WaterPoint))then

                                RelativeDiffFieldFacesU(i,j,k) =                        &
                                    (SecondPropertyFieldFacesU(i,j,k) -                 &
                                    FirstPropertyFieldFacesU(i,j,k)) /                  &
                                    FirstPropertyFieldFacesU(i,j,k)
                            end if

                        enddo
                        enddo
                        enddo
                    endif

                    if(DataSetWithNumber)then

                        call HDF5WriteData(OutputHDF5File%ObjHDF5,                      &
                                           "/AbsoluteDifference/"//trim(Property%Name), &
                                           trim(Property%OutputName), Property%Units,   &
                                           Array3D = AbsoluteDiffFieldFacesU,           &
                                           OutputNumber = instant,                      &
                                           STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR32'
                
                        if (RelativeDifferences) then
                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                           "/RelativeDifference/"//trim(Property%Name), &
                                               trim(Property%OutputName), "-",          & 
                                               Array3D = RelativeDiffFieldFacesU,       &
                                               OutputNumber = instant,                  &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR33'
                        endif

                    else

                        call HDF5WriteData(OutputHDF5File%ObjHDF5,                      &
                                           "/AbsoluteDifference/"//trim(Property%Name), &
                                           trim(Property%OutputName), Property%Units,   &
                                           Array3D = AbsoluteDiffFieldFacesU,           &
                                           STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR34'

                        if (RelativeDifferences) then
                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                           "/RelativeDifference/"//trim(Property%Name), &
                                               trim(Property%OutputName), "-",          &
                                               Array3D = RelativeDiffFieldFacesU,       &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR35'
                        endif

                    end if

                else !(TypeV_)

                    call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB + 1, JLB, JUB,   &
                                       KLB, KUB)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ComputeDifferences - HDF5Operator - ERR36'

                    call HDF5SetLimits(SecondHDF5File%ObjHDF5, ILB, IUB + 1, JLB, JUB,  &
                                       KLB, KUB)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ComputeDifferences - HDF5Operator - ERR37'

                    call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB + 1, JLB, JUB,  &
                                       KLB, KUB)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ComputeDifferences - HDF5Operator - ERR38'

                    if(DataSetWithNumber)then

                        call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path),   &
                                          trim(Property%Name),                          &
                                          Array3D      = FirstPropertyFieldFacesV,      &
                                          OutputNumber = instant,                       &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR39'
                    
                        call HDF5ReadData(SecondHDF5File%ObjHDF5, trim(Property%Path),  &
                                          trim(Property%Name),                          &
                                          Array3D      = SecondPropertyFieldFacesV,     &
                                          OutputNumber = instant,                       &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR40'

                    else

                        call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path),   &
                                          trim(Property%Name),                          &
                                          Array3D      = FirstPropertyFieldFacesV,      &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR41'
                    
                        call HDF5ReadData(SecondHDF5File%ObjHDF5, trim(Property%Path),  &
                                          trim(Property%Name),                          &
                                          Array3D      = SecondPropertyFieldFacesV,     &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR42'

                    end if

                    do j = JLB, JUB
                    do k = JLB, KUB

                        if (Mapping3D(ILB,j,k) == WaterPoint) then

                            AbsoluteDiffFieldFacesV(ILB,j,k) =                          &
                                FirstPropertyFieldFacesV(ILB,j,k) -                     &
                                SecondPropertyFieldFacesV(ILB,j,k)
                        end if

                        if (Mapping3D(IUB,j,k) == WaterPoint) then

                            AbsoluteDiffFieldFacesV(IUB+1,j,k) =                        &
                                FirstPropertyFieldFacesV(IUB+1,j,k) -                   &
                                SecondPropertyFieldFacesV(IUB+1,j,k)
                        end if

                    enddo
                    enddo

                    do j = JLB, JUB
                    do i = ILB + 1, IUB
                    do k = JLB, KUB

                        if((Mapping3D(i-1,j,k) == WaterPoint) .and.                     &
                           (Mapping3D(i,j,k) == WaterPoint))then

                            AbsoluteDiffFieldFacesV(i,j,k) =                            &
                                FirstPropertyFieldFacesV(i,j,k) -                       &
                                SecondPropertyFieldFacesV(i,j,k)

                        end if

                    enddo
                    enddo
                    enddo

                    if (RelativeDifferences) then

                        do j = JLB, JUB
                        do k = JLB, KUB

                            if (Mapping3D(ILB,j,k) == WaterPoint) then

                                AbsoluteDiffFieldFacesV(ILB,j,k) =                      &
                                    FirstPropertyFieldFacesV(ILB,j,k) -                 &
                                    SecondPropertyFieldFacesV(ILB,j,k)

                            end if

                            if (Mapping3D(IUB,j,k) == WaterPoint) then

                                AbsoluteDiffFieldFacesV(IUB + 1,j,k) =                  &
                                    FirstPropertyFieldFacesV(IUB + 1,j,k) -             &
                                    SecondPropertyFieldFacesV(IUB + 1,j,k)
                            end if

                        enddo
                        enddo


                        do j = JLB, JUB
                        do i = ILB + 1, IUB
                        do k = JLB, KUB

                            if((Mapping3D(i-1,j,k) == WaterPoint) .and.                 &
                           (Mapping3D(i,j,k) == WaterPoint))then

                                RelativeDiffFieldFacesV(i,j,k) =                        &
                                    (SecondPropertyFieldFacesV(i,j,k) -                 &
                                    FirstPropertyFieldFacesV(i,j,k)) /                  &
                                    FirstPropertyFieldFacesV(i,j,k)
                            end if

                        enddo
                        enddo
                        enddo
                    endif

                    if(DataSetWithNumber)then

                        call HDF5WriteData(OutputHDF5File%ObjHDF5,                      &
                                           "/AbsoluteDifference/"//trim(Property%Name), &
                                           trim(Property%OutputName), Property%Units,   &
                                           Array3D = AbsoluteDiffFieldFacesV,           &
                                           OutputNumber = instant,                      &
                                           STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'ComputeDifferences - HDF5Operator - ERR43'
                
                        if (RelativeDifferences) then
                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                          "/RelativeDifference/"//trim(Property%Name),  &
                                               trim(Property%OutputName), "-",          & 
                                               Array3D = RelativeDiffFieldFacesV,       &
                                               OutputNumber = instant,                  &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR44'
                        endif

                    else

                        call HDF5WriteData(OutputHDF5File%ObjHDF5,                      &
                                           "/AbsoluteDifference/"//trim(Property%Name), &
                                           trim(Property%OutputName), Property%Units,   &
                                           Array3D = AbsoluteDiffFieldFacesV,           &
                                           STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                    &
                            stop 'ComputeDifferences - HDF5Operator - ERR45'

                        if (RelativeDifferences) then
                            call HDF5WriteData(OutputHDF5File%ObjHDF5,                  &
                                           "/RelativeDifference/"//trim(Property%Name), &
                                               trim(Property%OutputName), "-",          &
                                               Array3D = RelativeDiffFieldFacesV,       &
                                               STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                &
                                stop 'ComputeDifferences - HDF5Operator - ERR46'
                        endif

                    end if

                endif

                Property => Property%Next

            enddo

        end do

        call HDF5FlushMemory(OutputHDF5File%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputeDifferences - HDF5Operator - ERR47'

        if (RelativeDifferences) then
            deallocate(RelativeDiffField  )
            if (PropertiesInFacesU) deallocate(RelativeDiffFieldFacesU  )
            if (PropertiesInFacesV) deallocate(RelativeDiffFieldFacesV  )
        endif

        deallocate(AbsoluteDiffField  )
        deallocate(FirstPropertyField )
        deallocate(SecondPropertyField)

        if (PropertiesInFacesU) then
            deallocate(AbsoluteDiffFieldFacesU  )
            deallocate(FirstPropertyFieldFacesU )
            deallocate(SecondPropertyFieldFacesU)
        endif

        if (PropertiesInFacesV) then
            deallocate(AbsoluteDiffFieldFacesV  )
            deallocate(FirstPropertyFieldFacesV )
            deallocate(SecondPropertyFieldFacesV)
        endif

        if (Properties2D) then
            deallocate(AbsoluteDiffField2D  )
            deallocate(FirstPropertyField2D )
            deallocate(SecondPropertyField2D)
        endif

    end subroutine ComputeDifferences
    
    !--------------------------------------------------------------------------

    subroutine ComputeSum
       
        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer          :: SumField, PropertyField
        integer                                     :: STAT_CALL
        integer                                     :: i, j, k, instant
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------
        
        ILB = FirstHDF5File%Size%ILB
        IUB = FirstHDF5File%Size%IUB
        JLB = FirstHDF5File%Size%JLB
        JUB = FirstHDF5File%Size%JUB
        KLB = FirstHDF5File%Size%KLB
        KUB = FirstHDF5File%Size%KUB

        !Auxiliar sum array
        allocate(SumField(ILB:IUB,JLB:JUB,KLB:KUB))
        SumField(:,:,:) = 0.

        allocate(PropertyField(ILB:IUB,JLB:JUB,KLB:KUB))
        PropertyField(:,:,:) = null_real

        call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR01'
        
        call HDF5ReadData (FirstHDF5File%ObjHDF5, trim(MappingGroup),                   &
                          trim(MappingDataSet),                                         &
                          Array3D      = Mapping3D,                                     &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR02'

        do instant = 1, FirstHDF5File%nInstants

            SumField(:,:,:) = 0.

            Property => FirstProperty

            do while(associated(Property))

                call HDF5SetLimits(FirstHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
                if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR01'

                if(DataSetWithNumber)then

                    call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path),       &
                                      trim(Property%Name),                              &
                                      Array3D      = PropertyField,                     &
                                      OutputNumber = instant,                           &
                                      STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR03'

                else

                    call HDF5ReadData(FirstHDF5File%ObjHDF5, trim(Property%Path),       &
                                      trim(Property%Name),                              &
                                      Array3D      = PropertyField,                     &
                                      STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR04'

                end if

                do j = JLB, JUB
                do i = ILB, IUB
                do k = JLB, KUB

                    if(Mapping3D(i,j,k) == WaterPoint)then

                        SumField(i,j,k) = SumField(i,j,k) + PropertyField(i,j,k) *      &
                                          Property%Factor

                    end if

                enddo
                enddo
                enddo

                Property => Property%Next

            enddo

            call HDF5SetLimits(OutputHDF5File%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB)
            if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR01'

            if(DataSetWithNumber)then

                call HDF5WriteData   (OutputHDF5File%ObjHDF5,                           &
                                      trim(OutputGroup),                                &
                                      trim(OutputDataSet), "-",                         &
                                      Array3D = SumField,                               &
                                      OutputNumber = instant,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR05'

            else

                call HDF5WriteData   (OutputHDF5File%ObjHDF5,                           &
                                      trim(OutputGroup),                                &
                                      trim(OutputDataSet), "-",                         &
                                      Array3D = SumField,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR06'

            end if

        end do

        call HDF5FlushMemory(OutputHDF5File%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ComputeSum - HDF5Operator - ERR07'

        deallocate(SumField, PropertyField)

    end subroutine ComputeSum

    !--------------------------------------------------------------------------

    subroutine Finish
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        deallocate(Mapping3D)

        call KillHDF5(FirstHDF5File%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Finish - HDF5Operator - ERR01'

        if(Operation .ne. Sum)then

            call KillHDF5(SecondHDF5File%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Finish - HDF5Operator - ERR02'

        end if

        call KillHDF5(OutputHDF5File%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Finish - HDF5Operator - ERR03'
    
        write(*,*)
        write(*,*)'Sucessfully terminated the program!'
        write(*,*)


    end subroutine Finish

end program HDF5Operator


