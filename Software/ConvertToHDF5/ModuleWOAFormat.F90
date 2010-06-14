!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid ConvertToHDF5
! MODULE        : WOAFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : March 2007
! REVISION      : Bartolomeu Bernardes - v4.0
! DESCRIPTION   : Module to convert World Ocean Atlas 2005 Format netcdf files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------

!DataFile
!
!   ACTION                      : char              -           !'CONVERT WOA FORMAT' to use this module
!   FILEPATH                    : char              -           !Path to WOA05 original netcdf files(n0112an1.nc for e.g.) 
!                                                               !with ending slash
!   OUTPUTFILENAME              : char              -           !Path to WOA05 file generated from NETCDF file
!
!   MOHID_UNITS                  : 0/1              -           !Convert uM to mg/l

!   <<begin_woafiles>>
!   (netcdf file)                 char              -           !WOA05 netcdf file in path defined in FILEPATH, one per line
!   (netcdf file2)                char                          ! e.g.: n0112an1.nc i0112an1.nc p0112an1.nc o0112an1.nc 
!    ...        
!   <<end_woafiles>>
!
!   CONVERT_TO_ASCII            : 0              [0]          !Flag to convert to ascii file: not implemented
!   CONVERT_TO_HDF5             : 1              [0]          !Flag to convert to hdf5 file
!
!   XX_VARIABLE                 : char              -           !'longitude'
!   YY_VARIABLE                 : char              -           !'latitude'
!   ZZ_VARIABLE                 : char              -           !'depth'
!   TIME_VARIABLE               : char              -           !'time'
! 
!   --------------------------------------------------------------------------------------------------
!   Will output into one hdf5 the vars on the netcdf files defined on block <<begin_woafiles>>
!
!   Refences:
!         ---WOA05 NAME---           ----Units-----   ---MOHID NAME---  ------Units-----------
!    
!       Temperature                      ºC            not imported
!       Salinity                        PPS            not imported
!       Dissolved oxygen                ml l-1         not imported           
!       Apparent oxygen utilization     ml l-1         not imported
!       Percentage oxygen saturation      %   dissolved oxygen percent saturation    %
!       Phosphate                       microM      inorganic phosphorus      microM or mg / l
!       Nitrate                         microM           nitrate              microM or mg / l
!       Silicate                        microM        silicate acid             microM or mg / l
!
!   MOHID names based on ModuleGlobalData
!   P,N,Si have 14 levels, 0-500m
!   O2 % has 24 levels, 0-1500m
!   Standart Levels: 0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1750,
!                    2000,2500,3000,3500,4000,4500,5000,5500
!   Conversion eq:
!   microM X = 10^-6 moles X/L == mmol X / m3
!   mg X / L = microM X * 10^-3 * PM mg / L
!
!   For clarifications,check: woa05_oxygen_final.pdf, woa05_nutrients_final.pdf on the kepler WOA05 data folder

Module ModuleWOAFormat

    use netcdf90
    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid

    implicit none

    private
    
    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertWOAFormat
    private ::      ReadOptions
    private ::      OpenReadAndConvertFile
    private ::          OpenASCIIFile
    private ::          OpenNetCDFFile
    private ::          Open_HDF5_OutPut_File
    private ::          ConstructGridAndMap
    private ::              ReadGridMapOptions
    private ::              ConstructGrid
    private ::                  GetVariable
    private ::              ConstructBathymetry
    private ::              ConstructTime
    private ::              ConstructMap
    private ::          HandleAttributtes
    private ::          IsGridVariable
    private ::          WriteDataType
    private ::          HandleDimensions
    private ::              WriteDimensions
    private ::          HandleChars
    private ::          HandleShorts
    private ::              ConvertToMohidUnits3D
    private ::          HandleIntegers
    private ::          HandleFloats
    private ::          HandleDoubles
    private ::      KillWOAFormat   

    !ERA40-----------------------------------------------------------------------

    !Parameters----------------------------------------------------------------
    integer, parameter                          :: DotGrid                  = 1 !Cell corners
    integer, parameter                          :: CrossGrid                = 2 !Cell center

    !Parameters-------------------------------------------------------------------------
    integer, parameter                          :: Cartesian                = 1
    integer, parameter                          :: Spherical                = 2
    integer, parameter                          :: Curvilinear              = 3
    integer, parameter                          :: Spherical_Curvilinear    = 4
    
    integer, parameter                          :: Sigma                    = 1
    integer, parameter                          :: Z_Coord                  = 2
    integer, parameter                          :: Generic                  = 3

    !Types---------------------------------------------------------------------
    type       T_Date
        type(T_Time)                                        :: Date
        type(T_Date), pointer                               :: Next
    end type  T_Date

    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        logical                                             :: Convert        = .false.
        integer                                             :: nDimensions
        integer                                             :: GridLocation
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        integer                                             :: OutputNumber   = 1
        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),              pointer                 :: Next
    end type  T_Field

    type     T_Value
        character(len=StringLength)                         :: Units
        !byte
        integer(1)                              :: Byte_
        !char
        character(len=80)                       :: Char     = char(32)
        character(len=80), dimension(:), pointer:: Char_1D
        !short
        integer                                 :: Short
        integer, dimension(:    ), pointer      :: Short_1D
        integer, dimension(:,:  ), pointer      :: Short_2D
        integer, dimension(:,:,:), pointer      :: Short_3D
        !int       
        integer                                 :: Int
        integer, dimension(:      ), pointer    :: Int_1D
        integer, dimension(:,:    ), pointer    :: Int_2D
        integer, dimension(:,:,:  ), pointer    :: Int_3D
        integer, dimension(:,:,:,:), pointer    :: Int_4D

        !float
        real                                    :: Float
        real,    dimension(:      ), pointer    :: FLOAT_1D
        real,    dimension(:,:    ), pointer    :: FLOAT_2D
        real,    dimension(:,:,:  ), pointer    :: FLOAT_3D
        real,    dimension(:,:,:,:), pointer    :: FLOAT_4D

        !double
        real(8)                                 :: Double
        real(8), dimension(:      ), pointer    :: Double_1D
        real(8), dimension(:,:    ), pointer    :: Double_2D
        real(8), dimension(:,:,:  ), pointer    :: Double_3D
        real(8), dimension(:,:,:,:), pointer    :: Double_4D

    end type T_Value

    type     T_ID
        character(len=80)                       :: Name
        integer                                 :: Number
        integer                                 :: Length
        integer                                 :: Type_
    end type T_ID

    type     T_Dimension
        character(len=80)                       :: Name
        type(T_Size3D)                          :: Size
    end type T_Dimension

    type     T_Attribute
        type(T_ID)                              :: ID
        type(T_Value)                           :: Value
    end type T_Attribute

    type     T_Variable
        type(T_ID)                              :: ID
        type(T_Attribute)                       :: Attribute
        type(T_Value)                           :: Value
        type(T_Dimension)                       :: Dim
        integer                                 :: nDimensions
        integer                                 :: nAttributes
        integer, dimension(4)                   :: DimensionsID
    end type T_Variable
    
    type       T_WOAFormat
        integer                                             :: ObjEnterData             = 0
        integer                                             :: ObjHDF5                  = 0
        integer                                             :: ObjHorizontalGrid        = 0
        integer                                             :: Unit
        integer                                             :: OutputUnit
        !implement 4 files
        !character(len=StringLength), dimension(:), pointer  :: InputFiles
        character(len=PathLength)                           :: InputFileName
        character(len=PathLength)                           :: InputFilePath
        character(len=PathLength)                           :: InputFile
        character(len=PathLength)                           :: InGridFileName
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: OutputFileName
        character(len=PathLength)                           :: HDF5FileName
        character(len=PathLength)                           :: PropName
        character(len=PathLength)                           :: PropUnits
        integer                                             :: Grid_Type
        real, dimension(:,:  ),       pointer               :: Bathymetry
        real, dimension(:,:  ),       pointer               :: CenterX
        real, dimension(:,:  ),       pointer               :: CenterY
        real, dimension(:,:  ),       pointer               :: ConnectionX
        real, dimension(:,:  ),       pointer               :: ConnectionY
        real, dimension(:    ),       pointer               :: Sigma
        real, dimension(:,:,:),       pointer               :: SZZ
        integer, dimension(:,:,:),       pointer               :: WaterPoints3D
        real, dimension(:),pointer                          :: Levels
        logical                                             :: InGridFileNotSpecified  = .false.
        logical                                             :: OutputGridNotSpecified  = .false.
        logical                                             :: ConvertToHDF5
        logical                                             :: ConvertToASCII
        logical                                             :: MOHID_UNITS

        real,        dimension(:    ), pointer              :: TimeValues
        integer                                             :: Time_LB, Time_UB
        integer                                             :: nDimensions, nVariables, nAttributes
        integer                                             :: HDFOutputNumber
        integer                                             :: NumDates = FillValueInt
        type(T_Time), dimension(:), pointer                 :: OutPutDates
        type(T_Time)                                        :: InitialDate
        type(T_Time)                                        :: CurrentDate
        type(T_Time)                                        :: PreviousDate
        type(T_Time),dimension(:    ), pointer              :: TimeArray
        integer, dimension(:,:  ), pointer                  :: Mapping2D
        character(len=80)                                   :: XX_Variable
        character(len=80)                                   :: YY_Variable
        character(len=80)                                   :: ZZ_Variable
        character(len=80)                                   :: Time_Variable
        character(len=80)                                   :: GridType_Variable
        character(len=80)                                   :: VertCoord_Variable
        character(len=80)                                   :: Batim_Variable
        type(T_Size3D)                                      :: Size
        type(T_Size3D)                                      :: WorkSize
        type(T_Variable), pointer                           :: Batim
        real                                                :: FillValue, Add_OffSet
        real                                                :: Missing_Value, Add_ScaleFactor
        integer                                             :: FileID
        integer                                             :: DoneGrid
        !Maybe the next are not needed:
        real                                                :: PTop                   = null_real
        type(T_Field),                pointer               :: FirstField
        type(T_Date),                 pointer               :: FirstDate
        character(len=StringLength), dimension(:), pointer  :: FieldsToConvert
        character(len=StringLength), dimension(:), pointer  :: FieldsToRead
        integer                                             :: NFiles = FillValueInt
        character(len=StringLength), dimension(:), pointer  :: Files

    end type  T_WOAFormat

    type(T_WOAFormat), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertWOAFormat(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions(ClientNumber)

        call OpenReadAndConvertFile

        call KillWOAFormat

        STAT = SUCCESS_


    end subroutine ConvertWOAFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions(ClientNumber)
        
        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag, i
        logical                                         :: BlockFound
        integer                                         :: FirstLine, LastLine

        !Begin-----------------------------------------------------------------
       
        call GetData(Me%InputFilePath,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'FILEPATH',                     &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWOAFormat - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleWOAFormat - ERR02'
        end if

        call GetData(Me%OutputFileName,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'OUTPUTFILENAME',               &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWOAFormat - ERR04'      

        call GetData(Me%ConvertToHDF5,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'CONVERT_TO_HDF5',              &
                     Default      = OFF,                            &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWOAFormat - ERR05'

        call GetData(Me%ConvertToASCII,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'CONVERT_TO_ASCII',             &
                     Default      = OFF,                            &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWOAFormat - ERR06'

        call GetData(Me%MOHID_UNITS,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'MOHID_UNITS',                  &
                     Default      = OFF,                            &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWOAFormat - ERR06'



!get file names
do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                    &
                                        '<<begin_woafiles>>', '<<end_woafiles>>', BlockFound,    &
                                        FirstLine = FirstLine, LastLine = LastLine,      &
                                        STAT = STAT_CALL)

if1 :       if(STAT_CALL .EQ. SUCCESS_) then    
if2 :           if (BlockFound) then

                    Me%nFiles = LastLine - FirstLine - 1

                    allocate (Me%Files(Me%nFiles))

                    do i = 1, Me%nFiles

                        call GetData(Me%Files(i), Me%ObjEnterData,  iflag,          & 
                                     Buffer_Line  = FirstLine + i,                       &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleWOAFormat - ERR02'

                    enddo

                    exit do1 
                        
                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'ReadOptions - ModuleWOAFormat - ERR03'
                    
            end if if1
        end do do1




    end subroutine ReadOptions
   
    !------------------------------------------------------------------------

     subroutine ConvertToMohidUnits3D(Variable)

        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable

        !Local-----------------------------------------------------------------
                
        !Begin-----------------------------------------------------------------

!data - Variable%Value%Double_4D
!mg X / L = microM X * 10-3 * PM X
!PM: peso molecular, P = 30.973, N = 14.0067, Si = 28.085 g/mol
        
        select case(trim(Variable%ID%Name))

            case('silicate acid')

                Variable%Value%Units     = 'mgSi/l'
                Variable%Value%Double_4D  = Variable%Value%Double_4D * 1E-3 * 28.085
                write(*,*)'     converting '//trim(Variable%ID%Name)//" units to "//trim(Variable%Value%Units)
            case('nitrate')

                Variable%Value%Units     = 'mgN/l'
               Variable%Value%Double_4D  = Variable%Value%Double_4D * 1E-3 * 14.0067
                write(*,*)'     converting '//trim(Variable%ID%Name)//" units to "//trim(Variable%Value%Units)
            case('inorganic phosphorus')

                Variable%Value%Units     = 'mgP/l'
               Variable%Value%Double_4D  = Variable%Value%Double_4D * 1E-3 * 30.973
                write(*,*)'     converting '//trim(Variable%ID%Name)//" units to "//trim(Variable%Value%Units)
            case('dissolved oxygen percent saturation')

                Variable%Value%Units     = '%'
               write(*,*)'      no conversion '
            case default

        end select

    end subroutine ConvertToMohidUnits3D

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%HDF5FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleWOAFormat - ERR01'
        

    end subroutine Open_HDF5_OutPut_File


    !----------------------------------------------------------------------
   
    subroutine KillWOAFormat
        
        !Local-----------------------------------------------------------------
        integer                                         :: nUsers
        
        !Begin-----------------------------------------------------------------
       


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillWOAFormat - ModuleWOAFormat - ERR03'

        deallocate(Me%Mapping2D)
        nullify   (Me%Mapping2D)
        deallocate(Me%ConnectionX, Me%ConnectionY)
        nullify   (Me%ConnectionX, Me%ConnectionY)
        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)


    
    end subroutine KillWOAFormat

    !--------------------------------------------------------------------------

    subroutine OpenReadAndConvertFile

        !Local-----------------------------------------------------------------
        integer                                         :: stat
        integer                                         :: VariableID, STAT_CALL,i
        type(T_Variable), pointer                       :: Variable
        character(len=PathLength)                       :: AuxFileName
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)"Converting netcdf file. Please wait..."
        write(*,*)
        Me%DoneGrid = 0
        

        Me%HDF5FileName = Me%OutputFileName

        call Open_HDF5_OutPut_File

        Do i=1, Me%NFiles

            AuxFileName = trim(Me%InputFilePath)//trim(Me%Files(i))
            !InputFile = Me%InputFilePath//"i0112an1.nc"
            call OpenNetCDFFile(AuxFileName)


            if(Me%nAttributes > 0)then

                !Handle global attributes
                call HandleAttributtes(NF90_GLOBAL, Me%nAttributes)

            endif


            if(Me%DoneGrid == 0)then ! do grid thingies
    
                    call ConstructGridAndMap
                    !Set grid
                    call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,  &
                                       Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR02'

                    !ConnectionX is not really a connection X but for could be useful
                    !if HDF5Extractor use is required (here is equal to Longitude)
                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX",           &
                                          "degrees",                                    &
                                          Array2D = Me%ConnectionX,                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR03a'

                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Longitude",             &
                                          "degrees",                                    &
                                          Array2D = Me%ConnectionX,                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR03b'

                    !ConnectionY is not really a connection Y but for could be useful
                    !if HDF5Extractor use is required (here is equal to Latitude)
                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY",           & 
                                          "degrees", Array2D = Me%ConnectionY,          &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR04a'

                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Latitude",              &
                                          "degrees", Array2D = Me%ConnectionY,          &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR04b'

                    !Set bathymetry          
                    call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB,Me%WorkSize%IUB,     &
                               Me%WorkSize%JLB,Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR05'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                               "/Grid",                                                 &
                               "Bathymetry",                                            &
                               'm',                                                     &
                               Array2D      = Me%Batim%Value%Float_2D,                  &
                               STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR06'

                    !Set map
                    call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,    & 
                                       Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR07'
                    
                    call HDF5FlushMemory(Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR10'


                endif

                !load vars into hdf5

            do VariableID = 1, Me%nVariables

                nullify(Variable); allocate(Variable)

                Variable%ID%Number = VariableID

                stat = nf90_inquire_variable(Me%FileID, Variable%ID%Number, Variable%ID%Name,   &
                                             Variable%ID%Type_, Variable%nDimensions,           &
                                             Variable%DimensionsID, Variable%nAttributes)
                if (stat /= nf90_noerr) stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR01'


                if (Variable%ID%Name == "i0112an1") Variable%ID%Name = "silicate acid"
                if (Variable%ID%Name == "O0112an1") Variable%ID%Name = "dissolved oxygen percent saturation"
                if (Variable%ID%Name == "n0112an1") Variable%ID%Name = "nitrate"
                if (Variable%ID%Name == "p0112an1") Variable%ID%Name = "inorganic phosphorus"
                

          
                call HandleDimensions (Variable)

                call HandleAttributtes(Variable%ID%Number,  Variable%nAttributes)

                !if(Me%ConvertToHDF5)then
                if( ((Variable%ID%Name == "silicate acid").or.(Variable%ID%Name == "nitrate") &
                     .or.(Variable%ID%Name == "dissolved oxygen percent saturation").or.      &
                     (Variable%ID%Name == "inorganic phosphorus")))then
                    write(*,*)trim(Variable%ID%Name)
                    select case(Variable%ID%Type_)

                        case(NF90_BYTE)
                
                            write(*,*)'Cannot convert information in bytes.'

                        case(NF90_CHAR)

                            call HandleChars(Variable)


                        case(NF90_SHORT)

                            call HandleShorts(Variable)


                        case(NF90_INT)

                            call HandleIntegers(Variable)

                        case(NF90_FLOAT)

                            call HandleFloats(Variable)


                        case(NF90_DOUBLE)
                
                            call HandleDoubles(Variable)

                    end select

                end if

                nullify(Variable); deallocate(Variable)

        
            end do !var cycle
        end do !file cycle
        call HDF5FlushMemory(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR09'

        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
        stop 'OpenReadAndConvertFile - ModuleWOAFormat - ERR11'

        write(*,*)'Finished converting'


    end subroutine OpenReadAndConvertFile

    !------------------------------------------------------------------
    
    subroutine OpenASCIIFile

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenASCIIFile - ModuleWOAFormat - ERR01'

        open(Unit   = Me%OutputUnit,                                            &
             File   = Me%OutputFileName,                                        &
             STATUS = 'UNKNOWN',                                                &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenASCIIFile - ModuleWOAFormat - ERR02'

    end subroutine OpenASCIIFile

    !------------------------------------------------------------------

    logical function IsGridVariable(Variable)
        
        !Arguments-----------------------------------------------------          
        type(T_Variable), pointer                       :: Variable
        
        !Local---------------------------------------------------------          
        character(len=80)                               :: Name

        !Begin---------------------------------------------------------          


        Name = trim(Variable%ID%Name)

        IsGridVariable = .false.

        if    (Name == Me%XX_Variable)then

            IsGridVariable = .true.

        elseif(Name == Me%YY_Variable)then

            IsGridVariable = .true.

        elseif(Name == Me%ZZ_Variable)then

            IsGridVariable = .true.

        elseif(Name == Me%Time_Variable)then

            IsGridVariable = .true.

        elseif(Name == Me%GridType_Variable)then

            IsGridVariable = .true.

        elseif(Name == Me%VertCoord_Variable)then

            IsGridVariable = .true.

        elseif(Name == Me%Batim_Variable)then

            IsGridVariable = .true.

        end if
    
    
    end function IsGridVariable

    !------------------------------------------------------------------
    
    subroutine OpenNetCDFFile(InputFile)

        !Arguments-------------------------------------------------------------
        character(len=PathLength)                           :: InputFile
        !full path to file
        !Local-----------------------------------------------------------------
        integer                                         :: stat
        logical                                         :: exist

        !----------------------------------------------------------------------
        
        
        !Verifies if file exists
        inquire(file = InputFile, exist = exist)
        if (.not. exist) then
            write(*,*)'NETCDF file does not exist'
            stop 'OpenNetCDFFile - ModuleWOAFormat - ERR01'
        endif

        stat = nf90_open(InputFile, NF90_NOWRITE, Me%FileID)
        if (stat /= nf90_noerr) stop 'OpenNetCDFFile - ModuleWOAFormat - ERR02'

        stat = nf90_inquire(Me%FileID, Me%nDimensions, Me%nVariables, Me%nAttributes)
        if (stat /= nf90_noerr) stop 'OpenNetCDFFile - ModuleWOAFormat - ERR03'

    end subroutine OpenNetCDFFile

    !------------------------------------------------------------------
  subroutine ConstructGridAndMap

        call ReadGridMapOptions
        
        call ConstructGrid

        call ConstructBathymetry

        call ConstructTime

        call ConstructMap

    end subroutine ConstructGridAndMap


    !------------------------------------------------------------------
    subroutine ConstructVerticalCoordinate(Variable,t)

        !Local---------------------------------------------------------          
!        type(T_Size3D)                                  :: Size
         integer                                         :: i, j, k, t
         integer                                         :: STAT_CALL
         type(T_Variable), pointer                       :: Variable
        !Begin---------------------------------------------------------          
        
        allocate(Me%SZZ(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB,0:24))

allocate(Me%Levels(0:25))
        Me%Levels(0) = 1600.
        Me%Levels(1) = 1500.
        Me%Levels(2) = 1400.
        Me%Levels(3) = 1300.
        Me%Levels(4) = 1200.
        Me%Levels(5) = 1100.
        Me%Levels(6) = 1000.
        Me%Levels(7) = 900.
        Me%Levels(8) = 800.
        Me%Levels(9) = 700.
        Me%Levels(10) = 600.
        Me%Levels(11) = 500.
        Me%Levels(12) = 400.
        Me%Levels(13) = 300.
        Me%Levels(14) = 250.
        Me%Levels(15) = 200.
        Me%Levels(16) = 150.
        Me%Levels(17) = 125.
        Me%Levels(18) = 100.
        Me%Levels(19) = 75.
        Me%Levels(20) = 50.
        Me%Levels(21) = 30.
        Me%Levels(22) = 20.
        Me%Levels(23) = 10.
        Me%Levels(24) = 0.
        Me%Levels(25) = -10.
        !reverse batims

        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        !do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        !do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (trim(Variable%ID%Name) == "dissolved oxygen percent saturation")then
                do k = 0,24
                    Me%SZZ(i,j,k) = Me%Levels(k) + (Me%Levels(k+1) - Me%Levels(k))/2
                    Me%SZZ(i,j,24) = -5.
                enddo
                call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                Me%WorkSize%JLB, Me%WorkSize%JUB, &
                0, 24, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ConstructVerticalCoordinate - GETMToHDF5 - ERR07'

            else

                do k = 10,23
                    Me%SZZ(i,j,k-10) = Me%Levels(k) + (Me%Levels(k+1) - Me%Levels(k))/2
                    Me%SZZ(i,j,14) = -5.
                enddo

                call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                Me%WorkSize%JLB, Me%WorkSize%JUB, &
                0, 14, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ConstructVerticalCoordinate - GETMToHDF5 - ERR07'

            end if
            !Me%SZZ(i,j,0:14) = Me%Levels(10:24)

          !  Me%SZZ(i,j,Me%WorkSize%KLB-1              ) = MaxDepth

        enddo
        enddo
        !writes stuff from !woa
        

        call HDF5WriteData(Me%ObjHDF5,                                  &
                           "/Grid/VerticalZ",                               &
                           "Vertical",                                     &
                           "m",                                             &
                           Array3D      = Me%SZZ,                           &
                           OutputNumber = t,                                &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructVerticalCoordinate - GETMToHDF5 - ERR08'
     



    end subroutine ConstructVerticalCoordinate
    !------------------------------------------------------------------

    subroutine HandleAttributtes(TypeOfAttribute, nAttributes, WriteToFile)

        !Arguments-------------------------------------------------------------
        integer                                         :: TypeOfAttribute, nAttributes
        logical, optional                               :: WriteToFile

        !Local-----------------------------------------------------------------
        type(T_Attribute), pointer                      :: Attribute
        integer                                         :: stat, n
        logical                                         :: WriteToFile_

        !Begin-----------------------------------------------------------------

        if(present(WriteToFile))then

            WriteToFile_ = WriteToFile

        else

            WriteToFile_ = .true.

        end if

        do n = 1, nAttributes

            nullify(Attribute); allocate(Attribute)

            stat = nf90_inq_attname(Me%FileID, TypeOfAttribute, n, Attribute%ID%Name)
            if (stat /= nf90_noerr) stop 'HandleAttributtes - ModuleWOAFormat - ERR01'

            stat = nf90_inquire_attribute(Me%FileID, TypeOfAttribute, Attribute%ID%Name, &
                                          Attribute%ID%Type_,                            &
                                          Attribute%ID%Length, Attribute%ID%Number)            
            if (stat /= nf90_noerr) stop 'HandleAttributtes - ModuleWOAFormat - ERR02'
            
            !Only requires the pre-existence of the ASCII file, if it is the case
            if(WriteToFile_ .and. Me%ConvertToASCII)then

                call WriteDataLine(Me%OutputUnit, "<begin_attribute>")
                call WriteDataLine(Me%OutputUnit, "Name", trim(Attribute%ID%Name))
                call WriteDataType(Attribute%ID%Type_)

            endif


            select case(Attribute%ID%Type_)

                case(NF90_BYTE)

                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Byte_)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleWOAFormat - ERR03'

                case(NF90_CHAR)

                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Char)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleWOAFormat - ERR04'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        call WriteDataLine(Me%OutputUnit,"Value", Attribute%Value%Char)
                   

                case(NF90_SHORT)
                    
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Short)
                    if (stat /= nf90_noerr)                                              &
                        stop 'HandleAttributtes - ModuleWOAFormat - ERR05'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        write(Me%OutputUnit,*)"Value            :", Attribute%Value%Short

                case(NF90_INT)
                   
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Int)
                    if (stat /= nf90_noerr)                                              &
                        stop 'HandleAttributtes - ModuleWOAFormat - ERR06'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             &
                        call WriteDataLine(Me%OutputUnit, "Value", Attribute%Value%Int)


                case(NF90_FLOAT)
                   
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Float)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleWOAFormat - ERR07'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        call WriteDataLine(Me%OutputUnit,"Value", Attribute%Value%Float)


                case(NF90_DOUBLE)
                    
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Double)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleWOAFormat - ERR08'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        write(Me%OutputUnit,*)"Value            :", Attribute%Value%Double

            end select

            if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<end_attribute>")

            if (trim(Attribute%ID%Name)=='add_offset'  )                                 & 
                Me%Add_OffSet      = Attribute%Value%Double

            if (trim(Attribute%ID%Name)=='scale_factor')                                 & 
                Me%Add_ScaleFactor = Attribute%Value%Double

            if (trim(Attribute%ID%Name)=='missing_value') then
                Me%Missing_Value = Me%Add_ScaleFactor*real(Attribute%Value%Short) + Me%Add_OffSet
            endif

            if ((trim(Attribute%ID%Name)=='units') .AND.                                 &
                (TypeOfAttribute > 3)) Me%PropUnits = trim(Attribute%Value%Char)

            if (trim(Attribute%ID%Name)=='long_name'    ) then
                Me%PropName        = trim(Attribute%Value%Char(1:33))

           endif

            deallocate(Attribute)

        end do

    end subroutine HandleAttributtes

    !------------------------------------------------------------------

    subroutine WriteDataType(DataType)
        
        !Arguments-------------------------------------------------------------
        integer                                         :: DataType
        
        !Begin-----------------------------------------------------------------
         
        if(Me%ConvertToASCII)then
          
            select case(DataType)

                case(NF90_BYTE)
                
                    call WriteDataLine(Me%OutputUnit,"Data type ID", NF90_BYTE)
                    call WriteDataLine(Me%OutputUnit,"Data type char", "byte")

                case(NF90_CHAR)

                    call WriteDataLine(Me%OutputUnit,"Data type ID", NF90_CHAR)
                    call WriteDataLine(Me%OutputUnit,"Data type char", "char")


                case(NF90_SHORT)

                    call WriteDataLine(Me%OutputUnit,"Data type ID", NF90_SHORT)
                    call WriteDataLine(Me%OutputUnit,"Data type char", "short")

                case(NF90_INT)
                    call WriteDataLine(Me%OutputUnit,"Data type ID", NF90_INT)
                    call WriteDataLine(Me%OutputUnit,"Data type char", "integer")


                case(NF90_FLOAT)
               
                    call WriteDataLine(Me%OutputUnit,"Data type ID", NF90_FLOAT)
                    call WriteDataLine(Me%OutputUnit,"Data type char", "float")

                case(NF90_DOUBLE)
                
                    call WriteDataLine(Me%OutputUnit,"Data type ID", NF90_DOUBLE)
                    call WriteDataLine(Me%OutputUnit,"Data type char", "double")

            end select

        end if

    end subroutine WriteDataType

    !--------------------------------------------------------------------------

    subroutine HandleDimensions(Variable, WriteDims)
        
        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable
        type(T_Size3D)                                  :: Size
        logical, optional                               :: WriteDims

        !Local-----------------------------------------------------------------
        integer                                         :: dim, i
        logical                                         :: WriteDimensions_

        !Begin-----------------------------------------------------------------

        Size%ILB = 1;        Size%JLB = 1;        Size%KLB = 1
        Size%IUB = null_int; Size%JUB = null_int; Size%KUB = null_int

        if(present(WriteDims))then
            WriteDimensions_ = WriteDims
        else
            WriteDimensions_ = .true.
        end if

        do i = 1, Variable%nDimensions

            call WriteDimensions(i, Variable, dim, WriteDimensions_)

            select case(i)

                case(1)

                    Size%IUB = dim

                case(2)

                    Size%JUB = dim

                case(3)

                    Size%KUB = dim

            end select

        end do

        Variable%Dim%Size = Size

    end subroutine HandleDimensions

    !--------------------------------------------------------------------------

    subroutine HandleFloats(Variable)
        
        !Arguments-------------------------------------------------------------
        type(T_Variable),   pointer                     :: Variable

        !Local-----------------------------------------------------------------
        integer                                         :: stat, i, j, k, t, STAT_CALL
        type(T_Size3D)                                  :: Size

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size

        select case(Variable%nDimensions)
            
            case(0)
                
                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Float)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleWOAFormat - ERR01'

                if(Me%ConvertToASCII)then
                    call WriteDataLine(Me%OutputUnit,"<begin_data>")
                    write(Me%OutputUnit,*)Variable%Value%Float
                    call WriteDataLine(Me%OutputUnit,"<end_data>")
                endif


            case(1)

                nullify (Variable%Value%FLOAT_1D)
                allocate(Variable%Value%FLOAT_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_1D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleWOAFormat - ERR02'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                        write(Me%OutputUnit,*)Variable%Value%FLOAT_1D(i)
                    enddo
                
                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                end if

                deallocate(Variable%Value%FLOAT_1D)

            case(2)

                nullify (Variable%Value%FLOAT_2D)
                allocate(Variable%Value%FLOAT_2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_2D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleWOAFormat - ERR03'
                
                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")
                
                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                        write(Me%OutputUnit,*)Variable%Value%FLOAT_2D(i,j)
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                if(Me%ConvertToHDF5)then

!                    call ConvertToMohidUnits2D(Field)
                    
                    call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,        & 
                                       Size%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleFloats - ModuleWOAFormat - ERR04'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array2D      = Variable%Value%FLOAT_2D,          &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleFloats - ModuleWOAFormat - ERR05'

                end if

                deallocate(Variable%Value%FLOAT_2D)

            
            case(3)

                nullify (Variable%Value%FLOAT_3D)
                allocate(Variable%Value%FLOAT_3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,  &
                         Size%KLB:Size%KUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_3D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleWOAFormat - ERR06'


                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                    do k = Size%KLB, Size%KUB
                        write(Me%OutputUnit,*)Variable%Value%FLOAT_3D(i,j,k)
                    enddo
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


                if(Me%ConvertToHDF5)then

                    call ConvertToMohidUnits3D(Variable)

                    call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,        & 
                                       Size%JUB, Size%KLB, Size%KUB, STAT = STAT_CALL)
                                       
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleFloats - ModuleWOAFormat - ERR07'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array3D      = Variable%Value%FLOAT_3D,          &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleFloats - ModuleWOAFormat - ERR08'
                
                end if

                deallocate(Variable%Value%FLOAT_3D)


            case(4)

                nullify (Variable%Value%FLOAT_4D)
                allocate(Variable%Value%FLOAT_4D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,  &
                                                 Size%KLB:Size%KUB, Me%Time_LB:Me%Time_UB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_4D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleWOAFormat - ERR09'


                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<begin_data>")

                do t = Me%Time_LB, Me%Time_UB
                    
                    if(Me%ConvertToASCII) then
                    
                        call WriteDataLine(Me%OutputUnit,"<begin_field>")

                        do i = Size%ILB, Size%IUB
                        do j = Size%JLB, Size%JUB
                        do k = Size%KLB, Size%KUB
                            write(Me%OutputUnit,*)Variable%Value%FLOAT_4D(i,j,k,t)
                        enddo
                        enddo
                        enddo

                        call WriteDataLine(Me%OutputUnit,"<end_field>")

                    end if

                    if(Me%ConvertToHDF5)then


                        call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,    & 
                                           Size%JUB, Size%KLB, Size%KUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleFloats - ModuleWOAFormat - ERR10'

                        Variable%Value%FLOAT_3D => Variable%Value%FLOAT_4D(:,:,:,t)

                        call HDF5WriteData(Me%ObjHDF5,                                  &
                                           "/Results/"//trim(Variable%ID%Name),         &
                                           trim(Variable%ID%Name),                      &
                                           Variable%Value%Units,                        &
                                           Array3D      = Variable%Value%FLOAT_3D,      &
                                           OutputNumber = t,                            &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleFloats - ModuleWOAFormat - ERR10'
                        
                        nullify(Variable%Value%FLOAT_3D)

                    end if

                    

                end do
                
                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<end_data>")

                deallocate(Variable%Value%FLOAT_4D)

            case default


        end select


    end subroutine HandleFloats


    !-------------------------------------------------------------------------


    subroutine HandleDoubles(Variable)
        
        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable
       !Local-----------------------------------------------------------------
        integer                                         :: stat, i, j, k, t
        integer                                         :: STAT_CALL
        type(T_Size3D)                                  :: Size
        real, dimension(:,:,:), pointer                 :: Temp, aux
        real(8), dimension(:,:  ), pointer             :: temp4

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size


        select case(Variable%nDimensions)
            
            case(4)

                nullify (Variable%Value%Double_4D)
                allocate(Variable%Value%Double_4D(Size%ILB:Size%IUB, Size%JLB:Size%JUB, &
                                                  Size%KLB:Size%KUB,                    & 
                                                  Me%Time_LB:Me%Time_UB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number,                      & 
                                    Variable%Value%Double_4D)
                if (stat /= nf90_noerr) stop 'HandleDoubles - ModuleWOAFormat - ERR09'

                !convert to mohid units
                if (Me%MOHID_UNITS) call ConvertToMohidUnits3D(Variable)

                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<begin_data>")



                do t = Me%Time_LB, Me%Time_UB
                    
                    if(Me%ConvertToASCII)then


                    
                        call WriteDataLine(Me%OutputUnit,"<begin_field>")

                        do i = Size%ILB, Size%IUB
                        do j = Size%JLB, Size%JUB
                        do k = Size%KLB, Size%KUB
                            write(Me%OutputUnit,*)Variable%Value%Double_4D(i,j,k,t)
                        enddo
                        enddo
                        enddo

                    endif


                    call HDF5SetLimits(Me%ObjHDF5, Size%JLB, Size%JUB, Size%ILB,    & 
                                       Size%IUB, Size%KLB, Size%KUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'HandleDoubles - ModuleWOAFormat - ERR10'
            
            
                    allocate(Temp(Size%JLB:Size%JUB, Size%ILB:Size%IUB, Size%KLB:Size%KUB))
                    allocate(aux(Size%JLB:Size%JUB, Size%ILB:Size%IUB, Size%KLB:Size%KUB))


                    
                    Variable%Value%Double_3D => Variable%Value%Double_4D(:,:,:,t)
                    

                    do k = Size%KLB, Size%KUB
                        temp4 => Variable%Value%Double_3D(:,:,k)

                        Temp(:,:,-(k-Size%KUB)+1) = transpose(temp4)
                    enddo
                  
                  
                     
                    aux(:,:,:) =  temp(:,:,:)

                    temp(:,   1:180, :) = aux(:, 181:360, :)
                    temp(:, 181:360, :) = aux(:,   1:180, :)

                    deallocate(aux); nullify(aux)
                    
                    

                    !woa
                    call HDF5WriteData(Me%ObjHDF5,                                  &
                                       "/Results/"//trim(Variable%ID%Name),         &
                                       trim(Variable%ID%Name),                      &
                                       Variable%Value%Units,                        &
                                       Array3D      = Temp,                         &
                                       OutputNumber = t,                            &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'HandleDoubles - ModuleWOAFormat - ERR11'
                    

!write verticalz
                    if (Me%DoneGrid == 0) then

                        call ConstructVerticalCoordinate(Variable,t)

!writes waterpoints3d only one time
                         if (t==1)then 

                            nullify(Me%WaterPoints3D)
                            allocate(Me%WaterPoints3D(Size%JLB:Size%JUB, Size%ILB:Size%IUB, Size%KLB:Size%KUB))
                            
                            
                            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            do j = Me%WorkSize%JLB, Me%WorkSize%JUB 
                            do k = Me%WorkSize%KLB, Me%WorkSize%KUB 
                                    if(Temp(i,j,k) > 1E+10)then
                                        Me%WaterPoints3D(i,j,k) = 0
                                    else
                                        Me%WaterPoints3D(i,j,k) = 1
                                    endif
                            enddo
                            enddo
                            enddo
                            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, &
                            Me%WorkSize%JLB, Me%WorkSize%JUB, &
                            1, 14, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'ConstructVerticalCoordinate - GETMToHDF5 - ERR07'


                            call HDF5WriteData(Me%ObjHDF5,                                  &
                                                "/Grid",                           &
                                                "WaterPoints3D",                                      &
                                                "m",                                             &
                                                Array3D      = Me%WaterPoints3D,                           &
                                                STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'ConstructVerticalCoordinate - GETMToHDF5 - ERR08'
                        endif
                   
                    end if !enddonegrid if
                    
                    nullify(Variable%Value%Double_3D)
                
                
                end do !end t do

                Me%DoneGrid = 1
                
                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<end_data>")

                deallocate(Variable%Value%Double_4D)
                deallocate(Temp); nullify(temp)

            case default


        end select


    end subroutine HandleDoubles

    !--------------------------------------------------------------------------

    subroutine HandleIntegers(Variable)
        
        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable

        !Local-----------------------------------------------------------------
        integer                                         :: stat, i, j, k, t
        integer                                         :: STAT_CALL
        type(T_Size3D)                                  :: Size

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size

        select case(Variable%nDimensions)
            
            case(0)
                
                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleWOAFormat - ERR01'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Int

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

            case(1)

                nullify (Variable%Value%Int_1D)
                allocate(Variable%Value%Int_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_1D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleWOAFormat - ERR02'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                        write(Me%OutputUnit,*)Variable%Value%Int_1D(i)
                    enddo
                
                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                deallocate(Variable%Value%Int_1D)

            case(2)

                nullify (Variable%Value%Int_2D)
                allocate(Variable%Value%Int_2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_2D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleWOAFormat - ERR03'
                
                if(Me%ConvertToASCII)then
                
                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                        write(Me%OutputUnit,*)Variable%Value%Int_2D(i,j)
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")
                
                endif

                if(Me%ConvertToHDF5)then

                    call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,        & 
                                       Size%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleIntegers - ModuleWOAFormat - ERR04'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array2D      = Variable%Value%Int_2D,            &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleIntegers - ModuleWOAFormat - ERR05'

                end if

                deallocate(Variable%Value%Int_2D)

            
            case(3)

                nullify (Variable%Value%Int_3D)
                allocate(Variable%Value%Int_3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,    & 
                         Size%KLB:Size%KUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_3D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleWOAFormat - ERR06'


                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                    do k = Size%KLB, Size%KUB
                        write(Me%OutputUnit,*)Variable%Value%Int_3D(i,j,k)
                    enddo
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                if(Me%ConvertToHDF5)then

                    call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,        & 
                                       Size%JUB, Size%KLB, Size%KUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleIntegers - ModuleWOAFormat - ERR07'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array3D      = Variable%Value%Int_3D,            &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleIntegers - ModuleWOAFormat - ERR08'
                
                end if


                deallocate(Variable%Value%Int_3D)

            case(4)

                nullify (Variable%Value%Int_4D)
                allocate(Variable%Value%Int_4D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,    &
                                               Size%KLB:Size%KUB, Me%Time_LB:Me%Time_UB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_4D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleWOAFormat - ERR09'


                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<begin_data>")

                do t = Me%Time_LB, Me%Time_UB
                    
                    if(Me%ConvertToASCII)then
                        
                        call WriteDataLine(Me%OutputUnit,"<begin_field>")

                        do i = Size%ILB, Size%IUB
                        do j = Size%JLB, Size%JUB
                        do k = Size%KLB, Size%KUB
                            write(Me%OutputUnit,*)Variable%Value%Int_4D(i,j,k,t)
                        enddo
                        enddo
                        enddo

                    end if

                    if(Me%ConvertToHDF5)then


                        call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,    & 
                                           Size%JUB, Size%KLB, Size%KUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleIntegers - ModuleWOAFormat - ERR10'

                        Variable%Value%Int_3D => Variable%Value%Int_4D(:,:,:,t)

                        call HDF5WriteData(Me%ObjHDF5,                                  &
                                           "/Results/"//trim(Variable%ID%Name),         &
                                           trim(Variable%ID%Name),                      &
                                           Variable%Value%Units,                        &
                                           Array3D      = Variable%Value%Int_3D,        &
                                           OutputNumber = t,                            &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleIntegers - ModuleWOAFormat - ERR10'
                        
                        nullify(Variable%Value%Int_3D)

                    end if

                    if(Me%ConvertToASCII)  call WriteDataLine(Me%OutputUnit,"<end_field>")

                end do
                
                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<end_data>")

                deallocate(Variable%Value%Int_4D)

            case default


        end select


    end subroutine HandleIntegers

    !--------------------------------------------------------------------------

    subroutine HandleShorts(Variable)
        
        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable

        !Local-----------------------------------------------------------------
        integer                                         :: stat, i, j, k
        type(T_Size3D)                                  :: Size
        integer                                         :: STAT_CALL
        integer                                         :: NVCount
        logical                                         :: GoodValue
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        real                                            :: Year, Month, Day, Hour 
        real                                            :: Minute, Second

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size

        Me%PreviousDate = Me%InitialDate

        select case(Variable%nDimensions)
            
            case(0)
                
                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short)
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleWOAFormat - ERR01'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Short

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


            case(1)

                nullify (Variable%Value%Short_1D)
                allocate(Variable%Value%Short_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short_1D)
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleWOAFormat - ERR02'


                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                        write(Me%OutputUnit,*)Variable%Value%Short_1D(i)
                    enddo
                
                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                deallocate(Variable%Value%Short_1D)

            case(2)

                nullify (Variable%Value%Short_2D)
                allocate(Variable%Value%Short_2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short_2D)
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleWOAFormat - ERR03'
                

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                        write(Me%OutputUnit,*)Variable%Value%Short_2D(i,j)
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                deallocate(Variable%Value%Short_2D)

            
            case(3)

                nullify (Variable%Value%Short_3D)
                allocate(Variable%Value%Short_3D(Me%WorkSize%JLB:Me%WorkSize%JUB,       & 
                         Me%WorkSize%ILB:Me%WorkSize%IUB, Size%KLB:Size%KUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short_3D)
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleWOAFormat - ERR04'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                    do k = Size%KLB, Size%KUB
                        write(Me%OutputUnit,*)Variable%Value%Short_3D(i,j,k)
                    enddo
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


                if(Me%ConvertToHDF5)then

                    Me%HDFOutputNumber = 1

                    do k = Size%KLB, Size%KUB !K is understood as time instant
                        
                        !Get units of this variable:
                        Variable%Value%Units = Me%PropUnits

                        GoodValue = .true.

                        NVCount = 0
                        
                        allocate(Variable%Value%Float_2D(Me%WorkSize%ILB:Me%WorkSize%IUB, &
                                 Me%WorkSize%JLB:Me%WorkSize%JUB))
                       
doloop1 :               do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        !if (.not. GoodValue) exit doloop1
doloop2 :               do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                            Variable%Value%Float_2D(i, j)=                              & 
                                Me%Add_ScaleFactor *                                    &
                                real(Variable%Value%Short_3D(j,Me%WorkSize%IUB-i+1,k))  &
                                + Me%Add_OffSet
                            !Check if value is valid: if not then exit loop and not write
                            if (Variable%Value%Float_2D(i, j) == Me%Missing_Value) then
                               GoodValue = .false.
                               NVCount = NVCount + 1
                            end if     
                        enddo doloop2 
                        enddo doloop1
                        
                        !Get time for this instant
                        call ConstructInstantTime(k, AuxTime)
                        
                        if (.not. GoodValue) then
889                         format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0,   & 
                                    1x, f3.0)
                            call ExtractDate(Me%CurrentDate, Year, Month, Day, Hour,    &  
                                             Minute, Second)                        
                            write(*,*)'Current time:'
                            write(*,fmt=889) Year, Month, Day, Hour, Minute, Second
                            write(*,*)'Number of non valid values:', NVCount

                         else
                            !Check if time is after previous time
                            if (Me%CurrentDate <= Me%PreviousDate) then
                                write(*,*)'Current time before or equal previous!'
                            else

                            call ConvertToMohidUnits3D(Variable)

                            !Write time
                            TimePtr => AuxTime

                            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleWOAFormat - ERR05'

                            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                     Array1D = TimePtr,                                 &
                                     OutputNumber = Me%HDFOutputNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleWOAFormat - ERR06'

                            !Write values
                            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB,             & 
                                               Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                                               Me%WorkSize%JUB,                         &
                                           STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleWOAFormat - ERR07'

                            call HDF5WriteData(Me%ObjHDF5,                              &
                                           "/Results/"//trim(Me%PropName),              &
                                           trim(Me%PropName), Variable%Value%Units,     &
                                           Array2D      = Variable%Value%Float_2D,      &
                                           OutputNumber = Me%HDFOutputNumber,           &
                                           STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleWOAFormat - ERR08'

                            deallocate(Variable%Value%Float_2D)

                            !Write instant time
                            write(*,*)'writing instant ', k

                            Me%PreviousDate = Me%CurrentDate
                            Me%HDFOutputNumber = Me%HDFOutputNumber + 1

                            endif

                        endif

                    enddo

                endif

                deallocate(Variable%Value%Short_3D)

            case default

        end select

    end subroutine HandleShorts

    !--------------------------------------------------------------------------

    subroutine HandleChars(Variable)
        
        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable

        !Local-----------------------------------------------------------------
        integer                                         :: stat, i
        type(T_Size3D)                                  :: Size

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size

        select case(Variable%nDimensions)
            
            case(0)
                
                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Char)
                if (stat /= nf90_noerr) stop 'HandleChars - ModuleWOAFormat - ERR01'


                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Char

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


            case(1)

                nullify (Variable%Value%Char_1D)
                allocate(Variable%Value%Char_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Char_1D)
                if (stat /= nf90_noerr) stop 'HandleChars - ModuleWOAFormat - ERR02'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                        write(Me%OutputUnit,*)Variable%Value%Char_1D(i)
                    enddo
                
                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                deallocate(Variable%Value%Char_1D)

            case(2,3)

                write(*,*) 'Could not convert 2D/3D array of characters.'
                write(*,*)
                
            case default


        end select


    end subroutine HandleChars

    !------------------------------------------------------------------
    
    subroutine ReadGridMapOptions

        !Local---------------------------------------------------------          
        integer                                         :: iflag, STAT_CALL

        !Begin---------------------------------------------------------

        call GetData(Me%XX_Variable,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'XX_VARIABLE',                  &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR01'


        call GetData(Me%YY_Variable,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'YY_VARIABLE',                  &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR02'


        call GetData(Me%ZZ_Variable,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'ZZ_VARIABLE',                  &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR03'


        call GetData(Me%Time_Variable,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'TIME_VARIABLE',                &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR04'


        call GetData(Me%GridType_Variable,                          &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'GRID_TYPE_VARIABLE',           &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR05'

        call GetData(Me%VertCoord_Variable,                         &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'VERT_COORD_VARIABLE',          &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR06'

        call GetData(Me%Batim_Variable,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'BATIM_VARIABLE',               &
                     ClientModule = 'ModuleWOAFormat',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleWOAFormat - ERR07'

    end subroutine ReadGridMapOptions

    !------------------------------------------------------------------
    
    subroutine ConstructTime

        !Local---------------------------------------------------------          
        type(T_Variable), pointer           :: Time
        integer                             :: stat_call,l
        integer                             :: stat
        integer                             :: t
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr

        !Begin---------------------------------------------------------          

        nullify(Time); allocate(Time)

        Time%ID%Name = Me%Time_Variable

        stat = nf90_inq_varid(Me%FileID, trim(Time%ID%Name), Time%ID%Number)
        if (stat /= nf90_noerr) return

        stat = nf90_inquire_variable(Me%FileID, Time%ID%Number, Time%ID%Name,       &
                                     Time%ID%Type_, Time%nDimensions,               &
                                     Time%DimensionsID, Time%nAttributes)
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleWOAFormat - ERR01'

        !if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<begin_time>")

        call HandleDimensions (Time, .true.)

        call HandleAttributtes(Time%ID%Number,  Time%nAttributes)

        nullify(Me%TimeValues);allocate(Me%TimeValues(Time%Dim%Size%ILB:Time%Dim%Size%IUB))
        nullify(Me%TimeArray );allocate(Me%TimeArray (Time%Dim%Size%ILB:Time%Dim%Size%IUB))


        stat = nf90_get_var(Me%FileID, Time%ID%Number, Me%TimeValues)
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleWOAFormat - ERR02'


        Me%Time_LB = Time%Dim%Size%ILB
        Me%Time_UB = Time%Dim%Size%IUB







    do t = 1,12
        Me%NumDates = 12
        allocate(Me%OutPutDates(Me%NumDates))  
        do l = 1, 12
            call SetDate(Me%OutPutDates(l),                                         &
                        Year    = CyclicTime,                                       &
                        Month   = real(l),                                          & 
                        Day     = 1.,                                               &
                        Hour    = 0.,                                               &
                        Minute  = 0.,                                               &
                        Second  = 0.) 
        enddo

                call ExtractDate   (Me%OutPutDates(t),                                      &
                AuxTime(1), AuxTime(2), AuxTime(3),                     &
                AuxTime(4), AuxTime(5), AuxTime(6))



            !Write time
            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                  &
                stop 'HandleShorts - ModuleWOAFormat - ERR05'

            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                     Array1D = TimePtr,                                 &
                     OutputNumber = t, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                  &
                stop 'HandleShorts - ModuleWOAFormat - ERR06'

            !Write instant time


            Me%PreviousDate = Me%CurrentDate
            Me%HDFOutputNumber = Me%HDFOutputNumber + 1


            deallocate(Time)
        enddo
       write(*,*)'Writing Time '

    end subroutine ConstructTime

    !------------------------------------------------------------------
    
    subroutine ConstructBathymetry

        !Local---------------------------------------------------------          
        integer                                         :: stat
        real, dimension(:,:), pointer                   :: tempBatim

        !Begin---------------------------------------------------------          

        nullify(Me%Batim); allocate(Me%Batim)

        stat = nf90_inq_varid(Me%FileID, trim(Me%Batim_Variable), Me%Batim%ID%Number)
        if (stat /= nf90_noerr) then
           
            nullify (Me%Batim%Value%Float_2D)
            allocate(Me%Batim%Value%Float_2D(Me%Size%ILB:Me%Size%IUB,                   & 
                     Me%Size%JLB:Me%Size%JUB))
            Me%Batim%Value%Float_2D = 0.


        else

            stat = nf90_inq_varid(Me%FileID, trim(Me%Batim_Variable), Me%Batim%ID%Number)

            stat = nf90_get_var(Me%FileID, Me%Batim%ID%Number, Me%Batim%Value%Int)
            if (stat /= nf90_noerr) stop 'ConstructBathymetry - ModuleWOAFormat - ERR01'

            stat = nf90_inquire_variable(Me%FileID, Me%Batim%ID%Number,                 & 
                                         Me%Batim%ID%Name,                              &
                                         Me%Batim%ID%Type_, Me%Batim%nDimensions,       &
                                         Me%Batim%DimensionsID, Me%Batim%nAttributes)
            if (stat /= nf90_noerr) stop 'ConstructBathymetry - ModuleWOAFormat - ERR02'

            call HandleDimensions (Me%Batim, .false.)
            call HandleAttributtes(Me%Batim%ID%Number,  Me%Batim%nAttributes)

            nullify (Me%Batim%Value%Float_2D, tempBatim)
            allocate(Me%Batim%Value%Float_2D(Me%Size%ILB:Me%Size%IUB,                   & 
                     Me%Size%JLB:Me%Size%JUB))
            allocate(tempBatim(Me%Size%JLB:Me%Size%JUB, Me%Size%ILB:Me%Size%IUB))

            !switch i,j
            stat = nf90_get_var(Me%FileID, Me%Batim%ID%Number,                          &    
                                tempBatim(Me%WorkSize%JLB:Me%WorkSize%JUB,              &
                                          Me%WorkSize%ILB:Me%WorkSize%IUB))
            if (stat /= nf90_noerr) stop 'ConstructBathymetry - ModuleWOAFormat - ERR03'


            Me%Batim%Value%Float_2D = TRANSPOSE(tempBatim)

        endif

    end subroutine ConstructBathymetry

    !------------------------------------------------------------------

    subroutine ConstructMap

        !Local---------------------------------------------------------          
        integer                                         :: i, j

        !Begin---------------------------------------------------------          

        nullify (Me%Mapping2D)
        
        if(Me%Size%ILB     < 0) Me%Size%ILB     = 1
        if(Me%Size%IUB     < 0) Me%Size%IUB     = 1
        if(Me%Size%JLB     < 0) Me%Size%JLB     = 1
        if(Me%Size%JUB     < 0) Me%Size%JUB     = 1


        allocate(Me%Mapping2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB

            if(Me%Batim%Value%Float_2D(i,j) <= -10.)then
                Me%Mapping2D(i,j) = 0
            else
                Me%Mapping2D(i,j) = 1
            endif

        enddo
        enddo

    end subroutine ConstructMap

    !--------------------------------------------------------------------------

    subroutine WriteDimensions(Index, Variable, Size, WriteDims)

        !Arguments-------------------------------------------------------------
        integer                                         :: Index
        type(T_Variable), pointer                       :: Variable
        integer                                         :: Size
        logical                                         :: WriteDims
        
        !Local-----------------------------------------------------------------
        integer                                         :: stat

        !Begin-----------------------------------------------------------------
        
        stat = nf90_inquire_dimension(Me%FileID, Variable%DimensionsID(Index),      & 
                                      Variable%Dim%Name, Size)


        if(WriteDims)then
            if(Me%ConvertToASCII)then
                call WriteDataLine(Me%OutputUnit,"<begin_dimension>")
                call WriteDataLine(Me%OutputUnit,"Name", trim(Variable%Dim%Name))
                call WriteDataLine(Me%OutputUnit,"Size", Size)
                call WriteDataLine(Me%OutputUnit,"<end_dimension>")
            endif
        end if

    end subroutine WriteDimensions

    !------------------------------------------------------------------

    subroutine ConstructGrid

        !Local---------------------------------------------------------          
        type(T_Variable), pointer                       :: XX, YY, ZZ, GridType
        integer                                         :: stat
        integer                                         :: i, j
        real                                            :: dx, dy

        !Begin---------------------------------------------------------          

        nullify(XX, YY, ZZ, GridType); allocate(XX, YY, ZZ, GridType)

        call GetVariable(XX, trim(Me%XX_Variable))
        call GetVariable(YY, trim(Me%YY_Variable))
        call GetVariable(ZZ, trim(Me%ZZ_Variable))


        stat = nf90_inq_varid(Me%FileID, trim(Me%GridType_Variable), GridType%ID%Number)
        if (stat /= nf90_noerr) then
            write(*,*)'No grid type defined. Assumed cartesian...'
            write(*,*)
            Me%Grid_Type = Cartesian
        else
            stat = nf90_get_var(Me%FileID, GridType%ID%Number, GridType%Value%Int)
            if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleWOAFormat - ERR01'

            Me%Grid_Type = GridType%Value%Int
        end if
        
        select case(Me%Grid_Type)

            case(Cartesian)

                write(*,*)'Grid type: Cartesian'

                Me%WorkSize%ILB = 1
                Me%WorkSize%JLB = 1
                Me%WorkSize%KLB = 1
                Me%WorkSize%IUB = YY%Dim%Size%IUB
                Me%WorkSize%JUB = XX%Dim%Size%IUB
                Me%WorkSize%KUB = ZZ%Dim%Size%IUB




                Me%Size%ILB = Me%WorkSize%ILB - 1
                Me%Size%JLB = Me%WorkSize%JLB - 1
                Me%Size%IUB = Me%WorkSize%IUB + 1
                Me%Size%JUB = Me%WorkSize%JUB + 1
                Me%Size%KUB = ZZ%Dim%Size%IUB


                nullify(Me%ConnectionX, Me%ConnectionY)

                allocate(Me%ConnectionX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                allocate(Me%ConnectionY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                
                nullify (XX%Value%FLOAT_1D, YY%Value%FLOAT_1D)
                allocate(XX%Value%FLOAT_1D(Me%Size%JLB:Me%Size%JUB))
                allocate(YY%Value%FLOAT_1D(Me%Size%ILB:Me%Size%IUB))

                stat = nf90_get_var(Me%FileID, XX%ID%Number,                            & 
                                    XX%Value%Float_1D(Me%WorkSize%JLB:Me%WorkSize%JUB))
                if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleWOAFormat - ERR02'

                stat = nf90_get_var(Me%FileID, YY%ID%Number,                            & 
                                    YY%Value%Float_1D(Me%WorkSize%ILB:Me%WorkSize%IUB))
                if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleWOAFormat - ERR03'


                dx = XX%Value%Float_1D(2) - XX%Value%Float_1D(1)
                dy = YY%Value%Float_1D(1) - YY%Value%Float_1D(2)

                !XOrigin = XX%Value%Float_1D(1) - dx/2.
                !YOrigin = YY%Value%Float_1D(YY%Dim%Size%IUB) - dy/2.

                !do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                !do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    !if (XX%Value%Float_1D(j  ) > 180.)then
                    !    Me%ConnectionX(i,j) = XX%Value%Float_1D(j  )   -360                &                     
                    !                      - dx/2.
                    !else
                    !    Me%ConnectionX(i,j) = XX%Value%Float_1D(j  )                        &                     
                    !                      - dx/2.
                    !end if
                    !Me%ConnectionY(i,j) = YY%Value%Float_1D(YY%Dim%Size%IUB - i + 1)    &
                    !                      - dy/2.
                    !Me%ConnectionY(i,j) = YY%Value%Float_1D(i)    &
                    !                      - dy/2.


                !enddo
                !enddo

                !Me%ConnectionY(Me%WorkSize%IUB + 1, :) =                                &
                !                          Me%ConnectionY(Me%WorkSize%IUB, :) + dy
                !Me%ConnectionX(Me%WorkSize%IUB + 1, :) =                                &
                !                          Me%ConnectionX(Me%WorkSize%IUB, :)

                !Me%ConnectionX(:, Me%WorkSize%JUB + 1) = Me%ConnectionX(:, Me%WorkSize%JUB) + dx
                !Me%ConnectionY(:, Me%WorkSize%JUB + 1) = Me%ConnectionY(:, Me%WorkSize%JUB)
                
                Me%ConnectionX(:,1) = -180.

                do j = 2, 361
                    Me%ConnectionX(:,j) = Me%ConnectionX(:,j-1) + 1.
                end do


                Me%ConnectionY(1,:) = -90.

                do i = 2, 181
                    Me%ConnectionY(i,:) = Me%ConnectionY(i-1,:) + 1.
                end do



                if(Me%ConvertToASCII)then
                    deallocate(Me%ConnectionX, Me%ConnectionY)
                    nullify   (Me%ConnectionX, Me%ConnectionY)
                endif

            case(Spherical)

                write(*,*)'grid type: Spherical'
                stop 'Spherical grids not yet supported'

            case(Curvilinear)

                write(*,*)'grid type: Curvilinear'
                stop 'Curvilinear grids not yet supported'

            case(Spherical_Curvilinear)

                write(*,*)'grid type: Spherical_Curvilinear'
                stop 'Spherical_Curvilinear grids not yet supported'


        end select

  

    end subroutine ConstructGrid

    !------------------------------------------------------------------

    subroutine GetVariable(Variable, Name)

        !Local---------------------------------------------------------          
        type(T_Variable), pointer                       :: Variable
        character(len=*)                                :: Name
        integer                                         :: stat

        !Begin---------------------------------------------------------          

        stat = nf90_inq_varid(Me%FileID, trim(Name), Variable%ID%Number)
        if (stat /= nf90_noerr) return

        stat = nf90_inquire_variable(Me%FileID, Variable%ID%Number, Variable%ID%Name,   &
                                     Variable%ID%Type_, Variable%nDimensions,           &
                                     Variable%DimensionsID, Variable%nAttributes)
        if (stat /= nf90_noerr) stop 'GetVariable - ModuleWOAFormat - ERR01'

        call HandleDimensions (Variable, .false.)

        call HandleAttributtes(Variable%ID%Number,  Variable%nAttributes, .false.)


    end subroutine GetVariable

    !--------------------------------------------------------------------------

    subroutine ConstructInstantTime(t, AuxTime)

        !Local---------------------------------------------------------          
        type(T_Variable), pointer           :: Time
        integer                             :: t
        real,    dimension(6), target       :: AuxTime
        real                                :: Year, Month, Day, Hour, Minute, Second

        !Begin---------------------------------------------------------          

        !Next part can be performed for each time instant together with instant value 
        1001 format(f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)

        if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<begin_time>")

        if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<begin_dates>")

        call ExtractDate(Me%TimeArray(t), Year, Month, Day, Hour, Minute, Second) 

        if(Me%ConvertToASCII) write(Me%OutputUnit, 1001) Year, Month, Day, Hour,    &
                                                         Minute, Second

        if(Me%ConvertToHDF5)then

            call ExtractDate   (Me%TimeArray(t),                                    &
                                AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                AuxTime(4), AuxTime(5), AuxTime(6))

                call SetDate(Me%CurrentDate, AuxTime(1),AuxTime(2),AuxTime(3),      &
                             AuxTime(4),AuxTime(5),AuxTime(6))

        end if

        if(Me%ConvertToASCII)call WriteDataLine(Me%OutputUnit, "<end_dates>")

        if(Me%ConvertToASCII)call WriteDataLine(Me%OutputUnit, "<end_time>")

        deallocate(Time)

    end subroutine ConstructInstantTime



    !------------------------------------------------------------------

end module ModuleWOAFormat








