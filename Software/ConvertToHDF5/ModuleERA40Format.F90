!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ERA40Format
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : Module to convert ERA40Format (NETCDF) files into HDF5 format.
!
!------------------------------------------------------------------------------
!DataFile
!
!   ACTION                      : char              -           !'CONVERT ERA40 FORMAT' to use this module
!   FILENAME                    : char              -           !Path to ERA40 original file
!   OUTPUTFILENAME              : char              -           !Path to ERA40 file generated from NETCDF file
!                                                               !(root of name for all files produced)
!
!   CONVERT_TO_ASCII            : 0/1              [0]          !Flag to convert to ascii file
!   CONVERT_TO_HDF5             : 0/1              [0]          !Flag to convert to hdf5 file
!   GRIDTO180                   : 0/1              [0]          !Flag to put grid from [0 360] to [-180 180]
!
!   XX_VARIABLE                 : char              -           !'longitude'
!   YY_VARIABLE                 : char              -           !'latitude'
!   TIME_VARIABLE               : char              -           !'time'
!
! In future version should be allowed to choose the variables to convert (file can become very large)

! In case CONVERT_TO_HDF5 is chosen then a hdf5 file is produced for each data variable.
! Appended to the OUTPUTFILENAME is included the short name for the variable:
!
!         ---ERA NAME---    ---MOHID NAME---
!
!            sshf           sensible heat
!            slhf           latent heat
!            msl            atmospheric pressure
!            tcc            cloud cover
!            p10u           wind velocity X
!            p10v           wind velocity Y
!            p2t            air temperature
!            ewss           wind stress X
!            nsss           wind stress Y
!            ssrd           solar radiation
!            r              relative humidity

Module ModuleERA40Format

#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif    
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
    public  :: ConvertERA40Format
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
    private ::      KillERA40Format   

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
    
    type       T_ERA40Format
        integer                                             :: ObjEnterData             = 0
        integer                                             :: ObjHDF5                  = 0
        integer                                             :: ObjHorizontalGrid        = 0
        integer                                             :: Unit
        integer                                             :: OutputUnit
        character(len=PathLength)                           :: InputFileName
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
        logical                                             :: InGridFileNotSpecified  = .false.
        logical                                             :: OutputGridNotSpecified  = .false.
        logical                                             :: ConvertToHDF5
        logical                                             :: ConvertToASCII
        logical                                             :: GridTo180 = .false.
        real,        dimension(:    ), pointer              :: TimeValues
        integer                                             :: Time_LB, Time_UB
        integer                                             :: nDimensions, nVariables, nAttributes
        integer                                             :: HDFOutputNumber
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
        !Maybe the next are not needed:
        real                                                :: PTop                   = null_real
        type(T_Field),                pointer               :: FirstField
        type(T_Date),                 pointer               :: FirstDate
        character(len=StringLength), dimension(:), pointer  :: FieldsToConvert
        character(len=StringLength), dimension(:), pointer  :: FieldsToRead

    end type  T_ERA40Format

    type(T_ERA40Format), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertERA40Format(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions

        call OpenReadAndConvertFile

        call KillERA40Format

        STAT = SUCCESS_


    end subroutine ConvertERA40Format

    !------------------------------------------------------------------------

    subroutine ReadOptions
        
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag

        !Begin-----------------------------------------------------------------
       
        call GetData(Me%InputFileName,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'FILENAME',                     &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleERA40Format - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleERA40Format - ERR02'
        end if

        call GetData(Me%OutputFileName,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'OUTPUTFILENAME',               &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleERA40Format - ERR04'      

        call GetData(Me%ConvertToHDF5,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'CONVERT_TO_HDF5',              &
                     Default      = OFF,                            &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleERA40Format - ERR05'

        call GetData(Me%ConvertToASCII,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'CONVERT_TO_ASCII',             &
                     Default      = OFF,                            &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleERA40Format - ERR06'

        call GetData(Me%GridTo180,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'GRIDTO180',                    &
                     Default      = .false.,                        &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)      
        !If not present, don't stop. 
        !if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleERA40Format - ERR06.5'


    end subroutine ReadOptions
   
    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits2D(Field)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                         :: Field
        
        !Begin-----------------------------------------------------------------

        select case(trim(Field%Units))


            case('K')

                Field%Units     = 'ºC'
                Field%Values2D  = Field%Values2D - AbsoluteZero
            
            case('m s{-1}')

                Field%Units     = 'm/s'

            case('cm')

                Field%Units     = 'm/s'
                Field%Values2D  = Field%Values2D / 3600. / 100.

            case('W/m^2')

                Field%Units     = 'W/m2'

            case default

        end select

    end subroutine ConvertToMohidUnits2D

    !------------------------------------------------------------------------

    subroutine ConvertToMohidUnits3D(Variable)

        !Arguments-------------------------------------------------------------
        type(T_Variable), pointer                       :: Variable

        !Local-----------------------------------------------------------------
                
        !Begin-----------------------------------------------------------------

        select case(trim(Variable%Value%Units))

            case('K')

                Variable%Value%Units     = 'ºC'
                Variable%Value%Float_2D  = Variable%Value%Float_2D - AbsoluteZero
            
            case('m s**-1')

                Variable%Value%Units     = 'm/s'

            case('N m**-2 s')

                Variable%Value%Units     = 'N/m2'
                Variable%Value%Float_2D  = Variable%Value%Float_2D/(3600.*6.) !just guessing this is right

            case('W m**-2 s')

                Variable%Value%Units     = 'W/m2'
                Variable%Value%Float_2D  = Variable%Value%Float_2D/(3600.*6.) !just guessing this is right

            case default

        end select

    end subroutine ConvertToMohidUnits3D

    !------------------------------------------------------------------------

    ! This is not considered for now, although may be relevant in future versions   
    !subroutine ComputeVerticalCoordinate
        
    !    !Local-----------------------------------------------------------------
    !    type(T_Field), pointer                  :: VerticalZ
    !    real, dimension(:,:  ), pointer         :: ReferenceSurfacePressure
    !    real, dimension(:,:,:), pointer         :: Pressure3D
    !    real, dimension(:,:,:), pointer         :: Temperature
    !    logical                                 :: Pref_OK , Pressure_OK, Temperature_OK            
    !    integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
    !    integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
    !    integer                                 :: i,j,k
    !    type(T_Field), pointer                  :: Field
    !    type(T_Date), pointer                   :: CurrentDate

    !    !Begin-----------------------------------------------------------------

    !    write(*,*)
    !    write(*,*)'Computing Vertical Coordinate...'

    !    ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
    !    IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
    !    JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
    !    JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
    !    KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
    !    KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 

    !    Temperature_OK  = .false.
    !    Pref_OK         = .false.
    !    Pressure_OK     = .false.


    !    CurrentDate => Me%FirstDate

    !    do while(associated(CurrentDate))

    !        Field => Me%FirstField

    !        do while(associated(Field))
            
    !            if(Field%Name == GetPropertyName(AtmosphericPressure_)              .and. &
    !               Field%Date == CurrentDate%Date)then
                    
    !                ReferenceSurfacePressure    => Field%Values2D
    !                Pref_OK                     = .true.

    !            end if

    !            if(Field%Name == trim(GetPropertyName(AtmosphericPressure_))//"_3D" .and. &
    !               Field%Date == CurrentDate%Date)then
                    
    !                Pressure3D                  => Field%Values3D
    !                Pressure_OK                 = .true.

    !            end if

    !            if(Field%Name == trim(GetPropertyName(AirTemperature_))//"_3D"      .and. &
    !               Field%Date == CurrentDate%Date)then
                    
    !                Temperature                 => Field%Values3D
    !                Temperature_OK              = .true.

    !            end if

    !            if(Pref_OK .and. Pressure_OK .and. Temperature_OK)then

    !                call AddField(Me%FirstField, VerticalZ)

    !                call SetNewFieldAttributes(Field         =  VerticalZ,          &
    !                                           Name          = 'VerticalZ',         &
    !                                           Units         = 'm',                 &
    !                                           Date          = CurrentDate%Date,    &
    !                                           WorkSize      = Me%WorkSize,         &
    !                                           nDimensions   = 3,                   &
    !                                           Convert       = .true.)

    !                allocate(VerticalZ%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))
    !                VerticalZ%Values3D = 0.
                    
    !                do j = WJLB, WJUB
    !                do i = WILB, WIUB
    !                do k = WKLB, WKUB
                        
    !                    !Pressure(i,j,k) = Sigma(k) * PSTARCRS(i,j) + PTOP + PP(i,j,k)
    !                    Pressure3D(i,j,k) = Me%Sigma(WKUB + 1 - k) * ReferenceSurfacePressure(i,j) + &
    !                                        Me%PTop + Pressure3D(i,j,k)

    !                enddo
    !                enddo
    !                enddo

    !                ReferenceSurfacePressure(WILB:WIUB, WJLB:WJUB) = Pressure3D(WILB:WIUB, WJLB:WJUB, WKLB) 

    !                Pressure3D(ILB:IUB, JLB:JUB,WKUB + 1) = Me%Ptop

    !                VerticalZ%Values3D(WILB:WIUB, WJLB:WJUB, 0) = - Me%Bathymetry(WILB:WIUB, WJLB:WJUB)

    !                do j = WJLB, WJUB
    !                do i = WILB, WIUB
    !                do k = WKLB+1, WKUB+1

    !                    VerticalZ%Values3D(i,j,k-1) = VerticalZ%Values3D(i,j,k-2) + &
    !                                                  (Pressure3D(i,j,k) - Pressure3D(i,j,k-1))/ &
    !                                                  (-Gravity * Pressure3D(i,j,k-1)     / &
    !                                                  (Perfect_gas_R * Temperature(i,j,k-1)))

    !                enddo
    !                enddo
    !                enddo

    !                Temperature_OK  = .false.
    !                Pref_OK         = .false.
    !                Pressure_OK     = .false.

    !                nullify(Pressure3D, ReferenceSurfacePressure) 

    !            end if

    !            Field => Field%Next

    !        end do


    !        CurrentDate => CurrentDate%Next

    !    end do 


    !end subroutine ComputeVerticalCoordinate

    !----------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%HDF5FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleERA40Format - ERR01'
        

    end subroutine Open_HDF5_OutPut_File


    !----------------------------------------------------------------------
   
    subroutine KillERA40Format
        
        !Local-----------------------------------------------------------------
        integer                                         :: nUsers
        
        !Begin-----------------------------------------------------------------
       


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillERA40Format - ModuleERA40Format - ERR03'

        deallocate(Me%Mapping2D)
        nullify   (Me%Mapping2D)
        deallocate(Me%ConnectionX, Me%ConnectionY)
        nullify   (Me%ConnectionX, Me%ConnectionY)
        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillERA40Format

    !--------------------------------------------------------------------------

    subroutine OpenReadAndConvertFile

        !Local-----------------------------------------------------------------
        integer                                         :: stat
        integer                                         :: VariableID, STAT_CALL
        type(T_Variable), pointer                       :: Variable
        character(len=PathLength)                       :: AuxFileName

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)"Converting netcdf file. Please wait..."
        write(*,*)

        if(Me%ConvertToASCII)then
            !Only one ASCII file is made: caution it can become very large
            call OpenASCIIFile

        endif

        call OpenNetCDFFile


        call ConstructGridAndMap


        if(Me%nAttributes > 0)then

            !Handle global attributes
            call HandleAttributtes(NF90_GLOBAL, Me%nAttributes)

        endif


        do VariableID = 1, Me%nVariables

            nullify(Variable); allocate(Variable)

            if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<begin_variable>")


            Variable%ID%Number = VariableID

            stat = nf90_inquire_variable(Me%FileID, Variable%ID%Number, Variable%ID%Name,   &
                                         Variable%ID%Type_, Variable%nDimensions,           &
                                         Variable%DimensionsID, Variable%nAttributes)
            if (stat /= nf90_noerr) stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR01'

            write(*,*)trim(Variable%ID%Name)

            !/Previous is only to get to know the variable: no calculations are made
            if(.not. IsGridVariable(Variable))then

                if(Me%ConvertToHDF5)then
                    !Arrange name for the HDF5 output file for this variable
                    AuxFileName = trim(Me%OutputFileName)//trim(Variable%ID%Name)
                    Me%HDF5FileName = trim(AuxFileName)//".hdf5"

                    !Arrange a new HDF5 output file for this variable
                    call Open_HDF5_OutPut_File

                    !Set grid
                    call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,  &
                                       Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR02'

                    !ConnectionX is not really a connection X but for could be useful
                    !if HDF5Extractor use is required (here is equal to Longitude)
                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX",           &
                                          "degrees",                                    &
                                          Array2D = Me%ConnectionX,                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR03a'

                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Longitude",             &
                                          "degrees",                                    &
                                          Array2D = Me%ConnectionX,                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR03b'

                    !ConnectionY is not really a connection Y but for could be useful
                    !if HDF5Extractor use is required (here is equal to Latitude)
                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY",           & 
                                          "degrees", Array2D = Me%ConnectionY,          &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR04a'

                    call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Latitude",              &
                                          "degrees", Array2D = Me%ConnectionY,          &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR04b'

                    !Set bathymetry          
                    call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB,Me%WorkSize%IUB,     &
                               Me%WorkSize%JLB,Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR05'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                               "/Grid",                                                 &
                               "Bathymetry",                                            &
                               'm',                                                     &
                               Array2D      = Me%Batim%Value%Float_2D,                  &
                               STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR06'

                    !Set map
                    call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,    & 
                                       Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR07'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                   "/Grid",                                             &
                                   "Mapping2D",                                         &
                                   '-',                                                 &
                                   Array2D      = Me%Mapping2D,                         &
                                   STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR08'

                    call HDF5FlushMemory(Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR09'

                    !flush mem
                    !deallocate(Me%ConnectionX, Me%ConnectionY,Me%Mapping2D,Me%Batim%Value%Float_2D)
                    !nullify   (Me%ConnectionX, Me%ConnectionY,Me%Mapping2D,Me%Batim%Value%Float_2D)


                endif

                if(Me%ConvertToASCII)then
                    call WriteDataLine(Me%OutputUnit,"Name", trim(Variable%ID%Name))
                    call WriteDataLine(Me%OutputUnit,"Dimensions", Variable%nDimensions)
                    call WriteDataLine(Me%OutputUnit,"Attributes", Variable%nAttributes)
                    call WriteDataType(Variable%ID%Type_)
                end if

                call HandleDimensions (Variable)

                call HandleAttributtes(Variable%ID%Number,  Variable%nAttributes)
                
                if(Me%ConvertToHDF5)then

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

                call HDF5FlushMemory(Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR10'

                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<end_variable>")

                call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenReadAndConvertFile - ModuleERA40Format - ERR11'

            endif
        
        end do

        write(*,*)'Finished converting'


    end subroutine OpenReadAndConvertFile

    !------------------------------------------------------------------
    
    subroutine OpenASCIIFile

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !----------------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenASCIIFile - ModuleERA40Format - ERR01'

        open(Unit   = Me%OutputUnit,                                            &
             File   = Me%OutputFileName,                                        &
             STATUS = 'UNKNOWN',                                                &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenASCIIFile - ModuleERA40Format - ERR02'

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
    
    subroutine OpenNetCDFFile

        !Local-----------------------------------------------------------------
        integer                                         :: stat
        logical                                         :: exist

        !----------------------------------------------------------------------
        

        !Verifies if file exists
        inquire(file = Me%InputFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'NETCDF file does not exist'
            stop 'OpenNetCDFFile - ModuleERA40Format - ERR01'
        endif

        stat = nf90_open(Me%InputFileName, NF90_NOWRITE, Me%FileID)
        if (stat /= nf90_noerr) stop 'OpenNetCDFFile - ModuleERA40Format - ERR02'

        stat = nf90_inquire(Me%FileID, Me%nDimensions, Me%nVariables, Me%nAttributes)
        if (stat /= nf90_noerr) stop 'OpenNetCDFFile - ModuleERA40Format - ERR03'


        if(Me%ConvertToASCII)then

            call WriteDataLine(Me%OutputUnit,"Number of variables", Me%nVariables)
            call WriteDataLine(Me%OutputUnit,"Dimensions", Me%nDimensions)
            call WriteDataLine(Me%OutputUnit,"Global Attributes", Me%nAttributes)

        endif

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
            if (stat /= nf90_noerr) stop 'HandleAttributtes - ModuleERA40Format - ERR01'

            stat = nf90_inquire_attribute(Me%FileID, TypeOfAttribute, Attribute%ID%Name, &
                                          Attribute%ID%Type_,                            &
                                          Attribute%ID%Length, Attribute%ID%Number)            
            if (stat /= nf90_noerr) stop 'HandleAttributtes - ModuleERA40Format - ERR02'
            
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
                        stop 'HandleAttributtes - ModuleERA40Format - ERR03'

                case(NF90_CHAR)

                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Char)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleERA40Format - ERR04'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        call WriteDataLine(Me%OutputUnit,"Value", Attribute%Value%Char)
                   

                case(NF90_SHORT)
                    
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Short)
                    if (stat /= nf90_noerr)                                              &
                        stop 'HandleAttributtes - ModuleERA40Format - ERR05'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        write(Me%OutputUnit,*)"Value            :", Attribute%Value%Short

                case(NF90_INT)
                   
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Int)
                    if (stat /= nf90_noerr)                                              &
                        stop 'HandleAttributtes - ModuleERA40Format - ERR06'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             &
                        call WriteDataLine(Me%OutputUnit, "Value", Attribute%Value%Int)


                case(NF90_FLOAT)
                   
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Float)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleERA40Format - ERR07'

                    if(WriteToFile_ .and. Me%ConvertToASCII)                             & 
                        call WriteDataLine(Me%OutputUnit,"Value", Attribute%Value%Float)


                case(NF90_DOUBLE)
                    
                    stat =  nf90_get_att(Me%FileID, TypeOfAttribute, Attribute%ID%Name,  & 
                                         Attribute%Value%Double)
                    if (stat /= nf90_noerr)                                              & 
                        stop 'HandleAttributtes - ModuleERA40Format - ERR08'

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

                if (Me%PropName == 'Surface sensible heat flux' ) Me%PropName=GetPropertyName(SensibleHeat_)
                if (Me%PropName == 'Surface latent heat flux'   ) Me%PropName=GetPropertyName(LatentHeat_)
                if (Me%PropName == 'Mean sea level pressure'    ) Me%PropName=GetPropertyName(AtmosphericPressure_)
                if (Me%PropName == 'Total cloud cover'          ) Me%PropName=GetPropertyName(CloudCover_)
                if (Me%PropName == '10 metre U wind component'  ) Me%PropName=GetPropertyName(WindVelocityX_)
                if (Me%PropName == '10 metre V wind component'  ) Me%PropName=GetPropertyName(WindVelocityY_)
                if (Me%PropName == '2 metre temperature'        ) Me%PropName=GetPropertyName(AirTemperature_)
                if (Me%PropName == 'East-West surface stress'   ) Me%PropName=GetPropertyName(WindStressX_)
                if (Me%PropName == 'North-South surface stress' ) Me%PropName=GetPropertyName(WindStressY_)
                if (Me%PropName == 'Surface solar radiation downwards') Me%PropName=GetPropertyName(SolarRadiation_)
                if (Me%PropName == 'Relative humidity') Me%PropName=GetPropertyName(RelativeHumidity_)
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
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleERA40Format - ERR01'

                if(Me%ConvertToASCII)then
                    call WriteDataLine(Me%OutputUnit,"<begin_data>")
                    write(Me%OutputUnit,*)Variable%Value%Float
                    call WriteDataLine(Me%OutputUnit,"<end_data>")
                endif


            case(1)

                nullify (Variable%Value%FLOAT_1D)
                allocate(Variable%Value%FLOAT_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_1D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleERA40Format - ERR02'

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
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleERA40Format - ERR03'
                
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
                        stop 'HandleFloats - ModuleERA40Format - ERR04'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array2D      = Variable%Value%FLOAT_2D,          &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleFloats - ModuleERA40Format - ERR05'

                end if

                deallocate(Variable%Value%FLOAT_2D)

            
            case(3)

                nullify (Variable%Value%FLOAT_3D)
                allocate(Variable%Value%FLOAT_3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,  &
                         Size%KLB:Size%KUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_3D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleERA40Format - ERR06'


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
                        stop 'HandleFloats - ModuleERA40Format - ERR07'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array3D      = Variable%Value%FLOAT_3D,          &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleFloats - ModuleERA40Format - ERR08'
                
                end if

                deallocate(Variable%Value%FLOAT_3D)


            case(4)

                nullify (Variable%Value%FLOAT_4D)
                allocate(Variable%Value%FLOAT_4D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,  &
                                                 Size%KLB:Size%KUB, Me%Time_LB:Me%Time_UB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%FLOAT_4D)
                if (stat /= nf90_noerr) stop 'HandleFloats - ModuleERA40Format - ERR09'


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
                            stop 'HandleFloats - ModuleERA40Format - ERR10'

                        Variable%Value%FLOAT_3D => Variable%Value%FLOAT_4D(:,:,:,t)

                        call HDF5WriteData(Me%ObjHDF5,                                  &
                                           "/Results/"//trim(Variable%ID%Name),         &
                                           trim(Variable%ID%Name),                      &
                                           Variable%Value%Units,                        &
                                           Array3D      = Variable%Value%FLOAT_3D,      &
                                           OutputNumber = t,                            &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleFloats - ModuleERA40Format - ERR10'
                        
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

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size


        select case(Variable%nDimensions)
            
            case(0)
                
                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Double)
                if (stat /= nf90_noerr) stop 'HandleDoubles - ModuleERA40Format - ERR01'

                if(Me%ConvertToASCII)then
                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Double

                    call WriteDataLine(Me%OutputUnit,"<end_data>")
                endif


            case(1)

                nullify (Variable%Value%Double_1D)
                allocate(Variable%Value%Double_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number,                      & 
                                    Variable%Value%Double_1D)
                if (stat /= nf90_noerr) stop 'HandleDoubles - ModuleERA40Format - ERR02'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                        write(Me%OutputUnit,*)Variable%Value%Double_1D(i)
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

                deallocate(Variable%Value%Double_1D)

            case(2)

                nullify (Variable%Value%Double_2D)
                allocate(Variable%Value%Double_2D(Size%ILB:Size%IUB, Size%JLB:Size%JUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number,                      & 
                                    Variable%Value%Double_2D)
                if (stat /= nf90_noerr) stop 'HandleDoubles - ModuleERA40Format - ERR03'
                
                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                        write(Me%OutputUnit,*)Variable%Value%Double_2D(i,j)
                    enddo
                    enddo

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


                if(Me%ConvertToHDF5)then

!                    call ConvertToMohidUnits2D(Field)

                    call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,        & 
                                       Size%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleDoubles - ModuleERA40Format - ERR04'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array2D      = Variable%Value%Double_2D,         &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleDoubles - ModuleERA40Format - ERR05'

                end if

                deallocate(Variable%Value%Double_2D)

            
            case(3)

                nullify (Variable%Value%Double_3D)
                allocate(Variable%Value%Double_3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB, & 
                         Size%KLB:Size%KUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number,                      & 
                                    Variable%Value%Double_3D)
                if (stat /= nf90_noerr) stop 'HandleDoubles - ModuleERA40Format - ERR06'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    do i = Size%ILB, Size%IUB
                    do j = Size%JLB, Size%JUB
                    do k = Size%KLB, Size%KUB
                        write(Me%OutputUnit,*)Variable%Value%Double_3D(i,j,k)
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
                        stop 'HandleDoubles - ModuleERA40Format - ERR07'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array3D      = Variable%Value%Double_3D,         &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleDoubles - ModuleERA40Format - ERR08'
                
                end if

                deallocate(Variable%Value%Double_3D)
            
            
            case(4)

                nullify (Variable%Value%Double_4D)
                allocate(Variable%Value%Double_4D(Size%ILB:Size%IUB, Size%JLB:Size%JUB, &
                                                  Size%KLB:Size%KUB,                    & 
                                                  Me%Time_LB:Me%Time_UB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number,                      & 
                                    Variable%Value%Double_4D)
                if (stat /= nf90_noerr) stop 'HandleDoubles - ModuleERA40Format - ERR09'


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



                    if(Me%ConvertToHDF5)then

                        call HDF5SetLimits(Me%ObjHDF5, Size%ILB, Size%IUB, Size%JLB,    & 
                                           Size%JUB, Size%KLB, Size%KUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleDoubles - ModuleERA40Format - ERR10'

                        Variable%Value%Double_3D => Variable%Value%Double_4D(:,:,:,t)

                        call HDF5WriteData(Me%ObjHDF5,                                  &
                                           "/Results/"//trim(Variable%ID%Name),         &
                                           trim(Variable%ID%Name),                      &
                                           Variable%Value%Units,                        &
                                           Array3D      = Variable%Value%Double_3D,     &
                                           OutputNumber = t,                            &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleDoubles - ModuleERA40Format - ERR11'
                        
                        nullify(Variable%Value%Double_4D)

                    end if

                    if(Me%ConvertToASCII)  call WriteDataLine(Me%OutputUnit,"<end_field>")

                end do
                
                if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit,"<end_data>")

                deallocate(Variable%Value%Double_4D)


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
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleERA40Format - ERR01'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Int

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif

            case(1)

                nullify (Variable%Value%Int_1D)
                allocate(Variable%Value%Int_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_1D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleERA40Format - ERR02'

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
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleERA40Format - ERR03'
                
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
                        stop 'HandleIntegers - ModuleERA40Format - ERR04'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array2D      = Variable%Value%Int_2D,            &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleIntegers - ModuleERA40Format - ERR05'

                end if

                deallocate(Variable%Value%Int_2D)

            
            case(3)

                nullify (Variable%Value%Int_3D)
                allocate(Variable%Value%Int_3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,    & 
                         Size%KLB:Size%KUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_3D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleERA40Format - ERR06'


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
                        stop 'HandleIntegers - ModuleERA40Format - ERR07'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//trim(Variable%ID%Name),             &
                                       trim(Variable%ID%Name),                          &
                                       Variable%Value%Units,                            &
                                       Array3D      = Variable%Value%Int_3D,            &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HandleIntegers - ModuleERA40Format - ERR08'
                
                end if


                deallocate(Variable%Value%Int_3D)

            case(4)

                nullify (Variable%Value%Int_4D)
                allocate(Variable%Value%Int_4D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,    &
                                               Size%KLB:Size%KUB, Me%Time_LB:Me%Time_UB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Int_4D)
                if (stat /= nf90_noerr) stop 'HandleIntegers - ModuleERA40Format - ERR09'


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
                            stop 'HandleIntegers - ModuleERA40Format - ERR10'

                        Variable%Value%Int_3D => Variable%Value%Int_4D(:,:,:,t)

                        call HDF5WriteData(Me%ObjHDF5,                                  &
                                           "/Results/"//trim(Variable%ID%Name),         &
                                           trim(Variable%ID%Name),                      &
                                           Variable%Value%Units,                        &
                                           Array3D      = Variable%Value%Int_3D,        &
                                           OutputNumber = t,                            &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'HandleIntegers - ModuleERA40Format - ERR10'
                        
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
        real, dimension(:,:), pointer                 :: aux

        !Begin-----------------------------------------------------------------

        Size = Variable%Dim%Size

        Me%PreviousDate = Me%InitialDate

        select case(Variable%nDimensions)
            
            case(0)
                
                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short)
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleERA40Format - ERR01'

                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Short

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


            case(1)

                nullify (Variable%Value%Short_1D)
                allocate(Variable%Value%Short_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short_1D)
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleERA40Format - ERR02'


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
                if (stat /= nf90_noerr) stop 'HandleShorts - ModuleERA40Format - ERR03'
                

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
                
                if(Me%ConvertToASCII)then

                    !get var
                    nullify (Variable%Value%Short_3D)
                    allocate(Variable%Value%Short_3D(Me%WorkSize%JLB:Me%WorkSize%JUB,       & 
                             Me%WorkSize%ILB:Me%WorkSize%IUB, Size%KLB:Size%KUB))
                    
                    stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short_3D)
                    if (stat /= nf90_noerr) stop 'HandleShorts - ModuleERA40Format - ERR04'
 
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

                        !get var for each k

                        nullify (Variable%Value%Short_2D)
                        allocate(Variable%Value%Short_2D(Me%WorkSize%JLB:Me%WorkSize%JUB,       & 
                                 Me%WorkSize%ILB:Me%WorkSize%IUB))

                        stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Short_2D, start=(/1,1,k/))
                        if (stat /= nf90_noerr) stop 'HandleShorts - ModuleERA40Format - ERR04'


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
                                real(Variable%Value%Short_2D(j,Me%WorkSize%IUB-i+1))  &
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

                            !change grid -180 180
                            if(Me%GridTo180)then

                                allocate(aux(Me%WorkSize%ILB:Me%WorkSize%IUB,  &
                                         Me%WorkSize%JLB:Me%WorkSize%JUB))

                                !write(*,*) 'Converting grid longitude from [0 360] to [-180 180]'
                                aux(:,:) =  Variable%Value%Float_2D(:,:)
                    
                                Variable%Value%Float_2D(:, 1:72) = aux(:, 73:144)
                                Variable%Value%Float_2D(:, 73:144) = aux(:, 1:72)

                                deallocate(aux);nullify(aux)

                            endif

                            !Write time
                            TimePtr => AuxTime

                            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleERA40Format - ERR05'

                            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                     Array1D = TimePtr,                                 &
                                     OutputNumber = Me%HDFOutputNumber, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleERA40Format - ERR06'

                            !Write values
                            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB,             & 
                                               Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                                               Me%WorkSize%JUB,                         &
                                           STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleERA40Format - ERR07'

                            call HDF5WriteData(Me%ObjHDF5,                              &
                                           "/Results/"//trim(Me%PropName),              &
                                           trim(Me%PropName), Variable%Value%Units,     &
                                           Array2D      = Variable%Value%Float_2D,      &
                                           OutputNumber = Me%HDFOutputNumber,           &
                                           STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                  &
                                stop 'HandleShorts - ModuleERA40Format - ERR08'

                            deallocate(Variable%Value%Float_2D)
                            deallocate(Variable%Value%Short_2D)

                            !Write instant time
                            write(*,*)'writing instant ', k

                            Me%PreviousDate = Me%CurrentDate
                            Me%HDFOutputNumber = Me%HDFOutputNumber + 1

                            endif

                        endif

                    enddo

                endif


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
                if (stat /= nf90_noerr) stop 'HandleChars - ModuleERA40Format - ERR01'


                if(Me%ConvertToASCII)then

                    call WriteDataLine(Me%OutputUnit,"<begin_data>")

                    write(Me%OutputUnit,*)Variable%Value%Char

                    call WriteDataLine(Me%OutputUnit,"<end_data>")

                endif


            case(1)

                nullify (Variable%Value%Char_1D)
                allocate(Variable%Value%Char_1D(Size%ILB:Size%IUB))

                stat = nf90_get_var(Me%FileID, Variable%ID%Number, Variable%Value%Char_1D)
                if (stat /= nf90_noerr) stop 'HandleChars - ModuleERA40Format - ERR02'

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
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR01'


        call GetData(Me%YY_Variable,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'YY_VARIABLE',                  &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR02'


        call GetData(Me%ZZ_Variable,                                &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'ZZ_VARIABLE',                  &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR03'


        call GetData(Me%Time_Variable,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'TIME_VARIABLE',                &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR04'


        call GetData(Me%GridType_Variable,                          &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'GRID_TYPE_VARIABLE',           &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR05'

        call GetData(Me%VertCoord_Variable,                         &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'VERT_COORD_VARIABLE',          &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR06'

        call GetData(Me%Batim_Variable,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromFile,                       &
                     keyword      = 'BATIM_VARIABLE',               &
                     ClientModule = 'ModuleERA40Format',            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGridMapOptions - ModuleERA40Format - ERR07'

    end subroutine ReadGridMapOptions

    !------------------------------------------------------------------
    
    subroutine ConstructTime

        !Local---------------------------------------------------------          
        type(T_Variable), pointer           :: Time
        integer                             :: stat, t, t1, t2
        integer                             :: ChangeDateOption 
        real                                :: Year, Month, Day, Hour, Minute, Second
        real                                :: NYear, NMonth, NDay, NHour 
        real                                :: NMinute, NSecond

        !Begin---------------------------------------------------------          

        nullify(Time); allocate(Time)

        Time%ID%Name = Me%Time_Variable

        stat = nf90_inq_varid(Me%FileID, trim(Time%ID%Name), Time%ID%Number)
        if (stat /= nf90_noerr) return

        stat = nf90_inquire_variable(Me%FileID, Time%ID%Number, Time%ID%Name,       &
                                     Time%ID%Type_, Time%nDimensions,               &
                                     Time%DimensionsID, Time%nAttributes)
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleERA40Format - ERR01'

        !if(Me%ConvertToASCII) call WriteDataLine(Me%OutputUnit, "<begin_time>")

        call HandleDimensions (Time, .true.)

        call HandleAttributtes(Time%ID%Number,  Time%nAttributes)

        nullify(Me%TimeValues);allocate(Me%TimeValues(Time%Dim%Size%ILB:Time%Dim%Size%IUB))
        nullify(Me%TimeArray );allocate(Me%TimeArray (Time%Dim%Size%ILB:Time%Dim%Size%IUB))


        stat = nf90_get_var(Me%FileID, Time%ID%Number, Me%TimeValues)
        if (stat /= nf90_noerr) stop 'ConstructTime - ModuleERA40Format - ERR02'


        Me%Time_LB = Time%Dim%Size%ILB
        Me%Time_UB = Time%Dim%Size%IUB

        call SetDate(Me%InitialDate, 1900.,1.,1.,0.,0.,0.) !This previously known

        do t = Me%Time_LB, Me%Time_UB
            
            Me%TimeArray(t) = Me%InitialDate + Me%TimeValues(t) * 3600.

        enddo

        !Check if there are equal times
        do t1 = Me%Time_LB, Me%Time_UB
        do t2 = Me%Time_LB, Me%Time_UB
            if ((Me%TimeArray(t1) == Me%TimeArray(t2)) .AND. (t1 /= t2)) then
890             format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
                write(*,*)'Repeated date found.' 
                write(*,*)'If not altered, first detected data will be used'
                write(*,*)'and second detected data will not be used.'
                write(*,*)'Check the dates:'                               
                write(*,*)'Last date - 1:'
                if (t1 /= Me%Time_LB) then
                    call ExtractDate(Me%TimeArray(t1-1), Year, Month, Day, Hour,    &  
                                     Minute, Second)
                    write(*,fmt=890) Year, Month, Day, Hour, Minute, Second
                else
                    write(*,*)'No last date. Repeated date is first date.'
                endif               
                write(*,*)'Repeated date - 2:'
                call ExtractDate(Me%TimeArray(t1), Year, Month, Day, Hour,          &  
                                 Minute, Second)
                write(*,fmt=890) Year, Month, Day, Hour, Minute, Second
                write(*,*)'Repeated date - 3:'
                call ExtractDate(Me%TimeArray(t2), Year, Month, Day, Hour,          &  
                                 Minute, Second)
                write(*,fmt=890) Year, Month, Day, Hour, Minute, Second
                write(*,*)'Next date (after date 2) - 4:'
                if (t1 /= Me%Time_UB) then
                    call ExtractDate(Me%TimeArray(t1+1), Year, Month, Day, Hour,    &  
                                     Minute, Second)
                    write(*,fmt=890) Year, Month, Day, Hour, Minute, Second
                else
                    write(*,*)'No next date. Repeated date is last date.'
                endif
                write(*,*)'Do you want to change any of the dates?'
                write(*,*)'(Answer 0 for not change OR answer number of date)'
                read(*, '(i1)') ChangeDateOption                                         
                if (ChangeDateOption /= 0) then
                    write(*,*)'Year?'
                    read(*, '(f5.0)') NYear
                    write(*,*)'Month?'
                    read(*, '(f3.0)') NMonth
                    write(*,*)'Day?'
                    read(*, '(f3.0)') NDay
                    write(*,*)'Hour?'
                    read(*, '(f3.0)') NHour
                    write(*,*)'Minute?'
                    read(*, '(f3.0)') NMinute
                    write(*,*)'Second?'
                    read(*, '(f3.0)') NSecond                                                           
                    select case(ChangeDateOption)

                        case(1)
                            call SetDate(Me%TimeArray(t1-1),NYear,NMonth,NDay,      &
                                         NHour,NMinute,NSecond)

                        case(2)
                            call SetDate(Me%TimeArray(t1),NYear,NMonth,NDay,NHour,  &
                                         NMinute,NSecond)

                        case(3)
                            call SetDate(Me%TimeArray(t2),NYear,NMonth,NDay,NHour,  &
                                         NMinute,NSecond)

                        case(4)
                            call SetDate(Me%TimeArray(t1+1),NYear,NMonth,NDay,      &
                                         NHour,NMinute,NSecond)

                        case default

                            write(*,*)'Not a valid answer. No date changed.'

                    end select
                else
                    write(*,*)'No date changed.'
                endif

            endif
        enddo
        enddo

        deallocate(Time)

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
            if (stat /= nf90_noerr) stop 'ConstructBathymetry - ModuleERA40Format - ERR01'

            stat = nf90_inquire_variable(Me%FileID, Me%Batim%ID%Number,                 & 
                                         Me%Batim%ID%Name,                              &
                                         Me%Batim%ID%Type_, Me%Batim%nDimensions,       &
                                         Me%Batim%DimensionsID, Me%Batim%nAttributes)
            if (stat /= nf90_noerr) stop 'ConstructBathymetry - ModuleERA40Format - ERR02'

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
            if (stat /= nf90_noerr) stop 'ConstructBathymetry - ModuleERA40Format - ERR03'


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
        type(T_Variable), pointer                       :: XX, YY, GridType
        integer                                         :: stat
        integer                                         :: i, j
        real                                            :: dx, dy, XOrigin, YOrigin

        !Begin---------------------------------------------------------          

        nullify(XX, YY, GridType); allocate(XX, YY, GridType)

        call GetVariable(XX, trim(Me%XX_Variable))
        call GetVariable(YY, trim(Me%YY_Variable))


        stat = nf90_inq_varid(Me%FileID, trim(Me%GridType_Variable), GridType%ID%Number)
        if (stat /= nf90_noerr) then
            write(*,*)'No grid type defined. Assumed cartesian...'
            write(*,*)
            Me%Grid_Type = Cartesian
        else
            stat = nf90_get_var(Me%FileID, GridType%ID%Number, GridType%Value%Int)
            if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleERA40Format - ERR01'

            Me%Grid_Type = GridType%Value%Int
        end if
        
        select case(Me%Grid_Type)

            case(Cartesian)

                write(*,*)'Grid type: Cartesian'

                Me%WorkSize%ILB = 1
                Me%WorkSize%JLB = 1
                Me%WorkSize%IUB = YY%Dim%Size%IUB
                Me%WorkSize%JUB = XX%Dim%Size%IUB


                Me%Size%ILB = Me%WorkSize%ILB - 1
                Me%Size%JLB = Me%WorkSize%JLB - 1
                Me%Size%IUB = Me%WorkSize%IUB + 1
                Me%Size%JUB = Me%WorkSize%JUB + 1

                nullify(Me%ConnectionX, Me%ConnectionY)

                allocate(Me%ConnectionX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                allocate(Me%ConnectionY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                
                nullify (XX%Value%FLOAT_1D, YY%Value%FLOAT_1D)
                allocate(XX%Value%FLOAT_1D(Me%Size%JLB:Me%Size%JUB))
                allocate(YY%Value%FLOAT_1D(Me%Size%ILB:Me%Size%IUB))

                stat = nf90_get_var(Me%FileID, XX%ID%Number,                            & 
                                    XX%Value%Float_1D(Me%WorkSize%JLB:Me%WorkSize%JUB))
                if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleERA40Format - ERR02'

                stat = nf90_get_var(Me%FileID, YY%ID%Number,                            & 
                                    YY%Value%Float_1D(Me%WorkSize%ILB:Me%WorkSize%IUB))
                if (stat /= nf90_noerr) stop 'ConstructGrid - ModuleERA40Format - ERR03'


                dx = XX%Value%Float_1D(2) - XX%Value%Float_1D(1)
                dy = YY%Value%Float_1D(1) - YY%Value%Float_1D(2)

                if  (Me%GridTo180)then
                !put -180 to 180
                    XOrigin = XX%Value%Float_1D(1) - 180.
                    else
                    XOrigin = XX%Value%Float_1D(1) - dx/2.
                end if
                YOrigin = YY%Value%Float_1D(YY%Dim%Size%IUB) - dy/2.

                if  (Me%GridTo180)then
                !put -180 to 180
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        Me%ConnectionX(i,j) = XX%Value%Float_1D(j  )                        &                     
                                              - dx/2. + XOrigin
                        Me%ConnectionY(i,j) = YY%Value%Float_1D(YY%Dim%Size%IUB - i + 1)    &
                                              - dy/2.

                    enddo
                    enddo
                else
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        Me%ConnectionX(i,j) = XX%Value%Float_1D(j  )                        &                     
                                              - dx/2.
                        Me%ConnectionY(i,j) = YY%Value%Float_1D(YY%Dim%Size%IUB - i + 1)    &
                                              - dy/2.

                    enddo
                    enddo
                 end if

                Me%ConnectionY(Me%WorkSize%IUB + 1, :) =                                &
                                          Me%ConnectionY(Me%WorkSize%IUB, :) + dy
                Me%ConnectionX(Me%WorkSize%IUB + 1, :) =                                &
                                          Me%ConnectionX(Me%WorkSize%IUB, :)

                Me%ConnectionX(:, Me%WorkSize%JUB + 1) = Me%ConnectionX(:, Me%WorkSize%JUB) + dx
                Me%ConnectionY(:, Me%WorkSize%JUB + 1) = Me%ConnectionY(:, Me%WorkSize%JUB)

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

        !By default only one layer is read
        Me%WorkSize%KLB = 1
        Me%Size%KLB     = 0
        Me%WorkSize%KUB = 1
        Me%Size%KUB     = 2


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
        if (stat /= nf90_noerr) stop 'GetVariable - ModuleERA40Format - ERR01'

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

end module ModuleERA40Format









