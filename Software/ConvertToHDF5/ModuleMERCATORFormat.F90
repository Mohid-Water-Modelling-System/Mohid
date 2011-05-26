!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : MERCATORFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert MERCATORFormat files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleMERCATORFormat

    use ModuleTime
    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap

    ! Manages NetCDF files
#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif
    implicit none

    private 
    
    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertMERCATORFormat
    private ::      ReadOptions
    private ::      OpenAndReadMERCATORFileV1
    private ::          AddField
    private ::          AddDate
    private ::          CheckName
    private ::      OpenAndReadMERCATORFileV2
    private ::          OpenAndReadBathymMERCATOR
    private ::          ConstructGridV2
    private ::          WriteMERCATORGeometry
    private ::          Open_HDF5_OutPut_File
    private ::          OpenAndReadMERCATORFields
    private ::      OpenAndReadMERCATORFileV3
    private ::          OpenAndReadBathymMERCATORV3
    private ::          OpenAndReadMERCATORFieldsV3
    private ::              ReadMercatorFileV3
    private ::                  CheckNameV3
    private ::      OutputFields
    private ::          ReadMercatorFile
    private ::          MapFromMercatorFile
    private ::      WriteVelocityModulus
    private ::      KillMERCATORFormat

    !Parameters----------------------------------------------------------------
    integer, parameter                          :: Version1     = 1
    integer, parameter                          :: Version2     = 2
    integer, parameter                          :: Version3     = 3
    integer, parameter                          :: Version4     = 4
    integer, parameter                          :: PSY2V4       = 5
    integer, parameter                          :: Version6     = 6
    
    
    integer, parameter                          :: MercatorLayers = 43

    character(LEN = StringLength), parameter    :: input_files_begin   = '<<begin_input_files>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<end_input_files>>'

    real,    parameter                          :: MissingValue = 99999.

    !Types---------------------------------------------------------------------
    
    private :: T_Date
    type       T_Date
        type(T_Time)                            :: Date
        type(T_Date), pointer                   :: Next
    end type  T_Date


    private :: T_Field
    type       T_Field
        character(len=StringLength)             :: Name
        character(len=StringLength)             :: Units
        integer                                 :: GridLocation
        type(T_Time)                            :: Date
        integer                                 :: nDimensions
        real, dimension(:,:),       pointer     :: Values2D
        real, dimension(:,:,:),     pointer     :: Values3D
        integer                                 :: OutputNumber         = 1
        type(T_Size3D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: Next
    end type  T_Field

    
    private :: T_MERCATORFormat
    type       T_MERCATORFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjGridData          = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjGeometry          = 0
        integer                                 :: ObjMap               = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit, ClientNumber
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName, InputGridFile
        character(len=PathLength)               :: InputGridFileU, InputGridFileV
        character(len=PathLength)               :: InputBathymetryFile
        character(len=PathLength)               :: InputMeshZGridFile
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: BaseBulletin
        character(len=PathLength)               :: DatesFile
        character(len=PathLength)               :: GeometryFileName
        integer                                 :: NumDates, LayerNumber
        integer                                 :: ReadOptionType
        integer                                 :: imax, jmax, kmax
        real, dimension(:,:),       pointer     :: Bathymetry, BathymetryMax
        real, dimension(:,:),       pointer     :: CenterX, CenterY
        real, dimension(:,:),       pointer     :: XX_IE, YY_IE
        real, dimension(:  ),       pointer     :: XX, YY, LayersThickness
        real, dimension(:,:,:),     pointer     :: SZZ
        real, dimension(:, :, :),   pointer     :: AreaU
        real, dimension(:, :, :),   pointer     :: AreaV
        real, dimension(:, :) ,     pointer     :: SurfaceElevation
        real, dimension(:, :, :),   pointer     :: ScaFactorU, ScaFactorV
        integer, dimension(:,:,:),  pointer     :: WaterPoints3D 
        integer, dimension(:,:,:),  pointer     :: ComputeFaces3DU, ComputeFaces3DV 
        integer, dimension(:,:,:),  pointer     :: ComputeFacesU3D
        integer, dimension(:,:,:),  pointer     :: ComputeFacesV3D
        integer, dimension(:,:,:),  pointer     :: ComputeFacesW3D
        type(T_Size3D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: FirstField         
        type(T_Date),               pointer     :: FirstDate 
        type(T_Time)                            :: RefDateTime    
        integer, dimension(12)                  :: Instants(1:12) = 0
        logical                                 :: ComputeBarotropicVel = .false.
    end type  T_MERCATORFormat

    type(T_MERCATORFormat), pointer             :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertMERCATORFormat(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        integer                                 STAT_CALL

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

!        !The time in Mercator is compute in days from 1950/1/1 : 0h:0m:0s
!        call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

!        call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0,    &
!                                 VariableDT = .false., STAT = STAT_CALL)   
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMERCATORFormat - ERR02a'


        call ReadOptions

        if      (Me%ReadOptionType == Version1) then

            !The time in Mercator is compute in days from 1950/1/1 : 0h:0m:0s
            call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

            call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0, &
                                     VariableDT = .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConvertMERCATORFormat - ModuleMERCATORFormat - ERR01'


            call OpenAndReadMERCATORFileV1

            call ConstructGrid

            call OutputFields

        else if (Me%ReadOptionType == Version2) then

            !The time in Mercator is compute in days from 1950/1/1 : 0h:0m:0s
            call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

            call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0, &
                                     VariableDT = .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConvertMERCATORFormat - ModuleMERCATORFormat - ERR02'


            call OpenAndReadMERCATORFileV2

        else if (Me%ReadOptionType == Version3) then

            !The time in Mercator is compute in seconds from 2006/10/11 : 0h:0m:0s
            call SetDate (Me%RefDateTime, Year=2006, Month=10, Day=11, Hour=0, Minute=0, Second=0) 
            !call SetDate (Me%RefDateTime, Year=2006, Month=10, Day=10, Hour=0, Minute=0, Second=0) 

            call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0, &
                                     VariableDT = .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConvertMERCATORFormat - ModuleMERCATORFormat - ERR03'


            call OpenAndReadMERCATORFileV3

        else if (Me%ReadOptionType == Version4) then

            !The time in Mercator is compute in days from 1950/1/1 : 0h:0m:0s
            call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

            call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0, &
                                  VariableDT = .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConvertMERCATORFormat - ModuleMERCATORFormat - ERR04'

            call OpenAndReadMERCATORFileV4
            
!MJ ************************************************************************
        else if (Me%ReadOptionType == Version6) then

            !The time in Mercator is compute in days from 1950/1/1 : 0h:0m:0s
            call SetDate (Me%RefDateTime, Year=1950, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

            call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0, &
                                  VariableDT = .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConvertMERCATORFormat - ModuleMERCATORFormat - ERR04'

            call OpenAndReadMERCATORFileV6

        else if (Me%ReadOptionType == PSY2V4) then

            !The time in Mercator is compute in seconds from 2006/10/11 : 0h:0m:0s
            call SetDate (Me%RefDateTime, Year=2006, Month=10, Day=11, Hour=0, Minute=0, Second=0) 

            call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0, &
                                  VariableDT = .false., STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConvertMERCATORFormat - ModuleMERCATORFormat - ERR05'

            call OpenAndReadMERCATORFileV5

        endif

        call WriteVelocityModulus

        call KillMERCATORFormat

        STAT = SUCCESS_

    end subroutine ConvertMERCATORFormat

    !------------------------------------------------------------------------

    subroutine WriteVelocityModulus

        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:  ), pointer         :: Aux3D, Aux3D_U, Aux3D_V
        integer                                     :: STAT_CALL, HDF5_READWRITE
        integer                                     :: iOut, nItems
        character(Len=StringLength)                 :: PropUnits, MohidName
        integer                                     :: i, j, k
        !Begin-----------------------------------------------------------------

        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR60'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR01'

        allocate(Aux3D  (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
        allocate(Aux3D_U(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
        allocate(Aux3D_V(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))        

        MohidName = GetPropertyName(VelocityU_)

        call GetHDF5GroupNumberOfItems (Me%ObjHDF5, "/Time", nItems, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR10'

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                             Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR15'


        do iOut = 1, nItems

            MohidName = GetPropertyName(VelocityU_)

            call HDF5ReadData  (Me%ObjHDF5, "/Results/"//trim(MohidName),               &
                                 trim(MohidName),Array3D = Aux3D_U,    &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR20'

            MohidName = GetPropertyName(VelocityV_)

            call HDF5ReadData  (Me%ObjHDF5, "/Results/"//trim(MohidName),               &
                                 trim(MohidName),Array3D = Aux3D_V,    &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR30'

            Aux3D(:,:,:) = 0.0

            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if(Me%WaterPoints3D(i,j,k) == 1)then
                    Aux3D(i,j,k) = sqrt(Aux3D_U(i,j,k)**2 + Aux3D_V(i,j,k)**2)
                endif

            enddo
            enddo
            enddo


            MohidName = GetPropertyName(VelocityModulus_)
            PropUnits = MohidUnits(MohidName)


            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array3D = Aux3D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR40'
        enddo

        deallocate(Aux3D, Aux3D_U, Aux3D_V)
        nullify   (Aux3D, Aux3D_U, Aux3D_V)

    end subroutine WriteVelocityModulus
    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                             &
                     ClientModule = 'ModuleMERCATORFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR10'


        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleMERCATORFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR20'

        call GetData(Me%ReadOptionType,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'READ_OPTION',                                      &
                     default      = Version3,                                           &
                     ClientModule = 'ModuleMERCATORFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR30'

        if (Me%ReadOptionType == Version1) then

            call GetData(Me%BaseBulletin,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BASE_BULLETIN',                                &
                         ClientModule = 'ModuleMERCATORFormat',                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR40'


            call GetData(Me%DatesFile,                                                  &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'DATES_FILE',                                   &
                         ClientModule = 'ModuleMERCATORFormat',                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR50'

            call GetData(Me%NumDates,                                                   &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NUM_DATES',                                    &
                         ClientModule = 'ModuleMERCATORFormat',                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR60'

        else

            call GetData(Me%GeometryFilename,                                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'OUTPUT_GEOMETRY_FILENAME',                     &
                         ClientModule = 'ModuleMERCATORFormat',                         &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR70'


            if (Me%ReadOptionType /= Version4 .and. Me%ReadOptionType /= Version6 ) then            

                !Version 2 and Version 3
                call GetData(Me%InputGridFile,                                              &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'INPUT_GRID_FILENAME',                          &
                             ClientModule = 'ModuleMERCATORFormat',                         &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR80'

                if (Me%ReadOptionType == Version2) then
            
                    call GetData(Me%InputGridFileU,                                         &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'INPUT_GRID_FILENAME_U',                    &
                                 ClientModule = 'ModuleMERCATORFormat',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR90'

                    call GetData(Me%InputGridFileV,                                         &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'INPUT_GRID_FILENAME_V',                    &
                                 ClientModule = 'ModuleMERCATORFormat',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR100'

                else

                    call GetData(Me%InputBathymetryFile,                                    &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'INPUT_BATHY_FILENAME',                     &
                                 ClientModule = 'ModuleMERCATORFormat',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR110'

                    !Version 3 and Version 4
                    call GetData(Me%ComputeBarotropicVel,                                   &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      = 'CALC_BAROTROPIC_VEL',                              &
                         default      = .false.,                                              &
                         ClientModule = 'ModuleMERCATORFormat',                             &
                         STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR120'

                    if (Me%ComputeBarotropicVel) then

                        call GetData(Me%InputMeshZGridFile,                                 &
                             Me%ObjEnterData, iflag,                                        &
                             SearchType   = FromBlock,                                      &
                             keyword      = 'INPUT_MESH_ZGRID_FILENAME',                    &
                             ClientModule = 'ModuleMERCATORFormat',                         &
                             STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR130'

                    endif

                endif

            else

                call GetData(Me%ComputeBarotropicVel,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'CALC_BAROTROPIC_VEL',                              &
                     default      = .false.,                                              &
                     ClientModule = 'ModuleMERCATORFormat',                             &
                     STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMERCATORFormat - ERR120'

            endif

        endif

    end subroutine ReadOptions


    !--------------------------------------------------------------------------

    subroutine OpenAndReadMERCATORFileV1

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist

        !Others----------------------------------------------------------------
        integer                                     :: nDimensions
        character (len=19)                          :: current_date
        character (len=25)                          :: units
        real, pointer,     dimension(:,:,:,:)       :: DataAux
        real,              dimension(6)             :: TimeReal
        integer                                     :: i, j
        logical                                     :: newtime = .TRUE.
        character(len=StringLength)                 :: MohidName
        type(T_Field), pointer                      :: NewField
        type(T_Date ), pointer                      :: NewDate
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB

        !Variables Pablo
        integer                                     :: n, k, nDates
        integer                                     :: ncid, dimid, status
        integer, dimension(4)                       :: dimids, lengths
        integer                                     :: nDims, nVars, nAtrr, xtype
        character (len=80)                          :: nameAux, nameAux2
        character (len=12)                          :: prediction
        real, pointer,     dimension(:,:,:)         :: tmp

        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading MERCATOR output file...'
        
        nullify(NewField        )
        nullify(Me%FirstField   )
        nullify(Me%FirstDate    )
        nullify(Me%Bathymetry   )
        nullify(Me%BathymetryMax)
        nullify(Me%CenterX      )
        nullify(Me%CenterY      )
        nullify(Me%XX           )
        nullify(Me%YY           )
        nullify(Me%SZZ          )
        nullify(Me%WaterPoints3D)


        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR00a'

        open(Unit   = Me%Unit,          &
             File   = Me%DatesFile,     &
             Form   = 'FORMATTED',      &
             STATUS = 'OLD',            &
             Action = 'READ',           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR00b'
        
        rewind(Me%Unit)

main:   do nDates=1,Me%NumDates

        read(Me%Unit,*) prediction

        !Verifies if file exists
        Me%FileName = "mercator_psy2v1_prestige_grid1o12_dc_galice_tsuv_b"//trim(Me%BaseBulletin)//"_f"//trim(prediction)//"_h.nc"
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'MERCATOR File does not exist'
            stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR01'
        endif

        status=NF90_OPEN(Me%FileName,NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR02'

                if (newtime) then

                    status=NF90_GET_ATT(ncid,NF90_GLOBAL,"field_date",current_date)
                    if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR02b'

                    do i=1,len_trim(current_date)

                        if (current_date(i:i) =='_'.or.current_date(i:i) ==':'.or. current_date(i:i) =='-') then
                            current_date(i:i) = ' '
                        endif

                    enddo

                    read(current_date,*) TimeReal
                    call AddDate(Me%FirstDate, NewDate)
                    call SetDate(NewDate%Date, Year = TimeReal(1), Month  = TimeReal(2), Day    = TimeReal(3), &
                                               Hour = TimeReal(4), Minute = TimeReal(5), Second = TimeReal(6))

                endif


        Me%WorkSize%ILB = 1

        status=NF90_INQ_DIMID(ncid,"latitude",dimid)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR03a'

        status=NF90_INQUIRE_DIMENSION(ncid,dimid,nameAux,Me%WorkSize%IUB)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR03a'

        Me%WorkSize%JLB = 1

        status=NF90_INQ_DIMID(ncid,"longitude",dimid)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR03b'

        status=NF90_INQUIRE_DIMENSION(ncid,dimid,nameAux,Me%WorkSize%JUB)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR03b'

        Me%Size%ILB     = Me%WorkSize%ILB - 1
        Me%Size%IUB     = Me%WorkSize%IUB + 1
        Me%Size%JLB     = Me%WorkSize%JLB - 1
        Me%Size%JUB     = Me%WorkSize%JUB + 1

        WILB            = Me%WorkSize%ILB 
        WIUB            = Me%WorkSize%IUB 
        WJLB            = Me%WorkSize%JLB 
        WJUB            = Me%WorkSize%JUB 

        ILB             = Me%Size%ILB
        IUB             = Me%Size%IUB
        JLB             = Me%Size%JLB
        JUB             = Me%Size%JUB


        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04'

do0:    do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions, dimids)
            if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04b'

            do j=1,nDimensions
                status=NF90_INQUIRE_DIMENSION(ncid,dimids(j),nameAux2,lengths(j))
                if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04c'
            end do

            status=NF90_GET_ATT(ncid,n,"units",units)
            if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04d'

            if (nDimensions == 1) then
                allocate(DataAux(lengths(1), 1, 1, 1))
            elseif (nDimensions == 2) then
                allocate(DataAux(lengths(1), lengths(2), 1, 1))
            elseif (nDimensions == 3) then
                allocate(DataAux(lengths(1), lengths(2), lengths(3), 1))
            endif

            status=NF90_GET_VAR(ncid,n,DataAux)
            if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04e'

            select case(nameAux)

                        case('depth' )
                            
                            if(.not. associated(Me%Bathymetry))then

!------------------------------------------------------------------------------------------------
!           Revisar esto. Es un procedimiento bastante chapucero para obtener la bati.
!------------------------------------------------------------------------------------------------
                                allocate(Me%WaterPoints3D(WIUB, WJUB, MercatorLayers))
                                allocate(tmp(WJUB, WIUB, MercatorLayers)             )

                                status=NF90_INQ_VARID(ncid,"temperature",dimid)
                                if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR03a'

                                status=NF90_GET_VAR(ncid,dimid,tmp)
                                if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04d'

                                do k=1,MercatorLayers
                                do i=WILB,WIUB
                                do j=WJLB,WJUB
                                if (tmp(j,i,k)==1e35) then
                                    Me%WaterPoints3D(i,j,44-k) = 0
                                    tmp(j,i,k) = 0
                                else
                                    tmp(j,i,k) = 1
                                    Me%WaterPoints3D(i,j,44-k) = 1
                                end if
                                enddo
                                enddo
                                enddo


                                do i=WILB,WIUB
                                do j=WJLB,WJUB
                                do k=1,42
                                    tmp(j, i, k)=tmp(j, i, k)*(DataAux(k,1,1,1)+DataAux(k+1,1,1,1))*0.5
                                enddo
                                tmp(j, i, MercatorLayers)=(tmp(j, i, 42)+302.1)*tmp(j, i, MercatorLayers)
                                enddo
                                enddo
                                

                                allocate(Me%Bathymetry(ILB:IUB, JLB:JUB))
                            
                                do i=WILB,WIUB
                                do j=WJLB,WJUB
                                if(maxval(tmp(j, i, 1:MercatorLayers))==0) then
                                    Me%Bathymetry(i,j) = -99.
                                else
                                    Me%Bathymetry(i,j) = maxval(tmp(j, i, 1:MercatorLayers))
                                    !One centimer precision for the bathymetry
                                    Me%Bathymetry(i,j) = real(int(Me%Bathymetry(i,j)*100.))/100. 
                                endif
                                enddo
                                enddo

                                deallocate(tmp)
!------------------------------------------------------------------------------------------------

                            end if

                        case('longitude')

                            if(.not. associated(Me%CenterX))then
                                
                                allocate(Me%CenterX(ILB:IUB, JLB:JUB))
                                allocate(Me%XX              (JLB:JUB))

                                do i = WILB, WIUB
                                    Me%CenterX(i, WJLB:WJUB) = DataAux(WJLB:WJUB, 1, 1, 1)
                                enddo

                                Me%XX(WJLB) = 0.
                                
                                do j = WJLB, WJUB - 1
                                    
                                    Me%XX(j + 1) = Me%XX(j) + (Me%CenterX(1, j + 1) - Me%CenterX(1, j))
                                
                                end do

                                Me%XX(WJUB + 1)  = Me%XX(WJUB) + (Me%CenterX(1, WJUB) - Me%CenterX(1, WJUB - 1))

                            end if

                        case('latitude')
                            
                            if(.not. associated(Me%CenterY))then
                                
                                allocate(Me%CenterY(ILB:IUB, JLB:JUB))
                                allocate(Me%YY     (ILB:IUB         ))

                                do j=WILB,WIUB
                                Me%CenterY(WILB:WIUB, j) = DataAux(WILB:WIUB, 1, 1, 1)
                                enddo
                                               
                                Me%YY(WILB) = 0.
                                
                                do i = WILB, WIUB - 1
                                    
                                    Me%YY(i + 1) = Me%YY(i) + (Me%CenterY(i + 1, 1) - Me%CenterY(i, 1))
                                
                                end do
                                
                                Me%YY(WIUB + 1)  = Me%YY(WIUB) + (Me%CenterY(WIUB, 1) - Me%CenterY(WIUB - 1, 1))

                            end if

                        case default

                if(CheckName(nameAux, MohidName))then

                    call AddField(Me%FirstField, NewField)

                    NewField%WorkSize%ILB   = 1
                    NewField%WorkSize%IUB   = lengths(2)
                    NewField%WorkSize%JLB   = 1
                    NewField%WorkSize%JUB   = lengths(1)
                    NewField%WorkSize%KLB   = 1
                    NewField%WorkSize%KUB   = lengths(3)

                    NewField%Size%ILB       = NewField%WorkSize%ILB - 1
                    NewField%Size%IUB       = NewField%WorkSize%IUB + 1
                    NewField%Size%JLB       = NewField%WorkSize%JLB - 1
                    NewField%Size%JUB       = NewField%WorkSize%JUB + 1
                    NewField%Size%KLB       = NewField%WorkSize%KLB - 1
                    NewField%Size%KUB       = NewField%WorkSize%KUB + 1

                    WILB                    = NewField%WorkSize%ILB 
                    WIUB                    = NewField%WorkSize%IUB 
                    WJLB                    = NewField%WorkSize%JLB 
                    WJUB                    = NewField%WorkSize%JUB 
                    WKLB                    = NewField%WorkSize%KLB 
                    WKUB                    = NewField%WorkSize%KUB 
                                            
                    ILB                     = NewField%Size%ILB
                    IUB                     = NewField%Size%IUB
                    JLB                     = NewField%Size%JLB
                    JUB                     = NewField%Size%JUB
                    KLB                     = NewField%Size%KLB
                    KUB                     = NewField%Size%KUB

                    NewField%Name        =   trim(MohidName)

                    ! Data para cada campo:
                    status=NF90_GET_ATT(ncid,NF90_GLOBAL,"field_date",current_date)
                    if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR02b'

                    do i=1,len_trim(current_date)

                        if (current_date(i:i) =='_'.or.current_date(i:i) ==':'.or. current_date(i:i) =='-') then
                            current_date(i:i) = ' '
                        endif

                    enddo

                    read(current_date,*) TimeReal
                    call SetDate(NewField%Date, Year = TimeReal(1), Month  = TimeReal(2), Day    = TimeReal(3), &
                                                Hour = TimeReal(4), Minute = TimeReal(5), Second = TimeReal(6))

                    NewField%Units       =   trim(Units)
                    NewField%nDimensions =   nDimensions
                    
                    select case(nDimensions)

                        case(2)
                                allocate(NewField%Values2D(ILB:IUB, JLB:JUB))

                                ! Mudamos os indices (por causa do NetCDF):
                                do i=WILB,WIUB
                                do j=WJLB,WJUB
                                    NewField%Values2D(i,j) = DataAux(j,i,1, 1)
                                enddo
                                enddo
                    
                        case(3)

                                allocate(NewField%Values3D(ILB:IUB, JLB:JUB,KLB:KUB))

                                ! Mudamos os indices (por causa do NetCDF):
                                do i=WILB,WIUB
                                do j=WJLB,WJUB
                                do k=WKLB,WKUB
                                    NewField%Values3D(i,j,WKUB-k+1) = DataAux(j,i,k, 1)
                                enddo
                                enddo
                                enddo

                        case default

                    end select

                    call ConvertToMohidUnits(NewField)

                end if

            end select

            deallocate(DataAux)

        end do do0

        status=NF90_CLOSE(ncid)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04'

        end do main

100     write(*,*)'Finished reading MERCATOR output file.'

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR05'

    end subroutine OpenAndReadMERCATORFileV1
    
    
    !--------------------------------------------------------------------------

    subroutine OpenAndReadMERCATORFileV2

        !Local-----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        call OpenAndReadBathymMERCATOR

        call ConstructGridV2

        call WriteMERCATORGeometry

        call Open_HDF5_OutPut_File

        call OpenAndReadMERCATORFields

    end subroutine OpenAndReadMERCATORFileV2
    
    
    !--------------------------------------------------------------------------

    subroutine OpenAndReadMERCATORFileV3

        !Local-----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        call OpenAndReadBathymMERCATORV3

        call ConstructGridV2

        call WriteMERCATORGeometry

        call Open_HDF5_OutPut_File

        call OpenAndReadMERCATORFieldsV3

    end subroutine OpenAndReadMERCATORFileV3
    
        !--------------------------------------------------------------------------

    subroutine OpenAndReadMERCATORFileV4

        !Local-----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        call OpenAndReadBathymMERCATORV4

        call ConstructGridV2

        call WriteMERCATORGeometry

        call Open_HDF5_OutPut_File

        call OpenAndReadMERCATORFieldsV4

    end subroutine OpenAndReadMERCATORFileV4
    
    !------------------------------------------------------------------------
!Mj***********************************************************************************

    subroutine OpenAndReadMERCATORFileV6

        !Local-----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        call OpenAndReadBathymMERCATORV6

        call ConstructGridV2

        call WriteMERCATORGeometry

        call Open_HDF5_OutPut_File

        call OpenAndReadMERCATORFieldsV6

    end subroutine OpenAndReadMERCATORFileV6
    


!***************************************************************************************
    subroutine OpenAndReadMERCATORFileV5
        
        !PSY2V4(version5) and PSYV3(version3) have the same format for grid/bathymetry 
        call OpenAndReadBathymMERCATORV3  

        call ConstructGridV2

        call WriteMERCATORGeometry

        call Open_HDF5_OutPut_File

        call OpenAndReadMERCATORFieldsV5
    
    end subroutine OpenAndReadMERCATORFileV5
    
    
    !------------------------------------------------------------------------

    subroutine OpenAndReadBathymMERCATOR

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), allocatable     :: Aux3D
        real, dimension(:,:  ), allocatable     :: Aux2D
        real, dimension(:    ), allocatable     :: Depth, LayersInterface
        logical                                 :: exist, FirstTime
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: i, j, k, n
        integer                                 :: ncid, status, dimid


        !Begin----------------------------------------------------------------

        !Verifies if file exists
        inquire(file = Me%InputGridFile, exist = exist)
        
i1:     if (exist) then

            status=NF90_OPEN(trim(Me%InputGridFile),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR10'

            status=NF90_INQ_DIMID(ncid,"y",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR20'

            status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Me%imax)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR30'


            status=NF90_INQ_DIMID(ncid,"x",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR40'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%jmax)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR50'

            status=NF90_INQ_DIMID(ncid,"depth",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR60'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%kmax)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR70'

            allocate(Depth(1:Me%kmax))

            status = nf90_inq_varid(ncid, 'depth', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR100'

            status = NF90_GET_VAR(ncid,n,Depth)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR110'

            
            !The border cells are not considered because with the avialable data
            !is not possible to compute the horizontal position of the cells corners
            Me%WorkSize%ILB = 1
            Me%WorkSize%IUB = Me%imax - 2

            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = Me%jmax - 2

            Me%WorkSize%KLB = 1
            Me%WorkSize%KUB = Me%kmax


            Me%Size%ILB     = Me%WorkSize%ILB - 1
            Me%Size%IUB     = Me%WorkSize%IUB + 1
            Me%Size%JLB     = Me%WorkSize%JLB - 1
            Me%Size%JUB     = Me%WorkSize%JUB + 1
            Me%Size%KLB     = Me%WorkSize%KLB - 1
            Me%Size%KUB     = Me%WorkSize%KUB + 1

            WILB            = Me%WorkSize%ILB 
            WIUB            = Me%WorkSize%IUB 
            WJLB            = Me%WorkSize%JLB 
            WJUB            = Me%WorkSize%JUB 
            WKLB            = Me%WorkSize%KLB 
            WKUB            = Me%WorkSize%KUB 


            ILB             = Me%Size%ILB
            IUB             = Me%Size%IUB
            JLB             = Me%Size%JLB
            JUB             = Me%Size%JUB
            KLB             = Me%Size%KLB
            KUB             = Me%Size%KUB

            allocate(Me%Bathymetry     (ILB:IUB, JLB:JUB))
            allocate(Me%BathymetryMax  (ILB:IUB, JLB:JUB))
            allocate(Me%WaterPoints3D  (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%ComputeFaces3DU(ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%ComputeFaces3DV(ILB:IUB, JLB:JUB, KLB:KUB))

            allocate(Me%SZZ            (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%XX_IE          (ILB:IUB, JLB:JUB))
            allocate(Me%YY_IE          (ILB:IUB, JLB:JUB))
            allocate(LayersInterface   (KLB:KUB))
            allocate(Me%LayersThickness(KLB:KUB))

            allocate(Aux3D(Me%jmax, Me%imax, Me%kmax))
            allocate(Aux2D(Me%jmax, Me%imax      ))

            status = nf90_inq_varid(ncid, 'mask', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR60'

            status = NF90_GET_VAR(ncid,n,Aux3D)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR70'


            do k=WKLB,WKUB
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%WaterPoints3D(i,j,WKUB+1-k) = Aux3D(j + 1, i + 1, k)
            enddo
            enddo
            enddo

            Me%SZZ(:,:,:) = FillValueReal

            LayersInterface(WKUB) = 0.
            
            do k=WKUB, WKLB, -1
                !Bathymetry with one centimetry precision
                LayersInterface(k-1) = LayersInterface(k) + (Depth(WKUB + 1 - k) - LayersInterface(k)) * 2
                LayersInterface(k-1) = real(int(LayersInterface(k-1)*100))/100.
            enddo

            status = nf90_inq_varid(ncid, 'bathy', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR80'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR90'

            do i=WILB,WIUB
            do j=WJLB,WJUB
                
                if (Me%WaterPoints3D(i,j,WKUB) == WaterPoint) then
                    Me%Bathymetry(i,j) = Aux2D(j + 1, i + 1)
                else
                    Me%Bathymetry(i,j) = -99.
                endif
            enddo
            enddo


            do i=WILB,WIUB
            do j=WJLB,WJUB

                if (Me%WaterPoints3D(i,j,WKUB) == WaterPoint) then
                    Me%SZZ(i, j,WKUB) = 0.
                    FirstTime = .true.
                    do k=WKUB, WKLB, -1
                        if (Me%WaterPoints3D(i,j,k) == WaterPoint) then
                            Me%SZZ(i, j, k-1) = LayersInterface(k-1)
                        else
                            if (FirstTime) then
                                Me%Bathymetry(i,j) = LayersInterface(k+1)
                                FirstTime          = .false.
                            endif
                            Me%SZZ(i, j, k-1) = FillValueReal
                        endif
                    enddo
                else
                    Me%SZZ(i, j, : ) = FillValueReal
                endif
            enddo
            enddo

            Me%LayerNumber = Me%kmax
            do k=WKLB, WKUB
                Me%LayersThickness(k) = LayersInterface(k-1) - LayersInterface(k)
            enddo



            Me%BathymetryMax(:,:) = LayersInterface(0)

            status = nf90_inq_varid(ncid, 'longitude', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR100'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR110'

            do i=WILB,WIUB+1
            do j=WJLB,WJUB+1
               Me%XX_IE(i,j) = (Aux2D(j + 1, i    ) + Aux2D(j + 1, i + 1) +              &
                                Aux2D(j    , i + 1) + Aux2D(j    , i   )) / 4.
            enddo
            enddo
  
            status = nf90_inq_varid(ncid, 'latitude', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR120'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR130'

            do i=WILB,WIUB+1
            do j=WJLB,WJUB+1
               Me%YY_IE(i,j) = (Aux2D(j + 1, i    ) + Aux2D(j + 1, i + 1) +              &
                                Aux2D(j    , i + 1) + Aux2D(j    , i   )) / 4.
            enddo
            enddo

            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR140'

            deallocate(Depth)
            deallocate(LayersInterface)


        else i1

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFile)
            stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR150'

        endif i1

        !Verifies if file exists
        inquire(file = Me%InputGridFileU, exist = exist)
        
i2:     if (exist) then

            status=NF90_OPEN(trim(Me%InputGridFileU),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR160'

            status = nf90_inq_varid(ncid, 'mask', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR170'

            status = NF90_GET_VAR(ncid,n,Aux3D)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR180'


            do k=WKLB,WKUB
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%ComputeFaces3DU(i,j,WKUB+1-k) = Aux3D(j + 1, i + 1, k)
            enddo
            enddo
            enddo

            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR190'

        else i2

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFileU)
            stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR210'
        
        end if i2

        !Verifies if file exists
        inquire(file = Me%InputGridFileV, exist = exist)
        
i3:     if (exist) then

            status=NF90_OPEN(trim(Me%InputGridFileV),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR220'

            status = nf90_inq_varid(ncid, 'mask', n)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR230'

            status = NF90_GET_VAR(ncid,n,Aux3D)
            if (status /= nf90_noerr) stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR240'


            do k=WKLB,WKUB
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%ComputeFaces3DV(i,j,WKUB+1-k) = Aux3D(j + 1, i + 1, k)
            enddo
            enddo
            enddo

            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR250'

        else i3

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFileV)
            stop 'OpenAndReadBathymMERCATOR - ModuleMERCATORFormat - ERR270'
        
        end if i3


        deallocate(Aux3D)
        deallocate(Aux2D)

    end subroutine OpenAndReadBathymMERCATOR
    
    !------------------------------------------------------------------------ 

    subroutine OpenAndReadBathymMERCATORV3

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), allocatable     :: Aux3D
        real, dimension(:,:  ), allocatable     :: Aux2D
        real, dimension(:    ), allocatable     :: Depth, LayersInterface
        logical                                 :: exist, FirstTime, BlockFound
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: i, j, k, n
        integer                                 :: ncid, status, dimid
        character(len=PathLength)               :: InputDepthFile
        integer                                 :: iflag, FirstLine, STAT_CALL


        !Begin----------------------------------------------------------------

        !Open Input Grid File:
        !Verifies if file exists
        inquire(file = Me%InputGridFile, exist = exist)
        
i1:     if (exist) then

            status=NF90_OPEN(trim(Me%InputGridFile),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR10'

            status=NF90_INQ_DIMID(ncid,"y",dimid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR20'

            status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Me%imax)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR30'


            status=NF90_INQ_DIMID(ncid,"x",dimid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR40'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%jmax)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR50'

            status=NF90_INQ_DIMID(ncid,"z",dimid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR60'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%kmax)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR70'

            allocate(Depth(1:Me%kmax))
           
            !The border cells are not considered because with the avialable data
            !is not possible to compute the horizontal position of the cells corners
            Me%WorkSize%ILB = 1
            Me%WorkSize%IUB = Me%imax - 2

            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = Me%jmax - 2

            Me%WorkSize%KLB = 1
            Me%WorkSize%KUB = Me%kmax


            Me%Size%ILB     = Me%WorkSize%ILB - 1
            Me%Size%IUB     = Me%WorkSize%IUB + 1
            Me%Size%JLB     = Me%WorkSize%JLB - 1
            Me%Size%JUB     = Me%WorkSize%JUB + 1
            Me%Size%KLB     = Me%WorkSize%KLB - 1
            Me%Size%KUB     = Me%WorkSize%KUB + 1

            WILB            = Me%WorkSize%ILB 
            WIUB            = Me%WorkSize%IUB 
            WJLB            = Me%WorkSize%JLB 
            WJUB            = Me%WorkSize%JUB 
            WKLB            = Me%WorkSize%KLB 
            WKUB            = Me%WorkSize%KUB 


            ILB             = Me%Size%ILB
            IUB             = Me%Size%IUB
            JLB             = Me%Size%JLB
            JUB             = Me%Size%JUB
            KLB             = Me%Size%KLB
            KUB             = Me%Size%KUB

            allocate(Me%Bathymetry     (ILB:IUB, JLB:JUB))
            allocate(Me%BathymetryMax  (ILB:IUB, JLB:JUB))
            allocate(Me%WaterPoints3D  (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%ComputeFaces3DU(ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%ComputeFaces3DV(ILB:IUB, JLB:JUB, KLB:KUB))

            allocate(Me%SZZ            (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(Me%XX_IE          (ILB:IUB, JLB:JUB))
            allocate(Me%YY_IE          (ILB:IUB, JLB:JUB))
            allocate(LayersInterface   (KLB:KUB))
            allocate(Me%LayersThickness(KLB:KUB))

            allocate(Aux3D(Me%jmax, Me%imax, Me%kmax))
            allocate(Aux2D(Me%jmax, Me%imax      ))

            status = nf90_inq_varid(ncid, 'tmask', n)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR80'

            status = NF90_GET_VAR(ncid,n,Aux3D)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR90'


            do k=WKLB,WKUB
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%WaterPoints3D(i,j,WKUB+1-k) = Aux3D(j + 1, i + 1, k)
            enddo
            enddo
            enddo

            status = nf90_inq_varid(ncid, 'umask', n)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR100'

            status = NF90_GET_VAR(ncid,n,Aux3D)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR110'


            do k=WKLB,WKUB
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%ComputeFaces3DU(i,j,WKUB+1-k) = Aux3D(j + 1, i + 1, k)
            enddo
            enddo
            enddo

            status = nf90_inq_varid(ncid, 'vmask', n)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR120'

            status = NF90_GET_VAR(ncid,n,Aux3D)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR130'


            do k=WKLB,WKUB
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%ComputeFaces3DV(i,j,WKUB+1-k) = Aux3D(j + 1, i + 1, k)
            enddo
            enddo
            enddo

            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR140'

        else i1

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFile)
            stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR150'
        
        end if i1

        !Open Input Vertical Levels File:
        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

BF:         if (BlockFound) then

                call GetData(InputDepthFile, EnterDataID = Me%ObjEnterData, flag = iflag, &
                             Buffer_Line = FirstLine + 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR160'

                inquire(file = InputDepthFile, exist = exist)

i2:             if (exist) then
                                                    
                    status=NF90_OPEN(trim(InputDepthFile),NF90_NOWRITE,ncid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR170'

                    status = nf90_inq_varid(ncid, 'deptht', n)
                    if (status /= nf90_noerr) then
                        
                        status = nf90_inq_varid(ncid, 'depthu', n)

                        if (status /= nf90_noerr) then

                            status = nf90_inq_varid(ncid, 'depthv', n)

                            if (status /= nf90_noerr)                                   &
                            stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR180'

                        endif

                    endif

                    status = NF90_GET_VAR(ncid,n,Depth)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR190'

                    status=NF90_CLOSE(ncid)
                    if (status /= nf90_noerr)                                                   &
                        stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR200'

                else i2

                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR210'

                endif i2

            else BF

                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR220'

            end if BF

            call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR230'

        else   IS

            stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR240'

        end if IS

        Me%SZZ(:,:,:) = FillValueReal

        LayersInterface(WKUB) = 0.
        
        do k=WKUB, WKLB, -1
            !Bathymetry with one centimetry precision
            LayersInterface(k-1) = LayersInterface(k) +                                 &
                                   (Depth(WKUB + 1 - k) - LayersInterface(k)) * 2
            LayersInterface(k-1) = real(int(LayersInterface(k-1)*100))/100.
        enddo

        !Open Input Bathymetry File:
        !Verifies if file exists
        inquire(file = Me%InputBathymetryFile, exist = exist)
        
i3:     if (exist) then

            status=NF90_OPEN(trim(Me%InputBathymetryFile),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR250'

            status = nf90_inq_varid(ncid, 'Bathymetry', n)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR260'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR270'

            do i=WILB,WIUB
            do j=WJLB,WJUB
                
                if (Me%WaterPoints3D(i,j,WKUB) == WaterPoint) then
                    Me%Bathymetry(i,j) = Aux2D(j + 1, i + 1)
                else
                    Me%Bathymetry(i,j) = -99.
                endif
            enddo
            enddo


            do i=WILB,WIUB
            do j=WJLB,WJUB

                if (Me%WaterPoints3D(i,j,WKUB) == WaterPoint) then
                    Me%SZZ(i, j,WKUB) = 0.
                    FirstTime = .true.
                    do k=WKUB, WKLB, -1
                        if (Me%WaterPoints3D(i,j,k) == WaterPoint) then
                            !Me%SZZ(i, j, k) = LayersInterface(k)
                            Me%SZZ(i, j, k-1) = LayersInterface(k-1)
                            !if (k==WKLB) then
                            !    Me%SZZ(i, j, k-1) = Me%Bathymetry(i,j)
                            !endif
                        else
                            if (FirstTime) then
                                Me%Bathymetry(i,j) = LayersInterface(k+1)
                                FirstTime          = .false.
                            endif
                            !Me%SZZ(i, j, k) = FillValueReal
                            Me%SZZ(i, j, k-1) = FillValueReal
                            !if (k==WKLB) then
                            !    Me%SZZ(i, j, k-1) = Me%Bathymetry(i,j)
                            !endif
                        endif
                    enddo
                else
                    Me%SZZ(i, j, : ) = FillValueReal
                endif
            enddo
            enddo

            Me%LayerNumber = Me%kmax
            do k=WKLB, WKUB
                Me%LayersThickness(k) = LayersInterface(k-1) - LayersInterface(k)
            enddo

            Me%BathymetryMax(:,:) = LayersInterface(0)

            status = nf90_inq_varid(ncid, 'nav_lon', n)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR280'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR290'

            do i=WILB,WIUB+1
            do j=WJLB,WJUB+1
               Me%XX_IE(i,j) = (Aux2D(j + 1, i    ) + Aux2D(j + 1, i + 1) +             &
                                Aux2D(j    , i + 1) + Aux2D(j    , i   )) / 4.
            enddo
            enddo
  
            status = nf90_inq_varid(ncid, 'nav_lat', n)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR300'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR310'

            do i=WILB,WIUB+1
            do j=WJLB,WJUB+1
               Me%YY_IE(i,j) = (Aux2D(j + 1, i    ) + Aux2D(j + 1, i + 1) +             &
                                Aux2D(j    , i + 1) + Aux2D(j    , i   )) / 4.
            enddo
            enddo

            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)                                                   &
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR320'

            deallocate(Depth)
            deallocate(LayersInterface)


        else i3

            write (*,*) "The input bathymetry file do not exist : ",trim(Me%InputBathymetryFile)
            stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR330'

        endif i3

        if (Me%ComputeBarotropicVel) then
            
            !Open Input Barotropic Scaling Factor File:
            !Verifies if file exists
            inquire(file = Me%InputMeshZGridFile, exist = exist)
        
i4:         if (exist) then

                status=NF90_OPEN(trim(Me%InputMeshZGridFile),NF90_NOWRITE,ncid)
                if (status /= nf90_noerr)                                               &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR340'

                status = nf90_inq_varid(ncid, 'e3u_ps', n)
                if (status /= nf90_noerr)                                               &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR350'

                status = NF90_GET_VAR(ncid,n,Aux3D)
                if (status /= nf90_noerr)                                               &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR360'

                allocate(Me%ScaFactorU (ILB:IUB, JLB:JUB, KLB:KUB))

                do k=WKLB,WKUB
                do i=WILB,WIUB
                do j=WJLB,WJUB
                    Me%ScaFactorU(i,j,WKUB+1-k) = Aux3D(j+1,i+1,WKUB+1-k)
                enddo
                enddo
                enddo

                status = nf90_inq_varid(ncid, 'e3v_ps', n)
                if (status /= nf90_noerr)                                               &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR370'

                status = NF90_GET_VAR(ncid,n,Aux3D)
                if (status /= nf90_noerr)                                               &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR380'

                allocate(Me%ScaFactorV (ILB:IUB, JLB:JUB, KLB:KUB))

                do k=WKLB,WKUB
                do i=WILB,WIUB
                do j=WJLB,WJUB
                    Me%ScaFactorV(i,j,WKUB+1-k) = Aux3D(j+1,i+1,WKUB+1-k)
                enddo
                enddo
                enddo

                status=NF90_CLOSE(ncid)
                if (status /= nf90_noerr)                                               &
                    stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR390'

            else i4

                write (*,*) "The input scaling file do not exist : ",trim(Me%InputMeshZGridFile)
                stop 'OpenAndReadBathymMERCATORV3 - ModuleMERCATORFormat - ERR400'

            endif i4

        endif

        deallocate(Aux3D)
        deallocate(Aux2D)

    end subroutine OpenAndReadBathymMERCATORV3
    
    !------------------------------------------------------------------------

    subroutine OpenAndReadBathymMERCATORV4

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), allocatable     :: Aux3D
        real, dimension(:,:  ), allocatable     :: Aux2D
        real, dimension(:    ), allocatable     :: AuxLong, AuxLat 
        real, dimension(:    ), allocatable     :: Depth, LayersInterface
        logical                                 :: exist, FirstTime
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: i, j, k, n !, d
        integer                                 :: ncid, status, dimid
        character(len=PathLength)               :: InputGridFile
        logical                                 :: BlockFound
        integer                                 :: iflag, FirstLine, STAT_CALL

        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: nDimensions
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName

        real, parameter :: MERCATORFillVReal    = 1.0e35

        !Begin----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

BF:         if (BlockFound) then

                call GetData(InputGridFile, EnterDataID = Me%ObjEnterData, flag = iflag, &
                             Buffer_Line = FirstLine + 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR10'

                inquire(file = InputGridFile, exist = exist)

i1:             if (exist) then

!                    status=NF90_OPEN(trim(InputGridFile),NF90_NOWRITE,ncid)
                    status=NF90_OPEN("teste.nc",NF90_NOWRITE, ncid)
                    if (status /= nf90_noerr)                                           &
                       !stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR20'

                    status=NF90_INQ_DIMID(ncid,"latitude",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR30'

                    status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Me%imax)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR40'


                    status=NF90_INQ_DIMID(ncid,"longitude",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4- ModuleMERCATORFormat - ERR50'

                    status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%jmax)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR60'

                    status=NF90_INQ_DIMID(ncid,"depth",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR70'

                    status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%kmax)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR80'

                    allocate(Depth(1:Me%kmax))

                    status = nf90_inq_varid(ncid, 'depth', n)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR90'

                    status = NF90_GET_VAR(ncid,n,Depth)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR100'
           
                    !The border cells are not considered because with the avialable data
                    !is not possible to compute the horizontal position of the cells corners
                    Me%WorkSize%ILB = 1
                    Me%WorkSize%IUB = Me%imax - 2

                    Me%WorkSize%JLB = 1
                    Me%WorkSize%JUB = Me%jmax - 2

                    Me%WorkSize%KLB = 1
                    Me%WorkSize%KUB = Me%kmax


                    Me%Size%ILB     = Me%WorkSize%ILB - 1
                    Me%Size%IUB     = Me%WorkSize%IUB + 1
                    Me%Size%JLB     = Me%WorkSize%JLB - 1
                    Me%Size%JUB     = Me%WorkSize%JUB + 1
                    Me%Size%KLB     = Me%WorkSize%KLB - 1
                    Me%Size%KUB     = Me%WorkSize%KUB + 1

                    WILB            = Me%WorkSize%ILB 
                    WIUB            = Me%WorkSize%IUB 
                    WJLB            = Me%WorkSize%JLB 
                    WJUB            = Me%WorkSize%JUB 
                    WKLB            = Me%WorkSize%KLB 
                    WKUB            = Me%WorkSize%KUB 


                    ILB             = Me%Size%ILB
                    IUB             = Me%Size%IUB
                    JLB             = Me%Size%JLB
                    JUB             = Me%Size%JUB
                    KLB             = Me%Size%KLB
                    KUB             = Me%Size%KUB

                    allocate(Me%Bathymetry     (ILB:IUB, JLB:JUB))
                    allocate(Me%BathymetryMax  (ILB:IUB, JLB:JUB))
                    allocate(Me%WaterPoints3D  (ILB:IUB, JLB:JUB, KLB:KUB))

                    allocate(Me%SZZ            (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%XX_IE          (ILB:IUB, JLB:JUB))
                    allocate(Me%YY_IE          (ILB:IUB, JLB:JUB))
                    allocate(LayersInterface   (KLB:KUB))
                    allocate(Me%LayersThickness(KLB:KUB))

                    allocate(Aux3D(Me%jmax, Me%imax, Me%kmax))
                    allocate(Aux2D(Me%jmax, Me%imax      ))
                    allocate(AuxLong(Me%jmax))
                    allocate(AuxLat(Me%imax))

                    status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR110'

d0:                 do n=1,nVars

                        status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
                        if (status /= nf90_noerr)                                       &
                            stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR120'

i2:                     if (CheckName(nameAux, MohidName)) then
                
                            if      (nDimensions == 3) then
                        
                                if (MohidName == GetPropertyName(Temperature_)) then 

                                    status = NF90_GET_VAR(ncid,n,Aux3D)
                                    if (status /= nf90_noerr)                           &
                                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR130'

                                    !The boundary cells are not read
                                    do k=WKLB,WKUB
                                    do i=WILB,WIUB
                                    do j=WJLB,WJUB

                                        if (Aux3D(j+1,i+1,WKUB+1-k) .lt.                &
                                            MERCATORFillVReal) then

                                            Me%WaterPoints3D(i,j,k) = 1

                                        else

                                            Me%WaterPoints3D(i,j,k) = 0

                                        endif

                                    enddo
                                    enddo
                                    enddo

                                    exit d0

                                endif

                            endif

                        endif i2 

                    enddo d0

                    Me%SZZ(:,:,:) = FillValueReal

                    LayersInterface(WKUB) = 0.
            
                    do k=WKUB, WKLB, -1
                        !Bathymetry with one centimetry precision
                        LayersInterface(k-1) = LayersInterface(k) +                     &
                                               (Depth(WKUB + 1 - k) -                   &
                                               LayersInterface(k)) * 2
                        LayersInterface(k-1) = real(int(LayersInterface(k-1)*100))/100.
                    enddo

                    do i=WILB,WIUB
                    do j=WJLB,WJUB

                        if (Me%WaterPoints3D(i,j,WKUB) == WaterPoint) then
                            Me%SZZ(i, j,WKUB) = 0.
                            FirstTime = .true.
                            do k=WKUB, WKLB, -1
                                if (Me%WaterPoints3D(i,j,k) == WaterPoint) then
                                    Me%SZZ(i, j, k-1) = LayersInterface(k-1)
                                else
                                    if (FirstTime) then
                                        Me%Bathymetry(i,j) = LayersInterface(k)
                                        FirstTime          = .false.
                                    endif
                                    Me%SZZ(i, j, k-1) = FillValueReal
                                endif
                            enddo
                        else
                            Me%SZZ(i, j, : ) = FillValueReal
                            
                            Me%Bathymetry(i,j) = -99.
                        endif
                    enddo
                    enddo

                    Me%LayerNumber = Me%kmax
                    do k=WKLB, WKUB
                        Me%LayersThickness(k) = LayersInterface(k-1) - LayersInterface(k)
                    enddo

                    Me%BathymetryMax(:,:) = LayersInterface(0)

                    status = nf90_inq_varid(ncid, 'longitude', n)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR140'

                    status = NF90_GET_VAR(ncid,n,AuxLong)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR150'

                    do i=WILB,WIUB+1
                    do j=WJLB,WJUB+1
                       Me%XX_IE(i,j) = (AuxLong(j + 1) + AuxLong(j)) / 2.
                    enddo
                    enddo
  
                    status = nf90_inq_varid(ncid, 'latitude', n)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR160'

                    status = NF90_GET_VAR(ncid,n,AuxLat)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR170'

                    do i=WILB,WIUB+1
                    do j=WJLB,WJUB+1
                       Me%YY_IE(i,j) = (AuxLat(i) + AuxLat(i + 1)) / 2.
                    enddo
                    enddo

                    status=NF90_CLOSE(ncid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR180'

                    deallocate(Depth)
                    deallocate(LayersInterface)

                else i1

                    write (*,*) "The input grid file do not exist : ",                  &
                                trim(Me%InputGridFile)
                    stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR190'

                endif i1

            else BF

                stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR200'

            end if BF

            call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR210'

        else   IS

            stop 'OpenAndReadBathymMERCATORV4 - ModuleMERCATORFormat - ERR220'

        end if IS

        deallocate(Aux3D)
        deallocate(Aux2D)

        deallocate(AuxLong)
        deallocate(AuxLat)

    end subroutine OpenAndReadBathymMERCATORV4
    
    !--------------------------------------------------------------------------

    subroutine OpenAndReadBathymMERCATORV6

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer, dimension(:,:,:,:), allocatable     :: Aux4D
        real, dimension(:,:  ), allocatable     :: Aux2D
        real, dimension(:    ), allocatable     :: AuxLong, AuxLat 
        real, dimension(:    ), allocatable     :: Depth, LayersInterface
        logical                                 :: exist
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: i, j, k, n
        integer                                 :: ncid, status, dimid
        character(len=PathLength)               :: InputGridFile
        logical                                 :: BlockFound
        integer                                 :: iflag, FirstLine, STAT_CALL

        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: nDimensions
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName
        integer                                 :: ninst, MERCATORFillVInteger
        real, parameter                         :: MERCATORFillVReal    = 1.0e35

        !Begin----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, STAT = STAT_CALL)


IS:     if(STAT_CALL .EQ. SUCCESS_) then

            
BF:         if (BlockFound) then

                call GetData(InputGridFile, EnterDataID = Me%ObjEnterData, flag = iflag, &
                             Buffer_Line = FirstLine + 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR10'

                inquire(file = InputGridFile, exist = exist)

i1:             if (exist) then

                    status=NF90_OPEN(trim(InputGridFile),NF90_NOWRITE,ncid)
                    if (status /= nf90_noerr)                                           &
                    stop 'OpenAndReadBathymMERCATORV6- ModuleMERCATORFormat - ERR20'

                    status=NF90_INQ_DIMID(ncid,"time",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR30'

                    status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Ninst)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR40'
                        
                    status=NF90_INQ_DIMID(ncid,"latitude",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR30'

                    status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Me%imax)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR40'

                    status=NF90_INQ_DIMID(ncid,"longitude",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6- ModuleMERCATORFormat - ERR50'

                    status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%jmax)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR60'

                    status=NF90_INQ_DIMID(ncid,"depth",dimid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR70'

                    status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%kmax)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR80'

                    allocate(Depth(1:Me%kmax))

                    status = nf90_inq_varid(ncid, 'depth', n)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR90'

                    status = NF90_GET_VAR(ncid,n,Depth)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR100'
           
                    !The border cells are not considered because with the avialable data
                    !is not possible to compute the horizontal position of the cells corners
                    Me%WorkSize%ILB = 1
                    Me%WorkSize%IUB = Me%imax - 2

                    Me%WorkSize%JLB = 1
                    Me%WorkSize%JUB = Me%jmax - 2

                    Me%WorkSize%KLB = 1
                    Me%WorkSize%KUB = Me%kmax


                    Me%Size%ILB     = Me%WorkSize%ILB - 1
                    Me%Size%IUB     = Me%WorkSize%IUB + 1
                    Me%Size%JLB     = Me%WorkSize%JLB - 1
                    Me%Size%JUB     = Me%WorkSize%JUB + 1
                    Me%Size%KLB     = Me%WorkSize%KLB - 1
                    Me%Size%KUB     = Me%WorkSize%KUB + 1

                    WILB            = Me%WorkSize%ILB 
                    WIUB            = Me%WorkSize%IUB 
                    WJLB            = Me%WorkSize%JLB 
                    WJUB            = Me%WorkSize%JUB 
                    WKLB            = Me%WorkSize%KLB 
                    WKUB            = Me%WorkSize%KUB 


                    ILB             = Me%Size%ILB
                    IUB             = Me%Size%IUB
                    JLB             = Me%Size%JLB
                    JUB             = Me%Size%JUB
                    KLB             = Me%Size%KLB
                    KUB             = Me%Size%KUB

                    allocate(Me%Bathymetry     (ILB:IUB, JLB:JUB))
                    allocate(Me%BathymetryMax  (ILB:IUB, JLB:JUB))
                    allocate(Me%WaterPoints3D  (ILB:IUB, JLB:JUB, KLB:KUB))

                    allocate(Me%SZZ            (ILB:IUB, JLB:JUB, KLB:KUB))
                    allocate(Me%XX_IE          (ILB:IUB, JLB:JUB))
                    allocate(Me%YY_IE          (ILB:IUB, JLB:JUB))
                    allocate(LayersInterface   (KLB:KUB))
                    allocate(Me%LayersThickness(KLB:KUB))


                    allocate(Aux4D(Me%jmax, Me%imax, Me%kmax,ninst))
                    allocate(Aux2D(Me%jmax, Me%imax      ))
                    allocate(AuxLong(Me%jmax))
                    allocate(AuxLat(Me%imax))
                       
                    status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR110'

d0:                 do n=1,nVars

                        status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
                        if (status /= nf90_noerr)                                       &
                            stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR120'


i2:                     if (CheckName(nameAux, MohidName)) then
                
                            if      (nDimensions == 4) then
                            
                                if (MohidName == GetPropertyName(Temperature_)) then 

                                    status = NF90_GET_VAR(ncid,n,Aux4D)
                                    if (status /= nf90_noerr)                           &
                                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR130'

                                    status=NF90_GET_ATT(ncid,n,"_FillValue",MERCATORFillVInteger)
                                    if (status /= nf90_noerr)                                                   &
                                        stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR81'
                                        
                                    Me%WaterPoints3D(:,:,:)= 0

                                    do k=WKLB,WKUB
                                    do i=WILB,WIUB
                                    do j=WJLB,WJUB
                                        
                                        if (Aux4D(j+1,i+1,k,1)/= MERCATORFillVInteger) then
                                            Me%WaterPoints3D(i,j,WKUB+1-k) = 1
                                        end if
                                    enddo
                                    enddo
                                    enddo

                                    exit d0

                                endif

                            endif

                        endif i2 

                    enddo d0

                    Me%SZZ(:,:,:) = FillValueReal

   !MJ ***********************************************
                    do k=WKUB,WKLB,-1
                        LayersInterface(k)=Depth(WKUB - WKLB - k + 2)
                    enddo
                    
                    LayersInterface(WKLB-1) = LayersInterface(WKLB) + LayersInterface(WKLB) - LayersInterface(WKLB+1)

                    Me%BathymetryMax(:,:) = LayersInterface(WKLB-1)                        


                    do i=WILB,WIUB
                    do j=WJLB,WJUB
                    
                        Me%Bathymetry(i,j)=-99.
                            
                        if(Me%WaterPoints3D(i,j,WKUB) == 1) then
                            Me%Bathymetry(i,j) = LayersInterface(WKLB-1)
                            do k=WKUB-1,WKLB,-1
                                if (Me%WaterPoints3D(i,j,k)==0) then
                                    Me%Bathymetry(i,j) = LayersInterface(k)
                                    exit
                                endif
                            enddo
                        endif
                    enddo
                    enddo


                    Me%SZZ(:,:,:) = FillValueReal


                    do i=WILB,WIUB
                    do j=WJLB,WJUB
                        if (Me%WaterPoints3D(i,j,WKUB) == WaterPoint) then
                            Me%SZZ(i, j, WKUB) = LayersInterface(WKUB)
                            do k=WKUB, WKLB, -1
                                if (Me%WaterPoints3D(i,j,k) == WaterPoint)              &                           
                                    Me%SZZ(i, j, k-1) = LayersInterface(k-1)
                            enddo
                        endif
                    enddo
                    enddo

                    Me%LayerNumber = Me%kmax

                    do k=WKUB, WKLB, -1
                        Me%LayersThickness(k) = LayersInterface(k-1) - LayersInterface(k)
                    enddo

                    status = nf90_inq_varid(ncid, 'longitude', n)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR140'

                    status = NF90_GET_VAR(ncid,n,AuxLong)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR150'

                    do i=WILB,WIUB+1
                    do j=WJLB,WJUB+1
                       Me%XX_IE(i,j) = (AuxLong(j + 1) + AuxLong(j)) / 2.
                    enddo
                    enddo
  
                    status = nf90_inq_varid(ncid, 'latitude', n)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR160'

                    status = NF90_GET_VAR(ncid,n,AuxLat)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR170'

                    do i=WILB,WIUB+1
                    do j=WJLB,WJUB+1
                       Me%YY_IE(i,j) = (AuxLat(i) + AuxLat(i + 1)) / 2.
                    enddo
                    enddo

                    status=NF90_CLOSE(ncid)
                    if (status /= nf90_noerr)                                           &
                        stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR180'

                    deallocate(Depth)
                    deallocate(LayersInterface)
 
                else i1

                    write (*,*) "The input grid file do not exist : ",                  &
                                trim(Me%InputGridFile)
                    stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR190'

                endif i1

            else BF

                stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR200'

            end if BF

            call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR210'

        else   IS

            stop 'OpenAndReadBathymMERCATORV6 - ModuleMERCATORFormat - ERR220'

        end if IS

        deallocate(Aux4D)
        deallocate(Aux2D)

        deallocate(AuxLong)
        deallocate(AuxLat)

    end subroutine OpenAndReadBathymMERCATORV6
    
    !--------------------------------------------------------------------------    
    
    subroutine OpenAndReadMERCATORFields

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: exist, BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL


        !Begin----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

BF:         if (BlockFound) then

                do line = FirstLine + 1, LastLine - 1

                    call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMERCATORFields - ModuleMERCATORFormat - ERR10'

                    inquire(file = InputFile, exist = exist)
       
i1:                 if (exist) then
                                                    
                        call ReadMercatorFile (InputFile) 

                    endif i1
                enddo

            else BF

                stop 'OpenAndReadMERCATORFields - ModuleMERCATORFormat - ERR20'

            end if BF

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadMERCATORFields - ModuleMERCATORFormat - ERR30'

        else   IS

            stop 'OpenAndReadMERCATORFields - ModuleMERCATORFormat - ERR40'

        end if IS



    end subroutine OpenAndReadMERCATORFields


    !------------------------------------------------------------------------

    subroutine OpenAndReadMERCATORFieldsV3

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: exist, BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL


        !Begin----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

            !The block is found to exist before when reading depth

            do line = FirstLine + 1, LastLine - 1

                call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                             Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadMERCATORFieldsV3 - ModuleMERCATORFormat - ERR10'

                inquire(file = InputFile, exist = exist)
   
i1:             if (exist) then
                                                
                    call ReadMercatorFileV3 (InputFile) 

                endif i1
            enddo

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadMERCATORFieldsV3 - ModuleMERCATORFormat - ERR20'

        else   IS

            stop 'OpenAndReadMERCATORFieldsV3 - ModuleMERCATORFormat - ERR30'

        end if IS


    end subroutine OpenAndReadMERCATORFieldsV3

    !------------------------------------------------------------------------

    subroutine OpenAndReadMERCATORFieldsV4

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: exist, BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL


        !Begin----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

            !The block is found to exist before when reading depth

            do line = FirstLine + 1, LastLine - 1

                call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                             Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadMERCATORFieldsV4 - ModuleMERCATORFormat - ERR10'

                inquire(file = InputFile, exist = exist)
   
i1:             if (exist) then
                                                
                    call ReadMercatorFileV4 (InputFile) 

                endif i1
            enddo

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadMERCATORFieldsV4 - ModuleMERCATORFormat - ERR20'

        else   IS

            stop 'OpenAndReadMERCATORFieldsV4 - ModuleMERCATORFormat - ERR30'

        end if IS


    end subroutine OpenAndReadMERCATORFieldsV4

    !------------------------------------------------------------------------

!MJ***********************
    subroutine OpenAndReadMERCATORFieldsV6

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: exist, BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL


        !Begin----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

            !The block is found to exist before when reading depth

            do line = FirstLine + 1, LastLine - 1

                call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                             Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadMERCATORFieldsV6 - ModuleMERCATORFormat - ERR10'

                inquire(file = InputFile, exist = exist)
   
i1:             if (exist) then
                                                
                    call ReadMercatorFileV6 (InputFile) 

                endif i1
            enddo

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadMERCATORFieldsV6 - ModuleMERCATORFormat - ERR20'

        else   IS

            stop 'OpenAndReadMERCATORFieldsV6 - ModuleMERCATORFormat - ERR30'

        end if IS


    end subroutine OpenAndReadMERCATORFieldsV6

    !------------------------------------------------------------------------




    
    subroutine OpenAndReadMERCATORFieldsV5

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: exist, BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL


        !Begin----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

            !The block is found to exist before when reading depth

            do line = FirstLine + 1, LastLine - 1

                call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                             Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OpenAndReadMERCATORFieldsV5 - ModuleMERCATORFormat - ERR10'

                inquire(file = InputFile, exist = exist)
   
i1:             if (exist) then
                                                
                    call ReadMercatorFileV5 (InputFile) 

                endif i1
            enddo

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadMERCATORFieldsV5 - ModuleMERCATORFormat - ERR20'

        else   IS

            stop 'OpenAndReadMERCATORFieldsV5 - ModuleMERCATORFormat - ERR30'

        end if IS


    end subroutine OpenAndReadMERCATORFieldsV5

    !------------------------------------------------------------------------


    subroutine ReadMercatorFile(InputFile)

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:,:), pointer       :: Aux4D
        real, dimension(:,:,:  ), pointer       :: Aux3D, Aux3DModulus
        real, dimension(:,:    ), pointer       :: Aux2D
        real(8), dimension(:    ), allocatable     :: AuxDays
        integer                                 :: ni, n, nInst, iOut, i, j, k
        integer                                 :: ncid, status, dimid
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: WestON, EastON, SouthON, NorthON
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName
        type (T_Time)                           :: FieldTime


        !Begin----------------------------------------------------------------

        !Verifies if file exists

        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR10'

        status=NF90_INQ_DIMID(ncid,"time",dimid)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR20'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = nInst)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR30'

        allocate(AuxDays(1:nInst))

        status = nf90_inq_varid(ncid, 'time', n)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR40'

        status = NF90_GET_VAR(ncid,n,AuxDays)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR50'

        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR60'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR70'

            if (nDimensions == 4) then
                !Temperature, Salinity, Meridional vel., Zonal vel.
                allocate(Aux4D(Me%jmax, Me%imax, Me%kmax, nInst))
                allocate(Aux3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                Aux3D(:,:,:) = FillValueReal

                allocate(Aux3DModulus(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                Aux3DModulus(:,:,:) = FillValueReal

            elseif (nDimensions == 3) then
                !Sea surface height
                allocate(Aux4D(Me%jmax, Me%imax, 1      , nInst))
                allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                Aux2D(:,:) = FillValueReal
            elseif (nDimensions == 2 .or. nDimensions == 1) then
                !Grid properties already written
                cycle
            endif

            status = NF90_GET_VAR(ncid,n,Aux4D)
            if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR80'

i2:         if (CheckName(nameAux, MohidName)) then
                
d1:             do ni = 1, nInst

                    iOut = OutputInstants(MohidName)

                    FieldTime = Me%RefDateTime + AuxDays(ni) * 86400.

                    if      (nDimensions == 4) then
                        
                        !The boundary cells are not read
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then
                                
                                if      (MohidName == GetPropertyName(VelocityU_)) then

                                    WestON = Me%ComputeFaces3DU(i, j-1, k) 
                                    EastON = Me%ComputeFaces3DU(i, j  , k) 

                                    Aux3D(i, j, k) = (Aux4D(j  ,i+1,Me%WorkSize%KUB+1-k,ni) * WestON  + &
                                                      Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) * EastON) / 2.

                                    if (abs(Aux3D(i, j, k))> 200.) Aux3D(i, j, k) = 0.

                                                    
                                else if (MohidName == GetPropertyName(VelocityV_)) then
                                    
                                    SouthON = Me%ComputeFaces3DV(i-1, j, k) 
                                    NorthON = Me%ComputeFaces3DV(i,   j, k) 
                                    
                                    Aux3D(i, j, k) = (Aux4D(j+1,i  ,Me%WorkSize%KUB+1-k,ni) * SouthON  + &
                                                      Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) * NorthON) / 2.

                                   if (abs(Aux3D(i, j, k))> 200.) Aux3D(i, j, k) = 0.

                                 
                                else if (MohidName == GetPropertyName(Temperature_) .or.  &
                                         MohidName == GetPropertyName(Salinity_   )) then

                                    Aux3D(i, j, k) = Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni)
                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        
                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux3D = Aux3D)


                    else if (nDimensions == 3) then

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, Me%WorkSize%KUB) == WaterPoint) then
                                
                                if      (MohidName == GetPropertyName(BarotropicVelocityU_)) then

                                    WestON = Me%ComputeFaces3DU(i, j-1, Me%WorkSize%KUB) 
                                    EastON = Me%ComputeFaces3DU(i, j  , Me%WorkSize%KUB) 

                                    Aux2D(i, j) = (Aux4D(j  ,i+1,1,ni) * WestON  + &
                                                   Aux4D(j+1,i+1,1,ni) * EastON) / 2.

                                    if (abs(Aux2D(i, j))> 200.) Aux2D(i, j) = 0.

                                else if (MohidName == GetPropertyName(BarotropicVelocityV_)) then
                                    
                                    SouthON = Me%ComputeFaces3DV(i-1, j, Me%WorkSize%KUB) 
                                    NorthON = Me%ComputeFaces3DV(i,   j, Me%WorkSize%KUB) 
                                    
                                
                                    Aux2D(i, j) = (Aux4D(j+1,  i  ,1,ni) * SouthON  + &
                                                   Aux4D(j+1,i+1  ,1,ni) * NorthON) / 2.

                                    if (abs(Aux2D(i, j))> 200.) Aux2D(i, j) = 0.
                                

                                else if (MohidName == GetPropertyName(WaterLevel_ )) then

                                    Aux2D(i, j) = Aux4D(j+1,i+1,1, ni)
                                endif
                            endif

                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D = Aux2D)

                    endif


                enddo d1

            else i2

                call SetError (WARNING_, INTERNAL_, &
                    'This property name was not defined yet - ReadMercatorFile - ModuleMERCATORFormat')

            endif i2

            deallocate(Aux4D)

            if     (nDimensions == 4) then
                deallocate(Aux3D       )
                deallocate(Aux3DModulus)
            elseif (nDimensions == 3) then
                deallocate(Aux2D)
            endif

        enddo d0

        deallocate(AuxDays)


    end subroutine ReadMercatorFile
    

    !------------------------------------------------------------------------

    subroutine ReadMercatorFileV3(InputFile)

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:,:), pointer       :: Aux4D, BaroAux4D
        real, dimension(:,:,:  ), pointer       :: Aux3D, Aux3DModulus
        real, dimension(:,:,:  ), pointer       :: BaroAux3D, SumScaFactor
        real, dimension(:,:    ), pointer       :: Aux2D, BaroAux2D
        real(8), dimension(:    ), allocatable     :: AuxDays
        integer                                 :: ni, n, nInst, iOut, i, j, k
        integer                                 :: ncid, status, dimid
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: WestON, EastON, SouthON, NorthON
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName
        type (T_Time)                           :: FieldTime


        !Begin----------------------------------------------------------------

        !Verifies if file exists

        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR10'

        status=NF90_INQ_DIMID(ncid,"time_counter",dimid)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR20'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = nInst)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR30'

        allocate(AuxDays(1:nInst))

        status = nf90_inq_varid(ncid, 'time_counter', n)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR40'

        status = NF90_GET_VAR(ncid,n,AuxDays)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR50'

        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR60'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr)                                                   &
                stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR70'

            if (nDimensions == 4) then

!                if (CheckNameV3(nameAux, MohidName)) then
                !Temperature, Salinity, Meridional vel., Zonal vel.
                allocate(Aux4D(Me%jmax, Me%imax, Me%kmax, nInst))
                allocate(Aux3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,        &
                        Me%Size%KLB:Me%Size%KUB))
                Aux3D(:,:,:) = FillValueReal

                allocate(Aux3DModulus(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, &
                        Me%Size%KLB:Me%Size%KUB))
                Aux3DModulus(:,:,:) = FillValueReal
!                endif

            elseif (nDimensions == 3) then

!                if (CheckNameV3(nameAux, MohidName)) then
                    !Sea surface height
                    allocate(Aux4D(Me%jmax, Me%imax, 1      , nInst))
                    allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    Aux2D(:,:) = FillValueReal
!                endif

            elseif (nDimensions == 2 .or. nDimensions == 1) then
                !Grid properties already written
                cycle
            endif

            status = NF90_GET_VAR(ncid,n,Aux4D)
            if (status /= nf90_noerr)                                                   &
                stop 'ReadMercatorFileV3 - ModuleMERCATORFormat - ERR80'

i2:         if (CheckNameV3(nameAux, MohidName)) then
                
d1:             do ni = 1, nInst

                    iOut = OutputInstants(MohidName)

                    FieldTime = Me%RefDateTime + (AuxDays(ni) - 43200.)

                    if      (nDimensions == 4) then
                        
                        !The boundary cells are not read
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then
                                
                                if      (MohidName == GetPropertyName(VelocityU_)) then

                                    WestON = Me%ComputeFaces3DU(i, j-1, k) 
                                    EastON = Me%ComputeFaces3DU(i, j  , k) 

                                    Aux3D(i, j, k) = (Aux4D(j  ,i+1,                    &
                                                      Me%WorkSize%KUB+1-k,ni) * WestON + &
                                                      Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k, &
                                                      ni) * EastON) / 2.

                                    if (abs(Aux3D(i, j, k))> 200.) Aux3D(i, j, k) = 0.

                                                    
                                else if (MohidName == GetPropertyName(VelocityV_)) then
                                    
                                    SouthON = Me%ComputeFaces3DV(i-1, j, k) 
                                    NorthON = Me%ComputeFaces3DV(i,   j, k) 
                                    
                                    Aux3D(i, j, k) = (Aux4D(j+1,i  ,                    &
                                                      Me%WorkSize%KUB+1-k,ni) * SouthON  + &
                                                      Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k, &
                                                      ni) * NorthON) / 2.

                                   if (abs(Aux3D(i, j, k))> 200.) Aux3D(i, j, k) = 0.

                                 
                                else if (MohidName == GetPropertyName(Temperature_) .or. &
                                         MohidName == GetPropertyName(Salinity_   )) then

                                    Aux3D(i, j, k) = Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni)
                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        
                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux3D = Aux3D)


                    else if (nDimensions == 3) then

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, Me%WorkSize%KUB) == WaterPoint) then
                                
                                if (MohidName == GetPropertyName(BarotropicVelocityU_)) then

                                    WestON = Me%ComputeFaces3DU(i, j-1, Me%WorkSize%KUB) 
                                    EastON = Me%ComputeFaces3DU(i, j  , Me%WorkSize%KUB) 

                                    Aux2D(i, j) = (Aux4D(j  ,i+1,1,ni) * WestON  + &
                                                   Aux4D(j+1,i+1,1,ni) * EastON) / 2.

                                    if (abs(Aux2D(i, j))> 200.) Aux2D(i, j) = 0.

                                else if (MohidName == GetPropertyName(BarotropicVelocityV_)) then
                                    
                                    SouthON = Me%ComputeFaces3DV(i-1, j, Me%WorkSize%KUB) 
                                    NorthON = Me%ComputeFaces3DV(i,   j, Me%WorkSize%KUB) 
                                    
                                
                                    Aux2D(i, j) = (Aux4D(j+1,  i  ,1,ni) * SouthON  + &
                                                   Aux4D(j+1,i+1  ,1,ni) * NorthON) / 2.

                                    if (abs(Aux2D(i, j))> 200.) Aux2D(i, j) = 0.
                                

                                else if (MohidName == GetPropertyName(WaterLevel_ )) then

                                    Aux2D(i, j) = Aux4D(j+1,i+1,1, ni)
                                endif
                            endif

                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D = Aux2D)

                    endif


                enddo d1

                if (MohidName == GetPropertyName(VelocityU_) .and.                      &
                    Me%ComputeBarotropicVel) then

                    allocate(BaroAux4D(Me%jmax, Me%imax, Me%kmax , nInst))
                    allocate(BaroAux2D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB))
                    allocate(BaroAux3D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB,                         &
                                       nInst))
                    allocate(SumScaFactor(Me%Size%ILB:Me%Size%IUB,                      &
                                          Me%Size%JLB:Me%Size%JUB,                      &
                                          nInst))
                    BaroAux2D(:,:) = FillValueReal

d2:                 do ni = 1, nInst

                        SumScaFactor(:,:,ni) = 0.

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                            if (k == Me%WorkSize%KLB) then

                                BaroAux3D(i+1,j, ni) = 0.
                                BaroAux3D(i+1,j+1, ni) = 0.

                                SumScaFactor(i+1,j, ni) = 0.
                                SumScaFactor(i+1,j+1, ni) = 0.

                            endif

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then

                                BaroAux2D(i, j) = 0.

                                WestON = Me%ComputeFaces3DU(i, j-1, k) 
                                EastON = Me%ComputeFaces3DU(i, j  , k) 

                                if ((WestON <= 1) .and. (EastON <= 1)) then

                                    BaroAux4D(j  ,i+1, Me%WorkSize%KUB+1-k,ni) =            &
                                        Aux4D(j  ,i+1, Me%WorkSize%KUB+1-k,ni) * WestON*    &
                                        Me%ScaFactorU(i+1,j, Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i+1,j, ni) = SumScaFactor(i+1,j, ni) +     &
                                        Me%ScaFactorU(i+1,j, Me%WorkSize%KUB+1-k)* WestON   

                                    BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) =            &
                                        Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) * EastON*    &
                                        Me%ScaFactorU(i+1,j+1,Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i+1,j+1, ni) = SumScaFactor(i+1,j+1,ni) +  &
                                        Me%ScaFactorU(i+1,j+1, Me%WorkSize%KUB+1-k)* EastON   

                                    BaroAux3D(i+1,j, ni) = BaroAux3D(i+1,j, ni) +           &
                                        BaroAux4D(j,i+1, Me%WorkSize%KUB+1-k, ni)

                                    BaroAux3D(i+1,j+1, ni) = BaroAux3D(i+1,j+1, ni) +       &
                                        BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni)

                                    if ((k == Me%WorkSize%KUB) .and.                        &
                                       (SumScaFactor(i+1,j, ni) > 0. .and.                 &
                                       SumScaFactor(i+1,j+1, ni) > 0.)) then

                                        BaroAux2D(i,j) =                                    &
                                        (BaroAux3D(i+1,j, ni)/SumScaFactor(i+1,j, ni) +     &
                                        BaroAux3D(i+1,j+1, ni)/SumScaFactor(i+1,j+1, ni)) / 2.

                                        if (abs(BaroAux2D(i, j))> 200.) BaroAux2D(i, j) = 0.
                                    endif

                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime,                                  &
                                            GetPropertyName(BarotropicVelocityU_),      &
                                            iOut, Aux2D = BaroAux2D)

                    enddo d2

                    deallocate(BaroAux4D)
                    deallocate(BaroAux3D)
                    deallocate(BaroAux2D)
                    deallocate(SumScaFactor)

                endif

                if (MohidName == GetPropertyName(VelocityV_) .and.                      &
                    Me%ComputeBarotropicVel) then

                    allocate(BaroAux4D(Me%jmax, Me%imax, Me%kmax , nInst))
                    allocate(BaroAux2D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB))
                    allocate(BaroAux3D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB,                         &
                                       nInst))
                    allocate(SumScaFactor(Me%Size%ILB:Me%Size%IUB,                      &
                                          Me%Size%JLB:Me%Size%JUB,                      &
                                          nInst))
                    BaroAux2D(:,:) = FillValueReal

d3:                 do ni = 1, nInst

                        SumScaFactor(:,:,ni) = 0.

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                            if (k == Me%WorkSize%KLB) then

                                BaroAux3D(i,j+1,ni) = 0.
                                BaroAux3D(i+1,j+1, ni) = 0.

                                SumScaFactor(i, j+1, ni) = 0.
                                SumScaFactor(i+1,j+1, ni) = 0.

                            endif

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then

                                BaroAux2D(i, j) = 0.

                                SouthON = Me%ComputeFaces3DV(i-1, j, k) 
                                NorthON = Me%ComputeFaces3DV(i,   j, k) 

                                if ((SouthON <= 1) .and. (NorthON <= 1)) then

                                    BaroAux4D(j+1,i ,Me%WorkSize%KUB+1-k,ni) =              &
                                        Aux4D(j+1,i ,Me%WorkSize%KUB+1-k,ni) * SouthON*     &
                                        Me%ScaFactorV(i, j+1,Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i, j+1, ni) = SumScaFactor(i, j+1, ni) +   &
                                        Me%ScaFactorV(i, j+1, Me%WorkSize%KUB+1-k)* SouthON   

                                    BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) =            &
                                        Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) * NorthON*   &
                                        Me%ScaFactorV(i+1,j+1,Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i+1,j+1, ni) = SumScaFactor(i+1,j+1, ni) + &
                                        Me%ScaFactorV(i+1,j+1, Me%WorkSize%KUB+1-k)* NorthON   

                                    BaroAux3D(i,j+1,ni) = BaroAux3D(i,j+1,ni) +             &
                                        BaroAux4D(j+1,i ,Me%WorkSize%KUB+1-k,ni)

                                    BaroAux3D(i+1,j+1, ni) = BaroAux3D(i+1,j+1, ni) +       &
                                        BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni)

                                    if ((k == Me%WorkSize%KUB) .and.                        &
                                        (SumScaFactor(i, j+1, ni) > 0. .and.                 &
                                        SumScaFactor(i+1,j+1, ni) > 0.)) then

                                        BaroAux2D(i,j) =                                    &
                                        (BaroAux3D(i,j+1,ni)/SumScaFactor(i, j+1, ni) +     &
                                        BaroAux3D(i+1,j+1, ni)/SumScaFactor(i+1,j+1, ni)) / 2.
                                
                                        if (abs(BaroAux2D(i, j))> 200.) BaroAux2D(i, j) = 0.
                                    endif

                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, GetPropertyName(BarotropicVelocityV_), &
                                            iOut, Aux2D = BaroAux2D)

                    enddo d3

                    deallocate(BaroAux4D)
                    deallocate(BaroAux3D)
                    deallocate(BaroAux2D)
                    deallocate(SumScaFactor)

                endif

!            else i2

!                call SetError (WARNING_, INTERNAL_, &
!                    'This property name was not defined yet - ReadMercatorFile - ModuleMERCATORFormat')

            endif i2

            deallocate(Aux4D)

            if     (nDimensions == 4) then
                deallocate(Aux3D       )
                deallocate(Aux3DModulus)
            elseif (nDimensions == 3) then
                deallocate(Aux2D)
            endif

        enddo d0

        deallocate(AuxDays)


    end subroutine ReadMercatorFileV3
    

    !------------------------------------------------------------------------

    subroutine ReadMercatorFileV4(InputFile)

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:  ), pointer       :: Aux3D, Field3D
        real, dimension(:,:,:  ), pointer       :: AuxU3D, AuxV3D
        real, dimension(:,:    ), pointer       :: Aux2D, Field2D
        real, dimension(:,:    ), pointer       :: AuxBaroU2D, AuxBaroV2D
        real, dimension(:,:    ), pointer       :: SumLayerTickness
        real(8)                                 :: AuxDays
        integer                                 :: n, iOut, i, j, k
        integer                                 :: ncid, status
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: VelU_ID = 0
        integer                                 :: VelV_ID = 0 
        integer                                 :: WL_ID = 0
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName
        type (T_Time)                           :: FieldTime
        real, parameter                         :: MERCATORFillVReal = 1.0e35

        !Begin----------------------------------------------------------------

        !Verifies if file exists
        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR10'

        status = nf90_get_att(ncid,NF90_Global,"field_julian_date",AuxDays)
        if (status /= nf90_noerr) stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR20'

        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR30'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr)                                                   &
                stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR40'

            if (nDimensions == 3) then
                !Temperature, Salinity, Meridional vel., Zonal vel.
                allocate(Aux3D(Me%jmax, Me%imax, Me%kmax))
                allocate(Field3D(Me%Size%ILB:Me%Size%IUB,                               &
                         Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                Field3D(:,:,:) = FillValueReal

            elseif (nDimensions == 2) then
                !Sea surface height
                allocate(Aux2D(Me%jmax, Me%imax))
                allocate(Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                Field2D(:,:) = FillValueReal

            elseif (nDimensions == 1) then
                !Grid properties already written
                cycle
            endif


i1:         if (CheckName(nameAux, MohidName)) then
                
                iOut = OutputInstants(MohidName)

                FieldTime = Me%RefDateTime + (AuxDays * 86400. + 43200.)

                if      (nDimensions == 3) then
                    
                    if (MohidName == GetPropertyName(VelocityU_) .or.                   &
                        MohidName == GetPropertyName(VelocityV_) .or.                   &
                        MohidName == GetPropertyName(Temperature_) .or.                 &
                        MohidName == GetPropertyName(Salinity_)) then

                        status = NF90_GET_VAR(ncid,n,Aux3D)
                        if (status /= nf90_noerr)                                       &
                            stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR50'

                        !The boundary cells are not read
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then

                                Field3D(i, j, k) = Aux3D(j+1,i+1,Me%WorkSize%KUB+1-k)

                                if (Field3D(i, j, k) >= MERCATORFillVReal/2.) then

                                    Field3D(i, j, k) = 0.

                                endif

                            endif

                        enddo
                        enddo
                        enddo
                    
                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux3D = Field3D)

                        if (MohidName == GetPropertyName(VelocityU_)) then 
                            VelU_ID = n
                        endif

                        if (MohidName == GetPropertyName(VelocityV_)) then 
                            VelV_ID = n
                        endif

                    endif

                else if (nDimensions == 2) then

                    if (MohidName == GetPropertyName(WaterLevel_ ) .or.                 &
                        MohidName == GetPropertyName(BarotropicVelocityU_) .or.         &
                        MohidName == GetPropertyName(BarotropicVelocityV_)) then

                        status = NF90_GET_VAR(ncid,n,Aux2D)
                        if (status /= nf90_noerr)                                       &
                            stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR60'

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, Me%WorkSize%KUB) == WaterPoint)  &
                                then
                            
                                Field2D(i, j) = Aux2D(j+1,i+1)

                            endif

                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D = Field2D)

                        if (MohidName == GetPropertyName(WaterLevel_)) WL_ID = n

                    endif

                endif

            endif i1

            if     (nDimensions == 3) then
                deallocate(Aux3D       )
                deallocate(Field3D     )
            elseif (nDimensions == 2) then
                deallocate(Aux2D)
                deallocate(Field2D)
            endif

        enddo d0

        !calculate barotropic velocities
        if (Me%ComputeBarotropicVel) then

            if ((WL_ID > 0) .and. (VelU_ID > 0) .and. (VelV_ID > 0)) then

                allocate(Aux2D(Me%jmax, Me%imax))
                allocate(AuxU3D(Me%jmax, Me%imax, Me%kmax))
                allocate(AuxV3D(Me%jmax, Me%imax, Me%kmax))
               
                allocate(AuxBaroU2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                allocate(AuxBaroV2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                allocate(SumLayerTickness(Me%Size%ILB:Me%Size%IUB,                      &
                         Me%Size%JLB:Me%Size%JUB))

                AuxBaroU2D(:,:) = FillValueReal
                AuxBaroV2D(:,:) = FillValueReal

                status = NF90_GET_VAR(ncid,WL_ID,Aux2D)
                if (status /= nf90_noerr)                                               &
                    stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR70'

                status = NF90_GET_VAR(ncid,VelU_ID,AuxU3D)
                if (status /= nf90_noerr)                                               &
                    stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR80'

                status = NF90_GET_VAR(ncid,VelV_ID,AuxV3D)  
                if (status /= nf90_noerr)                                               &
                    stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR90'

                SumLayerTickness(:,:) = 0.

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                    if (k == Me%WorkSize%KLB) then

                        SumLayerTickness(i, j) = 0.
                        AuxBaroU2D(i,j) = 0.
                        AuxBaroV2D(i,j) = 0.

                    endif

                    if (Me%WaterPoints3D(i,j,k) == WaterPoint) then

                        if (AuxU3D(j+1,i+1,Me%WorkSize%KUB+1-k) >= MERCATORFillVReal/2.)    &
                            then

                            AuxU3D(j+1,i+1, Me%WorkSize%KUB+1-k) = 0.

                        endif

                        if (AuxV3D(j+1,i+1,Me%WorkSize%KUB+1-k) >= MERCATORFillVReal/2.)    &
                            then

                            AuxV3D(j+1,i+1, Me%WorkSize%KUB+1-k) = 0.

                        endif

                        if (k /= Me%WorkSize%KUB) then

                            AuxBaroU2D(i,j) = AuxBaroU2D(i,j) +                         &
                                              AuxU3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                              *Me%LayersThickness(k)
                            AuxBaroV2D(i,j) = AuxBaroV2D(i,j) +                         &
                                              AuxV3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                              *Me%LayersThickness(k)

                            SumLayerTickness(i, j) = SumLayerTickness(i, j) +           &
                                Me%LayersThickness(k)

                        else

                            AuxBaroU2D(i,j) = AuxBaroU2D(i,j) +                         &
                                              AuxU3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                              *(Me%LayersThickness(k) + Aux2D(j+1,i+1)) 
                            AuxBaroV2D(i,j) = AuxBaroV2D(i,j) +                         &
                                              AuxV3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                              *(Me%LayersThickness(k) + Aux2D(j+1,i+1)) 

                            SumLayerTickness(i, j) = SumLayerTickness(i, j) +           &
                                              (Me%LayersThickness(k) + Aux2D(j+1, i+1))

                            AuxBaroU2D(i,j) = AuxBaroU2D(i,j)/SumLayerTickness(i,j)
                            AuxBaroV2D(i,j) = AuxBaroV2D(i,j)/SumLayerTickness(i,j)

                        endif

                    endif

                enddo
                enddo
                enddo

                call WriteHDF5Field(FieldTime,                                          &
                                    GetPropertyName(BarotropicVelocityU_),              &
                                    iOut, Aux2D = AuxBaroU2D)

                call WriteHDF5Field(FieldTime,                                          &
                                    GetPropertyName(BarotropicVelocityV_),              &
                                    iOut, Aux2D = AuxBaroV2D)

                deallocate(Aux2D)
                deallocate(AuxU3D)
                deallocate(AuxV3D)
                deallocate(AuxBaroU2D)
                deallocate(AuxBaroV2D)
                deallocate(SumLayerTickness)

            else

                write(*,*)'Unable to calculate batropic velocity.'
                write(*,*)'Check if velocity components and water level are available.'
                stop 'ReadMercatorFileV4 - ModuleMERCATORFormat - ERR100'
            
            endif

        endif

    end subroutine ReadMercatorFileV4
    

    subroutine ReadMercatorFileV6(InputFile)

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        
        
        integer, dimension(:,:,:,:  ), pointer  :: Aux4D,Aux4D_V
        integer, dimension(:,:,:    ), pointer  :: Aux3D
        
        real, dimension(:,:,:  ), pointer       :: Field3D,AuxU3D, AuxV3D
        real, dimension(:,:    ), pointer       :: Aux2D, Field2D
        real, dimension(:,:    ), pointer       :: AuxBaroU2D, AuxBaroV2D
        real, dimension(:,:    ), pointer       :: SumLayerTickness
        real(8), dimension(:   ), allocatable   :: auxdays
        real(8)                                 :: AddOffSet,ScaleFactor
        real(8)                                 :: Addoffset_U, ScaleFactor_U
        real(8)                                 :: Addoffset_V, ScaleFactor_V
        real(8)                                 :: Addoffset_ssh, ScaleFactor_ssh
        integer                                 :: n, iOut, i, j, k,ni
        integer                                 :: ncid, status
        integer                                 :: nDimensions,ninst,dimid
        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: VelU_ID = 0
        integer                                 :: VelV_ID = 0 
        integer                                 :: WL_ID = 0
        character (len=80)                      :: nameAux, Unity
        character(Len=StringLength)             :: MohidName
        type (T_Time)                           :: FieldTime
        real, parameter                         :: MERCATORFillVReal    = 1.0e35
        integer                                 :: MERCATORFillVInteger, MERCATORFillVInteger_ssh
        integer                                 :: MERCATORFillVInteger_U, MERCATORFillVInteger_V
        
       
        !Begin----------------------------------------------------------------

        !Verifies if file exists
        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR10'

        status=NF90_INQ_DIMID(ncid,"time",dimid)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR20'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = nInst)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR30'

        allocate(AuxDays(1:nInst))

        status = nf90_inq_varid(ncid, 'time', n)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR40'

        status = NF90_GET_VAR(ncid,n,AuxDays)
        if (status /= nf90_noerr) stop 'ReadMercatorFile - ModuleMERCATORFormat - ERR50'
        
        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR30'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr)                                                   &
                stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR40'

            if (nDimensions == 4) then
                !Temperature, Salinity, Meridional vel., Zonal vel.
                allocate(Aux4D(Me%jmax, Me%imax, Me%kmax,Ninst))
                allocate(Field3D(Me%Size%ILB:Me%Size%IUB,                               &
                         Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                Field3D(:,:,:) = FillValueReal

            elseif (nDimensions == 3) then
                !Sea surface height
                allocate(Aux3D(Me%jmax, Me%imax,Ninst))
                allocate(Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                Field2D(:,:) = FillValueReal
 
            elseif (nDimensions == 2) then
                !Grid properties already written
                cycle

            elseif (nDimensions == 1) then
                !Grid properties already written
                cycle
            endif


i1:         if (CheckName(nameAux, MohidName)) then


                status=NF90_GET_ATT(ncid,n,"_FillValue",MERCATORFillVInteger)
                if (status /= nf90_noerr)                                                   &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR81'

                
                status=NF90_GET_ATT(ncid,n,"add_offset",AddOffSet)
                if (status /= nf90_noerr)                                                   &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR81'
                
                status=NF90_GET_ATT(ncid,n,"scale_factor",ScaleFactor)
                if (status /= nf90_noerr)                                                   &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR82'

                unity=' '
               
                status=NF90_GET_ATT(ncid,n,"unit_long",unity)
                if (status /= nf90_noerr) status=NF90_GET_ATT(ncid,n,"units_long",unity)
                    if (status /= nf90_noerr)                                                   &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR82a'
 
                if( trim(unity) .eq. 'Kelvin' .and. trim(nameAux) .eq. 'temperature')  AddOffSet = AddOffSet - 273.15 

                do ni = 1, nInst                
                
                    iOut = OutputInstants(MohidName)

                    !FieldTime = Me%RefDateTime + (AuxDays * 86400. + 43200.)
                    FieldTime = Me%RefDateTime + (AuxDays(ni)/24. * 86400. + 43200.) !time is in hours

                    if      (nDimensions == 4) then
                        
                        if (MohidName == GetPropertyName(VelocityU_) .or.                   &
                            MohidName == GetPropertyName(VelocityV_) .or.                   &
                            MohidName == GetPropertyName(Temperature_) .or.                 &
                            MohidName == GetPropertyName(Salinity_)) then

                            status = NF90_GET_VAR(ncid,n,Aux4D)
                            if (status /= nf90_noerr)                                       &
                                stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR50'

                            !The boundary cells are not read
                            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                                if (Me%WaterPoints3D(i, j, k) == WaterPoint) then

                                    if (Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) == MERCATORFillVInteger) then
                                        
                                        !Field3D(i, j, k) = MERCATORFillVReal
                                        Field3D(i, j, k) = FillValueReal
                                    else
                                    
                                        Field3D(i, j, k) = float(Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni))  * ScaleFactor + AddOffSet
                                    
                                    endif
                                    
                                    
!                                    if (Field3D(i, j, k) >= MERCATORFillVReal/2.) then
!                                    
!
!                                        !Field3D(i, j, k) = 0.
!                                        Field3D(i, j, k) = FillValueReal
!                                        
!                                    endif

                                endif

                            enddo
                            enddo
                            enddo
                        
                            call WriteHDF5Field(FieldTime, MohidName, iOut, Aux3D = Field3D)

                            if (MohidName == GetPropertyName(VelocityU_)) then 
                                VelU_ID = n
                                Addoffset_U = Addoffset
                                ScaleFactor_U = ScaleFactor
                                MERCATORFillVInteger_U = MERCATORFillVInteger
                            endif

                            if (MohidName == GetPropertyName(VelocityV_)) then 
                                VelV_ID = n
                                Addoffset_V = Addoffset
                                ScaleFactor_V = ScaleFactor
                                MERCATORFillVInteger_V = MERCATORFillVInteger
                            endif

                        endif

                    else if (nDimensions == 3) then

                        if (MohidName == GetPropertyName(WaterLevel_ ) .or.                 &
                            MohidName == GetPropertyName(BarotropicVelocityU_) .or.         &
                            MohidName == GetPropertyName(BarotropicVelocityV_)) then

                                status = NF90_GET_VAR(ncid,n,Aux3D)

                            if (status /= nf90_noerr)                                       &
                                stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR60'

                                !The boundary cells are not read
                            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                                if (Me%WaterPoints3D(i, j, Me%WorkSize%KUB) == WaterPoint) then

                                    if (Aux3D(j+1,i+1,ni) == MERCATORFillVInteger) then
                                        
                                        !Field2D(i, j) = MERCATORFillVReal
                                        Field2D(i, j) = FillValueReal
                                    else
                                    
                                        Field2D(i, j) = float(Aux3D(j+1,i+1,ni)) * ScaleFactor + AddOffSet
                                    
                                    endif
                                    
                                    !if (Field2D(i, j) >= MERCATORFillVReal/2.) then
                                    if (Field2D(i, j) <= FillValueReal/2.) then

                                        Field2D(i, j) = 0.
                                        
                                    endif

                                endif

                            enddo
                            enddo

                            call WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D = Field2D)

                            if (MohidName == GetPropertyName(WaterLevel_))then
                            
                                WL_ID = n
                                Addoffset_ssh = Addoffset
                                ScaleFactor_ssh = ScaleFactor
                                MERCATORFillVInteger_ssh = MERCATORFillVInteger
                            
                            endif

                        endif

                    endif
                
                enddo

            endif i1

            if     (nDimensions == 4) then
                deallocate(Aux4D       )
                deallocate(Field3D     )

            else if (nDimensions == 3) then
                deallocate(Aux3D       )
                deallocate(Field2D     )
            endif

        enddo d0

        deallocate(auxdays)


        !calculate barotropic velocities
        if (Me%ComputeBarotropicVel) then

            if ((WL_ID > 0) .and. (VelU_ID > 0) .and. (VelV_ID > 0)) then


                allocate(AuxU3D(Me%jmax, Me%imax, Me%kmax))
                allocate(AuxV3D(Me%jmax, Me%imax, Me%kmax))
                allocate(Aux4D  (Me%jmax, Me%imax, Me%kmax,Ninst))
                allocate(Aux4D_V(Me%jmax, Me%imax, Me%kmax,Ninst))
                allocate(Aux3D(Me%jmax, Me%imax,ninst))
                allocate(Aux2D(Me%jmax, Me%imax)) 
               
                allocate(AuxBaroU2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                allocate(AuxBaroV2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                allocate(SumLayerTickness(Me%Size%ILB:Me%Size%IUB,                      &
                         Me%Size%JLB:Me%Size%JUB))

                AuxBaroU2D(:,:) = FillValueReal
                AuxBaroV2D(:,:) = FillValueReal

                status = NF90_GET_VAR(ncid,VelU_ID,Aux4D) !u
                if (status /= nf90_noerr)                                               &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR80'

                status = NF90_GET_VAR(ncid,VelV_ID,Aux4D_V) !v  
                if (status /= nf90_noerr)                                               &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR90'

                status = NF90_GET_VAR(ncid,WL_ID,Aux3D)  !ssh
                if (status /= nf90_noerr)                                               &
                    stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR70'

                SumLayerTickness(:,:) = 0.
            
                do ni=1, ninst
                
                    iOut = ni
                    
                    do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Aux3D(j+1,i+1,ni) == MERCATORFillVInteger_ssh) then
                                            
                                Aux2D(j+1,i+1) = 0.
                                    
                            else
                                    
                                Aux2D(j+1,i+1) = float(Aux3D(j+1,i+1,ni))  * ScaleFactor_ssh + AddOffSet_ssh     
                                        
                        endif
                        
                    enddo
                    enddo

                  
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                        if (k == Me%WorkSize%KLB) then

                            SumLayerTickness(i, j) = 0.
                            AuxBaroU2D(i,j) = 0.
                            AuxBaroV2D(i,j) = 0.

                        endif

    !                    if (Me%WaterPoints3D(i,j,k) == WaterPoint) then
    !
    !                        if (AuxU4D(j+1,i+1,Me%WorkSize%KUB+1-k) >= MERCATORFillVReal/2.)    &
    !                            then
    !
    !                            AuxU4D(j+1,i+1, Me%WorkSize%KUB+1-k) = 0.
    !
    !                        endif
    !
    !                        if (AuxV4D(j+1,i+1,Me%WorkSize%KUB+1-k) >= MERCATORFillVReal/2.)    &
    !                            then
    !
    !                            AuxV4D(j+1,i+1, Me%WorkSize%KUB+1-k) = 0.
    !
    !                        endif

                        if (Me%WaterPoints3D(i, j, k) == WaterPoint) then
                                
                            if (Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) == MERCATORFillVInteger_U) then
                                            
                                AuxU3D(j+1,i+1, Me%WorkSize%KUB+1-k) = 0.
                                    
                            else
                                    
                                AuxU3D(j+1,i+1, Me%WorkSize%KUB+1-k)=float(Aux4D(j+1,i+1, Me%WorkSize%KUB+1-k,ni))  *  &
                                                ScaleFactor_U + AddOffSet_U     
                                        
                            endif
                                        
                            if (Aux4D_V(j+1,i+1,Me%WorkSize%KUB+1-k,ni) == MERCATORFillVInteger_V) then
                                            
                                AuxV3D(j+1,i+1, Me%WorkSize%KUB+1-k) = 0.
                                    
                            else
                                    
                                AuxV3D(j+1,i+1, Me%WorkSize%KUB+1-k)=float(Aux4D_V(j+1,i+1, Me%WorkSize%KUB+1-k,ni))  * &
                                                ScaleFactor_V + AddOffSet_V     
                                        
                            endif
                                        
                            if (k /= Me%WorkSize%KUB) then

                                AuxBaroU2D(i,j) = AuxBaroU2D(i,j) +                         &
                                                  AuxU3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                                  *Me%LayersThickness(k)
                                AuxBaroV2D(i,j) = AuxBaroV2D(i,j) +                         &
                                                  AuxV3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                                  *Me%LayersThickness(k)

                                SumLayerTickness(i, j) = SumLayerTickness(i, j) +           &
                                    Me%LayersThickness(k)

                            else

                                AuxBaroU2D(i,j) = AuxBaroU2D(i,j) +                         &
                                                  AuxU3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                                  *(Me%LayersThickness(k) + Aux2D(j+1,i+1)) 
                                AuxBaroV2D(i,j) = AuxBaroV2D(i,j) +                         &
                                                  AuxV3D(j+1,i+1,Me%WorkSize%KUB+1-k)       &
                                                  *(Me%LayersThickness(k) + Aux2D(j+1,i+1)) 

                                SumLayerTickness(i, j) = SumLayerTickness(i, j) +           &
                                                  (Me%LayersThickness(k) + Aux2D(j+1, i+1))

                                AuxBaroU2D(i,j) = AuxBaroU2D(i,j)/SumLayerTickness(i,j)
                                AuxBaroV2D(i,j) = AuxBaroV2D(i,j)/SumLayerTickness(i,j)

                            endif

                        endif

                    enddo
                    enddo
                    enddo

                    call WriteHDF5Field(FieldTime,                                          &
                                        GetPropertyName(BarotropicVelocityU_),              &
                                        iOut, Aux2D = AuxBaroU2D)

                    call WriteHDF5Field(FieldTime,                                          &
                                        GetPropertyName(BarotropicVelocityV_),              &
                                        iOut, Aux2D = AuxBaroV2D)

                enddo                

                deallocate(Aux2D)
                deallocate(AuxU3D)
                deallocate(AuxV3D)
                deallocate(AuxBaroU2D)
                deallocate(AuxBaroV2D)
                deallocate(SumLayerTickness)
                deallocate(Aux4D)
                deallocate(Aux4D_V)
                deallocate(Aux3D)

           else

                write(*,*)'Unable to calculate batropic velocity.'
                write(*,*)'Check if velocity components and water level are available.'
                stop 'ReadMercatorFileV6 - ModuleMERCATORFormat - ERR100'
            
            endif

        endif
        
    end subroutine ReadMercatorFileV6

    subroutine ReadMercatorFileV5(InputFile)

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:,:), pointer       :: Aux4D, BaroAux4D
        real, dimension(:,:,:  ), pointer       :: Aux3D, Aux3DModulus
        real, dimension(:,:,:  ), pointer       :: BaroAux3D, SumScaFactor
        real, dimension(:,:    ), pointer       :: Aux2D, BaroAux2D
        real(8), dimension(:    ), allocatable  :: AuxDays
        integer                                 :: ni, n, nInst, iOut, i, j, k
        integer                                 :: ncid, status, dimid
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: WestON, EastON, SouthON, NorthON
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName
        type (T_Time)                           :: FieldTime
        real                                    :: AddOffSet    = 0.
        real                                    :: ScaleFactor  = 1.
        real                                    :: WestValue, EastValue
        real                                    :: SouthValue, NorthValue
    
        !Begin----------------------------------------------------------------

        !Verifies if file exists

        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR10'

        status=NF90_INQ_DIMID(ncid,"time_counter",dimid)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR20'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = nInst)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR30'

        allocate(AuxDays(1:nInst))

        status = nf90_inq_varid(ncid, 'time_counter', n)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR40'

        status = NF90_GET_VAR(ncid,n,AuxDays)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR50'

        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr)                                                       &
            stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR60'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr)                                                   &
                stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR70'

            if (nDimensions == 4) then

!                if (CheckNameV3(nameAux, MohidName)) then
                !Temperature, Salinity, Meridional vel., Zonal vel.
                allocate(Aux4D(Me%jmax, Me%imax, Me%kmax, nInst))
                allocate(Aux3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,        &
                        Me%Size%KLB:Me%Size%KUB))
                Aux3D(:,:,:) = FillValueReal

                allocate(Aux3DModulus(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, &
                        Me%Size%KLB:Me%Size%KUB))
                Aux3DModulus(:,:,:) = FillValueReal
!                endif

            elseif (nDimensions == 3) then

!                if (CheckNameV3(nameAux, MohidName)) then
                    !Sea surface height
                    allocate(Aux4D(Me%jmax, Me%imax, 1      , nInst))
                    allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                    Aux2D(:,:) = FillValueReal
!                endif

            elseif (nDimensions == 2 .or. nDimensions == 1) then
                !Grid properties already written
                cycle
            endif

            status = NF90_GET_VAR(ncid,n,Aux4D)
            if (status /= nf90_noerr)                                                   &
                stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR80'

i2:         if (CheckNameV3(nameAux, MohidName)) then

                status=NF90_GET_ATT(ncid,n,"add_offset",AddOffSet)
                if (status /= nf90_noerr)                                                   &
                    stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR81'
                
                status=NF90_GET_ATT(ncid,n,"scale_factor",ScaleFactor)
                if (status /= nf90_noerr)                                                   &
                    stop 'ReadMercatorFileV5 - ModuleMERCATORFormat - ERR82'

d1:             do ni = 1, nInst

                    iOut = OutputInstants(MohidName)

                    FieldTime = Me%RefDateTime + (AuxDays(ni) - 43200.)

                    if      (nDimensions == 4) then
                        
                        !The boundary cells are not read
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then
                                
                                if      (MohidName == GetPropertyName(VelocityU_)) then

                                    WestON = Me%ComputeFaces3DU(i, j-1, k) 
                                    EastON = Me%ComputeFaces3DU(i, j  , k) 
                                    
                                    WestValue = Aux4D(j,  i+1,Me%WorkSize%KUB+1-k,ni) * ScaleFactor + AddOffSet
                                    EastValue = Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) * ScaleFactor + AddOffSet

                                    Aux3D(i, j, k) = (WestValue * WestON  + &
                                                      EastValue * EastON) / 2.

                                    if (abs(Aux3D(i, j, k))> 200.) Aux3D(i, j, k) = 0.

                                                    
                                else if (MohidName == GetPropertyName(VelocityV_)) then
                                    
                                    SouthON = Me%ComputeFaces3DV(i-1, j, k) 
                                    NorthON = Me%ComputeFaces3DV(i,   j, k) 
                                    
                                    SouthValue = Aux4D(j+1,i,  Me%WorkSize%KUB+1-k,ni) * ScaleFactor + AddOffSet
                                    NorthValue = Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) * ScaleFactor + AddOffSet
                                    
                                    Aux3D(i, j, k) = (SouthValue * SouthON  + &
                                                      NorthValue * NorthON) / 2.

                                   if (abs(Aux3D(i, j, k))> 200.) Aux3D(i, j, k) = 0.

                                 
                                else if (MohidName == GetPropertyName(Temperature_) .or. &
                                         MohidName == GetPropertyName(Salinity_   )) then

                                    Aux3D(i, j, k) = Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni)* ScaleFactor + AddOffSet
                                    
                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        
                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux3D = Aux3D)


                    else if (nDimensions == 3) then

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            if (Me%WaterPoints3D(i, j, Me%WorkSize%KUB) == WaterPoint) then
                                
                                if (MohidName == GetPropertyName(BarotropicVelocityU_)) then

                                    WestON = Me%ComputeFaces3DU(i, j-1, Me%WorkSize%KUB) 
                                    EastON = Me%ComputeFaces3DU(i, j  , Me%WorkSize%KUB) 
                                    WestValue = Aux4D(j  ,i+1,1,ni) * ScaleFactor + AddOffSet
                                    EastValue = Aux4D(j+1,i+1,1,ni) * ScaleFactor + AddOffSet

                                    Aux2D(i, j) = (WestValue * WestON  + &
                                                   EastValue * EastON) / 2.

                                    if (abs(Aux2D(i, j))> 200.) Aux2D(i, j) = 0.

                                else if (MohidName == GetPropertyName(BarotropicVelocityV_)) then
                                    
                                    SouthON = Me%ComputeFaces3DV(i-1, j, Me%WorkSize%KUB) 
                                    NorthON = Me%ComputeFaces3DV(i,   j, Me%WorkSize%KUB) 
                                    SouthValue = Aux4D(j+1,  i  ,1,ni) * ScaleFactor + AddOffSet
                                    NorthValue = Aux4D(j+1,i+1  ,1,ni) * ScaleFactor + AddOffSet
                                
                                    Aux2D(i, j) = (SouthValue * SouthON  + &
                                                   NorthValue * NorthON) / 2.

                                    if (abs(Aux2D(i, j))> 200.) Aux2D(i, j) = 0.
                                

                                else if (MohidName == GetPropertyName(WaterLevel_ )) then

                                    Aux2D(i, j) = Aux4D(j+1,i+1,1, ni) * ScaleFactor + AddOffSet
                                    
                                endif
                            endif

                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D = Aux2D)

                    endif


                enddo d1

                if (MohidName == GetPropertyName(VelocityU_) .and.                      &
                    Me%ComputeBarotropicVel) then

                    allocate(BaroAux4D(Me%jmax, Me%imax, Me%kmax , nInst))
                    allocate(BaroAux2D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB))
                    allocate(BaroAux3D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB,                         &
                                       nInst))
                    allocate(SumScaFactor(Me%Size%ILB:Me%Size%IUB,                      &
                                          Me%Size%JLB:Me%Size%JUB,                      &
                                          nInst))
                    BaroAux2D(:,:) = FillValueReal

d2:                 do ni = 1, nInst

                        SumScaFactor(:,:,ni) = 0.

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                            if (k == Me%WorkSize%KLB) then

                                BaroAux3D(i+1,j, ni) = 0.
                                BaroAux3D(i+1,j+1, ni) = 0.

                                SumScaFactor(i+1,j, ni) = 0.
                                SumScaFactor(i+1,j+1, ni) = 0.

                            endif

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then

                                BaroAux2D(i, j) = 0.

                                WestON = Me%ComputeFaces3DU(i, j-1, k) 
                                EastON = Me%ComputeFaces3DU(i, j  , k) 

                                if ((WestON <= 1) .and. (EastON <= 1)) then

                                    BaroAux4D(j  ,i+1, Me%WorkSize%KUB+1-k,ni) =            &
                                        (Aux4D(j  ,i+1, Me%WorkSize%KUB+1-k,ni) *           & 
                                        ScaleFactor + AddOffSet) * WestON *    &
                                        Me%ScaFactorU(i+1,j, Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i+1,j, ni) = SumScaFactor(i+1,j, ni) +     &
                                        Me%ScaFactorU(i+1,j, Me%WorkSize%KUB+1-k)* WestON   

                                    BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) =            &
                                        (Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) *           &
                                        ScaleFactor + AddOffSet) * EastON*                  &
                                        Me%ScaFactorU(i+1,j+1,Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i+1,j+1, ni) = SumScaFactor(i+1,j+1,ni) +  &
                                        Me%ScaFactorU(i+1,j+1, Me%WorkSize%KUB+1-k)* EastON   

                                    BaroAux3D(i+1,j, ni) = BaroAux3D(i+1,j, ni) +           &
                                        BaroAux4D(j,i+1, Me%WorkSize%KUB+1-k, ni)

                                    BaroAux3D(i+1,j+1, ni) = BaroAux3D(i+1,j+1, ni) +       &
                                        BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni)

                                    if ((k == Me%WorkSize%KUB) .and.                        &
                                       (SumScaFactor(i+1,j, ni) > 0. .and.                 &
                                       SumScaFactor(i+1,j+1, ni) > 0.)) then

                                        BaroAux2D(i,j) =                                    &
                                        (BaroAux3D(i+1,j, ni)/SumScaFactor(i+1,j, ni) +     &
                                        BaroAux3D(i+1,j+1, ni)/SumScaFactor(i+1,j+1, ni)) / 2.

                                        if (abs(BaroAux2D(i, j))> 200.) BaroAux2D(i, j) = 0.
                                    endif

                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime,                                  &
                                            GetPropertyName(BarotropicVelocityU_),      &
                                            iOut, Aux2D = BaroAux2D)

                    enddo d2

                    deallocate(BaroAux4D)
                    deallocate(BaroAux3D)
                    deallocate(BaroAux2D)
                    deallocate(SumScaFactor)

                endif

                if (MohidName == GetPropertyName(VelocityV_) .and.                      &
                    Me%ComputeBarotropicVel) then

                    allocate(BaroAux4D(Me%jmax, Me%imax, Me%kmax , nInst))
                    allocate(BaroAux2D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB))
                    allocate(BaroAux3D(Me%Size%ILB:Me%Size%IUB,                         &
                                       Me%Size%JLB:Me%Size%JUB,                         &
                                       nInst))
                    allocate(SumScaFactor(Me%Size%ILB:Me%Size%IUB,                      &
                                          Me%Size%JLB:Me%Size%JUB,                      &
                                          nInst))
                    BaroAux2D(:,:) = FillValueReal

d3:                 do ni = 1, nInst

                        SumScaFactor(:,:,ni) = 0.

                        !The boundary cells are not read
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                            if (k == Me%WorkSize%KLB) then

                                BaroAux3D(i,j+1,ni) = 0.
                                BaroAux3D(i+1,j+1, ni) = 0.

                                SumScaFactor(i, j+1, ni) = 0.
                                SumScaFactor(i+1,j+1, ni) = 0.

                            endif

                            if (Me%WaterPoints3D(i, j, k) == WaterPoint) then

                                BaroAux2D(i, j) = 0.

                                SouthON = Me%ComputeFaces3DV(i-1, j, k) 
                                NorthON = Me%ComputeFaces3DV(i,   j, k) 

                                if ((SouthON <= 1) .and. (NorthON <= 1)) then

                                    BaroAux4D(j+1,i ,Me%WorkSize%KUB+1-k,ni) =              &
                                        (Aux4D(j+1,i ,Me%WorkSize%KUB+1-k,ni)               & 
                                        * ScaleFactor + AddOffSet) * SouthON*               &
                                        Me%ScaFactorV(i, j+1,Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i, j+1, ni) = SumScaFactor(i, j+1, ni) +   &
                                        Me%ScaFactorV(i, j+1, Me%WorkSize%KUB+1-k)* SouthON   

                                    BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) =            &
                                        (Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni) *           &
                                        ScaleFactor + AddOffSet) * NorthON *                &
                                        Me%ScaFactorV(i+1,j+1,Me%WorkSize%KUB+1-k)

                                    SumScaFactor(i+1,j+1, ni) = SumScaFactor(i+1,j+1, ni) + &
                                        Me%ScaFactorV(i+1,j+1, Me%WorkSize%KUB+1-k)* NorthON   

                                    BaroAux3D(i,j+1,ni) = BaroAux3D(i,j+1,ni) +             &
                                        BaroAux4D(j+1,i ,Me%WorkSize%KUB+1-k,ni)

                                    BaroAux3D(i+1,j+1, ni) = BaroAux3D(i+1,j+1, ni) +       &
                                        BaroAux4D(j+1,i+1,Me%WorkSize%KUB+1-k, ni)

                                    if ((k == Me%WorkSize%KUB) .and.                        &
                                        (SumScaFactor(i, j+1, ni) > 0. .and.                 &
                                        SumScaFactor(i+1,j+1, ni) > 0.)) then

                                        BaroAux2D(i,j) =                                    &
                                        (BaroAux3D(i,j+1,ni)/SumScaFactor(i, j+1, ni) +     &
                                        BaroAux3D(i+1,j+1, ni)/SumScaFactor(i+1,j+1, ni)) / 2.
                                
                                        if (abs(BaroAux2D(i, j))> 200.) BaroAux2D(i, j) = 0.
                                    endif

                                endif

                            endif

                        enddo
                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, GetPropertyName(BarotropicVelocityV_), &
                                            iOut, Aux2D = BaroAux2D)

                    enddo d3

                    deallocate(BaroAux4D)
                    deallocate(BaroAux3D)
                    deallocate(BaroAux2D)
                    deallocate(SumScaFactor)

                endif

!            else i2

!                call SetError (WARNING_, INTERNAL_, &
!                    'This property name was not defined yet - ReadMercatorFile - ModuleMERCATORFormat')

            endif i2

            deallocate(Aux4D)

            if     (nDimensions == 4) then
                deallocate(Aux3D       )
                deallocate(Aux3DModulus)
            elseif (nDimensions == 3) then
                deallocate(Aux2D)
            endif

        enddo d0

        deallocate(AuxDays)


    end subroutine ReadMercatorFileV5
    

    !------------------------------------------------------------------------


    subroutine MapFromMercatorFile(InputFile, MapDone)

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        logical                                 :: MapDone
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:,:), pointer       :: Aux4D
        real(8), dimension(:    ), allocatable  :: AuxDays
        integer                                 :: ni, n, nInst, i, j, k
        integer                                 :: ncid, status, dimid
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: MohidName


        !Begin----------------------------------------------------------------

        !Verifies if file exists

        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR10'

        status=NF90_INQ_DIMID(ncid,"time",dimid)
        if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR60'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = nInst)
        if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR70'

        allocate(AuxDays(1:nInst))

        status = nf90_inq_varid(ncid, 'time', n)
        if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR100'

        status = NF90_GET_VAR(ncid,n,AuxDays)
        if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR110'

        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr) stop 'OpenAndReadMERCATORFileV1 - ModuleMERCATORFormat - ERR04'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR04b'

i2:         if (CheckName(nameAux, MohidName)) then
                
                if (MohidName /= GetPropertyName(Temperature_)) cycle

                !First instant 
                ni = 1 
                !Temperature
                allocate(Aux4D(Me%jmax, Me%imax, Me%kmax, ni))

                status = NF90_GET_VAR(ncid,n,Aux4D)
                if (status /= nf90_noerr) stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR70'



                !The boundary cells are not read
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Aux4D(j+1,i+1,Me%WorkSize%KUB+1-k,ni) /= MissingValue) then
                        Me%WaterPoints3D  (i, j, k) = 1
                    endif

                enddo
                enddo
                enddo


                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%WaterPoints3D  (i, j, k) == 0) then
                        
                        Me%ComputeFaces3DV(i-1,   j,   k) = 0
                        Me%ComputeFaces3DV(i    , j,   k) = 0
                    endif

                    if (Me%WaterPoints3D  (i, j, k) == 0) then
                        
                        Me%ComputeFaces3DU(i,   j-1,   k) = 0
                        Me%ComputeFaces3DU(i,   j    , k) = 0

                    endif

                enddo
                enddo
                enddo

                deallocate(Aux4D)

                MapDone = .true.

            endif i2

        enddo d0

        status=NF90_CLOSE(ncid)
        if (status /= nf90_noerr)  stop 'MapFromMercatorFile - ModuleMERCATORFormat - ERR250'

        deallocate(AuxDays)

    end subroutine MapFromMercatorFile
   
    !------------------------------------------------------------------------

    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        type(T_Size2D)                              :: WorkSize2D
        real                                        :: Xorig, Yorig, Latitude, Longitude
        logical                                     :: ContinuesCompute
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'

        WorkSize2D%ILB = Me%WorkSize%ILB
        WorkSize2D%JLB = Me%WorkSize%JLB

        WorkSize2D%IUB = Me%WorkSize%IUB
        WorkSize2D%JUB = Me%WorkSize%JUB

        ContinuesCompute = .FALSE.

!       Os dados de MERCATOR nao tenhen informaciao de niveis:
        allocate(Me%SurfaceElevation(Me%Size%ILB:Me%Size%IUB,&
                                  Me%Size%JLB:Me%Size%JUB))
        Me%SurfaceElevation = 0.0

        Xorig = Me%CenterX(1, Me%WorkSize%JLB) - &
               (Me%XX(Me%WorkSize%JLB + 1) - Me%XX(Me%WorkSize%JLB)) / 2.
        
        Yorig = Me%CenterY(Me%WorkSize%ILB, 1) - &
               (Me%YY(Me%WorkSize%ILB + 1) - Me%YY(Me%WorkSize%ILB)) / 2.

        call WriteGridData (FileName        = Me%GridFileName,              &
                            XX              = Me%XX,                        &
                            YY              = Me%YY,                        &
                            COMENT1         = " MERCATOR Grid based on file :",  &
                            COMENT2         = Me%FileName,                  &
                            WorkSize        = WorkSize2D,                   & 
                            CoordType       = 4,                            &
                            Xorig           = Xorig,                        &
                            Yorig           = Yorig,                        &
                            Zone            = 0,                            &
                            Grid_Angle      = 0.,                           &
                            Latitude        = Latitude,                     &
                            Longitude       = Longitude,                    &
                            GridData2D_Real = Me%Bathymetry,                &
                            Overwrite       = ON,                           &
                            FillValue       = -99.,                            &
                            STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMERCATORFormat - ERR01'


        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMERCATORFormat - ERR02b'


        call ConstructGridData(GridDataID       = Me%ObjGridData,           &
                               HorizontalGridID = Me%ObjHorizontalGrid,     &
                               Filename         = Me%GridFileName,          &
                               STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR03'

        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjGridData, Me%ObjHorizontalGrid,   &
                                      Me%FirstField%Date, STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR04'

        call ConstructGeometry(Me%ObjGeometry, Me%ObjGridData, Me%ObjHorizontalGrid,    &
                                    Me%ObjHorizontalMap, Me%FirstField%Date,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR05'

        call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,                  &
                                    WaterPoints3D    = Me%WaterPoints3D,                &
                                    SurfaceElevation = Me%SurfaceElevation,             &
                                    ContinuesCompute = ContinuesCompute,                &
                                    NonHydrostatic   = .false.,                         &
                                    ActualTime       = Me%FirstField%Date,              &
                                    STAT             = STAT_CALL )
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenAndReadFatherHDF5File -  ModuleInterpolateGrids - ERR300'


        !SZZ
        call GetGeometryDistances (Me%ObjGeometry, SZZ = Me%SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR07'

        call GetGeometryAreas(Me%ObjGeometry, Me%AreaU, Me%AreaV, Me%FirstField%Date, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR10'

        call ConstructMap(Map_ID           = Me%ObjMap,             &   
                          GeometryID       = Me%ObjGeometry,        &
                          HorizontalMapID  = Me%ObjHorizontalMap,   & 
                          TimeID           = Me%ObjTime,            &
                          STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR08'

        call UpdateComputeFaces3D(Me%ObjMap, Me%SurfaceElevation, Me%FirstField%Date, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR09'

        call GetComputeFaces3D(Me%ObjMap, Me%ComputeFacesU3D, Me%ComputeFacesV3D, Me%ComputeFacesW3D, &
                                 Me%FirstField%Date, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridData - ModuleMERCATORFormat - ERR10'

    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------

    subroutine ConstructGridV2
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL !, UnitID, i, j
        type(T_Size2D)              :: WorkSize2D
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'


        WorkSize2D%ILB = Me%WorkSize%ILB
        WorkSize2D%JLB = Me%WorkSize%JLB

        WorkSize2D%IUB = Me%WorkSize%IUB
        WorkSize2D%JUB = Me%WorkSize%JUB

        call WriteGridData  (FileName       = Me%GridFileName,                              &
                             ConnectionX    = Me%XX_IE,                                     &
                             ConnectionY    = Me%YY_IE,                                     &
                             COMENT1        = 'Grid Data created from MERCATOR NetCDF file',&
                             COMENT2        = trim(Me%InputGridFile),                       &
                             WorkSize       = WorkSize2D,                                   &
                             CoordType      = 4,                                            &
                             Xorig          = 0.,                                           &
                             Yorig          = 0.,                                           &
                             Zone           = 29,                                           &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Me%XX_IE(1,1),                                &
                             Longitude      = Me%YY_IE(1,1),                                &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%Bathymetry,                                &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGridV2 - ModuleMERCATORFormat - ERR10'


        call WriteGridData  (FileName       = trim(Me%GridFileName)//"_MaxDepth",           &
                             ConnectionX    = Me%XX_IE,                                     &
                             ConnectionY    = Me%YY_IE,                                     &
                             COMENT1        = 'Grid Data created from MERCATOR NetCDF file',&
                             COMENT2        = trim(Me%InputGridFile)//"_MaxDepth",          &
                             WorkSize       = WorkSize2D,                                   &
                             CoordType      = 4,                                            &
                             Xorig          = 0.,                                           &
                             Yorig          = 0.,                                           &
                             Zone           = 29,                                           &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Me%XX_IE(1,1),                                &
                             Longitude      = Me%YY_IE(1,1),                                &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%BathymetryMax,                             &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGridV2 - ModuleMERCATORFormat - ERR10'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGridV2 - ModuleMERCATORFormat - ERR20'


    end subroutine ConstructGridV2

    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    
    subroutine WriteMERCATORGeometry
        
        !Local-----------------------------------------------------------------
        integer                 :: STAT_CALL, k
        
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
        write(Me%Unit,*) 'LAYERS : ', Me%LayerNumber
        write(Me%Unit,*) '<<beginlayers>>'
        do k = Me%LayerNumber, 1, -1
            write(Me%Unit,'(f9.2)') Me%LayersThickness(k)
        enddo
        write(Me%Unit,*) '<<endlayers>>'
        write(Me%Unit,*) 'MININITIALLAYERTHICKNESS : 1'
        write(Me%Unit,*) '<enddomain>'
        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteMERCATORGeometry - ModuleMERCATORFormat - ERR30'


    end subroutine WriteMERCATORGeometry

    
    !------------------------------------------------------------------------


    subroutine ConvertToMohidUnits(Field)

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: Field
        
        !Begin-----------------------------------------------------------------

        select case(trim(Field%Units))

            case('degC')

                Field%Units     = 'ºC'

            case default

        end select



    end subroutine ConvertToMohidUnits

    !------------------------------------------------------------------------

    
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


    subroutine AddDate (FirstDate, ObjDate)

        !Arguments-------------------------------------------------------------
        type (T_Date), pointer                   :: FirstDate
        type (T_Date), pointer                   :: ObjDate

        !Local-----------------------------------------------------------------
        type (T_Date), pointer                   :: NewDate
        type (T_Date), pointer                   :: PreviousDate
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewDate)
        nullify  (NewDate%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstDate)) then
            FirstDate         => NewDate
            ObjDate           => NewDate
        else
            PreviousDate      => FirstDate
            ObjDate           => FirstDate%Next
            do while (associated(ObjDate))
                PreviousDate  => ObjDate
                ObjDate       => ObjDate%Next
            enddo
            ObjDate           => NewDate
            PreviousDate%Next => NewDate
        endif


    end subroutine AddDate
    
    
    !------------------------------------------------------------------------
    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        real,    dimension(:,:,:), pointer              :: tmp
        integer                                         :: i, j, k
        integer                                         :: STAT_CALL, OutputNumber
        integer                                         :: ILB, JLB, KLB
        integer                                         :: IUB, JUB, KUB
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing HDF5 file...'

        call Open_HDF5_OutPut_File

        OutputNumber = 1
        CurrentDate  => Me%FirstDate

        do while(associated(CurrentDate))

!           Dados para escriver uma soa vez cada date:
            call ExtractDate   (CurrentDate%Date,                           &
                                AuxTime(1), AuxTime(2), AuxTime(3),         &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR01'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR02'

            !Writes SZZ
            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                                 Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR04'
    
            call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",           &
                                "m", Array3D = Me%SZZ, OutputNumber = OutPutNumber,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'

            !Writes Openpoints = WaterPoints por quanto mas pode trocar nun futuro
            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                                 Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR04'
    
            call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",                  &
                                "m", Array3D = Me%WaterPoints3D, OutputNumber = OutPutNumber,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'

            !Start: Writes ComputeComputeFacesU/V/W3D
            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                                 Me%WorkSize%JUB +1, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR04'
    
            call HDF5WriteData  (Me%ObjHDF5, "/Grid/ComputeFacesU3D", "ComputeFacesU3D",                  &
                                "m", Array3D = Me%ComputeFacesU3D, OutputNumber = OutPutNumber,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB + 1, Me%WorkSize%JLB,    &
                                 Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR04'
    
            call HDF5WriteData  (Me%ObjHDF5, "/Grid/ComputeFacesV3D", "ComputeFacesV3D",                  &
                                "m", Array3D = Me%ComputeFacesV3D, OutputNumber = OutPutNumber,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB,        &
                                 Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR04'
    
            call HDF5WriteData  (Me%ObjHDF5, "/Grid/ComputeFacesW3D", "ComputeFacesW3D",                  &
                                "m", Array3D = Me%ComputeFacesW3D, OutputNumber = OutPutNumber,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaterProperties - ERR05'
            !End: Writes ComputeComputeFacesU/V/W3D

            !Ojo aquí!!!! Esto tiene que cambiar. Tiene que pasar al bucle de abajo.
            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                                 Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR02'

        
            call HDF5WriteData   (Me%ObjHDF5, "/Results/Waterlevel", "Waterlevel", "m",       &
                                  Array2D =  Me%SurfaceElevation, OutputNumber = OutPutNumber, &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR03'            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

            Field => Me%FirstField

            do while(associated(Field))
            
                Field%OutputNumber = OutputNumber

                if(Field%Date == CurrentDate%Date)then

                ILB=Field%WorkSize%ILB
                IUB=Field%WorkSize%IUB
                JLB=Field%WorkSize%JLB
                JUB=Field%WorkSize%JUB
                KLB=Field%WorkSize%KLB
                KUB=Field%WorkSize%KUB

                if(trim(Field%Name) == 'velocity U') then

                allocate(tmp(ILB:IUB,JLB:JUB,KLB:KUB))

                
                do k=Field%WorkSize%KLB , Field%WorkSize%KUB
                do i=Field%WorkSize%ILB , Field%WorkSize%IUB
                do j=Field%WorkSize%JLB +1 , Field%WorkSize%JUB
                if(Me%ComputeFacesU3D(i,j,k)/=0) then
                    tmp(i,j,k)=(Field%Values3D(i,j,k)+Field%Values3D(i,j-1,k))*0.5
                endif
                enddo
                enddo
                enddo

                Field%Values3D(ILB:IUB,JLB:JUB,KLB:KUB)=tmp*Me%AreaU(ILB:IUB,JLB:JUB,KLB:KUB)
                Field%Values3D=Field%Values3D*Me%ComputeFacesU3D
                Field%Name='WaterFluxX'
                Field%WorkSize%JUB = Field%WorkSize%JUB +1

                deallocate(tmp)

                elseif(trim(Field%Name) == 'velocity V') then

                allocate(tmp(ILB:IUB,JLB:JUB,KLB:KUB))

                do k=Field%WorkSize%KLB , Field%WorkSize%KUB
                do i=Field%WorkSize%ILB + 1 , Field%WorkSize%IUB
                do j=Field%WorkSize%JLB , Field%WorkSize%JUB
                if(Me%ComputeFacesV3D(i,j,k)/=0) then
                    tmp(i,j,k)=(Field%Values3D(i,j,k)+Field%Values3D(i-1,j,k))*0.5
                endif
                enddo
                enddo
                enddo

                Field%Values3D(ILB:IUB,JLB:JUB,KLB:KUB)=tmp*Me%AreaV(ILB:IUB,JLB:JUB,KLB:KUB)
                Field%Values3D=Field%Values3D*Me%ComputeFacesV3D
                Field%Name='WaterFluxY'
                Field%WorkSize%IUB = Field%WorkSize%IUB + 1

                deallocate(tmp)

                endif
                
            
if0:                if(Field%nDimensions == 2) then

                    call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                       Field%WorkSize%JLB, Field%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR05'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//Field%Name,                         &
                                       Field%Name,                                      &
                                       Field%Units,                                     &
                                       Array2D      = Field%Values2D,                   &
                                       OutputNumber = Field%OutputNumber,               &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR06'

                elseif(Field%nDimensions == 3) then

                    call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                       Field%WorkSize%JLB, Field%WorkSize%JUB, &
                                       Field%WorkSize%KLB, Field%WorkSize%KUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR05'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//Field%Name,                         &
                                       Field%Name,                                      &
                                       Field%Units,                                     &
                                       Array3D      = Field%Values3D,                   &
                                       OutputNumber = Field%OutputNumber,               &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR06'

                endif if0

                end if

                Field => Field%Next

            end do

            OutputNumber = OutputNumber + 1

            CurrentDate => CurrentDate%Next

        end do

        write(*,*)
        write(*,*)'Closing HDF5 file...'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMecatorFormat - ERR07'


    end subroutine OutputFields


    !------------------------------------------------------------------------
    subroutine WriteHDF5Field(FieldTime, MohidName, iOut, Aux2D, Aux3D)

        !Arguments-------------------------------------------------------------
        type (T_Time)                                   :: FieldTime
        character(Len=StringLength)                     :: MohidName
        integer                                         :: iOut

        real, dimension(:,:  ), pointer, optional       :: Aux2D
        real, dimension(:,:,:), pointer, optional       :: Aux3D

        !Local-----------------------------------------------------------------
        type (T_Time)                                   :: AuxFieldTime
        character(Len=StringLength)                     :: PropUnits, AuxName
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, nItems
        integer                                         :: WorkILB, WorkJLB, WorkKLB
        integer                                         :: WorkIUB, WorkJUB, WorkKUB
        logical                                         :: Exist

        !Begin-----------------------------------------------------------------
        
        !Bounds

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        PropUnits = MohidUnits(MohidName)

        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR10'

        if      (present(Aux2D)) then

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array2D = Aux2D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR20'

        else if (present(Aux3D)) then

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array3D = Aux3D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR30'

        else 

            stop 'WriteHDF5Field - ModuleMecatorFormat - ERR40'

        endif

        call GetHDF5GroupExist (Me%ObjHDF5, "/Time", Exist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR35'

        if (Exist) then
        
            call GetHDF5GroupNumberOfItems (Me%ObjHDF5,  "/Time", nItems, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR37'
        
        endif    


        if (.not. Exist .or. iOut > nItems) then

            !Writes current time
            call ExtractDate   (FieldTime, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                           AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR50'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                                 Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR60'

            !Writes SZZ
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                                 WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR70'

            call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",                 &
                                 "m", Array3D = Me%SZZ, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR80'

            AuxName = MohidName

        else if (iOut == nItems) then

            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR100'

            call HDF5ReadData  (Me%ObjHDF5, "/Time", "Time",                               &
                                 Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR110'


            call SetDate   (AuxFieldTime, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                          AuxTime(4), AuxTime(5), AuxTime(6))


            if (FieldTime /= AuxFieldTime) then
!                write(*,*) 'The time instants of property ',trim(MohidName)
!                write(*,*) 'are not consistent with property ',trim(AuxName)
                stop 'WriteHDF5Field - ModuleMecatorFormat - ERR120'   
            endif             
            
        endif

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMecatorFormat - ERR90'


    end subroutine WriteHDF5Field


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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR02'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR04'            
   
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = Me%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR08'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMecatorFormat - ERR09'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    logical function CheckName(MERCATORName, MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: MERCATORName
        character(Len=StringLength)     :: MohidName
        
        !Begin-----------------------------------------------------------------


        select case(trim(MERCATORName))

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

            case('ssh')

                MohidName = GetPropertyName(WaterLevel_)
                CheckName = .true.


            case('ubar')

                MohidName = GetPropertyName(BarotropicVelocityU_)
                CheckName = .true.


            case('vbar')

                MohidName = GetPropertyName(BarotropicVelocityV_)
                CheckName = .true.


            case default
                
                CheckName = .false.

        end select


    end function CheckName
    
    
    !--------------------------------------------------------------------------

    logical function CheckNameV3(MERCATORName, MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: MERCATORName
        character(Len=StringLength)     :: MohidName
        
        !Begin-----------------------------------------------------------------


        select case(trim(MERCATORName))

            case('vozocrtx')

                MohidName = GetPropertyName(VelocityU_)
                CheckNameV3 = .true.

            case('vomecrty')

                MohidName = GetPropertyName(VelocityV_)
                CheckNameV3 = .true.

            case('votemper')

                MohidName = GetPropertyName(Temperature_)
                CheckNameV3 = .true.

            case('vosaline')

                MohidName = GetPropertyName(Salinity_)
                CheckNameV3 = .true.

            case('sossheig')

                MohidName = GetPropertyName(WaterLevel_)
                CheckNameV3 = .true.


            case('ubar')

                MohidName = GetPropertyName(BarotropicVelocityU_)
                CheckNameV3 = .true.


            case('vbar')

                MohidName = GetPropertyName(BarotropicVelocityV_)
                CheckNameV3 = .true.


            case default
                
                CheckNameV3 = .false.

        end select


    end function CheckNameV3
    
    
    !--------------------------------------------------------------------------

    integer function OutputInstants(MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=StringLength)     :: MohidName
        !Local-----------------------------------------------------------------
        integer                         :: MohidID
        
        !Begin-----------------------------------------------------------------


        MohidID = GetPropertyIDNumber(MohidName)

        select case(MohidID)

            case(VelocityU_)
                
                Me%Instants(1) = Me%Instants(1) + 1
                OutputInstants = Me%Instants(1)

            case(VelocityV_)
                
                Me%Instants(2) = Me%Instants(2) + 1
                OutputInstants = Me%Instants(2)

            case(Temperature_)
                
                Me%Instants(3) = Me%Instants(3) + 1
                OutputInstants = Me%Instants(3)

            case(Salinity_)
                
                Me%Instants(4) = Me%Instants(4) + 1
                OutputInstants = Me%Instants(4)

            case(WaterLevel_)
                
                Me%Instants(5) = Me%Instants(5) + 1
                OutputInstants = Me%Instants(5)

            case(BarotropicVelocityU_)

                Me%Instants(6) = Me%Instants(6) + 1
                OutputInstants = Me%Instants(6)

            case(BarotropicVelocityV_)

                Me%Instants(7) = Me%Instants(7) + 1
                OutputInstants = Me%Instants(7)

        end select


    end function OutputInstants
    
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    character(Len=StringLength) function MohidUnits(MohidName)
        
        !Arguments-----------------------------------------------------------
        character(Len=StringLength)     :: MohidName
        !Local-----------------------------------------------------------------
        integer                         :: MohidID
        
        !Begin-----------------------------------------------------------------


        MohidID = GetPropertyIDNumber(MohidName)

        select case(MohidID)

            case(VelocityU_, VelocityV_, VelocityModulus_, BarotropicVelocityU_, BarotropicVelocityV_)
               
                MohidUnits = 'm/s'

            case(Temperature_)
                
                MohidUnits = 'ºC'

            case(Salinity_)
                
                MohidUnits = 'psu'

            case(WaterLevel_)
                
                MohidUnits = 'm'

        end select


    end function MohidUnits
    
    
    !--------------------------------------------------------------------------

    
    subroutine KillMERCATORFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------


        if      (Me%ReadOptionType == Version1) then

            call KillMap(Me%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR10'

            call KillGeometry(Me%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR20'

            call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR30'

            call KillGridData(Me%ObjGridData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR40'

            deallocate(Me%CenterX)
            deallocate(Me%XX)
            deallocate(Me%CenterY)
            deallocate(Me%YY)
            deallocate(Me%FirstField)
            deallocate(Me%Bathymetry)
            deallocate(Me%WaterPoints3D)

        else if (Me%ReadOptionType == Version2) then


            deallocate(Me%Bathymetry)
            deallocate(Me%BathymetryMax)
            deallocate(Me%XX_IE)
            deallocate(Me%YY_IE)
            deallocate(Me%WaterPoints3D)
            deallocate(Me%ComputeFaces3DU)
            deallocate(Me%ComputeFaces3DV)
            deallocate(Me%SZZ)
            deallocate(Me%LayersThickness)


        else if (Me%ReadOptionType == Version3) then


        endif

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR50'
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR60'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillMERCATORFormat - ModuleMERCATORFormat - ERR70'


        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillMERCATORFormat

    !--------------------------------------------------------------------------

end module ModuleMERCATORFormat
