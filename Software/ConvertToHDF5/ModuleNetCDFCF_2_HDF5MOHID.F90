!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : NetCDFCF_2_HDF5MOHID
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert NetCDFCF_2_HDF5MOHID files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleNetCDFCF_2_HDF5MOHID

    use ModuleTime
    use ModuleGlobalData
    use ModuleFunctions
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
    public  :: ConvertNetCDFCF_2_HDF5MOHID
    private ::      ReadOptions
    private ::      KillNetCDFCF_2_HDF5MOHID

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: input_files_begin   = '<<begin_input_files>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<end_input_files>>'

    character(LEN = StringLength), parameter    :: begin_time          = '<<begin_time>>'
    character(LEN = StringLength), parameter    :: end_time            = '<<end_time>>'
    
    character(LEN = StringLength), parameter    :: begin_grid          = '<<begin_grid>>'
    character(LEN = StringLength), parameter    :: end_grid            = '<<end_grid>>'

    character(LEN = StringLength), parameter    :: begin_field         = '<<begin_field>>'
    character(LEN = StringLength), parameter    :: end_field           = '<<end_field>>'


    integer                      , parameter    :: Real4_ = 1, Real8_ =2, Integer4_ = 3


    !Types---------------------------------------------------------------------
    
    private :: T_ValueIn
    type       T_ValueIn    
        integer                                       :: DataType = Real8_                
        integer                                       :: Dim
        integer,    dimension(:),       allocatable   :: CountDim
        integer                                       :: x = -99, y = -99, z = -99, t = -99
        real(4),    dimension(:),       allocatable   :: R41D
        real(4),    dimension(:,:  ),   allocatable   :: R42D        
        real(4),    dimension(:,:,:),   allocatable   :: R43D
        real(4),    dimension(:,:,:,:), allocatable   :: R44D        
        real(8),    dimension(:),       allocatable   :: R81D
        real(8),    dimension(:,:  ),   allocatable   :: R82D        
        real(8),    dimension(:,:,:),   allocatable   :: R83D
        real(8),    dimension(:,:,:,:), allocatable   :: R84D        
        integer(4), dimension(:),       allocatable   :: I41D
        integer(4), dimension(:,:  ),   allocatable   :: I42D        
        integer(4), dimension(:,:,:),   allocatable   :: I43D
        integer(4), dimension(:,:,:,:), allocatable   :: I44D        
    end type T_ValueIn    
    
    private :: T_Date
    type       T_Date
        character(len=StringLength)             :: NetCDFName    
        type (T_ValueIn)                        :: ValueIn        
        type(T_Time), dimension(:), pointer     :: Value1DOut
        integer                                 :: NumberInst
        type(T_Time)                            :: RefDateTimeIn, RefDateTimeOut            
        real                                    :: UnitsFactor
        logical                                 :: RefAttribute
        character(len=StringLength)             :: RefAttributeName
        integer                                 :: NetCDFvar, NetCDFdim
    end type  T_Date

    private :: T_LongLat
    type       T_LongLat
        character(len=StringLength)             :: NetCDFNameLat
        character(len=StringLength)             :: NetCDFNameLong        
        type (T_ValueIn)                        :: LongIn,  LatIn        
        real, dimension(:,:),     pointer       :: LongOut, LatOut
        integer                                 :: imax, jmax
    end type  T_LongLat


    private :: T_Depth
    type       T_Depth
        character(len=StringLength)             :: NetCDFName
        logical                                 :: Dim3D = .false.
        type (T_ValueIn)                        :: ValueIn        
        real, dimension(:,:,:),   pointer       :: Value3DOut
        logical                                 :: ON
        logical                                 :: Sigma !true - sigma false - z-level                
        logical                                 :: SigmaDistortion            = .false.
        real                                    :: theta_s, theta_b, Hc
        integer                                 :: kmax
    end type  T_Depth

    private :: T_Bathym
    type       T_Bathym
        character(len=StringLength)             :: NetCDFName
        logical                                 :: ON = .false.
        real                                    :: Default
        type (T_ValueIn)                        :: ValueIn        
        real, dimension(:,:),     pointer       :: Value2DOut
    end type  T_Bathym

    private :: T_Mapping
    type       T_Mapping
        character(len=StringLength)             :: NetCDFName
        logical                                 :: ON = .true.
        real                                    :: Limit
        type (T_ValueIn)                        :: ValueIn        
        integer, dimension(:,:,:),   pointer    :: Value3DOut
        integer, dimension(:,:  ),   pointer    :: Value2DOut
    end type  T_Mapping    

    private :: T_Field
    type       T_Field
        type(T_PropertyID)                      :: ID
        integer                                 :: Dim
        character(len=StringLength)             :: NetCDFName
        real                                    :: Add, Multiply
        type (T_ValueIn)                        :: ValueIn        
        real, dimension(:,:),     pointer       :: Value2DOut
        real, dimension(:,:,:),   pointer       :: Value3DOut
    end type  T_Field

    
    private :: T_NetCDFCF_2_HDF5MOHID
    type       T_NetCDFCF_2_HDF5MOHID
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjGridData          = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjGeometry          = 0
        integer                                 :: ObjMap               = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: ObjNetCDF            = 0
        integer                                 :: Unit, ClientNumber
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: GeometryFileName
        character(len=PathLength)               :: OutputNetCDFFile
        type(T_Size3D)                          :: Size, WorkSize
        type(T_Size2D)                          :: Size2D, WorkSize2D
        type(T_Field),  dimension(:), allocatable :: Field 
        integer                                 :: PropNumber
        type(T_Date)                            :: Date  
        type(T_LongLat)                         :: LongLat
        type(T_Depth)                           :: Depth
        type(T_Mapping)                         :: Mapping
        type(T_Bathym)                          :: Bathym
        logical                                 :: OutHDF5, OutNetcdf
    end type  T_NetCDFCF_2_HDF5MOHID

    type(T_NetCDFCF_2_HDF5MOHID), pointer             :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertNetCDFCF_2_HDF5MOHID(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadOptions


        call ReadNetCDFCF_WriteHDF5MOHID

        !call WriteVelocityModulus

        !call KillNetCDFCF_2_HDF5MOHID

        STAT = SUCCESS_

    end subroutine ConvertNetCDFCF_2_HDF5MOHID

    !------------------------------------------------------------------------

    subroutine WriteVectorModulus(Field_U, Field_V, VectorModulus_, iOut, PropUnits)
    
        !Arguments-------------------------------------------------------------
        type(T_Field)                               :: Field_U, Field_V
        integer                                     :: VectorModulus_, iOut        
        character(Len=StringLength)                 :: PropUnits 
        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:  ), pointer         :: Aux3D
        real,   dimension(:,:    ), pointer         :: Aux2D        
        integer                                     :: STAT_CALL
        character(Len=StringLength)                 :: MohidName
        integer                                     :: i, j, k
        !Begin-----------------------------------------------------------------

        MohidName = GetPropertyName(VectorModulus_)

        if (Field_U%Dim==3) then
        
            allocate(Aux3D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))

            Aux3D(:,:,:) = 0.0

            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if(Me%Mapping%Value3DOut(i,j,k) == 1)then
                    Aux3D(i,j,k) = sqrt(Field_U%Value3DOut(i,j,k)**2 + Field_V%Value3DOut(i,j,k)**2)
                endif

            enddo
            enddo
            enddo

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,          &
                                             Me%WorkSize%JLB, Me%WorkSize%JUB,          &
                                             Me%WorkSize%KLB, Me%WorkSize%KUB,          &
                                             STAT = STAT_CALL)

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array3D = Aux3D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR40'
        
        endif
        
        if (Field_U%Dim==2) then
        
            allocate(Aux2D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

            Aux2D(:,:) = 0.0

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if(Me%Mapping%Value2DOut(i,j) == 1)then
                    Aux2D(i,j) = sqrt(Field_U%Value2DOut(i,j)**2 + Field_V%Value2DOut(i,j)**2)
                endif

            enddo
            enddo

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,          &
                                 Me%WorkSize%JLB, Me%WorkSize%JUB,                      &
                                 STAT = STAT_CALL)

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(MohidName),              &
                                 trim(MohidName),trim(PropUnits), Array2D = Aux2D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteVelocityModulus - ModuleMecatorFormat - ERR40'
        
        endif        

        if (Field_U%Dim==3) then
            deallocate(Aux3D)
            nullify   (Aux3D)
        endif
        
        if (Field_U%Dim==2) then
            deallocate(Aux2D)
            nullify   (Aux2D)
        endif
        

    end subroutine WriteVectorModulus
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
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


        call GetData(Me%GeometryFilename,                                           &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'OUTPUT_GEOMETRY_FILENAME',                     &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
        

        call GetData(Me%PropNumber,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PROPERTIES_NUMBER',                                &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
        if (iflag == 0)            stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
        
        call GetData(Me%OutHDF5,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'HDF5_OUT',                                         &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
        
        call GetData(Me%OutNetCDF,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NETCDF_OUT',                                       &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
        
        if (Me%OutNetCDF) then
        
            call GetData(Me%OutputNetCDFFile,                                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'OUTPUT_NETCDF_FILE',                           &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR120'
        
        endif

        allocate(Me%Field(1:Me%PropNumber))

        call ReadTimeOptions
        
        call ReadGridOptions
        
        call ReadFieldOptions


    end subroutine ReadOptions

    !------------------------------------------------------------------------
    
    subroutine ReadTimeOptions

        !Local-----------------------------------------------------------------
        real, dimension(6)                          :: AuxTime
        logical                                     :: BlockFound
        integer                                     :: iflag, STAT_CALL
        !begin-----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   begin_time, end_time,                                &
                                   BlockInBlockFound = BlockFound,  STAT = STAT_CALL)


IS:     if(STAT_CALL .EQ. SUCCESS_) then
            
BF:         if (BlockFound) then


                call GetData(Me%Date%NetCDFName,                                        &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME',                              &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "time",                                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

               !The time in netcdf is assumed in days from 1950/1/1 : 0h:0m:0s (Mercator Ocean convention)
                AuxTime(:) = 0.
                AuxTime(1) = 1950
                AuxTime(2) = 1
                AuxTime(3) = 1

                call GetData(AuxTime,                                                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'REFERENCE_DATE_IN',                        &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'    
    

                call SetDate (Me%Date%RefDateTimeIn, Year    = AuxTime(1),              &
                                                     Month   = AuxTime(2),              &
                                                     Day     = AuxTime(3),              &
                                                     Hour    = AuxTime(4),              &
                                                     Minute  = AuxTime(5),              &
                                                     Second  = AuxTime(6))

                call StartComputeTime(Me%ObjTime, Me%Date%RefDateTimeIn, Me%Date%RefDateTimeIn, &
                                      Me%Date%RefDateTimeIn, DT = 0.0,                        &      
                                      VariableDT = .false., STAT = STAT_CALL)   
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
                                     
                call GetData(AuxTime,                                                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'REFERENCE_DATE_OUT',                       &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'  
                
                if (iflag==0) then
                   !The time in netcdf is assumed in days from 1950/1/1 : 0h:0m:0s (Mercator Ocean convention)
                    AuxTime(:) = 0.
                    AuxTime(1) = 1950
                    AuxTime(2) = 1
                    AuxTime(3) = 1                
                endif  

                call SetDate (Me%Date%RefDateTimeOut, Year    = AuxTime(1),             &
                                                      Month   = AuxTime(2),             &
                                                      Day     = AuxTime(3),             &
                                                      Hour    = AuxTime(4),             &
                                                      Minute  = AuxTime(5),             &
                                                      Second  = AuxTime(6))                                      

                call GetData(Me%Date%ValueIn%DataType,                                  &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'DATA_TYPE_IN',                             &
                             default      = Real8_,                                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'    
                
                Me%Date%ValueIn%Dim = 1

                !call GetData(Me%Date%UnitsFactor,                                       &
                !             Me%ObjEnterData, iflag,                                    &
                !             SearchType   = FromBlock,                                  &
                !             keyword      = 'UNITS_FACTOR',                             &
                !             default      = 1.,                                         &
                !             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                !             STAT         = STAT_CALL)        
                !if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'    

                call GetData(Me%Date%RefAttribute,                                      &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'REF_DATE_ATTRIBUTE',                       &
                             default      = .true.,                                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'    
                
                if (Me%Date%RefAttribute) then
        
                    call GetData(Me%Date%RefAttributeName,                                  &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'REF_DATE_ATTRIBUTE_NAME',                  &
                                 !CF convention
                                 default      = "units",                                    &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'    
                    
                endif
                
                !Index 1 is time
                Me%Date%ValueIn%t = 1
                
                !call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                !if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'    
                
 
        
            else BF
            
                stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'    
                
            endif BF
        else IS
        
            stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR200'    
        
        endif IS                             
                              
                              
    end subroutine ReadTimeOptions                              


    !------------------------------------------------------------------------
    
    subroutine ReadGridOptions

        !Local-----------------------------------------------------------------
        logical                                     :: BlockFound
        integer                                     :: iflag, STAT_CALL
        !begin-----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   begin_grid , end_grid ,                              &
                                   BlockInBlockFound = BlockFound,  STAT = STAT_CALL)


IS:     if(STAT_CALL .EQ. SUCCESS_) then
            
BF:         if (BlockFound) then


                call GetData(Me%LongLat%NetCDFNameLat,                                  &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_LAT',                          &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "latitude",                                 &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                call GetData(Me%LongLat%NetCDFNameLong,                                 &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_LONG',                         &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "longitude",                                &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                
   
                
                Me%LongLat%LatIn%DataType  = Real8_
                Me%LongLat%LongIn%DataType = Real8_
                
                Me%LongLat%LatIn%Dim  = 2
                Me%LongLat%LongIn%Dim = 2

                call GetData(Me%Depth%NetCDFName,                                       &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_DEPTH',                        &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
                if  (iflag == 0) then 
                    Me%Depth%Dim3D = .false.
                    Me%Depth%kmax = -99
                else 
                    Me%Depth%Dim3D = .true.
                endif
                
                Me%Depth%ValueIn%DataType = Real8_
                Me%Depth%ValueIn%Dim      = 1
                
                call GetData(Me%Bathym%NetCDFName,                                      &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_BATHYM',                       &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                
                if  (iflag == 0) then 
                    Me%Bathym%ON = .false.
                else 
                    Me%Bathym%ON = .true.
                endif  
                
                Me%Bathym%ValueIn%DataType = Real8_
                Me%Bathym%ValueIn%Dim      = 2
                
                if (.not. Me%Bathym%ON) then
                
                    call GetData(Me%Bathym%Default,                                         &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlockInBlock,                           &
                                 keyword      = 'BATHYM_DEFAULT',                           &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                                 default      = 0.,                                         &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                
                endif
                
         
                
                call GetData(Me%Mapping%NetCDFName,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_MAPPING',                      &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'   
                
                if (iflag == 0) then
                    Me%Mapping%ON = .false.
                else
                    Me%Mapping%ON = .true.
                endif
                
                call GetData(Me%Mapping%Limit,                                          &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'MAPPING_LIMIT',                            &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = 0.5,                                        &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
                
                Me%Mapping%ValueIn%DataType = Real8_
 
                
                !call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                !if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'    
                
            else BF
            
                stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'    
                
            endif BF
        else IS
        
            stop 'ReadGrid2DOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR200'    
        
        endif IS                             
                              
                              
    end subroutine ReadGridOptions                              

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------
    
    subroutine ReadFieldOptions

        !Local-----------------------------------------------------------------
        logical                                     :: BlockFound
        integer                                     :: iflag, STAT_CALL, ip
        !begin-----------------------------------------------------------------

        do ip = 1, Me%PropNumber 
                        
            call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                &
                                       begin_field, end_field,                          &
                                       BlockInBlockFound = BlockFound,  STAT = STAT_CALL)


    IS:     if(STAT_CALL .EQ. SUCCESS_) then
                    
    BF:         if (BlockFound) then

                    call GetData(Me%Field(ip)%NetCDFName,                               &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'NETCDF_NAME',                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
                    if (iflag==0)              stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

                    call GetData(Me%Field(ip)%Add,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ADD_FACTOR',                           &
                                 default      = 0.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'    
                    
                    call GetData(Me%Field(ip)%Multiply,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'MULTIPLY_FACTOR',                      &
                                 default      = 1.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'    

                    call GetData(Me%Field(ip)%Dim,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DIM',                                  &
                                 default      = 2,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'    
                    
                    if (Me%Field(ip)%Dim /= 2 .and. Me%Field(ip)%Dim /= 3)              &
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'    

                    call ConstructPropertyID (Me%Field(ip)%ID,  Me%ObjEnterData, FromBlockInBlock)
                    
                    !Index 1 is time
                    Me%Field(ip)%ValueIn%DataType = Real8_
                    
                    if      (Me%Field(ip)%Dim == 2) then
                        Me%Field(ip)%ValueIn%Dim = 3
                    else if (Me%Field(ip)%Dim == 3) then
                        Me%Field(ip)%ValueIn%Dim = 4
                    endif
                    
                    !call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                    !if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
                    
     
            
                else BF
                
                    stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'    
                    
                endif BF
            else IS
            
                stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR200'    
            
            endif IS            
        
        enddo                 
                              
                              
    end subroutine ReadFieldOptions                              


    !------------------------------------------------------------------------

    subroutine ReadNetCDFCF_WriteHDF5MOHID

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: BlockFound
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL, HDF5_CREATE, iOut, status, ncid


        !Begin----------------------------------------------------------------
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        if (Me%OutHDF5) then
            !Opens HDF5 File
            call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Me%OutNetCDF) then
            !Opens Netcdf File
            STAT_CALL = nf90_create(path = trim(Me%OutputNetCDFFile), cmode = nf90_clobber, ncid = Me%ObjNetCDF)
            if (STAT_CALL /= SUCCESS_)stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif

        
        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   input_files_begin, input_files_end,                  &
                                   BlockInBlockFound = BlockFound,                      &
                                   FirstLine = FirstLine, LastLine = LastLine,          &
                                   STAT = STAT_CALL)

IS:     if(STAT_CALL .EQ. SUCCESS_) then

            !The block is found to exist before when reading depth
BF:         if (BlockFound) then            
            
                iOut = 0

                do line = FirstLine + 1, LastLine - 1

                    call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                                 Buffer_Line = line, STAT = STAT_CALL)

                    
                    !Verifies if file exists
                    status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
                    if (status /= nf90_noerr) stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                    
                    write(*,*) 'Reading ', trim(InputFile)

                    !Time
                    call ReadWriteTime          (ncid, iOut)                    
                        
                    if (line == FirstLine + 1) then
                        !Grid/Bathym and Grid/Longitude and Grid/Latitude
                         call ReadWriteGrid2D  (ncid) 
                         !Bathym, mapping
                         call ReadWriteGrid3D  (ncid)
                    endif 
                    !Grid/Vertical Z
                    if  (Me%Depth%On) then
    !                        call ReadWriteVerticalZ    (ncid, iOut)
                    endif
                    !Results/XXX    
                    call ReadWriteFields    (ncid, iOut) 
                    
                    iOut = iOut + Me%Date%NumberInst
                    
                    status=NF90_CLOSE(ncid)
                    if (status /= nf90_noerr)  stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

                
                enddo
                
            else  BF
                stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            endif BF

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        else   IS

            stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR60'

        end if IS


    end subroutine ReadNetCDFCF_WriteHDF5MOHID

    !------------------------------------------------------------------------


    !------------------------------------------------------------------------

    
  !  subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
 !       type (T_Field), pointer                   :: FirstField
 !       type (T_Field), pointer                   :: ObjField

        !Local-----------------------------------------------------------------
  !      type (T_Field), pointer                   :: NewField
  !      type (T_Field), pointer                   :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
  !      allocate (NewField)
  !      nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
  !      if (.not. associated(FirstField)) then
  !          FirstField         => NewField
  !          ObjField           => NewField
  !      else
  !          PreviousField      => FirstField
  !          ObjField           => FirstField%Next
  !          do while (associated(ObjField))
  !              PreviousField  => ObjField
  !              ObjField       => ObjField%Next
  !          enddo
  !          ObjField           => NewField
  !          PreviousField%Next => NewField
  !      endif


  !  end subroutine AddField
    
    
    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine ReadWriteTime(ncid, iOut)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut, ncid
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        call ReadTimeNetCDF(ncid)
        
        if (Me%OutHDF5  ) call WriteTimeHDF5  (iOut)
        
        if (Me%OutNetCDF) call WriteTimeNetCDF(iOut)
        
        call DeAllocateValueIn(Me%Date%ValueIn)            

    end subroutine ReadWriteTime

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine ReadWriteGrid2D(ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: ncid
        
        !Local-----------------------------------------------------------------
        real,   dimension(:),   pointer                 :: Dummy
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call ReadGrid2DNetCDF(ncid)
        
        if (Me%OutHDF5) then
        
            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%LongLat%LatOut, Me%LongLat%LongOut, &
                                         XX  = Dummy, YY = Dummy, Latitude = 45., Longitude = -8.,    &
                                         ILB = Me%WorkSize%ILB, IUB = Me%WorkSize%IUB,                &
                                         JLB = Me%WorkSize%JLB, JUB = Me%WorkSize%JUB,                &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadWriteGrid2D - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
            
            Me%WorkSize2D%ILB = Me%WorkSize%ILB
            Me%WorkSize2D%IUB = Me%WorkSize%IUB        
            Me%WorkSize2D%JLB = Me%WorkSize%JLB
            Me%WorkSize2D%JUB = Me%WorkSize%JUB        
                                         

            call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,                         &
                                     WorkSize = Me%WorkSize2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadWriteGrid2D - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        

    end subroutine ReadWriteGrid2D

    !------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine ReadWriteGrid3D(ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: ncid
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        call ReadGrid3DNetCDF(ncid)
        
        if (Me%OutHDF5) call WriteGrid3DHDF5

    end subroutine ReadWriteGrid3D

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine ReadWriteFields(ncid, iOut)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut, ncid
        
        !Local-----------------------------------------------------------------
        integer                                         :: iP , iFinal, iT, i, j, k, mask
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            call ReadFieldNetCDF(ncid, iP)
            
            if (Me%OutHDF5) then

                do iT = 1, Me%Date%NumberInst
                
                    if      (Me%Field(iP)%Dim==3 .and. .not. associated(Me%Field(iP)%Value3DOut)) then 
                        allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB,           &
                                                         Me%Size%KLB:Me%Size%KUB))
                    else if (Me%Field(iP)%Dim==2 .and. .not. associated(Me%Field(iP)%Value2DOut)) then 
                        allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB))
                    endif                            
                    
                    
                    if       (Me%Field(iP)%Dim==3) then
                    
                        do k= Me%WorkSize%KLB, Me%WorkSize%KUB                    
                        do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                        
                            if (Me%Mapping%Value3DOut(i,j,k) == 1) then
                        
                                Me%Field(iP)%Value3DOut(i, j, k) = GetNetCDFValue(Me%Field(iP)%ValueIn,  Dim1 = j+1, &
                                                                                  Dim2 = i+1, Dim3 = k, Dim4 = iT)
                                Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%Value3DOut(i, j, k) * &
                                                                   Me%Field(iP)%Multiply + Me%Field(iP)%Add

                            else 
                                Me%Field(iP)%Value3DOut(i, j, k) = FillValueReal
                            endif                
                        enddo
                        enddo
                        enddo
                        
                    else if  (Me%Field(iP)%Dim==2) then
                    
                        do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                        
                            if (Me%Depth%Dim3D) then
                                mask = Me%Mapping%Value3DOut(i,j,Me%WorkSize%KUB)
                            else
                                mask = Me%Mapping%Value2DOut(i,j)
                            endif
                            
                            if (mask == 1) then
                                Me%Field(iP)%Value2DOut(i, j) = GetNetCDFValue(Me%Field(iP)%ValueIn,  &
                                                                    Dim1 = j+1, Dim2 = i+1, Dim3 = iT)
                                Me%Field(iP)%Value2DOut(i, j) = Me%Field(iP)%Value2DOut(i, j) * &
                                                                Me%Field(iP)%Multiply + Me%Field(iP)%Add                            
                            else 
                                Me%Field(iP)%Value2DOut(i, j) = FillValueReal
                            endif                
                            
                        enddo
                        enddo
                        
                    endif

                    iFinal = iT + iOut
                    
                    call WriteFieldHDF5(iP, iFinal)    
        
                enddo
                
            endif
        
            call DeAllocateValueIn(Me%Field(iP)%ValueIn)

        enddo  

    end subroutine ReadWriteFields        
                  
   !------------------------------------------------------------------------
    subroutine WriteTimeHDF5(iOut)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut
        
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real(8)                                         :: Aux        
        real,    dimension(:), pointer                  :: TimePtr
        type(T_Time)                                    :: CurrentTime
        integer                                         :: STAT_CALL, i

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Time HDF5 file...'

        !allocate(Me%Date%Value1DOut(1:Me%Date%NumberInst))
        

        do i=1, Me%Date%NumberInst

            Aux = GetNetCDFValue(Me%Date%ValueIn, Dim1 = i)
            
            Aux = Aux * dble(Me%Date%UnitsFactor)
            
            CurrentTime = Me%Date%RefDateTimeIn + Aux
       

!           Dados para escriver uma soa vez cada date:
            call ExtractDate   (CurrentTime,                                            &
                                AuxTime(1), AuxTime(2), AuxTime(3),                     &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                                   &
                                 "Time", "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D = TimePtr,                                     &
                                 OutputNumber = i+iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
           
        enddo
        
       

    end subroutine WriteTimeHDF5


   !---------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine WriteTimeNetCDF(iOut)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut
        
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real(8)                                         :: Aux, SetValue        
        type(T_Time)                                    :: CurrentTime
        integer                                         :: STAT_CALL, i, iaux
        character(len=Stringlength)                     :: AuxChar

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Time NetCDF file...'

        !allocate(Me%Date%Value1DOut(1:Me%Date%NumberInst))

        if (iOut==0) then
            
            STAT_CALL = nf90_def_dim(ncid=Me%ObjNetCDF, name="time", len= nf90_unlimited, dimid= Me%Date%NetCDFdim)
            STAT_CALL = nf90_def_var(ncid=Me%ObjNetCDF, name="time", xtype=nf90_double, dimids= Me%Date%NetCDFdim, varid=Me%Date%NetCDFvar)
            if (STAT_CALL /= SUCCESS_)then
                write(*,*) trim(nf90_strerror(STAT_CALL))
                stop 'WriteTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
            endif

            call ExtractDate (Me%Date%RefDateTimeOut, Year    = AuxTime(1),             &
                                                      Month   = AuxTime(2),             &
                                                      Day     = AuxTime(3),             &
                                                      Hour    = AuxTime(4),             &
                                                      Minute  = AuxTime(5),             &
                                                      Second  = AuxTime(6))              
            
            AuxChar(1:13)="second since "
            write(AuxChar(14:17),'(I4)') int(AuxTime(1))
            write(AuxChar(18:18),'(A1)') "-"
            write(AuxChar(21:21),'(A1)') "-"            
            write(AuxChar(24:24),'(A1)') " "                        
            write(AuxChar(27:27),'(A1)') ":"                                    
            write(AuxChar(31:31),'(A1)') ":"   
            do i=2,6                    
                iaux = (i-2)*3                         
                if (AuxTime(i)>=10) then
                    write(AuxChar(19+iaux:20+iaux),'(I2)') int(AuxTime(i))                    
                else
                    write(AuxChar(19+iaux:19+iaux),'(A1)') "0"
                    write(AuxChar(20+iaux:20+iaux),'(I1)') int(AuxTime(i))                    
                endif
            enddo

            STAT_CALL = nf90_put_att(Me%ObjNetCDF, Me%Date%NetCDFvar, "units",trim(AuxChar))
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
            
            STAT_CALL = nf90_enddef(ncid=Me%ObjNetCDF)
            if (STAT_CALL /= SUCCESS_)then
                write(*,*) trim(nf90_strerror(STAT_CALL))
                stop 'WriteTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            endif
                    
        endif        

        do i=1, Me%Date%NumberInst

            Aux = GetNetCDFValue(Me%Date%ValueIn, Dim1 = i)
            
            Aux = Aux * dble(Me%Date%UnitsFactor)
            
            CurrentTime = Me%Date%RefDateTimeIn + Aux
            
            SetValue = CurrentTime - Me%Date%RefDateTimeOut
            
            call SetNetCDFValue(Me%Date%ValueIn, SetValue, Dim1 = i)
       
        enddo
        
 
        call PutNetCDFMatrix(Me%ObjNetCDF, Me%Date%NetCDFvar, Me%Date%ValueIn, StartTime = iOut+1)        
        
    end subroutine WriteTimeNetCDF


   !---------------------------------------------------------------------------
   
   !------------------------------------------------------------------------
    subroutine WriteGrid3DHDF5
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Grid 3D HDF5 file...'

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB, Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


        call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                       &
                             "Bathymetry", "m",                                         &
                             Array2D = Me%Bathym%Value2DOut,                            &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
        
        if (Me%Depth%Dim3D) then

            call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,          &
                                             Me%WorkSize%JLB, Me%WorkSize%JUB,          &
                                             Me%WorkSize%KLB, Me%WorkSize%KUB,          &
                                             STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'


            call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                   &
                                 "WaterPoints3D", "-",                                  &
                                 Array3D = Me%Mapping%Value3DOut,                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            
       
        else

            call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                   &
                                 "WaterPoints2D", "-",                                  &
                                 Array2D = Me%Mapping%Value2DOut,                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        
           
        endif         
        

    end subroutine WriteGrid3DHDF5


   !---------------------------------------------------------------------------
   
   
    subroutine ReadTimeNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        character(Len=StringLength)             :: ref_date
        real, dimension(6)                      :: AuxTime
        integer                                 :: n, status, dimid, i, tmax
        logical                                 :: ReadTime
        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Time NetCDF file...'


        status=NF90_INQ_DIMID(ncid,trim(Me%Date%NetCDFName),dimid)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%Date%NumberInst)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        call AllocateValueIn(Me%Date%ValueIn, Dim1 = Me%Date%NumberInst)

        status = nf90_inq_varid(ncid, trim(Me%Date%NetCDFName), n)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

        call GetNetCDFMatrix(ncid, n, Me%Date%ValueIn) 
        
        if (Me%Date%RefAttribute) then
        
            status=NF90_GET_ATT(ncid,n,trim(Me%Date%RefAttributeName), ref_date)
            if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            
            
            
            tmax = len_trim(ref_date)

            ReadTime =.false.
            
            Me%Date%UnitsFactor = 3600.

            do i=1,tmax-5
                if (ref_date(i:i+5)== "second") then
                    ReadTime =.true.
                    Me%Date%UnitsFactor = 1.
                    exit
                endif
            enddo

            do i=1,tmax-4            
                if (ref_date(i:i+4)== "since") then
                    ref_date = ref_date(i+5:tmax)
                    exit
                endif
                
            enddo

            do i=1,len_trim(ref_date)

                if (ref_date(i:i) =='_'.or.ref_date(i:i) ==':'.or. ref_date(i:i) =='-'&
                    .or. ref_date(i:i) =='Z'.or. ref_date(i:i) =='T') then
                    ref_date(i:i) = ' '
                endif

            enddo  
            
            ref_date(1:19) = trim(adjustl(ref_date))
            
            AuxTime(:) = 0.
            
            read(ref_date(1:4),*) AuxTime(1)    
            read(ref_date(6:7),*) AuxTime(2)                
            read(ref_date(9:10),*) AuxTime(3)
            if (ReadTime) then                            
                read(ref_date(12:13),*) AuxTime(4) 
                read(ref_date(15:16),*) AuxTime(5)                                                                   
                read(ref_date(18:19),*) AuxTime(6)                                                                               
            endif

                        
            call SetDate (Me%Date%RefDateTimeIn, Year    = AuxTime(1),                &
                                               Month   = AuxTime(2),                &
                                               Day     = AuxTime(3),                &
                                               Hour    = AuxTime(4),                &
                                               Minute  = AuxTime(5),                &
                                               Second  = AuxTime(6))
               
            
        endif
        
 
    end subroutine ReadTimeNetCDF

    !------------------------------------------------------------------------
   !---------------------------------------------------------------------------
    subroutine ReadGrid2DNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        integer                                 :: status, numDims, i, j
        integer                                 :: RhVarIdLat, RhVarIdLong
        integer, dimension(nf90_max_var_dims)   :: rhDimIdsLat, rhDimIdsLong
        real(8), dimension(:), allocatable      :: Long1D, Lat1D
        real(8)                                 :: X1, X2, X3, X4, Y1, Y2, Y3, Y4
        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Grid NetCDF file...'


        status=nf90_inq_varid(ncid,trim(Me%LongLat%NetCDFNameLong),RhVarIdLong)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        status = nf90_inquire_variable(ncid, RhVarIdLong, ndims = numDims)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        status = nf90_inquire_variable(ncid, RhVarIdLong, dimids = rhDimIdsLong(:numDims))
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        
        status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(1), len = Me%LongLat%jmax)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'        
        
        status=nf90_inq_varid(ncid,trim(Me%LongLat%NetCDFNameLat),RhVarIdLat)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        
        if      (numDims == 2) then     

            status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(2), len = Me%LongLat%imax)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50' 
                    
        else if (numDims == 1) then
            status = nf90_inquire_variable(ncid, RhVarIdLat, dimids = rhDimIdsLat(:numDims))
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'      
            
            status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLat(1), len = Me%LongLat%imax)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR70' 
            
        else
            stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
        endif
        
        call AllocateValueIn(Me%LongLat%LongIn, Dim1 =  Me%LongLat%jmax, Dim2 = Me%LongLat%imax)
        call AllocateValueIn(Me%LongLat%LatIn,  Dim1 =  Me%LongLat%jmax, Dim2 = Me%LongLat%imax)        
        
        if      (numDims == 2) then    
            call GetNetCDFMatrix(ncid, RhVarIdLong, Me%LongLat%LongIn) 
            call GetNetCDFMatrix(ncid, RhVarIdLat,  Me%LongLat%LatIn) 
        else if (numDims == 1) then
            allocate(Long1D(1:Me%LongLat%jmax))
            allocate(Lat1D (1:Me%LongLat%imax))
            
            status = nf90_get_var(ncid, RhVarIdLong, Long1D)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR90'      
            
            status = nf90_get_var(ncid, RhVarIdLat, Lat1D)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR100'     
             
            do j=1, Me%LongLat%jmax
            do i=1, Me%LongLat%imax
                call SetNetCDFValue(Me%LongLat%LatIn,  Lat1D (i), Dim1 = j,   Dim2 = i  )        
                call SetNetCDFValue(Me%LongLat%LongIn, Long1D(j), Dim1 = j,   Dim2 = i  )
            enddo
            enddo    
                   
        endif
        
        !Build HDF5 MOHID Grid
        Me%WorkSize%ILB = 1
        Me%WorkSize%IUB = Me%LongLat%imax-2
        Me%WorkSize%JLB = 1
        Me%WorkSize%JUB = Me%LongLat%jmax-2

        Me%Size%ILB     = 0
        Me%Size%IUB     = Me%WorkSize%IUB + 1
        Me%Size%JLB     = 0
        Me%Size%JUB     = Me%WorkSize%JUB + 1
       
        allocate(Me%LongLat%LongOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%LongLat%LatOut (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        
        do j=1, Me%Size%JUB
        do i=1, Me%Size%IUB
        
            X1 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j,   Dim2 = i  )
            X2 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j+1, Dim2 = i  )
            X3 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j,   Dim2 = i+1)
            X4 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j+1, Dim2 = i+1)

            Y1 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j,   Dim2 = i  )
            Y2 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j+1, Dim2 = i  )
            Y3 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j,   Dim2 = i+1)
            Y4 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j+1, Dim2 = i+1)
            
            Me%LongLat%LongOut(i, j) = (X1 + X2 + X3 + X4) / 4.
            Me%LongLat%LatOut (i, j) = (Y1 + Y2 + Y3 + Y4) / 4.
            
        enddo
        enddo
        
        call DeAllocateValueIn(Me%LongLat%LongIn)
        call DeAllocateValueIn(Me%LongLat%LatIn )
        
        if (numDims == 1) then
            deallocate(Long1D)
            deallocate(Lat1D )
        endif
        
    end subroutine ReadGrid2DNetCDF

   !---------------------------------------------------------------------------
    subroutine ReadGrid3DNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        integer                                 :: status, i, j, k, dn, mn, bn, dimid, numDims
        real(8)                                 :: Aux
        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Grid 3D NetCDF file...'

        !Read number of layers and their depth
        if (Me%Depth%Dim3D) then

            status=NF90_INQ_DIMID(ncid,trim(Me%Depth%NetCDFName),dimid)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%Depth%kmax)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
            call AllocateValueIn(Me%Depth%ValueIn, Dim1 = Me%Depth%kmax)

            status = nf90_inq_varid(ncid, trim(Me%Depth%NetCDFName), dn)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'            
            
            call GetNetCDFMatrix(ncid, dn, Me%Depth%ValueIn) 
            
        else
        
            Me%Depth%kmax = -99
        
        endif
        
        Me%WorkSize%KLB = 1
        Me%WorkSize%KUB = Me%Depth%kmax        

        Me%Size%KLB = Me%WorkSize%KLB - 1
        Me%Size%KUB = Me%WorkSize%KUB + 1      

        if (Me%Mapping%ON) then

            status=nf90_inq_varid(ncid,trim(Me%Mapping%NetCDFName),mn)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

            status = nf90_inquire_variable(ncid, mn, ndims = numDims)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
                   
            Me%Mapping%ValueIn%Dim = numDims

            
            if (Me%Depth%Dim3D) then 
                if      (numDims == 3) then
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax,        &
                                                             Dim3 = Me%Depth%kmax)
                else if (numDims == 4) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax,        &
                                                             Dim3 = Me%Depth%kmax,          &
                                                             Dim4 = Me%Date%NumberInst)
                else if (numDims == 2) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax)
                endif
            else
                if      (numDims == 2) then
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax)
                else if (numDims == 3) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax,        &
                                                             Dim3 = Me%Date%NumberInst)
                endif        
            endif     

            call GetNetCDFMatrix(ncid, mn, Me%Mapping%ValueIn)         
        
        endif
        
        if (Me%Depth%Dim3D) then 

            allocate(Me%Mapping%Value3DOut(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
            
            Me%Mapping%Value3DOut(:,:,:) = 1

            if (Me%Mapping%ON) then

                do k= Me%WorkSize%KLB, Me%WorkSize%KUB                    
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    if      (Me%Mapping%ValueIn%Dim == 4) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = k, Dim4 = 1)
                    else if (Me%Mapping%ValueIn%Dim == 3) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = k)
                    else if (Me%Mapping%ValueIn%Dim == 2) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1)
                    endif
                    
                    if (Me%Mapping%Limit <= 1) then
                        if (Aux > Me%Mapping%Limit) then
                            Me%Mapping%Value3DOut(i,j,k) = 1
                        else
                            Me%Mapping%Value3DOut(i,j,k) = 0
                        endif
                    else
                        if (Aux < Me%Mapping%Limit) then
                            Me%Mapping%Value3DOut(i,j,k) = 1
                        else
                            Me%Mapping%Value3DOut(i,j,k) = 0
                        endif
                    endif
                enddo
                enddo
                enddo
            endif
        else
            allocate(Me%Mapping%Value2DOut(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            
            Me%Mapping%Value2DOut(:,:) = 1
            
            if (Me%Mapping%ON) then

                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB

                    if      (Me%Mapping%ValueIn%Dim == 3) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1,   Dim2 = i+1, Dim3 = 1)
                    else if (Me%Mapping%ValueIn%Dim == 2) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1,   Dim2 = i+1)
                    endif

                    if (Aux > Me%Mapping%Limit) then
                        Me%Mapping%Value2DOut(i,j) = 1
                    else
                        Me%Mapping%Value2DOut(i,j) = 0
                    endif
                    
                    
                enddo
                enddo
                
            endif            
        endif        
        
        !Read bathymetry
        
        if (Me%Bathym%ON) then 
            
            call AllocateValueIn(Me%Bathym%ValueIn, Dim1 = Me%LongLat%jmax, Dim2 = Me%LongLat%imax)

            status = nf90_inq_varid(ncid, trim(Me%Bathym%NetCDFName), bn)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'            
            
            call GetNetCDFMatrix(ncid, bn, Me%Bathym%ValueIn) 
            
        endif
        
        
        allocate(Me%Bathym%Value2DOut(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        
        if (Me%Bathym%ON) then 
        
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                Me%Bathym%Value2DOut(i,j) = GetNetCDFValue(Me%Bathym%ValueIn,  Dim1 = j+1, Dim2 = i+1)
            enddo
            enddo
        
        else

            Me%Bathym%Value2DOut(:,:) = Me%Bathym%Default
        
        endif
        
        if (Me%Depth%Dim3D)                                                             &
            call DeAllocateValueIn(Me%Depth%ValueIn)

        if (Me%Mapping%ON)                                                              &
            call DeAllocateValueIn(Me%Mapping%ValueIn)
        
        if (Me%Bathym%ON)                                                               &
            call DeAllocateValueIn(Me%Bathym%ValueIn)
        
        
    end subroutine ReadGrid3DNetCDF

    !------------------------------------------------------------------------
   !---------------------------------------------------------------------------
    subroutine ReadFieldNetCDF(ncid, iP)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid, iP
        
        !Local-----------------------------------------------------------------
        integer                                 :: status, pn 
        !Begin-----------------------------------------------------------------
        
        !Read number of layers and their depth
        if      (Me%Field(iP)%Dim == 3) then

            call AllocateValueIn(Me%Field(iP)%ValueIn, Dim1 = Me%LongLat%jmax,          &
                                                       Dim2 = Me%LongLat%imax,          &
                                                       Dim3 = Me%Depth%kmax,            &
                                                       Dim4 = Me%Date%NumberInst)
        else if (Me%Field(iP)%Dim == 2) then
            call AllocateValueIn(Me%Field(iP)%ValueIn, Dim1 = Me%LongLat%jmax,          &
                                                       Dim2 = Me%LongLat%imax,          &
                                                       Dim3 = Me%Date%NumberInst)
        endif     

        status = nf90_inq_varid(ncid, trim(Me%Field(iP)%NetCDFName), pn)
        if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'            

        call GetNetCDFMatrix(ncid, pn, Me%Field(iP)%ValueIn)         
        

        
    end subroutine ReadFieldNetCDF

    !------------------------------------------------------------------------
    
    subroutine AllocateValueIn(ValueIn, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------        
        type(T_ValueIn)     :: ValueIn
        integer             :: Dim1
        integer, optional   :: Dim2, Dim3, Dim4
        !Local-----------------------------------------------------------------                
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Dim ==2) then
            if (.not.(present(Dim2))) stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        
        if (Dim ==3) then
            if (.not.(present(Dim2).and.present(Dim3))) stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        endif
        
        if (Dim ==4) then
            if (.not.(present(Dim2).and.present(Dim3).and.present(Dim4))) &
                stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        endif
        
        if      (Dim==1) then
            allocate(ValueIn%CountDim(1))
            ValueIn%CountDim(1) = Dim1
        else if (Dim==2) then
            allocate(ValueIn%CountDim(2))
            ValueIn%CountDim(1) = Dim1
            ValueIn%CountDim(2) = Dim2
        else if (Dim==3) then
            allocate(ValueIn%CountDim(3))
            ValueIn%CountDim(1) = Dim1
            ValueIn%CountDim(2) = Dim2
            ValueIn%CountDim(3) = Dim3
        else if (Dim==4) then
            allocate(ValueIn%CountDim(4))
            ValueIn%CountDim(1) = Dim1
            ValueIn%CountDim(2) = Dim2
            ValueIn%CountDim(3) = Dim3        
            ValueIn%CountDim(4) = Dim4
        else
            stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R81D(1:Dim1))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R41D(1:Dim1))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I41D(1:Dim1))                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R82D(1:Dim1,1:Dim2))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R42D(1:Dim1,1:Dim2))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I42D(1:Dim1,1:Dim2))                        
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R83D(1:Dim1,1:Dim2,1:Dim3))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R43D(1:Dim1,1:Dim2,1:Dim3))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I43D(1:Dim1,1:Dim2,1:Dim3))                        
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R84D(1:Dim1,1:Dim2,1:Dim3,1:Dim4))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R44D(1:Dim1,1:Dim2,1:Dim3,1:Dim4))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I44D(1:Dim1,1:Dim2,1:Dim3,1:Dim4))                        
            
            endif    

        else
            stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
        endif
        
    end subroutine AllocateValueIn

    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    
    subroutine DeAllocateValueIn(ValueIn)
        !Arguments-------------------------------------------------------------        
        type(T_ValueIn)     :: ValueIn
        !Local-----------------------------------------------------------------                
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        deallocate(ValueIn%CountDim)

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R81D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R41D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I41D)                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R82D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R42D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I42D)                        
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R83D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R43D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I43D)                        
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R84D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R44D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I44D)                        
            
            endif    

        else
            stop 'DeAllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        
    end subroutine DeAllocateValueIn

    !------------------------------------------------------------------------
    
    subroutine GetNetCDFMatrix(ncid, n, ValueIn)
        !Arguments-------------------------------------------------------------        
        integer             :: ncid, n
        type(T_ValueIn)     :: ValueIn
        !Local-----------------------------------------------------------------                
        integer             :: Dim, DataTypeIn, status
        character(len=StringLength) :: Error
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then

                status = NF90_GET_VAR(ncid,n,ValueIn%R81D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%R41D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%I41D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then

                status = NF90_GET_VAR(ncid,n,ValueIn%R82D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%R42D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%I42D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then

                status = NF90_GET_VAR(ncid,n,ValueIn%R83D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%R43D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%I43D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
            
            endif

        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then

                status = NF90_GET_VAR(ncid,n,ValueIn%R84D)
                Error = nf90_strerror(status)
                write(*,*) Error
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
                
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%R44D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR120'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%I44D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
            
            endif

        else
            stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR140'
        endif
        
    end subroutine GetNetCDFMatrix

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    
    subroutine PutNetCDFMatrix(ncid, n, ValueIn, StartTime)
        !Arguments-------------------------------------------------------------        
        integer             :: ncid, n
        type(T_ValueIn)     :: ValueIn
        integer, optional   :: StartTime
        !Local-----------------------------------------------------------------                
        integer, dimension(:), allocatable :: Start
        integer             :: Dim, DataTypeIn, status
        character(len=StringLength) :: Error
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif


        allocate(Start(1:Dim))
        Start(:)= 1
        if (present(StartTime)) Start(1) = StartTime
        

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then

                status = NF90_PUT_VAR(ncid,n,ValueIn%R81D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) then
                    write(*,*) trim(nf90_strerror(status))
                    stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                endif                    
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%R41D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%I41D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then

                status = NF90_PUT_VAR(ncid,n,ValueIn%R82D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%R42D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%I42D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then

                status = NF90_PUT_VAR(ncid,n,ValueIn%R83D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%R43D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%I43D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
            
            endif

        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then

                status = NF90_PUT_VAR(ncid,n,ValueIn%R84D, start=Start, count=ValueIn%CountDim)
                Error = nf90_strerror(status)
                write(*,*) Error
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
                
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%R44D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR120'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_PUT_VAR(ncid,n,ValueIn%I44D, start=Start, count=ValueIn%CountDim)
                if (status /= nf90_noerr) stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
            
            endif

        else
            stop 'PutNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR140'
        endif
        
        deallocate(Start)


        
    end subroutine PutNetCDFMatrix

    !------------------------------------------------------------------------


    function GetNetCDFValue(ValueIn, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------  
        real                :: GetNetCDFValue      
        type(T_ValueIn)     :: ValueIn
        integer             :: Dim1
        integer, optional   :: Dim2, Dim3, Dim4
        !Local-----------------------------------------------------------------                
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Dim ==2) then
            if (.not.(present(Dim2))) stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        
        if (Dim ==3) then
            if (.not.(present(Dim2).and.present(Dim3))) stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        endif
        
        if (Dim ==4) then
            if (.not.(present(Dim2).and.present(Dim3).and.present(Dim4))) &
                stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        endif

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R81D(Dim1)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R41D(Dim1)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I41D(Dim1)                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R82D(Dim1,Dim2)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R42D(Dim1,Dim2)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I42D(Dim1,Dim2)                        
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R83D(Dim1,Dim2,Dim3)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R43D(Dim1,Dim2,Dim3)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I43D(Dim1,Dim2,Dim3)                        
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R84D(Dim1,Dim2,Dim3,Dim4)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R44D(Dim1,Dim2,Dim3,Dim4)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I44D(Dim1,Dim2,Dim3,Dim4)                        
            
            endif    

        else
            stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        
    end function GetNetCDFValue

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------

    subroutine SetNetCDFValue(ValueIn, SetValue, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------  
        type(T_ValueIn)     :: ValueIn
        real(8)             :: SetValue                
        integer             :: Dim1
        integer, optional   :: Dim2, Dim3, Dim4
        !Local-----------------------------------------------------------------                
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Dim ==2) then
            if (.not.(present(Dim2))) stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        
        if (Dim ==3) then
            if (.not.(present(Dim2).and.present(Dim3))) stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        endif
        
        if (Dim ==4) then
            if (.not.(present(Dim2).and.present(Dim3).and.present(Dim4))) &
                stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        endif

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R81D(Dim1) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
               ValueIn%R41D(Dim1)  = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I41D(Dim1) = SetValue                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R82D(Dim1,Dim2) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
                ValueIn%R42D(Dim1,Dim2) = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I42D(Dim1,Dim2) = SetValue
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R83D(Dim1,Dim2,Dim3) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
                ValueIn%R43D(Dim1,Dim2,Dim3) = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I43D(Dim1,Dim2,Dim3) = SetValue
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R84D(Dim1,Dim2,Dim3,Dim4) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
                ValueIn%R44D(Dim1,Dim2,Dim3,Dim4) = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I44D(Dim1,Dim2,Dim3,Dim4) = SetValue
            
            endif    

        else
            stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        
    end subroutine SetNetCDFValue

    !------------------------------------------------------------------------
        
    subroutine WriteFieldHDF5(iP, iFinal)

        !Arguments-------------------------------------------------------------
        integer                                         :: iP, iFinal

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: WorkILB, WorkJLB, WorkKLB
        integer                                         :: WorkIUB, WorkJUB, WorkKUB

        !Begin-----------------------------------------------------------------
        
        !Bounds



        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 



        if      (Me%Field(iP)%Dim==2) then
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                 &
                                 WorkJUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR10'        

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array2D = Me%Field(iP)%Value2DOut,                     &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR20'

        else if (Me%Field(iP)%Dim==3) then
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                 &
                                 WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR30'        

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array3D = Me%Field(iP)%Value3DOut,                     &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR40'

        else 

            stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR50'

        endif


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR90'


    end subroutine WriteFieldHDF5


    !----------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    
    subroutine KillNetCDFCF_2_HDF5MOHID
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers, ip
        
        !Begin-----------------------------------------------------------------

            call KillMap(Me%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

            call KillGeometry(Me%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            call KillGridData(Me%ObjGridData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

            call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
            call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR60'

            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR70'

            if (associated(Me%Date%Value1DOut    )) deallocate(Me%Date%Value1DOut   )
            
            if (associated(Me%LongLat%LongOut    )) deallocate(Me%LongLat%LongOut   )
            if (associated(Me%LongLat%LatOut     )) deallocate(Me%LongLat%LatOut    )
            
            if (associated(Me%Bathym%Value2DOut  )) deallocate(Me%Bathym%Value2DOut )
            
            if (associated(Me%Mapping%Value2DOut )) deallocate(Me%Mapping%Value2DOut)
            if (associated(Me%Mapping%Value3DOut )) deallocate(Me%Mapping%Value3DOut)

            do ip = 1, Me%PropNumber   
                if (associated(Me%Field(iP)%Value2DOut)) deallocate(Me%Field(iP)%Value2DOut)
                if (associated(Me%Field(iP)%Value3DOut)) deallocate(Me%Field(iP)%Value3DOut)
            enddo
            deallocate(Me%Field)

            deallocate(Me)
            nullify   (Me)

    
    end subroutine KillNetCDFCF_2_HDF5MOHID

    !--------------------------------------------------------------------------

end module ModuleNetCDFCF_2_HDF5MOHID
