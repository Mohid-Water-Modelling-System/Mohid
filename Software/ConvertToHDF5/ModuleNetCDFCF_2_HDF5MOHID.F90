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

    use ModuleNetCDF
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

    type T_NetCDF_Out                                         
        character(len=PathLength)                           :: Name
        integer                                             :: ObjNETCDF    = 0
        character(len=StringLength)                         :: Title
        character(len=StringLength)                         :: Convention
        character(len=StringLength)                         :: Version
        character(len=StringLength)                         :: History
        character(len=StringLength)                         :: Source
        character(len=StringLength)                         :: Institution
        character(len=StringLength)                         :: References
        integer                                             :: iDate
    end type T_NetCDF_Out

    
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
        integer                                 :: Unit, ClientNumber
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: GeometryFileName
        type(T_NetCDF_Out)                      :: NetCDF_Out
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

        call KillNetCDFCF_2_HDF5MOHID

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
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


        call GetData(Me%GeometryFilename,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GEOMETRY_FILENAME',                         &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
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
        
            call ReadOutNetCDFOptions
        
        endif

        allocate(Me%Field(1:Me%PropNumber))

        call ReadTimeOptions
        
        call ReadGridOptions
        
        call ReadFieldOptions


    end subroutine ReadOptions

    !------------------------------------------------------------------------

    subroutine ReadOutNetCDFOptions
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

    
            call GetData(Me%NetCDF_Out%Name,                                            &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'OUTPUT_NETCDF_FILE',                           &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

            call GetData(Me%NetCDF_Out%Title,                                           &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_TITLE',                                 &
                         Default      = 'EASYCO MODEL SOLUTIONS',                       &                         
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            call GetData(Me%NetCDF_Out%Convention,                                      &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_CONVENTION',                            &
                         Default      = 'CF-1.0',                                       &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            call GetData(Me%NetCDF_Out%Version,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_VERSION',                               &
                         Default      = '3.6.1',                                        &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

            call GetData(Me%NetCDF_Out%History,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_HISTORY',                               &
                         Default      = '-',                                            &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
            call GetData(Me%NetCDF_Out%Source,                                          &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_SOURCE',                                &
                         default      = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'

            
            call GetData(Me%NetCDF_Out%Institution,                                     &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_INSTITUTION',                           &
                         Default      = 'Instituto Superior Tecnico',                   &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            call GetData(Me%NetCDF_Out%References,                                      &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_REFERENCES',                            &
                         Default      = 'http://www.mohid.com',                         &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'

            call GetData(Me%NetCDF_Out%iDate,                                           &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_DATE',                                  &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'    
            
    end subroutine ReadOutNetCDFOptions            
    
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
                                                   STAT_CALL, HDF5_CREATE, iOut, status,&
                                                   ncid


        !Begin----------------------------------------------------------------
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        if (Me%OutHDF5) then
            !Opens HDF5 File
            call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Me%OutNetCDF) then
            
            call OpenNCDFFile

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
       
    subroutine OpenNCDFFile

        !Local-----------------------------------------------------------------
        integer                                     :: NCDF_CREATE, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetNCDFFileAccess(NCDF_CREATE = NCDF_CREATE)
        
        call ConstructNETCDF(ObjNCDFID = Me%NetCDF_Out%ObjNetCDF,                       &
                             FileName  = trim(Me%NetCDF_Out%Name),                      &
                             Access    = NCDF_CREATE,                                   &
                             STAT      = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'OpenNCDFFile - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        

        call NETCDFWriteHeader(NCDFID         = Me%NetCDF_Out%ObjNetCDF     ,           &
                               Title          = Me%NetCDF_Out%Title         ,           &
                               Convention     = Me%NetCDF_Out%Convention    ,           &
                               Version        = Me%NetCDF_Out%Version       ,           &
                               History        = Me%NetCDF_Out%History       ,           &
                               iDate          = Me%NetCDF_Out%iDate         ,           &
                               Source         = Me%NetCDF_Out%Source        ,           &
                               Institution    = Me%NetCDF_Out%Institution   ,           &
                               References     = Me%NetCDF_Out%References    ,           &
                               STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenNCDFFile - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        
        write(*,*)
        write(*,*)'Opened ncdf file                : ', trim(Me%NetCDF_Out%Name)

    end subroutine OpenNCDFFile
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

        Me%WorkSize2D%ILB = Me%WorkSize%ILB
        Me%WorkSize2D%IUB = Me%WorkSize%IUB        
        Me%WorkSize2D%JLB = Me%WorkSize%JLB
        Me%WorkSize2D%JUB = Me%WorkSize%JUB         

        if (Me%OutHDF5) then

            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%LongLat%LatOut, Me%LongLat%LongOut, &
                                         XX  = Dummy, YY = Dummy, Latitude = 45., Longitude = -8.,    &
                                         ILB = Me%WorkSize%ILB, IUB = Me%WorkSize%IUB,                &
                                         JLB = Me%WorkSize%JLB, JUB = Me%WorkSize%JUB,                &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadWriteGrid2D - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
      

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
        
        if (Me%OutNetCDF) call WriteGridNetCDF
            
  

    end subroutine ReadWriteGrid3D

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine WriteGridNetCDF
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer    :: lon, lat, lon_stag, lat_stag, bathym
        real(8), dimension(:,:), pointer    :: SphericMercatorX_stag, SphericMercatorY_stag
        integer, dimension(:,:,:), pointer  :: mask3D
        integer, dimension(:,:  ), pointer  :: mask2D        
        character(len=StringLength)         :: NCDFName, LongName, StandardName, Units, Positive
        real                                :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                             :: i, j, k, STAT_CALL
     
        !Begin-----------------------------------------------------------------
        
        
        allocate(lon     (Me%WorkSize%IUB,  Me%WorkSize%JUB  ), lat     (Me%WorkSize%IUB,Me%WorkSize%JUB    ))
        allocate(lon_stag(Me%WorkSize%IUB+1,Me%WorkSize%JUB+1), lat_stag(Me%WorkSize%IUB+1,Me%WorkSize%JUB+1))

        allocate(SphericMercatorX_stag(Me%WorkSize%IUB+1,Me%WorkSize%JUB+1),            &
                 SphericMercatorY_stag(Me%WorkSize%IUB+1,Me%WorkSize%JUB+1))
        
        do j=1,  Me%WorkSize%JUB+1
        do i=1,  Me%WorkSize%IUB+1            
        
            lon_stag(i,j) = Me%LongLat%LongOut(i,j)
            lat_stag(i,j) = Me%LongLat%LatOut (i,j)
        
        enddo
        enddo
        
        do j=1,  Me%WorkSize%JUB
        do i=1,  Me%WorkSize%IUB            
        
            lon(i,j) = (lon_stag(i,j) + lon_stag(i+1,j) + lon_stag(i,j+1) + lon_stag(i+1,j+1))/4.
            lat(i,j) = (lat_stag(i,j) + lat_stag(i+1,j) + lat_stag(i,j+1) + lat_stag(i+1,j+1))/4.
        
        enddo
        enddo

        if (Me%Depth%Dim3D) then
            call NETCDFSetDimensions (Me%NetCDF_Out%ObjNetCDF,                          &
                                                    IUB  = Me%WorkSize%IUB,             &
                                                    JUB  = Me%WorkSize%JUB,             &
                                                    KUB  = Me%WorkSize%KUB,             &
                                                    STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        else
            call NETCDFSetDimensions (Me%NetCDF_Out%ObjNetCDF,                          &
                                                    IUB  = Me%WorkSize%IUB,             &
                                                    JUB  = Me%WorkSize%JUB,             &
                                                    KUB  = 0,                           &
                                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        endif

        call WGS84toGoogleMaps(lon_stag, lat_stag, Me%WorkSize2D,  SphericMercatorX_stag, SphericMercatorY_stag)

                
        call NETCDFWriteLatLon(Me%NetCDF_Out%ObjNetCDF, Lat, Lon, Lat_Stag, Lon_Stag,   &
                               SphericMercatorX_stag, SphericMercatorY_stag,            &
                               GeoCoordinates = .true., STAT = STAT_CALL)
                                
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

        deallocate(lon     , lat     )
        deallocate(lon_stag, lat_stag)
        deallocate(SphericMercatorX_stag, SphericMercatorY_stag)

        
        allocate(bathym(Me%WorkSize%IUB,Me%WorkSize%JUB))
        
        do j=1,  Me%WorkSize%JUB
        do i=1,  Me%WorkSize%IUB
            bathym(i,j) = Me%Bathym%Value2DOut(i,j)
        enddo
        enddo        
        
        call BuildAttributes("Bathymetry", NCDFName, LongName, StandardName,            &
                                           Units, ValidMin, ValidMax,                   &
                                           MinValue, MaxValue, MissingValue, Positive, bathym)
        
        call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,                  &
                              Name          = trim(NCDFName),                           &
                              LongName      = trim(LongName),                           &
                              StandardName  = trim(StandardName),                       & 
                              Units         = trim(Units),                              &
                              ValidMin      = ValidMin,                                 &
                              ValidMax      = ValidMax,                                 &
                              MinValue      = MinValue,                                 &
                              MaxValue      = MaxValue,                                 &
                              MissingValue  = MissingValue,                             &
                              Array2D       = bathym,                                   &
                              STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        deallocate(bathym)
        
        
        if (Me%Depth%Dim3D) then
        
            allocate(mask3D(1:Me%WorkSize%IUB, 1:Me%WorkSize%JUB, 1:Me%WorkSize%KUB))
            
            do k=1,  Me%WorkSize%KUB
            do j=1,  Me%WorkSize%JUB
            do i=1,  Me%WorkSize%IUB            
                mask3D(i,j,k) = Me%Mapping%Value3DOut(i,j,k)
            enddo
            enddo
            enddo
        
            call BuildAttributes("WaterPoints3D", NCDFName,                             &
                                 LongName, StandardName,                                &
                                 Units, ValidMin, ValidMax,                             &
                                 MinValue, MaxValue, MissingValue,                      &
                                 Int3D = mask3D)
    
            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,              &
                                  Name          = trim(NCDFName),                       &
                                  LongName      = trim(LongName),                       &
                                  StandardName  = trim(StandardName),                   & 
                                  Units         = trim(Units),                          &
                                  ValidMin      = ValidMin,                             &
                                  ValidMax      = ValidMax,                             &
                                  MinValue      = MinValue,                             &
                                  MaxValue      = MaxValue,                             &
                                  MissingValue  = MissingValue,                         &
                                  Array3D       = mask3D,                               &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
            
            deallocate(mask3D)

        else


            allocate(mask2D(1:Me%WorkSize%IUB, 1:Me%WorkSize%JUB))
            
            do j=1,  Me%WorkSize%JUB
            do i=1,  Me%WorkSize%IUB            
                mask2D(i,j) = Me%Mapping%Value2DOut(i,j)
            enddo
            enddo



            call BuildAttributes("WaterPoints2D", NCDFName, LongName, StandardName,     &
                                 Units, ValidMin, ValidMax,                             &
                                 MinValue, MaxValue, MissingValue,                      &
                                 Int2D = mask2d)

    
            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,              &
                                  Name          = trim(NCDFName),                       &
                                  LongName      = trim(LongName),                       &
                                  StandardName  = trim(StandardName),                   & 
                                  Units         = trim(Units),                          &
                                  ValidMin      = ValidMin,                             &
                                  ValidMax      = ValidMax,                             &
                                  MinValue      = MinValue,                             &
                                  MaxValue      = MaxValue,                             &
                                  MissingValue  = MissingValue,                         &
                                  Array2D       = mask2D,                               &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            deallocate(mask2D)

        endif
    
    
    end subroutine WriteGridNetCDF

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine ReadWriteFields(ncid, iOut)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut, ncid
        
        !Local-----------------------------------------------------------------
        integer                                         :: iP 
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            call ReadFieldNetCDF(ncid, iP)
            
            call WriteFieldAllInst(iOut, iP)

            call DeAllocateValueIn(Me%Field(iP)%ValueIn)

        enddo
        
          

    end subroutine ReadWriteFields        
                  
   !------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine WriteFieldAllInst(iOut, iP)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut, iP
        
        !Local-----------------------------------------------------------------
        integer                                         :: iFinal, iT, i, j, k, mask
        !Begin-----------------------------------------------------------------
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
                        Me%Field(iP)%Value3DOut(i, j, k) = (Me%Field(iP)%Value3DOut(i, j, k) + Me%Field(iP)%Add) &
                                                            * Me%Field(iP)%Multiply 

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
                        Me%Field(iP)%Value2DOut(i, j) = (Me%Field(iP)%Value2DOut(i, j) + Me%Field(iP)%Add) &
                                                        * Me%Field(iP)%Multiply 
                    else 
                        Me%Field(iP)%Value2DOut(i, j) = FillValueReal
                    endif                
                    
                enddo
                enddo
                
            endif

            iFinal = iT + iOut
            
            if (Me%OutHDF5) then
                call WriteFieldHDF5  (iP, iFinal)    
            endif

            if (Me%OutNetCDF) then
                call WriteFieldNetCDF(iP, iFinal)    
            endif
        enddo    
    
    end subroutine WriteFieldAllInst
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
        real(8), dimension(:), pointer                  :: Times        
        real(8)                                         :: Aux
        type(T_Time)                                    :: CurrentTime
        integer                                         :: STAT_CALL, i
        character(len=Stringlength)                     :: AuxChar

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Time NetCDF file...'
        
        if (iOut == 0) then

            AuxChar = TimeToStringV2(Me%Date%RefDateTimeOut)
           
        endif        
        
        allocate(Times(Me%Date%NumberInst))

        do i=1, Me%Date%NumberInst

            Aux = GetNetCDFValue(Me%Date%ValueIn, Dim1 = i)
            
            Aux = Aux * dble(Me%Date%UnitsFactor)
            
            CurrentTime = Me%Date%RefDateTimeIn + Aux
            
            Aux = CurrentTime - Me%Date%RefDateTimeOut
            
            Times(i) = Aux
            
        enddo
        
        call NETCDFWriteTime(NCDFID       = Me%NetCDF_Out%ObjNetCDF,                    &
                             InitialDate  = AuxChar,                                    &
                             nInstants    = Me%Date%NumberInst,                         &
                             Times        = Times,                                      &
                             StartInstant = iOut+1,                                     &
                             STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        deallocate(Times)

        
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

                        
            call SetDate (Me%Date%RefDateTimeIn, Year    = AuxTime(1),              &
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

    !-------------------------------------------------------------------------
    
    subroutine BuildAttributes(NameIn, NCDFName, LongName, StandardName, Units,           &
                               ValidMin, ValidMax, Min, Max, MissingValue, Positive,    &
                               Float2D, Float3D, Int2D, Int3D)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in )                       :: NameIn
        character(len=*), intent(out)                       :: NCDFName
        character(len=*), intent(out)                       :: LongName
        character(len=*), intent(out)                       :: StandardName
        character(len=*), intent(out)                       :: Units
        character(len=*), intent(out), optional             :: Positive
        real,             intent(out)                       :: ValidMin, Min
        real,             intent(out)                       :: ValidMax, Max
        real,             intent(out)                       :: MissingValue
        real,    dimension(:,:  ), pointer, optional        :: Float2D 
        real,    dimension(:,:,:), pointer, optional        :: Float3D 
        integer, dimension(:,:  ), pointer, optional        :: Int2D   
        integer, dimension(:,:,:), pointer, optional        :: Int3D   

        !Local-----------------------------------------------------------------
        character(len=StringLength)                         :: Name
        integer                                             :: i, j, k
        

        !Begin-----------------------------------------------------------------
        
        call CheckAndCorrectVarName(NameIn, Name)

        select case(trim(adjustl(Name)))

            case("Bathymetry")
                NCDFName        = "bathymetry"
                LongName        = "bathymetry"
                StandardName    = "sea_floor_depth_below_geoid"
                Units           = "m"
                ValidMin        = -50.
                ValidMax        = 11000.
                Positive        = "down"
                MissingValue    = -99.0

            case("WaterPoints2D", "WaterPoints3D", "MappingPoints2D", "WaterPoints")
                NCDFName        = "mask"
                LongName        = "mask of potential water points"
                StandardName    = "land_binary_mask"
                Units           = ""
                ValidMin        = 0
                ValidMax        = 1
            
            case("OpenPoints2D", "OpenPoints3D")
                NCDFName        = "mask"
                LongName        = "mask of effective water points at one given instant"
                StandardName    = "mask"
                Units           = ""
                ValidMin        = 0
                ValidMax        = 1

            case("temperature")
                NCDFName        = "temperature"
                LongName        = "temperature"
                StandardName    = "sea_water_temperature"
                Units           = "degC"
                ValidMin        = 0.
                ValidMax        = 50.
                MissingValue    = FillValueReal

            case("salinity")
                NCDFName        = "salinity"
                LongName        = "salinity"
                StandardName    = "sea_water_salinity"
                Units           = "1e-3"
                ValidMin        = 0.
                ValidMax        = 40.
                MissingValue    = FillValueReal

            case("density")
                NCDFName        = "sea_water_density"
                LongName        = "sea water density"
                StandardName    = "sea_water_density"
                Units           = "kg m-3"
                ValidMin        = 900.
                ValidMax        = 1200.
                MissingValue    = FillValueReal

            case("oxygen")
                NCDFName        = "dissolved_oxygen"
                LongName        = "dissolved oxygen"
                StandardName    = "oxygen"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 30.
                MissingValue    = FillValueReal

            case("dissolved_oxygen_percent_saturation")
                NCDFName        = "dissolved_oxygen_percent_saturation"
                LongName        = "dissolved oxygen percent saturation"
                StandardName    = "dissolved_oxygen_percent_saturation"
                Units           = "%"
                ValidMin        = 0.
                ValidMax        = 200.
                MissingValue    = FillValueReal

            case("velocity_U")
                NCDFName        = "u"
                LongName        = "east-west current velocity"
                StandardName    = "eastward_sea_water_velocity"
                Units           = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = FillValueReal

            case("velocity_V")
                NCDFName        = "v"
                LongName        = "north-south current velocity"
                StandardName    = "northward_sea_water_velocity"
                Units           = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = FillValueReal

            case("velocity_W")
                NCDFName        = "w"
                LongName        = "vertical current velocity"
                StandardName    = "upward_sea_water_velocity"
                Units           = "m s-1"
                ValidMin        = -2.
                ValidMax        = 2.
                MissingValue    = FillValueReal

            case("velocity_modulus")
                NCDFName        = "vm"
                LongName        = "sea water speed"
                StandardName    = "sea_water_speed"
                Units           = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = FillValueReal

            case("water_level")
                NCDFName        = "ssh"
                LongName        = "sea water level"
                StandardName    = "sea_surface_height"
                Units           = "m"
                ValidMin        = -20.
                ValidMax        = 20.
                MissingValue    = FillValueReal

            case("wind_velocity_X")
                NCDFName        = "x_wind"
                LongName        = "east-west wind velocity"
                StandardName    = "east-west_wind_velocity"
                Units           = "m s-1"
                ValidMin        = -100.
                ValidMax        = 100.
                MissingValue    = FillValueReal

            case("wind_velocity_Y")
                NCDFName        = "y_wind"
                LongName        = "north-south wind velocity"
                StandardName    = "north-south_wind_velocity"
                Units           = "m s-1"
                ValidMin        = -100.
                ValidMax        = 100.
                MissingValue    = FillValueReal

            case("air_temperature")
                NCDFName        = "air_temperature"
                LongName        = "air temperature"
                StandardName    = "air_temperature"
                Units           = "degC"
                ValidMin        = -90.
                ValidMax        = 60.
                MissingValue    = FillValueReal

            case("atmospheric_pressure")
                NCDFName        = "air_pressure"
                LongName        = "atmospheric pressure"
                StandardName    = "atmospheric_pressure"
                Units           = "Pa"
                ValidMin        = -90.
                ValidMax        = 60.
                MissingValue    = FillValueReal

            case("short_wave_solar_radiation_extinction")
                NCDFName        = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
                LongName        = "short wave solar radiation light extinction coefficient"
                StandardName    = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
                Units           = "m-1"
                ValidMin        = 0.
                ValidMax        = 100.
                MissingValue    = FillValueReal

            case("short_wave_solar_radiation")
                NCDFName        = "short_wave_solar_radiation"
                LongName        = "short wave solar radiation"
                StandardName    = "omnidirectional_photosynthetic_spherical_irradiance_in_sea_water"
                Units           = "W m-2"
                ValidMin        = 0.
                ValidMax        = 1400.
                MissingValue    = FillValueReal

            case("phytoplankton")
                NCDFName        = "phytoplankton"
                LongName        = "phytoplankton"
                StandardName    = "phytoplankton"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("zooplankton")
                NCDFName        = "zooplankton"
                LongName        = "zooplankton"
                StandardName    = "zooplankton"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("nitrate")
                NCDFName        = "nitrate"
                LongName        = "nitrate"
                StandardName    = "nitrate"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("ammonia")
                NCDFName        = "ammonia"
                LongName        = "ammonia"
                StandardName    = "ammonia"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case default

                NCDFName        = trim(adjustl(Name))
                LongName        = trim(adjustl(Name))
                StandardName    = trim(adjustl(Name))
                Units           = "unknown"
                ValidMin        =   null_real
                ValidMax        = - null_real
                MissingValue    = FillValueReal

        end select

        Min = ValidMax
        Max = ValidMin

if1:   if(present(Int2D) .or. present(Int3D))then
           
            Min = 0
            Max = 1

        elseif(present(Float2D))then if1

            do j = 1, Me%WorkSize%JUB
            do i = 1, Me%WorkSize%IUB

                if(Float2D(i,j) .gt. FillValueReal/2. .and. Float2D(i,j) .lt.  Min .and. &
                                                            Float2D(i,j) .ge. ValidMin)then

                    Min = Float2D(i,j)

                end if

                if(Float2D(i,j) > Max .and. Float2D(i,j) .le. ValidMax) then

                    Max = Float2D(i,j) 

                endif

            enddo
            enddo

        elseif(present(Float3D))then if1

            do j = 1, Me%WorkSize%JUB
            do i = 1, Me%WorkSize%IUB
            do k = 1, Me%WorkSize%KUB

                if(Float3D(i,j,k) .gt. FillValueReal/2. .and. Float3D(i,j,k) .lt.  Min .and. &
                                                              Float3D(i,j,k) .ge. ValidMin)then

                    Min = Float3D(i,j,k)

                end if

                if(Float3D(i,j,k) > Max .and. Float3D(i,j,k) .le. ValidMax) then

                    Max = Float3D(i,j,k) 

                endif

            enddo
            enddo
            enddo

        endif if1

    end subroutine BuildAttributes
    
    
    !--------------------------------------------------------------------------
    
    subroutine CheckAndCorrectVarName(obj_name, Name)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: obj_name

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: Name
        integer                                     :: i

        !Begin-----------------------------------------------------------------

        if(scan(obj_name, "_") .ne. 0 .and. scan(obj_name, "0") .ne. 0)then
            Name = obj_name(1:len_trim(obj_name)-6)
        else
            Name = trim(obj_name)
        endif

        do i = 1, len_trim(Name)
            if (Name(i:i) == ' ') then
                Name(i:i) =  '_'
            endif
        enddo

    end subroutine CheckAndCorrectVarName    
!--------------------------------------------------------------------------    
    
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


    function GetNetCDFValue(ValueIn, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------  
        real(8)             :: GetNetCDFValue      
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
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'        

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array2D = Me%Field(iP)%Value2DOut,                     &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        else if (Me%Field(iP)%Dim==3) then
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                 &
                                 WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'        

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array3D = Me%Field(iP)%Value3DOut,                     &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

        else 

            stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        endif


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleMecatorFormat - ERR90'


    end subroutine WriteFieldHDF5


    !----------------------------------------------------------------------

    !------------------------------------------------------------------------
        
    subroutine WriteFieldNetCDF(iP, iFinal)

        !Arguments-------------------------------------------------------------
        integer                                         :: iP, iFinal

        !Local-----------------------------------------------------------------
        real,  dimension(:,:,:), pointer                :: Field3D
        real,  dimension(:,:  ), pointer                :: Field2D        
        integer                                         :: STAT_CALL, i, j, k 
        integer                                         :: WorkIUB, WorkJUB, WorkKUB
        character(len=StringLength)                     :: NCDFName, LongName, StandardName, Units
        real                                            :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        !Begin-----------------------------------------------------------------
        
        !Bounds
        WorkIUB = Me%WorkSize%IUB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKUB = Me%WorkSize%KUB 



        if      (Me%Field(iP)%Dim==2) then
        
            allocate(Field2D(1:WorkIUB,1:WorkJUB))

            do j = 1, WorkJUB
            do i = 1, WorkIUB
                Field2D(i,j)=Me%Field(iP)%Value2DOut(i,j)
            enddo
            enddo                
             

            call BuildAttributes(trim(Me%Field(iP)%ID%Name), NCDFName, LongName, StandardName, &
                                       Units, ValidMin, ValidMax,        &
                                       MinValue, MaxValue, MissingValue, Float2D = Field2D)


            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,  &
                                  Name          = trim(NCDFName),           &
                                  LongName      = trim(LongName),           &
                                  StandardName  = trim(StandardName),       & 
                                  Units         = trim(Units),              &
                                  ValidMin      = ValidMin,                 &
                                  ValidMax      = ValidMax,                 &
                                  MinValue      = MinValue,                 &
                                  MaxValue      = MaxValue,                 &
                                  OutputNumber  = iFinal,                   &
                                  Array2D       = Field2D,                  &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'        
            
            deallocate(Field2D)

        else if (Me%Field(iP)%Dim==3) then

            allocate(Field3D(1:WorkIUB,1:WorkJUB,1:WorkKUB))   

            do k = 1, WorkKUB
            do j = 1, WorkJUB
            do i = 1, WorkIUB
                Field3D(i,j,k)=Me%Field(iP)%Value3DOut(i,j,k)
            enddo
            enddo                
            enddo

        
            call BuildAttributes(trim(Me%Field(iP)%ID%Name), NCDFName,                  &
                                      LongName, StandardName,                           &
                                      Units, ValidMin, ValidMax,                        &
                                      MinValue, MaxValue, MissingValue,                 &
                                      Float3D = Field3D)


            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,              &
                                  Name          = trim(NCDFName),                       &
                                  LongName      = trim(LongName),                       &
                                  StandardName  = trim(StandardName),                   & 
                                  Units         = trim(Units),                          &
                                  ValidMin      = ValidMin,                             &
                                  ValidMax      = ValidMax,                             &
                                  MinValue      = MinValue,                             &
                                  MaxValue      = MaxValue,                             &
                                  OutputNumber  = iFinal,                               &
                                  Array3D       = Field3D,                              &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
            deallocate(Field3D)        

        else 

            stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        endif

    end subroutine WriteFieldNetCDF


    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    

    character(len=19) function TimeToStringV2(Date)

        !Arguments-------------------------------------------------------------
        type(T_Time)                            :: Date
        real,    dimension(6)                   :: AuxTime
        character(len=4)                        :: CharYear
        character(len=2)                        :: CharMonth
        character(len=2)                        :: CharDay
        character(len=2)                        :: CharHour
        character(len=2)                        :: CharMinute
        character(len=2)                        :: CharSecond

        !Begin-----------------------------------------------------------------

        call ExtractDate(Date, Year     = AuxTime(1), Month  = AuxTime(2), &
                               Day      = AuxTime(3), Hour   = AuxTime(4), &
                               Minute   = AuxTime(5), Second = AuxTime(6))
        
        write(CharYear,  '(i4)')int(AuxTime(1))
        write(CharMonth, '(i2)')int(AuxTime(2))
        write(CharDay,   '(i2)')int(AuxTime(3))
        write(CharHour,  '(i2)')int(AuxTime(4))
        write(CharMinute,'(i2)')int(AuxTime(5))
        write(CharSecond,'(i2)')int(AuxTime(6))

        if(len_trim(trim(adjustl(CharMonth)))   < 2)then 
            CharMonth = "0"//trim(adjustl(CharMonth))
        endif
        
        if(len_trim(trim(adjustl(CharDay)))     < 2)then 
            CharDay = "0"//trim(adjustl(CharDay))
        endif

        if(len_trim(trim(adjustl(CharHour)))    < 2)then 
            CharHour = "0"//trim(adjustl(CharHour))
        endif

        if(len_trim(trim(adjustl(CharMinute)))  < 2)then 
            CharMinute = "0"//trim(adjustl(CharMinute))
        endif

        if(len_trim(trim(adjustl(CharSecond)))  < 2)then 
            CharSecond = "0"//trim(adjustl(CharSecond))
        endif

        TimeToStringV2 = CharYear//"-"//CharMonth//"-"//CharDay//" "//&
                         CharHour//":"//CharMinute//":"//CharSecond

    end function    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    
    subroutine KillNetCDFCF_2_HDF5MOHID
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers, ip
        
        !Begin-----------------------------------------------------------------

            if (Me%OutHDF5) then
                !Close HDF5 File
                call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
            endif
            
            if (Me%OutNetCDF) then
                
                !Close NetCDF File
                call KillNETCDF(Me%NetCDF_Out%ObjNetCDF, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
            endif

            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

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
