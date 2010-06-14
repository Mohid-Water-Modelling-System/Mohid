!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToHDF5
! MODULE        : MOG2D Format
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Januray 2008
! REVISION      : Rodrigo Fernandes
! DESCRIPTION   : Module to convert MOG2DFormat (NETCDF) files into HDF5 format.
!
!------------------------------------------------------------------------------
!DataFile
!
!   OUTPUT_GRID_FILENAME        : char              -           !Path to grid data file generated from MOG2D file
!   OUTPUTFILENAME              : char              -           !Path to HDF5 file generated from MOG2D file
!   INPUT_GRID_FILENAME         : char              -           !Path to netcdf griddata file
!
!   <<begin_input_files>>
!   ALADIN_BULKCLOUD_OPASYMP_19723_20088.nc
!   ALADIN_BULKHUMI_OPASYMP_19723_20088.nc
!   ...
!   ... (see below for available fields)
!   <<end_input_files>>


!
!         ---MOG2D---    ---MOHID NAME---
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




!Dimensions:
!X = 251;
!Y = 245;
!
!Variables:
!float lon(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "longitude"
!           units = "degree_east"
!           valid_min = ...
!           valid_max = ...
!float lat(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!     Attributes:
!           long_name = "latitude"
!           units = "degree_north"
!           valid_min = ...
!           valid_max = ...
!float x_elev(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!     Attributes:
!           long_name = "x_elev"
!           units = "m"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = 999.0
!float y_elev(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "y_elev"
!           units = "m"
!       associate = "time lat lon"
!       axis = "YX"
!       missing_value = -9999.0
!       _FillValue = -9999.0
!       valid_min = ...
!       valid_max = 999.0

!float x_charge(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "x_charge"
!           units = "m"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = ...

!float y_charge(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "y_charge"
!           units = "m"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = ...

!float x1_oe(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "x1_oe"
!           units = "m2/s"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = 999.9

!float x2_oe(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "x2_oe"
!           units = "m2/s"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = 999.9

!float y1_oe(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "y1_sn"
!           units = "m2/s"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = 999.9

!float y2_oe(Y,X)
!      Dimensions:
!           Y = 245
!           X = 251
!      Attributes:
!           long_name = "y2_sn"
!           units = "m2/s"
!           associate = "time lat lon"
!           axis = "YX"
!           missing_value = -9999.0
!           _FillValue = -9999.0
!           valid_min = ...
!           valid_max = 999.9


Module ModuleMOG2DFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid
    use ModuleHorizontalMap

    ! Manages NetCDF files
    use netcdf90
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertMOG2DFormat
    private ::      ReadOptions
    private ::      OpenAndReadMOG2DFile
    private ::          OpenAndReadGridMOG2D
    private ::          GetNamesInHDF
    private ::          Open_HDF5_OutPut_File
    private ::          ReadMOG2DFile
    private ::      KillMOG2DFormat
    
    !Parameters---------------------------------------------------------------
    character(LEN = StringLength), parameter    :: input_files_begin   = '<<begin_input_files>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<end_input_files>>'

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
        integer                                 :: OutputNumber         = 1
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: Next
    end type  T_Field

    private :: T_MOG2DFormat
    type T_MOG2DFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
 !       integer                                 :: ObjGridData          = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit, ClientNumber
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: InputGridFile
        character(len=PathLength)               :: OutputFileName
        integer                                 :: NumDates
        integer                                 :: ReadOptionType
        integer                                 :: imax, jmax
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:,:),       pointer     :: XX, YY
        real, dimension(:,:),       pointer     :: XX_IE, YY_IE
        real, dimension(:, :) ,     pointer     :: SurfaceElevation
        integer, dimension(:,:),  pointer       :: WaterPoints
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Date),               pointer     :: FirstDate 
        type(T_Time)                            :: RefDateTime    
        integer, dimension(12)                  :: Instants(1:12) = 0
    end type

    type(T_MOG2DFormat), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ConvertMOG2DFormat(EnterDataID, ClientNumber, STAT)

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

        call SetDate (Me%RefDateTime, Year=2004, Month=1, Day=1, Hour=0, Minute=0, Second=0) 

        call StartComputeTime(Me%ObjTime, Me%RefDateTime, Me%RefDateTime, DT = 0.0,    &
                                 VariableDT = .false., STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMOG2DFormat - ERR02a'

        call ReadOptions

        call OpenAndReadMOG2DFile

        call KillMOG2DFormat


        STAT = SUCCESS_


    end subroutine ConvertMOG2DFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        !Read output griddata filename
        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                             &
                     ClientModule = 'ConvertMOG2DFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertMOG2DFormat - ERR10'

        !Read output filename
        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ConvertMOG2DFormat',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertMOG2DFormat - ERR20'

        !Read input MOG2D netcdf gridded data file to generate the griddata
        call GetData(Me%InputGridFile,                                              &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'INPUT_GRID_FILENAME',                          &
                     ClientModule = 'ConvertMOG2DFormat',                         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ConvertMOG2DFormat - ERR80'

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine OpenAndReadMOG2DFile

        !Local-----------------------------------------------------------------

        !Begin----------------------------------------------------------------

        call OpenAndReadGridMOG2D

        call ConstructGrid

        call Open_HDF5_OutPut_File

        call OpenAndReadMOG2DFields

    end subroutine OpenAndReadMOG2DFile
    
    
    !------------------------------------------------------------------------

    !Here we're fooling converttohdf5 in thinking there's an all-water bathym
    subroutine OpenAndReadGridMOG2D

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:  ), allocatable     :: Aux2D
        real, dimension(:    ), allocatable     :: Aux1DLon, Aux1DLat
        logical                                 :: exist
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: WILB, WIUB, WJLB,WJUB
        integer                                 :: i, j, n
        integer                                 :: ncid, status, dimid


        !Begin----------------------------------------------------------------

        !Verifies if file exists
        inquire(file = Me%InputGridFile, exist = exist)
i1:     if (exist) then

            status=NF90_OPEN(trim(Me%InputGridFile),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR10'

            status=NF90_INQ_DIMID(ncid,"X",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR40'

            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%jmax)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR50'

            status=NF90_INQ_DIMID(ncid,"Y",dimid)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR20'

            status=NF90_INQUIRE_DIMENSION(ncid,dimid,len = Me%imax)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR30'

          
            !The border cells are not considered because with the available data
            !is not possible to compute the horizontal position of the cells corners
            Me%WorkSize%ILB = 1
            Me%WorkSize%IUB = Me%imax - 2
                
            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = Me%jmax - 2

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

            allocate(Me%WaterPoints (ILB:IUB, JLB:JUB))
            allocate(Me%Bathymetry  (ILB:IUB, JLB:JUB))
            
            !Longitude == XX == j
            !Latitude == YY == i
            allocate(Me%XX          (ILB:IUB, JLB:JUB))
            allocate(Me%YY          (ILB:IUB, JLB:JUB))

            allocate(Aux2D(Me%jmax, Me%imax))

            !There's no mask in the MOG2D netcdf variables so it's all water!
            do i=WILB,WIUB
            do j=WJLB,WJUB
               Me%WaterPoints(i,j) = WaterPoint
            enddo
            enddo


            !There's no orography in the MOG2D outputs so there's no real bathymetry
            do i=WILB,WIUB
            do j=WJLB,WJUB
                Me%Bathymetry(i,j) = 100.
            enddo
            enddo

            status = nf90_inq_varid(ncid, 'lat', n)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR100'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR110'



            !Calculating Latitude on the cell faces

            do i=WILB,WIUB +1
            do j=WJLB,WJUB +1
               Me%YY(i,j) = (Aux2D(j,i) + Aux2D(j, i+1) )/2.
            enddo
            enddo




            status = nf90_inq_varid(ncid, 'lon', n)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR120'

            status = NF90_GET_VAR(ncid,n,Aux2D)
            if (status /= nf90_noerr) stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR130'

            !Calculate Longitude on the cell faces

            do i=WILB,WIUB +1
            do j=WJLB,WJUB +1
               Me%XX(i,j) = (Aux2D(j,i) + Aux2D(j+1, i))/2.
            enddo
            enddo

            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR140'

            deallocate(Aux2D)

        else i1

            write (*,*) "The input grid file do not exist : ",trim(Me%InputGridFile)
            stop 'OpenAndReadGridMOG2D - ModuleMOG2DFormat - ERR150'

        endif i1

    end subroutine OpenAndReadGridMOG2D
    
    !------------------------------------------------------------------------

    
    !------------------------------------------------------------------------

    subroutine ConstructGrid
        
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

        call WriteGridData (FileName       = Me%GridFileName,                              &
                             ConnectionX             = Me%XX,                                        &
                             ConnectionY             = Me%YY,                                        &
                             COMENT1        = 'Grid Data created from MOG2D NetCDF file',   &
                             COMENT2        = trim(Me%InputGridFile),                       &
                             WorkSize       = WorkSize2D,                                   &
                             CoordType      = 4,                                            &
                             Xorig          = 0.,                                           &
                             Yorig          = 0.,                                           &
                             Zone           = 29,                                           &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Me%YY(1,1),                                   &
                             Longitude      = Me%XX(1,1),                                   &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%Bathymetry,                                &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleMOG2DFormat - ERR10'


        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMOG2DFormat - ERR20'


    end subroutine ConstructGrid

    
    !------------------------------------------------------------------------


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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR02'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR04'            
   
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints", "-",    &
                              Array2D = Me%WaterPoints,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR08'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMOG2DFormat - ERR09'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------


    !------------------------------------------------------------------------

    subroutine OpenAndReadMOG2DFields

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

                    if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMOG2DFields - ModuleMOG2DFormat - ERR10'

                    inquire(file = InputFile, exist = exist)
       
i1:                 if (exist) then
                                                    
                        call ReadMOG2DFile (InputFile) 

                    endif i1
                enddo

            else BF

                stop 'OpenAndReadMOG2DFields - ModuleMOG2DFormat - ERR20'

            end if BF

            call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OpenAndReadMOG2DFields - ModuleMOG2DFormat - ERR30'

        else   IS

            stop 'OpenAndReadMOG2DFields - ModuleMOG2DFormat - ERR40'

        end if IS



    end subroutine OpenAndReadMOG2DFields


    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    subroutine ReadMOG2DFile(InputFile)
    !
    !Local subsubroutines : 
    !   CheckName, check
    !   OutputInstants, check
    !   GetPropertyName, 
    !   GetPropertyIDNumber
    !   WriteHDF5Field

        !Arguments-------------------------------------------------------------
        character (Len=*)                       :: InputFile
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:    ), pointer       :: Aux2D, Aux2DRead
        integer                                 :: ni, n, nInst, iOut, i, j
        integer                                 :: ncid, status, dimid
        integer                                 :: nDimensions
        integer                                 :: nDims, nVars, nAtrr, xtype
        integer                                 :: WestON, EastON, SouthON, NorthON
        character (len=80)                      :: nameAux
        character(Len=StringLength)             :: NameInHDF, AuxNameInHDF
        type (T_Time)                           :: FieldTime
        real, dimension(:,:    ), pointer       :: X_Elev, Y_Elev, Amplitude, Phase
        real, dimension(:,:    ), pointer       :: x1_oe, x2_oe, y1_sn, y2_sn


        !Begin----------------------------------------------------------------

        !verifies if file exists
        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
        if (status /= nf90_noerr) stop 'ReadMOG2DFile - ModuleMOG2DFormat - ERR10'


        status=NF90_INQUIRE(ncid, nDims, nVars, nAtrr)
        if (status /= nf90_noerr) stop 'ReadMOG2DFile - ModuleMOG2DFormat - ERR60'

d0:     do n=1,nVars

            status=NF90_INQUIRE_VARIABLE(ncid, n, nameAux, xtype, nDimensions)
            if (status /= nf90_noerr) stop 'OpenAndReadMOG2DFileV1 - ModuleMOG2DFormat - ERR70'

            if (nDimensions == 2) then

                !Allocate memory
                allocate(Aux2DRead(Me%jmax, Me%imax))
                allocate(Aux2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                Aux2D(:,:) = FillValueReal
                
                !Fill memory
                status = NF90_GET_VAR(ncid,n,Aux2DRead)
                if (status /= nf90_noerr) stop 'ReadMOG2DFile - ModuleMOG2DFormat - ERR80'

            elseif (nDimensions == 1) then
                !Grid properties already written
                cycle
            endif


                    Call GetNamesInHDF(nameAux, NameInHDF)

                    iOut = 1

                    FieldTime = Me%RefDateTime

                    if (nDimensions == 2) then
                           
                                        
                        !The boundary cells are not read
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB

!                            if (Me%WaterPoints(i, j) == WaterPoint) then
                             if (Aux2DRead(j+1,i+1) < 999.9 .or. Aux2DRead(j+1,i+1) > 1000 ) then
                                

                                Aux2D(i, j) = Aux2DRead(j+1,i+1)
                             else
                                Aux2D(i, j) = null_real
                            endif

                        enddo
                        enddo

                        call WriteHDF5Field(FieldTime, NameInHDF, iOut, Aux2D = Aux2D)

                        select case (NameInHDF)
                        
                            case('x_elev')
                                allocate(X_Elev(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                                X_Elev(:,:) = Aux2D(:,:)
                                                             
                            case('y_elev')
                                allocate(Y_Elev(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                                Y_Elev(:,:) = Aux2D(:,:)

                        end select

                       allocate(Amplitude(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                       allocate(Phase(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                       !                        Phase  = SQRT(X_Elev*X_Elev + Y_Elev*Y_Elev)


                    endif


            if (nDimensions == 2) then
                deallocate(Aux2DRead)
                deallocate(Aux2D)
            endif

        enddo d0


!-------Write water level amplitude and phase properties
                
        !-------Amplitude        
        iOut = 1
        FieldTime = Me%RefDateTime
        AuxNameInHDF = 'Amplitude'

        !The boundary cells are not read
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
!            if (Me%WaterPoints(i, j) == WaterPoint) then
            if (X_Elev(i,j) /= null_real .and. Y_Elev(i,j) /= null_real) then
                                    
                Amplitude(i, j) = SQRT(X_Elev(i,j)*X_Elev(i,j) + Y_Elev(i,j)* Y_Elev(i,j))

            else
                Amplitude(i, j) = null_real
            endif

        enddo
        enddo

        call WriteHDF5Field(FieldTime, AuxNameInHDF, iOut, Aux2D = Amplitude)



        !-------Phase
        iOut = 1
        FieldTime = Me%RefDateTime
        AuxNameInHDF = 'Phase'

        !The boundary cells are not read
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            
!            if (Me%WaterPoints(i, j) == WaterPoint) then
            if (X_Elev(i,j) /= null_real .and. Y_Elev(i,j) /= null_real) then
                                    
                Phase(i, j) = ATAN2(-Y_Elev(i,j),X_Elev(i,j))*180 / Pi

            else
                Phase(i, j) = null_real

            endif

        enddo
        enddo

        call WriteHDF5Field(FieldTime, AuxNameInHDF, iOut, Aux2D = Phase)


!---------------------------


        !Closes Netcdf file handle
        status=NF90_CLOSE(ncid)
        if (status /= nf90_noerr) stop 'ReadMOG2DFile - ModuleMOG2DFormat - ERR91'


    end subroutine ReadMOG2DFile
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetNamesInHDF(MOG2DName, NameInHDF)
        
        !Arguments-----------------------------------------------------------
        character(Len=*)                :: MOG2DName
        character(Len=StringLength)     :: NameInHDF
        
        !Begin-----------------------------------------------------------------


!provisorio, enqanto nao seencontra correspondencia entre o mohid e o mog2d
NameInHDF = trim(MOG2DName)

!        select case(trim(MOG2DName))
!
!            case('lon')
!                NameInHDF = GetPropertyName(WindVelocityX_)
!                CheckName = .true.
!
!            case('lat')
!
!                NameInHDF = GetPropertyName(WindVelocityX_)
!                CheckName = .true.
!
!            case('x_elev')
!
!                NameInHDF = GetPropertyName(CloudCover_)
!                CheckName = .true.
!
!            case('y_elev')
!
!                NameInHDF = GetPropertyName(RelativeHumidity_)
!                CheckName = .true.
!
!            case('x_charge')
!
!                NameInHDF = GetPropertyName(NonSolarFlux_)
!                CheckName = .true.
!
!            case('y_charge')
!
!                NameInHDF = GetPropertyName(AtmosphericPressure_)
!                CheckName = .true.
!
!            case('x1_oe')
!
!               NameInHDF = GetPropertyName(SolarRadiation_)
!               CheckName = .true.
!
!
!            case('x2_oe')
!
!                NameInHDF = GetPropertyName(AirTemperature_)
!                CheckName = .true.
!
!
!            case('y1_sn')
!
!               NameInHDF = GetPropertyName(WindModulus_)
!               CheckName = .true.
!
!            case('y2_sn')
!
!                NameInHDF = GetPropertyName(Precipitation_)
!                CheckName = .true.
!
!            case default
!                
!               CheckName = .false.
!
!        end select


    end subroutine GetNamesInHDF
    
    
    
    !------------------------------------------------------------------------
    subroutine WriteHDF5Field(FieldTime, NameInHDF, iOut, Aux2D, Aux3D)

        !Arguments-------------------------------------------------------------
        type (T_Time)                                   :: FieldTime
        character(Len=StringLength)                     :: NameInHDF
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

        PropUnits = Units(NameInHDF)

        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR10'

        if      (present(Aux2D)) then

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(NameInHDF),              &
                                 trim(NameInHDF),trim(PropUnits), Array2D = Aux2D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR20'

        else if (present(Aux3D)) then

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(NameInHDF),              &
                                 trim(NameInHDF),trim(PropUnits), Array3D = Aux3D,      &
                                 OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR30'

        else 

            stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR40'

        endif

        call GetHDF5GroupExist (Me%ObjHDF5, "/Time", Exist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR35'

        if (Exist) then
        
            call GetHDF5GroupNumberOfItems (Me%ObjHDF5,  "/Time", nItems, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR37'
        
        endif    


        if (.not. Exist .or. iOut > nItems) then

            !Writes current time
            call ExtractDate   (FieldTime, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                           AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR50'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                                 Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR60'

            AuxName = NameInHDF

        else if (iOut == nItems) then

            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR100'

            call HDF5ReadData  (Me%ObjHDF5, "/Time", "Time",                               &
                                 Array1D = TimePtr, OutputNumber = iOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR110'


            call SetDate   (AuxFieldTime, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                          AuxTime(4), AuxTime(5), AuxTime(6))


            if (FieldTime /= AuxFieldTime) then
!                write(*,*) 'The time instants of property ',trim(NameInHDF)
!                write(*,*) 'are not consistent with property ',trim(AuxName)
                stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR120'   
            endif             
            
        endif

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5Field - ModuleMOG2DFormat - ERR90'


    end subroutine WriteHDF5Field
 
    !------------------------------------------------------------------------
 

    !--------------------------------------------------------------------------

    character(Len=StringLength) function Units(NameInHDF)
        
        !Arguments-----------------------------------------------------------
        character(Len=StringLength)     :: NameInHDF, AuxNameInHDF
        !Local-----------------------------------------------------------------
        integer                         :: MohidID
        
        !Begin-----------------------------------------------------------------


        !MohidID = GetPropertyIDNumber(NameInHDF)

AuxNameInHDF = trim(NameInHDF)
        select case(AuxNameInHDF)


            case('Amplitude')
               
                Units = 'm'

            case('Phase')
                
                Units = 'º'

            case('lon')

                Units = 'º'

            case('lat')

                Units = 'º'

            case('x_elev')

                Units = 'm'

            case('y_elev')

                Units = 'm'

            case('x_charge')

                Units = 'm'

            case('y_charge')

                Units = 'm'

            case('x1_oe')

                Units = 'm2/s'

            case('x2_oe')

                Units = 'm2/s'

            case('y1_sn')

                Units = 'm2/s'

            case('y2_sn')

                Units = 'm2/s'

        end select


    end function Units
    
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine KillMOG2DFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMOG2DFormat - KillMOG2DFormat - ERR50'
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMOG2DFormat - KillMOG2DFormat - ERR60'

        deallocate(Me%XX)
        deallocate(Me%YY)
        deallocate(Me%Bathymetry)
        deallocate(Me%WaterPoints)

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillMOG2DFormat - KillMOG2DFormat - ERR70'

        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillMOG2DFormat

    !--------------------------------------------------------------------------


end module ModuleMOG2DFormat









