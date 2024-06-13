!****************************************************************************!
!  RasterFunctions.f90                                                       !
!                                                                            !
!****************************************************************************!
!                                                                            !
!  PROGRAM    : RGBComposite                                                 !
!  PURPOSE    : Provide Global Raster Functions for multiple solutions       !
!  AUTHOR     : Fl�vio Santos                                                !
!  AFFILIATION: HIDROMOD                                                     !
!                                                                            !
!****************************************************************************!
    
Module RasterFunctions

    use ModuleGlobalData
    use gdal
    use proj
    use fortranc
    use iso_c_binding
    
    implicit none
        
    type T_Matrix
        real(4)   , dimension(:,:), pointer :: DataR4 => null()
        real(8)   , dimension(:,:), pointer :: DataR8 => null()
        integer(2), dimension(:,:), pointer :: DataI2 => null()
        integer(4), dimension(:,:), pointer :: DataI4 => null()
    end type
    
    type T_Band
        type(gdalrasterbandh)                :: BandObj
        integer                              :: DataType
        type(T_Matrix)                       :: DataMatrix
    end type
    
    type T_Raster
        character(len=PathLength)            :: FileName
        type (gdaldriverh)                   :: Driver
        type (gdaldataseth)                  :: Dataset
        real(4)                              :: OriginTL(2)
        real(kind=c_double)                  :: GeoTransform(6)
        character(len=800)                   :: ProjRef
        integer                              :: ii
        integer                              :: jj
        real(4), dimension(:), pointer       :: Origin => null()
        type (T_Band), dimension(:), pointer :: Bands => null()
    end type
    
    public
    public :: OpenRaster
    
    public :: CreateRaster
    interface CreateRaster
        module procedure CreateRasterI2
        module procedure CreateRasterI4
        module procedure CreateRasterR4
        module procedure CreateRasterR8
    end interface
    
    !public :: AverageRasters
    
    contains
    
    
    !----------------------------------------------------------------------
    subroutine OpenRaster(FileName, DriverName, Raster)
        !Arguments------------------------------------------------
        character(len=*), intent(IN)    :: FileName
        character(len=*),          intent(IN)    :: DriverName
        type(T_Raster), pointer,   intent(OUT)   :: Raster
        !Local----------------------------------------------------    
        type(gdaldriverh)                        :: driver
        type(gdaldataseth)                       :: ds
        integer                                  :: nbands, b
        integer                                  :: ierr
        integer                                  :: datatype
        real(kind=c_double)                      :: geotrans(6)
        real(kind=c_double)                      :: x1, y1
        !Body-----------------------------------------------------
        
        allocate(Raster)
        
        Raster%FileName = trim(adjustl(FileName))
    
        call gdalallregister() ! register all available gdal drivers "https:\\gdal.org\doxygen\gdal_8h.html#a9d40bc998bd6ed07ccde96028e85ae26"
        driver = gdalgetdriverbyname(DriverName//char(0))
        
        if (.NOT.gdalassociated(driver)) then
            write(*,*) 'Error getting GeoTIFF driver from gdal'
            stop 'OpenRaster - RasterFunctions - ERR01'
        endif
        Raster%Driver = driver
                
        ds = gdalopen(trim(adjustl(FileName))//char(0), GA_ReadOnly) !read a gdal dataset from a file, simplified Fortran interface
        if (.NOT.gdalassociated(ds)) then
            write(*,*) 'Error opening dataset on file ',trim(adjustl(FileName))
            stop 'OpenRaster - RasterFunctions - ERR02'
        endif
        Raster%Dataset = ds
        
        Raster%ii = gdalgetrasterxsize(ds) !Fetch raster width in pixels.
        Raster%jj = gdalgetrasterysize(ds) !Fetch raster height in pixels.        
        if (Raster%ii <= 0 .or. Raster%jj <= 0) then
            write(*,*) 'Error reading size from GeoTIFF dataset on file ', trim(adjustl(FileName))
            stop 'OpenRaster - RasterFunctions - ERR03'
        endif
        
        Raster%ProjRef = strtofchar(gdalgetprojectionref(ds), 800) !Fetch the projection definition string.
        ierr = gdalgetgeotransform(ds, geotrans) !Fetch the affine transformation coefficients for transforming between pixel/line (P,L) raster space, and projection coordinates (Xp,Yp) space.
        if (ierr /= 0) then
            write(*,*) 'Error reading Projection Reference from GeoTIFF on file ',trim(adjustl(FileName))
            stop 'OpenRaster - RasterFunctions - ERR04'
        end if
        Raster%GeoTransform = geotrans
        
        call gdalapplygeotransform_f(geotrans, 0.5_c_double, 0.5_c_double, x1, y1) !Apply GeoTransform to x/y coordinate. -> center of cell in top left corner
        Raster%OriginTL = [x1, y1]
                
        nbands = gdalgetrastercount(ds)
        if (nbands <= 0) then
            write(*,*) 'Error counting number of bands from GeoTIFF dataset on file ', trim(adjustl(FileName))
            stop 'OpenRaster - RasterFunctions - ERR04'
        endif
        
        allocate(Raster%Bands(1:nbands))
        
        do b=1, size(Raster%Bands)
            Raster%Bands(b)%BandObj = gdalgetrasterband(ds, b) !Fetch a band object for a dataset.
            if (.NOT.gdalassociated(Raster%Bands(b)%BandObj)) then
                write(*,*) 'Error getting raster band from GeoTIFF dataset on file ',trim(adjustl(FileName))
                stop 'OpenRaster - RasterFunctions - ERR10'
            endif
            
            !Get name of data type. -> https:\\naturalatlas.github.io\node-gdal\classes\Constants%20(GDT).html
            datatype = gdalgetrasterdatatype(Raster%Bands(b)%BandObj)
            Raster%Bands(b)%DataType = datatype
            if (datatype .eq. GDT_Float32) then
                allocate(Raster%Bands(b)%DataMatrix%DataR4(1:Raster%ii, 1:Raster%jj))
                
                ierr = gdalrasterio_f(Raster%Bands(b)%BandObj,                  & !Band object
                                      erwflag = GF_Read,                        & !Read/Write flag -> GF_Read or GF_Write
                                      ndsxoff = 0,                              & !pixel offset to the top left corner of the region of the band to be accessed. This would be zero to start from the left side.
                                      ndsyoff = 0,                              & !line offset to the top left corner of the region of the band to be accessed. This would be zero to start from the top.
                                      pbuffer = Raster%Bands(b)%DataMatrix%DataR4)
            elseif (datatype .eq. GDT_Float64) then
                allocate(Raster%Bands(b)%DataMatrix%DataR8(1:Raster%ii, 1:Raster%jj))
                
                ierr = gdalrasterio_f(Raster%Bands(b)%BandObj,                  & !Band object
                                      erwflag = GF_Read,                        & !Read/Write flag -> GF_Read or GF_Write
                                      ndsxoff = 0,                              & !pixel offset to the top left corner of the region of the band to be accessed. This would be zero to start from the left side.
                                      ndsyoff = 0,                              & !line offset to the top left corner of the region of the band to be accessed. This would be zero to start from the top.
                                      pbuffer = Raster%Bands(b)%DataMatrix%DataR8)
            elseif (datatype .eq. GDT_Int16) then
                allocate(Raster%Bands(b)%DataMatrix%DataI2(1:Raster%ii, 1:Raster%jj))
                
                ierr = gdalrasterio_f(Raster%Bands(b)%BandObj,                  & !Band object
                                      erwflag = GF_Read,                        & !Read/Write flag -> GF_Read or GF_Write
                                      ndsxoff = 0,                              & !pixel offset to the top left corner of the region of the band to be accessed. This would be zero to start from the left side.
                                      ndsyoff = 0,                              & !line offset to the top left corner of the region of the band to be accessed. This would be zero to start from the top.
                                      pbuffer = Raster%Bands(b)%DataMatrix%DataI2)
            elseif (datatype .eq. GDT_Int32) then
                allocate(Raster%Bands(b)%DataMatrix%DataI4(1:Raster%ii, 1:Raster%jj))
                
                ierr = gdalrasterio_f(Raster%Bands(b)%BandObj,                  & !Band object
                                      erwflag = GF_Read,                        & !Read/Write flag -> GF_Read or GF_Write
                                      ndsxoff = 0,                              & !pixel offset to the top left corner of the region of the band to be accessed. This would be zero to start from the left side.
                                      ndsyoff = 0,                              & !line offset to the top left corner of the region of the band to be accessed. This would be zero to start from the top.
                                      pbuffer = Raster%Bands(b)%DataMatrix%DataI4)
            else
                write(*,*) 'Error reading datatype of ',trim(adjustl(FileName))
                write(*,*) strtofchar(gdalgetdatatypename(datatype), StringLength)
                stop 'OpenRaster - RasterFunctions - ERR11'
            end if
            
            if (ierr /= 0) then
                write(*,*) 'Error reading data from GeoTIFF dataset on file ',trim(adjustl(FileName))
                stop 'OpenRaster - RasterFunctions - ERR20'
            end if
            
        end do
        
        call gdalclose(ds)
        
    end subroutine OpenRaster
    !----------------------------------------------------------------------
    
    !----------------------------------------------------------------------
    subroutine CreateRasterI2(FileName, DriverName, RasterWidth, RasterHeight, DataType, NumberBands, Projection, GeoTrans, DataMatrix3D)
        !Arguments------------------------------------------------
        character(len=*),            intent(IN)    :: FileName
        character(len=*),            intent(IN)    :: DriverName
        integer,                     intent(IN)    :: RasterWidth
        integer,                     intent(IN)    :: RasterHeight
        integer,                     intent(IN)    :: DataType
        integer,                     intent(IN)    :: NumberBands
        character(len=*),            intent(IN)    :: Projection
        real(kind=c_double),         intent(IN)    :: GeoTrans(6)
        integer(2), dimension(:,:,:),intent(INOUT)    :: DataMatrix3D
        !Local----------------------------------------------------
        type(gdaldriverh)                        :: driver
        type(gdaldataseth)                       :: ds
        integer                                  :: i, ierr
        !integer                                  :: x_dim, y_dim
        !integer(2), dimension(:,:,:), allocatable:: DataMatrix3D
        logical                                  :: exists
        !Body-----------------------------------------------------    
    
        call gdalallregister()
        driver = gdalgetdriverbyname(DriverName//char(0))
        if (.NOT.gdalassociated(driver)) then
            write(*,*) 'Error getting GeoTIFF driver from gdal'
            stop 'OpenRaster - CreateRaster - ERR01'
        endif
        
        !delete filename if exists
        !inquire(file=FileName, exist=exists)
        !if (exists) then
        !    call system("del "//trim(FileName))
        !endif
        
        ds = gdalcreate(driver, trim(adjustl(FileName))//char(0), &
                        RasterWidth, RasterHeight, NumberBands, DataType, c_ptr_ptr_getobject(c_ptr_ptr_new((/('',i=1,0)/))))
        if (.NOT.gdalassociated(ds)) then
            write(*,*) 'Error creating dataset on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR02'
        endif
        
        !projection
        ierr = gdalsetprojection(ds, trim(adjustl(Projection)))
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL projection on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR03'
        end if
        
        !geotransformation
        ierr = gdalsetgeotransform(ds, GeoTrans)
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL geotransformation properties on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR04'
        end if
        
        ierr = gdaldatasetrasterio_f(ds, GF_Write, 0, 0, DataMatrix3D)
        if (ierr /= 0) then
            write(*,*) 'Error setting matrix on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR04'
        end if
        
        !deallocate(DataMatrix3D)
        call gdalclose(ds)
    
    end subroutine CreateRasterI2
    !----------------------------------------------------------------------
    
    !----------------------------------------------------------------------
    subroutine CreateRasterI4(FileName, DriverName, RasterWidth, RasterHeight, DataType, NumberBands, Projection, GeoTrans, DataMatrix3D)
        !Arguments------------------------------------------------
        character(len=*),            intent(IN)    :: FileName
        character(len=*),            intent(IN)    :: DriverName
        integer,                     intent(IN)    :: RasterWidth
        integer,                     intent(IN)    :: RasterHeight
        integer,                     intent(IN)    :: DataType
        integer,                     intent(IN)    :: NumberBands
        character(len=*),            intent(IN)    :: Projection
        real(kind=c_double),         intent(IN)    :: GeoTrans(6)
        integer(4), dimension(:,:,:),intent(INOUT)    :: DataMatrix3D
        !Local----------------------------------------------------
        type(gdaldriverh)                        :: driver
        type(gdaldataseth)                       :: ds
        integer                                  :: i, ierr
        logical                                  :: exists
        !Body-----------------------------------------------------    
    
        call gdalallregister()
        driver = gdalgetdriverbyname(DriverName//char(0))
        if (.NOT.gdalassociated(driver)) then
            write(*,*) 'Error getting GeoTIFF driver from gdal'
            stop 'OpenRaster - CreateRaster - ERR01'
        endif
        
        !delete filename if exists
        !inquire(file=FileName, exist=exists)
        !if (exists) then
        !    call system("del "//trim(FileName))
        !endif
        
        ds = gdalcreate(driver, trim(adjustl(FileName))//char(0), &
                        RasterWidth, RasterHeight, NumberBands, DataType, c_ptr_ptr_getobject(c_ptr_ptr_new((/('',i=1,0)/))))
        if (.NOT.gdalassociated(ds)) then
            write(*,*) 'Error creating dataset on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR02'
        endif
        
        !projection
        ierr = gdalsetprojection(ds, trim(adjustl(Projection)))
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL projection on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR03'
        end if
        
        !geotransformation
        ierr = gdalsetgeotransform(ds, GeoTrans)
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL geotransformation properties on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR04'
        end if
                
        ierr = gdaldatasetrasterio_f(ds, GF_Write, 0, 0, DataMatrix3D)
        if (ierr /= 0) then
            write(*,*) 'Error setting matrix on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR04'
        end if
        
        call gdalclose(ds)
    
    end subroutine CreateRasterI4
    !----------------------------------------------------------------------
    
    !----------------------------------------------------------------------
    subroutine CreateRasterR4(FileName, DriverName, RasterWidth, RasterHeight, DataType, NumberBands, Projection, GeoTrans, DataMatrix3D)
        !Arguments------------------------------------------------
        character(len=*),          intent(IN)    :: FileName
        character(len=*),          intent(IN)    :: DriverName
        integer,                   intent(IN)    :: RasterWidth
        integer,                   intent(IN)    :: RasterHeight
        integer,                   intent(IN)    :: DataType
        integer,                   intent(IN)    :: NumberBands
        character(len=*),          intent(IN)    :: Projection
        real(kind=c_double),       intent(IN)    :: GeoTrans(6)
        real(4), dimension(:,:,:), intent(INOUT)    :: DataMatrix3D
        !Local----------------------------------------------------
        type(gdaldriverh)                        :: driver
        type(gdaldataseth)                       :: ds
        integer                                  :: i, ierr
        !integer                                  :: x_dim, y_dim
        !real(4), dimension(:,:,:), allocatable   :: DataMatrix3D
        logical                                  :: exists
        !Body-----------------------------------------------------    
    
        call gdalallregister()
        driver = gdalgetdriverbyname(DriverName//char(0))
        if (.NOT.gdalassociated(driver)) then
            write(*,*) 'Error getting GeoTIFF driver from gdal'
            stop 'OpenRaster - CreateRaster - ERR01'
        endif
        
        !delete filename if exists - does it work?
        !inquire(file=FileName, exist=exists)
        !if (exists) then
        !    call system("del "//trim(FileName))
        !endif
        
        ds = gdalcreate(driver, trim(adjustl(FileName))//char(0), &
                        RasterWidth, RasterHeight, NumberBands, DataType, c_ptr_ptr_getobject(c_ptr_ptr_new((/('',i=1,0)/))))
        if (.NOT.gdalassociated(ds)) then
            write(*,*) 'Error creating dataset on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR02'
        endif
        
        !projection
        ierr = gdalsetprojection(ds, trim(adjustl(Projection)))
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL projection on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR03'
        end if
        
        !geotransformation
        ierr = gdalsetgeotransform(ds, GeoTrans)
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL geotransformation properties on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR04'
        end if
        
        !x_dim = size(DataMatrix, dim=1)
        !y_dim = size(DataMatrix, dim=2)
        !allocate(DataMatrix3D(x_dim, y_dim, 1))
        
        !DataMatrix3D = reshape(DataMatrix, [x_dim, y_dim, 1])
        
        ierr = gdaldatasetrasterio_f(ds, GF_Write, 0, 0, DataMatrix3D)
        if (ierr /= 0) then
            write(*,*) 'Error setting matrix on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR05'
        end if
        
        !deallocate(DataMatrix3D)
        call gdalclose(ds)
    
    end subroutine CreateRasterR4
    !----------------------------------------------------------------------
    
    !----------------------------------------------------------------------
    subroutine CreateRasterR8(FileName, DriverName, RasterWidth, RasterHeight, DataType, NumberBands, Projection, GeoTrans, DataMatrix3D)
        !Arguments------------------------------------------------
        character(len=*),          intent(IN)    :: FileName
        character(len=*),          intent(IN)    :: DriverName
        integer,                   intent(IN)    :: RasterWidth
        integer,                   intent(IN)    :: RasterHeight
        integer,                   intent(IN)    :: DataType
        integer,                   intent(IN)    :: NumberBands
        character(len=*),          intent(IN)    :: Projection
        real(kind=c_double),       intent(IN)    :: GeoTrans(6)
        real(8), dimension(:,:,:), intent(INOUT)    :: DataMatrix3D
        !Local----------------------------------------------------
        type(gdaldriverh)                        :: driver
        type(gdaldataseth)                       :: ds
        integer                                  :: i, ierr
        logical                                  :: exists
        !Body-----------------------------------------------------    
    
        call gdalallregister()
        driver = gdalgetdriverbyname(DriverName//char(0))
        if (.NOT.gdalassociated(driver)) then
            write(*,*) 'Error getting GeoTIFF driver from gdal'
            stop 'OpenRaster - CreateRaster - ERR01'
        endif
        
        !delete filename if exists - does it work?
        !inquire(file=FileName, exist=exists)
        !if (exists) then
        !    call system("del "//trim(FileName))
        !endif
        
        ds = gdalcreate(driver, trim(adjustl(FileName))//char(0), &
                        RasterWidth, RasterHeight, NumberBands, DataType, c_ptr_ptr_getobject(c_ptr_ptr_new((/('',i=1,0)/))))
        if (.NOT.gdalassociated(ds)) then
            write(*,*) 'Error creating dataset on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR02'
        endif
        
        !projection
        ierr = gdalsetprojection(ds, trim(adjustl(Projection)))
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL projection on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR03'
        end if
        
        !geotransformation
        ierr = gdalsetgeotransform(ds, GeoTrans)
        if (ierr /= 0) then
            write(*,*) 'Error setting GDAL geotransformation properties on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR04'
        end if
                
        ierr = gdaldatasetrasterio_f(ds, GF_Write, 0, 0, DataMatrix3D)
        if (ierr /= 0) then
            write(*,*) 'Error setting matrix on file ',trim(adjustl(FileName))
            stop 'OpenRaster - CreateRaster - ERR05'
        end if
        
        !deallocate(DataMatrix3D)
        call gdalclose(ds)
    
    end subroutine CreateRasterR8
    !----------------------------------------------------------------------
    
    !----------------------------------------------------------------------
    !subroutine AverageRasters(file1, file2, file3)
    !    !Arguments------------------------------------------------
    !    character(len=PathLength)           :: file1
    !    character(len=PathLength)           :: file2
    !    character(len=PathLength)           :: file3
    !    !Local----------------------------------------------------   
    !    type(T_Raster), pointer             :: Raster1
    !    type(T_Raster), pointer             :: Raster2
    !    type(T_Raster), pointer             :: Raster3
    !    type(T_Raster), pointer             :: OutRaster
    !    real(4), dimension(:,:), pointer    :: OutMatrix
    !    character(len=8)                    :: current_date
    !    character(len=10)                   :: current_time
    !    !Body----------------------------------------------------- 
    !    
    !    call OpenRaster(file1, 'GTiff', Raster1)
    !    call OpenRaster(file2, 'GTiff', Raster2)
    !    call OpenRaster(file3, 'GTiff', Raster3)
    !    
    !    call DATE_AND_TIME(current_date, current_time)
    !    write(*,*) current_date(1:4)//"/"//current_date(5:6)//"/"//current_date(7:8)//" "//current_time(1:2)//":"//current_time(3:4)//":"//current_time(5:)//" : Performing Average"
    !    
    !    allocate(OutRaster)
    !    allocate(OutMatrix(1:Raster1%ii, 1:Raster1%jj))
    !    
    !    OutMatrix = (Raster1%Bands(1)%DataMatrix%DataR4 + Raster2%Bands(1)%DataMatrix%DataR4 + Raster3%Bands(1)%DataMatrix%DataR4)/3
    !    
    !    call DATE_AND_TIME(current_date, current_time)
    !    write(*,*) current_date(1:4)//"/"//current_date(5:6)//"/"//current_date(7:8)//" "//current_time(1:2)//":"//current_time(3:4)//":"//current_time(5:)//" : Done."
    !    
    !    call CreateRaster('out.tiff', 'GTiff',                    &
    !                      RasterWidth = Raster1%ii,               &
    !                      RasterHeight = Raster1%jj,              &
    !                      DataType = Raster1%Bands(1)%DataType,   &
    !                      NumberBands = size(Raster1%Bands),      &
    !                      Projection = Raster1%ProjRef,           &
    !                      GeoTrans = Raster1%GeoTransform,        &
    !                      DataMatrix = OutMatrix)
    !
    !end subroutine AverageRasters
    !----------------------------------------------------------------------
    
end module RasterFunctions