!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Shell
! PROGRAM       : MainShell
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig /Luis Fernandes - v4.0
! DESCRIPTION   : Shell to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program CreateHDF5BoxesFluxes

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleTimeSerie
    use ModuleHDF5
    use HDF5

    implicit none

    !private :: ConstructCreateHDF5BoxesFluxes
    !private ::     ReadKeyWords
    !private ::     ReadBxm         ! Water volume, Salinity(time)
    !private ::     ReadBxf         ! Water volume flux(time)
    !private ::     ReadBoxesHdf5   ! Boxes numbers
    !private :: ModifyCreateHDF5BoxesFluxes
    !private ::     CalculateBoxesBarycenter
    !private ::     CalculateBoxesBarycenterDistances
    !private ::     CalculateBoxesBarycenterAngles
    !private ::     CalculateBoxesFluxesComponents
    !private ::     WriteBoxesFluxesHdf5
    !private :: KillCreateHDF5BoxesFluxes

    interface WriteGridVar
        module procedure WriteGridVar2D_Real
        module procedure WriteGridVar2D_Integer
        module procedure WriteGridVar3D_Real
        module procedure WriteGridVar3D_Integer
    end interface WriteGridVar

    interface WriteVar
        module procedure WriteVar2D_Real
        module procedure WriteVar2D_Integer
        module procedure WriteVar3D_Real
        module procedure WriteVar3D_Integer
    end interface WriteVar

    interface Killocate
        module procedure Killocate1D_real
        module procedure Killocate1D_real8
        module procedure Killocate1D_integer
        module procedure Killocate2D_real
        module procedure Killocate2D_integer
        module procedure Killocate3D_real
        module procedure Killocate3D_integer
    end interface Killocate

    type T_HDFGrid

        real, dimension(:,:), pointer       :: Lat, Lon, Lat_Stag, Lon_Stag
        real, dimension(:,:), pointer       :: Bathymetry
        integer, dimension(:,:), pointer    :: WaterPoints2D
        integer, dimension(:,:,:), pointer  :: WaterPoints3D
        real, dimension(:,:,:), pointer     :: Vert3D
        real, dimension(:,:), pointer       :: Float2D
        real, dimension(:,:,:), pointer     :: Float3D
        integer, dimension(:,:), pointer    :: Integer2D
        integer, dimension(:,:,:), pointer  :: Integer3D

    end type T_HDFGrid

    type T_HDFFile

        integer                             :: FileID
        character(len=PathLength)           :: Name
        integer                             :: ObjHDF5      = 0
        integer                             :: nInstants
        real(8), dimension(:), pointer      :: Times
        type(T_Time)                        :: InitialDate
        type(T_Size3D)                      :: Size
        character(len=StringLength)         :: SizeGroup, SizeDataSet
        character(len=StringLength)         :: TimeVar
        logical                             :: ReadLatLon   = .true.
        logical                             :: Sigma   = .false.
        type(T_HDFGrid)                     :: Grid

    end type T_HDFFile

    type T_HDFBox

        type(T_HDFFile)                         :: HDFFile
        integer, dimension(:,:), pointer        :: Box2D
        integer, dimension(:,:,:), pointer      :: Box3D

    end type T_HDFBox

    type T_2Dsection
        
        real, dimension(:,:,:), pointer         :: U
        real, dimension(:,:,:), pointer         :: V

    end type T_2Dsection

    type T_Box
        
        integer                                 :: ID
        character(len=PathLength)               :: File
        character(len=line_length)              :: Header
        real, dimension(:,:), pointer           :: DataMatrix
        integer                                 :: DataColumns
        integer                                 :: DataRows
        type(T_Time)                            :: InitialData
        integer                                 :: Dt
        real, dimension(:), pointer             :: BoxData
        real, dimension(:,:), pointer           :: FluxData

    end type T_Box

    type T_FluxIndexes
        
        character(len=PathLength)               :: File
        integer                                 :: Length
        integer, dimension(:), pointer          :: Columns
        integer, dimension(:), pointer          :: FFrom
        integer, dimension(:), pointer          :: TTo

    end type T_FluxIndexes

    type T_BoxesFluxes

        character(len=PathLength)   :: HDF5BoxesFile
        type(T_HDFBox)              :: InFile
        type(T_HDFBox)              :: OuFile
        type(T_Box)                 :: BoxesWaterFlux
        type(T_Box)                 :: BoxesWater
        type(T_Box)                 :: BoxesSalinity
        type(T_FluxIndexes)         :: FluxIndexes
        real, dimension(:,:,:), pointer  :: Salinity3D
        type (T_2Dsection)                  :: XY
        type (T_2Dsection)                  :: XZ
        type (T_2Dsection)                  :: YZ
        integer                     :: BoxesNumber
        integer                     :: ObjEnterData = 0
        integer                     :: ObjTime      = 0

    end type T_BoxesFluxes

    integer                     :: STAT_CALL

    type (T_BoxesFluxes)        :: Me

    call ConstructCreateHDF5BoxesFluxes
    call ModifyCreateHDF5BoxesFluxes
    call KillCreateHDF5BoxesFluxes

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructCreateHDF5BoxesFluxes

        call ReadKeywords

        !Gets the boxes distribution grid-data in hdf5 format
        call ReadHDF5Box(Me%InFile)
        
        !Gets the pair of boxes number representing the fluxes
        call ReadBoxFluxIndexes(Me%FluxIndexes)

        !Reads the boxes time series for the water mass
        call LoadBox(Me%BoxesWater)

        !Reads the boxes time series for the salinity
        call LoadBox(Me%BoxesSalinity)

        !Set the available number of boxes
        Me%BoxesNumber = Me%BoxesWater%Datacolumns - 8

        !Reads the box fluxes time series for the water mass
        call LoadFlux(Me%BoxesWaterFlux, Me%FluxIndexes, Me%BoxesNumber)
        
    end subroutine ConstructCreateHDF5BoxesFluxes
   
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local-------------------------------------------------
        character(PathLength)                       :: DataFile
        integer                                     :: FromFile
        integer                                     :: iflag
        logical                                     :: exist

        call ReadFileName('IN_MODEL', DataFile, "CreateHDF5BoxesFluxes", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - CreateHDF5BoxesFluxes - ERR01'

        call ConstructEnterData (Me%ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - CreateHDF5BoxesFluxes - ERR02'

        call GetExtractType     (FromFile = FromFile)

        !Read the boxes numbers file
        call GetData(Me%InFile%HDFFile%Name,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'HDF5BOXES_INFILE',                                      &
                     default      = 'Boxes.hdf5',                                           &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR40'

        !Verifies if file exists
        inquire(FILE = trim(Me%InFile%HDFFile%Name), EXIST = exist)
        if (.not. exist) then
            write(*,*)'Water HDF5 file does not exist'
            stop 'ReadKeywords - ReadWriteHDF5 - ERR02'
        endif

        call GetData(Me%OuFile%HDFFile%Name,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'HDF5BOXES_OUFILE',                                      &
                     default      = 'OutBoxes.hdf5',                                           &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR40'

        call GetData(Me%InFile%HDFFile%TimeVar,                                                 &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_TIME_VAR',                                     &
                     Default      = 'Time',                                             &
                     ClientModule = 'ReadWriteHDF5',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR12'


        call GetData(Me%InFile%HDFFile%SizeGroup,                                               &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_SIZE_GROUP',                                   &
                     Default      = '/Grid',                                            &
                     ClientModule = 'ReadWriteHDF5',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR13'

        call GetData(Me%InFile%HDFFile%SizeDataSet,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF_SIZE_DATASET',                                 &
                     Default      = 'WaterPoints3D',                                    &
                     ClientModule = 'ReadWriteHDF5',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR14'

        !Read the boxes fluxes indexes file
        call GetData(Me%FluxIndexes%File,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'INDEXESFILE',                                      &
                     default      = 'BoxesFluxesIndexes.txt',                                           &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR41'

        !Read boxes fluxes indexes file lines size
        call GetData(Me%FluxIndexes%Length,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'INDEXESLENGTH',                                      &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR41'

        !Read the volume fluxes file
        call GetData(Me%BoxesWaterFlux%File,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'BOXESWATERFLUXFILE',                                      &
                     default      = 'water.bxf',                                           &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR50'

        !Read the boxes volume file
        call GetData(Me%BoxesWater%File,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'BOXESWATERFILE',                                      &
                     default      = 'water.bxm',                                           &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR50'

        !Read the boxes salinity
        call GetData(Me%BoxesSalinity%File,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'BOXESSALINITYFILE',                                      &
                     default      = 'salinity.bxm',                                           &
                     ClientModule = 'MainCreateHDF5BoxesFluxes',                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MainCreateHDF5BoxesFluxes - ERR50'

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - CreateHDF5BoxesFluxes - ERR03'

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine LoadBox(Box)
        
        !Arguments-------------------------------------------------------
        type(T_Box)                             :: Box

        call ReadBox(Box)
        call FillBox(Box)

    end subroutine LoadBox
    
    !--------------------------------------------------------------------------

    subroutine LoadFlux(Flux, Indexes, Datacolumns)
        
        !Arguments-------------------------------------------------------
        type(T_Box)                             :: Flux
        type(T_FluxIndexes)                     :: Indexes
        integer                                 :: Datacolumns

        call ReadBox(Flux)
        !The number of boxes is the number of columns in water.bxm  minus 
        !eight columns (due to the date).
        call FillFlux(Flux, Indexes, Datacolumns)

    end subroutine LoadFlux
    
    !--------------------------------------------------------------------------

    subroutine ReadBoxFluxIndexes(FluxIndexes)
        
        !Arguments------------------------------------------------
        type(T_FluxIndexes)                         :: FluxIndexes

        !Local-----------------------------------------------
        integer                                     :: Length
        integer                                     :: i

        Length = Me%FluxIndexes%Length

        allocate(Me%FluxIndexes%Columns(1:Length))
        allocate(Me%FluxIndexes%FFrom(1:Length))
        allocate(Me%FluxIndexes%TTo(1:Length))

        open(20, FILE=trim(FluxIndexes%File), STATUS='OLD')

        do i = 1,Length

            READ (20,*) FluxIndexes%Columns(i), FluxIndexes%FFrom(i), FluxIndexes%TTo(i)

        enddo

        close(20)        

    end subroutine ReadBoxFluxIndexes
    
    !--------------------------------------------------------------------------

    subroutine ReadHDF5Box(Box)
        
        !Arguments----------------------------------------------
        type(T_HDFBox)                             :: Box

        call OpenHDF5_File(Box%HDFFile)

        call ReadGrid(Box%HDFFile)

        call ReadBoxes(Box)

    end subroutine ReadHDF5Box
    
    !--------------------------------------------------------------------------
   
    subroutine OpenHDF5_File(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Local-----------------------------------------------------------------
        integer                             :: HDF5_READ

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (HDF_IO_File%ObjHDF5, HDF_IO_File%Name, Access = HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5_IO_File - ReadWriteHDF5 - ERR01'

        write(*,*)'Opened hdf5 file                : ', trim(HDF_IO_File%Name)


    end subroutine OpenHDF5_File

    !--------------------------------------------------------------------------
   
    subroutine OpenHDF5_FileW(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Local-----------------------------------------------------------------
        integer                             :: HDF5_CREATE

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        call ConstructHDF5 (HDF_IO_File%ObjHDF5, HDF_IO_File%Name, Access = HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5_FileW - ReadWriteHDF5 - ERR01'

        write(*,*)'Opened hdf5 file                : ', trim(HDF_IO_File%Name)


    end subroutine OpenHDF5_FileW

    !--------------------------------------------------------------------------

    subroutine ReadBoxes(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFBox)                     :: HDF_IO_File

        !Begin-----------------------------------------------------------------
       
        write(*,*)
        write(*,*)"Reading boxes..."
        write(*,*)

        allocate(HDF_IO_File%Box2D (1:HDF_IO_File%HDFFile%Size%IUB, &
                                1:HDF_IO_File%HDFFile%Size%JUB))
        allocate(HDF_IO_File%Box3D (1:HDF_IO_File%HDFFile%Size%IUB, &
                                1:HDF_IO_File%HDFFile%Size%JUB, &
                                1:HDF_IO_File%HDFFile%Size%KUB))

        call HDF5SetLimits(HDF_IO_File%HDFFile%ObjHDF5, ILB = 1, IUB = HDF_IO_File%HDFFile%Size%IUB, &
                               JLB = 1, JUB = HDF_IO_File%HDFFile%Size%JUB, &
                               KLB = 1, KUB = HDF_IO_File%HDFFile%Size%KUB, &
                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBoxes - ReadWriteHDF5 - ERR01'

        call HDF5ReadData(HDF5ID       = HDF_IO_File%HDFFile%ObjHDF5,            &
                          GroupName    = "/Boxes",                       &
                          Name         = "Boxes2D",                  &
                          Array2D      = HDF_IO_File%Box2D,                      &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBoxes - ReadWriteHDF5 - ERR01'

        call HDF5ReadData(HDF5ID       = HDF_IO_File%HDFFile%ObjHDF5,            &
                          GroupName    = "/Boxes",                       &
                          Name         = "Boxes3D",                 &
                          Array3D      = HDF_IO_File%Box3D,                        &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBoxes - ReadWriteHDF5 - ERR02'

    end subroutine ReadBoxes

    !--------------------------------------------------------------------------

    subroutine ReadGrid(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Begin-----------------------------------------------------------------
       
        write(*,*)
        write(*,*)"Reading grid..."
        write(*,*)

        call GetHDF5FileID (HDF_IO_File%ObjHDF5, HDF_IO_File%FileID, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteGrid - ReadWriteHDF5 - ERR01'

        !Sets up the Float2D and Float3D arrays dimensions
        call ReadSetDimensions(HDF_IO_FILE)

        call ReadLatLon(HDF_IO_FILE)
        call ReadBathymetry(HDF_IO_FILE)

        if (HDF_IO_File%Size%KUB .gt. 0) then

            call ReadVertical(HDF_IO_File)
    
        endif

    end subroutine ReadGrid

    !--------------------------------------------------------------------------

    subroutine ReadVertical(HDF_IO_File)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Local-----------------------------------------------------------------
        integer                             :: i, j, k
        real                                :: Vert, area
        character(len=StringLength)         :: VertVar
        character(len=StringLength)         :: PointsVar
        integer                             :: CurrentInstant

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading Vertical Coordinate..."
        
        !Same as SZZ
        allocate(HDF_IO_File%Grid%Vert3D (1:HDF_IO_File%Size%IUB, &
                              1:HDF_IO_File%Size%JUB, &
                              1:HDF_IO_File%Size%KUB+1))

        allocate(HDF_IO_File%Grid%WaterPoints3D (1:HDF_IO_File%Size%IUB, &
                                1:HDF_IO_File%Size%JUB, &
                                1:HDF_IO_File%Size%KUB))

        !ReadWrite VerticalZ code

        VertVar = "Vertical"

        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                               KLB = 1, KUB = HDF_IO_File%Size%KUB+1, &
                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR01'

        call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid/VerticalZ",                       &
                          Name         = trim(VertVar),                 &
                          Array3D      = HDF_IO_File%Grid%Vert3D,                        &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR02'

        !Loop through instants
        do CurrentInstant = 1, HDF_IO_File%nInstants

            !Read VerticalZ
            call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB+1, &
                                               STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR01'
        
            call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid/VerticalZ",                       &
                          Name         = trim(VertVar),                 &
                          Array3D      = HDF_IO_File%Grid%Vert3D,                        &
                          OutputNumber = CurrentInstant,                             &
                          STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR02'

        enddo
        
        !Read Waterpoints code
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR03'

        PointsVar = "WaterPoints3D"

        call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(PointsVar),               &
                          Array3D      = HDF_IO_File%Grid%WaterPoints3D,        &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR04'

        !Read Openpoints code
        PointsVar = "OpenPoints"

        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR03'

        !Loop through instants    
        
        do CurrentInstant = 1, HDF_IO_File%nInstants

                call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR03'
    
                call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid/"//trim(PointsVar),                       &
                          Name         = trim(PointsVar),               &
                          Array3D      = HDF_IO_File%Grid%WaterPoints3D,                 &
                          OutputNumber = CurrentInstant,                             &
                          STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteVertical - ReadWriteHDF5 - ERR04'

        enddo

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadVertical

    !--------------------------------------------------------------------------

    subroutine ReadLatLon (HDF_IO_File)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        character(len=StringLength)         :: LatVar, LonVar

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading latitude and longitude..."
        

        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB+1, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB+1, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - ReadWriteHDF5 - ERR01'

        if(HDF_IO_File%ReadLatLon)then
            LatVar = "Latitude"
            LonVar = "Longitude"
        else
            LatVar = "ConnectionY"
            LonVar = "ConnectionX"
        end if

        !Read
        call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(LatVar),                  &
                          Array2D      = HDF_IO_File%Grid%Lat_Stag,                      &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - ReadWriteHDF5 - ERR01'
        
        call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(LonVar),                  &
                          Array2D      = HDF_IO_File%Grid%Lon_Stag,                      &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteLatLon - ReadWriteHDF5 - ERR01'

        do i = HDF_IO_File%Size%ILB, HDF_IO_File%Size%IUB
        do j = HDF_IO_File%Size%JLB, HDF_IO_File%Size%JUB

            HDF_IO_File%Grid%Lon(i,j) = ( HDF_IO_File%Grid%Lon_Stag(i,j) + HDF_IO_File%Grid%Lon_Stag(i+1 ,j) ) * 0.5
            HDF_IO_File%Grid%Lat(i,j) = ( HDF_IO_File%Grid%Lat_Stag(i,j) + HDF_IO_File%Grid%Lat_Stag(i,j+1) ) * 0.5

        enddo
        enddo

        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadLatLon

    !--------------------------------------------------------------------------

    subroutine ReadBathymetry (HDF_IO_File)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Local------------------------------------------------------------------------
        character(len=StringLength)         :: NCDFName, LongName, StandardName, Units
        real                                :: MinValue, MaxValue, ValidMin, ValidMax

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading bathymetry..."
        
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteBathymetry - ReadWriteHDF5 - ERR01'

        call HDF5ReadData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = "Bathymetry",                  &
                          Array2D      = HDF_IO_File%Grid%Float2D,                    &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteBathymetry - ReadWriteHDF5 - ERR01'
        
       
        write(*,*)"Done!"
        write(*,*)

    end subroutine ReadBathymetry

    !--------------------------------------------------------------------------

    subroutine ReadSetDimensions(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: gr_id, dset_id, class_id, size
        integer(HID_T)                              :: space_id, datatype_id
        integer(HID_T)                              :: rank
        integer(HSIZE_T), dimension(7)              :: dims, maxdims

        !Begin-----------------------------------------------------------------

        write(*,*)"Reading sizes..."

        call h5gopen_f(HDF_IO_FILE%FileID, HDF_IO_FILE%SizeGroup, gr_id, STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadSetDimensions - ReadWriteHDF5 - ERR01'

        !Opens data set
        call h5dopen_f(gr_id, trim(adjustl(HDF_IO_FILE%SizeDataSet)), dset_id, STAT_CALL)
        call h5dget_space_f                 (dset_id, space_id, STAT_CALL)
        call h5sget_simple_extent_ndims_f   (space_id, rank, STAT_CALL)
        call h5dget_type_f                  (dset_id, datatype_id,   STAT_CALL)
        call h5tget_size_f                  (datatype_id, size,      STAT_CALL)
        call h5tget_class_f                 (datatype_id, class_id,  STAT_CALL) 
        call h5tclose_f                     (datatype_id, STAT_CALL) 
        call h5sget_simple_extent_dims_f    (space_id, dims, maxdims, STAT_CALL) 
        call h5sclose_f                     (space_id, STAT_CALL)
        call h5dclose_f                     (dset_id, STAT_CALL)
        call h5gclose_f                     (gr_id, STAT_CALL)

        allocate(HDF_IO_File%Grid%WaterPoints2D        (1:dims(1), 1:dims(2)))
        allocate(HDF_IO_File%Grid%WaterPoints3D        (1:dims(1), 1:dims(2), 1:dims(3)))

        HDF_IO_File%Grid%WaterPoints2D    = null_int
        HDF_IO_File%Grid%WaterPoints3D    = null_int  

        allocate(HDF_IO_File%Grid%Lat        (1:dims(1), 1:dims(2)))
        allocate(HDF_IO_File%Grid%Lon        (1:dims(1), 1:dims(2)))
        allocate(HDF_IO_File%Grid%Lat_Stag   (1:dims(1)+1, 1:dims(2)+1))
        allocate(HDF_IO_File%Grid%Lon_Stag   (1:dims(1)+1, 1:dims(2)+1))

        HDF_IO_File%Grid%Lat    = null_int
        HDF_IO_File%Grid%Lon    = null_int  
        HDF_IO_File%Grid%Lat_Stag  = null_real
        HDF_IO_File%Grid%Lon_Stag  = null_real  

        allocate(HDF_IO_File%Grid%Integer2D   (1:dims(1), 1:dims(2)))
        allocate(HDF_IO_File%Grid%Integer3D   (1:dims(1), 1:dims(2), 1:dims(3)))
        allocate(HDF_IO_File%Grid%Float2D (1:dims(1), 1:dims(2)))
        allocate(HDF_IO_File%Grid%Float3D (1:dims(1), 1:dims(2), 1:dims(3)))

        HDF_IO_File%Grid%Integer2D    = null_int
        HDF_IO_File%Grid%Integer3D    = null_int
        HDF_IO_File%Grid%Float2D  = null_real
        HDF_IO_File%Grid%Float3D  = null_real
       
        HDF_IO_FILE%Size%ILB = 1
        HDF_IO_FILE%Size%JLB = 1
        HDF_IO_FILE%Size%KLB = 1

        HDF_IO_FILE%Size%IUB = dims(1)
        HDF_IO_FILE%Size%JUB = dims(2)
        HDF_IO_FILE%Size%KUB = dims(3)

        write(*,*)
        write(*,*)"IUB", HDF_IO_FILE%Size%IUB
        write(*,*)"JUB", HDF_IO_FILE%Size%JUB
        write(*,*)"KUB", HDF_IO_FILE%Size%KUB
        write(*,*)

    end subroutine ReadSetDimensions

    !--------------------------------------------------------------------------

    subroutine ReadBox(Box)
        
        !Arguments----------------------------------------
        type(T_Box)                                 :: Box

        !Local-----------------------------------------------------------------

        !Open the timeserie file
        call StartTimeSerieInput(   TimeSerieID = Box%ID,                &
                                    TimeSerieDataFile = Box%File,          &
                                    ReadAll = .true.,                       &
                                    STAT = STAT_CALL                        &
                                    )
        if (STAT_CALL /= SUCCESS_) stop 'ReadBxm - MainCreateHDF5BoxesFluxes - ERR10'

        !Get the header for the names (optional)
        !call GetTimeSerieHeader(Box%ID, Box%Header, STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'ReadBxm - MainCreateHDF5BoxesFluxes - ERR20'

        !Get the data in a row/col matrix format
        call GetTimeSerieDataMatrix(Box%ID, Box%DataMatrix, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadBxm - MainCreateHDF5BoxesFluxes - ERR30'

        !Get number of columns
        call GetTimeSerieDataColumns(Box%ID, Box%DataColumns, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadBxm - MainCreateHDF5BoxesFluxes - ERR40'

        !Get number of rows
        call GetTimeSerieDataValues (Box%ID, Box%DataRows, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadBxm - MainCreateHDF5BoxesFluxes - ERR50'

    end subroutine ReadBox
    
    !--------------------------------------------------------------------------

    subroutine FillBox(Box)
        
        !Arguments----------------------------------------
        type(T_Box)                                 :: Box

        !Local------------------------------------------
        integer                                     :: i

        allocate( Box%BoxData(1:Box%Datacolumns) )

        do i = 1, Box%Datacolumns

            Box%BoxData(i) = Box%DataMatrix(Box%DataRows,i)

        enddo

    end subroutine FillBox

    !--------------------------------------------------------------------------

    subroutine FillFlux(Box, Indexes, Boxesnumber)
        
        !Arguments-------------------------------------------------------
        type(T_Box)                                 :: Box
        type(T_FluxIndexes)                         :: Indexes
        integer                                     :: Boxesnumber

        !Local-----------------------------------------------------------------
        integer                                     :: i, j

        allocate( Box%FluxData(1:Boxesnumber, 1:Boxesnumber) )

        Box%FluxData(:,:) = 0.0

        do i = 1, Indexes%Length           

            Box%FluxData(Indexes%FFrom(i),Indexes%TTo(i)) = Box%DataMatrix(Box%DataRows, Indexes%Columns(i)-1)

        enddo

    end subroutine FillFlux

    !--------------------------------------------------------------------------

    subroutine ModifyCreateHDF5BoxesFluxes

        !Local------------------------------------------------
        logical                                 :: Running        
        integer                                 :: a, b, i, j, k, x    
        Integer, dimension(:), pointer          :: CentreI, CentreJ, CentreK
        Integer, dimension(:), pointer          :: Lower, Upper
        real, dimension(:,:), pointer           :: DistU, DistV, DistK
        real, dimension(:,:), pointer           :: Cs, Sn
        real, dimension(:,:), pointer           :: FluxData
        real                                    :: volume
        type(T_Size3D)                          :: Size

        !Copy common data pointers from InFile to OuFile

        Me%OuFile%HDFFile%Grid = Me%InFile%HDFFile%Grid
        Me%OuFile%HDFFile%Size = Me%InFile%HDFFile%Size
        Me%OuFile%Box2D => Me%InFile%Box2D
        Me%OuFile%Box3D => Me%InFile%Box3D
        Size = Me%InFile%HDFFile%Size
        FluxData => Me%BoxesWaterFlux%FluxData

        !Let's calculate Vectors between boxes. These vectors are proportional
        !to the flux and are directed between the boxes geometric center.

        !Calculate centers

        allocate (CentreI(1:Me%BoxesNumber))
        allocate (CentreJ(1:Me%BoxesNumber))
        allocate (CentreK(1:Me%BoxesNumber))

        do a = 1, Me%BoxesNumber

            CentreI(a) = 0.0
            CentreJ(a) = 0.0
            CentreK(a) = 0.0

        enddo

        do a = 1, Me%BoxesNumber

            !Must start with a unitary volume.
            volume = 1.0

            do k = Size%KLB, Size%KUB
            do i = Size%ILB, Size%IUB
            do j = Size%JLB, Size%JUB

                if ( Me%InFile%Box3D(i,j,k) == a-1 ) then

                    CentreI(a) = CentreI(a) + i
                    CentreJ(a) = CentreJ(a) + j
                    CentreK(a) = Centrek(a) + k
                    volume = volume + 1.0
                    
                endif

            enddo
            enddo
            enddo

            CentreI(a) = CentreI(a) / volume
            CentreJ(a) = CentreJ(a) / volume
            CentreK(a) = CentreK(a) / volume

        enddo

        !Calculate distances
        
        allocate (DistU(1:Me%BoxesNumber, 1:Me%BoxesNumber))
        allocate (DistV(1:Me%BoxesNumber, 1:Me%BoxesNumber))
        allocate (DistK(1:Me%BoxesNumber, 1:Me%BoxesNumber))

        !a or b == 1 is the Zero box number. The zero box number doesn't interests us.
        do a = 1, Me%BoxesNumber
        do b = 1, Me%BoxesNumber

            if ( (CentreI(b) .ne. 0) .and. (CentreI(a) .ne. 0) ) then

                DistU(a,b) = Me%InFile%HDFFile%Grid%Lon(2,CentreJ(b)) - Me%InFile%HDFFile%Grid%Lon(2,CentreJ(a))

            else

                DistU(a,b) = 0.0

            endif

            if ( (CentreJ(b) .ne. 0) .and. (CentreJ(a) .ne. 0) ) then

                DistV(a,b) = Me%InFile%HDFFile%Grid%Lat(CentreI(b),2) - Me%InFile%HDFFile%Grid%Lat(CentreI(a),2)

            else

                DistV(a,b) = 0.0

            endif

            if ( (CentreK(b) .ne. 0) .and. (CentreK(a) .ne. 0) ) then
            
                DistK(a,b) = Me%InFile%HDFFile%Grid%Vert3D(2,2,CentreK(b)) - Me%InFile%HDFFile%Grid%Vert3D(2,2,CentreK(a))

            else

                DistK(a,b) = 0.0

            endif

        enddo
        enddo

        !Allocate variables
        allocate (Cs(1:Me%BoxesNumber, 1:Me%BoxesNumber))
        allocate (Sn(1:Me%BoxesNumber, 1:Me%BoxesNumber))            

        allocate( Lower(1:Me%BoxesNumber) )
        allocate( Upper(1:Me%BoxesNumber) )

        ! XY SURFACE -----------------------------------------------------

        !Calculate angles (Cos and Sin)
        call CalcAngles(Cs, Sn, DistU, DistV)

        !Calculate limits of each box along Z coordinate axis
        call CalcLimits(Lower, Upper, Size, 3)

        !Calculate fluxes
        call alloc2Dsection(Me%XY, Size)

        do a = 1, Me%BoxesNumber
        do b = 1, Me%BoxesNumber

            !If there's a non-null flow between the boxes
            if ( FluxData(a,b) /= 0.0 .and.         &
                .not. ( (b > a) .and. ( FluxData(a,b) == -FluxData(b,a) ) ) ) then

                !We also make sure that the boxes share a non-null vertical surface.
                if (    ( .not. Lower(a) >= Upper(b) ) .and.             &
                        ( .not. Lower(b) >= Upper(a) ) ) then

                    !If the flux is positive then it's driven from box a to box b
                    if ( FluxData(a,b) > 0.0 ) then

                        x = a

                    !If the flux is negative, then it's driven from box b to box a
                    else

                        x = b

                    endif

                    !Write down the fluxes, but first make sure that the cell isn't
                    !already occupied with another flux ...
                    do k = Lower(x), Upper(x)

                        if ( sqrt(Me%XY%U( CentreI(x), CentreJ(x), k)**2 &
                                + Me%XY%V( CentreI(x), CentreJ(x), k)**2 ) == 0.0 ) then
                                    
                                i = 0
                                j = 0

                        elseif ( sqrt(Me%XY%U( CentreI(x) + 1, CentreJ(x), k)**2  &
                                    + Me%XY%V( CentreI(x) + 1, CentreJ(x), k)**2 ) == 0.0 ) then
    
                                i = 1
                                j = 0

                        elseif ( sqrt(Me%XY%U( CentreI(x) - 1, CentreJ(x), k)**2  &
                                    + Me%XY%V( CentreI(x) - 1, CentreJ(x), k)**2 ) == 0.0 )then
    
                                i = -1
                                j = 0

                        elseif ( sqrt(Me%XY%U( CentreI(x), CentreJ(x) + 1, k)**2 &
                                    + Me%XY%V( CentreI(x), CentreJ(x) + 1, k)**2 ) == 0.0 )then
    
                                i = 0
                                j = 1

                        elseif   ( sqrt(Me%XY%U( CentreI(x), CentreJ(x) - 1, k)**2 &
                                      + Me%XY%V( CentreI(x), CentreJ(x) - 1, k)**2 ) == 0.0 ) then
    
                                i = 0
                                j = -1

                        endif

                        Me%XY%U( CentreI(x) + i, CentreJ(x) + j, k) = FluxData(a,b) * Cs(a,b)
                        Me%XY%V( CentreI(x) + i, CentreJ(x) + j, k) = FluxData(a,b) * Sn(a,b)

                    enddo

                endif

            endif

        enddo
        enddo

        ! YZ SURFACE -----------------------------------------------------

        !Calculate angles (Cos and Sin)
        call CalcAngles(Cs, Sn, DistU, DistK)

        !Calculate limits of each box along Y coordinate axis
        call CalcLimits(Lower, Upper, Size, 2)

        !Calculate fluxes
        call alloc2Dsection(Me%YZ, Size)

        do a = 1, Me%BoxesNumber
        do b = 1, Me%BoxesNumber

            !If there's a non-null flow between the boxes
            if ( FluxData(a,b) /= 0.0 .and.         &
                .not. ( (b > a) .and. ( FluxData(a,b) == -FluxData(b,a) ) ) ) then

                !We also make sure that the boxes share a non-null vertical surface.
                if (    ( .not. Lower(a) >= Upper(b) ) .and.             &
                        ( .not. Lower(b) >= Upper(a) ) ) then

                    !If the flux is positive then it's driven from box a to box b
                    if ( FluxData(a,b) > 0.0 ) then

                        x = a

                    !If the flux is negative, then it's driven from box b to box a
                    else

                        x = b

                    endif

                    do j = Lower(x), Upper(x)

                        if ( sqrt(Me%YZ%U( CentreI(x), j, CentreK(x))**2 &
                                + Me%YZ%V( CentreI(x), j, CentreK(x))**2 ) == 0.0 ) then
                                    
                                i = 0
                                k = 0

                        elseif ( sqrt(Me%YZ%U( CentreI(x) + 1, j, CentreK(x))**2  &
                                    + Me%YZ%V( CentreI(x) + 1, j, CentreK(x))**2 ) == 0.0 ) then
    
                                i = 1
                                k = 0

                        elseif ( sqrt(Me%YZ%U( CentreI(x) - 1, j, CentreK(x))**2  &
                                    + Me%YZ%V( CentreI(x) - 1, j, CentreK(x))**2 ) == 0.0 ) then
    
                                i = -1
                                k = 0

                        elseif ( sqrt(Me%YZ%U( CentreI(x), j, CentreK(x) + 1)**2 &
                                    + Me%YZ%V( CentreI(x), j, CentreK(x) + 1)**2 ) == 0.0 ) then
    
                                i = 0
                                k = 1

                        elseif   ( sqrt(Me%YZ%U( CentreI(x), j, CentreK(x) - 1)**2 &
                                      + Me%YZ%V( CentreI(x), j, CentreK(x) - 1)**2 ) == 0.0 ) then
    
                                i = 0
                                k = -1

                        endif

                        Me%YZ%U( CentreI(x) + i, j, CentreK(x) + k ) = FluxData(a,b) * Cs(a,b)
                        Me%YZ%V( CentreI(x) + i, j, CentreK(x) + k ) = FluxData(a,b) * Sn(a,b)

                    enddo

                endif

            endif

        enddo
        enddo

        ! XZ SURFACE -----------------------------------------------------

        !Calculate angles (Cos and Sin)
        call CalcAngles(Cs, Sn, DistV, DistK)

        !Calculate limits of each box along X coordinate axis
        call CalcLimits(Lower, Upper, Size, 1)

        !Calculate fluxes
        call alloc2Dsection(Me%XZ, Size)

        do a = 1, Me%BoxesNumber
        do b = 1, Me%BoxesNumber

            !If there's a non-null flow between the boxes
            if ( FluxData(a,b) /= 0.0 .and.         &
                .not. ( (b > a) .and. ( FluxData(a,b) == -FluxData(b,a) ) ) ) then

                !We also make sure that the boxes share a non-null vertical surface.
                if (    ( .not. Lower(a) >= Upper(b) ) .and.             &
                        ( .not. Lower(b) >= Upper(a) ) ) then

                    !If the flux is positive then it's driven from box a to box b
                    if ( FluxData(a,b) > 0.0 ) then

                        x = a

                    !If the flux is negative, then it's driven from box b to box a
                    else

                        x = b

                    endif

                    do i = Lower(x), Upper(x)

                        if ( sqrt(Me%XZ%U( i, CentreJ(x), CentreK(x))**2 &
                                + Me%XZ%V( i, CentreJ(x), CentreK(x))**2 ) == 0.0 ) then
                                    
                                j = 0
                                k = 0

                        elseif ( sqrt(Me%XZ%U( i, CentreJ(x) + 1, CentreK(x))**2  &
                                    + Me%XZ%V( i, CentreJ(x) + 1, CentreK(x))**2 ) == 0.0 ) then
    
                                j = 1
                                k = 0

                        elseif ( sqrt(Me%XZ%U( i, CentreJ(x) - 1, CentreK(x))**2  &
                                    + Me%XZ%V( i, CentreJ(x) - 1, CentreK(x))**2 ) == 0.0 ) then
    
                                j = -1
                                k = 0

                        elseif ( sqrt(Me%XZ%U( i, CentreJ(x), CentreK(x) + 1)**2 &
                                    + Me%XZ%V( i, CentreJ(x), CentreK(x) + 1)**2 ) == 0.0 ) then
    
                                j = 0
                                k = 1

                        elseif   (  sqrt(Me%XZ%U( i, CentreJ(x), CentreK(x) - 1)**2 &
                                       + Me%XZ%V( i, CentreJ(x), CentreK(x) - 1)**2 ) == 0.0 ) then
    
                                j = 0
                                k = -1

                        endif

                        Me%XZ%U( i, CentreJ(a) + j, CentreK(a) + k ) = FluxData(a,b) * Cs(a,b)
                        Me%XZ%V( i, CentreJ(a) + j, CentreK(a) + k ) = FluxData(a,b) * Sn(a,b)

                    enddo

                endif

            endif

        enddo
        enddo

        !Free memory
        call Killocate(CentreI)
        call Killocate(CentreJ)
        call Killocate(CentreK)

        call Killocate(DistU)
        call Killocate(DistV)
        call Killocate(DistK)

        call Killocate(Cs)
        call Killocate(Sn)

        call Killocate(Lower)
        call Killocate(Upper)

        !Fill salinity3D
        allocate ( Me%Salinity3D(Size%ILB:Size%IUB, Size%JLB:Size%JUB,  Size%KLB:Size%KUB) )

        do k = 1, Size%KUB
        do i = 1, Size%IUB
        do j = 1, Size%JUB

            do a = 1, Me%BoxesNumber
            
                if ( Me%InFile%Box3D(i,j,k) == a-1 ) then

                    Me%Salinity3D(i,j,k) = Me%BoxesSalinity%BoxData(a+8) / Me%BoxesWater%BoxData(a+8)

                endif

            enddo

        enddo
        enddo
        enddo

    end subroutine ModifyCreateHDF5BoxesFluxes
    
    !--------------------------------------------------------------------------

    subroutine CalcAngles(Cxy, Sxy, DistU, DistV)

        !Arguments--------------------------------------
        real, dimension(:,:), pointer           :: Cxy, Sxy
        real, dimension(:,:), pointer           :: DistU, DistV

        !Locals-----------------------------------------
        integer                                 :: a, b
        real                                    :: Hypot

        do a = 1, Me%BoxesNumber
        do b = 1, Me%BoxesNumber

            Hypot = sqrt( DistU(a,b)**2 + DistV(a,b)**2 )

            if (Hypot /= 0) then

                Cxy(a,b) = DistU(a,b) / Hypot
                Sxy(a,b) = DistV(a,b) / Hypot

            else

                Cxy(a,b) = 0.0
                Sxy(a,b) = 0.0

            endif
        
        enddo
        enddo

    end subroutine CalcAngles
    
    !--------------------------------------------------------------------------

    subroutine CalcLimits(Lower, Upper, Size, SizType)

        !Arguments----------------------------------------------------
        integer, dimension(:), pointer              :: Lower, Upper
        type(T_Size3D)                              :: Size
        integer                                     :: SizType    
        integer                                     :: a, i, j, k, x

        do a = 1, Me%BoxesNumber

            !Initialize with a very high number
            Lower(a) = 999999

            !Initialize with a lowest natural number
            Upper(a) = 0

            do k = 1, Size%KUB
            do i = 1, Size%IUB
            do j = 1, Size%JUB

                select case(SizType)

                    case(1)
                        x = i

                    case(2)
                        x = j

                    case(3)
                        x = k

                end select
                 
                if ( Me%InFile%Box3D(i,j,k) == a-1 ) then

                    if ( x < Lower(a) ) then
                        Lower(a) = x
                    endif

                    if ( x > Upper(a) ) then
                        Upper(a) = x
                    endif

                endif

            enddo
            enddo
            enddo

        enddo

    end subroutine CalcLimits

    !--------------------------------------------------------------------------

    subroutine alloc2Dsection(Surf, Size)

        !Arguments---------------------------------------------------
        type(T_2Dsection)               :: Surf
        type(T_Size3D)                  :: Size

        !Locals------------------------------------------------------
        integer                         :: i, j, k

        allocate ( Surf%U(Size%ILB:Size%IUB, Size%JLB:Size%JUB, Size%KLB:Size%KUB) )
        allocate ( Surf%V(Size%ILB:Size%IUB, Size%JLB:Size%JUB, Size%KLB:Size%KUB) )

        do k = Size%KLB, Size%KUB
        do i = Size%ILB, Size%IUB
        do j = Size%JLB, Size%JUB

            Surf%U(i, j, k) = 0.0
            Surf%V(i, j, k) = 0.0

        enddo
        enddo
        enddo

    end subroutine

    !--------------------------------------------------------------------------

    subroutine Kill2Dsection(Surf)

        !Arguments---------------------------------------------------
        type(T_2Dsection)               :: Surf

        call Killocate ( Surf%U )
        call Killocate ( Surf%V )

    end subroutine

    !--------------------------------------------------------------------------

    subroutine KillCreateHDF5BoxesFluxes

        !Local-----------------------------------------------
        integer                                     :: nUsers

        call WriteHDF5Boxes(Me%OuFile)
        
        call KillBox(Me%BoxesWater)
        call KillFluxIndexes(Me%FluxIndexes)

        !Common memory for InFile and OuFile
        call KillSetDimensions(Me%InFile%HDFFile)
        call KillHDF5Boxes(Me%InFile)

        call Kill2Dsection(Me%XY)
        call Kill2Dsection(Me%XZ)
        call Kill2Dsection(Me%YZ)

    end subroutine KillCreateHDF5BoxesFluxes
                

    !--------------------------------------------------------------------------

    subroutine WriteHDF5Boxes(Box)
        
        !Arguments----------------------------------------------
        type(T_HDFBox)                             :: Box

        call OpenHDF5_FileW(Box%HDFFile)

        call WriteGrid(Box%HDFFile)

        call WriteBoxes(Box)

        call WriteBoxesVectors(Box)

    end subroutine WriteHDF5Boxes
    
    !--------------------------------------------------------------------------

    subroutine WriteBoxes(HDF_Box)

        !Arguments-------------------------------------------------------------
        type(T_HDFBox)                     :: HDF_Box

        !Begin-----------------------------------------------------------------
       
        write(*,*)
        write(*,*)"Writing boxes..."
        write(*,*)

        call GetHDF5FileID (HDF_Box%HDFFile%ObjHDF5, HDF_Box%HDFFile%FileID, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteGrid - ReadWriteHDF5 - ERR01'

        !Write water Box
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Boxes2D", "m3 s-1", HDF_Box%Box2D)
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Boxes3D", "m3 s-1", HDF_Box%Box3D)

        !Write Salinity Box
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Salinity", "-", Me%Salinity3D)        

        write(*,*)
        write(*,*)"Done."
        write(*,*)

    end subroutine WriteBoxes

    !--------------------------------------------------------------------------

    subroutine WriteBoxesVectors(HDF_Box)

        !Arguments-------------------------------------------------------------
        type(T_HDFBox)                     :: HDF_Box

        !Begin------------------------------------------------------------
        write(*,*)
        write(*,*)"Writing boxes vectors..."
        write(*,*)

        call GetHDF5FileID (HDF_Box%HDFFile%ObjHDF5, HDF_Box%HDFFile%FileID, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteGrid - ReadWriteHDF5 - ERR01'

        !XY SURFACE ------------------------------------------------------
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Uxy", "m3 s-1", Me%XY%U)
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Vxy", "m3 s-1", Me%XY%V)

        !XZ SURFACE ------------------------------------------------------
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Uxz", "m3 s-1", Me%XZ%U)
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Vxz", "m3 s-1", Me%XZ%V)

        !YZ SURFACE ------------------------------------------------------
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Uyz", "m3 s-1", Me%YZ%U)
        call WriteVar(HDF_Box%HDFFile, "/Boxes", "Vyz", "m3 s-1", Me%YZ%V)

        write(*,*)
        write(*,*)"Done."
        write(*,*)

    end subroutine WriteBoxesVectors

    !--------------------------------------------------------------------------

    subroutine WriteVar2D_Real(HDF_IO_FILE, vargroup, varname, varunits, var2D)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: vargroup, varname, varunits
        real, dimension(:,:), pointer       :: var2D

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR10'

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,              &
                          GroupName    = trim(vargroup),                    &
                          Name         = trim(varname),                     &
                          Units        = trim(varunits),                    &
                          Array2D      = var2D,                             &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteVar2D_Real

    !--------------------------------------------------------------------------

    subroutine WriteVar3D_Real(HDF_IO_FILE, vargroup, varname, varunits, var3D, staggered)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: vargroup, varname, varunits
        real, dimension(:,:,:), pointer     :: var3D
        logical, optional                   :: staggered

        !Begin-----------------------------------------------------------------
        
        if (present(staggered)) then

            call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB+1, &
                                               STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR10'

        else

            call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR10'

        endif

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,              &
                          GroupName    = trim(vargroup),                    &
                          Name         = trim(varname),                     &
                          Units        = trim(varunits),                    &
                          Array3D      = var3D,                             &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteVar3D_Real

    !--------------------------------------------------------------------------

    subroutine WriteVar2D_Integer(HDF_IO_FILE, vargroup, varname, varunits, var2D)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: vargroup, varname, varunits
        integer, dimension(:,:), pointer       :: var2D

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR10'

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,              &
                          GroupName    = trim(vargroup),                    &
                          Name         = trim(varname),                     &
                          Units        = trim(varunits),                    &
                          Array2D      = var2D,                             &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteVar2D_Integer

    !--------------------------------------------------------------------------

    subroutine WriteVar3D_Integer(HDF_IO_FILE, vargroup, varname, varunits, var3D)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: vargroup, varname, varunits
        integer, dimension(:,:,:), pointer       :: var3D

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR10'

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,              &
                          GroupName    = trim(vargroup),                    &
                          Name         = trim(varname),                     &
                          Units        = trim(varunits),                    &
                          Array3D      = var3D,                             &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteVar - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteVar3D_Integer

    !--------------------------------------------------------------------------

    subroutine WriteGrid(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)"Writing grid..."
        write(*,*)

        call GetHDF5FileID (HDF_IO_File%ObjHDF5, HDF_IO_File%FileID, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadWriteGrid - ReadWriteHDF5 - ERR01'

        !Write lon
        call WriteGridVar(HDF_IO_FILE, "Longitude", "º", HDF_IO_FILE%Grid%Lon_Stag, .true.)

        !Write lat
        call WriteGridVar(HDF_IO_FILE, "Latitude", "º", HDF_IO_FILE%Grid%Lat_Stag, .true.)

        !Write Bath
        call WriteGridVar(HDF_IO_FILE, "Bathymetry", "m", HDF_IO_FILE%Grid%Float2D)

        !Write WaterPoints2D
        !call WriteGridVar(HDF_IO_FILE, "WaterPoints2D", "-", HDF_IO_FILE%Grid%WaterPoints2D)

        !Write WaterPoints3D
        call WriteGridVar(HDF_IO_FILE, "WaterPoints3D", "-", HDF_IO_FILE%Grid%WaterPoints3D)

        !Write VerticalZ
        call WriteVar(HDF_IO_FILE, "/Grid/Vertical","VerticalZ", "m", HDF_IO_FILE%Grid%Vert3D, .true.)

    end subroutine WriteGrid

    !--------------------------------------------------------------------------

    subroutine WriteGridVar2D_Real(HDF_IO_FILE, varname, varunits, var2D, staggered)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: varname, varunits
        real, dimension(:,:), pointer       :: var2D
        logical, optional                   :: staggered

        !Begin-----------------------------------------------------------------

        if (present(staggered)) then

            call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB+1, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB+1, &
                                               STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR10'
  
        else

            call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR10'

        endif

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(varname),                  &
                          Units        = trim(varunits),                        &
                          Array2D      = var2D,                    &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteGridVar2D_Real

    !--------------------------------------------------------------------------

    subroutine WriteGridVar3D_Real(HDF_IO_FILE, varname, varunits, var3D)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: varname, varunits
        real, dimension(:,:,:), pointer       :: var3D

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR10'

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(varname),                  &
                          Units        = trim(varunits),                        &
                          Array3D      = var3D,                    &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteGridVar3D_Real

    !--------------------------------------------------------------------------

    subroutine WriteGridVar2D_Integer(HDF_IO_FILE, varname, varunits, var2D)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: varname, varunits
        integer, dimension(:,:), pointer       :: var2D

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR10'

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(varname),                  &
                          Units        = trim(varunits),                        &
                          Array2D      = var2D,                    &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteGridVar2D_Integer

    !--------------------------------------------------------------------------

    subroutine WriteGridVar3D_Integer(HDF_IO_FILE, varname, varunits, var3D)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        character(len=stringlength)         :: varname, varunits
        integer, dimension(:,:,:), pointer       :: var3D

        !Begin-----------------------------------------------------------------
        call HDF5SetLimits(HDF_IO_File%ObjHDF5, ILB = 1, IUB = HDF_IO_File%Size%IUB, &
                                               JLB = 1, JUB = HDF_IO_File%Size%JUB, &
                                               KLB = 1, KUB = HDF_IO_File%Size%KUB, &
                                               STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR10'

        call HDF5WriteData(HDF5ID       = HDF_IO_File%ObjHDF5,            &
                          GroupName    = "/Grid",                       &
                          Name         = trim(varname),                  &
                          Units        = trim(varunits),                        &
                          Array3D      = var3D,                    &
                          STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'WriteGrid - ReadWriteHDF5 - ERR11'

        call HDF5FlushMemory (HDF_IO_File%ObjHDF5, STAT = STAT_CALL)                         
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteGrid - ModulePatchHDF5Files - ERR110'

    end subroutine WriteGridVar3D_Integer

    !--------------------------------------------------------------------------

    subroutine KillSetDimensions(HDF_IO_FILE)

        !Arguments-------------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File

        call Killocate(HDF_IO_File%Grid%WaterPoints2D)
        call Killocate(HDF_IO_File%Grid%WaterPoints2D)

        call Killocate(HDF_IO_File%Grid%Lat)
        call Killocate(HDF_IO_File%Grid%Lon)
        call Killocate(HDF_IO_File%Grid%Lat_Stag)
        call Killocate(HDF_IO_File%Grid%Lon_Stag)

        call Killocate(HDF_IO_File%Grid%Integer2D)
        call Killocate(HDF_IO_File%Grid%Integer3D)
        call Killocate(HDF_IO_File%Grid%Float2D)
        call Killocate(HDF_IO_File%Grid%Float3D)

    end subroutine KillSetDimensions
                
    !--------------------------------------------------------------------------

    subroutine KillHDF5Boxes(HDF_IO_File)

        !Arguments-------------------------------------------------------
        type(T_HDFBox)                     :: HDF_IO_File
        
        call Killocate(HDF_IO_File%Box2D)
        call Killocate(HDF_IO_File%Box3D)

        call KillHDF5File(HDF_IO_File%HDFFile)

    end subroutine KillHDF5Boxes

    !--------------------------------------------------------------------------

    subroutine KillHDF5File(HDF_IO_File)

        !Arguments-------------------------------------------------------
        type(T_HDFFile)                     :: HDF_IO_File
        
        call Killocate(HDF_IO_File%Times)         

    end subroutine KillHDF5File

    !--------------------------------------------------------------------------

    subroutine KillFluxIndexes(FluxIndexes)
        
        !Arguments-------------------------------------------------------
        type(T_FluxIndexes)                                :: FluxIndexes

        call Killocate(FluxIndexes%Columns)
        call Killocate(FluxIndexes%FFrom)
        call Killocate(FluxIndexes%TTo)

    end subroutine KillFluxIndexes

    !--------------------------------------------------------------------------

    subroutine KillBox(Box)
        
        !Arguments----------------------------------------
        type(T_Box)                                 :: Box

        call Killocate(Box%BoxData)

    end subroutine KillBox

    !--------------------------------------------------------------------------

    subroutine KillFlux(Box)
        
        !Arguments----------------------------------------
        type(T_Box)                                 :: Box

        call Killocate(Box%FluxData)

    end subroutine KillFlux

    !--------------------------------------------------------------------------

    subroutine Killocate1D_real(Var)
        
        !Arguments----------------------------------------
        real, dimension(:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate1D_real

    !--------------------------------------------------------------------------

    subroutine Killocate1D_real8(Var)
        
        !Arguments----------------------------------------
        real(8), dimension(:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate1D_real8

    !--------------------------------------------------------------------------

    subroutine Killocate2D_real(Var)
        
        !Arguments----------------------------------------
        real, dimension(:,:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate2D_real

    !--------------------------------------------------------------------------

    subroutine Killocate3D_real(Var)
        
        !Arguments----------------------------------------
        real, dimension(:,:,:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate3D_real

    !--------------------------------------------------------------------------

    subroutine Killocate1D_integer(Var)
        
        !Arguments----------------------------------------
        integer, dimension(:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate1D_integer

    !--------------------------------------------------------------------------

    subroutine Killocate2D_integer(Var)
        
        !Arguments----------------------------------------
        integer, dimension(:,:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate2D_integer

    !--------------------------------------------------------------------------

    subroutine Killocate3D_integer(Var)
        
        !Arguments----------------------------------------
        integer, dimension(:,:,:), pointer              :: Var

        if(associated(Var)) then
            deallocate(Var)
            nullify(Var)
        end if

    end subroutine Killocate3D_integer

end program CreateHDF5BoxesFluxes
