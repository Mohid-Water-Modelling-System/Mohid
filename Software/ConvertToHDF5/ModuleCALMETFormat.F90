!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : CALMETFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : December 2012
! REVISION      : Rosa Trancoso
! DESCRIPTION   : Module to convert CALMET ASCII files into HDF5 format.
!
!------------------------------------------------------------------------------
!DataFile
!
!   FILENAME                    : char              -           !Path to MM5 original file
!   TERRAIN_FILENAME            : char              -           !Path to MM5 TERRAIN file
!   OUTPUTFILENAME              : char              -           !Path to HDF5 file generated from MM5 file
!
!   START                       : YYYY MM DD HH MM SS   [-]     !Start date of new file
!   END                         : YYYY MM DD HH MM SS   [-]     !End date of new file

Module ModuleCALMETFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid    
    use proj4
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertCALMETFormat
    private ::      ReadOptions
    private ::      OpenAndReadTerrainFile
    private ::      WriteGridInformation 
    private ::      Open_HDF5_OutPut_File
    private ::      OpenAndReadCALMETFile
    private ::      BeginEndTime
    private ::      ComputeVerticalCoordinate
    private ::      ComputeWindSurface
    private ::      OutputFields
    private ::          WriteGridToHDF5File            
    private ::      KillCALMETFormat
    
    !Perfect gas constant in Pa m3 kg-1 K-1
    integer, parameter                                      :: Perfect_gas_R        = 287.04 !Pa m3 kg-1 K-1


    !Types---------------------------------------------------------------------
    type       T_Date
        type(T_Time)                                        :: Date
        type(T_Date), pointer                               :: Next
    end type  T_Date

    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        logical                                             :: Convert                  = .false.
        integer                                             :: nDimensions
        logical                                             :: StaggeredX               = .false.
        logical                                             :: StaggeredY               = .false.
        logical                                             :: StaggeredZ               = .false.        
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        integer, dimension(:,:  ),  pointer                 :: IValues2D
        logical                                             :: IsInteger                = .false.
        integer                                             :: OutputNumber             = 1
        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),              pointer                 :: Next
    end type  T_Field

 
    type       T_CALMETFormat
        integer                                             :: ObjEnterData             = 0
        integer                                             :: ObjHDF5                  = 0
        integer                                             :: ObjHorizontalGrid        = 0
        integer                                             :: Unit
        character(len=PathLength)                           :: FileName
        character(len=PathLength)                           :: TerrainFileName
        character(len=PathLength)                           :: OutputFileName
        character(len=PathLength)                           :: GridFileName    
        real, dimension(:,:  ),       pointer               :: Bathymetry
        real, dimension(:,:  ),       pointer               :: CenterX
        real, dimension(:,:  ),       pointer               :: CenterY
        real, dimension(:,:  ),       pointer               :: ConnectionX
        real, dimension(:,:  ),       pointer               :: ConnectionY
        real, dimension(:,:  ),       pointer               :: LandUse
        real, dimension(:    ),       pointer               :: ZFace
        
        
        logical                                             :: ConvertCALMET            = .false.
        logical                                             :: ConvertTERRAIN           = .false.

        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),                pointer               :: FirstField
        type(T_Date),                 pointer               :: FirstDate
        character(len=StringLength), dimension(:), pointer  :: FieldsToConvert
        character(len=StringLength), dimension(:), pointer  :: FieldsToRead
        logical                                             :: TimeWindow
        type(T_Time)                                        :: StartTime, EndTime
        real                                                :: OutputDTInterval = null_real
        
        integer                                             :: ProjType                 = null_int
        real                                                :: CenterLat                = null_real
        real                                                :: CenterLon                = null_real
        real                                                :: TrueLatUpper             = null_real
        real                                                :: TrueLatLower             = null_real
        real(8)                                             :: FalseEast                = null_real
        real(8)                                             :: FalseNorth               = null_real
        real                                                :: XOri, YOri               = null_real
        real                                                :: DX, DY

        type(prj90_projection)                              :: Proj
        character(len=20), dimension(8)                     :: Params

        

    end type  T_CALMETFormat

    type(T_CALMETFormat), pointer                              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertCALMETFormat(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions
        
        if (Me%ConvertTERRAIN) then
                call OpenAndReadTerrainFile
                call WriteGridInformation 
        endif
        
        if (Me%ConvertCALMET) then
            
            call Open_HDF5_OutPut_File        
            
            call OpenAndReadCALMETFile
            
            call WriteGridInformation ! terrain from CALMET file
        
            call BeginEndTime
            
            call ComputeVerticalCoordinate

            call ComputeWindSurface

            call OutputFields
            
        endif
        
        call KillCALMETFormat

        STAT = SUCCESS_


    end subroutine ConvertCALMETFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1, iflag2

        !Begin-----------------------------------------------------------------
       
        call GetData(Me%FileName,                                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'ModuleCALMETFormat',               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR01'
                
        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ModuleCALMETFormat',               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR02'

        call GetData(Me%GridFileName,                                   &
                     Me%ObjEnterData, iflag2,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ModuleCALMETFormat',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR03'

        
        if (iflag == 1) then

            Me%ConvertCALMET = .true.

            if (iflag1 == 0) then
                write(*,*)'Must specify OUTPUTFILENAME'
                stop 'ReadOptions - ModuleCALMETFormat - ERR04'            
            endif
            if (iflag2 == 0) then
                write(*,*)'Must specify OUTPUT_GRID_FILENAME'
                stop 'ReadOptions - ModuleCALMETFormat - ERR05'            
            endif
        endif
        
        !We may only want to covert TERREL ----------------
       
        !Terrain from TERREL
        call GetData(Me%TerrainFileName,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'TERRAIN_FILENAME',                 &
                     ClientModule = 'ModuleCALMETFormat',               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR06'
                
        if (iflag == 1) Me%ConvertTERRAIN = .true.

        write(*,*) 'Convert CALMET Output  = ', Me%ConvertCALMET
        write(*,*) 'Convert TERRAIN Output = ', Me%ConvertTERRAIN
        write(*,*) 
        
        call GetData(Me%StartTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ModuleCALMETFormat',               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR09'


        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ModuleCALMETFormat',               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR10'

        if (iflag==1 .AND. iflag1==1) Me%TimeWindow = .TRUE.
        
        if (Me%TimeWindow) then

            if (Me%StartTime .GE. Me%EndTime) then
                write (*,*)
                write (*,*) 'START greater or equal than END'
                stop 'ReadOptions - ModuleCALMETFormat - ERR11'
            endif

        endif

        call GetData(Me%OutputDTInterval,                               &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_DT',                        & 
                     Default      = 0.0,                                &
                     ClientModule = 'ModuleCALMETFormat',               &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCALMETFormat - ERR12'
        
        
    end subroutine ReadOptions
    
    !--------------------------------------------------------------------------
    
    subroutine OpenAndReadTerrainFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist
        
        !Others----------------------------------------------------------------
        integer                                     :: TerrainUnit        
        character(len=StringLength)                 :: AuxString
        character(len=StringLength)                 :: ProjName, Clat0, Clon0, Csp1, Csp2
        real                                        :: x1, y1
        integer                                     :: nx, ny, i, j
        integer                                     :: ILB,IUB, WILB, WIUB
        integer                                     :: JLB,JUB, WJLB, WJUB        
        real, pointer, dimension(:,:)               :: DataAux        
        
        !Begin-----------------------------------------------------------------
        

        !Verifies if file exists
        inquire(file = Me%TerrainFileName, exist = exist)
        if (.not. exist) then
            write(*,*)'TERRAIN file does not exist'
            stop 'OpenAndReadCALMETTerrainFile - ModuleCALMETFormat - ERR01'
        endif

        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading Terrain file...'
        write(*,*)

        call UnitsManager(TerrainUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenReadWriteTerrainFile - ModuleCALMETFormat - ERR02'

        open(Unit   = TerrainUnit,          &
             File   = Me%TerrainFileName,   &             
             STATUS = 'OLD',                &
             Action = 'READ',               &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenReadWriteTerrainFile - ModuleCALMETFormat - ERR03'

        rewind(TerrainUnit)

        read(TerrainUnit, *) AuxString   !TERREL.DAT      2.0             Header structure with coordinate parameters
        read(TerrainUnit, *) AuxString   !   2
        read(TerrainUnit, *) AuxString   !Produced by TERREL Version: 3.69  Level: 110330
        read(TerrainUnit, *) AuxString   !Internal Coordinate Transformations  ---  COORDLIB   Version: 1.99   Level: 070921
    
        read(TerrainUnit, *) ProjName                   !LCC     
        read(TerrainUnit, *) Clat0, Clon0, Csp1, Csp2   !38.78N          9.395W          30N             60N             
        read(TerrainUnit, *) Me%FalseEast, Me%FalseNorth  !   0.0000000       0.0000000    
        read(TerrainUnit, *)                            !WGS-84  02-21-2003  
        read(TerrainUnit, *) nx, ny, x1, y1, Me%DX, Me%DY     !      99      99       0.000       0.000       0.200       0.200
        read(TerrainUnit, *)                            !KM  M   
        read(TerrainUnit, *)                            !W_E N_S 

        !Allocate Size ------------------------------------
        
        Me%Size%ILB = 0             ; Me%WorkSize%ILB = Me%Size%ILB + 1
        Me%Size%IUB = nx + 1        ; Me%WorkSize%IUB = Me%Size%IUB - 1
        Me%Size%JLB = 0             ; Me%WorkSize%JLB = Me%Size%JLB + 1
        Me%Size%JUB = ny + 1        ; Me%WorkSize%JUB = Me%Size%JUB - 1

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        
        Me%DX = Me%DX * 1000.
        Me%DY = Me%DY * 1000.
        Me%FalseEast = Me%FalseEast * 1000.
        Me%FalseNorth = Me%FalseNorth * 1000.
        
        Me%XOri = x1
        Me%YOri = y1
                
        !Initialize Projection ----------------------------
        
        if (trim(adjustl(ProjName)) == 'LCC') then
            Me%ProjType = LAMB_CONF_CONIC_
    
            call ConvertLatLonCharToReal(Clat0, Me%CenterLat)
            call ConvertLatLonCharToReal(Clon0, Me%CenterLon)
            call ConvertLatLonCharToReal(Csp1, Me%TrueLatLower)
            call ConvertLatLonCharToReal(Csp2, Me%TrueLatUpper)
        

            Me%Params(1) = 'proj=lcc'
            Me%Params(2) = 'ellps=WGS84'
            write(Me%Params(3),'(a6,f8.3)') 'lat_1=', Me%TrueLatLower
            write(Me%Params(4),'(a6,f8.3)') 'lat_2=', Me%TrueLatUpper
            write(Me%Params(5),'(a6,f8.3)') 'lon_0=', Me%CenterLon
            write(Me%Params(6),'(a6,f8.3)') 'lat_0=', Me%CenterLat
            write(Me%Params(7),'(a4,f12.3)') 'x_0=' , Me%FalseEast
            write(Me%Params(8),'(a4,f12.3)') 'y_0=' , Me%FalseNorth
            
            
        else
            stop 'Not ready for projections other than LCC'
        endif


        write(*,*) Me%Params
        
        STAT_CALL = prj90_init(Me%Proj,Me%Params)
        call handle_proj_error(STAT_CALL)
        if (STAT_CALL /= PRJ90_NOERR) stop 'OpenReadWriteTerrainFile - ModuleCALMETFormat - ERR04'

        !Read terrain -------------------------------------
        
        allocate(DataAux(nx,ny))
            
        do j=ny,1,-1
            read(TerrainUnit, *) DataAux(:,j)
        enddo
                
        print*, 'DataAux(50,50) = ', DataAux(50,50)
        print*, 'DataAux(1,1)   = ', DataAux(1,1)
        print*, 'DataAux(99,99) = ', DataAux(99,99)
        print*, 'DataAux(99,98) = ', DataAux(99,98)
        
        allocate(Me%Bathymetry(ILB:IUB, JLB:JUB))               
        Me%Bathymetry(WILB:WIUB, WJLB:WJUB) = DataAux(WILB:WIUB, WJLB:WJUB)        

        do j = WJLB, WJUB
        do i = WILB, WIUB
            if (Me%Bathymetry(i,j) < -99.) Me%Bathymetry(i,j) = -99.
        enddo
        enddo

        
    end subroutine OpenAndReadTerrainFile
    
    !--------------------------------------------------------------------------
        
    subroutine WriteGridInformation

        !Local-----------------------------------------------------------------
        integer                                     :: ILB,IUB, WILB, WIUB
        integer                                     :: JLB,JUB, WJLB, WJUB   
        integer                                     :: i, j, STAT_CALL
        real(8)                                     :: x, y, lon, lat
        character(StringLength)                     :: COMENT1, COMENT2
        type(T_Size2D)                              :: Size2D
        
        nullify(Me%ConnectionX)
        nullify(Me%ConnectionY)
        nullify(Me%CenterX    )
        nullify(Me%CenterY    )
        
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB         
            
        allocate(Me%CenterX    (ILB:IUB, JLB:JUB))
        allocate(Me%CenterY    (ILB:IUB, JLB:JUB))
        allocate(Me%ConnectionX(ILB:IUB, JLB:JUB))
        allocate(Me%ConnectionY(ILB:IUB, JLB:JUB))
        
        do j = 1, WJUB + 1
        do i = 1, WIUB + 1
            
            x = dble(Me%XOri + (i-1)*Me%DX)
            y = dble(Me%YOri + (j-1)*Me%DY)
        
            STAT_CALL = prj90_inv(Me%Proj, x, y, lon, lat)
            call handle_proj_error(STAT_CALL)
            if (STAT_CALL /= PRJ90_NOERR) stop 'WriteGridInformation - ModuleCALMETFormat - ERR01'
        
            Me%ConnectionX(i,j) = lon
            Me%ConnectionY(i,j) = lat
            
        enddo
        enddo

        
        do j = 1, WJUB
        do i = 1, WIUB

            x = dble(Me%XOri + Me%DX/2. + (i-1)*Me%DX)
            y = dble(Me%YOri + Me%DY/2. + (j-1)*Me%DY)

            STAT_CALL = prj90_inv(Me%Proj, x, y, lon, lat)
            call handle_proj_error(STAT_CALL)
            if (STAT_CALL /= PRJ90_NOERR) stop 'WriteGridInformation - ModuleCALMETFormat - ERR02'

            Me%CenterX(i,j) = lon
            Me%CenterY(i,j) = lat

        enddo
        enddo
        
        
        write(*,*)
        write(*,*)'Constructing grid...'

        COMENT1 = 'Grid Data from '//trim(Me%GridFileName)
        write(COMENT2,'(A,F8.3,A,F8.3)') 'DX,DY =', Me%DX, ',', Me%DY

        Size2D%ILB = Me%WorkSize%ILB
        Size2D%IUB = Me%WorkSize%IUB
        Size2D%JLB = Me%WorkSize%JLB
        Size2D%JUB = Me%WorkSize%JUB
        
        call WriteGridData  (FileName          = Me%GridFileName,                              &
                                ConnectionX    = Me%ConnectionX,                               &
                                ConnectionY    = Me%ConnectionY,                               &
                                COMENT1        = COMENT1,                                      &
                                COMENT2        = COMENT2,                                      &
                                WorkSize       = Size2D,                                       &
                                CoordType      = SIMPLE_GEOG_,                                 &
                                Xorig          = Me%ConnectionX(1,1),                          &
                                Yorig          = Me%ConnectionY(1,1),                          &
                                Zone           = -99,                                          &
                                GRID_ANGLE     = 0.,                                           &
                                Latitude       = Me%CenterLat,                                 &
                                Longitude      = Me%CenterLon,                                 &
                                FillValue      = -99.,                                         &
                                Overwrite      = .true.,                                       &
                                GridData2D_Real= Me%Bathymetry,                                &
                                Datum          = WGS_84_DATUM,                                 &
                                ProjType       = Me%ProjType,                                  & 
                                SP1            = Me%TrueLatLower,                              &
                                SP2            = Me%TrueLatUpper,                              &
                                STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'WriteGridInformation - ModuleCALMETFormat - ERR06'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteGridInformation - ModuleCALMETFormat - ERR07'
        
    end subroutine WriteGridInformation               

    !----------------------------------------------------------------------
    
    subroutine ConvertLatLonCharToReal(AuxString, AuxReal)
    
        !Arguments ----------------------------------------
        character(len=StringLength), intent(in) :: AuxString
        real, intent(out)                       :: AuxReal
        !Local---------------------------------------------
        character(len=StringLength)             :: AuxString2
        integer                                 :: i, n
        logical                                 :: Found
        
        
        n = len(trim(adjustl(AuxString)))
        Found = .false.
        
        !Find North -----------------
        
        if (.not.Found) then            
            do i = 1, n
                if (AuxString(i:i) == 'N') then
                    AuxString2 = AuxString(1:(i-1))                
                    read(AuxString2,*) AuxReal
                    Found = .true.
                    exit
                endif
            enddo
        endif

        !Find South -----------------
        
        if (.not.Found) then            
            do i = 1, n
                if (AuxString(i:i) == 'S') then
                    AuxString2 = AuxString(1:(i-1))                
                    read(AuxString2,*) AuxReal
                    AuxReal = -1.*AuxReal
                    Found = .true.
                    exit
                endif
            enddo
        endif

        !Find East -----------------
        
        if (.not.Found) then            
            do i = 1, n
                if (AuxString(i:i) == 'E') then
                    AuxString2 = AuxString(1:(i-1))                
                    read(AuxString2,*) AuxReal
                    Found = .true.
                    exit
                endif
            enddo
        endif

        !Find North -----------------
        
        if (.not.Found) then            
            do i = 1, n
                if (AuxString(i:i) == 'W') then
                    AuxString2 = AuxString(1:(i-1))                
                    read(AuxString2,*) AuxReal
                    AuxReal = -1.*AuxReal
                    Found = .true.
                    exit
                endif
            enddo
        endif
        
        if (.not.Found) then
            write(*,*) 
            write(*,*) trim(adjustl(AuxString)), " didn't have N, S, E or W."
            stop 'GetStringUpToChar - ModuleCALMETFormat - ERR01'
        endif
        
    end subroutine ConvertLatLonCharToReal

    !------------------------------------------------------------------------
          
    subroutine SetNewFieldAttributes(Field, Name, Units, Date, WorkSize, nDimensions, Convert)

        !Arguments-------------------------------------------------------------
        type(T_Field),    pointer                 :: Field
        character(len=*),           intent(IN)    :: Name
        character(len=*),           intent(IN)    :: Units        
        type(T_Time)   , optional,  intent(IN)    :: Date
        type(T_Size3D) , optional,  intent(IN)    :: WorkSize
        integer        , optional,  intent(IN)    :: nDimensions
        logical,                    intent(IN)    :: Convert

        !Begin-----------------------------------------------------------------

        Field%Name        = Name        
        Field%Units       = Units
        Field%Convert     = Convert

        if(present(Date       ))    Field%Date          = Date
        if(present(WorkSize   ))    Field%WorkSize      = WorkSize
        if(present(nDimensions))    Field%nDimensions   = nDimensions


    end subroutine SetNewFieldAttributes

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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleCALMETFormat - ERR01'
        

    end subroutine Open_HDF5_OutPut_File
    
    !----------------------------------------------------------------------
    
    subroutine OpenAndReadCALMETFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j, k, it
        logical                                     :: exist
        type(T_Field), pointer                      :: NewField
        type(T_Date ), pointer                      :: CurrentDate
        type(T_Time)                                :: AuxTime
        
        character(len=16)                           :: dataset  ! dataset name (calmet.dat)
        character(len=16)                           :: dataver  ! dataset version
        character(len=64)                           :: datamod  ! dataset message field
        integer                                     :: ncom     ! number of comment lines
        character(len=StringLength)                 :: AuxString

        ! see description in gen.met, grid.met, map.met files 
        integer                                     :: ibyrn,ibmon,ibdyn,ibhrn,ibsecn   
        integer                                     :: ieyrn,iemon,iedyn,iehrn,iesecn
        character(len=8)                            :: axtz             ! UTC base time zone (utc+hhmm)
        integer                                     :: irlg             !Run length (hours) --> nt
        integer                                     :: irtype           !Run type flag (0 - compute wind only, 1 - all)
        integer                                     :: nx, ny, nz
        real                                        :: dgrid            ! Grid size (m)
        real                                        :: xorigr, yorigr   ! Reference grid coordinate of 
                                                                        ! SOUTHWEST corner of grid cell (1,1)
        integer                                     :: iwfcod           ! wind field code (0=objective analysis, 
                                                                        ! 1=diagnostic model, 2=single station wind model)
        integer                                     :: nssta            ! No. surface wind stations
        integer                                     :: nusta            ! No. upper air wind stations
        integer                                     :: npsta            ! No. precipitation stations
        integer                                     :: nowsta           ! No. of overwater met stations
        integer                                     :: nlu              ! Number of land use categories
        integer                                     :: iwat1, iwat2     ! Range of land use categories defining water
        logical                                     :: lcalgrd          ! flag controlling output of special data fields needed by 
                                                                        ! CALGRID (3-D fields of W and temperature)
        character(len=8)                            :: pmap             ! Map projection
                                                                        !UTM :  Universal Transverse Mercator
                                                                        !LCC :  Lambert Conformal Conic
                                                                        !PS  :  Polar Stereographic
                                                                        !EM  :  Equatorial Mercator
                                                                        !LAZA:  Lambert Azimuthal Equal Area
                                                                        !TTM :  Tangential Transverse Mercator
        character(len=8)                            :: datum            ! Datum
        character(len=12)                           :: daten            ! NIMA date (mm-dd-yyy) for datum definitions
        real                                        :: feast,fnorth     ! False Easting and Northing at projection origin
        character(len=4)                            :: utmhem           ! Base hemisphere for UTM projection (N or S)
        integer                                     :: iutmzn           ! UTM zone for UTM projection
        real                                        :: rnlat0,relon0    ! lat & long of x=0 and y=0 of map projection 
                                                                        ! (Used only if PMAP= LCC, PS, EM, TTM or LAZA)
        real                                        :: xlat1,xlat2      ! Matching N. latitude(s) for projection 
                                                                        ! (Used only if PMAP= LCC, PS, or EM)
                                                                        ! LCC :  Projection cone slices through Earth's surface at XLAT1 and XLAT2
                                                                        ! PS  :  Projection plane slices through Earth at XLAT1
                                                                        ! EM  :  Projection cylinder slices through Earth at [+/-] XLAT1        

        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB

        character(len=8)                            :: clabel, CALMETName   ! Variable name
        integer                                     :: ndathrb  ! Beginning Date and time of data (YYYYJJJHH) (explicit time)
        integer                                     :: ibsec    ! Beginning Second of data (SSSS)
        integer                                     :: ndathre  ! Ending Date and time of data (YYYYJJJHH) (explicit time)
        integer                                     :: iesec    ! Ending Second of data (SSSS)
        real,pointer, dimension(:)                  :: DataAux1D
        real,pointer, dimension(:,:,:)              :: DataAux3D1, DataAux3D2, DataAux3D3
        character(len=StringLength)                 :: MohidName, Units
        
        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading CALMET output file...'

        nullify(NewField      )
        nullify(Me%FirstField )
        nullify(Me%FirstDate  )

        !Verifies if file exists
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'CALMET output file does not exist'
            stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR01'
        endif

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR02'

        open(Unit   = Me%Unit,          &
             File   = Me%FileName,      &
             Form   = 'UNFORMATTED',    &
             STATUS = 'OLD',            &
             Action = 'READ',           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR03'

        
        rewind(Me%Unit)
        read  (Me%Unit) dataset,dataver,datamod
        print*, 'dataset = ', trim(adjustl(dataset)), '.'
        print*, 'dataver = ', trim(adjustl(dataver)), '.'
        print*, 'datamod = ', trim(adjustl(datamod)), '.'
        
        read  (Me%Unit) ncom
        print*, 'ncom = ', ncom
        
        do i=1, ncom
            read  (Me%Unit) AuxString
            !print*, 'i, AuxString = ', i, trim(adjustl(AuxString))            
        enddo

        !Only for 2.1 version??
        read (Me%Unit)  ibyrn,ibmon,ibdyn,ibhrn,ibsecn,                     &
                        ieyrn,iemon,iedyn,iehrn,iesecn,                     &
                        axtz,irlg,irtype,                                   &
                        nx, ny, nz, dgrid,                                  &
                        xorigr, yorigr, iwfcod,                             &
                        nssta, nusta, npsta, nowsta,                        &
                        nlu, iwat1, iwat2, lcalgrd,                         &
                        pmap,datum,daten,feast,fnorth,utmhem,iutmzn,        &
                        rnlat0,relon0,xlat1,xlat2                           
        
        print*, 'ibyrn,ibmon,ibdyn,ibhrn,ibsecn = ', ibyrn,ibmon,ibdyn,ibhrn,ibsecn
        print*, 'ieyrn,iemon,iedyn,iehrn,iesecn = ', ieyrn,iemon,iedyn,iehrn,iesecn
        print*, 'axtz,irlg,irtype               = ', axtz,irlg,irtype
        print*, 'nx, ny, nz, dgrid              = ', nx, ny, nz, dgrid
        print*, 'xorigr, yorigr, iwfcod         = ', xorigr, yorigr, iwfcod
        print*, 'nssta, nusta, npsta, nowsta    = ', nssta, nusta, npsta, nowsta
        print*, 'nlu, iwat1, iwat2, lcalgrd     = ', nlu, iwat1, iwat2, lcalgrd
        print*, 'pmap,datum,daten               = ', pmap,datum,daten
        print*, 'feast,fnorth                   = ', feast,fnorth
        print*, 'utmhem,iutmzn                  = ', utmhem,iutmzn
        print*, 'rnlat0,relon0,xlat1,xlat2      = ', rnlat0,relon0,xlat1,xlat2 
        
        !call SetDate(AuxTime, ibyrn, ibmon, ibdyn, ibhrn, ibsecn, 0)
        !call SetDate(Me%FileEndTime, ieyrn, iemon, iedyn, iehrn, iesecn, 0)
        
        if (.not. Me%ConvertTERRAIN) then

            Me%Size%ILB = 0             ; Me%WorkSize%ILB = Me%Size%ILB + 1
            Me%Size%IUB = nx + 1        ; Me%WorkSize%IUB = Me%Size%IUB - 1
            Me%Size%JLB = 0             ; Me%WorkSize%JLB = Me%Size%JLB + 1
            Me%Size%JUB = ny + 1        ; Me%WorkSize%JUB = Me%Size%JUB - 1
    
            Me%DX = dgrid
            Me%DY = dgrid
            Me%FalseEast = feast
            Me%FalseNorth = fnorth
                
            Me%XOri = xorigr
            Me%YOri = yorigr

            if (trim(adjustl(pmap)) == 'LCC') then
                
                Me%ProjType = LAMB_CONF_CONIC_    
                Me%CenterLat    = rnlat0
                Me%CenterLon    = relon0
                Me%TrueLatLower = xlat1
                Me%TrueLatUpper = xlat2
        
                Me%Params(1) = 'proj=lcc'
                Me%Params(2) = 'ellps=WGS84'
                write(Me%Params(3),'(a6,f8.3)') 'lat_1=', Me%TrueLatLower
                write(Me%Params(4),'(a6,f8.3)') 'lat_2=', Me%TrueLatUpper
                write(Me%Params(5),'(a6,f8.3)') 'lon_0=', Me%CenterLon
                write(Me%Params(6),'(a6,f8.3)') 'lat_0=', Me%CenterLat
                write(Me%Params(7),'(a4,f12.3)') 'x_0=' , Me%FalseEast
                write(Me%Params(8),'(a4,f12.3)') 'y_0=' , Me%FalseNorth
            
            else
                write(*,*) 'Map Projection = ', trim(adjustl(pmap))
                write(*,*) 'Not ready for projections other than LCC'
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR04'
            endif

            write(*,*) Me%Params
        
            STAT_CALL = prj90_init(Me%Proj,Me%Params)
            call handle_proj_error(STAT_CALL); 
            if (STAT_CALL /= PRJ90_NOERR) stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR05'
        
        else ! ConvertTERRAIN = TRUE, check that parameters are the same
            
            if (Me%WorkSize%IUB /= nx) then
                write(*,*) 'TERRAIN nx differs from CALMET nx = ', Me%WorkSize%IUB, nx
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR06'
            endif
            
            if (Me%WorkSize%JUB /= ny) then
                write(*,*) 'TERRAIN ny differs from CALMET ny = ', Me%WorkSize%JUB, ny
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR07'
            endif

            if (Me%DX /= dgrid) then
                write(*,*) 'TERRAIN dx differs from CALMET dgrid = ', Me%DX, dgrid
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR08'
            endif

            if (Me%DY /= dgrid) then
                write(*,*) 'TERRAIN dy differs from CALMET dgrid = ', Me%DY, dgrid
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR09'
            endif

            if (Me%FalseEast /= feast) then
                write(*,*) 'TERRAIN False Easting differs from CALMET = ', Me%FalseEast, feast
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR10'
            endif
                
            if (Me%FalseNorth /= fnorth) then
                write(*,*) 'TERRAIN False Northing differs from CALMET = ', Me%FalseNorth, fnorth
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR11'
            endif
                
            if (pmap /= 'LCC') then
                write(*,*) 'TERRAIN projection differs from CALMET = ',  Me%ProjType, pmap
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR12'
            endif

            if (Me%CenterLat /= rnlat0) then
                write(*,*) 'TERRAIN CenterLat differs from CALMET = ', Me%CenterLat, rnlat0
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR13'
            endif
                
            if (Me%CenterLon /= relon0) then
                write(*,*) 'TERRAIN CenterLon differs from CALMET = ', Me%CenterLon, relon0
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR14'
            endif

            if (Me%TrueLatLower /= xlat1) then
                write(*,*) 'TERRAIN TrueLatLower differs from CALMET = ', Me%TrueLatLower, xlat1
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR15'
            endif
                                
            if (Me%TrueLatUpper /= xlat2) then
                write(*,*) 'TERRAIN TrueLatUpper differs from CALMET = ', Me%TrueLatUpper, xlat2
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR16'
            endif

            if (Me%XOri /= xorigr) then
                write(*,*) 'TERRAIN XOri differs from CALMET = ', Me%XOri, xorigr
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR17'
            endif
            
            if (Me%YOri /= yorigr) then
                write(*,*) 'TERRAIN YOri differs from CALMET = ', Me%YOri, yorigr
                stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR18'
            endif

        endif !ConvertTERRAIN
        
        Me%Size%KLB = 0             ; Me%WorkSize%KLB = Me%Size%KLB + 1
        Me%Size%KUB = nz + 1        ; Me%WorkSize%KUB = Me%Size%KUB - 1
            
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 
                
        !----------------------------------------------------------------------
        ! Static Information 
        !----------------------------------------------------------------------
        
        ! ZFACE (m)            
        allocate(Me%ZFace(nz+1))                      
        read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, Me%ZFace
        
        if (clabel /= 'ZFACE') then
            write(*,*) 'Property should be ', trim(adjustl(CALMETName))
            stop 'OpenAndReadCALMETFile - ModuleCALMETFormat - ERR17'
        endif
        
        ! if stations 
        
        if (nssta > 1) then                    
            allocate(DataAux1D(nssta))              
            !xssta
            read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux1D
            !yssta
            read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux1D
            deallocate(DataAux1D)
        endif
        
        if (nusta > 1) then                    
            allocate(DataAux1D(nusta))              
            !XUSTA
            read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux1D
            !YUSTA
            read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux1D
            deallocate(DataAux1D)
        endif

        if (npsta > 1) then                    
            allocate(DataAux1D(npsta))              
            !XPSTA
            read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux1D
            !YPSTA
            read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux1D
            deallocate(DataAux1D)
        endif
        
        !Z0 - surface roughness lengths (NX * NY words)     
        CALMETName  = 'Z0'
        MohidName   = CALMETName
        Units       = 'm'
        call AddFieldXY (CALMETName, MohidName, Units, NewField)
        
        !ILANDU - land use categories (NX * NY words)
        CALMETName  = 'ILANDU'
        MohidName   = 'LandUse'
        Units       = '-'
        call AddIntegerFieldXY (CALMETName, MohidName, Units, NewField, ConvertToReal = .TRUE.)
                
        !ELEV -  elevations (NX * NY words)
        CALMETName  = 'ELEV'
        MohidName   = CALMETName
        Units       = 'm'
        call AddFieldXY (CALMETName, MohidName, Units, NewField)

        nullify(Me%Bathymetry)        
        allocate(Me%Bathymetry(ILB:IUB,JLB:JUB))
        Me%Bathymetry(WILB:WIUB,WJLB:WJUB) = NewField%Values2D(WILB:WIUB,WJLB:WJUB)    

        do j = WJLB, WJUB
        do i = WILB, WIUB
            if (Me%Bathymetry(i,j) < -99.) Me%Bathymetry(i,j) = -99.
        enddo
        enddo
        
        !XLAI - leaf area index (NX * NY words)        
        CALMETName  = 'XLAI'
        MohidName   = CALMETName
        Units       = '-'
        call AddFieldXY (CALMETName, MohidName, Units, NewField)

        !NEARS -- Number of the closest surface met station        
        if (nssta > 1) then
            CALMETName  = 'NEARS'
            MohidName   = CALMETName
            Units       = '-'
            call AddFieldXY (CALMETName, MohidName, Units, NewField)
        endif

        !----------------------------------------------------------------------
        ! Time Dependent Variables
        !----------------------------------------------------------------------
        
        do it = 1, irlg

            nullify (CurrentDate)
            
            !Wind Components ------------------------------------------------
            
            allocate(DataAux3D1(WILB:WIUB, WJLB:WJUB, WKLB:WKUB))            
            allocate(DataAux3D2(WILB:WIUB, WJLB:WJUB, WKLB:WKUB))            
            if (lcalgrd) allocate(DataAux3D3(WILB:WIUB, WJLB:WJUB, WKLB:WKUB+1))            
                        
            do k = 1, nz

                !U
                read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux3D1(:,:,k)

                !First time date and time are read from file
                if (k==1) then
                    call SetNewDate(CurrentDate, ndathr = ndathrb, isec=ibsec)            
                    AuxTime = CurrentDate%Date
                endif
                                
                !V
                read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux3D2(:,:,k)

                if(lcalgrd)then

                    !WFACE - W velocities at TOP cell face are written (NZ fields in all)
                    read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux3D3(:,:,k+1)
                    
                endif

            enddo ! do k

            !Velocity U            
            nullify(NewField)            
            call AddField(Me%FirstField, NewField)            
        
            MohidName = trim(adjustl(GetPropertyName(WindVelocityX_)))//'_3D'
            
            call SetNewFieldAttributes(Field         = NewField,                                &  
                                       Name          = MohidName,                               &
                                       Units         = 'm/s',                                   &
                                       WorkSize      = Me%WorkSize,                             &
                                       nDimensions   = 3,                                       &
                                       Date          = CurrentDate%Date,                        &
                                       Convert       = .true.) !FieldIsToConvert(MohidName))
        
            allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))        
            NewField%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) = DataAux3D1(WILB:WIUB, WJLB:WJUB, WKLB:WKUB)
            deallocate(DataAux3D1)
        
            !Velocity V
            nullify(NewField)            
            call AddField(Me%FirstField, NewField)            

            MohidName = trim(adjustl(GetPropertyName(WindVelocityY_)))//'_3D'
            
            call SetNewFieldAttributes(Field         = NewField,                                &  
                                       Name          = MohidName,                               &
                                       Units         = 'm/s',                                   &
                                       WorkSize      = Me%WorkSize,                             &
                                       nDimensions   = 3,                                       &
                                       Date          = CurrentDate%Date,                        &
                                       Convert       = .true.) !FieldIsToConvert(MohidName))
            
            allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))        
            NewField%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) = DataAux3D2(WILB:WIUB, WJLB:WJUB, WKLB:WKUB)
            deallocate(DataAux3D2)

            if (lcalgrd) then
                !Velocity W
                nullify(NewField)            
                call AddField(Me%FirstField, NewField)            
        
                call SetNewFieldAttributes(Field         = NewField,                            &  
                                           Name          = 'wind velocity Z_3D',                &
                                           Units         = 'm/s',                               &
                                           WorkSize      = Me%WorkSize,                         &
                                           nDimensions   = 3,                                   &
                                           Date          = CurrentDate%Date,                    &
                                           Convert       = .true.) !FieldIsToConvert(MohidName))
        
                allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))        
                DataAux3D3(:,:,WKLB) = 0.0
                NewField%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:(WKUB+1)) = DataAux3D3(WILB:WIUB, WJLB:WJUB, WKLB:(WKUB+1))
                deallocate(DataAux3D3)
        
            endif
   
if1:        if (irtype > 0) then

                !T-3D                
                allocate(DataAux3D1(WILB:WIUB, WJLB:WJUB, WKLB:WKUB))            
    
                do k = 1, nz
                    read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux3D1(:,:,k)
                enddo

                nullify(NewField)            
                call AddField(Me%FirstField, NewField)            

                MohidName = trim(adjustl(GetPropertyName(AirTemperature_)))//'_3D'
                
                call SetNewFieldAttributes(Field         = NewField,                                &  
                                           Name          = MohidName,                               &
                                           Units         = 'ºC',                                    &
                                           WorkSize      = Me%WorkSize,                             &
                                           nDimensions   = 3,                                       &
                                           Date          = CurrentDate%Date,                        &
                                           Convert       = .true.) !FieldIsToConvert(MohidName))
                
        
                allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))        
                NewField%Values3D(WILB:WIUB, WJLB:WJUB, WKLB:WKUB) = DataAux3D1(WILB:WIUB, WJLB:WJUB, WKLB:WKUB)
                deallocate(DataAux3D1)


                ! other 2d meteo variables ------------------------------------

                ! PGT stability class (1-6)
                CALMETName  = 'IPGT'
                MohidName   = CALMETName
                Units       = '-'
                call AddIntegerFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)

                ! Friction velocity (m/s)
                CALMETName  = 'USTAR'
                MohidName   = CALMETName
                Units       = 'm/s'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)

                ! Mixing height (m) from MIXHT
                CALMETName  = 'ZI'
                MohidName   = CALMETName
                Units       = 'm'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)
                
                ! Monin-Obukhov length (m)
                CALMETName  = 'EL'
                MohidName   = CALMETName
                Units       = 'm'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)

                ! Convective velocity scale (m/s)
                CALMETName  = 'WSTAR'
                MohidName   = CALMETName
                Units       = 'm/s'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)

                !! Precip. rate (mm/hr)
                !CALMETName  = 'RMM'
                !MohidName   = CALMETName
                !Units       = 'mm/h'
                !call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)

                !Air Temperature (K)
                CALMETName  = 'TEMPK'
                MohidName   = GetPropertyName(AirTemperature_)
                Units       = 'degC'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)                
                NewField%Values2D = NewField%Values2D - 273.15
                
                !Air density (kg/m**3)
                CALMETName  = 'RHO'
                MohidName   = 'air density'
                Units       = 'kg/m3'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)
 
                !Short-wave solar radiation (W/m**2)
                CALMETName  = 'QSW'
                MohidName   = 'solar radiation'
                Units       = 'W/m2'
                call AddFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)                

                !Relative humidity (percent)
                CALMETName  = 'IRH'
                MohidName   = 'relative humidity'
                Units       = '%'
                call AddIntegerFieldXY (CALMETName, MohidName, Units, NewField, AuxTime, ConvertToReal=.true.)
                
                !!Precipitation Code
                !CALMETName  = 'IPCODE'
                !MohidName   = CALMETName
                !Units       = '-'
                !call AddIntegerFieldXY (CALMETName, MohidName, Units, NewField, AuxTime)
                !print*, CALMETName, NewField%IValues2D(50,50)
                
                
            endif if1

        enddo ! do times
        
        !Close CALMET File
        close(Me%Unit)
        
    end subroutine OpenAndReadCALMETFile

    !--------------------------------------------------------------------------
    
    subroutine AddFieldXY (CALMETName, MohidName, Units, NewField, CurrentTime)
        
        !Arguments --------------------------------------------------------
        character(8), intent(in)                :: CALMETName
        character(StringLength), intent(in)     :: MohidName, Units
        type(T_Field), pointer                  :: NewField
        type(T_Time), optional                  :: CurrentTime
        
        !Local ------------------------------------------------------------
        real, pointer, dimension(:,:)           :: DataAux2D
        character(len=8)                        :: clabel  ! Variable name
        integer                                 :: ndathrb, ibsec, ndathre, iesec
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: WILB, WIUB, WJLB, WJUB
        
        !Begin ------------------------------------------------------------
            
        nullify(DataAux2D)
                       
        allocate(DataAux2D(Me%WorkSize%IUB, Me%WorkSize%JUB))
            
        read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux2D
        
        if (clabel /= trim(adjustl(CALMETName))) then
            write(*,*) 'Property should be ', trim(adjustl(CALMETName))
            write(*,*) 'instead of ', clabel
            stop 'ReadStaticVarXY - ModuleCALMETFormat - ERR17'
        endif
       
        call AddField(Me%FirstField, NewField)
    
        call SetNewFieldAttributes(Field         = NewField,         &
                                   Name          = MohidName,        &
                                   Units         = Units,            &
                                   WorkSize      = Me%WorkSize,      &
                                   nDimensions   = 2,                &
                                   Convert       = .true.)

        if (present(CurrentTime)) NewField%Date = CurrentTime
        
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        
        allocate(NewField%Values2D(ILB:IUB, JLB:JUB))        
        NewField%Values2D(WILB:WIUB, WJLB:WJUB) = DataAux2D(WILB:WIUB, WJLB:WJUB)
        deallocate(DataAux2D)

    end subroutine AddFieldXY
    
    !--------------------------------------------------------------------------
    
    subroutine AddIntegerFieldXY (CALMETName, MohidName, Units, NewField, CurrentTime, ConvertToReal)
        
        !Arguments --------------------------------------------------------
        character(8), intent(in)                :: CALMETName
        character(StringLength), intent(in)     :: MohidName, Units
        type(T_Field), pointer                  :: NewField
        type(T_Time), optional                  :: CurrentTime
        logical, optional                       :: ConvertToReal 
        
        !Local ------------------------------------------------------------
        integer, pointer, dimension(:,:)        :: DataAux2D
        character(len=8)                        :: clabel  ! Variable name
        integer                                 :: ndathrb, ibsec, ndathre, iesec
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: WILB, WIUB, WJLB, WJUB
        
        !Begin ------------------------------------------------------------
            
        nullify(DataAux2D)
                       
        allocate(DataAux2D(Me%WorkSize%IUB, Me%WorkSize%JUB))
            
        read(Me%Unit) clabel, ndathrb, ibsec, ndathre, iesec, DataAux2D
        
        if (clabel /= trim(adjustl(CALMETName))) then
            write(*,*) 'Property should be ', trim(adjustl(CALMETName))
            write(*,*) 'instead of ', clabel
            stop 'ReadStaticVarXY - ModuleCALMETFormat - ERR17'
        endif
       
        call AddField(Me%FirstField, NewField)

        call SetNewFieldAttributes(Field         = NewField,         &
                                   Name          = MohidName,        &
                                   Units         = Units,            &
                                   WorkSize      = Me%WorkSize,      &
                                   nDimensions   = 2,                &
                                   Convert       = .true.) !!FieldIsToConvert(MohidName))

        if (present(CurrentTime)) NewField%Date = CurrentTime
        
        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        
        if (present(ConvertToReal)) then

            allocate(NewField%Values2D(ILB:IUB, JLB:JUB))        
            NewField%Values2D(WILB:WIUB, WJLB:WJUB) = real(DataAux2D(WILB:WIUB, WJLB:WJUB))
            
        else 
            allocate(NewField%IValues2D(ILB:IUB, JLB:JUB))        
            NewField%IValues2D(WILB:WIUB, WJLB:WJUB) = DataAux2D(WILB:WIUB, WJLB:WJUB)
            NewField%IsInteger  = .true.

        endif
        
        deallocate(DataAux2D)

        
    end subroutine AddIntegerFieldXY
    
    !------------------------------------------------------------------------

    
    subroutine SetNewDate(NewDate, ndathr, isec, Year, Month, Day, Hour, Minute, Second)

        !Arguments-------------------------------------------------------------
        type(T_Date), pointer, intent(OUT)          :: NewDate        
        integer, optional, intent(IN)               :: ndathr  ! Date and time of data (YYYYJJJHH) (explicit time)
        integer, optional, intent(IN)               :: isec    ! Second of data (SSSS)
        integer, optional, intent(IN)               :: Year, Month, Day, Hour, Minute, Second
        
        !Local-----------------------------------------------------------------
        integer                                     :: IYear2, JulianDay, IHour2
        character(9)                                :: YYYYJJJHH
        type(T_Time)                                :: AuxTime

        !Begin-----------------------------------------------------------------

        if (.not. present(ndathr) .and. .not. present(Year)) then
            write(*,*) 'Optional arguments missing. Specify:'
            write(*,*) 'ndathr, isec or'
            write(*,*) 'integers: Year, Month, Day, Hour, Minute, Second'            
            stop 'SetNewDate - ModuleCALMETFormat -ERR01'
        endif
        
        call AddDate(Me%FirstDate, NewDate)
                
        if (present(ndathr)) then
            
            write(YYYYJJJHH, '(I9)') ndathr                
            read(YYYYJJJHH(1:4),*) IYear2
            read(YYYYJJJHH(5:7),*) JulianDay
            read(YYYYJJJHH(8:9),*) IHour2
        
            call JulianDayToMonthDay(IYear2, JulianDay, AuxTime)
            
            NewDate%Date = AuxTime + Ihour2*60. + isec
        
        else

            call SetDate(NewDate%Date, Year, Month, Day, Hour, Minute, Second)
        
        endif
        

    end subroutine SetNewDate    
    
    !------------------------------------------------------------------------    
    
    subroutine BeginEndTime

        !Local-------------------------------------------------------------
        type(T_Date), pointer                   :: CurrentDate
        type(T_Time)                            :: StartTime, EndTime
        real                                    :: Year, Month, Day
        real                                    :: Hour, Minute, Second

100 format (1x, f5.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0, 1x, f3.0)
        
        CurrentDate => Me%FirstDate        
        StartTime = CurrentDate%Date

        do while(associated(CurrentDate))                

            EndTime = CurrentDate%Date
            CurrentDate => CurrentDate%Next

        end do 
        
if1:    if (Me%TimeWindow) then

            if (Me%StartTime .GT. EndTime .OR. Me%EndTime .LT. StartTime) then
                
                write (*,*) 
                write (*,*) 'Time Window not in file'
                write (*,*) 'File START at:'
                call ExtractDate(StartTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,*) 'File END at:'
                call ExtractDate(EndTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,fmt=100) Year, Month, Day, Hour, Minute, Second
                stop 'BeginEndTime - ModuleCALMETFormat - ERR01'                

            endif


            if (Me%StartTime .LT. StartTime) then
                
                write (*,*) 
                write (*,*) 'START is before the start time of file'
                write (*,*) 'File begins at:'
                call ExtractDate(StartTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write (*,*) 'BeginEndTime - ModuleCALMETFormat - WARN01'                

            endif

            if (Me%EndTime .GT. EndTime) then                
                
                write (*,*) 
                write (*,*) 'END is after the end time of file'
                write (*,*) 'File ends at:'
                call ExtractDate(EndTime, Year, Month, Day, Hour, Minute, Second)        
                write (*,fmt=100) Year, Month, Day, Hour, Minute, Second
                write (*,*) 'BeginEndTime - ModuleCALMETFormat - WARN02'                
            
            endif
        
        else !if1

            Me%StartTime = StartTime
            Me%EndTime   = EndTime

        endif if1


    end subroutine BeginEndTime    

    
        !------------------------------------------------------------------------

    subroutine ComputeVerticalCoordinate
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                  :: NewField
        type(T_Date), pointer                   :: CurrentDate        
        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: k

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Computing Vertical Coordinate...'

        ILB = Me%Size%ILB; WILB = Me%WorkSize%ILB 
        IUB = Me%Size%IUB; WIUB = Me%WorkSize%IUB 
        JLB = Me%Size%JLB; WJLB = Me%WorkSize%JLB 
        JUB = Me%Size%JUB; WJUB = Me%WorkSize%JUB 
        KLB = Me%Size%KLB; WKLB = Me%WorkSize%KLB 
        KUB = Me%Size%KUB; WKUB = Me%WorkSize%KUB 

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            call AddField(Me%FirstField, NewField)

            call SetNewFieldAttributes(Field         = NewField,         &
                                       Name          = 'VerticalZ',      &
                                       Units         = 'm',              &
                                       Date          = CurrentDate%Date, &
                                       WorkSize      = Me%WorkSize,      &
                                       nDimensions   = 3,                &
                                       Convert       = .true.)
                
            allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))        
            do k = WKLB, WKUB + 1
                NewField%Values3D(:,:,k) = Me%ZFace(k)
            enddo

            CurrentDate => CurrentDate%Next

        end do 


    end subroutine ComputeVerticalCoordinate

    !------------------------------------------------------------------------

    subroutine ComputeWindSurface
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate
        real, dimension(:,:,:), pointer                 :: WindVelocityX
        real, dimension(:,:,:), pointer                 :: WindVelocityY
        logical                                         :: VelocityX_OK     = .false.
        logical                                         :: VelocityY_OK     = .false.
        type(T_Field), pointer                          :: WindModulus, WindDirection, WindX, WindY
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: WILB, WIUB, WJLB, WJUB
        integer                                         :: i, j
        real                                            :: wd

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 

        ILB  = Me%Size%ILB
        IUB  = Me%Size%IUB
        JLB  = Me%Size%JLB
        JUB  = Me%Size%JUB

        write(*,*)
        write(*,*)'Computing wind modulus ...'

        CurrentDate => Me%FirstDate

        do while(associated(CurrentDate))

            Field => Me%FirstField

            do while(associated(Field))
            
                if(Field%Name == trim(adjustl(GetPropertyName(WindVelocityX_)))//'_3D'    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    WindVelocityX       => Field%Values3D
                    VelocityX_OK        = .true.

                end if

                if(Field%Name == trim(adjustl(GetPropertyName(WindVelocityY_)))//'_3D'    .and. &
                   Field%Date == CurrentDate%Date)then
                    
                    WindVelocityY       => Field%Values3D
                    VelocityY_OK        = .true.

                end if


                if(VelocityX_OK .and. VelocityY_OK) then

                    
                    call AddField(Me%FirstField, WindModulus)
                    call AddField(Me%FirstField, WindDirection)                    
                    call AddField(Me%FirstField, WindX)
                    call AddField(Me%FirstField, WindY)                    
                    
                    call SetNewFieldAttributes(Field         = WindModulus,                         &
                                               Name          = GetPropertyName(WindModulus_),       &
                                               Units         = 'm/s',                               &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)

                    call SetNewFieldAttributes(Field         = WindDirection,                       &
                                               Name          = GetPropertyName(WindDirection_),     &
                                               Units         = 'degrees',                           &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)

                    call SetNewFieldAttributes(Field         = WindX,                               &
                                               Name          = GetPropertyName(WindVelocityX_),     &
                                               Units         = 'm/s',                               &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)

                    call SetNewFieldAttributes(Field         = WindY,                               &
                                               Name          = GetPropertyName(WindVelocityY_),     &
                                               Units         = 'm/s',                               &
                                               Date          = CurrentDate%Date,                    &
                                               WorkSize      = Me%WorkSize,                         &
                                               nDimensions   = 2,                                   &
                                               Convert       = .true.)
                    
                    allocate(WindModulus%Values2D   (ILB:IUB, JLB:JUB))
                    allocate(WindDirection%Values2D (ILB:IUB, JLB:JUB))
                    allocate(WindX%Values2D         (ILB:IUB, JLB:JUB))
                    allocate(WindY%Values2D         (ILB:IUB, JLB:JUB))

                    do i = WILB, WIUB
                    do j = WJLB, WJUB

                        WindModulus%Values2D(i,j)= sqrt(WindVelocityX(i,j,1)**2+WindVelocityY(i,j,1)**2.)

                        wd = 3.* Pi / 2. - atan2(WindVelocityY(i,j,1), WindVelocityX(i,j,1))
                        wd = wd * 180 / Pi

                        if (wd .gt. 360.) wd = wd - 360.

                        WindDirection%Values2D(i,j) = wd
                        
                    enddo
                    enddo
                    
                    WindX%Values2D(WILB:WIUB,WJLB:WJUB) = WindVelocityX(WILB:WIUB,WJLB:WJUB,1)
                    WindY%Values2D(WILB:WIUB,WJLB:WJUB) = WindVelocityY(WILB:WIUB,WJLB:WJUB,1)

                    nullify(WindVelocityX    )
                    nullify(WindVelocityY    )
                    
                    VelocityX_OK     = .false.
                    VelocityY_OK     = .false.
                    

                end if

                Field => Field%Next

            end do


            CurrentDate => CurrentDate%Next

        end do           
            
    end subroutine ComputeWindSurface

    !----------------------------------------------------------------------

    subroutine WriteGridToHDF5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer,    dimension(:,:), pointer         :: WaterPoints2D
        type(T_Field), pointer                      :: Field
        logical                                     :: LandUse_OK, Z0_OK, XLAI_OK


        !----------------------------------------------------------------------

        !call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB+1,&
        !                        Me%WorkSize%JLB, Me%WorkSize%JUB+1, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR01'
        !
        !call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX", "-",       &
        !                        Array2D =  Me%ConnectionX, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR02'
        !
        !call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY", "-",       &
        !                        Array2D =  Me%ConnectionY, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR03'

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR04'

        allocate(WaterPoints2D(Me%Size%ILB:Me%Size%IUB,&
                               Me%Size%JLB:Me%Size%JUB))
     
        WaterPoints2D = 1

        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, &
                             Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR05'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",    &
                              Array2D = WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR06'
        
        !Bathymetry
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR07'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR08'

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR09'

        nullify(Field)
        Field => Me%FirstField
        LandUse_OK   = .false.
        Z0_OK       = .false.
        XLAI_OK     = .false.
        
        do while (associated(Field))
            
            select case (Field%Name)
                
            case ('LandUse')
                call HDF5WriteData   (Me%ObjHDF5, "/Grid", "LandUse", "-",       &
                                        Array2D =  Field%Values2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR10'

                LandUse_OK = .true.
                                    
            case ('Z0')
                call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Z0", "-",       &
                                        Array2D =  Field%Values2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR11'

                Z0_OK = .true.

            case ('XLAI')
                call HDF5WriteData   (Me%ObjHDF5, "/Grid", "XLAI", "-",       &
                                        Array2D =  Field%Values2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR12'
                
                XLAI_OK = .true.
                
            end select

            if (LandUse_OK .and. Z0_OK .and. XLAI_OK) exit
            
            Field => Field%Next
            
        enddo
        
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridToHDF5File - ModuleCALMETFormat - ERR13'

        deallocate(WaterPoints2D)
        nullify   (WaterPoints2D)

    end subroutine WriteGridToHDF5File
    
    !------------------------------------------------------------------------
    
    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, OutputNumber
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate, PreviousOutputDate
        real                                            :: dt

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Writing HDF5 file...'

        call WriteGridToHDF5File
        
        OutputNumber = 1
        CurrentDate        => Me%FirstDate
        PreviousOutputDate => CurrentDate
        dt = 0.0

        do while(associated(CurrentDate))

            call ExtractDate   (CurrentDate%Date,                           &
                                AuxTime(1), AuxTime(2), AuxTime(3),         &
                                AuxTime(4), AuxTime(5), AuxTime(6))


ifT:        if (CurrentDate%Date .GE. Me%StartTime .AND. CurrentDate%Date .LE. Me%EndTime) then            
            dt = CurrentDate%Date - PreviousOutputDate%Date
ifDT:       if (CurrentDate%Date .EQ. Me%FirstDate%Date .OR. dt >= Me%OutputDTInterval) then

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR30'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR40'


            Field => Me%FirstField

            do while(associated(Field))

if1:            if (Field%Convert .and. (Field%Date == CurrentDate%Date)) then
                
                    Field%OutputNumber = OutputNumber

if2:                if(Field%nDimensions == 2)then

                        call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                            Field%WorkSize%JLB, Field%WorkSize%JUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR50'

                        if (Field%IsInteger) then
                                
                            call HDF5WriteData(Me%ObjHDF5,                                       &
                                                "/Results/"//Field%Name,                         &
                                                Field%Name,                                      &
                                                Field%Units,                                     &
                                                Array2D      = Field%IValues2D,                  &
                                                OutputNumber = Field%OutputNumber,               &
                                                STAT         = STAT_CALL)
                        else
                            call HDF5WriteData(Me%ObjHDF5,                                       &
                                                "/Results/"//Field%Name,                         &
                                                Field%Name,                                      &
                                                Field%Units,                                     &
                                                Array2D      = Field%Values2D,                   &
                                                OutputNumber = Field%OutputNumber,               &
                                                STAT         = STAT_CALL)
                        endif
                            
                        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR60'


                    elseif(Field%nDimensions == 3) then !if2

                        if(Field%Name == 'VerticalZ') then

                            call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,  &
                                                Field%WorkSize%JLB, Field%WorkSize%JUB,             &
                                                Field%WorkSize%KLB - 1, Field%WorkSize%KUB, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR70'

                            call HDF5WriteData(Me%ObjHDF5,                                      &
                                            "/Grid/"//Field%Name,                               &
                                            'Vertical',                                         &
                                            Field%Units,                                        &
                                            Array3D      = Field%Values3D,                      &
                                            OutputNumber = Field%OutputNumber,                  &
                                            STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR80'

                        else
                    
                            call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                                Field%WorkSize%JLB, Field%WorkSize%JUB,            &
                                                Field%WorkSize%KLB, Field%WorkSize%KUB, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR90'

                            call HDF5WriteData(Me%ObjHDF5,                                      &
                                                "/Results/"//Field%Name,                         &
                                                Field%Name,                                      &
                                                Field%Units,                                     &
                                                Array3D      = Field%Values3D,                   &
                                                OutputNumber = Field%OutputNumber,               &
                                                STAT         = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR100'

                        end if
                    end if if2
                end if if1

                Field => Field%Next

            end do

            OutputNumber = OutputNumber + 1
            PreviousOutputDate => CurrentDate

            endif ifDT            
            endif ifT

            CurrentDate  => CurrentDate%Next

        end do

        write(*,*)
        write(*,*)'Closing HDF5 file...'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCALMETFormat - ERR110'



    end subroutine OutputFields
    
    !--------------------------------------------------------------------------
    
    subroutine KillCALMETFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCALMETFormat - ModuleCALMETFormat - ERR01'
        
        deallocate(Me%Bathymetry )
        deallocate(Me%CenterX    )
        deallocate(Me%ConnectionX)
        deallocate(Me%CenterY    )
        deallocate(Me%ConnectionY)

        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCALMETFormat - ModuleCALMETFormat - ERR02'
        
        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillCALMETFormat - ModuleCALMETFormat - ERR03'

        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillCALMETFormat

    !--------------------------------------------------------------------------

    subroutine handle_proj_error (status)

        integer, intent(in)         :: status

        if (status /= PRJ90_NOERR) write(*,*) trim(prj90_strerrno(status))

    end subroutine handle_proj_error    

    !--------------------------------------------------------------------------    
 
end module ModuleCALMETFormat
