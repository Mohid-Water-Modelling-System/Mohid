!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Tools
! PROGRAM       : HDF5Operator
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : August 2013
! REVISION      : Pedro Chambel Leitão
! DESCRIPTION   : To make arithmetic operations with HDF5
!
!------------------------------------------------------------------------------

program ETCalculator

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleTime
    use ModuleEnterData
    use ModuleFillMatrix
    use ModuleFunctions

    implicit none

    !Parameter
    integer, parameter                   :: SELECTION_ALL = 1
    integer, parameter                   :: SELECTION_TIME_INSTANT = 2
    integer, parameter                   :: SELECTION_INSTANT = 3
    !Sun constant (W / m**2)
    real,    parameter                              :: KSun                 = 1367.0         

    real, dimension(:, :), pointer        :: ConnectionX         => null()
    real, dimension(:, :), pointer        :: ConnectionY         => null()
    real, dimension(:, :), pointer        :: aux_ConnectionX
    real(8), dimension(:, :), pointer        :: PropertyField
    real, dimension(:, :), pointer        :: PropertyFieldOld
    real, dimension(:, :), pointer        :: PropertyField2
    real, dimension(:, :), pointer        :: PropertyFieldOld2
    integer, dimension(:, :), pointer     :: WaterPoints2D
    integer, dimension(:, :), pointer     :: LatitudeGrid, LongitudeGrid
    real                                  :: Latitude
    real                                  :: Longitude
    integer                               :: CoordType

     type (T_Time)                        :: InitialSystemTime, FinalSystemTime, Now
     integer, dimension(8)                :: F95Time
     real                                 :: ElapsedSeconds, TotalCPUTime

     character(PathLength)                :: DataFile  = 'HDF5Operator.dat'

     character(PathLength)                :: HDFFile
     character(PathLength)                :: HDFFileB
     character(PathLength)                :: HDFFileOut
     logical                              :: AllInstants
     logical                              :: OperateTwoHDF5
     integer                              :: NumberOfHDFInstants
     integer                              :: Instant
     integer                              :: GroupRank
     
     real                                 :: TimeInterval
     real                                 :: WindHeight !FLAVIO SANTOS
     logical                              :: XYWind     !FLAVIO SANTOS
     
     type (T_Size2D)                      :: Size


     integer                              :: ObjHDF5           = 0
     integer                              :: ObjHDF5B          = 0
     integer                              :: ObjHDF5_Out       = 0
     integer                              :: ObjHorizontalGrid = 0
     integer                              :: ObjEnterData      = 0
  
    
    !Begin----------------------------------------------------------------------
    
    call StartHDF5Operator


    call ReadDataFile

    call ConvertHDF5Operator
    
    call KillHDF5Operator
    
    
    call ShutHDF5Operator

    
    contains
    
    
    !----------------------------------------------------------------------------
    
    subroutine StartHDF5Operator
    
        call date_and_time(Values = F95Time)
        call StartUpMohid("HDF5Operator")
        call SetDate(InitialSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                     float(F95Time(3)), float(F95Time(5)), &
                                     float(F95Time(6)), float(F95Time(7))+ &
                                     F95Time(8)/1000.)
    end subroutine StartHDF5Operator
    
    !----------------------------------------------------------------------------
    
    subroutine ReadDataFile

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag   
        logical                                     :: Exist   
   
        !Begin-------------------------------------------------------------------

        write(*,*)
        write(*,*)'Reading options...'
        write(*,*)
        
         call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
         if (STAT_CALL /= SUCCESS_) stop 'HDF5Operator - ERR000'

         call GetData(HDFFile,                             &
                      ObjEnterData, iflag,                 &
                      SearchType   = FromFile,             &
                      keyword      = 'HDF_FILE',           &
                      ClientModule = 'HDF5Operator',     &
                      STAT         = STAT_CALL)
         if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR010'

        !Verifies if file exists
        inquire(FILE = HDFFile, EXIST = Exist)
        if (.not. Exist) then
            write(*,*)'HDF5 file does not exist:'//trim(HDFFile)
            stop 'ReadDataFile - HDF5Operator - ERR20'
        endif

        call GetData(OperateTwoHDF5,                         &
                  ObjEnterData, iflag,                    &
                  SearchType   = FromFile,                &
                  keyword      = 'OPERATE_TWO_HDF5',          &
                  Default      = .false.,                 &
                  ClientModule = 'HDF5Operator',        &
                  STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR050'
        
        if (OperateTwoHDF5) then
             call GetData(HDFFileB,                             &
                          ObjEnterData, iflag,                 &
                          SearchType   = FromFile,             &
                          keyword      = 'HDF_FILE2',           &
                          ClientModule = 'HDF5Operator',     &
                          STAT         = STAT_CALL)
             if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR010'

            !Verifies if file exists
            inquire(FILE = HDFFileB, EXIST = Exist)
            if (.not. Exist) then
                write(*,*)'HDF5 file does not exist:'//trim(HDFFileB)
                stop 'ReadDataFile - HDF5Operator - ERR20'
            endif
        endif

        call GetData(HDFFileOut,                             &
                      ObjEnterData, iflag,                 &
                      SearchType   = FromFile,             &
                      keyword      = 'HDF_FILE_OUT',           &
                      ClientModule = 'HDF5Operator',     &
                      STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR010'

        call GetData(AllInstants,                         &
                  ObjEnterData, iflag,                    &
                  SearchType   = FromFile,                &
                  keyword      = 'ALL_INSTANTS',          &
                  Default      = .false.,                 &
                  ClientModule = 'HDF5Operator',        &
                  STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR050'

        if (.not. AllInstants) then
        
             call GetData(Instant,                                &
                          ObjEnterData, iflag,                    &
                          SearchType   = FromFile,                &
                          keyword      = 'INSTANT',               &
                          Default      = 1,                       &
                          ClientModule = 'HDF5Operator',        &
                          STAT         = STAT_CALL)
             if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR060'            
            
        endif   

         call GetData(TimeInterval,                                &
                      ObjEnterData, iflag,                    &
                      SearchType   = FromFile,                &
                      keyword      = 'TIME_INTERVAL',               &
                      Default      = 3.0,                       &
                      ClientModule = 'HDF5Operator',        &
                      STAT         = STAT_CALL)
         if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR060'       
         
         !wind modulus new implementation (FLAVIO SANTOS)
         call GetData(WindHeight,                             &
                      ObjEnterData, iflag,                    &
                      SearchType   = FromFile,                &
                      keyword      = 'WIND_HEIGHT',           &
                      Default      = 2.0,                     &
                      ClientModule = 'HDF5Operator',          &
                      STAT         = STAT_CALL)
         if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR060'       
         
         
         call GetData(XYWind,                                 &
                      ObjEnterData, iflag,                    &
                      SearchType   = FromFile,                &
                      keyword      = 'XY_WIND',               &
                      Default      = .false.,                 &
                      ClientModule = 'HDF5Operator',          &
                      STAT         = STAT_CALL)
         if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5Operator - ERR060'       
              
        
    end subroutine ReadDataFile
    
    !---------------------------------------------------------------------------

    subroutine ConvertHDF5Operator

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: HDF5_READ
        integer                                     :: HDF5_READWRITE
        integer                                     :: HDF5_CREATE            
        logical                                     :: Exist   
        character(len=PathLength)                   :: HDFGroup
        integer  , dimension(7)                     :: Dimensions
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: AuxPosition
        character(len=StringLength)                 :: Property
        integer                                     :: i,j
        integer                                     :: GroupRank
        
        !Begin-------------------------------------------------------------------

        write(*,*)
        write(*,*)'Opening HDF...'
        write(*,*)


        
        !Create HDF to read it 
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ, HDF5_READWRITE = HDF5_READWRITE , HDF5_CREATE = HDF5_CREATE)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5, trim(HDFFile), HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR00'

        if(OperateTwoHDF5) then
            !Open HDF5 file
            call ConstructHDF5 (ObjHDF5B, trim(HDFFileB), HDF5_READWRITE, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR00'
        end if

        call GetHDF5GroupExist (ObjHDF5, GroupName = "/Grid",    &
                                  Exist = Exist, STAT= STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR90'

        if (Exist) then
            call GetHDF5ArrayDimensions (HDF5ID = ObjHDF5,              &
                                        GroupName = trim("/Grid"),      &
                                        ItemName = trim("Bathymetry"),  &
                                        Imax = Size%IUB,                &
                                        Jmax = Size%JUB,                &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR90'
            
            Size%ILB = 1
            Size%JLB = 1
        end if

        nullify(WaterPoints2D)
        allocate (WaterPoints2D(1:Size%IUB, 1:Size%JUB))
        WaterPoints2D = FillValueReal

        nullify(PropertyField)
        allocate (PropertyField(1:Size%IUB, 1:Size%JUB))
        PropertyField = FillValueReal

        nullify(LatitudeGrid)
        allocate (LatitudeGrid(1:Size%IUB+1, 1:Size%JUB+1))
        LatitudeGrid = FillValueReal

        nullify(LongitudeGrid)
        allocate (LongitudeGrid(1:Size%IUB+1, 1:Size%JUB+1))
        LongitudeGrid = FillValueReal

                                   
        !Set grid
         call HDF5SetLimits (ObjHDF5, Size%ILB,             &
                            Size%IUB, Size%JLB,Size%JUB,   &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    & 
            stop 'ConvertHDF5Operator - HDF5Operator - ERR00'  
            
        call GetHDF5GroupNumberOfItems(ObjHDF5, "/Time", &
                                       NumberOfHDFInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR010'

        !verify instant given
        if (.not. AllInstants) then        
            if (Instant > NumberOfHDFInstants) then
                write(*,*)'HDF5 file instant required ("INSTANT" is a integer) does not exist in file :'//trim(HDFFile)
                stop 'ConvertHDF5Operator - HDF5Operator - ERR20'
            endif
        
        endif

        !Get HDF grid 
        call ConstructHorizontalGrid(ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR30'
        
        !GetConnectionX and ConnectionY
        !Load the grid corners coordinates XX and YY
        call GetCornersCoordinates(ObjHorizontalGrid, ConnectionX, ConnectionY, STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ConvertHDF5Operator - HDF5Operator - ERR31'

        !Gets Grid Latitude Longitude
        call GetLatitudeLongitude(ObjHorizontalGrid,                           &
                                  Latitude  = Latitude,                        &
                                  Longitude = Longitude,                       &
                                  STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR32'

        call GetGridCoordType(ObjHorizontalGrid, CoordType = CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'



        call HDF5SetLimits (ObjHDF5, 1,             &
                            Size%IUB, 1,Size%JUB,   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'


        
        call HDF5ReadData  (ObjHDF5, "/Grid",                        &
                             "Bathymetry",                               &
                             Array2D = PropertyField,                      &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'

        call HDF5SetLimits (ObjHDF5, 1,             &
                            Size%IUB+1, 1,Size%JUB+1,   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'


        
        call HDF5ReadData  (ObjHDF5, "/Grid",                        &
                             "Latitude",                               &
                             Array2D = LatitudeGrid,                      &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'

        call HDF5SetLimits (ObjHDF5, 1,             &
                            Size%IUB+1, 1,Size%JUB+1,   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'


        
        call HDF5ReadData  (ObjHDF5, "/Grid",                        &
                             "Longitude",                               &
                             Array2D = LongitudeGrid,                      &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'        
                                           
        call ConstructHDF5Out()
                              
        !search for each Groups from user list in HDF
        !Read block of parameters to extract
        
do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData,                           &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginParameter>',   &
                                        block_end       = '<EndParameter>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then        

                     call GetData(Property,                               &
                                  ObjEnterData, iflag,                    &
                                  SearchType   = FromBlock,               &
                                  keyword      = 'PROPERTY',              &
                                  ClientModule = 'HDF5Operator',        &
                                  STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0)          &
                            stop 'ReadDataFile - HDF5Operator - ERR060'  
                    
                        if (.not.CheckPropertyName(Property)) then
                            write(*,*)
                            write(*,*) 'The property name is not recognised by the model.'
                            write(*,*) 'ReadDataFile - HDF5Operator - WRN70'
                        endif           

                    ! Obtain parameter group
                    call GetData(HDFGroup,                                  &
                                 ObjEnterData, iflag,                        &
                                 SearchType   = FromBlock,                  &
                                 keyword      = 'HDF_GROUP',                &
                                 default      = "/Results/"//trim(Property),&
                                 ClientModule = 'HDF5Operator',           &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)           &
                    stop 'ConvertHDF5Operator - HDF5Operator - ERR080'
                            
                    !verify if group exists (not GetHDF5DataSetExist)
                    call GetHDF5GroupExist (ObjHDF5, GroupName = trim(HDFGroup),    &
                                              Exist = Exist, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR90'

                    if (Exist) then                    
                        AuxPosition = 1
                        !get field rank and dimensions
                        call GetHDF5GroupID(ObjHDF5, trim(HDFGroup),                            &
                                            AuxPosition, trim(HDFGroup),                        &
                                            Rank = GroupRank,                                   & 
                                            Dimensions = Dimensions,                            &
                                            STAT = STAT_CALL)                                
                        if (STAT_CALL .NE. SUCCESS_)                                            & 
                            stop 'ConvertHDF5Operator - HDF5Operator - ERR50'
                        
                        select case (GroupRank)
                            
                            case (2)
                                call Convert2DField(Dimensions, HDFGroup, Property)
                            case (3)
                                write(*,*)'Olny 2D cases for now. Use 2D properties groups'
                                stop 'ConvertHDF5Operator - HDF5Operator - ERR60'
                                !call Convert3DField(Dimensions, HDFGroup)
                            case default 
                            
                                write(*,*)'HDF5 group',  trim(HDFGroup) 
                                write(*,*) 'must be 2D or 3D.'
                                stop 'ConvertHDF5Operator - HDF5Operator - ERR70'
                            
                        end select 

                    else
                        write(*,*) 'The HDF group', trim(HDFGroup)
                        write(*,*) 'does not exist in the HDF file', trim(HDFFile)
                        stop 'ConvertHDF5Operator - HDF5Operator - ERR80'
                    endif
                 
                else cd2
                
                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR90' 
                    
                    exit do1                   
                
                endif cd2
                
            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConvertHDF5Operator - HDF5Operator - ERR100'
            else cd1
                stop 'ConvertHDF5Operator - HDF5Operator - ERR110'
            end if cd1
                
        enddo do1
    
    end subroutine ConvertHDF5Operator
    
    !---------------------------------------------------------------------------

    subroutine  ConstructHDF5Out ()
        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        integer                                     :: HDF5_READWRITE
        integer                                     :: HDF5_CREATE


        
        !Begin-------------------------------------------------------------------
 
         !Create HDF to read it 
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ, HDF5_READWRITE = HDF5_READWRITE , HDF5_CREATE = HDF5_CREATE)
       
        !Open HDF5 ouput file
        call ConstructHDF5 (ObjHDF5_Out, HDFFileOut, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR00'
                
        !Write the Horizontal Grid
        call WriteHorizontalGrid (ObjHorizontalGrid, ObjHDF5_Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'
          
        call HDF5SetLimits (ObjHDF5_Out, 1,             &
                            Size%IUB, 1,Size%JUB,   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'  
   
        call HDF5WriteData   (ObjHDF5_Out, "/Grid", "Bathymetry",             &
                              "m",                                     &
                              Array2D = PropertyField,                    &
                              STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'

        WaterPoints2D = 1                                       

        call HDF5WriteData   (ObjHDF5_Out, "/Grid", "WaterPoints2D",             &
                              "-",                                     &
                              Array2D = WaterPoints2D,                    &
                              STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConvertHDF5Operator - HDF5Operator - ERR33'
        
    end subroutine  ConstructHDF5Out

    
    !----------------------------------------------------------------------------
    
    subroutine Convert2DField (Dimensions, HDFGroup, Property)
        
        !Arguments---------------------------------------------------------------
        integer  , dimension(7)                     :: Dimensions
        character(len=PathLength)                   :: HDFGroup
        character(5)                                :: char_instant
       
        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL    
        integer                                     :: HDFInstant
        character(len=PathLength)                   :: Output
        character(len=StringLength)                 :: Property   
        character(StringLength)                     :: AuxChar    
        integer                                     :: Instant
        integer                                     :: NumberOfHDFInstants
        integer                                     :: Julday
        real                                        :: Hour,Minute,Second        
        real, dimension(:), pointer                 :: TimePtr
        real                                        :: RacingWithTheSun, Declination
        real                                        :: LatitudePI, LongitudePI
        real                                        :: HourAngle, SunHighAngle, GmtReference
        real                                        :: QSO, QSO10
        real                                        :: CloudCoverMinDay, CloudCoverNight
        
        integer                                     :: i,j
        real                                        :: AUX1, AUX2
        real                                        :: Topography = 100
        real                                        :: ATMTransmitivity = 1.0
        
        real, dimension(:,:), pointer               :: CloudCover
        real, dimension(:,:), pointer               :: LastRadiation
        real, dimension(:,:), pointer               :: SolarRadiation
        real, dimension(:,:), pointer               :: WindModulus
        real, dimension(:,:), pointer               :: WindVelocityX !FLAVIO SANTOS
        real, dimension(:,:), pointer               :: WindVelocityY !FLAVIO SANTOS
        real, dimension(:,:), pointer               :: AirTemperature
        real, dimension(:,:), pointer               :: RelativeHumidity

        real                                        :: NetRadiation
        real                                        :: SSVPC, psiconst
        real                                        :: LwradCorrection, Lwrad
        real                                        :: SVP, VP, SoilHeatFluxDensity
                 
        !Begin------------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid(s)...'
        write(*,*)
        
        !calculate Size2D
        Size%ILB = 1
        Size%IUB = Dimensions(1)
        Size%JLB = 1
        Size%JUB = Dimensions(2)

        call GetHDF5GroupNumberOfItems(ObjHDF5, "/Time", &
                                       NumberOfHDFInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR010'
        
        nullify(PropertyField)
        allocate (PropertyField(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        PropertyField = FillValueReal

        nullify(CloudCover)
        allocate (CloudCover(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        CloudCover = FillValueReal
        nullify(LastRadiation)
        allocate (LastRadiation(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        LastRadiation = FillValueReal
        nullify(SolarRadiation)
        allocate (SolarRadiation(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        SolarRadiation = FillValueReal
        nullify(WindModulus)
        allocate (WindModulus(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        WindModulus = FillValueReal
        nullify(AirTemperature)
        allocate (AirTemperature(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        AirTemperature = FillValueReal
        nullify(RelativeHumidity)
        allocate (RelativeHumidity(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        RelativeHumidity = FillValueReal


        nullify(PropertyFieldOld)
        allocate (PropertyFieldOld(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        PropertyFieldOld = 0

        if(OperateTwoHDF5) then
            nullify(PropertyField2)
            allocate (PropertyField2(1:Size%IUB, 1:Size%JUB))
            PropertyField2 = FillValueReal
            nullify(PropertyFieldOld2)
            allocate (PropertyFieldOld2(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
            PropertyFieldOld2 = 0
        end if

        nullify(TimePtr)
        allocate (TimePtr(1:6))
        TimePtr = FillValueReal

        do HDFInstant = 1, NumberOfHDFInstants 

            call HDF5SetLimits (ObjHDF5, Size%ILB,             &
                                Size%IUB, Size%JLB,Size%JUB,   &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'Convert2DField - HDF5Operator - ERR00'  

            call HDF5ReadData  (ObjHDF5, trim(HDFGroup),                       &
                                 trim(Property),                               &
                                 Array2D = PropertyField,                      &
                                 OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
            
            if(HDFInstant==1)then
                LastRadiation = 0
            else
                LastRadiation = SolarRadiation
            end if
            
            call HDF5ReadData  (ObjHDF5, '/Results/solar radiation',                       &
                                 'solar radiation',                               &
                                 Array2D = SolarRadiation,                      &
                                 OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
            
            if (.not. XYWind) then                
                
                call HDF5ReadData  (ObjHDF5, '/Results/wind modulus',            &
                                     'wind modulus',                             &
                                     Array2D = WindModulus,                      &
                                     OutputNumber = HDFInstant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
                
            else
                allocate (WindVelocityX(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
                allocate (WindVelocityY(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
                
                WindVelocityX = FillValueReal
                WindVelocityY = FillValueReal
                
                call HDF5ReadData  (ObjHDF5, '/Results/wind velocity X',           &
                                     'wind velocity X',                            &
                                     Array2D = WindVelocityX,                      &
                                     OutputNumber = HDFInstant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
                
                call HDF5ReadData  (ObjHDF5, '/Results/wind velocity Y',           &
                                     'wind velocity Y',                            &
                                     Array2D = WindVelocityY,                      &
                                     OutputNumber = HDFInstant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
                
                do i = 1, Size%IUB
                do j = 1, Size%JUB
                    WindModulus (i,j) = sqrt(WindVelocityX(i,j)**2 + WindVelocityY(i,j)**2)
                enddo
                enddo
                
            end if
            
            
            call HDF5ReadData  (ObjHDF5, '/Results/air temperature',                       &
                                 'air temperature',                               &
                                 Array2D = AirTemperature,                      &
                                 OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'

            call HDF5ReadData  (ObjHDF5, '/Results/relative humidity',                       &
                                 'relative humidity',                               &
                                 Array2D = RelativeHumidity,                      &
                                 OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'

            
            if(OperateTwoHDF5) then
                call HDF5SetLimits (ObjHDF5B, Size%ILB,             &
                                    Size%IUB, Size%JLB,Size%JUB,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                    stop 'Convert2DField - HDF5Operator - ERR00'  

                call HDF5ReadData  (ObjHDF5B, trim(HDFGroup),                       &
                                     trim(Property),                               &
                                     Array2D = PropertyField2,                      &
                                     OutputNumber = HDFInstant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
            endif

            call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
            
            call HDF5ReadData  (ObjHDF5, "Time",                       &
                                 "Time",                               &
                                 Array1D = TimePtr,                      &
                                 OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
            
            Hour = TimePtr(4) 
            Minute = TimePtr(5) 
            Second = TimePtr(6)
            
            call SetDate (Now,TimePtr(1),TimePtr(2),TimePtr(3),Hour,Minute,Second)

            call JulianDay  (Now, Julday)

            Hour = Hour + Minute/60. + Second/3600.

            !Sun "Racing"
            RacingWithTheSun = RacingWithTheSun_(JulDay)

            !Declination of the Sun
            Declination = SunDeclination_(JulDay)

            call HDF5SetLimits (ObjHDF5_Out, 1,             &
                                Size%IUB, 1,Size%JUB,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR33'

            if(OperateTwoHDF5) then
                do i = 1, Size%IUB
                    do j = 1, Size%JUB
                        AUX1 = PropertyField (i,j) - PropertyFieldOld (i,j)
                        AUX2 = PropertyField2(i,j) - PropertyFieldOld2(i,j)
                        if(AUX1<0) AUX1=0
                        if(AUX2<0) AUX2=0
                        PropertyField(i,j) = (AUX1+AUX2) * 10
                    end do                 
                end do 
            Else
                !PropertyField = 7.0012*PropertyField + 0.2978  
                do i = 1, Size%IUB
                    do j = 1, Size%JUB
                    
                    
                    LatitudePI  = LatitudeGrid (i, j) * PI / 180.0
                    LongitudePI = LongitudeGrid(i, j) 
                    GmtReference = 0
                    !Hour angle 
                    HourAngle = HourAngle_ (Hour, LongitudePI, GmtReference, RacingWithTheSun)
                 
                    !use a sunset angle of 10 as safety factor
                    SunHighAngle =  Asin( sin(LatitudePI) * sin(Declination)            +   &
                                   cos(LatitudePI) * cos(Declination) * cos(HourAngle)) *   &
                                   180./PI                    

                    !New implementation because cloud cover should be limited at night (and in sunrise and sunset)
                    !Because with the old formulation gives very low cloud cover values at this times (order of 0.001) 
                    !that are not consistent with bibliography (cloud cover values from 0.3 to 1 and night around 0.5).
                    !In the above old implementation it seems that code exists to correct it ("if (QSO==0.0)") 
                    !but will never correct it because the situation SunHighAngles > 10 is not at night.
                    !So compute normal radiation, start with night, than low angles, than normal day
                    
                    !Clear day radiation
                    QSO = TOARadiation (LatitudePI, Declination, HourAngle, Julday)
                    CloudCoverNight          = 0.595
                    CloudCoverMinDay         = 0.3
                    
                    !night - impose cloud cover
                    if (QSO == 0.0) then
                        CloudCover (i,j) = CloudCoverNight
                   
                    !day
                    else
                    
                        !low positive angles - sunrise or sunset - use last observed positive radiation 
                        !and 10º TOA radiation.
                        !Testing the formulation this part still creates low cloud cover at sunrise because   
                        !the last radiation saved was the low value before night.
                        if (SunHighAngle <= 10.) then
                            
                            !It's almost nigth use sunset angle of 10º                      
                            QSO10    = TOARadiation (LatitudePI  = LatitudePI ,    &
                                                     Declination = Declination,    &
                                                     Julday      = Julday )     
                                                                   
                            CloudCover (i,j) = LastRadiation(i,j) / QSO10 
                            
                        !normal day - use observed radiation and TOA radiation
                        else
                            CloudCover (i,j) = CloudCover (i,j) / QSO
                        endif

                        !Correct daily values only
                        !usually near sunrise values computed are very low (near 10º sun angle or more)
                        !so min value for day is important to correct it
                        if (CloudCover (i,j) < CloudCoverMinDay) then

                            CloudCover (i,j) = CloudCoverMinDay

                        elseif (CloudCover (i,j) > 1) then

                            CloudCover (i,j) = 1.0
                        
                        endif

                    endif
                    
                    ATMTransmitivity = CloudCover (i,j)
                    
                    !Calculate Psicrometric constant
                    !Calculation of the atmospheric pressure based on the heigth simplification of the ideal gas law
                    psiconst    = 0.665E-3 * (101.3 * ( (293.-0.0065 * Topography ) / 293.) **5.26) ![kPa / ºC]
                
                    !Saturation Vapour pressura [kPa - Eqn 11]
                    SVP         = 0.6108 * exp (17.27 * AirTemperature(i,j) / (AirTemperature(i,j) + 237.3))
                
                    ![kPa]
                    VP          = SVP * RelativeHumidity(i, j)  

                    !Calculates the Slope [kPa/ºC]
                    SSVPC       =  4098.* SVP / (AirTemperature(i,j) + 237.3)**2.0 ![kPa / ºC]
                
                    !StefanBoltzmann          = 5.669e-08     ![W/m2/K4]
                    LwradCorrection =   (0.34 - 0.14 * VP **(0.5)) * (1.35 * ATMTransmitivity  - 0.35)
                    Lwrad           =   5.669e-08 * (AirTemperature(i,j) + 273.15)** 4. * LwradCorrection   ![W / m2]
                    
                    !Calculation of net radiation (0.23 is the reference albedo)
                    NetRadiation    = (1-0.23) * SolarRadiation(i, j) - Lwrad            ![W / m2]   
                
                    !Converts Netradiation into MJ/m2/hour (requested by the formular below)
                    !1W = J / s => 1MJ/hour
                    NetRadiation    = NetRadiation /1.e6 * 3600.
                                    
                    if (NetRadiation .GE. 0) then
                        SoilHeatFluxDensity = 0.1 * NetRadiation
                    else
                        SoilHeatFluxDensity = 0.5 * NetRadiation
                    end if
                    
                    !wind modulus modification - new implementation (FLAVIO SANTOS)
                    !most wind modulus datasets present wind modulus at a height of 10m
                    !According to FAO: "Wind speeds measured at different heights above the soil surface are different. 
                    !Surface friction tends to slow down wind passing over it. Wind speed is slowest at the surface and 
                    !increases with height. For this reason anemometers are placed at a chosen standard height, i.e.
                    !, 10 m in meteorology and 2 or 3 m in agrometeorology. For the calculation of evapotranspiration, 
                    !wind speed measured at 2 m above the surface is required. To adjust wind speed data obtained from 
                    !instruments placed at elevations other than the standard height of 2m, a logarithmic wind speed profile 
                    !may be used for measurements above a short grassed surface"
                    !z0 (roughness coef) roughly = 0.017
                    !u2 = u1 * ln(z2/z0)/ln(z1/z0) if wind final height is to be implemented too
                    WindModulus(i,j) = WindModulus(i,j) * (4.87 / LOG((67.8 * WindHeight) - 5.42))

                
                    !FAO ET0 - mm/hour
                    !http://www.fao.org/docrep/X0490E/x0490e08.htm#calculation%20procedure (Hourly time step)
                    PropertyField(i, j) = (0.408 * SSVPC * (NetRadiation - SoilHeatFluxDensity) +  &
                                            psiconst * 37. /(AirTemperature(i, j) + 273.)  *        &
                                            WindModulus(i,j) * (SVP- VP)) /                         & 
                                            (SSVPC + psiconst * (1. + 0.34 * WindModulus(i,j)))
                    !mm/3hours
                    PropertyField(i, j) = PropertyField(i, j) * TimeInterval
                    end do                 
                end do                 
            end if

            call HDF5WriteData   (ObjHDF5_Out,                            &
                                  '/Results/reference evapotranspiration',          &
                                  'reference evapotranspiration',                       &
                                  "-",                      &
                                  Array2D = PropertyField,                     &
                                  OutputNumber = HDFInstant,                   &
                                  STAT = STAT_CALL)

            call HDF5SetLimits  (ObjHDF5_Out, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'

            call HDF5WriteData  (ObjHDF5_Out, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                 Array1D = TimePtr, OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'

            call HDF5SetLimits (ObjHDF5, Size%ILB,             &
                                Size%IUB, Size%JLB,Size%JUB,   &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'Convert2DField - HDF5Operator - ERR00'  

            call HDF5ReadData  (ObjHDF5, trim(HDFGroup),                       &
                                 trim(Property),                               &
                                 Array2D = PropertyFieldOld,                      &
                                 OutputNumber = HDFInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
 
            if(OperateTwoHDF5) then
                call HDF5SetLimits (ObjHDF5B, Size%ILB,             &
                                    Size%IUB, Size%JLB,Size%JUB,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    & 
                    stop 'Convert2DField - HDF5Operator - ERR00'  

                call HDF5ReadData  (ObjHDF5B, trim(HDFGroup),                       &
                                     trim(Property),                               &
                                     Array2D = PropertyFieldOld2,                      &
                                     OutputNumber = HDFInstant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Convert2DField - HDF5Operator - ERR10'
            endif
                                        
        enddo

        

        
        deallocate (PropertyField)
        deallocate (TimePtr)
        deallocate(CloudCover)
        deallocate(LastRadiation)
        deallocate(SolarRadiation)
        deallocate(WindModulus)
        deallocate(AirTemperature)
        deallocate(RelativeHumidity)
        deallocate(PropertyFieldOld)

    end subroutine Convert2DField
    
    !--------------------------------------------------------------------------

    subroutine WriteGrid2D(PropertyField, OutPut)
        
        !Arguments-------------------------------------------------------------
         real, dimension(:,:  ), pointer             :: PropertyField 
         character(len=PathLength)                   :: Output

        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL !, UnitID, i, j
        !type(T_Size2D)              :: Size2D
        
        !Begin-----------------------------------------------------------------

        
        !Size2D%ILB = Size%ILB
        !Size2D%IUB = Size%IUB
        !Size2D%JLB = Size%JLB
        !Size2D%JUB = Size%JUB


        call WriteGridData  (FileName       = Output,                                       &
                             ConnectionX    = ConnectionX,                                  &
                             ConnectionY    = ConnectionY,                                  &
                             COMENT1        = 'Grid Data created from HDF file',            &
                             COMENT2        = trim(HDFFile),                                &
                             WorkSize       = Size,                                         &
                             CoordType      = CoordType,                                    &
                             Xorig          = ConnectionX(1,1),                             &
                             Yorig          = ConnectionY(1,1),                             &
                             Zone           = -99,                                          &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Latitude,                                     &
                             Longitude      = Longitude,                                    &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= PropertyField,                                &
                             STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleWRFFormat - ERR01'



    end subroutine WriteGrid2D

    !------------------------------------------------------------------------

    subroutine ConstructDSName (Name, OutputNumber, AuxChar)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Name
        integer                                     :: OutputNumber
        character(len=*)                            :: AuxChar

        !Local-----------------------------------------------------------------
        character(StringLength)                     :: AuxNum

        write(AuxNum, fmt=*)OutputNumber

        if     (OutputNumber < 10     ) then
            AuxChar = trim(adjustl(Name))//"_0000"//trim(adjustl(AuxNum))
        elseif (OutputNumber < 100    ) then
            AuxChar = trim(adjustl(Name))//"_000" //trim(adjustl(AuxNum))
        elseif (OutputNumber < 1000   ) then
            AuxChar = trim(adjustl(Name))//"_00"  //trim(adjustl(AuxNum))
        elseif (OutputNumber < 10000  ) then
            AuxChar = trim(adjustl(Name))//"_0"   //trim(adjustl(AuxNum))
        else
            AuxChar = trim(adjustl(Name))//"_"    //trim(adjustl(AuxNum))
        endif

    end subroutine ConstructDSName

    !--------------------------------------------------------------------------

    !This function accounts for disturbances in the earth’s rotation rate that affect 
    !the time the suns takes to go through the longitude differences. 
    !Taken from "Evapotranspiration Technical Manual" eq. 53 and eq. 54 but reference not found.
    !It adds maximum 10 min or takes maximum 15 min to the hour depending on day of year.
    real function RacingWithTheSun_ (JulDay)

        !Arguments-------------------------------------------------------------
        integer                                     :: JulDay

        !Local-----------------------------------------------------------------
        real                                        :: Aux

        Aux = 2 * PI * (JulDay - 81)/364. 
        RacingWithTheSun_ = 0.1645 * sin(2*Aux) - 0.1255 * cos(Aux) - 0.025 * Aux 

    end function RacingWithTheSun_

    !--------------------------------------------------------------------------
    
    !This function computes the declination of the Sun which is the angle between &
    !the rays of the sun and the plane of the earth's equator. 
    !Possibly taken from Deas, M.L. and Lowney C.L. - "Water Temperature Modelling Review" Sept. 2000.
    real function SunDeclination_(JulDay)

        !Arguments-------------------------------------------------------------
        integer                                     :: JulDay

        !Local-----------------------------------------------------------------

        !Sun Declination
        SunDeclination_ = 23.45 * cos((172.0 - JulDay) * 2.0 * PI / 365.0) * PI / 180.0

    end function SunDeclination_

    !This function computes radian angles according to the (solar) hour for sun height computation.
    !At 12h, angle is zero (maximum co-sine in sun height computation) and increases towards sunrise & 
    !and sunset (co-sine in sun height computation decreases).
    !Reference? Possibly adapted from Deas, M.L. and Lowney C.L. - "Water Temperature Modelling Review" Sept. 2000.
    real function HourAngle_ (Hour, LongitudePI, GmtReference, RacingWithTheSun)

        !Arguments-------------------------------------------------------------
        real                                        :: LongitudePI, GmtReference
        real                                        :: RacingWithTheSun, Hour

        !Local-----------------------------------------------------------------
        real                                        :: HourC            
        
        !Corrected Hour for longitude, timezone and small disturbances (Racing with the sun)
        ! h   =  h   +    degrees    / deg/h -       h       +       h
        HourC = Hour + (LongitudePI / 15.0) - GmtReference + RacingWithTheSun 
        
        !Hour angle (to_change passar para SCT)
        if (HourC .LT. 12) then
            HourAngle_ = (HourC + 12.0) * PI / 12.0 
        else
            HourAngle_ = (HourC - 12.0) * PI / 12.0 
        end if

    end function HourAngle_
    !--------------------------------------------------------------------------

    function TOARadiation(LatitudePI, Declination, HourAngle, Julday)
        real TOARadiation

        !Arguments-------------------------------------------------------------
        integer , intent(IN)                :: Julday   !Julday = 1 -> 1st of January
        real    , intent(IN)                :: LatitudePI
        real    , intent(IN)                :: Declination
        real    , intent(IN), optional      :: HourAngle
        
        !Local-----------------------------------------------------------------
        real    :: SunHigh
        real    :: ROrbit 

        !--------------------------------------------------------------------------

        !Sun high  
        if (present (HourAngle))  then    
            SunHigh = sin(LatitudePI) * sin(Declination) + cos(LatitudePI) * cos(Declination) * cos(HourAngle)        
        else
            SunHigh = sin (10 * PI / 180.)
        endif
        !Atmosphere radiation
        if (SunHigh .LT. 0.0) then     !night
            
            TOARadiation = 0.0

        else                           !day
            
            ROrbit    = 1.0 + 0.017 * cos((186.0 - JulDay) * 2.0 * PI / 365.0)
            !Top atmosphere sun radiation
            TOARadiation       = KSun * SunHigh / ROrbit ** 2.0
        end if

        !----------------------------------------------------------------------

    end function TOARadiation    
    !--------------------------------------------------------------------------
    
    subroutine KillHDF5Operator
    
        !Local---------------------------------------------------------------
        integer                                  :: STAT_CALL
        !Begin---------------------------------------------------------------

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillHDF5Operator - HDF5Operator - ERR00'

        call UngetHorizontalGrid(ObjHorizontalGrid, ConnectionX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillHDF5Operator - HDF5Operator - ERR00'
        
        call UngetHorizontalGrid(ObjHorizontalGrid, ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillHDF5Operator - HDF5Operator - ERR00'
        
        call KillHorizontalGrid(ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillHDF5Operator - HDF5Operator - ERR00'

        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillHDF5Operator - HDF5Operator - ERR00'

        call KillHDF5 (ObjHDF5_Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillHDF5Operator - HDF5Operator - ERR00'
    
    end subroutine KillHDF5Operator
    
    !----------------------------------------------------------------------------
    
    subroutine ShutHDF5Operator
    
        call date_and_time(Values = F95Time)
        call SetDate (FinalSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                    float(F95Time(3)), float(F95Time(5)), &
                                    float(F95Time(6)), float(F95Time(7))+ &
                                    F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = FinalSystemTime - InitialSystemTime
        call ShutdownMohid ("HDF5Operator", ElapsedSeconds, TotalCPUTime)
    
    end subroutine ShutHDF5Operator
    
    !----------------------------------------------------------------------------
    
end program ETCalculator

