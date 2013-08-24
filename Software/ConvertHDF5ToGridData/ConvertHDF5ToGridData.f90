!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Tools
! PROGRAM       : HDF5ToGridData
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : August 2013
! REVISION      : David Brito
! DESCRIPTION   : To create grid data files from HDF5
!
!------------------------------------------------------------------------------

program HDF5ToGridData

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
    

    real, dimension(:, :), pointer        :: ConnectionX         => null()
    real, dimension(:, :), pointer        :: ConnectionY         => null()
    real, dimension(:, :), pointer        :: PropertyField
    real                                  :: Latitude
    real                                  :: Longitude
    integer                               :: CoordType

     type (T_Time)                        :: InitialSystemTime, FinalSystemTime
     integer, dimension(8)                :: F95Time
     real                                 :: ElapsedSeconds, TotalCPUTime

     character(PathLength)                :: DataFile  = 'HDF5ToGridData.dat'

     character(PathLength)                :: HDFFile
!     logical                              :: InstantOnName
     logical                              :: AllInstants
     integer                              :: NumberOfHDFInstants
!     type (T_Time)                        :: TimeInstant
     integer                              :: Instant
!     integer                              :: NumberOfInstants
!     character(1024)                      :: BaseGroup
     integer                              :: GroupRank
!     type (T_PropertyID)                  :: Property
     type (T_Size2D)                      :: Size
!     type (T_Size2D)                      :: WorkSize
!     real, dimension(:,:), pointer        :: DataSet
!     integer, dimension(:,:), pointer     :: PointsToFill
!     type (T_Time)                        :: BeginTime, EndTime, Current
!     real                                 :: FillValue
!     integer                              :: SelectionType
!     logical                              :: FindStartInHDF
!     logical                              :: FindEndInHDF


     integer                              :: ObjHDF5           = 0
     integer                              :: ObjHorizontalGrid = 0
!     integer                              :: ObjGridData       = 0
!     integer                              :: ObjFillMatrix     = 0
!     integer                              :: ObjTime           = 0
     integer                              :: ObjEnterData      = 0
  
    
    !Begin----------------------------------------------------------------------
    
    call StartHDF5ToGridData


    call ReadDataFile

    call ConvertHDF5ToGridData
    
    call KillHDF5ToGridData
    
    
    call ShutHDF5ToGridData

    
    contains
    
    
    !----------------------------------------------------------------------------
    
    subroutine StartHDF5ToGridData
    
        call date_and_time(Values = F95Time)
        call StartUpMohid("HDF5ToGridData")
        call SetDate(InitialSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                     float(F95Time(3)), float(F95Time(5)), &
                                     float(F95Time(6)), float(F95Time(7))+ &
                                     F95Time(8)/1000.)
    end subroutine StartHDF5ToGridData
    
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
         if (STAT_CALL /= SUCCESS_) stop 'HDF5ToGridData - ERR000'

         call GetData(HDFFile,                             &
                      ObjEnterData, iflag,                 &
                      SearchType   = FromFile,             &
                      keyword      = 'HDF_FILE',           &
                      ClientModule = 'HDF5ToGridData',     &
                      STAT         = STAT_CALL)
         if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5ToGridData - ERR010'

        !Verifies if file exists
        inquire(FILE = HDFFile, EXIST = Exist)
        if (.not. Exist) then
            write(*,*)'HDF5 file does not exist:'//trim(HDFFile)
            stop 'ReadDataFile - HDF5ToGridData - ERR20'
        endif

        call GetData(AllInstants,                         &
                  ObjEnterData, iflag,                    &
                  SearchType   = FromFile,                &
                  keyword      = 'ALL_INSTANTS',          &
                  Default      = .false.,                 &
                  ClientModule = 'HDF5ToGridData',        &
                  STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5ToGridData - ERR050'

        if (.not. AllInstants) then
        
             call GetData(Instant,                                &
                          ObjEnterData, iflag,                    &
                          SearchType   = FromFile,                &
                          keyword      = 'INSTANT',               &
                          Default      = 1,                       &
                          ClientModule = 'HDF5ToGridData',        &
                          STAT         = STAT_CALL)
             if (STAT_CALL .NE. SUCCESS_) stop 'ReadDataFile - HDF5ToGridData - ERR060'            
            
        endif   
        
    end subroutine ReadDataFile
    
    !---------------------------------------------------------------------------

    subroutine ConvertHDF5ToGridData

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: HDF5_READ    
        logical                                     :: Exist   
        character(len=PathLength)                   :: HDFGroup
        integer  , dimension(7)                     :: Dimensions
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: AuxPosition
        character(len=StringLength)                 :: Property
        
        !Begin-------------------------------------------------------------------

        write(*,*)
        write(*,*)'Opening HDF...'
        write(*,*)

        
        !Create HDF to read it 
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (ObjHDF5, trim(HDFFile), HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR00'

        call GetHDF5GroupNumberOfItems(ObjHDF5, "/Time", &
                                       NumberOfHDFInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR010'

        !verify instant given
        if (.not. AllInstants) then        
            if (Instant > NumberOfHDFInstants) then
                write(*,*)'HDF5 file instant required ("INSTANT" is a integer) does not exist in file :'//trim(HDFFile)
                stop 'ConvertToGridData - HDF5ToGridData - ERR20'
            endif
        
        endif

        !Get HDF grid 
        call ConstructHorizontalGrid(ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR30'
        
        !GetConnectionX and ConnectionY
        !Load the grid corners coordinates XX and YY
        call GetCornersCoordinates(ObjHorizontalGrid, ConnectionX, ConnectionY, STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ConvertToGridData - HDF5ToGridData - ERR31'

        !Gets Grid Latitude Longitude
        call GetLatitudeLongitude(ObjHorizontalGrid,                           &
                                  Latitude  = Latitude,                        &
                                  Longitude = Longitude,                       &
                                  STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR32'

        call GetGridCoordType(ObjHorizontalGrid, CoordType = CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR33'

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
                                  ClientModule = 'HDF5ToGridData',        &
                                  STAT         = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0)          &
                            stop 'ReadDataFile - HDF5ToGridData - ERR060'  
                    
                        if (.not.CheckPropertyName(Property)) then
                            write(*,*)
                            write(*,*) 'The property name is not recognised by the model.'
                            write(*,*) 'ReadDataFile - HDF5ToGridData - WRN70'
                        endif           

                    ! Obtain parameter group
                    call GetData(HDFGroup,                                  &
                                 ObjEnterData, iflag,                        &
                                 SearchType   = FromBlock,                  &
                                 keyword      = 'HDF_GROUP',                &
                                 default      = "/Results/"//trim(Property),&
                                 ClientModule = 'HDF5ToGridData',           &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)           &
                    stop 'ReadDataFile - HDF5ToGridData - ERR080'
                            
                    !verify if group exists (not GetHDF5DataSetExist)
                    call GetHDF5GroupExist (ObjHDF5, GroupName = trim(HDFGroup),    &
                                              Exist = Exist, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR90'

                    if (Exist) then                    
                        AuxPosition = 1
                        !get field rank and dimensions
                        call GetHDF5GroupID(ObjHDF5, trim(HDFGroup),                            &
                                            AuxPosition, trim(HDFGroup),                        &
                                            Rank = GroupRank,                                   & 
                                            Dimensions = Dimensions,                            &
                                            STAT = STAT_CALL)                                
                        if (STAT_CALL .NE. SUCCESS_)                                            & 
                            stop 'ConvertToGridData - HDF5ToGridData - ERR50'
                        
                        select case (GroupRank)
                            
                            case (2)
                                call Convert2DField(Dimensions, HDFGroup, Property)
                            case (3)
                                write(*,*)'Olny 2D cases for now. Use 2D properties groups'
                                stop 'ConvertToGridData - HDF5ToGridData - ERR60'
                                !call Convert3DField(Dimensions, HDFGroup)
                            case default 
                            
                                write(*,*)'HDF5 group',  trim(HDFGroup) 
                                write(*,*) 'must be 2D or 3D.'
                                stop 'ConvertToGridData - HDF5ToGridData - ERR70'
                            
                        end select 

                    else
                        write(*,*) 'The HDF group', trim(HDFGroup)
                        write(*,*) 'does not exist in the HDF file', trim(HDFFile)
                        stop 'ConvertToGridData - HDF5ToGridData - ERR80'
                    endif
                 
                else cd2
                
                    call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR90' 
                    
                    exit do1                   
                
                endif cd2
                
            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConvertToGridData - HDF5ToGridData - ERR100'
            else cd1
                stop 'ConvertToGridData - HDF5ToGridData - ERR110'
            end if cd1
                
        enddo do1
    
    end subroutine ConvertHDF5ToGridData
    
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
        !Begin------------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid(s)...'
        write(*,*)
        
        !calculate Size2D
        Size%ILB = 1
        Size%IUB = Dimensions(1)
        Size%JLB = 1
        Size%JUB = Dimensions(2)
        
        nullify(PropertyField)
        allocate (PropertyField(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
        PropertyField = FillValueReal

        call HDF5SetLimits (ObjHDF5, Size%ILB,             &
                            Size%IUB, Size%JLB,Size%JUB,   &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                    & 
            stop 'ConvertToGridData - HDF5ToGridData - ERR00'  
        
        if (AllInstants) then
            
            do HDFInstant = 1, NumberOfHDFInstants 

                call HDF5ReadData  (ObjHDF5, trim(HDFGroup),                       &
                                     trim(Property),                               &
                                     Array2D = PropertyField,                      &
                                     OutputNumber = HDFInstant, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR10'
                
                call ConstructDSName(Property,HDFInstant,AuxChar)              
                Output = trim(AuxChar) // ".dat" 
                
                call WriteGrid2D(PropertyField, Output)
            
            enddo
        
        else
            !one defined instant
            call HDF5ReadData  (ObjHDF5, trim(HDFGroup),                       &
                                 trim(Property),                               &
                                 Array2D = PropertyField,                      &
                                 OutputNumber = Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConvertToGridData - HDF5ToGridData - ERR90'                   
            
            call ConstructDSName(Property,Instant,AuxChar)              
            Output = trim(AuxChar) // ".dat" 
            
            call WriteGrid2D(PropertyField, Output)
            
        endif
        
        deallocate (PropertyField)
        
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

    subroutine KillHDF5ToGridData
    
        !Local---------------------------------------------------------------
        integer                                  :: STAT_CALL
        !Begin---------------------------------------------------------------

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillIHRadarFormat - HDF5ToGridData - ERR000'

        call UngetHorizontalGrid(ObjHorizontalGrid, ConnectionX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR10'
        
        call UngetHorizontalGrid(ObjHorizontalGrid, ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR20'
        
        call KillHorizontalGrid(ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR50'

        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillIHRadarFormat - ModuleIHRadarFormat - ERR60'
    
    end subroutine KillHDF5ToGridData
    
    !----------------------------------------------------------------------------
    
    subroutine ShutHDF5ToGridData
    
        call date_and_time(Values = F95Time)
        call SetDate (FinalSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                    float(F95Time(3)), float(F95Time(5)), &
                                    float(F95Time(6)), float(F95Time(7))+ &
                                    F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = FinalSystemTime - InitialSystemTime
        call ShutdownMohid ("HDF5ToGridData", ElapsedSeconds, TotalCPUTime)
    
    end subroutine ShutHDF5ToGridData
    
    !----------------------------------------------------------------------------
    
end program HDF5ToGridData

