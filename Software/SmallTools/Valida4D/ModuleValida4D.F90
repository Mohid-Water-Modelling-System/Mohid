!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Valida4D
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : June 2010
! REVISION      : Paulo Leitão - v4.0
! DESCRIPTION   : Module to serve as Valida4D to create new modules
!
!------------------------------------------------------------------------------


Module ModuleValida4D

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5
    use ModuleStopWatch,      only : CreateWatchGroup, KillWatchGroup    
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap
!    use ModuleValida4D
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructValida4D


    !Selector
                     
    
    !Modifier
    public  :: ModifyValida4D

    !Destructor
    public  :: KillValida4D                                                     

    !Management
    
    !Interfaces----------------------------------------------------------------


    !Types---------------------------------------------------------------------

    type (T_Time)                           :: InitialSystemTime, FinalSystemTime
    real                                    :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                   :: F95Time

    type     T_InterpolTime
        real, dimension(:, :   ),  pointer          :: Values2D, NextValues2D, PrevValues2D
        real, dimension(:, :, :),  pointer          :: Values3D, NextValues3D, PrevValues3D       
    end type T_InterpolTime

    type      T_External_Var
        integer, dimension(:, :   ), pointer        :: KFloor_Z, WaterPoints2D
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real   , dimension(:, :   ), pointer        :: DZX, DZY, DUX, DVY
    end type  T_External_Var
    

    
    private :: T_HDF5Files
    type       T_HDF5Files
        character (len = PathLength)                :: Name
        type (T_Time)                               :: StartTime, EndTime
        integer                                     :: NumberOfInstants
        integer                                     :: ObjHDF5 = 0
    end type   T_HDF5Files    
    
    type      T_VerticalZ
        type(T_PropertyID)                          :: ID      
!        integer, dimension(:      ), pointer        :: ObjFillMatrix
        integer, dimension(:, :, :), pointer        :: Mapping
        type (T_InterpolTime)                       :: InterpolTime     
    end type  T_VerticalZ

    private :: T_Properties
    type       T_Properties
        type(T_PropertyID)                          :: ID
        integer                                     :: Column
        real, dimension(:), pointer                 :: ValueTable
        real, dimension(:), pointer                 :: ValueHDF5
        logical                                     :: ThreeD
        type (T_InterpolTime)                       :: InterpolTime
!        integer, dimension(:), pointer              :: ObjFillMatrix                 
    end type   T_Properties        
    
    private :: T_Valida4D
    type       T_Valida4D
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: WorkSize2D
        character (len = PathLength)                :: BathymetryFile
        character (len = PathLength)                :: GeometryFile
        character (len = PathLength)                :: InputTable
        character (len = PathLength)                :: OutputTable        
        integer                                     :: TableColumns, PropNumber, HDF5Number, TableValues
        integer                                     :: Tcolumn, Xcolumn, Ycolumn, Zcolumn
        type(T_Time)                                :: InitialDate
        logical, dimension(:), pointer              :: NullValue        
        real                                        :: Zmax
        
        real,    dimension(:), pointer              :: T, X, Y, Z, PercI, PercJ
        integer, dimension(:), pointer              :: i, j
        character (len = StringLength), dimension(:), pointer :: StationName
      
        type (T_HDF5Files),  dimension(:), pointer  :: HDF5Files
        integer                                     :: NextInstant, PrevInstant
        type (T_Time)                               :: NextTime,PrevTime

        type (T_Properties), dimension(:), pointer  :: Properties
        type (T_VerticalZ)                          :: VerticalZ
        
        real,    dimension(:, :   ),  pointer       :: LayerValues
        integer, dimension(:, :   ),  pointer       :: FullMatrix2D
        
        integer, dimension(:, :, :),  pointer       :: FullMatrix3D        
        
        type (T_External_Var)                       :: External_Var
        
        integer                                     :: ObjTime           = 0
        integer                                     :: ObjEnterData      = 0
        integer                                     :: ObjHorizontalGrid = 0
        integer                                     :: ObjBathymetry     = 0
        integer                                     :: ObjHorizontalMap  = 0
        integer                                     :: ObjMap            = 0 
        integer                                     :: ObjGeometry       = 0
        
        type(T_Valida4D), pointer                   :: Next
    end type  T_Valida4D

    !Global Module Variables
    type (T_Valida4D), pointer                      :: FirstObjValida4D
    type (T_Valida4D), pointer                      :: Me

    integer                                         :: mValida4D_ = 0 !just to compile

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructValida4D()

        !Arguments---------------------------------------------------------------

        !External----------------------------------------------------------------


        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------
        write(*,*)
        write(*,*)'Running Valida4D...'
        write(*,*)
        
        call StartUpMohid("Valida4D")

        call StartCPUTime
        
        allocate(Me)
      
        call ReadOptionsFile
        
        call ReadInputTable
        
        call ConstructSubModules
        
        call AllocateMatrixes
                
 

        !----------------------------------------------------------------------

    end subroutine ConstructValida4D
 
    !--------------------------------------------------------------------------
    

    subroutine ReadOptionsFile

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        character(PathLength)                       :: WatchFile        
        
        !Begin-----------------------------------------------------------------
        call ConstructEnterData (Me%ObjEnterData, "InputValida4D.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptionsFile - ModuleValida4D - ERR10'
        
        call GetData(WatchFile,                                             &
                     Me%ObjEnterData, iflag,                                &
                     SearchType   = FromFile,                               &
                     keyword      = 'OUTWATCH',                             &
                     ClientModule = 'Valida4D',                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptionsFile - ModuleValida4D - ERR02'
        if (iflag == 0)then
            MonitorPerformance  = .false.
        else
            STAT_CALL = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        endif
                

        call ReadOptions
        
        call ReadProperties
        
        call ReadHDF5FileName
        
        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptionsFile - ModuleValida4D - ERR20'
        
        
    end subroutine ReadOptionsFile

    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        logical                                     :: exist

        !Begin-----------------------------------------------------------------
        
        call GetData(Me%BathymetryFile,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BATHYM_FILENAME',                                  &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR10'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR20'

        inquire(FILE = Me%BathymetryFile, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Bathym file does not exist'
            stop 'ReadOptions - ModuleValida4D - ERR30'
        endif
    
        call GetData(Me%InputTable,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'INPUT_TABLE',                                      &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR50'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR60'

        call GetData(Me%Tcolumn,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'T_COLUMN',                                         &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR70'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR80'

        call GetData(Me%Xcolumn,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'X_COLUMN',                                         &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR90'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR100'

        call GetData(Me%Ycolumn,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'Y_COLUMN',                                         &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR100'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR110'

        call GetData(Me%Zcolumn,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'Z_COLUMN',                                         &
                     default      = FillValueInt,                                       &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR130'
            

        call GetData(Me%GeometryFile,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GEOMETRY_FILENAME',                                &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR150'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR160'
        

        inquire(FILE = Me%GeometryFile, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Geometry file does not exist'
            stop 'ReadOptions - ModuleValida4D - ERR170'
        endif
        
        call GetData(Me%OutputTable,                                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUT_TABLE',                                     &
                     ClientModule = 'ModuleValida4D',                                   &
                     STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleValida4D - ERR180'
        if (iflag     == 0       ) stop 'ReadOptions - ModuleValida4D - ERR190'        

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ReadProperties

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber, iP, iflag
        logical                                     :: BlockFound
        logical                                     :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        Me%PropNumber = 0
        
        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleValida4D - ERR10'
        

        ! Check number of Properties
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = '<beginproperty>',            & 
                                        block_end       = '<endproperty>',              &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then  

                    Me%PropNumber = Me%PropNumber + 1                                                        

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                                  & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleValida4D - ERR20'

                    exit do1
                end if cd2

            end if cd1

        end do do1

        if (Me%PropNumber == 0) then                                            
            write(*,*) 'No property block is indicated in input file. '
            stop 'ReadProperties - ModuleValida4D - ERR30'
        end if
        
        allocate(Me%Properties(Me%PropNumber))

        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleValida4D - ERR40'


        ! Read Properties
do2 :   do  iP = 1, Me%PropNumber
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = '<beginproperty>',            &
                                        block_end       = '<endproperty>',              &
                                        BlockFound      = BlockFound)

            call ConstructPropertyID  (Me%Properties(iP)%ID, Me%ObjEnterData, FromBlock)

            call GetData(Me%Properties(iP)%Column,                                      &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'COLUMN',                                       &
                         ClientModule = 'ModuleValida4D',                               &
                         STAT         = STAT_CALL)                                          
            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleValida4D - ERR50'
            if (iflag     == 0       ) stop 'ReadProperties - ModuleValida4D - ERR60'
        

            call GetData(Me%Properties(iP)%ThreeD,                                      &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = '3D',                                           &
                         default      = .true.,                                         &
                         ClientModule = 'ModuleValida4D',                               &
                         STAT         = STAT_CALL)                                          
            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleValida4D - ERR70'

        

        enddo do2



        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleValida4D - ERR80'

        
        
        

    end subroutine ReadProperties

    !--------------------------------------------------------------------------

    subroutine ReadHDF5FileName

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, ClientNumber, line , iflag
        integer                                     :: LastLine, FirstLine, iH, HDF5_READ
        logical                                     :: BlockFound, exist

        !Begin-----------------------------------------------------------------


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5FileName - ModuleValida4D - ERR10'


        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginHDF5>',                    &
                                    block_end       = '<EndHDF5>',                      &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)


IS:     if(STAT_CALL .EQ. SUCCESS_) then

BF:         if (BlockFound) then

                Me%HDF5Number = LastLine - FirstLine - 1
                
                allocate(Me%HDF5Files(Me%HDF5Number))
                
                do iH = 1, Me%HDF5Number
                
                    Me%HDF5Files(iH)%ObjHDF5          = 0
                    
                    Me%HDF5Files(iH)%NumberOfInstants = FillValueInt
                    
                    Me%HDF5Files(iH)%Name             = ' ' 
                    
                    call null_time(Me%HDF5Files(iH)%StartTime )
                    call null_time(Me%HDF5Files(iH)%EndTime  )

                enddo                
                
                iH = 0

                do line = FirstLine + 1, LastLine - 1
                
                    iH = iH + 1

                    call GetData(Me%HDF5Files(iH)%Name, EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5FileName - ModuleValida4D - ERR20'

                    inquire(file = Me%HDF5Files(iH)%Name, exist = exist)
       
i1:                 if (exist) then
                                                    
                        call ConstructHDF5 (Me%HDF5Files(iH)%ObjHDF5, Me%HDF5Files(iH)%Name, HDF5_READ, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5FileName - ModuleValida4D - ERR30'
                        
                        call GetHDF5GroupNumberOfItems(Me%HDF5Files(iH)%ObjHDF5, "/Time", &
                                                       Me%HDF5Files(iH)%NumberOfInstants, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5FileName - ModuleValida4D - ERR40'
                        
                        Me%HDF5Files(iH)%StartTime = HDF5TimeInstant(iH, 1)
                        Me%HDF5Files(iH)%EndTime   = HDF5TimeInstant(iH, Me%HDF5Files(iH)%NumberOfInstants)                        
                        
                    else i1

                        stop 'ReadHDF5FileName - ModuleValida4D - ERR50'

                    endif i1
                enddo

            else BF

                stop 'ReadHDF5FileName - ModuleValida4D - ERR60'

            end if BF

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ReadHDF5FileName - ModuleValida4D - ERR70'

        else   IS

            stop 'ReadHDF5FileName - ModuleValida4D - ERR80'

        end if IS


    end subroutine ReadHDF5FileName

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine ReadInputTable

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        integer                                     :: iLength, iAux
        integer                                     :: FirstLine, LastLine, line
        integer                                     :: iP, iV, i, ClientNumber
        logical                                     :: exist, BlockFound
        character (len = PathLength)                :: AuxChar, AuxCharV
        real, dimension(100)                        :: AuxVector
        
        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData (Me%ObjEnterData, Me%InputTable, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadInputTable - ModuleValida4D - ERR10'


        !Gets start time
        call GetData(Me%InitialDate,                                                    &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      ='SERIE_INITIAL_DATA',                                &
                     ClientModule ='ModuleTimeSerie',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputTable - ModuleValida4D - ERR20'
        if (iflag     /= 1       ) stop 'ReadInputTable - ModuleValida4D - ERR30'


        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInputTable - ModuleValida4D - ERR40'


        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginTable>',                   &
                                    block_end       = '<EndTable>',                     &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS:     if (STAT_CALL == SUCCESS_) then

BF:         if (BlockFound) then
 
                Me%TableValues = LastLine - FirstLine - 1
                
                allocate(Me%T(Me%TableValues))
                allocate(Me%X(Me%TableValues),Me%Y(Me%TableValues))
                allocate(Me%i(Me%TableValues),Me%j(Me%TableValues))
                allocate(Me%PercI(Me%TableValues),Me%PercJ(Me%TableValues))      
                allocate(Me%NullValue(Me%TableValues))          
                
                Me%NullValue(:) = .false.
                
                allocate(Me%StationName(Me%TableValues))
                
                if (Me%Zcolumn > FillValueInt) allocate(Me%Z(Me%TableValues))
                
                do iP=1, Me%PropNumber
                    allocate(Me%Properties(iP)%ValueTable(Me%TableValues))
                    allocate(Me%Properties(iP)%ValueHDF5 (Me%TableValues))
                    
                    Me%Properties(iP)%ValueTable(:) = FillValueReal
                    Me%Properties(iP)%ValueHDF5 (:) = FillValueReal
                    
                enddo
                
               
                line = FirstLine + 1

                
                call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                                 
                if (STAT_CALL /= SUCCESS_) stop 'ReadInputTable - ModuleValida4D - ERR50'
                                 
                iAux = scan(AuxChar,',')
                
                if (iAux == 0) then 

                    iLength = len_trim(AuxChar)
                
                else
                
                    iLength = len_trim(AuxChar(1:iAux - 1))

                endif
                
                AuxCharV = adjustl(AuxChar(1:iLength))                

                
cd1 :           if (iLength > 0) then ! counts the number of data
                    Me%TableColumns = 1
do1 :                   do i = 2, iLength
                            if((AuxCharV (i  :i  ) == " ") .AND.                        &
                               (AuxCharV (i-1:i-1) /= " "))                             &
                                Me%TableColumns = Me%TableColumns + 1 
                        end do do1
                end if cd1                
                
                read(AuxCharV(1:iLength),*) (AuxVector(i),i=1,Me%TableColumns)
                
                if (Me%Tcolumn > Me%TableColumns ) stop 'ReadInputTable - ModuleValida4D - ERR60'
                
                if (Me%Zcolumn > Me%TableColumns ) stop 'ReadInputTable - ModuleValida4D - ERR70'
                
                if (Me%Xcolumn > Me%TableColumns ) stop 'ReadInputTable - ModuleValida4D - ERR80'
                
                if (Me%Ycolumn > Me%TableColumns ) stop 'ReadInputTable - ModuleValida4D - ERR90'
                
                do iP=1, Me%PropNumber
                    if (Me%Properties(iP)%Column >  Me%TableColumns) stop 'ReadInputTable - ModuleValida4D - ERR100'
                enddo
                
                iV = 0 
               
                do line = FirstLine + 1, LastLine - 1
                
                    iV = iV + 1

                    call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop 'ReadInputTable - ModuleValida4D - ERR110'
                    
                    iAux = scan(AuxChar,',')
                    
                    if (iAux == 0) then 

                        iLength = len_trim(AuxChar) 
                
                    else
                
                        iLength = iAux - 1
                
                    endif
    
                    
                    
                    read(AuxChar(1:iLength),*) AuxVector(1:Me%TableColumns)
                    
                    Me%T(iV) = AuxVector(Me%Tcolumn)
                    
                    if (iV > 1) then
                    
                        if (Me%T(iV) < Me%T(iV-1)) then
                            write(*,*) 'In MOHID Table format time must be order in a ascending way'
                            write(*,*) 'see table lines - ', iV-1, iV
                            stop 'ReadInputTable - ModuleValida4D - ERR120'
                        endif
                    
                    endif 

                    if (iAux == 0) then
                        Me%StationName(iV) = "NoName"
                    else
                        Me%StationName(iV) = trim(AuxChar(iAux+1:PathLength))
                    endif
                    
                    Me%X(iV) = AuxVector(Me%Xcolumn)
                    Me%Y(iV) = AuxVector(Me%Ycolumn)
                    
                    if (Me%Zcolumn > FillValueInt) then
                        Me%Z(iV) = AuxVector(Me%Zcolumn)
                    endif
                    
                    do iP=1, Me%PropNumber
                        Me%Properties(iP)%ValueTable(iV) = AuxVector(Me%Properties(iP)%Column)
                    enddo                    

                enddo

            else BF

                stop 'ReadInputTable - ModuleValida4D - ERR130'

            end if BF

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ReadInputTable - ModuleValida4D - ERR140'

        else   IS

            stop 'ReadInputTable - ModuleValida4D - ERR150'

        end if IS
        
        Me%Zmax = FillValueReal
        
        do iV =1, Me%TableValues
            if (Me%Z(iV)>Me%Zmax) Me%Zmax = Me%Z(iV)
        enddo
        
        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadInputTable - ModuleValida4D - ERR160'


              
    end subroutine ReadInputTable

    !--------------------------------------------------------------------------
    
    subroutine ConstructSubModules

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_Time)                               :: BeginTime, EndTime
        
        !Begin-----------------------------------------------------------------
        
        BeginTime = Me%InitialDate + Me%T(1)
        EndTime   = Me%InitialDate + Me%T(Me%TableValues)

        call StartComputeTime       (TimeID           = Me%ObjTime,                     &
                                     InitialSystemTime= BeginTime,                      &
                                     BeginTime        = BeginTime,                      &
                                     EndTime          = EndTime,                        &
                                     DT               = 1.,                             & 
                                     VariableDT       = .false.,                        &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSubModules - ModuleValida4D - ERR10'


        call ConstructHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     DataFile         = Me%BathymetryFile,              &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSubModules - ModuleValida4D - ERR20'


        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     FileName         = Me%BathymetryFile,              &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSubModules - ModuleValida4D - ERR30'

        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSubModules - ModuleValida4D - ERR40'

        !Geometry - Water Column
        call ConstructGeometry      (GeometryID       = Me%ObjGeometry,                 &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     NewDomain        = Me%GeometryFile,                &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSubModules - ModuleValida4D - ERR50'


        call ConstructMap           (Map_ID           = Me%ObjMap,                      &
                                     GeometryID       = Me%ObjGeometry,                 &
                                     HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     TimeID           = Me%ObjTime,                     &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSubModules - ModuleValida4D - ERR60'
        
    end subroutine ConstructSubModules        


    !--------------------------------------------------------------------------
    
    subroutine AllocateMatrixes

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                     :: iP, i, j, k, iH, STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                     :: HaveTimeLimits

        
        !Begin-----------------------------------------------------------------    
 
         call ConstructEnterData (Me%ObjEnterData, Me%InputTable, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR10'

        
        !Gets the size from the Geometry
        call GetGeometrySize(Me%ObjGeometry,                                            &
                             Size        = Me%Size,                                     &
                             WorkSize    = Me%WorkSize,                                 &
                             STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR20'
        
        Me%WorkSize2D%ILB = Me%WorkSize%ILB
        Me%WorkSize2D%IUB = Me%WorkSize%IUB
        Me%WorkSize2D%JLB = Me%WorkSize%JLB
        Me%WorkSize2D%JUB = Me%WorkSize%JUB


        !Auxiliar variables for the do loops
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB      
        
        allocate(Me%LayerValues (ILB:IUB,JLB:JUB))
        allocate(Me%FullMatrix2D(ILB:IUB,JLB:JUB))
        allocate(Me%FullMatrix3D(ILB:IUB,JLB:JUB,KLB:KUB))        
        
        Me%FullMatrix2D(:,:  ) = 1
        Me%FullMatrix3D(:,:,:) = 1        
        
        call GetWaterPoints3D(Me%ObjMap, Me%External_Var%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR30'

        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%External_Var%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR40'
        
         call GetHorizontalGrid(Me%ObjHorizontalGrid,                                   &
                                DZX  = Me%External_Var%DZX,                             &
                                DZY  = Me%External_Var%DZY,                             &
                                DUX  = Me%External_Var%DUX,                             &
                                DVY  = Me%External_Var%DVY,                             &
                                STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR50'                                
       
        do iP=1, Me%PropNumber
        
            if (Me%Properties(iP)%ThreeD) then
            
                allocate(Me%Properties(iP)%InterpolTime%Values3D    (ILB:IUB,JLB:JUB,KLB:KUB))
                allocate(Me%Properties(iP)%InterpolTime%PrevValues3D(ILB:IUB,JLB:JUB,KLB:KUB))
                allocate(Me%Properties(iP)%InterpolTime%NextValues3D(ILB:IUB,JLB:JUB,KLB:KUB))
                
                
                Me%Properties(iP)%InterpolTime%Values3D    (:,:,:) = FillValueReal
                Me%Properties(iP)%InterpolTime%PrevValues3D(:,:,:) = FillValueReal
                Me%Properties(iP)%InterpolTime%NextValues3D(:,:,:) = FillValueReal
                
            else
            
                allocate(Me%Properties(iP)%InterpolTime%Values2D    (ILB:IUB,JLB:JUB))
                allocate(Me%Properties(iP)%InterpolTime%PrevValues2D(ILB:IUB,JLB:JUB))
                allocate(Me%Properties(iP)%InterpolTime%NextValues2D(ILB:IUB,JLB:JUB))
                
                Me%Properties(iP)%InterpolTime%Values2D    (:,:) = FillValueReal
                Me%Properties(iP)%InterpolTime%PrevValues2D(:,:) = FillValueReal
                Me%Properties(iP)%InterpolTime%NextValues2D(:,:) = FillValueReal
                
            endif
            
        enddo
        
        allocate(Me%VerticalZ%Mapping                  (ILB:IUB,JLB:JUB,KLB:KUB))
        
        allocate(Me%VerticalZ%InterpolTime%Values3D    (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%VerticalZ%InterpolTime%PrevValues3D(ILB:IUB,JLB:JUB,KLB:KUB))
        allocate(Me%VerticalZ%InterpolTime%NextValues3D(ILB:IUB,JLB:JUB,KLB:KUB))
       
        Me%VerticalZ%InterpolTime%Values3D    (:,:,:) = FillValueReal
        Me%VerticalZ%InterpolTime%PrevValues3D(:,:,:) = FillValueReal
        Me%VerticalZ%InterpolTime%NextValues3D(:,:,:) = FillValueReal
        
        Me%VerticalZ%Mapping                  (:,:,:) = Me%External_Var%WaterPoints3D(:,:,:)
        
        Me%VerticalZ%ID%Name        = GetPropertyName(VerticalZ_)
        Me%VerticalZ%ID%Units       = 'm'
        Me%VerticalZ%ID%Description = 'vertical geometry of mohid'
        
        call GetGeometryKFloor(Me%ObjGeometry, Z = Me%External_Var%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR70'   
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%External_Var%WaterPoints3D(i, j, Me%WorkSize%KUB) == WaterPoint) then
                k = Me%External_Var%KFloor_Z(i, j) - 1
                Me%VerticalZ%Mapping(i, j, k) = WaterPoint
            endif
            
        enddo
        enddo

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AllocateMatrixes - ModuleValida4D - ERR110'

    end subroutine AllocateMatrixes
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyValida4D()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

        call ConvertHDF5_to_Table

    end subroutine ModifyValida4D

    !--------------------------------------------------------------------------

    subroutine ConvertHDF5_to_Table()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------

       
        !Interpolate the input table values from hdf5 files
        call GenerateNewTable

        !write new table
        call WriteNewTable
        
        !Compare the input table values with the hdf5 values
        !call ComputeStatistics
        
        
    end subroutine ConvertHDF5_to_Table


    !--------------------------------------------------------------------------

    subroutine GenerateNewTable()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension (:,:,:), pointer :: Values3D
        real, dimension (:,:  ), pointer :: Values2D
        real(8), dimension (:), pointer  :: DepthsSW, DepthsSE, DepthsNW, DepthsNE, ValuesAux
        real                             :: X_W, X_E, Xv, Y_S, Y_N, Yv        
        real                             :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        integer                          :: iV, iH, iP, STAT_CALL, Ndepths
        integer                          :: i, j, k, kb, iI, kbb
        integer                          :: jW, jE, iS, iN, imin, jmin, imax, jmax
        integer                          :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                          :: kmin_prev, kmin_next
        logical                          :: NewFields, NewPosition, InsideDomain, FoundBottom, FoundSurface
        type (T_Time)                    :: CurrentTime
        type (T_Size3D)                  :: SizeZ    
        !----------------------------------------------------------------------
        
        allocate(DepthsSW(Me%Size%KLB:Me%Size%KUB), DepthsNE(Me%Size%KLB:Me%Size%KUB),  &
                 DepthsSE(Me%Size%KLB:Me%Size%KUB), DepthsNW(Me%Size%KLB:Me%Size%KUB),  &
                 ValuesAux(Me%Size%KLB:Me%Size%KUB))
        
        DepthsSW(:)  = FillValueReal
        DepthsSE(:)  = FillValueReal
        DepthsNW(:)  = FillValueReal
        DepthsNE(:)  = FillValueReal
        
        ValuesAux(:)  = FillValueReal
        
        call null_time(Me%PrevTime)
        call null_time(Me%NextTime)
        
        SizeZ     = Me%WorkSize
        SizeZ%KLB = Me%WorkSize%KLB - 1

        ILB       = Me%WorkSize%ILB
        IUB       = Me%WorkSize%IUB

        JLB       = Me%WorkSize%JLB
        JUB       = Me%WorkSize%JUB

        KLB       = Me%WorkSize%KLB
        KUB       = Me%WorkSize%KUB        
        
        NewFields = .true. 
        
        imin = - FillValueInt
        imax =   FillValueInt        

        jmin = - FillValueInt
        jmax =   FillValueInt        

diV:    do iV = 1, Me%TableValues       

            InsideDomain = GetXYInsideDomain(Me%ObjHorizontalGrid, Me%X(iV), Me%Y(iV), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GenerateNewTable - ModuleValida4D - ERR10'
            
            if (.not. InsideDomain) then
                write(*,*) 'Point ',iV, ' not inside the model domain'
                !stop 'GenerateNewTable - ModuleValida4D - ERR20'
                Me%NullValue(iV) = .true.
                cycle
            endif
            
            if (iV == 1) NewPosition = .true.
            
            if (iV > 1) then
                if (Me%X(iV) /= Me%X(iV-1) .or. Me%Y(iV) /= Me%Y(iV-1)) then
                    NewPosition = .true.
                else
                    NewPosition = .false.
                endif
            endif
            
            if (NewPosition) then
            
                call GetXYCellZ(Me%ObjHorizontalGrid, Me%X(iV), Me%Y(iV), Me%i(iV), Me%j(iV), Me%PercI(iV), Me%PercJ(iV), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'GenerateNewTable - ModuleValida4D - ERR30'
                
                imin =min(Me%i(iV), imin)
                jmin =min(Me%j(iV), jmin)

                imax =max(Me%i(iV), imax)
                jmax =max(Me%j(iV), jmax)
            
            else 
                Me%X    (iV) = Me%X    (iV-1)
                Me%Y    (iV) = Me%Y    (iV-1)
                Me%i    (iV) = Me%i    (iV-1)
                Me%j    (iV) = Me%j    (iV-1)
                Me%PercI(iV) = Me%PercI(iV-1)
                Me%PercJ(iV) = Me%PercJ(iV-1)
            endif
            
        enddo diV
        
        imin =max(ILB, imin-1)
        jmin =max(JLB, jmin-1)

        imax =min(IUB, imax+1)
        jmax =min(JUB, jmax+1)        
        
diV1:   do iV = 1, Me%TableValues  
      
            if (Me%PercJ(iV) > 0.5) then
                jW = Me%j(iV)
                jE = Me%j(iV)+1
                Xv = Me%PercJ(iV) - 0.5
            else
                jW = Me%j(iV)-1
                jE = Me%j(iV)
                Xv = Me%PercJ(iV) + 0.5
            endif
            
            if (Me%PercI(iV) > 0.5) then
                iS = Me%i(iV)
                iN = Me%i(iV)+1
                Yv = Me%PercI(iV) - 0.5
            else
                iS = Me%i(iV)-1
                iN = Me%i(iV)
                Yv = Me%PercI(iV) + 0.5
            endif            
            
            X_W = 0.
            X_E = 1
            Y_S = 0.                
            Y_N = 1.


            !if (Me%External_Var%WaterPoints3D(i, j, Me%WorkSize%KUB) /= WaterPoint) cycle 

            CurrentTime = Me%InitialDate + Me%T(iV)
            
            call ActualizeCurrentTime(TimeID    = Me%ObjTime,                           &
                                      DT_Global = 0.,                                   &
                                      Current   = CurrentTime,                          &
                                      STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GenerateNewTable - ModuleValida4D - ERR40'

            if (CurrentTime > Me%NextTime) then 
                NewFields = .true. 
            else
                NewFields = .false. 
            endif
            
            if (NewFields) then
                    
                !Search for the correct hdf file
                do iH = 1, Me%HDF5Number       
                    if     (CurrentTime .GE. Me%HDF5Files(iH)%StartTime) then
                        if (CurrentTime .LE. Me%HDF5Files(iH)%EndTime  ) exit
                    endif
                enddo
                
                if (iH > Me%HDF5Number) then
                    write(*,*) 'data out the hdf time period'
                    Me%NullValue(iV) = .true.
                    cycle
                endif
                
            endif            
            
            if (NewFields) then
                    
                !Search for the correct hdf file
                do iH = 1, Me%HDF5Number       
                    if     (CurrentTime .GE. Me%HDF5Files(iH)%StartTime) then
                        if (CurrentTime .LE. Me%HDF5Files(iH)%EndTime  ) exit
                    endif
                enddo
                
                do iI = 1, Me%HDF5Files(iH)%NumberOfInstants
                
                    Me%PrevTime = HDF5TimeInstant(iH, iI)
                    Me%NextTime = HDF5TimeInstant(iH, iI+1)                        
                    if     (CurrentTime .GE. Me%PrevTime) then
                        if (CurrentTime .LE. Me%NextTime) exit
                    endif
                enddo
                
            endif
            

            !interpolate VerticalZ
            if (Me%WorkSize%KUB > 1) then
            
                if (NewFields) then
                    
                    call HDF5SetLimits  (Me%HDF5Files(iH)%ObjHDF5, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleValida4D - ERR50'
                    
                         
                    call HDF5ReadData(HDF5ID      = Me%HDF5Files(iH)%ObjHDF5,               &
                                      GroupName   = "/Grid/VerticalZ",                      &
                                      Name        = "Vertical",                             &
                                      Array3D     = Me%VerticalZ%InterpolTime%PrevValues3D, &
                                      OutputNumber= iI,                                     &
                                      STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_ )stop 'ReadHDF5Values3D - ModuleValida4D - ERR60' 
                    
                    kmin_prev = - FillValueInt
                    
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if (Me%External_Var%WaterPoints3D(i, j, KUB) == WaterPoint) then

                            do k = Me%External_Var%KFloor_Z(i, j), KUB
                        
                                if (Me%VerticalZ%InterpolTime%PrevValues3D(i,j,k  ) >= Me%Zmax .and. &
                                    Me%VerticalZ%InterpolTime%PrevValues3D(i,j,k+1) <= Me%Zmax) then
                                    if (k < kmin_prev) kmin_prev = k
                                endif 
                                
                            enddo
                            
                        endif
 
                    enddo  
                    enddo                  

                    call FillMatrix3DNearestCell (imin, imax, jmin, jmax, kmin_prev-1, KUB,     &
                                                  Me%VerticalZ%Mapping,                     &
                                                  Me%VerticalZ%InterpolTime%PrevValues3D)
                   

                    call HDF5ReadData(HDF5ID      = Me%HDF5Files(iH)%ObjHDF5,               &
                                      GroupName   = "/Grid/VerticalZ",                      &
                                      Name        = "Vertical",                             &
                                      Array3D     = Me%VerticalZ%InterpolTime%NextValues3D, &
                                      OutputNumber= iI + 1,                                 &
                                      STAT         = STAT_CALL)
                                      
                    if (STAT_CALL /= SUCCESS_ )stop 'ReadHDF5Values3D - ModuleValida4D - ERR70'
                    
                    kmin_next = - FillValueInt
                    
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if (Me%External_Var%WaterPoints3D(i, j, KUB) == WaterPoint) then

                            do k = Me%External_Var%KFloor_Z(i, j), KUB
                                                
                                if (Me%VerticalZ%InterpolTime%NextValues3D(i,j,k  ) >= Me%Zmax .and. &
                                    Me%VerticalZ%InterpolTime%NextValues3D(i,j,k+1) <= Me%Zmax) then
                                    if (k < kmin_next) kmin_next = k
                                endif 
                            
                            enddo
                            
                        endif
 
                    enddo
                    enddo                    
                    
                    call FillMatrix3DNearestCell (imin, imax, jmin, jmax, kmin_next-1, KUB, &
                                                  Me%VerticalZ%Mapping,                 &
                                                  Me%VerticalZ%InterpolTime%NextValues3D)                        
                
                
                endif
                
                nullify(Values3D)
                
                if      (CurrentTime .EQ. Me%PrevTime) then
                
                    Values3D => Me%VerticalZ%InterpolTime%PrevValues3D
                    
                else if (CurrentTime .EQ. Me%NextTime) then
                
                    Values3D => Me%VerticalZ%InterpolTime%NextValues3D
                    
                else
                
                    call InterpolateMatrix3DInTime(ActualTime = CurrentTime,            &
                                                   Size       = SizeZ,                  &
                                                   Time1      = Me%PrevTime,            &
                                                   Matrix1    = Me%VerticalZ%InterpolTime%PrevValues3D, &
                                                   Time2      = Me%NextTime,            &
                                                   Matrix2    = Me%VerticalZ%InterpolTime%NextValues3D, &
                                                   MatrixOUT  = Me%VerticalZ%InterpolTime%Values3D, &
                                                   PointsToFill3D = Me%FullMatrix3D)
                    
                    Values3D => Me%VerticalZ%InterpolTime%Values3D
                    
                endif
                
                kb = max(kmin_prev, kmin_next)
                
                Ndepths = Me%WorkSize%KUB - kb + 1

                do k = kb, Me%WorkSize%KUB

                    DepthsSW(k) = (Values3D(iS, jW, k) + Values3D(iS, jW, k-1)) / 2. - Values3D(iS, jW, Me%WorkSize%KUB)
                    DepthsSE(k) = (Values3D(iS, jE, k) + Values3D(iS, jE, k-1)) / 2. - Values3D(iS, jE, Me%WorkSize%KUB)
                    DepthsNW(k) = (Values3D(iN, jW, k) + Values3D(iN, jW, k-1)) / 2. - Values3D(iN, jW, Me%WorkSize%KUB)
                    DepthsNE(k) = (Values3D(iN, jE, k) + Values3D(iN, jE, k-1)) / 2. - Values3D(iN, jE, Me%WorkSize%KUB)                                        

                enddo                
            
            endif

diP:        do iP = 1, Me%PropNumber

i3D:            if (Me%Properties(iP)%ThreeD) then
                                 
i12:                if (NewFields) then
                        
                        call HDF5SetLimits  (Me%HDF5Files(iH)%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleValida4D - ERR80'
                        
                             
                        call HDF5ReadData(HDF5ID      = Me%HDF5Files(iH)%ObjHDF5,               &
                                          GroupName   = "/Results/"//trim(Me%Properties(iP)%ID%Name),&
                                          Name        = trim(Me%Properties(iP)%ID%Name),             &
                                          Array3D     = Me%Properties(iP)%InterpolTime%PrevValues3D, &
                                          OutputNumber= iI,                                     &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_ )stop 'ReadHDF5Values3D - ModuleValida4D - ERR90'          
                        
                        call FillMatrix3DNearestCell (imin, imax, jmin, jmax, kmin_prev, KUB,&
                                                      Me%External_Var%Waterpoints3D,     &
                                                      Me%Properties(iP)%InterpolTime%PrevValues3D)
                                  

                        call HDF5ReadData(HDF5ID      = Me%HDF5Files(iH)%ObjHDF5,               &
                                          GroupName   = "/Results/"//trim(Me%Properties(iP)%ID%Name),&
                                          Name        = trim(Me%Properties(iP)%ID%Name),             &
                                          Array3D     = Me%Properties(iP)%InterpolTime%NextValues3D, &
                                          OutputNumber= iI + 1,                                 &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_ )stop 'ReadHDF5Values3D - ModuleValida4D - ERR100'                    
                        
                        call FillMatrix3DNearestCell (imin, imax, jmin, jmax, kmin_next, KUB,&
                                                      Me%External_Var%Waterpoints3D,     &
                                                      Me%Properties(iP)%InterpolTime%NextValues3D)                        
                    
                    
                    endif i12
                    
                    nullify(Values3D)
                    
                    if      (CurrentTime .EQ. Me%PrevTime) then
                    
                        Values3D => Me%Properties(iP)%InterpolTime%PrevValues3D
                        
                    else if (CurrentTime .EQ. Me%NextTime) then
                    
                        Values3D => Me%Properties(iP)%InterpolTime%NextValues3D
                        
                    else
                    
                        call InterpolateMatrix3DInTime(ActualTime     = CurrentTime,                                 &
                                                       Size           = Me%WorkSize,                                 &
                                                       Time1          = Me%PrevTime,                                 &
                                                       Matrix1        = Me%Properties(iP)%InterpolTime%PrevValues3D, &
                                                       Time2          = Me%NextTime,                                 &
                                                       Matrix2        = Me%Properties(iP)%InterpolTime%NextValues3D, &
                                                       MatrixOUT      = Me%Properties(iP)%InterpolTime%Values3D,     &
                                                       PointsToFill3D = Me%FullMatrix3D)
                         
                        Values3D => Me%Properties(iP)%InterpolTime%Values3D
                        
                    endif                
                    
!                    do k = kb, Me%WorkSize%KUB
                    
!                         Me%LayerValues   (:,:) = Values3D(:,:, k)
               
!                         Values(k) = InterpolXYPoint(Me%ObjHorizontalGrid, Me%LayerValues,  &
!                                                     Me%FullMatrix2D, Me%X(iV), Me%Y(iV),   &
!                                                     STAT = STAT_CALL)                    
!                        if (STAT_CALL /= SUCCESS_) stop 'GenerateNewTable - ModuleValida4D - ERR110'
                        
                        
!                    enddo

                    if (Me%WorkSize%KUB > 1) then
                    
                        ValuesAux (kb:Me%WorkSize%KUB) = Values3D(iS, jW, kb:Me%WorkSize%KUB)
                        
                        do kbb = kb, Me%WorkSize%KUB
                            if (ValuesAux(kbb)> FillValueReal/1000.) exit
                        enddo 
                        
                        Ndepths = Me%WorkSize%KUB - kbb + 1

                        ValueSW = InterpolateProfileR8(dble(Me%Z(iV)), Ndepths,             &
                                                       DepthsSW  (kbb:Me%WorkSize%KUB),     &
                                                       ValuesAux (kbb:Me%WorkSize%KUB),     &
                                                       FoundBottom, FoundSurface)

                        ValuesAux (kb:Me%WorkSize%KUB) = Values3D(iS, jE, kb:Me%WorkSize%KUB)
                        
                        do kbb = kb, Me%WorkSize%KUB
                            if (ValuesAux(kbb)> FillValueReal/1000.) exit
                        enddo 
                        
                        Ndepths = Me%WorkSize%KUB - kbb + 1
                        

                        ValueSE = InterpolateProfileR8(dble(Me%Z(iV)), Ndepths,             &
                                                       DepthsSE  (kbb:Me%WorkSize%KUB),     &
                                                       ValuesAux (kbb:Me%WorkSize%KUB),     &
                                                       FoundBottom, FoundSurface)
                        ValuesAux (kb:Me%WorkSize%KUB) = Values3D(iN, jW, kb:Me%WorkSize%KUB)

                        do kbb = kb, Me%WorkSize%KUB
                            if (ValuesAux(kbb)> FillValueReal/1000.) exit
                        enddo 
                        
                        Ndepths = Me%WorkSize%KUB - kbb + 1

                        ValueNW = InterpolateProfileR8(dble(Me%Z(iV)), Ndepths,             &
                                                       DepthsNW(kbb:Me%WorkSize%KUB),       &
                                                       ValuesAux (kbb:Me%WorkSize%KUB),     &
                                                       FoundBottom, FoundSurface)

                        ValuesAux (kb:Me%WorkSize%KUB) = Values3D(iN, jE, kb:Me%WorkSize%KUB)
                        
                        do kbb = kb, Me%WorkSize%KUB
                            if (ValuesAux(kbb)> FillValueReal/1000.) exit
                        enddo 
                        
                        Ndepths = Me%WorkSize%KUB - kbb + 1
                        

                        ValueNE = InterpolateProfileR8(dble(Me%Z(iV)), Ndepths,             &
                                                       DepthsNE  (kbb:Me%WorkSize%KUB),     &
                                                       ValuesAux (kbb:Me%WorkSize%KUB),     &
                                                       FoundBottom, FoundSurface)
                                                       
                        ValueN = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                        ValueS = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                        
                        Me%Properties(iP)%ValueHDF5(iV) = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)
                        
 
                    
                    else
                    
                        Me%LayerValues   (:,:) = Values3D(:,:, Me%WorkSize%KUB)
                    
                        Me%Properties(iP)%ValueHDF5(iV) = InterpolXYPoint(Me%ObjHorizontalGrid, Me%LayerValues,  &
                                                          Me%FullMatrix2D, Me%X(iV), Me%Y(iV),   &
                                                          STAT = STAT_CALL)                    
                        if (STAT_CALL /= SUCCESS_) stop 'GenerateNewTable - ModuleValida4D - ERR120'

                    
                    endif
                
                else i3D
                
                    
iN2:                if (NewFields) then
                        
                        call HDF5SetLimits  (Me%HDF5Files(iH)%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5Values3D - ModuleValida4D - ERR130'
                        
                             
                        call HDF5ReadData(HDF5ID      = Me%HDF5Files(iH)%ObjHDF5,               &
                                          GroupName   = "/Results/"//trim(Me%Properties(iP)%ID%Name),&
                                          Name        = trim(Me%Properties(iP)%ID%Name),             &
                                          Array2D     = Me%Properties(iP)%InterpolTime%PrevValues2D, &
                                          OutputNumber= iI,                                     &
                                          STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_ )stop 'ReadHDF5Values3D - ModuleValida4D - ERR140'
                        
                        
                        call FillMatrix2DNearestCell (imin, imax, jmin, jmax,               &
                                                      Me%External_Var%Waterpoints2D,    &
                                                      Me%Properties(iP)%InterpolTime%PrevValues2D)

                        call HDF5ReadData(HDF5ID      = Me%HDF5Files(iH)%ObjHDF5,               &
                                          GroupName   = "/Results/"//trim(Me%Properties(iP)%ID%Name),&
                                          Name        = trim(Me%Properties(iP)%ID%Name),             &
                                          Array2D     = Me%Properties(iP)%InterpolTime%NextValues2D, &
                                          OutputNumber= iI + 1,                                 &
                                          STAT         = STAT_CALL)
                                          
                        if (STAT_CALL /= SUCCESS_ )stop 'ReadHDF5Values3D - ModuleValida4D - ERR150'                    

                        call FillMatrix2DNearestCell(imin, imax, jmin, jmax,                &
                                                     Me%External_Var%Waterpoints2D,     &
                                                     Me%Properties(iP)%InterpolTime%NextValues2D)                    
                    
                    endif iN2
                    
                    nullify(Values2D)
                    
                    if      (CurrentTime .EQ. Me%PrevTime) then
                    
                        Values2D => Me%Properties(iP)%InterpolTime%PrevValues2D
                        
                    else if (CurrentTime .EQ. Me%NextTime) then
                    
                        Values2D => Me%Properties(iP)%InterpolTime%NextValues2D
                        
                    else
                    
                        call InterpolateMatrix2DInTime(ActualTime     = CurrentTime,                                 &
                                                       Size           = Me%WorkSize2D,                               &
                                                       Time1          = Me%PrevTime,                                 &
                                                       Matrix1        = Me%Properties(iP)%InterpolTime%PrevValues2D, &
                                                       Time2          = Me%NextTime,                                 &
                                                       Matrix2        = Me%Properties(iP)%InterpolTime%NextValues2D, &
                                                       MatrixOUT      = Me%Properties(iP)%InterpolTime%Values2D,     &
                                                       PointsToFill2D = Me%FullMatrix2D)
                         
                        Values2D => Me%Properties(iP)%InterpolTime%Values2D
                        
                    endif                
                    
!                    Me%Properties(iP)%ValueHDF5(iV) = InterpolXYPoint(Me%ObjHorizontalGrid, Values2D,  &
!                                                     Me%FullMatrix2D, Me%X(iV), Me%Y(iV),              &
!                                                     STAT = STAT_CALL)                    
!                    if (STAT_CALL /= SUCCESS_) stop 'GenerateNewTable - ModuleValida4D - ERR160'

                    ValueSW = Values2D(iS, jW)

                    ValueSE = Values2D(iS, jE)

                    ValueNW = Values2D(iN, jW)

                    ValueNE = Values2D(iN, jE)
                                                   
                    ValueN = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                    ValueS = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                    
                    Me%Properties(iP)%ValueHDF5(iV) = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)
                        
 
                
                endif i3D
            enddo diP
        
        enddo diV1
        
    end subroutine GenerateNewTable
    
    !--------------------------------------------------------------------------

    subroutine     WriteNewTable()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                      :: iV, iP, STAT_CALL, unit
        real                         :: year, month, day, hour, minutes, seconds
        character(len=1000)          :: AuxC
        !----------------------------------------------------------------------
        
        call UnitsManager(unit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewTable - ModuleValida4D - ERR10'


        open(UNIT = unit, FILE = trim(adjustl(Me%OutputTable)),                         &
             FORM = "FORMATTED",   STATUS = "UNKNOWN", ACTION = "WRITE", IOSTAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewTable - ModuleValida4D - ERR20'

        !Gets Date in format YY MM DD hh mm ss
        call ExtractDate(Me%InitialDate, year, month, day, hour, minutes, seconds)        
        write(AuxC,*) "SERIE_INITIAL_DATA : ", year, month, day, hour, minutes, seconds 
        write(unit,'(A)') trim(adjustl(AuxC))
        
        AuxC = ' '
        
        write(unit,*) "Seconds   X     Y     Z   ", (trim(Me%Properties(iP)%ID%Name)," ", iP=1,Me%PropNumber), " StationName"
        
        write(unit,*) "<BeginTable>"
        
diV:    do iV = 1, Me%TableValues       
            if (Me%NullValue(iV)) then
                write(AuxC,*) Me%T(iV), Me%X(iV), Me%Y(iV), Me%Z(iV), (" null ", iP=1,Me%PropNumber), " , ", Me%StationName(iV)
            else
                write(AuxC,*) Me%T(iV), Me%X(iV), Me%Y(iV), Me%Z(iV), (Me%Properties(iP)%ValueHDF5(iV), iP=1,Me%PropNumber), " , ", Me%StationName(iV)
            endif
            
            write(unit,'(A)') trim(adjustl(AuxC))
            
        enddo diV
        
        write(unit,*) "<EndTable>"        
        
        call UnitsManager(unit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteNewTable - ModuleValida4D - ERR30'
        
        
    end subroutine WriteNewTable


    type(T_Time) function HDF5TimeInstant(iH, Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: iH, Instant
        

        !Local-----------------------------------------------------------------
!        type(T_Time)                            :: TimeInstant
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%HDF5Files(iH)%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%HDF5Files(iH)%ObjHDF5,                 &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleValida4D - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

                                     
        deallocate(TimeVector)

    end function HDF5TimeInstant
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillValida4D()

        !Arguments---------------------------------------------------------------



        !Local-------------------------------------------------------------------
        integer                             :: STAT_CALL           
        

        !Begin-------------------------------------------------------------------

        !Kills Map
        call KillMap            (Me%ObjMap,           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR10'

        !Kills Geometry
        call KillGeometry       (Me%ObjGeometry,      STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR20'

        !Kills HorizontalMap
        call KillHorizontalMap  (Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR30'

        !Kills Bathymetry
        call KillGridData       (Me%ObjBathymetry,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR40'

        !Kills HorizontalGrid
        call KillHorizontalGrid (Me%ObjHorizontalGrid,      STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR50'

        !Kills Compute Time
        call KillComputeTime    (Me%ObjTime,                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR60'            
          

        if (MonitorPerformance) then
            call KillWatchGroup (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillValida4D - ModuleValida4D - ERR10'
        endif

        call StopCPUTime

        call ShutdownMohid ("Valida4D", ElapsedSeconds, TotalCPUTime)
        !------------------------------------------------------------------------

    end subroutine KillValida4D
        

    !------------------------------------------------------------------------
    
    

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        
        call cpu_time(TotalCPUTime)

        ElapsedSeconds = FinalSystemTime - InitialSystemTime

    end subroutine StopCPUTime


end module ModuleValida4D









