!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToHDF5
! MODULE        : InterpolateTime
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod 
! DATE          : November 2010
! REVISION      : Paulo Leitão
! DESCRIPTION   : Module to Interpolate in time with data written in HDF5 format files
!
!------------------------------------------------------------------------------


Module ModuleInterpolateTime

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap
    use ModuleFunctions
    use ModuleFillMatrix

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartInterpolateTime
    private ::      ReadOptions
    private ::      ConstructGrid
    private ::      ConstructOutput
    private ::          Open_HDF5_OutPut_File
    private ::      ReadFieldsToInterpolate
    private ::          ReadNumberOfDataSources
    private ::      WriteHDF5File
    private ::          OutputInstants
    private ::          OutputGrid3D
    private ::          OutputFields2D
    private ::          OutputFields3D
    private ::      KillInterpolateTime


    
    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: prop_block_begin     = '<<beginproperty>>'
    character(LEN = StringLength), parameter    :: prop_block_end       = '<<endproperty>>'

    !Types---------------------------------------------------------------------

   
    type      T_OutW
        type(T_OutPutTime), dimension(:), pointer :: OutPutWindows
        logical                                   :: OutPutWindowsON            
        integer                                   :: WindowsNumber        
        integer,            dimension(:), pointer :: ObjHDF5 
    end type  T_OutW
    
   
    type       T_Property
        type (T_PropertyID)                                 :: ID
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        logical                                             :: Dim3D         
        integer,  dimension(:),     pointer                 :: IDdataSource
        integer                                             :: NdataSources = 0
    end type  T_Property

    type       T_InterpolateTime
        integer                                             :: ObjEnterData         = 0
        integer                                             :: ObjTime              = 0
        integer                                             :: ObjHorizontalGrid    = 0
        integer                                             :: ObjHorizontalMap     = 0
        integer                                             :: Objbathymetry        = 0
        integer                                             :: ObjGeometry          = 0
        integer                                             :: ObjMap               = 0
        character(len=PathLength)                           :: FileNameOut
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: GeometryFileName
        
        type (T_Time)                                       :: BeginTime
        type (T_Time)                                       :: EndTime        
        
        type(T_Size2D)                                      :: Size2D
        type(T_Size3D)                                      :: Size3D
        
        type(T_Size2D)                                      :: WorkSize2D
        type(T_Size3D)                                      :: WorkSize3D
        
        integer, dimension(:,:),    pointer                 :: WaterPoints2D
        integer, dimension(:,:, :), pointer                 :: WaterPoints3D
        
        integer                                             :: NumberDataSources
        type (T_Property), dimension (:), pointer           :: DataSource             
        integer,           dimension (:), pointer           :: IDDataSource                


        integer                                             :: NumberProperties
        type (T_Property), dimension (:), pointer           :: Property  
        integer,           dimension (:), pointer           :: NbyDataSource

        logical                                             :: NeedGeometry
        type (T_OutW)                                       :: OutW

    end type  T_InterpolateTime

    type(T_InterpolateTime), pointer                       :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartInterpolateTime(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)
        
        call ReadOptions(ClientNumber)

        call ConstructGrid
        
        call ReadFieldsToInterpolate(ClientNumber)
        
        call ConstructOutput

        call WriteHDF5File

        call KillInterpolateTime
        
        deallocate(Me)
        nullify(Me)

        STAT = SUCCESS_

    end subroutine StartInterpolateTime

    !------------------------------------------------------------------------
    
    subroutine ConstructGrid

        !Local-------------------------------------------------------------
        real,       dimension(:,:  ), pointer   :: SurfaceElevation
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------


        write(*,*)'Constructing grid...'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid - ModuleInterpolateTime - ERR10'


        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                                &
                                   WorkSize = Me%WorkSize2D,                            &
                                   Size     = Me%Size2D,                                &
                                   STAT     = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR20'

        call ConstructGridData      (GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     FileName         = Me%GridFileName,                &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR40'


        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     ActualTime       = Me%BeginTime,                   & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR50'
        
        if (Me%NeedGeometry) then

            call ConstructGeometry      (GeometryID       = Me%ObjGeometry,             &
                                         GridDataID       = Me%ObjBathymetry,           &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         HorizontalMapID  = Me%ObjHorizontalMap,        &
                                         ActualTime       = Me%BeginTime,               &
                                         NewDomain        = Me%GeometryFileName,        &
                                         STAT             = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR60'

            call GetGeometrySize(GeometryID     = Me%ObjGeometry,                       &
                                 Size           = Me%Size3D,                            &
                                 WorkSize       = Me%WorkSize3D,                        &
                                 STAT           = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR70'

            call ConstructMap ( Map_ID          = Me%ObjMap,                            &
                                GeometryID      = Me%ObjGeometry,                       &
                                HorizontalMapID = Me%ObjHorizontalMap,                  &
                                TimeID          = Me%ObjTime,                           &
                                STAT            = STAT_CALL)  
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR80'

            allocate(SurfaceElevation (Me%Size3D%ILB:Me%Size3D%IUB,                     &
                                       Me%Size3D%JLB:Me%Size3D%JUB))
            SurfaceElevation(:,:) = 0.

            call GetWaterPoints3D(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleInterpolateTime - ERR90'

            call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,              &
                                        WaterPoints3D    = Me%WaterPoints3D,            &
                                        SurfaceElevation = SurfaceElevation,            &
                                        ActualTime       = Me%BeginTime,                &
                                        STAT             = STAT_CALL )
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleInterpolateTime - ERR100'

            deallocate(SurfaceElevation)
            
        endif
        
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleInterpolateTime - ERR110'

        

    end subroutine ConstructGrid
    

    subroutine ConstructOutput
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iW

        !Begin-----------------------------------------------------------------


        
        Me%OutW%WindowsNumber = 1
        
        allocate(Me%OutW%ObjHDF5      (Me%OutW%WindowsNumber))
        allocate(Me%OutW%OutPutWindows(Me%OutW%WindowsNumber))
        
        Me%OutW%ObjHDF5(:) = 0
        
        do iW = 1, Me%OutW%WindowsNumber
        
            Me%OutW%OutPutWindows(iW)%ILB = Me%WorkSize2D%ILB
            Me%OutW%OutPutWindows(iW)%JLB = Me%WorkSize2D%JLB
            Me%OutW%OutPutWindows(iW)%IUB = Me%WorkSize2D%IUB
            Me%OutW%OutPutWindows(iW)%JUB = Me%WorkSize2D%JUB

            if (Me%NeedGeometry) then

                Me%OutW%OutPutWindows(iW)%KLB = Me%WorkSize3D%KLB
                Me%OutW%OutPutWindows(iW)%KUB = Me%WorkSize3D%KUB
                
            endif
            
            call GetOutPutTime(EnterDataID      = Me%ObjEnterData,                    & 
                                CurrentTime     = Me%BeginTime,                       &
                                EndTime         = Me%EndTime,                         &
                                keyword         = 'OUTPUT_TIME',                      &
                                SearchType      = FromBlock,                          &
                                OutPutsTime     = Me%OutW%OutPutWindows(1)%OutTime,   & 
                                OutPutsON       = Me%OutW%OutPutWindows(1)%ON,        & 
                                OutPutsNumber   = Me%OutW%OutPutWindows(1)%Number,    &                                         
                                STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutput - ModuleInterpolateTime - ERR10'
                     
            
            Me%OutW%OutPutWindows%NextOutPut = 1
            call Open_HDF5_OutPut_File(iW)

        enddo         
        
        
!        call GetOutPutTimeWindows(EnterDataID     = Me%ObjEnterData,                    & 
!                                  CurrentTime     = Me%BeginTime,                       &
!                                  EndTime         = Me%EndTime,                         &
!                                  OutPutWindows   = Me%OutW%OutPutWindows,              &   
!                                  OutPutWindowsON = Me%OutW%OutPutWindowsON,            & 
!                                  WindowsNumber   = Me%OutW%WindowsNumber,              &   
!                                  STAT            = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructOutput - ModuleInterpolateTime - ERR10'

            
!         if(Me%OutW%OutPutWindowsON)then
!
!            allocate(Me%OutW%ObjHDF5(Me%OutW%WindowsNumber))
!            
!            Me%OutW%ObjHDF5(:) = 0
!            
!            do iW = 1, Me%OutW%WindowsNumber
!            
!                if (Me%OutW%OutPutWindows(iW)%ILB < Me%WorkSize2D%ILB)                      &
!                    stop 'ConstructOutput - ModuleInterpolateTime - ERR20'
!
!                if (Me%OutW%OutPutWindows(iW)%IUB > Me%WorkSize2D%IUB)                      &
!                    stop 'ConstructOutput - ModuleInterpolateTime - ERR30'
!
!
!                if (Me%OutW%OutPutWindows(iW)%JLB < Me%WorkSize2D%JLB)                      &
!                    stop 'ConstructOutput - ModuleInterpolateTime - ERR40'
!
!                if (Me%OutW%OutPutWindows(iW)%JUB > Me%WorkSize2D%JUB)                      &
!                    stop 'ConstructOutput - ModuleInterpolateTime - ERR50'
!                    
!                if (Me%NeedGeometry) then
!
!                    if (Me%OutW%OutPutWindows(iW)%KLB < Me%WorkSize3D%KLB)                      &
!                        stop 'ConstructOutput - ModuleInterpolateTime - ERR60'
!
!                    if (Me%OutW%OutPutWindows(iW)%KUB > Me%WorkSize3D%KUB)                      &
!                        stop 'ConstructOutput - ModuleInterpolateTime - ERR70'
!                endif
!            
!                Me%OutW%OutPutWindows%NextOutPut = 1
!                call Open_HDF5_OutPut_File(iW)
!            enddo 
!        
!        else
!            call SetError(FATAL_, INTERNAL_, 'ConstructOutput - ModuleInterpolateTime - ERR20')
!        end if
!   
!   
!        
    end subroutine ConstructOutput
    
    !------------------------------------------------------------------------
    
    subroutine ReadOptions(ClientNumber)


        !Arguments---------------------------------------------------------------
        integer                                     :: ClientNumber
        !Local-----------------------------------------------------------------
        integer,  allocatable, dimension(:)         :: AuxInt, AuxInt2, AuxInt3, i
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1
        integer                                     :: ik, iP, iD, iPD

        !Begin-----------------------------------------------------------------
        
        write(*,*)'Reading instructions...'


        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'GRID_FILENAME',                                    &
                     ClientModule = 'ModuleInterpolateTime',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateTime - ERR10'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleInterpolateTime - ERR20'
        end if

        call GetData(Me%FileNameOut,                                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleInterpolateTime',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateTime - ERR30'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleInterpolateTime - ERR40'
        end if


        call GetData(Me%BeginTime,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'START',                                            &
                     ClientModule = 'ModuleInterpolateTime',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateTime - ERR50'

        call GetData(Me%EndTime,                                                        &
                     Me%ObjEnterData, iflag1,                                           &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'END',                                              &
                     ClientModule = 'ModuleInterpolateTime',                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateTime - ERR60'

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%BeginTime, Me%EndTime, 60., .false., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateTime - ERR70'

        Me%NumberDataSources = ReadNumberOfDataSources(ClientNumber)
        
        if (Me%NumberDataSources == -1) stop 'ReadOptions  - ModuleInterpolateTime - ERR80'
        
       
        allocate(Me%DataSource(Me%NumberDataSources))

        call ReadDataSourcesID(ClientNumber)

!------------------------------------------------------------------------------------------------------
        if(Me%NeedGeometry) then

            call GetData(Me%GeometryFileName,                                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'GEOMETRY_FILE',                                &
                         ClientModule = 'ModuleInterpolateTime',                        &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateTime - ERR90'


        endif
        

        if (Me%NumberDataSources > 1) then
        
            allocate(Me%IDDataSource(Me%NumberDataSources))
            
            allocate(AuxInt(Me%NumberDataSources), AuxInt2(Me%NumberDataSources), AuxInt3(Me%NumberDataSources)) 
            
            Me%IDDataSource  (:) = 0
                   
            !ID number of each DataSource (or data source)
            do iP = 1, Me%NumberDataSources
                Me%IDDataSource(iP) = Me%DataSource(iP)%ID%IDNumber
            enddo

            AuxInt3 = Me%IDDataSource
            
            allocate(i(Me%NumberDataSources)) 
            
            i(1:Me%NumberDataSources) = FillValueInt
            
            !Order ID properties in a descending way
            do iD = 1, Me%NumberDataSources
                i             = maxloc(AuxInt3)
                AuxInt2(iD)   = Me%IDDataSource(i(1))
                AuxInt (iD)   = i(1)
                AuxInt3(i(1)) = FillValueInt 
                i(:)          = FillValueInt                
            enddo
            
            Me%NumberProperties = 1
                    
            !Compute the number of different properties
            do iP = 2, Me%NumberDataSources
                if (AuxInt2(iP)/= AuxInt2(iP-1)) Me%NumberProperties= Me%NumberProperties+1
            enddo
            
        else
            Me%NumberProperties = Me%NumberDataSources
            
        
        endif
        
        allocate(Me%Property(Me%NumberProperties))
        
        
        ik=2
        do iP= 1, Me%NumberProperties 
            Me%Property(iP)%NdataSources = 1
            if (ik > Me%NumberDataSources) exit            
            do iD = ik, Me%NumberDataSources
                 if (AuxInt2(id)== AuxInt2(id-1)) then
                    Me%Property(iP)%NdataSources = Me%Property(iP)%NdataSources + 1
                 else
                    ik = iD + 1
                    exit
                 endif
            enddo
        enddo  
        
        id = 1
        do iP= 1, Me%NumberProperties 
            allocate(Me%Property(iP)%IDdataSource(Me%Property(iP)%NdataSources))
            do iPD = 1, Me%Property(iP)%NdataSources
                Me%Property(iP)%IDdataSource(iPD) = AuxInt(id)
                id = id + 1
            enddo
        enddo               

    end subroutine ReadOptions
    
    !--------------------------------------------------------------------------

    integer function ReadNumberOfDataSources (ClientNumber)
    
        !Arguments---------------------------------------------------------------
        integer                             :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, NumberOfDataSources
        logical                             :: BlockFound
        integer                             :: iflag
        

        !Sees how many  exists
        NumberOfDataSources = 0
        BlockFound = .true.
        Me%NeedGeometry = .false.
do1:    do while (BlockFound)
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                  &
                                       prop_block_begin, prop_block_end, BlockFound,   &
                                       STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) then
                NumberOfDataSources = -1
                exit
            endif
            if (BlockFound) then
            
                if (.not. Me%NeedGeometry) then
            
                    call GetData(Me%NeedGeometry, Me%ObjEnterData,  iflag,                      &
                                 SearchType     = FromBlockInBlock,                             &
                                 keyword        = 'DIM3D',                                      &
                                 default        = .true.,                                       &
                                 ClientModule   = 'ModuleInterpolateTime',                      &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadNumberOfDataSources - ModuleInterpolateTime - ERR10'
                    
                endif
            
                NumberOfDataSources = NumberOfDataSources + 1
            else
                call RewindBlock  (Me%ObjEnterData, ClientNumber) 
                exit
            endif
        enddo do1

        ReadNumberOfDataSources = NumberOfDataSources

    end function ReadNumberOfDataSources

    !--------------------------------------------------------------------------    
    !--------------------------------------------------------------------------

    subroutine ReadDataSourcesID (ClientNumber)
    
        !Arguments---------------------------------------------------------------
        integer                             :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iD
        logical                             :: BlockFound
        


        
 do1:    do iD = 1, Me%NumberDataSources
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                  &
                                       prop_block_begin, prop_block_end, BlockFound,   &
                                       STAT = STAT_CALL)
                                       
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataSourcesID - ModuleInterpolateTime - ERR10'
            
            if (.not. BlockFound     ) stop 'ReadDataSourcesID - ModuleInterpolateTime - ERR20'
                                       
            call ConstructPropertyID (Me%DataSource(iD)%ID, Me%ObjEnterData, FromBlockInBlock)                                      

        enddo do1
        
        call RewindBlock  (Me%ObjEnterData, ClientNumber) 
        
        

    end subroutine ReadDataSourcesID

    !--------------------------------------------------------------------------    

    !------------------------------------------------------------------------------------------------------

    subroutine ReadFieldsToInterpolate(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                             :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: iflag, STAT_CALL, iP, iD, iPD, iD1
        logical                             :: BlockFound

        !Begin-----------------------------------------------------------------
        
        
        if (Me%NumberDataSources < 1) stop 'ReadFieldsToInterpolate - ModuleInterpolateTime - ERR10'
        
        iP  = 1
        iPD = 1
    
        do iD = 1, Me%NumberDataSources 
        
            call ExtractBlockFromBlock(Me%ObjEnterData,                                 &
                                        ClientNumber,                                   &
                                        prop_block_begin,                               &
                                        prop_block_end,                                 &
                                        BlockFound,                                     &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToInterpolate - ModuleInterpolateTime - ERR20'
            
            if (.not.BlockFound) stop 'ReadFieldsToInterpolate - ModuleInterpolateTime - ERR30'

            
            call GetData(Me%DataSource(iD)%Dim3D, Me%ObjEnterData,  iflag,              &
                         SearchType     = FromBlockInBlock,                             &
                         keyword        = 'DIM3D',                                      &
                         default        = .true.,                                       &
                         ClientModule   = 'ModuleInterpolateTime',                      &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToInterpolate - ModuleInterpolateTime - ERR40'
            
            do iP = 1, Me%NumberProperties
                do iPD = 1, Me%Property(iP)%NdataSources
                    iD1 = Me%Property(iP)%IDdataSource(iPD) 
            
                    if (iD1==iD) then
                        Me%Property(iP)%ID    = Me%DataSource(iD)%ID
                        Me%Property(iP)%Dim3D = Me%DataSource(iD)%Dim3D
                        exit
                    endif
                enddo
                if (iD1==iD) then
                    exit
                endif                
            enddo

            if (Me%DataSource(iD)%Dim3D) then
            
                allocate(Me%DataSource(iD)%Values3D(Me%Size3D%ILB:Me%Size3D%IUB,        &
                                                    Me%Size3D%JLB:Me%Size3D%JUB,        &
                                                    Me%Size3D%KLB:Me%Size3D%KUB))
                if (iPD==1) then                                                          
                    allocate(Me%Property(iP)%Values3D  (Me%Size3D%ILB:Me%Size3D%IUB,    &
                                                        Me%Size3D%JLB:Me%Size3D%JUB,    &
                                                        Me%Size3D%KLB:Me%Size3D%KUB))
                    
                    Me%Property(iP)%Values3D(:,:,:) = FillValueReal
                
                endif


                !Uncovered cells are set to zero
                call ConstructFillMatrix  (PropertyID           = Me%DataSource(iD)%ID,  &
                                           EnterDataID          = Me%ObjEnterData,       &
                                           TimeID               = Me%ObjTime,            &
                                           HorizontalGridID     = Me%ObjHorizontalGrid,  &
                                           GeometryID           = Me%ObjGeometry,        &
                                           ExtractType          = FromBlockInBlock,      &
                                           PointsToFill3D       = Me%WaterPoints3D,      &
                                           Matrix3D             = Me%DataSource(iD)%Values3D,&
                                           TypeZUV              = TypeZ_,                &
                                           FillMatrix           = 0.,                    &             
                                           STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToInterpolate - ModuleInterpolateTime - ERR50'
                
            else                    


                allocate(Me%DataSource(iD)%Values2D(Me%Size2D%ILB:Me%Size2D%IUB,        &
                                                    Me%Size2D%JLB:Me%Size2D%JUB))

                if (iPD==1) then
                    allocate(Me%Property(iP)%Values2D  (Me%Size2D%ILB:Me%Size2D%IUB,    &
                                                        Me%Size2D%JLB:Me%Size2D%JUB))
                    Me%Property(iP)%Values2D(:,:) = FillValueReal
                endif

                !Uncovered cells are set to zero
                call ConstructFillMatrix  (PropertyID           = Me%DataSource(iD)%ID,  &
                                           EnterDataID          = Me%ObjEnterData,       &
                                           TimeID               = Me%ObjTime,            &
                                           HorizontalGridID     = Me%ObjHorizontalGrid,  &
                                           ExtractType          = FromBlockInBlock,      &
                                           PointsToFill2D       = Me%WaterPoints2D,      &
                                           Matrix2D             = Me%DataSource(iD)%Values2D,&
                                           TypeZUV              = TypeZ_,                &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToInterpolate - ModuleInterpolateTime - ERR60'
                                        
            endif
       
        enddo    
        
        
    end subroutine ReadFieldsToInterpolate

    !--------------------------------------------------------------------------
    
    subroutine WriteHDF5File
   
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Time)                       :: OutTime
        integer                             :: iW, io, iP, iD, iPD, STAT_CALL
        integer                             :: WorkILB, WorkIUB
        integer                             :: WorkJLB, WorkJUB
        integer                             :: WorkKLB, WorkKUB
        integer                             :: i, j, k

        !----------------------------------------------------------------------
        !Bounds



        do iW = 1, Me%OutW%WindowsNumber
        
            WorkILB = Me%OutW%OutPutWindows(iW)%ILB
            WorkIUB = Me%OutW%OutPutWindows(iW)%IUB

            WorkJLB = Me%OutW%OutPutWindows(iW)%JLB
            WorkJUB = Me%OutW%OutPutWindows(iW)%JUB

            WorkKLB = Me%OutW%OutPutWindows(iW)%KLB
            WorkKUB = Me%OutW%OutPutWindows(iW)%KUB        
            
            do io = 1, Me%OutW%OutPutWindows(iW)%Number

                OutTime = Me%OutW%OutPutWindows(iW)%OutTime(io)
                
                call ActualizeCurrentTime(TimeID = Me%ObjTime, DT_Global = 60.,             &
                                          Current = OutTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5File - ModuleInterpolateTime - ERR10'
                
                call OutputInstants(OutTime, iW, io)
                        
                if (Me%NeedGeometry) then
                    call OutputGrid3D(iw, io)
                endif
                
                do iP = 1, Me%NumberProperties
                    do iPD = 1, Me%Property(iP)%NdataSources
                        iD = Me%Property(iP)%IDdataSource(iPD)

                        if (Me%DataSource(iD)%Dim3D) then
                        
                            call ModifyFillMatrix (FillMatrixID   = Me%DataSource(iD)%ID%ObjFillMatrix,&
                                                   Matrix3D       = Me%DataSource(iD)%Values3D,         &
                                                   PointsToFill3D = Me%WaterPoints3D,                &
                                                   STAT           = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) &
                                stop 'WriteHDF5File - ModuleInterpolateTime - ERR10'

                            do k = WorkKLB, WorkKUB
                            do j = WorkJLB, WorkJUB
                            do i = WorkILB, WorkIUB
                                if (Me%WaterPoints3D(i, j, k) == WaterPoint) then
                                    if (Me%Property(iP)%Values3D(i, j, k) == FillValueReal) then
                                        if (Me%DataSource(iD)%Values3D(i, j, k) > FillValueReal/1000.) then
                                            Me%Property(iP)%Values3D(i,j,k) = Me%DataSource(iD)%Values3D(i,j,k)
                                        endif
                                    endif                                                            
                                endif
                            enddo
                            enddo
                            enddo
                            

                            if (iPD == Me%Property(iP)%NdataSources)  then
                                call OutputFields3D("/Results", iP, io, iW)
                                Me%Property(iP)%Values3D(:,:,:) = FillValueReal                                                  
                            endif                   
                    
                        else
                    
                            call ModifyFillMatrix (FillMatrixID   = Me%DataSource(iD)%ID%ObjFillMatrix,&
                                                   Matrix2D       = Me%DataSource(iD)%Values2D,         &
                                                   PointsToFill2D = Me%WaterPoints2D,                &
                                                   STAT           = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) &
                                stop 'WriteHDF5File - ModuleInterpolateTime - ERR20'
                                
                            do j = WorkJLB, WorkJUB
                            do i = WorkILB, WorkIUB
                                if (Me%WaterPoints2D(i, j) == WaterPoint) then
                                    if (Me%Property(iP)%Values2D(i, j) == FillValueReal) then
                                        if (Me%DataSource(iD)%Values2D(i, j) > FillValueReal/1000.) then
                                            Me%Property(iP)%Values2D(i,j) = Me%DataSource(iD)%Values2D(i,j)
                                        endif
                                    endif                                                            
                                endif
                            enddo
                            enddo
                                
                                
                            if (iPD == Me%Property(iP)%NdataSources) then                  
                                call OutputFields2D("/Results", iP, io, iW)  
                                Me%Property(iP)%Values2D(:,:) = FillValueReal                  
                            endif
                        endif
                    enddo                
                enddo
            enddo
        enddo
        
  
   
    end subroutine WriteHDF5File
   
    !--------------------------------------------------------------------------   

   subroutine Open_HDF5_OutPut_File(iW)
        
        !Arguments-------------------------------------------------------------
        integer, optional                           :: iW

        !Local-----------------------------------------------------------------
        real,    pointer, dimension(:, :   )        :: Bathymetry
        character (Len = PathLength)                :: FileName
        character (Len = StringLength)              :: AuxChar
        type(T_Size2D)                              :: WorkSize2D
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE, ObjHDF5, i, n, j

        !----------------------------------------------------------------------
        !Bounds

        FileName = Me%FileNameOut

        WorkILB = Me%OutW%OutPutWindows(iW)%ILB
        WorkIUB = Me%OutW%OutPutWindows(iW)%IUB

        WorkJLB = Me%OutW%OutPutWindows(iW)%JLB
        WorkJUB = Me%OutW%OutPutWindows(iW)%JUB

        WorkKLB = Me%OutW%OutPutWindows(iW)%KLB
        WorkKUB = Me%OutW%OutPutWindows(iW)%KUB
        
        WorkSize2D%ILB = Me%OutW%OutPutWindows(iW)%ILB            
        WorkSize2D%IUB = Me%OutW%OutPutWindows(iW)%IUB 
        
        WorkSize2D%JLB = Me%OutW%OutPutWindows(iW)%JLB 
        WorkSize2D%JUB = Me%OutW%OutPutWindows(iW)%JUB 
        
        if (Me%OutW%WindowsNumber > 1) then
        
            write(AuxChar,fmt='(i5)') iW
            Auxchar           = "_w"//trim(adjustl(Auxchar))//".hdf5"
            n                 = len_trim(Auxchar) 
            do j=1,len_trim(Filename)  
                if(FileName(j:j+4)==".hdf5") then
                    i = j
                    exit
                endif
            enddo
            FileName(i:i+n-1) = trim(Auxchar)
            
        endif
        
        call GetGridData      (Me%Objbathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR10'


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        ObjHDF5 = 0

        !Opens HDF5 File
        call ConstructHDF5      (ObjHDF5, trim(FileName),                               &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR20'
        
        Me%OutW%ObjHDF5(iW) = ObjHDF5
       
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5,                         &
                                 WorkSize = WorkSize2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR30'
        

        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                       &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR40'


        !Writes the Grid
        call HDF5WriteData   (ObjHDF5, "/Grid", "Bathymetry", "m",                      &
                              Array2D = Bathymetry,                                     &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR50'


        !Ungets the Bathymetry
        call UngetGridData (Me%Objbathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR60'

        if (Me%NeedGeometry) then

            call HDF5WriteData   (ObjHDF5, "/Grid", "WaterPoints3D", "-",               &
                                  Array3D = Me%WaterPoints3D,                           &
                                  STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR70'

        else

            call HDF5WriteData   (ObjHDF5, "/Grid", "WaterPoints2D", "-",               &
                                  Array2D = Me%WaterPoints2D,                           &
                                  STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR80'

        endif

        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleInterpolateTime - ERR90'


        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File
    !------------------------------------------------------------------------


    !------------------------------------------------------------------------
    
    subroutine OutputFields2D(RootGroup, iP, OutputNumber, iW)

        !Arguments-------------------------------------------------------------
        character(len=*)                                :: RootGroup
        integer                                         :: iP, OutputNumber, iW
        
        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, ObjHDF5
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        ILB = Me%OutW%OutPutWindows(iW)%ILB
        IUB = Me%OutW%OutPutWindows(iW)%IUB

        JLB = Me%OutW%OutPutWindows(iW)%JLB
        JUB = Me%OutW%OutPutWindows(iW)%JUB
       
        ObjHDF5 = Me%OutW%ObjHDF5(iW)

        call HDF5SetLimits  (ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields2D - ModuleInterpolateTime - ERR10'


        call HDF5WriteData(  ObjHDF5,                                                   &
                             trim(RootGroup)//"/"//trim(Me%Property(iP)%ID%Name),       &   
                             trim(Me%Property(iP)%ID%Name),                             &
                             trim(Me%Property(iP)%ID%Units),                            &
                             Array2D      = Me%Property(iP)%Values2D,                   &
                             OutputNumber = OutputNumber,                               &
                             STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields2D - ModuleInterpolateTime - ERR30'


        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields2D - ModuleInterpolateTime - ERR40'

    end subroutine OutputFields2D

    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
    
    subroutine OutputFields3D(RootGroup, iP, OutputNumber, iW)

        !Arguments-------------------------------------------------------------
        character(len=*)                                :: RootGroup
        integer                                         :: iP, OutputNumber, iW
        
        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: STAT_CALL, ObjHDF5

        !Begin-----------------------------------------------------------------
        
        ILB = Me%OutW%OutPutWindows(iW)%ILB
        IUB = Me%OutW%OutPutWindows(iW)%IUB

        JLB = Me%OutW%OutPutWindows(iW)%JLB
        JUB = Me%OutW%OutPutWindows(iW)%JUB

        KLB = Me%OutW%OutPutWindows(iW)%KLB
        KUB = Me%OutW%OutPutWindows(iW)%KUB
        
        ObjHDF5 = Me%OutW%ObjHDF5(iW)

        call HDF5SetLimits  (ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateTime - ERR10'


        call HDF5WriteData(  ObjHDF5,                                                   &
                             trim(RootGroup)//"/"//trim(Me%Property(iP)%ID%Name),       & 
                             trim(Me%Property(iP)%ID%Name),                             &
                             trim(Me%Property(iP)%ID%Units),                            &
                             Array3D      = Me%Property(iP)%Values3D,                   &
                             OutputNumber = OutputNumber,                               &
                             STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateTime - ERR30'


        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateTime - ERR40'

    end subroutine OutputFields3D

    !------------------------------------------------------------------------
    
    subroutine OutputGrid3D(iW, OutputNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: iW, OutputNumber
        
        !Local-----------------------------------------------------------------
        real   , dimension(:,:,:), pointer              :: SZZ
        integer                                         :: STAT_CALL, ObjHDF5
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB        

        !Begin-----------------------------------------------------------------

        ILB = Me%OutW%OutPutWindows(iW)%ILB
        IUB = Me%OutW%OutPutWindows(iW)%IUB

        JLB = Me%OutW%OutPutWindows(iW)%JLB
        JUB = Me%OutW%OutPutWindows(iW)%JUB

        KLB = Me%OutW%OutPutWindows(iW)%KLB
        KUB = Me%OutW%OutPutWindows(iW)%KUB
        
        ObjHDF5 = Me%OutW%ObjHDF5(iW)


        call HDF5SetLimits  (ObjHDF5, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateTime - ERR10'

        call GetGeometryDistances(Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateTime - ERR20'


        call HDF5WriteData(ObjHDF5,                                                     &
                           "/Grid/VerticalZ",                                           &
                           "Vertical",                                                  &
                           "m",                                                         &
                           Array3D      = SZZ,                                          &
                           OutputNumber = OutputNumber,                                 &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateTime - ERR30'


        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateTime - ERR40'

        call UnGetGeometry(Me%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateTime - ERR50'


    end subroutine OutputGrid3D


    !--------------------------------------------------------------------------



    subroutine OutputInstants(CurrentDate, iW, CurrentInstant)
        !Arguments-------------------------------------------------------------
        type(T_Time)                                    :: CurrentDate
        integer                                         :: iW, CurrentInstant
        
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, ObjHDF5



        !Begin-----------------------------------------------------------------

        ObjHDF5 = Me%OutW%ObjHDF5(iW)        
        
        call ExtractDate   (CurrentDate,                                                &
                            AuxTime(1), AuxTime(2), AuxTime(3),                         &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleInterpolateTime - ERR10'


        call HDF5WriteData  (ObjHDF5, "/Time",                                          &
                             "Time", "YYYY/MM/DD HH:MM:SS",                             &
                             Array1D = TimePtr,                                         &
                             OutputNumber = CurrentInstant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleInterpolateTime - ERR20'


    end subroutine OutputInstants

    
    !------------------------------------------------------------------------
   
    subroutine KillInterpolateTime
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, iW, iP, nUsers

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Killing InterpolateTime...'
        
        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillInterpolateTime - ModuleInterpolateGrids - ERR10' 
        

        if (Me%OutW%OutPutWindowsON)  then
        
            do iW = 1, Me%OutW%WindowsNumber
                
                call KillHDF5 (Me%OutW%ObjHDF5(iW), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR20'
                
                deallocate(Me%OutW%OutPutWindows(iW)%OutTime)
                
            enddo
            
            deallocate(Me%OutW%ObjHDF5      )
            deallocate(Me%OutW%OutPutWindows)
            
        endif              
        
        do iP = 1, Me%NumberDataSources
        
            if (Me%DataSource(iP)%Dim3D) then
                deallocate(Me%DataSource(iP)%Values3D)
            else
                deallocate(Me%DataSource(iP)%Values2D)
            endif
            
            call KillFillMatrix(Me%DataSource(iP)%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR30'

        
        enddo        
        
        deallocate(Me%DataSource)

        if(Me%NeedGeometry) then

            call UnGetMap(Me%ObjMap, Me%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleInterpolateTime - ERR40'

            call KillMap( Me%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR50'

            call KillGeometry( Me%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR60'
            
        endif
        
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR70'


        call KillHorizontalMap( Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR80'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR90'

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR100'
        
        call KillComputeTime(Me%ObjTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillInterpolateTime - ModuleInterpolateTime - ERR110'


    end subroutine KillInterpolateTime
    !------------------------------------------------------------------------

End Module ModuleInterpolateTime

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

