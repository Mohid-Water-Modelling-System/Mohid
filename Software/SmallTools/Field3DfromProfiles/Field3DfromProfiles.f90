!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Critical Erosion Shear 3D Field 
! PROGRAM       : Critical Erosion Shear 3D Field 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Paulo Chambel - v4.0
! DESCRIPTION   : Program to create 3D fields from profiles 
!
!------------------------------------------------------------------------------

!DataFile
!
!   BATIM_FILE                  : char              -           !Path to the bathymetry file to be created
!   OUTPUT_FILE                 : char              -           !Path to the outpufile 
!

program Field3DfromProfiles
    
    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap


    implicit none

    !Parameters----------------------------------------------------------------
    character(LEN = StringLength), parameter :: BeginProfile      = '<BeginProfile>'
    character(LEN = StringLength), parameter :: EndProfile        = '<EndProfile>'
    character(LEN = StringLength), parameter :: BeginLayers       = '<BeginLayers>'
    character(LEN = StringLength), parameter :: EndLayers         = '<EndLayers>'
    character(LEN = StringLength), parameter :: BeginProfileID    = '<BeginProfile_ID>'
    character(LEN = StringLength), parameter :: EndProfileID      = '<EndProfile_ID>'

    integer,                       parameter :: Profile           = 1
    integer,                       parameter :: Layers            = 2


    integer,                       parameter :: NoProfile         = -99

    !Globals-------------------------------------------------------------------
            
    type     T_ExternalVar
        real,           dimension(:,:,:),  pointer  :: SZZ      
        integer,        dimension(:,:,:),  pointer  :: PointsToFill3D
        real,           dimension(:,:),    pointer  :: Bathymetry, AuxGrid
        type (T_Size3D)                             :: WorkSize
        type (T_Size3D)                             :: Size
    end type T_ExternalVar
            
    type     T_Profile
        integer                                     :: Hmin, Hmax, Ndepths
        real,  dimension(:), pointer                :: Values, Depths
    end type T_Profile

    type     T_ProfileList
        integer                                     :: Number
        type(T_Profile),  dimension(:), pointer     :: Profiles
    end type T_ProfileList


    type T_Global
        integer                                     :: ObjHorizontalGrid            = 0
        integer                                     :: ObjHorizontalMap             = 0
        integer                                     :: ObjGeometry                  = 0
        integer                                     :: ObjMap                       = 0
        integer                                     :: ObjBathymetry                = 0
        integer                                     :: ObjAuxGridData               = 0
        integer                                     :: ObjEnterData                 = 0
        integer,        dimension(:,:  ),  pointer  :: ProfileID
        real,           dimension(:,:,:),  pointer  :: Field3D

        character(len=PathLength)                   :: BatimFilePath
        character(len=PathLength)                   :: AuxGridDataFilePath
        character(len=PathLength)                   :: GeometryFilePath
        character(len=PathLength)                   :: Field3DFilePath

        integer                                     :: InputDataType

        type(T_ProfileList)                         :: ProfileList
        type(T_ExternalVar)                         :: ExtVar

        type (T_Time)                               :: InitialSystemTime, FinalSystemTime
        integer, dimension(8)                       :: F95Time
    end type T_Global                                 
                                                     
    type (T_Global)                                 :: Me

    !Begin---------------------------------------------------------------------
    call OpenProject
    call RunProject  
    call CloseProject

    contains
    
    !--------------------------------------------------------------------------
    
    subroutine OpenProject
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(PathLength), parameter            :: ProjectFilePath  = "Field3DfromProfiles.dat"

        !Begin-----------------------------------------------------------------


        call StartupMohid ("Construct Field3D from Profiles")

        !Gets the actual time
        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%InitialSystemTime, float(Me%F95Time(1)), float(Me%F95Time(2)), &
                                                 float(Me%F95Time(3)), float(Me%F95Time(5)), &
                                                 float(Me%F95Time(6)), float(Me%F95Time(7))+ &
                                                 Me%F95Time(8)/1000.)


        call ConstructEnterData(Me%ObjEnterData, ProjectFilePath, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - Field3DfromProfiles - ERR10'

        call ReadGridFilesNames

        call ConstructAuxiliarClasses

        
        allocate(Me%Field3D (Me%ExtVar%Size%ILB:Me%ExtVar%Size%IUB,               &
                             Me%ExtVar%Size%JLB:Me%ExtVar%Size%JUB,               &
                             Me%ExtVar%Size%KLB:Me%ExtVar%Size%KUB))

        allocate(Me%ProfileID(Me%ExtVar%Size%ILB:Me%ExtVar%Size%IUB,              &
                              Me%ExtVar%Size%JLB:Me%ExtVar%Size%JUB))

        Me%Field3D  (:,:,:) = FillValueReal
        Me%ProfileID(:,:  ) = FillValueReal

        call ConstructProfileList

        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - Field3DfromProfiles - ERR20'
        

    end subroutine OpenProject
  
 
    !--------------------------------------------------------------------------

    subroutine ReadGridFilesNames
        
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL

        !Begin-----------------------------------------------------------------


        !Open bathymetry file
        call GetData(Me%BatimFilePath,                                                  &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='BATIM_FILE',                                        &
                     ClientModule ='Field3DfromProfiles',                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - Field3DfromProfiles - ERR10'


        !Open bathymetry file
        call GetData(Me%AuxGridDataFilePath,                                            &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='AUX_FILE',                                          &
                     Default      = trim(Me%BatimFilePath),                             &
                     ClientModule ='Field3DfromProfiles',                               &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - Field3DfromProfiles - ERR20'


        !Geometry 3D file
        call GetData(Me%GeometryFilePath,                                               &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     ClientModule ='Field3DfromProfiles',                               &
                     keyword      ='GEOMETRY_FILE',                                     &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - Field3DfromProfiles - ERR30'

 
        !Field 3D file
        call GetData(Me%Field3DFilePath,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     ClientModule ='Field3DfromProfiles',                               &
                     keyword      ='FIED3D_FILE',                                       &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - Field3DfromProfiles - ERR30'


    end subroutine ReadGridFilesNames

    !--------------------------------------------------------------------------

    subroutine ConstructAuxiliarClasses

        !Local-----------------------------------------------------------------
        real,   pointer,  dimension(:,:)            :: SurfaceElevation
        type (T_Time)                               :: ActualTime        
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Construct grid
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, trim(Me%BatimFilePath), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - Field3DfromProfiles - ERR10'


        call ConstructGridData      (GridDataID       = Me%ObjBathymetry,        &
                                     HorizontalGridID = Me%ObjHorizontalGrid,    &
                                     FileName         = Me%BatimFilePath,        &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR20'

        call ConstructGridData      (GridDataID       = Me%ObjAuxGridData,       &
                                     HorizontalGridID = Me%ObjHorizontalGrid,    &
                                     FileName         = Me%AuxGridDataFilePath,  &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR25'


        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,     &
                                     GridDataID       = Me%ObjBathymetry,        &
                                     HorizontalGridID = Me%ObjHorizontalGrid,    &
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR30'


        call ConstructGeometry      (GeometryID       = Me%ObjGeometry,         &
                                     GridDataID       = Me%ObjBathymetry,       &
                                     HorizontalGridID = Me%ObjHorizontalGrid,   &
                                     HorizontalMapID  = Me%ObjHorizontalMap,    &
                                     NewDomain        = Me%GeometryFilePath,    &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR40'

        call GetGeometrySize(GeometryID     = Me%ObjGeometry,                   &
                             Size           = Me%ExtVar%Size,                   &
                             WorkSize       = Me%ExtVar%WorkSize,               &
                             STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR50'

        call ConstructMap ( Map_ID          = Me%ObjMap,                        &
                            GeometryID      = Me%ObjGeometry,                   &
                            HorizontalMapID = Me%ObjHorizontalMap,              &
                            STAT            = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR60'

        allocate(SurfaceElevation (Me%ExtVar%Size%ILB:Me%ExtVar%Size%IUB,               &
                                   Me%ExtVar%Size%JLB:Me%ExtVar%Size%JUB))
        SurfaceElevation(:,:) = 0

        call GetWaterPoints3D(Me%ObjMap, Me%ExtVar%PointsToFill3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - ModuleInterpolateGrids - ERR70'

        call ComputeInitialGeometry(GeometryID      = Me%ObjGeometry,                   &
                                    WaterPoints3D   = Me%ExtVar%PointsToFill3D,         &
                                    SurfaceElevation= SurfaceElevation,                 &
                                    STAT            = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames -  ModuleInterpolateGrids - ERR80'

        deallocate(SurfaceElevation)

        call UnGetMap(Me%ObjMap, Me%ExtVar%PointsToFill3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - ModuleInterpolateGrids - ERR90'

       
    end subroutine ConstructAuxiliarClasses

    !--------------------------------------------------------------------------

    subroutine ConstructProfileList

        !Arguments-------------------------------------------------------------
        integer                                     :: NPoints
        real,    dimension(: ), pointer             :: Values, Depths

        !Local----------------------------------------------------------------
        character(PathLength)                       :: FileName
        character(StringLength)                     :: AuxString
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: NDEPTHS
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: n
        real, dimension(:), allocatable             :: Aux 

        !----------------------------------------------------------------------
        !Begin----------------------------------------------------------------

        n = 0

        do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        BeginProfileID, EndProfileID, BlockFound,       &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR10'

BF:         if (BlockFound) then

                n = n + 1

            else BF

                Me%ProfileList%Number = n

                if (Me%ProfileList%Number == 0) then
                    write (*,*) 'Zero profiles defined'
                    stop 'ConstructProfileList - Field3DfromProfiles - ERR20'
                endif

                allocate(Me%ProfileList%Profiles(1:Me%ProfileList%Number))

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR30'

                exit

            endif BF

        enddo

        do n = 1, Me%ProfileList%Number

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        BeginProfileID, EndProfileID, BlockFound,       &
                                        STAT = STAT_CALL)
        
            if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR40'

            call GetData(AuxString,                                                     &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'INPUT_DATA_TYPE',                              &
                         ClientModule = 'Field3DfromProfiles',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR50'

            select case (trim(adjustl(AuxString)))
                case ("Layers",     "LAYERS",     "layers")
                    Me%InputDataType = Layers
                case ("Profile",    "PROFILE",    "profile")
                    Me%InputDataType = Profile
                case default
                    write(*,*)'Invalid option for keyword INITIALIZATION_METHOD'
                    stop 'ReadOptions - ModuleFillMatrix - ERR04'
            end select


            call GetData(FileName,                                                      &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'FILENAME',                                     &
                         ClientModule = 'Field3DfromProfiles',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR50'

            call GetData(Me%ProfileList%Profiles(n)%Hmin,                               &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'HMIN',                                         &
                         ClientModule = 'Field3DfromProfiles',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR60'

            call GetData(Me%ProfileList%Profiles(n)%Hmax,                               &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'HMAX',                                         &
                         ClientModule = 'Field3DfromProfiles',                          &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR70'

            if      (Me%InputDataType == Layers ) then

                call ConstructSpaceLayers  (FileName, Me%ProfileList%Profiles(n)%Values)

            else if (Me%InputDataType == Profile) then 

                call ConstructSpaceProfile (FileName, Me%ProfileList%Profiles(n)%NDepths, &
                                                      Me%ProfileList%Profiles(n)%Values,  &
                                                      Me%ProfileList%Profiles(n)%Depths)
            endif

        enddo 

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProfileList - Field3DfromProfiles - ERR80'



    end subroutine ConstructProfileList

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine ConstructSpaceProfile (FileName, NDepths, Values, Depths)

        !Arguments-------------------------------------------------------------
        character(PathLength)                       :: FileName
        integer                                     :: NDepths
        real,    dimension(: ), pointer             :: Values, Depths

        !Local----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: line, l
        real, dimension(:), allocatable             :: Aux 

        !----------------------------------------------------------------------
        !Begin----------------------------------------------------------------


        !Opens File
        call ConstructEnterData(ObjEnterData, FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR10'


        call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                          &
                                    BeginProfile, EndProfile, BlockFound,                &
                                    FirstLine = FirstLine, LastLine = LastLine,          &
                                    STAT = STAT_CALL)
    
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR20'

BF:     if (BlockFound) then

            NDEPTHS =  LastLine - FirstLine - 1

            !Allocates auxiliar variables
            allocate (Values        (NDEPTHS))
            allocate (Depths        (NDEPTHS))

            Values(:) = FillValueReal
            Depths(:) = FillValueReal

            allocate (Aux(2))
        
            l = 1
            do line = FirstLine + 1, LastLine - 1

                call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                             SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR30'

                Depths(l) = Aux(1)
                Values(l) = Aux(2)
                l = l + 1

            enddo

            deallocate(Aux)

            call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR40'

            call KillEnterData(ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR50'


        else 

            write(*,*) 'Block <BeginProfile>, <EndProfile> not found'
            write(*,*) 'FileName = ', trim(FileName)
            stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR60'

        endif BF


    end subroutine ConstructSpaceProfile

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ConstructSpaceLayers (FileName, Values)

        !Arguments-------------------------------------------------------------
        character(PathLength)                       :: FileName
        real,    dimension(: ), pointer             :: Values

        !Local----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber
        integer                                     :: FirstLine, LastLine
        logical                                     :: BlockFound
        integer                                     :: line, l
        real, dimension(:), allocatable             :: Aux 

        !----------------------------------------------------------------------
        !Begin----------------------------------------------------------------


        !Opens File
        call ConstructEnterData(ObjEnterData, FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR10'


        call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                         &
                                    BeginLayers, EndLayers, BlockFound,                 &
                                    FirstLine = FirstLine, LastLine = LastLine,         &
                                    STAT = STAT_CALL)

BF:     if (BlockFound) then

            allocate (Values(Me%ExtVar%Size%KLB : Me%ExtVar%Size%KUB))

            Values(:) = FillValueReal

            allocate (Aux(2))
        
            l = 1
            do line = FirstLine + 1, LastLine - 1

                call GetData(Aux, EnterDataID = ObjEnterData, flag = iflag, &
                             SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR30'

                if (Aux(1) < Me%ExtVar%Worksize%KLB .or. Aux(1) > Me%ExtVar%Worksize%KUB) then
                    stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR35'
                endif 
                Values(Aux(1)) = Aux(2)
                l = l + 1

            enddo

            deallocate(Aux)

            call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR40'

            call KillEnterData(ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR50'


        else  BF

            write(*,*) 'Block <BeginLayers>, <EndLayers> not found'
            write(*,*) 'FileName = ', trim(FileName)
            stop 'ConstructSpaceProfile - Field3DfromProfiles - ERR60'

        endif BF


    end subroutine ConstructSpaceLayers

    !--------------------------------------------------------------------------

    
    subroutine RunProject

    !Begin-----------------------------------------------------------------
       
        write(*,*)"Running..."

        call Read_Lock_External_Var

        call BuildField3D

        call WriteField3D          

        call Read_UnLock_External_Var

    end subroutine RunProject

    !--------------------------------------------------------------------------
    
    subroutine BuildField3D

        !Local-----------------------------------------------------------------
        real        :: CellDepth
        integer     :: k, KLB, KUB
        integer     :: j, JLB, JUB
        integer     :: i, ILB, IUB, n


        !Begin-----------------------------------------------------------------
       
        Me%ProfileID(:,:) = NoProfile

        do n  = 1, Me%ProfileList%Number

            do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
            do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                if (Me%ExtVar%PointsToFill3D(i, j, Me%ExtVar%WorkSize%KUB) == WaterPoint) then

                    if (Me%ExtVar%AuxGrid(i,j) >= Me%ProfileList%Profiles(n)%Hmin .and. &
                        Me%ExtVar%AuxGrid(i,j) <= Me%ProfileList%Profiles(n)%Hmax) then

                        Me%ProfileID(i,j) = n

                    endif

                endif

            enddo
            enddo

        enddo

        if (Me%InputDataType == Profile) then

            do k = Me%ExtVar%WorkSize%KLB, Me%ExtVar%WorkSize%KUB
            do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
            do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                if (Me%ExtVar%PointsToFill3D(i, j, k) == WaterPoint) then

                    CellDepth          = (Me%ExtVar%SZZ(i, j, k) + Me%ExtVar%SZZ(i, j, k-1))&
                                          / 2. - Me%ExtVar%SZZ(i, j, Me%ExtVar%WorkSize%KUB)
                                                     
                    n = Me%ProfileID(i,j)

                    if (n /= NoProfile) then
                                       
                        Me%Field3D(i,j,k)  = InterpolateProfile (CellDepth,                     &
                                             Me%ProfileList%Profiles(n)%NDepths,                &
                                             Me%ProfileList%Profiles(n)%Values,                 &
                                             Me%ProfileList%Profiles(n)%Depths)
                    else

                        Me%Field3D(i,j,k)  =  FillValueReal

                    endif

                endif

            enddo
            enddo
            enddo

        else if (Me%InputDataType == Layers) then

            do k = Me%ExtVar%WorkSize%KLB, Me%ExtVar%WorkSize%KUB
            do j = Me%ExtVar%WorkSize%JLB, Me%ExtVar%WorkSize%JUB
            do i = Me%ExtVar%WorkSize%ILB, Me%ExtVar%WorkSize%IUB

                if (Me%ExtVar%PointsToFill3D(i, j, k) == WaterPoint) then

                    n = Me%ProfileID(i,j)

                    if (n /= NoProfile) then
                                       
                        Me%Field3D(i,j,k) = Me%ProfileList%Profiles(n)%Values(k)

                    else

                        Me%Field3D(i,j,k) = FillValueReal

                    endif
                                         
                endif

            enddo
            enddo
            enddo

        endif

   
    end subroutine BuildField3D

    !--------------------------------------------------------------------------
    
    subroutine WriteField3D
        
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len=StringLength)         :: Coment1, Coment2 
        
        !Begin-----------------------------------------------------------------

        write(*,*)"Writing Field 3D..."

        Coment1 = 'File generated by'
        Coment2 = 'Mohid - Construct Field3D from Profiles'


        call WriteGridData(FileName         = trim(Me%Field3DFilePath),         &
                           COMENT1          = Coment1,                          &
                           COMENT2          = Coment2,                          &
                           HorizontalGridID = Me%ObjHorizontalGrid,             &
                           FillValue        = -99.,                             &
                           Overwrite        = .true.,                           &
                           GridData3D_Real  = Me%Field3D,                       &
                           KLB              = Me%ExtVar%WorkSize%KLB,           &
                           KUB              = Me%ExtVar%WorkSize%KUB,           &
                           STAT             = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_) stop 'WriteField3D - Field3DfromProfiles - ERR10'

    end subroutine WriteField3D
    
    
    !--------------------------------------------------------------------------

 
    subroutine Read_Lock_External_Var
        
        !External-----------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetGridData (Me%ObjBathymetry, Me%ExtVar%Bathymetry,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_Lock_External_Var - Field3DfromProfiles - ERR10'

        call GetGridData (Me%ObjAuxGridData, Me%ExtVar%AuxGrid,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_Lock_External_Var - Field3DfromProfiles - ERR20'

        !Gets Center of the cells
        call GetGeometryDistances(Me%ObjGeometry, SZZ = Me%ExtVar%SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_Lock_External_Var - Field3DfromProfiles - ERR30'

        call GetWaterPoints3D(Me%ObjMap, Me%ExtVar%PointsToFill3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Read_Lock_External_Var - Field3DfromProfiles - ERR40'


    end subroutine Read_Lock_External_Var

    !--------------------------------------------------------------------------


    subroutine Read_UnLock_External_Var

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call UnGetGridData (Me%ObjBathymetry, Me%ExtVar%Bathymetry,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_UnLock_External_Var - Field3DfromProfiles - ERR10'

        call UnGetGridData (Me%ObjAuxGridData, Me%ExtVar%AuxGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_UnLock_External_Var - Field3DfromProfiles - ERR20'

        call UnGetGeometry (Me%ObjGeometry, Me%ExtVar%SZZ,        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_UnLock_External_Var - Field3DfromProfiles - ERR30'

        call UnGetMap      (Me%ObjMap, Me%ExtVar%PointsToFill3D,   STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Read_UnLock_External_Var - Field3DfromProfiles - ERR40'

    end subroutine Read_UnLock_External_Var

    
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    
    subroutine CloseProject

        !Local-----------------------------------------------------------------
        real                                        :: ElapsedSeconds
        real                                        :: TotalCPUTime
        integer                                     :: STAT_CALL, n

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Killing father grid...'

        call KillMap( Me%ObjMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR10'

        call KillGeometry( Me%ObjGeometry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR20'

        call KillHorizontalMap( Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR30'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR40'

        call KillGridData(Me%ObjAuxGridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR50'

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR60'


        deallocate(Me%Field3D)

        deallocate(Me%ProfileID)

        deallocate(Me%ProfileList%Profiles)

        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%FinalSystemTime,float(Me%F95Time(1)), float(Me%F95Time(2)),      &
                                              float(Me%F95Time(3)), float(Me%F95Time(5)),      &
                                              float(Me%F95Time(6)), float(Me%F95Time(7))+      &
                                              Me%F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = Me%FinalSystemTime - Me%InitialSystemTime

        call ShutdownMohid ("Construct Field3D from Profiles", ElapsedSeconds, TotalCPUTime)


    end subroutine CloseProject

    !--------------------------------------------------------------------------
    
end program Field3DfromProfiles

