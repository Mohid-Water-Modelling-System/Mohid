!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Smooth bathymetry for nesting
! PROGRAM       : Smooth bathymetry for nesting
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Paulo Chambel 
! DESCRIPTION   : Program to create a new sub-model batim to allow a smooth transition in the boundary between
!                 the coarser model and the high resolution model                  
!
!------------------------------------------------------------------------------

!DataFile  : "SmoothBatimNesting.dat"
!
!
!   FATHER_BATIM                : char              -           !Path to the bathymetry file of the father (or coarser) model
!   SON_BATIM                   : char              -           !Path to the bathymetry file of the son (or higher resolution) model
!   SMOOTH_COEF                 : char              -           !Path to the coefficient file use to smooth the bathymetry (0 - son model, 1 - Father model)
!   NEW_SON_BATIM               : char              -           !Path to the bathymetry file to be create

program SmoothBatimNesting
    
    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleDrawing
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleFillMatrix


    implicit none

    !Parameters----------------------------------------------------------------

    !Globals-------------------------------------------------------------------
    type     T_ExternalVar
        real, pointer, dimension(:,:)               :: XX_IE, YY_IE 
        real, pointer, dimension(:  )               :: XX, YY
        integer, pointer, dimension(:,:)            :: DefineCellsMap
        type (T_Size2D)                             :: WorkSize
        type (T_Size2D)                             :: Size
        logical                                     :: GridDistortion
        real                                        :: OriginX, OriginY, Rotation
        real,           dimension(:,:),    pointer  :: Batim
    end type T_ExternalVar

       
    type T_Global
        integer                                     :: ObjEnterData                 = 0

        integer                                     :: ObjGridFather                = 0
        integer                                     :: ObjGridSon                   = 0

        integer                                     :: ObjGridDataSon               = 0
        integer                                     :: ObjGridDataFather            = 0
        integer                                     :: ObjTime                      = 0
        integer                                     :: ObjHorizontalMap
        character(LEN=StringLength)                 :: SmoothBatimFile, FahterFile, SonFile, CoefFile
        real                                        :: LandPoint
        real,           dimension(:,:),    pointer  :: NewDepth, SmoothCoef
        type(T_PointF),  dimension(:,:),   pointer  :: GridPoint
        type(T_Polygon),                   pointer  :: Rect
        type(T_ExternalVar)                         :: ExtVarFather
        type(T_ExternalVar)                         :: ExtVarSon
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
        integer                             :: STAT_CALL
        character(len=*), parameter         :: ProjectFilePath  = "SmoothBatimNesting.dat"

        !Begin-----------------------------------------------------------------


        call StartupMohid ("Smooth Bathymetries for nesting")

        !Gets the actual time
        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%InitialSystemTime, float(Me%F95Time(1)), float(Me%F95Time(2)), &
                                                 float(Me%F95Time(3)), float(Me%F95Time(5)), &
                                                 float(Me%F95Time(6)), float(Me%F95Time(7))+ &
                                                 Me%F95Time(8)/1000.)


        call ConstructEnterData(Me%ObjEnterData, ProjectFilePath, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - SmoothBatimNesting - ERR01'

        call ReadGridFilesNames
        
        call Read_Lock_External_Var(Me%ObjGridFather, Me%ExtVarFather)

        call Read_Lock_External_Var(Me%ObjGridSon,    Me%ExtVarSon) 
       
        call ConstructPropertyCoefficients

        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - SmoothBatimNesting - ERR190'
        

    end subroutine OpenProject

    !----------------------------------------------------------------------    

    subroutine ReadGridFilesNames
        
        !Local-----------------------------------------------------------------
        type (T_Time)                       :: BeginTime, EndTime
        real                                :: DT
        integer                             :: flag, STAT_CALL
        
        
        !Begin-----------------------------------------------------------------

        !Get grid file path
        call GetData(Me%FahterFile,                                                     &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='FATHER_BATIM',                                      &
                     ClientModule ='SmoothBatimNesting',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR10'
        
        
        !Get grid file path
        call GetData(Me%SonFile,                                                        &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='SON_BATIM',                                         &
                     ClientModule ='SmoothBatimNesting',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR20'

        !Get grid file path
        call GetData(Me%CoefFile,                                                       &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='SMOOTH_COEF',                                       &
                     ClientModule ='SmoothBatimNesting',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR30'

        call GetData(Me%SmoothBatimFile,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='NEW_SON_BATIM',                                     &
                     ClientModule ='SmoothBatimNesting',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR40'

        call GetData(Me%LandPoint,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='LAND_POINT',                                        &
                     default      = -99.,                                               &
                     ClientModule ='SmoothBatimNesting',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR45'
        
        call StartComputeTime(Me%ObjTime, BeginTime,BeginTime, EndTime, DT, VariableDT = .false., STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR48'
        
        !Construct grids
        call ConstructHorizontalGrid(Me%ObjGridFather, trim(Me%FahterFile), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR50'

        call ConstructHorizontalGrid(Me%ObjGridSon, trim(Me%SonFile), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimNesting - ERR60'

        call ConstructGridData(Me%ObjGridDataFather, Me%ObjGridFather, FileName = trim(Me%FahterFile), STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - SmoothBatimNesting - ERR70'

        call ConstructGridData(Me%ObjGridDataSon, Me%ObjGridSon, FileName = trim(Me%SonFile), STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - SmoothBatimNesting - ERR80'
        
        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjGridDataSon, Me%ObjGridSon, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - SmoothBatimNesting - ERR90'
        

!        call ConstructGridData(Me%ObjGridDataCoef, Me%ObjGridSon, FileName = trim(Me%CoefFile), STAT = STAT_CALL)
!        if(STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - SmoothBatimNesting - ERR90'


    end subroutine ReadGridFilesNames

    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW

        !Begin-----------------------------------------------------------------


        allocate(Me%GridPoint(Me%ExtVarSon%Size%ILB:Me%ExtVarSon%Size%IUB,&
                              Me%ExtVarSon%Size%JLB:Me%ExtVarSon%Size%JUB))

        allocate(Me%NewDepth (Me%ExtVarSon%Size%ILB:Me%ExtVarSon%Size%IUB,&
                              Me%ExtVarSon%Size%JLB:Me%ExtVarSon%Size%JUB))



    end subroutine AllocateVariables
    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyCoefficients

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer, dimension(:,:  ), pointer      :: WaterPoints2D
        type (T_PropertyID)                     :: SmoothCoefID
        logical                                 :: BlockFound
        integer                                 :: ClientNumber
        integer                                 :: SizeILB, SizeIUB, SizeJLB, SizeJUB
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------
 
        SizeILB = Me%ExtVarSon%Size%ILB
        SizeIUB = Me%ExtVarSon%Size%IUB
        SizeJLB = Me%ExtVarSon%Size%JLB
        SizeJUB = Me%ExtVarSon%Size%JUB

        call GetWaterPoints2D(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - SmoothBatimNesting - ERR10'
        
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &
                                   "<begin_coef>", "<end_coef>", BlockFound,            &
                                   STAT = STAT_CALL)

        if(STAT_CALL .EQ. SUCCESS_)then

cd0:        if (BlockFound) then

                allocate(Me%SmoothCoef (SizeILB:SizeIUB, SizeJLB:SizeJUB))
                

                Me%SmoothCoef(:,:) = FillValueReal


                call ConstructFillMatrix  (PropertyID           = SmoothCoefID,         &
                                           EnterDataID          = Me%ObjEnterData,      &
                                           TimeID               = Me%ObjTime,           &
                                           HorizontalGridID     = Me%ObjGridSon,        &
                                           ExtractType          = FromBlock_,           &
                                           PointsToFill2D       = WaterPoints2D,        &
                                           Matrix2D             = Me%SmoothCoef,        &
                                           TypeZUV              = TypeZ_,               &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - SmoothBatimNesting - ERR20'


                if(.not. SmoothCoefID%SolutionFromFile)then

                    call KillFillMatrix(SmoothCoefID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - SmoothBatimNesting - ERR30'
                
                end if
                

            endif cd0
            
        endif          

        call UnGetHorizontalMap(Me%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyCoefficients - SmoothBatimNesting - ERR40'


    end subroutine ConstructPropertyCoefficients
    
    !----------------------------------------------------------------
    
    subroutine RunProject


        !Local-----------------------------------------------------------------
        integer                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
       
        write(*,*)"Running..."

        call GetGridData(Me%ObjGridDataFather, Me%ExtVarFather%Batim, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimNesting - ERR10'

        call GetGridData(Me%ObjGridDataSon,    Me%ExtVarSon%Batim,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimNesting - ERR20'

!        call GetGridData(Me%ObjGridDataCoef,   Me%SmoothCoef,         STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimNesting - ERR30'

        call AllocateVariables

        call DefineGridPoints

        call FillingCells

        call WriteSmoothBathym
        
        call Read_UnLock_External_Var(Me%ObjGridFather, Me%ExtVarFather)

        call Read_UnLock_External_Var(Me%ObjGridSon,    Me%ExtVarSon) 

        call UnGetGridData(Me%ObjGridDataFather, Me%ExtVarFather%Batim, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimNesting - ERR40'

        call UnGetGridData(Me%ObjGridDataSon,    Me%ExtVarSon%Batim,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimNesting - ERR50'

!        call UnGetGridData(Me%ObjGridDataCoef,   Me%SmoothCoef,         STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimNesting - ERR60'

    end subroutine RunProject
    
    !--------------------------------------------------------------------------
    
    subroutine WriteSmoothBathym
        
        !Local--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len=StringLength)         :: Coment1, Coment2 
        
        !Begin-----------------------------------------------------------------

        write(*,*)"Writing bathymetry..."

        Coment1 = 'File generated by'
        Coment2 = 'Mohid Smooth bathymetry for nesting'


        call WriteGridData(FileName         = trim(Me%SmoothBatimFile),         &
                           COMENT1          = Coment1,                          &
                           COMENT2          = Coment2,                          &
                           HorizontalGridID = Me%ObjGridSon,                    &
                           FillValue        = Me%LandPoint,                     &
                           Overwrite        = .true.,                           &
                           GridData2D_Real  = Me%NewDepth,                      &
                           STAT             = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_) stop 'WriteSmoothBathym - SmoothBatimNesting - ERR10'

    end subroutine WriteSmoothBathym
    
    
    !--------------------------------------------------------------------------



    subroutine Read_Lock_External_Var(ObjGrid, ExtVar) 

        !Arguments-------------------------------------------------------------
        integer                             :: ObjGrid
        type(T_ExternalVar)                 :: ExtVar
        
        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: GEOG, UTM, MIL_PORT, SIMPLE_GEOG
        integer                             :: GRID_COORD, CoordType

        !Begin-----------------------------------------------------------------

        call GetGridOrigin    (ObjGrid, ExtVar%OriginX, ExtVar%OriginY, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR00'

        call GetGridAngle     (ObjGrid, ExtVar%Rotation, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR01'

        call GetCheckDistortion(ObjGrid, ExtVar%GridDistortion, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR02'

        call GetHorizontalGridSize(ObjGrid, WorkSize = ExtVar%WorkSize, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR10'

        call GetHorizontalGridSize(ObjGrid, Size = ExtVar%Size, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR15'


        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,              &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD)

        !Gets Coordinates in use
        call GetGridCoordType(ObjGrid, CoordType, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR11'
        

        if    (CoordType == SIMPLE_GEOG)then
            
            call GetGridLatitudeLongitude(ObjGrid,                           &
                                          GridLatitudeConn  = ExtVar%YY_IE,  &
                                          GridLongitudeConn = ExtVar%XX_IE,  &
                                          STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR12'

        elseif(CoordType == UTM .or. CoordType == MIL_PORT .or. CoordType == GRID_COORD)then

            call GetHorizontalGrid(ObjGrid,                                  &
                                   XX_IE = ExtVar%XX_IE,                     &
                                   YY_IE = ExtVar%YY_IE,                     &
                                   STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR13'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR13'

        end if

        call GetHorizontalGrid(ObjGrid, XX = ExtVar%XX, YY = ExtVar%YY,  STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR40'

        call GetDefineCellsMap(ObjGrid, ExtVar%DefineCellsMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR50'

   
    end subroutine Read_Lock_External_Var

    !--------------------------------------------------------------------------



    subroutine Read_UnLock_External_Var(ObjGrid, ExtVar) 

        !Arguments-------------------------------------------------------------
        integer                             :: ObjGrid
        type(T_ExternalVar)                 :: ExtVar

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin-----------------------------------------------------------------


        call UngetHorizontalGrid(ObjGrid, ExtVar%XX_IE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR10'

        call UngetHorizontalGrid(ObjGrid, ExtVar%YY_IE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR20'

        call UngetHorizontalGrid(ObjGrid, ExtVar%XX, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR30'

        call UngetHorizontalGrid(ObjGrid, ExtVar%YY, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR40'

        call UngetHorizontalGrid(ObjGrid, ExtVar%DefineCellsMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimNesting - ERR50'

    end subroutine Read_UnLock_External_Var

    
    !--------------------------------------------------------------------------
    
    
    subroutine DefineGridPoints

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW

        !Begin-----------------------------------------------------------------

        write(*,*)"Defining grid points..."

        allocate(Me%GridPoint(Me%ExtVarSon%Size%ILB:Me%ExtVarSon%Size%IUB,&
                              Me%ExtVarSon%Size%JLB:Me%ExtVarSon%Size%JUB))

        do i = Me%ExtVarSon%WorkSize%ILB,  Me%ExtVarSon%WorkSize%IUB
        do j = Me%ExtVarSon%WorkSize%JLB , Me%ExtVarSon%WorkSize%JUB
                       
            XSW = Me%ExtVarSon%XX_IE(i, j)
            YSW = Me%ExtVarSon%YY_IE(i, j)
            XSE = Me%ExtVarSon%XX_IE(i, j + 1)
            YSE = Me%ExtVarSon%YY_IE(i, j + 1)
            XNE = Me%ExtVarSon%XX_IE(i + 1, j + 1)
            YNE = Me%ExtVarSon%YY_IE(i + 1, j + 1)
            XNW = Me%ExtVarSon%XX_IE(i + 1, j)
            YNW = Me%ExtVarSon%YY_IE(i + 1, j)

            Me%GridPoint(i,j)%X = (XSW+XNW+XNE+XSE) / 4.
            Me%GridPoint(i,j)%Y = (YSW+YNW+YNE+YSE) / 4.


        end do
        end do

    end subroutine DefineGridPoints

    
    !--------------------------------------------------------------------------


    subroutine FillingCells

        !Local-----------------------------------------------------------------
        type (T_PointF),   pointer          :: GridPoint
        real                                :: SumOfDepths
        integer                             :: nPointsInside, CurrentPoint
        integer                             :: i, j, ii, jj
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW
        logical                             :: flag

        !Begin-----------------------------------------------------------------

        flag = .false.
        
        allocate(Me%Rect)
        Me%Rect%Count = 5

        allocate(Me%Rect%VerticesF(1:Me%Rect%Count))

        write(*,*)"Filling cells with known data..."

        do i = Me%ExtVarSon%WorkSize%ILB,  Me%ExtVarSon%WorkSize%IUB
        do j = Me%ExtVarSon%WorkSize%JLB , Me%ExtVarSon%WorkSize%JUB
                       
idef:       if (Me%ExtVarSon%DefineCellsMap(i, j)==1 .and. Me%ExtVarSon%Batim(i, j) /= Me%LandPoint) then

                GridPoint => Me%GridPoint(i, j)

                flag = .false.

                do ii = Me%ExtVarFather%WorkSize%ILB,  Me%ExtVarFather%WorkSize%IUB
                do jj = Me%ExtVarFather%WorkSize%JLB , Me%ExtVarFather%WorkSize%JUB
                    
                    

id1:                if (Me%ExtVarFather%DefineCellsMap(ii, jj)==1 .and. Me%ExtVarFather%Batim(ii, jj) /= Me%LandPoint) then
                    
                        XSW = Me%ExtVarFather%XX_IE(ii,     jj)
                        YSW = Me%ExtVarFather%YY_IE(ii,     jj)
                        XSE = Me%ExtVarFather%XX_IE(ii,     jj + 1)
                        YSE = Me%ExtVarFather%YY_IE(ii,     jj + 1)
                        XNE = Me%ExtVarFather%XX_IE(ii + 1, jj + 1)
                        YNE = Me%ExtVarFather%YY_IE(ii + 1, jj + 1)
                        XNW = Me%ExtVarFather%XX_IE(ii + 1, jj)
                        YNW = Me%ExtVarFather%YY_IE(ii + 1, jj)
                    
                        Me%Rect%VerticesF(1)%X = XSW
                        Me%Rect%VerticesF(1)%Y = YSW
                        Me%Rect%VerticesF(2)%X = XSE
                        Me%Rect%VerticesF(2)%Y = YSE
                        Me%Rect%VerticesF(3)%X = XNE
                        Me%Rect%VerticesF(3)%Y = YNE
                        Me%Rect%VerticesF(4)%X = XNW
                        Me%Rect%VerticesF(4)%Y = YNW
                    
                        Me%Rect%VerticesF(5)%X = Me%Rect%VerticesF(1)%X
                        Me%Rect%VerticesF(5)%Y = Me%Rect%VerticesF(1)%Y
                    
                        call SetLimits(Me%Rect)

                        if(IsPointInsidePolygon(GridPoint, Me%Rect))then
                            Me%NewDepth(i, j) = Me%ExtVarFather%Batim(ii, jj) * (1. -  Me%SmoothCoef(i, j)) + &
                                                Me%ExtVarSon%Batim   (i,  j ) *        Me%SmoothCoef(i, j)
                            flag = .true.  
                            exit
                        end if

                    end if id1

                enddo
                    if (flag) exit
                enddo

                !If the ij son grid point lies outside the father's    &
                !grid then the son bathymetry value is assumed
                if(.not.flag) then
                    Me%NewDepth(i,j) = Me%ExtVarSon%Batim   (i,  j )
                end if

            else
            
                Me%NewDepth(i,j) = Me%ExtVarSon%Batim   (i,  j )
                
            endif idef     

        end do
        end do

        deallocate(Me%Rect%VerticesF)

    end subroutine FillingCells


    !--------------------------------------------------------------------------

    
    subroutine CloseProject

        !Local-----------------------------------------------------------------
        real                                        :: ElapsedSeconds
        real                                        :: TotalCPUTime
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

!        call KillGridData(Me%ObjGridDataCoef, STAT = STAT_CALL)
!        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimNesting - ERR25'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimNesting - ERR30'


        call KillGridData(Me%ObjGridDataSon, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimNesting - ERR30'

        call KillGridData(Me%ObjGridDataFather, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimNesting - ERR40'

        call KillHorizontalGrid(HorizontalGridID= Me%ObjGridSon, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimNesting - ERR10'

        call KillHorizontalGrid(HorizontalGridID= Me%ObjGridFather, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimNesting - ERR20'


        deallocate(Me%NewDepth, Me%GridPoint)

        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%FinalSystemTime,float(Me%F95Time(1)), float(Me%F95Time(2)),      &
                                              float(Me%F95Time(3)), float(Me%F95Time(5)),      &
                                              float(Me%F95Time(6)), float(Me%F95Time(7))+      &
                                              Me%F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = Me%FinalSystemTime - Me%InitialSystemTime

        call ShutdownMohid ("SmoothBathymNesting", ElapsedSeconds, TotalCPUTime)


    end subroutine CloseProject

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

end program SmoothBatimNesting

