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

!DataFile  : "SmoothBatimCoef.dat"
!
!
!   SON_BATIM                   : char              -           !Path to the bathymetry file of the son (or higher resolution) model
!   NEW_COEFS                   : char              -           !Path to the bathymetry file to be create
!   NUMBER_CELLS                : int               -           !Number of cells contemplated from the border in.
!   NUMBER_CELLS                : int               -           !Number of cells contemplated from the border in.

program SmoothBatimCoef
    
    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleDrawing
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleGridData

    implicit none

    !Parameters----------------------------------------------------------------

    !Globals-------------------------------------------------------------------
    type     T_ExternalVar
        real, pointer, dimension(:,:)               :: XX_IE, YY_IE 
        real, pointer, dimension(:  )               :: XX, YY
        integer, pointer, dimension(:,:)            :: DefineCellsMap
        type (T_Size2D)                             :: WindowSize
        type (T_Size2D)                             :: WorkSize
        type (T_Size2D)                             :: Size
        logical                                     :: GridDistortion
        real                                        :: OriginX, OriginY, Rotation
        real,           dimension(:,:),    pointer  :: Batim
    end type T_ExternalVar
       
    type T_Global
        integer                                     :: ObjEnterData                 = 0

        integer                                     :: ObjGridSon                   = 0

        integer                                     :: ObjGridDataSon               = 0
        character(LEN=StringLength)                 :: SonFile, CoefFile
        real                                        :: LandPoint
        real,           dimension(:,:),    pointer  :: SmoothCoef
        integer                                     :: ncells
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
        character(len=*), parameter         :: ProjectFilePath  = "SmoothBatimCoef.dat"

        !Begin-----------------------------------------------------------------


        call StartupMohid ("Smooth Bathymetries Ceofficients")

        !Gets the actual time
        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%InitialSystemTime, float(Me%F95Time(1)), float(Me%F95Time(2)), &
                                                 float(Me%F95Time(3)), float(Me%F95Time(5)), &
                                                 float(Me%F95Time(6)), float(Me%F95Time(7))+ &
                                                 Me%F95Time(8)/1000.)


        call ConstructEnterData(Me%ObjEnterData, ProjectFilePath, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - SmoothBatimCoef - ERR01'

        call ReadGridFilesNames

        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'OpenProject - SmoothBatimCoef - ERR190'
        

    end subroutine OpenProject

    !----------------------------------------------------------------------    

    subroutine ReadGridFilesNames
        
        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        integer, dimension(:), pointer      :: aux
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------
                
        !Get grid file path
        call GetData(Me%SonFile,                                                        &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='SON_BATIM',                                         &
                     ClientModule ='SmoothBatimCoef',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimCoef - ERR20'

        call GetData(Me%CoefFile,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='NEW_COEF',                                     &
                     ClientModule ='SmoothBatimCoef',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimCoef - ERR40'

        call GetData(Me%Ncells,                                                &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='NUMBER_CELLS',                               &
                     ClientModule ='SmoothBatimCoef',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimCoef - ERR40'

        call GetData(Me%LandPoint,                                                      &
                     Me%ObjEnterData, flag,                                             &
                     SearchType   = FromFile_,                                          &
                     keyword      ='LAND_POINT',                                        &
                     default      = -99.,                                               &
                     ClientModule ='SmoothBatimCoef',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimCoef - ERR45'

        !Read Smooth_Window ILB, IUB, JLB, JUB data
        allocate (aux(4))

        call GetData(aux,                                                               &
                     Me%ObjEnterData, flag,                                             &
                     SearchType = FromFile,                                             &
                     keyword    = 'SMOOTH_WINDOW',                                      &
                     Default    = FillValueInt,                                         &                                           
                     ClientModule ='SmoothBathimCoefs',                                 &
                     STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, KEYWORD_, 'ReadGridFilesNames - SmoothBathimCoefs - ERR02')

        if   (flag == 0) then

            Me%ExtVarSon%WindowSize%ILB = 0
            Me%ExtVarSon%WindowSize%IUB = 0

            Me%ExtVarSon%WindowSize%JLB = 0
            Me%ExtVarSon%WindowSize%JUB = 0

        else if (flag == 4) then

            Me%ExtVarSon%WindowSize%ILB = aux(1)
            Me%ExtVarSon%WindowSize%IUB = aux(2)
            Me%ExtVarSon%WindowSize%JLB = aux(3)
            Me%ExtVarSon%WindowSize%JUB = aux(4)

        else

            call SetError(FATAL_, KEYWORD_, 'ReadGridFilesNames - SmoothBathimCoefs - ERR03')

        endif

        deallocate (aux)
        
        !Construct grids
        call ConstructHorizontalGrid(Me%ObjGridSon, trim(Me%SonFile), STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ReadGridFilesNames - SmoothBatimCoef - ERR60'

        call ConstructGridData(Me%ObjGridDataSon, Me%ObjGridSon, FileName = trim(Me%SonFile), STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'ReadGridFilesNames - SmoothBatimCoef - ERR80'

    end subroutine ReadGridFilesNames


    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Local-----------------------------------------------------------------
        integer                             :: i, j
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW

        !Begin-----------------------------------------------------------------


        allocate(Me%SmoothCoef (Me%ExtVarSon%Size%ILB:Me%ExtVarSon%Size%IUB,&
                                Me%ExtVarSon%Size%JLB:Me%ExtVarSon%Size%JUB))



    end subroutine AllocateVariables
    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------

    subroutine RunProject


        !Local-----------------------------------------------------------------
        integer                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
       
        write(*,*)"Running..."

        call Read_Lock_External_Var(Me%ObjGridSon,    Me%ExtVarSon) 

        call GetGridData(Me%ObjGridDataSon,    Me%ExtVarSon%Batim,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimCoef - ERR20'

        call AllocateVariables

        call FillingCells

        call WriteSmoothCoef
        
        call Read_UnLock_External_Var(Me%ObjGridSon,    Me%ExtVarSon) 

        call UnGetGridData(Me%ObjGridDataSon,    Me%ExtVarSon%Batim,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunProject - SmoothBatimCoef - ERR50'

    end subroutine RunProject

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine WriteSmoothCoef
        
        !Local--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len=StringLength)         :: Coment1, Coment2 
        
        !Begin-----------------------------------------------------------------

        write(*,*)"Writing coefficients..."

        Coment1 = 'File generated by'
        Coment2 = 'Mohid Smooth bathymetry for nesting'


        call WriteGridData(FileName         = trim(Me%CoefFile),         &
                           COMENT1          = Coment1,                          &
                           COMENT2          = Coment2,                          &
                           HorizontalGridID = Me%ObjGridSon,                    &
                           FillValue        = Me%LandPoint,                     &
                           Overwrite        = .true.,                           &
                           GridData2D_Real  = Me%SmoothCoef,                      &
                           STAT             = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_) stop 'WriteSmoothBathym - SmoothBatimCoef - ERR10'

    end subroutine WriteSmoothCoef
    
    
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
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR00'

        call GetGridAngle     (ObjGrid, ExtVar%Rotation, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR01'

        call GetCheckDistortion(ObjGrid, ExtVar%GridDistortion, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR02'

        call GetHorizontalGridSize(ObjGrid, WorkSize = ExtVar%WorkSize, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR10'

        call GetHorizontalGridSize(ObjGrid, Size = ExtVar%Size, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR15'


        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,              &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD)

        !Gets Coordinates in use
        call GetGridCoordType(ObjGrid, CoordType, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR11'
        

        if    (CoordType == SIMPLE_GEOG)then
            
            call GetGridLatitudeLongitude(ObjGrid,                           &
                                          GridLatitudeConn  = ExtVar%YY_IE,  &
                                          GridLongitudeConn = ExtVar%XX_IE,  &
                                          STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR12'

        elseif(CoordType == UTM .or. CoordType == MIL_PORT .or. CoordType == GRID_COORD)then

            call GetHorizontalGrid(ObjGrid,                                  &
                                   XX_IE = ExtVar%XX_IE,                     &
                                   YY_IE = ExtVar%YY_IE,                     &
                                   STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR13'

        else

            write(*,*)'GEOG coordinate type cannot be used in digital terrain generation'
            stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR13'

        end if

        call GetHorizontalGrid(ObjGrid, XX = ExtVar%XX, YY = ExtVar%YY,  STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR40'

        call GetDefineCellsMap(ObjGrid, ExtVar%DefineCellsMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR50'

   
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
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR10'

        call UngetHorizontalGrid(ObjGrid, ExtVar%YY_IE, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR20'

        call UngetHorizontalGrid(ObjGrid, ExtVar%XX, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR30'

        call UngetHorizontalGrid(ObjGrid, ExtVar%YY, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR40'

        call UngetHorizontalGrid(ObjGrid, ExtVar%DefineCellsMap, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'Read_Lock_External_Var - SmoothBatimCoef - ERR50'

    end subroutine Read_UnLock_External_Var

    
    !--------------------------------------------------------------------------

    subroutine FillingCells

        !Local-----------------------------------------------------------------
        integer                             :: d, dd
        integer                             :: i, j
        integer                             :: ILB, IUB
        integer                             :: JLB, JUB
        integer                             :: DLI, DUI
        integer                             :: DLJ, DUJ
        integer                             :: DDLI, DDUI
        integer                             :: DDLJ, DDUJ
        integer                             :: WILB, WIUB
        integer                             :: WJLB, WJUB

        !Begin-----------------------------------------------------------------
        ILB = Me%ExtVarSon%WorkSize%ILB
        IUB = Me%ExtVarSon%WorkSize%IUB
        JLB = Me%ExtVarSon%WorkSize%JLB
        JUB = Me%ExtVarSon%WorkSize%JUB
        WILB = Me%ExtVarSon%WindowSize%ILB
        WIUB = Me%ExtVarSon%WindowSize%IUB
        WJLB = Me%ExtVarSon%WindowSize%JLB
        WJUB = Me%ExtVarSon%WindowSize%JUB
        
        if (WILB == 0 .or. WILB < ILB ) WILB = ILB
        if (WIUB == 0 .or. WIUB > IUB ) WIUB = IUB
        if (WJLB == 0 .or. WJLB < JLB ) WJLB = JLB
        if (WJUB == 0 .or. WJUB > JUB ) WJUB = JUB

        d = Me%Ncells;      
        dd = d*d;
        
        write(*,*)"Filling cells coefs values..."
        
        do i = ILB,  IUB
        do j = JLB , JUB
                       
idef:       if (Me%ExtVarSon%DefineCellsMap(i, j)==1 .and. Me%ExtVarSon%Batim(i, j) /= Me%LandPoint) then
                
                !Put main code here ...
                !In a picture frame there are 5 areas:
                !DLI area
                !DLJ area
                !DUI area
                !DUJ area
                !Inner and outer area
                !          DUJ
                !    ________________
                !   |\ ____________ /|
                !   | |            | |
                !   | |            | |
                !   | |            | |
                !DLI| |            | | DUI
                !   | |    Inner   | |
                !   | |     area   | |
                !   | |            | |
                !   | |____________| |
                !   |/______________\|
                !
                !          DLJ
                !
                DLI = WILB-i
                DDLI = DLI*DLI
                DUI = WIUB-i
                DDUI = DUI*DUI
                DLJ = WJLB-J
                DDLJ = DLJ*DLJ
                DUJ = WJUB-j
                DDUJ = DUJ*DUJ

cd3:            if (j > WJLB - d .and. j < WJUB + d                     &
                    .and. i > WILB - d .and. i < WIUB + d ) then
                
    cd2:            if (j > WJLB .and. j < WJUB                     &
                        .and. i > WILB .and. i < WIUB ) then

        cd1:            if( DDLI<dd .and. DDLI<DDLJ .and. DDLI<DDUJ) then
                    
                            Me%SmoothCoef(i,j) = real(dd-DDLI)/real(dd)
                               
                        elseif( DDLJ<dd .and. DDLJ<DDUI ) then

                            Me%SmoothCoef(i,j) = real(dd-DDLJ)/real(dd)
                
                        elseif( DDUI<dd .and. DDUI<DDUJ ) then
                
                            Me%SmoothCoef(i,j) = real(dd-DDUI)/real(dd)

                        elseif( DDUJ<dd ) then

                            Me%SmoothCoef(i,j) = real(dd-DDUJ)/real(dd)
                
                        else
                
                            Me%SmoothCoef(i,j) = 0.

                        endif cd1

                    else

        cd4:            if( DDLI < dd .and. i .le. WILB                 &
                            .and. .not.( j < WJLB .and. DDLI < DDLJ .or. &
                            j > WJUB .and. DDLI < DDUJ ) ) then
                    
                            Me%SmoothCoef(i,j) = real(dd-DDLI)/real(dd)
                               
                        elseif( DDLJ < dd .and. j .le. WJLB                      &
                                .and. .not.( i < WILB .and. DDLJ < DDLI .or.     &
                                i > WIUB .and. DDLJ < DDUI ) ) then

                            Me%SmoothCoef(i,j) = real(dd-DDLJ)/real(dd)
                
                        elseif( DDUI < dd .and. i .ge. WIUB                     &
                                .and. .not.( j < WJLB .and. DDUI < DDLJ .or.    &
                                j > WJUB .and. DDUI < DDUJ ) ) then
                
                            Me%SmoothCoef(i,j) = real(dd-DDUI)/real(dd)

                        elseif( DDUJ < dd ) then

                            Me%SmoothCoef(i,j) = real(dd-DDUJ)/real(dd)
                
                        else
                
                            Me%SmoothCoef(i,j) = 0.

                        endif cd4

                    endif cd2

                else

                    Me%SmoothCoef(i,j) = 0.

                endif cd3

            endif idef     

        end do
        end do

        !Uncomment this loop if the coefs are the other way around
        do i = ILB,  IUB
        do j = JLB , JUB
            Me%SmoothCoef(i,j) = 1. - Me%SmoothCoef(i,j)
        end do
        end do

    end subroutine FillingCells


    !--------------------------------------------------------------------------

    
    subroutine CloseProject

        !Local-----------------------------------------------------------------
        real                                        :: ElapsedSeconds
        real                                        :: TotalCPUTime
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call KillGridData(Me%ObjGridDataSon, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimCoef - ERR30'

        call KillHorizontalGrid(HorizontalGridID= Me%ObjGridSon, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'CloseProject - SmoothBatimCoef - ERR10'

        deallocate(Me%SmoothCoef)

        call date_and_time(Values = Me%F95Time)
        call SetDate      (Me%FinalSystemTime,float(Me%F95Time(1)), float(Me%F95Time(2)),      &
                                              float(Me%F95Time(3)), float(Me%F95Time(5)),      &
                                              float(Me%F95Time(6)), float(Me%F95Time(7))+      &
                                              Me%F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = Me%FinalSystemTime - Me%InitialSystemTime

        call ShutdownMohid ("SmoothBathimCoefs", ElapsedSeconds, TotalCPUTime)


    end subroutine CloseProject

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

end program SmoothBatimCoef


