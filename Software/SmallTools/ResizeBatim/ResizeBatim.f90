!------------------------------------------------------------------------------
!        INTECMAR/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Resize Batim
! PROJECT       : MohidResizeBatim
! MODULE        : -
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : December 2009
! REVISION      : Pedro Montero - v4.0
! DESCRIPTION   : Program to resize grid data batim filling with land rows or columns
!------------------------------------------------------------------------------

!Data file - default name 'Resize.dat' (must be placed in working directory)

!   IN_BATIM                    : char              -           !Original bathymetry file name
!   OUT_BATIM                   : char              -           !Filtered bathymetry file name
!   ADD_ROWS                    : int               -           !Number of rows to add
!   ADD_COLUMNS                 : int               -           !Number of columns to add


program MohidResizeBatim

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData     
    use ModuleHorizontalGrid                                     
    use ModuleGridData  
    use ModuleHorizontalMap
    use ModuleFunctions
             
    implicit none

    !Instances
    integer                                     :: ObjEnterData         = 0
    integer                                     :: ObjBathymetry        = 0
    integer                                     :: ObjHorizontalGrid    = 0
    integer                                     :: ObjHorizontalMap     = 0

    !Dummy variable time
    type (T_Time  )                             :: CurrentTime

    !Size
    type (T_Size2D)                             :: Size
    type (T_Size2D)                             :: WorkSize

    !Files
    character(len=PathLength)                   :: BathymetryFile
    character(len=PathLength)                   :: NewBathymetryFile
    character(len=PathLength)                   :: DataFile             = 'Resize.dat'

    !Input variables
    integer                                     :: AddRows
    integer                                     :: AddColumns

    !Working variables
    integer, dimension(:,:), pointer            :: WaterPoints2D
    real,    dimension(:,:), pointer            :: Bathymetry
   

    
    !Begin---------------------------------------------------------------------

    call StartUpMohid("Mohid Bathymetry Resizer")

    call ReadOptions

    call ConstructDomain

    call WriteNewBathymetry

    call KillBathymetryResizer

    contains

    !--------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Reading options...'
        write(*,*)

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidResizeBatim - ERR10'

        call GetData(BathymetryFile,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'IN_BATIM',                         &
                     ClientModule = 'MohidResizeBatim',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidResizeBatim - ERR20'

        call GetData(NewBathymetryFile,                                 &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'OUT_BATIM',                        &
                     ClientModule = 'MohidResizeBatim',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidResizeBatim - ERR30'

        call GetData(AddRows,                                           &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'ADD_ROWS',                         &
                     Default      = 0,                                  &
                     ClientModule = 'MohidResizeBatim',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidResizeBatim - ERR40'

        call GetData(AddColumns,                                        &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'ADD_COLUMNS',                      &
                     Default      = 0,                                  &
                     ClientModule = 'MohidResizeBatim',                 &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidResizeBatim - ERR50'

        if(AddRows < 0 .or. AddColumns < 0)then
            write(*,*)'Rows or Columns to add must be defined 0 or positive'
            stop 'ReadOptions - MohidResizeBatim - ERR60'
        end if




        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - MohidResizeBatim - ERR60'

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ConstructDomain

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
    
        write(*,*)
        write(*,*)"Constructing horizontal grid..."
        write(*,*)

        call SetDate(CurrentTime, 2004, 1, 1, 0, 0, 0)

        !Horizontal Grid
        call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGrid,          &
                                     DataFile         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidResizeBatim - ERR01'


        write(*,*)
        write(*,*)"Constructing bathymetry..."
        write(*,*)


        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     FileName         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidResizeBatim - ERR02'


        write(*,*)
        write(*,*)"Constructing mapping..."
        write(*,*)

        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     ActualTime       = CurrentTime,                &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidResizeBatim - ERR03'

        write(*,*)
        write(*,*)"Compiling information..."
        write(*,*)

        call GetWaterPoints2D(HorizontalMapID  = ObjHorizontalMap,                  &
                              WaterPoints2D    = WaterPoints2D,                     &
                              STAT             = STAT_CALL)  

        call GetHorizontalGridSize(ObjHorizontalGrid, Size, WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidResizeBatim - ERR04'

        call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidResizeBatim - ERR05'

        write(*,*)
        write(*,*)"Domain successfully constructed..."
        write(*,*)


    end subroutine ConstructDomain

    !--------------------------------------------------------------------------  

    subroutine WriteNewBathymetry

        !Local-----------------------------------------------------------------
        integer                                :: STAT_CALL
        real,    dimension(:  ), pointer    :: XX, YY, XX_Aux, YY_Aux
        type (T_Size2D)                     :: NewWorkSize,NewSize  
        integer                             :: GRID_COORD, CoordType, Zone
        real                                :: Xorig, Yorig, GRID_ANGLE, Latitude, Longitude
        character(len=StringLength)         :: Coment1, Coment2 
        integer                             :: i,j
        real,    dimension(:,:), pointer   :: NewBathymetry

        !Begin-----------------------------------------------------------------

        
        
        
        
        
        write(*,*)
        write(*,*)"Writing new bathymetry..."
        write(*,*)



        ! Resize grid limits
        NewWorkSize%ILB = WorkSize%ILB 
        NewWorkSize%JLB = WorkSize%JLB 
        NewWorkSize%IUB = WorkSize%IUB+AddRows
        NewWorkSize%JUB = WorkSize%JUB+AddColumns
        NewSize%ILB = Size%ILB 
        NewSize%JLB = Size%JLB 
        NewSize%IUB = Size%IUB+AddRows
        NewSize%JUB = Size%JUB+AddColumns

        Coment1 = 'File generated by Mohid Bathymetry Resizer'
        Coment2 = 'Original Bathymetry '//trim(BathymetryFile)

        !Gets Coordinates in use
        call GetGridCoordType(ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR10'

        !Gets Origin in use
        Call GetGridOrigin(ObjHorizontalGrid, Xorig, Yorig, STAT= STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'WriteNewBathymetry - MohidResizeBatim - ERR20'

        !Gets Grid in use
        call GetHorizontalGrid(ObjHorizontalGrid,                                       &
                                XX = XX,                                                &
                                YY = YY,                                                &
                                STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR30'
       
        !Gets Angle in use
        call GetGridAngle(ObjHorizontalGrid, GRID_ANGLE, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR40'
        
        !Gets Zone in use
        call GetGridZone(ObjHorizontalGrid, Zone, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR50'
        
        !Gets LatitudeLongitude in use
        call GetLatitudeLongitude(ObjHorizontalGrid, Latitude, Longitude, STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR60'


        nullify (NewBathymetry)
        allocate(NewBathymetry(NewWorkSize%ILB:NewWorkSize%IUB, NewWorkSize%JLB:NewWorkSize%JUB))

        NewBathymetry = -99
        
        

        do i=Worksize%ILB,Worksize%IUB
        do j=Worksize%JLB,Worksize%JUB
            NewBathymetry(i,j)=Bathymetry(i,j)
        enddo
        enddo

        

       
        allocate(XX_aux(NewWorkSize%JLB:NewWorkSize%JUB+1))
        allocate(YY_aux(NewWorkSize%ILB:NewWorkSize%IUB+1))

        do i=NewWorkSize%ILB,NewWorkSize%IUB+1

            YY_aux(i)=(YY(Worksize%ILB+1)-YY(Worksize%ILB))*(i-1)
        enddo
        do j=NewWorkSize%JLB,NewWorkSize%JUB+1
            XX_aux(j)=(XX(Worksize%JLB+1)-XX(Worksize%JLB))*(j-1)
        enddo

            

        call WriteGridData(FileName             = trim(NewBathymetryFile),         &
                          XX                   = XX_aux,                           &
                          YY                   = YY_aux,                           &
                          COMENT1              = COMENT1,                          &
                          COMENT2              = COMENT2,                          &
                          WorkSize             = NewWorkSize,                      &  
                          CoordType            = CoordType,                        &
                          Xorig                = Xorig,                            &
                          Yorig                = Yorig,                            &
                          Zone                 = Zone,                             &
                          GRID_ANGLE           = GRID_ANGLE,                       &
                          Latitude             = Latitude,                         &
                          Longitude            = Longitude,                        &
                          FillValue            = -99.,                             &
                          Overwrite            = .True.,                           &
                          GridData2D_Real      = NewBathymetry,                    &
                          STAT                 = STAT_CALL) 

           write(*,*)
           write(*,*)"Deallocating variables..."
           write(*,*)

           deallocate(NewBathymetry); 
           nullify (NewBathymetry)
           deallocate(XX_aux)
           deallocate(YY_aux)



        

            call UnGetHorizontalGrid(ObjHorizontalGrid, XX, STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR70'

            call UnGetHorizontalGrid(ObjHorizontalGrid, YY, STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_)stop 'WriteNewBathymetry - MohidResizeBatim - ERR80'
       





        
        write(*,*)
        write(*,*)"Done..."
        write(*,*)




    end subroutine WriteNewBathymetry

    !--------------------------------------------------------------------------

    

    subroutine KillBathymetryResizer

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)"Finishing freeing memory..."
        write(*,*)

        call KillHorizontalMap  (ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillBathymetryFilter - MohidResizeBatim - ERR01'

        call KillGridData       (ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillBathymetryFilter - MohidResizeBatim - ERR02'

        call KillHorizontalGrid (ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillBathymetryFilter - MohidResizeBatim - ERR03'

    end subroutine KillBathymetryResizer


end program MohidResizeBatim
