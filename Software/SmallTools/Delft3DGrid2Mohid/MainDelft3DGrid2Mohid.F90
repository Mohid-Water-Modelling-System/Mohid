!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Delft3DGrid2Mohid
! PROGRAM       : MainDelft3DGrid2Mohid
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig /Luis Fernandes - v4.0
! DESCRIPTION   : Delft3DGrid2Mohid to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MohidDelft3DGrid2Mohid

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleHorizontalGrid

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: STAT_CALL
    character(len=PathLength)   :: Input_Delft3D_Grid, Output_Mohid_Grid
    real                        :: Delft3D_Grid_FillValue
    
    integer, parameter          :: FileOpen = 1, FileClose = 0    


    call ConstructMohidDelft3DGrid2Mohid
    call ModifyMohidDelft3DGrid2Mohid
    call KillMohidDelft3DGrid2Mohid

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidDelft3DGrid2Mohid
        
        call StartUpMohid("MohidDelft3DGrid2Mohid")

        call StartCPUTime

        call ReadKeywords


    end subroutine ConstructMohidDelft3DGrid2Mohid
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidDelft3DGrid2Mohid
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        call ReadWriteGrids
    
    end subroutine ModifyMohidDelft3DGrid2Mohid
    
    !--------------------------------------------------------------------------

    subroutine KillMohidDelft3DGrid2Mohid

        call StopCPUTime

        call ShutdownMohid ("MohidDelft3DGrid2Mohid", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidDelft3DGrid2Mohid
    
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
    
    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: FromFile, iflag
        
        !Begin-----------------------------------------------------------------        

        call ReadFileName('IN_MODEL', DataFile, "MohidDelft3DGrid2Mohid", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDelft3DGrid2Mohid - ERR10'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDelft3DGrid2Mohid - ERR20'

        call GetExtractType     (FromFile = FromFile)
        
        call GetData(Input_Delft3D_Grid,                                                &
                     ObjEnterData, iflag,                                               &
                     keyword    = 'INPUT_DELFT3D_GRID',                                 &
                     SearchType = FromFile,                                             &
                     ClientModule ='MohidDelft3DGrid2Mohid',                            &
                     STAT       = STAT_CALL)

        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
            call SetError(FATAL_, INTERNAL_, 'ReadKeywords - MohidDelft3DGrid2Mohid - ERR30')
        endif
        
        call GetData(Delft3D_Grid_FillValue,                                            &
                     ObjEnterData, iflag,                                               &
                     keyword    = 'DELFT3D_GRID_FILLVALUE',                             &
                     SearchType = FromFile,                                             &
                     ClientModule ='MohidDelft3DGrid2Mohid',                            &
                     STAT       = STAT_CALL)

        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
            call SetError(FATAL_, INTERNAL_, 'ReadKeywords - MohidDelft3DGrid2Mohid - ERR40')
        endif        
        
        call GetData(Output_Mohid_Grid,                                                 &
                     ObjEnterData, iflag,                                               &
                     keyword    = 'OUTPUT_MOHID_GRID',                                  &
                     SearchType = FromFile,                                             &
                     ClientModule ='MohidDelft3DGrid2Mohid',                            &
                     STAT       = STAT_CALL)

        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
            call SetError(FATAL_, INTERNAL_, 'ReadKeywords - MohidDelft3DGrid2Mohid - ERR50')
        endif        
            

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDelft3DGrid2Mohid - ERR60'

    end subroutine ReadKeywords
    
    !--------------------------------------------------------------------------
    
    subroutine ReadWriteGrids

        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer     :: XX, YY
        integer                             :: Imax, Jmax, InputID, OutputID, i, j
        integer                             :: STAT_CALL, DummyInt
        character(13)                       :: DummyChar
        
        !Begin-----------------------------------------------------------------
        
        call UnitsManager(InputID, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SetError(FATAL_, INTERNAL_, 'ReadWriteGrids - MohidDelft3DGrid2Mohid - ERR10')
        endif

        open(unit = InputID, file = Input_Delft3D_Grid, form = 'FORMATTED', status = 'UNKNOWN')
        
        do i =1, 5
            read(InputID,*)
        enddo
        
        read(InputID,*)  jmax, imax
        read(InputID,*)
        
        allocate(XX(0:imax,0:jmax),YY(0:imax,0:jmax))
        
        do i=1,imax
            read(InputID,*) DummyChar,DummyInt,(XX(i,j),j=1,jmax)
        enddo
        
        do i=1,imax
            read(InputID,*) DummyChar,DummyInt,(YY(i,j),j=1,jmax)
        enddo        
        
        do i=1,imax
        do j=1,jmax
            if (XX(i,j) == Delft3D_Grid_FillValue) then
                XX(i,j) = FillValueReal
                YY(i,j) = FillValueReal
            endif
        enddo
        enddo
        
        call UnitsManager(InputID, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SetError(FATAL_, INTERNAL_, 'ReadWriteGrids - MohidDelft3DGrid2Mohid - ERR20')
        endif        
        
        
        call UnitsManager(OutputID, FileOpen, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SetError(FATAL_, INTERNAL_, 'ReadWriteGrids - MohidDelft3DGrid2Mohid - ERR30')
        endif

        open(unit = OutputID, file = Output_Mohid_Grid, form = 'FORMATTED', status = 'UNKNOWN')
        

        write(OutputID,*) "COMENT1    : From Delft3D Grid"
        write(OutputID,*) "COMENT2    : to MOHID Grid"
        write(OutputID,*) "ILB_IUB    : 1", Imax-1
        write(OutputID,*) "JLB_JUB    : 1", Jmax-1
        write(OutputID,*) "COORD_TIP  :  5"
        write(OutputID,*) "ORIGIN     :  0 0"
        write(OutputID,*) "GRID_ANGLE :  0"
        write(OutputID,*) "LATITUDE   :  42"
        write(OutputID,*) "LONGITUDE  :  -9"
        
        write(OutputID,*) "FILL_VALUE :  -99.0"
 
 
        write(OutputID,*) "<CornersXY>"
 
        do i=1,imax
        do j=1,jmax
            write(OutputID,*) XX(i,j), YY(i,j)
        enddo
        enddo
 
 
        write(OutputID,*) "<\CornersXY>"        
        
        call UnitsManager(OutputID, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            call SetError(FATAL_, INTERNAL_, 'ReadWriteGrids - MohidDelft3DGrid2Mohid - ERR50')
        endif               
        
        
    end subroutine ReadWriteGrids
    

end program MohidDelft3DGrid2Mohid
