!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Basin Delineator
! PROGRAM       : Basin Delineator
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Top Level of the Module Basin Geometry
!
!------------------------------------------------------------------------------
program BasinDelineator

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleBasinGeometry

    implicit none

    integer                         :: ObjGridData        = 0
    integer                         :: ObjHorizontalGrid    = 0
    integer                         :: ObjBasinGeomentry    = 0
    integer                         :: ObjEnterData         = 0
    integer                         :: STAT_CALL, flag
    character(PathLength)           :: ProjectFile          = "Basin.dat"
    character(PathLength)           :: TopographicFile

    !Other Stuff
    type (T_Time)                   :: InitialSystemTime, FinalSystemTime
    integer, dimension(8)           :: F95Time

    call ConstructBasinDelineator
    call KillBasinDelineator

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructBasinDelineator

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        call StartupMohid ("Basin Delineator")

        !Gets the actual time
        call date_and_time(Values = F95Time)
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

        call ConstructEnterData (ObjEnterData, ProjectFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasinDelineator - BasinDelineator - ERR01'

        !Get grid file path
        call GetData(TopographicFile,                                                        &
                     ObjEnterData, flag,                                                     &
                     SearchType   = FromFile,                                                &
                     keyword      ='TOPOGRAPHIC_FILE',                                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_ .or. flag == 0) stop 'ConstructBasinDelineator - BasinDelineator - ERR02'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasinDelineator - BasinDelineator - ERR07'

        call ConstructHorizontalGrid (ObjHorizontalGrid, TopographicFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasinDelineator - BasinDelineator - ERR03'

        call ConstructGridData      (ObjGridData, ObjHorizontalGrid, FileName = TopographicFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasinDelineator - BasinDelineator - ERR04'

        call ConstructBasinGeometry (ObjBasinGeomentry, ObjGridData, ObjHorizontalGrid,    &
                                     ProjectFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBasinDelineator - BasinDelineator - ERR05'


    end subroutine ConstructBasinDelineator

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillBasinDelineator


        !Local-----------------------------------------------------------------
        real                                        :: ElapsedSeconds, TotalCPUTime

        call KillBasinGeometry (ObjBasinGeomentry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BasinDelineator - BasinDelineator - ERR06'

        call KillGridData     (ObjGridData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BasinDelineator - BasinDelineator - ERR09'

        call KillHorizontalGrid (ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'BasinDelineator - BasinDelineator - ERR08'

        call date_and_time(Values = F95Time)
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)
        ElapsedSeconds = FinalSystemTime - InitialSystemTime

        call ShutdownMohid ("Basin Delineator", ElapsedSeconds, TotalCPUTime)

    end subroutine KillBasinDelineator

end program
