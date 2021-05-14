!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : MeshGlue
! PROGRAM       : MainMeshGlue
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod
! DATE          : December 2020
! REVISION      : Paulo Leitao - v1.0
! DESCRIPTION   : MeshGlue to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MohidMeshGlue

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions
    use ModuleTime
    use ModuleMeshGlue

    implicit none

    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: STAT_CALL
    integer                     :: MeshGlueID = 0


    call ConstructMohidMeshGlue
    call ModifyMohidMeshGlue
    call KillMohidMeshGlue

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidMeshGlue
        
        call StartUpMohid("MohidMeshGlue")

        call StartCPUTime

        !call ReadKeywords
        
        call ConstructMeshGlue(ObjMeshGlueID = MeshGlueID, STAT = STAT_CALL)   


    end subroutine ConstructMohidMeshGlue
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidMeshGlue
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running


        call ModifyMeshGlue(ObjMeshGlueID = MeshGlueID, STAT = STAT_CALL)              
    
    
    
    end subroutine ModifyMohidMeshGlue
    
    !--------------------------------------------------------------------------

    subroutine KillMohidMeshGlue
    
        call KillMeshGlue(ObjMeshGlueID = MeshGlueID, STAT = STAT_CALL)   

        call StopCPUTime

        call ShutdownMohid ("MohidMeshGlue", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidMeshGlue
    
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
        integer                                     :: FromFile

        call ReadFileName('IN_MODEL', DataFile, "MohidMeshGlue", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidMeshGlue - ERR01'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidMeshGlue - ERR02'

        call GetExtractType     (FromFile = FromFile)

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidMeshGlue - ERR03'

    end subroutine ReadKeywords

end program MohidMeshGlue
