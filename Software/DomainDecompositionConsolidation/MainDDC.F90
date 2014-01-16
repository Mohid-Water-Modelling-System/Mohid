!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model - Results consolidation
! PROJECT       : DDC - Domain Decomposition Consolidation
! PROGRAM       : MainDDC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          :  2013
! REVISION      : Ricardo Miranda
! DESCRIPTION   : Program to consolidate results from MPI run with domain decomposition
!
!------------------------------------------------------------------------------

program MohidDDC

    use ModuleDDC
    use ModuleGlobalData

    implicit none

    type       T_MohidDDC
        real                        :: start            = NULL_REAL
        real                        :: finish           = NULL_REAL
        real                        :: TotalCPUTime     = NULL_REAL
        real                        :: ElapsedSeconds   = NULL_REAL
        real(8)                     :: t1               = NULL_REAL
    
        type (T_DDC), pointer       :: ObjDDC
    end type  T_MohidDDC
    type(T_MohidDDC), pointer       :: Me

    call ConstructMohidDDC
    call ModifyMohidDDC
    call KillMohidDDC

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidDDC
    
        !----------------------------------------------------------------------

        call AllocateInstance

        call StartCPUTime

        Me%ObjDDC => ConstructDDC()

        !----------------------------------------------------------------------

    end subroutine ConstructMohidDDC
    
    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type(T_MohidDDC), pointer                           :: NewMohidDDC
        
        !----------------------------------------------------------------------

        allocate(NewMohidDDC)
        nullify(NewMohidDDC%ObjDDC)
        Me => NewMohidDDC

        !----------------------------------------------------------------------

    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidDDC
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL        = NULL_INT

        !----------------------------------------------------------------------

        STAT_CALL = ModifyDDC(Me%ObjDDC) 
        if (STAT_CALL /= SUCCESS_) stop 'MohidDDC - ModifyMohidDDC - ERR01'

        !----------------------------------------------------------------------

    end subroutine ModifyMohidDDC
    
    !--------------------------------------------------------------------------

    subroutine KillMohidDDC

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL        = NULL_INT

        !----------------------------------------------------------------------

        STAT_CALL = KillDDC(Me%ObjDDC)
        if (STAT_CALL /= SUCCESS_) stop 'MohidDDC - KillMohidDDC - ERR01'
        
        call StopCPUTime  
        call ShutdownDDC(ModelName      = "MohidDDC",                           &
                         ElapsedSeconds = Me%ElapsedSeconds,           &
                         TotalCPUTime   = Me%TotalCPUTime)
        call DeallocateInstance

        !----------------------------------------------------------------------

    end subroutine KillMohidDDC
    
    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !----------------------------------------------------------------------

        deallocate(Me)

        !----------------------------------------------------------------------

    end subroutine DeallocateInstance
        
    !--------------------------------------------------------------------------

    subroutine ShutdownDDC(ModelName, ElapsedSeconds, TotalCPUTime)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(IN)    :: ModelName
        real, intent(IN)                :: ElapsedSeconds
        real, intent(IN)                :: TotalCPUTime
        integer                         :: ElapsedHours, ElapsedMinutes, ElapsedSecremain

        !----------------------------------------------------------------------

        ElapsedHours        = INT(Elapsedseconds/3600)
        ElapsedMinutes      = INT((ElapsedSeconds-ElapsedHours*3600)/60)
        ElapsedSecremain    = INT((ElapsedSeconds-ElapsedMinutes*60-ElapsedHours*3600))

        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Program "//ModelName//" successfully terminated"
        write(*, *)                    
        write(*, *)
        write(*, 110)ElapsedSeconds, ElapsedHours, ElapsedMinutes,ElapsedSecremain
        write(*, 120)TotalCPUTime
        write(*, 130)100.*TotalCPUTime/ElapsedSeconds
        write(*, *)
        write(*, *)"----------------------------------------------------------"


    110 format(1x, "Total Elapsed Time     : ",f14.2," ",i3,"h ",i2,"min ",i2,"s",/)
    120 format(1x, "Total CPU time         : ",f14.2,/)
    130 format(1x, "CPU usage (%)          : ",f14.2,/)

        !----------------------------------------------------------------------

    end subroutine ShutdownDDC

    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        !Local-----------------------------------------------------------------
        real                            :: secnds
        
        !----------------------------------------------------------------------

        Me%t1 = secnds(0.0)
        call cpu_time(Me%start)
        
        !----------------------------------------------------------------------

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        !Local-----------------------------------------------------------------
        real                            :: secnds

        !----------------------------------------------------------------------

        Me%ElapsedSeconds = secnds(Me%t1)
        call cpu_time(Me%finish)
        Me%TotalCPUTime = Me%finish - Me%start

        !----------------------------------------------------------------------

    end subroutine StopCPUTime
    
    !--------------------------------------------------------------------------
    
end program MohidDDC
