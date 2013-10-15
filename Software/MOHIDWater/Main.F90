!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Water
! PROGRAM       : MohidWater
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Main Program which run Mohid Water
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------
!            
! Modelling system Mohid a general description
!              MOHID is a full 3D-baroclinic model and has been developed 
!              using an object oriented programming philosophy and using all 
!              the FORTRAN 95 potential. The system has two main classes: the 
!              first one manages the hydrodynamic properties (e.g. velocity, 
!              elevation, water fluxes, turbulent viscosity) and the second 
!              one the water properties (e.g. salinity, temperature, density, 
!              SPM, nutrients, phytoplankton, coliforms). The model is based 
!              on a finite volume concept. In this approach the discrete form 
!              of the governing equations are applied macroscopically to the 
!              cell control volume in the form of flux divergence. As a consequence 
!              this method automatically guarantees the conservation of transported 
!              properties .
!              The model computes the evolution of several flow and water properties
!              in a finite-volume.
! For further information please take a look at:
!              http://www.mohid.com
! or contact
!              e-mail: general@mohid.com
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!------------------------------------------------------------------------------
#ifdef _OPENMI_
module MohidWater
#else
program MohidWater
#endif

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData,        only : ReadFileName 
    use ModuleStopWatch,        only : CreateWatchGroup, KillWatchGroup, &
                                        StartWatch, StopWatch
    use ModuleHorizontalGrid,   only : ConstructHorizontalGrid,                          &
                                       ConstructFatherGridLocation,                      &
                                       GetGridFileName, GetHorizontalGridSize,           &
                                       GetNotDefinedCells
#ifdef OVERLAP
    use ModuleModel,            only : ConstructModel, UpdateTimeAndMapping,             &
                                       RunModel, KillModel, GetModelTimeStep,            &
                                       GetModelTimeLimits,  GetModelInstanceIDs,         &
                                       GetModelCurrentTime, GetSubModelWindow,           &
                                       GetSubModelWindowON, GetModelOverlap,             &
                                       GetModelOverlapInfo
#else  OVERLAP
    use ModuleModel,            only : ConstructModel, UpdateTimeAndMapping,             &
                                       RunModel, KillModel, GetModelTimeStep,            &
                                       GetModelTimeLimits,  GetModelInstanceIDs,         &
                                       GetModelCurrentTime, GetSubModelWindow,           &
                                       GetSubModelWindowON
#endif OVERLAP


#ifdef _USE_MPI
    use ModuleHydrodynamic,     only : GetHydroNeedsFather, SetHydroFather,              &
                                       SendHydrodynamicMPI, RecvHydrodynamicMPI,         &
                                       UpdateHydroMPI
    use ModuleWaterproperties,  only : GetWaterNeedsFather, GetPropListNeedsFather,      &
                                       SetWaterPropFather, SendWaterPropertiesMPI,       &
                                       RecvWaterPropertiesMPI, UpdateWaterMPI
    use ModuleFunctions,        only : MPIKind
    use mpi
#else _USE_MPI
                                       
#ifdef OVERLAP
    use ModuleHydrodynamic,     only : GetHydroNeedsFather, SetHydroFather,              &
                                       GetHydroOverlap, SetModelOverlapHydro
    use ModuleWaterproperties,  only : GetWaterNeedsFather, GetPropListNeedsFather,      &
                                       SetWaterPropFather,  GetWaterOverlap,             &
                                       SetModelOverlapWater
#else  OVERLAP
    use ModuleHydrodynamic,     only : GetHydroNeedsFather, SetHydroFather
    use ModuleWaterproperties,  only : GetWaterNeedsFather, GetPropListNeedsFather,      &
                                       SetWaterPropFather
#endif  OVERLAP

#endif _USE_MPI
    

    implicit none


    type T_ModelLink
        logical                                             :: Hydro   = .false.
        logical                                             :: Water   = .false.
        logical                                             :: Nesting = .false.
        integer                                             :: nProp
        integer, dimension(:), pointer                      :: PropertyIDNumbers
        type (T_Size2D)                                     :: Window
    end type 


    type T_MohidWater
        integer                                             :: ModelID
        integer                                             :: ModelLevel
        character(StringLength)                             :: ModelName
        character(PathLength)                               :: ModelPath
        integer                                             :: HorizontalGridID
        integer                                             :: HydrodynamicID
        integer                                             :: WaterpropertiesID
        type (T_Time)                                       :: CurrentTime
        type (T_Time)                                       :: InfoTime         !Time of the last father information
        logical                                             :: SubOn
        integer                                             :: nSubModels
        type (T_ModelLink)                                  :: FatherLink
        type (T_ModelLink), dimension(:), pointer           :: SubmodelLink
        integer                                             :: MPI_ID
        integer, dimension(:), pointer                      :: SubMPIID

#ifdef OVERLAP
        logical                                             :: Overlap
        integer                                             :: OverlapModelID
        integer, dimension(:,:), pointer                    :: OverlapCells
        type (T_ModelLink)                                  :: OverlapLink
        integer                                             :: OverlapHydrodynamicID
        integer                                             :: OverlapWaterPropertiesID
#endif OVERLAP
        
        integer                                             :: FatherGridID  = 0    !MPI only
        type (T_MohidWater), pointer                        :: FatherModel
        type (T_MohidWater), pointer                        :: Next    
    end type T_MohidWater


    !Variables
    type (T_MohidWater), pointer                            :: FirstModel
    type(T_Time)                                            :: GlobalBeginTime
    type(T_Time)                                            :: GlobalEndTime
    type(T_Time)                                            :: GlobalCurrentTime

    integer                                                 :: NumberOfModels       = 0

    !MPI Stuff        
    integer                                                 :: myMPI_ID = null_int
    integer                                                 :: nMPI_Processes

#ifdef _USE_MPI
    character(MPI_MAX_PROCESSOR_NAME)                       :: myMPI_Processor
    integer                                                 :: Precision
#else  _USE_MPI
    character(StringLength)                                 :: myMPI_Processor
#endif _USE_MPI

    logical                                                 :: RunInParallel
                                                            
    !CPU / System Time
    type (T_Time)                                           :: InitialSystemTime, FinalSystemTime
    real                                                    :: TotalCPUTime,     ElapsedSeconds

    real                                                    :: CPUTimeConstructor, CPUTimeModifier
    real                                                    :: WorkCycleElapsed, WorkCycleCPUTime
    type (T_Time)                                           :: WorkCycleInitial, WorkCycleFinal

    integer                                                 :: ObjLagrangianGlobal        = 0
    integer,                      dimension(:,:), pointer   :: LagInstance
    character(len=StringLength),  dimension(:  ), pointer   :: ModelNames 
    real, dimension(:), allocatable                         :: ModelDTs             !MPI only


#ifndef _OPENMI_


#ifdef _USE_MPI
    RunInParallel = .true.
#else  _USE_MPI
    RunInParallel = .false.

    !Disable unused variable warnings
    myMPI_Processor = 'My Processor'
#endif _USE_MPI
    
    if (RunInParallel) then
        call ConstructMohidWaterMPI
    else
        call ConstructMohidWater
    endif

    call ModifyMohidWater

    if (RunInParallel) then
        call KillMohidWaterMPI
    else
        call KillMohidWater
    endif
    
#endif

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidWater

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(PathLength)                       :: WatchFile, DTLogFile
        type (T_MohidWater), pointer                :: CurrentModel
        integer                                     :: iProp
        logical                                     :: SubModelWindowON, NotDefinedCells
        type(T_Time)                                :: SubModelBeginTime, SubModelEndTime
        !Begin-----------------------------------------------------------------

        !Common Startup tasks
        call StartUpMohid   ("Mohid Water")
        call GetSystemTime  (InitialSystemTime)

        !Monitor Performance of the model execution?
        call ReadFileName('OUTWATCH', WatchFile, Message = 'Stop Watch File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            STAT_CALL   = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        else
            MonitorPerformance  = .false.
        endif

        !Monitors DT
        call ReadFileName('DT_LOG', DTLogFile, Message = 'Start DTLog File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then            
            MonitorDT  = .true.
        else
            MonitorDT  = .false.
        endif

        if (MonitorDT) then
            call UnitsManager (UnitDT, OPEN_FILE)      
            open(UNIT   = UnitDT, FILE   = DTLogFile, STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR10'
        end if

        !Constucts the list of the model
        call ConstructModelList

        !Constructs the Models
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            call SetFilesName  (CurrentModel%ModelPath)
            call ConstructModel(LagInstance, ModelNames, NumberOfModels, ObjLagrangianGlobal, &
                                CurrentModel%ModelID, InitialSystemTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR20'

            !Get the Instance IDs of the Objects which are necessary for Model/SubModel
            !comunication
            call GetModelInstanceIDs      (CurrentModel%ModelID,                              &
                                           HorizontalGridID  = CurrentModel%HorizontalGridID, &
                                           HydrodynamicID    = CurrentModel%HydrodynamicID,   &
                                           WaterpropertiesID = CurrentModel%WaterpropertiesID,&
                                           STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR30'

            !Gets Current Time of the model
            call GetModelCurrentTime      (CurrentModel%ModelID, CurrentModel%CurrentTime,    &
                                           STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR40'

            CurrentModel => CurrentModel%Next

        enddo

        !Checks the limits of the model (always model one which controls it)
        call GetModelTimeLimits (FirstModel%ModelID, GlobalBeginTime, GlobalEndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR50'

        !Checks if the time limits for each sub-model is 
        CurrentModel => FirstModel%Next
        do while (associated(CurrentModel))

            call GetModelTimeLimits (CurrentModel%ModelID, SubModelBeginTime, SubModelEndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR51'

            if(SubModelBeginTime .ne. GlobalBeginTime .or. &
               SubModelEndTime   .ne. GlobalEndTime)then
                
                write(*, *)"---------------- CHECKING SUB-MODEL DATES -----------------"
                write(*,*)
                write(*,*)"The initial date and/or the final date of the model:"
                write(*,*)
                write(*,*)trim(CurrentModel%ModelName)
                write(*,*)
                write(*,*)"Differs from the initial/final date of the father model"
                stop 'ConstructMohidWater - MohidWater - ERR52'

            end if
            CurrentModel => CurrentModel%Next
        enddo


        !Construct the Father Grid
        CurrentModel => FirstModel
        do while (associated(CurrentModel))
            if (associated(CurrentModel%FatherModel)) then

                call GetSubModelWindowON (CurrentModel%ModelID, SubModelWindowON,       &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR60'

                call GetNotDefinedCells(CurrentModel%FatherModel%HorizontalGridID,&
                                        NotDefinedCells, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR70'

                if (.not. SubModelWindowON .and. NotDefinedCells) stop 'ConstructMohidWater - MohidWater - ERR80'

                if (SubModelWindowON) then

                    call GetSubModelWindow (CurrentModel%ModelID, CurrentModel%FatherLink%Window,&
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR90'

                else

                    call GetHorizontalGridSize(CurrentModel%FatherModel%HorizontalGridID,&
                                               WorkSize = CurrentModel%FatherLink%Window,&
                                                 STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR100'

                endif 


                call ConstructFatherGridLocation(CurrentModel%HorizontalGridID,             &
                                                 CurrentModel%FatherModel%HorizontalGridID, &
                                                 Window = CurrentModel%FatherLink%Window,   &
                                                 STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR110'

            endif
            CurrentModel => CurrentModel%Next
        enddo

        !Verifies which Model needs Hydrodynamic / Waterproperties conditions from father
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            if (associated(CurrentModel%FatherModel)) then

                call GetHydroNeedsFather (CurrentModel%HydrodynamicID,                   &
                                          CurrentModel%FatherLink%Hydro,                 &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR120'


                call GetWaterNeedsFather (CurrentModel%WaterpropertiesID,                &
                                          CurrentModel%FatherLink%Water,                 &
                                          CurrentModel%FatherLink%nProp,                 &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR130'


                if (CurrentModel%FatherLink%Hydro) then
                    call SetHydroFather (CurrentModel%HydrodynamicID,                     &
                                         CurrentModel%FatherModel%HydrodynamicID,         &
                                         InitialField = .true., STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR140'
                endif

                if (CurrentModel%FatherLink%Water) then

                    allocate (CurrentModel%FatherLink%PropertyIDNumbers(CurrentModel%FatherLink%nProp))
                    call GetPropListNeedsFather(CurrentModel%WaterPropertiesID,              &
                                                CurrentModel%FatherLink%PropertyIDNumbers,   &
                                                STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR150'


                    do iProp = 1, CurrentModel%FatherLink%nProp

                        call SetWaterPropFather (CurrentModel%WaterPropertiesID,                    &
                                                 CurrentModel%FatherModel%WaterPropertiesID,        &
                                                 CurrentModel%FatherLink%PropertyIDNumbers(iProp),  &
                                                 InitialField= .true.,                              &
                                                 STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR160'
                    enddo

                endif


            endif
            CurrentModel => CurrentModel%Next
        enddo

#ifdef OVERLAP

        call ConstructModelOverlapping

#endif OVERLAP
    

        !Gets System Time
        call GetSystemTime (WorkCycleInitial)
        call cpu_time      (CPUTimeConstructor)
        
        !Update Time - OpenMI needs this in the constructor
        GlobalCurrentTime  = GlobalBeginTime

    end subroutine ConstructMohidWater

    !--------------------------------------------------------------------------

#ifdef OVERLAP
    subroutine ConstructModelOverlapping
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_MohidWater), pointer                :: CurrentModel
        
        !Begin-----------------------------------------------------------------


        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            call GetModelOverlap    (CurrentModel%ModelID, CurrentModel%Overlap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModelOverlapping - MohidWater - ERR10'
                
            if(CurrentModel%Overlap)then
                
                call GetModelOverlapInfo(CurrentModel%ModelID,                              &
                                         CurrentModel%OverlapModelID,                       &
                                         CurrentModel%OverlapCells,                         &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModelOverlapping - MohidWater - ERR20'


                call GetModelInstanceIDs(CurrentModel%ModelID,                              &
                                         HorizontalGridID  = CurrentModel%HorizontalGridID, &
                                         HydrodynamicID    = CurrentModel%HydrodynamicID,   &
                                         WaterpropertiesID = CurrentModel%WaterpropertiesID,&
                            STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModelOverlapping - MohidWater - ERR30'


                call GetModelInstanceIDs(CurrentModel%OverlapModelID,                       &
                                 HydrodynamicID    = CurrentModel%OverlapHydrodynamicID,    &
                                 WaterpropertiesID = CurrentModel%OverlapWaterpropertiesID, &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModelOverlapping - MohidWater - ERR40'

                call GetWaterOverlap      (CurrentModel%WaterpropertiesID,                  &
                                           CurrentModel%OverlapLink%Water,                  &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModelOverlapping - MohidWater - ERR50'

                call GetHydroOverlap      (CurrentModel%HydrodynamicID,                     &
                                           CurrentModel%OverlapLink%Hydro,                  &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModelOverlapping - MohidWater - ERR60'





            end if
            
            CurrentModel => CurrentModel%Next
        enddo

    end subroutine ConstructModelOverlapping

#endif OVERLAP
   
    !--------------------------------------------------------------------------

    subroutine ConstructMohidWaterMPI

#ifdef _USE_MPI
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(PathLength)                       :: WatchFile, DTLogFile
        type (T_MohidWater), pointer                :: CurrentModel
        integer                                     :: tag = 0
        integer                                     :: status(MPI_STATUS_SIZE)
        integer                                     :: iSub, length, iProp
        character(PathLength)                       :: FileNameStr
        integer, dimension(PathLength)              :: FileNameInt
        real, dimension(6)                          :: AuxBeg, AuxEnd
        integer                                     :: iProc
        logical                                     :: SubModelWindowON, NotDefinedCells
        type(T_Time)                                :: SubModelBeginTime, SubModelEndTime

        !Begin-----------------------------------------------------------------

        call MPI_INIT               (STAT_CALL )
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR10'
        
        call MPI_COMM_RANK          (MPI_COMM_WORLD, myMPI_ID, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR20'
        
        call MPI_COMM_SIZE          (MPI_COMM_WORLD, nMPI_Processes, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR30'
        
        call MPI_GET_PROCESSOR_NAME (myMPI_Processor, length, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR40'
        
        call MPI_Barrier            (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR50'

        !This must be shifted elsewhere
        allocate (ModelDTs          (nMPI_Processes))

        !Reads the file Tree.dat and constructs the list of models
        !Like this every model has access to the whole model list
        !Do this sequentially, so there will be no conflict accessing Tree.dat
        do iProc = 0, nMPI_Processes - 1

            if (iProc == myMPI_ID) then

                !Common Startup tasks
                call StartUpMohid   ("Mohid Water on "//trim(adjustl(myMPI_Processor)))
                call GetSystemTime  (InitialSystemTime)
                call ConstructModelList

            endif

            !Wait for all processes
            call MPI_Barrier   (MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR60'
        enddo

        !Construct all modules. Each module is constructed one after the other. Between each 
        !construction, the other modules wait at MPI_Barrier below.
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            if (myMPI_ID == CurrentModel%MPI_ID) then

                call SetFilesName  (CurrentModel%ModelPath)

                !Monitor Performance of the model execution?
                call ReadFileName('OUTWATCH', WatchFile, Message = 'Stop Watch File', STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    STAT_CALL   = CreateWatchGroup (WatchFile)
                    MonitorPerformance  = .true.
                else
                    MonitorPerformance  = .false.
                endif

                !Monitors DT
                call ReadFileName('DT_LOG', DTLogFile, Message = 'Start DTLog File', STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then            
                    MonitorDT  = .true.
                else
                    MonitorDT  = .false.
                endif

                if (MonitorDT) then
                    call UnitsManager (UnitDT, OPEN_FILE)      
                    open(UNIT   = UnitDT, FILE   = DTLogFile, STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR70'
                end if
                
                !PCL 
                write(*,*) "Construct modelo MPI ID =", CurrentModel%MPI_ID

                call ConstructModel(LagInstance, ModelNames, NumberOfModels,            &
                                    ObjLagrangianGlobal, CurrentModel%ModelID,          &
                                    InitialSystemTime, CurrentModel%MPI_ID,             &
                                    STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR80'
                
                !Get the Instance IDs of the Objects which are necessary for Model/SubModel
                !comunication
                call GetModelInstanceIDs      (CurrentModel%ModelID,                              &
                                               HorizontalGridID  = CurrentModel%HorizontalGridID, &
                                               HydrodynamicID    = CurrentModel%HydrodynamicID,   &
                                               WaterpropertiesID = CurrentModel%WaterpropertiesID,&
                                               STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR90'

                !Gets Current Time of the model
                call GetModelCurrentTime      (CurrentModel%ModelID, CurrentModel%CurrentTime,    &
                                               STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR100'

                !Sets the last info time 
                CurrentModel%InfoTime = CurrentModel%CurrentTime

            endif

            CurrentModel => CurrentModel%Next

            !Waits for all processes
            call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR110'

        enddo

        !Waits for all processes
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR120'

        !Sends the information about the Grid to the submodels
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            !Send the information to the submodels
            if (CurrentModel%MPI_ID == myMPI_ID) then
           
                do iSub = 1, CurrentModel%nSubModels
                        
                    !Gets the name of the file
                    call GetGridFileName          (CurrentModel%HorizontalGridID, FileNameStr, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR130'

                    !Converts a String to an integer Array
                    call ConvertStrToInt (FileNameStr, FileNameInt, PathLength)

                    !Send name to submodel
                    call MPI_Send (FileNameInt, PathLength, MPI_INTEGER,                 &
                                   CurrentModel%SubMPIID(iSub), tag, MPI_COMM_WORLD,     &
                                   STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR140'
                    
                enddo


                if (associated(CurrentModel%FatherModel)) then

                    call MPI_Recv (FileNameInt, PathLength, MPI_INTEGER, CurrentModel%FatherModel%MPI_ID, &
                                   tag, MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR150'
                    
                    !Converts an integer Array to a String
                    call ConvertIntToStr (FileNameInt, FileNameStr, PathLength)

                    !Constructs the Father Grid
                    call ConstructHorizontalGrid (HorizontalGridID = CurrentModel%FatherGridID, &
                                                  DataFile         = FileNameStr,               &
                                                  STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR160'


                    call GetNotDefinedCells(CurrentModel%FatherGridID,                          &
                                            NotDefinedCells, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR170'


                    call GetSubModelWindowON (CurrentModel%ModelID, SubModelWindowON,           &
                                              STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR180'

                    if (SubModelWindowON) then

                        call GetSubModelWindow (CurrentModel%ModelID, CurrentModel%FatherLink%Window,&
                                                STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR190'

                    else

                        call GetHorizontalGridSize(CurrentModel%FatherGridID,                   &
                                                   WorkSize = CurrentModel%FatherLink%Window,   &
                                                   STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR200'

                    endif 
                    
                    if (.not. SubModelWindowON .and. NotDefinedCells) stop 'ConstructMohidWaterMPI - MohidWater - ERR210'

                    !Constructs the location of the father grid
                    call ConstructFatherGridLocation(CurrentModel%HorizontalGridID,             &
                                                     CurrentModel%FatherGridID,                 &
                                                     Window = CurrentModel%FatherLink%Window,   &
                                                     STAT   = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR220'

                endif

            endif

            CurrentModel => CurrentModel%Next

        enddo  
             
        !Wait for all
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR230'
           
        !Sends to all processes the information about GlobalBeginTime end GlobalEndTime
        if (myMPI_ID == 0) then

            !Checks the limits of the model (always model one which controls it)
            call GetModelTimeLimits (FirstModel%ModelID, GlobalBeginTime, GlobalEndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR240'

            call ExtractDate (GlobalBeginTime, AuxBeg(1), AuxBeg(2), AuxBeg(3),          &
                              AuxBeg(4), AuxBeg(5), AuxBeg(6))

            call ExtractDate (GlobalEndTime,   AuxEnd(1), AuxEnd(2), AuxEnd(3),          &
                              AuxEnd(4), AuxEnd(5), AuxEnd(6))

        endif

        Precision = MPIKind(AuxBeg(1))

        call MPI_Bcast(AuxBeg, 6, Precision, 0, MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR250'
        call MPI_Bcast(AuxEnd, 6, Precision, 0, MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR260'
        
if1 :   if (myMPI_ID /= 0) then
            
            CurrentModel => FirstModel%Next
do1 :       do while (associated(CurrentModel))
            if (myMPI_ID == CurrentModel%MPI_ID) then
                
                
            
            call SetDate (GlobalBeginTime, AuxBeg(1), AuxBeg(2), AuxBeg(3),              &
                          AuxBeg(4), AuxBeg(5), AuxBeg(6))

            call SetDate (GlobalEndTime, AuxEnd(1), AuxEnd(2), AuxEnd(3),                &
                          AuxEnd(4), AuxEnd(5), AuxEnd(6))





            write (*,*) 'MPIProcess: ', CurrentModel%MPI_ID
            call GetModelTimeLimits (CurrentModel%ModelID, SubModelBeginTime, SubModelEndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then 
                write (*,*) 'MPIProcess: ', CurrentModel%MPI_ID
                stop        'ConstructMohidWaterMPI - Main - ERR270'
            endif

if2 :       if(SubModelBeginTime .ne. GlobalBeginTime .or. &
               SubModelEndTime   .ne. GlobalEndTime)then
                
                write(*, *)"---------------- CHECKING SUB-MODEL DATES -----------------"
                write(*,*)
                write(*,*)"The initial date and/or the final date of the model:"
                write(*,*)
                write(*,*)trim(CurrentModel%ModelName)
                write(*,*)
                write(*,*)"Differs from the initial/final date of the father model"
                stop 'ConstructMohidWaterMPI - MohidWater - ERR280'

            end if if2
                    
            
            endif
                CurrentModel => CurrentModel%Next
            enddo do1
        endif if1

        !Wait for all
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR290'
        

        !Fills in FatherLink
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            !Send the information to the submodels
            if (CurrentModel%MPI_ID == myMPI_ID) then

                if (associated(CurrentModel%FatherModel)) then

                    call GetHydroNeedsFather (CurrentModel%HydrodynamicID,               &
                                              CurrentModel%FatherLink%Hydro,             &
                                              STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR300'
                    

                    call GetWaterNeedsFather (CurrentModel%WaterpropertiesID,            &
                                              CurrentModel%FatherLink%Water,             &
                                              CurrentModel%FatherLink%nProp,             &
                                              STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWater - MohidWater - ERR310'
                    
                    if (CurrentModel%FatherLink%Hydro .or. CurrentModel%FatherLink%Water) then
                        CurrentModel%FatherLink%Nesting = .true.
                    else
                        CurrentModel%FatherLink%Nesting = .false.
                    endif


                    if (CurrentModel%FatherLink%Water) then

                        allocate (CurrentModel%FatherLink%PropertyIDNumbers(CurrentModel%FatherLink%nProp))
                        call GetPropListNeedsFather(CurrentModel%WaterPropertiesID,                 &
                                                    CurrentModel%FatherLink%PropertyIDNumbers,      &
                                                    STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR320'

                    endif
                endif

            endif

            CurrentModel => CurrentModel%Next

        enddo  

        !Waits for all processes
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR330'

        !Fills in SubModelLink
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            !Send the information to the submodels
            if (CurrentModel%MPI_ID == myMPI_ID) then

                do iSub = 1, CurrentModel%nSubModels

                    !Recieves information if submodel needs 
                    call MPI_Recv (CurrentModel%SubmodelLink(iSub)%Hydro, 1, MPI_LOGICAL, &
                                   CurrentModel%SubMPIID(iSub), tag, MPI_COMM_WORLD,      &
                                   status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR340'

                    call MPI_Recv (CurrentModel%SubmodelLink(iSub)%Water, 1, MPI_LOGICAL, &
                                   CurrentModel%SubMPIID(iSub), tag, MPI_COMM_WORLD,      &
                                   status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR350'
                    
                    if (CurrentModel%SubmodelLink(iSub)%Water) then

                        call MPI_Recv (CurrentModel%SubmodelLink(iSub)%nProp, 1, MPI_INTEGER, &
                                       CurrentModel%SubMPIID(iSub), tag, MPI_COMM_WORLD,      &
                                       status, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR360'

                        allocate (CurrentModel%SubmodelLink(iSub)%PropertyIDNumbers(CurrentModel%SubmodelLink(iSub)%nProp))

                        call MPI_Recv (CurrentModel%SubmodelLink(iSub)%PropertyIDNumbers,  &
                                       CurrentModel%SubmodelLink(iSub)%nProp, MPI_INTEGER, &
                                       CurrentModel%SubMPIID(iSub), tag, MPI_COMM_WORLD,    &
                                       status, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR370'

                    endif
                
                    call MPI_Recv (CurrentModel%SubmodelLink(iSub)%Window%ILB, 1,         &
                                   MPI_INTEGER, CurrentModel%SubMPIID(iSub), tag,         &
                                   MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR380'
                    
                    call MPI_Recv (CurrentModel%SubmodelLink(iSub)%Window%IUB, 1,         &
                                   MPI_INTEGER, CurrentModel%SubMPIID(iSub), tag,         &
                                   MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR390'
                                   
                    call MPI_Recv (CurrentModel%SubmodelLink(iSub)%Window%JLB, 1,         &
                                   MPI_INTEGER, CurrentModel%SubMPIID(iSub), tag,         &
                                   MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR400'

                    call MPI_Recv (CurrentModel%SubmodelLink(iSub)%Window%JUB, 1,         &
                                   MPI_INTEGER, CurrentModel%SubMPIID(iSub), tag,         &
                                   MPI_COMM_WORLD, status, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR410'
                    

                enddo


                if (associated(CurrentModel%FatherModel)) then

                    !Send to father information about Father link
                    call MPI_Send (CurrentModel%FatherLink%Hydro, 1, MPI_LOGICAL,        &
                                   CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR420'

                    call MPI_Send (CurrentModel%FatherLink%Water, 1, MPI_LOGICAL,        &
                                   CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR430'
                    
                    if (CurrentModel%FatherLink%Water) then
                    
                        call MPI_Send (CurrentModel%FatherLink%nProp, 1, MPI_INTEGER,        &
                                       CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR440'
                        
                        call MPI_Send (CurrentModel%FatherLink%PropertyIDNumbers,            &
                                       CurrentModel%FatherLink%nProp, MPI_INTEGER,           & 
                                       CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR450'

                    endif 
                    
                    call MPI_Send (CurrentModel%FatherLink%Window%ILB, 1, MPI_INTEGER,   &
                                   CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR460'
                    
                    call MPI_Send (CurrentModel%FatherLink%Window%IUB, 1, MPI_INTEGER,   &
                                   CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR470'
                    
                    call MPI_Send (CurrentModel%FatherLink%Window%JLB, 1, MPI_INTEGER,   &
                                   CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR480'
                    
                    call MPI_Send (CurrentModel%FatherLink%Window%JUB, 1, MPI_INTEGER,   &
                                   CurrentModel%FatherModel%MPI_ID, tag, MPI_COMM_WORLD, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR490'

                endif

            endif

            CurrentModel => CurrentModel%Next

        enddo  

        !Waits for all processes
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR500'

        !Sends Initial Hydro Information.
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            !Send the information to the submodels
            if (CurrentModel%MPI_ID == myMPI_ID) then

                do iSub = 1, CurrentModel%nSubModels
                    
                    if (CurrentModel%SubmodelLink(iSub)%Hydro) then
                        call SendHydrodynamicMPI (CurrentModel%HydrodynamicID,           &
                                                  CurrentModel%SubMPIID(iSub),           &
                                                  CurrentModel%SubmodelLink(iSub)%Window,&
                                                  InitialField = .true.,                 &
                                                  STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR510'
                    end if

                enddo

                if (associated(CurrentModel%FatherModel)) then
                    if (CurrentModel%FatherLink%Hydro) then
                        call RecvHydrodynamicMPI (CurrentModel%HydrodynamicID,              &
                                                  CurrentModel%FatherModel%MPI_ID,          &
                                                  CurrentModel%FatherLink%Window,           &
                                                  InitialField = .true.,                    &
                                                  FatherGridID = CurrentModel%FatherGridID, &
                                                  STAT = STAT_CALL)                
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR520'
                    endif
                endif

            endif            

            CurrentModel => CurrentModel%Next

        enddo

        !Waits for all processes
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR530'

        !Sends Initial Water Information.
        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            !Send the information to the submodels
            if (CurrentModel%MPI_ID == myMPI_ID) then

                do iSub = 1, CurrentModel%nSubModels

                    if (CurrentModel%SubmodelLink(iSub)%Water) then

                        do iProp = 1, CurrentModel%SubmodelLink(iSub)%nProp
                            
                            call SendWaterPropertiesMPI (CurrentModel%WaterPropertiesID,                                    &
                                             CurrentModel%SubMPIID(iSub),                                       &
                                             CurrentModel%SubmodelLink(iSub)%Window,                            &
                                             InitialField = .true.,                                             &
                                             PropIDNumber = CurrentModel%SubmodelLink(iSub)%PropertyIDNumbers(iProp),   &
                                             STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR540'
                        enddo

                    endif
                enddo

                if (associated(CurrentModel%FatherModel)) then

                    if (CurrentModel%FatherLink%Water) then
                        
                        do iProp = 1, CurrentModel%FatherLink%nProp

                            call RecvWaterPropertiesMPI (CurrentModel%WaterPropertiesID,           &
                                                         CurrentModel%FatherModel%MPI_ID,          &
                                                         CurrentModel%FatherLink%Window,           &
                                                         InitialField = .true.,                    &
                                                         FatherGridID = CurrentModel%FatherGridID, &
                                                         PropIDNumber = CurrentModel%FatherLink%PropertyIDNumbers(iProp), &
                                                         STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR550'
                        end do

                    end if

                endif


            endif            

            CurrentModel => CurrentModel%Next

        enddo

        !Waits for all processes
        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidWaterMPI - MohidWater - ERR560'


        !Gets System Time
        call GetSystemTime (WorkCycleInitial)
        call cpu_time      (CPUTimeConstructor)

        !Update Time - OpenMI needs this in the constructor
        GlobalCurrentTime  = GlobalBeginTime

#endif

    end subroutine ConstructMohidWaterMPI

    !--------------------------------------------------------------------------

    subroutine ConstructModelList
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: TreeExists
        integer                                     :: iTree, STAT_CALL, i
        integer                                     :: MPI_ID
        type (T_MohidWater), pointer                :: CurrentModel, NextModel
        character(StringLength)                     :: Coment1, Coment2
        character(PathLength)                       :: AuxChar


        !Nullifies First Pointer
        nullify(FirstModel)

        !Verifies if Tree file exists and allocates the list of models
        inquire(file='tree.dat', EXIST = TreeExists)
      
        if (TreeExists) then

            call UnitsManager(iTree, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModelList - MohidWater - ERR01'

            open(UNIT = iTree, FILE = 'tree.dat', status = 'OLD')

            read(unit=iTree, fmt=*) Coment1
            read(unit=iTree, fmt=*) Coment2

            MPI_ID  = 0
            Do 
                read(unit=iTree, fmt='(a256)', END=100) AuxChar
                AuxChar = trim(adjustl(AuxChar))
                if (AuxChar(1:1) == '+') then
                    if (RunInParallel) then
                        call AddNewModel(trim(adjustl(AuxChar)), MPI_ID)
                    else
                        call AddNewModel(trim(adjustl(AuxChar)))
                    endif
                    MPI_ID = MPI_ID + 1
                else
                    exit
                endif
            enddo 

100         call UnitsManager(iTree, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModelList - MohidWater - ERR02'

            NumberOfModels = MPI_ID

        else

            call AddNewModel("*")

            NumberOfModels = 1

        endif

        allocate(ModelNames(1:NumberOfModels))
        allocate(LagInstance(1:TotalLagInst_, 1:NumberOfModels))
        LagInstance(:,:) = 0


        
        !Determines which models have submodels
        CurrentModel => FirstModel
        i = 0
        do while (associated(CurrentModel))
            i = i + 1
            ModelNames(i) = trim(CurrentModel%ModelName)

            NextModel => CurrentModel%Next
doNext:     do while (associated(NextModel))
                
                if (NextModel%ModelLevel  <= CurrentModel%ModelLevel) exit doNext
                
                if (NextModel%ModelLevel  == CurrentModel%ModelLevel+1) then
                    CurrentModel%SubOn     = .true.
                    CurrentModel%nSubModels= CurrentModel%nSubModels + 1
                    NextModel%FatherModel  => CurrentModel
                    if (RunInParallel) then
                        CurrentModel%SubMPIID(CurrentModel%nSubModels) = NextModel%MPI_ID
                    endif
                endif

                NextModel => NextModel%Next
            enddo doNext

            MPI_ID       =  MPI_ID + 1
            CurrentModel => CurrentModel%Next

        enddo
                
    end subroutine ConstructModelList

    !--------------------------------------------------------------------------

    subroutine AddNewModel(AuxString, MPI_ID)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: AuxString
        integer, optional                           :: MPI_ID

        !Local-----------------------------------------------------------------
        type (T_MohidWater), pointer                :: NewModel
        type (T_MohidWater), pointer                :: PreviousModel
        type (T_MohidWater), pointer                :: CurrentModel

        !Allocates new model
        allocate (NewModel)
        nullify  (NewModel%FatherModel       )
        nullify  (NewModel%Next              )
        nullify  (NewModel%FatherLink%PropertyIDNumbers)

        !Sets data
        NewModel%ModelID      = 0
        if (present(MPI_ID)) then
            NewModel%MPI_ID       = MPI_ID
        else
            NewModel%MPI_ID       = null_int
        endif
!        NewModel%FatherMPI_ID = null_int
        NewModel%ModelLevel   = ModelLevel(AuxString)
        NewModel%ModelPath    = ModelPath (AuxString, NewModel%ModelLevel)
        NewModel%ModelName    = ModelName (NewModel%ModelPath)
        NewModel%SubOn        = .false.
        NewModel%nSubModels   = 0
        
        if(RunInParallel) then
            allocate (NewModel%SubMPIID    (nMPI_Processes))
            allocate (NewModel%SubmodelLink(nMPI_Processes))
        else
            nullify  (NewModel%SubMPIID    )
            nullify  (NewModel%SubmodelLink)
        endif

        !Insert into model list
        if (.not. associated(FirstModel)) then
            FirstModel    => NewModel
        else
            PreviousModel => FirstModel
            CurrentModel  => FirstModel%Next
            do while (associated(CurrentModel))
                PreviousModel  => CurrentModel
                CurrentModel   => CurrentModel%Next
            enddo
            CurrentModel       => NewModel
            PreviousModel%Next => NewModel
        endif

    end subroutine AddNewModel
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidWater
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: Running
        !integer                                     :: STAT_CALL
        !type (T_MohidWater), pointer                :: CurrentModel
        !logical                                     :: DoNextStep
        !real                                        :: DT_Father
        real                                        :: DTmin, DTmax

#ifndef _OUTPUT_OFF_
        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Running MOHID, please wait..."
        write(*, *)                    
#endif
        if (MonitorPerformance) then
            call StartWatch ("Main", "ModifyMohidWater")
        endif

        Running            = .true.

        !Search for initial Min and Max Time Step
        DTmin   = - FillValueReal
        DTmax   =   FillValueReal
        call SearchMinMaxTimeStep (DTmin, DTmax)

        do while (Running)
            
            !The which was here has been refactored into a function, 
            !so it can be called from here and from OpenMP
            Running = DoOneTimeStep (DTmin)
            
        enddo
        
        if (MonitorPerformance) then
            call StopWatch ("Main", "ModifyMohidWater")
        endif

    end subroutine ModifyMohidWater
    
    !--------------------------------------------------------------------------

    logical function DoOneTimeStep (DTMinParallel)

        !Arguments-------------------------------------------------------------
        real                                        :: DTMinParallel
        !Local-----------------------------------------------------------------
        type (T_MohidWater), pointer                :: CurrentModel
        logical                                     :: DoNextStep
        integer                                     :: STAT_CALL
        real                                        :: DTmin, DTmax, DT_Father


        !Search for initial Min and Max Time Step
        
        if (RunInParallel) then
            DTmin   =   DTMinParallel
        else
            DTmin   = - FillValueReal
            DTmax   =   FillValueReal
            call SearchMinMaxTimeStep (DTmin, DTmax)
        endif
        
        GlobalCurrentTime = GlobalCurrentTime + DTmin
        
        if (DTmin == 0.) then
            write(*,*) 'Time step equal to zero dt =', dtmin
            DoOneTimeStep = .false.
            return 
        endif
        
        if (RunInParallel) then
        
            write(*,*) 'Start run MPI'
        
            CurrentModel => FirstModel
            
            do while (associated(CurrentModel))
            
                write(*,*) 'Start CurrentModel%MPI_ID == myMPI_ID'
            
                if (CurrentModel%MPI_ID == myMPI_ID) then
                
                    write(*,*) 'UpdateTimeAndMapping Main'    
                    call UpdateTimeAndMapping    (CurrentModel%ModelID, GlobalCurrentTime, DoNextStep, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidWater - MohidWater - ERR10'
                    
                    if (associated(CurrentModel%FatherModel))  then
                    
                    
                        !
                        !THESE NEXT LINES HAVE TO BE RECODED. The used approach only works if nesting is done like
                        !Main
                        !   Sub
                        !       SubSub
                        !Frank Fev 2011                        
                    
                        DT_Father = ModelDTs(CurrentModel%FatherModel%MPI_ID+1)
                        !write(*,*)"CurrentModel%FatherModel%MPI_ID: ", CurrentModel%FatherModel%MPI_ID
                        !write(*,*)"DT_Father                      : ", DT_Father
                        !call GetModelTimeStep (CurrentModel%FatherModel%ModelID, DT_Father, STAT = STAT_CALL)
                        !if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidWater - MohidWater - ERR20'
                    else
                        DT_Father = - FillValueReal
                    endif

                    if (DoNextStep) then
                        !Waits for information from father
                        if (associated(CurrentModel%FatherModel)) then
                            
                            if (CurrentModel%FatherLink%Nesting) then
                        
                                !Post next recieve...
                                if (CurrentModel%CurrentTime == CurrentModel%InfoTime) then
                                    write(*,*) 'ReceiveInformationMPI Main'
                                    call ReceiveInformationMPI (CurrentModel)
                                else if (CurrentModel%CurrentTime < CurrentModel%InfoTime) then 
                                    write(*,*) 'UpdateSubModelValues Main'
                                    call UpdateSubModelValues (CurrentModel)
                                else if (CurrentModel%CurrentTime > CurrentModel%InfoTime) then

                                    stop 'ModifyMohidWater - MohidWater - ERR30'

                                endif
                                
                            endif                                

                        endif
                        write(*,*) 'RunModel Main'
                        call RunModel(CurrentModel%ModelID, DT_Father, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidWater - MohidWater - ERR40'
                        write(*,*) 'End RunModel Main'
                        if (CurrentModel%FatherLink%Nesting) then
                            write(*,*) 'SendInformationMPI Main'
                            call SendInformationMPI (CurrentModel)
                        endif
                        write(*,*) 'end if 1' 
                    endif
                    write(*,*) 'end if 2'
                endif
                write(*,*) 'end if 3'
                
                CurrentModel => CurrentModel%Next
                
            enddo
            
            write(*,*) ' Finish run cycle'
            
        else

            CurrentModel => FirstModel
            do while (associated(CurrentModel))

                call UpdateTimeAndMapping (CurrentModel%ModelID, GlobalCurrentTime,  &
                                           DoNextStep, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidWater - MohidWater - ERR50'

                if (DoNextStep) then    
                    call SubModelComunication     (CurrentModel)

#ifdef OVERLAP
                    call OverlapModelCommunication(CurrentModel)
#endif OVERLAP


                    if (associated(CurrentModel%FatherModel))  then
                        call GetModelTimeStep (CurrentModel%FatherModel%ModelID, DT_Father, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidWater - MohidWater - ERR60'
                    else
                        DT_Father = - FillValueReal
                    endif

                    call RunModel             (CurrentModel%ModelID, DT_Father, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidWater - MohidWater - ERR70'
                    
                endif

                CurrentModel => CurrentModel%Next
            enddo

            !Search again MinMax, so the test DTmin / 10.0 can be safely done
            !
            DTmin   = - FillValueReal
            DTmax   =   FillValueReal
            call SearchMinMaxTimeStep (DTmin, DTmax)
        endif
        
        if (abs(GlobalCurrentTime - GlobalEndTime) > DTmin / 10.0) then
            DoOneTimeStep = .true.
        else
            DoOneTimeStep = .false.
        endif
        
        
    end function DoOneTimeStep
    
    !--------------------------------------------------------------------------

    subroutine KillMohidWater

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_MohidWater), pointer                :: CurrentModel

        !Gets System Time
        call GetSystemTime (WorkCycleFinal )
        call cpu_time      (CPUTimeModifier)

#ifndef _OUTPUT_OFF_
        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Shuting down MOHID, please wait..."
        write(*, *)                    
#endif

        CurrentModel => FirstModel
        do while (associated(CurrentModel))

            call KillModel(CurrentModel%ModelID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillMohidWater - MohidWater - ERR01'

            CurrentModel => CurrentModel%Next
        enddo

        if (MonitorPerformance) then
            call KillWatchGroup (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillMohidWater - MohidWater - ERR02'
        endif

        if (MonitorDT) call UnitsManager (UnitDT, CLOSE_FILE)

        !Total Time
        call GetSystemTime (FinalSystemTime)
        call cpu_time      (TotalCPUTime   )
        ElapsedSeconds = FinalSystemTime - InitialSystemTime

        !WorkCycle Time
        WorkCycleCPUTime = CPUTimeModifier - CPUTimeConstructor
        WorkCycleElapsed = WorkCycleFinal  - WorkCycleInitial

        call ShutdownMohid ("Mohid Water", ElapsedSeconds, TotalCPUTime,                 &
                            WorkCycleElapsed, WorkCycleCPUTime)

    end subroutine KillMohidWater
    
    !--------------------------------------------------------------------------

    subroutine KillMohidWaterMPI

        !Arguments-------------------------------------------------------------
#ifdef _USE_MPI

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_MohidWater), pointer                :: CurrentModel

        !Gets System Time
        call GetSystemTime (WorkCycleFinal)
        call cpu_time      (CPUTimeModifier)


        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Shuting down MOHID, please wait..."
        write(*, *)                    


        CurrentModel => FirstModel
        do while (associated(CurrentModel))
            !The next if statement is important due to the MPI
            if (CurrentModel%MPI_ID == myMPI_ID) then
                call KillModel(CurrentModel%ModelID, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillMohidWaterMPI - MohidWater - ERR01'
            endif

            !Wait for all processes
            call MPI_Barrier   (MPI_COMM_WORLD, STAT_CALL)

            CurrentModel => CurrentModel%Next
        enddo

        if (MonitorPerformance) then
            call KillWatchGroup (STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillMohidWaterMPI - MohidWater - ERR02'
        endif

        if (MonitorDT) call UnitsManager (UnitDT, CLOSE_FILE)

        call MPI_Barrier  (MPI_COMM_WORLD, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillMohidWaterMPI - MohidWater - ERR03'

        call MPI_Finalize (STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillMohidWaterMPI - MohidWater - ERR04'

        !Total Time
        call GetSystemTime (FinalSystemTime)
        call cpu_time      (TotalCPUTime   )
        ElapsedSeconds = FinalSystemTime - InitialSystemTime

        !WorkCycle Time
        WorkCycleCPUTime = CPUTimeModifier - CPUTimeConstructor
        WorkCycleElapsed = WorkCycleFinal  - WorkCycleInitial

        call ShutdownMohid ("Mohid Water", ElapsedSeconds, TotalCPUTime,                 &
                            WorkCycleElapsed, WorkCycleCPUTime)

#endif _USE_MPI

    end subroutine KillMohidWaterMPI
    
    !--------------------------------------------------------------------------

    subroutine GetSystemTime (Time)

        !Arguments-------------------------------------------------------------
        type (T_Time)                               :: Time

        !Local-----------------------------------------------------------------
        integer, dimension(8)                       :: F95Time

        call date_and_time(Values = F95Time)
        call SetDate      (Time, float(F95Time(1)), float(F95Time(2)),      &
                                 float(F95Time(3)), float(F95Time(5)),      &
                                 float(F95Time(6)), float(F95Time(7))+      &
                                 F95Time(8)/1000.)

    end subroutine GetSystemTime

    !--------------------------------------------------------------------------

    integer function ModelLevel (String)

        !Arguments-------------------------------------------------------------
        Character(len=*), intent(in)                :: String

        !Local-----------------------------------------------------------------
        integer                                     :: Level, i

        if(String(1:1)/='+') then
            ModelLevel=0
        else
            Level=1
do1:        do i=2,StringLength
                if (String(i:i)/='+') exit
                Level=Level+1
            enddo do1
            ModelLevel = Level
        endif

    end function ModelLevel

    !--------------------------------------------------------------------------

    character(len=PathLength) function ModelName (ModelPath)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in)                :: ModelPath

        !Local-----------------------------------------------------------------
        integer                                     :: LenName, i 
        !Begin-----------------------------------------------------------------
        ModelName = ModelPath

        LenName          = len(trim(ModelName))
        !Remove from the name "/exe"
        LenName          = LenName-4
        ModelName        = ModelName(1:LenName)
        do i =  LenName, 1, -1
            if (ModelName(i:i) =='\') then
                ModelName = ModelName(i+1:LenName)
                exit
            endif
        enddo

    end function ModelName

    !--------------------------------------------------------------------------

        !--------------------------------------------------------------------------

    character(len=PathLength) function ModelPath (String, Level)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in)                :: String
        integer,          intent(in)                :: Level

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        ModelPath = ''
        ModelPath = String(Level+1:)

    end function ModelPath

    !--------------------------------------------------------------------------

    subroutine SetFilesName (Folder)

    !Arguments------------------------------------
    character(LEN=*) :: Folder

    !Begin------------------------------------

    if (trim(Folder) /= '*') FilesName = trim(Folder)//'/nomfich.dat'

    end subroutine SetFilesName

    !--------------------------------------------------------------------------

#ifdef _USE_MPI
    subroutine ConvertIntToStr (FileNameInt, FileNameStr, Length)

        !Arguments-------------------------------------------------------------
        integer                                 :: Length
        integer, dimension(Length)              :: FileNameInt
        character(len = Length)                 :: FileNameStr

        !Local-----------------------------------------------------------------
        integer                                 :: i

        do i = 1, Length
            FileNameStr(i:i) = char(FileNameInt(i))
        enddo

    end subroutine ConvertIntToStr 

    !--------------------------------------------------------------------------

    subroutine ConvertStrToInt (FileNameStr, FileNameInt, Length)

        !Arguments-------------------------------------------------------------
        integer                                 :: Length
        character(len = Length)                 :: FileNameStr
        integer, dimension(Length)              :: FileNameInt

        !Local-----------------------------------------------------------------
        integer                                 :: i
        do i = 1, Length
             FileNameInt(i) = ichar(FileNameStr(i:i))
        enddo

    end subroutine ConvertStrToInt 
#endif _USE_MPI
    !--------------------------------------------------------------------------

    subroutine SearchMinMaxTimeStep (DTmin, DTmax)

        !Arguments-------------------------------------------------------------
        real, intent(INOUT)                         :: DTmin, DTmax
        


#ifdef _USE_MPI

        !Begin-----------------------------------------------------------------

        call UpdateTimeStepsByMPI()

        !Determines the minimum time step
        DTmin = minval(ModelDTs)
        DTmax = maxval(ModelDTs)
        
#else _USE_MPI

        !Local-----------------------------------------------------------------
        real                                        :: DT
        type (T_MohidWater), pointer                :: CurrentModel
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------        
                
        CurrentModel => FirstModel
        do while (associated(CurrentModel))
        
            call GetModelTimeStep (CurrentModel%ModelID, DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SearchMinMaxTimeStep - MohidWater - ERR01'
        
            DTmin = min(DT, DTmin)
            DTmax = max(DT, DTmax)

            CurrentModel => CurrentModel%Next
        enddo

#endif _USE_MPI

    end subroutine SearchMinMaxTimeStep

    !--------------------------------------------------------------------------

#ifdef _USE_MPI
    subroutine UpdateTimeStepsByMPI()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real                                        :: DT
        integer                                     :: STAT_CALL
        type (T_MohidWater), pointer                :: CurrentModel
        integer                                     :: status(MPI_STATUS_SIZE)
        logical                                     :: Variable

        !Get the DT of the current model
        CurrentModel => FirstModel
        do while (associated(CurrentModel))
            if (CurrentModel%MPI_ID == myMPI_ID) then
                call GetModelTimeStep (CurrentModel%ModelID, DT, Variable, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeStepsByMPI - MohidWater - ERR10'
                
                if (associated(CurrentModel%FatherModel))  then
                    if (Variable) then
                        write(*,*) 'Parallel processing and nesting can not have variable DT'
                        stop 'UpdateTimeStepsByMPI - MohidWater - ERR20'
                    endif
                endif
            endif
            CurrentModel => CurrentModel%Next
        enddo

        !Sends the information to all the others 
        call MPI_Allgather(DT, 1, Precision, ModelDTs, 1, Precision, MPI_COMM_WORLD, status)
    
    end subroutine UpdateTimeStepsByMPI
#endif _USE_MPI

    !--------------------------------------------------------------------------

#ifdef OVERLAP

    subroutine OverlapModelCommunication(CurrentModel)

        !Arguments-------------------------------------------------------------
        type (T_MohidWater), pointer                :: CurrentModel

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL 

        !Begin-----------------------------------------------------------------

        if(CurrentModel%OverlapLink%Water)then

            if(CurrentModel%WaterPropertiesID < CurrentModel%OverlapWaterpropertiesID)then

                call SetModelOverlapWater(CurrentModel%WaterPropertiesID,                   &
                                          CurrentModel%OverlapWaterpropertiesID,            &
                                          CurrentModel%OverlapCells,                        &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OverlapModelCommunication - MohidWater - ERR01'
            
            endif
        endif

        if(CurrentModel%OverlapLink%Hydro)then

            if(CurrentModel%HydrodynamicID < CurrentModel%OverlapHydrodynamicID)then

                call SetModelOverlapHydro(CurrentModel%HydrodynamicID,                      &
                                          CurrentModel%OverlapHydrodynamicID,               &
                                          CurrentModel%OverlapCells,                        &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OverlapModelCommunication - MohidWater - ERR01'
            
            endif
        endif


    end subroutine OverlapModelCommunication

#endif OVERLAP
  
    !--------------------------------------------------------------------------

    subroutine SubModelComunication (CurrentModel)

        !Arguments-------------------------------------------------------------
        type (T_MohidWater), pointer                :: CurrentModel

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i 

        !----------------------------------------------------------------------

        if (CurrentModel%FatherLink%Hydro) then

             call SetHydroFather (CurrentModel%HydrodynamicID,                           &
                                  CurrentModel%FatherModel%HydrodynamicID,               &
                                  InitialField = .false., STAT = STAT_CALL)
             if (STAT_CALL /= SUCCESS_) stop 'SubModelComunication - MohidWater - ERR01'

        endif

        if (CurrentModel%FatherLink%Water) then

            do i = 1, CurrentModel%FatherLink%nProp

                call SetWaterPropFather (CurrentModel%WaterPropertiesID,                 &
                                         CurrentModel%FatherModel%WaterPropertiesID,     &
                                         CurrentModel%FatherLink%PropertyIDNumbers(i),   &
                                         InitialField= .false., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'SubModelComunication - MohidWater - ERR02'

            enddo

        endif

    end subroutine SubModelComunication

    !--------------------------------------------------------------------------

    subroutine ReceiveInformationMPI (CurrentModel)

        !Arguments-------------------------------------------------------------
        type (T_MohidWater), pointer                :: CurrentModel

#ifdef _USE_MPI
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iProp
        integer                                     :: status(MPI_STATUS_SIZE)
        real, dimension(6)                          :: AuxTime
        
        !Begin-----------------------------------------------------------------

        call MPI_Recv (AuxTime, 6, Precision, CurrentModel%FatherModel%MPI_ID, 999, MPI_COMM_WORLD, status,  &
                       STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReceiveInformationMPI - MohidWater - ERR01'

        call SetDate     (CurrentModel%InfoTime, AuxTime(1), AuxTime(2), AuxTime(3),     &
                          AuxTime(4), AuxTime(5), AuxTime(6))

        if (CurrentModel%FatherLink%Hydro) then
            call RecvHydrodynamicMPI (CurrentModel%HydrodynamicID,             &
                                      CurrentModel%FatherModel%MPI_ID,         &
                                      CurrentModel%FatherLink%Window,          &
                                      InitialField = .false.,                  &
                                      FatherGridID = CurrentModel%FatherGridID,&
                                      STAT = STAT_CALL)
        endif

        if (CurrentModel%FatherLink%Water) then
            
            do iProp = 1, CurrentModel%FatherLink%nProp
                
                call RecvWaterPropertiesMPI (CurrentModel%WaterPropertiesID,           &
                                             CurrentModel%FatherModel%MPI_ID,          &
                                             CurrentModel%FatherLink%Window,           &
                                             InitialField = .false.,                   &
                                             FatherGridID = CurrentModel%FatherGridID, &
                                             PropIDNumber = CurrentModel%FatherLink%PropertyIDNumbers(iProp), &
                                             STAT = STAT_CALL)
            end do

        end if

#else _USE_MPI

    !Disable unused variable warnings
    CurrentModel = CurrentModel

#endif _USE_MPI

    end subroutine ReceiveInformationMPI

    !--------------------------------------------------------------------------

    subroutine SendInformationMPI (CurrentModel)

        !Arguments-------------------------------------------------------------
        type (T_MohidWater), pointer                :: CurrentModel

#ifdef _USE_MPI
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iSub, iProp
        real, dimension(6)                          :: AuxTime

        !Gets the Current Time
        call GetModelCurrentTime (CurrentModel%ModelID, CurrentModel%CurrentTime,        &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'SendInformationMPI - MohidWater - ERR01'

        call ExtractDate (CurrentModel%CurrentTime, AuxTime(1), AuxTime(2), AuxTime(3),  &
                          AuxTime(4), AuxTime(5), AuxTime(6))

        do iSub = 1, CurrentModel%nSubModels
        
            write(*,*) 'Run mpi', AuxTime

            call MPI_Send (AuxTime, 6, Precision, CurrentModel%SubMPIID(iSub), 999,       &
                           MPI_COMM_WORLD, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SendInformationMPI - MohidWater - ERR02'

            if (CurrentModel%SubmodelLink(iSub)%Hydro) then
                call SendHydrodynamicMPI (CurrentModel%HydrodynamicID,                   &
                                          CurrentModel%SubMPIID(iSub),                   &
                                          CurrentModel%SubmodelLink(iSub)%Window,        &
                                          InitialField = .false., STAT = STAT_CALL)
            endif

            if(CurrentModel%SubmodelLink(iSub)%Water)then

                do iProp = 1, CurrentModel%SubmodelLink(iSub)%nProp

                    call SendWaterPropertiesMPI (CurrentModel%WaterPropertiesID,                            &
                                         CurrentModel%SubMPIID(iSub),                                       &
                                         CurrentModel%SubmodelLink(iSub)%Window,                            &
                                         InitialField = .false.,                                            &
                                         PropIDNumber = CurrentModel%SubmodelLink(iSub)%PropertyIDNumbers(iProp),   &
                                         STAT         = STAT_CALL)

                enddo

            end if
            
        enddo

#else _USE_MPI

    !Disable unused variable warnings
    CurrentModel = CurrentModel

#endif _USE_MPI
    end subroutine SendInformationMPI 

    !--------------------------------------------------------------------------

    subroutine UpdateSubModelValues (CurrentModel)

        !Arguments-------------------------------------------------------------
        type (T_MohidWater), pointer                :: CurrentModel

#ifdef _USE_MPI
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iProp

        if (CurrentModel%FatherLink%Hydro) then
            call UpdateHydroMPI (CurrentModel%HydrodynamicID,             &
                                 InitialField = .false.,                  &
                                 STAT = STAT_CALL)
        endif

        if (CurrentModel%FatherLink%Water) then
            
            do iProp = 1, CurrentModel%FatherLink%nProp

                call UpdateWaterMPI (CurrentModel%WaterPropertiesID,           &
                                     InitialField = .false.,                   &
                                     PropIDNumber = CurrentModel%FatherLink%PropertyIDNumbers(iProp), &
                                     STAT = STAT_CALL)

            end do

        end if
#else _USE_MPI

    !Disable unused variable warnings
    CurrentModel = CurrentModel

#endif _USE_MPI

    end subroutine UpdateSubModelValues


#ifdef _OPENMI_

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::Initialize
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_INITIALIZE"::Initialize
    !DEC$ ENDIF
    logical function Initialize(workingDirectory)
                     
        !Arguments-------------------------------------------------------------
        character(*)                                :: workingDirectory
        
        !Local-----------------------------------------------------------------
        
        FilesName = workingDirectory
        
        !Disable unused variable warnings
        RunInParallel = .false.
        myMPI_Processor = 'My Processor'
    
        call ConstructMohidWater

        Initialize = .true.

        return
    
    end function Initialize
    
    !--------------------------------------------------------------------------

    !Perform a single time step
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::PerformTimeStep
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_PERFORMTIMESTEP"::PerformTimeStep
    !DEC$ ENDIF
    logical function PerformTimeStep()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        logical                                     :: dummy
        
        dummy = DoOneTimeStep()
        PerformTimeStep = .true.

    end function PerformTimeStep
    
    !--------------------------------------------------------------------------

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::Finish
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_FINISH"::Finish
    !DEC$ ENDIF
    logical function Finish()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        call KillMohidWater()
        Finish = .true.

    end function Finish

    !--------------------------------------------------------------------------

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::Dispose
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_DISPOSE"::Dispose
    !DEC$ ENDIF
    !The dispose function does not do anything. All Clean up is done by de Finish function
    logical function Dispose()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        Dispose = .true.

    end function Dispose
    
    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetModelID
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETMODELID"::GetModelID
    !DEC$ ENDIF
    logical function GetModelID(id)
    
        !Arguments-------------------------------------------------------------
        character(*)                                :: id       
    
        id = FirstModel%ModelName
        GetModelID = .true.
        return
    
    end function GetModelID

    !--------------------------------------------------------------------------

    !Test Function - Runs the whole model
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::RunSimulation
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_RUNSIMULATION"::RunSimulation
    !DEC$ ENDIF
    !Test method to run the whole simulation once
    logical function RunSimulation()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        call ModifyMohidWater
        
        RunSimulation = .true.
    
    end function RunSimulation

    !--------------------------------------------------------------------------

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetNumberOfMessages
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETNUMBEROFMESSAGES"::GetNumberOfMessages
    !DEC$ ENDIF
    !Return the number of Error Messages
    integer function GetNumberOfMessages()
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------

        GetNumberOfMessages = NumberOfErrorMessages
        
        return
    
    end function GetNumberOfMessages

    !--------------------------------------------------------------------------
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetMessage
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETMESSAGE"::GetMessage
    !DEC$ ENDIF
    logical function GetMessage(Number, Message)

        !Arguments-------------------------------------------------------------
        integer                                     :: Number
        character(len=*)                            :: Message        
        !Local-----------------------------------------------------------------


        if(Number .ge. 1 .and. Number .le. MaxErrorMessages)then
            Message=ErrorMessagesStack(Number)
            GetMessage=.true.
        else
            Message=' '
            GetMessage=.false.
        endif

      end function GetMessage
      
    
    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStartInstant
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTARTINSTANT"::GetStartInstant
    !DEC$ ENDIF
    logical function GetStartInstant(Instant)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Instant        

        !Local-----------------------------------------------------------------

        Instant = ConvertTimeToString(GlobalBeginTime)
        
        GetStartInstant = .true.

      end function GetStartInstant      

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetStopInstant
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETSTOPINSTANT"::GetStopInstant
    !DEC$ ENDIF
    logical function GetStopInstant(Instant)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Instant        

        !Local-----------------------------------------------------------------

        Instant = ConvertTimeToString(GlobalEndTime)
        
        GetStopInstant = .true.

      end function GetStopInstant      

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetCurrentInstant
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETCURRENTINSTANT"::GetCurrentInstant
    !DEC$ ENDIF
    logical function GetCurrentInstant(Instant)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Instant        

        !Local-----------------------------------------------------------------

        Instant = ConvertTimeToString(GlobalCurrentTime)
        
        GetCurrentInstant = .true.

      end function GetCurrentInstant      


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetCurrentTimeStep
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETCURRENTTIMESTEP"::GetCurrentTimeStep
    !DEC$ ENDIF
    real(8) function GetCurrentTimeStep()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: DTmin, DTmax
        
        
        DTmin   = - FillValueReal
        DTmax   =   FillValueReal
        call SearchMinMaxTimeStep (DTmin, DTmax)

        GetCurrentTimeStep = dble(DTmin)
        

      end function GetCurrentTimeStep      


    
#endif




#ifdef _OPENMI_
end module MohidWater
#else
end program MohidWater
#endif


!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

