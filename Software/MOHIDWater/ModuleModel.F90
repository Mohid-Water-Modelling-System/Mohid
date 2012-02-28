!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Model
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which coordinates the evolution of one Model (Main or SubModel)
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

Module ModuleModel

    use ModuleGlobalData
    use ModuleFunctions,                only : ReadTimeKeywords, GetDataOnlineString
    use ModuleEnterData           
    use ModuleTime
    use ModuleHorizontalGrid,           only : ConstructHorizontalGrid, KillHorizontalGrid, &
                                               GetHorizontalGridSize
    use ModuleGridData,                 only : ConstructGridData, KillGridData
    use ModuleHorizontalMap,            only : ConstructHorizontalMap,  KillHorizontalMap,  &
                                               GetWaterPoints2D, UnGetHorizontalMap
    use ModuleGeometry,                 only : ConstructGeometry,       KillGeometry,    &
                                               GetGeometrySize
    use ModuleMap,                      only : ConstructMap,            KillMap,         &
                                               UpdateComputeFaces3D
    use ModuleTurbulence,               only : ConstructTurbulence, Turbulence,          &
                                               KillTurbulence, GetTurbulenceOptions
#ifndef _USE_SEQASSIMILATION
    use ModuleHydrodynamic,             only : StartHydrodynamic, GetWaterLevel,         &
                                               GetHorizontalVelocity, GetChezy,          &
                                               GetVerticalVelocity,                      &
                                               Modify_Hydrodynamic, UnGetHydrodynamic,   &
                                               KillHydrodynamic          
#else
    use ModuleHydrodynamic,             only : StartHydrodynamic, GetWaterLevel,         &
                                               GetHorizontalVelocity, GetChezy,          &
                                               GetVerticalVelocity,                      &
                                               Modify_Hydrodynamic, UnGetHydrodynamic,   &
                                               KillHydrodynamic,                         &
                                               GetHydroSeqAssimilation,                  &
                                               SetHydroVirtualRun          
#endif _USE_SEQASSIMILATION
#ifndef _USE_SEQASSIMILATION
    use ModuleWaterProperties,          only : Construct_WaterProperties,                &
                                               GetDensity, GetConcentration, GetSigma,   &
                                               UnGetWaterProperties,                     &
                                               WaterProperties_Evolution,                &
                                               KillWaterProperties
#else
    use ModuleWaterProperties,          only : Construct_WaterProperties,                &
                                               GetDensity, GetConcentration, GetSigma,   &
                                               UnGetWaterProperties,                     &
                                               WaterProperties_Evolution,                &
                                               KillWaterProperties,                      &
                                               GetWaterSeqAssimilation,                  &
                                               SetWaterPropVirtualRun
#endif

#ifndef _LAGRANGIAN_
#ifdef _LAGRANGIAN_GLOBAL_  
    use ModuleLagrangianGlobal,         only : AllocateLagrangianGlobal, ConstructLagrangianGlobal, &
                                               ModifyLagrangianGlobal, DeallocateLagrangianGlobal,  &
                                               KillLagrangianGlobal
#else
    use ModuleLagrangian,               only : ConstructLagrangian, ModifyLagrangian,    &
                                               KillLagrangian
#endif
#endif


#ifndef _SEDIMENT_
    use ModuleConsolidation,            only : ConstructConsolidation, KillConsolidation,&
                                               ModifyConsolidation

    use ModuleSedimentProperties,       only : Construct_SedimentProperties,             &
                                               SedimentProperties_Evolution,             &
                                               KillSedimentProperties
#endif

    use ModuleInterfaceSedimentWater,   only : StartInterfaceSedimentWater,              &
                                               ModifyInterfaceSedimentWater,             &
                                               KillInterfaceSedimentWater
#ifndef _AIR_
    use ModuleInterfaceWaterAir,        only : StartInterfaceWaterAir,                   &
                                               ModifyInterfaceWaterAir,                  &
                                               KillInterfaceWaterAir

    use ModuleAtmosphere,               only : StartAtmosphere,                          &
                                               ModifyAtmosphere,                         &
                                               KillAtmosphere
#endif

#ifndef _WAVES_
    use ModuleWaves,                    only : StartWaves, ModifyWaves, KillWaves,       &
                                               ComputeWaveParameters, GetWavesOptions
#endif

#ifdef _USE_SEQASSIMILATION
    use ModuleSequentialAssimilation,   only : StartSequentialAssimilation,              &
                                               SetModelInitialState, GetModelResult,     &
                                               CovarianceCalculationSEEK,                &
                                               GetSeqAssimilationTime,                   &
                                               ModifySequentialAssimilation,             &
                                               KillSequentialAssimilation,               &
                                               GetSeqAssimilationOptions
#endif _USE_SEQASSIMILATION
    use ModuleStopWatch,            only: StartWatch, StopWatch
    
#ifdef _ENABLE_CUDA
    ! JPW, CUDA support
    use ModuleCuda
#endif

    !$ use omp_lib


    implicit none

    private

    !Subroutines--------------------------------------------

    !Constructor
    public  :: ConstructModel

    !Modifier
    public  :: UpdateTimeAndMapping
    public  :: RunModel

    !Destructor
    public  :: KillModel

    !Selector
    public  :: GetModelTimeLimits 
    public  :: GetModelTimeStep
    public  :: GetModelCurrentTime
    public  :: GetModelInstanceIDs
    public  :: GetSubModelWindow
    public  :: GetSubModelWindowON

#ifdef  OVERLAP
    public  :: GetModelOverlap
    public  :: GetModelOverlapInfo
#endif OVERLAP

    !Management
    private :: AllocateInstance
    private :: DeallocateInstance
    private :: Read_Lock
    private :: Read_Unlock
    private :: Ready


    integer, parameter                          :: ReadQueryString_         = 1
    integer, parameter                          :: ReadCommandLine_         = 2




    !Type----------------------------------------------------------------------
    
  
    type T_ExternalVar
        real, dimension(:,:,:), pointer         :: Density, SigmaDens       => null()
        real, dimension(:,:,:), pointer         :: Salinity                 => null()
        real, dimension(:,:,:), pointer         :: Temperature              => null()
        real, dimension(:,:,:), pointer         :: VelocityX                => null()
        real, dimension(:,:,:), pointer         :: VelocityY                => null()
        real, dimension(:,:,:), pointer         :: VelocityZ                => null()
        real, dimension(:,:  ), pointer         :: Chezy, WaterLevel        => null()
        logical                                 :: NeedsTempSalinity
    end type T_ExternalVar

    !Groups modules of the WaterColumn
    type T_Column
        integer                                 :: ObjBathymetry            = 0
        integer                                 :: ObjHorizontalMap         = 0
        integer                                 :: ObjGeometry              = 0
        integer                                 :: ObjMap                   = 0
    end type T_Column

#ifdef OVERLAP
    type     T_Overlap
        logical                                 :: Yes                      = .false.
        integer                                 :: OverlapModelID           = 0
        character(len=PathLength)               :: FileName
        integer, dimension(:,:), pointer        :: Cells                    => null()
    end type T_Overlap
#endif OVERLAP

    !public :: T_Model
    type T_Model 

        integer                                 :: InstanceID

        !integer, dimension(:,:), pointer        :: LagInstance

        !character(StringLength), dimension(:), pointer :: ModelNames

        character(StringLength)                   :: ModelName

        integer                                 :: NumberOfModels
        integer                                 :: ObjLagrangian
!        integer                                 :: ObjLagrangianX
        integer                                 :: ObjLagrangianGlobal

        type (T_ExternalVar)                    :: ExternalVar

        type (T_Time)                           :: InitialSystemTime 
        type (T_Time)                           :: CurrentTime
        type (T_Time)                           :: BeginTime
        type (T_Time)                           :: EndTime
        integer                                 :: Iteration
        real                                    :: DT, MaxDT
        logical                                 :: VariableDT
        real                                    :: InitialDT
        real                                    :: GmtReference
        real                                    :: DTPredictionInterval

        logical                                 :: RunSediments
        logical                                 :: RunLagrangian
        logical                                 :: RunWaves
        logical                                 :: NoIsolatedCells
        integer                                 :: OnLineType

        type (T_Size2D)                         :: SubModelWindow
        logical                                 :: SubModelWindowON
        
        logical                                 :: BackTracking

#ifdef OVERLAP
        type (T_Overlap)                        :: Overlap
#endif OVERLAP
        
#ifdef _USE_SEQASSIMILATION
        !Sequential data assimilation variables:
        logical                                 :: RunSeqAssimilation
        logical                                 :: StateCovEvolution
        integer                                 :: StateCovRank
        type (T_Time)                           :: SeqAssimilationTime
#endif _USE_SEQASSIMILATION

        !$ logical                              :: FatherModelFlag = .false.

        !Instance of other Modules
        integer                                 :: ObjTime                  = 0
        integer                                 :: ObjHorizontalGrid        = 0
        type (T_Column)                         :: Water
        type (T_Column)                         :: Sediment
        integer                                 :: ObjAssimilation          = 0
        integer                                 :: ObjDischarges            = 0
        integer                                 :: ObjTurbGOTM              = 0
        integer                                 :: ObjTurbulence            = 0
        integer                                 :: ObjHydrodynamic          = 0
        integer                                 :: ObjWaterProperties       = 0
        integer                                 :: ObjAtmosphere            = 0
        integer                                 :: ObjWaves                 = 0
        integer                                 :: ObjInterfaceWaterAir     = 0
        integer                                 :: ObjInterfaceSedimentWater= 0
!        integer                                 :: ObjInterfaceSedAir      = 0
        integer                                 :: ObjSedimentProperties    = 0
        integer                                 :: ObjConsolidation         = 0
#ifdef _USE_SEQASSIMILATION
        integer                                 :: ObjSeqAssimilation       = 0
#endif _USE_SEQASSIMILATION

#ifdef _ENABLE_CUDA
        ! JPW: CUDA support
        integer                                 :: ObjCuda                  = 0
#endif _ENABLE_CUDA

        !Linked list of Instances
        type(T_Model), pointer                  :: Next

    End Type T_Model

    !Global Module Variables
    type (T_Model), pointer                     :: FirstModel
    type (T_Model), pointer                     :: Me

    Contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !In this subroutine is read the model structure (number of parallel models and their sub-models)
    ! and the necessary memory is allocate and store in a structure of pointer lists

    subroutine ConstructModel (LagInstance, ModelNames, NumberOfModels, ObjLagrangianGlobal, ModelID, InitialSystemTime, STAT)

        !Arguments-------------------------------------------------------------
        integer         , dimension(:,:), pointer   :: LagInstance
        character(len=*), dimension(:  ), pointer   :: ModelNames
        integer                                     :: NumberOfModels, ObjLagrangianGlobal, ModelID
        type (T_Time)                               :: InitialSystemTime
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer, allocatable, dimension(:)          :: AuxInt4
        integer                                     :: STAT_
        integer                                     :: STAT_CALL, ready_
        integer                                     :: FromFile, flag
        character(PathLength)                       :: BathymetryFile, DataFile
        character(PathLength), save                 :: LagNomfich

        !character(StringLength), dimension(:), pointer :: AuxModelNames

#ifndef _SEDIMENT_
        character(PathLength)                       :: SedimentFile, SedGeometryFile
#endif
        integer                                     :: ObjEnterData = 0
#ifndef _AIR_
        integer, dimension(:, :), pointer           :: WaterPoints2D
#endif
#ifndef _WAVES_
        logical                                     :: WaveParametersON
#endif

#ifdef _USE_SEQASSIMILATION
        logical                                     :: HydroSeqAssim = .false.
        logical                                     :: WaterSeqAssim = .false.
#endif _USE_SEQASSIMILATION
        !$ integer                                     :: openmp_num_threads
        
#ifdef _ENABLE_CUDA
        type(T_Size3D)                              :: GeometrySize
#endif _ENABLE_CUDA

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mModel_)) then
            nullify (FirstModel)
            call RegisterModule (mModel_) 
        endif

        call Ready(ModelID, ready_)    

if0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance 

            Me%NumberOfModels   = NumberOfModels            
            !Me%ModelNames       => ModelNames
            !Me%LagInstance      => LagInstance

            !Stores name
            Me%ModelName        = trim(ModelNames(Me%InstanceID))

#ifndef _OUTPUT_OFF_
            write(*, *)"-------------------------- MODEL -------------------------"
            write(*, *)
            write(*, *)"Constructing      : ", trim(Me%ModelName)
            write(*, *)"ID                : ", Me%InstanceID
#endif

            Me%InitialSystemTime = InitialSystemTime

            !Gets the name of the data file
            call ReadFileName('IN_MODEL', DataFile, "Compute Time Data File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR10'

            !Constructs EnterData
            call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR20'


            !Read the Time related keywords
            call GetExtractType  (FromFile = FromFile)

            call ReadTimeKeyWords(ObjEnterData, FromFile, Me%BeginTime, Me%EndTime,      &
                                  Me%DT, Me%VariableDT, "Model", Me%MaxDT,               &
                                  Me%GmtReference, Me%DTPredictionInterval)
            Me%InitialDT = Me%DT
            Me%Iteration = 0

#ifdef      _ONLINE_

            !Check mohid reads the online data
            call GetData         (Me%OnLineType, ObjEnterData, flag,                     &
                                  SearchType    =  FromFile,                             &
                                  keyword       = 'ONLINE_TYPE',                         &
                                  default       = ReadCommandLine_,                      &
                                  ClientModule  = 'ModuleModel',                         &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR30'

            call ReadStringOnline
#endif

            !Use Sediment columns
            call GetData         (Me%RunSediments, ObjEnterData, flag,                   &
                                  SearchType    =  FromFile,                             &
                                  keyword       = 'SEDIMENT_COLUMN',                     &
                                  default       = .false.,                               &
                                  ClientModule  = 'ModuleModel',                         &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR40'


            !Use Lagrangian
            call GetData         (Me%RunLagrangian, ObjEnterData, flag,                  &
                                  SearchType    =  FromFile,                             &
                                  keyword       = 'LAGRANGIAN',                          &
                                  default       = .false.,                               &
                                  ClientModule  = 'ModuleModel',                         &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR50'

            !Use Waves
            call GetData         (Me%RunWaves, ObjEnterData, flag,                       &
                                  SearchType    =  FromFile,                             &
                                  keyword       = 'WAVES',                               &
                                  default       = .false.,                               &
                                  ClientModule  = 'ModuleModel',                         &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR60'

#ifdef _USE_SEQASSIMILATION
            !Use SequentialAssimilation (Ang)
            call GetData         (Me%RunSeqAssimilation, ObjEnterData, flag,             &
                                  SearchType    =  FromFile,                             &
                                  keyword       = 'SEQUENTIAL_ASSIMILATION',             &
                                  default       = .false.,                               &
                                  ClientModule  = 'ModuleModel',                         &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR70'
#endif _USE_SEQASSIMILATION

            call GetData         (Me%NoIsolatedCells, ObjEnterData, flag,                &
                                  SearchType    =  FromFile,                             &
                                  keyword       = 'NO_ISOLATED_CELLS',                   &
                                  default       = .true.,                                &
                                  ClientModule  = 'ModuleModel',                         &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR80'



            call GetData         (Me%BackTracking, ObjEnterData, flag,                  &
                                  SearchType    =  FromFile,                            &
                                  keyword       = 'BACKTRACKING',                       &
                                  default       = .false.,                              &
                                  ClientModule  = 'ModuleModel',                        &
                                  STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR90'


            !Start Module Time
            Me%CurrentTime = Me%BeginTime
            call StartComputeTime(Me%ObjTime, Me%InitialSystemTime, Me%BeginTime, Me%EndTime, Me%DT,     &
                                  Me%VariableDT, Me%MaxDT, Me%GmtReference, Me%BackTracking, STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR100'

            !Gets the file name of the Bathymetry
            call ReadFileName('IN_BATIM', BathymetryFile, "Bathymetry File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR110'

            !Horizontal Grid
            call ConstructHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,        &
                                         DataFile         = BathymetryFile,              &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR120'
          
            allocate(AuxInt4(1:4))

            AuxInt4(1) = FillValueInt

            call GetData         (AuxInt4, ObjEnterData, flag,                          &
                                  SearchType    =  FromFile,                            &
                                  keyword       = 'SUBMODEL_WINDOW',                    &
                                  ClientModule  = 'ModuleModel',                        &
                                  STAT          = STAT_CALL)

            if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= KEYWORD_NOT_FOUND_ERR_)        &
                stop 'ConstructModel - ModuleModel - ERR47'

            if (flag > 0) then
                if (flag /= 4) stop 'ConstructModel - ModuleModel - ERR130'
                Me%SubModelWindow%ILB = AuxInt4(1)
                Me%SubModelWindow%JLB = AuxInt4(2)
                Me%SubModelWindow%IUB = AuxInt4(3)
                Me%SubModelWindow%JUB = AuxInt4(4)
                Me%SubModelWindowON = .true.
            else
                Me%SubModelWindowON = .false.
            endif

            deallocate(AuxInt4)

#ifdef OVERLAP
            
            call GetData(Me%Overlap%Yes,                                                    & 
                         ObjEnterData, flag,                                           &
                         SearchType    =  FromFile,                                         &
                         keyword       = 'OVERLAP',                                         &
                         ClientModule  = 'ModuleModel',                                     &
                         STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR140'

            if(Me%Overlap%Yes)then
                call ConstructOverlap(ObjEnterData)
            endif

#endif OVERLAP

            !griflet: Adding a new keyword in model.dat to define the number of threads to use.
            !Reads OPENMP_NUM_THREADS (number of threads to use with openmp)
            !$ call GetData(openmp_num_threads, ObjEnterData, flag, keyword = 'OPENMP_NUM_THREADS',  &
            !$         SearchType   = FromFile,                                                      &
            !$         ClientModule = 'ModuleModel',                                                 &
            !$         default      = 0,                                                             &
            !$         STAT         = STAT_CALL)
            !$ if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR94'            
            !$ if ( .not. FirstModel%FatherModelFlag ) then
            !$    FirstModel%FatherModelFlag = .true.
            !$    write(*,*)
            !$    write(*,*)"OPENMP: Max number of threads available is ", omp_get_max_threads()
            !$    if ( openmp_num_threads .gt. 0 ) then
            !$       write(*,*)"OPENMP: Number of threads requested is ", openmp_num_threads
            !$       if (openmp_num_threads .gt. omp_get_max_threads()) then
            !$        openmp_num_threads = omp_get_max_threads()
            !$        write(*,*)"<Compilation Options Warning>"
            !$       endif
            !$       call omp_set_num_threads(openmp_num_threads)
            !$       write(*,*)"OPENMP: Number of threads implemented is ", openmp_num_threads
            !$    else
            !$       write(*,*)"OPENMP: Using the minimum value between max number of threads available"
            !$       write(*,*)"OPENMP: and OMP_NUM_THREADS environment variable value."
            !$    endif
            !$ else
            !$    if ( openmp_num_threads .gt. 0 ) then
            !$       write(*,*) "OPENMP: WARNING, OPENMP_NUM_THREADS should be defined in the father model only!"
            !$    endif
            !$ endif

            call KillEnterData    (ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR150'

            !Horizontal Grid Data - Water Column (Bathymetry)
            call ConstructGridData      (GridDataID       = Me%Water%ObjBathymetry,      &
                                         HorizontalGridID = Me%ObjHorizontalGrid,        &
                                         TimeID           = Me%ObjTime,                  &
                                         FileName         = BathymetryFile,              &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR160'

            !Horizontal Map
            call ConstructHorizontalMap (HorizontalMapID  = Me%Water%ObjHorizontalMap,   &
                                         GridDataID       = Me%Water%ObjBathymetry,      &
                                         HorizontalGridID = Me%ObjHorizontalGrid,        &
                                         ActualTime       = Me%CurrentTime,              &
                                         STAT             = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR170'

            !Geometry - Water Column
            call ConstructGeometry      (GeometryID       = Me%Water%ObjGeometry,        &
                                         GridDataID       = Me%Water%ObjBathymetry,      &
                                         HorizontalGridID = Me%ObjHorizontalGrid,        &
                                         HorizontalMapID  = Me%Water%ObjHorizontalMap,   &
                                         ActualTime       = Me%CurrentTime,              &
                                         STAT             = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR180'

            !Map - Water Column            
            if (Me%NoIsolatedCells) then
                call ConstructMap           (Map_ID           = Me%Water%ObjMap,          &
                                             GeometryID       = Me%Water%ObjGeometry,     &
                                             HorizontalMapID  = Me%Water%ObjHorizontalMap,&
                                             TimeID           = Me%ObjTime,               &
                                             GridDataID       = Me%Water%ObjBathymetry,   &
                                             HorizontalGridID = Me%ObjHorizontalGrid,     &
                                             STAT             = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR190'
            else

                call ConstructMap           (Map_ID           = Me%Water%ObjMap,          &
                                             GeometryID       = Me%Water%ObjGeometry,     &
                                             HorizontalMapID  = Me%Water%ObjHorizontalMap,&
                                             TimeID           = Me%ObjTime,               &
                                             HorizontalGridID = Me%ObjHorizontalGrid,     &
                                             STAT             = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR200'

            endif
            
#ifdef _ENABLE_CUDA
        ! Construct a ModuleCuda instance.
        call ConstructCuda(Me%ObjCuda, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR136'

        call GetGeometrySize(Me%Water%ObjGeometry, Size = GeometrySize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR137'
        
        ! JPW: Initialize a C++ Thomas instance for CUDA. 
        ! This will allocate device memory for all Thomas variables (D, E, F, TI, Res)
        ! Do this seperately from ConstructCuda, because later on ModuleCuda might be used for things other than Thomas algorithm
        call InitializeThomas(Me%ObjCuda, GeometrySize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR210'
#endif _ENABLE_CUDA

#ifndef _WAVES_
            if(Me%RunWaves)then

                call StartWaves(ModelName           = trim(Me%ModelName),               &  
                                WavesID             = Me%ObjWaves,                      &
                                TimeID              = Me%ObjTime,                       &
                                HorizontalMapID     = Me%Water%ObjHorizontalMap,        &
                                HorizontalGridID    = Me%ObjHorizontalGrid,             &
                                GridDataID          = Me%Water%ObjBathymetry,           &
                                GeometryID          = Me%Water%ObjGeometry,             &
                                STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop "Sub. ConstructModel - ModuleModel - ERR220"

            end if
#else
            write(*,*)
            write(*,*)"<Compilation Options Warning>"
            write(*,*)"This executable was not compiled with the Waves module."
            if (Me%RunWaves) then
                stop 'ConstructModel - ModuleModel - ERR230'
            endif
            write(*,*)"<Compilation Options Warning>"
            write(*,*)

#endif

            !Turbulence
            call ConstructTurbulence    (ModelName        = trim(Me%ModelName),         &  
                                         TurbulenceID     = Me%ObjTurbulence,           &
                                         TurbGOTMID       = Me%ObjTurbGOTM,             &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         GeometryID       = Me%Water%ObjGeometry,       &
                                         MapID            = Me%Water%ObjMap,            &
                                         GridDataID       = Me%Water%ObjBathymetry,     &
                                         HorizontalMapID  = Me%Water%ObjHorizontalMap,  &
                                         TimeID           = Me%ObjTime,                 &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR240'

            call GetTurbulenceOptions   (TurbulenceID     = Me%ObjTurbulence,            &
                                         NeedsTempSalinity= Me%ExternalVar%NeedsTempSalinity, &
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR250'

            !Constructs hydrodynamic
            call StartHydrodynamic      (ModelName        = trim(Me%ModelName),         &  
                                         HydrodynamicID   = Me%ObjHydrodynamic,         &
                                         GridDataID       = Me%Water%ObjBathymetry,     &
                                         HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         GeometryID       = Me%Water%ObjGeometry,       &
                                         HorizontalMapID  = Me%Water%ObjHorizontalMap,  &
                                         MapID            = Me%Water%ObjMap,            &
                                         AssimilationID   = Me%ObjAssimilation,         &
                                         TimeID           = Me%ObjTime,                 &
                                         TurbulenceID     = Me%ObjTurbulence,           &
                                         DischargesID     = Me%ObjDischarges,           &
                                         WavesID          = Me%ObjWaves,                &
#ifdef _ENABLE_CUDA
                                         CudaID           = Me%ObjCuda,                 &
#endif
                                         STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR260'

#ifndef _WAVES_
            !Need the watercolumn thickness that is compute the first time in the StartHydrodynamic
            if(Me%RunWaves)then
                call GetWavesOptions(Me%ObjWaves, WaveParametersON = WaveParametersON, STAT  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR270'
                
                if (WaveParametersON) then
                    call ComputeWaveParameters(Me%ObjWaves, STAT  = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR280'
                endif
            endif
#endif

            !Constructs waterProperties
            call Construct_WaterProperties(ModelName        = trim(Me%ModelName),       &  
                                           WaterPropertiesID= Me%ObjWaterProperties,    &
                                           GridDataID       = Me%Water%ObjBathymetry,   &
                                           HorizontalGridID = Me%ObjHorizontalGrid,     &
                                           GeometryID       = Me%Water%ObjGeometry,     &
                                           HorizontalMapID  = Me%Water%ObjHorizontalMap,&
                                           MapID            = Me%Water%ObjMap,          &
                                           HydrodynamicID   = Me%ObjHydrodynamic,       &
                                           TimeID           = Me%ObjTime,               &     
                                           TurbulenceID     = Me%ObjTurbulence,         &
                                           AssimilationID   = Me%ObjAssimilation,       &
                                           DischargesID     = Me%ObjDischarges,         &
#ifdef _ENABLE_CUDA
                                           CudaID           = Me%ObjCuda,                 &
#endif _ENABLE_CUDA
                                           STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "Sub. ConstructModel - ModuleModel - ERR290"

#ifndef _LAGRANGIAN_

il:         if (Me%RunLagrangian) then
            
#ifdef  _LAGRANGIAN_GLOBAL_       
                                  
                LagInstance(1,  Me%InstanceID) = Me%ObjTime                                    
                LagInstance(2,  Me%InstanceID) = Me%Water%ObjBathymetry
                LagInstance(3,  Me%InstanceID) = Me%ObjHorizontalGrid
                LagInstance(4,  Me%InstanceID) = Me%Water%ObjHorizontalMap
                LagInstance(5,  Me%InstanceID) = Me%Water%ObjGeometry
                LagInstance(6,  Me%InstanceID) = Me%Water%ObjMap
                LagInstance(7,  Me%InstanceID) = Me%ObjAssimilation
                LagInstance(8,  Me%InstanceID) = Me%ObjHydrodynamic
                LagInstance(9,  Me%InstanceID) = Me%ObjTurbulence
                LagInstance(10, Me%InstanceID) = Me%ObjWaves
                LagInstance(11, Me%InstanceID) = Me%ObjWaterProperties

                if (Me%InstanceID == 1) then
                    LagNomfich = FilesName

                    call AllocateLagrangianGlobal(LagrangianID = ObjLagrangianGlobal,   &
                                                  STAT         = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop "Sub. ConstructModel - ModuleModel - ERR300"
                    
                endif
                
                Me%ObjLagrangianGlobal = ObjLagrangianGlobal                

                if (Me%InstanceID == Me%NumberOfModels) then

                    call ConstructLagrangianGlobal(LagrangianID = Me%ObjLagrangianGlobal,&
                                                   Nmodels      = Me%NumberOfModels,    &
                                                   ModelNames   = ModelNames,           &
                                                   FileNomfich  = LagNomfich,           &
                                                   LagInstance  = LagInstance,          &             
                                                   STAT         = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_) stop "Sub. ConstructModel - ModuleModel - ERR310"

                endif

                !Me%ObjLagrangianX = Me%ObjLagrangianGlobal


#else
                 
                call ConstructLagrangian(Me%ObjLagrangian,                              &
                                         Me%ObjTime,                                    &
                                         Me%Water%ObjBathymetry,                        &
                                         Me%ObjHorizontalGrid,                          &
                                         Me%Water%ObjHorizontalMap,                     &
                                         Me%Water%ObjGeometry,                          &
                                         Me%Water%ObjMap,                               &
                                         Me%ObjAssimilation,                            &
                                         Me%ObjHydrodynamic,                            &
                                         Me%ObjTurbulence,                              &
                                         Me%ObjWaves,                                   &
                                         Me%ObjWaterProperties,                         &
                                         STAT )
                if (STAT_CALL /= SUCCESS_) stop "Sub. ConstructModel - ModuleModel - ERR320"

                !Removed warning for unused variable
                ObjLagrangianGlobal =  null_int
                LagNomfich          =  null_str
                nullify(LagInstance)
                !Me%ObjLagrangianX = Me%ObjLagrangian
    
#endif                            
            end if il
#else
            write(*,*)
            write(*,*)"<Compilation Options Warning>"
            write(*,*)"This executable was not compiled with the Lagrangian module."
            if (Me%RunLagrangian) then
                stop 'ConstructModel - ModuleModel - ERR330'
            endif
            write(*,*)"<Compilation Options Warning>"
            write(*,*)

#endif

#ifndef _SEDIMENT_

            if (Me%RunSediments) then

                !Gets the file name of the Bathymetry
                call ReadFileName('IN_SEDIMENT', SedimentFile, "Sediment File", STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR340'

                !Gets the file name of the Bathymetry
                call ReadFileName('SED_GEOM', SedGeometryFile, "Sediment Geometry File", STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR350'

                !Horizontal Grid Data - Sediment Column (Bathymetry)
                call ConstructGridData      (GridDataID       = Me%Sediment%ObjBathymetry,   &
                                             HorizontalGridID = Me%ObjHorizontalGrid,        &
                                             FileName         = SedimentFile,                &
                                             STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR360'

                !Horizontal Map
                call ConstructHorizontalMap (HorizontalMapID  = Me%Sediment%ObjHorizontalMap,&
                                             GridDataID       = Me%Sediment%ObjBathymetry,   &
                                             HorizontalGridID = Me%ObjHorizontalGrid,        &
                                             ActualTime       = Me%CurrentTime,              &
                                             STAT             = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR370'

                !Geometry - Sediment Column
                call ConstructGeometry      (GeometryID       = Me%Sediment%ObjGeometry,     &
                                             GridDataID       = Me%Sediment%ObjBathymetry,   &
                                             HorizontalGridID = Me%ObjHorizontalGrid,        &
                                             HorizontalMapID  = Me%Sediment%ObjHorizontalMap,&
                                             ActualTime       = Me%CurrentTime,              &
                                             NewDomain        = SedGeometryFile,             &
                                             STAT             = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR380'

                !Map - Sediment Column            
                call ConstructMap           (Map_ID           = Me%Sediment%ObjMap,          &
                                             GeometryID       = Me%Sediment%ObjGeometry,     &
                                             HorizontalMapID  = Me%Sediment%ObjHorizontalMap,&
                                             TimeID           = Me%ObjTime,                  &
                                             STAT             = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR390'


                call ConstructConsolidation(ConsolidationID     = Me%ObjConsolidation,          &
                                            TimeID              = Me%ObjTime,                   &
                                            GridDataID          = Me%Sediment%ObjBathymetry,    &
                                            HorizontalMapID     = Me%Sediment%ObjHorizontalMap, &
                                            HorizontalGridID    = Me%ObjHorizontalGrid,         &
                                            GeometryID          = Me%Sediment%ObjGeometry,      &
                                            MapID               = Me%Sediment%ObjMap,           &
                                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR400'


                !Constructs SedimentProperties
                call Construct_SedimentProperties(SedimentPropertiesID  = Me%ObjSedimentProperties,     &
                                                  TimeID                = Me%ObjTime,                   &     
                                                  HorizontalGridID      = Me%ObjHorizontalGrid,         &
                                                  GridDataID            = Me%Sediment%ObjBathymetry,    &
                                                  HorizontalMapID       = Me%Sediment%ObjHorizontalMap, &
                                                  MapID                 = Me%Sediment%ObjMap,           &
                                                  GeometryID            = Me%Sediment%ObjGeometry,      &
                                                  ConsolidationID       = Me%ObjConsolidation,          &
                                                  STAT                  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR410'

            endif
#else
            write(*,*)
            write(*,*)"<Compilation Options Warning>"
            write(*,*)"This executable was not compiled with the Sediment modules."
            if(Me%RunSediments)then
                stop 'ConstructModel - ModuleModel - ERR420'
            endif
            write(*,*)"<Compilation Options Warning>"
            write(*,*)

#endif

            call StartInterfaceSedimentWater(ModelName                   = trim(Me%ModelName),&  
                                             ObjInterfaceSedimentWaterID = Me%ObjInterfaceSedimentWater, &
                                             TimeID                      = Me%ObjTime,                   &     
                                             HorizontalGridID            = Me%ObjHorizontalGrid,         &
                                             WaterGridDataID             = Me%Water%ObjBathymetry,       &
                                             WaterHorizontalMapID        = Me%ObjHorizontalGrid,         &
                                             WaterMapID                  = Me%Water%ObjMap,              &
                                             WaterGeometryID             = Me%Water%ObjGeometry,         &
                                             SedimentGridDataID          = Me%Sediment%ObjBathymetry,    &
                                             SedimentHorizontalMapID     = Me%Sediment%ObjHorizontalMap, &
                                             SedimentMapID               = Me%Sediment%ObjMap,           &
                                             SedimentGeometryID          = Me%Sediment%ObjGeometry,      &
                                             HydrodynamicID              = Me%ObjHydrodynamic,           &
                                             TurbGOTMID                  = Me%ObjTurbGOTM,               &
                                             TurbulenceID                = Me%ObjTurbulence,             &
                                             WaterPropertiesID           = Me%ObjWaterProperties,        &
#ifdef  _LAGRANGIAN_GLOBAL_          
                                             LagrangianID                = Me%ObjLagrangianGlobal,       &
#else
                                             LagrangianID                = Me%ObjLagrangian,             &
#endif
                                             WavesID                     = Me%ObjWaves,                  &
                                             SedimentPropertiesID        = Me%ObjSedimentProperties,     &
                                             ConsolidationID             = Me%ObjConsolidation,          &
                                             DischargesID                = Me%ObjDischarges,             &
                                             RunsSediments               = Me%RunSediments,              &
                                             STAT                        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR430'

#ifndef _AIR_

            !Starts Atmosphere
            call GetWaterPoints2D (Me%Water%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR440'

            call StartAtmosphere(ModelName          = trim(Me%ModelName),&
                                 AtmosphereID       = Me%ObjAtmosphere,                 &
                                 TimeID             = Me%ObjTime,                       &
                                 GridDataID         = Me%Water%ObjBathymetry,           &
                                 HorizontalGridID   = Me%ObjHorizontalGrid,             &
                                 MappingPoints      = WaterPoints2D,                    &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR450'

            call UnGetHorizontalMap (Me%Water%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR460'

            call StartInterfaceWaterAir(ModelName                   = trim(Me%ModelName),&
                                        ObjInterfaceWaterAirID      = Me%ObjInterfaceWaterAir,&
                                        TimeID                      = Me%ObjTime,             &     
                                        HorizontalGridID            = Me%ObjHorizontalGrid,   &
                                        WaterGridDataID             = Me%Water%ObjBathymetry, &
                                        WaterHorizontalMapID        = Me%ObjHorizontalGrid,   &
                                        WaterMapID                  = Me%Water%ObjMap,        &
                                        WaterGeometryID             = Me%Water%ObjGeometry,   &
                                        HydrodynamicID              = Me%ObjHydrodynamic,     &
                                        TurbGOTMID                  = Me%ObjTurbGOTM,         &
                                        WavesID                     = Me%ObjWaves,            &  
                                        WaterPropertiesID           = Me%ObjWaterProperties,  &
                                        AtmosphereID                = Me%ObjAtmosphere,       &
#ifdef  _LAGRANGIAN_GLOBAL_          
                                        LagrangianID                = Me%ObjLagrangianGlobal, &
#else
                                        LagrangianID                = Me%ObjLagrangian,       &
#endif
                                        STAT                        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - 470'
#else
            write(*,*)
            write(*,*)"<Compilation Options Warning>"
            write(*,*)"This executable was not compiled with the Interface Water-Air module."
            write(*,*)"<Compilation Options Warning>"
            write(*,*)

#endif

#ifdef _USE_SEQASSIMILATION
            !Sequential data assimilation (Ang)
            if (Me%RunSeqAssimilation) then

                !Check if sequential data assimilation is commanded in modules
                !Hydrodynamic
                call GetHydroSeqAssimilation(Me%ObjHydrodynamic, HydroSeqAssim, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR480'

                !WaterProperties
                call GetWaterSeqAssimilation(Me%ObjWaterProperties, WaterSeqAssim, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR490'

                if (.not. HydroSeqAssim .and. .not. WaterSeqAssim) then
                    write(*,*)
                    write(*,*)"Sequential assimilation commanded in Model input file"
                    write(*,*)"but absent from Hydrodynamic and WaterProperties input files."
                    stop 'ConstructModel - ModuleModel - ERR500'
                endif

                call StartSequentialAssimilation(Me%ObjSeqAssimilation,                     &
                                         GridDataID         = Me%Water%ObjBathymetry,       &
                                         HorizontalGridID = Me%ObjHorizontalGrid,           &
                                         TimeID             = Me%ObjTime,                   &
                                         HydrodynamicID     = Me%ObjHydrodynamic,           &
                                         WaterPropertiesID  = Me%ObjWaterProperties,        &
                                         GeometryID         = Me%Water%ObjGeometry,         &
                                         MapID              = Me%Water%ObjMap,              &
                                         STAT               = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR510'

                !Get essential sequential assimilation variables
                !(Me%StateCovEvolution .and. Me%StateCovRank)
                call GetSeqAssimilationOptions (Me%ObjSeqAssimilation,                      &
                                                StateCovEvolution = Me%StateCovEvolution,   &
                                                StateCovRank      = Me%StateCovRank,        &
                                                STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR520'

                !Actualizes Me%NextSeqAssimilationTime
                call GetSeqAssimilationTime(Me%ObjSeqAssimilation,                          &
                                            Me%SeqAssimilationTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructModel - ModuleModel - ERR530'
            endif
#endif _USE_SEQASSIMILATION

            !nullify(Me%ModelNames )
            !nullify(Me%LagInstance)
            
            !Returns ID
            ModelID     = Me%InstanceID
            STAT_       = SUCCESS_

        else 

            stop 'ModuleModel - ConstructModel - ERR540' 

        end if if0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    End Subroutine ConstructModel

    !--------------------------------------------------------------------------

#ifdef OVERLAP

    subroutine ConstructOverlap(ObjModelEnterData)

        !Arguments-------------------------------------------------------------
        integer,intent(IN)                  :: ObjModelEnterData

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, flag, ClientNumber
        integer                             :: ObjEnterData                = 0
        logical                             :: BlockFound                  = .false.
        integer                             :: FirstLine, LastLine
        integer                             :: nOverlapCells, iCell, line
        integer, dimension(:),pointer       :: Aux
        !Begin-----------------------------------------------------------------
        
        call GetData(Me%Overlap%OverlapModelID,                                         & 
                     ObjModelEnterData, flag,                                           &
                     SearchType    =  FromFile,                                         &
                     keyword       = 'OVERLAP_MODEL_ID',                                &
                     ClientModule  = 'ModuleModel',                                     &
                     STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOverlap - ModuleModel - 05'

        call GetData(Me%Overlap%FileName,                                               & 
                     ObjModelEnterData, flag,                                           &
                     SearchType    =  FromFile,                                         &
                     keyword       = 'OVERLAP_FILE',                                    &
                     ClientModule  = 'ModuleModel',                                     &
                     STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOverlap - ModuleModel - 10'

        
        call ConstructEnterData(ObjEnterData, Me%Overlap%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOverlap - ModuleModel - 20'

        call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                         &
                                    '<begin_cells>', '<end_cells>', BlockFound,         &
                                    FirstLine = FirstLine, LastLine = LastLine,         &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOverlap - ModuleModel - 30'
            
        if(BlockFound)then

            nOverlapCells = LastLine - FirstLine - 1

            allocate(Me%Overlap%Cells(1:nOverlapCells, 1:4))
            allocate(Aux             (                 1:4))
            
            iCell = 1
            do line = FirstLine + 1, LastLine - 1

                call GetData(Aux, EnterDataID = ObjEnterData, flag = flag,              &
                             SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructOverlap - ModuleModel - 40'

                !if(Aux(1) < Me%WorkSize3D%ILB .or. Aux(1) > Me%WorkSize3D%IUB)then

                !endif

                !if(Aux(2) < Me%WorkSize3D%JLB .or. Aux(2) > Me%WorkSize3D%JUB)then

                !endif

                !if(Aux(3) < Me%WorkSize3D%KLB .or. Aux(3) > Me%WorkSize3D%KUB)then

                !endif

                if(flag .ne. 4)then
                    stop 'ConstructOverlap - ModuleModel - 50'
                end if

                Me%Overlap%Cells(iCell,1) = Aux(1)
                Me%Overlap%Cells(iCell,2) = Aux(2)
                Me%Overlap%Cells(iCell,3) = Aux(3)
                Me%Overlap%Cells(iCell,4) = Aux(4)

                iCell = iCell + 1

            enddo

            deallocate(Aux)

        else
            write(*,*)'<begin_cells>...<end_cells> block was not found in file:'
            write(*,*)trim(Me%Overlap%FileName)
            stop 'ConstructOverlap - ModuleModel - 60'
        end if

        call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOverlap - ModuleModel - 90'


    end subroutine

#endif OVERLAP

    !--------------------------------------------------------------------------
#ifdef _ONLINE_

    Subroutine ReadStringOnline

        
        !Local-----------------------------------------------------------------
        real,  dimension(:), pointer     :: TimeReal
        integer                          :: i, STAT_CALL
        !Begin-----------------------------------------------------------------

        if      (Me%OnLineType == ReadQueryString_) then
    
            CALL getenv('QUERY_STRING',OnlineString)
        
        else if (Me%OnLineType == ReadCommandLine_) then

!            OnlineString="%LAG_NORIGINS=1%EMISSION=1%LAG_XY=-3.7722_36.73037%LAG_START=2009_06_07_18_00_00%"
!            OnlineString=trim(OnlineString)//"WIND_COEF=0.00%FLOW=200.0%CONCENTRATION=1.e7%T90=43200.%"
!            OnlineString=trim(OnlineString)//"START=2009_06_07_04_00_00%END=2009_06_08_04_00_00%"

!            OnlineString="%LAG_NORIGINS=1%EMISSION=1%LAG_XY=-9.285285_38.69494%LAG_START=2010_09_27_10_00_00%"
!            OnlineString=trim(OnlineString)//"WIND_COEF=0.00%FLOW=1.0%CONCENTRATION=1.e7%T90=10800.%START=2010_09_27_10_00_00%END=2010_09_27_13_00_00%"


            
            call GETARG (1, OnlineString, STAT_CALL)
            write(*,*) 'On line String=',trim(OnlineString)
            if (STAT_CALL == -1) stop 'ConstructModel - ModuleModel - ERR18'
            do i=1, len(trim(OnlineString))
                if (OnlineString(i:i)=='%') OnlineString(i:i) = '&'
            enddo

        endif

        allocate(TimeReal(1:6))

        call GetDataOnlineString ('START', ArrayData = TimeReal)

        call SetDate(Me%BeginTime, TimeReal(1), TimeReal(2), TimeReal(3),               &
                                   TimeReal(4), TimeReal(5), TimeReal(6))

        call GetDataOnlineString ('END'  , ArrayData = TimeReal)

        call SetDate(Me%EndTime,   TimeReal(1), TimeReal(2), TimeReal(3),               &
                                   TimeReal(4), TimeReal(5), TimeReal(6))

        deallocate(TimeReal)

    end subroutine ReadStringOnline
    
#endif _ONLINE_


    !--------------------------------------------------------------------------

    subroutine AllocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Model), pointer                     :: NewModel
        type (T_Model), pointer                     :: PreviousModel


        !Allocates new instance
        allocate (NewModel)
        nullify  (NewModel%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstModel)) then
            FirstModel    => NewModel
            Me            => NewModel
        else
            PreviousModel => FirstModel
            Me            => FirstModel%Next
            do while (associated(Me))
                PreviousModel  => Me
                Me             => Me%Next
            enddo
            Me                 => NewModel
            PreviousModel%Next => NewModel
        endif

        Me%InstanceID = RegisterNewInstance (mMODEL_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    
    subroutine GetModelTimeLimits (ModelID, BeginTime, EndTime, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        type (T_Time), intent(OUT)                  :: BeginTime, EndTime
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            BeginTime = Me%BeginTime
            EndTime   = Me%EndTime
                
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetModelTimeLimits

    !--------------------------------------------------------------------------

    subroutine GetModelTimeStep (ModelID, DT, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        real, intent(OUT)                           :: DT
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call GetComputeTimeStep  (Me%ObjTime, DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetModelTimeStep - ModuleModel - ERR01'
                
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetModelTimeStep

    !--------------------------------------------------------------------------

    subroutine GetSubModelWindow (ModelID, SubModelWindow, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        type (T_Size2D)                             :: SubModelWindow
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            SubModelWindow = Me%SubModelWindow
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetSubModelWindow

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetSubModelWindowON (ModelID, SubModelWindowON, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        logical                                     :: SubModelWindowON
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            SubModelWindowON = Me%SubModelWindowON
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetSubModelWindowON

    !--------------------------------------------------------------------------

#ifdef OVERLAP

    subroutine GetModelOverlap (ModelID, Overlap, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        logical                                     :: Overlap
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            Overlap = Me%Overlap%Yes
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetModelOverlap
    
    !--------------------------------------------------------------------------

    subroutine GetModelOverlapInfo (ModelID, OverlapModelID, OverlapCells, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        integer                                     :: OverlapModelID
        integer, dimension(:,:), pointer            :: OverlapCells
        integer, intent(OUT),   optional            :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            OverlapModelID  =  Me%Overlap%OverlapModelID
            OverlapCells    => Me%Overlap%Cells
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetModelOverlapInfo

#endif


    !--------------------------------------------------------------------------

    subroutine GetModelCurrentTime (ModelID, CurrentTime, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        type (T_Time)                               :: CurrentTime
        integer,       intent(OUT),   optional      :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
    
            CurrentTime = Me%CurrentTime
                           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetModelCurrentTime

    !--------------------------------------------------------------------------

    subroutine GetModelInstanceIDs (ModelID, HorizontalGridID, HydrodynamicID,           &
                                    WaterPropertiesID, STAT)

        !Arguments------------------------------------------------------------ 
        integer                                     :: ModelID
        integer, intent(OUT), optional              :: HorizontalGridID
        integer, intent(OUT), optional              :: HydrodynamicID
        integer, intent(OUT), optional              :: WaterPropertiesID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_) 
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(HorizontalGridID )) HorizontalGridID  = Me%ObjHorizontalGrid
            if (present(HydrodynamicID   )) HydrodynamicID    = Me%ObjHydrodynamic
            if (present(WaterPropertiesID)) WaterPropertiesID = Me%ObjWaterProperties
                
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetModelInstanceIDs

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine UpdateTimeAndMapping (ModelID, Global_CurrentTime, DoNextStep, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ModelID
        type (T_Time),     intent(IN )              :: Global_CurrentTime
        logical                                     :: DoNextStep
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_)

if1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Gives the actualized Current time to the module time
            call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR01'

            DoNextStep = .false.

if9:        if (Me%CurrentTime .LT. Me%EndTime) then      

if2:            if (Global_CurrentTime .GE. Me%CurrentTime) then

                    call GetComputeTimeStep (Me%ObjTime, Me%DT, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR03'

                    !Actualize the CurrentTime with Model time interval AppTime%DT
                    call ActualizeCurrentTime(Me%ObjTime, Me%DT, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR04'
                   
                    !Gives the actualized Current time to the module time
                    call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR05'

                    call GetWaterLevel(Me%ObjHydrodynamic, Me%ExternalVar%WaterLevel,    &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR06'

                    !Update the moving boundary  (boundary of the tidal areas covered)
                    call UpdateComputeFaces3D(Me%Water%ObjMap, Me%ExternalVar%WaterLevel, &
                                              Me%CurrentTime, STAT = STAT_CALL)      
                    if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR07'

                    call UnGetHydrodynamic(Me%ObjHydrodynamic,                     &
                                           Me%ExternalVar%WaterLevel, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'UpdateTimeAndMapping - ModuleModel - ERR08'

                    if (Me%CurrentTime .LE. Me%EndTime) then
                        DoNextStep = .true.
                    endif

                end if if2

            end if if9

            STAT_ = SUCCESS_

        else               

            STAT_ = ready_

        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine UpdateTimeAndMapping

    !--------------------------------------------------------------------------

    subroutine RunModel(ModelID, DT_Father, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ModelID
        real                                        :: DT_Father
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        real                                        :: CPUTime, LastCPUTime = 0.
        real                                        :: NewDT
        type(T_NewDT)                               :: PredictedDT
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
#ifdef _USE_SEQASSIMILATION
        integer                                     :: Cyclenumber
        logical                                     :: VirtualRun = .false.
#endif _USE_SEQASSIMILATION
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_)

if1 :   if (ready_ .EQ. IDLE_ERR_) then

        if (MonitorPerformance) then
            call StartWatch ("ModuleModel", "RunModel")
        endif

#ifdef _USE_SEQASSIMILATION
           if (Me%RunSeqAssimilation) then

                if (Me%StateCovEvolution) then
                    
                    !Evolution of covariance structure and model state
                    
                    !   Run model for disturbed states
                    VirtualRun = .true.
                    call SetHydroVirtualRun(Me%ObjHydrodynamic, VirtualRun,     &
                                            STAT = STAT_CALL)
                    call SetWaterPropVirtualRun(Me%ObjWaterProperties,          &
                                                VirtualRun, STAT = STAT_CALL)

                    do cyclenumber = 1, Me%StateCovRank

                        call SetModelInitialState(cyclenumber)

                        call RunOneModel(PredictedDT, DT_Father)

                        !call GetModelResult(cyclenumber)
                    enddo

                    !   Run model for undisturbed state
                    VirtualRun = .false.
                    call SetHydroVirtualRun(Me%ObjHydrodynamic, VirtualRun,     &
                                            STAT = STAT_CALL)
                    call SetWaterPropVirtualRun(Me%ObjWaterProperties,          &
                                                VirtualRun, STAT = STAT_CALL)

                    call SetModelInitialState(Me%StateCovRank + 1)

                    call RunOneModel(PredictedDT, DT_Father)

                    !!   Calculate the new EOF set
                    !call CovarianceCalculationSEEK
                else

                    call RunOneModel(PredictedDT, DT_Father)
                endif

                if (Me%CurrentTime == Me%SeqAssimilationTime) then

                    !Calculate the new EOF set
                    call CovarianceCalculationSEEK

                    !Calculate analysis and correct model prediction
                    call ModifySequentialAssimilation(Me%ObjSeqAssimilation)

                    !Actualizes Me%SeqAssimilationTime
                    call GetSeqAssimilationTime(Me%ObjSeqAssimilation,          &
                                                Me%SeqAssimilationTime,         &
                                                STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunModel - ModuleModel - ERR01'
                endif
           else
#endif _USE_SEQASSIMILATION
                call RunOneModel(PredictedDT, DT_Father) 

#ifdef _USE_SEQASSIMILATION
           endif
#endif _USE_SEQASSIMILATION

           if (Me%VariableDT) then

                NewDT = PredictedDT%DT
                
                Me%Iteration = Me%Iteration + 1

!                if (NewDT > Me%DT) then
!                    !At maximum increase DT by 5%
!                    NewDT = min(Me%DT * 1.05, NewDT)
!                endif

                !Dont allow DT superior to Maximum DT
                if (NewDT > Me%MaxDT) then
                    NewDT = Me%MaxDT
                endif

                !Dont allow DT below Initial DT
                if (NewDT < Me%InitialDT) then
                    NewDT = Me%InitialDT
                endif

                !Fit DT so model will stop at the right time
                if (Me%CurrentTime + NewDT > Me%EndTime .and. Me%EndTime > Me%CurrentTime) then
                    NewDT = Me%EndTime - Me%CurrentTime
                endif

                !WriteDTLog: com i, j e k opcionais. O nome do módulo deve conter
                !a propriedade limitadora. O -99 deve conter o número da iteração.
                if (MonitorDT)                                                          &
                    call WriteDTLog (trim(Me%ModelName), Me%Iteration, NewDT,           &
                                     PredictedDT%i, PredictedDT%j, predictedDT%k,       &
                                     trim(predictedDT%property))

                !Actualize the Time Step
                call ActualizeDT(Me%ObjTime, NewDT, STAT = STAT_CALL)     
                if (STAT_CALL /= SUCCESS_) stop 'RunModel - ModuleModel - ERR26'
                
           endif 


            call CPU_TIME(CPUTime)
            if (CPUTime - LastCPUTime > Me%DTPredictionInterval) then
                LastCPUTime = CPUTime
                call PrintProgress(Me%ObjTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunModel - ModuleModel - ERR27'
            endif
                

            STAT_ = SUCCESS_
        else               

            STAT_ = ready_

        end if if1

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) then
            call StopWatch ("ModuleModel", "RunModel")
        endif

        !----------------------------------------------------------------------

    end subroutine RunModel

    !--------------------------------------------------------------------------

    subroutine RunOneModel(PredictedDT, DT_Father)

        !Arguments-------------------------------------------------------------
        real                                        :: DT_Father
        type(T_NewDT)                               :: PredictedDT

        !Local-----------------------------------------------------------------
        type(T_NewDT)                               :: PredHydroDT, PredWaterDT
        real                                        :: AuxReal
        integer                                     :: STAT_CALL, AuxInt
#ifndef _AIR_
        integer, dimension(:, :), pointer           :: WaterPoints2D
#endif
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleModel", "RunOneModel")

        !Updates the Interface between the sediment and the water column
        call ModifyInterfaceSedimentWater (Me%ObjInterfaceSedimentWater,                &
#ifdef  _LAGRANGIAN_GLOBAL_          
                                           LagrangianID = Me%ObjLagrangianGlobal,       &
#else
                                           LagrangianID = Me%ObjLagrangian,             &
#endif               
                                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR00'

#ifndef _AIR_
            
        !Gets Waterpoints
        call GetWaterPoints2D (Me%Water%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR01'

        !Starts Atmosphere
        call ModifyAtmosphere  (Me%ObjAtmosphere, MappingPoints = WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR02'

        !Gets Unget Waterpoints
        call UnGetHorizontalMap (Me%Water%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR03'

        call ModifyInterfaceWaterAir (Me%ObjInterfaceWaterAir,                          &
#ifdef  _LAGRANGIAN_GLOBAL_          
                                      LagrangianID = Me%ObjLagrangianGlobal,            &
#else
                                      LagrangianID = Me%ObjLagrangian,                  &
#endif               
                                      STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR04'
#endif
#ifndef _WAVES_
        if(Me%RunWaves)then
                            
            call ModifyWaves(Me%ObjWaves, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR05'
            
        end if
#endif
        call GetDensity(Me%ObjWaterProperties, Me%ExternalVar%Density,              &
                        Me%CurrentTime,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR06'

        call GetSigma(Me%ObjWaterProperties, Me%ExternalVar%SigmaDens,              &
                        Me%CurrentTime,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR06a'


        if(Me%ExternalVar%NeedsTempSalinity)then

            call GetConcentration(Me%ObjWaterProperties, Me%ExternalVar%Salinity,   &
                                  Salinity_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'RunOneModel - ModuleModel - ERR07'

            call GetConcentration(Me%ObjWaterProperties, Me%ExternalVar%Temperature,&
                                  Temperature_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'RunOneModel - ModuleModel - ERR08'

        endif

        call GetHorizontalVelocity(Me%ObjHydrodynamic,                              &
                                   Me%ExternalVar%VelocityX,                        &
                                   Me%ExternalVar%VelocityY,                        &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR10'

        call GetVerticalVelocity(HydrodynamicID = Me%ObjHydrodynamic,               &
                                   Velocity_W = Me%ExternalVar%VelocityZ,           &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunModel - ModuleModel - ERR11'

        call GetChezy(Me%ObjHydrodynamic, Me%ExternalVar%Chezy, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR12'
                    
        call Turbulence(Me%ObjTurbulence,                                           &
                        Me%ExternalVar%VelocityX,                                   &
                        Me%ExternalVar%VelocityY,                                   &
                        Me%ExternalVar%VelocityZ,                                   &
                        Me%ExternalVar%Chezy,                                       &
                        Me%ExternalVar%Salinity,                                    &
                        Me%ExternalVar%Temperature,                                 &
                        STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR13'

        call UnGetHydrodynamic(Me%ObjHydrodynamic,                                  &
                               Me%ExternalVar%VelocityX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR14'

        call UnGetHydrodynamic(Me%ObjHydrodynamic,                                  &
                               Me%ExternalVar%VelocityY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR15'
        
        call UnGetHydrodynamic(Me%ObjHydrodynamic,                                  &
                               Me%ExternalVar%VelocityZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR16'

        call UnGetHydrodynamic(Me%ObjHydrodynamic,                                  &
                               Me%ExternalVar%Chezy, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR17'

        PredictedDT%DT = Me%DT             
        call Modify_Hydrodynamic(Me%ObjHydrodynamic,                                &
                                 Me%ExternalVar%Density,                            &
                                 Me%ExternalVar%SigmaDens,                          &
                                 PredHydroDT,                                       &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR18'
                
        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExternalVar%Density,    &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR19'

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExternalVar%SigmaDens,  &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR19a'


        if(Me%ExternalVar%NeedsTempSalinity)then

            call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExternalVar%Salinity, &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR20'

            call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExternalVar%Temperature, &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR21'

        end if

        call WaterProperties_Evolution(Me%ObjWaterProperties, PredWaterDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR22'

        if (Me%VariableDT) then

            AuxReal = min(PredHydroDT%DT, PredWaterDT%DT, DT_Father)
            
            if (PredHydroDT%DT .eq. AuxReal) then
                PredictedDT = PredHydroDT
            elseif (PredWaterDT%DT .eq. AuxReal) then
                PredictedDT = PredWaterDT
            else
                PredictedDT%property = 'DT_Father'
                PredictedDT%DT = DT_Father
            end if
            
            if (DT_Father < - FillValueReal/100.) then
                AuxReal      = DT_Father/PredictedDT%DT
                AuxInt       = int(AuxReal) 
                if (AuxReal/= real(AuxInt)) AuxInt = AuxInt + 1
                PredictedDT%DT =  DT_Father / real(AuxInt)
            endif
            
            !write(998,*)  PredHydroDT, PredWaterDT, DT_Father, PredictedDT

        endif

#ifndef _SEDIMENT_
        if(Me%RunSediments)then
            
            call ModifyConsolidation(Me%ObjConsolidation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR23'

            call SedimentProperties_Evolution(Me%ObjSedimentProperties, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR24'

        end if
#endif

 
#ifndef _LAGRANGIAN_

#ifdef  _LAGRANGIAN_GLOBAL_          
            if (Me%InstanceID == Me%NumberOfModels .and. Me%ObjLagrangianGlobal /= 0) then
                call ModifyLagrangianGlobal   (Me%ObjLagrangianGlobal, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR30'
            endif
#else
         if (Me%ObjLagrangian /= 0) then           
            call ModifyLagrangian   (Me%ObjLagrangian, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOneModel - ModuleModel - ERR40'
        endif      
#endif
#endif 

        if (MonitorPerformance) call StopWatch ("ModuleModel", "RunOneModel")

        !----------------------------------------------------------------------

    end subroutine RunOneModel

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillModel (ModelID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ModelID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: nUsers
        integer                                     :: STAT_CALL
        real                                        :: DT_Error
               

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ModelID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mMODEL_,  Me%InstanceID)

            if (nUsers == 0) then

                call GetComputeCurrentTime(Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR01'

                DT_error = Me%EndTime - Me%CurrentTime
#ifndef _OUTPUT_OFF_
if9 :           if (DT_error /= 0.) then
                    write(*,*)  
                    write(*,*) 'Warning: The model = ',trim(Me%ModelName)
                    write(*,*) 'Didnt finish in the foreseen data'                   

if7 :               if     (DT_error > 0) then
                        write(*,*) 'The model stoped ',DT_error,' seconds before the foreseen data!'                  
                    elseif (DT_error < 0) then if7
                        write(*,*) 'The model stoped ',DT_error,' seconds after the foreseen data!'                             
                    end if if7
                    write(*,*) '                        '
                    write(*,*) 'SUBROUTINE KillModel; ModuleModel. WRN01.'
                    write(*,*)  
                end if if9
#endif
                !Last Progress message
                call PrintProgress(Me%ObjTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR10'

#ifdef _USE_SEQASSIMILATION
                if (Me%RunSeqAssimilation) then
                    call KillSequentialAssimilation (Me%ObjSeqAssimilation,         &
                                                     STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR15'
                endif
#endif _USE_SEQASSIMILATION
                call KillInterfaceSedimentWater (Me%ObjInterfaceSedimentWater, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR20'
#ifndef _AIR_
                call KillInterfaceWaterAir (Me%ObjInterfaceWaterAir, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR30'
#endif

#ifndef _LAGRANGIAN_
#ifdef  _LAGRANGIAN_GLOBAL_                  
                if (Me%ObjLagrangianGlobal /= 0) then

                    if (Me%InstanceID == 1) then
                        call DeallocateLagrangianGlobal (Me%ObjLagrangianGlobal, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR35'
                    endif

                    if (Me%InstanceID == Me%NumberOfModels) then
                        call KillLagrangianGlobal       (Me%ObjLagrangianGlobal, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR40'
                    endif
                endif
#else
                if (Me%ObjLagrangian /= 0) then
                    call KillLagrangian (Me%ObjLagrangian, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR45'
                endif
#endif
#endif

#ifndef _AIR_
                call KillAtmosphere (Me%ObjAtmosphere, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR50'
#endif
#ifndef _SEDIMENT_
                if (Me%RunSediments) then

                    call KillSedimentProperties (Me%ObjSedimentProperties,      STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR95'         

                    call KillConsolidation      (Me%ObjConsolidation,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR95'

                    !Kills Map
                    call KillMap                (Me%Sediment%ObjMap,           STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR60'

                    !Kills Geometry
                    call KillGeometry           (Me%Sediment%ObjGeometry,      STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR70'

                    !Kills HorizontalMap
                    call KillHorizontalMap      (Me%Sediment%ObjHorizontalMap, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR80'

                    !Kills Bathymetry
                    call KillGridData           (Me%Sediment%ObjBathymetry,    STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR90'


                endif
#endif

                !Kills water properties
                call KillWaterProperties(Me%ObjWaterProperties,     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR100'

                !Kills hydrodynamic properties
                call KillHydrodynamic(Me%ObjHydrodynamic,           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR110'

                !Kill Turbulence
                call KillTurbulence     (Me%ObjTurbulence,          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR120'
#ifndef _WAVES_
                if (Me%ObjWaves /= 0) then
                    call KillWaves (Me%ObjWaves, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR130'
                endif
#endif
                !Kills Map
                call KillMap            (Me%Water%ObjMap,           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR140'

                !Kills Geometry
                call KillGeometry       (Me%Water%ObjGeometry,      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR150'

                !Kills HorizontalMap
                call KillHorizontalMap  (Me%Water%ObjHorizontalMap, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR160'

                !Kills Bathymetry
                call KillGridData       (Me%Water%ObjBathymetry,    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR170'

                !Kills HorizontalGrid
                call KillHorizontalGrid (Me%ObjHorizontalGrid,      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR180'

                !Kills Compute Time
                call KillComputeTime    (Me%ObjTime,                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR190'
                
#ifdef _ENABLE_CUDA                
                !Kills ModuleCuda
                call KillCuda           (Me%ObjCuda,                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR200'
#endif _ENABLE_CUDA                
                
                !ObjHorizontalGrid
                call DeallocateInstance
                
                ModelID = 0
                STAT_   = SUCCESS_

            end if

        else 

            STAT_ = ready_

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    End Subroutine KillModel
    
    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Model), pointer          :: AuxObjModel
        type (T_Model), pointer          :: PreviousObjModel

        !Updates pointers
        if (Me%InstanceID == FirstModel%InstanceID) then
            FirstModel => FirstModel%Next
        else
            PreviousObjModel => FirstModel
            AuxObjModel      => FirstModel%Next
            do while (AuxObjModel%InstanceID /= Me%InstanceID)
                PreviousObjModel => AuxObjModel
                AuxObjModel      => AuxObjModel%Next
            enddo

            !Now update linked list
            PreviousObjModel%Next => AuxObjModel%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 
            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (ModelID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ModelID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ModelID > 0) then

            Me => FirstModel
            do while (associated (Me))
                if (Me%InstanceID == ModelID) exit
                Me => Me%Next
            enddo

            if (.not. associated(Me)) stop 'ModuleModel - Ready - ERR01'

            ready_ = VerifyReadLock (mMODEL_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

end module ModuleModel

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
