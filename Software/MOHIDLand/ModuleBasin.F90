!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! MODULE        : Basin
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank - v4.0
! DESCRIPTION   : Module Basin is the top level of RunOff and Infiltration 
!
!------------------------------------------------------------------------------
! ATMOSPHERE                    : 0/1              [1]          !Use Module Atmosphere 
! POROUS_MEDIA                  : 0/1              [1]          !Use Module Porous Media
! RUN_OFF                       : 0/1              [1]          !Use Module RunOff
! DRAINAGE_NET                  : 0/1              [1]          !Use Module DrainageNetork
! OUTPUT_TIME                   : sec. sec. sec.    -           !Output Time
! TIME_SERIE_LOCATION           : char              -           !Path to time serie location file
! INITIAL_LOSS                  : real             [0.0]        !Coeficient of initial rainfall losses
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

Module ModuleBasin

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleTimeSerie
    use ModuleHDF5
    use ModuleFunctions,      only : ReadTimeKeyWords, LatentHeat, ConstructPropertyID,  &
                                     TimeToString, ChangeSuffix, CHUNK_J
    use ModuleFillMatrix,     only : ConstructFillMatrix, ModifyFillMatrix,              &
                                     KillFillMatrix,GetIfMatrixRemainsConstant
    use ModuleHorizontalGrid, only : ConstructHorizontalGrid, KillHorizontalGrid,        &
                                     WriteHorizontalGrid, GetHorizontalGridSize,         &
                                     GetGridCellArea, UnGetHorizontalGrid,               &
                                     GetHorizontalGrid, GetXYCellZ
    use ModuleHorizontalMap,  only : ConstructHorizontalMap, KillHorizontalMap,          &
                                     UpdateComputeFaces2D, GetOpenPoints2D,              &
                                     GetBoundaries, UngetHorizontalMap 
    use ModuleGridData,       only : ConstructGridData,       KillGridData,              &
                                     GetGridData, UngetGridData
    use ModuleBasinGeometry,  only : ConstructBasinGeometry,  KillBasinGeometry,         &
                                     GetBasinPoints, GetRiverPoints, UngetBasin
    use ModuleAtmosphere,     only : StartAtmosphere, ModifyAtmosphere,                  &
                                     GetAtmosphereProperty, GetAtmosphereDTPrediction,   &
                                     GetAtmospherePropertiesIDByIdx,                     &
                                     GetAtmospherenProperties, AtmospherePropertyExists, &
                                     UnGetAtmosphere, KillAtmosphere
    use ModuleRunOff,         only : ConstructRunOff, ModifyRunOff, GetOverLandFlow,     &
                                     GetFlowToChannels, GetNextRunOffDT,                 &
                                     GetFlowAtBoundary, UnGetRunOff, KillRunOff
    use ModuleRunoffProperties,                                                          &
                              only : ConstructRunoffProperties,                          &
                                     ModifyRunoffProperties,                             &
                                     KillRunoffProperties, GetRPConcentration,           &
                                     GetRPnProperties, GetRPOptions,                     &
                                     GetRPPropertiesIDByIdx, SetDNConcRP,                &
                                     UngetRunoffProperties, SetBasinConcRP,              &
                                     SetBasinToRPSplash                                                                         
    use ModuleDrainageNetwork,only : ConstructDrainageNetwork, FillOutPutMatrix,         &
                                     ModifyDrainageNetwork,                              &
                                     GetHasProperties, GetDNnProperties,                 &
                                     GetHasToxicity,GetNeedsRadiation,                   &
                                     GetNeedsAtmosphere, SetAtmosphereDrainageNet,       &
                                     GetDNPropertiesIDByIdx, GetNextDrainageNetDT,       &
                                     GetVolumes, GetPropHasBottomFluxes,                 &
                                     GetChannelsNodeLength,GetChannelsID,GetDNConcentration,&
                                     SetPMPConcDN,SetRPConcDN, UnGetDrainageNetwork,     &
                                     KillDrainageNetwork, SetGWFlowLayersToDN
    use ModulePorousMedia,    only : ConstructPorousMedia, ModifyPorousMedia,            &
                                     KillPorousMedia, GetGWFlowToChannels,               &
                                     GetInfiltration, GetEfectiveEVTP,                   &
                                     GetTotalStoredVolume, GetEvaporation,               &
                                     GetNextPorousMediaDT, UngetPorousMedia, GetGWLayer, &
                                     GetGWFlowOption, GetGWFlowToChannelsByLayer,        &
                                     GetGWToChannelsLayers
    use ModulePorousMediaProperties,                                                     &
                              only : ConstructPorousMediaProperties,                     &
                                     ModifyPorousMediaProperties,                        &
                                     KillPorousMediaProperties, GetPMPConcentration,     &
                                     SetVegetationPMProperties, GetPMPCoupled,           &
                                     SetWindVelocity, GetPMPnProperties,                 &
                                     GetPMPPropertiesIDByIdx, SetDNConcPMP,              &
                                     UngetPorousMediaProperties,                         &
                                     SetInfColConcPMP
    use ModuleVegetation,     only : ConstructVegetation, ModifyVegetation,              &
                                     KillVegetation, GetLeafAreaIndex,                   &
                                     GetSpecificLeafStorage, GetEVTPCropCoefficient,     &
                                     GetTranspiration, GetVegetationSoilFluxes,          &
                                     SetSoilConcVegetation, GetVegetationOptions,        &
                                     GetVegetationDT, GetRootDepth, GetNutrientFraction, &
                                     UnGetVegetation, UnGetVegetationSoilFluxes,         &
                                     GetCanopyHeight
    use ModuleStopWatch      ,only : StartWatch, StopWatch
    use ModuleGeometry,       only : GetGeometrySize
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructBasin
    private ::      AllocateInstance
    private ::      ReadFileNames
    private ::      ReadDataFile
    private ::      VerifyOptions
    private ::      AllocateVariables
    private ::      ConstructCoupledModules
    private ::          ConstructHDF5Output
    private ::      ConstructTimeSeries
    private ::      ReadInitialFile

    !Selector
                    
    
    !Modifier
    public  :: ModifyBasin
    private ::      AtmosphereProcesses
    private ::          DividePrecipitation
    private ::          CalcPotEvapoTranspiration
    private ::      VegetationProcesses
    private ::      OverLandProcesses
    private ::      DrainageNetworkProcesses
    private ::      PorousMediaProcesses
    private ::      PorousMediaPropertiesProcesses
    private ::      HDF5Output
    private ::      TimeSerieOutput
    private ::      GlobalMassBalance

    !Destructor
    public  :: KillBasin                                                     
    private ::      DeAllocateInstance
    private ::      WriteFinalFile

    !Management
    private ::      Ready
    private ::          LocateObjBasin 
    private :: ReadLockExternalVar
    private :: ReadUnLockExternalVar

    
    !Interfaces----------------------------------------------------------------

    !Parameters----------------------------------------------------------------
    character(LEN = StringLength), parameter        :: block_begin          = '<beginproperty>'
    character(LEN = StringLength), parameter        :: block_end            = '<endproperty>'
    !Separate evapotranspiration
    integer, parameter                              :: SingleEvapoTranspiration   = 1
    integer, parameter                              :: SeparateEvapoTranspiration = 2
    !Gw link between porous media and drainage network
    integer, parameter                              :: Layer_ = 3

    !Types---------------------------------------------------------------------
    type T_OutPut
        type (T_Time), dimension(:), pointer        :: OutTime
        type (T_Time), dimension(:), pointer        :: RestartOutTime
        integer                                     :: NextOutPut           = 1
        real   , dimension(:,:), pointer            :: OutputChannels
        logical                                     :: Yes                  = .false.
        logical                                     :: WriteRestartFile     = .false.
        logical                                     :: RestartOverwrite     = .false.
        integer                                     :: NextRestartOutput    = 1
    end type T_OutPut

    type T_Coupling
        logical                                     :: Atmosphere
        logical                                     :: Evapotranspiration
        logical                                     :: RunOff
        logical                                     :: RunOffProperties
        logical                                     :: DrainageNetwork
        logical                                     :: PorousMedia
        logical                                     :: PorousMediaProperties
        logical                                     :: Vegetation
        logical                                     :: SimpleInfiltration
        logical                                     :: Snow
    end type T_Coupling

    type T_ExtVar
        integer, dimension(:,:), pointer            :: BasinPoints
        integer, dimension(:,:), pointer            :: RiverPoints
        integer, dimension(:,:), pointer            :: OpenPoints2D
        integer, dimension(:,:), pointer            :: BoundaryPoints2D
        real   , dimension(:,:), pointer            :: GridCellArea
        real   , dimension(:,:), pointer            :: Topography
        real   , dimension(:,:), pointer            :: LeafAreaIndex 
        real   , dimension(:,:), pointer            :: SpecificLeafStorage
        real   , dimension(:,:), pointer            :: CropCoefficient
        real   , dimension(:,:,:), pointer          :: ActualTranspiration
!        real   , dimension(:,:), pointer            :: PermeableFraction
    end type T_ExtVar

    type T_Files
        character(len=PathLength)                   :: ConstructData
        character(len=PathLength)                   :: TopographicFile
        character(len=PathLength)                   :: HDFFile
        character(len=PathLength)                   :: InitialFile
        character(len=PathLength)                   :: FinalFile
        character(len=PathLength)                   :: TimeSerieLocation
    end type T_Files

    type T_IntegratedFlow
        logical                                     :: On
        real                                        :: Flow
        integer                                     :: CurrentIndex
        integer                                     :: ObjTimeSerie         = 0
    end type T_IntegratedFlow

    type       T_BasinProperty
        type (T_PropertyID)                         :: ID
        real, dimension(:,:), pointer               :: Field                => null()
        logical                                     :: Constant             = .false.
        type (T_BasinProperty), pointer             :: Next                 => null()
    end type T_BasinProperty    

    type T_PropertyB
        type (T_PropertyID)                         :: ID
        real, dimension(:,:), pointer               :: Field                => null()
    end type T_PropertyB

    type T_SimpleInfiltration
        type (T_PropertyB)                          :: Ks                   
        type (T_PropertyB)                          :: MP                   ! MatricPotential
        type (T_PropertyB)                          :: ThetaS               
        type (T_PropertyB)                          :: ThetaI               
        type (T_PropertyB)                          :: InfRate              
        type (T_PropertyB)                          :: AccInf               
    end type T_SimpleInfiltration

    type T_MassBalance
        !Initial Mass
        real(8)                                     :: IniVolumeBasin       = 0.
        real(8)                                     :: IniVolumeVegetation  = 0.
        real(8)                                     :: IniVolumePorousMedia = 0.
        real(8)                                     :: IniVolumeChannels    = 0.
        !Sink / Sources
        real(8)                                     :: EvapFromVegetation   = 0.
        real(8)                                     :: EvapFromGround       = 0.
        real(8)                                     :: EvapFromSoil         = 0.
        real(8)                                     :: RainIn               = 0.
        real(8)                                     :: DischargesIn         = 0.
        real(8)                                     :: OverTopOut           = 0.
        real(8)                                     :: OutVolumeChannel     = 0.
        real(8)                                     :: OutVolumeOverLand    = 0.
        !Final Mass
        real(8)                                     :: FinalVolumeBasin
        real(8)                                     :: FinalVolumeVegetation
        real(8)                                     :: FinalVolumePorousMedia
        real(8)                                     :: FinalVolumeChannels
        !Flow exchange between modules
        real(8)                                     :: OLFlowToRiver
        real(8)                                     :: GWFlowToRiver
        real(8)                                     :: Infiltration
    end type T_MassBalance

    type       T_Basin
        integer                                     :: InstanceID           = 0
        character(len=StringLength)                 :: ModelName
        type (T_Size2D)                             :: Size, WorkSize
        type (T_Coupling)                           :: Coupled
        type (T_Files)                              :: Files
        type (T_ExtVar)                             :: ExtVar
        type (T_SimpleInfiltration)                 :: SI                                       
        type (T_BasinProperty), pointer             :: FirstProperty        => null()
        logical                                     :: Continuous           = .false.
        logical                                     :: StopOnWrongDate      = .true.
        logical                                     :: VerifyGlobalMass     = .false.
        logical                                     :: Calibrating1D        = .false.
        logical                                     :: ConcentrateRain      = .false.
        real                                        :: RainAverageDuration  = 600.0
        real                                        :: WCRemovalTime        = 600.
        real                                        :: DTDuringRain
        logical                                     :: DiffuseWaterSource   = .false.
        real                                        :: FlowPerCapita        = 0.0
        character(PathLength)                       :: PopDensityFile
        real                                        :: ExtinctionCoef       = 0.6
        real                                        :: CurrentDT
        integer                                     :: EvapoTranspirationMethod
        logical                                     :: ConstructEvaporation
        logical                                     :: ConstructTranspiration
        real(8), dimension(:,:), pointer            :: WaterLevel             => null()
        real(8), dimension(:,:), pointer            :: WaterColumn            => null()
        real(8), dimension(:,:), pointer            :: WaterColumnOld         => null()
        real(8), dimension(:,:), pointer            :: WaterColumnRemoved     => null()
        real(8), dimension(:,:), pointer            :: CanopyStorageCapacity  => null()
        real(8), dimension(:,:), pointer            :: CanopyStorage          => null()
        real(8), dimension(:,:), pointer            :: CanopyDrainage         => null()
        real(8), dimension(:,:), pointer            :: SnowPack               => null()
        real(8), dimension(:,:), pointer            :: ThroughFall            => null()
        real   , dimension(:,:), pointer            :: CoveredFraction        => null()
        real   , dimension(:,:), pointer            :: CropEvapotrans         => null()
        real   , dimension(:,:), pointer            :: PotentialTranspiration => null() 
        real   , dimension(:,:), pointer            :: PotentialEvaporation   => null() 
        real(8), dimension(:,:), pointer            :: PotentialInfCol        => null()
        real(8), dimension(:,:), pointer            :: FlowProduction         => null()
        real(8), dimension(:,:), pointer            :: InfiltrationRate       => null()
        real(8), dimension(:,:), pointer            :: PrecipRate             => null()
        real(8), dimension(:,:), pointer            :: ThroughRate            => null()
        real(8), dimension(:,:), pointer            :: EVTPRate               => null()
        real(8), dimension(:,:), pointer            :: EVTPRate2              => null()
        real(8), dimension(:,:), pointer            :: PlantWaterStress       => null()
        real(8), dimension(:,:), pointer            :: AccInfiltration        => null()
        real(8), dimension(:,:), pointer            :: AccFlowProduction      => null()
        real(8), dimension(:,:), pointer            :: AccEVTP                => null()
        real(8), dimension(:,:), pointer            :: AccRainFall            => null()
        real(8), dimension(:,:), pointer            :: AccEVPCanopy           => null()
        real,    dimension(:,:), pointer            :: AccRainHour            => null()
        real,    dimension(:,:), pointer            :: RainStartTime          => null()
        real,    dimension(:,:), pointer            :: RainDuration           => null()
        real,    dimension(:,:), pointer            :: DiffuseFlow            => null()
        real                                        :: InitialWaterColumn
        real                                        :: WaterColumnCoef
        type (T_MassBalance)                        :: MB
        type (T_IntegratedFlow)                     :: DailyFlow            
        type (T_IntegratedFlow)                     :: MonthlyFlow            
        real, dimension(:), pointer                 :: TimeSeriesBuffer
        real, dimension(:), pointer                 :: TimeSeriesBuffer2


        type (T_Time)                               :: CurrentTime
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime       
        type (T_OutPut)                             :: OutPut
        
        !Instance IDs
        integer                                     :: ObjTime                  = 0
        integer                                     :: ObjEnterData             = 0
        integer                                     :: ObjHorizontalGrid        = 0
        integer                                     :: ObjHorizontalMap         = 0
        integer                                     :: ObjGridData              = 0
        integer                                     :: ObjBasinGeometry         = 0
        integer                                     :: ObjGeometry              = 0  !only for getting KUB
        integer                                     :: ObjAtmosphere            = 0
        integer                                     :: ObjRunOff                = 0
        integer                                     :: ObjRunOffProperties      = 0
        integer                                     :: ObjDrainageNetwork       = 0
        integer                                     :: ObjPorousMedia           = 0
        integer                                     :: ObjVegetation            = 0
        integer                                     :: ObjHDF5                  = 0
        integer                                     :: ObjTimeSerie             = 0
        integer                                     :: ObjTimeSerieBasin        = 0
        integer                                     :: ObjPorousMediaProperties = 0
        type (T_Basin), pointer                     :: Next

        !Used by PorousMediaProperties        
        real(8), dimension(:,:), pointer            :: WaterColumnEvaporated    => null() !in meters (m)
    end type  T_Basin

    !Global Module Variables
    type (T_Basin), pointer                         :: FirstObjBasin        => null()
    type (T_Basin), pointer                         :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructBasin(ObjBasinID, ObjTime, ModelName, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjBasinID 
        integer                                         :: ObjTime 
        character(len=*)                                :: ModelName
        integer, optional, intent(OUT)                  :: STAT     

        !Local-------------------------------------------------------------------
        integer                                         :: ready_      
        integer                                         :: STAT_, STAT_CALL
        logical                                         :: VariableDT
        character (Len = StringLength)                  :: WarningString
        character (Len = StringLength)                  :: LockToWhichModules
        character (Len = StringLength)                  :: UnLockToWhichModules

        !------------------------------------------------------------------------
!        integer                                     :: TID, OMP_GET_THREAD_NUM


        !!call OMP_SET_DYNAMIC(.false.)
        !!call OMP_SET_NUM_THREADS(2)

        !!!!$OMP PARALLEL SHARED(Me) PRIVATE(TID)
        !!!!TID = OMP_GET_THREAD_NUM()
        !!!!PRINT *, 'Hello World from thread = ', TID
        !!!!$OMP END PARALLEL         


        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBasin_)) then
            nullify (FirstObjBasin)
            call RegisterModule (mBasin_) 
        endif

        call Ready(ObjBasinID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Returns ID
            ObjBasinID          = Me%InstanceID
            
            Me%ModelName        = ModelName

            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,   ObjTime)

            call GetComputeCurrentTime  (Me%ObjTime, CurrentTime = Me%CurrentTime,       &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR01'

            call GetVariableDT          (Me%ObjTime, VariableDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR02'
  
            call GetComputeTimeLimits   (Me%ObjTime, BeginTime = Me%BeginTime,           &
                                         EndTime = Me%EndTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR03'

            call ReadFileNames        
        
            !Constructs Horizontal Grid
            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%Files%TopographicFile, &
                                         STAT = STAT_CALL)           
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR04'

            !Constructs GridData
            call ConstructGridData      (Me%ObjGridData, Me%ObjHorizontalGrid,           &
                                         FileName = Me%Files%TopographicFile,            &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR05'

            !Constructs BasinGeometry
            call ConstructBasinGeometry (Me%ObjBasinGeometry, Me%ObjGridData,            &
                                         Me%ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR06'

            !Gets BasinPoints
            call GetBasinPoints         (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR07'

            !Constructs Mapping
            call ConstructHorizontalMap (Me%ObjHorizontalMap, Me%ObjGridData,            &
                                        Me%ObjHorizontalGrid, Me%CurrentTime,            &
                                        Me%ExtVar%BasinPoints, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR08'

            !Ungets BasinPoints
            call UngetBasin             (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR09'

            !Gets the size of the grid
            call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                        Size     = Me%Size,                              &
                                        WorkSize = Me%WorkSize,                          &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR09'

            !Gets ExternalVars
            LockToWhichModules = 'AllModules'
            call ReadLockExternalVar (LockToWhichModules)

            !Reads Data File
            call ReadDataFile

            !Allocates Variables
            call AllocateVariables      ()

            !Verifies User Options
            call VerifyOptions

            !Constructs Coupled Modules
            call ConstructCoupledModules()

            !Constructs the property list 
            call ConstructPropertyList

            !Constructs Output Time Series
            call ConstructTimeSeries
            
            !Reads conditions from previous run
            if (Me%Continuous) call ReadInitialFile
          
            !Actualizes the WaterColumn
            WarningString = 'ConstructBasin'
            call ActualizeWaterColumn (WarningString)

            !Closes Data File
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR99'


            !UnGets ExternalVars
            UnLockToWhichModules = 'AllModules'
            call ReadUnLockExternalVar (UnLockToWhichModules)

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleBasin - ConstructBasin - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructBasin

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Basin), pointer                         :: NewObjBasin
        type (T_Basin), pointer                         :: PreviousObjBasin


        !Allocates new instance
        allocate (NewObjBasin)
        nullify  (NewObjBasin%Next)

        nullify  (NewObjBasin%WaterLevel)
        nullify  (NewObjBasin%WaterColumn)
        nullify  (NewObjBasin%Output%OutTime)
        nullify  (NewObjBasin%Output%OutputChannels)

        nullify  (NewObjBasin%ExtVar%BasinPoints)
        nullify  (NewObjBasin%ExtVar%RiverPoints)
        nullify  (NewObjBasin%ExtVar%OpenPoints2D)
        nullify  (NewObjBasin%ExtVar%BoundaryPoints2D)
        nullify  (NewObjBasin%ExtVar%GridCellArea)
        nullify  (NewObjBasin%ExtVar%Topography )

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjBasin)) then
            FirstObjBasin         => NewObjBasin
            Me                    => NewObjBasin
        else
            PreviousObjBasin      => FirstObjBasin
            Me                    => FirstObjBasin%Next
            do while (associated(Me))
                PreviousObjBasin  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjBasin
            PreviousObjBasin%Next => NewObjBasin
        endif

        Me%InstanceID = RegisterNewInstance (mBASIN_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadFileNames

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_CALL        

        !------------------------------------------------------------------------

        !Opens Basin data file
        call ReadFileName('BASIN_DATA', Me%Files%ConstructData,                         &
                           Message = "Basin Data", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR01'      
       
       !Reads file name of the topographic file
        call ReadFileName('IN_BASIN', Me%Files%TopographicFile,                         &
                           Message = "Topographic Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR03'

        !Reads file name of the hdf outupt
        call ReadFileName('BASIN_HDF', Me%Files%HDFFile,                                &
                           Message = "Basin HDF Output File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR04'

        !Reads file name of final output file
        call ReadFileName('BASIN_FIN', Me%Files%FinalFile,                              &
                           Message = "Basin Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR05'

    end subroutine ReadFileNames

    !--------------------------------------------------------------------------

    subroutine ReadDataFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag


        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR00'

        !Basin Initial Water Column
        call GetData(Me%InitialWaterColumn,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'INITIAL_WATER_COLUMN',                              &
                     default      = 0.0,                                                 & 
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR01'

        !Basin Initial Water Column
        call GetData(Me%WaterColumnCoef,                                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'WATER_COLUMN_COEF',                                 &
                     default      = 0.0,                                                 & 
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR01'


        !Continuous Computation
        call GetData(Me%Continuous,                                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONTINUOUS',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR02'

        call GetData(Me%StopOnWrongDate,                                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'STOP_ON_WRONG_DATE',                                &
                     default      = .true.,                                              &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR03'

        !Verify Global Mass
        call GetData(Me%VerifyGlobalMass,                                                &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'VERIFY_MASS',                                       &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR04'

        !Calibrating 1D column?
        call GetData(Me%Calibrating1D,                                                   &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CALIBRATING_1D',                                    &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR04a'

        !How much to remove for the watercolumn 
        if (Me%Calibrating1D) then
            call GetData(Me%WCRemovalTime,                                                &
                         Me%ObjEnterData, iflag,                                          &
                         SearchType   = FromFile,                                         &
                         keyword      = 'WC_REMOVE_TIME',                                 &
                         default      = 600.,                                             &
                         ClientModule = 'ModuleBasin',                                    &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR04b'
        endif

        !DT During Rain
        call GetData(Me%DTDuringRain,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DT_DURING_RAIN',                                    &
                     default      = 60.0,                                                &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR04a'

        !Concentrate Rain for subhourly steps
        call GetData(Me%ConcentrateRain,                                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONCENTRATE_RAIN',                                  &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR04b'
        
        !How long does rain last...
        if (Me%ConcentrateRain) then
            call GetData(Me%RainAverageDuration,                                         &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'RAIN_AVERAGE_DURATION',                         &
                         default      = 600.00,                                          &
                         ClientModule = 'ModuleBasin',                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR04c'
            
            if (1.5 * Me%RainAverageDuration + Me%DTDuringRain > 3600.0) then
                write(*,*)'1.5 * RainAverageDuration + DT During Rain cannot exceed 3600s'
                stop 'ReadDataFile - ModuleBasin - ERR04d'
            endif
            
        endif

        !Imposes difuse sewage from a grid data
        call GetData(Me%DiffuseWaterSource,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DIFFUSE_WATER_SOURCE',                              &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR95'
        
        if (Me%DiffuseWaterSource) then
            call GetData(Me%FlowPerCapita,                                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'FLOW_PER_CAPITA',                               &
                         default      = 0.0,                                             &
                         ClientModule = 'ModuleBasin',                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR96'

            call GetData(Me%PopDensityFile,                                              &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'POPULATION_DENSITY',                            &
                         ClientModule = 'ModuleBasin',                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR97'
            if (iflag == 0) then
                write(*,*)'Population Density file not given - keyword POPULATION_DENSITY'
                stop 'ReadDataFile - ModuleBasin - ERR98'
            endif

        endif



        !Extinction Coeficient (relates LAI to covered fraction)
        call GetData(Me%ExtinctionCoef,                                                  &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'EXTINCTION_COEF',                                   &
                     default      = 0.6,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR40'

        !Reads file name of initial condition file
        if (Me%Continuous) then
            call ReadFileName('BASIN_INI', Me%Files%InitialFile,                         &
                               Message = "Basin Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR05'
        endif

        !Verifies if the user wants to use the Atmosphere Condition
        call GetData(Me%Coupled%Atmosphere,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'ATMOSPHERE',                                        &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR06'

        !Verifies if the user wants to use the Atmosphere Condition
        call GetData(Me%Coupled%Evapotranspiration,                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'EVAPOTRANSPIRATION',                                &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR07'

        if (Me%Coupled%Evapotranspiration) then
            !Verifies which method user wants for evapotranspiration (1-EvapoTranspiration  
            !2-Evaporation and Transpiration separated)
            call GetData(Me%EvapoTranspirationMethod,                                        &
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'EVAPOTRANSPIRATION_METHOD',                         &
                         default      = SingleEvapotranspiration,                            &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR10b'
        endif



        !Verifies if the user wants to simulate Infiltration
        call GetData(Me%Coupled%PorousMedia,                                             &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'POROUS_MEDIA',                                      &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR08'

        if (Me%Coupled%PorousMedia) then
            ! verifies if the user wants to simulate transport of properties
            call GetData(Me%Coupled%PorousMediaProperties,                                   &
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'POROUS_MEDIA_PROPERTIES',                           &
                         default      = OFF,                                                 &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR08'
        endif

        !Verifies if the user wants to simulate OverLand RunOff
        call GetData(Me%Coupled%RunOff,                                                  & 
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'RUN_OFF',                                           &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR09'
        
        if (Me%Coupled%RunOff) then
            !Verifies if the user wants to simulate OverLand RunOff propertie transport
            call GetData(Me%Coupled%RunOffProperties,                                        & 
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'RUN_OFF_PROPERTIES',                                &
                         default      = OFF,                                                 &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR09.5'
        endif

        !A Drainage Network is coupled?
        call GetData(Me%Coupled%DrainageNetwork,                                         &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DRAINAGE_NET',                                      &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR10'

        !The Vegetation Module is coupled_
        call GetData(Me%Coupled%Vegetation,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'VEGETATION',                                        &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR10a'

        !Verifies if the user wants to use simple 
        call GetData(Me%Coupled%SimpleInfiltration,                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SIMPLE_INFILTRATION',                               &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR10b'


        !Gets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime = Me%CurrentTime,                                 &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR11'

        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%CurrentTime,                               &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR11a'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleBasin - ERR11a'


        !Gets TimeSerieLocationFile
        call GetData(Me%Files%TimeSerieLocation,                                         &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TIME_SERIE_LOCATION',                               &
                     ClientModule = 'ModuleBasin',                                       &
                     Default      = Me%Files%ConstructData,                              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR12'

        !Output daily flow values?
        call GetData(Me%DailyFlow%On,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DAILY_FLOW',                                        &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR13'

        !Output monthly flow values?
        call GetData(Me%MonthlyFlow%On,                                                  &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MONTHLY_FLOW',                                      &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR14'

        !Verifies if the user wants to have precipitations has snow
        call GetData(Me%Coupled%Snow,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SNOW',                                              &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR15'



    end subroutine ReadDataFile
    !--------------------------------------------------------------------------

    subroutine VerifyOptions

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        if (Me%Coupled%DrainageNetwork .and. .not. Me%Coupled%Runoff) then
            write(*,*)'You must enable module Runoff if you want to use module Drainage Network'
            stop 'VerifyOptions - ModuleBasin - ERR01'
        endif
        
!        if (Me%Coupled%PorousMedia .and. .not. Me%Coupled%Vegetation) then
!            write(*,*)'You must enable module Vegetation if you want to use module Porous Media'
!            stop 'VerifyOptions - ModuleBasin - ERR02'
!        endif

        if (Me%Coupled%PorousMedia .and. Me%Coupled%SimpleInfiltration) then
            write(*,*)'You can use SimpleInfiltration and PorousMedia at the same time'
            stop 'VerifyOptions - ModuleBasin - ERR03'
        endif

        if (Me%Coupled%Vegetation) then
            
            if (.not. Me%Coupled%PorousMedia .or. .not. Me%Coupled%Evapotranspiration) then
                write(*,*)'You can not use Vegetation without PorousMedia or Evapotranspiration'
                write(*,*)'Check basin data file'
                stop 'VerifyOptions - ModuleBasin - ERR04'
            endif
            
            Me%ConstructTranspiration   = .true.
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                Me%ConstructEvaporation = .true.
            else
                Me%ConstructEvaporation = .false.
            endif

        else
            
            if(Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                write(*,*)'If Vegetation is not used then the consistent Evapotranspiration_Method '
                write(*,*)'is the default - without implicit transpiration computation.'
                write(*,*)'Check basin data file'
                stop 'VerifyOptions - ModuleBasin - ERR05'  
            endif             
            
            Me%ConstructTranspiration   = .false.
            if (Me%Coupled%Evapotranspiration) then
                Me%ConstructEvaporation     = .true.
            else
                Me%ConstructEvaporation     = .false.
            endif
        endif

        if (Me%Coupled%PorousMediaProperties) then
            if(.not. Me%Coupled%RunoffProperties) then
                write(*,*)'If using porous media properties also need runoff properties active '
                write(*,*)'Check basin data file in RUN_OFF_PROPERTIES'
                stop 'VerifyOptions - ModuleBasin - ERR06'  
            endif
        endif
        
        if (Me%Coupled%PorousMediaProperties .and. .not. Me%Coupled%PorousMedia) then
            write(*,*)'If using porous media properties also need porous media active '
            write(*,*)'Check basin data file in POROUS_MEDIA'
            stop 'VerifyOptions - ModuleBasin - ERR07'  
        endif            

        if (Me%Coupled%RunoffProperties .and. .not. Me%Coupled%Runoff) then
            write(*,*)'If using runoff properties also need runoff active '
            write(*,*)'Check basin data file in RUN_OFF'
            stop 'VerifyOptions - ModuleBasin - ERR08'  
        endif            

    end subroutine VerifyOptions

    !--------------------------------------------------------------------------
    
    subroutine ConstructSimpleInfiltration    
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        
        !
        call ConstructOneProperty (Me%SI%Ks,        "Ks",       "<BeginKS>",        "<EndKS>")
        call ConstructOneProperty (Me%SI%MP,        "MP",       "<BeginMP>",        "<EndMP>")
        call ConstructOneProperty (Me%SI%ThetaS,    "ThetaS",   "<BeginThetaS>",    "<EndThetaS>")
        call ConstructOneProperty (Me%SI%ThetaI,    "ThetaI",   "<BeginThetaI>",    "<EndThetaI>")
        call ConstructOneProperty (Me%SI%InfRate,   "InfRate",  "<BeginInfRate>",   "<EndInfRate>")
        call ConstructOneProperty (Me%SI%AccInf,    "AccInf",   "<BeginAccInf>",    "<EndAccInf>")
    
    end subroutine ConstructSimpleInfiltration    

    !--------------------------------------------------------------------------

    subroutine ConstructOneProperty (NewProperty, PropertyName, BlockBegin, BlockEnd)
    
        !Arguments-------------------------------------------------------------
        type (T_PropertyB)                          :: NewProperty
        character(Len=*)                            :: PropertyName, BlockBegin, BlockEnd

        !Local-----------------------------------------------------------------
        integer                                     :: ClientNumber
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine
        integer                                     :: STAT_CALL

        !call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)       
        !if (STAT_CALL /= SUCCESS_) stop 'ConstructSimpleInfiltration - ModuleBasin - ERR01'
        
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                   &
                                    BlockBegin, BlockEnd, BlockFound,                &
                                    FirstLine, LastLine, STAT_CALL) 
                                    
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)BlockBegin, BlockEnd
            write(*,*)BlockFound
            write(*,*)STAT_CALL
            stop 'ConstructSimpleInfiltration - ModuleBasin - ERR02'
        endif
                                    
        if (.not. BlockFound) then
            write(*,*)"Block : ",BlockBegin, " ", BlockEnd, " missing."
            stop 'ConstructOneProperty - ModuleBasin - ERR03'
        endif

        !Construct property ID
        NewProperty%ID%Name = PropertyName
        !call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call ConstructFillMatrix(PropertyID         = NewProperty%ID,                    &
                                 EnterDataID        = Me%ObjEnterData,                   &
                                 TimeID             = Me%ObjTime,                        &
                                 HorizontalGridID   = Me%ObjHorizontalGrid,              &
                                 ExtractType        = FromBlock,                         &
                                 PointsToFill2D     = Me%ExtVar%BasinPoints,             &
                                 Matrix2D           = NewProperty%Field,                 &
                                 TypeZUV            = TypeZ_,                            &
                                 STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneProperty - ModuleBasin - ERR04'

        call KillFillMatrix (NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneProperty - ModuleBasin - ERR05'

    end subroutine ConstructOneProperty

    !--------------------------------------------------------------------------

    subroutine ConstructDiffuseWaterSource
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ObjGD, i, j
        real, dimension(:,:), pointer           :: PopDensity
        

        !Starts a GridData on PopDensity
        ObjGD = 0
        call ConstructGridData  (ObjGD, Me%ObjHorizontalGrid, FileName = Me%PopDensityFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDiffuseWaterSource - ModuleBasin - ERR10'
        
        call GetGridData (ObjGD, PopDensity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDiffuseWaterSource - ModuleBasin - ERR20'

        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExtVar%RiverPoints(i, j) == 1) then
                Me%DiffuseFlow(i, j) = PopDensity(i, j) * Me%FlowPerCapita * Me%ExtVar%GridCellArea(i, j)
            endif
        enddo
        enddo
        
       
        call UnGetGridData (ObjGD, PopDensity, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDiffuseWaterSource - ModuleBasin - ERR30'
       
        
        call KillGridData (ObjGD, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDiffuseWaterSource - ModuleBasin - ERR40'

    end subroutine ConstructDiffuseWaterSource

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB   
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE

        !------------------------------------------------------------------------
       
        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        if (Me%Coupled%DrainageNetwork) then
            allocate(Me%OutPut%OutputChannels(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        endif

        Me%OutPut%NextOutPut = 1  
        
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%HDFFile)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleBasin - ERR02'

      
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleBasin - ERR02'


        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleBasin - ERR05'
        
        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunOff - ERR05'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleBasin - ERR07'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "RiverPoints", "-",          &
                              Array2D = Me%ExtVar%RiverPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleBasin - ERR07'

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleBasin - ERR08'       


    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSeries

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList        
        integer                                             :: STAT_CALL
        real, dimension(6), target                          :: AuxTime
        integer                                             :: i, AuxInt
        character(PathLength)                               :: AuxChar
        integer                                             :: TimeSerieNumber, dn, Id, Jd
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        character(len=StringLength)                         :: TimeSerieName


        !Time Serie of properties variable in the Basin    
        if (Me%Coupled%Evapotranspiration) then
            if (Me%Coupled%Vegetation) then
                if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then  
                    allocate(PropertyList(15))
                else
                    allocate(PropertyList(11))
                endif
            else
                allocate(PropertyList(8))
            endif        
        else
            if (Me%Coupled%Vegetation) then
                if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then  
                    allocate(PropertyList(14))
                else
                    allocate(PropertyList(10))
                endif
            else
                allocate(PropertyList(7))
            endif
        endif

        PropertyList(1)  = 'water column'
        PropertyList(2)  = 'water level'
        PropertyList(3)  = 'Infil. Rate [mm/hour]'
        PropertyList(4)  = 'Precipitation Rate [mm/hour]'
        PropertyList(5)  = 'Throughfall Rate [mm/hour]'
        PropertyList(6)  = 'EvapoTranspiration Rate [mm/hour]'
        PropertyList(7)  = 'Water Column Removed'
        i = 8
        if (Me%Coupled%Vegetation) then
            PropertyList(i)  = 'Potential Crop EVTP'
            i = i + 1
            PropertyList(i) = 'Canopy Capacity [m]'
            i = i + 1
            PropertyList(i) = 'Canopy Storage [m]'
            i = i + 1
        
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then        
                PropertyList(i) = 'Potential Evaporation [mm/h]'
                i = i + 1
                PropertyList(i) = 'Potential Transpiration [mm/h]'
                i = i +1
                PropertyList(i) = 'Actual Evaporation [mm/h]'
                i = i + 1
                PropertyList(i) = 'Actual Transpiration [mm/h]'
                i = i +1
            endif

        endif 


        if (Me%Coupled%Evapotranspiration) then
            PropertyList(i) = 'Reference Evapotranspiration [mm/h]'
        endif

        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            Me%Files%TimeSerieLocation,                                 &
                            PropertyList, "srb",                                        &
                            WaterPoints2D = Me%ExtVar%BasinPoints,                      &
                            ModelName = Me%ModelName,                                   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR01'

        deallocate(PropertyList)


        !Time Serie of properties variable integrated for the basin
        allocate(PropertyList(21))
        PropertyList(1)     = "Initial_Volume_Basin"
        PropertyList(2)     = "Initial_Volume_Vegetation"
        PropertyList(3)     = "Initial_Volume_PorousMedia"
        PropertyList(4)     = "Initial_Volume_Channels"
        PropertyList(5)     = "Evap_From_Vegetation"
        PropertyList(6)     = "Evap_From_Ground"
        PropertyList(7)     = "Evap_From_Soil"
        PropertyList(8)     = "Rain_In"
        PropertyList(9)     = "Out_Volume_Channel"
        PropertyList(10)    = "Out_Volume_Overland"
        PropertyList(11)    = "Final_Volume_Basin"
        PropertyList(12)    = "Final_Volume_Vegetation"
        PropertyList(13)    = "Final_Volume_PorousMedia"
        PropertyList(14)    = "Final_Volume_Channels"
        PropertyList(15)    = "Mass_Error_Global"
        PropertyList(16)    = "Mass_Error_Bacia"
        PropertyList(17)    = "Mass_Error_Porous"
        PropertyList(18)    = "Mass_Error_DNet"
        PropertyList(19)    = "Infiltration"
        PropertyList(20)    = "OL_Volume_To_Channels"
        PropertyList(21)    = "GW_Volume_To_Channels"


        call StartTimeSerie(Me%ObjTimeSerieBasin, Me%ObjTime,                           &
                            Me%Files%TimeSerieLocation,                                 &
                            PropertyList, "srb",                                        &
                            ResultFileName = 'Basin Mass Balance',                      &
                            STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR02'


        deallocate(PropertyList)

        
        call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2),         &
                                            AuxTime(3), AuxTime(4),         &
                                            AuxTime(5), AuxTime(6))
        

        AuxChar = Me%Files%TimeSerieLocation

        !Ouput of daily values
        if (Me%DailyFlow%On) then
        
            allocate(PropertyList(1)) 
            PropertyList(1) = "Daily Flow [m3]"
        
            AuxInt = 0
            call StartTimeSerie(AuxInt, Me%ObjTime,                      &
                                AuxChar,                                 &
                                PropertyList, "srb",                                        &
                                ResultFileName = 'Daily Flow',                            &
                                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR04'

            Me%DailyFlow%ObjTimeSerie = AuxInt
            
            deallocate(PropertyList)

            Me%DailyFlow%CurrentIndex = AuxTime(3)

                     
        endif

        !Ouput of monthly values
        if (Me%MonthlyFlow%On) then
        
            allocate(PropertyList(1)) 
            PropertyList(1) = "Monthly Flow [m3]"
        
            AuxInt = 0
            call StartTimeSerie(AuxInt, Me%ObjTime,                      &
                                AuxChar,                                 &
                                PropertyList, "srb",                                        &
                                ResultFileName = 'Monthly Flow',                              &
                                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR03'

            deallocate(PropertyList)

            Me%MonthlyFlow%ObjTimeSerie = AuxInt

            Me%MonthlyFlow%CurrentIndex = AuxTime(2)
                     
                     
        endif


        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR03'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR04'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR04'
            
i1:         if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR05'

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR06'

                    if (IgnoreOK) then
                        write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName),' - ',trim(Me%ModelName)
                        cycle
                    else
                        stop 'ConstructTimeSerie - ModuleBasin - ERR07'
                    endif

                endif


                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR08'

            endif i1

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI   = Id,                             &
                                      LocalizationJ   = Jd,                             & 
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleBasin - ERR09'

            if (Me%ExtVar%BasinPoints(Id, Jd) /= WaterPoint) then
                 write(*,*) 'Time Serie in a land cell - ',trim(TimeSerieName),' - ',trim(Me%ModelName)
            endif


        enddo



    end subroutine ConstructTimeSeries

    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                             :: ILB, IUB, JLB, JUB    
        integer                                             :: i, j

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        allocate(Me%WaterLevel              (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%WaterColumn             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%WaterColumnOld          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%WaterColumnRemoved      (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%CanopyDrainage          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        if (Me%Coupled%Vegetation) then
            allocate(Me%CanopyStorageCapacity   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
            allocate(Me%CanopyStorage           (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%CoveredFraction         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
            allocate(Me%CropEvapotrans          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                allocate(Me%PotentialTranspiration  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(Me%PotentialEvaporation    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            endif  
        endif      
        allocate(Me%ThroughFall             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%PotentialInfCol         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%FlowProduction          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%InfiltrationRate        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%PrecipRate              (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%ThroughRate             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%EVTPRate                (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AccInfiltration         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AccFlowProduction       (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AccEVTP                 (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AccRainFall             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AccEVPCanopy            (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%AccRainHour             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%RainStartTime           (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%RainDuration            (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%WaterColumnEvaporated   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        
       
        if (Me%Coupled%Snow) allocate(Me%SnowPack (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 

        allocate(Me%TimeSeriesBuffer    (21))
        allocate(Me%TimeSeriesBuffer2   (1))

        Me%WaterLevel               = FillValueReal
        Me%WaterColumn              = FillValueReal
        Me%WaterColumnOld           = FillValueReal
        
        if (Me%Coupled%Vegetation) then
            Me%CanopyStorageCapacity    = FillValueReal
            Me%CanopyStorage            = 0.0
            Me%CoveredFraction          = 0.0  !This variable was used before had value definition
            Me%CropEvapotrans           = 0.0
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                Me%PotentialTranspiration   = 0.0
                Me%PotentialEvaporation     = 0.0
            endif 
        endif       
        Me%ThroughFall              = 0.0
        Me%WaterColumnRemoved       = 0.0
        Me%CanopyDrainage           = 0.0
        Me%InfiltrationRate         = FillValueReal
        Me%PrecipRate               = FillValueReal
        Me%ThroughRate              = FillValueReal
        Me%EVTPRate                 = FillValueReal
        Me%EVTPRate2                = FillValueReal
        Me%AccInfiltration          = 0.0
        Me%AccFlowProduction        = 0.0
        Me%AccEVTP                  = 0.0
        Me%AccRainFall              = 0.0
        Me%AccEVPCanopy             = 0.0
        Me%PotentialInfCol          = 0.0
        Me%FlowProduction           = null_real
        Me%WaterColumnEvaporated    = 0.0

        
        if (Me%Coupled%Snow) Me%SnowPack = FillValueReal
        
        if (Me%Coupled%SimpleInfiltration) then
            allocate(Me%SI%Ks%Field        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%MP%Field        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%ThetaS%Field    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%ThetaI%Field    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%InfRate%Field   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%AccInf%Field    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%SI%Ks%Field                  = FillValueReal
            Me%SI%MP%Field                  = FillValueReal
            Me%SI%ThetaS%Field              = FillValueReal
            Me%SI%ThetaI%Field              = FillValueReal
            Me%SI%InfRate%Field             = FillValueReal
            Me%SI%AccInf%Field              = AllmostZero
         endif
         
        if (Me%DiffuseWaterSource) then
            allocate(Me%DiffuseFlow         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%DiffuseFlow            = null_real
        endif
         

        !Default Values
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                Me%WaterLevel(i, j) = Me%InitialWaterColumn + Me%ExtVar%Topography(i, j)
            endif

        enddo
        enddo

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructCoupledModules()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: VariableDT
!        integer                                     :: GeometryID 
        integer                                     :: MapID
        !Begin-----------------------------------------------------------------

        !Constructs Atmosphere
        if (Me%Coupled%Atmosphere) then
            call StartAtmosphere        (ModelName          = Me%ModelName,              &
                                         AtmosphereID       = Me%ObjAtmosphere,          &
                                         TimeID             = Me%ObjTime,                &
                                         GridDataID         = Me%ObjGridData,            &
                                         HorizontalGridID   = Me%ObjHorizontalGrid,      &
                                         MappingPoints      = Me%ExtVar%BasinPoints,     &
                                         STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR01'
        endif

        !Constructs Drainage Network
        if (Me%Coupled%DrainageNetwork) then

            call GetVariableDT(Me%ObjTime, VariableDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR03'
            
            call ConstructDrainageNetwork (ModelName         = Me%ModelName,                     &
                                           DrainageNetworkID = Me%ObjDrainageNetwork,            &
                                           TimeID            = Me%ObjTime,                       &
                                           Size              = Me%Size,                          &
                                           CheckMass         = .true.,                           &
                                           CoupledPMP        = Me%Coupled%PorousMediaProperties, &
                                           CoupledRP         = Me%Coupled%RunoffProperties,      &
                                           STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR04'
        endif

!        !Constructs RunOff
!        if (Me%Coupled%RunOff) then
!            call ConstructRunOff        (ModelName          = Me%ModelName,              &                 
!                                         RunOffID           = Me%ObjRunOff,              &
!                                         ComputeTimeID      = Me%ObjTime,                &
!                                         HorizontalGridID   = Me%ObjHorizontalGrid,      &
!                                         HorizontalMapID    = Me%ObjHorizontalMap,       &
!                                         GridDataID         = Me%ObjGridData,            &
!                                         BasinGeometryID    = Me%ObjBasinGeometry,       &
!                                         DrainageNetworkID  = Me%ObjDrainageNetwork,     &
!                                         STAT               = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR05'
!
!            !Constructs RunoffProperties 
!            if (Me%Coupled%RunoffProperties) then
!                call ConstructRunoffProperties        (ObjRunoffPropertiesID      = Me%ObjRunoffProperties,       &
!                                                       ComputeTimeID              = Me%ObjTime,                   &
!                                                       HorizontalGridID           = Me%ObjHorizontalGrid,         &
!                                                       HorizontalMapID            = Me%ObjHorizontalMap,          &
!                                                       BasinGeometryID            = Me%ObjBasinGeometry,          &
!                                                       RunoffID                   = Me%ObjRunoff,                 &
!                                                       GeometryID                 = Me%ObjGeometryID,             &
!                                                       STAT                       = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR07'
!            endif
!            
!        endif

        
        !Constructs PorousMedia
        if (Me%Coupled%PorousMedia) then
            call ConstructPorousMedia   (ModelName              = Me%ModelName,             &
                                         ObjPorousMediaID       = Me%ObjPorousMedia,        &
                                         ComputeTimeID          = Me%ObjTime,               &
                                         HorizontalGridID       = Me%ObjHorizontalGrid,     &
                                         HorizontalMapID        = Me%ObjHorizontalMap,      &
                                         TopographyID           = Me%ObjGridData,           &
                                         BasinGeometryID        = Me%ObjBasinGeometry,      &
                                         DrainageNetworkID      = Me%ObjDrainageNetwork,    &
                                         CheckGlobalMass        = .true.,                   &
                                         ConstructEvaporation   = Me%ConstructEvaporation,  &
                                         ConstructTranspiration = Me%ConstructTranspiration,&
                                         GeometryID             = Me%ObjGeometry,           &
                                         MapID                  = MapID,                    &
                                         STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR06'

               
            !Constructs PorousMediaProperties 
            if (Me%Coupled%PorousMediaProperties) then
                call ConstructPorousMediaProperties   (ObjPorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                                       ComputeTimeID              = Me%ObjTime,                   &
                                                       HorizontalGridID           = Me%ObjHorizontalGrid,         &
                                                       HorizontalMapID            = Me%ObjHorizontalMap,          &
                                                       BasinGeometryID            = Me%ObjBasinGeometry,          &
                                                       PorousMediaID              = Me%ObjPorousMedia,            &
                                                       GeometryID                 = Me%ObjGeometry,               &
                                                       MapID                      = MapID,                        &
                                                       CoupledDN                  = Me%Coupled%DrainageNetwork,   &
                                                       CheckGlobalMass            = .true.,                       &
                                                       STAT                       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR07'
            endif
        endif

        !Constructs RunOff
        if (Me%Coupled%RunOff) then
            call ConstructRunOff        (ModelName          = Me%ModelName,              &                 
                                         RunOffID           = Me%ObjRunOff,              &
                                         ComputeTimeID      = Me%ObjTime,                &
                                         HorizontalGridID   = Me%ObjHorizontalGrid,      &
                                         HorizontalMapID    = Me%ObjHorizontalMap,       &
                                         GridDataID         = Me%ObjGridData,            &
                                         BasinGeometryID    = Me%ObjBasinGeometry,       &
                                         DrainageNetworkID  = Me%ObjDrainageNetwork,     &
                                         STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR05'

            !Constructs RunoffProperties 
            if (Me%Coupled%RunoffProperties) then
                call ConstructRunoffProperties        (ObjRunoffPropertiesID      = Me%ObjRunoffProperties,       &
                                                       ComputeTimeID              = Me%ObjTime,                   &
                                                       HorizontalGridID           = Me%ObjHorizontalGrid,         &
                                                       HorizontalMapID            = Me%ObjHorizontalMap,          &
                                                       BasinGeometryID            = Me%ObjBasinGeometry,          &
                                                       RunoffID                   = Me%ObjRunoff,                 &
                                                       GridDataID                 = Me%ObjGridData,               &
                                                       InitialWaterColumn         = Me%InitialWaterColumn,        &
 !                                                      GeometryID                 = Me%ObjGeometry,               &
 !                                                      CoupledPMP                 = Me%Coupled%PorousMediaProperties, &
                                                       CoupledDN                  = Me%Coupled%DrainageNetwork,   &
                                                       STAT                       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR07'
            endif
            
        endif

        
        !Constructs Vegetation
        if (Me%Coupled%Vegetation) then
            call ConstructVegetation    (ObjVegetationID    = Me%ObjVegetation,              &
                                         TimeID             = Me%ObjTime,                    &
                                         GridDataID         = Me%ObjGridData,                &
                                         HorizontalGridID   = Me%ObjHorizontalGrid,          &
                                         HorizontalMapID    = Me%ObjHorizontalMap,           &
                                         AtmosphereID       = Me%ObjAtmosphere,              &
                                         PorousMediaID      = Me%ObjPorousMedia,             &
                                         MappingPoints      = Me%ExtVar%BasinPoints,         &
                                         GeometryID         = Me%ObjGeometry,                &
                                         BasinGeometryID    = Me%ObjBasinGeometry,           &
                                         CoupledAtmosphere  = Me%Coupled%Atmosphere,         &
                                         STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR08'
        endif

       
        !Constructs Simple Infiltration
        if (Me%Coupled%SimpleInfiltration) then
            call ConstructSimpleInfiltration ()
        endif
        
        !Constructs Diffuse Water Source
        if (Me%DiffuseWaterSource) then
            call ConstructDiffuseWaterSource ()
        endif
        
        !Constructs Output
        if (Me%Output%Yes) then
            call ConstructHDF5Output    ()
        endif

    end subroutine ConstructCoupledModules

    !--------------------------------------------------------------------------
    subroutine ConstructPropertyList

        !Local-----------------------------------------------------------------
        type (T_BasinProperty), pointer             :: NewProperty  => null()
        integer                                     :: ClientNumber
        integer                                     :: STAT_CALL, ErrorCount
        logical                                     :: BlockFound
        type (T_BasinProperty), pointer             :: PropertyX    => null()

        !----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    call ConstructProperty(NewProperty)

                    ! Add new Property to the Vegetation List 
                    call AddProperty(NewProperty)
                else
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR01'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructPropertyList - ModuleBasin - ERR02'
            end if cd1
        end do do1
        
        !Verifies existence of required properties
        if (Me%Coupled%Evapotranspiration) then

            ErrorCount = 0
            call SearchProperty(PropertyX, RefEvapotrans_        , .true., STAT = STAT_CALL)
            if (Me%Coupled%Evapotranspiration) ErrorCount = ErrorCount + STAT_CALL
            
            if (ErrorCount /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR03'
            
        endif
        
    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------

    !This subroutine reads all the information needed to construct a new property.           
    subroutine ConstructProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_BasinProperty), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !----------------------------------------------------------------------
             
        !Allocates new property
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructProperty - ModuleBasin - ERR01'

        !Construct property ID
        call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        !Construct property values
        call Construct_PropertyValues   (NewProperty)

        !----------------------------------------------------------------------

    end subroutine ConstructProperty    
    !--------------------------------------------------------------------------
    
    !This subroutine reads all the information needed to construct the property values       
    ! in the domain and in the boundaries            
    subroutine Construct_PropertyValues (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_BasinProperty),   pointer                 :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !----------------------------------------------------------------------

        !Fills Matrix
        allocate (NewProperty%Field (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        call ConstructFillMatrix(PropertyID         = NewProperty%ID,                    &
                                 EnterDataID        = Me%ObjEnterData,                   &
                                 TimeID             = Me%ObjTime,                        &
                                 HorizontalGridID   = Me%ObjHorizontalGrid,              &
                                 ExtractType        = FromBlock,                         &
                                 PointsToFill2D     = Me%ExtVar%BasinPoints,        &
                                 Matrix2D           = NewProperty%Field,                 &
                                 TypeZUV            = TypeZ_,                            &
                                 STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR01'

        call GetIfMatrixRemainsConstant(FillMatrixID    = NewProperty%ID%ObjFillMatrix,     &
                                        RemainsConstant = NewProperty%Constant,             &
                                        STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR02'

        !GR : changed the if condition because no constant default value in the wind would then be possible
        if (.not. NewProperty%ID%SolutionFromFile .and. .not. NewProperty%Constant) then
            call KillFillMatrix (NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR05'
        endif

        !----------------------------------------------------------------------

    end subroutine Construct_PropertyValues

    !--------------------------------------------------------------------------
    subroutine AddProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_BasinProperty),   pointer     :: NewProperty

        !----------------------------------------------------------------------
        type(T_BasinProperty),   pointer     :: PreviousProperty => null()
        type(T_BasinProperty),   pointer     :: CurrentProperty  => null()

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%FirstProperty    => NewProperty
        else
            PreviousProperty => Me%FirstProperty
            CurrentProperty  => PreviousProperty%Next
            do while (associated(CurrentProperty))
                PreviousProperty => CurrentProperty
                CurrentProperty  => PreviousProperty%Next
            enddo
            PreviousProperty%Next => NewProperty

        end if 

        !----------------------------------------------------------------------

    end subroutine AddProperty 

    !--------------------------------------------------------------------------

    subroutine ReadInitialFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: InitialFile
        type (T_Time)                               :: BeginTime, EndTimeFile, EndTime
        real                                        :: DT_error
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleBasin - ERR01'

        open(Unit = InitialFile, File = Me%Files%InitialFile, Form = 'UNFORMATTED',     &
             status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleBasin - ERR02'

        !Reads Date
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleBasin - ERR03'
        
        DT_error = EndTimeFile - BeginTime

        !Avoid rounding erros - Frank 08-2001
        if (abs(DT_error) >= 0.01) then
            
            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'Date in the file'
            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error', DT_error
            if (Me%StopOnWrongDate) stop 'ReadInitialFile - ModuleBasin - ERR04'   

        endif

        read(InitialFile)Me%WaterLevel
        if (Me%Coupled%Vegetation) then
            read(InitialFile)Me%CanopyStorage
        endif

        !Just reads the following values if the current run is not after a spin up period
        if (Me%StopOnWrongDate) then
            read(InitialFile)Me%AccInfiltration
            read(InitialFile)Me%AccEVTP
            read(InitialFile)Me%AccRainFall
            read(InitialFile)Me%AccEVPCanopy
            read(InitialFile, end=10)Me%AccFlowProduction 
        endif
        
   10   continue


        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleBasin - ERR05'  
      
    end subroutine ReadInitialFile


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyBasin(ObjBasinID, NewDT, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBasinID
        real                                        :: NewDT
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL, i, j
        character (Len = StringLength)              :: WarningString
        character (Len = StringLength)              :: LockToWhichModules
        character (Len = StringLength)              :: UnLockToWhichModules
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "ModifyBasin")

        STAT_ = UNKNOWN_

        call Ready(ObjBasinID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            call GetComputeCurrentTime  (Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBasin - ModuleBasin - ERR01'

            !Gets ExternalVars
            LockToWhichModules = 'AllModules'
            call ReadLockExternalVar (LockToWhichModules)

            !Verifies Global Mass
            if (Me%VerifyGlobalMass) then
                call CalculateMass(Me%MB%IniVolumeBasin,                           &
                                   Me%MB%IniVolumeVegetation,                      &
                                   Me%MB%IniVolumePorousMedia,                     &
                                   Me%MB%IniVolumeChannels)
            endif

            !Atmospheric Processes 
            if (Me%Coupled%Atmosphere) then

                call AtmosphereProcesses 

                !Actualizes the WaterColumn
                WarningString = 'AtmosphereProcesses'
                call ActualizeWaterColumn (WarningString) 
                
                if (Me%Coupled%RunoffProperties) then
                    call ActualizeWaterColumnConc(WarningString) 
                endif
                
            endif

            !Updates Vegetation 
            if (Me%Coupled%Vegetation) then

                call VegetationProcesses
            
            endif            

            !Porous Media
            if (Me%Coupled%PorousMedia) then

                call PorousMediaProcesses

                !Actualizes the WaterColumn
                WarningString = 'PorousMediaProcesses'
                call ActualizeWaterColumn (WarningString)                 
                
                
                if (Me%Coupled%PorousMediaProperties) then

                    call PorousMediaPropertiesProcesses     
                    call ActualizeWaterColumnConc(WarningString) 
                    
                endif           
               
            else
            
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then
                    
                        !Potencial Infiltration Column = ThroughFall
                        Me%WaterLevel (i, j) = Me%WaterLevel (i, j) + Me%Throughfall(i, j)
                    end if
                end do
                end do

                !Actualizes the WaterColumn
                WarningString = 'PorousMediaProcesses 2'
                call ActualizeWaterColumn (WarningString) 

                if (Me%Coupled%RunoffProperties) then
                    call ActualizeWaterColumnConc(WarningString) 
                endif
                
            endif
            

            !Simplified Infiltration / Evapotranspiration model
            if (Me%Coupled%SimpleInfiltration) then
                call SimpleInfiltration

                !Actualizes the WaterColumn
                WarningString = 'SimpleInfiltration'
                call ActualizeWaterColumn (WarningString) 
                
                if (Me%Coupled%RunoffProperties) then
                    call ActualizeWaterColumnConc(WarningString) 
                endif                
            endif
            
            !Overland Flow
            if (Me%Coupled%RunOff) then
                call OverLandProcesses

                !Actualizes the WaterColumn
                WarningString = 'OverLandProcesses'
                call ActualizeWaterColumn (WarningString) 

                if (Me%Coupled%RunoffProperties) then

                    call RunoffPropertiesProcesses     

                endif

            endif

            !Drainage Network
            if (Me%Coupled%DrainageNetwork) then
                call DrainageNetworkProcesses
           endif
            
            !Calibrating 1D
            if (Me%Calibrating1D) then
                call RemoveWatercolumn
            endif

            !Verifies Global Mass
            if (Me%VerifyGlobalMass) then
                call CalculateMass(Me%MB%FinalVolumeBasin,                           &
                                   Me%MB%FinalVolumeVegetation,                      &
                                   Me%MB%FinalVolumePorousMedia,                     &
                                   Me%MB%FinalVolumeChannels)
                !Verifies Mass and makes global balance...
                call GlobalMassBalance
            endif

            call TimeSerieOutput

            !HDF 5 Output
            if (Me%Output%Yes) then
                call HDF5OutPut       
            endif
            
            !Restart Output
            if (Me%Output%WriteRestartFile .and. .not. (Me%CurrentTime == Me%EndTime)) then
                if(Me%CurrentTime >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then
                    call WriteFinalFile
                    Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1
                endif
            endif

            !UnGets ExternalVars
            UnLockToWhichModules = 'AllModules'
            call ReadUnLockExternalVar (UnLockToWhichModules)
            call PredictNewDT(NewDT)

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "ModifyBasin")

    end subroutine ModifyBasin

    !--------------------------------------------------------------------------

    subroutine AtmosphereProcesses

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real, dimension(:, :), pointer              :: PrecipitationFlux
        type (T_BasinProperty), pointer             :: RefEvapotrans
        character (Len = StringLength)              :: LockToWhichModules
        character (Len = StringLength)              :: UnLockToWhichModules
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "AtmosphereProcesses")

        !Updates Rainfall
        call ModifyAtmosphere       (Me%ObjAtmosphere, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR01'

        !Gets Rainfall [m3/s]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, PrecipitationFlux, ID = Precipitation_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR02'

        if (Me%Coupled%Vegetation) then 
            
            !Read Vegetation Properties
            LockToWhichModules = 'Atmosphere'
            call ReadLockExternalVar (LockToWhichModules)
        
        endif

        !Divides Precipitation into Snow, Rain, Throughfall and Canopy storage
        call DividePrecipitation(PrecipitationFlux)


        if (Me%Coupled%Evapotranspiration) then
            call SearchProperty(RefEvapotrans, RefEvapotrans_        , .true., STAT = STAT_CALL)        

            if(.not. RefEvapotrans%Constant) then
                call CalcPotEvapoTranspiration
            endif   
        endif

        if (Me%Coupled%Vegetation) then 
            
            !UnlockVegetation Properties so they can be run when vegetation processes are called
            UnLockToWhichModules = 'Atmosphere'
            call ReadUnLockExternalVar (UnLockToWhichModules)

        endif

        call UnGetAtmosphere    (Me%ObjAtmosphere, PrecipitationFlux, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR03'

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "AtmosphereProcesses")
    
    end subroutine AtmosphereProcesses

    !--------------------------------------------------------------------------

    subroutine DividePrecipitation (PrecipitationFlux)
    
        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer              :: PrecipitationFlux

        !Local-----------------------------------------------------------------     
        integer                                     :: i, j, STAT_CALL
        real, dimension(:,:), pointer               :: AirTemperature
        real                                        :: GrossPrecipitation
        real                                        :: SnowInput
        real                                        :: CanopyDrainage
        real                                        :: CurrentFlux, Rand, SecondsPassed
        real, dimension(6), target                  :: AuxTime
        integer, save                               :: LastHour = -99
        logical                                     :: ChangeRain
        real(8)                                     :: NewVolumeOnLeafs, OldVolumeOnLeafs

        if (Me%Coupled%Snow) then
            !Gets Air Temperature [ºC]
            call GetAtmosphereProperty  (Me%ObjAtmosphere, AirTemperature, ID = AirTemperature_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DividePrecipitation - ModuleBasin - ERR10'
        end if


        !Verifies if the hourly average rainfall is to be concentrated in a shorted time period...
        !Only do this if DT is already shorten by DTDuringRain variable.
        if (Me%CurrentDT <= Me%DTDuringRain .and. Me%ConcentrateRain) then
            ChangeRain = .true.
        else
            ChangeRain = .false.
        endif


        if (ChangeRain) then

            !Gets time
            call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2),         &
                                                AuxTime(3), AuxTime(4),         &
                                                AuxTime(5), AuxTime(6))
                                                
            SecondsPassed = AuxTime(5)*60.0 +  AuxTime(6)
                                                
            !If hour changed set accumulated rain in hour to zero
            if (nint(AuxTime(4)) /= LastHour .and. SecondsPassed > 0.0) then
                LastHour       = nint(AuxTime(4))
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then
                        Me%AccRainHour  (i, j) = 0.0

                        !RainDuration
                        call random_number(Rand)
                        Me%RainDuration(i, j)  = Rand * Me%RainAverageDuration + Me%RainAverageDuration / 2.0

                        !RainStartTime
                        call random_number(Rand)
                        Me%RainStartTime(i, j) = Rand * (3600. - Me%RainDuration(i, j) - Me%DTDuringRain)
                    endif
                enddo
                enddo
                
                
            endif
            
        endif
        

        Me%MB%RainIn = 0.0
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then


                if (ChangeRain) then
                
                    !At the beginning there is no rain
                    if (SecondsPassed < Me%RainStartTime(i, j)) then
                    
                        CurrentFlux         = 0.0
                
                    !Rain during the RainDuration
                    elseif (SecondsPassed >= Me%RainStartTime(i, j) .and. &
                            SecondsPassed + Me%CurrentDT <  Me%RainStartTime(i, j) + Me%RainDuration(i, j)) then
                        
                        !CurrentFlux
                        !m3/s               = m3/s                    * s      / s
                        CurrentFlux         = PrecipitationFlux(i, j) * 3600.0 / Me%RainDuration(i, j)
                        
                    else
                    
                        !Sets flux so total rain in hour will be correct
                        !m3/s               =   (m3/s * s - m3) / s
                        CurrentFlux         =  (PrecipitationFlux(i, j) * 3600. - Me%AccRainHour(i, j)) / &
                                                Me%CurrentDT
                    endif
                    
                    !Accumalted during period
                    !m3                 = m3                   + m3/s * s
                    Me%AccRainHour(i, j)= Me%AccRainHour(i, j) + CurrentFlux * Me%CurrentDT

                else
                
                    CurrentFlux = PrecipitationFlux(i, j)
                
                endif
                
                !Gross Rain 
                !m                 = m3/s                    * s            /  m2
                GrossPrecipitation = CurrentFlux * Me%CurrentDT / Me%ExtVar%GridCellArea(i, j)

                !Precipitation Rate for (output only)
                !mm/ hour            m3/s                    / m2                           * mm/m s/h
                Me%PrecipRate  (i,j) = CurrentFlux / Me%ExtVar%GridCellArea(i, j) * 1000.0 * 3600.0

                !Accumulated rainfall
                !m
                Me%AccRainfall(i, j) = Me%AccRainfall (i, j) + GrossPrecipitation

                !Integrates Total Input Volume
                if (Me%VerifyGlobalMass) then
                    Me%MB%RainIn = Me%MB%RainIn + GrossPrecipitation * Me%ExtVar%GridCellArea(i, j)
                endif

                if (Me%Coupled%Snow) then
                
                    !Divides Precipitation into rain / snow
                    !Taken from Daisy describtion
                    if      (AirTemperature(i, j) < -2.0) then
                        SnowInput = 0.0 !Smow module not yet active
                        !SnowInput = GrossPrecipitation
                    elseif  (AirTemperature(i, j) <  2.0) then
                        SnowInput = 0.0 !Smow module not yet active
                        !SnowInput = (2.0 - AirTemperature(i, j)) / 4.0 * GrossPrecipitation
                    else
                        SnowInput = 0.0
                    endif
                
                    !Rain is all what remains after snow input
                    GrossPrecipitation = GrossPrecipitation - SnowInput
                
                    !Increase SnowPack
                    Me%SnowPack (i, j) = Me%SnowPack (i, j) + SnowInput
                end if
                
                if (Me%Coupled%Vegetation) then
                
                    !Volume on Leaf (before leaf grow)
                    OldVolumeOnLeafs = Me%CoveredFraction(i, j) * Me%CanopyStorage(i, j) * Me%ExtVar%GridCellArea(i, j)
                    
                    !Calculates Covered Fraction
                    Me%CoveredFraction(i, j) = 1.0 - exp(-Me%ExtinctionCoef * Me%ExtVar%LeafAreaIndex(i, j))
                    
                    !Recalculates Canopy Storage Capacity (for the case that LAI is variable in time)
                    Me%CanopyStorageCapacity(i, j) = Me%ExtVar%SpecificLeafStorage(i, j) * Me%ExtVar%LeafAreaIndex(i, j)
                    
                    !Throughfall on uncovered area
                    !m                    
                    Me%ThroughFall(i, j)    = GrossPrecipitation * (1.0 - Me%CoveredFraction(i, j))
                    
                    !New volumes on leafs is the old volume plus the gross rain on the covered fraction of the cell
                    NewVolumeOnLeafs = OldVolumeOnLeafs + GrossPrecipitation * Me%CoveredFraction(i, j)          &
                                        * Me%ExtVar%GridCellArea(i, j)

                    !Put all into the canopy, excess will be removed 10 lines downwards...
                    if (Me%CoveredFraction(i, j) > 0.0) then
                        Me%CanopyStorage(i, j) = NewVolumeOnLeafs / (Me%CoveredFraction(i, j) * Me%ExtVar%GridCellArea(i, j))
                    else
                        !LAI got igual to zero... put everythink to throughfall
                        Me%ThroughFall  (i, j) = Me%ThroughFall(i, j) + NewVolumeOnLeafs / Me%ExtVar%GridCellArea(i, j)
                        Me%CanopyStorage(i, j) = 0.0
                    endif

                    !Calculates CanopyDrainage so that CanopyStorage as maximum is full
                    if (Me%CanopyStorage(i, j) > Me%CanopyStorageCapacity(i, j)) then
                        CanopyDrainage         = Me%CanopyStorage(i, j) - Me%CanopyStorageCapacity(i, j)
                        Me%CanopyStorage(i, j) = Me%CanopyStorageCapacity(i, j)
                    else
                        CanopyDrainage         = 0.0
                    endif

                    !Adds Canopy drainage to Throughfall considering the covered area
                    !m
                    Me%ThroughFall(i, j)   = Me%ThroughFall(i, j) + CanopyDrainage * Me%CoveredFraction(i, j)
                    Me%CanopyDrainage(i,j) = CanopyDrainage * Me%CoveredFraction(i, j)
            
                else
                    
                    Me%ThroughFall(i, j)    = GrossPrecipitation
                    Me%CanopyDrainage(i,j)  = 0.0
                endif
                    
                !mm/ hour                   m                       s         mm/m     s/h
                Me%ThroughRate(i, j) = Me%ThroughFall(i, j) / Me%CurrentDT * 1000.0 * 3600.0
                
                          
            endif
        enddo
        enddo

        if (Me%Coupled%Snow) then
            !UnGets Air Temperature [ºC]
            call UnGetAtmosphere  (Me%ObjAtmosphere, AirTemperature,    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DividePrecipitation - ModuleBasin - ERR90'
        end if


    end subroutine DividePrecipitation

    !--------------------------------------------------------------------------

    subroutine CalcPotEvapoTranspiration

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------     
        integer                                     :: i, j, STAT_CALL
        real                                        :: NetRadiation
        real                                        :: SSVPC, psiconst
        real                                        :: LwradCorrection, Lwrad
        real                                        :: SVP, VP, SoilHeatFluxDensity
        real, dimension(:,:), pointer               :: SolarRadiation
        real, dimension(:,:), pointer               :: WindModulus
        real, dimension(:,:), pointer               :: AirTemperature
        real, dimension(:,:), pointer               :: RelativeHumidity
        real, dimension(:,:), pointer               :: ATMTransmitivity
        real,    parameter                          :: LatentHeatOfVaporization = 2.5e6         ![J/kg]
        real,    parameter                          :: ReferenceDensity         = 1000.         ![kg/m3]
        real                                        :: LatentHeat_
        real(8)                                     :: EvaporationRate, dH
        type(T_BasinProperty),pointer               :: RefEvapotrans


       
        !Gets Horizontal Sun Radiation [W/m2]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, SolarRadiation, ID = SolarRadiation_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR01'

        !Gets ATMTransmitivity
        call GetAtmosphereProperty  (Me%ObjAtmosphere, ATMTransmitivity, ID = AtmTransmitivity_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR02'

        !Gets Wind Modulus [m/s]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, WindModulus, ID = WindModulus_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR03'

        !Gets Air Temperature [ºC]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, AirTemperature, ID = AirTemperature_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR04'

        !Gets Air Temperature [ºC]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, RelativeHumidity, ID = RelativeHumidity_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR05'

        call SearchProperty(RefEvapotrans, RefEvapotrans_        , .true., STAT = STAT_CALL)        

        if (RefEvapotrans%ID%SolutionFromFile) then
            call ModifyFillMatrix (FillMatrixID   = RefEvapotrans%ID%ObjFillMatrix,          &
                                   Matrix2D       = RefEvapotrans%Field,                     &
                                   PointsToFill2D = Me%ExtVar%BasinPoints,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR01'
        endif 
        
        !Calculates evaporation from canopy / Watercolumn on the ground
        Me%MB%EvapFromVegetation = 0.0
        Me%MB%EvapFromGround     = 0.0
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then

                ![J/m2/s]
                LatentHeat_ = -1.0 * LatentHeat (ReferenceDensity, AirTemperature(i, j), AirTemperature(i, j),     &
                                                 RelativeHumidity(i, j), WindModulus(i,j))
                                
                ![m/s]          = [J/m2/s] / [J/kg] / [kg/m3] 
                EvaporationRate = LatentHeat_ / LatentHeatOfVaporization / ReferenceDensity 
                
                if (Me%Coupled%Vegetation) then
                    !dH
                    dH              = min(dble(EvaporationRate * Me%CurrentDT), Me%CanopyStorage(i, j))
                
                    !Accumulated EVAP from Canopy
                    Me%AccEVPCanopy(i, j) = Me%AccEVPCanopy(i, j) + dH
                
                    !New Canopy Storage
                    Me%CanopyStorage(i, j) = Me%CanopyStorage(i, j) - dH
                
                    !Sistem loss
                    if (Me%VerifyGlobalMass) then
                        Me%MB%EvapFromVegetation = Me%MB%EvapFromVegetation +                         &
                                                   dH * Me%ExtVar%GridCellArea(i, j) * Me%CoveredFraction(i, j)
                    endif

                endif
                
                !Also EVAP from watercolumn - important for 1D cases to avoid accumulation of Water on the surface
                if (Me%WaterLevel (i, j) > Me%ExtVar%Topography(i, j)) then
                    !dH
                    dH              = min(dble(EvaporationRate * Me%CurrentDT), Me%WaterLevel (i, j) - Me%ExtVar%Topography(i, j))

                    Me%AccEVTP(i, j) = Me%AccEVTP(i, j) + dH
                    
                    !New Water Level
                    Me%WaterLevel(i, j) = Me%WaterLevel(i, j) - dH

                    !Evaporation from Ground
                    if (Me%VerifyGlobalMass) then
                        Me%MB%EvapFromGround = Me%MB%EvapFromGround + dH * Me%ExtVar%GridCellArea(i, j)
                    endif
                    
                    Me%WaterColumnEvaporated(i, j) = dH
                        
                else
                    
                    Me%WaterColumnEvaporated(i, j) = 0.0
                        
                endif
                
            endif
        enddo
        enddo
        
               

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then


etr_fao:        if (.not. RefEvapotrans%ID%SolutionFromFile) then
            
                    !Calculate Psicrometric constant
                    !Calculation of the atmospheric pressure based on the heigth simplification of the ideal gas law
                    psiconst    = 0.665E-3 * (101.3 * ( (293.-0.0065 * Me%ExtVar%Topography(i,j)) / 293.) **5.26) ![kPa / ºC]
                
                    !Saturation Vapour pressura [kPa - Eqn 11]
                    SVP         = 0.6108 * exp (17.27 * AirTemperature(i,j) / (AirTemperature(i,j) + 237.3))
                
                    ![kPa]
                    VP          = SVP * RelativeHumidity(i, j)  

                    !Calculates the Slope [kPa/ºC]
                    SSVPC       =  4098.* SVP / (AirTemperature(i,j) + 237.3)**2.0 ![kPa / ºC]
                
                    !StefanBoltzmann          = 5.669e-08     ![W/m2/K4]
                    LwradCorrection =   (0.34 - 0.14 * VP **(0.5)) * (1.35 * ATMTransmitivity(i,j)  - 0.35)
                    Lwrad           =   5.669e-08 * (AirTemperature(i,j) + 273.15)** 4. * LwradCorrection   ![W / m2]
                
                    !Calculation of net radiation (0.23 is the reference albedo)
                    NetRadiation    = (1-0.23) * SolarRadiation(i, j) - Lwrad            ![W / m2]   
                
                    !Converts Netradiation into MJ/m2/hour (requested by the formular below)
                    !1W = J / s => 1MJ/hour
                    NetRadiation    = NetRadiation /1.e6 * 3600.
                
                    if (NetRadiation .GE. 0) then
                        SoilHeatFluxDensity = 0.1 * NetRadiation
                    else
                        SoilHeatFluxDensity = 0.5 * NetRadiation
                    end if

                
                    !FAO ET0 - mm/hour
                    !http://www.fao.org/docrep/X0490E/x0490e08.htm#calculation%20procedure (Hourly time step)
                    RefEvapotrans%Field(i, j) = (0.408 * SSVPC * (NetRadiation - SoilHeatFluxDensity) +  &
                                            psiconst * 37. /(AirTemperature(i, j) + 273.)  *        &
                                            WindModulus(i,j) * (SVP- VP)) /                         & 
                                            (SSVPC + psiconst * (1. + 0.34 * WindModulus(i,j)))
                endif etr_fao
                                                            
                !m/s - Porous media consistency
                RefEvapotrans%Field(i, j)  = max(RefEvapotrans%Field(i, j) / 1000. / 3600., 0.0)

                if (Me%Coupled%Vegetation) then
                    
                    !m/s
                    Me%CropEvapotrans(i, j) = RefEvapotrans%Field(i, j) * Me%ExtVar%CropCoefficient(i, j)

                    if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                        
                        !m/s
                        Me%PotentialTranspiration(i, j) = Me%CropEvapotrans(i, j)                      &
                                                         * ( 1.0 - exp(-0.463 * Me%ExtVar%LeafAreaIndex(i, j)))
                
                        Me%PotentialEvaporation  (i, j) = Me%CropEvapotrans(i, j) - Me%PotentialTranspiration(i, j)

                    endif
                
                endif

            endif   

        enddo
        enddo
        
         !Gets Horizontal Sun Radiation [W/m2]
        call UnGetAtmosphere  (Me%ObjAtmosphere, SolarRadiation,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR06'

        !Gets Horizontal Sun Radiation [W/m2]
        call UnGetAtmosphere  (Me%ObjAtmosphere, ATMTransmitivity,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR07'

        !Gets Wind Modulus [m/s]
        call UnGetAtmosphere  (Me%ObjAtmosphere, WindModulus,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR08'

        !Gets Air Temperature [ºC]
        call UnGetAtmosphere  (Me%ObjAtmosphere, AirTemperature,    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR09'

        !Gets Air Temperature [ºC]
        call UnGetAtmosphere  (Me%ObjAtmosphere, RelativeHumidity,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR10'

    end subroutine CalcPotEvapoTranspiration

    !------------------------------------------------------------------------------
   
    real function CalcRadTerm(Radiation, SoilLoss, WindVelocity, SSVPC, psiconst)

        !Arguments-------------------------------------------------------------
        real                                        :: Radiation                ![W / m2]
        real                                        :: WindVelocity             ![m / s]
        real                                        :: SSVPC                    ![kPa / ºC]
        real                                        :: psiconst                 ![]
        real                                        :: SoilLoss
                
        !Local-----------------------------------------------------------------

        CalcRadTerm = 0.408 * SSVPC * 36./10000. * (Radiation - SoilLoss)                  & 
                      /(SSVPC + psiconst * (1. + 0.34 * WindVelocity))           ![mm / h]       

    end function CalcRadTerm
   

    !--------------------------------------------------------------------------

    real function CalcAeroTerm(Temperature, Wind, SSVPC, SVP, VP, psiconst)

        !Arguments-------------------------------------------------------------
        real                                        :: Temperature  ![ºC]
        real                                        :: Wind         ![m / s]
        real                                        :: SSVPC        ![kPa / ºC]
        real                                        :: psiconst     ![kPa / ºC]
        real                                        :: SVP          ![kPa]
        real                                        :: VP           ![kPa]      
        
        !Local-----------------------------------------------------------------
        real                                        :: R   = 0.287 ![kJ / kg / K]
        real                                        :: RMW = 0.622 !molecular weight ratio WaterVapor/DryAir  
        
        CalcAeroTerm =  60 * 60 * RMW /(1.01 * (Temperature+273) * R * 208) * psiconst * Wind * (SVP - VP )       &
                        / (SSVPC + psiconst * (1 + 0.34 * Wind) )  ![mm / h]         

    end function CalcAeroTerm
   
    !--------------------------------------------------------------------------

    subroutine VegetationProcesses

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real, dimension(:,:,:), pointer             :: ActualTranspiration
        real, dimension(:,:  ), pointer             :: PotentialTranspiration
        real, dimension(:,:,:), pointer             :: Nitrate
        real, dimension(:,:,:), pointer             :: InorganicPhosphorus
        logical                                     :: ModelNitrogen
        logical                                     :: ModelPhosphorus

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "VegetationProcesses")

        
        call GetVegetationOptions (Me%ObjVegetation,                                     &
                                   ModelNitrogen = ModelNitrogen,                        &
                                   ModelPhosphorus = ModelPhosphorus,                    &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'VegetationProcesses - ModuleBasin - ERR01'

        if (ModelNitrogen) then 
            
            if (Me%Coupled%PorousMediaProperties) then
            
                call GetPMPConcentration(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                         ConcentrationX          = Nitrate,                      &
                                         PropertyXIDNumber       = Nitrate_,                     &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                    write (*,*) 'Trying to model nitrogen in vegetation but nitrate property not'
                    write (*,*) 'defined in porousmediaproperties. Check options'                
                    stop 'VegetationProcesses - ModuleBasin - ERR010'
                endif

                call SetSoilConcVegetation (Me%ObjVegetation,                                 &
                                            Nitrate             = Nitrate,                    &
                                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR020' 


                call UnGetPorousMediaProperties (Me%ObjPorousMediaProperties, Nitrate, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR025'   
                  
            else
                write (*,*) 'Can not model nitrogen in vegetation if porous media properties'
                write (*,*) 'model not connected. Check basin keyword.'
                stop 'VegetationProcesses - ModuleBasin - ERR030'
            endif
        endif

        if (ModelPhosphorus) then 

            if (Me%Coupled%PorousMediaProperties) then

                call GetPMPConcentration(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                         ConcentrationX          = InorganicPhosphorus,          &
                                         PropertyXIDNumber       = Inorganic_Phosphorus_,        &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                    write (*,*) 'Trying to model phosphorus in vegetation but inorganic phosphorus '
                    write (*,*) 'property not defined in porousmediaproperties. Check options'                
                    stop 'VegetationProcesses - ModuleBasin - ERR040'
                endif

                call SetSoilConcVegetation (Me%ObjVegetation,                                 &
                                            InorganicPhosphorus = InorganicPhosphorus,        &
                                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR050'   

                call UnGetPorousMediaProperties (Me%ObjPorousMediaProperties, InorganicPhosphorus, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR055'   

            else
                write (*,*) 'Can not model phosphorus in vegetation if porous media properties'
                write (*,*) 'model not connected. Check basin keyword.'
                stop 'VegetationProcesses - ModuleBasin - ERR060'
            endif
        endif

        
        !Transpiration
        if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
            PotentialTranspiration => Me%PotentialTranspiration
        else
            PotentialTranspiration => Me%CropEvapoTrans
        endif                
        
        call ModifyVegetation(ObjVegetationID        = Me%ObjVegetation,      &
                              MappingPoints          = Me%ExtVar%BasinPoints, &
                              PotentialTranspiration = PotentialTranspiration,&
                              ActualTranspiration    = ActualTranspiration,   & 
                              STAT                   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR070'
    
        !Points to computed variables to be used in porous media (vegetation uses ModulePorousMedia)
        Me%ExtVar%ActualTranspiration => ActualTranspiration


        if (MonitorPerformance) call StopWatch ("ModuleBasin", "VegetationProcesses")
    
    end subroutine VegetationProcesses

    !--------------------------------------------------------------------------
    
    subroutine SimpleInfiltration
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: I, J
!        integer                                     :: TID, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

        integer, parameter                          :: ChunkSize = 10
        integer                                     :: Chunk
        
        if (MonitorPerformance) call StartWatch ("ModuleBasin", "SimpleInfiltration")
        CHUNK = Me%WorkSize%JUB / 10
        
!$OMP PARALLEL PRIVATE(I,J)
        
!        TID = OMP_GET_THREAD_NUM()
!        IF (TID .EQ. 0) THEN
!            NTHREADS = OMP_GET_NUM_THREADS()
!            PRINT *, 'Number of threads =', NTHREADS
!        END IF
!        PRINT *, 'Thread',TID,' starting...'
        
        !Updates Water column
        
!$OMP DO SCHEDULE(DYNAMIC)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                if (Me%Watercolumn(i, j) > 0.0) then

                    !Infiltration Rate
                    Me%SI%InfRate%Field(i, j) = Me%SI%Ks%Field(i, j) + Me%SI%Ks%Field(i, j) * (Me%SI%MP%Field(i, j) * &
                                                (Me%SI%ThetaS%Field(i, j) - Me%SI%ThetaI%Field(i, j))) /              &
                                                Me%SI%AccInf%Field(i, j)
                    
                    if (Me%SI%InfRate%Field(i, j) > Me%Watercolumn(i, j) / Me%CurrentDT) then
                        Me%SI%InfRate%Field(i, j) = Me%Watercolumn(i, j) / Me%CurrentDT
                    endif

                    Me%SI%AccInf%Field(i, j) = Me%SI%AccInf%Field(i, j) + Me%SI%InfRate%Field(i, j) * Me%CurrentDT

                else

                    !Sets Infiltration Rate to zero
                    Me%SI%InfRate%Field(i, j) = 0.0
                
                    !Resets accumulated infiltration
                    Me%SI%AccInf%Field(i, j) = AllmostZero
                endif
                
                
                Me%WaterLevel       (i, j) = Me%WaterLevel (i, j) - Me%SI%InfRate%Field(i, j) * Me%CurrentDT
                !mm /hour
                Me%InfiltrationRate (i, j) = Me%SI%InfRate%Field(i, j) * 1000.0 * 3600.0
                !mm /hour
                Me%EVTPRate         (i, j) = 0.0
                !m
                Me%AccInfiltration  (i, j) = Me%AccInfiltration  (i, j) + Me%SI%InfRate%Field(i, j) * Me%CurrentDT
                !m
                Me%AccEVTP          (i, j) = 0.0
                
            endif
            
!            WRITE(*,*) TID,I,J
            
        enddo
       enddo
 !$OMP END DO NOWAIT
      !PRINT *, 'Thread',TID,' done.'

!$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "SimpleInfiltration")

    
    end subroutine SimpleInfiltration

    !--------------------------------------------------------------------------

    subroutine OverLandProcesses

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :), pointer              :: FlowX, FlowY
        real, dimension(:, :), pointer              :: OLFlowToChannels
        real, dimension(:, :), pointer              :: FlowAtBoundary
        integer                                     :: STAT_CALL
        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        character (Len = StringLength)              :: UnLockToWhichModules
        character (Len = StringLength)              :: LockToWhichModules
        !Begin-----------------------------------------------------------------
            
        if (MonitorPerformance) call StartWatch ("ModuleBasin", "OverLandProcesses")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Unlocks External Var for Module Runoff to change Computefaces
        UnLockToWhichModules = 'AllModules'
        call ReadUnLockExternalVar (UnLockToWhichModules)

        call ModifyRunOff   (Me%ObjRunOff, Me%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR01'

        !Re-Locks External 
        LockToWhichModules = 'AllModules'
        call ReadLockExternalVar (LockToWhichModules)
        
        call GetOverLandFlow    (Me%ObjRunOff, FlowX, FlowY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR02'

        call GetFlowToChannels  (Me%ObjRunOff, OLFlowToChannels, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR03'

        call GetFlowAtBoundary  (Me%ObjRunOff, FlowAtBoundary, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR04'

        !Updates Water column
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then

                Me%WaterLevel (i, j) = Me%WaterLevel (i, j)                     &
                                     + ( FlowX(i, j) - FlowX(i, j+1)            &
                                     +   FlowY(i, j) - FlowY(i+1, j)            &
                                     -  OLFlowToChannels(i, j) -                &
                                        FlowAtBoundary(i, j)      )             &
                                     * Me%CurrentDT /                           &
                                     Me%ExtVar%GridCellArea(i, j)

            endif
        enddo
        enddo
        
        call UnGetRunOff (Me%ObjRunOff, FlowX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR05'

        call UnGetRunOff (Me%ObjRunOff, FlowY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR06'

        call UnGetRunOff            (Me%ObjRunOff, OLFlowToChannels, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR07' 

        call UnGetRunOff            (Me%ObjRunOff, FlowAtBoundary, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR08' 

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "OverLandProcesses")

    end subroutine OverLandProcesses    

    !--------------------------------------------------------------------------

    subroutine DrainageNetworkProcesses

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:, :), pointer              :: OLFlowToChannels
        real, dimension(:, :), pointer              :: GWFlowToChannels
        real, dimension(:, :, :), pointer           :: GWFlowToChannelsLayer
        integer, dimension(:, :), pointer           :: GWFlowBottomLayer
        integer, dimension(:, :), pointer           :: GWFlowTopLayer
        real, dimension(:, :), pointer              :: SolarRadiation
        real, dimension(:, :), pointer              :: AirTemperature
        real, dimension(:, :), pointer              :: CloudCover
        real, dimension(:, :), pointer              :: RelativeHumidity
        real, dimension(:, :), pointer              :: WindSpeed
        integer                                     :: ILB, IUB, JLB, JUB, i, j, k
        integer                                     :: STAT_CALL
        logical                                     :: NeedAtmosphere, NeedsRadiation, PropAdvDiff
        integer                                     :: nProperties, iProp 
        integer                                     :: PropID, GW_Link
        real, dimension(:, :   ), pointer           :: PMPConcentration2D
        real, dimension(:, :, :), pointer           :: PMPConcentration
        real, dimension(:, :   ), pointer           :: RPConcentration
        integer, dimension(:, :), pointer           :: GWLayer
        

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "DrainageNetworkProcesses")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        call GetFlowToChannels  (Me%ObjRunOff, OLFlowToChannels, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01'


        nullify(GWFlowToChannels)
        if (Me%Coupled%PorousMedia) then
        
            call GetGWFlowOption (ObjPorousMediaID    = Me%ObjPorousMedia,                  &
                                  DrainageNetworkLink = GW_Link,                            &
                                  STAT                = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01a'        
    

            call GetGWFlowToChannels    (Me%ObjPorousMedia, GWFlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01b'
            
            if (GW_Link == Layer_) then 

                call GetGWFlowToChannelsByLayer    (Me%ObjPorousMedia, GWFlowToChannelsLayer, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01c'
                
                call GetGWToChannelsLayers (ObjPorousMediaID        = Me%ObjPorousMedia,        &
                                            GWToChannelsBottomLayer = GWFlowBottomLayer,        &
                                            GWToChannelsTopLayer    = GWFlowTopLayer,           &
                                            STAT                    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01d' 
                
                call SetGWFlowLayersToDN  (DrainageNetworkID    = Me%ObjDrainageNetwork,        &
                                           GWFlowBottomLayer    = GWFlowBottomLayer,            &
                                           GWFlowTopLayer       = GWFlowTopLayer,               &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01e' 
                                           
            endif
            
        endif
        
        Me%MB%OLFlowToRiver = 0.0
        Me%MB%GWFlowToRiver = 0.0
        if (Me%VerifyGlobalMass) then
            do j = JLB, JUB
            do i = ILB, IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    Me%MB%OLFlowToRiver = Me%MB%OLFlowToRiver + OLFlowToChannels(i, j) * Me%CurrentDT
                    if (Me%Coupled%PorousMedia) then
                        Me%MB%GWFlowToRiver = Me%MB%GWFlowToRiver + GWFlowToChannels(i, j) * Me%CurrentDT
                    endif
                endif
            enddo
            enddo        
            
        endif
        
        if (Me%Coupled%PorousMediaProperties) then
        
            call GetPMPnProperties (Me%ObjPorousMediaProperties, nProperties, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR02'

            do iProp = 1, nProperties

                call GetPMPPropertiesIDByIdx(PorousMediaPropertiesID = Me%ObjPorousMediaProperties, &
                                             Idx                     = iProp,                       &
                                             ID                      = PropID,                      &
                                             PropAdvDiff             = PropAdvDiff,                 &
                                             STAT                    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR03' 
                
                !Only not particulate properties (in water phase)with advection diffusion may interact with drainage network
                !Particulate prop do not have advection diffusion in porous media so does not need double check
                if (PropAdvDiff) then
                    
                    !Get the conc from PMP
                    call GetPMPConcentration(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                             ConcentrationX          = PMPConcentration,             &
                                             PropertyXIDNumber       = PropID,                       &
                                             STAT                    = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR04'
                    
                    !check if flow computation by layers
                    if (GW_Link /= Layer_) then 
                        call GetGWLayer   (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR05'
                        
                        allocate(PMPConcentration2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                        PMPConcentration2D = FillValueReal                    
                        
                        !2D conc for drainage network
                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                                k = GWLayer(i,j)
                                PMPConcentration2D (i,j) = PMPConcentration(i,j,k)
                            endif
                        enddo
                        enddo          

                        call UnGetPorousMedia   (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR06'
                        
                        !Send conc from PMP to drainage network
                        call SetPMPConcDN             (DrainageNetworkID       = Me%ObjDrainageNetwork,        &
                                                       ConcentrationX2D        = PMPConcentration2D,           &
                                                       PropertyXIDNumber       = PropID,                       &
                                                       STAT                    = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR07'
                        
                        deallocate (PMPConcentration2D)
                        
                        call UngetPorousMediaProperties(Me%ObjPorousMediaProperties, PMPConcentration, STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR08'
                    
                    else

                        !Send conc from PMP to drainage network - 3D
                        call SetPMPConcDN             (DrainageNetworkID       = Me%ObjDrainageNetwork,        &
                                                       ConcentrationX3D        = PMPConcentration,             &
                                                       PropertyXIDNumber       = PropID,                       &
                                                       STAT                    = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR09'                    
 
                         call UngetPorousMediaProperties(Me%ObjPorousMediaProperties, PMPConcentration, STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR09a'
                        
                    endif
                endif
                
            enddo                      
        endif

        if (Me%Coupled%RunoffProperties) then
        
            call GetRPnProperties (Me%ObjRunoffProperties, nProperties, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR010'

            do iProp = 1, nProperties

                call GetRPPropertiesIDByIdx(RunoffPropertiesID = Me%ObjRunoffProperties,            &
                                             Idx                     = iProp,                       &
                                             ID                      = PropID,                      &
                                             PropAdvDiff             = PropAdvDiff,                 &
                                             STAT                    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR020' 
                
                !Particulate and not particulate properties (in water phase)with advection diffusion may interact 
                !with drainage network
                if (PropAdvDiff) then
                    
                    !Get the property conc from RP
                    call GetRPConcentration(RunoffPropertiesID       = Me%ObjRunoffProperties,       &
                                             ConcentrationX          = RPConcentration,              &
                                             PropertyXIDNumber       = PropID,                       &
                                             STAT                    = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR030'
                    
                    !And send it to drainage network
                    call SetRPConcDN             (DrainageNetworkID        = Me%ObjDrainagenetwork,        &
                                                   ConcentrationX          = RPConcentration,              &
                                                   PropertyXIDNumber       = PropID,                       &
                                                   STAT                    = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR040'
                    
                    call UngetRunoffProperties(Me%ObjRunoffProperties, RPConcentration, STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR050'
                
                endif
                
            enddo                      
        endif

        !Check if Radiation is needed
        call GetNeedsRadiation (Me%ObjDrainageNetwork, NeedsRadiation, STAT_CALL)

        !Check if Radiation is needed
        call GetNeedsAtmosphere(Me%ObjDrainageNetwork, NeedAtmosphere, STAT_CALL)
                
        if (NeedAtmosphere .or. NeedsRadiation) then

            !SolarRadiation
            call GetAtmosphereProperty  (Me%ObjAtmosphere, SolarRadiation, ID = SolarRadiation_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR060'

            call SetAtmosphereDrainageNet (Me%ObjDrainageNetwork, TopRadiation = SolarRadiation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR070'
        
        endif

        if (NeedAtmosphere) then

            !AirTemperature
            call GetAtmosphereProperty  (Me%ObjAtmosphere, AirTemperature, ID = AirTemperature_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR080'

            !CloudCover
            call GetAtmosphereProperty  (Me%ObjAtmosphere, CloudCover, ID = CloudCover_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR090'

            !RelativeHumidity
            call GetAtmosphereProperty  (Me%ObjAtmosphere, RelativeHumidity, ID = RelativeHumidity_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR0100'

            !WindSpeed
            call GetAtmosphereProperty  (Me%ObjAtmosphere, WindSpeed, ID = WindModulus_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR0110'

            call SetAtmosphereDrainageNet (Me%ObjDrainageNetwork,                           &
                                           AirTemperature      = AirTemperature,            &
                                           CloudCover          = CloudCover,                &
                                           RelativeHumidity    = RelativeHumidity,          &
                                           WindSpeed           = WindSpeed,                 &
                                           STAT                = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR0120'   

        endif
        
        if (GW_Link /= Layer_) then 
            !Runs DrainageNetwork with GW flow given by one value for each river cell
            call ModifyDrainageNetwork  ( Me%ObjDrainageNetwork, OLFlowToChannels = OLFlowToChannels,          &
                                          GWFlowToChannels = GWFlowToChannels, DiffuseFlow = Me%DiffuseFlow,   &
                                          STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR130'     
        else
            !Runs DrainageNetwork with GW flow given by layers for each river cell
            call ModifyDrainageNetwork  ( Me%ObjDrainageNetwork, OLFlowToChannels = OLFlowToChannels,                    &
                                          GWFlowToChannels = GWFlowToChannels,                                           &
                                          GWFlowToChannelsLayer = GWFlowToChannelsLayer, DiffuseFlow = Me%DiffuseFlow,   &
                                          STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR140'     
        endif    
                    
        !Ungets
        if (NeedAtmosphere .or. NeedsRadiation) then
            call UnGetAtmosphere  (Me%ObjAtmosphere, SolarRadiation,    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR150'
        endif

        !Ungets
        if (NeedAtmosphere) then
            
            !AirTemperature
            call UnGetAtmosphere  (Me%ObjAtmosphere, AirTemperature,    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR160'

            !CloudCover
            call UnGetAtmosphere  (Me%ObjAtmosphere, CloudCover,        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR170'

            !RelativeHumidity
            call UnGetAtmosphere  (Me%ObjAtmosphere, RelativeHumidity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR180'

            !WindSpeed
            call UnGetAtmosphere  (Me%ObjAtmosphere, WindSpeed,         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR190'
        endif

        call UnGetRunOff            (Me%ObjRunOff, OLFlowToChannels, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR200' 

        if (Me%Coupled%PorousMedia) then
            
            call UnGetPorousMedia       (Me%ObjPorousMedia, GWFlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR210' 
            
            if (GW_Link == Layer_) then 
                call UnGetPorousMedia       (Me%ObjPorousMedia, GWFlowToChannelsLayer, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR220' 

                call UnGetPorousMedia       (Me%ObjPorousMedia, GWFlowBottomLayer, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR230' 
            
                call UnGetPorousMedia       (Me%ObjPorousMedia, GWFlowTopLayer, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR240' 
            endif
        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "DrainageNetworkProcesses")

    end subroutine DrainageNetworkProcesses
    
    !--------------------------------------------------------------------------

    subroutine PorousMediaProcesses

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i,j 
        integer                                     :: STAT_CALL
        real(8), dimension(:, :), pointer           :: Infiltration
        real(8), dimension(:, :), pointer           :: EfectiveEVTP
        real,    dimension(:,: ), pointer           :: PotentialEvaporation
        type (T_BasinProperty), pointer             :: RefEvapotrans
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "PorousMediaProcesses")

        !Calculates Column which may infiltrate
        !Throughfall + Watercolumn over MIN_WC
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then
            
                !Potencial Infiltration Column = ThroughFall
                Me%PotentialInfCol (i, j) = Me%ThroughFall(i, j) + Me%WaterColumnCoef * Me%WaterColumn(i, j)
            end if
        end do
        end do

        if (Me%Coupled%Vegetation) then

            !if the user choosed to separate transpiration and evaporation,
            !evaporation is explicitly defined
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                
                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        PotentialEvaporation = Me%PotentialEvaporation,            &
                                        ActualTranspiration  = Me%ExtVar%ActualTranspiration,      &
                                        STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR01' 
            
            else !No evaporation

                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        ActualTranspiration  = Me%ExtVar%ActualTranspiration,      &
                                        STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR02' 


            endif
        
        else 

            if (Me%Coupled%Evapotranspiration) then
                !If no vegetation, there is no transpiration. 
                !And potential evapotranspiration is all in form of potential evaporation
                call SearchProperty(RefEvapotrans, RefEvapotrans_        , .true., STAT = STAT_CALL)        

                !m/s
                PotentialEvaporation => RefEvapotrans%Field
                
                
                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        PotentialEvaporation = PotentialEvaporation,               &
                                        STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR21' 
            
            else
                
                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR022' 
            
            endif                
        
        endif

        call GetInfiltration   (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR023'

        call GetEfectiveEVTP   (Me%ObjPorousMedia, EfectiveEVTP, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR024'

        !Updates Water column
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                !Flow Production = Potencial Infiltration - Efective infiltration - WC
                Me%FlowProduction   (i, j) = Me%PotentialInfCol (i, j)      - Infiltration      (i, j) - &
                                             Me%WaterColumnCoef * Me%WaterColumn(i, j)
            
                !Increase Waterlevel due to flow production
                Me%WaterLevel       (i, j) = Me%WaterLevel (i, j)           + Me%FlowProduction (i, j)
                
                !Output - mm /hour
                Me%InfiltrationRate (i, j) = Infiltration  (i, j) / Me%CurrentDT * 1000.0 * 3600.0

                !Output - mm /hour
                Me%EVTPRate         (i, j) = EfectiveEVTP  (i, j) / Me%CurrentDT * 1000.0 * 3600.0

                !m
                Me%AccInfiltration  (i, j) = Me%AccInfiltration  (i, j) + Infiltration (i,j)

                !m
                Me%AccEVTP          (i, j) = Me%AccEVTP  (i, j)         + EfectiveEVTP (i,j)
                
                !m
                Me%AccFlowProduction(i, j) = Me%AccFlowProduction(i, j) + Me%FlowProduction (i, j)
                
            endif
        enddo
        enddo

        !Overall Volume
        !EvapFromSoil -> System Loss
        !Infiltration -> FluxAmong Modules
        Me%MB%EvapFromSoil = 0.0
        Me%MB%Infiltration = 0.0
        if (Me%VerifyGlobalMass) then
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    Me%MB%EvapFromSoil = Me%MB%EvapFromSoil + EfectiveEVTP(i,j) * Me%ExtVar%GridCellArea(i, j)
                    Me%MB%Infiltration = Me%MB%Infiltration + Infiltration(i,j) * Me%ExtVar%GridCellArea(i, j)
                endif
            enddo
            enddo
        endif


        call UnGetPorousMedia     (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR03'
                
        call UnGetPorousMedia     (Me%ObjPorousMedia, EfectiveEVTP, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR03a'

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "PorousMediaProcesses")


    end subroutine PorousMediaProcesses

    !--------------------------------------------------------------------------

    subroutine PorousMediaPropertiesProcesses

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical, dimension(:,:), pointer            :: SoilFluxesActive
        real, dimension(:,:), pointer               :: GrazingBiomass
        real, dimension(:,:), pointer               :: GrazingNitrogen
        real, dimension(:,:), pointer               :: GrazingPhosphorus
        real, dimension(:,:), pointer               :: ManagementAerialBiomass
        real, dimension(:,:), pointer               :: ManagementNitrogen
        real, dimension(:,:), pointer               :: ManagementPhosphorus
        real, dimension(:,:), pointer               :: ManagementRootBiomass
        real, dimension(:,:), pointer               :: DormancyBiomass
        real, dimension(:,:), pointer               :: DormancyNitrogen
        real, dimension(:,:), pointer               :: DormancyPhosphorus
        real, dimension(:,:), pointer               :: FertilNitrateSurface
        real, dimension(:,:), pointer               :: FertilNitrateSubSurface
        real, dimension(:,:), pointer               :: FertilAmmoniaSurface
        real, dimension(:,:), pointer               :: FertilAmmoniaSubSurface
        real, dimension(:,:), pointer               :: FertilOrganicNSurface
        real, dimension(:,:), pointer               :: FertilOrganicNSubSurface
        real, dimension(:,:), pointer               :: FertilOrganicPSurface
        real, dimension(:,:), pointer               :: FertilOrganicPSubSurface
        real, dimension(:,:), pointer               :: FertilMineralPSurface
        real, dimension(:,:), pointer               :: FertilMineralPSubSurface
        real, dimension(:,:), pointer               :: RootDepth
        real, dimension(:,:), pointer               :: NitrogenFraction
        real, dimension(:,:), pointer               :: PhosphorusFraction
        real, dimension(:,:,:),pointer              :: NitrogenUptake
        real, dimension(:,:,:),pointer              :: PhosphorusUptake
        logical                                    :: Grazing
        logical                                    :: Management
        logical                                    :: Dormancy
        logical                                    :: Fertilization
        logical                                    :: NutrientFluxesWithSoil
        logical                                    :: CoupledSedimentQuality
        logical                                    :: ModelNitrogen
        logical                                    :: ModelPhosphorus
        logical                                    :: GrowthModel
        real, dimension(:,:), pointer               :: WindVelocity
        real                                        :: VegetationDT
        real, dimension (:), pointer               :: DNConcentration 
        real, dimension (:,:), pointer             :: RPConcentration
        real, dimension (:,:), pointer             :: AtmConcentration  , InfColConcentration
        integer, dimension(:, :), pointer          :: ChannelsID
        integer                                    :: nProperties, iProp, PropID
        logical                                    :: PropAdvDiff, PropParticulate !, PropRain, PropIrri

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "PorousMediaPropertiesProcesses")

        !Compute mass balance due to rain. 
        call GetRPnProperties (Me%ObjRunoffProperties, nProperties, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.2'

        do iProp = 1, nProperties

            call GetRPPropertiesIDByIdx(RunoffPropertiesID = Me%ObjRunoffProperties,            &
                                         Idx                     = iProp,                       &
                                         ID                      = PropID,                      &
                                         PropAdvDiff             = PropAdvDiff,                 &
                                         Particulate             = PropParticulate,             &
                                         STAT                    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.3' 
            
            if (PropAdvDiff .and. (.not. PropParticulate)) then

                allocate (InfColConcentration(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                InfColConcentration = FillValueReal

            
                call GetRPConcentration(RunoffPropertiesID       = Me%ObjRunoffProperties,       &
                                         ConcentrationX          = RPConcentration,              &
                                         PropertyXIDNumber       = PropID,                       &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.4'
                
                if ((Me%Coupled%Atmosphere) .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then

                    call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
                                               Scalar             = AtmConcentration,       &
                                               ID                 = PropID,                 &
                                               STAT               = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR10.5'
 
                    call ComputePropertyInfilColumn (AtmConcentration, RPConcentration, InfColConcentration)
                    
                    call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.7'   
                    
                else
                
                    InfColConcentration => RPConcentration
                
                endif
                
                
                call SetInfColConcPMP         (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                               ConcentrationX          = InfColConcentration,          &
                                               PropertyXIDNumber       = PropID,                       &
                                               STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.5'

                deallocate (InfColConcentration)        
                
               
                call UngetRunoffProperties (Me%ObjRunoffProperties, RPConcentration, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.6'
            
            endif
            
        enddo  
        
!        if (Me%Coupled%Atmosphere) then
!        
!            call GetAtmospherenProperties (Me%ObjAtmosphere, nProperties, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR00.1'
!
!            do iProp = 1, nProperties
!
!                call GetAtmospherePropertiesIDByIdx(AtmosphereID            = Me%ObjAtmosphere,            &
!                                                    Idx                     = iProp,                       &
!                                                    ID                      = PropID,                      &
!                                                    PropRain                = PropRain,                    &
!                                                    PropIrri                = PropIrri,                    &
!                                                    STAT                    = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.2' 
!                
!                if (PropRain .or. PropIrri) then
!                    
!                    allocate (InfColConcentration(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
!                    InfColConcentration = FillValueReal
!                    
!                    call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
!                                               Scalar             = AtmConcentration,       &
!                                               ID                 = PropID,                 &
!                                               STAT               = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.3'
!                    
!!                    if (Me%Coupled%RunoffProperties) then
!                        call GetRPConcentration   (RunoffPropertiesID = Me%ObjRunoffProperties, &
!                                                   ConcentrationX     = RPConcentration,        &
!                                                   PropertyXIDNumber  = PropID,                 &
!                                                   STAT               = STAT_CALL)        
!                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.4'
!!                    else
!!                        allocate (RPConcentration(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
!!                        RPConcentration = FillValueReal
!!                    
!!                        call ComputePropertyWaterColumn(AtmConcentration, RPConcentration)
!!                    endif
!                    
!                    call ComputePropertyInfilColumn (AtmConcentration, RPConcentration, InfColConcentration)
!                    
!                    call SetInfColConcPMP         (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
!                                                   ConcentrationX          = InfColConcentration,          &
!                                                   PropertyXIDNumber       = PropID,                       &
!                                                   STAT                    = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.5'
!
!!                    if (Me%Coupled%RunoffProperties) then
!                        call UngetRunoffProperties (Me%ObjRunoffProperties, RPConcentration, STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.6'
!!                    else
!!                        deallocate(RPConcentration)
!!                    endif
!
!                    
!                    call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.7'   
!                    
!                    deallocate (InfColConcentration)               
!                
!                endif
!                
!            enddo                  
!        
!        endif
        
        !Vegetation organic matter and nutrient fluxes to/from soil
        if (Me%Coupled%Vegetation) then

            call GetVegetationOptions (Me%ObjVegetation,                                          &
                                      NutrientFluxesWithSoil = NutrientFluxesWithSoil,           &
                                      Grazing                = Grazing,                          &
                                      Management             = Management,                       &
                                      Dormancy               = Dormancy,                         &
                                      Fertilization          = Fertilization,                    &
                                      ModelNitrogen          = ModelNitrogen,                    &
                                      ModelPhosphorus        = ModelPhosphorus,                  &
                                      GrowthModel            = GrowthModel,                      &
                                      STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR01'

            
            if (NutrientFluxesWithSoil) then

                call GetVegetationDT  (Me%ObjVegetation, VegetationDT, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR02'

                call GetRootDepth  (Me%ObjVegetation, RootDepth, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR03'
                
                call GetVegetationSoilFluxes  (VegetationID             = Me%ObjVegetation,                &
                                               NitrogenUptake           = NitrogenUptake,                  &
                                               PhosphorusUptake         = PhosphorusUptake,                &
                                               SoilFluxesActive         = SoilFluxesActive,                &  
                                               STAT                     = STAT_CALL)               
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR010'

                
                !Set to porous media properties the information needed from vegetation                
                call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties,     &
                                               NutrientFluxesWithSoil   = NutrientFluxesWithSoil,          &
                                               NitrogenUptake           = NitrogenUptake,                  &
                                               PhosphorusUptake         = PhosphorusUptake,                &
                                               SoilFluxesActive         = SoilFluxesActive,                &
                                               RootDepth                = RootDepth,                       &
                                               ModelNitrogen            = ModelNitrogen,                   &
                                               ModelPhosphorus          = ModelPhosphorus,                 &
                                               GrowthModel              = GrowthModel,                     &
                                               CoupledVegetation        = .true.,                          &
                                               VegetationDT             = VegetationDT,                    &
                                               STAT                     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR020'

                call UnGetVegetation  (Me%ObjVegetation, RootDepth, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR021'
                                


                if (GrowthModel) then

                   
                    call GetNutrientFraction (VegetationID             = Me%ObjVegetation,                &
                                              NitrogenFraction         = NitrogenFraction,                &
                                              PhosphorusFraction       = PhosphorusFraction,              &
                                              STAT                     = STAT_CALL)               
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR022'

                    call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                   NitrogenFraction         = NitrogenFraction,            &
                                                   PhosphorusFraction       = PhosphorusFraction,          &
                                                   STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR023'
 
                    call UnGetVegetation  (Me%ObjVegetation, NitrogenFraction, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR024'

                    call UnGetVegetation  (Me%ObjVegetation, PhosphorusFraction, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR025'

            
                
                    if (Grazing) then
                        call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                        GrazingBiomass           = GrazingBiomass,              &
                                                        GrazingNitrogen          = GrazingNitrogen,             &   
                                                        GrazingPhosphorus        = GrazingPhosphorus,           &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR030'

                        call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                       Grazing                  = Grazing,                     &
                                                       GrazingBiomass           = GrazingBiomass,              &
                                                       GrazingNitrogen          = GrazingNitrogen,             &
                                                       GrazingPhosphorus        = GrazingPhosphorus,           &
                                                       STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR040'
                    endif
                
                    if (Management) then
                        call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                        ManagementAerialBiomass  = ManagementAerialBiomass,     &
                                                        ManagementNitrogen       = ManagementNitrogen,          &
                                                        ManagementPhosphorus     = ManagementPhosphorus,        &
                                                        ManagementRootBiomass    = ManagementRootBiomass,       &
                                                       STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR050'

                        call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                       Management               = Management,                  &
                                                       ManagementAerialBiomass  = ManagementAerialBiomass,     &
                                                       ManagementNitrogen       = ManagementNitrogen,          &
                                                       ManagementPhosphorus     = ManagementPhosphorus,        &
                                                       ManagementRootBiomass    = ManagementRootBiomass,       &     
                                                      STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR060'
                    endif

                    if (Dormancy) then
                        call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                        DormancyBiomass          = DormancyBiomass,             &
                                                        DormancyNitrogen         = DormancyNitrogen,            &
                                                        DormancyPhosphorus       = DormancyPhosphorus,          &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR070'

                        call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                       Dormancy                 = Dormancy,                    &
                                                       DormancyBiomass          = DormancyBiomass,             &
                                                       DormancyNitrogen         = DormancyNitrogen,            &
                                                       DormancyPhosphorus       = DormancyPhosphorus,          &
                                                       STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR080'
                    endif

                    if (Fertilization) then
                        call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,           &
                                                        FertilNitrateSurface     = FertilNitrateSurface,       &
                                                        FertilNitrateSubSurface  = FertilNitrateSubSurface,    &
                                                        FertilAmmoniaSurface     = FertilAmmoniaSurface,       &
                                                        FertilAmmoniaSubSurface  = FertilAmmoniaSubSurface,    &
                                                        FertilOrganicNSurface    = FertilOrganicNSurface,      &
                                                        FertilOrganicNSubSurface = FertilOrganicNSubSurface,   &
                                                        FertilOrganicPSurface    = FertilOrganicPSurface,      &
                                                        FertilOrganicPSubSurface = FertilOrganicPSubSurface,   &
                                                        FertilMineralPSurface    = FertilMineralPSurface,      &
                                                        FertilMineralPSubSurface = FertilMineralPSubSurface,   &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR090'

                        call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                        Fertilization            = Fertilization,              &
                                                        FertilNitrateSurface     = FertilNitrateSurface,       &
                                                        FertilNitrateSubSurface  = FertilNitrateSubSurface,    &
                                                        FertilAmmoniaSurface     = FertilAmmoniaSurface,       &
                                                        FertilAmmoniaSubSurface  = FertilAmmoniaSubSurface,    &
                                                        FertilOrganicNSurface    = FertilOrganicNSurface,      &
                                                        FertilOrganicNSubSurface = FertilOrganicNSubSurface,   &
                                                        FertilOrganicPSurface    = FertilOrganicPSurface,      &
                                                        FertilOrganicPSubSurface = FertilOrganicPSubSurface,   &
                                                        FertilMineralPSurface    = FertilMineralPSurface,      &
                                                        FertilMineralPSubSurface = FertilMineralPSubSurface,   &       
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0100'
                    endif
                endif
            endif
        endif                                                    
        
        !sediment quality needs wind velocity
        call GetPMPCoupled (PorousMediaPropertiesID    = Me%ObjPorousMediaProperties,                    &
                            SoilQuality                = CoupledSedimentQuality,                         &       
                            STAT                       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0101'
        
        if (CoupledSedimentQuality) then

            !Wind Velocity
            call GetAtmosphereProperty  (Me%ObjAtmosphere, WindVelocity, ID = WindModulus_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0102'


            call SetWindVelocity (PorousMediaPropertiesID    = Me%ObjPorousMediaProperties,              &
                                  WindModulus                = WindVelocity,                             &       
                                  STAT                       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0103'

            call UngetAtmosphere (Me%ObjAtmosphere, WindVelocity, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.1'   

        endif


!        if (Me%Coupled%RunoffProperties) then
!        
!            call GetRPnProperties (Me%ObjRunoffProperties, nProperties, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.2'
!
!            do iProp = 1, nProperties
!
!                call GetRPPropertiesIDByIdx(RunoffPropertiesID = Me%ObjRunoffProperties,            &
!                                             Idx                     = iProp,                       &
!                                             ID                      = PropID,                      &
!                                             PropAdvDiff             = PropAdvDiff,                 &
!                                             STAT                    = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.3' 
!                
!                if (PropAdvDiff) then
!                
!                    call GetRPConcentration(RunoffPropertiesID       = Me%ObjRunoffProperties,       &
!                                             ConcentrationX          = RPConcentration,              &
!                                             PropertyXIDNumber       = PropID,                       &
!                                             STAT                    = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.4'
!                    
!                    call SetRPConcPMP             (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
!                                                   ConcentrationX          = RPConcentration,              &
!                                                   PropertyXIDNumber       = PropID,                       &
!                                                   STAT                    = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.6'
!
!                    call UngetRunoffProperties (Me%ObjRunoffProperties, RPConcentration, STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.7'
!                
!                endif
!                
!            enddo                      
!        endif
        
        if (Me%Coupled%DrainageNetwork) then
            
            call GetChannelsID   (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0104'              
            
            call GetDNnProperties (Me%ObjDrainageNetwork, nProperties, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0105'

            do iProp = 1, nProperties

                call GetDNPropertiesIDByIdx(DrainageNetworkID = Me%ObjDrainageNetwork,                 &
                                            Idx               = iProp,                                 &
                                            ID                = PropID,                                &
                                            PropAdvDiff       = PropAdvDiff,                           &
                                            STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR106' 
                
                !Only advection diffusion properties and not particulate may interact with soil
                if (PropAdvDiff .and. (.not. Check_Particulate_Property(PropID))) then
                    
                    !Get the property conc from Drainage Network
                    call GetDNConcentration   (DrainageNetworkID = Me%ObjDrainageNetwork,                  &
                                               ConcentrationX    = DNConcentration,                        &
                                               PropertyXIDNumber = PropID,                                 &
                                               STAT              = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0107'
                
                    !And send it to PorousMediaProperties
                    call SetDNConcPMP (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,               &
                                       PropertyID                 = PropID,                                 &
                                       DNConcentration            = DNConcentration,                        &
                                       ChannelsID                 = ChannelsID,                             &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0108'   

                    call UnGetDrainageNetwork (Me%ObjDrainageNetwork, DNConcentration, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR109'
                
                endif   

                
            enddo   
            
            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR109a'      
             
        endif
        
        call ModifyPorousMediaProperties(ObjPorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                         WaterColumn                = Me%WaterColumn,               &
!        								 ThroughFall                = Me%ThroughFall,               &
!        								 WCEvaporated               = Me%WaterColumnEvaporated,     &
                                         STAT                       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0110' 
        
        
        !Unlock
        if (Me%Coupled%Vegetation) then

            if (NutrientFluxesWithSoil) then

                call UngetVegetationSoilFluxes  (VegetationID             = Me%ObjVegetation,                &
                                               NitrogenUptake           = NitrogenUptake,                  &
                                               PhosphorusUptake         = PhosphorusUptake,                &
                                               SoilFluxesActive         = SoilFluxesActive,                &  
                                               STAT                     = STAT_CALL)               
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0120'
               
                if (Grazing) then
                    call UngetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                    GrazingBiomass           = GrazingBiomass,              &
                                                    GrazingNitrogen          = GrazingNitrogen,             &   
                                                    GrazingPhosphorus        = GrazingPhosphorus,           &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0130'
                endif
                
                if (Management) then
                    call UngetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                    ManagementAerialBiomass  = ManagementAerialBiomass,     &
                                                    ManagementNitrogen       = ManagementNitrogen,          &
                                                    ManagementPhosphorus     = ManagementPhosphorus,        &
                                                    ManagementRootBiomass    = ManagementRootBiomass,       &
                                                   STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0140'
                endif

                if (Dormancy) then
                    call UngetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                    DormancyBiomass          = DormancyBiomass,             &
                                                    DormancyNitrogen         = DormancyNitrogen,            &
                                                    DormancyPhosphorus       = DormancyPhosphorus,          &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0150'
                endif

                if (Fertilization) then
                    call UngetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,           &
                                                    FertilNitrateSurface     = FertilNitrateSurface,       &
                                                    FertilNitrateSubSurface  = FertilNitrateSubSurface,    &
                                                    FertilAmmoniaSurface     = FertilAmmoniaSurface,       &
                                                    FertilAmmoniaSubSurface  = FertilAmmoniaSubSurface,    &
                                                    FertilOrganicNSurface    = FertilOrganicNSurface,      &
                                                    FertilOrganicNSubSurface = FertilOrganicNSubSurface,   &
                                                    FertilOrganicPSurface    = FertilOrganicPSurface,      &
                                                    FertilOrganicPSubSurface = FertilOrganicPSubSurface,   &
                                                    FertilMineralPSurface    = FertilMineralPSurface,      &
                                                    FertilMineralPSubSurface = FertilMineralPSubSurface,   &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0160'
                endif
            endif
        endif

        if (CoupledSedimentQuality) then
            !Wind Velocity
            call UnGetAtmosphere  (Me%ObjAtmosphere, WindVelocity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0162'
        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "PorousMediaPropertiesProcesses")


    end subroutine PorousMediaPropertiesProcesses

    !--------------------------------------------------------------------------

    subroutine ComputePropertyInfilColumn (AtmConcentration, RPConcentration, InfColConcentration)
    
        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer            :: AtmConcentration    !IN
        real, dimension(:,:), pointer            :: RPConcentration     !IN
        real, dimension(:,:), pointer            :: InfColConcentration !OUT
        
        !Local-----------------------------------------------------------------
        real                                :: RainVolume, MassOnRain
        real                                :: WaterColumnVolume, MassOnWaterColumn
        real                                :: InfColumnVolume, MassOnInfColumn
        integer                             :: i,j !,CHUNK

        !----------------------------------------------------------------------       
        
        !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             

                !mass of the property on rain
                !m3        =  m * m2
                RainVolume = Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j)
                !g         =     g/m3       *           m                     m2 
                MassOnRain = AtmConcentration(i,j) * (Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j))
                
                !mass of the property in water column
                !m3 = m * m2
                WaterColumnVolume = (Me%PotentialInfCol(i,j) - Me%ThroughFall(i, j)) * Me%ExtVar%GridCellArea(i, j)
                !g = g/m3 * m3
                MassOnWaterColumn = RPConcentration(i,j) * WaterColumnVolume
                
                !mass of the property in infiltration column
                MassOnInfColumn   = MassOnRain + MassOnWaterColumn
                InfColumnVolume   = RainVolume + WaterColumnVolume
                
                !Compute infiltration column concetration 
                if (Me%PotentialInfCol(i,j) .gt. 0.0) then
                    !g/m3 = g  /  m3
                    InfColConcentration(i, j) = MassOnInfColumn / InfColumnVolume
                else
                    InfColConcentration(i,j) = 0.0
                endif
                
               
            endif
        enddo
        enddo
        !!$OMP END DO 
        
        !!$OMP END PARALLEL
    
    
    end subroutine ComputePropertyInfilColumn

    !--------------------------------------------------------------------------


    subroutine RunoffPropertiesProcesses

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                    :: PropAdvDiff
!        real, dimension(:,:), pointer               :: WindVelocity
!        real, dimension(:, :   ), pointer           :: PMPConcentration2D
!        real, dimension(:, :, :), pointer           :: PMPConcentration
        real, dimension (:), pointer               :: DNConcentration
        integer, dimension(:, :), pointer          :: ChannelsID
        integer                                    :: nProperties, iProp, PropID
        logical                                    :: SplashErosion, ModelCanopyHeight
        real, dimension(:,:), pointer              :: CanopyHeight
 !       type (T_Size3D)                            :: WorkSize3D
 !       real, dimension(:, :), pointer             :: InfiltrationFlux
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "RunoffPropertiesProcesses")
        
        !sediment quality needs wind velocity
!        call GetRPCoupled (RunoffPropertiesID          = Me%ObjRunoffProperties,                    &
!                            SoilQuality                = CoupledSedimentQuality,                         &       
!                            STAT                       = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0101'
!        
!        if (CoupledSedimentQuality) then
!
!            !Wind Velocity
!            call GetAtmosphereProperty  (Me%ObjAtmosphere, WindVelocity, ID = WindModulus_, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0102'
!
!
!            call SetWindVelocity (RunoffPropertiesID        = Me%ObjRunoffProperties,              &
!                                  WindModulus                = WindVelocity,                             &       
!                                  STAT                       = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0103'
!
!        endif
        
!        if (Me%Coupled%PorousMedia) then
!        
!            call GetInfiltration   (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0103.1'              
!            
!            allocate(InfiltrationFlux(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB        ))
!            InfiltrationFlux = null_real
!            
!            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!            do i = Me%WorkSize%ILB, Me%WorkSize%IUB 
!            
!                if (Me%ExtVar%BasinPoints (i, j) == 1) then
!                    !m3/s = m * m2 / s
!                    InfiltrationFlux(i,j) = - Infiltration(i, j) * Me%ExtVar%GridCellArea(i, j) / Me%CurrentDT
!                    
!                endif
!            
!            enddo
!            enddo           
!            
!            call SetInfiltrationRP(RunoffPropertiesID         = Me%ObjRunoffProperties,                  &
!                                   InfiltrationFlux           = InfiltrationFlux,                        &
!                                   STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0103.2'      
!            
!            deallocate (InfiltrationFlux)
!            
!            call UnGetPorousMedia   (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0103.3'    
!        
!        endif

!        if (Me%Coupled%PorousMediaProperties) then
!        
!            call GetPMPnProperties (Me%ObjPorousMediaProperties, nProperties, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR02'
!
!            do iProp = 1, nProperties
!
!                call GetPMPPropertiesIDByIdx(PorousMediaPropertiesID = Me%ObjPorousMediaProperties, &
!                                             Idx                     = iProp,                       &
!                                             ID                      = PropID,                      &
!                                             PropAdvDiff             = PropAdvDiff,                 &
!                                             STAT                    = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR03' 
!                
!                if (PropAdvDiff) then
!                
!                    call GetPMPConcentration(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
!                                             ConcentrationX          = PMPConcentration,             &
!                                             PropertyXIDNumber       = PropID,                       &
!                                             STAT                    = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR04'
!
!                    !Geometry Size
!                    call GetGeometrySize    (Me%ObjGeometry,             &    
!                                             WorkSize =    WorkSize3D,   &
!                                             STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR05'
!                    
!                    allocate(PMPConcentration2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
!                    PMPConcentration2D = FillValueReal
!                    
!                    !2D conc for runoff properties
!                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!                        if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
!                            k = WorkSize3D%KUB
!                            PMPConcentration2D (i,j) = PMPConcentration(i,j,k)
!                        endif
!                    enddo
!                    enddo          
!                    
!                    call SetPMPConcRP(RunoffPropertiesID      = Me%ObjRunoffProperties,       &
!                                      PMPConcentration        = PMPConcentration2D,           &
!                                      PropertyXIDNumber       = PropID,                       &
!                                      STAT                    = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR06'
!                    
!                    deallocate(PMPConcentration2D)
!                    
!                    call UngetPorousMediaProperties (Me%ObjPorousMediaProperties, PMPConcentration, STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR07'
!                
!                endif
!                
!            enddo                      
!        endif        
        
        if (Me%Coupled%DrainageNetwork) then
            
            call GetChannelsID   (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0104'              
            
            call GetDNnProperties (Me%ObjDrainageNetwork, nProperties, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0105'

            do iProp = 1, nProperties
            
                call GetDNPropertiesIDByIdx(DrainageNetworkID = Me%ObjDrainageNetwork,                 &
                                            Idx               = iProp,                                 &
                                            ID                = PropID,                                &
                                            PropAdvDiff       = PropAdvDiff,                           &
                                            STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR106' 
                
                !Particulate and not particulate properties (in water phase)with advection diffusion may interact with runoff
                if (PropAdvDiff) then                        
                    call GetDNConcentration   (DrainageNetworkID = Me%ObjDrainageNetwork,                  &
                                               ConcentrationX    = DNConcentration,                        &
                                               PropertyXIDNumber = PropID,                                 &
                                               STAT              = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0107'

                
                    call SetDNConcRP (RunoffPropertiesID          = Me%ObjRunoffProperties,                 &
                                       PropertyID                 = PropID,                                 &
                                       DNConcentration            = DNConcentration,                        &
                                       ChannelsID                 = ChannelsID,                             &
                                       STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0108'      
                endif

                call UnGetDrainageNetwork (Me%ObjDrainageNetwork, DNConcentration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - RunoffPropertiesProcesses - ERR109'
                
            enddo   
            
            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - RunoffPropertiesProcesses - ERR109a'      
             
        endif
        
        !See if splash erosion is computed
        call GetRPOptions (RunoffPropertiesID = Me%ObjRunoffProperties, SplashErosion = SplashErosion, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0110'
        
        !if true need troughfall from basin and canopy height from vegetation
        if (SplashErosion) then
            
            allocate(CanopyHeight(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
            CanopyHeight = 0.0
            
            if (Me%Coupled%Vegetation) then
                call GetVegetationOptions (VegetationID = Me%ObjVegetation, ModelCanopyHeight = ModelCanopyHeight, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0120'
                
                if (ModelCanopyHeight) then
                    call GetCanopyHeight (VegetationID = Me%ObjVegetation, Scalar = CanopyHeight, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0130'
                endif
            endif
        
            call SetBasinToRPSplash   (RunoffPropertiesID      = Me%ObjRunoffProperties,        &
                                       ThroughFall             = Me%ThroughFall,                &
                                       CanopyDrainage          = Me%CanopyDrainage,             &
                                       CanopyHeight            = CanopyHeight,                  &
                                       STAT                    = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR140'
            
            deallocate (CanopyHeight)
        
        endif
        
        call ModifyRunoffProperties(ObjRunoffPropertiesID           = Me%ObjRunoffProperties,  &
 !       								 ThroughFall                = Me%ThroughFall,               &
 !       								 WCEvaporated               = Me%WaterColumnEvaporated,     &
                                         STAT                       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0120' 
           
        
!        if (CoupledSedimentQuality) then
!            !Wind Velocity
!            call UnGetAtmosphere  (Me%ObjAtmosphere, WindVelocity, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0132'
!        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "RunoffPropertiesProcesses")


    end subroutine RunoffPropertiesProcesses

    !--------------------------------------------------------------------------

    subroutine ActualizeWaterColumn (WarningString)

        !Arguments-------------------------------------------------------------
        character (Len = *), intent(in)             :: WarningString

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, CHUNK

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "ActualizeWaterColumn")

        !$OMP PARALLEL PRIVATE(I,J, WarningString)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                Me%WaterColumnOld(i,j) = Me%WaterColumn(i,j)
                
                Me%WaterColumn(i, j) = Me%WaterLevel(i, j) - Me%ExtVar%Topography(i, j)                

                if (Me%WaterColumn(i, j) < 0.0) then
                    
                    !Rounding Error
                    if (Me%WaterColumn(i, j) < -1.e-10) then 
                        write (*,*) WarningString
                        write (*,*) 'Negative Water Column corrected', i, j, Me%WaterColumn(i, j)
                    endif
                    Me%WaterColumn(i, j) = 0.0
                    Me%WaterLevel (i, j) = Me%ExtVar%Topography(i, j)
                
                end if
            endif

        enddo
        enddo
        !$OMP END DO


        !$OMP END PARALLEL


        if (MonitorPerformance) call StopWatch ("ModuleBasin", "ActualizeWaterColumn")

    end subroutine ActualizeWaterColumn

    !--------------------------------------------------------------------------

    subroutine ActualizeWaterColumnConc (WarningString)

        !Arguments-------------------------------------------------------------
        character (Len = *), intent(in)             :: WarningString

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        integer                                     :: nProperties, iProp 
        integer                                     :: PropID
        real,    dimension(:, :   ), pointer        :: RPConcentration, NewRPConcentration
        logical                                     :: PropAdvDiff, PropParticulate
        real(8)                                     :: PropertyMassOld, PropertyMassNew
        real(8), dimension(:,:), pointer            :: MassInFlow, MassToBottom
        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleBasin", "ActualizeWaterColumnConcentration")
        
        
        call GetRPnProperties (Me%ObjRunoffProperties, nProperties, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR1'

        do iProp = 1, nProperties

            call GetRPPropertiesIDByIdx(RunoffPropertiesID       = Me%ObjRunoffProperties,      &
                                         Idx                     = iProp,                       &
                                         ID                      = PropID,                      &
                                         PropAdvDiff             = PropAdvDiff,                 &
                                         Particulate             = PropParticulate,             &
                                         STAT                    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR10' 
            
            if (PropAdvDiff) then
                
                allocate(MassInFlow(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                MassInFlow = 0.0
                allocate(MassToBottom(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                MassToBottom = 0.0
                allocate(NewRPConcentration(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                NewRPConcentration = null_real
                
                !Get the most recent conc from RP
                call GetRPConcentration(RunoffPropertiesID       = Me%ObjRunoffProperties,       &
                                         ConcentrationX          = RPConcentration,              &
                                         PropertyXIDNumber       = PropID,                       &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR20'
                
                !Compute mass flow matrix to update concentrations
                call ComputeMassInFlow (WarningString, RPConcentration, PropID, PropParticulate, MassInFlow)

                !Compute the new RP conc based on the new fluxes
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                        !g = g/m3 * (m * m2)
                        PropertyMassOld = RPConcentration(i,j) * (Me%WaterColumnOld(i,j) * Me%ExtVar%GridCellArea(i,j))
                        
                        PropertyMassNew = PropertyMassOld + MassInFlow(i,j)
                        
                        if (Me%WaterColumn(i,j) .gt. AlmostZero) then
                            !g/m3 = g / (m * m2)
                            NewRPConcentration(i,j) = PropertyMassNew / (Me%WaterColumn(i,j) * Me%ExtVar%GridCellArea(i,j))
                        else
                            NewRPConcentration(i,j) = 0.0
                            
                            !deposition driven by complete infiltration of water column (WC totally infiltrated in time step) 
                            if ((PropParticulate) .and. (Me%WaterColumnOld(i,j) .gt. AlmostZero)) then
                                MassToBottom(i,j) = PropertyMassOld
                            endif                            
                        endif                
                        

                       
                    endif

                enddo
                enddo
 
                call UngetRunoffProperties(Me%ObjRunoffProperties, RPConcentration, STAT_Call)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR40'
                
                !Send new conc to RP module
                call SetBasinConcRP (RunoffPropertiesID        = Me%ObjRunoffProperties,       &
                                       BasinConcentration      = NewRPConcentration,           &
                                       PropertyXIDNumber       = PropID,                       &
                                       MassToBottom            = MassToBottom,                 &
                                       WaterColumn             = Me%WaterColumn,               &
                                       STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR30'
                
                deallocate(MassInFlow)
                deallocate(NewRPConcentration)
                deallocate(MassToBottom)
                
            endif
            
        enddo

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "ActualizeWaterColumnConcentration")

    end subroutine ActualizeWaterColumnConc

    !--------------------------------------------------------------------------

    subroutine ComputeMassInFlow(WarningString, RPConcentration, PropID, Particulate, MassInFlow)
    
        !Arguments-------------------------------------------------------------
        character (Len = *), intent(in)             :: WarningString
        real(8), dimension(:,:), pointer            :: MassInFlow
        real, dimension(:, :   ), pointer           :: RPConcentration, InfColConcentration
        real, dimension(:, : ,:), pointer           :: PMPConcentration
        real, dimension(:, :   ), pointer           :: AtmConcentration
        integer                                     :: PropID,i,j, k, STAT_CALL
        type (T_Size3D)                             :: WorkSize3D   
        logical                                     :: Particulate  
        real(8), dimension(:, :), pointer           :: Infiltration   
        !Local-----------------------------------------------------------------
        
        MassInFlow = 0.0
        
        !Evaporation from WC does not take mass
        if(WarningString == 'AtmosphereProcesses') then
        
            !!Precipitation flux is computed in porous media processes (for consistence with water computation)

            !Evaporation flux
!            Do nothing to evaporation. it takes water but not mass from water column
!            MassInFlow = 0.0
        
        elseif(WarningString == 'PorousMediaProcesses') then

            !Flux with PorousMedia - the balance between rain and infiltration

            if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then
                
                call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
                                           Scalar             = AtmConcentration,       &
                                           ID                 = PropID,                 &
                                           STAT               = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR05.6'               
                
            endif            
            
            !Particulate properties do not enter the soil or exit (particulate do not have advection diffusion in soil)
            if (.not. Particulate) then

                allocate (InfColConcentration(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                InfColConcentration = FillValueReal
                
                if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then
                    
                    !Update conc in infiltration column - this sub is called twice, here and before PorousMediaProperties
                    call ComputePropertyInfilColumn (AtmConcentration, RPConcentration, InfColConcentration)
                    
                else
                    InfColConcentration => RPConcentration
                endif            
                
                call GetPMPConcentration(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                         ConcentrationX          = PMPConcentration,             &
                                         PropertyXIDNumber       = PropID,                       &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR04'
                
                !Geometry Size
                call GetGeometrySize    (Me%ObjGeometry,             &    
                                         WorkSize =    WorkSize3D,   &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR05'

                call GetInfiltration   (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR05.5'

            
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                        
                        !Mass flux with soil (infiltration/exfiltration)
                        if(Infiltration(i,j) .gt. 0.0) then ! positive infiltration - removing mass from WC
                            
                            !g = g/m3 * (m * m2)
                            MassInFlow(i,j) = InfColConcentration(i,j) * (Infiltration(i,j) * Me%ExtVar%GridCellArea(i,j))
                        
                        else ! negative infiltration (exfiltration) - adding mass to WC or zero
                            
                            k = WorkSize3D%KUB
                            MassInFlow(i,j) = PMPConcentration (i,j,k) * (Infiltration(i,j) * Me%ExtVar%GridCellArea(i,j))
                        endif
                        
                        !Mass product for water column (similar to balance in PorousMediaprocesses for water) 
                        !after real infiltration computed
                        !Mass production =  MassInRain (MassInInfilColumn - MassInWaterColumn) - RealInfilMass
                        if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then 
                            !g = g/m3 * (m*m2)
                            MassInflow(i,j) = (AtmConcentration(i,j) * Me%ThroughFall(i,j) * Me%ExtVar%GridCellArea(i,j))   &
                                               - MassInFlow(i,j)
                        else
                            !g = g/m3 * (m*m2)
                            MassInflow(i,j) = - MassInFlow(i,j)
                        endif
                        
                    
                    endif
                enddo
                enddo          
                
                deallocate (InfColConcentration)

                call UnGetPorousMedia     (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR05.7'

                call UngetPorousMediaProperties(Me%ObjPorousMediaProperties, PMPConcentration, STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR06'

            else

                if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then            
                    !Particulate properties do not enter soil but for now may come from rain - for testing
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
                            !mass of the property on rain
                            !g         =     g/m3       *           m                     m2 
                            MassInFlow(i,j) = AtmConcentration(i,j) * (Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j))

                        endif
                    enddo
                    enddo
                endif                
            endif

            if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then
                call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR06.5' 
            endif
        
        elseif (WarningString == 'PorousMediaProcesses 2') then
            
            !Precipitation Flux
            if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then
                
                call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
                                           Scalar             = AtmConcentration,       &
                                           ID                 = PropID,                 &
                                           STAT               = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR07'                    
    
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
                        !mass of the property on rain
                        !g         =     g/m3       *           m                     m2 
                        MassInFlow(i,j) = AtmConcentration(i,j) * (Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j))
                    endif
                enddo
                enddo

                call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR08' 
!            else
                !do nothing - if not atmosphere property mass is zero
            endif 

        elseif (WarningString == 'SimpleInfiltration') then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
                    !mass removed by infiltration
                    !g         =     g/m3       *           m/s   *    s  *     m2 
                    MassInFlow(i,j) = - RPConcentration(i,j) * Me%SI%InfRate%Field(i, j)   &
                                      * Me%CurrentDT * Me%ExtVar%GridCellArea(i,j)
                endif
            enddo
            enddo

        endif
           
             
    end subroutine ComputeMassInFlow
    
    !--------------------------------------------------------------------------
    

    subroutine RemoveWaterColumn

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                if (Me%WaterColumn(i, j) > 0.001) then
                    Me%WaterColumnRemoved(i, j) = Me%WaterColumn(i, j) / Me%WCRemovalTime 
                    Me%WaterColumn(i, j)        = Me%WaterColumn(i, j) - (Me%WaterColumnRemoved(i, j) * Me%CurrentDT)
                    if (Me%WaterColumn(i, j) < 0.0) then
                        Me%WaterColumn(i, j)    = 0.0
                        Me%WaterLevel (i, j)    = Me%ExtVar%Topography(i, j)
                    else
                        Me%WaterLevel (i, j)    = Me%ExtVar%Topography(i, j) + Me%WaterColumn(i, j)
                    endif
                else
                    Me%WaterColumnRemoved(i, j) = 0.0
                endif
                !Convert to m3/s
                Me%WaterColumnRemoved(i, j) = Me%WaterColumnRemoved(i, j) * 1.413E+08
            endif

        enddo
        enddo


    end subroutine RemoveWaterColumn

    !--------------------------------------------------------------------------
    
    subroutine CalculateMass(VolumeBasin, VolumeVegetation, VolumePorousMedia, VolumeChannels)

        !Arguments-------------------------------------------------------------
        real(8)                                     :: VolumeBasin, VolumeVegetation
        real(8)                                     :: VolumePorousMedia, VolumeChannels
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL

        !Calculates VolumeBasin / VolumeVegetation
        VolumeBasin       = 0.0
        VolumeVegetation  = 0.0
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                VolumeBasin     = VolumeBasin       + Me%WaterColumn    (i, j) * Me%ExtVar%GridCellArea(i, j)
                if (Me%Coupled%Vegetation) then
                    VolumeVegetation = VolumeVegetation  + Me%CanopyStorage  (i, j) * Me%ExtVar%GridCellArea(i, j) * &
                                       Me%CoveredFraction(i, j)
                else
                    VolumeVegetation = 0.0
                endif
            endif
        enddo
        enddo
        
        if (Me%Coupled%PorousMedia) then
            call GetTotalStoredVolume (Me%ObjPorousMedia, VolumePorousMedia, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR01'
        endif

        if (Me%Coupled%DrainageNetwork) then
            call GetVolumes(Me%ObjDrainageNetwork, TotalStoredVolume = VolumeChannels, STAT = STAT_CALL)
        endif

   
    end subroutine CalculateMass    

    !--------------------------------------------------------------------------

    subroutine TimeSerieOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i,j
        real(8)                                     :: TotalFlowVolume
        real, dimension(6), target                  :: AuxTime
        type (T_BasinProperty), pointer             :: RefEvapotrans
        real, dimension(:,:), pointer               :: PotentialTranspiration, PotentialEvaporation
        real, dimension(:,:), pointer               :: PotentialEvapoTranspiration
        real, dimension(:,:), pointer               :: ActualTranspiration, ActualEvaporation
        real, dimension(:,:), pointer               :: ActualTP, ActualEVAP
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "TimeSerieOutput")

        !Watercolumn
        call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                              &
                             Data2D_8    = Me%WaterColumn,                               &                             
                             STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR02'

        !Waterlevel
        call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                              &
                             Data2D_8    = Me%WaterLevel,                                &                             
                             STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR03'


        !Infiltration Rate
        call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                             Data2D_8 = Me%InfiltrationRate,                    &                             
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR04'

               !Rain Rate
        call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                             Data2D_8 = Me%PrecipRate,                          &                             
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR04a'
        
        !Through Fall Rate
        call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                             Data2D_8 = Me%ThroughRate,                         &                             
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR04b'

        !Efective EVTP Rate
        call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                             Data2D_8 = Me%EVTPRate,                            &                             
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR05'

        !Watercolumn Removed
        call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                     &
                             Data2D_8    = Me%WaterColumnRemoved,               &                             
                             STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR01'

        if (Me%Coupled%Vegetation) then
            
            !Potential Crop Evapotranspiration
            call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                 Data2D = Me%CropEvapotrans,            &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR07'            
            
            !Canopy Capacity 
            call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                     &
                                 Data2D_8    = Me%CanopyStorageCapacity,            &                             
                                 STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR08'

            !Canopy Storage
            call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                     &
                                 Data2D_8    = Me%CanopyStorage,                    &                             
                                 STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR09'


            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then

                !m3/s
                call GetTranspiration(Me%ObjVegetation, ActualTranspiration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR095'
                !m3/s
                call GetEvaporation(Me%ObjPorousMedia, ActualEvaporation, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR096'
                
                allocate(PotentialTranspiration(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(PotentialEvaporation(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(ActualTP(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(ActualEVAP(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            
                !Convert Units to mm/h
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                        !m/s * 1000 mm/m * 3600 s/h = mm/h
                        PotentialEvaporation(i,j)   = Me%PotentialEvaporation(i,j) * 1000. * 3600.
                        PotentialTranspiration(i,j) = Me%PotentialTranspiration(i,j) * 1000. * 3600.
                        
                        !m3/s / m2 * 1000 mm/m * 3600 s/h = mm/h
                        ActualEVAP(i,j) = ActualEvaporation(i,j) / Me%ExtVar%GridCellArea(i,j) * 1000. * 3600.
                        ActualTP(i,j)   = ActualTranspiration(i,j) / Me%ExtVar%GridCellArea(i,j) * 1000. * 3600.
                    endif
                enddo
                enddo

                call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                     Data2D = PotentialEvaporation,                     &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR10'


                call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                     Data2D = PotentialTranspiration,                   &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR11'

                call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                     Data2D = ActualEVAP,                               &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR115'


                call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                     Data2D = ActualTP,                                 &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR116'


                call UnGetVegetation(Me%ObjVegetation, ActualTranspiration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR0117'

                call UnGetPorousMedia(Me%ObjPorousMedia, ActualEvaporation, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR0118'            

                deallocate (PotentialTranspiration)
                deallocate (PotentialEvaporation)
                deallocate (ActualTP)
                deallocate (ActualEVAP)

            endif

        endif
  
        !Reference EVTP Rate
        if (Me%Coupled%Evapotranspiration) then
            call SearchProperty(RefEvapotrans, RefEvapotrans_        , .true., STAT = STAT_CALL)        

            allocate(PotentialEvapoTranspiration(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

            !Convert Units to mm/h
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    !m/s * 1000 mm/m * 3600 s/h = mm/h
                    PotentialEvapoTranspiration(i,j) = RefEvapotrans%Field(i,j) * 1000. * 3600.
                endif
            enddo
            enddo

            call WriteTimeSerie (Me%ObjTimeSerie,                                       &
                                 Data2D = PotentialEvapoTranspiration,                        &                             
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR13'

            deallocate(PotentialEvapoTranspiration)

        endif

        call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2),         &
                                            AuxTime(3), AuxTime(4),         &
                                            AuxTime(5), AuxTime(6))
        !Integrated Daily Values
        
        !Ouput of daily values
        if (Me%DailyFlow%On .and. Me%Coupled%DrainageNetwork) then
        
            call GetVolumes(Me%ObjDrainageNetwork, TotalFlowVolume = TotalFlowVolume, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR14'
        
            !Sum Volume
            Me%DailyFlow%Flow = Me%DailyFlow%Flow + TotalFlowVolume * Me%CurrentDT
        
            if (Me%DailyFlow%CurrentIndex /= AuxTime(3)) then
            
                !Write time series entry
                Me%TimeSeriesBuffer2(1) = Me%DailyFlow%Flow
                call WriteTimeSerieLine (Me%DailyFlow%ObjTimeSerie, Me%TimeSeriesBuffer2, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR15'
                
                !Resets integrated values
                Me%DailyFlow%Flow = 0.0
                
                !Stores New Index
                Me%DailyFlow%CurrentIndex = AuxTime(3)
                                
            endif
            
        endif
            
        !Ouput of montly values
        if (Me%MonthlyFlow%On .and. Me%Coupled%DrainageNetwork) then
        
            call GetVolumes(Me%ObjDrainageNetwork, TotalFlowVolume = TotalFlowVolume, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR16'
        
            !Sum Volume
            Me%MonthlyFlow%Flow = Me%MonthlyFlow%Flow + TotalFlowVolume * Me%CurrentDT
        
            if (Me%MonthlyFlow%CurrentIndex /= AuxTime(2)) then
            
                !Write time series entry
                Me%TimeSeriesBuffer2(1) = Me%MonthlyFlow%Flow
                call WriteTimeSerieLine (Me%MonthlyFlow%ObjTimeSerie, Me%TimeSeriesBuffer2, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR17'
                
                !Resets integrated values
                Me%MonthlyFlow%Flow = 0.0
                
                !Stores New Index
                Me%MonthlyFlow%CurrentIndex = AuxTime(2)
                                
            endif
            
        endif            

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "TimeSerieOutput")

    end subroutine TimeSerieOutput

    !--------------------------------------------------------------------------

    subroutine HDF5Output

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB    
        integer                                         :: STAT_CALL           
        real, dimension(6), target                      :: AuxTime
        real, dimension(:), pointer                     :: TimePointer

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "HDF5Output")

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB   

        if (Me%CurrentTime >= Me%OutPut%OutTime(Me%OutPut%NextOutPut)) then

            !Gets Open

            !Writes current time
            call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2),         &
                                                AuxTime(3), AuxTime(4),         &
                                                AuxTime(5), AuxTime(6))
            TimePointer => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR01'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = Me%OutPut%NextOutPut,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR02'

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR03'


            !Writes the Open Points
            call HDF5WriteData   (Me%ObjHDF5, "//Grid/OpenPoints",              &
                                  "OpenPoints", "-",                            &
                                  Array2D = Me%ExtVar%OpenPoints2D,             &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR04'


            !Writes the Water Column
            call HDF5WriteData   (Me%ObjHDF5, "//Results/water column",         &
                                  "water column", "m",                          &
                                  Array2D      = Me%WaterColumn,                &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR05'


            !Writes the Water Level
            call HDF5WriteData   (Me%ObjHDF5, "//Results/water level",          &
                                  "water level", "m",                           &
                                  Array2D      = Me%WaterLevel,                 &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR06'

            !Writes the Acc. Rain
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccRainFall",          &
                                  "Acc. Rainfall", "m",                         &
                                  Array2D      = Me%AccRainFall,                &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR06'

            !Writes the Acc Infil
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccInfiltration",      &
                                  "Acc. Infiltration", "m",                     &
                                  Array2D      = Me%AccInfiltration,            &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR06'

            !Writes the Acc Flow Production
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccFlowProduction",    &
                                  "Acc. Flow Production", "m",                  &
                                  Array2D      = Me%AccFlowProduction,          &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR06'

            !Writes the Acc EVTP
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccEVTP",              &
                                  "Acc. Evapotrans", "m",                       &
                                  Array2D      = Me%AccEVTP,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR06'

!            if (Me%Coupled%DrainageNetwork) then
!
!                !Writes Flow Into Channels
!
!                call GetFlowToChannels      (Me%ObjRunOff, FlowToChannels, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01'
!
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/FlowToChannels",       &
!                                      "FlowToChannels", "m3/s",                     &
!                                      Array2D      = FlowToChannels,                &
!                                      OutputNumber = Me%OutPut%NextOutPut,          &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR08'
!
!                call UnGetRunOff (Me%ObjRunOff, FlowToChannels, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01'
!
!
!                !Writes Flow Into Channels GW
!                if (Me%Coupled%PorousMedia) then
!                    call GetGWFlowToChannels    (Me%ObjPorousMedia, FlowToChannels, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01a'
!
!                    call HDF5WriteData   (Me%ObjHDF5, "//Results/GWFlowToChannels",     &
!                                          "GWFlowToChannels", "m3/s",                   &
!                                          Array2D      = FlowToChannels,                &
!                                          OutputNumber = Me%OutPut%NextOutPut,          &
!                                          STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR08a'
!
!                    call UnGetPorousMedia       (Me%ObjPorousMedia, FlowToChannels, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR01a'
!                endif
!
!                !Writes ChannelsWaterDepth 
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       WaterDepth_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8b'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/WaterDepth",       &
!                                      "ChanWaterDepth", "m",                             &
!                                      Array2D      = Me%OutPut%OutputChannels,           &
!                                      OutputNumber = Me%OutPut%NextOutPut,               &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8c'
!
!                !Writes ChannelsWaterLevel 
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       WaterLevel_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8b'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/WaterLevel",       &
!                                      "ChanWaterLevel", "m",                             &
!                                      Array2D      = Me%OutPut%OutputChannels,           &
!                                      OutputNumber = Me%OutPut%NextOutPut,               &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8c'
!
!                !Writes ChannelsFlowX 
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       WaterFluxX_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8d'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/Flow/X",           &
!                                      "ChanFlowX", "m3/s",                               &
!                                      Array2D      = Me%OutPut%OutputChannels,           &
!                                      OutputNumber = Me%OutPut%NextOutPut,               &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8e'
!
!                !Writes ChannelsFlowY
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       WaterFluxY_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/Flow/Y",           &
!                                      "ChanFlowY", "m3/s",                               &
!                                      Array2D      = Me%OutPut%OutputChannels,           &
!                                      OutputNumber = Me%OutPut%NextOutPut,               &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8g'
!                
!                !Writes ChannelsFlowModulus
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       FlowModulus_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/Flow/Modulus",     &
!                                      "ChanFlowMod", "m3/s",                             &
!                                      Array2D      = Me%OutPut%OutputChannels,           &
!                                      OutputNumber = Me%OutPut%NextOutPut,               &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8g'
!
!
!                !Writes ChannelsVelocityX 
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       VelocityU_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8d'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/Velocity/X",          &
!                                      "ChanVelX", "m/s",                                    &
!                                      Array2D      = Me%OutPut%OutputChannels,              &
!                                      OutputNumber = Me%OutPut%NextOutPut,                  &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8e'
!
!                !Writes ChannelsVelocityY
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,     &
!                                       VelocityV_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/Velocity/Y",          &
!                                      "ChanVelY", "m/s",                                    &
!                                      Array2D      = Me%OutPut%OutputChannels,              &
!                                      OutputNumber = Me%OutPut%NextOutPut,                  &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8g'
!                
!                !Writes ChannelsVelocity
!                call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                       VelocityModulus_, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!
!                
!                call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/Velocity/Modulus",    &
!                                      "ChanVelMod", "m/s",                                  &
!                                      Array2D      = Me%OutPut%OutputChannels,              &
!                                      OutputNumber = Me%OutPut%NextOutPut,                  &
!                                      STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8g'
!
!                call GetHasProperties (Me%ObjDrainageNetwork, HasProperties, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8h'
!
!                if (HasProperties) then
!
!                    call GetnProperties (Me%ObjDrainageNetwork, nProperties, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'
!
!                    do iProp = 1, nProperties
!
!                        call GetPropertiesIDByIdx(Me%ObjDrainageNetwork, iProp, PropID,             &
!                                                  OutputName, STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'
!                                                
!                        call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,     &
!                                               PropID%IDNumber, STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!
!                        call HDF5WriteData   (Me%ObjHDF5,                                           &
!                                              "//Results/Channels/Properties/"//trim(OutputName),   &
!                                              trim(OutputName), trim(PropID%Units),                 &
!                                              Array2D      = Me%OutPut%OutputChannels,              &
!                                              OutputNumber = Me%OutPut%NextOutPut,                  &
!                                              STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'
!
!                        call GetPropHasBottomFluxes (Me%ObjDrainageNetwork, PropID%IDNumber, HasBottomFluxes, STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR10'
!
!                        if (HasBottomFluxes) then
!
!                            WriteShearStress = .TRUE.
!
!                            call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,     &
!                                                   PropID%IDNumber, Deposited = HasBottomFluxes, STAT = STAT_CALL)
!                            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!            
!                            OutputName = trim(OutputName)//'_deposited'
!                        
!                            call HDF5WriteData   (Me%ObjHDF5,                                           &
!                                                  "//Results/Channels/Properties/"//trim(OutputName),   &
!                                                  trim(OutputName), trim('kg/m2'),                      &
!                                                  Array2D      = Me%OutPut%OutputChannels,              &
!                                                  OutputNumber = Me%OutPut%NextOutPut,                  &
!                                                  STAT = STAT_CALL)
!                            if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'
!
!                        end if                        
!                    enddo
!
!                    if (WriteShearStress) then
!
!                        call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,  &
!                                               ShearStress_, STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!
!                
!                        call HDF5WriteData   (Me%ObjHDF5, "//Results/Channels/ShearStress",         &
!                                              "ChanShearStress", "Pa",                              &
!                                              Array2D      = Me%OutPut%OutputChannels,              &
!                                              OutputNumber = Me%OutPut%NextOutPut,                  &
!                                              STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8g'
!
!                    end if
!
!                    call GetHasToxicity (Me%ObjDrainageNetwork, HasToxicity, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8i'
!
!                    if (HasToxicity) then
!
!                        call FillOutputMatrix (Me%ObjDrainageNetwork, Me%OutPut%OutputChannels,     &
!                                               GenericProperty_, STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR8f'
!                        
!                        OutputName = 'global_toxicity'
!                        
!                        call HDF5WriteData   (Me%ObjHDF5,                                           &
!                                              "//Results/Channels/Properties/"//trim(OutputName),   &
!                                              trim(OutputName), trim(PropID%Units),                 &
!                                              Array2D      = Me%OutPut%OutputChannels,              &
!                                              OutputNumber = Me%OutPut%NextOutPut,                  &
!                                              STAT = STAT_CALL)
!                        if (STAT_CALL /= SUCCESS_) stop 'RunOffOutput - ModuleRunOff - ERR09'
!    
!                    end if
!                end if
!
!            endif


            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR99'

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "HDF5Output")

    end subroutine HDF5Output

    !--------------------------------------------------------------------------

    subroutine GlobalMassBalance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL
        real, dimension(:,:), pointer               :: OL_FlowAtBoundary
        real(8)                                     :: InitialMass, FinalMass
        real(8)                                     :: Sinks, Sources, MassError
        real(8)                                     :: ME_Bacia, ME_Porous, ME_Drainage
        
        !Gets Mass Balance From Drainage Network
        if (Me%Coupled%DrainageNetwork) then
            call GetVolumes(Me%ObjDrainageNetwork,                                      &
                            TotalInputVolume  = Me%MB%DischargesIn,                     &
                            TotalOutputVolume = Me%MB%OutVolumeChannel,                 &
                            TotalOverTopVolume= Me%MB%OverTopOut,                       &
                            STAT = STAT_CALL)
        endif

       !Gets flow at open boundary of Module Runoff
       Me%MB%OutVolumeOverLand = 0.0
       if (Me%Coupled%Runoff) then
            call GetFlowAtBoundary (Me%ObjRunoff, OL_FlowAtBoundary, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR020'

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BoundaryPoints2D(i, j) == BasinPoint) then
                    Me%MB%OutVolumeOverLand = Me%MB%OutVolumeOverLand + OL_FlowAtBoundary(i, j) * Me%CurrentDT
                endif
                
            enddo
            enddo
            call UnGetRunoff(Me%ObjRunoff, OL_FlowAtBoundary, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR021'
        endif


        !Total
        InitialMass = Me%MB%IniVolumeBasin  + Me%MB%IniVolumeVegetation     + Me%MB%IniVolumePorousMedia  & 
                      + Me%MB%IniVolumeChannels
        FinalMass   = Me%MB%FinalVolumeBasin+ Me%MB%FinalVolumeVegetation   + Me%MB%FinalVolumePorousMedia &
                      + Me%MB%FinalVolumeChannels
        
        !Missing Losses to deep aquifer from PM - not yet implemented
        Sinks       = Me%MB%EvapFromVegetation + Me%MB%EvapFromGround + Me%MB%EvapFromSoil + Me%MB%OutVolumeChannel   &
                      + Me%MB%OutVolumeOverLand
        
        !
        Sources     = Me%MB%RainIn + Me%MB%DischargesIn
        MassError   = InitialMass + Sources - Sinks - FinalMass

        !DEBUG only
        !BACIA
        InitialMass = Me%MB%IniVolumeBasin  + Me%MB%IniVolumeVegetation
        FinalMass   = Me%MB%FinalVolumeBasin+ Me%MB%FinalVolumeVegetation
        Sinks       = Me%MB%EvapFromVegetation + Me%MB%EvapFromGround + Me%MB%OutVolumeOverLand + Me%MB%OLFlowToRiver &
                      + Me%MB%Infiltration
        Sources     = Me%MB%RainIn + Me%MB%OverTopOut
        ME_Bacia    = InitialMass + Sources - Sinks - FinalMass
        
        !POROUSMEDIA
        Sinks       = Me%MB%EvapFromSoil + Me%MB%GWFlowToRiver
        Sources     = Me%MB%Infiltration
        ME_Porous   = Me%MB%IniVolumePorousMedia + Sources - Sinks - Me%MB%FinalVolumePorousMedia
        
        !DRAINAGE NET
        InitialMass = Me%MB%IniVolumeChannels
        FinalMass   = Me%MB%FinalVolumeChannels
        Sinks       = Me%MB%OutVolumeChannel - Me%MB%OverTopOut
        Sources     = Me%MB%OLFlowToRiver + Me%MB%GWFlowToRiver + Me%MB%DischargesIn
        ME_Drainage = InitialMass + Sources - Sinks - FinalMass


        !Final Mass Balance
        Me%TimeSeriesBuffer(1)  = Me%MB%IniVolumeBasin
        Me%TimeSeriesBuffer(2)  = Me%MB%IniVolumeVegetation
        Me%TimeSeriesBuffer(3)  = Me%MB%IniVolumePorousMedia
        Me%TimeSeriesBuffer(4)  = Me%MB%IniVolumeChannels
        Me%TimeSeriesBuffer(5)  = Me%MB%EvapFromVegetation
        Me%TimeSeriesBuffer(6)  = Me%MB%EvapFromGround
        Me%TimeSeriesBuffer(7)  = Me%MB%EvapFromSoil
        Me%TimeSeriesBuffer(8)  = Me%MB%RainIn
        Me%TimeSeriesBuffer(9)  = Me%MB%OutVolumeChannel
        Me%TimeSeriesBuffer(10) = Me%MB%OutVolumeOverLand
        Me%TimeSeriesBuffer(11) = Me%MB%FinalVolumeBasin
        Me%TimeSeriesBuffer(12) = Me%MB%FinalVolumeVegetation
        Me%TimeSeriesBuffer(13) = Me%MB%FinalVolumePorousMedia
        Me%TimeSeriesBuffer(14) = Me%MB%FinalVolumeChannels
        Me%TimeSeriesBuffer(15) = MassError
        Me%TimeSeriesBuffer(16) = ME_Bacia
        Me%TimeSeriesBuffer(17) = ME_Porous
        Me%TimeSeriesBuffer(18) = ME_Drainage
        
        !Flux among modules
        Me%TimeSeriesBuffer(19) = Me%MB%Infiltration
        Me%TimeSeriesBuffer(20) = Me%MB%OLFlowToRiver
        Me%TimeSeriesBuffer(21) = Me%MB%GWFlowToRiver
       

        call WriteTimeSerieLine (Me%ObjTimeSerieBasin, Me%TimeSeriesBuffer, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR08'     
        
        
    end subroutine GlobalMassBalance

    !--------------------------------------------------------------------------

    subroutine PredictNewDT (NewDT)

        !Arguments---------------------------------------------------------------
        real                                        :: NewDT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL     
        real                                        :: AtmosfereDT, DTForNextEvent
        real                                        :: DNetDT, RunOffDT
        real                                        :: PorousMediaDT
        integer                                     :: ID_DT   

        if (Me%Coupled%Atmosphere) then
            call GetAtmosphereDTPrediction(Me%ObjAtmosphere, AtmosfereDT, DTForNextEvent, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PredictNewDT - ModuleBasin - ERR00'
        else
            AtmosfereDT     = -null_real
            DTForNextEvent  = -null_real
        endif

        if (Me%Coupled%DrainageNetwork) then
            call GetNextDrainageNetDT   (Me%ObjDrainageNetwork, DNetDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PredictNewDT - ModuleBasin - ERR01'
        else
            DNetDT = -null_real
        endif

        if (Me%Coupled%RunOff) then
            call GetNextRunOffDT        (Me%ObjRunOff, RunOffDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PredictNewDT - ModuleBasin - ERR02'
        else
            RunOffDT = -null_real
        end if

        if (Me%Coupled%PorousMedia) then
!            call GetNextPorousMediapropDT   (Me%ObjPorousMediaProperties, PorousMediaDT, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'PredictNewDT - ModuleBasin - ERR03'
            call GetNextPorousMediaDT   (Me%ObjPorousMedia, PorousMediaDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PredictNewDT - ModuleBasin - ERR03'

            if (Me%Coupled%PorousMediaProperties) then
!                call GetNextPorousMediaPropDT   (Me%ObjPorousMediaProperties, PorousMediaPropDT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'PredictNewDT - ModuleBasin - ERR04'
            endif

        else
            PorousMediaDT = -null_real
        end if

        ID_DT = 0
        if (AtmosfereDT < NewDT) then
            NewDT = AtmosfereDT
            ID_DT = 1
        endif
        
        if (DNetDT      < NewDT) then
            NewDT = DNetDT
            ID_DT = 2
        endif

        if (RunOffDT    < NewDT) then
            NewDT = RunOffDT
            ID_DT = 3
        endif

        if (PorousMediaDT < NewDT) then
            NewDT = PorousMediaDT 
            ID_DT = 4
        endif

        if (DTForNextEvent == 0.0) then
            if (Me%DTDuringRain < NewDT) then
                NewDT = Me%DTDuringRain
                ID_DT = 5
            endif
        else
            if (DTForNextEvent < NewDT) then
                NewDT = DTForNextEvent
                ID_DT = 6
            endif
        endif
        
        !if (ID_DT == 0) then
        !    call GetMaxComputeTimeStep(Me%ObjTime, MaxDT, STAT = STAT_CALL)
        !    NewDT = MaxDT
        !endif
        
        if (NewDT < 10.0) then
            call WriteDTLog ('ModuleBasin < 10', ID_DT, NewDT)
            NewDT = 10.0
            ID_DT = 7
        endif

        call WriteDTLog ('ModuleBasin', ID_DT, NewDT)

    end subroutine PredictNewDT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillBasin(ObjBasinID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjBasinID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           
        integer                             :: STAT_CALL     
        type(T_BasinProperty),  pointer     :: PropertyX => null()

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBasinID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBASIN_,  Me%InstanceID)

            if (nUsers == 0) then

                !Writes file with final condition
                call WriteFinalFile

                nUsers = DeassociateInstance(mTIME_,  Me%ObjTime)
                if (nUsers == 0)           stop 'KillBasin - ModuleBasin - ERR01'

                if (Me%Coupled%Runoff) then
                    if (Me%Coupled%RunoffProperties) then
                        call KillRunoffProperties (Me%ObjRunoffProperties, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR03d'
                    endif                    
                    
                    call KillRunOff         (Me%ObjRunOff,      STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR02'
                    
                endif
                
                if (Me%Coupled%Vegetation) then
                    call KillVegetation  (Me%ObjVegetation, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR03a'
                endif

                if (Me%Coupled%PorousMedia) then 
                    if (Me%Coupled%PorousMediaProperties) then
                        call KillPorousMediaProperties (Me%ObjPorousMediaProperties, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR03d'
                    endif
                    
                    call KillPorousMedia (Me%ObjPorousMedia, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR03'

                endif
                                
                if (Me%Coupled%DrainageNetwork) then
                    call KillDrainageNetwork(Me%ObjDrainageNetwork, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR04'
                endif


                if (Me%Coupled%Atmosphere) then
                    call KillAtmosphere     (Me%ObjAtmosphere,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR05'
                endif                
               
                PropertyX => Me%FirstProperty
                do while(associated(PropertyX))  
                    if (PropertyX%ID%ObjFillMatrix /= 0) then
                        call KillFillMatrix (PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR04'
                    endif
                    PropertyX => PropertyX%Next
                end do

                call KillHorizontalMap      (Me%ObjHorizontalMap,   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR06'

                call KillBasinGeometry      (Me%ObjBasinGeometry,   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR07'

                call KillGridData           (Me%ObjGridData,        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR08'

                call KillHorizontalGrid     (Me%ObjHorizontalGrid,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR09'

                if (Me%Output%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillModuleBasin - ModuleBasin - ERR10'
                endif            

                !Kills the TimeSerie
                if (Me%ObjTimeSerie > 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillModuleBasin - ModuleBasin - ERR11'
                endif

                if ((Me%ObjTimeSerieBasin > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerieBasin, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillModuleBasin - ModuleBasin - ERR12'
                endif
                
                if (Me%DailyFlow%ObjTimeSerie > 0) then
                    call KillTimeSerie(Me%DailyFlow%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillModuleBasin - ModuleBasin - ERR13'
                endif
                    
                if (Me%MonthlyFlow%ObjTimeSerie > 0) then
                    call KillTimeSerie(Me%MonthlyFlow%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillModuleBasin - ModuleBasin - ERR14'
                endif


                deallocate (Me%WaterLevel)

                !Deallocates Instance
                call DeallocateInstance ()

                ObjBasinID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillBasin

    !------------------------------------------------------------------------

    subroutine WriteFinalFile
        
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                                        :: Year_File, Month_File, Day_File
        real                                        :: Hour_File, Minute_File, Second_File
        integer                                     :: FinalFile
        integer                                     :: STAT_CALL
        character(LEN = PathLength)                 :: FileName
        
        !----------------------------------------------------------------------

        !Gets Date
        call ExtractDate(Me%CurrentTime, Year_File, Month_File, Day_File,               &
                         Hour_File, Minute_File, Second_File)
        
        
        if (Me%CurrentTime == Me%EndTime) then
            FileName = Me%Files%FinalFile
        else
            FileName = ChangeSuffix(Me%Files%FinalFile,                                 &
                            "_"//trim(TimeToString(Me%CurrentTime))//".fin")
        endif            
        
        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleBasin - ERR01'

        open(Unit = FinalFile, File = FileName, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleBasin - ERR02'

        !Writes Date
        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File,       &
                         Second_File

        write(FinalFile)Me%WaterLevel
        if (Me%Coupled%Vegetation) then
            write(FinalFile)Me%CanopyStorage
        endif

        write(FinalFile)Me%AccInfiltration
        write(FinalFile)Me%AccEVTP
        write(FinalFile)Me%AccRainFall
        write(FinalFile)Me%AccEVPCanopy
        write(FinalFile)Me%AccFlowProduction


        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalFile - ModuleBasin - ERR03'

    end subroutine WriteFinalFile

    !------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Basin), pointer          :: AuxObjBasin
        type (T_Basin), pointer          :: PreviousObjBasin

        !Updates pointers
        if (Me%InstanceID == FirstObjBasin%InstanceID) then
            FirstObjBasin => FirstObjBasin%Next
        else
            PreviousObjBasin => FirstObjBasin
            AuxObjBasin      => FirstObjBasin%Next
            do while (AuxObjBasin%InstanceID /= Me%InstanceID)
                PreviousObjBasin => AuxObjBasin
                AuxObjBasin      => AuxObjBasin%Next
            enddo

            !Now update linked list
            PreviousObjBasin%Next => AuxObjBasin%Next

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

    subroutine Ready (ObjBasin_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBasin_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBasin_ID > 0) then
            call LocateObjBasin (ObjBasin_ID)
            ready_ = VerifyReadLock (mBASIN_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjBasin (ObjBasinID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBasinID

        !Local-----------------------------------------------------------------

        Me => FirstObjBasin
        do while (associated (Me))
            if (Me%InstanceID == ObjBasinID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleBasin - LocateObjBasin - ERR01'

    end subroutine LocateObjBasin

    !----------------------------------------------------------------------
    
    subroutine SearchProperty(PropertyX, PropertyXIDNumber, PrintWarning, STAT)


        !Arguments-------------------------------------------------------------
        type(T_BasinProperty), optional, pointer    :: PropertyX
        integer         , optional, intent (IN)     :: PropertyXIDNumber
        logical,          optional, intent (IN)     :: PrintWarning
        integer         , optional, intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------

        integer                                     :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

do2 :   do while (associated(PropertyX)) 
if5 :       if (PropertyX%ID%IDNumber==PropertyXIDNumber) then
                exit do2 
            else
                PropertyX => PropertyX%Next                 
            end if if5
        end do do2

       !A PropertyX was found
       if (associated(PropertyX)) then
            STAT_ = SUCCESS_  
        else
            if (present(PrintWarning)) then
                if (PrintWarning) write (*,*)'Property Not Found in Module Basin ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchProperty

    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar(LockToWhichModules)
        
        !Arguments-------------------------------------------------------------
        character (Len = StringLength)              :: LockToWhichModules

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        if (LockToWhichModules == 'AllModules') then

            !Time Stuff
            call GetComputeCurrentTime  (Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR01'

            call GetComputeTimeStep     (Me%ObjTime, Me%CurrentDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR02'

            !Gets Basin Points
            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR03'

            !Gets River Points
            call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR04'

            !Gets Grid Cell Area
            call GetGridCellArea(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR05'

            !Gets a pointer to Topography
            call GetGridData        (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR06'

            !Gets a pointer to OpenPoints2D
            call GetOpenPoints2D  (Me%ObjHorizontalMap, Me%ExtVar%OpenPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR07'
        
            call GetBoundaries    (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR07a'

        elseif (LockToWhichModules == 'Atmosphere') then        
        
            !Leaf Area Index
            call GetLeafAreaIndex(Me%ObjVegetation, Me%ExtVar%LeafAreaIndex, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR08'
            
            !Specific Storage Capacity
            call GetSpecificLeafStorage(Me%ObjVegetation, Me%ExtVar%SpecificLeafStorage, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR09'
            
            !EVTP CropCoefficient
            call GetEVTPCropCoefficient(Me%ObjVegetation, Me%ExtVar%CropCoefficient, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR10'
            
        endif
        

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar (UnLockToWhichModules)
        
        !Arguments-------------------------------------------------------------
        character (Len = StringLength)              :: UnLockToWhichModules

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        if (UnLockToWhichModules == 'AllModules') then

            !UnGets Basin Points
            call UnGetBasin         (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR01'

            !UnGets River Points
            call UnGetBasin         (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR02'

            !UnGets Grid Cell Area
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR03'

            !UnGets Topography
            call UnGetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR04'

            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%OpenPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR05'
    
            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR05a'


        elseif (UnLockToWhichModules == 'Atmosphere') then        

            !Ungets Leaf Area Index
            call UngetVegetation(Me%ObjVegetation, Me%ExtVar%LeafAreaIndex, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR06'
            
            !Specific Storage Capacity
            call UngetVegetation(Me%ObjVegetation, Me%ExtVar%SpecificLeafStorage, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR07'
            
            !EVTP CropCoefficient
            call UngetVegetation(Me%ObjVegetation, Me%ExtVar%CropCoefficient, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR08'

        endif

    end subroutine ReadUnLockExternalVar

    !--------------------------------------------------------------------------

end module ModuleBasin

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








