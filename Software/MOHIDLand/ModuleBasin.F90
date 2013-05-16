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
                                     TimeToString, ChangeSuffix, CHUNK_J, SetMatrixValue
                                     
    use ModuleFillMatrix,     only : ConstructFillMatrix, ModifyFillMatrix,              &
                                     KillFillMatrix,GetIfMatrixRemainsConstant,          &
                                     GetDefaultValue
                                     
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
                                     GetFlowAtBoundary, UnGetRunOff, KillRunOff,         &
                                     SetBasinColumnToRunoff, GetRunoffWaterColumn,       &
                                     GetRunoffWaterColumnOld, GetRunoffWaterLevel,       &
                                     GetRunoffTotalStoredVolume, GetMassError,           &
                                     GetRunOffStoredVolumes, GetRunOffBoundaryFlowVolume
                                     
                                     
    use ModuleRunoffProperties,                                                          &
                              only : ConstructRunoffProperties,                          &
                                     ModifyRunoffProperties,                             &
                                     KillRunoffProperties, GetRPConcentration,           &
                                     GetRPnProperties, GetRPOptions,                     &
                                     GetRPPropertiesIDByIdx, SetDNConcRP,                &
                                     UngetRunoffProperties, SetBasinConcRP,              &
                                     SetBasinToRPSplash, GetRPMassBalance,               &
                                     CheckRPProperty, GetRPConcentrationOld,             &
                                     GetRPDecayRate, SetVegetationRP, GetRPConcentrationAT
                                     
    use ModuleDrainageNetwork,only : ConstructDrainageNetwork, FillOutPutMatrix,         &
                                     ModifyDrainageNetwork,                              &
                                     GetHasProperties, GetDNnProperties,                 &
                                     GetHasToxicity,GetNeedsRadiation,                   &
                                     GetNeedsAtmosphere, SetAtmosphereDrainageNet,       &
                                     GetDNPropertiesIDByIdx, GetNextDrainageNetDT,       &
                                     GetVolumes, GetPropHasBottomFluxes,                 &
                                     GetChannelsNodeLength,GetChannelsID,                &
                                     SetPMPConcDN,SetRPConcDN, UnGetDrainageNetwork,     &
                                     KillDrainageNetwork, SetGWFlowLayersToDN,           &
                                     GetDNMassBalance, CheckDNProperty,                  &
                                     GetDNConcentration, GetDNStoredVolume
                                     
    use ModulePorousMedia,    only : ConstructPorousMedia, ModifyPorousMedia,            &
                                     KillPorousMedia, GetGWFlowToChannels,               &
                                     GetInfiltration, GetEfectiveEVTP,                   &
                                     GetPorousMediaTotalStoredVolume, GetEvaporation,    &
                                     GetNextPorousMediaDT, UngetPorousMedia, GetGWLayer, &
                                     GetGWFlowOption, GetGWFlowToChannelsByLayer,        &
                                     GetGWToChannelsLayers, GetIgnoreWaterColumnOnEVAP,  &
                                     GetPMStoredVolume, GetEVTPVolumes,                  &
                                     GetPMBoundaryFlowVolume
                                     
    use ModulePorousMediaProperties,                                                     &
                              only : ConstructPorousMediaProperties,                     &
                                     ModifyPorousMediaProperties,                        &
                                     KillPorousMediaProperties, GetPMPConcentration,     &
                                     SetVegetationPMProperties, GetPMPCoupled,           &
                                     SetWindVelocity, GetPMPnProperties,                 &
                                     GetPMPPropertiesIDByIdx, SetDNConcPMP,              &
                                     UngetPorousMediaProperties, CheckPMPProperty,       &
                                     SetInfColConcPMP, GetPMPMassBalance, GetECw,        &
                                     GetPMPConcentrationOld
                                     
    use ModuleVegetation,     only : ConstructVegetation, ModifyVegetation,              &
                                     KillVegetation, GetLeafAreaIndex,                   &
                                     GetSpecificLeafStorage, GetEVTPCropCoefficient,     &
                                     GetTranspiration, GetVegetationSoilFluxes,          &
                                     SetSoilConcVegetation, GetVegetationOptions,        &
                                     GetVegetationDT, GetRootDepthOld, GetNutrientFraction, &
                                     UnGetVegetation, UnGetVegetationSoilFluxes,         &
                                     GetCanopyHeight, GetTranspirationBottomLayer,       &
                                     GetPotLeafAreaIndex, GetCanopyStorageType, SetECw,  &
                                     GetVegetationAerialFluxes, GetVegetationGrowing,    &
                                     UnGetVegetationAerialFluxes         
                                     
    use ModuleStopWatch,      only : StartWatch, StopWatch
    
    use ModuleGeometry,       only : GetGeometrySize
#ifdef _ENABLE_CUDA
    use ModuleCuda
#endif _ENABLE_CUDA    
    
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
    private ::          ConstructEVTPHDFOutput     
    private ::      ConstructTimeSeries
    private ::      ReadInitialFile

    !Selector
                    
    
    !Modifier
    public  :: ModifyBasin
    private ::      AtmosphereProcesses
    private ::          DividePrecipitation
    private ::          CalcPotEvapoTranspiration
    private ::              AdjustCropCoefficient
    private ::      VegetationProcesses
    private ::      OverLandProcesses
    private ::      DrainageNetworkProcesses
    private ::      PorousMediaProcesses
    private ::      PorousMediaPropertiesProcesses
    private ::      ComputeBasinWaterBalance
    private ::      HDF5Output
    private ::      EVTPHDFOutput
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
    !Surface evaporation method
    integer, parameter                              :: LatentHeatMethod         = 1
    integer, parameter                              :: ET0Method                = 2
    integer, parameter                              :: NoEvaporation            = 3
    !Initial or Final time instant, used on Basin Water Balance
    integer, parameter                              :: InitialInstant = 1
    integer, parameter                              :: FinalInstant   = 2
    !Gw link between porous media and drainage network
!    integer, parameter                              :: Layer_ = 2

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
        real   , dimension(:,:), pointer            :: PotLeafAreaIndex 
        real   , dimension(:,:), pointer            :: SpecificLeafStorage
        real   , dimension(:,:), pointer            :: CropCoefficient
        real   , dimension(:,:,:), pointer          :: ActualTranspiration
        real   , dimension(:,:), pointer            :: MassError
    end type T_ExtVar
    
    !External variables from Runoff but that will be updated and sent to Runoff
    type T_ExtUpdate
        real(8), dimension(:,:), pointer            :: WaterLevel             => null() 
        real(8), dimension(:,:), pointer            :: WaterColumn            => null() 
        real(8), dimension(:,:), pointer            :: WaterColumnOld         => null()
    end type T_ExtUpdate

    type T_Files
        character(len=PathLength)                   :: ConstructData
        character(len=PathLength)                   :: TopographicFile
        character(len=PathLength)                   :: HDFFile
        character(len=PathLength)                   :: EVTPHDFFile
        character(len=PathLength)                   :: EVTPHDFFile2
        character(len=PathLength)                   :: EVTPInstHDFFile
        character(len=PathLength)                   :: InitialFile
        character(len=PathLength)                   :: FinalFile
        character(len=PathLength)                   :: TimeSerieLocation
        character(len=PathLength)                   :: BWBTimeSeriesLocation
    end type T_Files

    type T_IntegratedFlow
        logical                                     :: On
        real                                        :: Flow
        integer                                     :: CurrentIndex
        integer                                     :: ObjTimeSerie         = 0
    end type T_IntegratedFlow

    type T_PropMassBalance
        !All masses in kg
        !Initial and final volumes
        real(8)                                     :: InitialRPStoredMass           = 0.
        real(8)                                     :: InitialVegetationStoredMass   = 0.
        real(8)                                     :: InitialPMPStoredMass          = 0.
        real(8)                                     :: InitialDNStoredMass           = 0.
        real(8)                                     :: FinalRPStoredMass             = 0.
        real(8)                                     :: FinalVegetationStoredMass     = 0.
        real(8)                                     :: FinalPMPStoredMass            = 0.
        real(8)                                     :: FinalDNStoredMass             = 0.
        
        !PorousMediaProperties Mass Fluxes
        real(8)                                     :: PMPTranspiredMass             = 0.
        real(8)                                     :: PMPExchangeMassToDN           = 0.
        real(8)                                     :: RPExchangeMassToPMP           = 0.
        
        !RunoffProperties Mass Fluxes
        real(8)                                     :: RPExchangeMassToDN            = 0.
        
        !DrainageNetwork Mass Fluxes
        real(8)                                     :: DNDischargeMass               = 0.
        real(8)                                     :: DNOutflowMass                 = 0.
        
        !Basin Mass Flux
        real(8)                                     :: TotalRainMass                 = 0. !above leafs (covered + uncovered)
        real(8)                                     :: CoveredRainMass               = 0. !rain on leafs (covered)
        real(8)                                     :: UncoveredRainMass             = 0. !direct rain (uncovered)
        real(8)                                     :: VegDrainedMass                = 0. !leaf leak
        real(8)                                     :: RunoffInputMass               = 0. !leaf leak + direct rain
        
    end type T_PropMassBalance

    type       T_BasinProperty
        type (T_PropertyID)                         :: ID
        real, dimension(:,:), pointer               :: Field                => null()
        logical                                     :: Constant             = .false. !in time
        logical                                     :: ConstantInSpace      = .false. !in space
        type (T_BasinProperty), pointer             :: Next                 => null()
        logical                                     :: AdvectionDiffusion
        logical                                     :: Particulate
        logical                                     :: Decay
        logical                                     :: DecayRate
        type (T_PropMassBalance)                    :: MB
        real(8), dimension(:,:), pointer            :: VegetationConc       => null()
        real(8), dimension(:,:), pointer            :: VegetationDrainage   => null()
        real(8), dimension(:,:), pointer            :: VegetationOldMass    => null()
        logical                                     :: Inherited            = .false. 
        real                                        :: VegTotalStoredMass   = 0.       
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
        type (T_PropertyB)                          :: ImpFrac
        type (T_PropertyB)                          :: TimeWithNoWC
    end type T_SimpleInfiltration
    
    type T_WaterMassBalance
        !Initial Mass
        real(8)                                     :: IniVolumeRunoff      = 0.
        real(8)                                     :: IniVolumeVegetation  = 0.
        real(8)                                     :: IniVolumePorousMedia = 0.
        real(8)                                     :: IniVolumeChannels    = 0.
        !Sink / Sources
        real(8)                                     :: EvapFromVegetation   = 0.
        real(8)                                     :: EvapFromGround       = 0.
        real(8)                                     :: EvapFromSoil         = 0.
        real(8)                                     :: TotalRainIn          = 0. !above leafs
        real(8)                                     :: RainDirect           = 0. !on uncovered area 
        real(8)                                     :: RainInVeg            = 0. !on covered area (leafs)
        real(8)                                     :: DrainageFromVeg      = 0. !leak from covered (leafs)
        real(8)                                     :: RainRunoff           = 0. !arriving at WC (direct + drainage)      
        real(8)                                     :: DischargesIn         = 0.
        real(8)                                     :: OverTopOut           = 0.
        real(8)                                     :: OutVolumeChannel     = 0.
        real(8)                                     :: OutVolumeOverLand    = 0.
        !Final Mass
        real(8)                                     :: FinalVolumeRunoff
        real(8)                                     :: FinalVolumeVegetation
        real(8)                                     :: FinalVolumePorousMedia
        real(8)                                     :: FinalVolumeChannels
        !Flow exchange between modules
        real(8)                                     :: OLFlowToRiver
        real(8)                                     :: GWFlowToRiver
        real(8)                                     :: Infiltration
    end type T_WaterMassBalance
    
    type T_BasinWaterBalance
        real(8)  :: Rain                    = 0.0 !m3
        real(8)  :: Irrigation              = 0.0 !m3
        real(8)  :: DischargesOnSoil        = 0.0 !m3
        real(8)  :: DischargesOnSurface     = 0.0 !m3
        real(8)  :: DischargesOnChannels    = 0.0 !m3
        real(8)  :: IniStoredInSoil         = 0.0 !m3
        real(8)  :: IniStoredInChannels     = 0.0 !m3
        real(8)  :: IniStoredInLeaves       = 0.0 !m3
        real(8)  :: IniStoredInSurface      = 0.0 !m3
        real(8)  :: IniStoredInStormWater   = 0.0 !m3
        real(8)  :: FinStoredInSoil         = 0.0 !m3
        real(8)  :: FinStoredInChannels     = 0.0 !m3
        real(8)  :: FinStoredInLeaves       = 0.0 !m3
        real(8)  :: FinStoredInSurface      = 0.0 !m3
        real(8)  :: FinStoredInStormWater   = 0.0 !m3
        real(8)  :: StoredInSoil            = 0.0 !m3
        real(8)  :: StoredInChannels        = 0.0 !m3
        real(8)  :: StoredInLeaves          = 0.0 !m3
        real(8)  :: StoredInSurface         = 0.0 !m3
        real(8)  :: StoredInStormWater      = 0.0 !m3        
        real(8)  :: OutletFlowVolume        = 0.0 !m3
        real(8)  :: EvapFromSoil            = 0.0 !m3
        real(8)  :: EvapFromLeaves          = 0.0 !m3
        real(8)  :: EvapFromChannels        = 0.0 !m3
        real(8)  :: EvapFromSurface         = 0.0 !m3
        real(8)  :: Transpiration           = 0.0 !m3
        real(8)  :: BoundaryFromSoil        = 0.0 !m3
        real(8)  :: BoundaryFromSurface     = 0.0 !m3
        real(8)  :: Input                   = 0.0 !m3 rain and irrigation
        real(8)  :: Stored                  = 0.0 !m3 on soil, on leaves, on channels, on surface
        real(8)  :: Output                  = 0.0 !m3 evap from soil, from leaves, from surface, from channels and transpiration + outlet flow
        real(8)  :: Discharges              = 0.0 !m3 on Soil, on surface, on channels (can be positive or negative)
        real(8)  :: Boundary                = 0.0 !m3 on Soil and Surface
        real(8)  :: ErrorInVolume           = 0.0 !m3  
        real(8)  :: AccErrorInVolume        = 0.0 !m3
        real(8)  :: ErrorInPercentage       = 0.0 !% of error from the actual stored (variation from initial to final) water content to the expected due the inputs/outputs/discharges
        logical  :: StoreInitial            = .true.
        !real(8)  :: BasinArea               = 0.0 !m2 -> this should be a line in the TS header. 
        !integer  :: NumberOfCells           = 0
    end type T_BasinWaterBalance

    type       T_Basin
        integer                                     :: InstanceID           = 0
        character(len=StringLength)                 :: ModelName
        type (T_Size2D)                             :: Size, WorkSize
        type (T_Coupling)                           :: Coupled
        type (T_Files)                              :: Files
        type (T_ExtVar)                             :: ExtVar
        type (T_ExtUpdate)                          :: ExtUpdate
        type (T_SimpleInfiltration)                 :: SI                                       
        type (T_BasinProperty), pointer             :: FirstProperty        => null()
        logical                                     :: Continuous           = .false.
        logical                                     :: StopOnWrongDate      = .true.
        logical                                     :: VerifyGlobalMass     = .false.
        logical                                     :: VerifyAtmosphereValues = .true.
        logical                                     :: Calibrating1D        = .false.
        logical                                     :: ConcentrateRain      = .false.
        logical                                     :: EvapFromWaterColumn
        logical                                     :: EvapFromCanopy
        integer                                     :: EvapMethod
        real                                        :: RefEvapotranspirationConstant
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
        real(8), dimension(:,:), pointer            :: WaterColumnRemoved     => null()
        real(8), dimension(:,:), pointer            :: CanopyStorageCapacity  => null()
        real(8), dimension(:,:), pointer            :: CanopyStorage          => null() !mH20 in m2 cov
        real(8), dimension(:,:), pointer            :: CanopyStorageOld       => null() !mH20 in m2 cov
        real(8), dimension(:,:), pointer            :: CanopyDrainage         => null() !mH20 in m2 cell
        real(8), dimension(:,:), pointer            :: SnowPack               => null() 
        real(8), dimension(:,:), pointer            :: ThroughFall            => null() !mH20 in m2 cell
        real(8), dimension(:,:), pointer            :: RainUncovered          => null() !mH20 in m2 cov
        real(8), dimension(:,:), pointer            :: RainCovered            => null() !mH20 in m2 uncov
        real   , dimension(:,:), pointer            :: CoveredFraction        => null()
        real   , dimension(:,:), pointer            :: CoveredFractionOld     => null()
        real   , dimension(:,:), pointer            :: CropEvapotrans         => null()
        real   , dimension(:,:), pointer            :: PotentialTranspiration => null() 
        real   , dimension(:,:), pointer            :: PotentialEvaporation   => null() 
        real(8), dimension(:,:), pointer            :: PotentialInfCol        => null()
        real(8), dimension(:,:), pointer            :: FlowProduction         => null()
        real(8), dimension(:,:), pointer            :: InfiltrationRate       => null()
        real(8), dimension(:,:), pointer            :: PrecipRate             => null()
        real(8), dimension(:,:), pointer            :: ThroughRate            => null()
        real(8), dimension(:,:), pointer            :: EVTPRate               => null()
        !real(8), dimension(:,:), pointer            :: EVTPRate2              => null()
        real(8), dimension(:,:), pointer            :: PlantWaterStress       => null()
        real(8), dimension(:,:), pointer            :: AccInfiltration        => null()
        real(8), dimension(:,:), pointer            :: AccFlowProduction      => null()
        real(8), dimension(:,:), pointer            :: AccEVTP                => null()
        real(8), dimension(:,:), pointer            :: PartialAccEVTP         => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVTPRef      => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVTPCrop     => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVPot        => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVAct        => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccETPot        => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccETAct        => null() !mm
        real(8), dimension(:,:), pointer            :: AccEVTP2               => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVTP2        => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVTPRef2     => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVTPCrop2    => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVPot2       => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccEVAct2       => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccETPot2       => null() !mm
        real(8), dimension(:,:), pointer            :: PartialAccETAct2       => null() !mm      
        real(8), dimension(:,:), pointer            :: AccRainFall            => null()
        real(8), dimension(:,:), pointer            :: AccEVPCanopy           => null()
        real,    dimension(:,:), pointer            :: AccRainHour            => null()
        real,    dimension(:,:), pointer            :: RainStartTime          => null()
        real,    dimension(:,:), pointer            :: RainDuration           => null()
        real,    dimension(:,:), pointer            :: DiffuseFlow            => null()
!        real                                        :: WaterColumnCoef
        real                                        :: ETConversionFactor       = 1
        type (T_WaterMassBalance)                   :: MB
        type (T_IntegratedFlow)                     :: DailyFlow            
        type (T_IntegratedFlow)                     :: MonthlyFlow            
        real, dimension(:), pointer                 :: TimeSeriesBuffer  !Water Balance
        real, dimension(:), pointer                 :: TimeSeriesBuffer2
        real, dimension(:), pointer                 :: TimeSeriesBuffer3 !Properties Balance
        real, dimension(:), pointer                 :: TimeSeriesBuffer4 !Water Fluxes
        real, dimension(:), pointer                 :: TimeSeriesBuffer5 !Properties Fluxes
        real, dimension(:), pointer                 :: BWBBuffer !buffer to be used for Basin Water Balance
        
        !Basin Water Balance
        type (T_BasinWaterBalance)                  :: BWB
        logical                                     :: ComputeBasinWaterBalance = .false.
               
        !Basin is responsable by Total vegetation volume
        real(8)                                     :: VolumeVegetation
        
        type (T_Time)                               :: CurrentTime
        type (T_Time)                               :: BeginTime
        type (T_Time)                               :: EndTime       
        type (T_OutPut)                             :: OutPut
        type (T_OutPut)                             :: EVTPOutPut
        type (T_OutPut)                             :: EVTPOutPut2
        type (T_OutPut)                             :: EVTPInstOutPut
        !type (T_OutPut)                             :: EVTPRefOutPut
        
        real                                        :: DefaultKcWhenLAI0        = 0.3
        logical                                     :: UseDefaultKcWhenLAI0     = .false.
        logical                                     :: UsePotLAI                = .false.
        real                                        :: KcMin                    = 0.3
        logical                                     :: UseKcMin                 = .false.
        logical                                     :: UseRefEVTPIfNeeded       = .true.
        
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
        integer                                     :: ObjEVTPHDF               = 0
        integer                                     :: ObjEVTPHDF2              = 0
        integer                                     :: ObjEVTPInstHDF           = 0
        !integer                                     :: ObjEVTPRefHDF            = 0
        integer                                     :: ObjTimeSerie             = 0
        integer                                     :: ObjTimeSerieBasin        = 0
        integer                                     :: ObjTimeSerieBasin2       = 0
        integer                                     :: ObjTimeSerieBasinMass    = 0
        integer                                     :: ObjTimeSerieBasinMass2   = 0
        integer                                     :: ObjPorousMediaProperties = 0
        integer                                     :: ObjBWB                   = 0 !Basin Water Balance Timeseries ID
#ifdef _ENABLE_CUDA
        integer                                     :: ObjCuda                  = 0
#endif _ENABLE_CUDA
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
        integer                                         :: i, j
        integer                                         :: STAT_, STAT_CALL
        logical                                         :: VariableDT
        character (Len = StringLength)                  :: WarningString
        character (Len = StringLength)                  :: LockToWhichModules
        character (Len = StringLength)                  :: UnLockToWhichModules
        character (Len = StringLength)                  :: OptionsType
#ifdef _ENABLE_CUDA        
        type (T_Size3D)                                 :: GeometrySize
#endif _ENABLE_CUDA        
        !------------------------------------------------------------------------
!        integer                                     :: TID, OMP_GET_THREAD_NUM


        !!call OMP_SET_DYNAMIC(.false.)
        !!call OMP_SET_NUM_THREADS(2)

        !!!! $OMP PARALLEL SHARED(Me) PRIVATE(TID)
        !!!!TID = OMP_GET_THREAD_NUM()
        !!!!PRINT *, 'Hello World from thread = ', TID
        !!!! $OMP END PARALLEL         


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
                                        Me%ExtVar%BasinPoints, 2, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR08'
            
            !Gets the size of the grid
            call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                        Size     = Me%Size,                              &
                                        WorkSize = Me%WorkSize,                          &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR09a'              

            call GetGridCellArea(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR09b'            
            
            !Me%BWB%BasinArea     = 0.0
            !Me%BWB%NumberOfCells = 0
            !do j = Me%Size%JLB, Me%Size%JUB
            !do i = Me%Size%ILB, Me%Size%IUB
            !    if (Me%ExtVar%BasinPoints(i, j) == 1) then
            !        Me%BWB%BasinArea     = Me%BWB%BasinArea + Me%ExtVar%GridCellArea(i, j)
            !        Me%BWB%NumberOfCells = Me%BWB%NumberOfCells + 1
            !    endif
            !enddo
            !enddo                        
            
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR09c'			            
            
            !Ungets BasinPoints
            call UngetBasin             (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR09d'                     

#ifdef _ENABLE_CUDA            
            ! Construct a ModuleCuda instance.
            call ConstructCuda(Me%ObjCuda, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR20'

            call GetGeometrySize (Me%ObjGeometry,             &    
                                  WorkSize = GeometrySize,    &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR21'
        
            ! JPW: Initialize a C++ Thomas instance for CUDA. 
            ! This will allocate device memory for all Thomas variables (D, E, F, TI, Res)
            ! Do this seperately from ConstructCuda, because later on ModuleCuda might be used for things other than Thomas algorithm
            call InitializeThomas(Me%ObjCuda, GeometrySize, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR22'
#endif _ENABLE_CUDA

            !Gets ExternalVars
            LockToWhichModules = 'AllModules'
            OptionsType = 'ConstructBasin'
            call ReadLockExternalVar (LockToWhichModules, OptionsType)

            !Reads Data File
            call ReadDataFile

            !Allocates Variables
            call AllocateVariables      ()

            !Verifies User Options
            OptionsType = "GlobalOptions"
            call VerifyOptions(OptionsType)

            !Constructs Coupled Modules
#ifdef _ENABLE_CUDA
            call ConstructCoupledModules(Me%ObjCuda)
#else
            call ConstructCoupledModules()
#endif _ENABLE_CUDA
            
            !Checks property related options
            if (Me%Coupled%RunoffProperties) then
                WarningString = "PropertyOptions"
                call VerifyOptions(WarningString)
            endif

            !Constructs the property list 
            call ConstructPropertyList

            !Constructs Output Time Series
            call ConstructTimeSeries
            
            !Reads conditions from previous run
            if (Me%Continuous) call ReadInitialFile
            
            
!            !Gets water level to control water column updates
!            call GetRunoffWaterLevel(Me%ObjRunOff, Me%WaterLevel, STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR00'            
!            
!            !Actualizes the WaterColumn
!            OptionsType = 'ConstructBasin'
!            call ActualizeWaterColumn (OptionsType)
!            
!            call UngetRunoff (Me%ObjRunOff, Me%WaterLevel, STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR10'     


            !Closes Data File
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructBasin - ModuleBasin - ERR99'

            !compute vegetation mass on top of foils
            if (Me%Coupled%Vegetation) then
                call CalculateVegTotalStoredMass
            endif

            !First HDF 5 Output
            if (Me%Output%Yes) then
                call HDF5OutPut       
            endif


            
            !UnGets ExternalVars
            UnLockToWhichModules = 'AllModules'
            OptionsType = 'ConstructBasin'
            call ReadUnLockExternalVar (UnLockToWhichModules, OptionsType)

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

        nullify  (NewObjBasin%ExtUpdate%WaterLevel)
        nullify  (NewObjBasin%ExtUpdate%WaterColumn)
        nullify  (NewObjBasin%Output%OutTime)
        nullify  (NewObjBasin%Output%OutputChannels)
        nullify  (NewObjBasin%EVTPOutput%OutTime)
        nullify  (NewObjBasin%EVTPOutput%OutputChannels)   
        nullify  (NewObjBasin%EVTPOutput2%OutTime)
        nullify  (NewObjBasin%EVTPOutput2%OutputChannels)   
        nullify  (NewObjBasin%EVTPInstOutput%OutTime)
        nullify  (NewObjBasin%EVTPInstOutput%OutputChannels) 

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
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR010'      
       
       !Reads file name of the topographic file
        call ReadFileName('IN_BASIN', Me%Files%TopographicFile,                         &
                           Message = "Topographic Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR020'

        !Reads file name of the hdf outupt
        call ReadFileName('BASIN_HDF', Me%Files%HDFFile,                                &
                           Message = "Basin HDF Output File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR030'             

        !Reads file name of final output file
        call ReadFileName('BASIN_FIN', Me%Files%FinalFile,                              &
                           Message = "Basin Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileNames - ModuleBasin - ERR050'

    end subroutine ReadFileNames

    !--------------------------------------------------------------------------

    subroutine ReadDataFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        real                                        :: dummy

        !Constructs the DataFile
        call ConstructEnterData (Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR010'

        !Basin Initial Water Column
        call GetData(dummy,                                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'INITIAL_WATER_COLUMN',                              &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR01'
        if (iflag /= 0) then
            write(*,*)'The keyword INITIAL_WATER_COLUMN was moved to Module Runoff'
            write(*,*)'Change the data files'
            stop 'ReadDataFile - ModuleBasin - ERR01.5'
        endif
!        !Basin Initial Water Column
!        call GetData(Me%WaterColumnCoef,                                                 &
!                     Me%ObjEnterData, iflag,                                             &
!                     SearchType   = FromFile,                                            &
!                     keyword      = 'WATER_COLUMN_COEF',                                 &
!                     default      = 0.0,                                                 & 
!                     ClientModule = 'ModuleBasin',                                       &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR01'


        !Continuous Computation
        call GetData(Me%Continuous,                                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONTINUOUS',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR040'

        call GetData(Me%StopOnWrongDate,                                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'STOP_ON_WRONG_DATE',                                &
                     default      = .true.,                                              &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR050'

        !Verify Global Mass
        call GetData(Me%VerifyGlobalMass,                                                &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'VERIFY_MASS',                                       &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR060'
        
        !Verify if do not exist negative values in atmospheric properties 
        !they may be causes of error that the user does not get aware
        call GetData(Me%VerifyAtmosphereValues,                                          &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'VERIFY_ATMOSPHERE_VALUES',                          &
                     default      = .true.,                                              &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR061'
        
        call GetData(Me%ComputeBasinWaterBalance,                                        &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'COMPUTE_WATER_BALANCE',                             &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR062'       
        
        if (Me%ComputeBasinWaterBalance) then
        !Gets TimeSerieLocationFile
            call GetData(Me%Files%BWBTimeSeriesLocation,                                 &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'TIME_SERIE_LOCATION_BWB',                       &
                         ClientModule = 'ModuleBasin',                                   &
                         Default      = Me%Files%ConstructData,                          &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR063'        
        
        endif

        !Calibrating 1D column?
        call GetData(Me%Calibrating1D,                                                   &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CALIBRATING_1D',                                    &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR070'

        !How much to remove for the watercolumn 
        if (Me%Calibrating1D) then
            call GetData(Me%WCRemovalTime,                                                &
                         Me%ObjEnterData, iflag,                                          &
                         SearchType   = FromFile,                                         &
                         keyword      = 'WC_REMOVE_TIME',                                 &
                         default      = 600.,                                             &
                         ClientModule = 'ModuleBasin',                                    &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR080'
        endif

        !DT During Rain
        call GetData(Me%DTDuringRain,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DT_DURING_RAIN',                                    &
                     default      = 60.0,                                                &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR090'

        !Concentrate Rain for subhourly steps
        call GetData(Me%ConcentrateRain,                                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONCENTRATE_RAIN',                                  &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR100'
        
        !How long does rain last...
        if (Me%ConcentrateRain) then
            call GetData(Me%RainAverageDuration,                                         &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'RAIN_AVERAGE_DURATION',                         &
                         default      = 600.00,                                          &
                         ClientModule = 'ModuleBasin',                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR110'
            
            if (1.5 * Me%RainAverageDuration + Me%DTDuringRain > 3600.0) then
                write(*,*)'1.5 * RainAverageDuration + DT During Rain cannot exceed 3600s'
                stop 'ReadDataFile - ModuleBasin - ERR0120'
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
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR130'

        call GetData(Me%EvapMethod,                                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'EVAP_METHOD',                                       &
                     default      = ET0Method,                                           & !ET0Method
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR140'

        if (Me%EvapMethod .NE. NoEvaporation) then
            !
            call GetData(Me%EvapFromWaterColumn,                                             &
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'EVAP_FROM_WATER_COLUMN',                            &
                         default      = .false.,                                              &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR150'
           
            !
            call GetData(Me%EvapFromCanopy,                                                  &
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'EVAP_FROM_CANOPY',                                  &
                         default      = .true.,                                              &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR160'
        
            if ((.NOT. Me%EvapFromWaterColumn) .AND. (.NOT. Me%EvapFromCanopy)) then
            
                Me%EvapMethod = NoEvaporation 
                                
            endif
        
        else
        
            Me%EvapFromWaterColumn = .false.
            Me%EvapFromCanopy      = .false.
        
        endif
        
        if (Me%DiffuseWaterSource) then
            call GetData(Me%FlowPerCapita,                                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'FLOW_PER_CAPITA',                               &
                         default      = 0.0,                                             &
                         ClientModule = 'ModuleBasin',                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR170'

            call GetData(Me%PopDensityFile,                                              &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromFile,                                        &
                         keyword      = 'POPULATION_DENSITY',                            &
                         ClientModule = 'ModuleBasin',                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR180'
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
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR190'

        !Reads file name of initial condition file
        if (Me%Continuous) then
            call ReadFileName('BASIN_INI', Me%Files%InitialFile,                         &
                               Message = "Basin Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR200'
        endif

        !Verifies if the user wants to use the Atmosphere Condition
        call GetData(Me%Coupled%Atmosphere,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'ATMOSPHERE',                                        &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR210'

        !Verifies if the user wants to use the Atmosphere Condition
        call GetData(Me%Coupled%Evapotranspiration,                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'EVAPOTRANSPIRATION',                                &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR220'

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
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR230'
        endif

        if ((Me%EvapoTranspirationMethod .EQ. 1) .AND. (Me%EvapMethod .NE. NoEvaporation)) then
                write(*,*)  
                write(*,*) 'If EVAPOTRANSPIRATION_METHOD = 1, then '
                write(*,*) 'EVAP_METHOD must be set to 3 (NoEvaporation)'
                stop 'ReadDataFile - ModuleBasin - ERR240'        
        endif

        !Verifies if the user wants to simulate Infiltration
        call GetData(Me%Coupled%PorousMedia,                                             &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'POROUS_MEDIA',                                      &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR250'

        if (Me%Coupled%PorousMedia) then
            ! verifies if the user wants to simulate transport of properties
            call GetData(Me%Coupled%PorousMediaProperties,                                   &
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'POROUS_MEDIA_PROPERTIES',                           &
                         default      = OFF,                                                 &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR260'
        endif

        !Verifies if the user wants to simulate OverLand RunOff
        call GetData(Me%Coupled%RunOff,                                                  & 
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'RUN_OFF',                                           &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR270'
        
        if (Me%Coupled%RunOff) then
            !Verifies if the user wants to simulate OverLand RunOff propertie transport
            call GetData(Me%Coupled%RunOffProperties,                                        & 
                         Me%ObjEnterData, iflag,                                             &
                         SearchType   = FromFile,                                            &
                         keyword      = 'RUN_OFF_PROPERTIES',                                &
                         default      = OFF,                                                 &
                         ClientModule = 'ModuleBasin',                                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR0280'
        endif

        !A Drainage Network is coupled?
        call GetData(Me%Coupled%DrainageNetwork,                                         &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DRAINAGE_NET',                                      &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR290'

        !The Vegetation Module is coupled_
        call GetData(Me%Coupled%Vegetation,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'VEGETATION',                                        &
                     default      = ON,                                                  &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR300'

        !Verifies if the user wants to use simple 
        call GetData(Me%Coupled%SimpleInfiltration,                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SIMPLE_INFILTRATION',                               &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR310'


        !Gets Output Time 
        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime = Me%CurrentTime,                                 &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR320'

        call GetOutPutTime(Me%ObjEnterData,                          &
                           CurrentTime   = Me%CurrentTime,           &
                           EndTime       = Me%EndTime,               &
                           keyword       = 'EVTP_OUTPUT_TIME',       &
                           SearchType    = FromFile,                 &
                           OutPutsTime   = Me%EVTPOutPut%OutTime,    &
                           OutPutsOn     = Me%EVTPOutPut%Yes,        &
                           STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR330'
        
        call GetOutPutTime(Me%ObjEnterData,                          &
                           CurrentTime   = Me%CurrentTime,           &
                           EndTime       = Me%EndTime,               &
                           keyword       = 'EVTP_OUTPUT_TIME2',       &
                           SearchType    = FromFile,                 &
                           OutPutsTime   = Me%EVTPOutPut2%OutTime,    &
                           OutPutsOn     = Me%EVTPOutPut2%Yes,        &
                           STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR331'
        
        call GetOutPutTime(Me%ObjEnterData,                          &
                           CurrentTime   = Me%CurrentTime,           &
                           EndTime       = Me%EndTime,               &
                           keyword       = 'EVTPINST_OUTPUT_TIME',       &
                           SearchType    = FromFile,                 &
                           OutPutsTime   = Me%EVTPInstOutPut%OutTime,    &
                           OutPutsOn     = Me%EVTPInstOutPut%Yes,        &
                           STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR332'        
        
        !Output for restart
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime  = Me%CurrentTime,                               &
                           EndTime      = Me%EndTime,                                   &
                           keyword      = 'RESTART_FILE_OUTPUT_TIME',                   &
                           SearchType   = FromFile,                                     &
                           OutPutsTime  = Me%OutPut%RestartOutTime,                     &
                           OutPutsOn    = Me%OutPut%WriteRestartFile,                   &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR340'

        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleBasin',                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ReadDataFile - ModuleBasin - ERR350'

        !Gets TimeSerieLocationFile
        call GetData(Me%Files%TimeSerieLocation,                                         &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TIME_SERIE_LOCATION',                               &
                     ClientModule = 'ModuleBasin',                                       &
                     Default      = Me%Files%ConstructData,                              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR360'

        !Output daily flow values?
        call GetData(Me%DailyFlow%On,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DAILY_FLOW',                                        &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR370'

        !Output monthly flow values?
        call GetData(Me%MonthlyFlow%On,                                                  &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MONTHLY_FLOW',                                      &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR380'

        !Verifies if the user wants to have precipitations has snow
        call GetData(Me%Coupled%Snow,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SNOW',                                              &
                     default      = OFF,                                                 &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR390'

        !Set a default value for Kc if LAI = 0
        call GetData(Me%DefaultKcWhenLAI0,                                               &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DEFAULT_KC_WHEN_LAI_ZERO',                          &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR400'
        
        if (iflag .EQ. 0) then
            Me%UseDefaultKcWhenLAI0 = .false.
        else
            Me%UseDefaultKcWhenLAI0 = .true.
        endif

        !Set a minimum value for Kc (for LAI = 0)
        call GetData(Me%KcMin,                                                           &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'KC_MIN',                                            &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR410'    
        if (iflag .EQ. 0) then
            Me%UseKcMin = .false.
        else
            Me%UseKcMin = .true.
        endif

        call GetData(Me%UseRefEVTPIfNeeded,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     keyword      = 'USE_REF_EVTP_IF_NEEDED',                            &
                     default      = .false.,                                              &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleBasin - ERR410'   

    end subroutine ReadDataFile
    !--------------------------------------------------------------------------

    subroutine VerifyOptions (WarningString)

        !Arguments-------------------------------------------------------------
        character (Len = StringLength)                   :: WarningString
        !Local-----------------------------------------------------------------
        integer                                          :: nProperties, STAT_CALL
        integer                                          :: iProp, PropID
        logical                                          :: PropAdvDiff, PropParticulate
        !Begin-----------------------------------------------------------------
        
        if (WarningString == "GlobalOptions") then
            if (Me%Coupled%DrainageNetwork .and. .not. Me%Coupled%Runoff) then
                write(*,*)'You must enable module Runoff if you want to use module Drainage Network'
                stop 'VerifyOptions - ModuleBasin - ERR01'
            endif
            
    !        if (Me%Coupled%PorousMedia .and. .not. Me%Coupled%Vegetation) then
    !            write(*,*)'You must enable module Vegetation if you want to use module Porous Media'
    !            stop 'VerifyOptions - ModuleBasin - ERR02'
    !        endif

            if (Me%Coupled%PorousMedia .and. Me%Coupled%SimpleInfiltration) then
                write(*,*)'You can not use SimpleInfiltration and PorousMedia at the same time'
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
            
            if (Me%Coupled%PorousMedia .or. Me%Coupled%SimpleInfiltration) then
                if (.not. Me%Coupled%Runoff) then
                    write(*,*)'If using Porous Media model than need to link RUN_OFF '
                    write(*,*)'Because Runoff is now the module responsible for the'
                    write(*,*)'water column. Module Basin is now only an updater'
                    stop 'VerifyOptions - ModuleBasin - ERR08' 
                endif                
            endif
            
        elseif (WarningString == "PropertyOptions") then
            !!!Check if properties that have advection diffusion have it in all modules
            !Checking with Runoff Properties
            if (Me%Coupled%RunoffProperties) then
                call GetRPnProperties (Me%ObjRunoffProperties, nProperties, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR10'

                do iProp = 1, nProperties

                    call GetRPPropertiesIDByIdx(RunoffPropertiesID = Me%ObjRunoffProperties,            &
                                                 Idx                     = iProp,                       &
                                                 ID                      = PropID,                      &
                                                 PropAdvDiff             = PropAdvDiff,                 &
                                                 Particulate             = PropParticulate,             &
                                                 STAT                    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR20' 
                    
                    !Check in PMP if the same dissolved properties also have advection diffusion connected
                    !Only dissolved properties can communicate with PMP
                    if (PropAdvDiff .and. (.not. PropParticulate)) then
                        if (Me%Coupled%PorousMediaProperties) then
                            call CheckPMPProperty (Me%ObjPorousMediaProperties, PropID, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR30' 
                        endif
                    endif
                    !Check in DN if the same properties have advection diffusion connected
                    !Dissolved and particulate properties can communicate with DN
                    if (PropAdvDiff) then
                        if (Me%Coupled%DrainageNetwork) then
                            call CheckDNProperty(Me%ObjDrainageNetwork, PropID, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR40' 
                        endif
                    endif
                enddo
            endif
            
            !Now checking with PMP
            if (Me%Coupled%PorousMediaProperties) then
                call GetPMPnProperties (Me%ObjPorousMediaProperties, nProperties, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR50'

                do iProp = 1, nProperties

                    call GetPMPPropertiesIDByIdx(PorousMediaPropertiesID = Me%ObjPorousMediaProperties, &
                                                 Idx                     = iProp,                       &
                                                 ID                      = PropID,                      &
                                                 PropAdvDiff             = PropAdvDiff,                 &
                                                 STAT                    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR60' 
                    
                    !Check in RP if the same properties also have advection diffusion connected
                    !Only Dissolved properties can communicate with RP but particulate props may not have adv-diff
                    if (PropAdvDiff) then
                        if (Me%Coupled%RunoffProperties) then
                            call CheckRPProperty (Me%ObjRunoffProperties, PropID, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR70' 
                        endif
                    endif
                    !Check in DN if the same properties have advection diffusion connected
                    !Only Dissolved properties can communicate with DN but particulate props may not have adv-diff
                    if (PropAdvDiff) then
                        if (Me%Coupled%DrainageNetwork) then
                            call CheckDNProperty(Me%ObjDrainageNetwork, PropID, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR80' 
                        endif
                    endif
                enddo
            endif

            !Now checking with DN
            if (Me%Coupled%DrainageNetwork) then
                call GetDNnProperties (Me%ObjDrainageNetwork, nProperties, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR90'

                do iProp = 1, nProperties

                    call GetDNPropertiesIDByIdx(DrainageNetworkID        = Me%ObjDrainageNetwork,       &
                                                 Idx                     = iProp,                       &
                                                 ID                      = PropID,                      &
                                                 PropAdvDiff             = PropAdvDiff,                 &
                                                 Particulate             = PropParticulate,             &
                                                 STAT                    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR100' 
                    
                    !Check in RP if the same properties also have advection diffusion connected
                    !Dissolved and particulate properties can communicate with RP
                    if (PropAdvDiff) then
                        if (Me%Coupled%RunoffProperties) then
                            call CheckRPProperty (Me%ObjRunoffProperties, PropID, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR110' 
                        endif
                    endif
                    !Check in PMP if the same properties have advection diffusion connected
                    !Only dissolved properties can communicate with PMP
                    if (PropAdvDiff .and. (.not. PropParticulate)) then
                        if (Me%Coupled%PorousMediaProperties) then
                            call CheckPMPProperty (Me%ObjPorousMediaProperties, PropID, STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'VerifyOptions - ModuleBasin - ERR120' 
                        endif
                    endif
                enddo
            endif
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
        !call ConstructOneProperty (Me%SI%InfRate,   "InfRate",  "<BeginInfRate>",   "<EndInfRate>")
        !call ConstructOneProperty (Me%SI%AccInf,    "AccInf",   "<BeginAccInf>",    "<EndAccInf>")
        call ConstructOneProperty (Me%SI%ImpFrac,   "ImpFrac",  "<BeginImpFrac>",   "<EndImpFrac>")
        
        Me%SI%InfRate%ID%Name = "InfRate"
        Me%SI%AccInf%ID%Name = "AccInf"
        Me%SI%TimeWithNoWC%ID%Name = "NotPondedTime"
    
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

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)       
        if (STAT_CALL /= SUCCESS_) stop 'ConstructSimpleInfiltration - ModuleBasin - ERR01'
        
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

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructOneProperty - ModuleBasin - ERR06'

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
        integer                                             :: nProperties, ColNumber
        real, dimension(6), target                          :: AuxTime
        integer                                             :: i, AuxInt
        character(PathLength)                               :: AuxChar
        integer                                             :: TimeSerieNumber, dn, Id, Jd
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        character(len=StringLength)                         :: TimeSerieName
        type (T_BasinProperty), pointer                     :: PropertyX
        !Begin------------------------------------------------------------------


        !Basin Water Balance====================================================
        if (Me%ComputeBasinWaterBalance) then
        
            !Time Serie for BASIN WATER BALANCE (BWB)
            allocate(PropertyList(36))
            allocate(Me%BWBBuffer(36))
            
            PropertyList(1)     = "Rain_m3"
            PropertyList(2)     = "Irrigation_m3"
            PropertyList(3)     = "DischargesOnSoil_m3"
            PropertyList(4)     = "DischargesOnSurface_m3"
            PropertyList(5)     = "DischargesOnChannels_m3"
            PropertyList(6)     = "InitialStoredInSoil_m3"
            PropertyList(7)     = "FinalStoredInSoil_m3"
            PropertyList(8)     = "StoredInSoil_m3"
            PropertyList(9)     = "InitialStoredInChannels_m3"
            PropertyList(10)    = "FinalStoredInChannels_m3"
            PropertyList(11)    = "StoredInChannels_m3"
            PropertyList(12)    = "InitialStoredInStormWater_m3"
            PropertyList(13)    = "FinalStoredInStormWater_m3"
            PropertyList(14)    = "StoredInStormWater_m3"
            PropertyList(15)    = "InitialStoredInLeaves_m3"
            PropertyList(16)    = "FinalStoredInLeaves_m3"
            PropertyList(17)    = "StoredInLeaves_m3"
            PropertyList(18)    = "InitialStoredInSurface_m3"
            PropertyList(19)    = "FinalStoredInSurface_m3"
            PropertyList(20)    = "StoredInSurface_m3"
            PropertyList(21)    = "OutletFlowVolume_m3"
            PropertyList(22)    = "EvaporationFromSoil_m3"
            PropertyList(23)    = "EvaporationFromLeaves_m3"
            PropertyList(24)    = "EvaporationFromChannels_m3"
            PropertyList(25)    = "EvaporationFromSurface_m3"
            PropertyList(26)    = "Transpiration_m3"
            PropertyList(27)    = "BoundaryFromSoil_m3"
            PropertyList(28)    = "BoundaryFromSurface_m3"
            PropertyList(29)    = "Input_m3"
            PropertyList(30)    = "Discharges_m3"
            PropertyList(31)    = "Output_m3"
            PropertyList(32)    = "Stored_m3"
            PropertyList(33)    = "Boundary_m3"            
            PropertyList(34)    = "Error_m3"
            PropertyList(35)    = "AccumulatedError_m3"
            PropertyList(36)    = "Error_%"
            !PropertyList(37)    = "NumberOfBasinCells"
            !PropertyList(38)    = "HDFAccEVTP_m3"
            
            call StartTimeSerie(Me%ObjBWB, Me%ObjTime,                                      &
                                Me%Files%BWBTimeSeriesLocation,                             &
                                PropertyList, "srbx",                                       &
                                ResultFileName = 'Basin Water Balance',                     &
                                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR01'

            deallocate(PropertyList)
        endif
        !End of Basin Water Balance==============================================
        
        !Time Serie of properties variable in the Basin 
        i = 7
        if (Me%Coupled%Vegetation) then
            i = i + 4
        
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then        
                i = i + 4
            endif
            
            PropertyX => Me%FirstProperty
            do while (associated (PropertyX))
                if (PropertyX%Inherited) then
                    i = i + 2
                endif
                PropertyX => PropertyX%Next
            enddo
            
        endif 

        if (Me%Coupled%Evapotranspiration) then
            i = i + 1
        endif

        allocate(PropertyList(i))             

        PropertyList(1)  = 'water column [m]'
        PropertyList(2)  = 'water level [m]'
        PropertyList(3)  = 'Infil. Rate [mm/hour]'
        PropertyList(4)  = 'Precipitation Rate [mm/hour]'
        PropertyList(5)  = 'Throughfall Rate [mm/hour]'
        PropertyList(6)  = 'EvapoTranspiration Rate [mm/hour]'
        PropertyList(7)  = 'Water Column Removed [m]'
        i = 8
        if (Me%Coupled%Vegetation) then
            PropertyList(i) = 'Canopy Capacity [m]'
            i = i + 1
            PropertyList(i) = 'Canopy Storage [m]'
            i = i + 1
            PropertyList(i) = 'Canopy Drainage [m]'
            i = i + 1
            PropertyList(i)  = 'Potential Crop EVTP [mm/h]'
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
            
            PropertyX => Me%FirstProperty
            do while (associated (PropertyX))
                if (PropertyX%Inherited) then
                    PropertyList(i) = 'Leaf Conc '//trim(PropertyX%ID%Name)//'[mg/l]'
                    i = i + 1
                    PropertyList(i) = 'Leaf Mass '//trim(PropertyX%ID%Name)//'[g]'
                    i = i + 1
                endif
                PropertyX => PropertyX%Next
            enddo
            
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


        !Time Serie of water balance integrated for the basin
        allocate(PropertyList(26))
        PropertyList(1)     = "Initial_Volume_Runoff_m3"
        PropertyList(2)     = "Initial_Volume_Vegetation_m3"
        PropertyList(3)     = "Initial_Volume_PorousMedia_m3"
        PropertyList(4)     = "Initial_Volume_Channels_m3"
        PropertyList(5)     = "Evap_From_Vegetation_m3"
        PropertyList(6)     = "Evap_From_Ground_m3"
        PropertyList(7)     = "Evap_From_Soil_m3"
        PropertyList(8)     = "Rain_Above_Leafs_m3"
        PropertyList(9)     = "Rain_Uncovered_m3"
        PropertyList(10)     = "Rain_Covered_m3"
        PropertyList(11)     = "Leaf_Drainage_m3"
        PropertyList(12)     = "Rain_Below_Leafs_m3"
        PropertyList(13)     = "Out_Volume_Channel_m3"
        PropertyList(14)    = "Out_Volume_Overland_m3"
        PropertyList(15)    = "Final_Volume_Runof_m3"
        PropertyList(16)    = "Final_Volume_Vegetation_m3"
        PropertyList(17)    = "Final_Volume_PorousMedia_m3"
        PropertyList(18)    = "Final_Volume_Channels_m3"
        PropertyList(19)    = "Volume_Error_Ratio_Global"
        PropertyList(20)    = "Volume_Error_Ratio_Runoff"
        PropertyList(21)    = "Volume_Error_Ratio_Porous"
        PropertyList(22)    = "Volume_Error_Ratio_DNet"
        PropertyList(23)    = "Volume_Error_Ratio_Veg"
        PropertyList(24)    = "Infiltration_m3"
        PropertyList(25)    = "OL_Volume_To_Channels_m3"
        PropertyList(26)    = "GW_Volume_To_Channels_m3"


        call StartTimeSerie(Me%ObjTimeSerieBasin, Me%ObjTime,                           &
                            Me%Files%TimeSerieLocation,                                 &
                            PropertyList, "srb",                                        &
                            ResultFileName = 'Basin Water Balance',                     &
                            STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR02'

        deallocate(PropertyList)


        !Time Serie of water balance integrated for the basin
        allocate(PropertyList(13))
        PropertyList(1)     = "Evap_From_Vegetation_m3/s"
        PropertyList(2)     = "Evap_From_Ground_m3/s"
        PropertyList(3)     = "Evap_From_Soil_m3/s"
        PropertyList(4)     = "Rain_Above_Leafs_m3/s"
        PropertyList(5)     = "Rain_Uncovered_m3/s"
        PropertyList(6)     = "Rain_Covered_m3/s"
        PropertyList(7)     = "Leaf_Drainage_m3/s"
        PropertyList(8)     = "Rain_Below_Leafs_m3/s"
        PropertyList(9)     = "Out_Volume_Channel_m3/s"
        PropertyList(10)    = "Out_Volume_Overland_m3/s"
        PropertyList(11)    = "Infiltration_m3/s"
        PropertyList(12)    = "OL_Volume_To_Channels_m3/s"
        PropertyList(13)    = "GW_Volume_To_Channels_m3/s"


        call StartTimeSerie(Me%ObjTimeSerieBasin2, Me%ObjTime,                           &
                            Me%Files%TimeSerieLocation,                                 &
                            PropertyList, "srb",                                        &
                            ResultFileName = 'Basin Water Balance Fluxes',              &
                            STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR020'

        deallocate(PropertyList)

        !Time Serie of properties mass balance integrated for the basin
        if (Me%Coupled%RunoffProperties) then
            
            nProperties = 0
            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))
                if (PropertyX%AdvectionDiffusion .and. PropertyX%Inherited) then
                    nProperties = nProperties + 1
                endif
                PropertyX => PropertyX%Next
            enddo
!            if (Me%Coupled%Evapotranspiration) then
!                nProperties = nProperties - 1
!            endif
            
            !23 outputs per property
            allocate(PropertyList(nProperties * 23))
            allocate (Me%TimeSeriesBuffer3(nProperties * 23))

            ColNumber = 1
            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))
                if (PropertyX%AdvectionDiffusion) then
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Initial_Mass_Runoff_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Initial_Mass_Vegetation_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Initial_Mass_PorousMedia_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Initial_Mass_Channels_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Transp_Mass_From_Soil_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Rain_Mass_In_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Rain_Mass_UnCovered_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Rain_Mass_Covered_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" LeafDrainage_Mass_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Drainage+Uncovered_Mass_kg"
                    ColNumber = ColNumber + 1                    
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Out_Mass_Channel_kg"
                    ColNumber = ColNumber + 1
!                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Out_Mass_Overland_kg"
!                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Final_Mass_Runoff_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Final_Mass_Vegetation_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Final_Mass_PorousMedia_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Final_Mass_Channels_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Mass_Error_Ratio_Global"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Mass_Error_Ratio_Runoff"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Mass_Error_Ratio_Porous"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Mass_Error_Ratio_DNet"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Mass_Error_Ratio_Veg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Infiltration_Mass_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" OL_Mass_To_Channels_kg"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" GW_Mass_To_Channels_kg"
                    ColNumber = ColNumber + 1
                endif
                PropertyX => PropertyX%Next
            enddo

            call StartTimeSerie(Me%ObjTimeSerieBasinMass, Me%ObjTime,                       &
                                Me%Files%TimeSerieLocation,                                 &
                                PropertyList, "srb",                                        &
                                ResultFileName = 'Basin Mass Balance',                      &
                                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR020'


            deallocate(PropertyList)
            

            !10 outputs per property
            allocate(PropertyList(nProperties * 10))
            allocate (Me%TimeSeriesBuffer5(nProperties * 10))

            ColNumber = 1
            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))
                if (PropertyX%AdvectionDiffusion) then
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Transp_Mass_From_Soil_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Rain_Mass_In_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Rain_Mass_UnCovered_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Rain_Mass_Covered_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" LeafDrainage_Mass_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Drainage+Uncovered_Mass_kg/s"
                    ColNumber = ColNumber + 1                    
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Out_Mass_Channel_kg/s"
                    ColNumber = ColNumber + 1
!                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Out_Mass_Overland_kg"
!                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" Infiltration_Mass_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" OL_Mass_To_Channels_kg/s"
                    ColNumber = ColNumber + 1
                    PropertyList(ColNumber)     = trim(PropertyX%ID%Name)//" GW_Mass_To_Channels_kg/s"
                    ColNumber = ColNumber + 1
                endif
                PropertyX => PropertyX%Next
            enddo

            call StartTimeSerie(Me%ObjTimeSerieBasinMass2, Me%ObjTime,                      &
                                Me%Files%TimeSerieLocation,                                 &
                                PropertyList, "srb",                                        &
                                ResultFileName = 'Basin Mass Balance Fluxes',               &
                                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSeries - ModuleBasin - ERR020'


            deallocate(PropertyList)
            
            
        endif
        
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

    subroutine ConstructEVTPHDFOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB   
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        
        !------------------------------------------------------------------------
       
       
        !Reads file name of the EVTP hdf outupt
        if (Me%EVTPOutput%Yes) then
            call ReadFileName('BASIN_EVTPHDF', Me%Files%EVTPHDFFile,                        &
                               Message = "Basin EVTPHDF Output File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR010'                    
            !Bounds
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB           
        
            Me%EVTPOutPut%NextOutPut = 1 
        
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

            !Opens HDF File
            call ConstructHDF5      (Me%ObjEVTPHDF, trim(Me%Files%EVTPHDFFile)//"5", HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR020'

      
            !Write the Horizontal Grid
            call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjEVTPHDF, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR030'


            !Sets limits for next write operations
            call HDF5SetLimits      (Me%ObjEVTPHDF, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR040'
        
            !Writes the Grid
            call HDF5WriteData   (Me%ObjEVTPHDF, "/Grid", "Bathymetry", "m",           &
                                    Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR050'

            !WriteBasinPoints
            call HDF5WriteData   (Me%ObjEVTPHDF, "/Grid", "BasinPoints", "-",          &
                                    Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR060'

            !Flushes All pending HDF5 commands
            call HDF5FlushMemory (Me%ObjEVTPHDF, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR070'              
        endif
        
        !Reads file name of the EVTP 2 hdf outupt
        if (Me%EVTPOutput2%Yes) then
            call ReadFileName('BASIN_EVTPHDF2', Me%Files%EVTPHDFFile2,                      &
                               Message = "Basin EVTPHDF2 Output File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR011'   
            !Bounds
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
    
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB           
        
            Me%EVTPOutPut2%NextOutPut = 1 
        
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    
            !Opens HDF File
            call ConstructHDF5      (Me%ObjEVTPHDF2, trim(Me%Files%EVTPHDFFile2)//"5", HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR021'
    
      
            !Write the Horizontal Grid
            call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjEVTPHDF2, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR031'
    
    
            !Sets limits for next write operations
            call HDF5SetLimits      (Me%ObjEVTPHDF2, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR041'
        
            !Writes the Grid
            call HDF5WriteData   (Me%ObjEVTPHDF2, "/Grid", "Bathymetry", "m",           &
                                    Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR051'
    
            !WriteBasinPoints
            call HDF5WriteData   (Me%ObjEVTPHDF2, "/Grid", "BasinPoints", "-",          &
                                    Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR061'
    
            !Flushes All pending HDF5 commands
            call HDF5FlushMemory (Me%ObjEVTPHDF2, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR071'         
         endif
        
        !Reads file name of the EVTP Instanteneous hdf outupt
        if (Me%EVTPInstOutput%Yes) then
            call ReadFileName('BASIN_EVTP_INST_HDF', Me%Files%EVTPInstHDFFile,                        &
                               Message = "Basin EVTP Instantaneous HDF Output File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR012'                    
            !Bounds
            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB

            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB           
        
            Me%EVTPInstOutPut%NextOutPut = 1 
        
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

            !Opens HDF File
            call ConstructHDF5      (Me%ObjEVTPInstHDF, trim(Me%Files%EVTPInstHDFFile)//"5", HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR022'

      
            !Write the Horizontal Grid
            call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjEVTPInstHDF, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR032'


            !Sets limits for next write operations
            call HDF5SetLimits      (Me%ObjEVTPInstHDF, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR042'
        
            !Writes the Grid
            call HDF5WriteData   (Me%ObjEVTPInstHDF, "/Grid", "Bathymetry", "m",           &
                                    Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR052'

            !WriteBasinPoints
            call HDF5WriteData   (Me%ObjEVTPInstHDF, "/Grid", "BasinPoints", "-",          &
                                    Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR062'

            !Flushes All pending HDF5 commands
            call HDF5FlushMemory (Me%ObjEVTPInstHDF, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPHDFOutput - ModuleBasin - ERR072'        
        endif
    end subroutine ConstructEVTPHDFOutput

   !--------------------------------------------------------------------------
   
    !subroutine ConstructEVTPRefHDFOutput
    !
    !    !Arguments-------------------------------------------------------------
    !
    !    !Local-----------------------------------------------------------------
    !    integer                                             :: ILB, IUB, JLB, JUB   
    !    integer                                             :: STAT_CALL
    !    integer                                             :: HDF5_CREATE
    !    
    !    !------------------------------------------------------------------------
    !   
    !   
    !    !Reads file name of the EVTP hdf outupt
    !    call ReadFileName('BASIN_EVTPREFHDF', Me%Files%EVTPRefHDFFile,                        &
    !                       Message = "Basin EVTPREFHDF Output File", STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR010'             
    !
    !   
    !    !Bounds
    !    ILB = Me%WorkSize%ILB
    !    IUB = Me%WorkSize%IUB
    !
    !    JLB = Me%WorkSize%JLB
    !    JUB = Me%WorkSize%JUB           
    !    
    !    Me%EVTPRefOutPut%NextOutPut = 1 
    !    
    !    call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
    !
    !    !Opens HDF File
    !    call ConstructHDF5      (Me%ObjEVTPRefHDF, trim(Me%Files%EVTPRefHDFFile)//"5", HDF5_CREATE, STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR020'
    !
    !  
    !    !Write the Horizontal Grid
    !    call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjEVTPRefHDF, STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR030'
    !
    !
    !    !Sets limits for next write operations
    !    call HDF5SetLimits      (Me%ObjEVTPRefHDF, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR040'
    !    
    !    !Writes the Grid
    !    call HDF5WriteData   (Me%ObjEVTPRefHDF, "/Grid", "Bathymetry", "m",           &
    !                          Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR050'
    !
    !    !WriteBasinPoints
    !    call HDF5WriteData   (Me%ObjEVTPRefHDF, "/Grid", "BasinPoints", "-",          &
    !                          Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR060'
    !
    !    !Flushes All pending HDF5 commands
    !    call HDF5FlushMemory (Me%ObjEVTPRefHDF, STAT = STAT_CALL)
    !    if (STAT_CALL /= SUCCESS_) stop 'ConstructEVTPRefHDFOutput - ModuleBasin - ERR070'       
    !
    !end subroutine ConstructEVTPRefHDFOutput
    !
   !--------------------------------------------------------------------------    
    
    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                             :: ILB, IUB, JLB, JUB    
!        integer                                             :: i, j

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        allocate(Me%ExtUpdate%WaterLevel     (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%ExtUpdate%WaterColumn    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%ExtUpdate%WaterColumnOld (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%WaterColumnRemoved      (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
        allocate(Me%CanopyDrainage          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 

        if (Me%Coupled%Vegetation) then
            allocate(Me%RainCovered             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%CanopyStorageCapacity   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
            allocate(Me%CanopyStorage           (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%CanopyStorageOld        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%CoveredFraction         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%CoveredFractionOld      (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 
            allocate(Me%CropEvapotrans          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                allocate(Me%PotentialTranspiration  (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(Me%PotentialEvaporation    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            endif  
        endif      
        allocate(Me%ThroughFall             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%RainUncovered           (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
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
        
        
        if (Me%Coupled%Evapotranspiration) then
            if (Me%EVTPOutput%Yes) then
                allocate(Me%PartialAccEVTP          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(Me%PartialAccEVTPRef       (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                Me%PartialAccEVTP     = 0.0
                Me%PartialAccEVTPRef  = 0.0
            
                if (Me%Coupled%Vegetation) then
                    allocate(Me%PartialAccEVTPCrop      (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    Me%PartialAccEVTPCrop = 0.0
                endif
                
                if (Me%EvapoTranspirationMethod .eq. SeparateEvapoTranspiration) then
                    allocate(Me%PartialAccEVPot         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%PartialAccEVAct         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%PartialAccETPot         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))            
                    allocate(Me%PartialAccETAct         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))                
                    Me%PartialAccEVPot = 0.0
                    Me%PartialAccEVAct = 0.0
                    Me%PartialAccETPot = 0.0
                    Me%PartialAccETAct = 0.0                
                endif
            endif
            
            if (Me%EVTPOutput2%Yes) then
                allocate(Me%PartialAccEVTP2          (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                allocate(Me%PartialAccEVTPRef2       (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                Me%PartialAccEVTP2     = 0.0
                Me%PartialAccEVTPRef2  = 0.0
            
                if (Me%Coupled%Vegetation) then
                    allocate(Me%PartialAccEVTPCrop2      (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    Me%PartialAccEVTPCrop2 = 0.0
                endif
                
                if (Me%EvapoTranspirationMethod .eq. SeparateEvapoTranspiration) then
                    allocate(Me%PartialAccEVPot2         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%PartialAccEVAct2         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%PartialAccETPot2         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))            
                    allocate(Me%PartialAccETAct2         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))                
                    Me%PartialAccEVPot2 = 0.0
                    Me%PartialAccEVAct2 = 0.0
                    Me%PartialAccETPot2 = 0.0
                    Me%PartialAccETAct2 = 0.0                
                endif
            endif                      
        endif                         
       
        if (Me%Coupled%Snow) allocate(Me%SnowPack (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB)) 

        allocate(Me%TimeSeriesBuffer    (26))
        allocate(Me%TimeSeriesBuffer4   (13))
        allocate(Me%TimeSeriesBuffer2   (1))

        Me%ExtUpdate%WaterLevel               = FillValueReal
        Me%ExtUpdate%WaterColumn              = FillValueReal
        Me%ExtUpdate%WaterColumnOld           = FillValueReal
        
        if (Me%Coupled%Vegetation) then
            Me%RainCovered              = 0.0
            Me%CanopyStorageCapacity    = FillValueReal
            Me%CanopyStorage            = 0.0
            Me%CanopyStorageOld         = 0.0
            Me%CoveredFraction          = 0.0  !This variable was used before had value definition
            Me%CoveredFractionOld       = 0.0
            Me%CropEvapotrans           = 0.0
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                Me%PotentialTranspiration   = 0.0
                Me%PotentialEvaporation     = 0.0
            endif 
        endif       
        Me%ThroughFall              = 0.0
        Me%RainUncovered            = 0.0
        Me%WaterColumnRemoved       = 0.0
        Me%CanopyDrainage           = 0.0
        Me%InfiltrationRate         = FillValueReal
        Me%PrecipRate               = FillValueReal
        Me%ThroughRate              = FillValueReal
        Me%EVTPRate                 = FillValueReal
        !Me%EVTPRate2                = FillValueReal
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
            allocate(Me%SI%Ks%Field             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%MP%Field             (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%ThetaS%Field         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%ThetaI%Field         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%InfRate%Field        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%AccInf%Field         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%ImpFrac%Field        (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            allocate(Me%SI%TimeWithNoWC%Field   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%SI%Ks%Field                  =   FillValueReal
            Me%SI%MP%Field                  =   FillValueReal
            Me%SI%ThetaS%Field              =   FillValueReal
            Me%SI%ThetaI%Field              =   FillValueReal
            Me%SI%InfRate%Field             =   FillValueReal
            Me%SI%AccInf%Field              =   AllmostZero
            Me%SI%ImpFrac%Field             =   FillValueReal
            Me%SI%TimeWithNoWC%Field        =   0.0
         endif
         
        if (Me%DiffuseWaterSource) then
            allocate(Me%DiffuseFlow         (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%DiffuseFlow            = null_real
        endif         
    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructCoupledModules(                 &
#ifdef _ENABLE_CUDA
                                        ObjCudaID       &
#endif _ENABLE_CUDA                                        
                                      )

        !Arguments-------------------------------------------------------------
#ifdef _ENABLE_CUDA
        integer                                     :: ObjCudaID
#endif
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: VariableDT
        integer                                     :: MapID
        logical                                     :: IgnoreWaterColumnOnEvap
        logical                                     :: IsConstant
        real                                        :: ConstantValue
        logical                                     :: VegParticFertilization
        logical                                     :: Pesticide
        !Begin-----------------------------------------------------------------

        !Constructs Atmosphere
        if (Me%Coupled%Atmosphere) then
            call StartAtmosphere        (ModelName          = Me%ModelName,              &
                                         AtmosphereID       = Me%ObjAtmosphere,          &
                                         TimeID             = Me%ObjTime,                &
                                         GridDataID         = Me%ObjGridData,            &
                                         HorizontalGridID   = Me%ObjHorizontalGrid,      &
                                         MappingPoints      = Me%ExtVar%BasinPoints,     &
                                         CheckValues        = Me%VerifyAtmosphereValues, &
                                         STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR010'
        endif

        !Constructs Drainage Network
        if (Me%Coupled%DrainageNetwork) then

            call GetVariableDT(Me%ObjTime, VariableDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR020'
            
            call ConstructDrainageNetwork (ModelName         = Me%ModelName,                     &
                                           DrainageNetworkID = Me%ObjDrainageNetwork,            &
                                           TimeID            = Me%ObjTime,                       &
                                           Size              = Me%Size,                          &
                                           CheckMass         = Me%VerifyGlobalMass,              &
                                           CoupledPMP        = Me%Coupled%PorousMediaProperties, &
                                           CoupledRP         = Me%Coupled%RunoffProperties,      &
                                           STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR030'
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
                                         CheckGlobalMass        = Me%VerifyGlobalMass,      &
                                         ConstructEvaporation   = Me%ConstructEvaporation,  &
                                         ConstructTranspiration = Me%ConstructTranspiration,&
                                         GeometryID             = Me%ObjGeometry,           &
                                         MapID                  = MapID,                    &
                                         STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR040'
            
            call GetIgnoreWaterColumnOnEVAP(Me%ObjPorousMedia, IgnoreWaterColumnOnEvap, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR050' 
            
            if ((.NOT. Me%EvapFromWaterColumn) .AND. (.NOT. IgnoreWaterColumnOnEvap)) then
                write (*,*)
                write (*,*) 'WARNING: Both EVAP_FROM_WATER_COLUMN keyword (Basin) and '
                write (*,*) '         IGNORE_WATER_COLUMN_ON_EVAP keyword (PorousMedia) are set to' 
                write (*,*) '         FALSE.'
                write (*,*)
            endif

               
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
                                                       CheckGlobalMass            = Me%VerifyGlobalMass,          &
#ifdef _ENABLE_CUDA
                                                       CudaID                     = ObjCudaID,                    &
#endif _ENABLE_CUDA
                                                       STAT                       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR060'
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
                                         UsePotLAI          = Me%UsePotLAI,                  &
                                         STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR090'
            
            call GetCanopyStorageType(Me%ObjVegetation, IsConstant, ConstantValue, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR100'            
            
            if (.not. Me%EvapFromCanopy .AND. ((.NOT. IsConstant) .OR. (ConstantValue /= 0.0))) then
                write (*,*)
                write (*,*) 'WARNING: EVAP_FROM_CANOPY keyword (Basin) is set to FALSE and '
                write (*,*) '         "specific leaf storage" value is not constant or have '
                write (*,*) '         a value other than 0.'
                write (*,*)
            endif
            
            !Send info about organic particulate material to runoff properties to allocate flux
            call GetVegetationOptions(Me%ObjVegetation,                                &
                                      VegParticFertilization = VegParticFertilization, &
                                      Pesticide = Pesticide, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR110' 
            
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
            if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR070'
            
           
            !Constructs RunoffProperties 
            if (Me%Coupled%RunoffProperties) then
                call ConstructRunoffProperties        (ObjRunoffPropertiesID      = Me%ObjRunoffProperties,       &
                                                       ComputeTimeID              = Me%ObjTime,                   &
                                                       HorizontalGridID           = Me%ObjHorizontalGrid,         &
                                                       HorizontalMapID            = Me%ObjHorizontalMap,          &
                                                       BasinGeometryID            = Me%ObjBasinGeometry,          &
                                                       RunoffID                   = Me%ObjRunoff,                 &
                                                       GridDataID                 = Me%ObjGridData,               &
 !                                                      GeometryID                 = Me%ObjGeometry,               &
 !                                                      CoupledPMP                 = Me%Coupled%PorousMediaProperties, &
                                                       CoupledDN                  = Me%Coupled%DrainageNetwork,   &
                                                       CheckGlobalMass            = Me%VerifyGlobalMass,          &
                                                       CoupledVegetation          = Me%Coupled%Vegetation,        &
                                                       VegParticFertilization     = VegParticFertilization,       &
                                                       Pesticide                  = Pesticide,                    &
                                                       STAT                       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructCoupledModules - ModuleBasin - ERR080'
            endif
            
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

        !Constructs Output
        if (Me%EVTPOutput%Yes .or. Me%EVTPOutput2%Yes .or. Me%EVTPInstOutput%Yes) then
            call ConstructEVTPHDFOutput()
        endif           
        
    end subroutine ConstructCoupledModules

    !--------------------------------------------------------------------------
    subroutine ConstructPropertyList

        !Local-----------------------------------------------------------------
        type (T_BasinProperty), pointer             :: NewProperty  => null()
        integer                                     :: ClientNumber, nProperties
        integer                                     :: PropID, iProp
        integer                                     :: STAT_CALL, ErrorCount
        logical                                     :: BlockFound
        type (T_BasinProperty), pointer             :: PropertyX    => null()
        logical                                     :: PropAdvDiff, Decay
        integer                                     :: i,j

        !----------------------------------------------------------------------

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)       
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR00'


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                        block_begin, block_end, BlockFound,     &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    call ConstructProperty(NewProperty)

                    ! Add new Property to the Basin List 
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
        
        !Properties inherited from Runoff Properties - needed for mass balance and for vegetation conc
        if (Me%Coupled%RunoffProperties) then
        
            if (Me%VerifyGloBalMass .or. Me%Coupled%Vegetation ) then
                call GetRPnProperties (Me%ObjRunoffProperties, nProperties, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR010'

                do iProp = 1, nProperties

                    call GetRPPropertiesIDByIdx(RunoffPropertiesID = Me%ObjRunoffProperties,            &
                                                 Idx                     = iProp,                       &
                                                 ID                      = PropID,                      &
                                                 PropAdvDiff             = PropAdvDiff,                 &
                                                 Decay                   = Decay,                       &
                                                 STAT                    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR020'

                    !For now only for properties with Advection diffusion 
!                    if (PropAdvDiff) then
                        
                        !Allocates new property
                        allocate (NewProperty, STAT = STAT_CALL)            
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR030'                    
                        
                        NewProperty%Inherited          = .true.
                        NewProperty%ID%IDNumber        = PropID
                        NewProperty%ID%Name            = trim(GetPropertyName(NewProperty%ID%IDNumber))                       
                        NewProperty%Particulate        = Check_Particulate_Property(NewProperty%ID%IDNumber)
                        NewProperty%AdvectionDiffusion = PropAdvDiff
                        NewProperty%Decay              = Decay
                        if (Me%Coupled%Vegetation) then
                            allocate (NewProperty%VegetationConc(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
                            NewProperty%VegetationConc = 0.0
                            allocate (NewProperty%VegetationDrainage(Me%WorkSize%ILB:Me%WorkSize%IUB,                         &
                                                                     Me%WorkSize%JLB:Me%WorkSize%JUB))
                            NewProperty%VegetationDrainage = 0.0
                            allocate (NewProperty%VegetationOldMass(Me%WorkSize%ILB:Me%WorkSize%IUB,                         &
                                                                     Me%WorkSize%JLB:Me%WorkSize%JUB))
                            NewProperty%VegetationOldMass = 0.0                            
                        endif
                    
                        ! Add new Property to the Basin List 
                        call AddProperty(NewProperty)
                        
!                    endif
                enddo
           endif
        endif
                    
        
        !Verifies existence of required properties
        if (Me%Coupled%Evapotranspiration) then

            ErrorCount = 0
            call SearchProperty(PropertyX, RefEvapotrans_        , .true., STAT = STAT_CALL)
            if (Me%Coupled%Evapotranspiration) ErrorCount = ErrorCount + STAT_CALL
            
            if (ErrorCount /= SUCCESS_) stop 'ConstructPropertyList - ModuleBasin - ERR03'
            
            !Define the conversion factor for ET0 based on the UNITS given by the user
            if (trim(adjustl(PropertyX%ID%Units)) .EQ. "mm/8d") then
                Me%ETConversionFactor = 1. / 691200000.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "mm/d") then
                Me%ETConversionFactor = 1. / 86400000.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "mm/h") then 
                Me%ETConversionFactor = 1. / 3600000.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "mm/s") then 
                Me%ETConversionFactor = 1. / 1000.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "cm/d") then 
                Me%ETConversionFactor = 1. / 8640000.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "cm/h") then 
                Me%ETConversionFactor = 1. / 360000.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "cm/s") then 
                Me%ETConversionFactor = 1. / 100.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "m/d") then 
                Me%ETConversionFactor = 1. / 86400.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "m/h") then 
                Me%ETConversionFactor = 1. / 3600.
            elseif (trim(adjustl(PropertyX%ID%Units)) .EQ. "m/s") then 
                Me%ETConversionFactor = 1.
            else 
                write(*,*)  
                write(*,*) 'Unknown unit for reference evapotranspiration property. '
                write(*,*) 'The available units (L/T) are: '
                write(*,*) 'L: mm, cm, m '
                write(*,*) 'T: d (for days), h (for hours), s (for seconds), 8d (for 8 days) '
                stop 'ConstructPropertyList - ModuleBasin - ERR04'
            endif
            
            if (PropertyX%Constant .and. (.not. PropertyX%ConstantInSpace)) then
            
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        
                    if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then
                        
                        PropertyX%Field(i, j) = PropertyX%Field(i, j) * Me%ETConversionFactor
                        !
                        !if (.NOT. PropertyX%ID%SolutionFromFile) then           
                        !    Me%RefEvapotranspirationConstant = PropertyX%Field(i, j)
                        !endif
                    endif
                    
                enddo
                enddo
                
            elseif (PropertyX%Constant) then
                
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        
                    if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then
                        
                        PropertyX%Field(i, j) = Me%RefEvapotranspirationConstant * Me%ETConversionFactor

                    endif
                    
                enddo
                enddo                
                
            endif            
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
        integer                                     :: STAT_CALL, iflag
        
        !----------------------------------------------------------------------

        !Fills Matrix
        allocate (NewProperty%Field (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        call ConstructFillMatrix(PropertyID         = NewProperty%ID,                    &
                                 EnterDataID        = Me%ObjEnterData,                   &
                                 TimeID             = Me%ObjTime,                        &
                                 HorizontalGridID   = Me%ObjHorizontalGrid,              &
                                 ExtractType        = FromBlock,                         &
                                 PointsToFill2D     = Me%ExtVar%BasinPoints,             &
                                 Matrix2D           = NewProperty%Field,                 &
                                 TypeZUV            = TypeZ_,                            &
                                 STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR010'

        call GetIfMatrixRemainsConstant(FillMatrixID    = NewProperty%ID%ObjFillMatrix,     &
                                        RemainsConstant = NewProperty%Constant,             &
                                        STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR020'

        call GetData(NewProperty%ConstantInSpace,                                        &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'CONSTANT_IN_SPACE',                                 &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleBasin',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR040'
        
        if (NewProperty%Constant .and. NewProperty%ConstantInSpace) then            
            call GetDefaultValue (FillMatrixID = NewProperty%ID%ObjFillMatrix,     &
                                  DefaultVAlue = Me%RefEvapotranspirationConstant, &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR030'
        endif
        
        !GR : changed the if condition because no constant default value in the wind would then be possible
        if (.not. NewProperty%ID%SolutionFromFile .and. .not. NewProperty%Constant) then
            call KillFillMatrix (NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleBasin - ERR040'
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

        if (Me%ComputeBasinWaterBalance) then
            read(InitialFile)Me%BWB%AccErrorInVolume
        endif
        
!        read(InitialFile)Me%WaterLevel
        if (Me%Coupled%Vegetation) then
            read(InitialFile)Me%CanopyStorage
            Me%CanopyStorageOld = Me%CanopyStorage
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
        integer                                     :: STAT_CALL !, i, j
        character (Len = StringLength)              :: WarningString
        character (Len = StringLength)              :: LockToWhichModules
        character (Len = StringLength)              :: UnLockToWhichModules
        character (Len = StringLength)              :: MassEvaluationTime
        character (Len = StringLength)              :: OptionsType
        type (T_BasinProperty), pointer             :: Property
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "ModifyBasin")

        STAT_ = UNKNOWN_

        call Ready(ObjBasinID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            call GetComputeCurrentTime  (Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyBasin - ModuleBasin - ERR01'

            !Gets ExternalVars
            LockToWhichModules = 'AllModules'
            OptionsType = 'ModifyBasin'
            call ReadLockExternalVar (LockToWhichModules, OptionsType)

            if (Me%ComputeBasinWaterBalance .and. Me%BWB%StoreInitial) then
                call GetStoredVolumes (InitialInstant)                
            endif
            
            !Verifies Global Mass
            if (Me%VerifyGlobalMass) then
!                call CalculateMass(Me%MB%IniVolumeBasin,                           &
!                                   Me%MB%IniVolumeVegetation,                      &
!                                   Me%MB%IniVolumePorousMedia,                     &
!                                   Me%MB%IniVolumeChannels)
                MassEvaluationTime = "Initial"
                call CalculateMass (MassEvaluationTime)
                
                Property => Me%FirstProperty
                do while (associated(Property)) 
                    
                    Property%MB%TotalRainMass      = 0.0  !total input (covered + uncovered)
                    Property%MB%UncoveredRainMass  = 0.0  !direct rain (uncovered)
                    Property%MB%CoveredRainMass    = 0.0  !on leaf (covered)
                    Property%MB%VegDrainedMass     = 0.0  !leaf leak
                    Property%MB%RunoffInputMass    = 0.0  !leaf leak + direct rain
                    Property%VegTotalStoredMass    = 0.0
                    
                    Property => Property%Next
                enddo               
            endif

            !Atmospheric Processes 
            if (Me%Coupled%Atmosphere) then

                call AtmosphereProcesses 

                !Actualizes the WaterColumn
                call ActualizeWaterColumn  
                
                if (Me%Coupled%RunoffProperties) then
                    WarningString = 'AtmosphereProcesses'
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
                call ActualizeWaterColumn                 
                
                
                if (Me%Coupled%PorousMediaProperties) then

                    call PorousMediaPropertiesProcesses     
                    WarningString = 'PorousMediaProcesses'
                    call ActualizeWaterColumnConc(WarningString) 
                    
                endif           
               
            else
                !do nothing - water column and concentration already updated in atmosphere processes
                
!                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!                    if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then
!                    
!                        !Potencial Infiltration Column = ThroughFall
!                        Me%ExtUpdate%WaterLevel (i, j) = Me%ExtUpdate%WaterLevel (i, j) + Me%Throughfall(i, j)
!                    end if
!                end do
!                end do
!
!                !Actualizes the WaterColumn
!                call ActualizeWaterColumn (WarningString) 
!
!                if (Me%Coupled%RunoffProperties) then
!                    WarningString = 'PorousMediaProcesses 2'
!                    call ActualizeWaterColumnConc(WarningString) 
!                endif
                
            endif
            

            !Simplified Infiltration / Evapotranspiration model
            if (Me%Coupled%SimpleInfiltration) then
                call SimpleInfiltration

                !Actualizes the WaterColumn
                call ActualizeWaterColumn  
                
                if (Me%Coupled%RunoffProperties) then
                    WarningString = 'SimpleInfiltration'
                    call ActualizeWaterColumnConc(WarningString) 
                endif                
            endif
            
            !Overland Flow
            if (Me%Coupled%RunOff) then
                call OverLandProcesses
                
                !No actualization - Runoff is the water column handler
                !Actualizes the WaterColumn
!                WarningString = 'OverLandProcesses'
!                call ActualizeWaterColumn (WarningString) 

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
            
            !compute vegetation mass on top of foils
            if (Me%Coupled%Vegetation) then
                call CalculateVegTotalStoredMass
            endif
            
            !Verifies Global Mass
            if (Me%VerifyGlobalMass) then
!                call CalculateMass(Me%MB%FinalVolumeBasin,                           &
!                                   Me%MB%FinalVolumeVegetation,                      &
!                                   Me%MB%FinalVolumePorousMedia,                     &
!                                   Me%MB%FinalVolumeChannels)
                MassEvaluationTime = "Final"
                call CalculateMass (MassEvaluationTime)
                !Verifies Mass and makes global balance...
                call GlobalMassBalance
            endif
            
            if (Me%ComputeBasinWaterBalance) then
                call GetStoredVolumes (FinalInstant)
                call ComputeBasinWaterBalance
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

            !EVTP HDF5 Output
            if (Me%EVTPOutput%Yes .or. Me%EVTPOutput2%Yes .or. Me%EVTPInstOutput%Yes) then
                call EVTPHDFOutPut       
            endif                    
            
            !UnGets ExternalVars
            UnLockToWhichModules = 'AllModules'
            OptionsType = 'ModifyBasin'
            call ReadUnLockExternalVar (UnLockToWhichModules, OptionsType)
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
        real, dimension(:, :), pointer              :: IrrigationFlux
        character (Len = StringLength)              :: LockToWhichModules
        character (Len = StringLength)              :: UnLockToWhichModules
        character (Len = StringLength)              :: WarningString
        logical                                     :: IrrigationExists = .false.
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "AtmosphereProcesses")

        !Updates Rainfall
        call ModifyAtmosphere       (Me%ObjAtmosphere, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR010'

        !Gets Rainfall [m3/s]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, PrecipitationFlux, ID = Precipitation_, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR020'

        !Gets Irrigation [m3/s]
        call GetAtmosphereProperty  (Me%ObjAtmosphere, IrrigationFlux, ID = Irrigation_, STAT = STAT_CALL, ShowWarning = .false.)
        if (STAT_CALL == SUCCESS_) then
            IrrigationExists = .true.
        elseif (STAT_CALL == NOT_FOUND_ERR_) then
            IrrigationExists = .false.
        else 
            stop 'AtmosphereProcesses - ModuleBasin - ERR030'
        endif
        
        if (Me%Coupled%Vegetation) then 
            
            !Read Vegetation Properties
            LockToWhichModules = 'Atmosphere'
            call ReadLockExternalVar (LockToWhichModules)
        
        endif

        !Divides Precipitation into Snow, Rain, Throughfall and Canopy storage
        if (IrrigationExists) then
            call DividePrecipitation(PrecipitationFlux, IrrigationFlux)
        else
            call DividePrecipitation(PrecipitationFlux)
        endif
        
        !Mass Balance in rain (covered and uncovered) and drained from veg (if active)
        if (Me%Coupled%RunoffProperties) then
            
            if (Me%VerifyGlobalMass) then
                call DividePrecipitationMassFluxes
            endif
        
            !Actualize Vegetation conc because of mix with rain and canopy drainage
            if (Me%Coupled%Vegetation) then
                WarningString = "WaterMix"
                call UpdateVegConcentration(WarningString)
            endif
        endif


        if (Me%Coupled%Evapotranspiration) then

                call CalcPotEvapoTranspiration
                
                !Actualize Vegetation conc because of evaporation
                if (Me%Coupled%Vegetation .and. Me%Coupled%RunoffProperties) then
                    WarningString = "Evaporation"
                    call UpdateVegConcentration(WarningString)
                endif
        endif

        if (Me%Coupled%Vegetation) then 
            
            !UnlockVegetation Properties so they can be run when vegetation processes are called
            UnLockToWhichModules = 'Atmosphere'
            call ReadUnLockExternalVar (UnLockToWhichModules)

        endif

        call UnGetAtmosphere    (Me%ObjAtmosphere, PrecipitationFlux, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR040'


        if (IrrigationExists) then
            call UnGetAtmosphere    (Me%ObjAtmosphere, IrrigationFlux, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'AtmosphereProcesses - ModuleBasin - ERR050'
        endif
        
        if (MonitorPerformance) call StopWatch ("ModuleBasin", "AtmosphereProcesses")
    
    end subroutine AtmosphereProcesses

    !--------------------------------------------------------------------------

    subroutine DividePrecipitation (PrecipitationFlux, IrrigationFlux)
!    subroutine DividePrecipitation (PrecipitationFlux)
    
        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer                       :: PrecipitationFlux
        real, dimension(:, :), pointer, optional, intent(IN) :: IrrigationFlux

        !Local-----------------------------------------------------------------     
        integer                                     :: i, j, STAT_CALL
        real, dimension(:,:), pointer               :: AirTemperature
        real                                        :: GrossPrecipitation
        real                                        :: SnowInput
!        real                                        :: CanopyDrainage
        real                                        :: CurrentFlux, Rand, SecondsPassed
        real, dimension(6), target                  :: AuxTime
        integer, save                               :: LastHour = -99
        logical                                     :: ChangeRain
        real(8)                                     :: NewVolumeOnLeafs, OldVolumeOnLeafs
        real                                        :: AreaFraction
        real, dimension(:, :), pointer              :: Irri
        logical                                     :: IsPresent
        !Begin-----------------------------------------------------------------
        
        NewVolumeOnLeafs = 0.0
        OldVolumeOnLeafs = 0.0
        
        if (present(IrrigationFlux)) then
            Irri => IrrigationFlux
            IsPresent = .true.
        else
            Irri => null()
            IsPresent = .false.    
        endif
        
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
        

        Me%MB%TotalRainIn     = 0.0  !total input
        Me%MB%RainDirect      = 0.0  !on uncovered area
        Me%MB%RainInVeg       = 0.0  !on covered area (leafs)
        Me%MB%DrainageFromVeg = 0.0  !leaf leak
        Me%MB%RainRunoff      = 0.0  !arriving to runoff (leaf drainage + direct)
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
                          
                if (Me%ComputeBasinWaterBalance) then
                    
                    !    m3     =     m3      +     m3/s    *      s
                    Me%BWB%Rain = Me%BWB%Rain + CurrentFlux * Me%CurrentDT 
                    
                    if (IsPresent) then
                        !      m3         =         m3        +    m3/s    *      s
                        Me%BWB%Irrigation = Me%BWB%Irrigation + Irri(i, j) * Me%CurrentDT
                    endif
                
                endif                
                
                if (IsPresent) then
                    CurrentFlux = CurrentFlux + Irri(i, j)
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
                    Me%MB%TotalRainIn = Me%MB%TotalRainIn + GrossPrecipitation * Me%ExtVar%GridCellArea(i, j)
                endif
                ! For now uncovered rain is total. it will be changed ir there are leafs
                Me%RainUncovered(i, j) = GrossPrecipitation
                
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
                    
                    Me%CoveredFractionOld(i,j) = Me%CoveredFraction(i, j)
                    
                    !Volume on Leaf (before leaf grow)
                    !m3 = m * m2plant
                    OldVolumeOnLeafs = Me%CanopyStorage(i, j) * Me%ExtVar%GridCellArea(i, j) * Me%CoveredFraction(i, j)
                    
                    !Calculates Covered Fraction
                    Me%CoveredFraction(i, j) = 1.0 - exp(-Me%ExtinctionCoef * Me%ExtVar%LeafAreaIndex(i, j))
                    
                    !Recalculates Canopy Storage Capacity (for the case that LAI is variable in time)
                    !mH20 = m3H20/m2leaf * m2leaf/m2soil
                    Me%CanopyStorageCapacity(i, j) = Me%ExtVar%SpecificLeafStorage(i, j) * Me%ExtVar%LeafAreaIndex(i, j)
                    
                    !Precipitations heights - heights are the same as gross (area different and volume different).
                    !height * (areacov + areauncov) = same volume as gross volume
                    !m         
                    Me%RainUncovered(i, j)    = GrossPrecipitation 
                    Me%RainCovered(i, j)      = GrossPrecipitation 
                    
                    !Throughfall on uncovered area - has to be transformed in height in total area
                    !m = (m * m2plant)/m2cell
                    Me%ThroughFall(i, j)      = Me%RainUncovered(i, j) * (1.0 - Me%CoveredFraction(i, j))               &
                                               * Me%ExtVar%GridCellArea(i, j) / Me%ExtVar%GridCellArea(i, j)
                    
                    !New volumes on leafs is the old volume plus the gross rain on the covered fraction of the cell
                    !m3 = m3 + (m*m2plant)
                    NewVolumeOnLeafs = OldVolumeOnLeafs + Me%RainCovered(i, j) * Me%ExtVar%GridCellArea(i, j)           &
                                       * Me%CoveredFraction(i, j)
                    
                    Me%CanopyStorageOld(i,j) = Me%CanopyStorage(i, j)
                    
                    !If LAI exists plant can drip from leafs
                    if (Me%CoveredFraction(i, j) > 0.0) then
                        
                        !m = m3 / m2plant
                        Me%CanopyStorage(i, j) = NewVolumeOnLeafs / (Me%CoveredFraction(i, j) * Me%ExtVar%GridCellArea(i, j))
                        
                        !Calculates CanopyDrainage so that CanopyStorage as maximum is full
                        if (Me%CanopyStorage(i, j) > Me%CanopyStorageCapacity(i, j)) then
                            !Canopy drainage is height in terms of cell area and not in terms of plant area
                            !m = (m * m2plant) / m2cell
                            Me%CanopyDrainage(i,j) = ((Me%CanopyStorage(i, j) - Me%CanopyStorageCapacity(i, j))               &
                                                      * (Me%CoveredFraction(i, j) * Me%ExtVar%GridCellArea(i, j)))            &
                                                      / Me%ExtVar%GridCellArea(i, j)
                            Me%CanopyStorage(i, j) = Me%CanopyStorageCapacity(i, j)
                        else
                            Me%CanopyDrainage(i,j) = 0.0
                        endif
                        
                    else !plant may drip if this was the instant that lost the leafs (old volume not zero)
                        
                        !m = m3 / m2cell - water equally distributed because leafs disappeared
                        Me%CanopyDrainage(i,j)     = NewVolumeOnLeafs / Me%ExtVar%GridCellArea(i, j)
!                        Me%ThroughFall  (i, j)     = Me%ThroughFall(i, j) + Me%CanopyDrainage(i,j)
                        Me%CanopyStorage(i, j)     = 0.0
                    endif

                    !Adds Canopy drainage to Throughfall - both have the same area associated (cell area)
                    !m
                    Me%ThroughFall(i, j)   = Me%ThroughFall(i, j) + Me%CanopyDrainage(i,j) 
                    
                    !Integrates Total Input Volume
                    if (Me%VerifyGlobalMass) then
                        !m3 = m3 + (m * m2plant) - rain covered is height in plant area
                        Me%MB%RainInVeg     = Me%MB%RainInVeg + (Me%RainCovered(i, j) * Me%ExtVar%GridCellArea(i, j)       &
                                            * Me%CoveredFraction(i, j))
                        !m3 = m3 + (m * m2cell) - Canopy drainage is height in cell area                   
                        Me%MB%DrainageFromVeg = Me%MB%DrainageFromVeg + (Me%CanopyDrainage(i,j) * Me%ExtVar%GridCellArea(i, j))
                    endif            
                else
                    
                    Me%ThroughFall(i, j)    = GrossPrecipitation
                    Me%CanopyDrainage(i,j)  = 0.0
                endif
                
                !Put rain on water column - it will facilitate the structure of the property mixture
                Me%ExtUpdate%WaterLevel(i,j) = Me%ExtUpdate%WaterLevel(i,j) + Me%ThroughFall(i, j) 
                    
                !mm/ hour                   m                       s         mm/m     s/h
                Me%ThroughRate(i, j) = Me%ThroughFall(i, j) / Me%CurrentDT * 1000.0 * 3600.0

                if (Me%VerifyGlobalMass) then
                    if (Me%Coupled%Vegetation) then
                        AreaFraction = 1 - Me%CoveredFraction(i,j)
                    else
                        AreaFraction = 1.0
                    endif
                    !Rain in uncovered area
                    Me%MB%RainDirect = Me%MB%RainDirect + Me%RainUncovered(i,j) * AreaFraction   &
                                          * Me%ExtVar%GridCellArea(i, j)
                    !Total rain arriving at the runoff -> uncovered + drainage
                    Me%MB%RainRunoff = Me%MB%RainRunoff + Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j)
                endif                
                    
                if (Me%Coupled%Vegetation .and. Me%ComputeBasinWaterBalance) then
                    !m3 = m3
                    Me%BWB%IniStoredInLeaves = OldVolumeOnLeafs
                    !The amount of water on leaves can change if there is evaporation from leaves.
                    !In this case, the value in StoredInLeaves must be updated.
                    !The new value will be the actual value calculated here minus the total evaporated from leaves.
                    ! m3 = m3 + m * m2
                    Me%BWB%FinStoredInLeaves = Me%BWB%FinStoredInLeaves + Me%CanopyStorage(i, j) * Me%ExtVar%GridCellArea(i, j)
                endif
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
        real                                        :: Kc
        real,    parameter                          :: LatentHeatOfVaporization = 2.5e6         ![J/kg]
        real,    parameter                          :: ReferenceDensity         = 1000.         ![kg/m3]
        real                                        :: LatentHeat_
        real(8)                                     :: EvaporationRate, dH
        real(8)                                     :: DTFactor
        type(T_BasinProperty), pointer              :: RefEvapotrans
        logical                                     :: EvaporateFromCanopy      = .false.
        logical                                     :: EvaporateFromWaterColumn = .false.
        logical                                     :: CalcET0                  = .false.
        logical                                     :: CalcET0OnlyForNoData     = .false.
        real(8)                                     :: Evaporation
        real, dimension(:,:), pointer               :: EvaporationMatrix => null()
!        real                                        :: LAI
!        real                                        :: Year, Month, Day, hour, minute, second
!        integer                                     :: days
        !Begin-----------------------------------------------------------------
      
        call SearchProperty(RefEvapotrans, RefEvapotrans_, .true., STAT = STAT_CALL)        

        if (RefEvapotrans%ID%SolutionFromFile) then
            call ModifyFillMatrix (FillMatrixID   = RefEvapotrans%ID%ObjFillMatrix,          &
                                   Matrix2D       = RefEvapotrans%Field,                     &
                                   PointsToFill2D = Me%ExtVar%BasinPoints,                   &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalcPotEvapoTranspiration - ModuleBasin - ERR01'
            
        elseif (RefEvapotrans%Constant .and. RefEvapotrans%ConstantInSpace) then
        
            RefEvapotrans%Field = Me%RefEvapotranspirationConstant    
        
        endif 
        
        DTFactor = 1000 * Me%CurrentDT
        
        !Calculates evaporation from canopy / Watercolumn on the ground
        Me%MB%EvapFromVegetation = 0.0
        Me%MB%EvapFromGround     = 0.0
                
        !CalcET0             = .NOT. (RefEvapotrans%ID%SolutionFromFile .OR. RefEvapotrans%Constant)
        CalcET0             = (.not. RefEvapotrans%ID%SolutionFromFile) .and. (.not. RefEvapotrans%Constant) .and. (.not. RefEvapotrans%ConstantInSpace)
        
        if (Me%UseRefEVTPIfNeeded) then
            CalcET0 = .true.
            CalcET0OnlyForNoData = .true.
        else
            CalcET0OnlyForNoData = .false.
        endif
                
        if ((CalcET0) .OR. (Me%EvapMethod .EQ. LatentHeatMethod)) then
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
        endif
                
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) .EQ. WaterPoint) then
               
                if (CalcET0 .and. ((.not. CalcET0OnlyForNoData ) .or. (RefEvapotrans%Field(i,j) < 0.0)))then
            
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
                
                    RefEvapotrans%Field(i, j)  = max(RefEvapotrans%Field(i, j) / 3600000., 0.0)
                
                elseif (RefEvapotrans%ID%SolutionFromFile) then
                
                    !m/s - Porous media consistency. If constant, already converted in the construction of the property
                    RefEvapotrans%Field(i, j)  = max(RefEvapotrans%Field(i, j) * Me%ETConversionFactor, 0.0)                    
                
                endif        
                    
                if (Me%EVTPOutput%yes) then
                    !mm                           = mm                         +  m/s                       * mm/m * s
                    Me%PartialAccEVTPRef   (i, j) = Me%PartialAccEVTPRef (i,j) + (RefEvapotrans%Field(i, j) * DTFactor)                 
                endif
                
                if (Me%EVTPOutput2%yes) then
                    !mm                           = mm                         +  m/s                       * mm/m * s
                    Me%PartialAccEVTPRef2   (i, j) = Me%PartialAccEVTPRef2 (i,j) + (RefEvapotrans%Field(i, j) * DTFactor)                 
                endif                
                
!                call ExtractDate(Me%CurrentTime, Year, Month, Day, hour, minute, second) 
!                call DateToGregorianDay(Me%CurrentTime, days)
!                                                                     
                if (Me%Coupled%Vegetation) then
                    
!                    if (Me%ExtVar%CropCoefficient(i, j) < 0.0) then
!                        write (*,*) 'Me%ExtVar%CropCoefficient(', i, ',', j, ') NEGATIVO em ', days, Year, &
!                                     Month, Day, hour, minute, second
!                        stop
!                    endif
                 
                    if (Me%UsePotLAI) then
                        Kc = AdjustCropCoefficient(i, j)
                    else
                        Kc = Me%ExtVar%CropCoefficient(i, j)
                    endif
                    
!                    if (Kc < 0.0) then
!                        write (*,*) 'Kc(', i, ',', j, ') NEGATIVO em ', days, Year, Month, Day, hour, minute, second
!                        stop
!                    endif
                    
                    !m/s
                    Me%CropEvapotrans(i, j) = RefEvapotrans%Field(i, j) * Kc
!                    if (Me%CropEvapotrans(i, j) < 0.0) then
!                        write (*,*) 'I/J: ', i, '/', j, 'CropEvapotrans NEGATIVO em ',Year, Month, Day, hour, minute, second
!                        stop
!                    endif

                    if (Me%EVTPOutput%yes) then
                        !mm                           = mm                         +  m/s                       * mm/m * s
                        Me%PartialAccEVTPCrop   (i, j) = Me%PartialAccEVTPCrop (i,j) + (Me%CropEvapotrans(i, j) * DTFactor)                 
                    endif
                
                    if (Me%EVTPOutput2%yes) then
                        !mm                           = mm                         +  m/s                       * mm/m * s
                        Me%PartialAccEVTPCrop2   (i, j) = Me%PartialAccEVTPCrop2 (i,j) + (Me%CropEvapotrans(i, j) * DTFactor)                 
                    endif  

                    if (Me%EvapoTranspirationMethod .EQ. SeparateEvapoTranspiration) then
                        
!                        if (Me%ExtVar%LeafAreaIndex(i, j) < 0.) then
!!                            write (*,*) 'LAI NEGATIVO em ', days, Year, Month, Day, hour, minute, second                            
!                            LAI = 0.
!                        else                            
!                            LAI = Me%ExtVar%LeafAreaIndex(i, j)                            
!                        endif
                        
                        !m/s
                        Me%PotentialTranspiration(i, j) = Me%CropEvapotrans(i, j) * &
                                                          (1.0 - exp(-0.463 * Me%ExtVar%LeafAreaIndex(i, j)))
                        if (Me%EVTPOutput%yes) then
                            !mm                           = mm                         +  m/s                       * mm/m * s
                            Me%PartialAccETPot   (i, j) = Me%PartialAccETPot (i,j) + (Me%PotentialTranspiration(i, j) * DTFactor)                 
                        endif
                
                        if (Me%EVTPOutput2%yes) then
                            !mm                           = mm                         +  m/s                       * mm/m * s
                            Me%PartialAccETPot2   (i, j) = Me%PartialAccETPot2 (i,j) + (Me%PotentialTranspiration(i, j) * DTFactor)                 
                        endif                          
                        
                        !m/s
                        Me%PotentialEvaporation  (i, j) = Me%CropEvapotrans(i, j) - Me%PotentialTranspiration(i, j)
                        if (Me%EVTPOutput%yes) then
                            !mm                           = mm                         +  m/s                       * mm/m * s
                            Me%PartialAccEVPot   (i, j) = Me%PartialAccEVPot (i,j) + (Me%PotentialEvaporation(i, j) * DTFactor)                 
                        endif
                
                        if (Me%EVTPOutput2%yes) then
                            !mm                           = mm                         +  m/s                       * mm/m * s
                            Me%PartialAccEVPot2   (i, j) = Me%PartialAccEVPot2 (i,j) + (Me%PotentialEvaporation(i, j) * DTFactor)                 
                        endif  
                    endif
                
                endif

            endif
            
        enddo
        enddo
        
        EvaporateFromCanopy      = Me%Coupled%Vegetation .AND. Me%EvapFromCanopy 
        EvaporateFromWaterColumn = Me%Coupled%Vegetation .AND. Me%EvapFromWaterColumn
        
        if (Me%Coupled%Vegetation) then
            !Pointer so that mass balance right even if there is no evaporation
            call SetMatrixValue (Me%CanopyStorageOld, Me%Size, Me%CanopyStorage)
        endif
        
        if (Me%EvapMethod .EQ. LatentHeatMethod) then
        
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExtVar%BasinPoints(i, j) .EQ. WaterPoint) then

                    LatentHeat_ = -1.0 * LatentHeat (ReferenceDensity, AirTemperature(i, j), AirTemperature(i, j), &
                                                     RelativeHumidity(i, j), WindModulus(i,j))                
                    ![m/s]          = [J/m2/s] / [J/kg] / [kg/m3] 
                    EvaporationRate = LatentHeat_ / LatentHeatOfVaporization / ReferenceDensity 
        
                    if (EvaporateFromCanopy) then
                                       
                        !dH
                        dH = min(dble(EvaporationRate * Me%CurrentDT), Me%CanopyStorage(i, j))
                    
                        !Accumulated EVAP from Canopy
                        Me%AccEVPCanopy(i, j) = Me%AccEVPCanopy(i, j) + dH
                        
                        !New Canopy Storage
                        Me%CanopyStorage(i, j) = Me%CanopyStorage(i, j) - dH
                    
                        !Sistem loss
                        if (Me%VerifyGlobalMass) then
                            Me%MB%EvapFromVegetation = Me%MB%EvapFromVegetation +                                   &
                                                       dH * Me%ExtVar%GridCellArea(i, j) * Me%CoveredFraction(i, j)
                        endif
                        
                        if (Me%ComputeBasinWaterBalance) then
                            ! m3 = m3 + m * m2 * [%]
                            Me%BWB%EvapFromLeaves = Me%BWB%EvapFromLeaves + dH * Me%ExtVar%GridCellArea(i, j) * Me%CoveredFraction(i, j)
                            
                            !Actualization of the water stored on leaves
                            ! m3 = m3 - m3
                            Me%BWB%FinStoredInLeaves = Me%BWB%FinStoredInLeaves - Me%BWB%EvapFromLeaves
                        endif
                        
                    endif 
                    
                    !Also EVAP from watercolumn - important for 1D cases to avoid accumulation of Water on the surface
                    if (EvaporateFromWaterColumn) then
                    
                        !dH
                        dH = min(dble(EvaporationRate * Me%CurrentDT), Me%ExtUpdate%WaterLevel (i, j) - Me%ExtVar%Topography(i, j))

                        Me%AccEVTP(i, j) = Me%AccEVTP(i, j) + dH
                        
                        !New Water Level
                        Me%ExtUpdate%WaterLevel(i, j) = Me%ExtUpdate%WaterLevel(i, j) - dH

                        !Evaporation from Ground
                        if (Me%VerifyGlobalMass) then
                            Me%MB%EvapFromGround = Me%MB%EvapFromGround + dH * Me%ExtVar%GridCellArea(i, j)
                        endif
                        
                        Me%WaterColumnEvaporated(i, j) = dH   
                                             
                    else                    
                    
                        Me%WaterColumnEvaporated(i, j) = 0.0   
                                             
                    endif 
                    
                    if (Me%ComputeBasinWaterBalance) then
                        ! m3 = m3 + m * m2
                        Me%BWB%EvapFromSurface = Me%BWB%EvapFromSurface + Me%WaterColumnEvaporated(i, j) * Me%ExtVar%GridCellArea(i, j)
                    endif                    
                        
                endif
                
            enddo
            enddo
        
        else if (Me%EvapMethod .EQ. ET0Method) then
               
               
            if (Me%Coupled%Vegetation) then                   
                if (Me%EvapoTranspirationMethod .EQ. SeparateEvapoTranspiration) then                        
                    !m/s
                    EvaporationMatrix => Me%PotentialEvaporation
                else                        
                    !m/s
                    EvaporationMatrix => Me%CropEvapotrans
                endif                
            else 
                !m/s               
                EvaporationMatrix => RefEvapotrans%Field
            endif              
               
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                !m
                Evaporation = EvaporationMatrix(i, j) * Me%CurrentDT

                if (Me%ExtVar%BasinPoints(i, j) .EQ. WaterPoint) then

                    if (EvaporateFromCanopy) then
                                       
                        !dH
                        dH = min(Evaporation, Me%CanopyStorage(i, j))
                    
                        !Accumulated EVAP from Canopy
                        Me%AccEVPCanopy(i, j) = Me%AccEVPCanopy(i, j) + dH
                    
                        !New Canopy Storage
                        Me%CanopyStorage(i, j) = Me%CanopyStorage(i, j) - dH
                    
                        !Sistem loss
                        if (Me%VerifyGlobalMass) then
                            Me%MB%EvapFromVegetation = Me%MB%EvapFromVegetation +                                   &
                                                       dH * Me%ExtVar%GridCellArea(i, j) * Me%CoveredFraction(i, j)
                        endif
                        
                        Evaporation = Evaporation - dH
                        
                        if (Me%ComputeBasinWaterBalance) then
                            ! m3 = m3 + m * m2 * [%]
                            Me%BWB%EvapFromLeaves = Me%BWB%EvapFromLeaves + dH * Me%ExtVar%GridCellArea(i, j) * Me%CoveredFraction(i, j)
                            
                            !Actualization of the water stored on leaves
                            ! m3 = m3 - m3
                            Me%BWB%FinStoredInLeaves = Me%BWB%FinStoredInLeaves - Me%BWB%EvapFromLeaves
                        endif                        
                        
                    endif 
                                      
                    !Also EVAP from watercolumn - important for 1D cases to avoid accumulation of Water on the surface
                    if (EvaporateFromWaterColumn) then
                    
                        !dH
                        dH = min(Evaporation, Me%ExtUpdate%WaterLevel (i, j) - Me%ExtVar%Topography(i, j))

                        Me%AccEVTP(i, j) = Me%AccEVTP(i, j) + dH
                        
                        !New Water Level
                        Me%ExtUpdate%WaterLevel(i, j) = Me%ExtUpdate%WaterLevel(i, j) - dH

                        !Evaporation from Ground
                        if (Me%VerifyGlobalMass) then
                            Me%MB%EvapFromGround = Me%MB%EvapFromGround + dH * Me%ExtVar%GridCellArea(i, j)
                        endif
                        
                        Me%WaterColumnEvaporated(i, j) = dH   
                        Evaporation = Evaporation - dH
                                             
                    else                    
                    
                        Me%WaterColumnEvaporated(i, j) = 0.0   
                                             
                    endif 

                    if (Me%ComputeBasinWaterBalance) then
                        ! m3 = m3 + m * m2
                        Me%BWB%EvapFromSurface = Me%BWB%EvapFromSurface + Me%WaterColumnEvaporated(i, j) * Me%ExtVar%GridCellArea(i, j)
                    endif
                    
                    !m/s              
                    EvaporationMatrix(i, j) = Evaporation / Me%CurrentDT
                                           
                endif
                                
            enddo
            enddo        
        
        endif
                
        if ((CalcET0) .OR. (Me%EvapMethod .EQ. LatentHeatMethod)) then                                
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
        endif

    end subroutine CalcPotEvapoTranspiration

    !------------------------------------------------------------------------------
   
    real function AdjustCropCoefficient (I, J)
    
        !Arguments-----------------------------------------------------------------
        integer :: I, J
    
        !Local---------------------------------------------------------------------
        real :: LAI, PotLAI, Kc
        
        !Begin---------------------------------------------------------------------
        LAI = max(Me%ExtVar%LeafAreaIndex(I, J), 0.0)
        PotLAI = Me%ExtVar%PotLeafAreaIndex(I, J)
        Kc = Me%ExtVar%CropCoefficient(I, J)
        
        if (LAI .EQ. 0) then
            if (Me%UseKcMin) then
                AdjustCropCoefficient = Me%KcMin                
            else
                AdjustCropCoefficient = 0.0
            endif
        else
        
            !From 'Necessidades de Água e Métodos de Rega', Luis Santos Pereira, 2004
            !Publicações Europa-América, pg 85-86
            if ((Me%UseDefaultKcWhenLAI0) .AND. (LAI .EQ. 0.0)) then
                if (Me%UseKcMin) then
                    AdjustCropCoefficient = max(Me%DefaultKcWhenLAI0, Me%KcMin)
                else
                    AdjustCropCoefficient = Me%DefaultKcWhenLAI0
                endif
            else if ((PotLAI > LAI) .AND. (PotLAI > 0.0)) then           
                if (Me%UseKcMin) then
                    AdjustCropCoefficient = max(Kc - (1 - (LAI / PotLAI)**0.5), Me%KcMin)                
                else
                    AdjustCropCoefficient = max(Kc - (1 - (LAI / PotLAI)**0.5), 0.0)
                endif
            else
                if (Me%UseKcMin) then
                    AdjustCropCoefficient = max(Kc, Me%KcMin)
                else
                    AdjustCropCoefficient = max(Kc, 0.0)
                endif
            endif
            
        endif
        
!        if (AdjustCropCoefficient < 0.0) then
!            write (*,*) 'I: ', I, ', J: ', J
!            write (*,*) 'LAI: ', LAI
!            write (*,*) 'PotLAI: ', PotLAI
!            write (*,*) 'CropCoefficient: ', Me%ExtVar%CropCoefficient(I, J) 
!            if (Me%UseKcMin) then
!                write (*,*) 'Uses Kc Min. Value: ', Me%KcMin
!            endif
!            write (*,*) 'Kc: ', kc
!            stop
!        endif
        
        !--------------------------------------------------------------------------
    
    end function AdjustCropCoefficient
       
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
    
    subroutine UpdateVegConcentration(WarningString)
    
        !Arguments-------------------------------------------------------------
        character (Len = *), intent(in)             :: WarningString

        !Local-----------------------------------------------------------------
        real(8)                                     :: OldVolumeOnLeafs
        type (T_BasinProperty), pointer             :: Property
        real, dimension(:,:), pointer               :: AtmConcentration
        real(8)                                     :: RainMassToVeg
        real(8)                                     :: DrainageMassFromVeg, VegetationNewMass
        real(8)                                     :: VegetationNewVolume, RainVolume
        real(8)                                     :: DrainageVolume
        integer                                     :: STAT_CALL, i,j
        integer                                     :: NumberOfPesticides, Pest, PesticideIDNumber
        real, dimension(:,:), pointer               :: PesticideVegetation
        real(8)                                     :: MassAdded, MassSink
        real                                        :: VegDT, DT
        real                                        :: DecayRate
        !Begin-----------------------------------------------------------------
                                  
        !Mass Balance to vegetation water on leafs
                
        !Mix after routine DividePrecipitation
        if (WarningString == "WaterMix") then      

            Property => Me%FirstProperty
            do while (associated(Property))
            
                !only for runoff properties and exclude for instance evapotranspiration
                if (Property%Inherited) then
                    
                    if (AtmospherePropertyExists (Me%ObjAtmosphere, Property%ID%IDNumber)) then
            
                        call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
                                                   Scalar             = AtmConcentration,       &
                                                   ID                 = Property%ID%IDNumber,   &
                                                   STAT               = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'UpdateVegConcentration - ModuleBasin - ERR10'  
                    else
                        allocate (AtmConcentration(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
                        AtmConcentration = 0.0
                    endif                
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then
                            
                            !At this point old volume has to be the one from previous time step - covered fraction old
                            !and canopystorage before the mix with rain and removal from drainage
                            OldVolumeOnLeafs = Me%CoveredFractionOld(i, j) * Me%CanopyStorageOld(i, j)                        &
                                               * Me%ExtVar%GridCellArea(i, j)
                            
                            !Mass from Rain
                            !m3 = m * m2plant
                            RainVolume          = Me%RainCovered(i,j) * Me%CoveredFraction(i, j)                              &
                                                   * Me%ExtVar%GridCellArea(i, j) 
                            !g = (m * m2plant) * g/m3
                            RainMassToVeg       = RainVolume * AtmConcentration(i,j)
                                                   
                            !m3
                            VegetationNewVolume = OldVolumeOnLeafs + RainVolume   
                            !g
                            VegetationNewMass   = Property%VegetationOldMass(i,j) + RainMassToVeg
                            
                            !Update veg conc with rain because the drainage flux needs a new conc
                            if (VegetationNewVolume .gt. 0.0) then
                                !g/m3 = g / (m * m2)
                                Property%VegetationConc(i,j) = VegetationNewMass / VegetationNewVolume
                            else
                                Property%VegetationConc(i,j) = 0.0
                            endif                                                   
                            
                            !Mass drained to soil - with the new veg conc.
                            !m3 = m * m2plant
                            DrainageVolume      = Me%CanopyDrainage(i,j) * Me%ExtVar%GridCellArea(i, j)
                            
                            !g = (m * m2cel) * g/m3 - Canopy drainage is height in cell area
                            DrainageMassFromVeg = DrainageVolume * Property%VegetationConc(i,j)
                            
                            Property%VegetationDrainage(i,j) = DrainageMassFromVeg

                            if (Me%VerifyGlobalMass) then
                                !output for vegetation leafs
                                !kg = kg + ((m*m2cell * g/m3 * 1e-3 kg/g) - Canopy drainage is height in cell area
                                Property%MB%VegDrainedMass   = Property%MB%VegDrainedMass + DrainageMassFromVeg * 1e-3
                            endif

                            !m3
                            VegetationNewVolume = VegetationNewVolume - DrainageVolume
                            
                            !g
                            VegetationNewMass   = VegetationNewMass - DrainageMassFromVeg
                            
                            !Update the new conc now with drainage flux
                            if (VegetationNewVolume .gt. 0.0) then
                                !g/m3 = g / (m * m2)
                                Property%VegetationConc(i,j) = VegetationNewMass / VegetationNewVolume
                            else
                                Property%VegetationConc(i,j) = 0.0
                            endif 
                            
                            Property%VegetationOldMass(i,j) = VegetationNewMass
                        
                        endif
                    enddo
                    enddo

                    if (AtmospherePropertyExists (Me%ObjAtmosphere, Property%ID%IDNumber)) then
                        call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'UpdateVegConcentration - ModuleBasin - ERR20' 
                    else
                        deallocate (AtmConcentration)
                    endif

                
                endif
                
                Property => Property%Next

           enddo
                        
        !compute new concentration after routine ComputePotentialEvapotranspiration
        elseif (WarningString == "Evaporation") then 

            Property => Me%FirstProperty
            do while (associated(Property))
            
                !only for runoff properties and exclude for instance evapotranspiration
                if (Property%Inherited) then

                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then
                        
                            !Mass does not change with evaporation
                            VegetationNewMass = Property%VegetationOldMass(i,j)
                        
                            if (Me%CanopyStorage(i, j) .gt. 0.0) then
                                !g/m3 = g / (m * m2)
                                Property%VegetationConc(i,j) = VegetationNewMass / (Me%CanopyStorage(i, j)                    &
                                                               * Me%CoveredFraction(i, j) * Me%ExtVar%GridCellArea(i, j))
                            else
                                Property%VegetationConc(i,j) = 0.0
                            endif                        

                        endif
                    enddo
                    enddo
                    
                endif
                
                Property => Property%Next
                
            enddo                        
                
        elseif (WarningString == "Pesticide") then
                
            call GetVegetationOptions (Me%ObjVegetation,                                     &
                                       NumberOfPesticides = NumberOfPesticides,              &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'UpdateVegConcentration - ModuleBasin - ERR30'        

            call GetVegetationDT  (Me%ObjVegetation, VegDT, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateVegConcentration - ModuleBasin - ERR40'
        
            do Pest = 1, NumberOfPesticides
                call GetVegetationAerialFluxes(Me%ObjVegetation,                            &
                                                Pest              = Pest,                   &
                                                PesticideIDNumber = PesticideIDNumber,      &
                                                PesticideVegetation = PesticideVegetation,  &
                                                STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'UpdateVegConcentration - ModuleBasin - ERR50'                 
            
                call SearchProperty(Property, PropertyXIDNumber = PesticideIDNumber, STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then

                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                        if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then

                            !g = kg/ha * dt/dtveg * 1e3 g/kg * (area) * 1ha/10000m2
                            !area is from all cell because leaf pesticide did account all area
                            MassAdded = PesticideVegetation(i,j) * Me%CurrentDT / VegDT * 1e3 * Me%ExtVar%GridCellArea(i,j)  &
                                         / 10000.
                            
                            VegetationNewMass = Property%VegetationOldMass(i,j) + MassAdded
                            
                            if (Me%CanopyStorage(i, j) .gt. 0.0) then
                                !g/m3 = g / (m * m2)
                                Property%VegetationConc(i,j) = VegetationNewMass / (Me%CanopyStorage(i, j)                    &
                                                               * Me%CoveredFraction(i, j) * Me%ExtVar%GridCellArea(i, j))
                            else
                                Property%VegetationConc(i,j) = 0.0
                            endif                        
                            
                            Property%VegetationOldMass(i,j) = VegetationNewMass               
                    
                        endif                    
                    
                    enddo 
                    enddo
                else
                    write(*,*)
                    write(*,*) 'Not found pesticide', GetPropertyName(PesticideIDNumber) 
                    write(*,*) 'in runoff property list'
                    stop ' UpdateVegConcentration - ModuleBasin - ERR60'
                endif                                
                
                call UngetVegetationAerialFluxes (Me%ObjVegetation,                            &
                                                Pest              = Pest,                   &
                                                PesticideVegetation = PesticideVegetation,  &
                                                STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'UpdateVegConcentration - ModuleBasin - ERR70' 
                
            enddo    

        elseif (WarningString == "Decay") then
            

            Property => Me%FirstProperty
            do while (associated(Property))
            
                !only for runoff properties and exclude for instance evapotranspiration
                !and properties with decay option
                if (Property%Inherited .and. Property%Decay) then
                    
                    call GetRPDecayRate (Me%ObjRunoffProperties, DecayRate, Property%ID%IDNumber, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'UpdateVegConcentration - ModuleBasin - ERR80' 
                    
                    DT = Me%CurrentDT / 86400.
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Me%ExtVar%BasinPoints(i, j) == WaterPoint) then
                            !WQ process without volume change. only mass change
                           
                            !P = P0*exp(-kt)
                            !g
                            MassSink = min (Property%VegetationOldMass(i,j)                                       &
                                            - Property%VegetationOldMass(i,j) * exp(-DecayRate * DT),             &
                                              Property%VegetationOldMass(i,j))
                            
                            VegetationNewMass = Property%VegetationOldMass(i,j) - MassSink
                            
                            if (Me%CanopyStorage(i, j) .gt. 0.0) then
                                !g/m3 = g / (m * m2)
                                Property%VegetationConc(i,j) = VegetationNewMass / (Me%CanopyStorage(i, j)                    &
                                                               * Me%CoveredFraction(i, j) * Me%ExtVar%GridCellArea(i, j))
                            else
                                Property%VegetationConc(i,j) = 0.0
                            endif                        
                            
                            Property%VegetationOldMass(i,j) = VegetationNewMass
                            
                        endif
                    enddo
                    enddo
                    
                endif
                
                Property => Property%Next
                
            enddo                        
            
        endif                   
    
    end subroutine UpdateVegConcentration

    !--------------------------------------------------------------------------

    subroutine DividePrecipitationMassFluxes
        !Local-----------------------------------------------------------------
        type (T_BasinProperty), pointer             :: Property
        real, dimension(:,:), pointer               :: AtmConcentration
        integer                                     :: STAT_CALL, i, j
        !Begin-----------------------------------------------------------------


        Property => Me%FirstProperty
        do while (associated(Property))
            
            !only for runoff properties, exclude evapotranspiration
            if (Property%Inherited) then
                if (AtmospherePropertyExists (Me%ObjAtmosphere, Property%ID%IDNumber)) then
        
                    call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
                                               Scalar             = AtmConcentration,       &
                                               ID                 = Property%ID%IDNumber,   &
                                               STAT               = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DividePrecipitation - ModuleBasin - ERR01'  
                else
                    allocate (AtmConcentration(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
                    AtmConcentration = 0.0
                endif

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then                
                        

                        if (Me%Coupled%Vegetation) then
                            !kg = kg + (m*m2 * g/m3 * 1e-3 kg/g)
!                            Property%MB%TotalRainMass = Property%MB%TotalRainMass + ((Me%RainUncovered(i,j)                  &
!                                                         + Me%RainCovered(i,j))                                              &
!                                                         * Me%ExtVar%GridCellArea(i, j) * AtmConcentration(i,j)              &
!                                                         * 1e-3)
                        
                            !kg = kg + ((m*m2unc * g/m3 * 1e-3 kg/g)
                            Property%MB%UncoveredRainMass = Property%MB%UncoveredRainMass + (Me%RainUncovered(i,j)           &
                                                         *  (1 - Me%CoveredFraction(i, j))                                   &
                                                         * Me%ExtVar%GridCellArea(i, j) * AtmConcentration(i,j)              &
                                                         * 1e-3)                                                     
                            !input for vegetation leafs
                            !kg = kg + ((m*m2plant * g/m3 * 1e-3 kg/g)
                            Property%MB%CoveredRainMass   = Property%MB%CoveredRainMass + (Me%RainCovered(i,j)               &
                                                         *  Me%CoveredFraction(i, j)                                         &
                                                         * Me%ExtVar%GridCellArea(i, j) * AtmConcentration(i,j)              &
                                                         * 1e-3)
                            
                            Property%MB%TotalRainMass = Property%MB%UncoveredRainMass + Property%MB%CoveredRainMass
                                                                                     
!                            !output for vegetation leafs
!                            !kg = kg + ((m*m2cell * g/m3 * 1e-3 kg/g) - Canopy drainage is height in cell area
!                            Property%MB%VegDrainedMass   = Property%MB%VegDrainedMass + (Me%CanopyDrainage(i,j)              &
!                                                           * Me%ExtVar%GridCellArea(i, j) * Property%VegetationConc(i,j)     &
!                                                           * 1e-3)
                                                           
                            Property%MB%RunoffInputMass   = Property%MB%VegDrainedMass + Property%MB%UncoveredRainMass
                                
                        else !there is no cover rain
                            !kg = kg + (m*m2 * g/m3 * 1e-3 kg/g)
                            Property%MB%TotalRainMass = Property%MB%TotalRainMass + ((Me%RainUncovered(i,j))                 &
                                                         * Me%ExtVar%GridCellArea(i, j) * AtmConcentration(i,j)              &
                                                         * 1e-3) 
                            Property%MB%UncoveredRainMass = Property%MB%TotalRainMass                                           
                        endif
                    endif
                enddo
                enddo
                
                if (AtmospherePropertyExists (Me%ObjAtmosphere, Property%ID%IDNumber)) then
                    call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'DividePrecipitation - ModuleBasin - ERR10' 
                else
                    deallocate (AtmConcentration)
                endif  
            endif
            
            Property => Property%Next
            
        enddo     
    
    end subroutine DividePrecipitationMassFluxes

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
        logical                                     :: Pesticide
        logical                                     :: Decay
        character (Len = StringLength)              :: WarningString
        real, dimension(:,:,:), pointer             :: ECw

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "VegetationProcesses")

        
        call GetVegetationOptions (Me%ObjVegetation,                                     &
                                   ModelNitrogen = ModelNitrogen,                        &
                                   ModelPhosphorus = ModelPhosphorus,                    &
                                   Pesticide       = Pesticide,                          &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'VegetationProcesses - ModuleBasin - ERR010'

        if (ModelNitrogen) then 
            
            if (Me%Coupled%PorousMediaProperties) then
            
                call GetPMPConcentration(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                         ConcentrationX          = Nitrate,                      &
                                         PropertyXIDNumber       = Nitrate_,                     &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) then
                    write (*,*) 'Trying to model nitrogen in vegetation but nitrate property not'
                    write (*,*) 'defined in porousmediaproperties. Check options'                
                    stop 'VegetationProcesses - ModuleBasin - ERR020'
                endif

                call SetSoilConcVegetation (Me%ObjVegetation,                                 &
                                            Nitrate             = Nitrate,                    &
                                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR030' 


                call UnGetPorousMediaProperties (Me%ObjPorousMediaProperties, Nitrate, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR040'   
                  
            else
                write (*,*) 'Can not model nitrogen in vegetation if porous media properties'
                write (*,*) 'model not connected. Check basin keyword.'
                stop 'VegetationProcesses - ModuleBasin - ERR050'
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
                    stop 'VegetationProcesses - ModuleBasin - ERR060'
                endif

                call SetSoilConcVegetation (Me%ObjVegetation,                                 &
                                            InorganicPhosphorus = InorganicPhosphorus,        &
                                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR070'   

                call UnGetPorousMediaProperties (Me%ObjPorousMediaProperties, InorganicPhosphorus, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR080'   

            else
                write (*,*) 'Can not model phosphorus in vegetation if porous media properties'
                write (*,*) 'model not connected. Check basin keyword.'
                stop 'VegetationProcesses - ModuleBasin - ERR090'
            endif
        endif

        if (Me%Coupled%PorousMediaProperties) then

            call GetECw (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                         ECw                     = ECw,                          &
                         STAT                    = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then
            
                call SetECw (Me%ObjVegetation, ECw, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'VegetationProcesses - ModuleBasin - ERR100'

            elseif (STAT_CALL /= NOT_ASSOCIATE_) then
            
                stop 'VegetationProcesses - ModuleBasin - ERR110'
                
            endif
                             
            call UnGetPorousMediaProperties (Me%ObjPorousMediaProperties, ECw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'VegetationProcesses - ModuleBasin - ERR120'

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
        if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR130'
    
        !Points to computed variables to be used in porous media (vegetation uses ModulePorousMedia)
        Me%ExtVar%ActualTranspiration => ActualTranspiration

        !update leaf pesticide conc with applications
        if (Pesticide) then
                
            WarningString = "Pesticide"
            call UpdateVegConcentration(WarningString)
                                                
        endif
        
        if (Me%Coupled%RunoffProperties) then
            !update leaf properties with decay
            call GetRPOptions (Me%ObjRunoffProperties, Decay = Decay, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'VegetationProcesses - ModuleBasin - ERR140'
            
            if (Decay) then

                WarningString = "Decay"
                call UpdateVegConcentration(WarningString)
            
            endif
        endif
        
        if (MonitorPerformance) call StopWatch ("ModuleBasin", "VegetationProcesses")
    
    end subroutine VegetationProcesses

    !--------------------------------------------------------------------------
    
    subroutine SimpleInfiltration
    
        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: I, J

        integer, parameter                          :: ChunkSize = 10
        integer                                     :: Chunk
        
        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleBasin", "SimpleInfiltration")
        
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
         
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
            
                if (Me%ExtUpdate%Watercolumn(i, j) > 0.0) then


                    !Infiltration Rate
                    Me%SI%InfRate%Field(i, j) = Me%SI%Ks%Field(i, j) + Me%SI%Ks%Field(i, j) * (Me%SI%MP%Field(i, j) * &
                                                (Me%SI%ThetaS%Field(i, j) - Me%SI%ThetaI%Field(i, j))) /              &
                                                max(Me%SI%AccInf%Field(i, j), 0.001)
                    
                    
                    if (Me%SI%InfRate%Field(i, j) > Me%ExtUpdate%Watercolumn(i, j) / Me%CurrentDT) then
                        Me%SI%InfRate%Field(i, j) = Me%ExtUpdate%Watercolumn(i, j) / Me%CurrentDT
                    endif

                    Me%SI%AccInf%Field(i, j)        = Me%SI%AccInf%Field(i, j) + Me%SI%InfRate%Field(i, j) * Me%CurrentDT

                    Me%SI%TimeWithNoWC%Field(i, j)  = 0.0

                else

                    !Sets Infiltration Rate to zero
                    Me%SI%InfRate%Field(i, j) = 0.0

                    Me%SI%TimeWithNoWC%Field(i, j)  = Me%SI%TimeWithNoWC%Field(i, j) + Me%CurrentDT
                    
                    if (Me%SI%TimeWithNoWC%Field(i, j) > 3600.0 * 3.0) then
                   
                        !Resets accumulated infiltration
                        Me%SI%AccInf%Field(i, j) = AllmostZero
                        
                    endif
                    
                endif
                
                !m                            = m - m/s * s
                Me%ExtUpdate%WaterLevel(i, j) = Me%ExtUpdate%WaterLevel (i, j) - Me%SI%InfRate%Field(i, j) * &
                                                (1.0- Me%SI%ImpFrac%Field(i, j)) * Me%CurrentDT
                
                !mm /hour
                Me%InfiltrationRate (i, j)    = Me%SI%InfRate%Field(i, j) * 1000.0 * 3600.0
                
                !mm /hour
                Me%EVTPRate         (i, j)    = 0.0
                
                !m - Accumulated Infiltration of the entite area
                Me%AccInfiltration  (i, j)    = Me%AccInfiltration  (i, j) + Me%SI%InfRate%Field(i, j) *     &
                                                (1.0- Me%SI%ImpFrac%Field(i, j)) * Me%CurrentDT
                
                !m
                Me%AccEVTP          (i, j)    = 0.0
                
            endif
            
           
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "SimpleInfiltration")

    
    end subroutine SimpleInfiltration

    !--------------------------------------------------------------------------

    subroutine OverLandProcesses

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
!        real, dimension(:, :), pointer              :: FlowX, FlowY
!        real, dimension(:, :), pointer              :: OLFlowToChannels
!        real, dimension(:, :), pointer              :: FlowAtBoundary
        integer                                     :: STAT_CALL
!        integer                                     :: i, j
        integer                                     :: ILB, IUB, JLB, JUB
        character (Len = StringLength)              :: UnLockToWhichModules
        character (Len = StringLength)              :: LockToWhichModules
        character (Len = StringLength)              :: OptionsType
        !Begin-----------------------------------------------------------------
            
        if (MonitorPerformance) call StartWatch ("ModuleBasin", "OverLandProcesses")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        !Unlocks External Var for Module Runoff to change Computefaces
        !and water level
        UnLockToWhichModules = 'AllModules'
        OptionsType          = 'ModifyBasin'
        call ReadUnLockExternalVar (UnLockToWhichModules, OptionsType)

!        call ModifyRunOff   (Me%ObjRunOff, Me%WaterColumn, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR01'
        
        !Runoff automatic has the most recent water column
        call ModifyRunOff   (Me%ObjRunOff, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR01'

        !Re-Locks External 
        LockToWhichModules = 'AllModules'
        OptionsType        = 'ModifyBasin'
        call ReadLockExternalVar (LockToWhichModules, OptionsType)
        
        !Level update is made inside the Runoff model. here it would produce duplication
        !Basin only handles water column changes in other modules to send them to Runoff
        
!        call GetOverLandFlow    (Me%ObjRunOff, FlowX, FlowY, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR02'
!
!        call GetFlowToChannels  (Me%ObjRunOff, OLFlowToChannels, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR03'
!
!        call GetFlowAtBoundary  (Me%ObjRunOff, FlowAtBoundary, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR04'
!
!        !Updates Water column
!        do j = JLB, JUB
!        do i = ILB, IUB
!
!            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
!
!                Me%ExtUpdate%WaterLevel (i, j) = Me%ExtUpdate%WaterLevel (i, j) &
!                                     + ( FlowX(i, j) - FlowX(i, j+1)            &
!                                     +   FlowY(i, j) - FlowY(i+1, j)            &
!                                     -  OLFlowToChannels(i, j) -                &
!                                        FlowAtBoundary(i, j)      )             &
!                                     * Me%CurrentDT /                           &
!                                     Me%ExtVar%GridCellArea(i, j)
!
!            endif
!        enddo
!        enddo
!        
!        call UnGetRunOff (Me%ObjRunOff, FlowX, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR05'
!
!        call UnGetRunOff (Me%ObjRunOff, FlowY, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR06'
!
!        call UnGetRunOff            (Me%ObjRunOff, OLFlowToChannels, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR07' 
!
!        call UnGetRunOff            (Me%ObjRunOff, FlowAtBoundary, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OverLandProcesses - ModuleBasin - ERR08' 

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
        
        GW_Link = GWFlowToChanByCell_
        
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
            
            if (GW_Link == GWFlowToChanByLayer_) then 

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
                    
                    !Get the conc from PMP Old because is explicit flux (is the same used in PMP)
                    call GetPMPConcentrationOld(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                                ConcentrationXOld       = PMPConcentration,             &
                                                PropertyXIDNumber       = PropID,                       &
                                                STAT                    = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR04'
                    
                    !check if flow computation by layers
                    if (GW_Link /= GWFlowToChanByLayer_) then 
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
                    
                    !Get the property conc from RP Old because is explicit flux (is the same used in RP fluxes)
!                    call GetRPConcentrationOld(RunoffPropertiesID    = Me%ObjRunoffProperties,       &
!                                               ConcentrationXOld     = RPConcentration,              &
!                                               PropertyXIDNumber     = PropID,                       &
!                                               STAT                  = STAT_CALL)        
!                    if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkProcesses - ModuleBasin - ERR030'
                    !The conc used for flux with DN is the one saved After Transport
                    call GetRPConcentrationAT(RunoffPropertiesID    = Me%ObjRunoffProperties,       &
                                               ConcentrationXAT      = RPConcentration,              &
                                               PropertyXIDNumber     = PropID,                       &
                                               STAT                  = STAT_CALL)        
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
        
        if (GW_Link /= GWFlowToChanByLayer_) then 
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
            
            if (GW_Link == GWFlowToChanByLayer_) then 
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

        !The column that infiltrates is already computed and now is water column (rain updated water column)

        !Calculates Column which may infiltrate 
!        !Throughfall + Watercolumn over MIN_WC
!        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if(Me%ExtVar%BasinPoints (i,j) == BasinPoint) then
!            
!                !Potencial Infiltration Column = ThroughFall
!                Me%PotentialInfCol (i, j) = Me%ThroughFall(i, j) + Me%WaterColumnCoef * Me%WaterColumn(i, j)
!            end if
!        end do
!        end do

!        !Get the most recent water col from Runof
!        call GetRunoffWaterColumn     (Me%ObjRunoff, WaterColumn, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR00'


        if (Me%Coupled%Vegetation) then

            !if the user choosed to separate transpiration and evaporation,
            !evaporation is explicitly defined
            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then
                
                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
!                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        InfiltrationColumn   = Me%ExtUpdate%Watercolumn,           &
                                        PotentialEvaporation = Me%PotentialEvaporation,            &
                                        ActualTranspiration  = Me%ExtVar%ActualTranspiration,      &
                                        STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR01' 
            
            else !No evaporation

                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
!                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        InfiltrationColumn   = Me%ExtUpdate%Watercolumn,           &
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
!                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        InfiltrationColumn   = Me%ExtUpdate%Watercolumn,           &
                                        PotentialEvaporation = PotentialEvaporation,               &
                                        STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaProcesses - ModuleBasin - ERR21' 
            
            else
                
                call ModifyPorousMedia (ObjPorousMediaID     = Me%ObjPorousMedia,                  &
!                                        InfiltrationColumn   = Me%PotentialInfCol,                 &
                                        InfiltrationColumn   = Me%ExtUpdate%Watercolumn,           &
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
                Me%FlowProduction   (i, j) = - Infiltration      (i, j)               
!                Me%FlowProduction   (i, j) = Me%PotentialInfCol (i, j)      - Infiltration      (i, j) - &
!                                             Me%WaterColumnCoef * Me%WaterColumn(i, j)
            
                !Increase Waterlevel due to flow production
                Me%ExtUpdate%WaterLevel       (i, j) = Me%ExtUpdate%WaterLevel (i, j) + Me%FlowProduction (i, j)
                
                !Output - mm /hour
                Me%InfiltrationRate (i, j) = Infiltration  (i, j) / Me%CurrentDT * 1000.0 * 3600.0

                !Output - mm /hour
                Me%EVTPRate         (i, j) = EfectiveEVTP  (i, j) / Me%CurrentDT * 1000.0 * 3600.0

                !m
                Me%AccInfiltration  (i, j) = Me%AccInfiltration  (i, j) + Infiltration (i,j)

                !m
                Me%AccEVTP          (i, j) = Me%AccEVTP  (i, j)         + EfectiveEVTP (i,j)
                
                if (Me%EVTPOutput%Yes) then
                    !mm         
                    Me%PartialAccEVTP   (i, j) = Me%PartialAccEVTP (i,j) + (EfectiveEVTP (i,j) * 1000) 
                    !Me%PartialAccEVTP   (i, j) = Me%PartialAccEVTP (i,j) + (EfectiveEVTP (i,j) * Me%ExtVar%GridCellArea(i, j)) 
                endif    
                if (Me%EVTPOutput2%Yes) then
                    !mm         
                    Me%PartialAccEVTP2   (i, j) = Me%PartialAccEVTP2 (i,j) + (EfectiveEVTP (i,j) * 1000) 
                endif                 
                  
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
        real, dimension(:,:), pointer               :: HarvestKillAerialBiomass
        real, dimension(:,:), pointer               :: HarvestKillNitrogen
        real, dimension(:,:), pointer               :: HarvestKillPhosphorus
        real, dimension(:,:), pointer               :: HarvestKillRootBiomass
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
        logical                                     :: Grazing
        logical                                     :: HarvestKill
        logical                                     :: Dormancy
        logical                                     :: Fertilization
        logical                                     :: NutrientFluxesWithSoil
        logical                                     :: CoupledSedimentQuality
        logical                                     :: ModelNitrogen, ModelWater
        logical                                     :: ModelPhosphorus
        logical                                     :: GrowthModel
        logical                                     :: Pesticide
        real, dimension(:,:), pointer               :: WindVelocity
        real                                        :: VegetationDT
        real, dimension (:), pointer                :: DNConcentration 
        real, dimension (:,:), pointer              :: RPConcentration
        integer, dimension(:, :), pointer           :: ChannelsID, TranspirationBottomLayer
        integer                                     :: nProperties, iProp, PropID
        integer                                     :: NutrientUptakeMethod
        logical                                     :: PropAdvDiff, PropParticulate !, PropRain, PropIrri
        integer                                     :: NumberOfPesticides
        integer                                     :: PesticideIDNumber, Pest
        real, dimension (:,:), pointer              :: PesticideSoil
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "PorousMediaPropertiesProcesses")

        if (Me%Coupled%RunoffProperties) then
            
            !Send concentrations in runoff (at this time, the potencial infiltration column) to PMP
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
                    
                    !Get the new RP concentration (after surface fluxes) because is the same used in PMP for infil flux
                    !and in MassInFlow (the routine that computes the actualization of Runoff Properties
                    call GetRPConcentration(RunoffPropertiesID       = Me%ObjRunoffProperties,       &
                                             ConcentrationX          = RPConcentration,              &
                                             PropertyXIDNumber       = PropID,                       &
                                             STAT                    = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.4'
                    
                    call SetInfColConcPMP         (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                                   ConcentrationX          = RPConcentration,          &
                                                   PropertyXIDNumber       = PropID,                       &
                                                   STAT                    = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0.5'

                   
                    call UngetRunoffProperties (Me%ObjRunoffProperties, RPConcentration, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR103.6'
                
                endif
                
            enddo  
         endif
        
       
        !Vegetation organic matter and nutrient fluxes to/from soil
        if (Me%Coupled%Vegetation) then

            call GetVegetationOptions (Me%ObjVegetation,                                          &
                                       NutrientFluxesWithSoil = NutrientFluxesWithSoil,           &
                                       Grazing                = Grazing,                          &
                                       HarvestKill            = HarvestKill,                      &
                                       Dormancy               = Dormancy,                         &
                                       Fertilization          = Fertilization,                    &
                                       ModelWater             = ModelWater,                       &
                                       ModelNitrogen          = ModelNitrogen,                    &
                                       ModelPhosphorus        = ModelPhosphorus,                  &
                                       GrowthModel            = GrowthModel,                      &
                                       NutrientUptakeMethod   = NutrientUptakeMethod,             &
                                       Pesticide              = Pesticide,                        &
                                       STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR01'

            call GetTranspirationBottomLayer(VegetationID            = Me%ObjVegetation,                &
                                            Scalar                   = TranspirationBottomLayer,        &
                                            STAT                     = STAT_CALL)               
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR01.2'

            call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties,     &
                                           CoupledVegetation        = .true.,                          &
                                           ModelWater               = ModelWater,                      &
                                           NutrientUptakeMethod     = NutrientUptakeMethod,            &
                                           TranspirationBottomLayer = TranspirationBottomLayer,        &
                                           STAT                     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR01.5'

            call UnGetVegetation  (Me%ObjVegetation, TranspirationBottomLayer, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR01.7'

            
            if (NutrientFluxesWithSoil) then

                call GetVegetationDT  (Me%ObjVegetation, VegetationDT, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR02'

!                call GetRootDepth  (Me%ObjVegetation, RootDepth, STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR03'
                
                !It is needed the same root depth that gave origin to the fluxes.
                !Specially needed in the case of kill where the value here would be zero but
                !fluxes were computed before cutting.
                call GetRootDepthOld  (Me%ObjVegetation, RootDepth, STAT_CALL)
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
                
                    if (HarvestKill) then
                        call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,            &
                                                        HarvestKillAerialBiomass = HarvestKillAerialBiomass,    &
                                                        HarvestKillNitrogen      = HarvestKillNitrogen,         &
                                                        HarvestKillPhosphorus    = HarvestKillPhosphorus,       &
                                                        HarvestKillRootBiomass   = HarvestKillRootBiomass,      &
                                                       STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR050'

                        call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                       HarvestKill              = HarvestKill,                 &
                                                       HarvestKillAerialBiomass = HarvestKillAerialBiomass,    &
                                                       HarvestKillNitrogen      = HarvestKillNitrogen,         &
                                                       HarvestKillPhosphorus    = HarvestKillPhosphorus,       &
                                                       HarvestKillRootBiomass   = HarvestKillRootBiomass,      &
                                                      STAT                      = STAT_CALL)
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

            !get vegetation pesticide fluxes and send them to PMP
            if (GrowthModel .and. Pesticide) then
                call GetVegetationOptions      (VegetationID             = Me%ObjVegetation,           &
                                                NumberOfPesticides       = NumberOfPesticides,         &
                                                STAT                     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0110'

                call GetVegetationDT  (Me%ObjVegetation, VegetationDT, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR112'

                
                do Pest = 1, NumberOfPesticides

                    call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,           &
                                                    Pest                     = Pest,                       &
                                                    PesticideIDNumber        = PesticideIDNumber,          &
                                                    PesticideSoil            = PesticideSoil,              &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0120'
                    

                    call SetVegetationPMProperties(PorousMediaPropertiesID  = Me%ObjPorousMediaProperties, &
                                                    PesticideIDNumber       = PesticideIDNumber,           &
                                                    PesticideSoilFlux       = PesticideSoil,               &
                                                    GrowthModel             = GrowthModel,                 &
                                                    VegetationDT            = VegetationDT,                &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0130'
                enddo                        
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
                    !The DN new concentration is computed at the end of timestep. So the conc used
                    !in PMP here is the new of the last time step or the old of the new that is the same
                    !used in DN for the water discharges                    
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
                                         WaterColumn                = Me%ExtUpdate%Watercolumn,     &
!                                         ThroughFall                = Me%ThroughFall,               &
!                                         WCEvaporated               = Me%WaterColumnEvaporated,     &
                                         WaterColumnOLd             = Me%ExtUpdate%WatercolumnOld,   &
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
               
                if (GrowthModel) then

                    if (Grazing) then
                        call UngetVegetationSoilFluxes (VegetationID             = Me%ObjVegetation,            &
                                                        GrazingBiomass           = GrazingBiomass,              &
                                                        GrazingNitrogen          = GrazingNitrogen,             &   
                                                        GrazingPhosphorus        = GrazingPhosphorus,           &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0130'
                    endif
                    
                    if (HarvestKill) then
                        call UngetVegetationSoilFluxes (VegetationID             = Me%ObjVegetation,            &
                                                        HarvestKillAerialBiomass = HarvestKillAerialBiomass,    &
                                                        HarvestKillNitrogen      = HarvestKillNitrogen,         &
                                                        HarvestKillPhosphorus    = HarvestKillPhosphorus,       &
                                                        HarvestKillRootBiomass   = HarvestKillRootBiomass,      &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0140'
                    endif

                    if (Dormancy) then
                        call UngetVegetationSoilFluxes (VegetationID             = Me%ObjVegetation,            &
                                                        DormancyBiomass          = DormancyBiomass,             &
                                                        DormancyNitrogen         = DormancyNitrogen,            &
                                                        DormancyPhosphorus       = DormancyPhosphorus,          &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0150'
                    endif

                    if (Fertilization) then
                        call UngetVegetationSoilFluxes (VegetationID             = Me%ObjVegetation,           &
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

            if (GrowthModel .and. Pesticide) then
                do Pest = 1, NumberOfPesticides
                    call UnGetVegetationSoilFluxes (VegetationID             = Me%ObjVegetation,         &
                                                    PesticideSoil            = PesticideSoil,            &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0161'                
                enddo
            endif

        endif

!        if (CoupledSedimentQuality) then
!            !Wind Velocity
!            call UnGetAtmosphere  (Me%ObjAtmosphere, WindVelocity, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'PorousMediaPropertiesProcesses - ModuleBasin - ERR0162'
!        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "PorousMediaPropertiesProcesses")


    end subroutine PorousMediaPropertiesProcesses

    !--------------------------------------------------------------------------

!    subroutine ComputePropertyInfilColumn (AtmConcentration, RPConcentration, InfColConcentration)
!    
!        !Arguments-------------------------------------------------------------
!        real, dimension(:,:), pointer            :: AtmConcentration    !IN
!        real, dimension(:,:), pointer            :: RPConcentration     !IN
!        real, dimension(:,:), pointer            :: InfColConcentration !OUT
!        
!        !Local-----------------------------------------------------------------
!        real                                :: RainVolume, MassOnRain
!        real                                :: WaterColumnVolume, MassOnWaterColumn
!        real                                :: InfColumnVolume, MassOnInfColumn
!        integer                             :: i,j !,CHUNK
!
!        !----------------------------------------------------------------------       
!        
!        !! $OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
!
!                !mass of the property on rain
!                !m3        =  m * m2
!                RainVolume = Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j)
!                !g         =     g/m3       *           m                     m2 
!                MassOnRain = AtmConcentration(i,j) * (Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j))
!                
!                !mass of the property in water column
!                !m3 = m * m2
!                WaterColumnVolume = (Me%PotentialInfCol(i,j) - Me%ThroughFall(i, j)) * Me%ExtVar%GridCellArea(i, j)
!                !g = g/m3 * m3
!                MassOnWaterColumn = RPConcentration(i,j) * WaterColumnVolume
!                
!                !mass of the property in infiltration column
!                MassOnInfColumn   = MassOnRain + MassOnWaterColumn
!                InfColumnVolume   = RainVolume + WaterColumnVolume
!                
!                !Compute infiltration column concetration 
!                if (Me%PotentialInfCol(i,j) .gt. 0.0) then
!                    !g/m3 = g  /  m3
!                    InfColConcentration(i, j) = MassOnInfColumn / InfColumnVolume
!                else
!                    InfColConcentration(i,j) = 0.0
!                endif
!                
!               
!            endif
!        enddo
!        enddo
!        !! $OMP END DO 
!        
!        !! $OMP END PARALLEL
!    
!    
!    end subroutine ComputePropertyInfilColumn
!
!    !--------------------------------------------------------------------------


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
        logical                                    :: SplashErosion !, ModelCanopyHeight
        real, dimension(:,:), pointer              :: CanopyHeight
        real, dimension(:,:), pointer              :: FertilOrganicNParticFluff
        real, dimension(:,:), pointer              :: FertilOrganicPParticFluff
        logical                                    :: NutrientFluxesWithSoil, Fertilization
        logical                                    :: ModelNitrogen, ModelPhosphorus
        logical                                    :: GrowthModel, VegParticFertilization, Pesticide
        real                                       :: VegetationDT
        integer                                     :: NumberOfPesticides        
        integer                                     :: PesticideIDNumber, Pest
        real, dimension (:,:), pointer              :: PesticideSoil       
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
                    !The DN new concentration is computed at the end of timestep. So the conc used
                    !in runoff here is the new of the last time step or the old of the new that is the same
                    !used in DN for the water discharges
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
                if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - RunoffPropertiesProcesses - ERR109'
                
            enddo   
            
            call UnGetDrainageNetwork (Me%ObjDrainageNetwork, ChannelsID, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - RunoffPropertiesProcesses - ERR109a'      
             
        endif
        
        !See if splash erosion is computed
        call GetRPOptions (RunoffPropertiesID = Me%ObjRunoffProperties, SplashErosion = SplashErosion, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0110'
        
        !if true need troughfall from basin and canopy height from vegetation
        if (SplashErosion) then
            
            if (Me%Coupled%Vegetation) then
!                call GetVegetationOptions (VegetationID = Me%ObjVegetation, ModelCanopyHeight = ModelCanopyHeight, &
!                                           STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0120'
!                
!                if (ModelCanopyHeight) then
                    call GetCanopyHeight (VegetationID = Me%ObjVegetation, Scalar = CanopyHeight, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0130'
!                endif
                if (.not. associated (CanopyHeight)) then
                    write(*,*)'For Splash Erosion need to be defined canopy height in ModuleVegetation'
                    stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0135'
                endif
            else
                allocate(CanopyHeight(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                CanopyHeight = 0.0
            endif
        
            call SetBasinToRPSplash   (RunoffPropertiesID      = Me%ObjRunoffProperties,        &
                                       ThroughFall             = Me%ThroughFall,                &
                                       CanopyDrainage          = Me%CanopyDrainage,             &
                                       CanopyHeight            = CanopyHeight,                  &
                                       STAT                    = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR140'

            if (Me%Coupled%Vegetation) then
                call UngetVegetation (Me%ObjVegetation, CanopyHeight, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0140'
            else        
                deallocate (CanopyHeight)
            endif
        
        endif
        
        !Vegetation organic matter fluxes to fluffly layer - particulate organic matter from fertilization
        if (Me%Coupled%Vegetation) then

            call GetVegetationOptions (Me%ObjVegetation,                                         &
                                      NutrientFluxesWithSoil = NutrientFluxesWithSoil,           &
                                      Fertilization          = Fertilization,                    &
                                      ModelNitrogen          = ModelNitrogen,                    &
                                      ModelPhosphorus        = ModelPhosphorus,                  &
                                      GrowthModel            = GrowthModel,                      &
                                      VegParticFertilization = VegParticFertilization,           &
                                      Pesticide              = Pesticide,                        &
                                      STAT                   = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0160'
            
            call GetVegetationDT (Me%ObjVegetation, VegetationDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0165'
            
            if (GrowthModel) then
                
                if (Fertilization .and. VegParticFertilization) then
                
                    if (ModelNitrogen) then
                        call GetVegetationSoilFluxes(VegetationID              = Me%ObjVegetation,           &
                                                     FertilOrganicNParticFluff = FertilOrganicNParticFluff,  &
                                                     STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0170'   
                        
                        call SetVegetationRP        (RunoffPropertiesID        = Me%ObjRunoffProperties,     &
                                                           FertilOrganicNParticFluff = FertilOrganicNParticFluff,  &
                                                           VegetationDT              = VegetationDT,               &
                                                           STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0172'   
                        
                        call UnGetVegetationSoilFluxes(VegetationID              = Me%ObjVegetation,           &
                                                       FertilOrganicNParticFluff = FertilOrganicNParticFluff,  &
                                                       STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0174'  
                    endif
                    
                    if (ModelPhosphorus) then
                         call GetVegetationSoilFluxes(VegetationID              = Me%ObjVegetation,           &
                                                     FertilOrganicPParticFluff = FertilOrganicPParticFluff,   &
                                                     STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0176'                       
                        
                        call SetVegetationRP        (RunoffPropertiesID        = Me%ObjRunoffProperties,     &
                                                           FertilOrganicPParticFluff = FertilOrganicPParticFluff,  &
                                                           VegetationDT              = VegetationDT,               &
                                                           STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0178'   
                        
                        call UnGetVegetationSoilFluxes(VegetationID              = Me%ObjVegetation,           &
                                                       FertilOrganicPParticFluff = FertilOrganicPParticFluff,  &
                                                       STAT                      = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0179'        
                    endif
                
                endif
            
                !get vegetation pesticide fluxes and send them to RP - computed if water column exists
                if (Pesticide) then
                    call GetVegetationOptions      (VegetationID             = Me%ObjVegetation,           &
                                                    NumberOfPesticides       = NumberOfPesticides,         &
                                                    STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0180'

                    call GetVegetationDT  (Me%ObjVegetation, VegetationDT, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR190'

                    
                    do Pest = 1, NumberOfPesticides

                        call GetVegetationSoilFluxes   (VegetationID             = Me%ObjVegetation,           &
                                                        Pest                     = Pest,                       &
                                                        PesticideIDNumber        = PesticideIDNumber,          &
                                                        PesticideSoil            = PesticideSoil,              &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0200'
                        

                        call SetVegetationRP(RunoffPropertiesID  = Me%ObjRunoffProperties, &
                                                        PesticideIDNumber       = PesticideIDNumber,           &
                                                        PesticideSoilFlux       = PesticideSoil,               &
                                                        VegetationDT            = VegetationDT,                &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0210'
                        
                        call UnGetVegetationSoilFluxes (VegetationID             = Me%ObjVegetation,         &
                                                        PesticideSoil            = PesticideSoil,            &
                                                        STAT                     = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0220'                
                        
                    enddo
                endif            
            
            endif
            
            
        endif
        
        
        call ModifyRunoffProperties(ObjRunoffPropertiesID           = Me%ObjRunoffProperties,  &
 !                                        ThroughFall                = Me%ThroughFall,               &
 !                                        WCEvaporated               = Me%WaterColumnEvaporated,     &
                                         STAT                       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0190' 
           
        
!        if (CoupledSedimentQuality) then
!            !Wind Velocity
!            call UnGetAtmosphere  (Me%ObjAtmosphere, WindVelocity, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'RunoffPropertiesProcesses - ModuleBasin - ERR0132'
!        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "RunoffPropertiesProcesses")


    end subroutine RunoffPropertiesProcesses

    !--------------------------------------------------------------------------

    subroutine ActualizeWaterColumn ()

        !Arguments-------------------------------------------------------------
        !character (Len = *), intent(in)             :: WarningString

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, CHUNK, STAT_CALL

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "ActualizeWaterColumn")

        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                Me%ExtUpdate%WatercolumnOld(i,j) = Me%ExtUpdate%Watercolumn(i,j)
                
                Me%ExtUpdate%Watercolumn(i, j) = Me%ExtUpdate%WaterLevel(i, j) - Me%ExtVar%Topography(i, j)                

                if (Me%ExtUpdate%Watercolumn(i, j) < 0.0) then
                    
                    !Mass Error
                    Me%ExtVar%MassError     (i, j) = Me%ExtVar%MassError (i, j) -        &
                                                     Me%ExtUpdate%Watercolumn(i, j) *    &
                                                     Me%ExtVar%GridCellArea(i, j)

                    Me%ExtUpdate%Watercolumn(i, j) = 0.0
                    Me%ExtUpdate%WaterLevel (i, j) = Me%ExtVar%Topography(i, j)
                
                end if
            endif

        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
        !Send water column to runoff - the beholder of the water column
        call SetBasinColumnToRunoff (ObjRunOffID             = Me%ObjRunoff,                 &
                                     WaterColumnOld          = Me%ExtUpdate%WatercolumnOld,  &
                                     WaterColumn             = Me%ExtUpdate%Watercolumn,     &
                                     STAT                    = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumn - ModuleBasin - ERR30'        

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
        real(8)                                     :: DebugMass
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
            
!            if (PropAdvDiff) then
                
                allocate(MassInFlow(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                MassInFlow = 0.0
                allocate(MassToBottom(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                MassToBottom = 0.0

                !Get the most recent conc. from RP (e.g. is the same used in PMP)
                call GetRPConcentration(RunoffPropertiesID       = Me%ObjRunoffProperties,       &
                                         ConcentrationX          = RPConcentration,              &
                                         PropertyXIDNumber       = PropID,                       &
                                         STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR20'
                
                allocate(NewRPConcentration(Me%WorkSize%ILB:Me%WorkSize%IUB,Me%WorkSize%JLB:Me%WorkSize%JUB))
                call SetMatrixValue (NewRPConcentration, Me%WorkSize, RPConcentration)

                !Compute mass flow matrix to update concentrations
                call ComputeMassInFlow (WarningString, RPConcentration, PropID, PropParticulate, MassInFlow)

                !Compute the new RP conc based on the new fluxes
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                        !g = g/m3 * (m * m2)
                        PropertyMassOld = RPConcentration(i,j) * (Me%ExtUpdate%WatercolumnOld(i,j) * Me%ExtVar%GridCellArea(i,j))
                        
                        PropertyMassNew = PropertyMassOld + MassInFlow(i,j)
                        
                        DebugMass = MassInFlow(i,j)
                        
                        if (Me%ExtUpdate%Watercolumn(i,j) .gt. AlmostZero) then
                            !g/m3 = g / (m * m2)
                            NewRPConcentration(i,j) = PropertyMassNew / (Me%ExtUpdate%Watercolumn(i,j)        &
                                                                         * Me%ExtVar%GridCellArea(i,j))
                        else
                            !Do not zero, leave unchanged (RP Conc)
                            !NewRPConcentration(i,j) = 0.0
                            
                            !deposition driven by complete infiltration of water column (WC totally infiltrated in time step) 
                            if ((PropParticulate) .and. (Me%ExtUpdate%WatercolumnOld(i,j) .gt. AlmostZero)) then
                                MassToBottom(i,j) = PropertyMassOld
                            endif  
                            
                            !IT IS NEEDED A FORMULATION SIMILAR FOR DISSOLVED PROPERTIES IF EVAPORATION REMOVES ALL WATER
                            !COLUMN IN ONE TIME STEP. THIS MAY BE DONE WITH MASS MATRIX HERE AND IN RUNOFFPROPERTIES
                            !The problem is that RunoffProperties has to update mass even if cell returned no water because
                            !the water may have been removed by infiltration.
                                                      
                        endif                
                        
                      
                    endif

                enddo
                enddo
 
                call UngetRunoffProperties(Me%ObjRunoffProperties, RPConcentration, STAT_Call)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR40'

                !Send new conc to RP module - the beholder of the concentrations
                call SetBasinConcRP (RunoffPropertiesID        = Me%ObjRunoffProperties,       &
                                       BasinConcentration      = NewRPConcentration,           &
                                       PropertyXIDNumber       = PropID,                       &
                                       MassToBottom            = MassToBottom,                 &
                                       STAT                    = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ActualizeWaterColumnConcentration - ModuleBasin - ERR70'
                
                deallocate(MassInFlow)
                deallocate(NewRPConcentration)
                deallocate(MassToBottom)
                
!            endif
            
        enddo

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "ActualizeWaterColumnConcentration")

    end subroutine ActualizeWaterColumnConc

    !--------------------------------------------------------------------------

    subroutine ComputeMassInFlow(WarningString, RPConcentration, PropID, Particulate, MassInFlow)
    
        !Arguments-------------------------------------------------------------
        character (Len = *), intent(in)             :: WarningString
        real(8), dimension(:,:), pointer            :: MassInFlow
        real, dimension(:, :   ), pointer           :: RPConcentration
        real, dimension(:, : ,:), pointer           :: PMPConcentration
        real, dimension(:, :   ), pointer           :: AtmConcentration
        integer                                     :: PropID,i,j, k, STAT_CALL
        type (T_Size3D)                             :: WorkSize3D   
        logical                                     :: Particulate  
        real(8), dimension(:, :), pointer           :: Infiltration
        real(8)                                     :: MassInRain, MassInDrainage
        type (T_BasinProperty), pointer             :: Property 
        !Local-----------------------------------------------------------------
        
        MassInFlow = 0.0
        
        if(WarningString == 'AtmosphereProcesses') then
        
            !!Precipitation flux
            if (Me%Coupled%Atmosphere) then

                if (AtmospherePropertyExists (Me%ObjAtmosphere, PropID)) then
                
                    call GetAtmosphereProperty(AtmosphereID       = Me%ObjAtmosphere,       &
                                               Scalar             = AtmConcentration,       &
                                               ID                 = PropID,                 &
                                               STAT               = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR07'

                else
                    allocate (AtmConcentration(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
                    AtmConcentration = 0.0
                endif

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then 
                        if (Me%Coupled%Vegetation) then
                            !mass of the property on rain in uncovered area + mass of property from leafs in covered area
                            !g         =     (g/m3) * m *  m2uncov 
                            MassInRain = AtmConcentration(i,j) * Me%RainUncovered(i, j) * (1 - Me%CoveredFraction(i,j))  &
                                          * Me%ExtVar%GridCellArea(i, j)
                            
                            call SearchProperty(Property, PropID, .true., STAT = STAT_CALL) 
                            if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR07.5'

                            !g         
                            MassInDrainage = Property%VegetationDrainage(i,j) 
                             
!                            !g         =     (g/m3) * m *  m2cell - - Canopy drainage is height in cell area 
!                            MassInDrainage = Property%VegetationConc(i,j) * Me%CanopyDrainage(i, j)                     &
!                                              * Me%ExtVar%GridCellArea(i, j) 
                                          
                            MassInFlow(i,j) = MassInRain +  MassInDrainage                 
                        
                        else            
                            !mass of the property on rain
                            !g         =     g/m3 * m *  m2 
                            MassInFlow(i,j) = (AtmConcentration(i,j) * Me%ThroughFall(i, j)                              &
                                              * Me%ExtVar%GridCellArea(i, j))
                        endif
                    endif
                enddo
                enddo

                if (AtmospherePropertyExists (Me%ObjAtmosphere, PropID)) then
                    call UngetAtmosphere (Me%ObjAtmosphere, AtmConcentration, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR08' 
                else
                    deallocate (AtmConcentration)
                endif
            else
                !do nothing - if not atmosphere property mass is zero
            endif 

            !Evaporation flux from water colum or plants
!            Do nothing to evaporation. it takes water but not mass from water column
!            MassInFlow = 0.0
        
        elseif(WarningString == 'PorousMediaProcesses') then

            !Flux with PorousMedia - runoff concentration is already the product of the mixing with rain

            !Particulate properties do not enter the soil or exit (particulate do not have advection diffusion in soil)
            if (.not. Particulate) then

                !Get PMP old conc because is explicit infiltration flux (is the same used in PMP fluxes)
                !RP concentration used is the new because is the same used in PMP
                call GetPMPConcentrationOld(PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                            ConcentrationXOld       = PMPConcentration,             &
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
                            MassInFlow(i,j) = RPConcentration(i,j) * (-Infiltration(i,j) * Me%ExtVar%GridCellArea(i,j))
                        
                        else ! negative infiltration (exfiltration) - adding mass to WC or zero
                            
                            k = WorkSize3D%KUB
                            MassInFlow(i,j) = PMPConcentration (i,j,k) * (-Infiltration(i,j) * Me%ExtVar%GridCellArea(i,j))
                        endif
                        
                    endif
                enddo
                enddo          
                
                call UnGetPorousMedia     (Me%ObjPorousMedia, Infiltration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR05.7'

                call UngetPorousMediaProperties(Me%ObjPorousMediaProperties, PMPConcentration, STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ComputeMassInFlow - ModuleBasin - ERR06'
                
!            else
!
!                if (Me%Coupled%Atmosphere .and. (AtmospherePropertyExists (Me%ObjAtmosphere, PropID))) then            
!                    !Particulate properties do not enter soil but for now may come from rain - for testing
!                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!                        if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
!                            !mass of the property on rain
!                            !g         =     g/m3       *           m                     m2 
!                            MassInFlow(i,j) = AtmConcentration(i,j) * (Me%ThroughFall(i, j) * Me%ExtVar%GridCellArea(i, j))
!
!                        endif
!                    enddo
!                    enddo
!                endif                
            endif

        
        elseif (WarningString == 'PorousMediaProcesses 2') then
            
            !Precipitation Flux already accounted in "Atmosphere Processes"

        elseif (WarningString == 'SimpleInfiltration') then
            
            !mass flux with "soil"
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
        integer                                     :: i, j, STAT_CALL
        
        !Me%ExtUpdate%WaterColumnOld => Me%ExtUpdate%WaterColumn
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                if (Me%ExtUpdate%WaterColumn(i, j) > 0.001) then
                    Me%WaterColumnRemoved(i, j) = Me%ExtUpdate%WaterColumn(i, j) / Me%WCRemovalTime 
                    Me%ExtUpdate%WaterColumnOld(i,j) = Me%ExtUpdate%WaterColumn(i, j) 
                    Me%ExtUpdate%WaterColumn(i, j)  = Me%ExtUpdate%WaterColumn(i, j) - (Me%WaterColumnRemoved(i, j) * Me%CurrentDT)
                    if (Me%ExtUpdate%WaterColumn(i, j) < 0.0) then
                        Me%ExtUpdate%WaterColumn(i, j)    = 0.0
                        Me%ExtUpdate%WaterLevel (i, j)    = Me%ExtVar%Topography(i, j)
                    else
                        Me%ExtUpdate%WaterLevel (i, j)    = Me%ExtVar%Topography(i, j) + Me%ExtUpdate%WaterColumn(i, j)
                    endif
                else
                    Me%WaterColumnRemoved(i, j) = 0.0
                endif
                !Convert to m3/s
                Me%WaterColumnRemoved(i, j) = Me%WaterColumnRemoved(i, j) * 1.413E+08
            endif

        enddo
        enddo

        !Send water column to runoff
        call SetBasinColumnToRunoff (ObjRunOffID             = Me%ObjRunoff,                    &
                                     WaterColumnOld          = Me%ExtUpdate%WatercolumnOld,     &
                                     WaterColumn             = Me%ExtUpdate%Watercolumn,        &
                                     STAT                    = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'RemoveWaterColumn - ModuleBasin - ERR30'     

    end subroutine RemoveWaterColumn

    !--------------------------------------------------------------------------
    
    subroutine GetStoredVolumes(instant)
    
        !Arguments-------------------------------------------------------------
        integer :: instant
        
        !Local-----------------------------------------------------------------
        integer :: stat
        
        !Begin-----------------------------------------------------------------
        
        !RunOff
        if (Me%Coupled%Runoff) then                   
            if (instant == InitialInstant) then
                call GetRunOffStoredVolumes (Me%ObjRunoff,                                  &
                                             Surface     = Me%BWB%IniStoredInSurface,       &
                                             StormSystem = Me%BWB%IniStoredInStormWater,    &
                                             STAT = stat)
                if (stat /= SUCCESS_) stop 'GetStoredVolumes - ModuleBasin - ERR010'            
            else !instant = FinalInstant
                call GetRunOffStoredVolumes (Me%ObjRunoff,                                  &
                                             Surface     = Me%BWB%FinStoredInSurface,       &
                                             StormSystem = Me%BWB%FinStoredInStormWater,    &
                                             STAT = stat)
                if (stat /= SUCCESS_) stop 'GetStoredVolumes - ModuleBasin - ERR020'  
            endif                
        endif        
        
        !PorousMedia
        if (Me%Coupled%PorousMedia) then            
            if (instant == InitialInstant) then
                call GetPMStoredVolume (Me%ObjPorousMedia, StoredVolume = Me%BWB%IniStoredInSoil, STAT = stat)
                if (stat /= SUCCESS_) stop 'GetStoredVolumes - ModuleBasin - ERR030'
            else !instant = FinalInstant
                call GetPMStoredVolume (Me%ObjPorousMedia, StoredVolume = Me%BWB%FinStoredInSoil, STAT = stat)
                if (stat /= SUCCESS_) stop 'GetStoredVolumes - ModuleBasin - ERR040'
            endif            
        endif

        !DrainageNetwork
        if (Me%Coupled%DrainageNetwork) then            
            if (instant == InitialInstant) then
                call GetDNStoredVolume(Me%ObjDrainageNetwork, Me%BWB%IniStoredInChannels, STAT = stat)
                if (stat /= SUCCESS_) stop 'GetStoredVolumes - ModuleBasin - ERR050'
            else !instant = FinalInstant
                call GetDNStoredVolume(Me%ObjDrainageNetwork, Me%BWB%FinStoredInChannels, STAT = stat)
                if (stat /= SUCCESS_) stop 'GetStoredVolumes - ModuleBasin - ERR060' 
            endif              
        endif       
        
        !Vegetation
        !It's already stored
        
        if (Me%BWB%StoreInitial) Me%BWB%StoreInitial = .false.
        !----------------------------------------------------------------------
        
    end subroutine GetStoredVolumes
    
    !--------------------------------------------------------------------------
    
 !   subroutine CalculateMass(VolumeBasin, VolumeVegetation, VolumePorousMedia, VolumeChannels)
    subroutine CalculateMass(Time)

        !Arguments-------------------------------------------------------------
        real(8)                                     :: VolumeRunoff
        real(8)                                     :: VolumePorousMedia, VolumeChannels
        character (Len = StringLength)               :: Time
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real(8)                                     :: PMPTotalStoredMass, RPTotalStoredMass
        real(8)                                     :: DNTotalStoredMass
        type (T_BasinProperty), pointer             :: PropertyX
        !Begin-----------------------------------------------------------------

        !Calculates VolumeRunoff / VolumeVegetation
        if (Me%Coupled%Runoff) then
            call GetRunoffTotalStoredVolume (Me%ObjRunoff, VolumeRunoff, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR00'        
            
            if (Time == "Initial") then
                Me%MB%IniVolumeRunoff     = VolumeRunoff
            elseif (Time == "Final") then
                Me%MB%FinalVolumeRunoff   = VolumeRunoff
            endif        
        
        endif
        
        !Vegetation is now the only thing to be cared by Basin
        if (Time == "Initial") then
            Me%MB%IniVolumeVegetation = Me%VolumeVegetation
        elseif (Time == "Final") then
            Me%MB%FinalVolumeVegetation   = Me%VolumeVegetation
        endif
        
        if (Me%Coupled%PorousMedia) then
            call GetPorousMediaTotalStoredVolume (Me%ObjPorousMedia, VolumePorousMedia, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR01'
            
            if (Time == "Initial") then
                Me%MB%IniVolumePorousMedia   = VolumePorousMedia
            elseif (Time == "Final") then
                Me%MB%FinalVolumePorousMedia     = VolumePorousMedia
            endif            
        endif

        if (Me%Coupled%DrainageNetwork) then
            call GetVolumes(Me%ObjDrainageNetwork, TotalStoredVolume = VolumeChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR02'
            
            if (Time == "Initial") then
                Me%MB%IniVolumeChannels      = VolumeChannels
            elseif (Time == "Final") then
                Me%MB%FinalVolumeChannels        = VolumeChannels 
            endif              
        endif

        PropertyX => Me%FirstProperty
        do while(associated(PropertyX))  
            
            if (PropertyX%AdvectionDiffusion) then
                
                if (Me%Coupled%Vegetation) then
                    if (Time == "Initial") then
                        PropertyX%MB%InitialVegetationStoredMass = PropertyX%VegTotalStoredMass
                    elseif (Time == "Final") then
                        PropertyX%MB%FinalVegetationStoredMass   = PropertyX%VegTotalStoredMass
                    endif                    
                endif
                
                if (Me%Coupled%PorousMediaProperties .and. .not. PropertyX%Particulate) then
                    call GetPMPMassBalance (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,  &
                                            PropertyID              = PropertyX%ID%IDNumber,        &
                                            TotalStoredMass         = PMPTotalStoredMass,           &
                                            STAT                    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR03'
                    
                    if (Time == "Initial") then
                        PropertyX%MB%InitialPMPStoredMass = PMPTotalStoredMass
                    elseif (Time == "Final") then
                        PropertyX%MB%FinalPMPStoredMass   = PMPTotalStoredMass
                    endif
                endif

                if (Me%Coupled%RunoffProperties) then
                    call GetRPMassBalance (RunoffPropertiesID = Me%ObjRunoffProperties,              &
                                           PropertyID         = PropertyX%ID%IDNumber,               &
                                           TotalStoredMass    = RPTotalStoredMass,                   &
                                           STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR04'
                    
                    if (Time == "Initial") then
                        PropertyX%MB%InitialRPStoredMass = RPTotalStoredMass
                    elseif (Time == "Final") then
                        PropertyX%MB%FinalRPStoredMass   = RPTotalStoredMass
                    endif
                endif
   
                if (Me%Coupled%DrainageNetwork) then
                    call GetDNMassBalance (DrainageNetworkID = Me%ObjDrainageNetwork,                &
                                           PropertyID        = PropertyX%ID%IDNumber,                &
                                           TotalStoredMass   = DNTotalStoredMass,                    &
                                           STAT              = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CalculateMass - ModuleBasin - ERR05'
                    
                    if (Time == "Initial") then
                        PropertyX%MB%InitialDNStoredMass = DNTotalStoredMass
                    elseif (Time == "Final") then
                        PropertyX%MB%FinalDNStoredMass   = DNTotalStoredMass
                    endif
                endif
            endif
            
            PropertyX => PropertyX%Next
        enddo
        
   
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
        real, dimension(:,:), pointer               :: PotentialEvapoTranspiration, CropEvapotrans
        real, dimension(:,:), pointer               :: ActualTranspiration
        real(8), dimension(:,:), pointer            :: ActualEvaporation
        real, dimension(:,:), pointer               :: ActualTP, ActualEVAP
        type (T_BasinProperty), pointer             :: Property

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "TimeSerieOutput")

        !Watercolumn - should be done in runoff
        call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                              &
                             Data2D_8    = Me%ExtUpdate%Watercolumn,                     &                             
                             STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR02'

        !Waterlevel
        call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                              &
                             Data2D_8    = Me%ExtUpdate%WaterLevel,                      &                             
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

            !Canopy Drainage
            call WriteTimeSerie (TimeSerieID = Me%ObjTimeSerie,                     &
                                 Data2D_8    = Me%CanopyDrainage,                   &                             
                                 STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR091'

            allocate(CropEvapotrans(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

            !Convert Units to mm/h
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                    !m/s * 1000 mm/m * 3600 s/h = mm/h
                    CropEvapotrans(i,j)   = Me%CropEvapotrans(i,j) * 1000. * 3600.
                endif
            enddo
            enddo
            
            !Potential Crop Evapotranspiration
            call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                 Data2D = CropEvapotrans,                           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR07'            

            deallocate (CropEvapotrans)

            if (Me%EvapoTranspirationMethod == SeparateEvapoTranspiration) then

                !m3/s
                call GetTranspiration(Me%ObjVegetation, ActualTranspiration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR095'
                !m
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
                        !mm/h                       =           m/s                * mm/m  * s/h
                        PotentialEvaporation(i,j)   = Me%PotentialEvaporation(i,j) * 1000. * 3600.                      
                        
                        !mm/h                       =           m/s                  * mm/m  * s/h 
                        PotentialTranspiration(i,j) = Me%PotentialTranspiration(i,j) * 1000. * 3600.
                        
                        !mm/h           =           m            /       s      * 1000 mm/m * 3600s/h
                        ActualEVAP(i,j) = ActualEvaporation(i,j) / Me%CurrentDT * 1000.0 * 3600.0
                        if (Me%EVTPOutput%Yes) then
                            !mm                         =         mm            +    (          m            * mm/m)
                            Me%PartialAccEVAct   (i, j) = Me%PartialAccEVAct (i,j) + (ActualEvaporation(i,j) * 1000) 
                        endif    
                        if (Me%EVTPOutput2%Yes) then
                            !mm                         =              mm            +    (          m          * mm/m)          
                            Me%PartialAccEVAct2   (i, j) = Me%PartialAccEVAct2 (i,j) + (ActualEvaporation (i,j) * 1000) 
                        endif                           
                        
                        !mm/h           =           m3/s           /              m2             * mm/m  * s/h
                        ActualTP(i,j)   = ActualTranspiration(i,j) / Me%ExtVar%GridCellArea(i,j) * 1000. * 3600.
                        if (Me%EVTPOutput%Yes) then
                            !mm                         =          mm              + (           m3/s          /            m2               * mm/m *     s       ) 
                            Me%PartialAccETAct   (i, j) = Me%PartialAccETAct (i,j) + (ActualTranspiration(i,j) / Me%ExtVar%GridCellArea(i,j) * 1000 * Me%CurrentDT) 
                        endif    
                        if (Me%EVTPOutput2%Yes) then
                            !mm                         =          mm                + (           m3/s           /         m2                  * mm/m *     s       )         
                            Me%PartialAccETAct2   (i, j) = Me%PartialAccETAct2 (i,j) + (ActualTranspiration (i,j) / Me%ExtVar%GridCellArea(i,j) * 1000 * Me%CurrentDT) 
                        endif                         
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
            
            Property => Me%FirstProperty
            do while (associated (Property))
                if (Property%Inherited) then
                    call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                         Data2D_8 = Property%VegetationConc,                &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR120'  
                    call WriteTimeSerie (Me%ObjTimeSerie,                                   &
                                         Data2D_8 = Property%VegetationOldMass,             &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'TimeSerieOutput - ModuleBasin - ERR130' 
                                      
                endif
                Property => Property%Next
            enddo            

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
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR010'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = Me%OutPut%NextOutPut,           &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR020'

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR030'


            !Writes the Open Points
            call HDF5WriteData   (Me%ObjHDF5, "//Grid/OpenPoints",              &
                                  "OpenPoints", "-",                            &
                                  Array2D = Me%ExtVar%OpenPoints2D,             &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR040'

            !Writes the Acc. Rain
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccRainFall",          &
                                  "Acc. Rainfall", "m",                         &
                                  Array2D      = Me%AccRainFall,                &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR070'

            !Writes the Acc Infil
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccInfiltration",      &
                                  "Acc. Infiltration", "m",                     &
                                  Array2D      = Me%AccInfiltration,            &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR080'

            !Writes the Acc Flow Production
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccFlowProduction",    &
                                  "Acc. Flow Production", "m",                  &
                                  Array2D      = Me%AccFlowProduction,          &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR090'

            !Writes the Acc EVTP
            call HDF5WriteData   (Me%ObjHDF5, "//Results/AccEVTP",              &
                                  "Acc. Evapotrans", "m",                       &
                                  Array2D      = Me%AccEVTP,                    &
                                  OutputNumber = Me%OutPut%NextOutPut,          &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR100'


            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'HDF5Output - ModuleBasin - ERR110'

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "HDF5Output")

    end subroutine HDF5Output

    !--------------------------------------------------------------------------

    subroutine EVTPHDFOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                           :: ILB, IUB, JLB, JUB    
        integer                           :: STAT_CALL, i,j           
        real, dimension(6), target        :: AuxTime
        real(4), dimension(6), target     :: AuxTimeSingle
        real(4), dimension(:), pointer       :: TimePointer
        type (T_BasinProperty), pointer   :: RefEvapotrans
        real(8), dimension(:,:), pointer  :: aux
        real, dimension(:,:), pointer     :: ActualTranspiration
        real(8), dimension(:,:), pointer  :: ActualEvaporation
        real(8)                           :: DTFactor

        if (MonitorPerformance) call StartWatch ("ModuleBasin", "EVTPHDFOutput")
        
        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB           
        call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2), &
                                            AuxTime(3), AuxTime(4), &
                                            AuxTime(5), AuxTime(6))
        AuxTimeSingle = AuxTime
        TimePointer => AuxTimeSingle        

        if (Me%EVTPOutput%Yes) then
            if (Me%CurrentTime >= Me%EVTPOutPut%OutTime(Me%EVTPOutPut%NextOutPut)) then
                call HDF5SetLimits  (Me%ObjEVTPHDF, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR010'

                !Writes current time                
                call HDF5WriteData  (Me%ObjEVTPHDF, "/Time", "Time",          &
                                     "YYYY/MM/DD HH:MM:SS",                   &
                                     Array1D      = TimePointer,              &
                                     OutputNumber = Me%EVTPOutPut%NextOutPut, &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR020'

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjEVTPHDF, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR030'


                !Writes the Open Points
                call HDF5WriteData   (Me%ObjEVTPHDF, "//Grid/OpenPoints",      &
                                      "OpenPoints", "-",                       &
                                      Array2D = Me%ExtVar%OpenPoints2D,        &
                                      OutputNumber = Me%EVTPOutPut%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR040'
                
                !Writes results
                call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccRefEVTP",   &
                                      "AccRefEVTP", "mm",                      &
                                      Array2D      = Me%PartialAccEVTPRef,     &
                                      OutputNumber = Me%EVTPOutPut%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR050' 
                Me%PartialAccEVTPRef = 0.0
                
                if (Me%Coupled%Vegetation) then
                    call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccCropEVTP",  &
                                          "AccCropEVTP", "mm",                     &
                                          Array2D      = Me%PartialAccEVTPCrop,    &
                                          OutputNumber = Me%EVTPOutPut%NextOutPut, &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR060'                  
                    Me%PartialAccEVTPCrop = 0.0
                endif                
                    
                call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccEVTP",      &
                                      "AccEVTP", "mm",                         &
                                      Array2D      = Me%PartialAccEVTP,        &
                                      OutputNumber = Me%EVTPOutPut%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR070'
                Me%PartialAccEVTP = 0.0
                
                if (Me%EvapotranspirationMethod == SeparateEvapoTranspiration) then                
                    call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccPotEV", &
                                          "AccPotEV", "mm",                    &
                                          Array2D      = Me%PartialAccEVPot,           &
                                          OutputNumber = Me%EVTPOutPut%NextOutPut,     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR080' 
                    Me%PartialAccEVPot = 0.0
                    
                    call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccPotTP", &
                                          "AccPotTP", "mm",                    &
                                          Array2D      = Me%PartialAccETPot,             &
                                          OutputNumber = Me%EVTPOutPut%NextOutPut,       &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR090'
                    Me%PartialAccETPot = 0.0
                    
                    call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccActualEV",  &
                                          "AccActualEV", "mm",                     &
                                          Array2D      = Me%PartialAccEVAct,         &
                                          OutputNumber = Me%EVTPOutPut%NextOutPut,   &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR100'      
                    Me%PartialAccEVAct = 0.0 
                    
                    call HDF5WriteData   (Me%ObjEVTPHDF, "//Results/AccActualTP", &
                                          "AccActualTP", "mm",                    &
                                          Array2D      = Me%PartialAccETAct,          &
                                          OutputNumber = Me%EVTPOutPut%NextOutPut,    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR110'  
                    Me%PartialAccETAct = 0.0
                endif
                
                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjEVTPHDF, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR120'
                                                
                Me%EVTPOutPut%NextOutPut = Me%EVTPOutPut%NextOutPut + 1
            endif
        endif
        
        if (Me%EVTPOutput2%Yes) then
            if (Me%CurrentTime >= Me%EVTPOutPut2%OutTime(Me%EVTPOutPut2%NextOutPut)) then

                call HDF5SetLimits  (Me%ObjEVTPHDF2, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR011'

                !Writes current time
                call HDF5WriteData  (Me%ObjEVTPHDF2, "/Time", "Time",          &
                                     "YYYY/MM/DD HH:MM:SS",                    &
                                     Array1D      = TimePointer,               &
                                     OutputNumber = Me%EVTPOutPut2%NextOutPut, &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR021'

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjEVTPHDF2, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR031'


                !Writes the Open Points
                call HDF5WriteData   (Me%ObjEVTPHDF2, "//Grid/OpenPoints",      &
                                      "OpenPoints", "-",                        &
                                      Array2D = Me%ExtVar%OpenPoints2D,         &
                                      OutputNumber = Me%EVTPOutPut2%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR041'

                !Writes results
                call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccRefEVTP",   &
                                      "AccRefEVTP", "mm",                       &
                                      Array2D      = Me%PartialAccEVTPRef2,     &
                                      OutputNumber = Me%EVTPOutPut2%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR051' 
                Me%PartialAccEVTPRef2 = 0.0
                
                if (Me%Coupled%Vegetation) then
                    call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccCropEVTP",  &
                                          "AccCropEVTP", "mm",                      &
                                          Array2D      = Me%PartialAccEVTPCrop2,    &
                                          OutputNumber = Me%EVTPOutPut2%NextOutPut, &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR061'                  
                    Me%PartialAccEVTPCrop2 = 0.0
                endif                
                    
                call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccEVTP",      &
                                      "AccEVTP", "mm",                          &
                                      Array2D      = Me%PartialAccEVTP2,        &
                                      OutputNumber = Me%EVTPOutPut2%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR071'
                Me%PartialAccEVTP2 = 0.0
                
                if (Me%EvapotranspirationMethod == SeparateEvapoTranspiration) then                
                    call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccPotEV", &
                                          "AccPotEV", "mm",                     &
                                          Array2D      = Me%PartialAccEVPot2,           &
                                          OutputNumber = Me%EVTPOutPut2%NextOutPut,     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR081' 
                    Me%PartialAccEVPot2 = 0.0
                    
                    call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccPotTP", &
                                          "AccPotTP", "mm",                     &
                                          Array2D      = Me%PartialAccETPot2,             &
                                          OutputNumber = Me%EVTPOutPut2%NextOutPut,       &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR091'
                    Me%PartialAccETPot2 = 0.0
                    
                    call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccActualEV",  &
                                          "AccActualEV", "mm",                      &
                                          Array2D      = Me%PartialAccEVAct2,         &
                                          OutputNumber = Me%EVTPOutPut2%NextOutPut,   &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR101'      
                    Me%PartialAccEVAct2 = 0.0 
                    
                    call HDF5WriteData   (Me%ObjEVTPHDF2, "//Results/AccActualTP", &
                                          "AccActualTP", "mm",                     &
                                          Array2D      = Me%PartialAccETAct2,          &
                                          OutputNumber = Me%EVTPOutPut2%NextOutPut,    &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR111'  
                    Me%PartialAccETAct2 = 0.0
                endif
      
                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjEVTPHDF2, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR121'

                Me%EVTPOutPut2%NextOutPut = Me%EVTPOutPut2%NextOutPut + 1

            endif
        endif
    
        if (Me%EVTPInstOutput%Yes) then
            if (Me%CurrentTime >= Me%EVTPInstOutPut%OutTime(Me%EVTPInstOutPut%NextOutPut)) then        
                call HDF5SetLimits  (Me%ObjEVTPInstHDF, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR012'

                !Writes current time
                call HDF5WriteData  (Me%ObjEVTPInstHDF, "/Time", "Time",          &
                                     "YYYY/MM/DD HH:MM:SS",                       &
                                     Array1D      = TimePointer,                  &
                                     OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR022'

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjEVTPInstHDF, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR032'


                !Writes the Open Points
                call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Grid/OpenPoints",      &
                                      "OpenPoints", "-",                           &
                                      Array2D = Me%ExtVar%OpenPoints2D,            &
                                      OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR042'
        
                !m3/s
                call GetTranspiration(Me%ObjVegetation, ActualTranspiration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR043'
            
                !m
                call GetEvaporation(Me%ObjPorousMedia, ActualEvaporation, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR044'  
            
                !m/s
                call SearchProperty(RefEvapotrans, RefEvapotrans_, .true., STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR045'  
                
                DTFactor = 1000 * Me%CurrentDT
                
                allocate(aux(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                aux = null_real
                
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB   
                    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                        aux (i, j) = RefEvapotrans%Field(i,j) * 1000. * 3600.
                    endif
                enddo
                enddo                
                call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/RefEVTP",        &
                                        "RefEVTP", "mm/h",                           &
                                        Array2D      = aux,                          &
                                        OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR052'                 
                
                if (Me%Coupled%Vegetation) then
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB         
                        if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                            aux (i, j) = Me%CropEvapotrans(i,j) * 1000. * 3600.
                        endif
                    enddo
                    enddo                
                    call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/CropEVTP",       &
                                            "CropEVTP", "mm/h",                          &
                                            Array2D      = aux,                          &
                                            OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR062'                                          
                endif
                
                call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/ActualEVTP",     &
                                        "ActualEVTP", "mm/h",                        &
                                        Array2D      = Me%EVTPRate,                  &
                                        OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR072'                 
                
                if (Me%EvapotranspirationMethod == SeparateEvapoTranspiration) then
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB                                                            
                        if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                            aux (i, j) = Me%PotentialEvaporation(i,j) * 1000. * 3600.
                        endif
                    enddo
                    enddo
                    call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/PotEV",        &
                                          "PotEV", "mm/h",                             &
                                          Array2D      = aux,                          &
                                          OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR082' 
                
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB                                                            
                        if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                            aux (i, j) = Me%PotentialTranspiration(i,j) * 1000. * 3600.
                        endif
                    enddo
                    enddo
                    call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/PotTP",        &
                                          "PotTP", "mm/h",                             &
                                          Array2D      = aux,                          &
                                          OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR092'                     

                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB                                                            
                        if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                            aux (i, j) = ActualEvaporation(i,j) / Me%CurrentDT * 1000.0 * 3600.0
                        endif
                    enddo
                    enddo
                    call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/ActualEV",     &
                                          "ActualEV", "mm/h",                          &
                                          Array2D      = aux,                          &
                                          OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR102'                     
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB                                                            
                        if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                            aux (i, j) = ActualTranspiration(i,j) / Me%ExtVar%GridCellArea(i,j) * 1000. * 3600.
                        endif
                    enddo
                    enddo     
                    call HDF5WriteData   (Me%ObjEVTPInstHDF, "//Results/ActualTP",     &
                                          "ActualTP", "mm/h",                          &
                                          Array2D      = aux,                          &
                                          OutputNumber = Me%EVTPInstOutPut%NextOutPut, &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR112'                     
                endif
                deallocate(aux)
                
                call UnGetVegetation(Me%ObjVegetation, ActualTranspiration, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR0113'

                call UnGetPorousMedia(Me%ObjPorousMedia, ActualEvaporation, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR0114'   
                
                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjEVTPInstHDF, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'EVTPHDFOutput - ModuleBasin - ERR122'

                Me%EVTPInstOutPut%NextOutPut = Me%EVTPInstOutPut%NextOutPut + 1                
            endif
        endif

        if (MonitorPerformance) call StopWatch ("ModuleBasin", "EVTPHDFOutput")
            
    end subroutine EVTPHDFOutput

    !--------------------------------------------------------------------------

    !subroutine EVTPRefHDFOutput
    !
    !    !Arguments-------------------------------------------------------------
    !
    !    !Local-----------------------------------------------------------------
    !    integer                                         :: ILB, IUB, JLB, JUB    
    !    integer                                         :: STAT_CALL           
    !    real, dimension(6), target                      :: AuxTime
    !    real, dimension(:), pointer                     :: TimePointer
    !
    !    if (MonitorPerformance) call StartWatch ("ModuleBasin", "EVTPRefHDFOutput")
    !
    !    !Bounds
    !    ILB = Me%WorkSize%ILB
    !    IUB = Me%WorkSize%IUB
    !
    !    JLB = Me%WorkSize%JLB
    !    JUB = Me%WorkSize%JUB   
    !
    !    if (Me%CurrentTime >= Me%EVTPRefOutPut%OutTime(Me%EVTPRefOutPut%NextOutPut)) then
    !
    !        !Writes current time
    !        call ExtractDate   (Me%CurrentTime, AuxTime(1), AuxTime(2), &
    !                                            AuxTime(3), AuxTime(4), &
    !                                            AuxTime(5), AuxTime(6))
    !        TimePointer => AuxTime
    !
    !        call HDF5SetLimits  (Me%ObjEVTPRefHDF, 1, 6, STAT = STAT_CALL)
    !        if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR010'
    !
    !        call HDF5WriteData  (Me%ObjEVTPRefHDF, "/Time", "Time",          &
    !                             "YYYY/MM/DD HH:MM:SS",                      &
    !                             Array1D      = TimePointer,                 &
    !                             OutputNumber = Me%EVTPRefOutPut%NextOutPut, &
    !                             STAT = STAT_CALL)
    !        if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR020'
    !
    !        !Sets limits for next write operations
    !        call HDF5SetLimits   (Me%ObjEVTPRefHDF, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
    !        if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR030'
    !
    !
    !        !Writes the Open Points
    !        call HDF5WriteData   (Me%ObjEVTPRefHDF, "//Grid/OpenPoints",      &
    !                              "OpenPoints", "-",                          &
    !                              Array2D = Me%ExtVar%OpenPoints2D,           &
    !                              OutputNumber = Me%EVTPRefOutPut%NextOutPut, &
    !                              STAT = STAT_CALL)
    !        if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR040'
    !
    !       
    !        !Writes the Water Column - should be on runoff
    !        call HDF5WriteData   (Me%ObjEVTPRefHDF, "//Results/ReferencePartialAccEVTP", &
    !                              "ReferencePartialAccEVTP", "m",                        &
    !                              Array2D      = Me%PartialAccEVTPRef,                   &
    !                              OutputNumber = Me%EVTPRefOutPut%NextOutPut,            &
    !                              STAT = STAT_CALL)
    !        if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR050'
    !  
    !        !Writes everything to disk
    !        call HDF5FlushMemory (Me%ObjEVTPRefHDF, STAT = STAT_CALL)
    !        if (STAT_CALL /= SUCCESS_) stop 'EVTPRefHDFOutput - ModuleBasin - ERR100'
    !
    !        Me%PartialAccEVTPRef = 0.0
    !        Me%EVTPRefOutPut%NextOutPut = Me%EVTPRefOutPut%NextOutPut + 1
    !
    !    endif
    !
    !    if (MonitorPerformance) call StopWatch ("ModuleBasin", "EVTPRefHDFOutput")
    !        
    !end subroutine EVTPRefHDFOutput

    !--------------------------------------------------------------------------
    
    
    subroutine GlobalMassBalance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL, ColNumber, FluxColNumber
        real, dimension(:,:), pointer               :: OL_FlowAtBoundary
        real(8)                                     :: InitialMass, FinalMass
        real(8)                                     :: Sinks, Sources, MassError
        real(8)                                     :: ME_Runoff, ME_Porous, ME_Drainage, ME_Veg
        type (T_BasinProperty), pointer             :: PropertyX
        !Begin-----------------------------------------------------------------
        
        !Gets Mass Balance From Drainage Network
        if (Me%Coupled%DrainageNetwork) then
            call GetVolumes(Me%ObjDrainageNetwork,                                      &
                            TotalInputVolume  = Me%MB%DischargesIn,                     &
                            TotalOutputVolume = Me%MB%OutVolumeChannel,                 &
!                            TotalOverTopVolume= Me%MB%OverTopOut,                       &
                            STAT = STAT_CALL)
        endif

       !Gets flow at open boundary of Module Runoff - should be done in module runoff
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
        
        !Misses get to the discharges in runoff.

        !Total
        InitialMass = Me%MB%IniVolumeRunoff  + Me%MB%IniVolumeVegetation     + Me%MB%IniVolumePorousMedia  & 
                      + Me%MB%IniVolumeChannels
        FinalMass   = Me%MB%FinalVolumeRunoff + Me%MB%FinalVolumeVegetation   + Me%MB%FinalVolumePorousMedia &
                      + Me%MB%FinalVolumeChannels
        
        !Missing Losses to deep aquifer from PM - not yet implemented
        Sinks       = Me%MB%EvapFromVegetation + Me%MB%EvapFromGround + Me%MB%EvapFromSoil + Me%MB%OutVolumeChannel   &
                      + Me%MB%OutVolumeOverLand
        
        !Misses discharges in runoff and diffuse discharge in drainage network
        Sources     = Me%MB%TotalRainIn + Me%MB%DischargesIn !Me%MB%RunoffDischarges + Me%MB%DrainageNetDiffuseDisch
        
        if (FinalMass /= 0.0) then
            MassError   = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
        else
            MassError = 0.0
        endif
        
        !DEBUG only
        !BACIA
        InitialMass = Me%MB%IniVolumeRunoff
        FinalMass   = Me%MB%FinalVolumeRunoff
        Sinks       = Me%MB%EvapFromGround + Me%MB%OutVolumeOverLand + Me%MB%OLFlowToRiver                            &
                      + Me%MB%Infiltration
        !Sources     = Me%MB%RainIn + Me%MB%OverTopOut
        !Direct rain and veg input
        Sources     = Me%MB%RainDirect + Me%MB%DrainageFromVeg  !+ Me%MB%RunoffDischarges
        if (FinalMass /= 0.0) then
            ME_Runoff   = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
        else
            ME_Runoff = 0.0
        endif
        
        !VEGETATION 
        InitialMass = Me%MB%IniVolumeVegetation
        FinalMass   = Me%MB%FinalVolumeVegetation
        Sinks       = Me%MB%DrainageFromVeg + Me%MB%EvapFromVegetation
        Sources     = Me%MB%RainInVeg
        if (FinalMass /= 0.0) then
            ME_Veg  = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
        else
            ME_Veg = 0.0
        endif
        
        !POROUSMEDIA
        InitialMass = Me%MB%IniVolumePorousMedia
        FinalMass   = Me%MB%FinalVolumePorousMedia
        Sinks       = Me%MB%EvapFromSoil + Me%MB%GWFlowToRiver
        Sources     = Me%MB%Infiltration
        if (FinalMass /= 0.0) then
            ME_Porous   = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
        else
            ME_Porous   = 0.0
        endif
        
        !DRAINAGE NET - misses diffuse discharge
        InitialMass = Me%MB%IniVolumeChannels
        FinalMass   = Me%MB%FinalVolumeChannels
        !Sinks       = Me%MB%OutVolumeChannel - Me%MB%OverTopOut
        Sinks       = Me%MB%OutVolumeChannel
        Sources     = Me%MB%OLFlowToRiver + Me%MB%GWFlowToRiver + Me%MB%DischargesIn !+ Me%MB%DrainageNetDiffuseDisch
        if (FinalMass /= 0.0) then
            ME_Drainage = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
        else
            ME_Drainage = 0.0
        endif

        !Final Mass Balance
        Me%TimeSeriesBuffer(1)  = Me%MB%IniVolumeRunoff
        Me%TimeSeriesBuffer(2)  = Me%MB%IniVolumeVegetation
        Me%TimeSeriesBuffer(3)  = Me%MB%IniVolumePorousMedia
        Me%TimeSeriesBuffer(4)  = Me%MB%IniVolumeChannels
        Me%TimeSeriesBuffer(5)  = Me%MB%EvapFromVegetation
        Me%TimeSeriesBuffer(6)  = Me%MB%EvapFromGround
        Me%TimeSeriesBuffer(7)  = Me%MB%EvapFromSoil
        Me%TimeSeriesBuffer(8)  = Me%MB%TotalRainIn      !rain total (direct + veg)
        Me%TimeSeriesBuffer(9)  = Me%MB%RainDirect       !rain direct (uncovered)
        Me%TimeSeriesBuffer(10)  = Me%MB%RainInVeg       !rain veg (covered)
        Me%TimeSeriesBuffer(11)  = Me%MB%DrainageFromVeg
        Me%TimeSeriesBuffer(12)  = Me%MB%RainRunoff      !arriving at runoff (direct + drainage)
        Me%TimeSeriesBuffer(13)  = Me%MB%OutVolumeChannel
        Me%TimeSeriesBuffer(14) = Me%MB%OutVolumeOverLand
        Me%TimeSeriesBuffer(15) = Me%MB%FinalVolumeRunoff
        Me%TimeSeriesBuffer(16) = Me%MB%FinalVolumeVegetation
        Me%TimeSeriesBuffer(17) = Me%MB%FinalVolumePorousMedia
        Me%TimeSeriesBuffer(18) = Me%MB%FinalVolumeChannels
        Me%TimeSeriesBuffer(19) = MassError
        Me%TimeSeriesBuffer(20) = ME_Runoff
        Me%TimeSeriesBuffer(21) = ME_Porous
        Me%TimeSeriesBuffer(22) = ME_Drainage
        Me%TimeSeriesBuffer(23) = ME_Veg
        
        !Flux among modules
        Me%TimeSeriesBuffer(24) = Me%MB%Infiltration
        Me%TimeSeriesBuffer(25) = Me%MB%OLFlowToRiver
        Me%TimeSeriesBuffer(26) = Me%MB%GWFlowToRiver
       

        call WriteTimeSerieLine (Me%ObjTimeSerieBasin, Me%TimeSeriesBuffer, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR30'     
        
        
        !Final Mass Balance
        Me%TimeSeriesBuffer4(1)  = Me%MB%EvapFromVegetation / Me%CurrentDT
        Me%TimeSeriesBuffer4(2)  = Me%MB%EvapFromGround    / Me%CurrentDT
        Me%TimeSeriesBuffer4(3)  = Me%MB%EvapFromSoil      / Me%CurrentDT
        Me%TimeSeriesBuffer4(4)  = Me%MB%TotalRainIn       / Me%CurrentDT    !rain total (direct + veg)
        Me%TimeSeriesBuffer4(5)  = Me%MB%RainDirect        / Me%CurrentDT    !rain direct (uncovered)
        Me%TimeSeriesBuffer4(6)  = Me%MB%RainInVeg         / Me%CurrentDT   !rain veg (covered)
        Me%TimeSeriesBuffer4(7)  = Me%MB%DrainageFromVeg   / Me%CurrentDT
        Me%TimeSeriesBuffer4(8)  = Me%MB%RainRunoff        / Me%CurrentDT !arriving at runoff (direct + drainage)
        Me%TimeSeriesBuffer4(9)  = Me%MB%OutVolumeChannel  / Me%CurrentDT
        Me%TimeSeriesBuffer4(10) = Me%MB%OutVolumeOverLand / Me%CurrentDT
        
        !Flux among modules
        Me%TimeSeriesBuffer4(11) = Me%MB%Infiltration   / Me%CurrentDT
        Me%TimeSeriesBuffer4(12) = Me%MB%OLFlowToRiver  / Me%CurrentDT
        Me%TimeSeriesBuffer4(13) = Me%MB%GWFlowToRiver  / Me%CurrentDT
       

        call WriteTimeSerieLine (Me%ObjTimeSerieBasin2, Me%TimeSeriesBuffer4, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR40'     
                
        
        !!!Mass Balance
        if (Me%Coupled%RunoffProperties) then
            ColNumber = 1
            FluxColNumber = 1
            PropertyX => Me%FirstProperty
            do while(associated(PropertyX))  
                
                PropertyX%MB%PMPTranspiredMass     = 0.0
                PropertyX%MB%PMPExchangeMassToDN   = 0.0
                PropertyX%MB%RPExchangeMassToPMP   = 0.0
                PropertyX%MB%RPExchangeMassToDN    = 0.0
                PropertyX%MB%DNDischargeMass       = 0.0
                PropertyX%MB%DNOutFlowMass         = 0.0
    !            PropertyX%MB%DNOverTopMass         = 0.0
                
                if (PropertyX%AdvectionDiffusion) then
                
                    if (Me%Coupled%PorousMediaProperties .and. .not. PropertyX%Particulate) then
                    
                        call GetPMPMassBalance (PorousMediaPropertiesID = Me%ObjPorousMediaProperties,     &
                                                PropertyID              = PropertyX%ID%IDNumber,           &
                                                TranspiredMass          = PropertyX%MB%PMPTranspiredMass,  &
                                                DNExchangeMass          = PropertyX%MB%PMPExchangeMassToDN,&
                                                RPExchangeMass          = PropertyX%MB%RPExchangeMassToPMP,&
                                                STAT                    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR40'     
                        
                    endif
                    if (Me%Coupled%RunoffProperties) then

                        call GetRPMassBalance (RunoffPropertiesID       = Me%ObjRunoffProperties,          &
                                                PropertyID              = PropertyX%ID%IDNumber,           &
                                                DNExchangeMass          = PropertyX%MB%RPExchangeMassToDN, &
                                                STAT                    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR50'     
                    
                    endif
                    if (Me%Coupled%DrainageNetwork) then
                        call GetDNMassBalance (DrainageNetworkID        = Me%ObjDrainageNetwork,           &
                                                PropertyID              = PropertyX%ID%IDNumber,           &
                                                TotalDischargeMass      = PropertyX%MB%DNDischargeMass,    &
                                                TotalOutFlowMass        = PropertyX%MB%DNOutFlowMass,      &
    !                                            TotalOverTopMass        = DNOverTopMass,                   &
                                                STAT                    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR60'  
                    endif                   
                    
                    !Global MASS BALANCE misses discharges and boundary flow in Runoff and diffuse discharge in DrainageNet  
                    InitialMass = PropertyX%MB%InitialRPStoredMass     + PropertyX%MB%InitialVegetationStoredMass    &
                                  + PropertyX%MB%InitialPMPStoredMass  + PropertyX%MB%InitialDNStoredMass
                    
                    FinalMass   = PropertyX%MB%FinalRPStoredMass     + PropertyX%MB%FinalVegetationStoredMass        &
                                  + PropertyX%MB%FinalPMPStoredMass  + PropertyX%MB%FinalDNStoredMass
                    
                    Sinks       = PropertyX%MB%PMPTranspiredMass + PropertyX%MB%DNOutFlowMass !+ PropertyX%MB%RPBoundaryOutFlowMass
                    
                    Sources     = PropertyX%MB%TotalRainMass + PropertyX%MB%DNDischargeMass     
                                  !+ PropertyX%MB%RPDischargeMass + PropertyX%MB%DNDiffuseDischMass
                    
                    if (FinalMass /= 0.0) then
                        MassError   = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
                    else
                        MassError = 0.0
                    endif
                    
                    !VEGETATION
                    InitialMass = PropertyX%MB%InitialVegetationStoredMass   
                    FinalMass   = PropertyX%MB%FinalVegetationStoredMass       
                    Sinks       = PropertyX%MB%VegDrainedMass
                    Sources     = PropertyX%MB%CoveredRainMass
                    if (FinalMass /= 0.0) then
                        ME_Veg  = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
                    else
                        ME_Veg = 0.0
                    endif
                    
                    !RUNOFF - misses boundary flow and discharge
                    InitialMass = PropertyX%MB%InitialRPStoredMass
                    FinalMass   = PropertyX%MB%FinalRPStoredMass
                    Sinks       = PropertyX%MB%RPExchangeMassToDN   + PropertyX%MB%RPExchangeMassToPMP   
                                  !+ PropertyX%MB%RPBoundaryOutFlowMass 
                    Sources     = PropertyX%MB%UncoveredRainMass + PropertyX%MB%VegDrainedMass !+ PropertyX%MB%RPDischargeMass
                    if (FinalMass /= 0.0) then
                        ME_Runoff   = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
                    else
                        ME_Runoff = 0.0
                    endif
                    
                    !POROUSMEDIA
                    InitialMass = PropertyX%MB%InitialPMPStoredMass
                    FinalMass   = PropertyX%MB%FinalPMPStoredMass
                    Sinks       = PropertyX%MB%PMPTranspiredMass + PropertyX%MB%PMPExchangeMassToDN
                    Sources     = PropertyX%MB%RPExchangeMassToPMP
                    if (FinalMass /= 0.0) then
                        ME_Porous   = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
                    else
                        ME_Porous = 0.0
                    endif
                    
                    !DRAINAGE NET - misses DN discharge mass
                    InitialMass = PropertyX%MB%InitialDNStoredMass
                    FinalMass   = PropertyX%MB%FinalDNStoredMass
                    Sinks       = PropertyX%MB%DNOutFlowMass
                    Sources     = PropertyX%MB%RPExchangeMassToDN + PropertyX%MB%PMPExchangeMassToDN       &
                                  + PropertyX%MB%DNDischargeMass !+ PropertyX%MB%DNDiffuseDischMass
                    if (FinalMass /= 0.0) then
                        ME_Drainage = (InitialMass + Sources - Sinks - FinalMass) / FinalMass
                    else
                        ME_drainage = 0.0
                    endif
                    
                    !Final Mass Balance
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%InitialRPStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%InitialVegetationStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%InitialPMPStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%InitialDNStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%PMPTranspiredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%TotalRainMass       !total input
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%UncoveredRainMass   !rain direct
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%CoveredRainMass     !rain leafs
                    ColNumber                       = ColNumber + 1       
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%VegDrainedMass      !leaf leak
                    ColNumber                       = ColNumber + 1   
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%RunoffInputMass     !leaf leak + direct rain
                    ColNumber                       = ColNumber + 1                                                       
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%DNOutFlowMass
                    ColNumber                       = ColNumber + 1 
    !                Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%RPBoundaryOutFlowMass
    !                ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%FinalRPStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%FinalVegetationStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%FinalPMPStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%FinalDNStoredMass
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = MassError
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = ME_Runoff
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = ME_Porous
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = ME_Drainage
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = ME_Veg
                    ColNumber                       = ColNumber + 1
                                        
                    !Flux among modules
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%RPExchangeMassToPMP
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%RPExchangeMassToDN
                    ColNumber                       = ColNumber + 1 
                    Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%PMPExchangeMassToDN
                    ColNumber                       = ColNumber + 1 
                    
                    
                    !Final Mass Balance fluxes
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%PMPTranspiredMass / Me%CurrentDT
                    FluxColNumber                       = FluxColNumber + 1 
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%TotalRainMass / Me%CurrentDT       !total input
                    FluxColNumber                       = FluxColNumber + 1 
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%UncoveredRainMass / Me%CurrentDT   !rain direct
                    FluxColNumber                       = FluxColNumber + 1 
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%CoveredRainMass / Me%CurrentDT     !rain leafs
                    FluxColNumber                       = FluxColNumber + 1       
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%VegDrainedMass / Me%CurrentDT      !leaf leak
                    FluxColNumber                       = FluxColNumber + 1   
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%RunoffInputMass  / Me%CurrentDT    !leaf leak + direct rain
                    FluxColNumber                       = FluxColNumber + 1                                                       
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%DNOutFlowMass / Me%CurrentDT
                    FluxColNumber                       = FluxColNumber + 1 
    !                Me%TimeSeriesBuffer3(ColNumber)  = PropertyX%MB%RPBoundaryOutFlowMass
    !                ColNumber                       = ColNumber + 1 
                                        
                    !Flux among modules
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%RPExchangeMassToPMP / Me%CurrentDT
                    FluxColNumber                       = FluxColNumber + 1 
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%RPExchangeMassToDN / Me%CurrentDT
                    FluxColNumber                       = FluxColNumber + 1 
                    Me%TimeSeriesBuffer5(FluxColNumber)  = PropertyX%MB%PMPExchangeMassToDN / Me%CurrentDT
                    FluxColNumber                       = FluxColNumber + 1 
               
                endif
                PropertyX => PropertyX%Next
            
            end do
            
            call WriteTimeSerieLine (Me%ObjTimeSerieBasinMass, Me%TimeSeriesBuffer3, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR40'                        
  
            call WriteTimeSerieLine (Me%ObjTimeSerieBasinMass2, Me%TimeSeriesBuffer5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR50'                        

        
        endif
        
        
    end subroutine GlobalMassBalance

    !--------------------------------------------------------------------------

    subroutine ComputeBasinWaterBalance

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer :: stat, i, j
        real(8) :: Volume, OldValue
        !real(8) :: sum
        real(8) :: InitialVol, FinalVol
        logical :: LineStored        
        
        !Begin-----------------------------------------------------------------
        
        if (Me%Coupled%RunOff) then
            call GetRunOffBoundaryFlowVolume (Me%ObjRunOff, Volume, STAT = stat)
            if (stat /= SUCCESS_) stop 'ComputeBasinWaterBalance - ModuleBasin - ERR005'
            Me%BWB%BoundaryFromSurface = Me%BWB%BoundaryFromSurface + Volume
        endif
        
        if (Me%Coupled%PorousMedia) then
        
            !Gets Evaporation from soil
            call GetEVTPVolumes(Me%ObjPorousMedia, Evaporation = Volume, STAT = stat)
            if (stat /= SUCCESS_) stop 'ComputeBasinWaterBalance - ModuleBasin - ERR010'
            Me%BWB%EvapFromSoil = Me%BWB%EvapFromSoil + Volume
            
            call GetPMBoundaryFlowVolume (Me%ObjPorousMedia, Volume, STAT = stat)
            if (stat /= SUCCESS_) stop 'ComputeBasinWaterBalance - ModuleBasin - ERR011'
            Me%BWB%BoundaryFromSoil = Me%BWB%BoundaryFromSoil + Volume
        
            if (Me%Coupled%Vegetation) then
            
                !Gets Transpiration
                call GetEVTPVolumes(Me%ObjPorousMedia, Transpiration = Volume, STAT = stat)
                if (stat /= SUCCESS_) stop 'ComputeBasinWaterBalance - ModuleBasin - ERR020'
                Me%BWB%Transpiration = Me%BWB%Transpiration + Volume
                
            endif
        
        endif
        
        if (Me%Coupled%DrainageNetwork) then
        
            !Gets Outlet Flow volume from DN
            call GetVolumes(Me%ObjDrainageNetwork, OutletFlowVolume = Volume, STAT = stat)
            if (stat /= SUCCESS_) stop 'ComputeBasinWaterBalance - ModuleBasin - ERR030'
            Me%BWB%OutletFlowVolume = Me%BWB%OutletFlowVolume + Volume
            
        endif
        
        Me%BWB%StoredInSoil        = Me%BWB%FinStoredInSoil       - Me%BWB%IniStoredInSoil
        Me%BWB%StoredInChannels    = Me%BWB%FinStoredInChannels   - Me%BWB%IniStoredInChannels
        Me%BWB%StoredInLeaves      = Me%BWB%FinStoredInLeaves     - Me%BWB%IniStoredInLeaves
        Me%BWB%StoredInSurface     = Me%BWB%FinStoredInSurface    - Me%BWB%IniStoredInSurface
        Me%BWB%StoredInStormWater  = Me%BWB%FinStoredInStormWater - Me%BWB%IniStoredInStormWater        
        
        Me%BWB%Input      = Me%BWB%Rain + Me%BWB%Irrigation
        Me%BWB%Stored     = Me%BWB%StoredInSoil + Me%BWB%StoredInChannels + Me%BWB%StoredInLeaves +     &
                            Me%BWB%StoredInSurface + Me%BWB%StoredInStormWater
        Me%BWB%Output     = Me%BWB%OutletFlowVolume + Me%BWB%EvapFromSoil + Me%BWB%EvapFromLeaves +     &
                            Me%BWB%EvapFromChannels + Me%BWB%EvapFromSurface + Me%BWB%Transpiration
        Me%BWB%Discharges = Me%BWB%DischargesOnSoil + Me%BWB%DischargesOnSurface +                      &  
                            Me%BWB%DischargesOnChannels    
        Me%BWB%Boundary   = Me%BWB%BoundaryFromSurface + Me%BWB%BoundaryFromSoil
        
        InitialVol = Me%BWB%IniStoredInSoil + Me%BWB%IniStoredInChannels + Me%BWB%IniStoredInLeaves +   &
                     Me%BWB%IniStoredInSurface + Me%BWB%IniStoredInStormWater
        
        FinalVol   = Me%BWB%FinStoredInSoil + Me%BWB%FinStoredInChannels + Me%BWB%FinStoredInLeaves +   &
                     Me%BWB%FinStoredInSurface + Me%BWB%FinStoredInStormWater
        
        Me%BWB%ErrorInVolume     = (InitialVol + Me%BWB%Input + Me%BWB%Discharges) -                    &
                                   (FinalVol + Me%BWB%Output + Me%BWB%Boundary)
        OldValue                 = Me%BWB%AccErrorInVolume
        Me%BWB%AccErrorInVolume  = Me%BWB%AccErrorInVolume + Me%BWB%ErrorInVolume
        Me%BWB%ErrorInPercentage = (Me%BWB%ErrorInVolume / FinalVol) * 100.0
        
        !sum = 0
        !do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        !do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        !    if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then 
        !        !m3 = m3 + mm * m/mm * m2
        !        sum = sum + Me%PartialAccEVTP(i,j) * 0.001 * Me%ExtVar%GridCellArea(i, j)
        !    endif
        !enddo
        !enddo
        
        !Writes to Timeseries
        
        Me%BWBBuffer(1)  = Me%BWB%Rain
        Me%BWBBuffer(2)  = Me%BWB%Irrigation
        Me%BWBBuffer(3)  = Me%BWB%DischargesOnSoil
        Me%BWBBuffer(4)  = Me%BWB%DischargesOnSurface
        Me%BWBBuffer(5)  = Me%BWB%DischargesOnChannels
        Me%BWBBuffer(6)  = Me%BWB%IniStoredInSoil
        Me%BWBBuffer(7)  = Me%BWB%FinStoredInSoil
        Me%BWBBuffer(8)  = Me%BWB%StoredInSoil
        Me%BWBBuffer(9)  = Me%BWB%IniStoredInChannels
        Me%BWBBuffer(10) = Me%BWB%FinStoredInChannels
        Me%BWBBuffer(11) = Me%BWB%StoredInChannels
        Me%BWBBuffer(12) = Me%BWB%IniStoredInStormWater
        Me%BWBBuffer(13) = Me%BWB%FinStoredInStormWater
        Me%BWBBuffer(14) = Me%BWB%StoredInStormWater        
        Me%BWBBuffer(15) = Me%BWB%IniStoredInLeaves
        Me%BWBBuffer(16) = Me%BWB%FinStoredInLeaves
        Me%BWBBuffer(17) = Me%BWB%StoredInLeaves
        Me%BWBBuffer(18) = Me%BWB%IniStoredInSurface
        Me%BWBBuffer(19) = Me%BWB%FinStoredInSurface
        Me%BWBBuffer(20) = Me%BWB%StoredInSurface
        Me%BWBBuffer(21) = Me%BWB%OutletFlowVolume
        Me%BWBBuffer(22) = Me%BWB%EvapFromSoil
        Me%BWBBuffer(23) = Me%BWB%EvapFromLeaves
        Me%BWBBuffer(24) = Me%BWB%EvapFromChannels
        Me%BWBBuffer(25) = Me%BWB%EvapFromSurface
        Me%BWBBuffer(26) = Me%BWB%Transpiration
        Me%BWBBuffer(27) = Me%BWB%BoundaryFromSoil
        Me%BWBBuffer(28) = Me%BWB%BoundaryFromSurface
        Me%BWBBuffer(29) = Me%BWB%Input
        Me%BWBBuffer(30) = Me%BWB%Discharges
        Me%BWBBuffer(31) = Me%BWB%Output
        Me%BWBBuffer(32) = Me%BWB%Stored
        Me%BWBBuffer(33) = Me%BWB%Boundary
        Me%BWBBuffer(34) = Me%BWB%ErrorInVolume 
        Me%BWBBuffer(35) = Me%BWB%AccErrorInVolume
        Me%BWBBuffer(36) = Me%BWB%ErrorInPercentage
        !Me%BWBBuffer(36) = Me%BWB%BasinArea
        !Me%BWBBuffer(37) = Me%BWB%NumberOfCells
        !Me%BWBBuffer(38) = sum
        
        call WriteTimeSerieLine (Me%ObjBWB, Me%BWBBuffer, LineStored = LineStored, STAT = stat)
        if (stat /= SUCCESS_) stop 'GlobalMassBalance - ModuleBasin - ERR40'         
        
        if (LineStored) then
        
            Me%BWB%StoreInitial = .true.
        
            Me%BWB%Rain                  = 0.0
            Me%BWB%Irrigation            = 0.0
            Me%BWB%DischargesOnSoil      = 0.0
            Me%BWB%DischargesOnSurface   = 0.0
            Me%BWB%DischargesOnChannels  = 0.0
            Me%BWB%IniStoredInSoil       = 0.0
            Me%BWB%FinStoredInSoil       = 0.0
            Me%BWB%StoredInSoil          = 0.0
            Me%BWB%IniStoredInChannels   = 0.0
            Me%BWB%FinStoredInChannels   = 0.0
            Me%BWB%StoredInChannels      = 0.0
            Me%BWB%IniStoredInStormWater = 0.0
            Me%BWB%FinStoredInStormWater = 0.0
            Me%BWB%StoredInStormWater    = 0.0
            Me%BWB%IniStoredInLeaves     = 0.0
            Me%BWB%FinStoredInLeaves     = 0.0
            Me%BWB%StoredInLeaves        = 0.0
            Me%BWB%IniStoredInSurface    = 0.0
            Me%BWB%FinStoredInSurface    = 0.0
            Me%BWB%StoredInSurface       = 0.0
            Me%BWB%OutletFlowVolume      = 0.0
            Me%BWB%EvapFromSoil          = 0.0
            Me%BWB%EvapFromLeaves        = 0.0
            Me%BWB%EvapFromChannels      = 0.0
            Me%BWB%EvapFromSurface       = 0.0
            Me%BWB%Transpiration         = 0.0
            Me%BWB%BoundaryFromSoil      = 0.0
            Me%BWB%BoundaryFromSurface   = 0.0
            Me%BWB%Input                 = 0.0
            Me%BWB%Discharges            = 0.0
            Me%BWB%Output                = 0.0
            Me%BWB%Boundary              = 0.0
            Me%BWB%Stored                = 0.0
            Me%BWB%ErrorInVolume         = 0.0
            Me%BWB%ErrorInPercentage     = 0.0            
        
        else
        
            Me%BWB%AccErrorInVolume = OldValue
            
        endif
        
        !----------------------------------------------------------------------        
        
    end subroutine ComputeBasinWaterBalance
    
    !--------------------------------------------------------------------------
    
    subroutine CalculateVegTotalStoredMass

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        type (T_BasinProperty), pointer             :: CurrProperty
        !Begin-----------------------------------------------------------------

        Me%VolumeVegetation = 0.0
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                Me%VolumeVegetation = Me%VolumeVegetation  + Me%CanopyStorage  (i, j) * Me%ExtVar%GridCellArea(i, j) * &
                                      Me%CoveredFraction(i, j)
            endif
        enddo
        enddo
        
        CurrProperty => Me%FirstProperty
        
        do while (associated(CurrProperty)) 

            if (CurrProperty%Inherited) then
                
                CurrProperty%VegTotalStoredMass = 0.0

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                    if (Me%ExtVar%BasinPoints(i, j) == 1) then
                        !kg = kg + (g/m3 * m * m2plant * 1E-3 kg/g)
                        CurrProperty%VegTotalStoredMass = CurrProperty%VegTotalStoredMass                                  &
                                                          + (CurrProperty%VegetationConc(i,j) * Me%CanopyStorage(i,j)      &
                                                          * Me%CoveredFraction(i,j) * Me%ExtVar%GridCellArea(i,j) * 1E-3)
                            
                    endif

                enddo
                enddo
            
            endif
                
            CurrProperty => CurrProperty%Next
        end do 
        


    end subroutine CalculateVegTotalStoredMass

    !--------------------------------------------------------------------------

    subroutine PredictNewDT (NewDT)

        !Arguments---------------------------------------------------------------
        real                                        :: NewDT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL     
        real                                        :: AtmosfereDT, DTForNextEvent
        real                                        :: DNetDT, RunOffDT
        real                                        :: PorousMediaDT, MaxDT
        integer                                     :: ID_DT
        character(len=132)                          :: AuxString   

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
!        if (AtmosfereDT < NewDT) then
!            NewDT = AtmosfereDT
!            ID_DT = 1
!        endif
        
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
        !    NewDT = MaxDT
        !endif
        
!        if (NewDT < 1.0) then
!            call WriteDTLog ('ModuleBasin < 1', ID_DT, NewDT)
!            NewDT = 1.0
!            ID_DT = 7
!        endif

        call GetMaxComputeTimeStep(Me%ObjTime, MaxDT, STAT = STAT_CALL)
        write(AuxString, fmt=10)min(DNetDT, MaxDT), min(RunOffDT, MaxDT), min(PorousMediaDT, MaxDT)
        10 format(f12.4, 1x, f12.4, 1x, f12.4)

        call WriteDTLog_ML ('ModuleBasin', ID_DT, NewDT, AuxString)



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
                if (nUsers == 0)           stop 'KillBasin - ModuleBasin - ERR010'

                if (Me%Coupled%Runoff) then
                    if (Me%Coupled%RunoffProperties) then
                        call KillRunoffProperties (Me%ObjRunoffProperties, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR020'
                    endif                    
                    
                    call KillRunOff         (Me%ObjRunOff,      STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR030'
                    
                endif
                
                if (Me%Coupled%Vegetation) then
                    call KillVegetation  (Me%ObjVegetation, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR040'
                endif

                if (Me%Coupled%PorousMedia) then 
                    if (Me%Coupled%PorousMediaProperties) then
                        call KillPorousMediaProperties (Me%ObjPorousMediaProperties, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR050'
                    endif
                    
                    call KillPorousMedia (Me%ObjPorousMedia, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR60'

                endif
                                
                if (Me%Coupled%DrainageNetwork) then
                    call KillDrainageNetwork(Me%ObjDrainageNetwork, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR070'
                endif


                if (Me%Coupled%Atmosphere) then
                    call KillAtmosphere     (Me%ObjAtmosphere,  STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR080'
                endif                
               
                PropertyX => Me%FirstProperty
                do while(associated(PropertyX))  
                    if (PropertyX%ID%ObjFillMatrix /= 0) then
                        call KillFillMatrix (PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR90'
                    endif
                    PropertyX => PropertyX%Next
                end do

                call KillHorizontalMap      (Me%ObjHorizontalMap,   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR100'

                call KillBasinGeometry      (Me%ObjBasinGeometry,   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR110'

                call KillGridData           (Me%ObjGridData,        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR120'

                call KillHorizontalGrid     (Me%ObjHorizontalGrid,  STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR130'

                if (Me%Output%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR140'
                endif            

                if (Me%EVTPOutput%Yes) then
                    call KillHDF5 (Me%ObjEVTPHDF, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR150'
                endif    
                if (Me%EVTPOutput2%Yes) then
                    call KillHDF5 (Me%ObjEVTPHDF2, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR151'
                endif                   
                if (Me%EVTPInstOutput%Yes) then
                    call KillHDF5 (Me%ObjEVTPInstHDF, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR152'
                endif                                  
                
                !Kills the TimeSerie
                if (Me%ObjTimeSerie > 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR160'
                endif

                if ((Me%ObjTimeSerieBasin > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerieBasin, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR170'
                endif

                if ((Me%ObjTimeSerieBasin2 > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerieBasin2, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR170'
                endif
               
                if ((Me%ObjTimeSerieBasinMass > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerieBasinMass, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR180'
                endif
                
                if (Me%DailyFlow%ObjTimeSerie > 0) then
                    call KillTimeSerie(Me%DailyFlow%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR190'
                endif
                    
                if (Me%MonthlyFlow%ObjTimeSerie > 0) then
                    call KillTimeSerie(Me%MonthlyFlow%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR200'
                endif
                
                if (Me%ComputeBasinWaterBalance) then
                    if (Me%ObjBWB > 0) then
                        call KillTimeSerie(Me%ObjBWB, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillBasin - ModuleBasin - ERR200a'                
                    endif
                    deallocate(Me%BWBBuffer)
                endif
                
#ifdef _ENABLE_CUDA                
                !Kills ModuleCuda
                call KillCuda           (Me%ObjCuda,                STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'KillModel - ModuleModel - ERR200'
#endif _ENABLE_CUDA


!                deallocate (Me%WaterLevel)

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

        !Writes AccErrorInVolume from the Basin Water Balance
        if (Me%ComputeBasinWaterBalance) then
            write(FinalFile)Me%BWB%AccErrorInVolume
        endif
        
!        write(FinalFile)Me%WaterLevel
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

    subroutine ReadLockExternalVar(LockToWhichModules, OptionsType)
        
        !Arguments-------------------------------------------------------------
        character (Len = StringLength)              :: LockToWhichModules
        character (Len = StringLength), optional    :: OptionsType

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        if (LockToWhichModules == 'AllModules') then
            
            !Matrixes from runoff that basin will update and send to runoff
            if (present(OptionsType)) then
                if (OptionsType == 'ModifyBasin') then
                    !Gets water level to control water column updates
                    call GetRunoffWaterLevel(Me%ObjRunOff, Me%ExtUpdate%WaterLevel, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR010'

                    !Get the most recent water col from Runof
                    call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtUpdate%WaterColumn, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR020'    

                    call GetRunoffWaterColumnOld  (Me%ObjRunoff, Me%ExtUpdate%WaterColumnOld, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR030'    
                    
                    call GetMassError   (Me%ObjRunoff, Me%ExtVar%MassError, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR040'    
                    
                endif
            endif

            !Time Stuff
            call GetComputeCurrentTime  (Me%ObjTime, Me%CurrentTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR040'

            call GetComputeTimeStep     (Me%ObjTime, Me%CurrentDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR050'

            !Gets Basin Points
            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR060'

            !Gets River Points
            call GetRiverPoints (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR070'

            !Gets Grid Cell Area
            call GetGridCellArea(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR080'

            !Gets a pointer to Topography
            call GetGridData        (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR090'

            !Gets a pointer to OpenPoints2D
            call GetOpenPoints2D  (Me%ObjHorizontalMap, Me%ExtVar%OpenPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR100'
        
            call GetBoundaries    (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR110'

        elseif (LockToWhichModules == 'Atmosphere') then        
        
            !Leaf Area Index
            call GetLeafAreaIndex(Me%ObjVegetation, Me%ExtVar%LeafAreaIndex, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR120'
            
            if (Me%UsePotLAI) then
                !Potencial Leaf Area Index            
                call GetPotLeafAreaIndex(Me%ObjVegetation, Me%ExtVar%PotLeafAreaIndex, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR130'
            endif
            
            !Specific Storage Capacity
            call GetSpecificLeafStorage(Me%ObjVegetation, Me%ExtVar%SpecificLeafStorage, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR140'
            
            !EVTP CropCoefficient
            call GetEVTPCropCoefficient(Me%ObjVegetation, Me%ExtVar%CropCoefficient, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleBasin - ERR150'
            
        endif
        

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar (UnLockToWhichModules, OptionsType)
        
        !Arguments-------------------------------------------------------------
        character (Len = StringLength)              :: UnLockToWhichModules
        character (Len = StringLength), optional    :: OptionsType

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        if (UnLockToWhichModules == 'AllModules') then
            
            !Matrixes from runoff that basin will update and send to runoff
            if (present(OptionsType)) then
                if (OptionsType == 'ModifyBasin') then
                    call UngetRunoff (Me%ObjRunOff, Me%ExtUpdate%WaterLevel, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR010'

                    call UngetRunoff (Me%ObjRunOff, Me%ExtUpdate%WaterColumn, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR020'  

                    call UngetRunoff (Me%ObjRunOff, Me%ExtUpdate%WaterColumnOld, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR030'  

                    call UngetRunoff (Me%ObjRunOff, Me%ExtVar%MassError, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR040'  
                      
                endif
            endif   

            !UnGets Basin Points
            call UnGetBasin         (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR040'

            !UnGets River Points
            call UnGetBasin         (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR050'

            !UnGets Grid Cell Area
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%GridCellArea, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR060'

            !UnGets Topography
            call UnGetGridData      (Me%ObjGridData, Me%ExtVar%Topography, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR070'

            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%OpenPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR080'
    
            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExtVar%BoundaryPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR090'


        elseif (UnLockToWhichModules == 'Atmosphere') then        

            !Ungets Leaf Area Index
            call UngetVegetation(Me%ObjVegetation, Me%ExtVar%LeafAreaIndex, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR100'
            
            if (Me%UsePotLAI) then           
                !Ungets Potencial Leaf Area Index
                call UngetVegetation(Me%ObjVegetation, Me%ExtVar%PotLeafAreaIndex, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR110'
            endif
                            
            !Specific Storage Capacity
            call UngetVegetation(Me%ObjVegetation, Me%ExtVar%SpecificLeafStorage, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR120'
            
            !EVTP CropCoefficient
            call UngetVegetation(Me%ObjVegetation, Me%ExtVar%CropCoefficient, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleBasin - ERR130'

        endif

    end subroutine ReadUnLockExternalVar

    !--------------------------------------------------------------------------

#ifdef _OPENMI_


    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::IsDrainageNetworkActive
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_ISDRAINAGENETWORKACTIVE"::IsDrainageNetworkActive
    !DEC$ ENDIF
    logical function IsDrainageNetworkActive(BasinID)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: BasinID
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ready_         

        call Ready(BasinID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            IsDrainageNetworkActive = Me%Coupled%DrainageNetwork
        else 
            IsDrainageNetworkActive = .false.
        end if
           
        return

    end function IsDrainageNetworkActive

#endif



end module ModuleBasin




!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








