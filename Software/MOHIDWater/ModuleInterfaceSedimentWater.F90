!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : InterfaceSedimentWater
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute water-sediment interface fluxes
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
!DataFile
!   <begin_manning>
!   See module FillMatrix       : -                -            !Initialization of Manning coefficients
!                                                               !NOTE: mandatory to be included in data file
!                                                               !if not found searches rugosity block
!   <end_manning>

!   <begin_rugosity>
!   See module FillMatrix       : -                [m]          !Initialization of rugosity coefficients
!                                                               !NOTE: mandatory to be included in data file
!   <end_rugosity>

!   <begin_waverugosity>
!   See module FillMatrix       : -                [m]          !Initialization of wave rugosity coefficients
!   <end_waverugosity>

!   <begin_consolidation_rate>
!   See module FillMatrix       : -               [s-1]         !Initialization of consolidation rate
!   <end_consolidation_rate>

!   <begin_critical_shear_deposition>
!   See module FillMatrix       : -               [N/m2]        !Initialization of critical shear stress for deposition
!   <end_critical_shear_deposition>
    
!   <begin_critical_shear_erosion>
!   See module FillMatrix       : -               [N/m2]        !Initialization of critical shear stress for erosion
!   <end_critical_shear_erosion>
!
!   <begin_erosion_rate>
!   See module FillMatrix       : -              [kg/m2s]       !Initialization of erosion rate
!   <end_erosion_rate>


!   CONSOLIDATION               : 0/1              [0]          !Compute consolidation
!   WAVETENSION                 : 0/1              [0]          !Compute wave induced shear stress
!   SAND_TRANSPORT              : 0/1              [0]          !Compute sand tranport
!   SHEAR_STRESS_LIMITATION     : 0/1              [0]          !Limit shear stress values in shallow zones
!       REFERENCE_DEPTH         : real           [0.2 m]        !Reference depth below which shear stress is
!                                                               !is limited
!       REFERENCE_SHEAR_STRESS  : real         [0.25 N/m2]      !Shear stress value assumed in limited cells
!   STATISTICS_SHEAR            : 0/1              [0]          !Perform statistics to shear velocity
!       STATISTICS_SHEAR_FILE   : char              -           !Path to statistics input data file
!   OUTPUT_SHEAR_STRESS         : 0/1              [0]          !Output shear stress in HDF format
!   OUTPUT_TIME                 : sec. sec. sec.    -           !Output Time
!   RESTART_FILE_OUTPUT_TIME    : sec. sec. sec.    -           !Output Time to write restart files
!   RESTART_FILE_OVERWRITE      : 0/1              [1]          !Overwrite intermediate restart files    
!   TIME_SERIE_LOCATION         : char              -           !Path to time serie location file
!   BOXFLUXES                   : char              -           !If specified computes box integration
!                                                               !based on boxes file defined by this keyword

!<beginproperty>
!   NAME                        : char              -           !Property name
!   UNITS                       : char              -           !Property units
!   DESCRIPTION                 : char              -           !Small description of the property
!   See module FillMatrix       : -                 -           !Initialization of concentration values
!                                                               !NOTE: Dissolved properties don´t have mass available
!   MASS_LIMITATION             : 0/1               -           !Property mass is finite
!   MASS_MIN                    : real         [1e-6 kg/m2]     !Minimum mass allowed
!   PARTICULATE                 : 0/1              [0]          !Property physical state: 0 - Dissolved ; 1 - Particulate
!   OLD                         : 0/1              [0]          !Initialization from previous run (overrides FillMatrix)
!   BENTHOS                     : 0/1              [0]          !Compute benthic ecological processes
!   BENTHICECOLOGY              : 0/1              [0]          !Compute benthic ecological processes with module BenthicEcology
!   CEQUALW2                    : 0/1              [0]          !Compute CEQUALW2 benthic ecological processes
!   DETRITUS                    : 0/1              [0]          !Computed as detritus
!   SEDIMENT_WATER_FLUXES       : 0/1              [0]          !Compute fluxes between sediment and water column
!       <<begin_diff_coef>>
!       See module FillMatrix   : -              [m2/s]         !Initialization of diffusion coefficient at the ISW
!       <<end_diff_coef>>
!   SEDIMENT_FLUXES             : 0/1              [0]          !Compute fluxes between ISW and sediment column
!   WATER_FLUXES                : 0/1              [0]          !Compute fluxes to/from water column
!   EROSION                     : 0/1              [0]          !Compute erosion
!   DEPOSITION                  : 0/1              [0]          !Compute deposition
!   DT_INTERVAL                 : real          [ModelDT]       !Property evolution time step (seconds)
!   TIME_SERIE                  : 0/1              [0]          !Ouputs results in time series  
!   BOX_TIME_SERIE              : 0/1              [0]          !Ouputs results in box time series
!   OUTPUT_HDF                  : 0/1              [0]          !Ouputs results in HDF5 format
!<endproperty>
!

!<beginbenthicrate>
!   NAME                        : char              -           !Name of the rate to perform output
!   DESCRIPTION                 : char              -           !Description of the rate to perform output
!   FIRSTPROP                   : char              -           !Name of the first property involved in the rate
!   SECONDPROP                  : char              -           !Name of the second property involved in the rate
!<endbenthicrate>


Module ModuleInterfaceSedimentWater

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleFunctions,            only: ConstructPropertyID, TimeToString, ChangeSuffix,  &
                                          CHUNK_J, CHUNK_I
    use ModuleHDF5,                 only: ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,  &
                                          HDF5WriteData, HDF5FlushMemory, HDF5ReadData, KillHDF5
    use ModuleGridData,             only: GetGridData, UngetGridData
    use ModuleHorizontalGrid,       only: GetHorizontalGrid, GetHorizontalGridSize,         &
                                          WriteHorizontalGrid, UnGetHorizontalGrid,         &
                                          GetGridCellArea, GetXYCellZ
    use ModuleHorizontalMap,        only: GetOpenPoints2D, GetWaterPoints2D, GetBoundaries, &
                                          UnGetHorizontalMap
    use ModuleGeometry,             only: GetGeometrySize, GetGeometryWaterColumn,          &
                                          UnGetGeometry, GetGeometryKFloor,                 &
                                          GetGeometryVolumes, GetGeometryDistances,         &
                                          GetGeometryMinWaterColumn, GetGeometryKtop
    use ModuleMap,                  only: GetWaterPoints3D, GetOpenPoints3D,                &
                                          GetLandPoints3D, UngetMap
    use ModuleBoxDif,               only: StartBoxDif, GetBoxes, GetNumberOfBoxes, BoxDif,  &
                                          UngetBoxDif, KillBoxDif
    use ModuleTimeSerie,            only: StartTimeSerie, WriteTimeSerie, KillTimeSerie,    &
                                          GetTimeSerieLocation, CorrectsCellsTimeSerie,     &
                                          GetNumberOfTimeSeries, TryIgnoreTimeSerie
    use ModuleStatistic,            only: ConstructStatistic, GetStatisticMethod,           &
                                          GetStatisticParameters, ModifyStatistic,          &
                                          KillStatistic
    use ModuleStopWatch,            only: StartWatch, StopWatch
    use ModuleWaterProperties,      only: GetWaterPropertiesSubModulesID, GetConcentration, &
                                          UnGetWaterProperties, GetWaterPropertyOptions,    &
                                          SetFluxToWaterProperties, WaterPropertyExists,    &
                                          GetWaterPropertiesBottomOptions,                  &
                                          SetMacroAlgaeParameters,                          &
                                          GetShortWaveRadiationAverage
    use ModuleHydrodynamic,         only: GetChezy,GetHorizontalVelocity, UnGetHydrodynamic,&
                                          SetHydrodynamicManning, SetBottomWaterFlux,       &
                                          SetHydrodynamicRugosityMatrix, GetWavesStressON,  &
                                          SetWaveChezyVel, SetHydrodynamicChezy
#ifndef _LAGRANGIAN_                                         
#ifdef  _LAGRANGIAN_GLOBAL_                                         
    use ModuleLagrangianGlobal,     only: SetLagrangianShearGlobal
#else
    use ModuleLagrangian,           only: SetLagrangianShear
#endif    
#endif



    use ModuleTurbGOTM,             only: SetTurbGOTMBottomRugosity,                        &
                                          SetTurbGOTMBottomShearVelocity
    use ModuleTurbulence,           only: SetTurbulenceBottomRugosity
    use ModuleFreeVerticalMovement, only: Get_FreeConvFlux, SetDepositionProbability,       &
                                          UngetFreeVerticalMovement, FreeVertPropertyExists,&
                                          FreeVertPropertyHasDeposition
#ifndef _SEDIMENT_ 
    use ModuleSedimentProperties,   only: SedimentPropertyExists,GetSedimentPropertyOptions,&
                                          GetSedimentConcentration, UnGetSedimentProperties,&
                                          SetFluxToSedimentProperties, SetSedimentWaterFlux,&
                                          GetSedimentDryDensity

    use ModuleConsolidation,        only: GetConsolidationWaterFluxes,                      &
                                          GetConsolidationOptions,                          &
                                          GetConsolidationPorosity,                         &
                                          GetConsolidationCriticalShear,                    &
                                          GetConsolidationDrySedVolume,                     &
                                          SetConsolidationFlux,                             &
                                          GetSedimentColumnFull,                            &
                                          GetConsolidationMinThickness, UngetConsolidation 
#endif
    use ModuleInterface,            only: ConstructInterface, SetSOD, Modify_Interface, GetRateFlux,&
                                          KillInterface
#ifndef _WAVES_
    use ModuleWaves,                only: GetWaves, UnGetWaves
#endif
    use ModuleFillMatrix,           only: ConstructFillMatrix, GetIfMatrixRemainsConstant,  &
                                          GetDefaultValue, KillFillMatrix
    use ModuleSand,                 only: StartSand, ModifySand, KillSand,                  &
                                          GetSandDiameters, GetSandDensity, UnGetSand 

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartInterfaceSedimentWater
    private ::      AllocateInstance
    private ::      ReadWaterSedimentFilesName
    private ::      ConstructGlobalVariables
    private ::      ConstructShearStress
    private ::          ConstructShearLimitation
    private ::          ConstructShearStatistics
    private ::          ConstructWaveShearStress
    private ::      ConstructRugosity
    private ::      Construct_PropertyList
    private ::          Construct_Property
    private ::              Construct_PropertyValues
    private ::                  Read_Property_2D
    private ::                  Read_Old_Properties_2D
    private ::              Construct_PropertyEvolution
    private ::              Construct_PropertyOutPut
    private ::          Add_Property
    private ::      Construct_BenthicRateList
    private ::          Construct_BenthicRate
    private ::              Construct_BenthicRateID
    private ::              Construct_BenthicRateValues
    private ::          Add_BenthicRate
    private ::      ConstructConsolidation
    private ::      Construct_Sub_Modules
    private ::          CoupleBenthos
    private ::          CoupleBenthicEcology
    private ::          CoupleCEQUALW2
    private ::          CheckOptionsWaterFluxes
#ifndef _SEDIMENT_
    private ::          CheckOptionsSedimentFluxes
    private ::          CheckOptionsSedimentWaterFluxes
#endif
    private ::          Construct_Time_Serie
    private ::          StartOutputBoxFluxes
    private ::      ConstructSandTransport


    private ::      SetSubModulesConstructor
    
    !Selector
    private :: Search_Property

                     
    !Modifier
    public  :: ModifyInterfaceSedimentWater
    private ::      ModifyShearStress
    private ::          ComputeWaveRugosity
    private ::          ComputeWaveTension
    private ::      Benthos_Processes
    private ::      BenthicEcology_Processes
    private ::      Detritus_Processes
    private ::      CEQUALW2_Processes
    private ::      ModifyWaterColumnFluxes
    private ::          InitializeFluxesToWaterColumn
    private ::          InitializeFluxesToWaterColumn_Benthic
    private ::          ModifyErosionFluxes
    private ::              ModifyErosionCoefficient
    private ::              ErosionFlux             !Function
    private ::              ShearStressLimitation   !Function
    private ::          ModifyDepositionFluxes
    private ::              DepositionProbability   !Function
    private ::          ModifyDissolvedFluxes
    private ::      ModifySedimentColumnFluxes
    private ::          InitializeFluxesToSediment
    private ::          ComputeConsolidation
    private ::          ModifyConsolidatedErosionFluxes
    private ::      ModifySedimentWaterFluxes
    private ::          ComputeWaterFlux
    private ::          DissolvedSedimentWaterFluxes
    private ::          ParticulateSedimentWaterFluxes
    private ::      ModifySandTransport
    private ::      Output_TimeSeries
    private ::      Output_BoxTimeSeries  
    private ::      OutPut_Results_HDF
    private ::      OutPut_Statistics
    private ::      TimeStepActualization
    private ::      Actualize_Time_Evolution
    private ::      OutputRestartFile

    private ::      SetSubModulesModifier
    private ::          SetFluxesToWaterColumn
#ifndef _SEDIMENT_
    private ::          SetFluxesToSedimentColumn
#endif

    !Destructor
    public  :: KillInterfaceSedimentWater                                                     
    private ::      DeAllocateInstance


    !Management
    private ::      Ready
    private ::          LocateObjInterfaceSedimentWater
    
    private ::              ReadLockExternalGlobal
#ifndef _SEDIMENT_
    private ::              ReadLockExternalSediment
    private ::              ReadUnlockExternalSediment
#endif
    private ::              ReadLockExternalWater
    private ::              ReadUnlockExternalGlobal
    private ::              ReadUnlockExternalWater

    
    !Interfaces----------------------------------------------------------------


    !Parameters----------------------------------------------------------------
    character(LEN = StringLength), parameter        :: cse_begin                = '<begin_critical_shear_erosion>'
    character(LEN = StringLength), parameter        :: cse_end                  = '<end_critical_shear_erosion>'
    character(LEN = StringLength), parameter        :: csd_begin                = '<begin_critical_shear_deposition>'
    character(LEN = StringLength), parameter        :: csd_end                  = '<end_critical_shear_deposition>'
    character(LEN = StringLength), parameter        :: erosion_begin            = '<begin_erosion_rate>'
    character(LEN = StringLength), parameter        :: erosion_end              = '<end_erosion_rate>'
    
    character(LEN = StringLength), parameter        :: manning_begin            = '<begin_manning>'
    character(LEN = StringLength), parameter        :: manning_end              = '<end_manning>'
    character(LEN = StringLength), parameter        :: rugosity_begin           = '<begin_rugosity>'
    character(LEN = StringLength), parameter        :: rugosity_end             = '<end_rugosity>'    
    character(LEN = StringLength), parameter        :: SOD_begin                = '<begin_SOD_Rate>'
    character(LEN = StringLength), parameter        :: SOD_end                  = '<end_SOD_Rate>'
    

    character(LEN = StringLength), parameter        :: waverugosity_begin       = '<begin_waverugosity>'
    character(LEN = StringLength), parameter        :: waverugosity_end         = '<end_waverugosity>'

    character(LEN = StringLength), parameter        :: consolidation_begin      = '<begin_consolidation_rate>'
    character(LEN = StringLength), parameter        :: consolidation_end        = '<end_consolidation_rate>'
    
    character(LEN = StringLength), parameter        :: prop_block_begin         = '<beginproperty>'
    character(LEN = StringLength), parameter        :: prop_block_end           = '<endproperty>'
    character(LEN = StringLength), parameter        :: rate_block_begin         = '<beginbenthicrate>'
    character(LEN = StringLength), parameter        :: rate_block_end           = '<endbenthicrate>'

    character(LEN = StringLength), parameter        :: diff_coef_begin          = '<<begin_diff_coef>>'
    character(LEN = StringLength), parameter        :: diff_coef_end            = '<<end_diff_coef>>'    

    !Types---------------------------------------------------------------------

    type     T_ID
        integer                                     :: IDNumber
        character(LEN = StringLength)               :: Name
        character(LEN = StringLength)               :: Description
        character(LEN = StringLength)               :: Units
    end type T_ID

    type       T_Files 
         character(len=PathLength)                  :: InputData
         character(len=PathLength)                  :: Final
         character(len=PathLength)                  :: Initial
         character(len=PathLength)                  :: Results
         character(len=PathLength)                  :: BoxesFile
    end type T_Files

    type       T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime, RestartOutTime
         integer                                    :: NextOutPut, NextRestartOutPut
         logical                                    :: WriteRestartFile     = .false. 
         logical                                    :: Yes
         logical                                    :: WriteFinalFile       = .false.
         logical                                    :: RestartOverwrite
    end type T_OutPut
    
    type       T_Ext_Global
        type(T_Time)                                :: Now
        real,    pointer, dimension(:,:  )          :: XX_IE, YY_IE
        real,    pointer, dimension(:,:  )          :: GridCellArea
    end type T_Ext_Global
    
    type       T_Ext_Water
        real,    pointer, dimension(:,:  )          :: Chezy
        real,    pointer, dimension(:,:,:)          :: DWZ
        real,    pointer, dimension(:,:,:)          :: SZZ
        real,    pointer, dimension(:,:,:)          :: Velocity_U
        real,    pointer, dimension(:,:,:)          :: Velocity_V
        real(8), pointer, dimension(:,:,:)          :: VolumeZ
        real(8), pointer, dimension(:,:,:)          :: VolumeZOld
        integer, pointer, dimension(:,:  )          :: WaterPoints2D
        integer, pointer, dimension(:,:  )          :: OpenPoints2D
        integer, pointer, dimension(:,:  )          :: BoundaryPoints2D
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        integer, pointer, dimension(:,:,:)          :: OpenPoints3D
        integer, pointer, dimension(:,:,:)          :: LandPoints3D
        real,    pointer, dimension(:,:  )          :: Bathymetry
        real,    pointer, dimension(:,:  )          :: WaterColumn
        real,    pointer, dimension(:,:)            :: WaterVolume
        real,    pointer, dimension(:,:)            :: Sediment
        real                                        :: MinWaterColumn        
        integer, pointer, dimension(:,:  )          :: KFloor_Z
        real,    pointer, dimension(:,:  )          :: WavePeriod, WaveHeight
        real,    pointer, dimension(:,:  )          :: Ubw, Abw
        type(T_Time)                                :: LastComputeWave
    end type T_Ext_Water

    type       T_Ext_Sed
        real                                        :: MinLayerThickness
        real,    pointer, dimension(:,:,:)          :: DWZ
        real,    pointer, dimension(:,:,:)          :: SZZ
        real(8), pointer, dimension(:,:,:)          :: VolumeZ
        real(8), pointer, dimension(:,:,:)          :: VolumeZOld
        integer, pointer, dimension(:,:  )          :: WaterPoints2D, KTop
        integer, pointer, dimension(:,:  )          :: OpenPoints2D
        integer, pointer, dimension(:,:  )          :: BoundaryPoints2D
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        integer, pointer, dimension(:,:,:)          :: OpenPoints3D
        integer, pointer, dimension(:,:,:)          :: LandPoints3D
        real,    pointer, dimension(:,:,:)          :: Porosity
        real(8), pointer, dimension(:,:,:)          :: WaterFluxZ, DrySedVolume
        real   , pointer, dimension(:,:  )          :: TopCriticalShear
        real   , pointer, dimension(:,:,:)          :: SedimentDryDensity
        integer, pointer, dimension(:,:  )          :: SedimentColumnFull
        logical                                     :: ComputeConsolidation = .false.
    end type T_Ext_Sed

    type       T_Property_2D
        type(T_PropertyID)                          :: ID
        real, pointer, dimension (:,:)              :: Field 
        real                                        :: Scalar
        logical                                     :: Constant
    end type T_Property_2D

    type       T_Evolution
        logical                                     :: Variable             = .false.
        logical                                     :: Benthos              = .false.
        logical                                     :: BenthicEcology       = .false.
        logical                                     :: BenthicOnly          = .false.
        logical                                     :: CEQUALW2             = .false.
        logical                                     :: Detritus             = .false.
        logical                                     :: SedimentWaterFluxes  = .false.
        logical                                     :: SedimentFluxes       = .false.
        logical                                     :: WaterFluxes          = .false.
        logical                                     :: Erosion              = .false.
        logical                                     :: Deposition           = .false.
        real                                        :: DTInterval           =  null_real
        type(T_Time)                                :: LastCompute
        type(T_Time)                                :: NextCompute
    end type  T_Evolution

    type       T_Property
         type(T_PropertyID)                         :: ID
         logical                                    :: Particulate          = .false.
         real, dimension(:,:), pointer              :: Mass_Available
         real, dimension(:,:), pointer              :: WaterConcentration
         real, dimension(:,:), pointer              :: MassInKg
         real, dimension(:,:), pointer              :: Mass_FromWater
         real, dimension(:,:), pointer              :: MassInKgFromWater
         real, dimension(:,:), pointer              :: SedimentConcentration
         real, dimension(:,:), pointer              :: FluxToWater
         real, dimension(:,:), pointer              :: FluxToSediment
         real, dimension(:,:), pointer              :: ErosionFlux
         real, dimension(:,:), pointer              :: DepositionFlux
         real, dimension(:,:), pointer              :: ErosionCoefficient
         real                                       :: Mass_Min             = FillValueReal
         logical                                    :: Old                  = .false.
         logical                                    :: Mass_Limitation      = .false.
         logical                                    :: TimeSerie            = .false.
         logical                                    :: BoxTimeSerie         = .false.
         logical                                    :: OutputHDF            = .false.
         type(T_Property_2D)                        :: MolecularDifCoef
         type(T_Evolution)                          :: Evolution
         type(T_Property), pointer                  :: Next
         type(T_Property), pointer                  :: Prev
    end type  T_Property 

    type       T_BenthicRate
        type (T_ID)                                 :: ID
        type (T_ID)                                 :: FirstProp
        type (T_ID)                                 :: SecondProp
        real, pointer, dimension(:,:)               :: Field 
        type(T_BenthicRate), pointer                :: Next, Prev
    end type   T_BenthicRate


    type      T_Coupling
         type(T_Time)                               :: NextCompute
         real                                       :: DT_Compute           = FillValueReal
         logical                                    :: Yes                  = .false.
         integer                                    :: NumberOfProperties   = 0
    end type T_Coupling  

    type       T_Coupled
         type(T_Coupling)                           :: Benthos
         type(T_Coupling)                           :: BenthicEcology
         type(T_Coupling)                           :: CEQUALW2
         type(T_Coupling)                           :: Detritus
         type(T_Coupling)                           :: SedimentFluxes   
         type(T_Coupling)                           :: WaterFluxes   
         type(T_Coupling)                           :: SedimentWaterFluxes
         type(T_Coupling)                           :: Erosion
         type(T_Coupling)                           :: Deposition
         type(T_Coupling)                           :: TimeSerie
         type(T_Coupling)                           :: BoxTimeSerie
         type(T_Coupling)                           :: OutputHDF
    end type T_Coupled

    type       T_Statistics
         integer                                    :: ID
         character(LEN = StringLength)              :: File
         logical                                    :: ON
    end type   T_Statistics
    
    type       T_Consolidation
         logical                                    :: Yes                  = .false.
         real, pointer, dimension(:,:)              :: Flux
         type(T_Property_2D)                        :: Rate
         type(T_Time)                               :: LastCompute
    end type  T_Consolidation

    type       T_Shear
         real, pointer, dimension (:,:)             :: CurrentVel, CurrentU, CurrentV
         real, pointer, dimension (:,:)             :: Velocity
         real, pointer, dimension (:,:)             :: Tension
         logical                                    :: Limitation           = .false.      
         real                                       :: ReferenceDepth
         real                                       :: ReferenceShearStress
         type (T_Time)                              :: LastCompute
         type (T_Statistics)                        :: Statistics
         logical                                    :: IntertidalRunOff     = .false.
    end type T_Shear

    type       T_WaveShear
         logical                                    :: Yes                  = .false.
         logical                                    :: NonLinear            = .false.         
         real, pointer, dimension (:,:)             :: Tension
         real, pointer, dimension (:,:)             :: TensionCurrents
         real, pointer, dimension (:,:)             :: ChezyVel
         type(T_Property_2D)                        :: Rugosity
         logical                                    :: RugosityRead
         type(T_Time)                               :: LastCompute
    end type   T_WaveShear

    private :: T_InterfaceSedimentWater
    type       T_InterfaceSedimentWater
        integer                                     :: InstanceID
        character(PathLength)                       :: ModelName
        type(T_Time       )                         :: BeginTime
        type(T_Time       )                         :: EndTime
        type(T_Time       )                         :: ActualTime
        type(T_Size2D     )                         :: Size2D
        type(T_Size2D     )                         :: WorkSize2D
        type(T_Size3D     )                         :: WaterSize3D
        type(T_Size3D     )                         :: WaterWorkSize3D
        type(T_Size3D     )                         :: SedimentSize3D
        type(T_Size3D     )                         :: SedimentWorkSize3D
        type(T_Files      )                         :: Files
        type(T_OutPut     )                         :: Output
        type(T_Coupled    )                         :: Coupled
        type(T_Shear      )                         :: Shear_Stress
        type(T_WaveShear  )                         :: WaveShear_Stress
        type(T_Property_2D)                         :: Critical_Shear_Erosion
        type(T_Property_2D)                         :: Critical_Shear_Deposition
        type(T_Property_2D)                         :: ErosionRate
        type(T_Ext_Global )                         :: ExternalVar
        type(T_Ext_Water  )                         :: ExtWater
        type(T_Ext_Sed    )                         :: ExtSed
        type(T_Property   ), pointer                :: FirstProperty
        type(T_Property   ), pointer                :: LastProperty
        type(T_BenthicRate), pointer                :: FirstBenthicRate
        type(T_BenthicRate), pointer                :: LastBenthicRate
        real,dimension(:,:), pointer                :: DepositionProbability
        logical                                     :: RunsSediments            = .false.
        logical                                     :: RunsSandTransport        = .false.
        logical                                     :: Manning                  = .false.
        logical                                     :: UseSOD                   = .false.
        logical                                     :: MacroAlgae               = .false.
        type(T_Property_2D)                         :: SOD
        type(T_Property_2D)                         :: Rugosity
        type(T_Property_2D)                         :: ManningCoef
        logical                                     :: Chezy
        real                                        :: ChezyCoef
        real(8), pointer, dimension(:,:)            :: Scalar2D
        real(8), pointer, dimension(:,:)            :: WaterFlux
        integer                                     :: PropertiesNumber         = 0
        integer                                     :: BenthicRatesNumber       = 0
        type(T_Consolidation)                       :: Consolidation

        !Instance of ModuleBoxDif                   
        integer                                     :: ObjBoxDif                = 0
                                                                            
        !Instance of ModuleEnterData                                        
        integer                                     :: ObjEnterData             = 0
                                                                            
        !Instance of ModuleTimeSerie                                        
        integer                                     :: ObjTimeSerie             = 0
                                                                                
        !Instance of ModuleHDF5                                                 
        integer                                     :: ObjHDF5                  = 0
                                                                                
        !Instance of ModuleTime                                                 
        integer                                     :: ObjTime                  = 0

        !Instance of ModuleDischarges           
        integer                                     :: ObjDischarges            = 0

        !Instance of ModuleHorizontalGrid                                       
        integer                                     :: ObjHorizontalGrid        = 0
                                                                                
        !Instance of ModuleGridData                                             
        integer                                     :: ObjSedimentGridData      = 0
                                                                                
        !Instance of ModuleGridData                                             
        integer                                     :: ObjWaterGridData         = 0
                                                                                
        !Instance of ModuleGeometry                                             
        integer                                     :: ObjSedimentGeometry      = 0
                                                                                
        !Instance of ModuleGeometry                                             
        integer                                     :: ObjWaterGeometry         = 0
        
        !Instance of ModuleHorizontalMap                                    
        integer                                     :: ObjSedimentHorizontalMap = 0
        
        !Instance of ModuleHorizontalMap                                    
        integer                                     :: ObjWaterHorizontalMap    = 0
                                                                                
        !Instance of ModuleMap                                                  
        integer                                     :: ObjSedimentMap           = 0
                                                                                
        !Instance of ModuleMap                                                  
        integer                                     :: ObjWaterMap              = 0

        !Instance of ModuleSand
        integer                                     :: ObjSand                  = 0
        
        !Instance of ModuleHydrodynamic
        integer                                     :: ObjHydrodynamic          = 0
        
        !Instance of ModuleTurbulence
        integer                                     :: ObjTurbulence            = 0
        
        !Instance of ModuleTurbGOTM
        integer                                     :: ObjTurbGOTM              = 0
        
        !Instance of ModuleWaves
        integer                                     :: ObjWaves                 = 0
        
        !Instance of ModuleInterfaceSedimentWater
        integer                                     :: ObjWaterProperties       = 0
        
        !Instance of ModuleFreeVerticalMovement
        integer                                     :: ObjFreeVerticalMovement  = 0
        
        !Instance of ModuleSedimentProperties
        integer                                     :: ObjSedimentProperties    = 0
        
        !Instance of ModuleConsolidation
        integer                                     :: ObjConsolidation         = 0

        !Instance of ModuleInterface
        integer                                     :: ObjInterface             = 0
        
        !Instance of ModuleLagrangian
        integer                                     :: ObjLagrangian            = 0

        type(T_InterfaceSedimentWater), pointer     :: Next
    end type  T_InterfaceSedimentWater

    !Global Module Variables
    type (T_InterfaceSedimentWater), pointer        :: FirstObjInterfaceSedimentWater
    type (T_InterfaceSedimentWater), pointer        :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartInterfaceSedimentWater( ModelName,                      &    
                                            ObjInterfaceSedimentWaterID,    &
                                            TimeID,                         &
                                            HorizontalGridID,               &
                                            WaterGridDataID,                &
                                            WaterHorizontalMapID,           &
                                            WaterMapID,                     &
                                            WaterGeometryID,                &
                                            SedimentGridDataID,             &
                                            SedimentHorizontalMapID,        &
                                            SedimentMapID,                  &
                                            SedimentGeometryID,             &
                                            HydrodynamicID,                 &
                                            TurbulenceID,                   &
                                            TurbGOTMID,                     &
                                            WavesID,                        &
                                            WaterPropertiesID,              &
                                            LagrangianID,                   &
                                            SedimentPropertiesID,           &
                                            ConsolidationID,                &
                                            DischargesID,                   &
                                            RunsSediments,                  &
                                            STAT)

        !Arguments---------------------------------------------------------------
        character(Len=*)                                :: ModelName
        integer                                         :: ObjInterfaceSedimentWaterID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: WaterGridDataID
        integer                                         :: WaterHorizontalMapID
        integer                                         :: WaterMapID
        integer                                         :: WaterGeometryID
        integer                                         :: SedimentGridDataID
        integer                                         :: SedimentHorizontalMapID
        integer                                         :: SedimentMapID
        integer                                         :: SedimentGeometryID
        integer                                         :: HydrodynamicID
        integer                                         :: LagrangianID    
        integer                                         :: TurbulenceID         
        integer                                         :: TurbGOTMID
        integer                                         :: WavesID
        integer                                         :: WaterPropertiesID     
        integer                                         :: SedimentPropertiesID 
        integer                                         :: ConsolidationID
        integer                                         :: DischargesID
        logical                                         :: RunsSediments     
        integer, optional, intent(OUT)                  :: STAT    

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL       

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterfaceSedimentWater_)) then
            nullify (FirstObjInterfaceSedimentWater)
            call RegisterModule (mInterfaceSedimentWater_) 
        endif

        call Ready(ObjInterfaceSedimentWaterID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            nullify(Me%FirstProperty    )
            nullify(Me%LastProperty     )
            nullify(Me%FirstBenthicRate )
            nullify(Me%LastBenthicRate  )

            !Associates External Instances
            Me%ObjTime                 = AssociateInstance(mTIME_,           TimeID                     )
            Me%ObjHorizontalGrid       = AssociateInstance(mHORIZONTALGRID_, HorizontalGridID           )
            Me%ObjDischarges           = DischargesID


            !Water column
            Me%ObjWaterGridData        = AssociateInstance(mGRIDDATA_,          WaterGridDataID         )
            Me%ObjWaterHorizontalMap   = AssociateInstance(mHORIZONTALMAP_,     WaterHorizontalMapID    )
            Me%ObjWaterGeometry        = AssociateInstance(mGEOMETRY_,          WaterGeometryID         )
            Me%ObjWaterMap             = AssociateInstance(mMAP_,               WaterMapID              )
            Me%ObjHydrodynamic         = AssociateInstance(mHYDRODYNAMIC_,      HydrodynamicID          )
            Me%ObjWaterProperties      = AssociateInstance(mWATERPROPERTIES_,   WaterPropertiesID       )
            Me%ObjTurbulence           = AssociateInstance(mTURBULENCE_,        TurbulenceID            )
            

            if(LagrangianID /= 0)then
                Me%ObjLagrangian       = AssociateInstance(mLAGRANGIAN_,        LagrangianID            )
            end if

            if(TurbGOTMID /= 0)then
                Me%ObjTurbGOTM         = AssociateInstance(mTURBGOTM_,          TurbGOTMID              )
            end if

            if(WavesID /= 0)then
                Me%ObjWaves            = AssociateInstance(mWAVES_,             WavesID                 )
            end if
            
            Me%RunsSediments = RunsSediments

            !Sediment column
            if(Me%RunsSediments)then
                Me%ObjSedimentGridData     = AssociateInstance(mGRIDDATA_,          SedimentGridDataID      )
                Me%ObjSedimentHorizontalMap= AssociateInstance(mHORIZONTALMAP_,     SedimentHorizontalMapID )    
                Me%ObjSedimentGeometry     = AssociateInstance(mGEOMETRY_,          SedimentGeometryID      )
                Me%ObjSedimentMap          = AssociateInstance(mMAP_,               SedimentMapID           )
                Me%ObjConsolidation        = AssociateInstance(mCONSOLIDATION_,     ConsolidationID         )
                Me%ObjSedimentProperties   = AssociateInstance(mSEDIMENTPROPERTIES_,SedimentPropertiesID    )
            end if

            call ReadLockExternalGlobal
            
            call ReadLockExternalWater

#ifndef _SEDIMENT_
            if(Me%RunsSediments) call ReadLockExternalSediment
#endif

            call ConstructGlobalVariables

            call ReadWaterSedimentFilesName

            call ConstructEnterData(Me%ObjEnterData, Me%Files%InPutData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                 &
                stop 'StartInterfaceSedimentWater - InterfaceSedimentWater - ERR02'
            

            call ConstructShearStress

            call ConstructRugosity

            call Construct_PropertyList

            call Construct_BenthicRateList

            call ConstructConsolidation

            call Construct_Sub_Modules
            
            call ConstructSandTransport

            call ConstructLog
            
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                 &
                stop 'StartInterfaceSedimentWater - InterfaceSedimentWater - ERR03'

#ifndef _SEDIMENT_
            if(Me%RunsSediments) call ReadUnlockExternalSediment
#endif

            call ReadUnlockExternalWater

            call ReadUnlockExternalGlobal

            call SetSubModulesConstructor

            !Returns ID
            ObjInterfaceSedimentWaterID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleInterfaceSedimentWater - StartInterfaceSedimentWater - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartInterfaceSedimentWater
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance
                                                    
        !Local-----------------------------------------------------------------
        type (T_InterfaceSedimentWater), pointer    :: NewObjWaterSedInterface
        type (T_InterfaceSedimentWater), pointer    :: PreviousObjWaterSedInterface


        !Allocates new instance
        allocate (NewObjWaterSedInterface)
        nullify  (NewObjWaterSedInterface%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjInterfaceSedimentWater)) then
            FirstObjInterfaceSedimentWater      => NewObjWaterSedInterface
            Me                                  => NewObjWaterSedInterface
        else
            PreviousObjWaterSedInterface        => FirstObjInterfaceSedimentWater
            Me                                  => FirstObjInterfaceSedimentWater%Next
            do while (associated(Me))
                PreviousObjWaterSedInterface    => Me
                Me                              => Me%Next
            enddo
            Me                                  => NewObjWaterSedInterface
            PreviousObjWaterSedInterface%Next   => NewObjWaterSedInterface
        endif

        Me%InstanceID = RegisterNewInstance (mINTERFACESEDIMENTWATER_)


    end subroutine AllocateInstance
    
    
    !--------------------------------------------------------------------------
    
    
    subroutine ReadWaterSedimentFilesName

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len = StringLength)       :: Message

        !----------------------------------------------------------------------

        Message ='ASCII file used to construct new sediment-water properties.'
        Message = trim(Message)

        call ReadFileName('BOT_DAT', Me%Files%InPutData, Message = Message, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'StartInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR01'

        Message   ='Instant fields of bottom properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('BOT_HDF', Me%Files%Results, Message = Message, TIME_END = Me%EndTime, &
                           Extension = 'bot', STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'StartInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR02'

        ! ---> Bottom properties final values in HDF format
        Message   ='Bottom properties final values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('BOT_FIN', Me%Files%Final, Message = Message, TIME_END = Me%EndTime, &
                           Extension = 'bof', STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'StartInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR03'



        ! ---> Bottom properties initial values in HDF format
        Message   ='Bottom properties initial values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('BOT_INI', Me%Files%Initial, Message = Message, TIME_END = Me%ActualTime, STAT = STAT_CALL)
                                                                                               
cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_   ) then
           
            write(*,*)'Initial file not found.'
            stop 'StartInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR04'

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then cd1
            
            Message   = 'Keyword BOT_INI not found - StartInterfaceSedimentWater - ModuleInterfaceSedimentWater'
            Message   = trim(Message)

            call SetError(WARNING_, KEYWORD_, Message, Screen = .false.)
                          
        else if (STAT_CALL .EQ. SUCCESS_              ) then cd1

            continue
        
        else cd1

            stop 'StartInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR06'

        end if cd1  
                                                                             
        !----------------------------------------------------------------------

    end subroutine ReadWaterSedimentFilesName
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructGlobalVariables
                                                    
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                        &
                                   Size        = Me%Size2D,                     &
                                   WorkSize    = Me%WorkSize2D,                 &
                                   STAT        = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructGlobalVariables - ModuleInterfaceSedimentWater - ERR01'

        call GetGeometrySize(Me%ObjWaterGeometry,                               &
                             Size       = Me%WaterSize3D,                       &
                             WorkSize   = Me%WaterWorkSize3D,                   &
                             STAT       = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructGlobalVariables - ModuleInterfaceSedimentWater - ERR02'

        if(Me%RunsSediments)then

            call GetGeometrySize(Me%ObjSedimentGeometry,                        &
                                 Size       = Me%SedimentSize3D,                &
                                 WorkSize   = Me%SedimentWorkSize3D,            &
                                 STAT       = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructGlobalVariables - ModuleInterfaceSedimentWater - ERR03'

        end if

        call GetComputeTimeLimits(Me%ObjTime,                                   &
                                  BeginTime = Me%BeginTime,                     &
                                  EndTime   = Me%EndTime,                       &
                                  STAT      = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructGlobalVariables - ModuleInterfaceSedimentWater - ERR04'


        !Actualize the time
        Me%ActualTime = Me%BeginTime


    end subroutine ConstructGlobalVariables
    
    
    !--------------------------------------------------------------------------

    
    subroutine ConstructShearStress

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, iflag
        logical                             :: OutputShearStress

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        nullify(Me%Shear_Stress%Tension, Me%Shear_Stress%Velocity)
        
        !Shear stress 
        allocate(Me%Shear_Stress%Tension(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructShearStress - ModuleInterfaceSedimentWater - ERR10'
        Me%Shear_Stress%Tension(:,:) = FillValueReal
        
        !Shear velocity 
        allocate(Me%Shear_Stress%Velocity(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructShearStress - ModuleInterfaceSedimentWater - ERR20'
        Me%Shear_Stress%Velocity(:,:) = FillValueReal


        call GetData(OutputShearStress,                                             &     
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     Keyword      = 'OUTPUT_SHEAR_STRESS',                          &
                     ClientModule = 'ModuleInterfaceSedimentWater',                 &
                     Default      = .false.,                                        &
                     STAT         = STAT_CALL)            
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructShearStress - ModuleInterfaceSedimentWater - ERR30'
        
        if(OutputShearStress)then
            Me%Coupled%OutputHDF%Yes = ON
        end if
        
        call ConstructShearLimitation

        call ConstructShearStatistics

        call ConstructWaveShearStress

        call GetData(Me%Shear_Stress%IntertidalRunOff,                              &     
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     Keyword      = 'INTERTIDAL_RUN_OFF_SHEAR',                     &
                     ClientModule = 'ModuleInterfaceSedimentWater',                 &
                     Default      = .false.,                                        &
                     STAT         = STAT_CALL)            
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructShearStress - ModuleInterfaceSedimentWater - ERR40'


    end subroutine ConstructShearStress

    !--------------------------------------------------------------------------

    subroutine ConstructWaveShearStress

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        call GetData(Me%WaveShear_Stress%Yes,                                       &     
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     Keyword      = 'WAVETENSION',                                  &
                     ClientModule = 'ModuleInterfaceSedimentWater',                 &
                     Default      = .false.,                                        &
                     STAT         = STAT_CALL)            
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructWaveShearStress - ModuleInterfaceSedimentWater - ERR01'


        call GetWavesStressON (Me%ObjHydrodynamic,                                  &
                               WavesStressON = Me%WaveShear_Stress%NonLinear, STAT = STAT_CALL)            

        if  (STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructWaveShearStress - ModuleInterfaceSedimentWater - ERR02'

        if(Me%WaveShear_Stress%Yes)then
            if(Me%ObjWaves == 0)then
                write(*,*)
                write(*,*)'Must activate Module Waves in order to'
                write(*,*)'compute wave shear stress'
                stop 'ConstructWaveShearStress - ModuleInterfaceSedimentWater - ERR03'
            end if
            !Shear stress 
            allocate(Me%WaveShear_Stress%Tension   (ILB:IUB, JLB:JUB)) 
            Me%WaveShear_Stress%Tension   (:,:) = FillValueReal

            allocate(Me%WaveShear_Stress%ChezyVel  (ILB:IUB, JLB:JUB)) 
            Me%WaveShear_Stress%ChezyVel  (:,:) = FillValueReal

        end if


    end subroutine ConstructWaveShearStress

    !--------------------------------------------------------------------------
    
    subroutine ReadSOD

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------

        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: ClientNumber
        logical                             :: BlockFound

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB
        
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,          &
                                    SOD_begin, SOD_end, BlockFound,         &
                                    STAT = STAT_CALL)
        if (STAT_CALL  /= SUCCESS_) stop 'ReadSOD - ModuleInterfaceSedimentWater - ERR01'
            
        if(BlockFound)then
            
            Me%UseSOD = ON

            nullify(Me%SOD%Field)
            allocate(Me%SOD%Field(ILB:IUB, JLB:JUB))

            call ConstructFillMatrix  (PropertyID           = Me%SOD%ID,                        &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExtWater%WaterPoints2D,        &
                                       Matrix2D             = Me%SOD%Field,                     &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadSOD - ModuleInterfaceSedimentWater - ERR02'


            call GetDefaultValue(Me%SOD%ID%ObjFillMatrix, Me%SOD%Scalar, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadSOD - ModuleInterfaceSedimentWater - ERR03'

            call KillFillMatrix(Me%SOD%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ReadSOD - ModuleInterfaceSedimentWater - ERR04'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ReadSOD - ModuleInterfaceSedimentWater - ERR05'        
        
        endif
    
        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ReadSOD - ModuleInterfaceSedimentWater - ERR06'
        
        
    end subroutine ReadSOD
    
    !--------------------------------------------------------------------------

    subroutine ConstructSandTransport

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------

        integer                             :: ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        call GetData(Me%RunsSandTransport,                                               &     
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromFile,                                            &
                     Keyword      = 'SAND_TRANSPORT',                                    &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)            
        if(STAT_CALL .ne. SUCCESS_)                                                      &
            stop 'ConstructSandTransport - ModuleInterfaceSedimentWater - ERR10.'

        if (Me%RunsSandTransport) then

            if (.not. associated(Me%Rugosity%Field))then

                write(*,*) 'When SAND TRANSPORT is ON is necessary to defined the RUGOSITY'
                stop 'ConstructSandTransport - ModuleInterfaceSedimentWater - ERR20.'

            endif

            call StartSand(ObjSandID            = Me%ObjSand,                            &
                           ObjGridDataID        = Me%ObjWaterGridData,                   &
                           ObjHorizontalGridID  = Me%ObjHorizontalGrid,                  &
                           ObjHorizontalMapID   = Me%ObjWaterHorizontalMap,              &
                           ObjTimeID            = Me%ObjTime,                            &
                           ObjWavesID           = Me%ObjWaves,                           &
                           ObjDischargesID      = Me%ObjDischarges,                      &
                           WaterDensity         = SigmaDensityReference,          &
                           WaveTensionON        = Me%WaveShear_Stress%Yes,               &          
                           STAT                 = STAT_CALL)

            if(STAT_CALL /= SUCCESS_)                                                    &
                stop 'ConstructSandTransport - ModuleInterfaceSedimentWater - ERR30.'



            if(Me%WaveShear_Stress%Yes)then
                allocate(Me%WaveShear_Stress%TensionCurrents(ILB:IUB, JLB:JUB)) 
                Me%WaveShear_Stress%TensionCurrents(:,:) = FillValueReal
            endif

            !Current velocity 
            allocate(Me%Shear_Stress%CurrentVel(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructSandTransport - ModuleInterfaceSedimentWater - ERR40'
            Me%Shear_Stress%CurrentVel(:,:) = FillValueReal

            allocate(Me%Shear_Stress%CurrentU(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructSandTransport - ModuleInterfaceSedimentWater - ERR50'
            Me%Shear_Stress%CurrentU(:,:) = FillValueReal

            allocate(Me%Shear_Stress%CurrentV(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructSandTransport - ModuleInterfaceSedimentWater - ERR60'
            Me%Shear_Stress%CurrentV(:,:) = FillValueReal

        endif

    end subroutine ConstructSandTransport

    !--------------------------------------------------------------------------

    subroutine ConstructConsolidation

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag
        integer                             :: ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        call GetData(Me%Consolidation%Yes,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     Keyword      = 'CONSOLIDATION',                                    &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)            
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructConsolidation - ModuleInterfaceSedimentWater - ERR10'


        if(Me%Consolidation%Yes)then

            Me%Coupled%SedimentFluxes%Yes      = ON  !review this
            Me%Coupled%SedimentWaterFluxes%Yes = ON  !review this

            call Read_Property_2D (Me%Consolidation%Rate, FromBlock, &
                                   consolidation_begin, consolidation_end)

            allocate(Me%Consolidation%Flux(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructConsolidation - ModuleInterfaceSedimentWater - ERR20'
            Me%Consolidation%Flux(:,:) = 0.

#ifndef _SEDIMENT_

            call GetConsolidationMinThickness(Me%ObjConsolidation, Me%ExtSed%MinLayerThickness, &
                                              STAT = STAT_CALL)            
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructConsolidation - ModuleInterfaceSedimentWater - ERR30'
#else 
            write(*,*)"Cannot compute consolidation because the model was"
            write(*,*)"not compiled with the SEDIMENT option"
            stop 'ConstructConsolidation - ModuleInterfaceSedimentWater - ERR40'

#endif
        end if

    end subroutine ConstructConsolidation
    
    !--------------------------------------------------------------------------


    subroutine ConstructShearLimitation

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        !Keyword so shear stress is limited in low depth waters
        call GetData(Me%Shear_Stress%Limitation,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     Keyword      = 'SHEAR_STRESS_LIMITATION',                          &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)            
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructShearLimitation - ModuleInterfaceSedimentWater - ERR01'
            
        if(Me%Shear_Stress%Limitation)then

            !Reference depth to compute shear stress limitation [m]
            call GetData(Me%Shear_Stress%ReferenceDepth,                                &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         Keyword      = 'REFERENCE_DEPTH',                              &
                         ClientModule = 'ModuleInterfaceSedimentWater',                 &
                         Default      = 0.2,                                            &
                         STAT         = STAT_CALL)            
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructShearLimitation - ModuleInterfaceSedimentWater - ERR02'

            !Reference shear stress in low low depth waters [N/m2]
            call GetData(Me%Shear_Stress%ReferenceShearStress,                          &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         Keyword      = 'REFERENCE_SHEAR_STRESS',                       &
                         ClientModule = 'ModuleInterfaceSedimentWater',                 &
                         Default      = 0.25,                                           &
                         STAT         = STAT_CALL)            
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructShearLimitation - ModuleInterfaceSedimentWater - ERR03'

        end if

    end subroutine ConstructShearLimitation

    !--------------------------------------------------------------------------
    
    subroutine ConstructShearStatistics

        !Local ------------------------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag

        !Begin----------------------------------------------------------------------------
        
        !Checks out if the user pretends the statistics of bottom shear velocity
        call GetData(Me%Shear_Stress%Statistics%ON,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'STATISTICS_SHEAR',                               &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleInterfaceSedimentWater',                   &
                     STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'ConstructShearStatistics - ModuleInterfaceSedimentWater - ERR01'

        
        if (Me%Shear_Stress%Statistics%ON) then

            call GetData(Me%Shear_Stress%Statistics%File,                               &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword        = 'STATISTICS_SHEAR_FILE',                      &
                         SearchType     = FromFile,                                     &
                         ClientModule   = 'ModuleInterfaceSedimentWater',               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag /= 1)                                  &
                stop 'ConstructShearStatistics - ModuleInterfaceSedimentWater - ERR02'


            call ConstructStatistic (StatisticID      = Me%Shear_Stress%Statistics%ID,  &
                                     ObjTime          = Me%ObjTime,                     &
                                     ObjHDF5          = Me%ObjHDF5,                     &
                                     DataFile         = Me%Shear_Stress%Statistics%File,&
                                     Size             = Me%WaterSize3D,                 &
                                     WorkSize         = Me%WaterWorkSize3D,             &
                                     Name             = GetPropertyName (ShearVelocity_),&
                                     STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructShearStatistics - ModuleInterfaceSedimentWater - ERR04'

        endif


    end subroutine ConstructShearStatistics
    

    !--------------------------------------------------------------------------

    subroutine ConstructRugosity

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: WILB, WIUB, WJLB, WJUB
        integer                             :: ClientNumber, iflag
        logical                             :: BlockFound

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        WILB = Me%WorkSize2D%ILB
        WIUB = Me%WorkSize2D%IUB
        WJLB = Me%WorkSize2D%JLB
        WJUB = Me%WorkSize2D%JUB

        Me%Manning = OFF
        Me%Chezy   = OFF

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                    manning_begin, manning_end, BlockFound,         &
                                    STAT = STAT_CALL)
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR01'
            
        if(BlockFound)then
            

            if(Me%WaterWorkSize3D%KUB > 1)then
                
                write(*,*)'Cannot use Manning coefficient in 3D.'
                stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR02'
            end if
            
            Me%Manning = ON

            nullify(Me%ManningCoef%Field)
            allocate(Me%ManningCoef%Field(ILB:IUB, JLB:JUB))

            call ConstructFillMatrix  (PropertyID           = Me%ManningCoef%ID,                &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExtWater%WaterPoints2D,        &
                                       Matrix2D             = Me%ManningCoef%Field,             &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR03'


            call GetDefaultValue(Me%ManningCoef%ID%ObjFillMatrix, Me%ManningCoef%Scalar, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR04'

            call KillFillMatrix(Me%ManningCoef%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR05'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR07'

        endif

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR06'

        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                          &
                                    rugosity_begin, rugosity_end, BlockFound,               &
                                    STAT = STAT_CALL)
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR08'

        if(BlockFound)then

            nullify(Me%Rugosity%Field)
            allocate(Me%Rugosity%Field(ILB:IUB, JLB:JUB))
            Me%Rugosity%Field(:,:) = FillValueReal

            call ConstructFillMatrix  (PropertyID           = Me%Rugosity%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExtWater%WaterPoints2D,        &
                                       Matrix2D             = Me%Rugosity%Field,                &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR09'


            call GetDefaultValue(Me%Rugosity%ID%ObjFillMatrix, Me%Rugosity%Scalar, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR10'

            call KillFillMatrix(Me%Rugosity%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR11'
           
            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR13'

        else if (.not. BlockFound .and. .not. Me%Manning) then

            call GetData(Me%ChezyCoef,                                                  &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword        = 'CHEZY_COEFFICIENT',                          &
                         SearchType     = FromFile,                                     &
                         ClientModule   = 'ModuleInterfaceSedimentWater',               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &         
                stop 'ConstructShearStatistics - ModuleInterfaceSedimentWater - ERR200'
            if (iflag == 1) then

                Me%Chezy = ON

            else
                write(*,*)'Must define rugosity in ModuleInterfaceSedimentWater' 
                stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR14'
            endif

        end if

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceSedimentWater - ERR12'


        if (BlockFound .and. Me%Manning) then

            write(*,*)'Both absolute rugosity and Manning coefficient were defined' 
            write(*,*)'by default the hydrodynamic module use the Manning coefficient' 
            write(*,*)'ConstructRugosity - ModuleInterfaceSedimentWater - WARN01'

        endif

        if(Me%WaveShear_Stress%Yes)then

            call Read_Property_2D(Me%WaveShear_Stress%Rugosity, FromBlock, &
                                  waverugosity_begin, waverugosity_end)

        endif

    end subroutine ConstructRugosity


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                                 :: ClientNumber
        integer                                 :: STAT_CALL
        logical                                 :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_Property), pointer              :: NewProperty

        !------------------------------------------------------------------------

        !Initialize the properties number   
        Me%PropertiesNumber = 0


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = prop_block_begin,     &
                                        block_end       = prop_block_end,       &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Property 
                    Call Construct_Property(NewProperty, ClientNumber)

                    ! Add new Property to the WaterProperties List 
                    Call Add_Property(NewProperty)
                else cd2
                    call Block_Unlock(Me%ObjEnterData , ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_PropertyList - ModuleInterfaceSedimentWater - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_PropertyList - ModuleInterfaceSedimentWater - ERR02'
            else cd1
                stop 'Construct_PropertyList - ModuleInterfaceSedimentWater - ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------


    subroutine Construct_Property(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty
        integer                         :: ClientNumber

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_Property - ModuleInterfaceSedimentWater - ERR01' 

        nullify(NewProperty%Mass_Available       )
        nullify(NewProperty%FluxToWater          )
        nullify(NewProperty%FluxToSediment       )
        nullify(NewProperty%ErosionCoefficient   )
        nullify(NewProperty%Prev,NewProperty%Next)

        !Construct property ID
        call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        !Construct property evolution parameters
        call Construct_PropertyEvolution(NewProperty, ClientNumber)

        if (NewProperty%Evolution%Variable) Me%OutPut%WriteFinalFile = .true.

        !Construct property values
        call Construct_PropertyValues   (NewProperty)


        !Defines the property output
        call Construct_PropertyOutPut   (NewProperty)


    end subroutine Construct_Property
        
    !--------------------------------------------------------------------------
    
    subroutine Construct_PropertyValues(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),        pointer    :: NewProperty

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: iflag, i, j
        integer                             :: ILB, IUB, JLB, JUB
        integer                             :: WILB, WIUB, WJLB, WJUB
        
        !----------------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        WILB = Me%WorkSize2D%ILB
        WIUB = Me%WorkSize2D%IUB
        WJLB = Me%WorkSize2D%JLB
        WJUB = Me%WorkSize2D%JUB

        !By default there's always mass limitation for every property
        call GetData(NewProperty%Mass_Limitation,                               & 
                     Me%ObjEnterData, iflag,                                    &
                     keyword      ='MASS_LIMITATION',                           &
                     SearchType   = FromBlock,                                  &
                     ClientModule = 'ModuleInterfaceSedimentWater',             &
                     Default      = .true.,                                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR01'
        
        !By default the minimum mass available is 1e-6 kg/m^2 
        if (NewProperty%Mass_Limitation) then 
            
            call GetData(NewProperty%Mass_Min,                                  &
                         Me%ObjEnterData, iflag,                                &
                         keyword      ='MASS_MIN',                              &
                         SearchType   = FromBlock,                              &
                         ClientModule = 'ModuleInterfaceSedimentWater',         &
                         Default      = 1e-6,                                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR02'

            call GetData(NewProperty%Old,                                       &
                         Me%ObjEnterData, iflag,                                &
                         keyword      ='OLD',                                   &
                         SearchType   = FromBlock,                              &
                         ClientModule = 'ModuleInterfaceSedimentWater',         &
                         Default      = .false.,                                &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR03'


            allocate(NewProperty%Mass_Available(ILB:IUB, JLB:JUB))
            NewProperty%Mass_Available  = null_real
            allocate(NewProperty%MassInKg      (ILB:IUB, JLB:JUB))
            NewProperty%MassInKg        = null_real
            allocate(NewProperty%MassInKgFromWater(ILB:IUB, JLB:JUB))
            allocate(NewProperty%Mass_FromWater(ILB:IUB, JLB:JUB))

            if(NewProperty%Old)then

                call Read_Old_Properties_2D(NewProperty%Mass_Available, NewProperty%ID%Name)
            
            else

                call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                           EnterDataID          = Me%ObjEnterData,                  &
                                           TimeID               = Me%ObjTime,                       &
                                           HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                           ExtractType          = FromBlock,                        &
                                           PointsToFill2D       = Me%ExtWater%WaterPoints2D,        &
                                           Matrix2D             = NewProperty%Mass_Available,       &
                                           TypeZUV              = TypeZ_,                           &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) &
                    stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR04'

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) &
                    stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR05'

                !Test to verify mass matrix consistence
                do j = WJLB, WJUB
                do i = WILB, WIUB

                    if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                        if(NewProperty%Mass_Available(i,j) < NewProperty%Mass_Min)then

                            NewProperty%Mass_Available(i,j) = NewProperty%Mass_Min

                        end if

                    end if

                enddo
                enddo
        
            end if

        end if

        if((.not. NewProperty%Particulate) .or. NewProperty%Evolution%WaterFluxes)then

            allocate(NewProperty%WaterConcentration(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR07'
            NewProperty%WaterConcentration(:,:) = null_real

        end if

       
       if(NewProperty%Evolution%BenthicOnly)then

            allocate(NewProperty%WaterConcentration(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR07.1'
            NewProperty%WaterConcentration(:,:) = 0.

        end if
       
       
        if(NewProperty%Evolution%WaterFluxes .or. NewProperty%Evolution%SedimentWaterFluxes)then
            
            allocate(NewProperty%FluxToWater(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR06'
            NewProperty%FluxToWater(:,:) = 0.


            if(NewProperty%Evolution%Erosion)then

                allocate(NewProperty%ErosionCoefficient(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR08'
                NewProperty%ErosionCoefficient(:,:) = null_real

                allocate(NewProperty%ErosionFlux(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR09'
                NewProperty%ErosionFlux(:,:) = null_real

            end if


            if(NewProperty%Evolution%Deposition)then

                allocate(NewProperty%DepositionFlux(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR10'
                NewProperty%DepositionFlux(:,:) = null_real

            end if

        end if

        if(NewProperty%Evolution%SedimentFluxes .or. NewProperty%Evolution%SedimentWaterFluxes)then

            allocate(NewProperty%FluxToSediment(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR11'
            NewProperty%FluxToSediment(:,:) = 0.

            allocate(NewProperty%SedimentConcentration(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_PropertyValues - ModuleInterfaceSedimentWater - ERR12'
            NewProperty%FluxToSediment(:,:) = null_real

            
        end if

    end subroutine Construct_PropertyValues

    
    !--------------------------------------------------------------------------


    subroutine Construct_PropertyEvolution(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty
        integer                                 :: ClientNumber

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        real                                    :: ModelDT

        !Local-----------------------------------------------------------------
        integer                                 :: iflag
        real                                    :: ErrorAux, auxFactor, DTaux
        logical                                 :: VariableDT

        !----------------------------------------------------------------------

        !Checks if the user wants this property to be particulate.
        !This property will be used to define particulated properties
        call GetData(NewProperty%Particulate,                                            &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'PARTICULATE',                                       &
                     Default      = .false.,                                             &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)&
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR00' 

        if (NewProperty%Particulate)then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// ' is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR10'
            end if
        endif

        call GetData(NewProperty%Evolution%Benthos,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BENTHOS',                                          &
                     Default      = OFF,                                                &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR20'

        if(NewProperty%evolution%Benthos) NewProperty%Evolution%Variable = .true.
       
       
        call GetData(NewProperty%Evolution%BenthicEcology,                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BENTHICECOLOGY',                                   &
                     Default      = OFF,                                                &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR25'
        
         if(NewProperty%evolution%BenthicEcology) NewProperty%Evolution%Variable = .true.
        
        
        call GetData(NewProperty%Evolution%BenthicOnly,                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BENTHICONLY',                                   &
                     Default      = OFF,                                                &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR25'
            
            
            
            
       

        !CEQUALW2 as a sink and source      
        call GetData(NewProperty%Evolution%CEQUALW2,                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'CEQUALW2',                                         &
                     Default      = OFF,                                                &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR30'

        if(NewProperty%Evolution%CEQUALW2) NewProperty%Evolution%Variable = .true.
    
        ! This property is a DETRITUS component
        call GetData(NewProperty%Evolution%Detritus,                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DETRITUS',                                         &
                     Default      = OFF,                                                &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR40'

        if(NewProperty%Evolution%Detritus) NewProperty%Evolution%Variable = .true.


        !Fluxes between water and sediment column
        call GetData(NewProperty%Evolution%SedimentWaterFluxes,                          &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SEDIMENT_WATER_FLUXES',                             &
                     Default      = .false.,                                             &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR60'

        if (NewProperty%Evolution%SedimentWaterFluxes)then

            NewProperty%Evolution%Variable = .true.

            if(.not. NewProperty%Particulate)then

                call Read_Property_2D(NewProperty%MolecularDifCoef, FromBlockInBlock,          &
                                      diff_coef_begin, diff_coef_end, ClientNumber = ClientNumber)

            end if

        end if


        !Fluxes between water and interface sediment-water
        call GetData(NewProperty%Evolution%WaterFluxes,                                  &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'WATER_FLUXES',                                      &
                     Default      = .false.,                                             &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR80'

        if (NewProperty%Evolution%WaterFluxes)then
            NewProperty%Evolution%Variable = .true.
        end if


        !Compute erosion fluxes  no - 0;  yes - 1      
        call GetData(NewProperty%Evolution%Erosion,                                      &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'EROSION',                                           &
                     Default      = .false.,                                             &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR90'

        if (NewProperty%Evolution%Erosion)then
            if(.not. NewProperty%Particulate)then
                write(*,*)
                write(*,*)'Cannot specify EROSION for a dissolved property.'
                write(*,*)'Property: '//trim(NewProperty%ID%Name)
                write(*,*)
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR100'
            end if
            if(.not. NewProperty%Evolution%WaterFluxes)then
                write(*,*)
                write(*,*)'Must specify WATER_FLUXES to property '//trim(NewProperty%ID%Name)
                write(*,*)'in order to compute its EROSION'
                write(*,*)
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR110'
            end if
            
            NewProperty%Evolution%Variable = .true.
                        
        end if


        !Compute deposition fluxes  no - 0;  yes - 1      
        call GetData(NewProperty%Evolution%Deposition,                                   &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'DEPOSITION',                                        &
                     Default      = .false.,                                             &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR120'

        if (NewProperty%Evolution%Deposition)then
            if(.not. NewProperty%Particulate)then
                write(*,*)'Cannot specify DEPOSITION for a dissolved property.'
                write(*,*)'Property: '//trim(NewProperty%ID%Name)
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR130'
            end if
            if(.not. NewProperty%Evolution%WaterFluxes)then
                write(*,*)
                write(*,*)'Must specify WATER_FLUXES to property '//trim(NewProperty%ID%Name)
                write(*,*)'in order to compute its DEPOSITION'
                write(*,*)
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR140'
            end if
            
            NewProperty%Evolution%Variable = .true.
                        
        end if


        ! This property has fluxes between sediment and 
        ! sediment-water interface? no - 0;  yes - 1
        call GetData(NewProperty%Evolution%SedimentFluxes,                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SEDIMENT_FLUXES',                                  &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR150'

        if (NewProperty%Evolution%SedimentFluxes)NewProperty%Evolution%Variable = .true.
            

        !Time Step if the property field is variable in time
        if ((NewProperty%Evolution%Variable         ).or. &
            (NewProperty%ID%IDNumber.eq.Temperature_).or. &
            (NewProperty%ID%IDNumber.eq.Salinity_   )) then

            call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR160'

            call GetVariableDT (Me%ObjTime, VariableDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR170'
    
            if (VariableDT) then
                
                NewProperty%Evolution%DTInterval = ModelDT

            else

                call GetData(NewProperty%Evolution%DTInterval,                          &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType     = FromBlock,                                &
                             keyword        = 'DT_INTERVAL',                            &
                             Default        = ModelDT,                                  &
                             ClientModule   = 'ModuleInterfaceSedimentWater',           &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR180'


                if (NewProperty%evolution%DTInterval < (ModelDT)) then
                    write(*,*) 
                    write(*,*) ' Time step error.'
                    stop       'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR190'

                elseif (NewProperty%evolution%DTInterval > (ModelDT)) then

                    !Property DT  must be a multiple of the ModelDT
                    auxFactor = NewProperty%evolution%DTInterval  / ModelDT

                    Erroraux = auxFactor - int(auxFactor)
                    
                    if (Erroraux /= 0) then
                        write(*,*) 
                        write(*,*) ' Time step error.'
                        stop       'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR200'
                    endif

                        ! Run period in seconds
                        DTaux = Me%EndTime - Me%ExternalVar%Now

                        !The run period must be a multiple of the Property DT
                        auxFactor = DTaux / NewProperty%evolution%DTInterval

                        ErrorAux = auxFactor - int(auxFactor)
                        
                        if (ErrorAux /= 0) then
                            write(*,*) 
                            write(*,*) ' Time step error.'
                            stop       'Construct_PropertyEvolution - ModuleInterfaceSedimentWater - ERR210'
                        endif

                    end if

                endif

            NewProperty%Evolution%NextCompute = Me%ExternalVar%Now + NewProperty%Evolution%DTInterval

        else 

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution
    
    
    !-------------------------------------------------------------------------
    
    
    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: iflag

        !Begin----------------------------------------------------------------

        !Checks out if the user pretends to write a HDF format file for this property
        call GetData(NewProperty%OutputHDF,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'OUTPUT_HDF',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceSedimentWater',                   &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceSedimentWater - ERR01'
           
        ! Checks out if the user pretends to write a time serie for this property
        call GetData(NewProperty%TimeSerie,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceSedimentWater',                   &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceSedimentWater - ERR02'

        ! Checks out if the user pretends to write a time serie inside each box for this property
        call GetData(NewProperty%BoxTimeSerie,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOX_TIME_SERIE',                                 &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceSedimentWater',                   &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceSedimentWater - ERR03'

    end subroutine Construct_PropertyOutPut
    
    !--------------------------------------------------------------------------
    
    subroutine Add_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer     :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber  = 1
            Me%FirstProperty     => NewProperty
            Me%LastProperty      => NewProperty
        else
            NewProperty%Prev     => Me%LastProperty
            Me%LastProperty%Next => NewProperty
            Me%LastProperty      => NewProperty
            Me%PropertiesNumber  = Me%PropertiesNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Property 

    !--------------------------------------------------------------------------

    
    subroutine Add_BenthicRate(NewBenthicRate)

        !Arguments-------------------------------------------------------------
        type(T_BenthicRate), pointer       :: NewBenthicRate

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstBenthicRate)) then
            Me%BenthicRatesNumber      = 1
            Me%FirstBenthicRate        => NewBenthicRate
            Me%LastBenthicRate         => NewBenthicRate
        else
            NewBenthicRate%Prev        => Me%LastBenthicRate
            Me%LastBenthicRate%Next    => NewBenthicRate
            Me%LastBenthicRate         => NewBenthicRate
            Me%BenthicRatesNumber      = Me%BenthicRatesNumber + 1
        end if 

    end subroutine Add_BenthicRate 

    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE
        
        !----------------------------------------------------------------------

        WorkILB = Me%WaterWorkSize3D%ILB 
        WorkIUB = Me%WaterWorkSize3D%IUB 
        WorkJLB = Me%WaterWorkSize3D%JLB 
        WorkJUB = Me%WaterWorkSize3D%JUB
        WorkKLB = Me%WaterWorkSize3D%KLB
        WorkKUB = Me%WaterWorkSize3D%KUB

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, trim(Me%Files%Results)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR02'
        
        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR03'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",               &
                              Array2D = Me%ExtWater%Bathymetry,                     &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR04'

        
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                          &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR05'

        call HDF5WriteData  (Me%ObjHDF5, "/Grid", "WaterPoints",                    &
                             "-", Array3D = Me%ExtWater%WaterPoints3D,              &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR06'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceSedimentWater - ERR07'
       
        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
        
        nullify(Me%OutPut%OutTime)

        call GetOutPutTime(Me%ObjEnterData,                              &
                           CurrentTime = Me%ExternalVar%Now,             &
                           EndTime     = Me%EndTime,                     &
                           keyword     = 'OUTPUT_TIME',                  &
                           SearchType  = FromFile,                       &
                           OutPutsTime = Me%OutPut%OutTime,              &
                           OutPutsOn   = Me%OutPut%Yes,                  &
                           STAT        = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                       &
            stop 'ConstructGlobalOutput - ModuleInterfaceSedimentWater - ERR01' 

        if (Me%OutPut%Yes) then

            Me%OutPut%NextOutPut = 1

            call Open_HDF5_OutPut_File

        else
            write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
            write(*,*)'one property has HDF format outputs.'
            stop 'ConstructGlobalOutput - ModuleInterfaceSedimentWater - ERR02'
        endif

        call GetOutPutTime(Me%ObjEnterData,                                         &
                           CurrentTime = Me%ExternalVar%Now,                        &
                           EndTime     = Me%EndTime,                                &
                           keyword     = 'RESTART_FILE_OUTPUT_TIME',                &
                           SearchType  = FromFile,                                  &
                           OutPutsTime = Me%OutPut%RestartOutTime,                  &
                           OutPutsOn   = Me%OutPut%WriteRestartFile,                &
                           STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                       &
            stop 'ConstructGlobalOutput - ModuleInterfaceSedimentWater - ERR03'

        if(Me%OutPut%WriteRestartFile)then

            Me%OutPut%NextRestartOutput = 1

        end if 

        !<BeginKeyword>
            !Keyword          : RESTART_FILE_OVERWRITE
            !<BeginDescription>       
               ! 
               ! This option checks wether the restart file is to be overwritten or not
               ! 
            !<EndDescription>
            !Type             : logical 
            !Default          : .true.
            !Multiple Options : Do not have
            !Search Type      : FromFile
        !<EndKeyword>
        
        call GetData(Me%OutPut%RestartOverwrite,                                        &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESTART_FILE_OVERWRITE',                           &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleWaterProperties',                            &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalOutput - ModuleInterfaceSedimentWater - ERR04'



    end subroutine ConstructGlobalOutput

    !------------------------------------------------------------------------

    subroutine Construct_BenthicRateList

        !External----------------------------------------------------------------
        integer                                :: ClientNumber
        integer                                :: STAT_CALL
        logical                                :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_BenthicRate),    pointer      :: NewBenthicRate

        !------------------------------------------------------------------------

 
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
                                        ClientNumber    = ClientNumber,     &
                                        block_begin     = rate_block_begin, &
                                        block_end       = rate_block_end,   &
                                        BlockFound      = BlockFound,       &
                                        STAT            = STAT_CALL)
            if(STAT_CALL .EQ. SUCCESS_)then    
                if (BlockFound) then                                                  
                    
                    !Construct a New Benthic Rate
                    Call Construct_BenthicRate  (NewBenthicRate)

                    !Add new Rate to the Benthic Rates List 
                    Call Add_BenthicRate    (NewBenthicRate)

                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Construct_BenthicRateList - ModuleInterfaceSedimentWater - ERR01'

                    exit do1    !No more blocks
                end if


            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_BenthicRateList - ModuleInterfaceSedimentWater - ERR02'
            else
                stop       'Construct_BenthicRateList - ModuleInterfaceSedimentWater - ERR03'
            end if
        end do do1
         
    end subroutine Construct_BenthicRateList

    !--------------------------------------------------------------------------
    
    subroutine Construct_BenthicRate(NewBenthicRate)

        !Arguments-------------------------------------------------------------
        type(T_BenthicRate), pointer        :: NewBenthicRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewBenthicRate, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                          &
            stop 'Construct_BenthicRate - ModuleInterfaceSedimentWater - ERR01' 

        nullify(NewBenthicRate%Field, NewBenthicRate%Prev, NewBenthicRate%Next)

        call Construct_BenthicRateID       (NewBenthicRate)

        call Construct_BenthicRateValues   (NewBenthicRate)


    end subroutine Construct_BenthicRate

    !--------------------------------------------------------------------------
    
    subroutine Construct_BenthicRateID(NewBenthicRate)

        !Arguments-------------------------------------------------------------
        type(T_BenthicRate), pointer       :: NewBenthicRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: iflag, PropNumber
        logical                             :: CheckName
        type (T_Property), pointer          :: PropertyX
      
        !----------------------------------------------------------------------

        !First Property defined in a rate relation
        call GetData(NewBenthicRate%FirstProp%name,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'FIRSTPROP',                                        &
                     ClientModule = 'ModuleInterfaceSedimentWater',                     &
                     SearchType   = FromBlock,                                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR01' 
        if (iflag==0)                                                                   &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR02' 

        !Check if the property name is valid
        CheckName = CheckPropertyName(NewBenthicRate%FirstProp%Name, Number = PropNumber)
        if (CheckName) then
            NewBenthicRate%FirstProp%IDnumber = PropNumber
        else
            write(*,*)
            write(*,*) 'The first property name is not recognised by the model.'
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR03' 
        end if 

        call Search_Property(PropertyX, PropertyXID = PropNumber, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                     &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR04' 
        
        !second Property defined in a rate relation
        call GetData(NewBenthicRate%SecondProp%Name,                                     &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'SECONDPROP',                                        &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR05' 
        if (iflag==0)                                                                    &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR06' 
      
        ! Check if the property name is valid OR not
        CheckName = CheckPropertyName(NewBenthicRate%SecondProp%name, Number = PropNumber)
        if (CheckName) then
            NewBenthicRate%SecondProp%IDnumber = PropNumber
        else
            write(*,*)
            write(*,*) 'The Second property name is not recognised by the model.'
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR07' 
        end if
        
        call Search_Property(PropertyX, PropertyXID = PropNumber, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR08' 
  
  
       !Rate description ex: zooplankton grazing over phytoplankton
        call GetData(NewBenthicRate%ID%Description,                                      &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'DESCRIPTION',                                       &
                     Default      = 'No description was given.',                         &
                     ClientModule = 'ModuleInterfaceSedimentWater',                      &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR09' 
         
        !Rate name (this is the name of the output boxes file)
        call GetData(NewBenthicRate%ID%Name,                                               &
                     Me%ObjEnterData, iflag,                                               &
                     keyword      = 'NAME',                                                &
                     ClientModule = 'ModuleInterfaceSedimentWater',                        &
                     SearchType   = FromBlock,                                             &
                     Default      = 'No name was given to sediment rate.',                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                         &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR10' 
        if (iflag==0)                                                                      &
            stop 'Construct_BenthicRateID - ModuleInterfaceSedimentWater - ERR11' 

    end subroutine Construct_BenthicRateID
    
    
    !--------------------------------------------------------------------------

    
    subroutine Construct_BenthicRateValues(NewBenthicRate)

        !Arguments-------------------------------------------------------------
        type(T_BenthicRate), pointer        :: NewBenthicRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB
          
        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB


        allocate(NewBenthicRate%Field(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_BenthicRateValues - ModuleInterfaceSedimentWater - ERR01' 
        NewBenthicRate%Field(:,:) = FillValueReal

    end subroutine Construct_BenthicRateValues
    
    
    !----------------------------------------------------------------------

    
    subroutine ConstructLog


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty


#ifndef _OUTPUT_OFF_
        write(*, *)"---------------- INTERFACE SEDIMENT-WATER -----------------"
        write(*, *)
        write(*, *)"Num of Properties : ", Me%PropertiesNumber
        write(*, *)

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))

            write(*, *)"Property          : ", trim(CurrentProperty%ID%Name)
            write(*, *)"---Benthos        : ", CurrentProperty%Evolution%Benthos
            write(*, *)"---BenthicEcology : ", CurrentProperty%Evolution%BenthicEcology
            write(*, *)"---BenthicCEQUALW2: ", CurrentProperty%Evolution%CEQUALW2
            write(*, *)"---Detritus       : ", CurrentProperty%Evolution%Detritus
            write(*, *)"---WaterFluxes    : ", CurrentProperty%Evolution%WaterFluxes
            write(*, *)"---SedimentFluxes : ", CurrentProperty%Evolution%SedimentFluxes
            write(*, *)"---Sed-Water Flux : ", CurrentProperty%Evolution%SedimentWaterFluxes
            write(*, *)"---Erosion        : ", CurrentProperty%Evolution%Erosion
            write(*, *)"---Deposition     : ", CurrentProperty%Evolution%Deposition
            write(*, *)

            CurrentProperty=>CurrentProperty%Next
        enddo
#endif

    end subroutine ConstructLog
    
    
    !--------------------------------------------------------------------------

    
    subroutine Construct_Sub_Modules

        !Local-----------------------------------------------------------------
        type (T_Property),           pointer                 :: PropertyX            
        integer                                              :: STAT_CALL

        !----------------------------------------------------------------------

        Me%Coupled%WaterFluxes%NumberOfProperties          = 0
        Me%Coupled%SedimentWaterFluxes%NumberOfProperties  = 0
        Me%Coupled%SedimentFluxes%NumberOfProperties       = 0
        Me%Coupled%Benthos%NumberOfProperties              = 0
        Me%Coupled%CEQUALW2%NumberOfProperties             = 0
        Me%Coupled%Detritus%NumberOfProperties             = 0
        Me%Coupled%TimeSerie%NumberOfProperties            = 0
        Me%Coupled%BoxTimeSerie%NumberOfProperties         = 0
        Me%Coupled%OutputHDF%NumberOfProperties            = 0
        Me%Coupled%BenthicEcology%NumberOfProperties       = 0

        PropertyX => Me%FirstProperty

do1 :   do while (associated(PropertyX))

            if (PropertyX%Evolution%WaterFluxes) then
                Me%Coupled%WaterFluxes%NumberOfProperties           = &
                Me%Coupled%WaterFluxes%NumberOfProperties           + 1
                Me%Coupled%WaterFluxes%Yes                          = ON
            endif


            if (PropertyX%Evolution%SedimentWaterFluxes) then
                Me%Coupled%SedimentWaterFluxes%NumberOfProperties   = &
                Me%Coupled%SedimentWaterFluxes%NumberOfProperties   + 1
                Me%Coupled%SedimentWaterFluxes%Yes                  = ON

            endif

            if (PropertyX%Evolution%SedimentFluxes) then
                Me%Coupled%SedimentFluxes%NumberOfProperties        = &
                Me%Coupled%SedimentFluxes%NumberOfProperties        + 1
                Me%Coupled%SedimentFluxes%Yes                       = ON
            endif

            if (PropertyX%Evolution%Benthos) then
                Me%Coupled%Benthos%NumberOfProperties               = &
                Me%Coupled%Benthos%NumberOfProperties               + 1
                Me%Coupled%Benthos%Yes                              = ON
            endif
            
             if (PropertyX%Evolution%BenthicEcology) then
                Me%Coupled%BenthicEcology%NumberOfProperties        = &
                Me%Coupled%BenthicEcology%NumberOfProperties        + 1
                Me%Coupled%BenthicEcology%Yes                       = ON
            endif

            if (PropertyX%Evolution%CEQUALW2) then
                Me%Coupled%CEQUALW2%NumberOfProperties              = &
                Me%Coupled%CEQUALW2%NumberOfProperties              + 1
                Me%Coupled%CEQUALW2%Yes                             = ON
            endif

            if (PropertyX%Evolution%Detritus) then
                Me%Coupled%Detritus%NumberOfProperties              = &
                Me%Coupled%Detritus%NumberOfProperties              + 1
                Me%Coupled%Detritus%Yes                             = ON
            endif

            if (PropertyX%Evolution%Erosion) then
                Me%Coupled%Erosion%NumberOfProperties               = &
                Me%Coupled%Erosion%NumberOfProperties               + 1
                Me%Coupled%Erosion%Yes                              = ON
            endif

            if (PropertyX%Evolution%Deposition) then
                Me%Coupled%Deposition%NumberOfProperties            = &
                Me%Coupled%Deposition%NumberOfProperties            + 1
                Me%Coupled%Deposition%Yes                           = ON
            endif


            if (PropertyX%TimeSerie) then
                Me%Coupled%TimeSerie%NumberOfProperties             = &
                Me%Coupled%TimeSerie%NumberOfProperties             + 1
                Me%Coupled%TimeSerie%Yes                            = ON
            endif

            if (PropertyX%BoxTimeSerie) then
                Me%Coupled%BoxTimeSerie%NumberOfProperties          = &
                Me%Coupled%BoxTimeSerie%NumberOfProperties          + 1
                Me%Coupled%BoxTimeSerie%Yes                         = ON
            endif
            
            if (PropertyX%OutputHDF) then
                Me%Coupled%OutputHDF%NumberOfProperties             = &
                Me%Coupled%OutputHDF%NumberOfProperties             + 1
                Me%Coupled%OutputHDF%Yes                            = ON
            endif

            PropertyX=>PropertyX%Next

        end do do1

        if(Me%Coupled%Benthos%Yes .and. Me%Coupled%CEQUALW2%Yes)then
            
            write(*,*)'Benthos and CEQUALW2 models cannot be simulated at the same time.'
            stop 'Construct_Sub_Modules - ModuleInterfaceSedimentWater - ERR01'

        end if

        if(Me%Coupled%Benthos%Yes .and. Me%Coupled%BenthicEcology%Yes)then
            
            write(*,*)'Benthos and BenthicEcology models cannot be simulated at the same time.'
            stop 'Construct_Sub_Modules - ModuleInterfaceSedimentWater - ERR01.1'

        end if
        
        if(Me%Coupled%BenthicEcology%Yes .and. Me%Coupled%CEQUALW2%Yes)then
            
            write(*,*)'BenthicEcology and CEQUALW2 models cannot be simulated at the same time.'
            stop 'Construct_Sub_Modules - ModuleInterfaceSedimentWater - ERR01.2'

        end if


        if(Me%Coupled%Detritus%Yes .and. .not. Me%Coupled%CEQUALW2%Yes)then
            
            write(*,*)'Detritus and CEQUALW2 must be simulated at the same time.'
            stop 'Construct_Sub_Modules - ModuleInterfaceSedimentWater - ERR02'

        end if


        if(Me%Coupled%WaterFluxes%Yes) then

            call CheckOptionsWaterFluxes

        end if

#ifndef _SEDIMENT_
        if(Me%Coupled%SedimentFluxes%Yes) then

            call CheckOptionsSedimentFluxes

        end if

        if(Me%Coupled%SedimentWaterFluxes%Yes) then

            call CheckOptionsSedimentWaterFluxes

        end if
#endif
 
        if(Me%Coupled%Benthos%Yes)then

            call CoupleBenthos

        end if
        
       if(Me%Coupled%BenthicEcology%Yes)then

            call CoupleBenthicEcology

        end if

        if(Me%Coupled%CEQUALW2%Yes)then

            call CoupleCEQUALW2

            call ReadSOD
                
            if (Me%UseSOD)then
                call SetSOD(Me%SOD%Field, Me%ExtWater%OpenPoints2D, Me%ExtWater%WaterPoints2D)
            end if            

        end if

   
        if(Me%Coupled%TimeSerie%Yes)then

            call Construct_Time_Serie
             
        end if
        
        if(Me%Coupled%BoxTimeSerie%Yes)then

            call StartOutputBoxFluxes
               
        end if

        if(Me%Coupled%OutputHDF%Yes)then

            call ConstructGlobalOutput

        end if

        call GetWaterPropertiesBottomOptions(WaterPropertiesID = Me%ObjWaterProperties, &
                                             MacroAlgae        = Me%MacroAlgae,         &
                                             STAT              = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Sub_Modules - ModuleInterfaceSedimentWater - ERR03' 


    end subroutine Construct_Sub_Modules


    !--------------------------------------------------------------------------


    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer                                             :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber  
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: nProperties
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !----------------------------------------------------------------------

        !First checks out how many properties will have time series
        PropertyX   => Me%FirstProperty
        nProperties =  1
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) nProperties = nProperties + 1
            PropertyX=>PropertyX%Next
        enddo


        !Allocates PropertyList
        allocate(PropertyList(nProperties), STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR10'

        PropertyList(1) = 'ShearStress'

        !Fills up PropertyList
        PropertyX   => Me%FirstProperty
        nProperties =  1
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                nProperties = nProperties + 1
                PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name))
            endif
            PropertyX=>PropertyX%Next
        enddo

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleWaterProperties',                            &
                     Default      = Me%Files%InputData,                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR20' 


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                    &
                            trim(TimeSerieLocationFile),                    &
                            PropertyList, "srb",                            &
                            WaterPoints3D = Me%ExtWater%WaterPoints3D,      &
                            ModelName     = Me%ModelName,                   & 
                            STAT          = STAT_CALL)
        if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR30'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR40'

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR50'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR60'

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR70'

                    if (IgnoreOK) then
                        cycle
                    else
                        stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR80'
                    endif

                endif


                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceSedimentWater - ERR90'
            endif

        enddo

        
    end subroutine Construct_Time_Serie

    !----------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        integer                                             :: ILB, IUB, JLB, JUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        type(T_Property),                           pointer :: PropertyX
        type(T_BenthicRate    ),                    pointer :: BenthicRateX
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList
        integer                                             :: nScalars, n

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        ! This keyword have two functions if exist fluxes between boxes are compute 
        ! and the value read is the name file where the boxes are defined
        call GetData(Me%Files%BoxesFile,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule = 'ModuleInterfaceSedimentWater',                 &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR01'
        if (iflag .EQ. 0)                                                           &
            stop 'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR02'    
        
        inquire(File = Me%Files%BoxesFile, Exist = exist)
        if (exist) then
            inquire(File = Me%Files%BoxesFile, Opened  = Opened)
            if (opened) then
                write(*,*    ) 
                write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Files%BoxesFile))
                write(*,*    ) 'Already opened.'
                stop           'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR03'    
            end if
        else
            write(*,*) 
            write(*,*)     'Could not find the boxes file.'
            write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
            stop           'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR04'    
        end if

        nScalars = Me%Coupled%BoxTimeSerie%NumberOfProperties + Me%BenthicRatesNumber

        allocate(ScalarOutputList(nScalars), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR06'

        PropertyX  => Me%FirstProperty
        n = 0
        do while (associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then
                n = n + 1
                ScalarOutputList(n) = "Bottom "//trim(PropertyX%ID%name)
            end if 

            PropertyX=>PropertyX%Next
        end do

        BenthicRateX => Me%FirstBenthicRate
        do while(associated(BenthicRateX))
            n = n + 1
            ScalarOutputList(n) = "Bottom "//trim(BenthicRateX%ID%name)
            BenthicRateX => BenthicRateX%Next
        end do

        call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &
                         TimeID             = Me%ObjTime,                   &
                         HorizontalGridID   = Me%ObjHorizontalGrid,         &
                         BoxesFilePath      = Me%Files%BoxesFile,           &
                         ScalarOutputList   = ScalarOutputList,             &
                         WaterPoints2D      = Me%ExtWater%WaterPoints2D,    &
                         STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR07'

        deallocate(ScalarOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR07'

        allocate(Me%Scalar2D(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceSedimentWater - ERR12'
        Me%Scalar2D(:,:) = 0.

    end subroutine StartOutputBoxFluxes


    !--------------------------------------------------------------------------

    subroutine CoupleBenthos

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX, Temperature
        integer, pointer, dimension(:)                      :: BenthosPropertyList
        integer                                             :: STAT_CALL
        real                                                :: BenthosDT
        integer                                             :: Index = 0

        !----------------------------------------------------------------------

        Index = 0

        nullify (BenthosPropertyList)
        allocate(BenthosPropertyList(1:Me%Coupled%Benthos%NumberOfProperties))


        call Search_Property(Temperature, PropertyXID = Temperature_, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Please define property temperature in the'
            write(*,*)'InterfaceSedimentWater data file.'
            stop 'CoupleBenthos - ModuleInterfaceSedimentWater - ERR00' 
        end if


        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%Benthos)then
                Index = Index + 1
                BenthosPropertyList(Index)  = PropertyX%ID%IDNumber
            end if

            PropertyX => PropertyX%Next

        enddo

        nullify(PropertyX)

        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = BenthosModel,                  &
                                DT                  = BenthosDT,                     &
                                PropertiesList      = BenthosPropertyList,           &
                                WaterPoints2D       = Me%ExtWater%WaterPoints2D,     &
                                Size2D              = Me%WorkSize2D,                 &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleBenthos - ModuleInterfaceSedimentWater - ERR01'

        Me%Coupled%Benthos%DT_Compute  = BenthosDT 
        Me%Coupled%Benthos%NextCompute = Me%ExternalVar%Now

        deallocate(BenthosPropertyList)
        nullify   (BenthosPropertyList)

    end subroutine CoupleBenthos


  
    !--------------------------------------------------------------------------

    subroutine CoupleBenthicEcology

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX, Temperature
        integer, pointer, dimension(:)                      :: BenthicEcologyPropertyList
        integer                                             :: STAT_CALL
        real                                                :: BenthicEcologyDT
        integer                                             :: Index = 0
        integer                                             :: ILB,IUB,JLB,JUB
        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB


        Index = 0

        nullify (BenthicEcologyPropertyList)
        allocate(BenthicEcologyPropertyList(1:Me%Coupled%BenthicEcology%NumberOfProperties))


        call Search_Property(Temperature, PropertyXID = Temperature_, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) then
            write(*,*)'Please define property temperature in the'
            write(*,*)'InterfaceSedimentWater data file.'
            stop 'CoupleBenthicEcology - ModuleInterfaceSedimentWater - ERR00' 
        end if


        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%BenthicEcology)then
                Index = Index + 1
                BenthicEcologyPropertyList(Index)  = PropertyX%ID%IDNumber
            end if

            PropertyX => PropertyX%Next

        enddo

        nullify(PropertyX)

        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = BenthicEcologyModel,                  &
                                DT                  = BenthicEcologyDT,                     &
                                PropertiesList      = BenthicEcologyPropertyList,           &
                                WaterPoints2D       = Me%ExtWater%WaterPoints2D,     &
                                Size2D              = Me%WorkSize2D,                 &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleBenthicEcology - ModuleInterfaceSedimentWater - ERR01'

        Me%Coupled%BenthicEcology%DT_Compute  = BenthicEcologyDT 
        Me%Coupled%BenthicEcology%NextCompute = Me%ExternalVar%Now

        allocate(Me%ExtWater%WaterVolume(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'CoupleBenthicEcology - ModuleInterfaceSedimentWater - ERR05'
     
       allocate(Me%ExtWater%Sediment(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'CoupleBenthicEcology - ModuleInterfaceSedimentWater - ERR06'
       
      

        deallocate(BenthicEcologyPropertyList)
        nullify   (BenthicEcologyPropertyList)

    end subroutine CoupleBenthicEcology
    !--------------------------------------------------------------------------

    

    subroutine CoupleCEQUALW2

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer, pointer, dimension(:)                      :: CEQUALW2PropertyList
        integer                                             :: STAT_CALL
        real                                                :: CEQUALW2DT
        integer                                             :: Index = 0

        !----------------------------------------------------------------------

        Index = 0

        nullify (CEQUALW2PropertyList)
        allocate(CEQUALW2PropertyList(1:Me%Coupled%CEQUALW2%NumberOfProperties))

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%CEQUALW2)then
                Index = Index + 1
                CEQUALW2PropertyList(Index)  = PropertyX%ID%IDNumber
            end if

            PropertyX => PropertyX%Next

        enddo

        nullify(PropertyX)
        
        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = BenthicCEQUALW2Model,          &
                                DT                  = CEQUALW2DT,                    &
                                PropertiesList      = CEQUALW2PropertyList,          &
                                WaterPoints2D       = Me%ExtWater%WaterPoints2D,     &
                                Size2D              = Me%WorkSize2D,                 &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleCEQUALW2 - ModuleInterfaceSedimentWater - ERR01'

        Me%Coupled%CEQUALW2%DT_Compute  = CEQUALW2DT 
        Me%Coupled%CEQUALW2%NextCompute = Me%ExternalVar%Now

        deallocate(CEQUALW2PropertyList)
        nullify   (CEQUALW2PropertyList)

    end subroutine CoupleCEQUALW2


    !--------------------------------------------------------------------------


    subroutine CheckOptionsWaterFluxes

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                             :: FreeVerticalMovementID
        integer                                             :: ILB, IUB, JLB, JUB
        type(T_Property), pointer                           :: PropertyX
        real                                                :: DTInterval
        logical                                             :: WaterFluxes
        logical                                             :: FreeVerticalMovement

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB


        if(Me%Coupled%Deposition%Yes)then

            call Read_Property_2D (Me%Critical_Shear_Deposition, FromBlock, csd_begin, csd_end)
            
            call GetWaterPropertiesSubModulesID(WaterPropertiesID      = Me%ObjWaterProperties, &
                                                FreeVerticalMovementID = FreeVerticalMovementID,&
                                                STAT                   = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR10'
            
            if(FreeVerticalMovementID /= 0)then
                Me%ObjFreeVerticalMovement = AssociateInstance(mFREEVERTICALMOVEMENT_,FreeVerticalMovementID)
            end if

            allocate(Me%DepositionProbability(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR20'
            Me%DepositionProbability(:,:) = FillValueReal

        end if

        if(Me%Coupled%Erosion%Yes)then

            call Read_Property_2D (Me%Critical_Shear_Erosion, FromBlock, cse_begin, cse_end)

            call Read_Property_2D (Me%ErosionRate, FromBlock, erosion_begin, erosion_end)

        end if

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if(PropertyX%Evolution%WaterFluxes)then

                if(.not. WaterPropertyExists(Me%ObjWaterProperties, PropertyX%ID%IDNumber))then
                    write(*,*)
                    write(*,*)'Property : '//trim(PropertyX%ID%Name)
                    write(*,*)'must be defined as a water property.'
                    stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR50'
                end if

                call GetWaterPropertyOptions(Me%ObjWaterProperties,                 &
                                             PropertyX%ID%IDNumber,                 &
                                             DTInterval,                            &
                                             WaterFluxes,                           &
                                             FreeVerticalMovement,                  &
                                             STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR60'

                if(.not. WaterFluxes)then
                    write(*,*)
                    write(*,*) 'Property '//trim(PropertyX%ID%Name)//' has option WATER_FLUXES'
                    write(*,*) 'activated, but does not have option BOTTOM_FLUXES activated'
                    write(*,*) 'in the Water Properties file. Please review your options.'
                    write(*,*)
                    stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR70'
                end if

                if(DTInterval /= PropertyX%Evolution%DTInterval)then
                    write(*,*)
                    write(*,*) ' Assumed time step for sediment-water interface '
                    write(*,*) ' property '//trim(PropertyX%ID%Name)
                    write(*,*) ' equal to the water property time step.'
                    write(*,*)
                    write(*,*) ' CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - WRN01'
                    
                    PropertyX%Evolution%DTInterval  = DTInterval

                    PropertyX%Evolution%NextCompute = Me%ExternalVar%Now + PropertyX%Evolution%DTInterval

                end if

                if(PropertyX%Evolution%Deposition)then

                    if(.not. FreeVerticalMovement)then
                        write(*,*)
                        write(*,*)'Property : '//trim(PropertyX%ID%Name)//' must be defined'
                        write(*,*)'with option VERTICAL_MOVEMENT in WaterProperties file.'
                        stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR80'
                    end if

                    if(.not. FreeVertPropertyExists(Me%ObjFreeVerticalMovement, PropertyX%ID%IDNumber))then
                        write(*,*)
                        write(*,*)'Property : '//trim(PropertyX%ID%Name)//' must be defined in'
                        write(*,*)'ModuleFreeVerticalMovement.'
                        stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR90'
                    end if

                    if(.not. FreeVertPropertyHasDeposition(Me%ObjFreeVerticalMovement, PropertyX%ID%IDNumber))then
                        write(*,*)
                        write(*,*)'Property : '//trim(PropertyX%ID%Name)//' must be defined in'
                        write(*,*)'ModuleFreeVerticalMovement with DEPOSITION option activated.'
                        stop 'CheckOptionsWaterFluxes - ModuleInterfaceSedimentWater - ERR100'
                    end if

                end if

            endif

            PropertyX=>PropertyX%Next
        enddo

        Me%Coupled%WaterFluxes%NextCompute = Me%ExternalVar%Now
   
    end subroutine CheckOptionsWaterFluxes

    
    !--------------------------------------------------------------------------

#ifndef _SEDIMENT_
    subroutine CheckOptionsSedimentFluxes

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        real                                                :: DTInterval
        logical                                             :: SedimentFluxes
        !----------------------------------------------------------------------
        
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if(PropertyX%Evolution%SedimentFluxes)then

                if(.not. SedimentPropertyExists(Me%ObjSedimentProperties, PropertyX%ID%IDNumber))then
                    write(*,*)
                    write(*,*)'Property : '//trim(PropertyX%ID%Name)
                    write(*,*)'must be defined as a sediment property.'
                    stop 'CheckOptionsSedimentFluxes - ModuleInterfaceSedimentWater - ERR10'
                end if

                call GetSedimentPropertyOptions(Me%ObjSedimentProperties,               &
                                                PropertyX%ID%IDNumber,                  &
                                                DTInterval,                             &
                                                SedimentFluxes,                         &
                                                STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                             &
                    stop 'CheckOptionsSedimentFluxes - ModuleInterfaceSedimentWater - ERR20'

                if(.not. SedimentFluxes)then
                    write(*,*)
                    write(*,*) 'Property '//trim(PropertyX%ID%Name)//' has option SEDIMENT_FLUXES'
                    write(*,*) 'activated, but does not have option SURFACE_FLUXES or            '
                    write(*,*) 'SEDIMENT_WATER_FLUXES activated in the Sediment Properties file.'
                    write(*,*) 'Please review your options.'
                    write(*,*)
                    stop 'CheckOptionsSedimentFluxes - ModuleInterfaceSedimentWater - ERR30'
                end if

                !if(DTInterval /= PropertyX%Evolution%DTInterval)then
                !    write(*,*)
                !    write(*,*) ' Property '//trim(PropertyX%ID%Name)//' must have the same'
                !    write(*,*) ' time step for sediment-water interface and the sediment column.'
                !    write(*,*)
                !    stop       'CheckOptionsSedimentFluxes - ModuleInterfaceSedimentWater - ERR40'
                !end if

                if(PropertyX%Mass_Min > 0.)then
                    write(*,*)
                    write(*,*) ' Property '//trim(PropertyX%ID%Name)//' must have MASS_MIN'
                    write(*,*) ' set to zero, as there are fluxes being computed from the'
                    write(*,*) ' sediment compartment.'
                    write(*,*)
                    stop       'CheckOptionsSedimentFluxes - ModuleInterfaceSedimentWater - ERR50'
                endif

            end if

            PropertyX=>PropertyX%Next

        enddo

        Me%Coupled%SedimentFluxes%NextCompute = Me%ExternalVar%Now

    end subroutine CheckOptionsSedimentFluxes

    !--------------------------------------------------------------------------

    subroutine CheckOptionsSedimentWaterFluxes

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: ILB, IUB, JLB, JUB
        integer                                             :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

                
        PropertyX   => Me%FirstProperty
        do while (associated(PropertyX))

            if(PropertyX%Evolution%SedimentWaterFluxes)then

                if((.not. PropertyX%Evolution%SedimentFluxes) .or. &
                   (.not. PropertyX%Evolution%WaterFluxes)  )then

                    write(*,*)
                    write(*,*)'Property : '//trim(PropertyX%ID%Name)
                    write(*,*)'must have both WATER_FLUXES and SEDIMENT_FLUXES.'
                    stop 'CheckOptionsSedimentWaterFluxes - ModuleInterfaceSedimentWater - ERR10'

                end if

            end if

            PropertyX=>PropertyX%Next
        enddo


        Me%Coupled%SedimentWaterFluxes%NextCompute = Me%ExternalVar%Now

        call GetConsolidationOptions (Me%ObjConsolidation, Me%ExtSed%ComputeConsolidation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'CheckOptionsSedimentWaterFluxes - ModuleInterfaceSedimentWater - ERR30'

        if(Me%ExtSed%ComputeConsolidation .and. .not. Me%Consolidation%Yes)then
            write(*,*)'Please activate option CONSOLIDATION' 
            stop 'CheckOptionsSedimentWaterFluxes - ModuleInterfaceSedimentWater - ERR35'
        endif

        if(.not. Me%ExtSed%ComputeConsolidation .and. Me%Consolidation%Yes)then
            write(*,*)'CONSOLIDATION option is active and ModuleConsolidation is not.' 
            stop 'CheckOptionsSedimentWaterFluxes - ModuleInterfaceSedimentWater - ERR36'
        endif

        allocate(Me%WaterFlux (ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)&
            stop 'CheckOptionsSedimentWaterFluxes - ModuleInterfaceSedimentWater - ERR40'
        Me%WaterFlux(:,:) = 0.

    end subroutine CheckOptionsSedimentWaterFluxes
#endif

    !--------------------------------------------------------------------------


    subroutine Read_Old_Properties_2D(Scalar_2D, PropertyName)

        !Arguments--------------------------------------------------------------
        real, dimension(:,:), pointer               :: Scalar_2D
        character (Len=*), Intent(IN)               :: PropertyName

        !Local-----------------------------------------------------------------
        integer                                     :: IUB, JUB, ILB, JLB 
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
        integer                                     :: ObjHDF5 = 0
        logical                                     :: exist

        !----------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        JLB = Me%WorkSize2D%JLB
        IUB = Me%WorkSize2D%IUB 
        JUB = Me%WorkSize2D%JUB

        ObjHDF5 = 0


        inquire(File = trim(Me%Files%Initial)//"5", Exist = exist)
        
        if(.not. exist)then
            write(*,*) 
            write(*,*)     'Could not find the final InterfaceSedimentWater file.'
            write(*,'(A)') 'BoxFileName = ', trim(Me%Files%Initial)//"5"
            stop           'Read_Old_Properties_2D - ModuleInterfaceSedimentWater - ERR00'    
        end if

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5, trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'read_Old_Properties_2D - ModuleInterfaceSedimentWater - ERR01'
            
        !Reads from HDF file the Property concentration and open boundary values
        call HDF5SetLimits  (ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'read_Old_Properties_2D - ModuleInterfaceSedimentWater - ERR02'
            
        call HDF5ReadData(ObjHDF5, "/Mass/"//PropertyName, &
                          PropertyName, Array2D = Scalar_2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'read_Old_Properties_2D - ModuleInterfaceSedimentWater - ERR03'            
        
        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'read_Old_Properties_2D - ModuleInterfaceSedimentWater - ERR04'            

    end subroutine read_Old_Properties_2D


    !--------------------------------------------------------------------------


    subroutine Read_Property_2D(Property, ExtractType, BeginBlock, EndBlock, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Property_2D)                :: Property
        character(len=*)                    :: BeginBlock, EndBlock
        integer                             :: ExtractType
        integer, optional                   :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: LocalClientNumber

        !Begin-----------------------------------------------------------------

        select case(ExtractType)

            case(FromBlock)
                
                call ExtractBlockFromBuffer (Me%ObjEnterData, LocalClientNumber,    &
                                             BeginBlock, EndBlock,                  &
                                             BlockFound, STAT = STAT_CALL)
            case(FromBlockInBlock)
                
                call ExtractBlockFromBlock (Me%ObjEnterData, ClientNumber,          &
                                            BeginBlock, EndBlock,                   &
                                            BlockFound, STAT = STAT_CALL)
            case default

        end select


        if(present(ClientNumber)) LocalClientNumber = ClientNumber
        
        if (BlockFound) then

            !Allocates Variables
            allocate (Property%Field (Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB))

            Property%Field (:,:) = null_real

            call ConstructFillMatrix  (PropertyID           = Property%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill2D       = Me%ExtWater%WaterPoints2D,        &
                                       Matrix2D             = Property%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR00'


            call GetDefaultValue(Property%ID%ObjFillMatrix, Property%Scalar, STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR10'

            call GetIfMatrixRemainsConstant(FillMatrixID    = Property%ID%ObjFillMatrix,        &
                                            RemainsConstant = Property%Constant,                &
                                            STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR20'


            if(.not. Property%ID%SolutionFromFile)then
                call KillFillMatrix(Property%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR30'
            end if

            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, LocalClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR40'
            end if

        else

            write (*,*)'Block ',trim(BeginBlock),' ',trim(EndBlock),' not found'
            stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR60'

        endif
        
        
        select case(ExtractType)

            case(FromBlock)

                call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR70'

            case(FromBlockInBlock)
                
                call RewindBlock  (Me%ObjEnterData, LocalClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Read_Property_2D - ModuleInterfaceSedimentWater - ERR80'

            case default

        end select


    end subroutine Read_Property_2D
    
    
    !--------------------------------------------------------------------------

    
    subroutine SetSubModulesConstructor 

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL 

        !Begin-----------------------------------------------------------------

        call SetHydrodynamicManning(HydrodynamicID = Me%ObjHydrodynamic,                &
                                    Manning        = Me%Manning,                        &
                                    STAT           = STAT_CALL)
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'SetSubModulesConstructor - ModuleInterfaceSedimentWater - ERR10'

        call SetHydrodynamicChezy  (HydrodynamicID = Me%ObjHydrodynamic,                &
                                    Chezy          = Me%Chezy,                          &
                                    ChezyCoef      = Me%ChezyCoef,                      &
                                    STAT           = STAT_CALL)
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'SetSubModulesConstructor - ModuleInterfaceSedimentWater - ERR20'


        if(Me%Manning)then
            
            call SetHydrodynamicRugosityMatrix(HydrodynamicID = Me%ObjHydrodynamic,     &
                                               RugosityMatrix = Me%ManningCoef%Field,   &
                                               STAT           = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_)                                                 &
                stop 'SetSubModulesConstructor - ModuleInterfaceSedimentWater - ERR30'


        else

            call SetHydrodynamicRugosityMatrix(HydrodynamicID = Me%ObjHydrodynamic,     &
                                               RugosityMatrix = Me%Rugosity%Field,      &
                                               STAT           = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_)                                                 &
                stop 'SetSubModulesConstructor - ModuleInterfaceSedimentWater - ERR40'

        end if

        call SetTurbulenceBottomRugosity(TurbulenceID   = Me%ObjTurbulence,             &
                                         BottomRugosity = Me%Rugosity%Scalar,           &
                                         STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'SetSubModulesConstructor - ModuleInterfaceSedimentWater - ERR50'

    end subroutine SetSubModulesConstructor


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !----------------------------------------------------------------------
    
    
    subroutine Search_Property(PropertyX, PropertyXID, STAT)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer         ,           intent (IN)         :: PropertyXID
        integer         , optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber==PropertyXID) then
                exit        
            else
                PropertyX => PropertyX%Next                 
            end if    
        end do    

       if (associated(PropertyX)) then

            STAT_ = SUCCESS_  

        else
            STAT_  = NOT_FOUND_ERR_  
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine Search_Property



    
    !--------------------------------------------------------------------------
    
    
    subroutine ReadLockExternalGlobal
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceSedimentWater - ERR01'

        !XX_IE and YY_IE
        call GetHorizontalGrid (Me%ObjHorizontalGrid,                                   &
                                XX_IE = Me%ExternalVar%XX_IE,                           &
                                YY_IE = Me%ExternalVar%YY_IE,                           &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceSedimentWater - ERR02'

        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceSedimentWater - ERR03'

    end subroutine ReadLockExternalGlobal

    
    !--------------------------------------------------------------------------

    
    subroutine ReadLockExternalWater
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjWaterHorizontalMap, Me%ExtWater%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR01'

        !OpenPoints2D
        call GetOpenPoints2D(Me%ObjWaterHorizontalMap, Me%ExtWater%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR02'
        
        !WaterPoints3D
        call GetWaterPoints3D(Me%ObjWaterMap, Me%ExtWater%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR03'

        !OpenPoints3D
        call GetOpenPoints3D(Me%ObjWaterMap, Me%ExtWater%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR04'

        !BoundaryPoints2D
        call GetBoundaries(Me%ObjWaterHorizontalMap, Me%ExtWater%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR05'

        !LandPoints3D
        call GetLandPoints3D(Me%ObjWaterMap, Me%ExtWater%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR06'

        !SZZ
        call GetGeometryDistances (Me%ObjWaterGeometry, SZZ = Me%ExtWater%SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR08'


        call GetGeometryVolumes(Me%ObjWaterGeometry,                        &
                                VolumeZ    = Me%ExtWater%VolumeZ,           &
                                VolumeZOld = Me%ExtWater%VolumeZOld,        &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR09'

        call GetGeometryDistances (Me%ObjWaterGeometry, DWZ = Me%ExtWater%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR10'
        
       
        call GetChezy(HydrodynamicID = Me%ObjHydrodynamic,                  &
                      Chezy          = Me%ExtWater%Chezy,                   &
                      STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR12'
        
        call GetHorizontalVelocity(HydrodynamicID = Me%ObjHydrodynamic,     &
                                   Velocity_U     = Me%ExtWater%Velocity_U, &
                                   Velocity_V     = Me%ExtWater%Velocity_V, &
                                   STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceSedimentWater - ERR13'

        call GetGridData(Me%ObjWaterGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleWaterProperties - ERR14'

        call GetGeometryWaterColumn(Me%ObjWaterGeometry,                    &
                                    WaterColumn = Me%ExtWater%WaterColumn,  &
                                    STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleWaterProperties - ERR15'

        call GetGeometryMinWaterColumn(Me%ObjWaterGeometry,                          &
                                       MinWaterColumn = Me%ExtWater%MinWaterColumn,  &
                                       STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleWaterProperties - ERR15a'

        call GetGeometryKFloor(Me%ObjWaterGeometry, Z = Me%ExtWater%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleWaterProperties - ERR16'

#ifndef _WAVES_
        if (Me%ObjWaves /=0) then        

            call GetWaves (WavesID       = Me%ObjWaves,                                  &
                           WavePeriod    = Me%ExtWater%WavePeriod,                       &
                           WaveHeight    = Me%ExtWater%WaveHeight,                       &
                           Abw           = Me%ExtWater%Abw,                              &
                           Ubw           = Me%ExtWater%Ubw,                              &
                           LastCompute   = Me%ExtWater%LastComputeWave,                  &
                           STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterfaceSedimentWater - ERR17'

        endif
#endif


    end subroutine ReadLockExternalWater


    !--------------------------------------------------------------------------
#ifndef _SEDIMENT_
    subroutine ReadLockExternalSediment
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjSedimentHorizontalMap, Me%ExtSed%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR01'

        !OpenPoints2D
        call GetOpenPoints2D(Me%ObjSedimentHorizontalMap, Me%ExtSed%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR02'

        !WaterPoints3D
        call GetWaterPoints3D(Me%ObjSedimentMap, Me%ExtSed%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR03'

        !OpenPoints3D
        call GetOpenPoints3D(Me%ObjSedimentMap, Me%ExtSed%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR04'


        !BoundaryPoints2D
        call GetBoundaries(Me%ObjSedimentHorizontalMap, Me%ExtSed%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR05'


        !SZZ
        call GetGeometryDistances (Me%ObjSedimentGeometry, SZZ = Me%ExtSed%SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR07'
        
        !DWZ
        call GetGeometryDistances (Me%ObjSedimentGeometry, DWZ = Me%ExtSed%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR08'
        
        !Volumes
        call GetGeometryVolumes(Me%ObjSedimentGeometry,                             &
                                VolumeZ     = Me%ExtSed%VolumeZ,                    & 
                                VolumeZOld  = Me%ExtSed%VolumeZOld,                 & 
                                STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR09'

        !LandPoints3D
        call GetLandPoints3D(Me%ObjSedimentMap, Me%ExtSed%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR10'


        !WaterFluxes
        call GetConsolidationWaterFluxes(Me%ObjConsolidation,                       & 
                                         WaterFluxZ = Me%ExtSed%WaterFluxZ,         &
                                         STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR11'

        !TopCriticalShear
        call GetConsolidationCriticalShear(Me%ObjConsolidation,                     & 
                                           TopCriticalShear = Me%ExtSed%TopCriticalShear,&
                                           STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR12'


        call GetConsolidationDrySedVolume(Me%ObjConsolidation,                          & 
                                          DrySedimentVolume = Me%ExtSed%DrySedVolume,   &
                                          STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR13'


        call GetSedimentDryDensity(Me%ObjSedimentProperties,                            &
                                   SedimentDryDensity = Me%ExtSed%SedimentDryDensity,   &
                                   STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR14'

        call GetGeometryKTop(Me%ObjSedimentGeometry,                                    &
                             KTopZ  = Me%ExtSed%KTop,                                   &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR15'

        call GetConsolidationPorosity(Me%ObjConsolidation,                              & 
                                      Porosity = Me%ExtSed%Porosity,                    &
                                      STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR16'

        call GetSedimentColumnFull(Me%ObjConsolidation,                                 & 
                                   SedimentColumnFull = Me%ExtSed%SedimentColumnFull,   &
                                   STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadLockExternalSediment - ModuleInterfaceSedimentWater - ERR17'
        

    end subroutine ReadLockExternalSediment

#endif
    !----------------------------------------------------------------------


    subroutine ReadUnlockExternalGlobal
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalGlobal - ModuleInterfaceSedimentWater - ERR01'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalGlobal - ModuleInterfaceSedimentWater - ERR03'
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalGlobal - ModuleInterfaceSedimentWater - ERR03'
        

        call null_time(Me%ExternalVar%Now)

    end subroutine ReadUnlockExternalGlobal

    
    !----------------------------------------------------------------------


    subroutine ReadUnlockExternalWater
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call UnGetHorizontalMap(Me%ObjWaterHorizontalMap, Me%ExtWater%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR01'

        call UnGetHorizontalMap(Me%ObjWaterHorizontalMap, Me%ExtWater%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR02'

        call UngetHorizontalMap (Me%ObjWaterHorizontalMap, Me%ExtWater%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR03'

        call UnGetMap(Me%ObjWaterMap, Me%ExtWater%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR04'
        
        call UnGetMap(Me%ObjWaterMap, Me%ExtWater%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR05'

        call UnGetMap(Me%ObjWaterMap, Me%ExtWater%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR06'

        call UnGetGeometry(Me%ObjWaterGeometry,Me%ExtWater%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR07'
        
        call UnGetGeometry(Me%ObjWaterGeometry,Me%ExtWater%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR08'

        call UnGetGeometry(Me%ObjWaterGeometry, Me%ExtWater%VolumeZ, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR09'

        call UnGetGeometry(Me%ObjWaterGeometry, Me%ExtWater%VolumeZOld, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR10'
                               
        call UnGetHydrodynamic(Me%ObjHydrodynamic, Me%ExtWater%Chezy, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR12'

        call UnGetHydrodynamic(Me%ObjHydrodynamic, Me%ExtWater%Velocity_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR13'

        call UnGetHydrodynamic(Me%ObjHydrodynamic, Me%ExtWater%Velocity_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR14'

        call UnGetGridData(Me%ObjWaterGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR15'

        call UnGetGeometry(Me%ObjWaterGeometry, Me%ExtWater%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR16'

        call UnGetGeometry(Me%ObjWaterGeometry, Me%ExtWater%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR17'
        

#ifndef _WAVES_
        if (Me%ObjWaves /=0) then

            call UnGetWaves(Me%ObjWaves, Me%ExtWater%WavePeriod, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR18'

            call UnGetWaves(Me%ObjWaves, Me%ExtWater%WaveHeight, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR19'

            call UnGetWaves(Me%ObjWaves, Me%ExtWater%Abw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR20'

            call UnGetWaves(Me%ObjWaves, Me%ExtWater%Ubw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceSedimentWater - ERR21'

        endif
#endif


    end subroutine ReadUnlockExternalWater


    !----------------------------------------------------------------------

#ifndef _SEDIMENT_
    subroutine ReadUnlockExternalSediment
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalMap(Me%ObjSedimentHorizontalMap, Me%ExtSed%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR01'

        call UnGetHorizontalMap(Me%ObjSedimentHorizontalMap, Me%ExtSed%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR02'

        call UnGetMap(Me%ObjSedimentMap, Me%ExtSed%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR03'
        
        call UnGetMap(Me%ObjSedimentMap, Me%ExtSed%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR04'

        call UnGetMap(Me%ObjSedimentMap, Me%ExtSed%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR05'

        call UnGetGeometry(Me%ObjSedimentGeometry,Me%ExtSed%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR06'
        
        call UnGetGeometry(Me%ObjSedimentGeometry,Me%ExtSed%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR07'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExtSed%WaterFluxZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR13'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExtSed%TopCriticalShear, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR14'

        call UnGetGeometry(Me%ObjSedimentGeometry, Me%ExtSed%VolumeZ, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR15'

        call UnGetGeometry(Me%ObjSedimentGeometry, Me%ExtSed%VolumeZOld, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR16'

        call UngetHorizontalMap (Me%ObjSedimentHorizontalMap, Me%ExtSed%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR17'
                 
        call UngetConsolidation(Me%ObjConsolidation, Me%ExtSed%DrySedVolume, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR18'  
         
        call UngetSedimentProperties(Me%ObjConsolidation, Me%ExtSed%SedimentDryDensity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR19'  

        call UnGetGeometry(Me%ObjSedimentGeometry, Me%ExtSed%KTop, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR20'  

        call UngetConsolidation(Me%ObjConsolidation, Me%ExtSed%Porosity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR21'  
                       
        call UngetConsolidation(Me%ObjConsolidation, Me%ExtSed%SedimentColumnFull, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalSediment - ModuleInterfaceSedimentWater - ERR22'  

    end subroutine ReadUnlockExternalSediment

#endif
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyInterfaceSedimentWater(ObjInterfaceSedimentWaterID, LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceSedimentWaterID, LagrangianID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterfaceSedimentWaterID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if(Me%ObjLagrangian == 0 .and. LagrangianID /= 0)then
                Me%ObjLagrangian  = AssociateInstance(mLAGRANGIAN_, LagrangianID)
            endif

            if (MonitorPerformance)                 &
                call StartWatch ("ModuleInterfaceSedimentWater", "ModifyInterfaceSedimentWater")

            call ReadLockExternalGlobal

            call ReadLockExternalWater

#ifndef _SEDIMENT_
            if(Me%RunsSediments) call ReadLockExternalSediment
#endif

            call TimeStepActualization

            call ModifyShearStress
            
            if(Me%Coupled%SedimentFluxes%Yes)       &
                call ModifySedimentColumnFluxes

            if(Me%Coupled%WaterFluxes%Yes)          &
                call ModifyWaterColumnFluxes

            if(Me%Coupled%Detritus%Yes)             &
                call Detritus_Processes

            if(Me%Coupled%Benthos%Yes)              &
                call Benthos_Processes
                
            if(Me%Coupled%BenthicEcology%Yes)       &
                call BenthicEcology_Processes  
                                   
            if(Me%Coupled%CEQUALW2%Yes)             &
                call CEQUALW2_Processes
            
            if(Me%Coupled%SedimentWaterFluxes%Yes)  &
                call ModifySedimentWaterFluxes

            if (Me%RunsSandTransport)               &
                call ModifySandTransport

            if(Me%Coupled%TimeSerie%Yes)            &
                call Output_TimeSeries

            if(Me%Coupled%BoxTimeSerie%Yes)         &
                call Output_BoxTimeSeries

            if(Me%Coupled%OutputHDF%Yes)            &
                call OutPut_Results_HDF

            if(Me%Output%WriteRestartFile)          &
                call OutputRestartFile

            call Actualize_Time_Evolution

#ifndef _SEDIMENT_
            if(Me%RunsSediments) call ReadUnlockExternalSediment
#endif
            call ReadUnlockExternalWater

            call ReadUnlockExternalGlobal

            call SetSubModulesModifier

            if (MonitorPerformance)                 &
                call StopWatch ("ModuleInterfaceSedimentWater", "ModifyInterfaceSedimentWater")



            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyInterfaceSedimentWater

    
    !--------------------------------------------------------------------------
    
    subroutine TimeStepActualization

        !Local-----------------------------------------------------------------
        type (T_Property), pointer       :: Property
        real                             :: NewDT
        integer                          :: STAT_CALL
        logical                          :: VariableDT

        !Begin-----------------------------------------------------------------


        call GetVariableDT (Me%ObjTime, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'TimeStepActualization - ModuleInterfaceSedimentWater - ERR01'
    
cd1:    if (VariableDT) then

            call GetComputeTimeStep(Me%ObjTime, NewDT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'TimeStepActualization - ModuleInterfaceSedimentWater - ERR02'

            Property => Me%FirstProperty

do1 :       do while (associated(Property))

                Property%Evolution%DTInterval = NewDT

                Property => Property%Next

            enddo do1

            nullify(Property)

        endif cd1


    end subroutine TimeStepActualization
   
    !------------------------------------------------------------------------
        
    
    subroutine ModifyShearStress

        !Local-----------------------------------------------------------------
        real,    pointer, dimension(:, :, :)    :: Velocity_U, Velocity_V
        real,    pointer, dimension(:, :   )    :: Chezy
        integer, pointer, dimension(:, :, :)    :: LandPoints3D, OpenPoints3D
        integer, pointer, dimension(:, :   )    :: KFloorZ
        real                                    :: VC,UC, UVC2, WaterDensity
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j, kbottom
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------

        !A if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ModifyShearStress")
                

        if (Me%Shear_Stress%LastCompute .LT. Me%ExternalVar%Now) then

            IUB = Me%WaterWorkSize3D%IUB
            JUB = Me%WaterWorkSize3D%JUB
            ILB = Me%WaterWorkSize3D%ILB
            JLB = Me%WaterWorkSize3D%JLB
            KUB = Me%WaterWorkSize3D%KUB

            Velocity_U      => Me%ExtWater%Velocity_U
            Velocity_V      => Me%ExtWater%Velocity_V
            Chezy           => Me%ExtWater%Chezy
            WaterDensity    =  SigmaDensityReference
            LandPoints3D    => Me%ExtWater%LandPoints3D
            OpenPoints3D    => Me%ExtWater%OpenPoints3D
            KFloorZ         => Me%ExtWater%KFloor_Z

            if (Me%WaveShear_Stress%Yes) then

                if (Me%ExtWater%LastComputeWave > Me%WaveShear_Stress%LastCompute) then


                    if (.not. Me%WaveShear_Stress%Rugosity%Constant) then
                        call ComputeWaveRugosity
                    endif
                    call ComputeWaveTension

                    Me%WaveShear_Stress%LastCompute = Me%ExternalVar%Now

                endif

            end if

            if (MonitorPerformance) then
                call StartWatch ("ModuleInterfaceSedimentWater", "ModifyShearStress")
            endif

            CHUNK = CHUNK_J(JLB, JUB)
            
            !$OMP PARALLEL PRIVATE(i,j,kbottom,VC,UC,UVC2)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            ! Se alisa cantos com terra. Aumenta artificialmente a tensão de corte 
            ! nos cantos com terra para evitar que estes se tornem em zonas 
            ! de deposição acentuada
            do j = JLB, JUB
            do i = ILB, IUB
                                        
                if (Me%ExtWater%OpenPoints3D(i, j, KUB) == OpenPoint) then

                    kbottom = KFloorZ(i, j)

                    if (((LandPoints3D(i,  j+1, kbottom)==1 .or. LandPoints3D(i,  j-1, kbottom)==1) .and.   &
                         (LandPoints3D(i+1,j,   kbottom)==1 .or. LandPoints3D(i-1,j,   kbottom)==1))) then

                        VC = max(abs(Velocity_V(i+1,j,kbottom)),abs(Velocity_V(i,j,  kbottom)))
                        UC = max(abs(Velocity_U(i,  j,kbottom)),abs(Velocity_U(i,j+1,kbottom)))

                    else

                        VC = (Velocity_V(i+1,j,kbottom)+Velocity_V(i,j,  kbottom))/2.
                        UC = (Velocity_U(i,  j,kbottom)+Velocity_U(i,j+1,kbottom))/2.

                    endif

                    UVC2 = UC*UC+VC*VC


                    Me%Shear_Stress%Tension (i,j) = Chezy(i,j) * UVC2 * WaterDensity

                    if (Me%RunsSandTransport) then

                        Me%Shear_Stress%CurrentVel(i, j) = sqrt(UVC2)
                        Me%Shear_Stress%CurrentU  (i, j) = UC
                        Me%Shear_Stress%CurrentV  (i, j) = VC

                    endif

                    if (Me%WaveShear_Stress%Yes) then

                        if (Me%WaveShear_Stress%NonLinear) then

                            Me%WaveShear_Stress%Tension (i,j) = Me%WaveShear_Stress%ChezyVel   (i,j) * &
                                                                sqrt(UVC2) * WaterDensity


                        else

                            Me%WaveShear_Stress%Tension (i,j) = Me%WaveShear_Stress%ChezyVel   (i,j) * &
                                                                Me%ExtWater%Ubw (i,j) * WaterDensity
                        endif

                        !Currents contribution 
                        if (Me%RunsSandTransport) then
                            Me%WaveShear_Stress%TensionCurrents(i, j) = Me%Shear_Stress%Tension (i,j)
                        endif

                        !Wave contribution
                        Me%Shear_Stress%Tension(i, j) = Me%Shear_Stress%Tension (i,j) + Me%WaveShear_Stress%Tension(i, j)

                    endif

                    if (Me%Shear_Stress%IntertidalRunOff) then

                        if (Me%ExtWater%WaterColumn(i,j) < 2. * Me%ExtWater%MinWaterColumn) then

                            Me%Shear_Stress%Tension(i, j) = Me%Shear_Stress%Tension (i,j) + 1.

                        endif

                    endif 

                    ! [m/s]                       = [N/m^2/ (kg/m^3)]^0.5 = 
                    Me%Shear_Stress%Velocity(i,j) = sqrt(Me%Shear_Stress%Tension (i,j) / WaterDensity)   

                endif

            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            if (MonitorPerformance) then
                call StopWatch ("ModuleInterfaceSedimentWater", "ModifyShearStress")
            endif

            if (Me%Shear_Stress%Statistics%ON) then

                call OutPut_Statistics (Me%Shear_Stress%Velocity, Me%Shear_Stress%Statistics%ID)
                                        
            endif

     
            Me%Shear_Stress%LastCompute = Me%ExternalVar%Now

        
        endif

        !A if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ModifyShearStress")

    
    end subroutine ModifyShearStress


    !--------------------------------------------------------------------------

    subroutine ComputeWaveRugosity

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer :: D50Field2D, D90Field2D
        real                          :: Roughness, WaterDensity, Abw, Ubw
        real                          :: D50,D90,SedimentDryDensity,RelativeDensity
        real                          :: Klsw, Kllsw, psi, RippleHeight, RippleLength
        integer                       :: IUB, JUB, ILB, JLB, i, j
        real                          :: WaveHeight, WaterColumn, WavePeriod, LimitMin
        integer                       :: STAT_CALL
        integer                       :: CHUNK
        
        !Begin-----------------------------------------------------------------


        if (Me%RunsSandTransport) then

            call GetSandDiameters(ObjSandID = Me%ObjSand,                                &
                                  D50       = D50Field2D,                                &
                                  D90       = D90Field2D,                                & 
                                  STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveRugosity - ModuleInterfaceSedimentWater - ERR02'

            call GetSandDensity (ObjSandID   = Me%ObjSand,                               &
                                 SandDensity = SedimentDryDensity,                       & 
                                 STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveRugosity - ModuleInterfaceSedimentWater - ERR03'

        else
               
            D50                 = 2E-3
            D90                 = 3E-3
            SedimentDryDensity  = 2300.

        endif


        IUB = Me%WorkSize2D%IUB
        JUB = Me%WorkSize2D%JUB
        ILB = Me%WorkSize2D%ILB
        JLB = Me%WorkSize2D%JLB

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceSedimentWater", "ComputeWaveRugosity")
        endif

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,WavePeriod,WaveHeight,Abw,Ubw,WaterColumn,WaterDensity) &
        !$OMP PRIVATE(RelativeDensity,Psi,LimitMin,Klsw,Kllsw,RippleHeight,RippleLength) &
        !$OMP PRIVATE(Roughness) &
        !$OMP FIRSTPRIVATE(D50,D90)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
                                    
cd0:        if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                WavePeriod      = Me%ExtWater%WavePeriod(i,j)
                WaveHeight      = Me%ExtWater%WaveHeight(i,j)
                Abw             = Me%ExtWater%Abw(i,j)
                Ubw             = Me%ExtWater%Ubw(i,j)
                WaterColumn     = Me%ExtWater%WaterColumn(i,j)

                WaterDensity    = SigmaDensityReference
        
cd1:            if (WavePeriod .ne. 0. .or. WaveHeight .ne. 0.) then
        
                    if (Me%RunsSandTransport) then

                        D50 = D50Field2D(i, j)
                        D90 = D90Field2D(i, j)

                    endif

       
                    RelativeDensity     = (SedimentDryDensity - WaterDensity) / &
                                           WaterDensity
        
cd6:                if(WaveHeight.LT.0.10)then
        
                        Me%WaveShear_Stress%Rugosity%Field(i,j) = 3. * D90
        
                    else cd6
        
! ---> Computes Ripple Height and Lenght (Eq. 6.3.7 and 6.3.8 VanRijn)
        
                        Psi = Ubw**2./(RelativeDensity*Gravity*D50)

                        LimitMin = 1e-6
        
cd4:                    if(Psi .LT. LimitMin**2)then

                            Klsw            = 3.  * D90 
                            Kllsw           = 0.

                        elseif(Psi.LE.10 .and. Psi .GT. LimitMin**2)then
        
                            RippleHeight    = Abw * 0.2
                            RippleLength    = RippleHeight / 0.18
                            Klsw            = 3.  * D90 
                            Kllsw           = 16. * RippleHeight * (RippleHeight / RippleLength)
        
                        elseif (Psi .GT. 10. .AND. Psi .LE. 250.)then cd4
        
                            RippleHeight    = Abw * 2.8E-13 * (250   - Psi)**5.
                            RippleLength    = RippleHeight  / (2.E-7 * (250 - Psi)**2.5)
                            Klsw            = 3.  * D90
                            Kllsw           = 16. * RippleHeight * (RippleHeight / RippleLength)
        
                        else cd4
        
                            RippleHeight    = 0.
                            RippleLength    = 0.
                            Klsw            = 3.  * (0.04 * Psi - 9.) * D90 
                            Kllsw           = 0.
        
                        end if cd4
        
! ---> Computes Wave Related bed BedRoughnesshness
        
                        Roughness                         = Klsw+Kllsw
                        Me%WaveShear_Stress%Rugosity%Field(i,j) = Roughness
                        !Me%WaveShear_Stress%Rugosity%Field(i,j) = min(Roughness, 0.1)
        
                    end if  cd6  ! WaveHeight>0

                endif cd1

            endif cd0
        

        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL        

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceSedimentWater", "ComputeWaveRugosity")
        endif

        if (Me%RunsSandTransport) then

            call UnGetSand(Me%ObjSand, D50Field2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveRugosity - ModuleInterfaceSedimentWater - ERR09'

            call UnGetSand(Me%ObjSand, D90Field2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveRugosity - ModuleInterfaceSedimentWater - ERR10'

        endif


    end subroutine ComputeWaveRugosity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ComputeWaveTension

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                          :: fw, WaterDensity, Abw, Ubw
        real                          :: C1
        integer                       :: IUB, JUB, ILB, JLB, i, j
        real                          :: WaveHeight, WavePeriod, LimitMin
        integer                       :: CHUNK
        
        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ComputeWaveTension")


        IUB = Me%WorkSize2D%IUB
        JUB = Me%WorkSize2D%JUB
        ILB = Me%WorkSize2D%ILB
        JLB = Me%WorkSize2D%JLB

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,WavePeriod,WaveHeight,Abw,Ubw,WaterDensity,LimitMin,C1) &
        !$OMP PRIVATE(Fw)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
                                    
cd0:        if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                WavePeriod      = Me%ExtWater%WavePeriod(i,j)
                WaveHeight      = Me%ExtWater%WaveHeight(i,j)
                Abw             = Me%ExtWater%Abw(i,j)
                Ubw             = Me%ExtWater%Ubw(i,j)

                WaterDensity    = SigmaDensityReference
        
cd1:            if (WavePeriod .ne. 0. .or. WaveHeight .ne. 0.) then
        
! ---> Finally WaveShearStress           
        
                    LimitMin = 1.e-6

                    C1       = -0.194
        
cd7:                if(WaveHeight .GT. 0.05 .and. Abw > LimitMin)then

                        Fw = exp(-5.977 + 5.213 * (Abw / Me%WaveShear_Stress%Rugosity%Field(i,j))**C1)
    
                        if (Fw.GT.0.3) Fw = 0.3
    
                    else cd7
    
                        Fw = 0.
    
                    end if cd7

                    Me%WaveShear_Stress%ChezyVel  (i,j)    = 0.25 * Fw * Ubw

    
                else cd1

                    Me%WaveShear_Stress%ChezyVel  (i,j)    = 0.
    
                endif cd1


            end if cd0

        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ComputeWaveTension")

    end subroutine ComputeWaveTension

    !--------------------------------------------------------------------------


    subroutine ModifyWaterColumnFluxes

        
        
        
        if (Me%Coupled%BenthicEcology%Yes)then
        call InitializeFluxesToWaterColumn_Benthic
        else
        call InitializeFluxesToWaterColumn
        endif
        
        
        if(Me%Coupled%Erosion%Yes)then
            
            !Proportional factor between property mass and cohesive sediment mass
            call ModifyErosionCoefficient

            call ModifyErosionFluxes

        end if

        if(Me%Coupled%Deposition%Yes)then

            call ModifyDepositionFluxes
        
        endif

        call ModifyDissolvedFluxes

    end subroutine ModifyWaterColumnFluxes

    
    !--------------------------------------------------------------------------
    
    
    subroutine ModifySedimentColumnFluxes

#ifndef _SEDIMENT_

        call InitializeFluxesToSediment

        if(Me%Consolidation%Yes)then

            call ComputeConsolidation

        end if

#endif


    end subroutine ModifySedimentColumnFluxes

    
    !--------------------------------------------------------------------------
    
    
    subroutine ModifySedimentWaterFluxes

#ifndef _SEDIMENT_

        if(Me%Consolidation%Yes)then

            call ModifyConsolidatedErosionFluxes

            call ComputeWaterFlux

        end if


        call ParticulateSedimentWaterFluxes

        call DissolvedSedimentWaterFluxes




#endif

    end subroutine ModifySedimentWaterFluxes


    !--------------------------------------------------------------------------

    subroutine InitializeFluxesToWaterColumn
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB
        integer                                 :: i, j, kbottom
        type(T_Property), pointer               :: PropertyX
        real, dimension(:,:,:), pointer         :: WaterPropertyConcentration
        character(len=StringLength)             :: WaterPropertyUnits
        real                                    :: WaterPropertyISCoef
        integer                                 :: CHUNK 

        !Begin-----------------------------------------------------------------

        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then

                if (MonitorPerformance) then
                    call StartWatch ("ModuleInterfaceSedimentWater", "InitializeFluxesToWaterColumn")
                endif

                if(.not. PropertyX%Particulate)then

                    call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                                          ConcentrationX    = WaterPropertyConcentration,       &
                                          PropertyXIDNumber = PropertyX%ID%IDNumber,            &
                                          PropertyXUnits    = WaterPropertyUnits,               &
                                          PropertyXISCoef   = WaterPropertyISCoef,              &
                                          STAT              = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'InitializeFluxesToWaterColumn - ModuleInterfaceSedimentWater - ERR01'

                    CHUNK = CHUNK_I(ILB, IUB)
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do i = ILB, IUB
                    do j = JLB, JUB

                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            PropertyX%WaterConcentration(i,j) = WaterPropertyConcentration(i,j,kbottom) * &
                                                                WaterPropertyISCoef

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                    call UnGetWaterProperties(Me%ObjWaterProperties, WaterPropertyConcentration, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'InitializeFluxesToWaterColumn - ModuleInterfaceSedimentWater - ERR02'
                end if


                if(PropertyX%Evolution%WaterFluxes .or. PropertyX%Evolution%SedimentWaterFluxes)then

                    CHUNK = CHUNK_I(ILB, IUB)
                    !$OMP PARALLEL PRIVATE(i,j)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do i = ILB, IUB
                    do j = JLB, JUB
                
                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            PropertyX%FluxToWater(i,j) = 0.

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                end if

                if (MonitorPerformance) then
                    call StopWatch ("ModuleInterfaceSedimentWater", "InitializeFluxesToWaterColumn")
                endif

            end if

            PropertyX => PropertyX%Next

        end do


    
    end subroutine InitializeFluxesToWaterColumn


    !--------------------------------------------------------------------------

    subroutine InitializeFluxesToWaterColumn_Benthic
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB
        integer                                 :: i, j, kbottom
        type(T_Property), pointer               :: PropertyX
        real, dimension(:,:,:), pointer         :: WaterPropertyConcentration
        character(len=StringLength)             :: WaterPropertyUnits
        real                                    :: WaterPropertyISCoef
        integer                                 :: CHUNK 

        !Begin-----------------------------------------------------------------

        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then

                if (MonitorPerformance) then
                    call StartWatch ("ModuleInterfaceSedimentWater", "InitializeFluxesToWaterColumn_Benthic")
                endif

                !if(.not. PropertyX%Particulate)then
                if(.not. PropertyX%Evolution%BenthicOnly)then
                    call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                                          ConcentrationX    = WaterPropertyConcentration,       &
                                          PropertyXIDNumber = PropertyX%ID%IDNumber,            &
                                          PropertyXUnits    = WaterPropertyUnits,               &
                                          PropertyXISCoef   = WaterPropertyISCoef,              &
                                          STAT              = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'InitializeFluxesToWaterColumn_Benthic - ModuleInterfaceSedimentWater - ERR01'

                    CHUNK = CHUNK_I(ILB, IUB)
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do i = ILB, IUB
                    do j = JLB, JUB

                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            PropertyX%WaterConcentration(i,j) = WaterPropertyConcentration(i,j,kbottom) * &
                                                                WaterPropertyISCoef
                                                                
                            
                            PropertyX%Mass_FromWater(i,j) =  PropertyX%WaterConcentration(i,j) * &
                                                                 Me%ExtWater%VolumeZ(i,j,kbottom)                                    

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                    call UnGetWaterProperties(Me%ObjWaterProperties, WaterPropertyConcentration, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'InitializeFluxesToWaterColumn_Benthic - ModuleInterfaceSedimentWater - ERR02'
                

                
                end if


                if(PropertyX%Evolution%WaterFluxes .or. PropertyX%Evolution%SedimentWaterFluxes)then

                    CHUNK = CHUNK_I(ILB, IUB)
                    !$OMP PARALLEL PRIVATE(i,j)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do i = ILB, IUB
                    do j = JLB, JUB
                
                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            PropertyX%FluxToWater(i,j) = 0.

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                end if

                if (MonitorPerformance) then
                    call StopWatch ("ModuleInterfaceSedimentWater", "InitializeFluxesToWaterColumn_Benthic")
                endif

            end if

            PropertyX => PropertyX%Next

        end do


    
    end subroutine InitializeFluxesToWaterColumn_Benthic

    !--------------------------------------------------------------------------

    subroutine InitializeFluxesToSediment
        
#ifndef _SEDIMENT_

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j
        type(T_Property), pointer               :: PropertyX
        real, dimension(:,:,:), pointer         :: SedimentPropertyConcentration
        character(len=StringLength)             :: SedimentPropertyUnits
        real                                    :: SedimentPropertyISCoef
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------

        IUB = Me%SedimentWorkSize3D%IUB
        JUB = Me%SedimentWorkSize3D%JUB
        ILB = Me%SedimentWorkSize3D%ILB
        JLB = Me%SedimentWorkSize3D%JLB

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then

                if(PropertyX%Evolution%SedimentFluxes .or. PropertyX%Evolution%SedimentWaterFluxes)then

                    call GetSedimentConcentration(SedimentPropertiesID = Me%ObjSedimentProperties,  &
                                                  ConcentrationX    = SedimentPropertyConcentration,&
                                                  PropertyXIDNumber = PropertyX%ID%IDNumber,        &
                                                  PropertyXUnits    = SedimentPropertyUnits,        &
                                                  PropertyXISCoef   = SedimentPropertyISCoef,       &
                                                  STAT              = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'InitializeFluxesToSediment - ModuleInterfaceSedimentWater - ERR01'

                    if (MonitorPerformance) then
                        call StartWatch ("ModuleInterfaceSedimentWater", "InitializeFluxesToSediment")
                    endif

                    CHUNK = CHUNK_J(JLB, JUB)
                    !$OMP PARALLEL PRIVATE(i,j,KUB)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)    
                    do j = JLB, JUB
                    do i = ILB, IUB
            
                        if(Me%ExtSed%WaterPoints2D(i,j) == WaterPoint)then

                            KUB = Me%ExtSed%KTop(i,j)

                            PropertyX%SedimentConcentration(i,j) = SedimentPropertyConcentration(i,j,KUB) * &
                                                                   SedimentPropertyISCoef

                        end if

                    enddo
                    enddo
                    !$OMP END DO

                    !$OMP MASTER
                    call UnGetSedimentProperties(Me%ObjSedimentProperties, SedimentPropertyConcentration, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                            stop 'InitializeFluxesToSediment - ModuleInterfaceSedimentWater - ERR02'
                    !$OMP END MASTER

                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)    
                    do j = JLB, JUB
                    do i = ILB, IUB
                
                        if(Me%ExtSed%WaterPoints2D(i,j) == WaterPoint)then

                            PropertyX%FluxToSediment(i,j) = 0.

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
                    
                    if (MonitorPerformance) then
                        call StopWatch ("ModuleInterfaceSedimentWater", "InitializeFluxesToSediment")
                    endif
                    
                end if

            end if

            PropertyX => PropertyX%Next

        end do

#endif
    
    end subroutine InitializeFluxesToSediment


    !--------------------------------------------------------------------------


    subroutine ModifyErosionCoefficient

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB
        integer                                 :: i, j
        integer                                 :: STAT_CALL
        type(T_Property), pointer               :: PropertyX, CohesiveSediment
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------

        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB

        call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ModifyErosionCoefficient - ModuleInterfaceSedimentWater - ERR01'

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%Erosion)then

                if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then
                    
                    if (MonitorPerformance) then
                        call StartWatch ("ModuleInterfaceSedimentWater", "ModifyErosionCoefficient")
                    endif
                    
                    CHUNK = CHUNK_J(JLB, JUB)
                    !$OMP PARALLEL PRIVATE(i,j)
                    if(PropertyX%Mass_Limitation)then
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = JLB, JUB
                        do i = ILB, IUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                if(CohesiveSediment%Mass_Available(i,j) > 0.) then

                                    PropertyX%ErosionCoefficient(i,j) = Me%ErosionRate%Field(i,j)      * &
                                                                        PropertyX%Mass_Available (i,j) / &
                                                                        CohesiveSediment%Mass_Available (i,j)
                                else

                                    PropertyX%ErosionCoefficient(i,j) = 0.

                                end if

                            end if
                        enddo
                        enddo
                        !$OMP END DO

                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = JLB, JUB
                        do i = ILB, IUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                PropertyX%ErosionCoefficient(i,j) = Me%ErosionRate%Field(i,j)
                            else

                                PropertyX%ErosionCoefficient(i,j) = 0.

                            end if

                        enddo
                        enddo
                        !$OMP END DO
                    end if
                    !$OMP END PARALLEL

                    if (MonitorPerformance) then
                        call StopWatch ("ModuleInterfaceSedimentWater", "ModifyErosionCoefficient")
                    endif

                end if
            end if

            PropertyX => PropertyX%Next

        end do

        nullify(CohesiveSediment)
    
    end subroutine ModifyErosionCoefficient

    
    !--------------------------------------------------------------------------
    
    
    subroutine ModifyErosionFluxes

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB
        integer                                 :: i, j
        type(T_Property), pointer               :: PropertyX
        real                                    :: ShearStress
        real                                    :: PotentialNewMass, ReallyErodedMass

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ModifyErosionFluxes")

        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB


        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%Erosion)then

                if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then
                    
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                            !Limit shear stress values in small depths
                            if(Me%Shear_Stress%Limitation)then

                                ShearStress = ShearStressLimitation(Me%ExtWater%WaterColumn(i,j),&
                                                                    Me%Shear_Stress%Tension(i,j))
                            else

                                ShearStress = Me%Shear_Stress%Tension(i,j)

                            end if

                            if(PropertyX%Mass_Limitation)then

    
                                ![kg/m2s]
                                PropertyX%ErosionFlux(i,j) =                                                 &
                                    ErosionFlux(CriticalShearErosion = Me%Critical_Shear_Erosion%Field(i,j), &
                                                ShearStress          = ShearStress,                          &
                                                ErosionRate          = PropertyX%ErosionCoefficient(i,j),    &
                                                Mass_Available       = PropertyX%Mass_Available(i,j),        &
                                                MinimumMass          = PropertyX%Mass_Min)


                                ![kg/m2] = [kg/m2] - [kg/m2/s * s]
                                PotentialNewMass                     = PropertyX%Mass_Available(i,j)       - &
                                                                       PropertyX%ErosionFlux(i,j)          * &
                                                                       PropertyX%Evolution%DTInterval
                                
                                !if erosion flux will erode more than it exists
                                if(PotentialNewMass <= PropertyX%Mass_Min)then
                                    
                                    !erosion will happen until mass minimum
                                    ReallyErodedMass                = PropertyX%Mass_Available(i,j) - PropertyX%Mass_Min
                                    
                                    !Re-compute real erosion flux
                                    PropertyX%ErosionFlux(i,j)      = ReallyErodedMass / PropertyX%Evolution%DTInterval

                                    PropertyX%Mass_Available(i,j)   = PropertyX%Mass_Min

                                else

                                    PropertyX%Mass_Available(i,j)   = PotentialNewMass

                                end if

                            else

                                ![kg/m2s]
                                PropertyX%ErosionFlux(i,j) =                                                 &
                                    ErosionFlux(CriticalShearErosion = Me%Critical_Shear_Erosion%Field(i,j), &
                                                ShearStress          = ShearStress,                          &
                                                ErosionRate          = PropertyX%ErosionCoefficient(i,j))
                            end if


                            !Sum erosion flux to the flux to water
                            PropertyX%FluxToWater(i,j) = PropertyX%FluxToWater(i,j) + PropertyX%ErosionFlux(i,j)

                        end if

                    enddo
                    enddo

                end if

            end if

            PropertyX => PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ModifyErosionFluxes")

    end subroutine ModifyErosionFluxes


    !--------------------------------------------------------------------------

    real function ErosionFlux (CriticalShearErosion, ShearStress, ErosionRate, &
                               Mass_Available, MinimumMass) 
        
        !Arguments-------------------------------------------------------------
        real, intent(IN)                        :: CriticalShearErosion
        real, intent(IN)                        :: ShearStress
        real, intent(IN)                        :: ErosionRate
        real, intent(IN), optional              :: Mass_Available
        real, intent(IN), optional              :: MinimumMass

        !Begin-----------------------------------------------------------------


        if (present(Mass_Available).and.present(MinimumMass)) then
            
            if (ShearStress >= CriticalShearErosion .and. Mass_Available > MinimumMass)then
                
                ErosionFlux = ErosionRate * (ShearStress / CriticalShearErosion - 1.)
            
            else
                
                ErosionFlux = 0.
            endif

        else

            if (ShearStress >= CriticalShearErosion) then

                ErosionFlux = ErosionRate * (ShearStress / CriticalShearErosion - 1.)

            else

                ErosionFlux = 0.

            endif

        endif

    end function ErosionFlux

    !--------------------------------------------------------------------------

    real function ShearStressLimitation(WaterColumn, ComputedShearStress)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                        :: WaterColumn
        real, intent(IN)                        :: ComputedShearStress
        
        !Local-----------------------------------------------------------------
        real                                    :: ReferenceDepth
        real                                    :: ReferenceShearStress

        !Begin-----------------------------------------------------------------

        ReferenceDepth       = Me%Shear_Stress%ReferenceDepth
        ReferenceShearStress = Me%Shear_Stress%ReferenceShearStress

        if(WaterColumn < ReferenceDepth) then

            ShearStressLimitation = min(ComputedShearStress, ReferenceShearStress)
        else
            ShearStressLimitation = ComputedShearStress
        
        end if

    end function ShearStressLimitation

    !--------------------------------------------------------------------------

    subroutine ModifyConsolidatedErosionFluxes
#ifndef _SEDIMENT_


        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j, STAT_CALL
        type(T_Property), pointer               :: CohesiveSediment
        real                                    :: MaximumFlux, MaxErosionDepth
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ModifyConsolidatedErosionFluxes")


        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB

        call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_)  &
            stop 'ModifyConsolidatedErosionCoef - ModuleInterfaceSedimentWater - ERR01'

        if(CohesiveSediment%Evolution%Erosion)then

            if(Me%ExternalVar%Now .ge. CohesiveSediment%Evolution%NextCompute)then

                CHUNK = CHUNK_J(JLB, JUB)
                !$OMP PARALLEL PRIVATE(i,j,KUB,MaxErosionDepth,MaximumFlux)
                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                do j = JLB, JUB
                do i = ILB, IUB

                    if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint .and. &
                       Me%ExtSed%WaterPoints2D (i,j) == WaterPoint)then

                        KUB = Me%ExtSed%KTop(i, j)

                        if(CohesiveSediment%Mass_Available(i,j) .le. CohesiveSediment%Mass_Min)then

                            MaxErosionDepth = Me%ExtSed%DWZ(i,j,KUB) - Me%ExtSed%MinLayerThickness

                            !kgsed/m2s = msed * m3sed / m3 / s * kgsed / m3sed 
                            MaximumFlux = MaxErosionDepth * (1.-Me%ExtSed%Porosity(i, j, KUB))/ &
                                          CohesiveSediment%Evolution%DTInterval               * &
                                          Me%ExtSed%SedimentDryDensity(i, j, KUB)
                            
                            if(MaximumFlux .lt. 0.)then
                                !$OMP CRITICAL (MCEF1_ERR02)
                                write(*,*)'Maximum erosion flux cannot negative.', 'i,j', i,j, 'KUB = ', KUB
                                stop 'ModifyConsolidatedErosionFluxes - ModuleInterfaceSedimentWater - ERR02'
                                !MaximumFlux = 0.
                                !$OMP END CRITICAL (MCEF1_ERR02)
                            end if       

                            if(Me%ExtSed%TopCriticalShear(i, j) < Me%Shear_Stress%Tension(i,j))then

                                !$OMP CRITICAL (MCEF2_FNC01)
                                !kgsed/m2s
                                CohesiveSediment%ErosionFlux(i,j) =                                         &
                                ErosionFlux(CriticalShearErosion = Me%ExtSed%TopCriticalShear(i, j),        &
                                            ShearStress          = Me%Shear_Stress%Tension(i,j),            &
                                            ErosionRate          = Me%ErosionRate%Field(i,j))
                                !$OMP END CRITICAL (MCEF2_FNC01)
                                            
                                if(CohesiveSediment%ErosionFlux(i,j) .gt. MaximumFlux)then

                                    CohesiveSediment%ErosionFlux(i,j) = MaximumFlux

                                end if

                            else

                                CohesiveSediment%ErosionFlux(i,j) = 0.

                            end if

                            Me%Consolidation%Flux(i,j) = Me%Consolidation%Flux(i,j) -                   &
                                                         CohesiveSediment%ErosionFlux(i,j)

                            CohesiveSediment%FluxToWater(i,j) = CohesiveSediment%FluxToWater(i,j)     + &
                                                                CohesiveSediment%ErosionFlux(i,j)

                        endif
                        
                    end if

                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL

            end if

        end if


        nullify(CohesiveSediment)


        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ModifyConsolidatedErosionFluxes")

#endif   
    end subroutine ModifyConsolidatedErosionFluxes

    !--------------------------------------------------------------------------

    subroutine ModifyDepositionFluxes
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j, kbottom
        type(T_Property),       pointer         :: PropertyX
        real, dimension(:,:,:), pointer         :: FreeConvFlux
        real, dimension(:,:,:), pointer         :: WaterPropertyConcentration
        character(len=StringLength)             :: WaterPropertyUnits
        real                                    :: WaterPropertyISCoef
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        !A if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ModifyDepositionFluxes")

        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB
        KUB = Me%WaterWorkSize3D%KUB

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%Deposition)then
                if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then

                    call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                                          ConcentrationX    = WaterPropertyConcentration,       &
                                          PropertyXIDNumber = PropertyX%ID%IDNumber,            &
                                          PropertyXUnits    = WaterPropertyUnits,               &
                                          PropertyXISCoef   = WaterPropertyISCoef,              &
                                          STAT              = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'ModifyDepositionFluxes - ModuleInterfaceSedimentWater - ERR01'
                
                    call Get_FreeConvFlux(FreeVerticalMovementID = Me%ObjFreeVerticalMovement,  &
                                          PropertyID             = PropertyX%ID%IDNumber,       &
                                          FreeConvFlux           = FreeConvFlux,                & 
                                          STAT                   = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'ModifyDepositionFluxes - ModuleInterfaceSedimentWater - ERR02'

                   
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            Me%DepositionProbability(i,j) =                                       &
                               DepositionProbability(Me%Critical_Shear_Deposition%Field(i,j),     &
                                                     Me%Shear_Stress%Tension(i,j))
                        
                            !   kg          kg        1
                            !--------- = ------- * ------- 
                            ! m2 * s        s         m2 
                            PropertyX%DepositionFlux(i,j) = -1. * FreeConvFlux(i,j,kbottom)     * &
                                                            WaterPropertyISCoef                 / &
                                                            Me%ExternalVar%GridCellArea(i,j)

                        end if

                    enddo
                    enddo

                    if(PropertyX%Mass_Limitation)then

                        if (MonitorPerformance) then
                            call StartWatch ("ModuleInterfaceSedimentWater", "ModifyDepositionFluxes")
                        endif

                        CHUNK = CHUNK_J(JLB, JUB)
                        !$OMP PARALLEL PRIVATE(i,j)
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = JLB, JUB
                        do i = ILB, IUB
                    
                            if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                                PropertyX%Mass_Available(i,j) = PropertyX%Mass_Available(i,j)       + &
                                                                PropertyX%DepositionFlux(i,j)       * &
                                                                PropertyX%Evolution%DTInterval

                            end if

                        enddo
                        enddo
                        !$OMP END DO
                        !$OMP END PARALLEL

                        if (MonitorPerformance) then
                            call StopWatch ("ModuleInterfaceSedimentWater", "ModifyDepositionFluxes")
                        endif

                    end if

                    call UnGetWaterProperties(Me%ObjWaterProperties, WaterPropertyConcentration, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'ModifyDepositionFluxes - ModuleInterfaceSedimentWater - ERR03'


                    call UngetFreeVerticalMovement(Me%ObjFreeVerticalMovement, FreeConvFlux, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'ModifyDepositionFluxes - ModuleInterfaceSedimentWater - ERR04'

                end if
            end if

            PropertyX => PropertyX%Next

        enddo

        !A if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ModifyDepositionFluxes")

    end subroutine ModifyDepositionFluxes

    
    !--------------------------------------------------------------------------
    

    real function DepositionProbability (CriticalShearDeposition, ShearStress)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                        :: CriticalShearDeposition
        real, intent(IN)                        :: ShearStress

        !Begin-----------------------------------------------------------------

        if (ShearStress <= CriticalShearDeposition .and. ShearStress >= 0.)then
            
            DepositionProbability = 1. - ShearStress / CriticalShearDeposition
                             

        elseif (ShearStress < 0.) then 
            
            stop 'DepositionProbability - ModuleInterfaceSedimentWater - ERR01'
        
        else
            
            DepositionProbability = 0.
        
        endif

    end function DepositionProbability

    
    !--------------------------------------------------------------------------
    
    
    subroutine ModifyDissolvedFluxes

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j, kbottom
        type(T_Property),       pointer         :: PropertyX
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ModifyDissolvedFluxes")

        IUB = Me%WaterWorkSize3D%IUB
        JUB = Me%WaterWorkSize3D%JUB
        ILB = Me%WaterWorkSize3D%ILB
        JLB = Me%WaterWorkSize3D%JLB
        KUB = Me%WaterWorkSize3D%KUB

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%WaterFluxes .and. (.not. PropertyX%Particulate))then
                
                if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then
                    
                    CHUNK = CHUNK_J(JLB, JUB)
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint) then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            PropertyX%FluxToWater(i,j) = PropertyX%FluxToWater(i,j)     + &
                                                         PropertyX%Mass_Available(i,j)  / &
                                                         PropertyX%Evolution%DTInterval

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                end if
            end if

            PropertyX => PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ModifyDissolvedFluxes")

    end subroutine ModifyDissolvedFluxes


    !--------------------------------------------------------------------------
    
    
    subroutine ComputeWaterFlux

#ifndef _SEDIMENT_

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ComputeWaterFlux")

        IUB = Me%WorkSize2D%IUB
        JUB = Me%WorkSize2D%JUB
        ILB = Me%WorkSize2D%ILB
        JLB = Me%WorkSize2D%JLB
        
        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,KUB)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if(Me%ExtSed%WaterPoints2D(i,j) == WaterPoint)then

                KUB = Me%ExtSed%KTop(i, j)

                !m3water/s  = kgsed.m-2.s-1 * m2 * m3water/m3 /(kgsed.m-3sed * m3sed.m3)
                Me%WaterFlux(i,j) = -1. * Me%Consolidation%Flux(i,j)            * &
                                    Me%ExternalVar%GridCellArea(i,j)            * &
                                    Me%ExtSed%Porosity(i, j, KUB)               / &
                                    (Me%ExtSed%SedimentDryDensity(i, j, KUB)    * &
                                    (1.- Me%ExtSed%Porosity(i, j, KUB)))

            end if

        end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ComputeWaterFlux")

#endif

    end subroutine ComputeWaterFlux


    !--------------------------------------------------------------------------
    
    
    subroutine ComputeConsolidation

#ifndef _SEDIMENT_

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB
        integer                                 :: i, j, STAT_CALL
        type(T_Property), pointer               :: CohesiveSediment
        real                                    :: Initial_Mass_Available
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ComputeConsolidation")
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

        call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_)stop 'ComputeConsolidation - ModuleInterfaceSedimentWater - ERR10'

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if(Me%ExtSed%WaterPoints2D(i,j) == WaterPoint)then
                Me%Consolidation%Flux(i,j) = 0.
            end if
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
              
        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ComputeConsolidation")
                
        !T !$OMP PARALLEL PRIVATE(i,j,Initial_Mass_Available)
        !T !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i,j) .eq. WaterPoint) then

                if(Me%Shear_Stress%Tension(i,j) &
                    < Me%Critical_Shear_Deposition%Field(i,j) .and. &
                   Me%ExtSed%SedimentColumnFull(i,j) == 0) then 

                    !Consolidation is solved implicitly to avoid instability
                    Initial_Mass_Available = CohesiveSediment%Mass_Available(i,j)

                    !kg/m2 = kg/m2 / (s * s-1)
                    CohesiveSediment%Mass_Available(i,j) = &
                        CohesiveSediment%Mass_Available(i,j)        / &
                        (1. + CohesiveSediment%Evolution%DTInterval * &
                        Me%Consolidation%Rate%Field(i,j))

                                                        
                    !kg/m2s = kg/m2 / s  (Positive if it entering consolidated sediment compartment)
                    Me%Consolidation%Flux(i,j) = Me%Consolidation%Flux(i,j)            + &
                                                (Initial_Mass_Available                - &
                                                 CohesiveSediment%Mass_Available(i,j)) / &
                                                 CohesiveSediment%Evolution%DTInterval
                else

                    Me%Consolidation%Flux(i,j) = 0.

                end if

            endif

        end do
        end do
        !T !$OMP END DO               
        !T !$OMP END PARALLEL

#endif

    end subroutine ComputeConsolidation

    
    !--------------------------------------------------------------------------
    
    
    subroutine DissolvedSedimentWaterFluxes

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB
        integer                                 :: i, j, kbottom
        type(T_Property),       pointer         :: PropertyX
        real                                    :: DiffusiveFlux
        real                                    :: AdvectiveFlux
        real                                    :: DT
        real                                    :: BoundaryLayerThickness
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "DissolvedSedimentWaterFluxes")

        IUB = Me%WorkSize2D%IUB
        JUB = Me%WorkSize2D%JUB
        ILB = Me%WorkSize2D%ILB
        JLB = Me%WorkSize2D%JLB

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%SedimentWaterFluxes .and. (.not. PropertyX%Particulate))then

                DT = PropertyX%Evolution%DTInterval

                if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then

                    CHUNK = CHUNK_J(JLB, JUB)
                    !$OMP PARALLEL PRIVATE(i,j,kbottom,BoundaryLayerThickness,AdvectiveFlux,DiffusiveFlux)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if(Me%ExtWater%OpenPoints2D(i,j) == WaterPoint .and. &
                           Me%ExtSed%WaterPoints2D (i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            if(Me%Shear_Stress%Velocity(i, j) > 0.)then
                                BoundaryLayerThickness = 2. * WaterCinematicVisc / Me%Shear_Stress%Velocity(i,j)
                            else
                                BoundaryLayerThickness = Me%ExtWater%DWZ(i,j,kbottom)/2.
                            end if


                            if    (Me%WaterFlux(i,j) > 0.)then

                                AdvectiveFlux = PropertyX%SedimentConcentration(i,j) * Me%WaterFlux(i,j) / &
                                                Me%ExternalVar%GridCellArea(i,j)

                                DiffusiveFlux = 0.

                            elseif(Me%WaterFlux(i,j) < 0.)then

                                AdvectiveFlux = PropertyX%WaterConcentration(i,j)    * Me%WaterFlux(i,j) / &
                                                Me%ExternalVar%GridCellArea(i,j)


                                !kg/m2s = (kg/m3 - kg/m3) * m2/s / m
                                DiffusiveFlux = (PropertyX%SedimentConcentration(i,j)               - &
                                                 PropertyX%WaterConcentration(i,j))                 * &
                                                 PropertyX%MolecularDifCoef%Field(i,j)              / &
                                                 BoundaryLayerThickness



                            elseif(Me%WaterFlux(i,j)== 0.)then

                                AdvectiveFlux = 0.

                                !kg/m2s = (kg/m3 - kg/m3) * m2/s / m
                                DiffusiveFlux = (PropertyX%SedimentConcentration(i,j)               - &
                                                 PropertyX%WaterConcentration(i,j))                 * &
                                                 PropertyX%MolecularDifCoef%Field(i,j)              / &
                                                 BoundaryLayerThickness


                            end if
                            
                            PropertyX%FluxToWater(i,j)    = PropertyX%FluxToWater(i,j)    + DiffusiveFlux + AdvectiveFlux

                            PropertyX%FluxToSediment(i,j) = PropertyX%FluxToSediment(i,j) - DiffusiveFlux - AdvectiveFlux

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                end if
            end if

            PropertyX => PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "DissolvedSedimentWaterFluxes")

    end subroutine DissolvedSedimentWaterFluxes



    !--------------------------------------------------------------------------
    
    
    subroutine ParticulateSedimentWaterFluxes

        !Local-----------------------------------------------------------------
        integer                                 :: IUB, JUB, ILB, JLB, KUB
        integer                                 :: i, j
        type(T_Property),       pointer         :: PropertyX, CohesiveSediment
        real                                    :: DT
        integer                                 :: STAT_CALL
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "ParticulateSedimentWaterFluxes")

        IUB = Me%WorkSize2D%IUB
        JUB = Me%WorkSize2D%JUB
        ILB = Me%WorkSize2D%ILB
        JLB = Me%WorkSize2D%JLB

        call Search_Property(CohesiveSediment, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_)stop 'ComputeConsolidation - ModuleInterfaceSedimentWater - ERR10'



        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%SedimentWaterFluxes .and. PropertyX%Particulate)then

                DT = PropertyX%Evolution%DTInterval

                if(Me%ExternalVar%Now .ge. PropertyX%Evolution%NextCompute)then

                    CHUNK = CHUNK_J(JLB, JUB)
                    !$OMP PARALLEL PRIVATE(i,j,KUB)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = JLB, JUB
                    do i = ILB, IUB
                    
                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            if    (Me%Consolidation%Flux(i,j) > 0.)then !consolidation

                                PropertyX%FluxToSediment(i,j) = PropertyX%FluxToSediment(i,j)        + &
                                                                PropertyX%Mass_Available(i,j)        / &
                                                                CohesiveSediment%Mass_Available(i,j) * &
                                                                Me%Consolidation%Flux(i,j)

                            elseif(Me%Consolidation%Flux(i,j) < 0.)then !erosion

                                KUB = Me%ExtSed%Ktop(i, j)

                                PropertyX%FluxToWater(i,j)    = PropertyX%FluxToWater(i,j)          - &
                                                                Me%Consolidation%Flux(i,j)          * &
                                                                PropertyX%SedimentConcentration(i,j)


                                PropertyX%FluxToSediment(i,j) = PropertyX%FluxToSediment(i,j)       - &
                                                                PropertyX%FluxToWater(i,j)



                            else

                                PropertyX%FluxToSediment(i,j) = 0. 

                            end if

                            PropertyX%Mass_Available(i,j) = PropertyX%Mass_Available(i,j)        - &
                                                            PropertyX%FluxToSediment(i,j)        * &
                                                            PropertyX%Evolution%DTInterval
                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                end if
            end if

            PropertyX => PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "ParticulateSedimentWaterFluxes")

    end subroutine ParticulateSedimentWaterFluxes
    
    
    !--------------------------------------------------------------------------

    
    subroutine SetFluxesToWaterColumn
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property),       pointer         :: PropertyX

        !Begin-----------------------------------------------------------------

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%WaterFluxes)then

                call SetFluxToWaterProperties(WaterPropertiesID = Me%ObjWaterProperties, &
                                              PropertyID        = PropertyX%ID%IDNumber, &
                                              Flux              = PropertyX%FluxToWater, &
                                              STAT              = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                             &
                    stop 'SetFluxesToWaterColumn - ModuleInterfaceSedimentWater - ERR01'

            end if

            PropertyX => PropertyX%Next

        enddo


    end subroutine SetFluxesToWaterColumn

    !--------------------------------------------------------------------------


#ifndef _SEDIMENT_
    subroutine SetFluxesToSedimentColumn
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property),       pointer         :: PropertyX

        !Begin-----------------------------------------------------------------

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%SedimentWaterFluxes)then

                call SetFluxToSedimentProperties(SedimentPropertiesID = Me%ObjSedimentProperties, &
                                                 PropertyID           = PropertyX%ID%IDNumber,    &
                                                 Flux                 = PropertyX%FluxToSediment, &
                                                 STAT                 = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                &
                    stop 'SetFluxesToSedimentColumn - ModuleInterfaceSedimentWater - ERR01'

            end if

            PropertyX => PropertyX%Next

        enddo


    end subroutine SetFluxesToSedimentColumn

    !--------------------------------------------------------------------------

#endif
    
    !--------------------------------------------------------------------------


    subroutine ModifySandTransport

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call UnGetGridData(Me%ObjWaterGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'ModifySandTransport - ModuleInterfaceSedimentWater - ERR01'

        call ModifySand(Me%ObjSand, Me%Shear_Stress%Tension,                     &
                        Me%Rugosity%Field,                                       &
                        Me%WaveShear_Stress%Rugosity%Field,                      &
                        Me%ExtWater%WaterColumn,                                 &
                        Me%Shear_Stress%CurrentU,                                &
                        Me%Shear_Stress%CurrentV,                                &
                        Me%Shear_Stress%CurrentVel,                              &
                        Me%WaveShear_Stress%Tension,                             &
                        Me%WaveShear_Stress%TensionCurrents,                     &
                        Me%Shear_Stress%Velocity,                                &
                        Me%ExtWater%MinWaterColumn,                              &
                        STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_)                                                &
            stop 'ModifySandTransport - ModuleInterfaceSedimentWater - ERR02'

        call GetGridData(Me%ObjWaterGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'ModifySandTransport - ModuleInterfaceSedimentWater - ERR03'


    end subroutine ModifySandTransport


    !--------------------------------------------------------------------------
    
    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX

        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D  = Me%Shear_Stress%Tension,                  &
                            STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleInterfaceSedimentWater - ERR01'
        
        
        PropertyX   => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                call WriteTimeSerie(Me%ObjTimeSerie,                            &
                                    Data2D  = PropertyX%Mass_Available,         &
                                    STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                      &
                    stop 'OutPut_TimeSeries - ModuleInterfaceSedimentWater - ERR02'
            endif
            PropertyX=>PropertyX%Next
        enddo
    
    end subroutine OutPut_TimeSeries


    !--------------------------------------------------------------------------

    subroutine OutPut_BoxTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: CHUNK

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB 
        IUB = Me%Size2D%IUB 
        JLB = Me%Size2D%JLB 
        JUB = Me%Size2D%JUB 

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then

                Me%Scalar2D(:,:) = 0.

                if (MonitorPerformance) then
                    call StartWatch ("ModuleInterfaceSedimentWater", "OutPut_BoxTimeSeries")
                endif

                CHUNK = CHUNK_J(JLB, JUB)
                !$OMP PARALLEL PRIVATE(I,J)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do J = JLB, JUB
                do I = ILB, IUB
                    if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then
                        Me%Scalar2D(i,j) = PropertyX%Mass_Available(i, j) * &
                                           Me%ExternalVar%GridCellArea(i,j)
                    endif
                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL
                
                if (MonitorPerformance) then
                    call StopWatch ("ModuleInterfaceSedimentWater", "OutPut_BoxTimeSeries")
                endif
                
                call BoxDif(Me%ObjBoxDif,                       &
                            Me%Scalar2D,                        &
                            "Bottom "//trim(PropertyX%ID%name), &
                            Me%ExtWater%WaterPoints2D,          &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  &
                    stop 'OutPut_BoxTimeSeries - ModuleInterfaceSedimentWater - ERR01'

                Me%Scalar2D(:,:) = null_real

            endif
            PropertyX=>PropertyX%Next
        enddo

    end subroutine OutPut_BoxTimeSeries

    !--------------------------------------------------------------------------

    subroutine OutPut_Statistics (Value2D, StatisticsID)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer :: Value2D
        integer                       :: StatisticsID

        !Local-----------------------------------------------------------------
        integer                       :: MethodStatistic, Value2DStat2D
        integer                       :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetStatisticMethod (StatisticsID, MethodStatistic, STAT = STAT_CALL)                                     
        if (STAT_CALL /= SUCCESS_)  &
            stop 'OutPut_Statistics - ModuleInterfaceSedimentWater - ERR01'
                                                                                    
        call GetStatisticParameters (StatisticsID,                                   &
                                     Value2DStat2D = Value2DStat2D,                  &
                                     STAT          = STAT_CALL)                        
        if (STAT_CALL /= SUCCESS_)  &
            stop 'OutPut_Statistics - ModuleInterfaceSedimentWater - ERR02'
                                                                                    
                                                                                    
        if (MethodStatistic /= Value2DStat2D)                                            &
            stop 'OutPut_Statistics - ModuleInterfaceSedimentWater - ERR03'
                                                                                    
                                                                                    
        call ModifyStatistic (StatisticsID,                                              &
                              Value2D       = Value2D,                                   &
                              WaterPoints2D = Me%ExtWater%WaterPoints2D,                 &
                              STAT          = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)  &
            stop 'OutPut_Statistics - ModuleInterfaceSedimentWater - ERR04'


    end subroutine OutPut_Statistics


    !--------------------------------------------------------------------------


    subroutine Benthos_Processes
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property   ),       pointer     :: PropertyX
        type (T_BenthicRate),       pointer     :: BenthicRateX
        integer                                 :: WILB, WIUB, WJLB, WJUB
        integer                                 :: i, j, kbottom
        real, dimension(:,:,:),     pointer     :: ConcentrationOld
        character(len=StringLength)             :: WaterPropertyUnits
        real                                    :: WaterPropertyISCoef
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        
        WIUB = Me%WorkSize2D%IUB
        WJUB = Me%WorkSize2D%JUB
        WILB = Me%WorkSize2D%ILB
        WJLB = Me%WorkSize2D%JLB

        CHUNK = CHUNK_J(WJLB, WJUB)
        
        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceSedimentWater", "Benthos_Processes")
        endif
        
        if (Me%ExternalVar%Now .GE. Me%Coupled%Benthos%NextCompute) then

            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))

                if(PropertyX%ID%IDNumber == Temperature_)then

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%WaterConcentration,     &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR01'

                elseif(PropertyX%ID%IDNumber == Oxygen_  )then
                    
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = WJLB, WJUB
                    do i = WILB, WIUB

                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                      Me%ExtWater%VolumeZ(i,j,kbottom)

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               &
                                          Oxygen2D      = PropertyX%WaterConcentration,     &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR04'


                else
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    if(PropertyX%Particulate)then
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%Mass_Available(i,j)     * &
                                                          Me%ExternalVar%GridCellArea(i,j)
                            end if

                        enddo
                        enddo
                        !$OMP END DO

                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)

                            end if

                        enddo
                        enddo
                        !$OMP END DO
                    
                    end if
                    !$OMP END PARALLEL

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR05'

                end if


                PropertyX => PropertyX%Next

            end do
            
            Me%Coupled%Benthos%NextCompute = Me%Coupled%Benthos%NextCompute +               &
                                             Me%Coupled%Benthos%DT_Compute
        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if (PropertyX%Evolution%Benthos) then

                if (Me%ExternalVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    if(PropertyX%Particulate)then
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%Mass_Available(i,j)     * &
                                                          Me%ExternalVar%GridCellArea(i,j)
                            end if

                        enddo
                        enddo
                        !$OMP END DO

                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)

                            end if

                        enddo
                        enddo
                        !$OMP END DO
                    end if
                    !$OMP END PARALLEL


                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               & 
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR06'

                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    if(.not. PropertyX%Particulate)then
                        !$OMP MASTER
                        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,    &
                                              ConcentrationX    = ConcentrationOld,         &
                                              PropertyXIDNumber = PropertyX%ID%IDNumber,    &
                                              PropertyXUnits    = WaterPropertyUnits,       &
                                              PropertyXISCoef   = WaterPropertyISCoef,      &
                                              STAT              = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)                                         &
                            stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR07'
                        !$OMP END MASTER
                        !$OMP BARRIER
                        
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB
                    
                            if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)
                                !kg m-2 s-1                = kg m-2 s-1 + (kg - kg/m3 * m3)/(m2 * s)
                                PropertyX%FluxToWater(i,j) = PropertyX%FluxToWater(i,j)                           +   &    
                                                             (PropertyX%MassInKg(i,j)                             -   &
                                                              ConcentrationOld(i,j,kbottom) * WaterPropertyISCoef *   &
                                                              Me%ExtWater%VolumeZ(i,j,kbottom))                   /   &
                                                             (Me%ExternalVar%GridCellArea(i,j)                    *   &
                                                              PropertyX%Evolution%DTInterval)
 
                            end if

                        enddo
                        enddo
                        !$OMP END DO

                        !$OMP MASTER
                        call UnGetWaterProperties(Me%ObjWaterProperties, ConcentrationOld, STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR08'
                        !$OMP END MASTER
                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%Mass_Available(i,j) =  PropertyX%MassInKg(i,j) / &
                                                                 Me%ExternalVar%GridCellArea(i,j)
                            end if

                        enddo
                        enddo
                        !$OMP END DO

                    end if
                    !$OMP END PARALLEL
                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceSedimentWater", "Benthos_Processes")
        endif
        
        BenthicRateX => Me%FirstBenthicRate

        do while (associated(BenthicRateX))

            call GetRateFlux(InterfaceID    = Me%ObjInterface,                          &
                             FirstProp      = BenthicRateX%FirstProp%IDNumber,          &
                             SecondProp     = BenthicRateX%SecondProp%IDNumber,         &
                             RateFlux2D     = BenthicRateX%Field,                       &
                             WaterPoints2D  = Me%ExtWater%WaterPoints2D,                &
                             STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR09'

            where (Me%ExtWater%WaterPoints2D == WaterPoint)                             &
                BenthicRateX%Field = BenthicRateX%Field * Me%ExternalVar%GridCellArea / &
                                     Me%Coupled%Benthos%DT_Compute


            call BoxDif(Me%ObjBoxDif,                                                   &
                        BenthicRateX%Field,                                             &
                        BenthicRateX%ID%Name,                                           &
                        Me%ExtWater%WaterPoints2D,                                      &
                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR10'

            BenthicRateX => BenthicRateX%Next
        enddo


    end subroutine Benthos_Processes


!--------------------------------------------------------------------------

subroutine BenthicEcology_Processes
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property   ),       pointer     :: PropertyX
        type (T_BenthicRate),       pointer     :: BenthicRateX
        integer                                 :: WILB, WIUB, WJLB, WJUB
        integer                                 :: i, j, kbottom
        real, dimension(:,:,:),     pointer     :: ConcentrationOld
        character(len=StringLength)             :: WaterPropertyUnits
        real                                    :: WaterPropertyISCoef
        real, dimension(:,:,:), pointer         :: WaterPropertyConcentration
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        
        WIUB = Me%WorkSize2D%IUB
        WJUB = Me%WorkSize2D%JUB
        WILB = Me%WorkSize2D%ILB
        WJLB = Me%WorkSize2D%JLB

        CHUNK = CHUNK_J(WJLB, WJUB)
        
        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceSedimentWater", "BenthicEcology_Processes")
        endif
        
       
        if (Me%ExternalVar%Now .GE. Me%Coupled%BenthicEcology%NextCompute) then
        
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = WJLB, WJUB
                    do i = WILB, WIUB

                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            Me%ExtWater%WaterVolume(i,j) =  Me%ExtWater%VolumeZ(i,j,kbottom)
                           

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL
        
        

            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))

                if(PropertyX%ID%IDNumber == Temperature_)  then

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%WaterConcentration,     &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          WaterVolume2D = Me%ExtWater%WaterVolume ,         &
                                          CellArea2D    = Me%ExternalVar%GridCellArea,      &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR01'

        
                elseif(PropertyX%ID%IDNumber == Cohesive_Sediment_  )then
                    
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = WJLB, WJUB
                    do i = WILB, WIUB

                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                      Me%ExtWater%VolumeZ(i,j,kbottom)

                        end if

                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               &
                                          Sediment      = PropertyX%WaterConcentration,     &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR04'


                else
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                    if(PropertyX%Particulate)then
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)
                       !For particulate properties, MassInkg is different from Mass_FromWater 
                       ! MassInKg is the mass of particulate material (In Kg) laying on the bottom;
                       ! Mass_FromWater is the mass of particulate material (In Kg) in the cell closest to the bottom
                                PropertyX%MassInKg(i,j) = PropertyX%Mass_Available(i,j)     * &   
                                                          Me%ExternalVar%GridCellArea(i,j)
                                PropertyX%WaterConcentration(i,j) = PropertyX%Mass_FromWater(i,j)     / &
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)                          
                                                                                        
                            end if

                        enddo
                        enddo
                        !$OMP END DO
                        
                        
                        
                        

                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                !PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &  !For particulate properties, MassInkg and Mass_FromWater are the same
                                !                          Me%ExtWater%VolumeZ(i,j,kbottom)
                                PropertyX%Mass_FromWater(i,j) = PropertyX%WaterConcentration(i,j) * &  !For particulate properties, MassInkg and Mass_FromWater are the same
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)


                            end if

                        enddo
                        enddo
                        !$OMP END DO
                    
                    end if
                    !$OMP END PARALLEL


                            
                           !if(PropertyX%Evolution%BenthicOnly)then
                           
                           !call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                           !                       PropertyID    = PropertyX%ID%IDNumber,            &
                           !                       Concentration = PropertyX%MassInKg,               &
                           !                       WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                           !                       OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                           !                       STAT          = STAT_CALL)
                           ! if (STAT_CALL .NE. SUCCESS_)                                            &
                           !     stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR05'
                           
                           !else
                            
                            call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                                  PropertyID    = PropertyX%ID%IDNumber,            &
                                                  Concentration = PropertyX%MassInKg,               &  ! Concentration is NOT a concentration but a mass
                                              MassInKgFromWater = PropertyX%Mass_FromWater,      &
                                                  WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                                  OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                                  STAT          = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                            &
                                stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR05'
                          
                          
                         ! endif
                           
                
                end if


                PropertyX => PropertyX%Next

            end do
            
            Me%Coupled%BenthicEcology%NextCompute = Me%Coupled%BenthicEcology%NextCompute +               &
                                             Me%Coupled%BenthicEcology%DT_Compute
        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

       if01:   if (PropertyX%Evolution%BenthicEcology) then

        if02:   if (Me%ExternalVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                   if1: if(PropertyX%Particulate)then
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%Mass_Available(i,j)     * &
                                                          Me%ExternalVar%GridCellArea(i,j)
                                                          
                                PropertyX%WaterConcentration(i,j) = PropertyX%Mass_FromWater(i,j)     / &
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)                        
                                                          
                            end if

                        enddo
                        enddo
                        !$OMP END DO

                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%Mass_FromWater(i,j) = PropertyX%WaterConcentration(i,j) * &  ! If the property is dissolved, the mass is retrieved from the water
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)
                                                          
                                PropertyX%MassInKg(i,j)=0. ! since the property is dissolved, it has no mass on the bottom 

                            end if

                        enddo
                        enddo
                        !$OMP END DO
                    end if if1
                    !$OMP END PARALLEL


                   !if(PropertyX%Evolution%BenthicOnly)then
                           
                           !call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                           !                       PropertyID    = PropertyX%ID%IDNumber,            &
                           !!                       Concentration = PropertyX%MassInKg,               & 
                           !!                       WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                           !                       OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                           !                       DTProp        = PropertyX%Evolution%DTInterval,   &
                           !                       STAT          = STAT_CALL)
                           ! if (STAT_CALL .NE. SUCCESS_)                                            &
                           !     stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR06'
                           
                     !else

                            call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                                  PropertyID    = PropertyX%ID%IDNumber,            &
                                                  Concentration = PropertyX%MassInKg,               & 
                                              MassInKgFromWater = PropertyX%Mass_FromWater,      &
                                                  WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                                  OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                                  DTProp        = PropertyX%Evolution%DTInterval,   &
                                                  STAT          = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_)                                            &
                                stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR06'
                     
                     
                     
                     !endif
                    !$OMP PARALLEL PRIVATE(i,j,kbottom)
                if2: if(.not. PropertyX%Particulate)then
                        !$OMP MASTER
                        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,    &
                                              ConcentrationX    = ConcentrationOld,         &
                                              PropertyXIDNumber = PropertyX%ID%IDNumber,    &
                                              PropertyXUnits    = WaterPropertyUnits,       &
                                              PropertyXISCoef   = WaterPropertyISCoef,      &
                                              STAT              = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)                                         &
                            stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR07'
                        !$OMP END MASTER
                        !$OMP BARRIER
                        
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB
                    
                            if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)
                                !kg m-2 s-1                = kg m-2 s-1 + (kg - kg/m3 * m3)/(m2 * s)
                               
                               
                                PropertyX%FluxToWater(i,j) = PropertyX%FluxToWater(i,j)                           +   &    
                                                            ! (PropertyX%MassInKg(i,j)                             -   &
                                                            (PropertyX%Mass_FromWater(i,j)                             -   &
                                                              ConcentrationOld(i,j,kbottom) * WaterPropertyISCoef *   &
                                                              Me%ExtWater%VolumeZ(i,j,kbottom))                   /   &
                                                             (Me%ExternalVar%GridCellArea(i,j)                    *   &
                                                              PropertyX%Evolution%DTInterval)
 
                            end if

                        enddo
                        enddo
                        !$OMP END DO

                        !$OMP MASTER
                        call UnGetWaterProperties(Me%ObjWaterProperties, ConcentrationOld, STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR08'
                        !$OMP END MASTER
                    else
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%Mass_Available(i,j) =  PropertyX%MassInKg(i,j) / &
                                                                 Me%ExternalVar%GridCellArea(i,j)
                                PropertyX%WaterConcentration(i,j) = PropertyX%Mass_FromWater(i,j)  / &
                                                                 Me%ExtWater%VolumeZ(i,j,kbottom)
                                
                                if(.not.PropertyX%Evolution%BenthicOnly)then
                                
                                        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,    &
                                                      ConcentrationX    = ConcentrationOld,         &
                                                      PropertyXIDNumber = PropertyX%ID%IDNumber,    &
                                                      PropertyXUnits    = WaterPropertyUnits,       &
                                                      PropertyXISCoef   = WaterPropertyISCoef,      &
                                                      STAT              = STAT_CALL)
                                if(STAT_CALL .ne. SUCCESS_)                                         &
                                    stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR03'
                                
                                                                 
                                    PropertyX%FluxToWater(i,j) = PropertyX%FluxToWater(i,j)                           +   &    
                                                                 (PropertyX%Mass_FromWater(i,j)                    -   &
                                                                  ConcentrationOld(i,j,kbottom) * WaterPropertyISCoef *   &
                                                                  Me%ExtWater%VolumeZ(i,j,kbottom))                   /   &
                                                                 (Me%ExternalVar%GridCellArea(i,j)                    *   &
                                                                  PropertyX%Evolution%DTInterval)                                 
                                   
                                    
                                    call UnGetWaterProperties(Me%ObjWaterProperties, ConcentrationOld, STAT = STAT_CALL)
                                    if(STAT_CALL .ne. SUCCESS_)&
                                     stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR08'
                                    
                        
                        
                                   
                                   endif                                  
                            end if
                            
                            

                        enddo
                        enddo
                        !$OMP END DO
                        
                       
                        

                    end if if2
                    !$OMP END PARALLEL
                end if if02

            end if if01

            PropertyX => PropertyX%Next
            
        end do

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceSedimentWater", "BenthicEcology_Processes")
        endif
        
        !BenthicRateX => Me%FirstBenthicRate

       ! do while (associated(BenthicRateX))

        !    call GetRateFlux(InterfaceID    = Me%ObjInterface,                          &
        !                     FirstProp      = BenthicRateX%FirstProp%IDNumber,          &
        !                     SecondProp     = BenthicRateX%SecondProp%IDNumber,         &
        !!                     RateFlux2D     = BenthicRateX%Field,                       &
        !                    WaterPoints2D  = Me%ExtWater%WaterPoints2D,                &
        !                     STAT           = STAT_CALL)
        !    if (STAT_CALL .NE. SUCCESS_)                                                &
        !        stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR09'

        !    where (Me%ExtWater%WaterPoints2D == WaterPoint)                             &
        !        BenthicRateX%Field = BenthicRateX%Field * Me%ExternalVar%GridCellArea / &
        !                             Me%Coupled%BenthicEcology%DT_Compute


        !    call BoxDif(Me%ObjBoxDif,                                                   &
        !                BenthicRateX%Field,                                             &
        !                BenthicRateX%ID%Name,                                           &
        !                Me%ExtWater%WaterPoints2D,                                      &
        !                STAT = STAT_CALL)
         !   if (STAT_CALL .NE. SUCCESS_)                                                &
         !       stop 'BenthicEcology_Processes - ModuleInterfaceSedimentWater - ERR10'

         !   BenthicRateX => BenthicRateX%Next
        !enddo


    end subroutine BenthicEcology_Processes
    !--------------------------------------------------------------------------
    

    subroutine Detritus_Processes
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property),       pointer        :: Detritus
        type (T_Property),       pointer        :: PropertyX
        integer                                 :: WILB, WIUB, WJLB, WJUB, i, j
        integer                                 :: CHUNK

        !Begin-----------------------------------------------------------------
        
        WIUB = Me%WorkSize2D%IUB
        WJUB = Me%WorkSize2D%JUB
        WILB = Me%WorkSize2D%ILB
        WJLB = Me%WorkSize2D%JLB

        call Search_Property(Detritus, PropertyXID = Detritus_, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_) stop 'Detritus_Processes - ModuleInterfaceSedimentWater - ERR01' 
            
        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if(PropertyX%Evolution%Detritus)then

                if (MonitorPerformance) then
                    call StartWatch ("ModuleInterfaceSedimentWater", "Detritus_Processes")
                endif

                CHUNK = CHUNK_J(WJLB, WJUB)
                !$OMP PARALLEL PRIVATE(i,j)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = WJLB, WJUB
                do i = WILB, WIUB

                    if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                        !DepositionFlux is positive if deposition occurs
                        Detritus%Mass_Available(i,j) = Detritus%Mass_Available(i,j)   + &
                                                       PropertyX%DepositionFlux(i,j)  * &
                                                       PropertyX%Evolution%DTInterval

                        PropertyX%Mass_Available(i,j) = PropertyX%Mass_Available(i,j) - &
                                                        PropertyX%DepositionFlux(i,j) * &
                                                        PropertyX%Evolution%DTInterval
                    end if

                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL

                if (MonitorPerformance) then
                    call StopWatch ("ModuleInterfaceSedimentWater", "Detritus_Processes")
                endif

            end if

            PropertyX => PropertyX%Next

        end do
                    

    end subroutine Detritus_Processes

    
    !--------------------------------------------------------------------------


    subroutine CEQUALW2_Processes
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property   ),       pointer     :: PropertyX
        type (T_BenthicRate),       pointer     :: BenthicRateX
        integer                                 :: WILB, WIUB, WJLB, WJUB
        integer                                 :: i, j, kbottom
        real, dimension(:,:,:),     pointer     :: ConcentrationOld
        character(len=StringLength)             :: WaterPropertyUnits
        real                                    :: WaterPropertyISCoef

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "CEQUALW2_Processes")

        WIUB = Me%WorkSize2D%IUB
        WJUB = Me%WorkSize2D%JUB
        WILB = Me%WorkSize2D%ILB
        WJLB = Me%WorkSize2D%JLB

        if (Me%ExternalVar%Now .GE. Me%Coupled%CEQUALW2%NextCompute) then
            
            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))

                if(PropertyX%ID%IDNumber == Temperature_)then

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%WaterConcentration,     &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR01'
                        
                elseif( PropertyX%ID%IDNumber == Oxygen_  )then
                    
                    do j = WJLB, WJUB
                    do i = WILB, WIUB

                        if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                            kbottom = Me%ExtWater%KFloor_Z(i, j)

                            PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                      Me%ExtWater%VolumeZ(i,j,kbottom)

                        end if

                    enddo
                    enddo

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               &
                                          Oxygen2D      = PropertyX%WaterConcentration,     &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Benthos_Processes - ModuleInterfaceSedimentWater - ERR04'

                else

                    if(PropertyX%Particulate)then
                
                        do i = WILB, WIUB
                        do j = WJLB, WJUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%Mass_Available(i,j)     * &
                                                          Me%ExternalVar%GridCellArea(i,j)
                            end if

                        enddo
                        enddo

                    else

                        do i = WILB, WIUB
                        do j = WJLB, WJUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)

                            end if

                        enddo
                        enddo
                    
                    end if

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR01'

                end if

                PropertyX => PropertyX%Next

            end do
            
            Me%Coupled%CEQUALW2%NextCompute = Me%Coupled%CEQUALW2%NextCompute +         &
                                              Me%Coupled%CEQUALW2%DT_Compute

        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if (PropertyX%Evolution%CEQUALW2) then

                if (Me%ExternalVar%Now .GE. PropertyX%Evolution%NextCompute) then
                    
                    if(PropertyX%Particulate)then
                        
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%Mass_Available(i,j)     * &
                                                          Me%ExternalVar%GridCellArea(i,j)
                            end if

                        enddo
                        enddo

                    else
                        
                        do j = WJLB, WJUB
                        do i = WILB, WIUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%MassInKg(i,j) = PropertyX%WaterConcentration(i,j) * &
                                                          Me%ExtWater%VolumeZ(i,j,kbottom)

                            end if

                        enddo
                        enddo
                    
                    end if


                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%MassInKg,               &
                                          WaterPoints2D = Me%ExtWater%WaterPoints2D,        &
                                          OpenPoints2D  = Me%ExtWater%OpenPoints2D,         &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR02'

                    if(.not. PropertyX%Particulate)then
                        
                        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,    &
                                              ConcentrationX    = ConcentrationOld,         &
                                              PropertyXIDNumber = PropertyX%ID%IDNumber,    &
                                              PropertyXUnits    = WaterPropertyUnits,       &
                                              PropertyXISCoef   = WaterPropertyISCoef,      &
                                              STAT              = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)                                         &
                            stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR03'

                        do j = WJLB, WJUB
                        do i = WILB, WIUB
                    
                            if(Me%ExtWater%OpenPoints2D(i,j) == OpenPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)
                                !kg m-2 s-1                = kg m-2 s-1 + (kg - kg/m3 * m3)/(m2 * s)
                                PropertyX%FluxToWater(i,j) = PropertyX%FluxToWater(i,j)                           +   &    
                                                             (PropertyX%MassInKg(i,j)                             -   &
                                                              ConcentrationOld(i,j,kbottom) * WaterPropertyISCoef *   &
                                                              Me%ExtWater%VolumeZ(i,j,kbottom))                   /   &
                                                             (Me%ExternalVar%GridCellArea(i,j)                    *   &
                                                              PropertyX%Evolution%DTInterval)

                            end if

                        enddo
                        enddo

                        call UnGetWaterProperties(Me%ObjWaterProperties, ConcentrationOld, STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR04'
                    else

                        do i = WILB, WIUB
                        do j = WJLB, WJUB

                            if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                                kbottom = Me%ExtWater%KFloor_Z(i, j)

                                PropertyX%Mass_Available(i,j) =  PropertyX%MassInKg(i,j) / &
                                                                 Me%ExternalVar%GridCellArea(i,j)
                            end if

                        enddo
                        enddo

                    end if

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        BenthicRateX => Me%FirstBenthicRate

        do while (associated(BenthicRateX))

            call GetRateFlux(InterfaceID    = Me%ObjInterface,                          &
                             FirstProp      = BenthicRateX%FirstProp%IDNumber,          &
                             SecondProp     = BenthicRateX%SecondProp%IDNumber,         &
                             RateFlux2D     = BenthicRateX%Field,                       &
                             WaterPoints2D  = Me%ExtWater%WaterPoints2D,                &
                             STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR05'

            where (Me%ExtWater%WaterPoints2D == WaterPoint)           &
                BenthicRateX%Field = BenthicRateX%Field             * &
                                     Me%ExternalVar%GridCellArea    / &
                                     Me%Coupled%CEQUALW2%DT_Compute


            call BoxDif(Me%ObjBoxDif,                                                   &
                        BenthicRateX%Field,                                             &
                        BenthicRateX%ID%Name,                                           &
                        Me%ExtWater%WaterPoints2D,                                      &
                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'CEQUALW2_Processes - ModuleInterfaceSedimentWater - ERR06'

            BenthicRateX => BenthicRateX%Next

        enddo
        
        
        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "CEQUALW2_Processes")


    end subroutine CEQUALW2_Processes

    !--------------------------------------------------------------------------

    subroutine OutPut_Results_HDF
        
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        real                                :: Year, Month, Day, Hour, Minute, Second
         
        !Local-----------------------------------------------------------------
        type (T_Property), pointer          :: PropertyX
        integer                             :: OutPutNumber
        integer, dimension(6)               :: TimeAux
        real,    dimension(6), target       :: AuxTime
        real,    dimension(:), pointer      :: TimePtr
        integer                             :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                             :: WorkKLB, WorkKUB
        !----------------------------------------------------------------------


        if (MonitorPerformance) call StartWatch ("ModuleInterfaceSedimentWater", "OutPut_Results_HDF")

        WorkILB = Me%WaterWorkSize3D%ILB 
        WorkIUB = Me%WaterWorkSize3D%IUB 
        WorkJLB = Me%WaterWorkSize3D%JLB 
        WorkJUB = Me%WaterWorkSize3D%JUB
        WorkKLB = Me%WaterWorkSize3D%KLB
        WorkKUB = Me%WaterWorkSize3D%KUB

        OutPutNumber = Me%OutPut%NextOutPut


TOut:   if (Me%ExternalVar%Now >= Me%OutPut%OutTime(OutPutNumber)) then
            
            call ExtractDate(Me%ExternalVar%Now,                         &
                             Year = Year, Month  = Month,  Day    = Day, &
                             Hour = Hour, Minute = Minute, Second = Second)

            TimeAux(1) = int(Year  )
            TimeAux(2) = int(Month )
            TimeAux(3) = int(Day   )
            TimeAux(4) = int(Hour  )
            TimeAux(5) = int(Minute)
            TimeAux(6) = int(Second)

            !Writes current time
            call ExtractDate   (Me%ExternalVar%Now,                  &
                                AuxTime(1), AuxTime(2), AuxTime(3),  &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR01'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                                 Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR02'

            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                 &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR03'

            call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints", &
                                 "-", Array3D = Me%ExtWater%OpenPoints3D,      &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR04'

            call HDF5WriteData  (Me%ObjHDF5, "/Results/ShearStress", "ShearStress",  &
                                 "N/m2", Array2D = Me%Shear_Stress%Tension,          &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR05'

            if(Me%Consolidation%Yes)then

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Consolidation",          &
                                     "consolidation flux",                          &
                                     "kg/m2s", Array2D = Me%Consolidation%Flux,     &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR50'

            end if

                        

            PropertyX => Me%FirstProperty

PropX:      do while (associated(PropertyX))

                if(PropertyX%OutputHDF)then

                    if(PropertyX%Mass_Limitation)then

                        call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyX%ID%Name,                &
                                             PropertyX%ID%Name, PropertyX%ID%Units,                     &
                                             Array2D = PropertyX%Mass_Available,                        &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR06'
                    end if

                    if(PropertyX%Evolution%WaterFluxes)then

                        call HDF5WriteData  (Me%ObjHDF5, "/Results/FluxToWater/"//PropertyX%ID%Name,    &
                                             PropertyX%ID%Name, 'kg/m2s',                               &
                                             Array2D = PropertyX%FluxToWater,                           &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR07'
                    end if
                    
                    if(PropertyX%Evolution%SedimentWaterFluxes)then
                       
                        call HDF5WriteData  (Me%ObjHDF5, "/Results/FluxToSediment/"//PropertyX%ID%Name, &
                                             PropertyX%ID%Name, 'kg/m2s',                               &
                                             Array2D = PropertyX%FluxToSediment,                        &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR08'

                    end if
                    
                    if(PropertyX%Evolution%Deposition)then
                       
                        call HDF5WriteData  (Me%ObjHDF5, "/Results/Deposition/"//PropertyX%ID%Name,     &
                                             PropertyX%ID%Name, 'kg/m2s',                               &
                                             Array2D = PropertyX%DepositionFlux,                        &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR09'

                    end if

                end if

                PropertyX => PropertyX%Next

            enddo PropX

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleInterfaceSedimentWater - ERR10'

            Me%OutPut%NextOutPut = OutPutNumber + 1


        endif  TOut    

        nullify(PropertyX)

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceSedimentWater", "OutPut_Results_HDF")

    end subroutine OutPut_Results_HDF

    !--------------------------------------------------------------------------

    subroutine OutputRestartFile
        
        !Local-----------------------------------------------------------------
        real                                :: Year, Month, Day, Hour, Minute, Second

        !----------------------------------------------------------------------


        if(Me%ExternalVar%Now >= Me%OutPut%RestartOutTime(Me%OutPut%NextRestartOutput))then


            call ExtractDate(Me%ExternalVar%Now,                         &
                             Year = Year, Month  = Month,  Day    = Day, &
                             Hour = Hour, Minute = Minute, Second = Second)

            call Write_Final_HDF

            Me%OutPut%NextRestartOutput = Me%OutPut%NextRestartOutput + 1

            call SetError(WARNING_, INTERNAL_, "Interface Sed.Water restart file saved : ", &
                                  Year, Month, Day, Hour, Minute, Second)

        end if


    end subroutine OutputRestartFile


    !--------------------------------------------------------------------------
    ! This subroutine is responsable for defining       
    ! the next time to actualize the value of each      
    ! property                                          
    subroutine Actualize_Time_Evolution

        !Local-----------------------------------------------------------------
        type (T_Property), pointer :: Property

        !----------------------------------------------------------------------

        Property => Me%FirstProperty  

        do while (associated(Property))

            if (Property%Evolution%Variable) then
                if (Me%ExternalVar%Now.GE.Property%Evolution%NextCompute) then
                        Property%Evolution%LastCompute = Property%Evolution%NextCompute
                        Property%Evolution%NextCompute = Property%Evolution%NextCompute &
                                                       + Property%Evolution%DTInterval
                end if
            end if


            Property => Property%Next
        end do   

        nullify(Property)

    end subroutine Actualize_Time_Evolution


    !--------------------------------------------------------------------------

    
    subroutine SetSubModulesModifier 

        !Local-----------------------------------------------------------------
        logical                                         :: WavesStressON
        integer                                         :: STAT_CALL 
        type(T_Property), pointer                       :: SPM

        !Begin-----------------------------------------------------------------

        
        call SetFluxesToWaterColumn

        if(Me%MacroAlgae)then

            call Search_Property(SPM, PropertyXID = Cohesive_Sediment_, STAT = STAT_CALL)  
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR000'
        
            call SetMacroAlgaeParameters(WaterPropertiesID = Me%ObjWaterProperties,     &
                                         ShearStress       = Me%Shear_Stress%Tension,   &
                                         SPMFlux           = SPM%DepositionFlux,        &
                                         STAT              = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR001' 
            
            nullify(SPM)

        endif

#ifndef _SEDIMENT_ 
        
        call SetFluxesToSedimentColumn

        if(Me%Consolidation%Yes)then

            call SetConsolidationFlux(Me%ObjConsolidation, Me%Consolidation%Flux,               &
                                      STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                        &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR01'

            call SetBottomWaterFlux(Me%ObjHydrodynamic,                                         &
                                    Me%WaterFlux,                                               &
                                    STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                         &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR02'

            call SetSedimentWaterFlux(Me%ObjSedimentProperties,                                 &
                                      Me%WaterFlux,                                             &
                                      STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                         &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR03'

        end if

#endif 
        
        if(Me%ObjTurbGOTM /= 0)then
            call SetTurbGOTMBottomShearVelocity(TurbGOTMID          = Me%ObjTurbGOTM,           &
                                                BottomShearVelocity = Me%Shear_Stress%Velocity, &
                                                STAT                = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR04'
        
            call SetTurbGOTMBottomRugosity(TurbGOTMID     = Me%ObjTurbGOTM,                     &
                                           BottomRugosity = Me%Rugosity%Field,                  &
                                           STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR05'
        end if

#ifndef _LAGRANGIAN_                                         
        if(Me%ObjLagrangian /= 0)then       
#ifdef  _LAGRANGIAN_GLOBAL_                                         
            call SetLagrangianShearGlobal                                               &
                                   (LagrangianID   = Me%ObjLagrangian,                  &
                                    ModelName      = Me%ModelName,                      &
                                    ShearStress    = Me%Shear_Stress%Tension,           &
                                    ShearVelocity  = Me%Shear_Stress%Velocity,          &
                                    STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR06a'
#else
            call SetLagrangianShear(LagrangianID   = Me%ObjLagrangian,                  &
                                    ShearStress    = Me%Shear_Stress%Tension,           &
                                    ShearVelocity  = Me%Shear_Stress%Velocity,          &
                                    STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR06b'
#endif    
        end if
#endif

        if(Me%ObjFreeVerticalMovement /= 0)then

            call SetDepositionProbability(Me%ObjFreeVerticalMovement,                   &
                                          Me%DepositionProbability,                     &
                                          STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)                                                 &
                stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR07'

        end if

        if (Me%ObjHydrodynamic /= 0) then
            if (Me%WaveShear_Stress%Yes) then

                call GetWavesStressON(Me%ObjHydrodynamic,                                &
                                     WavesStressON, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)                                              &
                    stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR08'
                
                if (WavesStressON) then

                    call SetWaveChezyVel(Me%ObjHydrodynamic,                             &
                                         Me%WaveShear_Stress%ChezyVel,                   &
                                         STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)                                          &
                        stop 'SetSubModulesModifier - ModuleInterfaceSedimentWater - ERR09'
                endif
            endif
        endif
            

    end subroutine SetSubModulesModifier

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillInterfaceSedimentWater(ObjInterfaceSedimentWaterID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjInterfaceSedimentWaterID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, STAT_CALL, nUsers              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_           
        type (T_Property),    pointer       :: PropertyX
        type (T_BenthicRate), pointer       :: BenthicRateX

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterfaceSedimentWaterID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mINTERFACESEDIMENTWATER_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%OutPut%Yes) then

                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR01'
                endif

                if (Me%Shear_Stress%Statistics%ON) then
                    call KillStatistic (Me%Shear_Stress%Statistics%ID, STAT = STAT_CALL)
                    stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR02'
                endif

                if (Me%RunsSandTransport) then

                    call KillSand(Me%ObjSand, STAT = STAT_CALL)

                    if(STAT_CALL /= SUCCESS_)                                            &
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR03.'

                endif



                call ReadLockExternalGlobal

                call ReadLockExternalWater
#ifndef _SEDIMENT_
                if(Me%RunsSediments) call ReadLockExternalSediment
#endif
                if (Me%OutPut%WriteFinalFile) call Write_Final_HDF( Final = .true.)

#ifndef _SEDIMENT_
                if(Me%RunsSediments) call ReadUnlockExternalSediment
#endif
                call ReadUnlockExternalWater

                call ReadUnlockExternalGlobal


                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR05'
                
                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR06'
                
                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjWaterGridData)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR07'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjWaterHorizontalMap)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR08'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjWaterGeometry)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR09'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjWaterMap)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR10'

                nUsers = DeassociateInstance(mTURBULENCE_,      Me%ObjTurbulence)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR11'

                nUsers = DeassociateInstance(mHYDRODYNAMIC_,    Me%ObjHydrodynamic)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR12'

                nUsers = DeassociateInstance(mWATERPROPERTIES_, Me%ObjWaterProperties)
                if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR13'

#ifndef _LAGRANGIAN_
                if(Me%ObjLagrangian /= 0)then
                    nUsers = DeassociateInstance(mLAGRANGIAN_, Me%ObjLagrangian)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR14'
                end if
#endif

                if(Me%ObjTurbGOTM /= 0)then
                    nUsers = DeassociateInstance(mTURBGOTM_, Me%ObjTurbGOTM)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR15'
                end if


                if(Me%ObjFreeVerticalMovement /= 0)then
                    nUsers = DeassociateInstance(mFREEVERTICALMOVEMENT_, Me%ObjFreeVerticalMovement)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR16'
                end if

                if(Me%RunsSediments)then
                                    
                    nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjSedimentGridData)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR17'

                    nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjSedimentHorizontalMap)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR18'

                    nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjSedimentGeometry)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR19'

                    nUsers = DeassociateInstance(mMAP_,             Me%ObjSedimentMap)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR20'
#ifndef _SEDIMENT_
                    nUsers = DeassociateInstance(mSEDIMENTPROPERTIES_, Me%ObjSedimentProperties)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR21'

                    nUsers = DeassociateInstance(mCONSOLIDATION_,   Me%ObjConsolidation)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR22'
#endif
                end if

                if (Me%Coupled%BoxTimeSerie%Yes) then
                    
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR23'


                    deallocate(Me%Scalar2D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR24'
                    nullify(Me%Scalar2D)

                end if

                !Kills the TimeSerie
                if (Me%Coupled%TimeSerie%Yes) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR25'
                endif

                deallocate(Me%Shear_Stress%Tension,    STAT = STAT_CALL) 
                if(STAT_CALL .ne. SUCCESS_)&
                    stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR26'
                nullify(Me%Shear_Stress%Tension)

                deallocate(Me%Shear_Stress%Velocity,   STAT = STAT_CALL) 
                if(STAT_CALL .ne. SUCCESS_)&
                    stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR27'
                nullify(Me%Shear_Stress%Velocity)

                if(Me%UseSOD)then
                    deallocate(Me%SOD%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR27a'
                    nullify(Me%SOD%Field)
                endif


                if (Me%RunsSandTransport) then

                    deallocate(Me%Shear_Stress%CurrentVel,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR27a'
                    nullify(Me%Shear_Stress%CurrentVel)

                    deallocate(Me%Shear_Stress%CurrentU,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR27b'
                    nullify(Me%Shear_Stress%CurrentU)

                    deallocate(Me%Shear_Stress%CurrentV,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR27c'
                    nullify(Me%Shear_Stress%CurrentV)


                endif


                if (Me%WaveShear_Stress%Yes) then

                    !Shear stress 
                    deallocate(Me%WaveShear_Stress%Tension        ) 

                    deallocate(Me%WaveShear_Stress%ChezyVel       ) 

                    deallocate(Me%WaveShear_Stress%Rugosity%Field ) 

                    if (Me%RunsSandTransport) then

                        deallocate(Me%WaveShear_Stress%TensionCurrents) 

                    endif

                endif
                
                if (.not. Me%Chezy) then

                    if(Me%Manning)then
                        deallocate(Me%ManningCoef%Field,    STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR28'
                        nullify(Me%ManningCoef%Field)
                    else
                        deallocate(Me%Rugosity%Field,       STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR29'
                        nullify(Me%Rugosity%Field)
                    endif
                
                endif


#ifndef _WAVES_

                if (Me%ObjWaves /= 0) then

                    nUsers = DeassociateInstance (mWAVES_,Me%ObjWaves)
                    if (nUsers == 0) stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR30'
                end if
#endif

                if(Me%Consolidation%Yes)then
                    
                    deallocate(Me%Consolidation%Flux,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR310'
                    nullify(Me%Consolidation%Flux)

                    deallocate(Me%Consolidation%Rate%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR311'
                    nullify(Me%Consolidation%Rate%Field)


                endif   

                if(Me%Coupled%Erosion%Yes)then

                    deallocate(Me%ErosionRate%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR31'
                    nullify(Me%ErosionRate%Field)

                    deallocate(Me%Critical_Shear_Deposition%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR32'
                    nullify(Me%Critical_Shear_Deposition%Field)

                    deallocate(Me%Critical_Shear_Erosion%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR33'
                    nullify(Me%Critical_Shear_Erosion%Field)

                end if
                
                if(Me%Coupled%Deposition%Yes)then
                    deallocate(Me%DepositionProbability, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ER41'
                    nullify(Me%DepositionProbability)
                end if

                if (Me%Coupled%Benthos%Yes) then
                    call KillInterface (Me%ObjInterface, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ER41a'
                endif


                if (Me%Coupled%CEQUALW2%Yes) then
                    call KillInterface (Me%ObjInterface, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ER41b'
                endif
                
                
                if (Me%Coupled%BenthicEcology%Yes) then
                    call KillInterface (Me%ObjInterface, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ER41c'
                
                 deallocate(Me%ExtWater%WaterVolume, STAT = STAT_CALL) 
                     if(STAT_CALL .ne. SUCCESS_)&
                     stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR41d'
                     
                deallocate(Me%ExtWater%Sediment, STAT = STAT_CALL) 
                     if(STAT_CALL .ne. SUCCESS_)&
                     stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR41d'
                
                endif

                PropertyX => Me%FirstProperty

                do while(associated(PropertyX))
                    
                    deallocate(PropertyX%Mass_Available,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR36'
                    nullify(PropertyX%Mass_Available)

                    deallocate(PropertyX%MassInKg,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR36a'
                    nullify(PropertyX%MassInKg)
                    
                    deallocate(PropertyX%Mass_FromWater,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR36b'
                    nullify(PropertyX%Mass_FromWater)

                    if(PropertyX%Evolution%Deposition)then

                        deallocate(PropertyX%DepositionFlux,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR37'
                        nullify(PropertyX%DepositionFlux)

                    end if

                    if(PropertyX%Evolution%Erosion)then

                        deallocate(PropertyX%ErosionCoefficient,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR38'
                        nullify(PropertyX%ErosionCoefficient)

                        deallocate(PropertyX%ErosionFlux,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR39'
                        nullify(PropertyX%ErosionFlux)

                    end if


                    if((.not. PropertyX%Particulate) .or. PropertyX%Evolution%WaterFluxes)then

                        deallocate(PropertyX%WaterConcentration,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR40'
                        nullify(PropertyX%WaterConcentration)
                                                
                    end if


                    if(PropertyX%Evolution%WaterFluxes .or. PropertyX%Evolution%SedimentWaterFluxes)then
                        
                        deallocate(PropertyX%FluxToWater,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR41'
                        nullify(PropertyX%FluxToWater)
                        
                    end if

                    if(PropertyX%Evolution%SedimentFluxes .or. PropertyX%Evolution%SedimentWaterFluxes)then
                        
                        deallocate(PropertyX%FluxToSediment,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR41a'
                        nullify(PropertyX%FluxToSediment)

                        deallocate(PropertyX%SedimentConcentration,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR41b'
                        nullify(PropertyX%SedimentConcentration)
                        
                    end if


                    PropertyX => PropertyX%Next

                end do
                
                BenthicRateX => Me%FirstBenthicRate

                do while(associated(BenthicRateX))

                    deallocate(BenthicRateX%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceSedimentWater - ModuleInterfaceSedimentWater - ERR42'
                    nullify(BenthicRateX%Field)

                    BenthicRateX => BenthicRateX%Next

                end do
                
                !Deallocates Instance
                call DeallocateInstance

                ObjInterfaceSedimentWaterID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

    end subroutine KillInterfaceSedimentWater
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_InterfaceSedimentWater), pointer          :: AuxObjInterfaceSedimentWater
        type (T_InterfaceSedimentWater), pointer          :: PreviousObjWaterSedInterface

        !Updates pointers
        if (Me%InstanceID == FirstObjInterfaceSedimentWater%InstanceID) then
            FirstObjInterfaceSedimentWater      => FirstObjInterfaceSedimentWater%Next
        else
            PreviousObjWaterSedInterface        => FirstObjInterfaceSedimentWater
            AuxObjInterfaceSedimentWater        => FirstObjInterfaceSedimentWater%Next
            do while (AuxObjInterfaceSedimentWater%InstanceID /= Me%InstanceID)
                PreviousObjWaterSedInterface    => AuxObjInterfaceSedimentWater
                AuxObjInterfaceSedimentWater    => AuxObjInterfaceSedimentWater%Next
            enddo

            !Now update linked list
            PreviousObjWaterSedInterface%Next   => AuxObjInterfaceSedimentWater%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine Write_Final_HDF( Final )

        !Arguments ---------------------------------------------------------
        logical, optional, intent(IN)           :: Final

        !Local--------------------------------------------------------------
        type(T_Property),           pointer     :: Property
        integer                                 :: ObjHDF5
        character (Len = StringLength)          :: PropertyName
        integer                                 :: WorkILB, WorkIUB
        integer                                 :: WorkJLB, WorkJUB
        integer                                 :: WorkKLB, WorkKUB
        integer                                 :: STAT_CALL
        integer(4)                              :: HDF5_CREATE
        logical                                 :: Final_ = .false.
        character (Len = Pathlength)            :: filename
        
        !----------------------------------------------------------------------

        WorkILB = Me%WaterWorkSize3D%ILB 
        WorkIUB = Me%WaterWorkSize3D%IUB 
        WorkJLB = Me%WaterWorkSize3D%JLB 
        WorkJUB = Me%WaterWorkSize3D%JUB 
        WorkKLB = Me%WaterWorkSize3D%KLB 
        WorkKUB = Me%WaterWorkSize3D%KUB 

        if(present(Final)) Final_ = Final

        !Checks if it's at the end of the run 
        !or !if it's supposed to overwrite the final HDF file
        if (Final_ .or. Me%Output%RestartOverwrite) then

            filename = trim(Me%Files%Final)

        else

            filename =  ChangeSuffix(Me%Files%Final,                         &
                            "_"//trim(TimeToString(Me%ExternalVar%Now))//".fin")

        endif

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        ObjHDF5 = 0

        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,trim(filename)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR01'

        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR02'


        !Writes the Grid
        call HDF5WriteData   (ObjHDF5, "/Grid", "Bathymetry", "m",                      &
                              Array2D = Me%ExtWater%Bathymetry,                         &
                              STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR03'

        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                                 &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR10'

        call HDF5WriteData   (ObjHDF5, "/Grid", "WaterPoints3D", "-",                   &
                              Array3D = Me%ExtWater%WaterPoints3D,                      &
                              STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR04'

        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB+1, WorkJLB, WorkJUB+1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR05'

        call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionX", "m",                     &
                              Array2D = Me%ExternalVar%XX_IE,                           &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR06'

        call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionY", "m",                     &
                              Array2D = Me%ExternalVar%YY_IE,                           &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR07'

        !Writes SZZ
        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                        &
                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR08'

        call HDF5WriteData  (ObjHDF5, "/Grid", "VerticalZ",                             &
                             "m", Array3D = Me%ExtWater%SZZ,                            &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR09'

        !Writes OpenPoints
        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                                 &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR10'

        call HDF5WriteData  (ObjHDF5, "/Grid", "OpenPoints",                            &
                             "-", Array3D = Me%ExtWater%OpenPoints3D,                   &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR11'


        Property => Me%FirstProperty

        do while (associated(Property))
      
            PropertyName = trim(Property%ID%name)
            
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR12'

            !Final concentration
            call HDF5WriteData  (ObjHDF5, "/Mass/"//Property%ID%Name,                   &
                                 Property%ID%Name, Property%ID%Units,                   &
                                 Array2D = Property%Mass_Available,                     &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR13'

            Property => Property%Next

        enddo

        nullify (Property)
   
        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR14'

        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_Final_HDF - ModuleInterfaceSedimentWater - ERR15'

    end subroutine Write_Final_HDF

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjInterfaceSedimentWater_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceSedimentWater_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjInterfaceSedimentWater_ID > 0) then
            call LocateObjInterfaceSedimentWater (ObjInterfaceSedimentWater_ID)
            ready_ = VerifyReadLock (mINTERFACESEDIMENTWATER_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjInterfaceSedimentWater (ObjInterfaceSedimentWaterID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceSedimentWaterID

        !Local-----------------------------------------------------------------

        Me => FirstObjInterfaceSedimentWater
        do while (associated (Me))
            if (Me%InstanceID == ObjInterfaceSedimentWaterID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleInterfaceSedimentWater - LocateObjInterfaceSedimentWater - ERR01'

    end subroutine LocateObjInterfaceSedimentWater

    !--------------------------------------------------------------------------

end module ModuleInterfaceSedimentWater

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------




