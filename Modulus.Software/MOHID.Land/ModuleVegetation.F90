!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Land
! MODULE        : Vegetation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Mar 2006
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module to simulate plant characteristics
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
!         [Keyword]                 [Format]  [Units]  [Default]  [Short Description]
! VEGETATION_ID_FILE               : string      -        [-]     !Vegetation distribution grid path
! VEGETATION_DT                    : real        s      [ModelDT] !Vegetation DT
! INTEGRATION_DT                   : real        s      [ModelDT] !DT to integrate external variables until vegetation is
!                                                                 ! is called (vegetation DT)
! TRANSPIRATION_METHOD             : integer     -        [1]     !Plant transpiration method: 1-MOHID 
!                                                                   !(previous approach); 2- SWAT based
!   ROOT_PROFILE                   : integer     -        [1]     !Read if TRANSPIRATION_METHOD == 1: 
!                                                                   !1-Triangular; 2-Constant; 3-Exponential(SWAT like)
!   WATER_UPTAKE_STRESS_METHOD     : integer     -        [1]     !Read if TRANSPIRATION_METHOD == 1: 1-Feddes; 2-VanGenuchten
! ATMOSPHERE_CO2                   : real       ppm      [330.]   !Atmosphere CO2 concetrations - should be atmosphere property               
! WATER_UPTAKE_COMPENSATION_FACTOR : real        -        [0.]    !Factor for uptake compensation from lower layers if computed  
!                                                                   !layer demand is not met
!                                                                   !If zero there will exist no compensation. If 1. total demand  
!                                                                   !no met may come from lower layers
! NITROGEN_DISTRIBUTION_PARAMETER  : real                [-20.]
! PHOSPHORUS_DISTRIBUTION_PARAMETER: real                [-20.]
!
! TEMPERATURE_STRESS               : 0/1         -        [1]
! GRAZING                          : 0/1         -        [0]
! MANAGEMENT                       : 0/1         -        [0]
! DORMANCY                         : 0/1         -        [0]
! CHANGE_LAI_SENESCENCE            : 0/1         -        [0] 
! CHANGE_CANOPY_HEIGHT             : 0/1         -        [0]
!    
! <beginproperty>
!  See module fillmatrix
! EVOLUTION                        : integer     -         1      !Property evolution: 1-Read from file
!                                                                 !2-vegetation growth model
! <endproperty>
!
! <beginvegetationtype>
!  ID
!  NAME        
!  HAS_LEAVES
!  FEDDES_H1
!  FEDDES_H2
!  FEDDES_H3
!  FEDDES_H4

!   <begintimingdatabase>
!    MATURITY_HU                       : 1700.
!    PLANTING_JULIANDAY                : -99.
!    PLANTING_HUBASE                   : 0.15
!   <endgtimingdatabase>
!
!   <begingrowthdatabase>
!    PLANT_TYPE                        : 5
!    OPTIMAL_NITROGENFRACTION_N1       : 0.0663
!    OPTIMAL_NITROGENFRACTION_N2       : 0.0255
!    OPTIMAL_NITROGENFRACTION_N3       : 0.0148
!    OPTIMAL_PHOSPHORUSFRACTION_P1     : 0.0053
!    OPTIMAL_PHOSPHORUSFRACTION_P2     : 0.0020
!    OPTIMAL_PHOSPHORUSFRACTION_P3     : 0.0012
!    BASE_TEMPERATURE                  : 0.
!    OPTIMAL_TEMPERATURE               : 18.0
!    RADIATION_EXTINCTION_COEF         : 0.65
!    BIOMASS_ENERGY_RATIO              : 30.0
!    CO2_HIGH                          : 660.0
!    BIOMASS_ENERGY_RATIO_HIGH         : 39.0
!    RUE_DECLINE_RATE                  : 6.0
!    LAI_MAX                           : 4.0
!    OPTIMAL_LAIMAXFRACTION_1          : 0.05
!    OPTIMAL_LAIMAXFRACTION_2          : 0.95
!    GROWFRACTION_1                    : 0.05
!    GROWFRACTION_2                    : 0.45
!    GROWFRACTION_LAIDECLINE           : 0.50
!    ROOT_DEPTH_MAX                    : 1.30
!    CANOPY_HEIGHT_MAX                 : 0.9
!    OPTIMAL_HARVEST_INDEX             : 0.4
!    MINIMUM_HARVEST_INDEX             : 0.2
!    YELD_NITROGENFRACTION             : 0.0250
!    YELD_PHOSPHORUSFRACTION           : 0.0022
!   <endgrowthdatabase>
!        
!   <beginmanagementandgrazedatabase>
!    GRAZING_START_JULIANDAY           : -99.
!    GRAZING_START_PLANTHU             : 0.5
!    GRAZING_DAYS                      : 10
!    MINIMUM_BIOMASS_FOR_GRAZING       : 10.
!    GRAZING_BIOMASS                   : 70.
!    TRAMPLING_BIOMASS                 : 30.
!    HARVESTKILL_JULIANDAY             : -99.
!    HARVESTKILL_PLANTHU               : 1.2
!    HARVEST_JULIANDAY                 : -99.
!    HARVEST_PLANTHU                   : -99.
!    HARVEST_EFFICIENCY                : 1.0
!    KILL_JULIANDAY                    : -99.
!    KILL_PLANTHU                      : -99.
!   <endmanagementandgrazedatabase>
!      
! <endvegetationtype>
            


Module ModuleVegetation

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5,           only : ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,             &
                                     HDF5WriteData, HDF5FlushMemory, HDF5ReadData, KillHDF5
    use ModuleAtmosphere,     only : GetAtmosphereProperty, AtmospherePropertyExists, UnGetAtmosphere
    use ModuleBasinGeometry,  only : GetBasinPoints, UngetBasin
    use ModulePorousMedia,    only : GetWaterContent, GetHead, GetThetaR, GetComputeSoilField,    &
                                     GetThetaField, GetLimitThetaLow, GetUnsatK, UngetPorousMedia
    use ModuleFunctions,      only : ConstructPropertyID, LinearInterpolation, InterpolateValueInTime
    use ModuleHorizontalGrid, only : GetHorizontalGridSize, GetGridCellArea, WriteHorizontalGrid, &
                                     UngetHorizontalGrid
    use ModuleHorizontalMap,  only : GetOpenPoints2D, UngetHorizontalMap
    use ModuleFillMatrix,     only : ConstructFillMatrix, GetDefaultValue, ModifyFillMatrix,      &
                                     KillFillMatrix, GetIfMatrixRemainsConstant
    use ModuleTimeSerie,      only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,               &
                                     StartTimeSerieInput, GetTimeSerieValue
    use ModuleGridData,       only : ConstructGridData, GetGridData, UngetGridData,               &
                                     KillGridData
    use ModuleGeometry,       only : GetGeometryDistances, GetGeometrySize,                       &
                                     GetGeometryVolumes, UngetGeometry           

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructVegetation
    private ::      AllocateInstance
    private ::      ReadVegetationFileNames           
    private ::      ConstructGlobalVariables
    private ::      ConstructPropertyList             
    private ::          ConstructProperty             
    private ::              ConstructPropertyID       
    private ::              ConstructPropertyValues
    private ::              ConstructPropertyEvolution   
    private ::              ConstructPropertyOutput  
    private ::          AddProperty
    private ::      CheckOptionsConsistence 
    private ::      ConstructVegetationGrid
    private ::      ConstructVegetationParameters
    private ::          ReadTimingDatabase
    private ::          ReadGrowthDatabase
    private ::          ReadManagementAndGrazeDatabase           !new
    private ::      ConstructTimeSerie
    private ::      ConstructHDF
    private ::          Open_HDF5_Output_File 
    private ::      ConstructLog                 

    !Selector
    public  :: GetLeafAreaIndex
    public  :: GetSpecificLeafStorage
    public  :: GetEVTPCropCoefficient
    public  :: GetRootDepth
    public  :: GetTranspiration
!    public  :: GetFeddesH
!    public  :: GetVegetationRootProfile
    public  :: UngetVegetation
    
    !Modifier
    public  :: ModifyVegetation
                    !Use model for property evolution
    public  ::      CheckPlantState
    public  ::      ModifyFluxes
    public  ::          PlantRootFluxes
    public  ::              WaterUptakeSWAT
    public  ::              NitrogenUptakeSWAT
    public  ::              NitrogenFixationSWAT
    public  ::              PhosphorusUptakeSWAT
    public  ::          PlantAerialFluxes
    public  ::              BiomassGrowthFromRadiationSWAT
    public  ::              LAIGrowthSWAT
    public  ::          GrazingFluxes
    public  ::          ManagementFluxes
    public  ::              HarvestOperation
    public  ::              HarvestKillOperation
    public  ::              KillOperation
    public  ::      ModifyModelledProperties
    public  ::          UpdateGlobalPlantProperties
    public  ::          UpdateRootProperties
    public  ::          UpdateStemProperties
    public  ::          UpdateLeafProperties
                    !Read property evolution
    public  ::      ModifyReadedProperties

    public  ::      Modify_OutputHDF
    public  ::      Modify_OutputTimeSeries

    !Destructor
    public  :: KillVegetation                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjVegetation 
    
    !Interfaces----------------------------------------------------------------
    private :: UngetVegetation2D
    interface  UngetVegetation
        module procedure UngetVegetation2D
    end interface UngetVegetation

    !Parameters----------------------------------------------------------------
    !Property Evolution
    integer, parameter                          :: ReadValue                     = 1
    integer, parameter                          :: SWAT                          = 2
    integer, parameter                          :: DEB                           = 3
    !Transpiration Method               
    integer, parameter                          :: TranspirationMOHID            = 1
    integer, parameter                          :: TranspirationSWAT             = 2
    !RootProfile
    integer, parameter                          :: RootTriangular                = 1
    integer, parameter                          :: RootConstant                  = 2
    integer, parameter                          :: RootExponential               = 3
    !MohidTranspiration Water Uptake Factor
    integer, parameter                          :: Feddes                        = 1
    integer, parameter                          :: Genuchten                     = 2

    character(StringLength), parameter          :: beginproperty                 = '<beginproperty>'
    character(StringLength), parameter          :: endproperty                   = '<endproperty>'

    !Types---------------------------------------------------------------------
    type       T_ID
        logical                                         :: SolutionFromFile
        integer                                         :: IDNumber
        character(LEN = StringLength)                   :: Name
        character(LEN = StringLength)                   :: Description
        character(LEN = StringLength)                   :: Units
    end type T_ID
    
    type       T_Integration
        real,    dimension(:,:  ), pointer              :: AveragePotTPDuringDT
        real,    dimension(:,:  ), pointer              :: AverageAirTempDuringDT
        real,    dimension(:,:  ), pointer              :: AverageAirHumidityDuringDT
        real,    dimension(:,:  ), pointer              :: AverageRadiationDuringDT
        real,    dimension(:,:  ), pointer              :: SumPotTP
        real,    dimension(:,:  ), pointer              :: SumTemperature
        real,    dimension(:,:  ), pointer              :: SumHumidity
        real,    dimension(:,:  ), pointer              :: SumRadiation
    end type   T_Integration
    
    type       T_External
        real                                            :: DT
        type(T_Time)                                    :: Now
        integer                                         :: JulianDay_Old
        logical                                         :: ComputeSoilField
        logical                                         :: CoupledAtmosphere
        real, dimension(:,:,:), pointer                 :: DWZ
        real, dimension(:,:,:), pointer                 :: SZZ
        integer, dimension(:,:  ), pointer              :: OpenPoints2D         => null()
        integer, dimension(:,:  ), pointer              :: MappingPoints2D      => null()
        real, dimension(:,:  ), pointer                 :: Topography           => null()
        real,    dimension(:,:  ), pointer              :: AirTemperature                 !ºC
        real,    dimension(:,:  ), pointer              :: SolarRadiation                 !W/m2
        real,    dimension(:,:  ), pointer              :: RelativeHumidity               !0-1
        real,    dimension(:,:  ), pointer              :: PotentialTranspiration         !m/s     
        real,    dimension(:,:,:), pointer              :: SoilWaterContent               !m3H2O/m3soil
        real,    dimension(:,:,:), pointer              :: Head
        real,    dimension(:,:,:), pointer              :: ResidualWaterContent                          
        real,    dimension(:,:,:), pointer              :: SoilNitrogenConcentration      !kg/m3H2O??    !Ver como fazer o Get
        real,    dimension(:,:,:), pointer              :: SoilPhosphorusConcentration    !kg/m3H2O??    !Ver como fazer o Get
        real,    dimension(:,:,:), pointer              :: FieldCapacity                  !m3H2O/m3soil  !criar no ModulePorousMedia
        real,    dimension(:,:  ), pointer              :: GridCellArea
        real(8), dimension(:,:,:), pointer              :: CellVolume
        type(T_Integration)                             :: Integration
    end type   T_External


    type       T_Files
        character(len=Pathlength)                       :: ConstructData
        character(len=Pathlength)                       :: Results
        character(len=Pathlength)                       :: Initial
        character(PathLength)                           :: VegetationIDFile

    end type T_Files

    type       T_OutPut
        type (T_Time), pointer, dimension(:)            :: OutTime
        logical                                         :: HDF_ON
        integer                                         :: NextHDF5
        integer                                         :: NextOutput
        integer                                         :: Number
        logical                                         :: TimeSerie_ON
    end type T_OutPut
    
    type T_Property
        type(T_PropertyID)                              :: ID
        real, dimension(:,:), pointer                   :: Field
        logical                                         :: Old
        logical                                         :: OutputHDF
        logical                                         :: TimeSerie
        logical                                         :: BoxTimeSerie
        integer                                         :: Evolution            = ReadValue
        logical                                         :: IsConstant
        real                                            :: ConstantValue        = null_real
        integer                                         :: TimeSeriesID         = 0
        integer                                         :: TimeSeriesColumn
        character(PathLength)                           :: TimeSeriesName
        type(T_Property), pointer                       :: Next, Prev
    end type T_Property

    type T_TimingDatabase
        real                                            :: PlantHUatMaturity
        real                                            :: PlantingJulianDay
        real                                            :: PlantingHUBase
    end type T_TimingDatabase

    type T_GrowthDatabase
        real                                            :: PlantType

        real                                            :: PlantFractionN1
        real                                            :: PlantFractionN2
        real                                            :: PlantFractionN3            
        real                                            :: PlantFractionP1
        real                                            :: PlantFractionP2
        real                                            :: PlantFractionP3  
        real                                            :: PlantBaseTemperature
        real                                            :: PlantOptimalTemperature
        real                                            :: ExtinctCoef
        real                                            :: BiomassEnergyRatio
        real                                            :: CO2ConcHigh
        real                                            :: BiomassEnergyRatioHigh
        real                                            :: RUEDeclineRate
        real                                            :: MaximumRootDepth
        real                                            :: FrLAIMax1
        real                                            :: FrLAIMax2
        real                                            :: FrGrow1
        real                                            :: FrGrow2
        real                                            :: FrGrowLAIDecline
        real                                            :: LAIMax
        real                                            :: MaxCanopyHeight
        real                                            :: OptimalHarvestIndex
        real                                            :: MinimumHarvestIndex
        real                                            :: NitrogenFractionInYeld
        real                                            :: PhosphorusFractionInYeld
        real                                            :: TreeYearsToMaturity
        real                                            :: TreeMaximumBiomass
    end type T_GrowthDatabase

    type T_ManagementDatabase
        real                                            :: GrazingStartJulianDay
        real, dimension (:), pointer                    :: GrazingStartPlantHU
        integer                                         :: GrazingDays
        real                                            :: GrazingMinimumBiomass
        real                                            :: GrazingBiomass
        real                                            :: TramplingBiomass
        real                                            :: HarvestKillJulianDay
        real                                            :: HarvestKillPlantHU
        real                                            :: HarvestJulianDay
        real, dimension (:), pointer                    :: HarvestPlantHU
        real                                            :: KillJulianDay
        real                                            :: KillPlantHU
        real                                            :: HarvestEfficiency
    end type T_ManagementDatabase

    type       T_TranspirationMOHID
        real, dimension(:,:), pointer                   :: RootFeddesH1
        real, dimension(:,:), pointer                   :: RootFeddesH2
        real, dimension(:,:), pointer                   :: RootFeddesH3
        real, dimension(:,:), pointer                   :: RootFeddesH4
    end type   T_TranspirationMOHID

    type T_Evolution
        logical                                         :: ReadNeeded     = .false.
        logical                                         :: GrowthModelNeeded    = .false.
        logical                                         :: ModelSWAT      = .false.
        logical                                         :: ModelDEB       = .false.
    end type T_Evolution 

    type       T_ComputeOptions
        logical                                         :: Grazing              = .false.
        logical                                         :: Management           = .false.
        logical                                         :: Dormancy             = .false.
        integer                                         :: RootProfile          
        integer                                         :: TranspirationMethod
        logical                                         :: LimitTPVel
        logical                                         :: WaterUptakeOld  
        real                                            :: VegetationDT
        real                                            :: IntegrationDT 
        logical                                         :: IsPlantGrowing
        logical                                         :: ChangeLAISenescence
        logical                                         :: ChangeCanopyHeight
        real                                            :: PlantingHUBase, PlantHUatMaturity
        real                                            :: PotentialHUBase, PotentialHUTotal
        real                                            :: WaterUptakeCompensationFactor
        real                                            :: NitrogenDistributionParameter
        real                                            :: PhosphorusDistributionParameter
        real                                            :: AtmosphereCO2
        integer                                         :: WaterUptakeStressMethod
        logical                                         :: ModelNitrogen, ModelPhosphorus
        logical                                         :: ModelRootBiomass, ModelPlantBiomass
        logical                                         :: ModelCanopyHeight, ModelTemperatureStress
        type(T_TranspirationMOHID)                      :: TranspirationMOHID
        type(T_Evolution)                               :: Evolution
    end type   T_ComputeOptions

    type       T_HeatUnits
        real, dimension(:,:), pointer                   :: PotentialHUTotal
        real, dimension(:,:), pointer                   :: PotentialHUBase
        real                                            :: PotentialHUBase_Old
        real, dimension(:,:), pointer                   :: PlantHUAccumulated
        real, dimension(:,:), pointer                   :: PlantHUAccumulated_Old
    end type   T_HeatUnits

    type       T_Growth
        real, dimension(:,:), pointer                   :: WaterStress
        real, dimension(:,:), pointer                   :: NitrogenStress
        real, dimension(:,:), pointer                   :: PhosphorusStress
        real, dimension(:,:), pointer                   :: TemperatureStress
        real, dimension(:,:), pointer                   :: GlobalStress
        real, dimension(:,:), pointer                   :: BiomassGrowthOld
        integer, dimension(:,:), pointer                :: TreeCurrentYear
        real, dimension(:,:), pointer                   :: TreeFractionToMaturity
        real, dimension(:,:), pointer                   :: TreeMaximumAnnualBiomass
    end type   T_Growth


    type T_VegetationType
        integer                                         :: ID
        character(StringLength)                         :: Name
        type(T_VegetationType), pointer                 :: Next, Prev
        type(T_VegetationType), pointer                 :: FirstVegetation
        type(T_VegetationType), pointer                 :: LastVegetation
        integer                                         :: VegetationsNumber     = 0
        real                                            :: RootFeddesH1
        real                                            :: RootFeddesH2
        real                                            :: RootFeddesH3
        real                                            :: RootFeddesH4
        type (T_TimingDatabase)                         :: TimingDatabase
        type (T_GrowthDatabase)                         :: GrowthDatabase
        type (T_ManagementDatabase)                     :: ManagementDatabase
        logical                                         :: ComputeRoot
        logical                                         :: ComputeStem
        logical                                         :: HasLeaves

    end type T_VegetationType

    type T_FluxesToSoil
        real, dimension(:,:), pointer                   :: GrazingBiomassToSoil
        real, dimension(:,:), pointer                   :: GrazingNitrogenToSoil
        real, dimension(:,:), pointer                   :: GrazingPhosphorusToSoil
        real, dimension(:,:), pointer                   :: ManagementBiomassToSoil
        real, dimension(:,:), pointer                   :: ManagementNitrogenToSoil
        real, dimension(:,:), pointer                   :: ManagementPhosphorusToSoil
        real, dimension(:,:), pointer                   :: ManagementRootBiomassLeftInSoil  !In Kill operation.
        real, dimension(:,:), pointer                   :: DormancyBiomassToSoil
        real, dimension(:,:), pointer                   :: DormancyNitrogenToSoil
        real, dimension(:,:), pointer                   :: DormancyPhosphorusToSoil
    end type T_FluxesToSoil
    
    type T_FluxesFromSoil
        real, dimension(:,:), pointer                   :: WaterUptakeFromSoil
        real, dimension(:,:), pointer                   :: NitrogenUptakeFromSoil
        real, dimension(:,:), pointer                   :: PhosphorusUptakeFromSoil
    end type T_FluxesFromSoil

    type T_Fluxes
        real, dimension(:,:  ), pointer                 :: BiomassGrazed
        real, dimension(:,:  ), pointer                 :: NitrogenGrazed
        real, dimension(:,:  ), pointer                 :: PhosphorusGrazed
        real, dimension(:,:  ), pointer                 :: BiomassGrazedFraction
        real, dimension(:,:,:), pointer                 :: WaterUptakeLayer     !m/s
        real, dimension(:,:  ), pointer                 :: WaterUptake          !m/s
        real, dimension(:,:,:), pointer                 :: NitrogenUptakeLayer
        real, dimension(:,:  ), pointer                 :: NitrogenUptake
        real, dimension(:,:,:), pointer                 :: PhosphorusUptakeLayer
        real, dimension(:,:  ), pointer                 :: PhosphorusUptake
        real, dimension(:,:  ), pointer                 :: BiomassGrowth
        real, dimension(:,:  ), pointer                 :: LAIGrowth
        real, dimension(:,:  ), pointer                 :: LAIDeclineFraction
        real, dimension(:,:  ), pointer                 :: BiomassRemovedInHarvest
        real, dimension(:,:  ), pointer                 :: NitrogenRemovedInHarvest
        real, dimension(:,:  ), pointer                 :: PhosphorusRemovedInHarvest
        real, dimension(:,:  ), pointer                 :: BiomassHarvestedFraction
        real, dimension(:,:  ), pointer                 :: BiomassRemovedInDormancy
        real, dimension(:,:  ), pointer                 :: NitrogenRemovedInDormancy
        real, dimension(:,:  ), pointer                 :: PhosphorusRemovedInDormancy
        type(T_FluxesToSoil)                            :: ToSoil
        type(T_FluxesFromSoil)                          :: FromSoil
    end type T_Fluxes

    type T_StateVariables
        type(T_Property), pointer                       :: TotalPlantNitrogen
        type(T_Property), pointer                       :: TotalPlantPhosphorus
        type(T_Property), pointer                       :: TotalPlantBiomass
        type(T_Property), pointer                       :: RootBiomass
        type(T_Property), pointer                       :: RootDepth
        type(T_Property), pointer                       :: LeafAreaIndex
        type(T_Property), pointer                       :: SpecificLeafStorage
        type(T_Property), pointer                       :: EVTPCropCoefficient
        
        type(T_Property), pointer                       :: CanopyHeight


    end type T_StateVariables

    private :: T_Vegetation
    type       T_Vegetation
        integer                                         :: InstanceID
        type(T_Size3D)                                  :: Size, WorkSize
        type(T_External)                                :: ExternalVar
        type(T_Files)                                   :: Files
        type(T_OutPut)                                  :: OutPut
        type(T_Time     )                               :: BeginTime
        type(T_Time     )                               :: EndTime
        type(T_Time     )                               :: ActualTime
        type(T_Time     )                               :: NextCompute
        type(T_Time     )                               :: NextIntegration
        type(T_Time     )                               :: LastOutPutHDF5
        type(T_Property    ), pointer                   :: FirstProperty
        type(T_Property    ), pointer                   :: LastProperty
        integer                                         :: PropertiesNumber     = 0

       
        !DataMatrixes
        integer, dimension(:,:), pointer                :: VegetationID
        integer                                         :: ObjEnterData         = 0
        integer                                         :: ObjTime              = 0
        integer                                         :: ObjGridData          = 0
        integer                                         :: ObjHorizontalGrid    = 0
        integer                                         :: ObjHorizontalMap     = 0
        integer                                         :: ObjGeometry          = 0
        integer                                         :: ObjBasinGeometry     = 0
        integer                                         :: ObjTimeSerie         = 0
        integer                                         :: ObjFillMatrix        = 0
        integer                                         :: ObjAtmosphere        = 0
        integer                                         :: ObjPorousMedia       = 0
        integer                                         :: ObjHDF5              = 0
        type(T_Vegetation), pointer                     :: Next
        
        logical, dimension(:,:), pointer                :: IsPlantGrowing       
        logical, dimension(:,:), pointer                :: PlantingOccurred      
        logical, dimension(:,:), pointer                :: KillOccurred        
        logical, dimension(:,:), pointer                :: HarvestOnlyOccurred   
        logical, dimension(:,:), pointer                :: HarvestKillOccurred   
        logical, dimension(:,:), pointer                :: IsPlantDormant       
        logical, dimension(:,:), pointer                :: PlantGoingDormant    
        logical, dimension(:,:), pointer                :: IsPlantBeingGrazed
        logical, dimension(:,:), pointer                :: LAISenescence        
        logical, dimension(:,:), pointer                :: HarvestFinished
        logical, dimension(:,:), pointer                :: GrazingFinished
        logical, dimension(:,:), pointer                :: ChangeCanopyEnabled      


        !global variables to be used in several routines
        real, dimension(:,:), pointer                   :: PlantNitrogenFraction
        real, dimension(:,:), pointer                   :: PlantPhosphorusFraction
        real, dimension(:,:), pointer                   :: PlantLAIMaxFraction
        real, dimension(:,:), pointer                   :: LAIBeforeSenescence
        !global counter
        real, dimension(:,:), pointer                   :: DaysOfGrazing                    !counter to days of grazing
        integer, dimension(:,:), pointer                :: HarvestOperations
        integer, dimension(:,:), pointer                :: GrazingOperations
        integer, dimension(:,:), pointer                :: HarvestKillOperations
        integer, dimension(:,:), pointer                :: KillOperations
        integer                                         :: nIterations                      !counter to atmosphere integration


        integer                                         :: VegetationsNumber     = 0
        type (T_VegetationType), dimension(:), pointer  :: VegetationTypes  => null()
        type(T_ComputeOptions)                          :: ComputeOptions
        type(T_HeatUnits)                               :: HeatUnits
        type(T_Growth)                                  :: Growth
        type(T_Fluxes)                                  :: Fluxes
        type(T_StateVariables)                          :: StateVariables
    end type  T_Vegetation

    !Global Module Variables
    type (T_Vegetation), pointer                        :: FirstObjVegetation  => null()
    type (T_Vegetation), pointer                        :: Me                  => null()


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructVegetation(ObjVegetationID,                 &
                                   TimeID,                          &
                                   GridDataID,                      &
                                   HorizontalGridID,                &
                                   HorizontalMapID,                 &
                                   AtmosphereID,                    &
                                   PorousMediaID,                   &
                                   MappingPoints,                   &
                                   GeometryID,                      &
                                   CoupledAtmosphere,               &
                                   STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjVegetationID 
        integer                                         :: TimeID         
        integer                                         :: GridDataID     
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: AtmosphereID
        integer                                         :: PorousMediaID
        integer, dimension(:, :), pointer               :: MappingPoints
        integer                                         :: GeometryID
        logical                                         :: CoupledAtmosphere
        integer, optional, intent(OUT)                  :: STAT     

        !Local-------------------------------------------------------------------
        integer                                         :: ready_         
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mVegetation_)) then
            nullify (FirstObjVegetation)
            call RegisterModule (mVegetation_) 
        endif

        call Ready(ObjVegetationID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjAtmosphere     = AssociateInstance (mATMOSPHERE_,     AtmosphereID    )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            
            Me%ExternalVar%MappingPoints2D    => MappingPoints
            Me%ExternalVar%CoupledAtmosphere  =  CoupledAtmosphere

            call ReadVegetationFileNames

            
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetation - ModuleVegetation - ERR01' 
            
            !Read Data file Options
            call ConstructGlobalVariables

            call ConstructPropertyList

            call ConstructVegetationParameters 
            
            call CheckOptionsConsistence

            
            !Grid operations
            call AllocateVariables
            
            call ConstructVegetationGrid

            if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then
                call ConstructFeddes
            endif

            !Output
            call ConstructTimeSerie

            call ConstructHDF


            call ConstructLog

            call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetation - ModuleVegetation - ERR02'


            !Returns ID
            ObjVegetationID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleVegetation - ConstructVegetation - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructVegetation
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Vegetation), pointer                         :: NewObjVegetation
        type (T_Vegetation), pointer                         :: PreviousObjVegetation


        !Allocates new instance
        allocate (NewObjVegetation)
        nullify  (NewObjVegetation%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjVegetation)) then
            FirstObjVegetation    => NewObjVegetation
            Me                    => NewObjVegetation
        else
            PreviousObjVegetation => FirstObjVegetation
            Me                    => FirstObjVegetation%Next
            do while (associated(Me))
                PreviousObjVegetation  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjVegetation
            PreviousObjVegetation%Next => NewObjVegetation
        endif

        Me%InstanceID = RegisterNewInstance (mVegetation_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ReadVegetationFileNames

        !External--------------------------------------------------------------
        character(len = StringLength)           :: Message
        integer                                 :: STAT_CALL

        !Begin------------------------------------------------------------------

        !Opens the Vegetation data file 
        ! ASCII file used to construct new properties
        call ReadFileName('VEGETATION_DATA', Me%Files%ConstructData, Message = 'Vegetation Data File', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadVegetationFileNames - ModuleVegetation - ERR01'

        ! ---> File in HDF format where is written instant fields of Vegetation properties
        Message   ='Instant fields of Vegetation properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('VEGETATION_HDF', Me%Files%Results, Message = 'Vegetation HDF File', &
                           TIME_END = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadVegetationFileNames - ModuleVegetation - ERR02'

    end subroutine ReadVegetationFileNames

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag
        real                                    :: auxFactor, Erroraux, DTaux, ModelDT

        !Begin------------------------------------------------------------------

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = Me%BeginTime,                 &
                                  EndTime = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR01'
        
        !Actualize the time
        Me%ActualTime  = Me%BeginTime
        Me%NextCompute = Me%ActualTime

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

        call GetGeometrySize    (Me%ObjGeometry,                                        &    
                                 Size     = Me%Size,                                    &
                                 WorkSize = Me%WorkSize,                                &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalDiscretization - ModulePorousMedia - ERR010'


        call GetData(Me%Files%VegetationIDFile,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'VEGETATION_ID_FILE',                             &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR020'

        call GetData(Me%ComputeOptions%TranspirationMethod,                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TRANSPIRATION_METHOD',                           &
                     Default        = TranspirationMOHID,                               &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR030'
        
        if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then        

            call GetData(Me%ComputeOptions%LimitTPVel,                                  &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword        = 'LIMIT_TRANSP_WATER_VEL',                     &
                         Default        = .false.,                                      &
                         SearchType     = FromFile,                                     &
                         ClientModule   = 'ModuleVegetation',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR040'
            
            
            call GetData(Me%ComputeOptions%RootProfile,                                 &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword        = 'ROOT_PROFILE',                               &
                         Default        = RootTriangular,                               &
                         SearchType     = FromFile,                                     &
                         ClientModule   = 'ModuleVegetation',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR050'

            call GetData(Me%ComputeOptions%WaterUptakeOld,                              &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword        = 'WATER_UPTAKE_OLD',                           &
                         Default        = .false.,                                      &
                         SearchType     = FromFile,                                     &
                         ClientModule   = 'ModuleVegetation',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR060'        

            call GetData(Me%ComputeOptions%WaterUptakeStressMethod,                     &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword        = 'WATER_UPTAKE_STRESS_METHOD',                 &
                         Default        = Feddes,                                       &
                         SearchType     = FromFile,                                     &
                         ClientModule   = 'ModuleVegetation',                           &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR070'
        endif

        call GetData(Me%ComputeOptions%AtmosphereCO2,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ATMOSPHERE_CO2',                                 &
                     Default        = 330.,                                             &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR080'

        call GetData(Me%ComputeOptions%WaterUptakeCompensationFactor,                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'WATER_UPTAKE_COMPENSATION_FACTOR',               &
                     Default        = 0.,                                               &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR090'

        call GetData(Me%ComputeOptions%NitrogenDistributionParameter,                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'NITROGEN_DISTRIBUTION_PARAMETER',                &
                     Default        = 20.,                                              &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0100'

        call GetData(Me%ComputeOptions%PhosphorusDistributionParameter,                 &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'PHOSPHORUS_DISTRIBUTION_PARAMETER',              &
                     Default        = 20.,                                              &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0110'
    
        call GetData(Me%ComputeOptions%ModelTemperatureStress,                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TEMPERATURE_STRESS',                             &
                     Default        = .true.,                                           &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0120'
        
        !Tests with initial conditions for one cell. Model in the future should try to contruct initial conditions for grid
        call GetData(Me%ComputeOptions%IsPlantGrowing,                                  &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ISPLANTGROWING',                                 &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0130'
        call GetData(Me%ComputeOptions%PotentialHUTotal,                                &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'POTENTIALHUTOTAL',                               &
                     Default        = 0.,                                               &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0140'
        call GetData(Me%ComputeOptions%PotentialHUBASE,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'POTENTIALHUBASE',                                &
                     Default        = 0.,                                               &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0150'


        call GetData(Me%ComputeOptions%Grazing,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'GRAZING',                                        &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0160'

        call GetData(Me%ComputeOptions%Management,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'MANAGEMENT',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0170'
  
        call GetData(Me%ComputeOptions%Dormancy,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'DORMANCY',                                       &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0180'

        call GetData(Me%ComputeOptions%ChangeLAISenescence,                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'CHANGE_LAI_SENESCENCE',                          &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0190'

        call GetData(Me%ComputeOptions%ChangeCanopyHeight,                              &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'CHANGE_CANOPY_HEIGHT',                           &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0200'


        call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR0210'

        Me%ExternalVar%DT = ModelDT

        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0220'

        call GetData(Me%ComputeOptions%VegetationDT,                                     &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'VEGETATION_DT',                                     &
                     Default      = ModelDT,                                             &
                     SearchType   = FromFile,                                            &
                     ClientModule = 'ModuleVegetation',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR0230'
                                   
        call GetData(Me%ComputeOptions%IntegrationDT,                                    &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'INTEGRATION_DT',                                    &
                     Default      = ModelDT,                                             &
                     SearchType   = FromFile,                                            &
                     ClientModule = 'ModuleVegetation',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR240'

        if (Me%ComputeOptions%VegetationDT .lt. Me%ComputeOptions%IntegrationDT) then
            write(*,*) 
            write(*,*) 'Vegetation DT time step is smaller then the external variables integration time step'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR250'
        endif
        
        if (Me%ComputeOptions%IntegrationDT .lt. Me%ExternalVar%DT) then
            write(*,*) 
            write(*,*) 'Integration DT time step is smaller then the model time step'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR260'  
        endif      
        
        if (Me%ComputeOptions%VegetationDT .lt. ModelDT) then
            write(*,*) 
            write(*,*) 'Vegetation DT time step is smaller then model time step'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR270'

        else

            !Vegetation time step must be a multiple of the model time step
            auxFactor = Me%ComputeOptions%VegetationDT  / ModelDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) 'Vegetation time step must be a multiple of model time step.'
                write(*,*) 'Please review your input data.'
                stop 'ConstructGlobalVariables - ModuleVegetation - ERR280'
            endif

            !Run period in seconds
            DTaux = Me%EndTime - Me%ExternalVar%Now

            !The run period   must be a multiple of the Vegetation DT
            auxFactor = DTaux / Me%ComputeOptions%VegetationDT

            ErrorAux = auxFactor - int(auxFactor)
            if (ErrorAux /= 0) then

                write(*,*) 
                write(*,*) 'Model run time is not a multiple of vegetation time step.'
                write(*,*) 'Please review your input data.'
                stop 'ConstructGlobalVariables - ModuleVegetation - ERR290'
            end if
        endif

        Me%NextCompute = Me%ExternalVar%Now + Me%ComputeOptions%VegetationDT

        Me%NextIntegration = Me%ExternalVar%Now + Me%ComputeOptions%IntegrationDT

    end subroutine ConstructGlobalVariables

    !--------------------------------------------------------------------------

    subroutine ConstructPropertyList

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
                                  
        integer :: STAT_CALL
        integer :: ClientNumber
        logical :: BlockFound
                                  
        !Local-----------------------------------------------------------------

        type (T_Property), pointer :: NewProperty

        !----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                &
                                        beginproperty, endproperty, BlockFound,       &
                                        STAT = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_      ) then    
cd2 :           if (BlockFound) then                                  
                
                    ! Construct a New Property 
                    call ConstructProperty(NewProperty)
                                             

                    ! Add new Property to the Property List 
                    call AddProperty(NewProperty)

                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                      &
                        stop 'Subroutine ConstructPropertyList; ModuleVegetation. ERR01.'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructPropertyList - ModuleVegetation - ERR02'
            else cd1
                write(*,*)  
                write(*,*) 'Vegetation Properties are now defined as MOHID'
                write(*,*) 'standards. For support verify options at MOHID'
                write(*,*) 'wiki in http://www.mohid.com/wiki'
                stop 'ConstructPropertyList - ModuleVegetation - ERR03'
            end if cd1
        end do do1


    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new property.           
    subroutine ConstructProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !----------------------------------------------------------------------
             
        allocate (NewProperty)

        nullify(NewProperty%Field                )
        nullify(NewProperty%Prev,NewProperty%Next)
        

        call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call ConstructPropertyValues    (NewProperty)

        call ConstructPropertyEvolution (NewProperty)

        call ConstructPropertyOutPut    (NewProperty)


    end subroutine ConstructProperty

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct the property values       
    ! in the domain and in the boundaries            
    subroutine ConstructPropertyValues(NewProperty)
        !Arguments-------------------------------------------------------------
        type(T_property),              pointer      :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        integer                                     :: ILB,IUB
        integer                                     :: JLB,JUB
        
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB


        allocate(NewProperty%Field(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleVegetation - ERR03'
        NewProperty%Field(:,:) = FillValueReal
        
        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleVegetation',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleVegetation - ERR05'


        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not.NewProperty%Old) then
            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExternalVar%MappingPoints2D,   &
                                       Matrix2D             = NewProperty%Field,                &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'ConstructPropertyValues - ModuleVegetation - ERR07'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructPropertyValues - ModuleVegetation - ERR08'
            end if

            call CheckFieldConsistence (NewProperty)

        else

            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Vegetation initial file written in HDF format  
            call ReadOldHDF(NewProperty)

        end if


    end subroutine ConstructPropertyValues

    !--------------------------------------------------------------------------
    
    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer               :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: i,j
        integer                                 :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
            
        !Verification if the values read are lower than zero in water points
        do I = ILB, IUB
        do J = JLB, JUB
            
            if (Me%ExternalVar%MappingPoints2D(i, j) == BasinPoint) then
                               
                if (NewProperty%Field(i, j) < 0.) then

                    NewProperty%Field(i, j) = 0.
                    
                end if

            else

                NewProperty%Field(i, j) = FillValueReal

            endif
        
        enddo
        enddo

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------

    !If the user wants to use the values of a previous   
    ! run, read the property values from the final      
    ! results file. By default this      
    ! file is in HDF format                                

    subroutine ReadOldHDF(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character (Len=StringLength)                :: PropertyName
        logical                                     :: EXIST
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        !----------------------------------------------------------------------


        inquire (FILE=trim(Me%Files%Initial)//"5", EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%Initial)//"5",&
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldHDF - ModuleVegetation - ERR01'


            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Field(:,:) = FillValueReal


            ! Reads from HDF file the Property value and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldHDF - ModuleVegetation - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Results/"//NewProperty%ID%Name,              &
                                 NewProperty%ID%Name,                                    &
                                 Array2D = NewProperty%Field,                            &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldHDF - ModuleVegetation - ERR03'


            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldHDF - ModuleVegetation - ERR06'

        else
            
            write(*,*)
            stop 'ReadOldHDF - ModuleVegetation - ERR07'

        end if cd0

    end subroutine ReadOldHDF

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct the property                          
    ! evolution parameters             
    subroutine ConstructPropertyEvolution(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        !----------------------------------------------------------------------

        call GetData(NewProperty%Evolution,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'EVOLUTION',                                         &
                     ClientModule = 'ModuleVegetation',                                  &
                     Default      = ReadValue,                                           &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR01'

        if (NewProperty%Evolution == DEB) then
            write(*,*) 
            write(*,*) 'DEB vegetation model is not yet implemented.'
            write(*,*) 'Check EVOLUTION keyword options'
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR010'
        
        endif
        
        if (NewProperty%Evolution == ReadValue) then
            Me%ComputeOptions%Evolution%ReadNeeded = .true.
        elseif (NewProperty%Evolution == SWAT) then
            Me%ComputeOptions%Evolution%GrowthModelNeeded = .true.
            Me%ComputeOptions%Evolution%ModelSWAT   = .true.
        elseif (NewProperty%Evolution == DEB) then
            Me%ComputeOptions%Evolution%GrowthModelNeeded = .true.
            Me%ComputeOptions%Evolution%ModelDEB    = .true.
        else
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Option for vegetation evolution' 
            write(*,*    ) 'not known. For now it is possible to run a Vegetation'
            write(*,*    ) 'model (SWAT) or read from file. Check EVOLUTION keyword '
            stop 'ConstructPropertyEvolution - ModuleVegetation - ERR020'  

        endif

        call GetData(NewProperty%IsConstant,                                             &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'REMAIN_CONSTANT',                                   &
                     ClientModule = 'ModuleVegetation',                                  &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR030'

        if (.not. NewProperty%Evolution == ReadValue .and. NewProperty%IsConstant) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Property ', trim(NewProperty%ID%Name) 
            write(*,*    ) 'in Vegetation file has model and constant value options active.'
            write(*,*    ) 'Check EVOLUTION and REMAIN_CONSTANT keywords. '
            stop 'ConstructPropertyEvolution - ModuleVegetation - ERR040'  
        endif

        if (.not. NewProperty%Evolution == ReadValue .and. NewProperty%ID%SolutionFromFile) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Property ', trim(NewProperty%ID%Name)
            write(*,*    ) 'in Vegetation file has model and read value options active.'
            write(*,*    ) 'Check EVOLUTION and FILE_IN_TIME keywords. '
            stop 'ConstructPropertyEvolution - ModuleVegetation - ERR050'  
        endif

        if (NewProperty%IsConstant) then
            call GetData(NewProperty%ConstantValue,                                      &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DEFAULTVALUE',                                  &
                         ClientModule = 'ModuleVegetation',                              &
                         Default      = 0.,                                              &
                         STAT         = STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_)                                                  &
                stop 'Construct_PropertyEvolution - ModuleVegetation - ERR060'
        endif
            

    end subroutine ConstructPropertyEvolution

    !-------------------------------------------------------------------------

    subroutine ConstructPropertyOutPut(NewProperty)

        !Arguments------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: iflag

        !Begin----------------------------------------------------------------

        call GetData(NewProperty%OutputHDF,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'OUTPUT_HDF',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPropertyOutPut - ModuleVegetation - ERR01'
           
        call GetData(NewProperty%TimeSerie,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPropertyOutPut - ModuleVegetation - ERR02'


        call GetData(NewProperty%BoxTimeSerie,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOX_TIME_SERIE',                                 &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleAtmosphere',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPropertyOutPut - ModuleVegetation - ERR03'

    end subroutine ConstructPropertyOutPut

    !--------------------------------------------------------------------------

    ! This subroutine adds a new property to the Water Property List  
    subroutine AddProperty(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),              pointer     :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstProperty)) then
            Me%PropertiesNumber     = 1
            Me%FirstProperty        => NewProperty
            Me%LastProperty         => NewProperty
        else
            NewProperty%Prev        => Me%LastProperty
            Me%LastProperty%Next    => NewProperty
            Me%LastProperty         => NewProperty
            Me%PropertiesNumber     = Me%PropertiesNumber + 1
        end if 


    end subroutine AddProperty 

   !--------------------------------------------------------------------------

    subroutine CheckOptionsConsistence

        !Local----------------------------------------------------------------
        integer                                 :: STAT_CALL
        type(T_Property), pointer               :: TotalPlantBiomass, TotalPlantNitrogen
        type(T_Property), pointer               :: TotalPlantPhosphorus, RootBiomass, RootDepth
        type(T_Property), pointer               :: LeafAreaIndex, SpecificLeafStorage, CanopyHeight
        type(T_Property), pointer               :: EVTPCropCoefficient
        logical                                 :: NitrogenFound, PhosphorusFound, RootBiomassFound
        logical                                 :: Grazing, Management, Dormancy
        !Begin----------------------------------------------------------------
        
        if(Me%ComputeOptions%Evolution%ModelSWAT .and. Me%ComputeOptions%Evolution%ModelDEB) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Properties Evolution can not be modeled' 
            write(*,*    ) 'by SWAT and DEB at the same time. Check EVOLUTION keyword in' 
            write(*,*    ) 'properties block.'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR000'
        endif

        call GetComputeSoilField(Me%ObjPorousMedia, Me%ExternalVar%ComputeSoilField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsConsistence - ModuleVegetation - ERR050'   
        
        if(Me%ComputeOptions%TranspirationMethod .eq. TranspirationSWAT .and. .not. Me%ExternalVar%ComputeSoilField) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! User defined method for water uptake implies soil field capacity' 
            write(*,*    ) 'computation. Check keyword COMPUTE_SOIL_FIELD in porous media file' 
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR010'
        endif
        
        !Not modeled mandatory properties. Needed in basin module
        call SearchProperty(SpecificLeafStorage, SpecificLeafStorage_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if(SpecificLeafStorage%Evolution /= ReadValue) then
                write(*,*    )
                write(*,*    ) 'specific leaf storage can not yet be modeled'
                write(*,*    ) 'Check property EVOLUTION keyword'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR021'
            endif                                                      
        elseif (STAT_CALL /= SUCCESS_)  then
            write(*,*    )
            write(*,*    ) 'specific leaf storage is a mandatory Property for'
            write(*,*    ) 'vegetation. Needed in basin processes'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR022'
        endif

        call SearchProperty(EVTPCropCoefficient, EVTPCropCoefficient_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if(EVTPCropCoefficient%Evolution /= ReadValue) then
                write(*,*    )
                write(*,*    ) 'crop coefficient can not yet be modeled'
                write(*,*    ) 'Check property EVOLUTION keyword'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR023'
            endif                                                      
        elseif (STAT_CALL /= SUCCESS_)  then
            write(*,*    )
            write(*,*    ) 'crop coefficent is a mandatory Property for'
            write(*,*    ) 'vegetation. Needed in basin processes'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR024'
        endif

        !Vegetation module properties
        !Must have this property - RootDepth for water uptake
        call SearchProperty(RootDepth, RootDepth_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if((RootDepth%Evolution == ReadValue .or. RootDepth%ID%SolutionFromFile)            &
                .and. Me%ComputeOptions%Evolution%ModelSWAT) then
                write(*,*    ) 
                write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                write(*,*    ) 'root depth must be modeled in this conditions.'
                write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR025'
            endif
        elseif (STAT_CALL /= SUCCESS_)  then
            write(*,*    )
            write(*,*    ) 'root depth is a mandatory Property for'
            write(*,*    ) 'vegetation. It is needed for water uptake'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR026'
        endif

        !Must have this property - LAI for potential transpiration
        call SearchProperty(LeafAreaIndex, LeafAreaIndex_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if((LeafAreaIndex%Evolution == ReadValue .or. LeafAreaIndex%ID%SolutionFromFile)        &
                .and. Me%ComputeOptions%Evolution%ModelSWAT) then
                write(*,*    ) 
                write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                write(*,*    ) 'leaf area index must be modeled in this conditions.'
                write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                 stop 'CheckOptionsConsistence - ModuleVegetation - ERR030'
            endif                                                      
        elseif (STAT_CALL /= SUCCESS_)  then
            write(*,*    )
            write(*,*    ) 'leaf area index is a mandatory Property for'
            write(*,*    ) 'vegetation. It is needed for potential transpiration'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR040'
        endif

        if (Me%ComputeOptions%Evolution%ModelSWAT) then

            NitrogenFound = .false.
            PhosphorusFound = .false.
            !Optional Property
            Me%ComputeOptions%ModelNitrogen = .false.
            call SearchProperty(TotalPlantNitrogen, TotalPlantNitrogen_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                Me%ComputeOptions%ModelNitrogen = .true.
                if(TotalPlantNitrogen%Evolution == ReadValue) then
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                    write(*,*    ) 'total plant nitrogen must be modeled in this conditions.'
                    write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR050'
                endif                                                      
            endif
        
            !Optional Property
            Me%ComputeOptions%ModelPhosphorus = .false.
            call SearchProperty(TotalPlantPhosphorus, TotalPlantPhosphorus_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                Me%ComputeOptions%ModelPhosphorus = .true.
                if(TotalPlantPhosphorus%Evolution == ReadValue) then
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                    write(*,*    ) 'total plant phosphorus must be modeled in this conditions.'
                    write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR060'
                endif                                                      
            endif
            
            !Optional
            Management       = Me%ComputeOptions%Management

            Me%ComputeOptions%ModelRootBiomass = .false.
            call SearchProperty(RootBiomass, RootBiomass_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                Me%ComputeOptions%ModelRootBiomass = .true.
                if(RootBiomass%Evolution == ReadValue) then
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                    write(*,*    ) 'root biomass must be modeled in this conditions.'
                    write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR070'
                elseif (STAT_CALL /= SUCCESS_ .and. Management) then
                    write(*,*    )
                    write(*,*    ) 'root biomass is a mandatory Property'
                    write(*,*    ) 'if you want to run vegetation model and want to model'
                    write(*,*    ) 'management practices (e.g. hatvest where aerial biomass must be known'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR071'

                endif                                                      
            endif
            
            NitrogenFound    = Me%ComputeOptions%ModelNitrogen
            PhosphorusFound  = Me%ComputeOptions%ModelPhosphorus
            RootBiomassFound = Me%ComputeOptions%ModelRootBiomass
            Grazing          = Me%ComputeOptions%Grazing
            Dormancy         = Me%ComputeOptions%Dormancy

            Me%ComputeOptions%ModelPlantBiomass = .false.
            call SearchProperty(TotalPlantBiomass, TotalPlantBiomass_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                Me%ComputeOptions%ModelPlantBiomass = .true.
                if(TotalPlantBiomass%Evolution == ReadValue) then
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                    write(*,*    ) 'total plant biomass must be modeled in this conditions.'
                    write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR080'
                endif                                                      
            elseif (STAT_CALL /= SUCCESS_ .and.(NitrogenFound .or. PhosphorusFound .or. RootBiomassFound .or. Grazing &
                                                .or. Management .or. Dormancy))  then
                write(*,*    )
                write(*,*    ) 'total plant biomass is a mandatory Property'
                write(*,*    ) 'if you want to run vegetation model and want to model'
                write(*,*    ) 'nitrogen, phosphorus, root biomass, grazing or management practices'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR090'
            endif

            !Optional
            Me%ComputeOptions%ModelCanopyHeight = .false.
            call SearchProperty(CanopyHeight, CanopyHeight_, .false., STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                Me%ComputeOptions%ModelCanopyHeight = .true.
                if(CanopyHeight%Evolution == ReadValue) then
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                    write(*,*    ) 'canopy height must be modeled in this conditions.'
                    write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR0100'
                endif                                                      
            endif
  
        endif

        !Atmosphere properties needed
        if (Me%ComputeOptions%Evolution%GrowthModelNeeded .and. .not. Me%ExternalVar%CoupledAtmosphere) then
            write(*,*    )
            write(*,*    ) 'Vegetation growth model needs atmosphere to be connected'
            write(*,*    ) 'Check atmosphere definition in basin data file'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR0101'
        endif

        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                write(*,*) 
                write(*,*) 'Please specify in Atmosphere, Property air temperature'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR110'
            endif

            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, RelativeHumidity_)) then
                write(*,*) 
                write(*,*) 'Please specify in Atmosphere, Property relative humidity'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR120'
            endif

            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_)) then
                write(*,*) 
                write(*,*) 'Please specify in Atmosphere, Property solar radiation'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR130'
            endif
        endif
        
        
    end subroutine CheckOptionsConsistence


   !--------------------------------------------------------------------------

    subroutine ConstructLog

        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------

        write(*, *)
        write(*, *)"------------------------ VEGETATION ----------------------"         
        write(*, *)                              

        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
            write(*,*    ) 'Using Vegetation Growth Model'
            write(*,*    ) '---Root Growth              : ON'
            write(*,*    ) '---LAI Growth               : ON'
            write(*,*    ) '---Water Uptake             : ON'
            
            if (Me%ComputeOptions%TranspirationMethod .eq. 1) then
                write(*,*    ) '   ---Uptake Method         : 1'
            else
                write(*,*    ) '   ---Uptake Method         : 2'
            endif
            if (Me%ComputeOptions%ModelPlantBiomass) then
                write(*,*    ) '---Biomass Growth           : ON'
            else
                write(*,*    ) '---Biomass Growth           : OFF'
            endif
            if (Me%ComputeOptions%ModelNitrogen) then
                write(*,*    ) '---Nitrogen Uptake          : ON'
            else
                write(*,*    ) '---Nitrogen Uptake          : OFF'
            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
                write(*,*    ) '---Phosphorus Uptake        : ON'
            else
                write(*,*    ) '---Phosphorus Uptake        : OFF'
            endif
            if (Me%ComputeOptions%ModelTemperatureStress) then
                write(*,*    ) '---Temperature Stress       : ON'
            else
                write(*,*    ) '---Temperature Stress       : OFF'
            endif
            if (Me%ComputeOptions%Management) then
                write(*,*    ) '---Management Operations    : ON'
            else
                write(*,*    ) '---Management Operations    : OFF'
            endif              
            if (Me%ComputeOptions%Grazing) then
                write(*,*    ) '---Grazing Operations       : ON'
            else
                write(*,*    ) '---Grazing Operations       : OFF'
            endif                                   
            if (Me%ComputeOptions%Dormancy) then
                write(*,*    ) '---Dormancy                 : ON'
                write(*,*    ) 
            else
                write(*,*    ) '---Dormancy                 : OFF'
                write(*,*    ) 
            endif                                   

        else
            write(*,*    ) 'Vegetation Growth Model not Used'
            write(*,*    ) '---Root readed from file'
            write(*,*    ) '---LAI  readed from file'
            write(*,*    ) '---Water Uptake             : ON'
            if (Me%ComputeOptions%TranspirationMethod .eq. 1) then
                write(*,*    ) '   ---UptakeMethod          : 1'
            else
                write(*,*    ) '   ---UptakeMethod          : 2'
            endif
            write(*,*    ) '---Biomass Growth           : OFF'
            write(*,*    ) '---Nitrogen Uptake          : OFF'
            write(*,*    ) '---Phosphorus Uptake        : OFF'
            write(*,*    ) '---Temperature Stress       : OFF'
            write(*,*    ) '---Grazing, Management'       
            write(*,*    ) '    and Dormancy            : OFF'
            write(*,*    ) 
        endif

    end subroutine ConstructLog

   !--------------------------------------------------------------------------


    subroutine SearchProperty(PropertyX, PropertyXIDNumber, PrintWarning, STAT)


        !Arguments-------------------------------------------------------------
        type(T_Property), optional, pointer         :: PropertyX
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
                if (PrintWarning) write (*,*)'Property Not Found in Module Vegetation: ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchProperty

    !--------------------------------------------------------------------------

    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB
        integer                                             :: KLB, KUB  
        integer                                             :: STAT_CALL
        !Begin-----------------------------------------------------------------


        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
        
        !Change this
        allocate(Me%ExternalVar%SoilNitrogenConcentration                 (ILB:IUB,JLB:JUB,KLB:KUB)) 
        allocate(Me%ExternalVar%SoilPhosphorusConcentration               (ILB:IUB,JLB:JUB,KLB:KUB)) 
        

        call SearchProperty(Me%StateVariables%TotalPlantBiomass   , TotalPlantBiomass_    , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%TotalPlantNitrogen  , TotalPlantNitrogen_   , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%TotalPlantPhosphorus, TotalPlantPhosphorus_ , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%RootBiomass         , RootBiomass_          , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%RootDepth           , RootDepth_            , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%LeafAreaIndex       , LeafAreaIndex_        , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%CanopyHeight        , CanopyHeight_         , .false., STAT = STAT_CALL)         
        call SearchProperty(Me%StateVariables%SpecificLeafStorage , SpecificLeafStorage_  , .false., STAT = STAT_CALL)        
        call SearchProperty(Me%StateVariables%EVTPCropCoefficient , EVTPCropCoefficient_  , .false., STAT = STAT_CALL)        

        allocate(Me%VegetationID                                          (ILB:IUB,JLB:JUB))
        Me%VegetationID (:,:) = FillValueInt

        !Soil Fluxes
        allocate(Me%Fluxes%WaterUptake                                    (ILB:IUB,JLB:JUB))  
        allocate(Me%Fluxes%WaterUptakeLayer                               (ILB:IUB,JLB:JUB,KLB:KUB))  
!        allocate(Me%Fluxes%FromSoil%WaterUptakeFromSoil                   (ILB:IUB,JLB:JUB))  
        Me%Fluxes%WaterUptake                                             (:,:  ) = 0.0 
        Me%Fluxes%WaterUptakeLayer                                        (:,:,:) = 0.0 
        
        allocate(Me%Growth%WaterStress                                    (ILB:IUB,JLB:JUB)) 
        Me%Growth%WaterStress                                             (:,:  ) = 1.0

        allocate(Me%ExternalVar%Integration%SumPotTP                      (ILB:IUB,JLB:JUB))
        allocate(Me%ExternalVar%Integration%AveragePotTPDuringDT          (ILB:IUB,JLB:JUB))
        Me%ExternalVar%Integration%SumPotTP         (:,:) = 0.0
        
        if (Me%ComputeOptions%Evolution%ModelSWAT) then
            allocate(Me%IsPlantGrowing                                    (ILB:IUB,JLB:JUB))
            allocate(Me%PlantingOccurred                                  (ILB:IUB,JLB:JUB))
            Me%ExternalVar%JulianDay_Old  = 0
            Me%IsPlantGrowing       (:,:) = Me%ComputeOptions%IsPlantGrowing
            
            allocate(Me%HeatUnits%PotentialHUTotal                        (ILB:IUB,JLB:JUB)) 
            allocate(Me%HeatUnits%PotentialHUBase                         (ILB:IUB,JLB:JUB)) 
            allocate(Me%HeatUnits%PlantHUAccumulated                      (ILB:IUB,JLB:JUB)) 
            allocate(Me%HeatUnits%PlantHUAccumulated_Old                  (ILB:IUB,JLB:JUB)) 
            allocate(Me%ExternalVar%Integration%AverageAirTempDuringDT    (ILB:IUB,JLB:JUB))  
            allocate(Me%ExternalVar%Integration%AverageAirHumidityDuringDT(ILB:IUB,JLB:JUB)) 
            allocate(Me%ExternalVar%Integration%AverageRadiationDuringDT  (ILB:IUB,JLB:JUB))
            allocate(Me%ExternalVar%Integration%SumTemperature            (ILB:IUB,JLB:JUB))  
            allocate(Me%ExternalVar%Integration%SumHumidity               (ILB:IUB,JLB:JUB)) 
            allocate(Me%ExternalVar%Integration%SumRadiation              (ILB:IUB,JLB:JUB))
            

            Me%ExternalVar%Integration%SumTemperature   (:,:) = 0.0
            Me%ExternalVar%Integration%SumHumidity      (:,:) = 0.0  
            Me%ExternalVar%Integration%SumRadiation     (:,:) = 0.0
            Me%HeatUnits%PotentialHUBase (:,:) = Me%ComputeOptions%PotentialHUBASE
            Me%HeatUnits%PotentialHUBase_Old = 0.0 
            Me%HeatUnits%PlantHUAccumulated (:,:) = 0.0             
            Me%HeatUnits%PlantHUAccumulated_Old (:,:) = 0.0             

            allocate(Me%Growth%GlobalStress                               (ILB:IUB,JLB:JUB))
            allocate(Me%Growth%TemperatureStress                          (ILB:IUB,JLB:JUB)) 
            allocate(Me%Growth%NitrogenStress                             (ILB:IUB,JLB:JUB)) 
            allocate(Me%Growth%PhosphorusStress                           (ILB:IUB,JLB:JUB)) 
            Me%Growth%GlobalStress                                        (:,:  ) = 1.0
            Me%Growth%TemperatureStress                                   (:,:  ) = 1.0
            Me%Growth%NitrogenStress                                      (:,:  ) = 1.0
            Me%Growth%PhosphorusStress                                    (:,:  ) = 1.0

            allocate(Me%Growth%TreeCurrentYear                            (ILB:IUB,JLB:JUB))        
            allocate(Me%Growth%TreeFractionToMaturity                     (ILB:IUB,JLB:JUB))        
            allocate(Me%Growth%TreeMaximumAnnualBiomass                   (ILB:IUB,JLB:JUB))        
            Me%Growth%TreeCurrentYear (:,:) = 0

            !Soil fluxes
            if(Me%ComputeOptions%ModelNitrogen) then
                allocate(Me%PlantNitrogenFraction                         (ILB:IUB,JLB:JUB)) 
                allocate(Me%Fluxes%NitrogenUptake                         (ILB:IUB,JLB:JUB))  
                allocate(Me%Fluxes%NitrogenUptakeLayer                    (ILB:IUB,JLB:JUB,KLB:KUB))  
!                allocate(Me%Fluxes%FromSoil%NitrogenUptakeFromSoil        (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%NitrogenUptake                                  (:,:  ) = 0.0 
                Me%Fluxes%NitrogenUptakeLayer                             (:,:,:) = 0.0 
                
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                allocate(Me%PlantPhosphorusFraction                       (ILB:IUB,JLB:JUB))
                allocate(Me%Fluxes%PhosphorusUptake                       (ILB:IUB,JLB:JUB)) 
                allocate(Me%Fluxes%PhosphorusUptakeLayer                  (ILB:IUB,JLB:JUB,KLB:KUB))  
!                allocate(Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil      (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%PhosphorusUptake                                (:,:  ) = 0.0 
                Me%Fluxes%PhosphorusUptakeLayer                           (:,:,:) = 0.0 

            endif
            
            !Aerial Fluxes
            allocate(Me%Fluxes%LAIGrowth                                  (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%LAIDeclineFraction                         (ILB:IUB,JLB:JUB)) 
            allocate(Me%LAISenescence                                     (ILB:IUB,JLB:JUB))             
            allocate(Me%PlantLAIMaxFraction                               (ILB:IUB,JLB:JUB)) 
            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                allocate(Me%LAIBeforeSenescence                           (ILB:IUB,JLB:JUB)) 
            endif



        endif

        if (Me%ComputeOptions%Grazing) then
            allocate(Me%DaysOfGrazing                                     (ILB:IUB,JLB:JUB)) 
            allocate(Me%IsPlantBeingGrazed                                (ILB:IUB,JLB:JUB)) 
            Me%DaysOfGrazing        (:,:) = 0.0 
            Me%IsPlantBeingGrazed   (:,:) = .false.
            allocate(Me%GrazingFinished                                   (ILB:IUB,JLB:JUB))
            allocate(Me%GrazingOperations                                 (ILB:IUB,JLB:JUB))
            Me%GrazingFinished      (:,:) = .false.
            Me%GrazingOperations    (:,:) = 1
            
            allocate(Me%Fluxes%BiomassGrazed                              (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%BiomassGrazedFraction                      (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%ToSoil%GrazingBiomassToSoil                (ILB:IUB,JLB:JUB)) 
            Me%Fluxes%BiomassGrazed                                       (:,:) = 0.0 
            Me%Fluxes%BiomassGrazedFraction                               (:,:) = 0.0 
            Me%Fluxes%ToSoil%GrazingBiomassToSoil                         (:,:) = 0.0 
            
            if (Me%ComputeOptions%ModelNitrogen) then
                allocate(Me%Fluxes%NitrogenGrazed                         (ILB:IUB,JLB:JUB))  
                allocate(Me%Fluxes%ToSoil%GrazingNitrogenToSoil           (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%NitrogenGrazed                                  (:,:) = 0.0
                Me%Fluxes%ToSoil%GrazingNitrogenToSoil                    (:,:) = 0.0 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                allocate(Me%Fluxes%PhosphorusGrazed                       (ILB:IUB,JLB:JUB)) 
                allocate(Me%Fluxes%ToSoil%GrazingPhosphorusToSoil         (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%PhosphorusGrazed                                (:,:) = 0.0 
                Me%Fluxes%ToSoil%GrazingPhosphorusToSoil                  (:,:) = 0.0 
            endif

        endif
        
        if (Me%ComputeOptions%Management) then
            allocate(Me%HarvestOnlyOccurred                               (ILB:IUB,JLB:JUB))
            allocate(Me%HarvestKillOccurred                               (ILB:IUB,JLB:JUB))
            allocate(Me%KillOccurred                                      (ILB:IUB,JLB:JUB))
            allocate(Me%HarvestFinished                                   (ILB:IUB,JLB:JUB))
            allocate(Me%HarvestOperations                                 (ILB:IUB,JLB:JUB))
            allocate(Me%HarvestKillOperations                             (ILB:IUB,JLB:JUB)) 
            allocate(Me%KillOperations                                    (ILB:IUB,JLB:JUB)) 
            Me%HarvestFinished      (:,:) = .false.
            Me%HarvestOperations    (:,:) = 1
            Me%HarvestKillOperations(:,:) = 1 
            Me%KillOperations       (:,:) = 1   
                      
            allocate(Me%Fluxes%BiomassRemovedInHarvest                    (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%BiomassHarvestedFraction                   (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%ToSoil%ManagementBiomassToSoil             (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%ToSoil%ManagementRootBiomassLeftInSoil     (ILB:IUB,JLB:JUB))
            Me%Fluxes%BiomassRemovedInHarvest                             (:,:) = 0.0  
            Me%Fluxes%BiomassHarvestedFraction                            (:,:) = 0.0
            Me%Fluxes%ToSoil%ManagementBiomassToSoil                      (:,:) = 0.0 
            Me%Fluxes%ToSoil%ManagementRootBiomassLeftInSoil              (:,:) = 0.0
            
            if (Me%ComputeOptions%ModelNitrogen) then            
                allocate(Me%Fluxes%NitrogenRemovedInHarvest               (ILB:IUB,JLB:JUB))  
                allocate(Me%Fluxes%ToSoil%ManagementNitrogenToSoil        (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%NitrogenRemovedInHarvest                        (:,:) = 0.0 
                Me%Fluxes%ToSoil%ManagementNitrogenToSoil                 (:,:) = 0.0 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                allocate(Me%Fluxes%PhosphorusRemovedInHarvest             (ILB:IUB,JLB:JUB)) 
                allocate(Me%Fluxes%ToSoil%ManagementPhosphorusToSoil      (ILB:IUB,JLB:JUB))
                Me%Fluxes%PhosphorusRemovedInHarvest                      (:,:) = 0.0
                Me%Fluxes%ToSoil%ManagementPhosphorusToSoil               (:,:) = 0.0 
            endif

        endif

        if (Me%ComputeOptions%Dormancy) then
            allocate(Me%IsPlantDormant                                    (ILB:IUB,JLB:JUB)) 
            allocate(Me%PlantGoingDormant                                 (ILB:IUB,JLB:JUB))
            Me%IsPlantDormant       (:,:) = .false.
            
            allocate(Me%Fluxes%BiomassRemovedInDormancy                   (ILB:IUB,JLB:JUB)) 
            allocate(Me%Fluxes%ToSoil%DormancyBiomassToSoil               (ILB:IUB,JLB:JUB)) 
            Me%Fluxes%BiomassRemovedInDormancy                            (:,:) = 0.0
            Me%Fluxes%ToSoil%DormancyBiomassToSoil                        (:,:) = 0.0
            
            if (Me%ComputeOptions%ModelNitrogen) then
                allocate(Me%Fluxes%NitrogenRemovedInDormancy              (ILB:IUB,JLB:JUB))  
                allocate(Me%Fluxes%ToSoil%DormancyNitrogenToSoil          (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%NitrogenRemovedInDormancy                       (:,:) = 0.0
                Me%Fluxes%ToSoil%DormancyNitrogenToSoil                   (:,:) = 0.0
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                allocate(Me%Fluxes%PhosphorusRemovedInDormancy            (ILB:IUB,JLB:JUB)) 
                allocate(Me%Fluxes%ToSoil%DormancyPhosphorusToSoil        (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%PhosphorusRemovedInDormancy                     (:,:) = 0.0
                Me%Fluxes%ToSoil%DormancyPhosphorusToSoil                 (:,:) = 0.0
            endif
 
       endif
        
        if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH1 (ILB:IUB,JLB:JUB))
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH2 (ILB:IUB,JLB:JUB))
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3 (ILB:IUB,JLB:JUB))
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH4 (ILB:IUB,JLB:JUB))
        endif
        


        if (Me%ComputeOptions%ModelPlantBiomass) then
            allocate(Me%Fluxes%BiomassGrowth                              (ILB:IUB,JLB:JUB)) 
            allocate(Me%Growth%BiomassGrowthOld                           (ILB:IUB,JLB:JUB))  
            
            Me%Fluxes%BiomassGrowth (:,:) = 0.0 
            Me%Growth%BiomassGrowthOld(:,:) = 0.0 
        endif

        if (Me%ComputeOptions%ModelCanopyHeight) then
            allocate(Me%ChangeCanopyEnabled                               (ILB:IUB,JLB:JUB)) 
            Me%ChangeCanopyEnabled  (:,:) = .false.
        endif
        

    end subroutine AllocateVariables

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: nProperties, STAT_CALL
        integer                                             :: iflag, AddProperties
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !First checks out how many properties will have time series
        PropertyX   => Me%FirstProperty
        nProperties =  0
        do while (associated(PropertyX))
            
            if (PropertyX%TimeSerie) then 
            
                nProperties = nProperties + 1
                
            end if
            
            PropertyX=>PropertyX%Next
        enddo


        if (nProperties > 0) then
            
            Me%OutPut%TimeSerie_ON = .true.
        
            !Allocates PropertyList
            AddProperties = 0
            !Always water uptake and water stress
            AddProperties = AddProperties + 2
            if (Me%ComputeOptions%Evolution%ModelSWAT) then
                AddProperties = AddProperties + 5
                if (Me%ComputeOptions%ModelNitrogen) then
                    AddProperties = AddProperties + 1
                endif  
                if (Me%ComputeOptions%ModelPhosphorus) then
                    AddProperties = AddProperties + 1
                endif 
            endif
            
            !Debug average computation
!            AddProperties = AddProperties + 1
            
            allocate(PropertyList(1:nProperties + AddProperties), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructTimeSerie - ModuleVegetation - ERR01'

            !Fills up PropertyList
            PropertyX   => Me%FirstProperty
            nProperties =  0
            do while (associated(PropertyX))
                if (PropertyX%TimeSerie) then
                    nProperties = nProperties + 1
                    PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name))
                endif
                PropertyX=>PropertyX%Next
            enddo
            
            nProperties = nProperties + 1
            PropertyList(nProperties) = "Water Uptake m3/s"
            nProperties = nProperties + 1
            PropertyList(nProperties) = "WaterStressFactor"
            
            if (Me%ComputeOptions%Evolution%ModelSWAT) then
                nProperties = nProperties + 1 
                PropertyList(nProperties) = "HU Accumulated"
                nProperties = nProperties + 1
                PropertyList(nProperties) = "Potential HU"
                nProperties = nProperties + 1
                PropertyList(nProperties) = "NitrogenStressFactor"  
                nProperties = nProperties + 1
                PropertyList(nProperties) = "PhosphorusStressFactor"          
                nProperties = nProperties + 1
                PropertyList(nProperties) = "TemperatureStressFactor"
                if (Me%ComputeOptions%ModelNitrogen) then
                    nProperties = nProperties + 1
                    PropertyList(nProperties) = "Nitrogen Uptake kg/ha"            
                endif
                if (Me%ComputeOptions%ModelPhosphorus) then
                    nProperties = nProperties + 1
                    PropertyList(nProperties) = "Phosphorus Uptake kg/ha"
                endif
            endif
 
!            nProperties = nProperties + 1
!            PropertyList(nProperties) = "Average Radiation W/m2"

            call GetData(TimeSerieLocationFile,                                 &
                         Me%ObjEnterData,iflag,                                 &
                         SearchType   = FromFile,                               &
                         keyword      = 'TIME_SERIE_LOCATION',                  &
                         ClientModule = 'ModuleVegetation',                     &
                         Default      = Me%Files%ConstructData,                 &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                        &
                stop 'Construct_Time_Serie - ModuleVegetation - ERR02' 


            !Constructs TimeSerie
            call StartTimeSerie(Me%ObjTimeSerie,                                &
                                Me%ObjTime,                                     &
                                TimeSerieLocationFile,                          &
                                PropertyList, "svg",                            &
                                WaterPoints2D = Me%ExternalVar%MappingPoints2D, &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_Time_Serie - ModuleVegetation - ERR03'

            !Deallocates PropertyList
            deallocate(PropertyList, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_Time_Serie - ModuleVegetation - ERR04'

        endif
       
    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------

    subroutine ConstructHDF
        
        !External-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        logical                                     :: OutputON
        integer :: STAT_CALL

        !Begin-----------------------------------------------------------------

        nullify(Me%OutPut%OutTime)

        OutputON        = OFF

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))
            
            if(CurrentProperty%OutputHDF) OutputON = ON
            CurrentProperty => CurrentProperty%Next

        enddo

        if(OutputON)then
  
            call GetOutPutTime(Me%ObjEnterData,                              &
                               CurrentTime = Me%BeginTime,                   &
                               EndTime     = Me%EndTime,                     &
                               keyword     = 'OUTPUT_TIME',                  &
                               SearchType  = FromFile,                       &
                               OutPutsTime = Me%OutPut%OutTime,              &
                               OutPutsOn   = Me%OutPut%HDF_ON,               &
                               OutPutsNumber = Me%OutPut%Number,             &
                               STAT        = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                       &
                stop 'ConstructHDF - ModuleVegetation - ERR01' 

            if (Me%OutPut%HDF_ON) then

                Me%OutPut%NextOutPut = 1

                call Open_HDF5_OutPut_File

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructHDF - ModuleVegetation - ERR02'
            endif 

        endif

    end subroutine ConstructHDF

    !--------------------------------------------------------------------------
    
    subroutine Open_HDF5_OutPut_File

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


        Me%OutPut%NextOutPut = 1  
        
        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR01'
        
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%Results)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR02'

      
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR03'


        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR05'
        
        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR06'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExternalVar%MappingPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR07'

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR08'       

        !UnGets Topography
        call UnGetGridData      (Me%ObjGridData, Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR09'


    end subroutine Open_HDF5_OutPut_File

   !--------------------------------------------------------------------------

    subroutine ConstructVegetationGrid

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ObjGD, i, j, ivt, k
        real, dimension(:,:), pointer           :: AuxID
        real                                    :: MaxRootDepth
        
        !Begin------------------------------------------------------------------


        !Gets Vegetation IDs
        ObjGD = 0
        call ConstructGridData  (ObjGD, Me%ObjHorizontalGrid, FileName = Me%Files%VegetationIDFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR03'
        
        !allocate (AuxID(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        call GetGridData (ObjGD, AuxID, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR04'

        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  SZZ         = Me%ExternalVar%SZZ,                     &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR05'
                                                     

        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                do ivt = 1, Me%VegetationsNumber
                    if (Me%VegetationTypes(ivt)%ID == NINT(AuxID(i, j))) then
                        !Me%VegetationID(i, j) = NINT(AuxID(i, j))
                        Me%VegetationID(i, j) = ivt
                    endif
                    if (ivt .eq. Me%VegetationsNumber .and. Me%VegetationID(i,j) < FillValueInt / 2) then
                        write(*,*)'Vegetation ID not found in vegetation type definition.'
                        write(*,*)'Check in vegetation grid or in vegetation types in data file, the ID:', NINT(AuxID(i,j))
                        write(*,*)'or Name:', trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ConstructVegetationGrid - ModuleVegetation - ERR06'
                    endif
                enddo

                !Check maximum root depth
                MaxRootDepth = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaximumRootDepth
                k = Me%WorkSize%KUB
                if(MaxRootDepth .gt. -1. * (Me%ExternalVar%SZZ(i,j,k))) then
                    Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaximumRootDepth = -1. * (Me%ExternalVar%SZZ(i,j,k))
                    write(*,*)'Maximum Root Depth is greater then soil depth. Maximum root depth was set to soil depth'
                    write(*,*)'Cell:', i, j
                endif
            endif
        enddo
        enddo
        
        call UnGetGeometry( Me%ObjGeometry, Me%ExternalVar%SZZ,          STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR07'


        !deallocate (AuxID)
        
        call UnGetGridData (ObjGD, AuxID, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR08'
       
        
        call KillGridData (ObjGD, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR09'

    end subroutine ConstructVegetationGrid

   !--------------------------------------------------------------------------

    subroutine ConstructFeddes

        !Local----------------------------------------------------------------
        integer                                      :: i,j
        !Begin----------------------------------------------------------------   
    

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH1(i,j) = &
                Me%VegetationTypes(Me%VegetationID(i, j))%RootFeddesH1 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH2(i,j) = &
                Me%VegetationTypes(Me%VegetationID(i, j))%RootFeddesH2 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH3(i,j) = &
                Me%VegetationTypes(Me%VegetationID(i, j))%RootFeddesH3 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH4(i,j) = &
                Me%VegetationTypes(Me%VegetationID(i, j))%RootFeddesH4 

            endif
        enddo
        enddo
    
   
    end subroutine ConstructFeddes

   !--------------------------------------------------------------------------
    
    subroutine ConstructVegetationParameters

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: VegetationTypeFound
        logical                                     :: BlockFound
        integer                                     :: ClientNumber, iflag
        integer                                     :: ivt
        !----------------------------------------------------------------------

        !Counts the number of VegetationTypes
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        Me%VegetationsNumber = 0
doS:    do
            !Gets Vegetations type
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginvegetationtype>',              &
                                        block_end       = '<endvegetationtype>',                &
                                        BlockFound      = VegetationTypeFound,                  &   
                                        STAT            = STAT_CALL)
SF:         if (STAT_CALL == SUCCESS_ .and. VegetationTypeFound) then
                Me%VegetationsNumber = Me%VegetationsNumber + 1
            else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR02'
               
                exit doS

            end if SF

            !old version of vegetation
            call ExtractBlockFromBlock(Me%ObjEnterData,                                         &
                                        ClientNumber        = ClientNumber,                     &
                                        block_begin         = '<<begin_property>>',             &
                                        block_end           = '<<end_property>>',               &
                                        BlockInBlockFound   = BlockFound,                       &   
                                        STAT                = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_ .and. BlockFound) then  
                write(*,*)  
                write(*,*) 'Using old formulation on vegetation' 
                write(*,*) 'Vegetation properties are now defined as MOHID'
                write(*,*) 'fill matrix standard with <beginproperty>, <endproperty> block.'
                write(*,*) 'For support verify options at MOHID'
                write(*,*) 'wiki in http://www.mohid.com/wiki'
                stop 'ConstructVegetationParameters - ModuleVegetation - ERR03'

            endif                   
        enddo doS            
        
        allocate(Me%VegetationTypes(Me%VegetationsNumber))
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)

        ivt = 0

doH:    do 
            !Gets Vegetations type
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<beginvegetationtype>',              &
                                        block_end       = '<endvegetationtype>',                &
                                        BlockFound      = VegetationTypeFound,                  &   
                                        STAT            = STAT_CALL)
HF:         if (STAT_CALL == SUCCESS_ .and. VegetationTypeFound) then

                ivt = ivt + 1
                
                !Reads ID
                call GetData(Me%VegetationTypes(ivt)%ID, Me%ObjEnterData,  iflag,               &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'ID',                                             &
                             ClientModule   = 'ModuleVegetation',                               &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR10'
                if (iflag /= 1) then
                    write(*,*)'Missing ID in Vegetation Type definition'
                    stop 'ConstructVegetationParameters - ModuleVegetation - ERR15'
                endif

                !Reads Name
                call GetData(Me%VegetationTypes(ivt)%Name, Me%ObjEnterData,  iflag,             &
                             SearchType     = FromBlock,                                        &
                             keyword        = 'NAME',                                           &
                             ClientModule   = 'ModuleVegetation',                               &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR20'
                if (iflag /= 1) then
                    write(*,*)'Missing NAME in Vegetation Type definition'
                    stop 'ConstructVegetationParameters - ModuleVegetation - ERR30'
                endif
                
                if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

                    call GetData(Me%VegetationTypes(ivt)%HasLeaves, Me%ObjEnterData,  iflag,        &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'HAS_LEAVES',                                     &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR31'
                    if (iflag /= 1) then
                        write(*,*)'Missing HAS_LEAVES in Vegetation Type definition'
                        stop 'ConstructVegetationParameters - ModuleVegetation - ERR32'
                    endif
                endif
                
                if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then
                    !Reads Feddes
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH1, Me%ObjEnterData,  iflag,     &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H1',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR40'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H1 in Vegetation Type definition'
                        stop 'ConstructVegetationParameters - ModuleVegetation - ERR50'
                    endif
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH2, Me%ObjEnterData,  iflag,     &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H2',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR60'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H2 in Vegetation Type definition'
                        stop 'ConstructVegetationParameters - ModuleVegetation - ERR70'
                    endif
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH3, Me%ObjEnterData,  iflag,     &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H3',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR80'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H3 in Vegetation Type definition'
                        stop 'ConstructVegetationParameters - ModuleVegetation - ERR90'
                    endif
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH4, Me%ObjEnterData,  iflag,     &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H4',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR100'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H4 in Vegetation Type definition'
                        stop 'ConstructVegetationParameters - ModuleVegetation - ERR110'
                    endif
  
                endif

                if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                    call ReadTimingDatabase(ivt, ClientNumber)
                    call ReadGrowthDatabase(ivt, ClientNumber)
                    if (Me%ComputeOptions%Management .or. Me%ComputeOptions%Grazing) then 
                        call ReadManagementAndGrazeDatabase(ivt, ClientNumber)
                    endif
                endif
            else

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR120'
               
                exit doH

            end if HF
        
        enddo doH

   
    end subroutine ConstructVegetationParameters

   !--------------------------------------------------------------------------
    
    subroutine ReadTimingDatabase(ivt, ClientNumber)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound
        integer                                     :: iflag, FailRead

        !Begin-----------------------------------------------------------------



        call ExtractBlockFromBlock(Me%ObjEnterData,                                           &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<begintimingdatabase>',              &
                                    block_end         = '<endtimingdatabase>',                &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

            !Reads timing Parameters

            call GetData(Me%VegetationTypes(ivt)%TimingDatabase%PlantHUAtMaturity, Me%ObjEnterData,  iflag,               &
                         SearchType     = FromBlockInBlock,                                                               &
                         Default        = -99.,                                                                           &
                         keyword        = 'MATURITY_HU',                                                                  &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadTimingDatabase - ModuleVegetation - ERR10'
            if (iflag /= 1) then
                write(*,*)'Missing HU to maturity in timing parameters definition -', trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadTimingDatabase - ModuleVegetation - ERR15'
            endif
                    
            call GetData(Me%VegetationTypes(ivt)%TimingDatabase%PlantingJulianDay, Me%ObjEnterData,  iflag,               &
                         SearchType     = FromBlockInBlock,                                                               &
                         Default        = -99.,                                                                           &
                         keyword        = 'PLANTING_JULIANDAY',                                                           &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadTimingDatabase - ModuleVegetation - ERR20'       
            if (iflag /= 1) then
                FailRead = 0
                FailRead = FailRead +1
            endif
                        
            call GetData(Me%VegetationTypes(ivt)%TimingDatabase%PlantingHUBase, Me%ObjEnterData,  iflag,                  &
                         SearchType     = FromBlockInBlock,                                                               &
                         Default        = -99.,                                                                           &
                         keyword        = 'PLANTING_HUBASE',                                                              &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadTimingDatabase - ModuleVegetation - ERR30'        
            if (iflag /= 1) then
                FailRead = FailRead + 1
                if (FailRead == 2) then
                    write(*,*)'Missing HU or julian day in planting timing parameters definition -', trim(Me%VegetationTypes(ivt)%Name)
                    stop 'ReadTimingDatabase - ModuleVegetation - ERR35'
                endif
            endif

        else

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadTimingDatabase - ModuleVegetation - ERR40'
           
        endif HF

    end subroutine ReadTimingDatabase

    !--------------------------------------------------------------------------
    
    subroutine ReadGrowthDatabase(ivt, ClientNumber)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound
        integer                                     :: iflag, PlantType

        !Begin-----------------------------------------------------------------



        call ExtractBlockFromBlock(Me%ObjEnterData,                                           &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<begingrowthdatabase>',              &
                                    block_end         = '<endgrowthdatabase>',                &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

            !Reads Growth Parameters
            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantType, Me%ObjEnterData,  iflag,                       &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'PLANT_TYPE',                                                                   &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR10'
            if (iflag /= 1) then
                write(*,*)'Missing plant type in growth parameters definition -', trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR15'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionN1, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_NITROGENFRACTION_N1',                                                  &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR20'
            if (iflag /= 1) then
                write(*,*)'Missing optimal nitrogen fraction parameter in growth parameters definition -',                &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR25'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionN2, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_NITROGENFRACTION_N2',                                                  &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR30'
            if (iflag /= 1) then
                write(*,*)'Missing optimal nitrogen fraction parameter in growth parameters definition -',                &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR35'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionN3, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               & 
                         keyword        = 'OPTIMAL_NITROGENFRACTION_N3',                                                  &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR40'
            if (iflag /= 1) then
                write(*,*)'Missing optimal nitrogen fraction parameter in growth parameters definition -',                &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR45'
            endif
            
            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionP1, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_PHOSPHORUSFRACTION_P1',                                                &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR50'
            if (iflag /= 1) then
                write(*,*)'Missing optimal phosphorus fraction parameter in growth parameters definition -',              &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR55'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionP2, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_PHOSPHORUSFRACTION_P2',                                                &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR60'
            if (iflag /= 1) then
                write(*,*)'Missing optimal phosphorus fraction parameter in growth parameters definition -',              &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR65'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionP3, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_PHOSPHORUSFRACTION_P3',                                                &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR70'
            if (iflag /= 1) then
                write(*,*)'Missing optimal phosphorus fraction parameter in growth parameters definition -',              &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR75'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantBaseTemperature, Me%ObjEnterData,  iflag,            &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'BASE_TEMPERATURE',                                                             &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR80'
            if (iflag /= 1) then
                write(*,*)'Missing base temperature in growth parameters definition -',                                   &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR85'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantOptimalTemperature, Me%ObjEnterData,  iflag,         &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_TEMPERATURE',                                                          &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR90'
            if (iflag /= 1) then
                write(*,*)'Missing optimal temperature in growth parameters definition -',                                &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR95'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%ExtinctCoef, Me%ObjEnterData,  iflag,                     &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'RADIATION_EXTINCTION_COEF',                                                    &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR100'
            if (iflag /= 1) then
                write(*,*)'Missing radiation extiction coefficient in growth parameters definition -',                    &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR105'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%BiomassEnergyRatio, Me%ObjEnterData,  iflag,              &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'BIOMASS_ENERGY_RATIO',                                                         &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR110'
            if (iflag /= 1) then
                write(*,*)'Missing biomass energy ratio in growth parameters definition -',                               &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR115'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%CO2ConcHigh, Me%ObjEnterData,  iflag,                     &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'CO2_HIGH',                                                                     &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR120'
            if (iflag /= 1) then
                write(*,*)'Missing high CO2 concentration in growth parameters definition -',                             &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR125'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%BiomassEnergyRatioHigh, Me%ObjEnterData,  iflag,          &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'BIOMASS_ENERGY_RATIO_HIGH',                                                    &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR130'
            if (iflag /= 1) then
                write(*,*)'Missing high biomass energy ratio in growth parameters definition -',                          &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR135'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%RUEDeclineRate, Me%ObjEnterData,  iflag,                  &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'RUE_DECLINE_RATE',                                                             &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR140'
            if (iflag /= 1) then
                write(*,*)'Missing RUE decline rate in growth parameters definition -',                                   &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR145'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrLAIMax1, Me%ObjEnterData,  iflag,                       &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_LAIMAXFRACTION_1',                                                     &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR150'
            if (iflag /= 1) then
                write(*,*)'Missing optimal fraction of maximum LAI in growth parameters definition -',                    &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR155'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrLAIMax2, Me%ObjEnterData,  iflag,                       &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_LAIMAXFRACTION_2',                                                     &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR160'
            if (iflag /= 1) then
                write(*,*)'Missing optimal fraction of maximum LAI in growth parameters definition -',                    &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR165'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrGrow1, Me%ObjEnterData,  iflag,                         &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'GROWFRACTION_1',                                                               &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR170'
            if (iflag /= 1) then
                write(*,*)'Missing grow stage fraction asscoiated with LAI curve in growth parameters definition -',      &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR175'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrGrow2, Me%ObjEnterData,  iflag,                         &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'GROWFRACTION_2',                                                               &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR180'
            if (iflag /= 1) then
                write(*,*)'Missing grow stage fraction asscoiated with LAI curve in growth parameters definition -',      &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR185'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrGrowLAIDecline, Me%ObjEnterData,  iflag,                &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'GROWFRACTION_LAIDECLINE',                                                      &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR190'
            if (iflag /= 1) then
                write(*,*)'Missing grow stage fraction asscoiated with LAI curve in growth parameters definition -',      &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR195'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%LAIMax, Me%ObjEnterData,  iflag,                          &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'LAI_MAX',                                                                      &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR200'
            if (iflag /= 1) then
                write(*,*)'Missing maximum LAI in growth parameters definition -',                                        &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR205'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%MaximumRootDepth, Me%ObjEnterData,  iflag,                &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'ROOT_DEPTH_MAX',                                                               &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR210'
            if (iflag /= 1) then
                write(*,*)'Missing maximum root depth in growth parameters definition -',                                 &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR215'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%MaxCanopyHeight, Me%ObjEnterData,  iflag,                 &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'CANOPY_HEIGHT_MAX',                                                            &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR220'
            if (iflag /= 1) then
                write(*,*)'Missing maximum canopy height in growth parameters definition -',                              &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR225'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%OptimalHarvestIndex, Me%ObjEnterData,  iflag,             &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'OPTIMAL_HARVEST_INDEX',                                                        &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR230'
            if (iflag /= 1) then
                write(*,*)'Missing optimal harvest index in growth parameters definition -',                              &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR235'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%MinimumHarvestIndex, Me%ObjEnterData,  iflag,             &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'MINIMUM_HARVEST_INDEX',                                                        &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR240'
            if (iflag /= 1) then
                write(*,*)'Missing minimum harvest index in growth parameters definition -',                              &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR245'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%NitrogenFractionInYeld, Me%ObjEnterData,  iflag,          &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'YELD_NITROGENFRACTION',                                                        &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR250'
            if (iflag /= 1) then
                write(*,*)'Missing nitrogen fraction in yeld in growth parameters definition -',                          &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR255'
            endif

            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PhosphorusFractionInYeld, Me%ObjEnterData,  iflag,        &
                         SearchType     = FromBlockInBlock,                                                               &
                         keyword        = 'YELD_PHOSPHORUSFRACTION',                                                      &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR260'
            if (iflag /= 1) then
                write(*,*)'Missing phosphorus fraction in yeld in growth parameters definition -',                        &
                           trim(Me%VegetationTypes(ivt)%Name)
                stop 'ReadGrowthDatabase - ModuleVegetation - ERR265'
            endif

            PlantType = Me%VegetationTypes(ivt)%GrowthDatabase%PlantType
            if (PlantType == 7) then
            
                call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%TreeYearsToMaturity, Me%ObjEnterData,  iflag,         &
                             SearchType     = FromBlockInBlock,                                                           &
                             keyword        = 'TREE_YEARSTOMATURITY',                                                     &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR270'
                if (iflag /= 1) then
                    write(*,*)'Missing tree years too maturity in growth parameters definition -',                        &
                               trim(Me%VegetationTypes(ivt)%Name)
                    stop 'ReadGrowthDatabase - ModuleVegetation - ERR275'
                endif

                call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%TreeMaximumBiomass, Me%ObjEnterData,  iflag,          &
                             SearchType     = FromBlockInBlock,                                                           &
                             keyword        = 'TREE_MAXIMUMBIOMASS',                                                      &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR280'
                if (iflag /= 1) then
                    write(*,*)'Missing tree maximum biomass in growth parameters definition -',                           &
                               trim(Me%VegetationTypes(ivt)%Name)
                    stop 'ReadGrowthDatabase - ModuleVegetation - ERR285'
                endif
            
            endif

        else

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR290'
           
        endif HF

    end subroutine ReadGrowthDatabase

    !--------------------------------------------------------------------------


    subroutine ReadManagementAndGrazeDatabase (ivt, ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock(Me%ObjEnterData,                                           &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<beginmanagementandgrazedatabase>',  &
                                    block_end         = '<endmanagementandgrazedatabase>',    &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

            !Reads Parameters
            
            if (Me%ComputeOptions%Grazing) then

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%GrazingStartJulianDay, Me%ObjEnterData,  iflag,   &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'GRAZING_START_JULIANDAY',                                                  &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR40'     
                
                allocate (Me%VegetationTypes(ivt)%ManagementDatabase%GrazingStartPlantHU(1))
                
                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%GrazingStartPlantHU(1), Me%ObjEnterData,  iflag,  &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'GRAZING_START_PLANTHU',                                                    &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR50'     

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%GrazingDays, Me%ObjEnterData,  iflag,             &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = 0,                                                                          &
                             keyword        = 'GRAZING_DAYS',                                                             &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR60'     
                
                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%GrazingMinimumBiomass, Me%ObjEnterData,  iflag,   &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = 0.,                                                                         &
                             keyword        = 'MINIMUM_BIOMASS_FOR_GRAZING',                                              &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR70'      

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%GrazingBiomass, Me%ObjEnterData,  iflag,          &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = 0.,                                                                         &
                             keyword        = 'GRAZING_BIOMASS',                                                          &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR80'                               

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%TramplingBiomass, Me%ObjEnterData,  iflag,        &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = 0.,                                                                         &
                             keyword        = 'TRAMPLING_BIOMASS',                                                        &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR90'     

            endif

            if (Me%ComputeOptions%Management) then
                
                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%HarvestKillJulianDay, Me%ObjEnterData,  iflag,    &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'HARVESTKILL_JULIANDAY',                                                    &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR100'     

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%HarvestKillPlantHU, Me%ObjEnterData,  iflag,      &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'HARVESTKILL_PLANTHU',                                                      &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR110'     

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%HarvestJulianDay, Me%ObjEnterData,  iflag,        &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'HARVEST_JULIANDAY',                                                        &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR120'     


                allocate (Me%VegetationTypes(ivt)%ManagementDatabase%HarvestPlantHU(1))

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%HarvestPlantHU(1), Me%ObjEnterData,  iflag,       &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'HARVEST_PLANTHU',                                                          &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR130'     

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%KillJulianDay, Me%ObjEnterData,  iflag,           &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'KILL_JULIANDAY',                                                           &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR140'     

                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%KillPlantHU, Me%ObjEnterData,  iflag,             &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = -99.,                                                                       &
                             keyword        = 'KILL_PLANTHU',                                                             &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR150'     
                
                call GetData(Me%VegetationTypes(ivt)%ManagementDatabase%HarvestEfficiency, Me%ObjEnterData,  iflag,       &
                             SearchType     = FromBlockInBlock,                                                           &
                             Default        = 1.,                                                                         &
                             keyword        = 'HARVEST_EFFICIENCY',                                                       &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadManagementAndGrazeDatabase - ModuleVegetation - ERR160'     
            endif
        else

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR170'
           
        endif HF


    end subroutine ReadManagementAndGrazeDatabase


    !--------------------------------------------------------------------------        
   
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetLeafAreaIndex(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer, optional     :: Scalar
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mVEGETATION_, Me%InstanceID)

            Scalar  => Me%StateVariables%LeafAreaIndex%Field

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetLeafAreaIndex
    
    !--------------------------------------------------------------------------

    subroutine GetSpecificLeafStorage(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer, optional     :: Scalar
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mVEGETATION_, Me%InstanceID)

            Scalar  => Me%StateVariables%SpecificLeafStorage%Field

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetSpecificLeafStorage
        
    !--------------------------------------------------------------------------

    subroutine GetEVTPCropCoefficient(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer, optional     :: Scalar
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mVEGETATION_, Me%InstanceID)

            Scalar  => Me%StateVariables%EVTPCropCoefficient%Field

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetEVTPCropCoefficient 
           
    !--------------------------------------------------------------------------

    subroutine GetRootDepth(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer, optional     :: Scalar
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mVEGETATION_, Me%InstanceID)

            Scalar  => Me%StateVariables%RootDepth%Field

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetRootDepth
    
    !--------------------------------------------------------------------------

    subroutine GetTranspiration(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer, optional     :: Scalar
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mVEGETATION_, Me%InstanceID)

            Scalar  => Me%Fluxes%WaterUptake

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine GetTranspiration

    !--------------------------------------------------------------------------

    subroutine UngetVegetation2D(VegetationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: VegetationID
        real, pointer, dimension(:,:)   :: Array
        integer, optional, intent (OUT) :: STAT
   
        !External--------------------------------------------------------------
        integer                         :: ready_   

        !Local-----------------------------------------------------------------
        integer                         :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mVegetation_, Me%InstanceID, "UngetVegetation2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetVegetation2D    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyVegetation(ObjVegetationID,                  &
                                MappingPoints,                    &
                                PotentialTranspiration,           &
                                ActualTranspiration,              &
                                STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                            :: ObjVegetationID         !IN
        integer, dimension(:,:  ), pointer             :: MappingPoints           !IN
        real,    dimension(:,:  ), pointer             :: PotentialTranspiration  !IN
        real,    dimension(:,:,:), pointer             :: ActualTranspiration     !OUT
        integer, intent(OUT),      optional            :: STAT                    !OUT

        !Local-----------------------------------------------------------------
        integer                                        :: STAT_, ready_
        integer                                        :: STAT_CALL
        integer                                        :: JulDay
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjVegetationID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Actualize the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyVegetation - ModuleVegetation - ERR01'

            !Sets External Variable
            Me%ExternalVar%MappingPoints2D => MappingPoints
            Me%Externalvar%PotentialTranspiration => PotentialTranspiration
        
            ! Transpiration is averaged during Integration DT until vegetation processes are called
            ! If vegetation growth model is used also atmosphere properties are averaged
            if (Me%ExternalVar%Now >= Me%NextIntegration) then     

                call AverageExtPropDuringDT
                Me%NextIntegration =  Me%NextIntegration + Me%ComputeOptions%IntegrationDT

            endif


            !Time when vegetation processes are called
            if (Me%ExternalVar%Now >= Me%NextCompute) then     
                

                call ReadLockExternalVar(ReadAtmosphere = .false.)

            
                if (Me%ComputeOptions%Evolution%ReadNeeded) then
                
                    call ModifyReadedProperties

                endif            
            
                if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                
                    call NullifyFluxes
                    
                    !Check growing, dormancy, grazing, counters...
                    call CheckPlantState

                    !Compute fluxes needed to update properties
                    call ModifyFluxes

                    !Compute new property values
                    call ModifyModelledProperties
                
                    !Compute nutrient applyed to soil in fertilization and organic resiude in debris.
                    call InterfaceWithSoil

                else 
                    
                    !Do not use vegetation growth model. Only water uptake is modeled
                    !Use readed vegetation properties (LAI, root depth) to uptake water
                    call WaterUptake
                
                endif
            

                !Output
                if (Me%OutPut%HDF_ON)          call Modify_OutputHDF                        
                if (Me%OutPut%TimeSerie_ON)     call Modify_OutPutTimeSeries


                !Get julian day from iteration to keep to next
                call JulianDay(Me%ExternalVar%Now, JulDay)
                Me%ExternalVar%JulianDay_Old = JulDay

                nullify (Me%ExternalVar%MappingPoints2D)
                
                
                call ReadUnLockExternalVar(ReadAtmosphere = .false.)

                Me%NextCompute =  Me%NextCompute + Me%ComputeOptions%VegetationDT

            
            endif

            !Outgoing Variables
            ActualTranspiration => Me%Fluxes%WaterUptakeLayer

           
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyVegetation
    
    !--------------------------------------------------------------------------

    subroutine AverageExtPropDuringDT
        
        !Local-----------------------------------------------------------------
        integer                                            :: i,j
        integer, dimension(:,:), pointer                   :: MappingPoints
        !Begin-----------------------------------------------------------------

        !If dt's the same use values without integration
        if (Me%ComputeOptions%VegetationDT .eq. Me%ComputeOptions%IntegrationDT) then

            Me%ExternalVar%Integration%AveragePotTPDuringDT => Me%ExternalVar%PotentialTranspiration


            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

                call ReadLockExternalVar(ReadAtmosphere = .true.)

                Me%ExternalVar%Integration%AverageAirTemPDuringDT     => Me%ExternalVar%AirTemperature
    
                Me%ExternalVar%Integration%AverageAirHumidityDuringDT => Me%ExternalVar%RelativeHumidity

                Me%ExternalVar%Integration%AverageRadiationDuringDT   => Me%ExternalVar%SolarRadiation
    
                call ReadUnLockExternalVar(ReadAtmosphere = .true.)

            endif
        

        else
        
            MappingPoints => Me%ExternalVar%MappingPoints2D
            call ReadLockExternalVar(ReadAtmosphere = .true.)

            Me%nIterations = Me%nIterations + 1
      
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (MappingPoints (i, j) == 1) then  
                

                    ! Atmosphere properties averaged if needed
                    if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

               
                        !Temperature
                        Me%ExternalVar%Integration%SumTemperature(i,j)   = Me%ExternalVar%Integration%SumTemperature(i,j)       &
                                                                           + Me%ExternalVar%AirTemperature(i,j)
                        !Humidity
                        Me%ExternalVar%Integration%SumHumidity(i,j)      = Me%ExternalVar%Integration%SumHumidity(i,j)          &
                                                                           + Me%ExternalVar%RelativeHumidity(i,j)
                        !Radiation
                        Me%ExternalVar%Integration%SumRadiation(i,j)     = Me%ExternalVar%Integration%SumRadiation(i,j)         &
                                                                           + Me%ExternalVar%SolarRadiation(i,j)
                    
               
                    
                        if (Me%ExternalVar%Now .ge. Me%NextCompute) then
                            Me%ExternalVar%Integration%AverageAirTemPDuringDT(i,j)      = Me%ExternalVar%Integration%SumTemperature(i,j) &
                                                                                          / Me%nIterations
                            Me%ExternalVar%Integration%AverageAirHumidityDuringDT(i,j)  = Me%ExternalVar%Integration%SumHumidity(i,j)    &
                                                                                          / Me%nIterations
                            Me%ExternalVar%Integration%AverageRadiationDuringDT(i,j)    = Me%ExternalVar%Integration%SumRadiation(i,j)   &
                                                                                          / Me%nIterations
                            Me%ExternalVar%Integration%SumTemperature(i,j) = 0.0
                            Me%ExternalVar%Integration%SumHumidity   (i,j) = 0.0
                            Me%ExternalVar%Integration%SumRadiation  (i,j) = 0.0                                                           
                        endif
                    endif

                    ! Transpiration always needed
                    Me%ExternalVar%Integration%SumPotTP(i,j) = Me%ExternalVar%Integration%SumPotTP(i,j)                                  &
                                                               + Me%ExternalVar%PotentialTranspiration(i,j)
                
                    if (Me%ExternalVar%Now .ge. Me%NextCompute) then
                        Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j)  = Me%ExternalVar%Integration%SumPotTP(i,j)                 &
                                                                                 / Me%nIterations      
                        Me%ExternalVar%Integration%SumPotTP(i,j) = 0.0
                    
                                                                             
                    endif

                endif
            enddo
            enddo
            
            if (Me%ExternalVar%Now .ge. Me%NextCompute) then
                !If end, nullify counter. 
                Me%nIterations = 0
            endif

            call ReadUnLockExternalVar(ReadAtmosphere = .true.)


        endif

    end subroutine AverageExtPropDuringDT

    !--------------------------------------------------------------------------

    subroutine NullifyFluxes

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------


        if (Me%ComputeOptions%Grazing) then

            Me%Fluxes%BiomassGrazed                                       (:,:) = 0.0 
            Me%Fluxes%NitrogenGrazed                                      (:,:) = 0.0
            Me%Fluxes%PhosphorusGrazed                                    (:,:) = 0.0 
            Me%Fluxes%BiomassGrazedFraction                               (:,:) = 0.0 
            Me%Fluxes%ToSoil%GrazingBiomassToSoil                         (:,:) = 0.0 
            Me%Fluxes%ToSoil%GrazingNitrogenToSoil                        (:,:) = 0.0 
            Me%Fluxes%ToSoil%GrazingPhosphorusToSoil                      (:,:) = 0.0 

        endif
        
        if (Me%ComputeOptions%Management) then

            Me%Fluxes%BiomassRemovedInHarvest                             (:,:) = 0.0  
            Me%Fluxes%NitrogenRemovedInHarvest                            (:,:) = 0.0 
            Me%Fluxes%PhosphorusRemovedInHarvest                          (:,:) = 0.0
            Me%Fluxes%BiomassHarvestedFraction                            (:,:) = 0.0
            Me%Fluxes%ToSoil%ManagementBiomassToSoil                      (:,:) = 0.0 
            Me%Fluxes%ToSoil%ManagementNitrogenToSoil                     (:,:) = 0.0 
            Me%Fluxes%ToSoil%ManagementPhosphorusToSoil                   (:,:) = 0.0 
            Me%Fluxes%ToSoil%ManagementRootBiomassLeftInSoil              (:,:) = 0.0

        endif

        if (Me%ComputeOptions%Dormancy) then
            
            Me%Fluxes%BiomassRemovedInDormancy                          (:,:) = 0.0
            Me%Fluxes%NitrogenRemovedInDormancy                         (:,:) = 0.0
            Me%Fluxes%PhosphorusRemovedInDormancy                       (:,:) = 0.0
            Me%Fluxes%ToSoil%DormancyBiomassToSoil                      (:,:) = 0.0
            Me%Fluxes%ToSoil%DormancyNitrogenToSoil                     (:,:) = 0.0
            Me%Fluxes%ToSoil%DormancyPhosphorusToSoil                   (:,:) = 0.0

        endif

        Me%Fluxes%WaterUptake           (:,:  ) = 0.0 
        Me%Fluxes%WaterUptakeLayer      (:,:,:) = 0.0 
        
        if (Me%ComputeOptions%ModelNitrogen) then
            Me%Fluxes%NitrogenUptake        (:,:  ) = 0.0 
            Me%Fluxes%NitrogenUptakeLayer   (:,:,:) = 0.0 
        endif
        if (Me%ComputeOptions%ModelNitrogen) then
            Me%Fluxes%PhosphorusUptake      (:,:  ) = 0.0 
            Me%Fluxes%PhosphorusUptakeLayer (:,:,:) = 0.0 
        endif
        if (Me%ComputeOptions%ModelPlantBiomass) then
            Me%Fluxes%BiomassGrowth (:,:) = 0.0
        endif

       
    end subroutine NullifyFluxes

    !--------------------------------------------------------------------------

    subroutine CheckPlantState
 
       !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                            :: i,j, PlantType
        logical                                            :: Dormant
       !Begin-----------------------------------------------------------------
        
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (MappingPoints (i, j) == 1) then
                
                !Reset global state variables
                if (Me%ComputeOptions%Management) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        Me%IsPlantGrowing(i,j)      = .false.
                        Me%KillOccurred(i,j)        = .false.        
                        Me%HarvestKillOccurred(i,j) = .false.
                        if (Me%ComputeOptions%Dormancy) then
                            Me%IsPlantDormant(i,j)      = .false.
                        endif
                    endif
                    if (Me%HarvestOnlyOccurred(i,j)) then
                        Me%HarvestOnlyOccurred(i,j) = .false.
                    endif                
                endif

                if (Me%PlantingOccurred(i,j)) then
                    Me%PlantingOccurred(i,j) = .false.
                endif

                !SWAT Base Heat Units counter (for planting schedule)
                if (Me%ComputeOptions%Evolution%ModelSWAT) then
                    
                    call ComputePlantAnnualStage(i,j)
                
                endif
                
                !Check if planting will occur
                if (.not. Me%IsPlantGrowing(i,j)) then

                    call CheckPlanting(i,j)
                
                endif
    
                if (Me%IsPlantGrowing(i,j)) then
                    
                    !Check tree age counter
                    PlantType  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                    if (PlantType == 7) then
                        
                        call CheckTreeAge(i,j)
                    
                    endif
                    
                    !Check if it is dormant period
                    if (Me%ComputeOptions%Dormancy) then
                        
!                        call CheckPlantDormancy(i,j)
                    
                    endif
                    
                    !Check if it is grazing period
                    if (Me%ComputeOptions%Grazing) then
                        
                        call CheckPlantGrazing(i,j)
                    
                    endif
 
                    !Check if management will occurr (harvest, kill... etc)
                    if (Me%ComputeOptions%Management) then
                        
                        call CheckPlantManagement(i,j)
                    
                    endif


                    !SWAT Growing HU counter
                    if (Me%ComputeOptions%Evolution%ModelSWAT) then
                        
                        Dormant = .false.
                        if (Me%ComputeOptions%Dormancy) then
                            if (Me%IsPlantDormant(i,j)) then
                                Dormant = .true.
                            endif
                        endif
                        
                        if (.not. Dormant) then
                            call ComputePlantGrowingStage(i,j)
                        endif
                    
                    endif
                
                endif
            
            endif
        
        enddo
        enddo


    end subroutine CheckPlantState

    !--------------------------------------------------------------------------

    subroutine ModifyFluxes
 
       !Local-----------------------------------------------------------------
       !Begin-----------------------------------------------------------------
 
        call PlantRootFluxes

        call PlantAerialFluxes
 
        !Harvest and kill operations
        if (Me%ComputeOptions%Management) then
            !Management fluxes are computed if plant is being managed (kill, harvest)            
            call ManagementFluxes
        endif

        if (Me%ComputeOptions%Dormancy) then
            !Dormancy fluxes are computed if plant is going dormant (leaf fall, etc)
            call DormancyFluxes
        endif
    
        if(Me%ComputeOptions%Grazing) then
            !Grazing fluxes are computed if plant is being grazed
            call GrazingFluxes
        endif
    end subroutine ModifyFluxes

    !--------------------------------------------------------------------------

    subroutine ModifyModelledProperties
 
       !Local-----------------------------------------------------------------
       !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%ModelNitrogen     .or.           &
            Me%ComputeOptions%ModelPhosphorus   .or.           &
            Me%ComputeOptions%ModelPlantBiomass       ) then
            
            call UpdateGlobalPlantProperties
        
        endif


        call UpdateRootProperties


        call UpdateLeafProperties

        
        if (Me%ComputeOptions%ModelCanopyHeight) then
            
            call UpdateStemProperties
        
        endif


    end subroutine ModifyModelledProperties

    !--------------------------------------------------------------------------

    subroutine ModifyReadedProperties
 
       !Local-----------------------------------------------------------------
       !Begin-----------------------------------------------------------------
 
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))

            if (PropertyX%IsConstant) then
                PropertyX%Field = PropertyX%ConstantValue

            elseif (PropertyX%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                       Matrix2D       = PropertyX%Field,            &
                                       PointsToFill2D = Me%ExternalVar%MappingPoints2D,      &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadSolution - ModuleAtmosphere - ERR01'
            
            endif
            
            PropertyX => PropertyX%Next
        
        enddo

    end subroutine ModifyReadedProperties

    !--------------------------------------------------------------------------

    subroutine ComputePlantAnnualStage(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        integer                                         :: JulDay, JulDay_Old 
        !SandBoxTest-----------------------------------------------------------
!        Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) = Me%ExternalVar%AirTemperature(i,j)
        Me%HeatUnits%PotentialHUTotal(i,j) = Me%ComputeOptions%PotentialHUTotal !5475. !15ºC * 365 days
        !----------------------------------------------------------------------

        Me%HeatUnits%PotentialHUBase_Old = Me%HeatUnits%PotentialHUBase(i,j)

        call JulianDay(Me%ExternalVar%Now, JulDay)
        JulDay_Old = Me%ExternalVar%JulianDay_Old

        if (JulDay_Old .gt. 364 .and. JulDay .ge. 1) then
            
            !Base heat units nullified
            Me%HeatUnits%PotentialHUBase(i,j) = 0.0
            
        else
            !! update base zero total heat units
            if (Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) .gt. 0. .and.                          &
                Me%HeatUnits%PotentialHUTotal(i,j) .gt. 0.01) then
    
                Me%HeatUnits%PotentialHUBase(i,j) = Me%HeatUnits%PotentialHUBase(i,j)                         &
                + Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) / Me%HeatUnits%PotentialHUTotal(i,j)
            endif

        end if

    end subroutine ComputePlantAnnualStage

    !--------------------------------------------------------------------------

    subroutine CheckTreeAge(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        integer                                         :: JulDay, JulDay_Old 
        integer                                         :: TreeYearsToMaturity
        real                                            :: TreeMaximumBiomass
        !SandBoxTest-----------------------------------------------------------
        !----------------------------------------------------------------------

        call JulianDay(Me%ExternalVar%Now, JulDay)
        JulDay_Old = Me%ExternalVar%JulianDay_Old

        
        !Update if tree was just planted
        if (Me%Growth%TreeCurrentYear(i,j) .eq. 0) then

            !Tree maximum annual biomass update
            TreeYearsToMaturity = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeYearsToMaturity
            
            if (TreeYearsToMaturity .gt. 0) then
                
                Me%Growth%TreeCurrentYear(i,j) = 1           
            
                TreeMaximumBiomass                      = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeMaximumBiomass
                Me%Growth%TreeFractionToMaturity(i,j)   = float(Me%Growth%TreeCurrentYear(i,j)) / float(TreeYearsToMaturity)
                Me%Growth%TreeMaximumAnnualBiomass(i,j) = Me%Growth%TreeFractionToMaturity(i,j) * TreeMaximumBiomass
            else
                Me%Growth%TreeFractionToMaturity(i,j) = 1.
            endif
       
       !Or update at beggining of year
        elseif (JulDay_Old .gt. 364 .and. JulDay .ge. 1) then
            
            !Tree maximum annual biomass update
            TreeYearsToMaturity = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeYearsToMaturity
            
            if (TreeYearsToMaturity .gt. 0) then
                
                Me%Growth%TreeCurrentYear(i,j)          = Me%Growth%TreeCurrentYear(i,j) + 1
                Me%Growth%TreeCurrentYear(i,j)          = min(Me%Growth%TreeCurrentYear(i,j), TreeYearsToMaturity)

                TreeMaximumBiomass                      = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeMaximumBiomass
                Me%Growth%TreeFractionToMaturity(i,j)   = float(Me%Growth%TreeCurrentYear(i,j)) / float(TreeYearsToMaturity)
                Me%Growth%TreeMaximumAnnualBiomass(i,j) = Me%Growth%TreeFractionToMaturity(i,j) * TreeMaximumBiomass
            
            else
                Me%Growth%TreeFractionToMaturity(i,j) = 1.
            end if

        end if

    end subroutine CheckTreeAge

    !--------------------------------------------------------------------------

    subroutine CheckPlanting(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        integer                                         ::         VegetationID
        real                                            ::    PlantingJulianDay
        real                                            ::       PlantingHUBase
        real                                            ::      PotentialHUBase
        real                                            ::  PotentialHUBase_Old
        integer                                         ::   JulDay, JulDay_Old
        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%PlantingJulianDay = -99
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%PlantingHUBase = Me%ComputeOptions%PlantingHUBase
        !Begin-----------------------------------------------------------------


        VegetationID      = Me%VegetationID (i,j)
        PlantingJulianDay = Me%VegetationTypes(VegetationID)%TimingDatabase%PlantingJulianDay
        PlantingHUBase    = Me%VegetationTypes(VegetationID)%TimingDatabase%PlantingHUBase
        PotentialHUBase   = Me%HeatUnits%PotentialHUBase(i,j)

        call JulianDay(Me%ExternalVar%Now, JulDay)
        JulDay_Old = Me%ExternalVar%JulianDay_Old

        !! Check if plant operation. Check if scheduled in days or HU
        if (PlantingJulianDay .gt. 0) then

            if (JulDay .ge. PlantingJulianDay .and. JulDay_Old .lt. PlantingJulianDay) then
                Me%PlantingOccurred(i,j) = .true.
                Me%IsPlantGrowing(i,j)  = .true.
                if (Me%ComputeOptions%Dormancy) then
                    Me%IsPlantDormant(i,j)  = .false.
                endif
            endif

        else if (PlantingHUBase .gt. 0.) then
            PotentialHUBase_Old = Me%HeatUnits%PotentialHUBase_Old
            if (PotentialHUBase .ge. PlantingHUBase .and. PotentialHUBase_Old .lt. PlantingHUBase) then
                Me%PlantingOccurred(i,j) = .true.
                Me%IsPlantGrowing(i,j)   = .true.
                if (Me%ComputeOptions%Dormancy) then
                    Me%IsPlantDormant       = .false.
                endif
            endif

        else 
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Julian day and Base HU for planting' 
            write(*,*    ) 'inconsistently defined in crop database.'
            stop 'CheckIfPlantWillStartGrowing - ModuleVegetation - ERR10'  

        endif

        !Nullify Variables for new growth cycle
        if (Me%PlantingOccurred(i,j)) then
            if (Me%ComputeOptions%ModelPlantBiomass) then
                Me%StateVariables%TotalPlantBiomass%Field      (i,j) = 0.0
            endif
            if (Me%ComputeOptions%ModelNitrogen) then
                Me%StateVariables%TotalPlantNitrogen%Field     (i,j) = 0.0
            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
                Me%StateVariables%TotalPlantPhosphorus%Field   (i,j) = 0.0
            endif
            if (Me%ComputeOptions%ModelRootBiomass) then
                Me%StateVariables%RootBiomass%Field            (i,j) = 0.0
            endif
            Me%StateVariables%RootDepth%Field              (i,j) = 0.0
            Me%StateVariables%LeafAreaIndex%Field          (i,j) = 0.0
!            Me%StateVariables%SpecificLeafStorage%Field    (i,j) = 0.0
            
            Me%HeatUnits%PlantHUAccumulated_Old      (i,j) = 0.0
            Me%HeatUnits%PlantHUAccumulated          (i,j) = 0.0

            Me%PlantLAIMaxFraction                   (i,j) = 0.0
            
            if (Me%ComputeOptions%Management) then
                Me%HarvestOperations       = 1
                Me%HarvestKillOperations   = 1
                Me%KillOperations          = 1
                Me%HarvestFinished(i,j)    = .false.     
            endif           
            if (Me%ComputeOptions%Grazing) then
                Me%GrazingOperations       = 1
                Me%GrazingFinished(i,j)    = .false.
            endif
            
        endif


!        if (Me%PlantingOccurred(i,j)) then
!REVIEW NULLIFY VARIABLES
!        Me%
!        Me%PlantingOccurred(i,j) = .false. !melhor aqui ou no fim?
!        endif

      

!          j = 0
!        j = ihru

!        igro(j) = 1
!        idorm(j) = 0
!        phuacc(j) = 0.
!        plantn(j) = 0.
!        plantp(j) = 0.
!        plt_et(j) = 0.
!        plt_pet(j) = 0.                                          
!        laimxfr(j) = 0.
!        hvstiadj(j) = 0.
!        olai(j) = 0.
!        rwt(j) = 0.
!        icr(j) = icr(j) + 1

        !! initialize transplant variables
!        if (lai_init(nro(j),icr(j),j) > 0.) then
!          laiday(j) = lai_init(nro(j),icr(j),j)
!          bio_ms(j) = bio_init(nro(j),icr(j),j)
!        endif

    end subroutine CheckPlanting

    !--------------------------------------------------------------------------

    subroutine ComputePlantGrowingStage(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        integer                                         ::         VegetationID
        real                                            ::     PlantHUVariation
        real                                            ::    PlantHUatMaturity
        real                                            ::  AverageTempDuringDT
        real                                            :: PlantBaseTemperature
        real                                            :: HUAcc, AccumulatedTemperature
        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%PlantHUatMaturity = Me%ComputeOptions%PlantHUatMaturity
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantBaseTemperature = 0.
        !Begin-----------------------------------------------------------------

        VegetationID         = Me%VegetationID(i,j)
        PlantHUatMaturity    = Me%VegetationTypes(VegetationID)%TimingDatabase%PlantHUatMaturity
        AverageTempDuringDT  = Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j)
        PlantBaseTemperature = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantBaseTemperature

        !! update accumulated heat units for the plant
        PlantHUVariation = 0.
        if (PlantHUatMaturity .gt. 0.1) then
            AccumulatedTemperature = AverageTempDuringDT * (Me%ComputeOptions%VegetationDT/86400.)
            PlantHUVariation = (AverageTempDuringDT - PlantBaseTemperature) / PlantHUatMaturity
        end if
        if (PlantHUVariation .lt. 0.) then
            PlantHUVariation = 0.
        endif
        Me%HeatUnits%PlantHUAccumulated_Old(i,j) = Me%HeatUnits%PlantHUAccumulated (i,j)
        Me%HeatUnits%PlantHUAccumulated    (i,j) = Me%HeatUnits%PlantHUAccumulated (i,j) + PlantHUVariation
        !Debug
        HUAcc = Me%HeatUnits%PlantHUAccumulated (i,j)


    end subroutine ComputePlantGrowingStage

    !--------------------------------------------------------------------------

    subroutine PlantRootFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%Evolution%ModelSWAT) then
    
            call WaterUptake
             
         !! if plant is not dormant nor hasn't reached maturity, nutrient assimilation is allowed. 
            if (Me%ComputeOptions%ModelNitrogen) then
                
                call NitrogenUptakeSWAT

                call NitrogenFixationSWAT

            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
            
                call PhosphorusUptakeSWAT

            endif
        endif            


    end subroutine PlantRootFluxes

    !--------------------------------------------------------------------------

    subroutine WaterUptake

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%TranspirationMethod .eq. TranspirationSWAT) then
            
            call WaterUptakeSWAT
        
        elseif (Me%ComputeOptions%TranspirationMethod .eq. TranspirationMOHID) then
            
            call WaterUptakeMOHID
        
        else
            write(*,*) 
            write(*,*) 'Transpiration method not known. Verify input'
            stop 'WaterUptake - ModuleVegetation - ERR100'
        endif

    end subroutine WaterUptake
    
    !--------------------------------------------------------------------------


    subroutine PlantAerialFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%Evolution%ModelSWAT) then

         !! if plant is not dormant nor hasn't reached maturity, this fluxes are allowed. 
            if (Me%ComputeOptions%ModelPlantBiomass) then
                
                call BiomassGrowthFromRadiationSWAT
           
            endif

            call LAIGrowthSWAT

        endif            


    end subroutine PlantAerialFluxes

    !--------------------------------------------------------------------------


!    subroutine CheckPlantDormancy(i,j)

        !Arguments-------------------------------------------------------------
!        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
!        !! check for beginning of dormant season
!        if (.not. Me%IsPlantDormant(i,j) .and. 

!! check for beginning of dormant season
!      if (idorm(j) == 0 .and. dayl(j)-dormhr(j) < daylmn(hru_sub(j)))   &
!     &                                                              then
!
!        select case (idc(idplt(nro(j),icr(j),j)))
!
!          !! beginning of forest dormant period
!          case (7)
!            idorm(j) = 1
!            resnew = 0.
!            resnew = bio_ms(j) * bio_leaf(idplt(nro(j),icr(j),j))
!            sol_rsd(1,j) = sol_rsd(1,j) + resnew
!            sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
!            sol_fon(1,j) = resnew * pltfr_n(j) + sol_fon(1,j)
!            sol_fop(1,j) = resnew * pltfr_p(j) + sol_fop(1,j)
!            bio_hv(nro(j),icr(j),j) = bio_ms(j) +                       &
!     &                                           bio_hv(nro(j),icr(j),j)
!            bio_yrms(j) = bio_yrms(j) + bio_ms(j) / 1000.
!            bio_ms(j) = .2 * bio_ms(j) *                                &
!     &                           (1. - bio_leaf(idplt(nro(j),icr(j),j)))
!            plantn(j) = 0.
!            plantp(j) = 0.
!            strsw(j) = 1.
!            laiday(j) = alai_min(idplt(nro(j),icr(j),j))
!            phuacc(j) = 0.
!
!          !! beginning of perennial (pasture/alfalfa) dormant period
!          case (3, 6)
!            idorm(j) = 1
!            resnew = 0.
!            bm_dieoff = 0.7
!            resnew = bm_dieoff * bio_ms(j)
!            sol_rsd(1,j) = sol_rsd(1,j) + resnew
!            sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
!            sol_fon(1,j) = sol_fon(1,j) + bm_dieoff * plantn(j)
!            sol_fop(1,j) = sol_fop(1,j) + bm_dieoff * plantp(j)
!            bio_hv(nro(j),icr(j),j) = bio_ms(j) * bm_dieoff +           &
!    &                                           bio_hv(nro(j),icr(j),j)
!            bio_yrms(j) = bio_yrms(j) + bio_ms(j) * bm_dieoff / 1000.
!            bio_ms(j) = (1. - bm_dieoff) * bio_ms(j)
!            plantn(j) = (1. - bm_dieoff) * plantn(j)
!            plantp(j) = (1. - bm_dieoff) * plantp(j)
!            strsw(j) = 1.
!            laiday(j) = alai_min(idplt(nro(j),icr(j),j))
!            phuacc(j) = 0.
!
!          !! beginning of cool season annual dormant period
!          case (2, 5)
!            if (phuacc(j) < 0.75) then
!              idorm(j) = 1
!              strsw(j) = 1.
!            end if
!
!        end select
!      end if


!! check if end of dormant period
!       if (idorm(j) == 1 .and. dayl(j)-dormhr(j) >= daylmn(hru_sub(j)))&
!     &                                                              then
!
!          select case (idc(idplt(nro(j),icr(j),j)))
!          
!            !! end of perennial dormant period
!            case (3, 6, 7)
!              idorm(j) = 0
!
!            !! end of cool season annual dormant period
!            case (2, 5)
!              idorm(j) = 0
!
!          end select
!
!        end if


!    end subroutine CheckPlantDormancy

    !--------------------------------------------------------------------------
    
    subroutine DormancyFluxes
       
        !Arguments-------------------------------------------------------------
!        integer, dimension(:,:), pointer                   :: MappingPoints
!        integer                                         :: i, j
        !Begin-----------------------------------------------------------------

!        MappingPoints => Me%ExternalVar%MappingPoints2D

!do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
!            if (MappingPoints (i, j) == 1 .and. Me%PlantGoingDormant(i,j)) then        
! IT IS IMPORTANT TO DO A PREDICTIVE BIOMASS AND NUTRIENT VALUE 
!TAKING INTO ACCOUNT MANAGEMENT TO ASSURE THAT VALUES NOT NEGATIVE.  IF KILL DONT REMOVE BIOMASS
! IF HARVEST PREDICT BIOMASS AND NUTRIENT
    
    end subroutine DormancyFluxes
   
    !--------------------------------------------------------------------------

    subroutine WaterUptakeSWAT

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i,j
        real                                            :: TopDepth, BottomDepth
        real                                            ::              PotTP
        logical                                         ::         FoundRoot
        real                                            :: NormalizationParameter
        real                                            ::            RootDepth
        real                                            :: TopUptake, BottomUptake
        real                                            :: PotentialWaterUptake
        real                                            :: EffectiveWaterUptake
        real                                            :: GridCellArea, CellVolume
        real                                            ::   LayerFieldCapacity
        real                                            ::    LayerWaterContent
        real                                            :: SumDemand, SumUptake
        real                                            :: DemandNotMetInUpperLayers
        real                                            :: IncreaseUptake, ReduceUptake
        integer                                         :: k, KUB, KLB
        logical                                         :: UptakeAllowed
        !SandBoxTest-----------------------------------------------------------
!        Me%ExternalVar%FieldCapacity(:,:,:) = 0.3 !m3/m3
!        Me%ExternalVar%ResidualWaterContent(:,:,:) = 0.0 
!        Me%ExternalVar%PotentialTranspiration(i,j) = 1e-8 !m/s - 1 mm/dia
        !Begin-----------------------------------------------------------------
        
!        if (ep_max <= 0.01) then
!        strsw(j) = 1.
!        else
        !! initialize variables
!        gx = 0.
!        ir = 0
!        sump = 0.
!        wuse = 0.
!        xx = 0.

        !!  compute aeration stress
!       if (sol_sw(j) > sol_sumfc(j)) then
!          satco = (sol_sw(j) - sol_sumfc(j)) / (sol_sumul(j) - 
!     &                                                 sol_sumfc(j))
!          strsa(j) = 1. - (satco / (satco + Exp(.176 - 4.544 *
!     &                                                      satco)))
!        else
!          strsa(j) = 1.
!        end if
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            UptakeAllowed = .true.
            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                !To by-pass no growing season periods. 
                if (.not. Me%IsPlantGrowing(i,j)) then
                    UptakeAllowed = .false.
                endif
            endif
               
            if (MappingPoints (i, j) == 1 .and. UptakeAllowed) then        
                
!                Me%ExternalVar%PotentialTranspiration(i,j) = 1e-8 !m/s - 1 mm/dia
                ! mm         = m/s * s * 1000mm/m
                PotTP        = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT * 1000.0 
        
cd1:            if (PotTP .le. 0.01) then

                    Me%Fluxes%WaterUptake(i,j) = 0.0
                    Me%Fluxes%FromSoil%WaterUptakeFromSoil => Me%Fluxes%WaterUptake
                    Me%Growth%WaterStress(i,j) = 1.0

                else
                
                    BottomDepth  = 0.
                    FoundRoot = .false.
                    
                    if (Me%StateVariables%RootDepth%Field(i,j) .lt. 0.01) then
                        Me%StateVariables%RootDepth%Field(i,j) = 0.01
                    endif
                    
                    !                               m
                    RootDepth    = Me%StateVariables%RootDepth%Field(i,j)
!                    RootDepth     = PropRootDepth%Field(i,j)

                    GridCellArea = Me%ExternalVar%GridCellArea(i,j)
                    SumDemand    = 0.
                    SumUptake    = 0.

                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3 :               do k = KUB, KLB, -1
            
                        if (FoundRoot) then
                            exit do3
                        endif
            
                        !Potential Water uptake is computed with the difference between top and bottom cell uptakes
                        TopDepth = BottomDepth
                        BottomDepth = BottomDepth + Me%ExternalVar%DWZ(i,j,k)
            
                        !except when root end is found. bottom depth used is root end
                        if (RootDepth .le. BottomDepth ) then
                            BottomDepth = RootDepth
                            FoundRoot = .true.
                        endif
                     
                        NormalizationParameter = 1. - exp(-10.0)
                        if (RootDepth .le. 1e-5) then
                            !       mm
                            PotentialWaterUptake = PotTP/NormalizationParameter
                        else
                            !      mm
                            BottomUptake = PotTP *(1.0 - exp(-10.0 * BottomDepth/RootDepth))/NormalizationParameter
                            TopUptake = PotTP *(1.0 - exp(-10.0 * TopDepth/RootDepth))/NormalizationParameter
                            PotentialWaterUptake = BottomUptake - TopUptake
                        endif

            !         !! don't allow compensation for aeration stress
            !          if (strsa(j) > .99) then
            !            yy = 0.
            !          else
            !            yy= sump - xx
            !          end if
            !          wuse(k) = sum - sump + yy * epco(j)
            !          wuse(k) = sum - sump + (sump - xx) * epco(j)
            !          sump = sum

                        !!Effective uptake will account for compensation from lower layers if demand is not met
                        !in upper layers. Increase to uptake is applied. WaterUptakeCompensationFactor ranges from zero to one
                        !to account for the deviation to the potential uptake allowed.
                        !Compare the sum of previous demand with previous uptake
                        DemandNotMetInUpperLayers = SumDemand - SumUptake
                        IncreaseUptake = DemandNotMetInUpperLayers * Me%ComputeOptions%WaterUptakeCompensationFactor

                        ! Demand for next iteration. mm
                        SumDemand = SumDemand + PotentialWaterUptake

                        !! Effective uptake will be updated if low soil water content occurs. Reduction to uptake is
                        !! applyed if occurs
                        !! adjust uptake if sw is less than 25% of plant available water

                        CellVolume = Me%ExternalVar%CellVolume(i,j,k)
                        !      mm          =        m3Water/m3cell               *       m3cell   /  m2      *  mm/m
                        LayerFieldCapacity = Me%ExternalVar%FieldCapacity(i,j,k) * (CellVolume/GridCellArea) * 1000
                        LayerWaterContent = (Me%ExternalVar%SoilWaterContent(i,j,k)                                 &
                                             - Me%ExternalVar%ResidualWaterContent(i,j,k)) * (CellVolume/GridCellArea) * 1000
            
                        ReduceUptake = 0.
                        if (LayerWaterContent .lt. LayerFieldCapacity/4.) then
                
                            ReduceUptake = exp(5. * (4. * LayerWaterContent / LayerFieldCapacity - 1.))
                        else
                            ReduceUptake = 1.
                        endif
            
                        !        mm
                        EffectiveWaterUptake = PotentialWaterUptake * ReduceUptake + IncreaseUptake
                        if (LayerWaterContent .lt. EffectiveWaterUptake) then
                            EffectiveWaterUptake = LayerWaterContent
                        end if
                        !  mm
                        SumUptake = SumUptake + EffectiveWaterUptake
                        !    m3/s                         =          mm          *    m/mm *   m2   *     1/DT
                        Me%Fluxes%WaterUptakeLayer(i,j,k) = EffectiveWaterUptake *  (1./1000.) * GridCellArea * (1./Me%ComputeOptions%VegetationDT)

                    enddo do3
        
                    !    m3/s
                    Me%Fluxes%WaterUptake(i,j) = SumUptake *  (1./1000.) * GridCellArea * (1./Me%ComputeOptions%VegetationDT)
                    Me%Fluxes%FromSoil%WaterUptakeFromSoil => Me%Fluxes%WaterUptake

            !       strsw(j) = strsa (j) * xx / ep_max
                    Me%Growth%WaterStress(i,j) = SumUptake/PotTP
        
                endif cd1

            endif

        enddo do2
        enddo do1


    end subroutine WaterUptakeSWAT

    !--------------------------------------------------------------------------

    subroutine WaterUptakeMOHID

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        real, dimension(:,:,:), pointer                 :: UnsatK
        integer                                         :: i, j
        real                                            :: SumUptake, RootDepth, BottomDepth
        real                                            :: TopDepth, CenterDepth, TotalCol
        real                                            :: FeddesH1, FeddesH2, FeddesH3, FeddesH4
        real                                            :: Head, Factor, h50, p1
        real                                            :: Aux, RootHeightInCell, CellBase
        real                                            :: ColToTransp, TranspirationDistribution
        real                                            :: SoilVolume, WaterVolume, TranspVolume
        real                                            :: ActualTranspiration, PotentialTranspiration
        logical                                         :: FoundRoot, UptakeAllowed
        integer                                         :: k, KLB, KUB, RootProfile
        real                                            :: VelocityVolume, NewTheta
        real                                            :: LimitThetaLow
        integer                                         :: STAT_CALL
        !Begin-----------------------------------------------------------------


!                if (Me%ExtVar%WaterPoints3D(i, j, Me%WorkSize%KUB) == 1  ) then   
! 
!  !first of all,it calculates the layer of the roots. 
!              
!   ! calculate root layer!
!       deep = 0
!       w    = Me%WorkSize%KUB
!       seguir = .True.
!
!       Do While( w>0.and. seguir==.True.) 
!
!       deep = deep + Me%ExtVar%DWZ(i,j,w)
!       y = Me%ExtVar%rootdepth(i,j)
!
!           write(*,*)'A profunidade ate onde descem as raizes é', deep            
!            If(y-deep<=0.001 .or. deep>y) then
!            
!               layerroot = (Me%WorkSize%KUB-w) +1
!               seguir = .False.
!                    
!           write(*,*)'Alayer das raizes e', layerroot         
!
!            Else
!
!            w = w-1
!
!           End If 
!       Enddo
        
        MappingPoints => Me%ExternalVar%MappingPoints2D
        
        !Fluxes Nullification
        Me%Fluxes%WaterUptakeLayer(:,:,:) = 0.0
        Me%Fluxes%WaterUptake     (:,:  ) = 0.0

        call GetLimitThetaLow(Me%ObjPorousMedia, LimitThetaLow, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterUptakeMOHID - ModuleVegetation - ERR000'   
        
        if (Me%ComputeOptions%LimitTPVel) then 

            call GetUnsatK(Me%ObjPorousMedia, UnsatK, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WaterUptakeMOHID - ModuleVegetation - ERR001'   
        
        endif

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            UptakeAllowed = .true.    
            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                !To by-pass no growing season periods. 
                if (.not. Me%IsPlantGrowing(i,j)) then
                    UptakeAllowed = .false.
                endif
            endif

            if (MappingPoints (i, j) == 1 .and. UptakeAllowed) then        

                !Total Water Column to Remove (Potentialy) 
                !m = m/s * s
                TotalCol    = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT    
                !  m
                RootDepth    = Me%StateVariables%RootDepth%Field(i,j)
                BottomDepth = 0.0                   
                
cd1:            if (TotalCol .eq. 0.0 .or. RootDepth .eq. 0.0) then

                    Me%Fluxes%WaterUptake(i,j) = 0.0
                    Me%Growth%WaterStress(i,j) = 1.0

                else
                 
                    SumUptake = 0.0
                    FoundRoot = .false.

!                    if (Me%StateVariables%RootDepth%Field(i,j) .lt. 0.01) then
!                        Me%StateVariables%RootDepth%Field(i,j) = 0.01
!                    endif

                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3 :               do k = KUB, KLB, -1
        
                        if (FoundRoot) then
                            exit do3
                        endif

                        !Layer depths to consider trnaspiration
                        TopDepth        = BottomDepth
                        BottomDepth     = BottomDepth + Me%ExternalVar%DWZ(i, j, k)

                        !Center of layer
                        CenterDepth    = BottomDepth - Me%ExternalVar%DWZ(i,j,k) / 2.0 


                        if (RootDepth .le. BottomDepth ) then
                            FoundRoot = .true.
                            BottomDepth = RootDepth
                        endif

             
                        !Reduce Factor of EVTP due to water stress
                        if (Me%ComputeOptions%WaterUptakeStressMethod .eq. Feddes) then   

                            ! stress functions by Feddes
                            FeddesH1 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH1(i,j)
                            FeddesH2 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH2(i,j)
                            FeddesH3 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH3(i,j)
                            FeddesH4 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH4(i,j)
                            Head     = Me%ExternalVar%Head(i, j, k)

                            if      (Head .gt. FeddesH1) then
                                Factor = 0.0
                            else if (Head .gt. FeddesH2) then
                                Factor = LinearInterpolation(FeddesH2, 1.0, FeddesH1, 0.0, Head)
                            else if (Head .gt. FeddesH3) then
                                Factor = 1.0
                            else if (Head .gt. FeddesH4) then
                                Factor = LinearInterpolation(FeddesH4, 0.0, FeddesH3, 1.0, Head)
                            else
                                Factor = 0.0
                            endif
              
                        elseif (Me%ComputeOptions%WaterUptakeStressMethod .eq. Genuchten) then  
            
                            !or stress functions by S-Shape (van Genuchten)
                            h50 =-8.0
                            p1 =3.0
                            Factor= 1/(1+(Head/h50)**p1)

                        endif
                    
                        !If transpiration occurs
                        if (Factor .gt. 0.0) then
                                
                            !Water lost by transpiration according to different root density distribution
                            RootProfile = Me%ComputeOptions%RootProfile

                            if (RootProfile == RootTriangular) then ! triangular (base on top)
           
                               !base of "transpiration triangle"
                                aux             = 2.0 * TotalCol / RootDepth
            
                                !This method computes succesive transpiration triangle bases (inside the "transpiration triangle") 
                                !at the center of each cell. the transpiration is the base * cell height
                                !The method was changed because underestimated transpiration in last cell when root end was 
                                !different from cell boundaries
                                !Error made: assumed a root height inside the cell equal to cell height. Also in root end, root 
                                !geometric form(inside the cell) is a triangle (base*height/2.0).
                                if (Me%ComputeOptions%WaterUptakeOld) then
                                    RootHeightInCell = Me%ExternalVar%DWZ(i,j,k)
                                    CellBase         = aux * (RootDepth - CenterDepth) / RootDepth
                                    ColToTransp      = Factor * (CellBase * RootHeightInCell)
                                else
                                    if (FoundRoot) then
                                        RootHeightInCell = RootDepth - TopDepth
                                        CellBase         = aux * (RootDepth - TopDepth) / RootDepth
                                        ColToTransp      = Factor * (CellBase * RootHeightInCell)/2.0
                                    else 
                                        RootHeightInCell = Me%ExternalVar%DWZ(i,j,k)
                                        CellBase         = aux * (RootDepth - CenterDepth) / RootDepth
                                        ColToTransp      = Factor * (CellBase * RootHeightInCell)
                                    endif
                                endif

                            elseif (RootProfile == RootConstant) then !constant root distribution
            
                                !This method was changed because overestimated transpiration in last cell when root end was 
                                !different from cell boundaries
                                !Error made: assumed a root height inside the cell equal to cell height. 
                                if (Me%ComputeOptions%WaterUptakeOld) then
                                    RootHeightInCell = Me%ExternalVar%DWZ(i,j,k)
                                else
                                    if (FoundRoot) then
                                        RootHeightInCell = RootDepth - TopDepth
                                    else
                                        RootHeightInCell = Me%ExternalVar%DWZ(i,j,k)
                                    endif
                                endif

                                ColToTransp          =  Factor * TotalCol * (RootHeightInCell / RootDepth)

                            !Same method as first but fraction computed with layer numbers - if root exists inside layer implicitily 
                            !occupies all layer. 
                            !Removed because it would need more variables to compute and the geometry is the same as first option.
                !           elseif (RootProfile == 3) then !linear root distribution
                !   
                !                ColToEvap      =  2 * Factor * TotalCol * (1. / layerroot - (prof/layerroot**2))              
                ! 
                            elseif (RootProfile == RootExponential) then  ! exponential (from SWAT)

                !                a =10
                !                ColToEvap = Factor*TotalCol*a*exp(-a*(Me%ExtVar%RootDepth(i,j)-deep2))

                                TranspirationDistribution = (1.0 - exp(-10.0 * BottomDepth/RootDepth))/(1.0 - exp(-10.0))                   &
                                                           - (1.0 - exp(-10.0 * TopDepth/RootDepth))/(1.0 - exp(-10.0))

                                ColToTransp = Factor * TotalCol * TranspirationDistribution
                            else
                               write(*,*)'Invalid Root Profile'
                               stop 'WaterUptakeMOHID - ModuleVegetation - ERR02'
                            endif


                            !Potential Water to transpirate from layer
                            ! m3 = m * m2
                            WaterVolume = ColToTransp * Me%ExternalVar%GridCellArea (i,j)
                            
                            !Velocity Volume
                            if (Me%ComputeOptions%LimitTPVel) then
                                VelocityVolume = UnSatK (i,j,k) * Me%ComputeOptions%VegetationDT * Me%ExternalVar%GridCellArea(i,j)
                                TranspVolume    = min(WaterVolume, VelocityVolume)
                            else
                                TranspVolume    = WaterVolume
                            endif        
                            
                            SoilVolume  = (Me%ExternalVar%SoilWaterContent(i,j,k) - Me%ExternalVar%ResidualWaterContent(i,j,k))           &
                                            * Me%ExternalVar%CellVolume(i,j,k)                

                            TranspVolume = min(TranspVolume, SoilVolume)
                            
                            !Estimates new Theta
                            NewTheta = Me%ExternalVar%SoilWaterContent(i,j,k) - TranspVolume / Me%ExternalVar%CellVolume(i,j,k)
                        
                            !Don't let dry... stability reasons
                            if (NewTheta .le. Me%ExternalVar%ResidualWaterContent(i,j,k) + LimitThetaLow) then
        
                                TranspVolume = 0.0
                            
                            endif

                    
                        else !No transpiration

                            TranspVolume = 0.0

                        endif
        
                        !m3
                        SumUptake = SumUptake + TranspVolume
        
                        !m3/s = m3  / s
                        Me%Fluxes%WaterUptakeLayer(i,j,k) = TranspVolume / Me%ComputeOptions%VegetationDT

                    enddo do3
    
                    !    m3/s
                    Me%Fluxes%WaterUptake(i,j) = SumUptake / Me%ComputeOptions%VegetationDT
!                    Me%Fluxes%FromSoil%WaterUptakeFromSoil => Me%Fluxes%WaterUptake
    
                    ! m  = m3/m2
                    ActualTranspiration = SumUptake / Me%ExternalVar%GridCellArea(i,j)
                    !m  = m/s * s
                    PotentialTranspiration = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT
                    Me%Growth%WaterStress(i,j) = ActualTranspiration/PotentialTranspiration
            
                endif cd1
            
            endif

        enddo do2
        enddo do1

        if (Me%ComputeOptions%LimitTPVel) then 

            call UnGetPorousMedia(Me%ObjPorousMedia, UnsatK, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WaterUptakeMOHID - ModuleVegetation - ERR003'   
        
        endif

    end subroutine WaterUptakeMOHID

    !--------------------------------------------------------------------------

    subroutine NitrogenUptakeSWAT

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        integer                                         ::         VegetationID
        real                                            :: TotalPlantNitrogen, TotalPlantBiomass
        real                                            :: PlantFractionN1, PlantFractionN2
        real                                            :: PlantFractionN3
        real                                            :: b1, b2, b3
        real                                            :: PlantShape1, PlantShape2
        real                                            :: NitrogenDemand, OptimalNContent
        real                                            :: Stress, HUAcc
        real                                            :: TopDepth, BottomDepth
        logical                                         ::         FoundRoot
        real                                            :: NormalizationParameter
        real                                            ::            RootDepth
        real                                            :: TopUptake, BottomUptake
        real                                            :: PotentialNitrogenUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: DistributionParameter
        real                                            ::   LayerNitrogenContent
        real                                            :: SumDemand, SumUptake
        real                                            :: DemandNotMetInUpperLayers
        integer                                         :: k, KUB, KLB, PlantType
        logical                                         :: Dormant
       !Begin-----------------------------------------------------------------

        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            Dormant = .false.
            if (Me%ComputeOptions%Dormancy) then
                if (Me%IsPlantDormant(i,j)) then
                    Dormant = .true.
                endif
            endif
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j) .and. .not. Dormant             &
                .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then  
                
         !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantFractionN1 = 0.0663 
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantFractionN2 = 0.0255 
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantFractionN3 = 0.0148 
        Me%ExternalVar%SoilNitrogenConcentration(:,:,:)                         = 1.0 !kg/m3H20??
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType       = 5      
        !------------------------------------------------------------------------

                VegetationID    = Me%VegetationID(i,j)
                PlantFractionN1 = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantFractionN1
                PlantFractionN2 = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantFractionN2
                PlantFractionN3 = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantFractionN3
    
                !Parameters for Nitrogen Shape curve
                b1 = PlantFractionN1 - PlantFractionN3       
                b2 = 1. - (PlantFractionN2 - PlantFractionN3) / b1
                b3 = 1. - .00001 / b1
                call ComputeShapeCoefficients(b2, b3, 0.5, 1.0, PlantShape1, PlantShape2)
    
                HUAcc = Me%HeatUnits%PlantHUAccumulated (i,j)
                Me%PlantNitrogenFraction(i,j) = (PlantFractionN1 - PlantFractionN3)                                &
                                                * (1. - HUAcc / (HUAcc + exp(PlantShape1 - PlantShape2 *  HUAcc)))  &
                                                + PlantFractionN3

                !   kgN/ha
                OptimalNContent = 0.
                TotalPlantBiomass  = Me%StateVariables%TotalPlantBiomass%Field(i,j)
                TotalPlantNitrogen = Me%StateVariables%TotalPlantNitrogen%Field(i,j)
                OptimalNContent = Me%PlantNitrogenFraction(i,j) * TotalPlantBiomass
                if (OptimalNContent .lt. TotalPlantNitrogen) then
                    OptimalNContent = TotalPlantNitrogen
                endif
    
                !Compute Nitrogen Stress
                PlantType = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                select case (PlantType)
                    !Legumes
                    case (1,2,3)
                        Me%Growth%NitrogenStress(i,j) = 1.
                    case default
                        call ComputeNutrientStress(TotalPlantNitrogen,OptimalNContent,Stress)
                        Me%Growth%NitrogenStress(i,j) = Stress
        !          if (uno3d > 1.e-5) then
        !            xx = nplnt(j) / uno3d
        !          else
        !            xx = 1.
        !          end if
        !          strsn = amax1(strsn, xx)
        !          strsn = amin1(strsn, 1.)
                end select

                !  kgN/ha
                NitrogenDemand = OptimalNContent - TotalPlantNitrogen
                NitrogenDemand = min(4. * PlantFractionN3 * Me%Growth%BiomassGrowthOld(i,j), NitrogenDemand)
    
cd1:            if(NitrogenDemand .lt. 1e-6) then
        
                    Me%Fluxes%NitrogenUptake(i,j) = 0.0
    
                else
    
                    BottomDepth = 0.0
                    FoundRoot = .false.
                    if (Me%StateVariables%RootDepth%Field(i,j) .lt. 0.01) then
                        Me%StateVariables%RootDepth%Field(i,j) = 0.01
                    endif
                    !                               m
                    RootDepth             = Me%StateVariables%RootDepth%Field(i,j)
!                    RootDepth             = PropRootDepth%Field(i,j)
                    SumDemand             = 0.
                    SumUptake             = 0.
                    GridCellArea          = Me%ExternalVar%GridCellArea(i,j)
                    DistributionParameter = Me%ComputeOptions%NitrogenDistributionParameter

                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3:                do k = KUB, KLB, -1
        
                        if (FoundRoot) then
                            exit do3
                        endif
        
                        !Potential Nitrogen uptake is computed with the difference between top and bottom cell uptakes
                        TopDepth = BottomDepth
                        BottomDepth = BottomDepth + Me%ExternalVar%DWZ(i,j,k)
        
                        !except when root end is found. bottom depth used is root end
                        if (RootDepth .le. BottomDepth ) then
                            BottomDepth = RootDepth
                            FoundRoot = .true.
                        endif
                 
                        NormalizationParameter = 1. - exp(-DistributionParameter)
                        !   kgN/ha
                        BottomUptake = NitrogenDemand *(1.0 - exp(-DistributionParameter * BottomDepth/RootDepth))     &
                                        /NormalizationParameter
                        TopUptake = NitrogenDemand *(1.0 - exp(-DistributionParameter * TopDepth/RootDepth))           &
                                        /NormalizationParameter
                        PotentialNitrogenUptake = BottomUptake - TopUptake

                        DemandNotMetInUpperLayers = SumDemand - SumUptake
                        !Demand in next iteration
                        SumDemand = SumDemand + PotentialNitrogenUptake

                        CellVolume = Me%ExternalVar%CellVolume(i,j,k)
        
                        !      kgN/ha        = kgN/m3H2O * m3H20/m3cell * m3cell / (m2) * 10000m2/ha                     
                        LayerNitrogenContent = Me%ExternalVar%SoilNitrogenConcentration(i,j,k)                &
                                                * Me%ExternalVar%SoilWaterContent(i,j,k)                      &
                                                * (CellVolume) / (GridCellArea) * 10000

                        Me%Fluxes%NitrogenUptakeLayer(i,j,k) = PotentialNitrogenUptake + DemandNotMetInUpperLayers
                        if (LayerNitrogenContent .lt. Me%Fluxes%NitrogenUptakeLayer(i,j,k)) then
                            Me%Fluxes%NitrogenUptakeLayer(i,j,k) = LayerNitrogenContent
                        end if
        
                        SumUptake = SumUptake + Me%Fluxes%NitrogenUptakeLayer(i,j,k)
                    enddo do3

                    Me%Fluxes%NitrogenUptake(i,j) = SumUptake
                    Me%Fluxes%FromSoil%NitrogenUptakeFromSoil => Me%Fluxes%NitrogenUptake
    
                endif cd1

            endif

        enddo do2
        enddo do1


    end subroutine NitrogenUptakeSWAT

    !--------------------------------------------------------------------------

    subroutine ComputeShapeCoefficients(Factor1, Factor2, Factor3, Factor4, PlantShape1, PlantShape2)

        !Arguments-------------------------------------------------------------
        real, intent (IN)                              :: Factor1
        real, intent (IN)                              :: Factor2
        real, intent (IN)                              :: Factor3
        real, intent (IN)                              :: Factor4
        real, intent (OUT)                             :: PlantShape1
        real, intent (OUT)                             :: PlantShape2
        !Local-----------------------------------------------------------------
        real                                           :: temp
        !Begin-----------------------------------------------------------------
        
        temp = 0.0
        PlantShape1 = 0.0
        PlantShape2 = 0.0

        temp = Log(Factor3/Factor1 - Factor3)
        PlantShape2 = (temp - Log(Factor4/Factor2 - Factor4)) / (Factor4 - Factor3)
        PlantShape1 = temp + (Factor3 * PlantShape2)

    end subroutine ComputeShapeCoefficients

    !--------------------------------------------------------------------------

    subroutine ComputeNutrientStress (PlantContent, OptimalContent, ComputedStress)
        
        !Arguments-------------------------------------------------------------
        real, intent (IN)                              :: PlantContent
        real, intent (IN)                              :: OptimalContent
        real, intent (OUT)                             :: ComputedStress
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------

        ComputedStress = 0.
        
        !Added this condition. SWAT produces error when content and optimal are small
        !Because of added .0001 in denominator for computed stress
        if (PlantContent .eq. OptimalContent) then
            ComputedStress = 1.0
        else
        
!Error            ComputedStress = 200. * (PlantContent / (OptimalContent + .0001) - .5)
            
            ComputedStress = 200. * (PlantContent / (OptimalContent) - .5)

            if (ComputedStress .lt. 0.) then
                ComputedStress = 0.
            elseif (ComputedStress < 99.) then
                ComputedStress = ComputedStress / (ComputedStress + Exp(3.535 - .02597 * ComputedStress))
            else
                ComputedStress = 1.
            endif

            if (OptimalContent .lt. 1.e-6) then
                ComputedStress = 1.
            endif
        endif

    end subroutine ComputeNutrientStress

    !--------------------------------------------------------------------------

    subroutine PhosphorusUptakeSWAT

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        integer                                         ::         VegetationID
        real                                            :: TotalPlantBiomass, TotalPlantPhosphorus
        real                                            :: PlantFractionP1, PlantFractionP2
        real                                            :: PlantFractionP3
        real                                            :: b1, b2, b3
        real                                            :: PlantShape1, PlantShape2
        real                                            :: PhosphorusDemand, OptimalPContent
        real                                            :: Stress, HUAcc
        real                                            :: TopDepth, BottomDepth
        logical                                         ::         FoundRoot
        real                                            :: NormalizationParameter
        real                                            ::            RootDepth
        real                                            :: TopUptake, BottomUptake
        real                                            :: PotentialPhosphorusUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: DistributionParameter
        real                                            ::   LayerPhosphorusContent
        real                                            :: SumDemand, SumUptake
        real                                            :: DemandNotMetInUpperLayers
        integer                                         :: k, KUB, KLB
        logical                                         :: Dormant
       !Begin-----------------------------------------------------------------
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            Dormant = .false.
            if (Me%ComputeOptions%Dormancy) then
                if (Me%IsPlantDormant(i,j)) then
                    Dormant = .true.
                endif
            endif
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j) .and. .not. Dormant                         &
                .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then    

         !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantFractionP1 = 0.0053 
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantFractionP2 = 0.0020 
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantFractionP3 = 0.0012 
!        Me%PhosphorusDistributionParameter                                      = -20.0
        Me%ExternalVar%SoilPhosphorusConcentration(:,:,:)                       = 0.1 !kg/m3H20??
        !----------------------------------------------------------------------- 
        
                VegetationID    = Me%VegetationID(i,j)
                PlantFractionP1 = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantFractionP1
                PlantFractionP2 = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantFractionP2
                PlantFractionP3 = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantFractionP3
        
                !Parameters for Nitrogen Shape curve
                b1 = PlantFractionP1 - PlantFractionP3       
                b2 = 1. - (PlantFractionP2 - PlantFractionP3) / b1
                b3 = 1. - .00001 / b1
                call ComputeShapeCoefficients(b2, b3, 0.5, 1.0, PlantShape1, PlantShape2)
        
                HUAcc = Me%HeatUnits%PlantHUAccumulated (i,j)
                Me%PlantPhosphorusFraction(i,j) = (PlantFractionP1 - PlantFractionP3)        &
                                                * (1. - HUAcc/ (HUAcc + exp(PlantShape1 - PlantShape2 *  HUAcc)))  &
                                                + PlantFractionP3

                !   kgP/ha
                OptimalPContent      = 0.
                TotalPlantBiomass    = Me%StateVariables%TotalPlantBiomass%Field(i,j)
                TotalPlantPhosphorus = Me%StateVariables%TotalPlantPhosphorus%Field(i,j)
                OptimalPContent = Me%PlantPhosphorusFraction(i,j) * TotalPlantBiomass
                if (OptimalPContent .lt. TotalPlantPhosphorus) then
                    OptimalPContent = TotalPlantPhosphorus
                endif
        
                !Compute phosphorus stress
                call ComputeNutrientStress(TotalPlantPhosphorus, OptimalPContent, Stress)
                Me%Growth%PhosphorusStress(i,j) = Stress

                !  kgP/ha
                PhosphorusDemand = OptimalPContent - TotalPlantPhosphorus
                PhosphorusDemand = min(4. * PlantFractionP3 * Me%Growth%BiomassGrowthOld(i,j), PhosphorusDemand)
                !Luxury P uptake
                PhosphorusDemand =  PhosphorusDemand * 1.5
        
cd1:            if(PhosphorusDemand .lt. 1e-6) then
            
                    Me%Fluxes%PhosphorusUptake(i,j) = 0.0
        
                else
        
                    BottomDepth  = 0.0
                    FoundRoot = .false.
                    if (Me%StateVariables%RootDepth%Field(i,j) .lt. 0.01) then
                        Me%StateVariables%RootDepth%Field(i,j) = 0.01
                    endif
                    !                               m
                    RootDepth    = Me%StateVariables%RootDepth%Field(i,j)
!                    RootDepth    = PropRootDepth%Field(i,j)
                    SumDemand    = 0.
                    SumUptake    = 0.
                    GridCellArea = Me%ExternalVar%GridCellArea(i,j)
                    DistributionParameter = Me%ComputeOptions%PhosphorusDistributionParameter
        
                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3 :               do k = KUB, KLB, -1
            
                        if (FoundRoot) then
                            exit do3
                        endif
            
                        !Potential Phosphorus uptake is computed with the difference between top and bottom cell uptakes
                        TopDepth = BottomDepth
                        BottomDepth = BottomDepth + Me%ExternalVar%DWZ(i,j,k)
            
                        !except when root end is found. bottom depth used is root end
                        if (RootDepth .le. BottomDepth ) then
                            BottomDepth = RootDepth
                            FoundRoot = .true.
                        endif
                     
                        NormalizationParameter = 1. - exp(-DistributionParameter)
                        !   kgP/ha
                        BottomUptake = PhosphorusDemand *(1.0 - exp(-DistributionParameter * BottomDepth/RootDepth))           &
                                       /NormalizationParameter
                        TopUptake = PhosphorusDemand *(1.0 - exp(-DistributionParameter * TopDepth/RootDepth))                 &
                                        /NormalizationParameter
                        PotentialPhosphorusUptake = BottomUptake - TopUptake
  
                        DemandNotMetInUpperLayers = SumDemand - SumUptake
                        !Demand in next iteration
                        SumDemand = SumDemand + PotentialPhosphorusUptake

                        CellVolume = Me%ExternalVar%CellVolume(i,j,k)
            
                        !      kgP/ha        = kgN/m3H2O * m3H20/m3cell * m3cell / (m2) * 10000m2/ha                     
                        LayerPhosphorusContent = Me%ExternalVar%SoilPhosphorusConcentration(i,j,k)                &
                                                * Me%ExternalVar%SoilWaterContent(i,j,k)                      &
                                                * (CellVolume) / (GridCellArea) * 10000

                        Me%Fluxes%PhosphorusUptakeLayer(i,j,k) = PotentialPhosphorusUptake + DemandNotMetInUpperLayers
                        if (LayerPhosphorusContent .lt. Me%Fluxes%PhosphorusUptakeLayer(i,j,k)) then
                            Me%Fluxes%PhosphorusUptakeLayer(i,j,k) = LayerPhosphorusContent
                        end if
            
                        SumUptake = SumUptake + Me%Fluxes%PhosphorusUptakeLayer(i,j,k)
                    enddo do3

                    Me%Fluxes%PhosphorusUptake(i,j) = SumUptake
                    Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil => Me%Fluxes%PhosphorusUptake

                endif cd1

            endif

        enddo do2
        enddo do1


    end subroutine PhosphorusUptakeSWAT

    !--------------------------------------------------------------------------

    subroutine NitrogenFixationSWAT

        !Local-----------------------------------------------------------------
!        integer, dimension(:,:), pointer                     :: MappingPoints
!        integer                                           :: i, j

!        MappingPoints => Me%ExternalVar%MappingPoints2D
!
!do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!                
!            if (MappingPoints (i, j) == 1 .and. .not. Me%IsPlantDormant(i,j)                              &
!                .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then     
 

!    !! if crop is a legume, call nitrogen fixation routine
!      select case (idc(idplt(nro(j),icr(j),j)))
!        case (1,2,3)
!         call nfix
!      end select
!
!      nplnt(j) = nplnt(j) + fixn
!      plantn(j) = plantn(j) + nplnt(j)
!
!
!            j = 0
!      j = ihru
! 
!! compute the difference between supply and demand
!      if (uno3d > nplnt(j)) then
!        uno3l = 0.
!        uno3l = uno3d - nplnt(j)
!      else
!        !! if supply is being met, fixation=0 and return
!        fixn = 0.
!        return
!      endif
!
!! compute fixation as a function of no3, soil water, and growth stage
!
!      !! compute soil water factor
!      fxw = 0.
!      fxw = sol_sw(j) / (.85 * sol_sumfc(j))
!
!      !! compute no3 factor
!      sumn = 0.
!      fxn = 0.
!      do l = 1, sol_nly(j)
!        sumn = sumn + sol_no3(l,j)
!      end do
!      if (sumn > 300.) fxn = 0.
!      if (sumn > 100. .and. sumn <= 300.) fxn = 1.5 - .0005 * sumn
!      if (sumn <= 100.) fxn = 1.
!
!      !! compute growth stage factor
!      fxg = 0.
!      if (phuacc(j) > .15 .and. phuacc(j) <= .30) then
!         fxg = 6.67 * phuacc(j) - 1.
!      endif
!      if (phuacc(j) > .30 .and. phuacc(j) <= .55) fxg = 1.
!      if (phuacc(j) > .55 .and. phuacc(j) <= .75) then
!         fxg = 3.75 - 5. * phuacc(j)
!      endif
!
!      fxr = Min(1., fxw, fxn) * fxg
!      fxr = Max(0., fxr)
!
!      fixn = Min(6., fxr * uno3l)
!      fixn = fixco * fixn + (1. - fixco) * uno3l
!             !! if fixco=0 then fix the entire demand
!      fixn = Min(fixn, uno3l)
!

    end subroutine NitrogenFixationSWAT

    !--------------------------------------------------------------------------

    subroutine BiomassGrowthFromRadiationSWAT

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                     :: MappingPoints
        integer                                           :: i, j
        real                                              ::        ExtinctCoef
        integer                                           ::       VegetationID
        real                                              :: LAI, PAR, RUE
        real                                              :: BiomassEnergyRatio
        real                                              :: BiomassEnergyRatioHigh
        real                                              :: CO2ConcHigh, RUEDeclineRate
        real                                              :: Decline, PotentialBiomassGrowth
        real                                              :: SolarRadiation  
        real                                              :: b1, b2, c1 
        real                                              :: PlantShape1, PlantShape2
        real                                              :: VapourPressureDeficit
        real                                              :: SaturatedVapourPressure
        real                                              :: ActualVapourPressure
        real                                              :: PlantType
        logical                                           :: Dormant
       !Begin-----------------------------------------------------------------

        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            Dormant = .false.
            if (Me%ComputeOptions%Dormancy) then
                if (Me%IsPlantDormant(i,j)) then
                    Dormant = .true.
                endif
            endif

            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j) .and. .not. Dormant                        &
                .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then     

        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%ExtinctCoef = 0.65
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%BiomassEnergyRatio = 30.0
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%CO2ConcHigh = 660.
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%BiomassEnergyRatioHigh = 39.0
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%RUEDeclineRate = 6.0
!        Me%ExternalVar%Integration%AverageRadiationDuringDT(i,j) = Me%ExternalVar%SolarRadiation(i,j)
!        Me%ExternalVar%Integration%AverageAirHumidityDuringDT(i, j) = Me%ExternalVar%RelativeHumidity(i,j)
        !----------------------------------------------------------------------- 
        
                VegetationID           =   Me%VegetationID(i,j)
                ExtinctCoef            =   Me%VegetationTypes(VegetationID)%GrowthDatabase%ExtinctCoef
                BiomassEnergyRatio     =   Me%VegetationTypes(VegetationID)%GrowthDatabase%BiomassEnergyRatio
                CO2ConcHigh            =   Me%VegetationTypes(VegetationID)%GrowthDatabase%CO2ConcHigh
                BiomassEnergyRatioHigh =   Me%VegetationTypes(VegetationID)%GrowthDatabase%BiomassEnergyRatioHigh
                RUEDeclineRate         =   Me%VegetationTypes(VegetationID)%GrowthDatabase%RUEDeclineRate
                !  m2/m2
                LAI                    =   Me%StateVariables%LeafAreaIndex%Field(i,j)
!                LAI                    =   PropLeafAreaIndex%Field(i,j)
                !  MJ/m2               =   W/m2 * s * 1e-6MJ/J
                SolarRadiation         =   Me%ExternalVar%Integration%AverageRadiationDuringDT(i,j) * Me%ComputeOptions%VegetationDT * 1e-6
                PlantType              =   Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
       
                !! calculate optimal biomass
                !Photossintetically Active Radiation (MJ/m2)
                PAR = 0.
                PAR = .5 * SolarRadiation * (1. - Exp(-ExtinctCoef * (LAI + .05)))

                !Parameters used to construct the radiation efficience use shape curve
                c1 = 330.                                 !! ambient CO2 concentration
                b1 = BiomassEnergyRatio * .01             !! "ambient" bio-e ratio/100
                b2 = BiomassEnergyRatioHigh * .01         !! "elevated" bio-e ratio/100
                call ComputeShapeCoefficients(b1, b2, c1, CO2ConcHigh, PlantShape1, PlantShape2)

                !Adjust radiation-use efficiency for CO2
                !RUE units in (kg/ha)/(MJ/m**2)
                !CO2 conc units in uL CO2/L air
                RUE = 0.
                if (Me%ComputeOptions%AtmosphereCO2 .gt. 330.) then
                    RUE = (100. * Me%ComputeOptions%AtmosphereCO2) / (Me%ComputeOptions%AtmosphereCO2 +         &
                                  Exp(PlantShape1 - Me%ComputeOptions%AtmosphereCO2) * PlantShape2)     
                else
                    RUE = BiomassEnergyRatio
                end if

                !! adjust radiation-use efficiency for vapor pressure deficit. 
                !Compute VPD (KPa). Equations in ModuleBasin adapted to average atmosphere values
                SaturatedVapourPressure  = 0.6108 * exp (17.27 * Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j)              &
                                            / (Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) + 237.3))
        
                ActualVapourPressure     = SaturatedVapourPressure * Me%ExternalVar%Integration%AverageAirHumidityDuringDT(i, j)  
        
                VapourPressureDeficit    = SaturatedVapourPressure - ActualVapourPressure
        
                !!assumes vapor pressure threshold of 1.0 kPa
                if (VapourPressureDeficit > 1.0) then
                    RUE = 0.
                    Decline = VapourPressureDeficit - 1.0
                    RUE = RUE - RUEDeclineRate * Decline
                    RUE = Max(RUE, 0.27 * BiomassEnergyRatio)
                end if

                PotentialBiomassGrowth = RUE * PAR
                if (PotentialBiomassGrowth .lt. 0.) then
                    PotentialBiomassGrowth = 0.
                endif

        !!MOVE TO MANAGEMENT
        !          !! auto fertilization-nitrogen demand (non-legumes only)
        !          select case (idc(idp))
        !            case (4, 5, 6, 7)
        !            if (auto_nstrs(j) > 0.) call anfert
        !          end select

                !! reduce predicted biomass due to stress on plant
                call ComputeGlobalStress(i,j)


        !BIOTARGET NOT INCLUDED - overrides stresses and goes to defined biomass
        !          if (bio_targ(nro(j),icr(j),j) > 1.e-2) then
        !            bioday = bioday * (bio_targ(nro(j),icr(j),j) - bio_ms(j)) / &
        !     &                                         bio_targ(nro(j),icr(j),j)
        !            reg = 1.
        !          end if

                Me%Fluxes%BiomassGrowth(i,j) = PotentialBiomassGrowth * Me%Growth%GlobalStress(i,j)
                !Duplication to be used in nutrient uptake. 
                Me%Growth%BiomassGrowthOld(i,j)     = Me%Fluxes%BiomassGrowth(i,j)

        !TREE CYCLE NOT YET USED - TO DO - COUNTER TO TREE AGE
        !        if (PlantType == 7) then
        !           if (mat_yrs(idp) > 0) then
        !              rto = float(curyr_mat(j)) / float(mat_yrs(idp))
        !              biomxyr = rto * bmx_trees(idp)
        !              bio_ms(j) = Min (bio_ms(j), biomxyr)
        !            else
        !              rto = 1.
        !            end if
        !          end if

            endif

        enddo do2
        enddo do1

    end subroutine BiomassGrowthFromRadiationSWAT

    !--------------------------------------------------------------------------

    subroutine ComputeGlobalStress(i,j)
        
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                 :: i, j
        !Local-----------------------------------------------------------------
        real                                              :: TemperatureStress, NitrogenStress
        real                                              :: PhosphorusStress, WaterStress
        real                                              :: GlobalStress
        !begin-----------------------------------------------------------------


        WaterStress = Me%Growth%WaterStress(i,j)
        
        if (Me%ComputeOptions%ModelTemperatureStress) then
            call ComputeTemperatureStress(i,j)
            TemperatureStress = Me%Growth%TemperatureStress(i,j)
        else
            TemperatureStress = 1.0
        endif
        if (Me%ComputeOptions%ModelNitrogen) then
            NitrogenStress = Me%Growth%NitrogenStress(i,j)
        else
            NitrogenStress = 1.0
        endif
        if (Me%ComputeOptions%ModelPhosphorus) then
            PhosphorusStress = Me%Growth%PhosphorusStress(i,j)
        else
            PhosphorusStress = 1.0
        endif

        GlobalStress = 0.

        GlobalStress = min(WaterStress, TemperatureStress, NitrogenStress, PhosphorusStress)
        if (GlobalStress .lt. 0.) then
            GlobalStress = 0.
        elseif (GlobalStress .gt. 1.) then
            GlobalStress = 1.
        endif

        Me%Growth%GlobalStress(i,j) = GlobalStress
    
    end subroutine ComputeGlobalStress

    !--------------------------------------------------------------------------

    subroutine ComputeTemperatureStress(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                                 :: i, j
        !Local-----------------------------------------------------------------
        real                                                :: BaseTemperature, OptimalTemperature
        real                                                :: AverageTempDuringDT
        real                                                :: Temperature1, Temperature2
        !SandTest--------------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantOptimalTemperature = 18.0
        !Begin-----------------------------------------------------------------
        
        BaseTemperature     = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantBaseTemperature
        OptimalTemperature  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantOptimalTemperature
        AverageTempDuringDT = Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j)

        Temperature1 = 0.
        Temperature1 = AverageTempDuringDT - BaseTemperature

        if (Temperature1 .lt. 0.) then
            Me%Growth%TemperatureStress(i,j) = 0.
        elseif (AverageTempDuringDT .gt. OptimalTemperature) then
            Temperature1 = 2. * OptimalTemperature - BaseTemperature - AverageTempDuringDT
        end if

        Temperature2 = 0.
        Temperature2 = ((OptimalTemperature - AverageTempDuringDT) / (Temperature1 + 1.e-6)) ** 2

        if (Temperature2 .le. 200. .and. Temperature1 .gt. 0.) then
            Me%Growth%TemperatureStress(i,j) = Exp(-0.1054 * Temperature2)
        else
            Me%Growth%TemperatureStress(i,j) = 0.
        end if

!        if (tmn(j) <= tmp_an(hru_sub(j)) - 15.) strstmp(j) = 0.

!        end if



    end subroutine ComputeTemperatureStress

    !--------------------------------------------------------------------------

    subroutine LAIGrowthSWAT

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        integer                                         :: VegetationID
        real                                            :: FrLAIMax1, FrLAIMax2
        real                                            :: FrGrow1, FrGrow2, HUAcc, HUAcc_Old
        real                                            :: LAIShape1, LAIShape2
        real                                            :: FractionLAIMax_Old
        real                                            :: FrGrowLAIDecline
        integer                                         :: PlantType
        real                                            :: DBLAIMax
        real                                            :: LAIMax, DeltaLAI, LAI
        real                                            :: DeltaFractionLAIMax
        real                                            :: GlobalStress
        logical                                         :: Dormant, HasLeaves
        !Begin----------------------------------------------------------
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            Dormant = .false.
            if (Me%ComputeOptions%Dormancy) then
                if (Me%IsPlantDormant(i,j)) then
                    Dormant = .true.
                endif
            endif
            HasLeaves = .false.
            if (MappingPoints (i, j) == 1 .and. Me%VegetationTypes(Me%VegetationID(i,j))%HasLeaves) then
                HasLeaves = .true.
            endif

            if (MappingPoints (i, j) == 1  .and. HasLeaves .and. Me%IsPlantGrowing(i,j) .and. .not. Dormant                         &
                .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then  

        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%FrLAIMax1 = 0.05
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%FrLAIMax2 = 0.95
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%FrGrow1   = 0.05
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%FrGrow2   = 0.45
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%FrGrowLAIDecline = 0.50
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%LAIMax    = 4.0
        !----------------------------------------------------------------------
            
                VegetationID = Me%VegetationID(i,j)
                FrLAIMax1 = Me%VegetationTypes(VegetationID)%GrowthDatabase%FrLAIMax1
                FrLAIMax2 = Me%VegetationTypes(VegetationID)%GrowthDatabase%FrLAIMax2
                FrGrow1   = Me%VegetationTypes(VegetationID)%GrowthDatabase%FrGrow1
                FrGrow2   = Me%VegetationTypes(VegetationID)%GrowthDatabase%FrGrow2
        
                call ComputeShapeCoefficients (FrLAIMax1, FrLAIMax2, FrGrow1, FrGrow2, LAIShape1, LAIShape2)
        
                HUAcc = Me%HeatUnits%PlantHUAccumulated(i,j)
                HUAcc_Old = Me%HeatUnits%PlantHUAccumulated_Old(i,j)

                !Save old value. It will be used in next step
                FractionLAIMax_Old = Me%PlantLAIMaxFraction(i,j)
                !Update
                Me%PlantLAIMaxFraction(i,j) = HUAcc / (HUAcc + exp(LAIShape1 - LAIShape2 * HUAcc))
        
                !Compute new Leaf Area Index
                FrGrowLAIDecline = Me%VegetationTypes(VegetationID)%GrowthDatabase%FrGrowLAIDecline
                PlantType        = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                DBLAIMax         = Me%VegetationTypes(VegetationID)%GrowthDatabase%LAIMax

                Me%LAISenescence(i,j) = .false.
         
                if (HUAcc .lt. FrGrowLAIDecline) then
            
                    LAIMax = 0.
                    DeltaLAI = 0.
            
                    if (PlantType .eq. 7) then
                        LAIMax = DBLAIMax * Me%Growth%TreeFractionToMaturity(i,j)
                    else
                        LAIMax = DBLAIMax
                    endif
            
                    DeltaFractionLAIMax = Me%PlantLAIMaxFraction(i,j) - FractionLAIMax_Old
                    
                    call ComputeGlobalStress(i,j)
                    GlobalStress = Me%Growth%GlobalStress(i,j)
                    
                    LAI = Me%StateVariables%LeafAreaIndex%Field(i,j)
!                    LAI = PropLeafAreaIndex%Field(i,j)

                    DeltaLAI = (DeltaFractionLAIMax * LAIMax * (1.0 - exp(5.0 * (LAI - LAIMax)))) * sqrt(GlobalStress)

                    Me%Fluxes%LAIGrowth(i,j) = DeltaLAI
        
                else
            
                    Me%LAISenescence(i,j) = .true.
            
                    if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                        Me%Fluxes%LAIDeclineFraction(i,j) = (1.0 - HUAcc) / (1.0 - FrGrowLAIDecline)
                    else
                        if (HUAcc .gt. HUAcc_Old) then
                            Me%Fluxes%LAIDeclineFraction(i,j) = (1.0 - HUAcc) / (1.0 - HUAcc_Old)
                        else
                            Me%Fluxes%LAIDeclineFraction(i,j) = 1.0
                        endif
                    endif
                endif

            endif

        enddo do2
        enddo do1


    end subroutine LAIGrowthSWAT

    !--------------------------------------------------------------------------

    subroutine CheckPlantGrazing(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            ::               DT_day
        real                                            :: GrazingStartJulianDate
!        real                                            ::   GrazingStartPlantHU
        integer                                         ::  GrazingDays, JulDay, JulDay_Old
        real                                            :: HUAcc, HUAcc_Old
        real, dimension (:), pointer                    :: GrazingStartPlantHU
        integer                                         :: Op

        !SandBoxTest-----------------------------------------------------------
!        allocate(Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartPlantHU(1))
       
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartJulianDay = -99.
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartPlantHU(1) = 0.5
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartPlantHU(2) = 0.8

        !Begin-----------------------------------------------------------------
        
        if (.not. Me%GrazingFinished(i,j)) then
            
            Op = Me%GrazingOperations(i,j)
            

            !! if cell currently not grazed, check to see if it is time
            !! to initialize grazing
            if (.not.Me%IsPlantBeingGrazed(i,j)) then

                call JulianDay(Me%ExternalVar%Now, JulDay)
                JulDay_Old = Me%ExternalVar%JulianDay_Old

                GrazingStartJulianDate = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartJulianDay 
                GrazingStartPlantHU    =>Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartPlantHU                 
                DT_day                 = Me%ComputeOptions%VegetationDT / (60 * 60 * 24)
                HUAcc                  = Me%HeatUnits%PlantHUAccumulated(i,j)
                HUAcc_Old              = Me%HeatUnits%PlantHUAccumulated_Old(i,j)

                !Check if grazing start is scheduled for day 
                if (GrazingStartJulianDate .gt. 0) then
                
                    if(JulDay .ge. GrazingStartJulianDate .and. JulDay_Old .lt. GrazingStartJulianDate) then
                
                        Me%IsPlantBeingGrazed(i,j) = .true.
                        Me%DaysOfGrazing(i,j) = 1
                
                    endif

                else if (GrazingStartPlantHU(Op) .gt. 0.) then
                
                    if(HUAcc .ge. GrazingStartPlantHU(Op) .and. HUAcc_Old .lt. GrazingStartPlantHU(Op)) then
                
                        Me%IsPlantBeingGrazed(i,j) = .true.
                        Me%DaysOfGrazing(i,j) = 1
                
                    endif
                else
                    Me%GrazingFinished(i,j) = .true. 

                endif
                
            else

                !! if not first day of grazing increment total days of grazing by one
                DT_day                 = Me%ComputeOptions%VegetationDT / (60 * 60 * 24)
                Me%DaysOfGrazing(i,j)  = Me%DaysOfGrazing(i,j) + DT_day
                GrazingDays            = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingDays 
                !! check to set if grazing period is over
                if (Me%DaysOfGrazing(i,j) .ge. GrazingDays) then
    
                    Me%IsPlantBeingGrazed(i,j) = .false.
    
                    Me%DaysOfGrazing(i,j) = 0.

                    Me%GrazingOperations(i,j) = Me%GrazingOperations(i,j) + 1

                    if (Me%GrazingOperations(i,j) .gt. size(GrazingStartPlantHU)) then
                        Me%GrazingFinished(i,j) = .true. 
                    endif

                end if

            end if

        endif

!        deallocate(Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingStartPlantHU)
    
    end subroutine CheckPlantGrazing

    !--------------------------------------------------------------------------

    subroutine GrazingFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer             :: MappingPoints
        integer                                   :: i, j
        real                                      :: TotalPlantBiomass, PredictedBiomass
        real                                      :: BiomassRemovedInDormancy, BiomassRemovedInHarvest
        real                                      :: TotalPlantNitrogen, PredictedNitrogen
        real                                      :: NitrogenRemovedInDormancy, NitrogenRemovedInHarvest
        real                                      :: TotalPlantPhosphorus, PredictedPhosphorus
        real                                      :: PhosphorusRemovedInDormancy, PhosphorusRemovedInHarvest
        real                                      ::GrazingMinimumBiomass
        real                                      ::       GrazingBiomass
        real                                      ::     TramplingBiomass
        real                                      ::        BiomassGrazed
        real                                      ::       NitrogenGrazed
        real                                      ::     PhosphorusGrazed
        real                                      ::      BiomassTrampled
        real                                      ::     NitrogenTrampled
        real                                      ::   PhosphorusTrampled
        real                                      :: BiomassGrazedFraction
        character (Len = StringLength)            :: WarningString
        logical                                   :: PlantKilled
        !Begin-----------------------------------------------------------------
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            PlantKilled = .false.
            if (Me%HarvestKillOccurred(i,j) .or. Me%KillOccurred(i,j)) then
                PlantKilled = .true.
            endif
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j) .and. Me%IsPlantBeingGrazed(i,j)      &
                .and. .not. PlantKilled) then     

        !SandBoxTest-----------------------------------------------------------
!         Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingDays = 10
!         Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingMinimumBiomass = 10.
!         Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingBiomass = 70.
!         Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%TramplingBiomass = 30.
        !----------------------------------------------------------------------
                
                GrazingMinimumBiomass = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingMinimumBiomass
                TotalPlantBiomass     = Me%StateVariables%TotalPlantBiomass%Field(i,j)
                
                BiomassRemovedInDormancy = 0.0
                if (Me%ComputeOptions%Dormancy .and. Me%PlantGoingDormant(i,j)) then
                    BiomassRemovedInDormancy = Me%Fluxes%BiomassRemovedInDormancy(i,j)
                endif
                BiomassRemovedInHarvest = 0.0
                if (Me%ComputeOptions%Management .and. Me%HarvestOnlyOccurred(i,j)) then
                    BiomassRemovedInHarvest = Me%Fluxes%BiomassRemovedInHarvest(i,j)
                endif
                PredictedBiomass = TotalPlantBiomass - BiomassRemovedInDormancy - BiomassRemovedInHarvest
                 
                !! graze only if adequate biomass in HRU
                if (PredictedBiomass .gt. GrazingMinimumBiomass) then

                    !Kg/ha
                    GrazingBiomass        = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%GrazingBiomass
                    TramplingBiomass      = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%TramplingBiomass
            
                    !! biomass grazing flux
                    BiomassGrazed = GrazingBiomass

                    !! Grazing only occurs until minimum biomass for grazing
                    if ((PredictedBiomass - GrazingBiomass) .lt. GrazingMinimumBiomass) then
                        BiomassGrazed = PredictedBiomass - GrazingMinimumBiomass
                    endif

                    !! biomass trampling flux - debris that is not eaten but taken from plant and deposited 
                    !in soil (related to grazing efficiency) 
                    BiomassTrampled = TramplingBiomass

                    !! Trampling only ocurrs until minimum biomass for grazing
                    if ((PredictedBiomass - BiomassGrazed - TramplingBiomass) &
                        .lt. GrazingMinimumBiomass) then
                        BiomassTrampled = PredictedBiomass - BiomassGrazed - GrazingMinimumBiomass
                    endif


                    !!Fill fluxes
                    !!Plant point of view - Sum grazing and trampling
                    Me%Fluxes%BiomassGrazed(i,j) = BiomassGrazed + BiomassTrampled

                    !Soil Point of View - Recieve only trampling fluxes. Manure in future
                    Me%Fluxes%ToSoil%GrazingBiomassToSoil(i,j) = BiomassTrampled

                    
                    if (Me%ComputeOptions%ModelNitrogen) then
                        
                        TotalPlantNitrogen    = Me%StateVariables%TotalPlantNitrogen%Field(i,j)

                        NitrogenRemovedInDormancy = 0.0
                        if (Me%ComputeOptions%Dormancy .and. Me%PlantGoingDormant(i,j)) then
                            NitrogenRemovedInDormancy = Me%Fluxes%NitrogenRemovedInDormancy(i,j)
                        endif
                        BiomassRemovedInHarvest = 0.0
                        if (Me%ComputeOptions%Management .and. Me%HarvestOnlyOccurred(i,j)) then
                            NitrogenRemovedInHarvest = Me%Fluxes%NitrogenRemovedInHarvest(i,j)
                        endif
                        
                        PredictedNitrogen = TotalPlantNitrogen - NitrogenRemovedInDormancy - NitrogenRemovedInHarvest

                        NitrogenGrazed = BiomassGrazed * Me%PlantNitrogenFraction (i,j)
                        if (NitrogenGrazed .gt. PredictedNitrogen) then
                            NitrogenGrazed  =   PredictedNitrogen
                        endif                        
                        
                        NitrogenTrampled = BiomassTrampled * Me%PlantNitrogenFraction(i,j)
                        if ((NitrogenGrazed + NitrogenTrampled) .gt. PredictedNitrogen) then
                            NitrogenTrampled  =   PredictedNitrogen - NitrogenGrazed
                        endif                               
                        
                        Me%Fluxes%NitrogenGrazed(i,j) = NitrogenGrazed + NitrogenTrampled
                        
                        Me%Fluxes%ToSoil%GrazingNitrogenToSoil(i,j) = NitrogenTrampled


                    endif
                    if (Me%ComputeOptions%ModelPhosphorus) then

                        TotalPlantPhosphorus  = Me%StateVariables%TotalPlantNitrogen%Field(i,j)

                        PhosphorusRemovedInDormancy = 0.0
                        if (Me%ComputeOptions%Dormancy .and. Me%PlantGoingDormant(i,j)) then
                            PhosphorusRemovedInDormancy = Me%Fluxes%PhosphorusRemovedInDormancy(i,j)
                        endif
                        PhosphorusRemovedInHarvest = 0.0
                        if (Me%ComputeOptions%Management .and. Me%HarvestOnlyOccurred(i,j)) then
                            PhosphorusRemovedInHarvest = Me%Fluxes%PhosphorusRemovedInHarvest(i,j)
                        endif
                        
                        PredictedPhosphorus = TotalPlantPhosphorus - PhosphorusRemovedInDormancy - PhosphorusRemovedInHarvest
                        
                        PhosphorusGrazed = BiomassGrazed * Me%PlantPhosphorusFraction (i,j)
                        if (PhosphorusGrazed .gt. PredictedPhosphorus) then
                            PhosphorusGrazed  =   PredictedPhosphorus
                        endif                        

                        PhosphorusTrampled = BiomassTrampled * Me%PlantPhosphorusFraction (i,j)
                        if ((PhosphorusGrazed + PhosphorusTrampled) .gt. PredictedPhosphorus) then
                            PhosphorusTrampled  =   PredictedPhosphorus - PhosphorusGrazed
                        endif                               
                         
                        Me%Fluxes%PhosphorusGrazed(i,j) = PhosphorusGrazed + PhosphorusTrampled

                        Me%Fluxes%ToSoil%GrazingPhosphorusToSoil(i,j) = PhosphorusTrampled

                    endif

                    

        ! SEND TO FLUXES TO SOIL - DEPENDS ON BIOMASS EVALUATION
        !                !! remove trampled biomass and add to residue
        !                dmii = 0.
        !                dmii = bio_ms(j)
        !                bio_ms(j) = bio_ms(j) - bio_trmp(nro(j),ngr(j),j)
        !                if (bio_ms(j) < bio_min(j))  then
        !                  sol_rsd(1,j) = sol_rsd(1,j) + dmii - bio_min(j)
        !                  bio_ms(j) = bio_min(j)
        !                else
        !                  sol_rsd(1,j) = sol_rsd(1,j) + bio_trmp(nro(j),ngr(j),j)
        !                endif
        !                sol_rsd(1,j) = Max(sol_rsd(1,j),0.)
        !                bio_ms(j) = Max(bio_ms(j),0.)

                    !! adjust nutrient content of residue and biomass for
                    !! trampling
        !                if (dmii - bio_ms(j) > 0.) then
        !                  sol_fon(1,j) = (dmii - bio_ms(j)) * pltfr_n(j) + sol_fon(1,j)
        !                 sol_fop(1,j) = (dmii - bio_ms(j)) * pltfr_p(j) + sol_fop(1,j) 
        !                end if


        !RECHECK AFTER TO INCLUDE MANURE PROCESSES
        !        !! apply manure
        !       it = 0
        !        it = manure_id(nro(j),ngr(j),j)
        !        if (manure_kg(nro(j),ngr(j),j) > 0.) then
        !          l = 1
        !
        !          sol_no3(l,j) = sol_no3(l,j) + manure_kg(nro(j),ngr(j),j) *    &
        !     &                 (1. - fnh3n(it)) * fminn(it)
        !          sol_fon(l,j) = sol_fon(l,j) + manure_kg(nro(j),ngr(j),j) *    &
        !     &                 forgn(it)
        !          sol_nh3(l,j) = sol_nh3(l,j) + manure_kg(nro(j),ngr(j),j) *    &
        !     &                 fnh3n(it) * fminn(it)
        !          sol_solp(l,j) = sol_solp(l,j) + manure_kg(nro(j),ngr(j),j) *  &
        !     &                 fminp(it)
        !          sol_fop(l,j) = sol_fop(l,j) + manure_kg(nro(j),ngr(j),j) *    &
        !     &                 forgp(it)
        !
        !! add bacteria - #cfu/g * t(manure)/ha * 1.e6 g/t * ha/10,000 m^2 = 100.
        !! calculate ground cover
        !          gc = 0.
        !          gc = (1.99532 - Erfc(1.333 * laiday(j) - 2.)) / 2.1
        !          if (gc < 0.) gc = 0.
        !
        !          gc1 = 0.
        !          gc1 = 1. - gc
        !
        !          swf = .15
        !
        !          frt_t = 0.
        !          frt_t = bact_swf * manure_kg(nro(j),ngr(j),j) / 1000.
        !
        !          bactp_plt(j) = gc * bactpdb(it) * frt_t * 100. + bactp_plt(j)
        !          bactlp_plt(j) = gc * bactlpdb(it) * frt_t * 100.+bactlp_plt(j)
        !
        !         bactpq(j) = gc1 * bactpdb(it)  * frt_t * 100. + bactpq(j)
        !          bactpq(j) = bactkddb(it) * bactpq(j)
        !
        !          bactps(j) = gc1 * bactpdb(it) * frt_t * 100. + bactps(j)
        !          bactps(j) = (1. - bactkddb(it)) * bactps(j)
        !
        !          bactlpq(j) = gc1 * bactlpdb(it) * frt_t * 100. + bactlpq(j)     
        !          bactlpq(j) = bactkddb(it) * bactlpq(j)
        !
        !          bactlps(j) = gc1 * bactlpdb(it) * frt_t * 100. + bactlps(j)
        !          bactlps(j) = (1. - bactkddb(it)) * bactlps(j)
        !
        !        endif

                    !! leaf area index and fraction of growing season fluxes
                    BiomassGrazedFraction = (BiomassGrazed + BiomassTrampled)/TotalPlantBiomass
        !                if (TotalPlantBiomass .gt. 1.) then
                    Me%Fluxes%BiomassGrazedFraction(i,j) = BiomassGrazedFraction
                    
                    WarningString = 'Graze'
                    call UpdatePlantGrowingStage(i,j, WarningString)
        !                else
        !                  laiday(j) = 0.05
        !                  phuacc(j) = 0.
        !                endif


                
                endif

            endif
        
        enddo do2
        enddo do1


    end subroutine GrazingFluxes

    !--------------------------------------------------------------------------

    subroutine CheckPlantManagement(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: HUAcc, HUAcc_Old
        real                                            :: HarvestKillJulianDay 
        real                                            :: HarvestJulianDay, KillJulianDay
        integer                                         :: PlantType, JulDay, JulDay_Old
        integer                                         :: Op
        real                                            :: HarvestKillPlantHU 
        real                                            :: KillPlantHU
        real, dimension (:), pointer                    :: HarvestHU
        logical                                         :: Dormant
        !SandBoxTest-----------------------------------------------------------
!        allocate(Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestPlantHU(1))
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestKillJulianDay = -99.
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestKillPlantHU   = 1.2
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestJulianDay     = -99.
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestPlantHU(1)    = -99.
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestPlantHU(2)    = -99.
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%KillJulianDay        = -99.
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%KillPlantHU          = -99.
        !Begin-----------------------------------------------------------------


        HUAcc                = Me%HeatUnits%PlantHUAccumulated(i,j)
        HUAcc_Old            = Me%HeatUnits%PlantHUAccumulated_Old(i,j)
        PlantType            = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
        HarvestKillJulianDay = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestKillJulianDay
        HarvestKillPlantHU   = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestKillPlantHU
        HarvestJulianDay     = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestJulianDay
        HarvestHU            =>Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestPlantHU
        KillJulianDay        = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%KillJulianDay
        KillPlantHU          = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%KillPlantHU

        !! check if end of annual growing season

        Dormant = .false.
        if (Me%ComputeOptions%Dormancy) then
            if (Me%IsPlantDormant(i,j)) then
                Dormant = .true.
            endif
        endif

        if (Dormant .and. HUAcc .gt. 0.75) then
            select case (PlantType)
                case (1,4,5)
                    Me%HarvestKillOccurred(i,j) = .true.
            end select
        endif
        
        call JulianDay(Me%ExternalVar%Now, JulDay)
        JulDay_Old = Me%ExternalVar%JulianDay_Old
        
            
        !! harvest and kill operation
        if (HarvestKillJulianDay .gt. 0.0) then
            if(JulDay .ge. HarvestKillJulianDay .and. JulDay_Old .lt. HarvestKillJulianDay) then
            
                Me%HarvestKillOccurred(i,j) = .true.
            
            endif
        elseif (HarvestKillPlantHU .gt. 0.0) then
            if(HUAcc .ge. HarvestKillPlantHU .and. HUAcc_Old .lt. HarvestKillPlantHU) then
            
                Me%HarvestKillOccurred(i,j) = .true.
            
            endif
        endif

        !! harvest operation (no kill)
        if (.not. Me%HarvestFinished(i,j)) then
            
            Op = Me%HarvestOperations(i,j)
            
            if (HarvestJulianDay .gt. 0.0) then
                if(JulDay .ge. HarvestJulianDay .and. JulDay_Old .lt. HarvestJulianDay) then
            
                    Me%HarvestOnlyOccurred(i,j) = .true.
            
                endif
            elseif (HarvestHU(Op) .gt. 0.0) then
                if(HUAcc .ge. HarvestHU(Op) .and. HUAcc_Old .lt. HarvestHU(Op)) then
            
                    Me%HarvestOnlyOccurred(i,j) = .true.
                
                    Me%HarvestOperations(i,j) = Me%HarvestOperations(i,j) + 1
                    if (Me%HarvestOperations(i,j) .gt. size(HarvestHU)) then
                        Me%HarvestFinished(i,j) = .true. 
                    endif

                endif
            else
                Me%HarvestFinished(i,j) = .true. 
            endif
        endif


        !! kill operation
        if (KillJulianDay .gt. 0.0) then
            if(JulDay .ge. KillJulianDay .and. JulDay_Old .lt. KillJulianDay) then
            
                Me%KillOccurred(i,j) = .true.

            endif
        elseif (KillPlantHU .gt. 0.0) then
            if(HUAcc .ge. KillPlantHU .and. HUAcc .lt. KillPlantHU) then
            
                Me%KillOccurred(i,j) = .true.
            endif
        endif
        
!        deallocate(Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestPlantHU)

   
    end subroutine CheckPlantManagement
    !--------------------------------------------------------------------------
    
    subroutine ManagementFluxes

        !Arguments-------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
                
        !Harvest, kill operations
        call PlantManagementFluxes
        
        !Fertilization, tillage
        call SoilManagementFluxes


    end subroutine ManagementFluxes

    !--------------------------------------------------------------------------

    subroutine PlantManagementFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        character (Len = StringLength)                  :: WarningString
        !SandBoxTest-----------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j)) then  

                if (Me%HarvestKillOccurred(i,j)) then
                    call HarvestKillOperation(i,j)
                    WarningString = 'Kill'
                    call UpdatePlantGrowingStage(i,j, WarningString)
                endif

                if (Me%HarvestOnlyOccurred(i,j)) then
                    call HarvestOperation(i,j)
                    WarningString = 'Harves'
                    call UpdatePlantGrowingStage(i,j, WarningString)
                endif

                if (Me%KillOccurred(i,j)) then
                    call KillOperation(i,j)
                    WarningString = 'Kill'
                    call UpdatePlantGrowingStage(i,j, WarningString)
                endif

            endif
        
        enddo do2
        enddo do1

   
    end subroutine PlantManagementFluxes
    !--------------------------------------------------------------------------

    subroutine HarvestOperation(i,j)
        
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: OptimalHarvestIndex, HarvestIndex, HUAcc
        real                                            :: MinimumHarvestIndex, ActualHarvestIndex
        real                                            :: PotTP, ActualTP, WaterDeficiencyFactor
        real                                            :: TotalPlantBiomass, AerialBiomass, Yeld, Residue
        real                                            :: NitrogenYeld, PhosphorusYeld
        real                                            :: NitrogenFractionInYeld, PhosphorusFractionInYeld
        real                                            :: TotalPlantNitrogen, TotalPlantPhosphorus
        real                                            :: NitrogenClipping, PhosphorusClipping
        real                                            :: Clip, HarvestEfficiency
        real                                            :: BiomassHarvestedFraction
        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%OptimalHarvestIndex      = 0.40
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MinimumHarvestIndex      = 0.20
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%NitrogenFractionInYeld   = 0.0250
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PhosphorusFractionInYeld = 0.0022
!        Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestEfficiency    = 1.0
        !Begin-----------------------------------------------------------------



!HARVEST INDEX OVERRIDE not yet implemented
!
!      hiad1 = 0.
!      if (hi_targ(nro(j),icr(j),j) > 0.) then
!        hiad1 = hi_targ(nro(j),icr(j),j)
!      else

        HUAcc                = Me%HeatUnits%PlantHUAccumulated(i,j)
        OptimalHarvestIndex  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%OptimalHarvestIndex
        MinimumHarvestIndex  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MinimumHarvestIndex
        
        HarvestIndex         = OptimalHarvestIndex * 100. * HUAcc / (100. * HUAcc + exp(11.1 - 10. * HUAcc))
        
        ! mm       =   m/s * mm/m * s
        PotTP      = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT * 1000.0
        !  mm      =   m/s * mm/m * s
        ActualTP   = Me%Fluxes%WaterUptake(i,j) * Me%ComputeOptions%VegetationDT * 1000. 

        !!Actualize Harvest Index according to water stress.
        if (PotTP .lt. 10.) then
            WaterDeficiencyFactor = 100.
        else
            WaterDeficiencyFactor = 100. * ActualTP / PotTP
        endif
        
        ActualHarvestIndex = (HarvestIndex - MinimumHarvestIndex) *                                                           &
                              (WaterDeficiencyFactor /  (WaterDeficiencyFactor + exp(6.13 - 0.883 * WaterDeficiencyFactor)))  &
                              + MinimumHarvestIndex
        
        if (ActualHarvestIndex .gt. OptimalHarvestIndex) then
            ActualHarvestIndex = OptimalHarvestIndex
        endif

        TotalPlantBiomass = Me%StateVariables%TotalPlantBiomass%Field(i,j)
        AerialBiomass     = TotalPlantBiomass - Me%StateVariables%RootBiomass%Field(i,j)
        Yeld              = 0.0
        Residue           = 0.0

        !! check if yield is from above or below ground. 
        !Yeld (profitable part of plant) is removed from watershed. The remainder plant is added as residue to soil
        if (OptimalHarvestIndex .gt. 1.001) then
            Yeld    = TotalPlantBiomass * (1. - 1. / (1. + ActualHarvestIndex))
            Residue = TotalPlantBiomass / (1. + ActualHarvestIndex)
        else
            Yeld    = AerialBiomass * ActualHarvestIndex
            Residue = AerialBiomass * (1. - ActualHarvestIndex)
        endif

        if (Yeld .lt. 0.0) then
            Yeld = 0.0
        endif
        if (Residue .lt. 0.0) then
            Residue = 0.0
        endif
        
        !! determine clippings (biomass left behind) and update yield
        HarvestEfficiency = Me%VegetationTypes(Me%VegetationID(i,j))%ManagementDatabase%HarvestEfficiency
        Clip = 0.0
        Clip = Yeld * (1. - HarvestEfficiency)
        Yeld = Yeld *  HarvestEfficiency
        
        if (Yeld .lt. 0.0) then
            Yeld = 0.0
        endif
        if (Clip .lt. 0.0) then
            Clip = 0.0
        endif

!! update residue on soil surface
!      sol_rsd(1,j) = resnew + sol_rsd(1,j)
!      sol_rsd(1,j) = Max(sol_rsd(1,j),0.)



!      if (hi_ovr(nro(j),ncut(j),j) > 0.) then
!        !! calculate nutrients removed with yield
!        yieldn = 0.
!        yieldp = 0.
!        yieldn = yield * pltfr_n(j)
!        yieldp = yield * pltfr_p(j)
!        yieldn = Min(yieldn, 0.9 * plantn(j))
!        yieldp = Min(yieldp, 0.9 * plantp(j))
!        !! calculate nutrients removed with clippings
!        clipn = 0.
!        clipp = 0.
!        clipn = clip * pltfr_n(j)
!        clipp = clip * pltfr_p(j)
!        clipn = Min(clipn,plantn(j)-yieldn)
!        clipp = Min(clipp,plantp(j)-yieldp)
!      else


        !!Biomass to soil. The fraction left by clippings
        Me%Fluxes%ToSoil%ManagementBiomassToSoil(i,j) = Clip
        
        !StateVariables Fluxes
        Me%Fluxes%BiomassRemovedInHarvest(i,j)    = Yeld + Clip
        
        TotalPlantBiomass        = Me%StateVariables%TotalPlantBiomass%Field(i,j)
        
        !For HU update
        BiomassHarvestedFraction = (Yeld + Clip) / TotalPlantBiomass
        Me%Fluxes%BiomassHarvestedFraction = BiomassHarvestedFraction

        

        if (Me%ComputeOptions%ModelNitrogen) then

            TotalPlantNitrogen       = Me%StateVariables%TotalPlantNitrogen%Field(i,j)
            NitrogenFractionInYeld   = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%NitrogenFractionInYeld
            
            !Nutrient yeld
            NitrogenYeld   = 0.0
            NitrogenYeld   = Yeld * NitrogenFractionInYeld
            NitrogenYeld   = min(NitrogenYeld,   0.9 * TotalPlantNitrogen)
            NitrogenYeld   = max(NitrogenYeld, 0.0)

            !Nutrient Clipping
            NitrogenClipping   = 0.0
            NitrogenClipping   = Clip * NitrogenFractionInYeld
            NitrogenClipping   = min(NitrogenClipping, TotalPlantNitrogen - NitrogenYeld)
            NitrogenClipping   = max(NitrogenClipping, 0.0)
            
            !!Nitrogen to soil. The fraction left by clippings
            Me%Fluxes%ToSoil%ManagementNitrogenToSoil(i,j) = NitrogenClipping
 
            !StateVariables Fluxes
            Me%Fluxes%NitrogenRemovedInHarvest(i,j)   = NitrogenYeld + NitrogenClipping
        
        endif

        if (Me%ComputeOptions%ModelPhosphorus) then
        
            TotalPlantPhosphorus     = Me%StateVariables%TotalPlantPhosphorus%Field(i,j)
            PhosphorusFractionInYeld = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PhosphorusFractionInYeld
            
            !Nutrient yeld
            PhosphorusYeld = 0.0
            PhosphorusYeld = Yeld * PhosphorusFractionInYeld
            PhosphorusYeld = min(PhosphorusYeld, 0.9 * TotalPlantPhosphorus)
            PhosphorusYeld = max(PhosphorusYeld, 0.0)

            !Nutrient Clipping
            PhosphorusClipping = 0.0
            PhosphorusClipping = Clip * PhosphorusFractionInYeld
            PhosphorusClipping = min(PhosphorusClipping, TotalPlantPhosphorus - PhosphorusYeld)
            PhosphorusClipping = max(PhosphorusClipping, 0.0)
            
            !!Phosphorus to soil. The fraction left by clippings
            Me%Fluxes%ToSoil%ManagementPhosphorusToSoil(i,j) = PhosphorusClipping

            !StateVariables Fluxes
            Me%Fluxes%PhosphorusRemovedInHarvest(i,j) = PhosphorusYeld + PhosphorusClipping
        
        endif





!        Me%HarvestOnlyOccurred(i,j) = .true.
!
!! adjust foliar pesticide for plant removal
!      if (hrupest(j) == 1) then
!        do k = 1, npmx
!          !! calculate amount of pesticide removed with yield and clippings
!          yldpst = 0.
!          clippst = 0.
!          if (hvsti(idplt(nro(j),icr(j),j)) > 1.001) then
!            yldpst = plt_pst(k,j)
!            plt_pst(k,j) = 0.
!          else
!            yldpst = hiad1 * plt_pst(k,j)
!            plt_pst(k,j) = plt_pst(k,j) - yldpst
!            if (plt_pst(k,j) < 0.) plt_pst(k,j) = 0.
!          endif
!          clippst = yldpst * (1. - harveff(nro(j),ncut(j),j))
!          if (clippst < 0.) clippst = 0.
!          !! add pesticide in clippings to soil surface
!          sol_pst(k,j,1) = sol_pst(k,j,1) + clippst
!        end do   
!      end if
      

!! calculate modifier for autofertilization target nitrogen content
!      tnyld(nro(j),icr(j),j) = 0.
!      tnyld(nro(j),icr(j),j) = (1. - rwt(j)) * bio_ms(j) * pltfr_n(j) * &
!     &                                                       auto_eff(j)
!      if (icr(j) > 1) then
!        tnyld(nro(j),icr(j)-1,j) = tnyld(nro(j),icr(j),j)
!      else
!        tnyld(nro(j),icr(j)+1,j) = tnyld(nro(j),icr(j),j)
!      end if

!! summary calculations
!      if (curyr > nyskip) then
!        wshd_yldn = wshd_yldn + yieldn * hru_dafr(j)
!        wshd_yldp = wshd_yldp + yieldp * hru_dafr(j)
!        yldkg(nro(j),icr(j),j) = yldkg(nro(j),icr(j),j) + yield + clip
!        yldanu(j) = yldanu(j) + (yield + clip) / 1000.
!
!       ! select case (idc(idplt(nro(j),icr(j),j)))
!       !   case (3, 6, 7)
!       !     bio_hv(nro(j),icr(j),j) = (yield + clip) + bio_hv(nro(j),icr(j),j)
!       !     bio_yrms(j) = bio_yrms(j) + (yield + clip) / 1000.
!       !   case default
!           bio_hv(nro(j),icr(j),j) = bio_ms(j) + bio_hv(nro(j),icr(j),j)
!            bio_yrms(j) = bio_yrms(j) + bio_ms(j) / 1000.
!       ! end select
!      endif


!! reset leaf area index and fraction of growing season
!      xx = 0.
!      xx = bio_ms(j)
!      if (xx > 0.001) then
!        bio_ms(j) = bio_ms(j) - yield - clip
!        if (bio_ms(j) < 0.) bio_ms(j) = 0.
!        laiday(j) = laiday(j) * bio_ms(j) / xx
!        phuacc(j) = phuacc(j) * bio_ms(j) / xx
!        rwt(j) = rwt(j) * xx / bio_ms(j)
!      else
!        bio_ms(j) = 0.
!        laiday(j) = 0.
!        phuacc(j) = 0.
!      endif

!! increment harvest sequence number
!      ncut(j) = ncut(j) + 1


    end subroutine HarvestOperation

    !--------------------------------------------------------------------------

    subroutine HarvestKillOperation(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: OptimalHarvestIndex, HarvestIndex, HUAcc
        real                                            :: MinimumHarvestIndex, ActualHarvestIndex
        real                                            :: PotTP, ActualTP, WaterDeficiencyFactor
        real                                            :: TotalPlantBiomass, AerialBiomass, Yeld, Residue
        real                                            :: NitrogenYeld, PhosphorusYeld
        real                                            :: NitrogenFractionInYeld, PhosphorusFractionInYeld
        real                                            :: TotalPlantNitrogen, TotalPlantPhosphorus
        real                                            :: NitrogenToSoil, PhosphorusToSoil
        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%OptimalHarvestIndex      = 0.40
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MinimumHarvestIndex      = 0.20
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%NitrogenFractionInYeld   = 0.0250
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PhosphorusFractionInYeld = 0.0022
        !Begin-----------------------------------------------------------------

!HARVEST INDEX OVERRIDE not yet implemented
!
!      hiad1 = 0.
!      if (hi_targ(nro(j),icr(j),j) > 0.) then
!        hiad1 = hi_targ(nro(j),icr(j),j)
!      else

        HUAcc                = Me%HeatUnits%PlantHUAccumulated(i,j)
        OptimalHarvestIndex  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%OptimalHarvestIndex
        MinimumHarvestIndex  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MinimumHarvestIndex
        
        HarvestIndex         = OptimalHarvestIndex * 100. * HUAcc / (100. * HUAcc + exp(11.1 - 10. * HUAcc))
        
        ! mm       =   m/s * mm/m * s
        PotTP      = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT * 1000.0
        !  mm      =   m/s * mm/m * s
        ActualTP   = Me%Fluxes%WaterUptake(i,j) * Me%ComputeOptions%VegetationDT * 1000. 

        !!Actualize Harvest Index according to water stress.
        if (PotTP .lt. 10.) then
            WaterDeficiencyFactor = 100.
        else
            WaterDeficiencyFactor = 100. * ActualTP / PotTP
        endif
        
        ActualHarvestIndex = (HarvestIndex - MinimumHarvestIndex) *                                                           &
                              (WaterDeficiencyFactor /  (WaterDeficiencyFactor + exp(6.13 - 0.883 * WaterDeficiencyFactor)))  &
                              + MinimumHarvestIndex
        
        if (ActualHarvestIndex .gt. OptimalHarvestIndex) then
            ActualHarvestIndex = OptimalHarvestIndex
        endif

        TotalPlantBiomass = Me%StateVariables%TotalPlantBiomass%Field(i,j)
        AerialBiomass     = TotalPlantBiomass - Me%StateVariables%RootBiomass%Field(i,j)
        Yeld              = 0.0
        Residue           = 0.0

        !! check if yield is from above or below ground. 
        !Yeld (profitable part of plant) is removed from watershed. The remainder plant is added as residue to soil
        if (OptimalHarvestIndex .gt. 1.001) then
            Yeld    = TotalPlantBiomass * (1. - 1. / (1. + ActualHarvestIndex))
            Residue = TotalPlantBiomass / (1. + ActualHarvestIndex)
        else
            Yeld    = AerialBiomass * ActualHarvestIndex
            Residue = AerialBiomass * (1. - ActualHarvestIndex)
        endif

        if (Yeld .lt. 0.0) then
            Yeld = 0.0
        endif
        if (Residue .lt. 0.0) then
            Residue = 0.0
        endif
        
!! update residue on soil surface
!      sol_rsd(1,j) = resnew + sol_rsd(1,j)
!      sol_rsd(1,j) = Max(sol_rsd(1,j),0.)

        
        !!Biomass to soil. The fraction not removed by yeld because plant will die
        Me%Fluxes%ToSoil%ManagementBiomassToSoil(i,j) = Residue
        !The remaining biomass in soil - Not accounted in SWAT because N and P to soil come from all plant
        Me%Fluxes%ToSoil%ManagementRootBiomassLeftInSoil(i,j) = Me%StateVariables%RootBiomass%Field(i,j)

        if (Me%ComputeOptions%ModelNitrogen) then        
            
            NitrogenYeld   = 0.0
            NitrogenFractionInYeld   = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%NitrogenFractionInYeld            
            TotalPlantNitrogen       = Me%StateVariables%TotalPlantNitrogen%Field(i,j)
            
            NitrogenYeld   = Yeld * NitrogenFractionInYeld
            NitrogenYeld   = min (NitrogenYeld,   0.9 * TotalPlantNitrogen)

            !!Nitrogen to soil. The fraction not removed by yeld because plant will die
            if (TotalPlantNitrogen .gt. NitrogenYeld) then
                NitrogenToSoil = TotalPlantNitrogen - NitrogenYeld
            else
                NitrogenToSoil = 0.0
            endif

            Me%Fluxes%ToSoil%ManagementNitrogenToSoil(i,j) = NitrogenToSoil

        endif

        if (Me%ComputeOptions%ModelPhosphorus) then        

            PhosphorusYeld = 0.0
            PhosphorusFractionInYeld = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PhosphorusFractionInYeld
            TotalPlantPhosphorus     = Me%StateVariables%TotalPlantPhosphorus%Field(i,j)

            PhosphorusYeld = Yeld * PhosphorusFractionInYeld
            PhosphorusYeld = min (PhosphorusYeld, 0.9 * TotalPlantPhosphorus)

            !!Phosphorus to soil. The fraction not removed by yeld because plant will die
            if (TotalPlantPhosphorus .gt. PhosphorusYeld) then
                PhosphorusToSoil = TotalPlantPhosphorus - PhosphorusYeld
            else
                PhosphorusToSoil = 0.0
            endif
            
            Me%Fluxes%ToSoil%ManagementPhosphorusToSoil(i,j) = PhosphorusToSoil

        endif

        !Biomass or nutrients fluxes for plant are not computed because plant dies.
        !Instead a warning variable is constructed and when state variables are updated, 
        !model knows that it needs to kill plant - biomass to zero.

        
        Me%Growth%WaterStress(i,j)       = 1.0
        Me%Growth%TemperatureStress(i,j) = 1.0
        Me%Growth%NitrogenStress(i,j)    = 1.0
        Me%Growth%PhosphorusStress(i,j)  = 1.0
        Me%Growth%TreeCurrentYear(i,j)   = 0


!! adjust foliar pesticide for plant removal
!      if (hrupest(j) == 1) then
!        do k = 1, npmx
!          !! calculate amount of pesticide removed with yield
!          yldpst = 0.
!          if (hvsti(idplt(nro(j),icr(j),j)) > 1.001) then
!            yldpst = plt_pst(k,j)
!            plt_pst(k,j) = 0.
!          else
!            yldpst = hiad1 * plt_pst(k,j)
!            plt_pst(k,j) = plt_pst(k,j) - yldpst
!            if (plt_pst(k,j) < 0.) plt_pst(k,j) = 0.
!          endif
!          !! add pesticide in residue to soil surface
!          sol_pst(k,j,1) = sol_pst(k,j,1) + plt_pst(k,j)
!          plt_pst(k,j) = 0.
!        end do
!      end if

!! calculate modifier for autofertilization target nitrogen content
!      tnyld(nro(j),icr(j),j) = 0.
!      tnyld(nro(j),icr(j),j) = (1. - rwt(j)) * bio_ms(j) * pltfr_n(j) * &
!     &                                                       auto_eff(j)
!      if (icr(j) > 1) then
!        tnyld(nro(j),icr(j)-1,j) = tnyld(nro(j),icr(j),j)
!      else
!        tnyld(nro(j),icr(j)+1,j) = tnyld(nro(j),icr(j),j)
!      end if


!! reset variables
!      igro(j) = 0
!      idorm(j) = 0
!      bio_ms(j) = 0.
!      plantn(j) = 0.
!      plantp(j) = 0.
!      strsw(j) = 1.
!      laiday(j) = 0.
!      hvstiadj(j) = 0.
!      phuacc(j) = 0.



    end subroutine HarvestKillOperation

    !--------------------------------------------------------------------------

    subroutine KillOperation(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: AerialBiomass, Residue
        !SandBoxTest-----------------------------------------------------------


!      if (curyr > nyskip) then
!        ncrops(nro(j),icr(j),j) = ncrops(nro(j),icr(j),j) + 1
!      endif

        AerialBiomass     = Me%StateVariables%TotalPlantBiomass%Field(i,j) - Me%StateVariables%RootBiomass%Field(i,j)
        Residue           = 0.0

        !Biomass to soil
        Residue = AerialBiomass
        Me%Fluxes%ToSoil%ManagementBiomassToSoil(i,j) = Residue
        !The remaining biomass in soil - Not accounted in SWAT because N and P to soil come from all plant
        Me%Fluxes%ToSoil%ManagementRootBiomassLeftInSoil(i,j) = Me%StateVariables%RootBiomass%Field(i,j)

!      sol_rsd(1,j) = sol_rsd(1,j) + resnew
!      sol_rsd(1,j) = Max(sol_rsd(1,j),0.)

        if (Me%ComputeOptions%ModelNitrogen) then
            !Nitrogen to soil
            Me%Fluxes%ToSoil%ManagementNitrogenToSoil(i,j) = Me%StateVariables%TotalPlantNitrogen%Field(i,j)
        endif
        
        if (Me%ComputeOptions%ModelPhosphorus) then
            !Phosphorus to soil
            Me%Fluxes%ToSoil%ManagementPhosphorusToSoil(i,j) = Me%StateVariables%TotalPlantPhosphorus%Field(i,j)
        endif

        !Biomass or nutrients fluxes for plant are not computed because plant dies.
        !Instead a warning variable is constructed and when state variables are updated, 
        !model knows that it needs to kill plant.


        Me%Growth%WaterStress(i,j)       = 1.0
        Me%Growth%TemperatureStress(i,j) = 1.0
        Me%Growth%NitrogenStress(i,j)    = 1.0
        Me%Growth%PhosphorusStress(i,j)  = 1.0
        Me%Growth%TreeCurrentYear(i,j)   = 0


!      if (hrupest(j) == 1) then
!        do k = 1, npmx
!          sol_pst(k,j,1) = sol_pst(k,j,1) + plt_pst(k,j)
!          plt_pst(k,j) = 0.
!        end do
!      end if


!! reset variables
!      igro(j) = 0
!      idorm(j) = 0
!      bio_ms(j) = 0.
!      plantn(j) = 0.
!      plantp(j) = 0.
!      strsw(j) = 1.
!      laiday(j) = 0.
!      hvstiadj(j) = 0.
!      phuacc(j) = 0.
!
    end subroutine KillOperation

    !--------------------------------------------------------------------------

    subroutine SoilManagementFluxes
        
        !Arguments-------------------------------------------------------------
!        integer, intent(IN)                             :: i, j


    end subroutine SoilManagementFluxes

    !--------------------------------------------------------------------------

    subroutine UpdatePlantGrowingStage (i,j, WarningString )
       
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        character (Len = StringLength), intent(IN)      :: WarningString
        !Local-----------------------------------------------------------------
        real                                            :: HUAcc, HarvestedFraction, GrazedFraction
        !Begin-----------------------------------------------------------------
        
        if (WarningString .eq. 'Kill') then
        
            Me%HeatUnits%PlantHUAccumulated_Old(i,j) = Me%HeatUnits%PlantHUAccumulated(i,j)            
            Me%HeatUnits%PlantHUAccumulated    (i,j) = 0.0

        elseif (WarningString .eq. 'Harvest') then

            HUAcc = Me%HeatUnits%PlantHUAccumulated(i,j)

            HarvestedFraction = Me%Fluxes%BiomassHarvestedFraction(i,j)
            HUAcc = HUAcc - (HUAcc * HarvestedFraction)
            Me%HeatUnits%PlantHUAccumulated_Old(i,j)  = Me%HeatUnits%PlantHUAccumulated(i,j) 
            Me%HeatUnits%PlantHUAccumulated(i,j) = max(HUAcc, 0.0)
        
        elseif (WarningString .eq. 'Graze') then
        
            HUAcc = Me%HeatUnits%PlantHUAccumulated(i,j)

            GrazedFraction = Me%Fluxes%BiomassGrazedFraction(i,j)
            HUAcc = HUAcc - (HUAcc * GrazedFraction)
            Me%HeatUnits%PlantHUAccumulated_Old(i,j)  = Me%HeatUnits%PlantHUAccumulated(i,j) 
            Me%HeatUnits%PlantHUAccumulated(i,j) = max(HUAcc, 0.0)
        endif

    end subroutine UpdatePlantGrowingStage

    !--------------------------------------------------------------------------

    subroutine UpdateGlobalPlantProperties

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j, PlantType, TreeYearsToMaturity
        logical                                         :: PlantKilled
        real                                            :: PlantBiomass, BiomassGrowth, BiomassGrazed
        real                                            :: BiomassRemovedInHarvest, BiomassRemovedInDormancy
        real                                            :: PlantNitrogen, NitrogenUptake, NitrogenGrazed
        real                                            :: NitrogenRemovedInHarvest, NitrogenRemovedInDormancy
        real                                            :: PlantPhosphorus, PhosphorusUptake, PhosphorusGrazed
        real                                            :: PhosphorusRemovedInHarvest, PhosphorusRemovedInDormancy
        !Begin------------------------------------------------------------------------------------------------
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j)) then    
                     
                !Grazing fluxes
                BiomassGrazed     = 0.0
                NitrogenGrazed    = 0.0
                PhosphorusGrazed  = 0.0   
                if (Me%ComputeOptions%Grazing) then
                    if (Me%IsPlantBeingGrazed(i,j)) then
                        
                        BiomassGrazed     = Me%Fluxes%BiomassGrazed(i,j)
                    
                        if (Me%ComputeOptions%ModelNitrogen) then
                            NitrogenGrazed    = Me%Fluxes%NitrogenGrazed(i,j)
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            PhosphorusGrazed  = Me%Fluxes%PhosphorusGrazed(i,j)
                        endif
                    endif
                endif

                !Harvesting only fluxes
                BiomassRemovedInHarvest     = 0.0
                NitrogenRemovedInHarvest    = 0.0
                PhosphorusRemovedInHarvest  = 0.0
                if (Me%ComputeOptions%Management) then
                    if (Me%HarvestOnlyOccurred(i,j)) then
                        
                        BiomassRemovedInHarvest     = Me%Fluxes%BiomassRemovedInHarvest(i,j)
                    
                        if (Me%ComputeOptions%ModelNitrogen) then
                            NitrogenRemovedInHarvest    = Me%Fluxes%NitrogenRemovedInHarvest(i,j)
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            PhosphorusRemovedInHarvest  = Me%Fluxes%PhosphorusRemovedInHarvest(i,j)
                        endif
                    endif
                endif

                !Dormancy fluxes
                BiomassRemovedInDormancy     = 0.0
                NitrogenRemovedInDormancy    = 0.0
                PhosphorusRemovedInDormancy  = 0.0
                if (Me%ComputeOptions%Dormancy) then
                    if (Me%PlantGoingDormant(i,j)) then
 
                        BiomassRemovedInDormancy     = Me%Fluxes%BiomassRemovedInDormancy(i,j)

                        if (Me%ComputeOptions%ModelNitrogen) then
                            NitrogenRemovedInDormancy    = Me%Fluxes%NitrogenRemovedInDormancy(i,j)
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            PhosphorusRemovedInDormancy  = Me%Fluxes%PhosphorusRemovedInDormancy(i,j)
                        endif
                    endif
                endif


                if (Me%ComputeOptions%ModelPlantBiomass) then
                    PlantBiomass      = Me%StateVariables%TotalPlantBiomass%Field(i,j)
                    BiomassGrowth     = Me%Fluxes%BiomassGrowth(i,j)
                endif

                if (Me%ComputeOptions%ModelNitrogen) then
                    PlantNitrogen     = Me%StateVariables%TotalPlantNitrogen%Field(i,j)
                    NitrogenUptake    = Me%Fluxes%NitrogenUptake(i,j)
                endif

                if (Me%ComputeOptions%ModelPhosphorus) then
                    PlantPhosphorus   = Me%StateVariables%TotalPlantPhosphorus%Field(i,j)
                    PhosphorusUptake  = Me%Fluxes%PhosphorusUptake(i,j)
                endif
                
                PlantKilled = .false.
                if (Me%ComputeOptions%Management) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif

                if (PlantKilled) then
            
                    PlantBiomass    = 0.0
                    PlantNitrogen   = 0.0
                    PlantPhosphorus = 0.0

                else

                    PlantBiomass    = PlantBiomass + BiomassGrowth - BiomassGrazed - BiomassRemovedInHarvest - BiomassRemovedInDormancy
                    
                    PlantNitrogen   = PlantNitrogen + NitrogenUptake - NitrogenGrazed - NitrogenRemovedInHarvest - NitrogenRemovedInDormancy

                    PlantPhosphorus = PlantPhosphorus + PhosphorusUptake - PhosphorusGrazed - PhosphorusRemovedInHarvest - PhosphorusRemovedInDormancy
                    
      
                endif

                PlantType            = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                TreeYearsToMaturity  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeYearsToMaturity

                !For trees, growth is slower, getting maximum biomass not in a single year. As so growth may be limited.
                if (PlantType == 7 .and. TreeYearsToMaturity .gt. 0 ) then
                    
                    PlantBiomass = min (PlantBiomass, Me%Growth%TreeMaximumAnnualBiomass(i,j))

                end if

                if (Me%ComputeOptions%ModelPlantBiomass) then
                    Me%StateVariables%TotalPlantBiomass%Field(i,j)    = PlantBiomass
                endif
                if (Me%ComputeOptions%ModelNitrogen) then
                    Me%StateVariables%TotalPlantNitrogen%Field(i,j)   = PlantNitrogen
                endif
                if (Me%ComputeOptions%ModelPhosphorus) then
                    Me%StateVariables%TotalPlantPhosphorus%Field(i,j) = PlantPhosphorus
                endif

            endif
        enddo do2
        enddo do1

    end subroutine UpdateGlobalPlantProperties

    !--------------------------------------------------------------------------

    subroutine UpdateRootProperties

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        logical                                         :: PlantKilled
        real                                            :: MaxRootDepth, RootBiomassFraction
        integer                                         :: VegetationID, PlantType
        !Begin-----------------------------------------------------------------


        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j)) then   

        !SandBox---------------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaximumRootDepth = 1.30
        !----------------------------------------------------------------------- 
                
                PlantKilled = .false.
                if (Me%ComputeOptions%Management) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif


                if (PlantKilled) then
            
                    if(Me%ComputeOptions%ModelRootBiomass) then
                        Me%StateVariables%RootBiomass%Field(i,j) = 0.0
                    endif
                    Me%StateVariables%RootDepth%Field(i,j)   = 0.0

                else

                    !!Root Biomass
                    if(Me%ComputeOptions%ModelRootBiomass) then
                            
                        RootBiomassFraction = 0.4 - 0.2 * Me%HeatUnits%PlantHUAccumulated(i,j)
        
                        if (RootBiomassFraction .lt. 0.0) then
                            RootBiomassFraction = 0.0
                        endif

                        Me%StateVariables%RootBiomass%Field(i,j) = RootBiomassFraction * Me%StateVariables%TotalPlantBiomass%Field(i,j)
                    
                    endif

                    !!Root Depth
                    VegetationID = Me%VegetationID(i,j)
                    PlantType    = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                    MaxRootDepth = Me%VegetationTypes(VegetationID)%GrowthDatabase%MaximumRootDepth
!                    if(MaxRootDepth .gt. Me%ExternalVar%Topography(i,j)) then
!                        MaxRootDepth = Me%ExternalVar%Topography(i,j)
!                    endif

                    select case (PlantType)
                        case (1, 2, 4, 5)
                
                            Me%StateVariables%RootDepth%Field(i,j) = 2.5 * Me%HeatUnits%PlantHUAccumulated (i,j) * MaxRootDepth
                
                            if (Me%StateVariables%RootDepth%Field(i,j) .gt. MaxRootDepth) then
                                Me%StateVariables%RootDepth%Field(i,j) = MaxRootDepth
                            endif

                            if (Me%StateVariables%RootDepth%Field(i,j) .lt. 0.01) then
                                Me%StateVariables%RootDepth%Field(i,j) = 0.01
                            endif
                        case default
                
                            Me%StateVariables%RootDepth%Field(i,j) = MaxRootDepth

                    end select
                
                endif
            
            endif

        enddo do2
        enddo do1 


    end subroutine UpdateRootProperties

    !--------------------------------------------------------------------------

    subroutine UpdateLeafProperties

        !Arguments-------------------------------------------------------------
        !Local----------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        logical                                         :: PlantKilled
        real                                            :: LAI, BiomassGrazedFraction, BiomassHarvestedFraction
        real                                            :: LAIGrowth, LAIDeclineFraction
        real                                            :: LAIMax, DBLAIMax
        integer                                         :: PlantType
        !Begin----------------------------------------------------------
        
        MappingPoints => Me%ExternalVar%MappingPoints2D

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j)) then            

        !SandBoxTest----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%ComputeLeaf = .true.
        !---------------------------------------------------------------------
        
                if (Me%VegetationTypes(Me%VegetationID(i,j))%HasLeaves) then

                    PlantKilled = .false.
                    if (Me%ComputeOptions%Management) then
                        if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                            PlantKilled = .true.
                        endif
                    endif

                    if (PlantKilled) then

                        Me%StateVariables%LeafAreaIndex%Field(i,j) = 0.0

                    else

                        LAI               = Me%StateVariables%LeafAreaIndex%Field(i,j)
                        DBLAIMax          = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%LAIMax    
                        PlantType         = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType

                        BiomassGrazedFraction = 0.0
                        if (Me%ComputeOptions%Grazing) then
                            if (Me%IsPlantBeingGrazed(i,j)) then
                                BiomassGrazedFraction = Me%Fluxes%BiomassGrazedFraction(i,j)
                            endif
                        endif

                        BiomassHarvestedFraction = 0.0
                        if (Me%ComputeOptions%Management) then
                            if (Me%HarvestOnlyOccurred(i,j)) then
                                BiomassHarvestedFraction = Me%Fluxes%BiomassHarvestedFraction(i,j)
                            endif
                        endif


                        if (.not. Me%LAISenescence(i,j)) then
                
                            LAIGrowth = Me%Fluxes%LAIGrowth(i,j)

                            Me%StateVariables%LeafAreaIndex%Field(i,j) = LAI + LAIGrowth                                       &
                                                                    - (LAI * BiomassGrazedFraction)                            &
                                                                    - (LAI * BiomassHarvestedFraction)
                    
                            if(.not. Me%ComputeOptions%ChangeLAISenescence) then

                                Me%LAIBeforeSenescence(i,j) = Me%StateVariables%LeafAreaIndex%Field(i,j)
                            endif

                        elseif (Me%LAISenescence(i,j) .and. LAI .ne. 0.0) then
                
                            LAIDeclineFraction = Me%Fluxes%LAIDeclineFraction(i,j)

                            if(.not. Me%ComputeOptions%ChangeLAISenescence) then

                                !LAI Computed from maximum value reached and not from last (SWAT theory). If grazing or 
                                !harvesting occurr during senescence LAI is increased due to the computation formulation.
                                Me%StateVariables%LeafAreaIndex%Field(i,j) = Me%LAIBeforeSenescence(i,j) * LAIDeclineFraction  &
                                                                        - (LAI * BiomassGrazedFraction)                        &
                                                                        - (LAI * BiomassHarvestedFraction)
                            else
                                !LAI computed from last value. It avoids erorrs when grazing or harvesting occurrs 
                                !during senescence.
                                Me%StateVariables%LeafAreaIndex%Field(i,j) = (LAI * LAIDeclineFraction)                        &
                                                                        - (LAI * BiomassGrazedFraction)                        &
                                                                        - (LAI * BiomassHarvestedFraction)
                            endif
                        endif
            
                        if (PlantType .eq. 7) then
                            LAIMax = Me%Growth%TreeFractionToMaturity(i,j) * DBLAIMax
                        else
                            LAIMax = DBLAIMax
                        endif


                        if (Me%StateVariables%LeafAreaIndex%Field(i,j) .gt. LAIMax) then
                
                            Me%StateVariables%LeafAreaIndex%Field(i,j) = LAIMax
            
                        endif
                        
                        !If in the future a leaf flux is computed explicitly than LAI can not be
                        !taken to zero here. If mass flux .gt. LAI then flux is LAI
                        if (Me%StateVariables%LeafAreaIndex%Field(i,j) .lt. 0.0) then
                
                            Me%StateVariables%LeafAreaIndex%Field(i,j) = 0.0
            
                        endif
            
                    endif

                 else
            
                    Me%StateVariables%LeafAreaIndex%Field(i,j) = 0.0

                endif

            endif

        enddo do2
        enddo do1

    end subroutine UpdateLeafProperties

    !--------------------------------------------------------------------------

    subroutine UpdateStemProperties

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer, dimension(:,:), pointer                   :: MappingPoints
        integer                                         :: i, j
        logical                                         :: PlantKilled
        real                                            :: MaxCanopyHeight, CanopyHeight, NearMaximum
        integer                                         :: PlantType
        real                                            :: BiomassHarvestedFraction
        !Begin-----------------------------------------------------------------
 
        MappingPoints => Me%ExternalVar%MappingPoints2D
       
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (MappingPoints (i, j) == 1  .and. Me%IsPlantGrowing(i,j)) then    

        !SandBoxTest-----------------------------------------------------------
!        Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaxCanopyHeight = 0.9
        !----------------------------------------------------------------------
        

                Me%ChangeCanopyEnabled(i,j) = .false.
                PlantKilled = .false.
                if (Me%ComputeOptions%Management) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif
    
                if (Me%ComputeOptions%Management .and. PlantKilled) then

                    Me%StateVariables%CanopyHeight%Field(i,j) = 0.0

                else

                    PlantType         = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                    MaxCanopyHeight   = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaxCanopyHeight

                    !!Compute new canopy height
                    if (PlantType == 7) then
                
                        Me%StateVariables%CanopyHeight%Field(i,j) = MaxCanopyHeight * Me%Growth%TreeFractionToMaturity(i,j)
            
                    else

                        NearMaximum = MaxCanopyHeight - 0.01
                    
                        !Change canopy was enabled because after reaching its maximum, crop canopy height was not affected  
                        !by harvest(or grazing but grazing can be consider only to remove leafs and not reduce stem height)
                        Me%ChangeCanopyEnabled(i,j) = .false.
                        if (Me%StateVariables%CanopyHeight%Field(i,j) .gt. NearMaximum .and. &
                            Me%ComputeOptions%ChangeCanopyHeight) then
                            Me%ChangeCanopyEnabled(i,j) = .true.
                        endif
            
                        if (Me%ChangeCanopyEnabled(i,j)) then

!                            BiomassGrazedFraction = 0.0
!                            if (Me%ComputeOptions%Grazing) then
!                                if (Me%IsPlantBeingGrazed(i,j)) then
!                                    BiomassGrazedFraction = Me%Fluxes%BiomassGrazedFraction(i,j)
!                                endif
!                            endif

                            BiomassHarvestedFraction = 0.0
                            if (Me%ComputeOptions%Management) then
                                if (Me%HarvestOnlyOccurred(i,j)) then
                                    BiomassHarvestedFraction = Me%Fluxes%BiomassHarvestedFraction(i,j)
                                endif
                            endif
                
                            CanopyHeight = Me%StateVariables%CanopyHeight%Field(i,j)
                
                            Me%StateVariables%CanopyHeight%Field(i,j) = CanopyHeight                                      &
!                                                                    - CanopyHeight * BiomassGrazedFraction         &
                                                                    - CanopyHeight * BiomassHarvestedFraction
                        
                            !If in the future a height reduction is computed explicitly than canopy height can not be
                            !taken to zero here. If mass flux .gt. canopy height then flux is canopy height                                
                            if (Me%StateVariables%CanopyHeight%Field(i,j) .lt. 0.0) then
        
                                Me%StateVariables%CanopyHeight%Field(i,j) = 0.0

                            endif    
                    
                        else                        
                    

                            Me%StateVariables%CanopyHeight%Field(i,j) = MaxCanopyHeight * sqrt(Me%PlantLAIMaxFraction(i,j))
                
                        endif 

                    endif

                endif
    

            endif

        enddo do2
        enddo do1


    end subroutine UpdateStemProperties

    !--------------------------------------------------------------------------

    subroutine InterfaceWithSoil

    end subroutine InterfaceWithSoil

    !--------------------------------------------------------------------------

    
    subroutine Modify_OutPutHDF

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        type(T_Time)                                :: Actual, LastTime, EndTime
        integer                                     :: OutPutNumber
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
             
        !Begin----------------------------------------------------------------

        Actual   = Me%ExternalVar%Now
        EndTime  = Me%EndTime
        LastTime = Me%LastOutPutHDF5
         
        OutPutNumber = Me%OutPut%NextOutput

TNum:   if (OutPutNumber <= Me%OutPut%Number)            then 
TOut:       if (Actual .GE. Me%OutPut%OutTime(OutPutNumber)) then 
                
                call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))

First:          if (LastTime.LT.Actual) then 
                    
                    !Writes Time
                    TimePtr => AuxTime
                    call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR00'

                    call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                         Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR01'
           
                    Me%LastOutPutHDF5 = Actual
       
                endif First

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                      Me%WorkSize%ILB,                           &
                                      Me%WorkSize%IUB,                           &
                                      Me%WorkSize%JLB,                           &
                                      Me%WorkSize%JUB,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR02'

                !Writes the Open Points
                call HDF5WriteData   (Me%ObjHDF5, "//Grid/OpenPoints",              &
                                      "OpenPoints", "-",                            &
                                      Array2D = Me%ExternalVar%OpenPoints2D,        &
                                      OutputNumber = OutPutNumber,                  &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR03'


                PropertyX => Me%FirstProperty
                do while (associated(PropertyX))

                    if (PropertyX%OutputHDF) then
 
                        call HDF5WriteData   (Me%ObjHDF5,                                    &
                                              "/Results/"//trim(PropertyX%ID%Name),          &
                                              trim(PropertyX%ID%Name),                       &
                                              trim(PropertyX%ID%Units),                      &
                                              Array2D = PropertyX%Field,                     &
                                              OutputNumber = OutPutNumber,                   &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR04'

                    endif

                    PropertyX => PropertyX%Next

                enddo

!                if (Me%ComputeOptions%Evolution%ModelSWAT) then
    !                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//"HUAccumulated",          &
    !                                     "HUAccumulated", "-",                              &
    !                                     Array2D      = Me%HeatUnits%PlantHUAccumulated,    &
    !                                     OutputNumber = OutPutNumber,                       &
    !                                     STAT         = STAT_CALL)                      
    !                if (STAT_CALL /= SUCCESS_)                                              &
    !                    stop 'OutPutHDF - ModuleVegetation - ERR05'                 
!                endif

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//"WaterUptake",&
                                     "WaterUptake", "m3/s",                                  &
                                     Array2D      = Me%Fluxes%WaterUptake,                   &
                                     OutputNumber = OutPutNumber,                            &
                                     STAT         = STAT_CALL)                      
                if (STAT_CALL /= SUCCESS_)                                                   &
                    stop 'OutPutHDF - ModuleVegetation - ERR05'                 


                Me%OutPut%NextOutput = OutPutNumber + 1

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR06'
            
            endif  TOut
        endif  TNum

    end subroutine Modify_OutPutHDF

    !--------------------------------------------------------------------------

    subroutine Modify_OutPutTimeSeries

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))

            if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie, Data2D = PropertyX%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR01'
                
            endif
            
            PropertyX => PropertyX%Next

        enddo
        
        call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Fluxes%WaterUptake, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR02'
        
        call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%WaterStress, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR03'

        if (Me%ComputeOptions%Evolution%ModelSWAT) then
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%HeatUnits%PlantHUAccumulated, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR04'

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%HeatUnits%PotentialHUBase, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR05'

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%NitrogenStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR06'

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%PhosphorusStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR07'
            
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%TemperatureStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR08'

        endif

        
        if (Me%ComputeOptions%ModelNitrogen) then

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Fluxes%NitrogenUptake, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR09'
            
        endif
        
        if (Me%ComputeOptions%ModelPhosphorus) then

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Fluxes%PhosphorusUptake, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR10'

        endif

!        call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%ExternalVar%Integration%AverageRadiationDuringDT, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR11'
        

    end subroutine Modify_OutPutTimeSeries



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillVegetation(ObjVegetationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: ObjVegetationID              
        integer, optional, intent(OUT)              :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_, nUsers          
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjVegetationID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mVegetation_,  Me%InstanceID)

            if (nUsers == 0) then

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALMAP_, Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR02'

                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR03'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjGridData)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR04'

                nUsers = DeassociateInstance (mGEOMETRY_, Me%ObjGeometry)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR05'

                nUsers = DeassociateInstance (mATMOSPHERE_, Me%ObjAtmosphere)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR06'

                nUsers = DeassociateInstance (mPOROUSMEDIA_, Me%ObjPorousMedia)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR07'

                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - ModuleVegetation  - ERR08'
                endif
                
                if ((Me%ObjTimeSerie > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_)
                    if (STAT_ .NE. SUCCESS_)                                                 &
                        stop 'KillVegetation - ModuleVegetation - ERR09'
                endif

                call DeAllocateVariables

                !Deallocates Instance
                call DeallocateInstance ()

                ObjVegetationID = 0
                STAT_           = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           
        !------------------------------------------------------------------------

    end subroutine KillVegetation

    !------------------------------------------------------------------------

   
    subroutine DeAllocateVariables

        !Local------------------------------------------------------------------
        type (T_Property), pointer :: PropertyX
        integer                    :: STAT_CALL
        !----------------------------------------------------------------------


        ! Deallocates all the vegetation properties 

        PropertyX => Me%FirstProperty

do1 :   do while(associated(PropertyX))  

            if (associated(PropertyX%Field)) then

                deallocate(PropertyX%Field    )
                nullify   (PropertyX%Field    )
            else
                deallocate(PropertyX%Field    )
                nullify   (PropertyX%Field    )
            endif


            if (PropertyX%ID%ObjFillMatrix /= 0) then
                call KillFillMatrix (PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'DeAllocateVariables - ModuleVegetation - ERR02'
            endif

           
            PropertyX => PropertyX%Next

        end do do1

        !Sets the number of properties equal to the FillValueInt
        Me%PropertiesNumber = FillValueInt

        Nullify   (Me%FirstProperty,Me%LastProperty)
        deallocate(Me%ExternalVar%SoilNitrogenConcentration                 ) 
        deallocate(Me%ExternalVar%SoilPhosphorusConcentration               ) 
        
        deallocate(Me%VegetationID                                          )

        !Soil Fluxes
        deallocate(Me%Fluxes%WaterUptake                                    )  
        deallocate(Me%Fluxes%WaterUptakeLayer                               )  
!        deallocate(Me%Fluxes%FromSoil%WaterUptakeFromSoil                   )  
        
        deallocate(Me%Growth%WaterStress                                    ) 
        deallocate(Me%ExternalVar%Integration%SumPotTP                      )
        deallocate(Me%ExternalVar%Integration%AveragePotTPDuringDT  )

        if (Me%ComputeOptions%Evolution%ModelSWAT) then
            deallocate(Me%IsPlantGrowing                                    )
            deallocate(Me%PlantingOccurred                                  )
            deallocate(Me%HeatUnits%PotentialHUTotal                        ) 
            deallocate(Me%HeatUnits%PotentialHUBase                         ) 
            deallocate(Me%HeatUnits%PlantHUAccumulated                      ) 
            deallocate(Me%HeatUnits%PlantHUAccumulated_Old                  ) 
            deallocate(Me%ExternalVar%Integration%AverageAirTempDuringDT    )  
            deallocate(Me%ExternalVar%Integration%AverageAirHumidityDuringDT) 
            deallocate(Me%ExternalVar%Integration%AverageRadiationDuringDT  )
            deallocate(Me%ExternalVar%Integration%SumTemperature            )  
            deallocate(Me%ExternalVar%Integration%SumHumidity               ) 
            deallocate(Me%ExternalVar%Integration%SumRadiation              )
             
            deallocate(Me%Growth%GlobalStress                               )
            deallocate(Me%Growth%TemperatureStress                          ) 
            deallocate(Me%Growth%NitrogenStress                             ) 
            deallocate(Me%Growth%PhosphorusStress                           ) 

            deallocate(Me%Growth%TreeCurrentYear                            )        
            deallocate(Me%Growth%TreeFractionToMaturity                     )        
            deallocate(Me%Growth%TreeMaximumAnnualBiomass                   )        


            !Soil fluxes
            if(Me%ComputeOptions%ModelNitrogen) then
                deallocate(Me%PlantNitrogenFraction                         ) 
                deallocate(Me%Fluxes%NitrogenUptake                         )  
                deallocate(Me%Fluxes%NitrogenUptakeLayer                    )  
!                deallocate(Me%Fluxes%FromSoil%NitrogenUptakeFromSoil        ) 
                
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                deallocate(Me%PlantPhosphorusFraction                       )
                deallocate(Me%Fluxes%PhosphorusUptake                       ) 
                deallocate(Me%Fluxes%PhosphorusUptakeLayer                  )  
!                deallocate(Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil      ) 

            endif
            
            !Aerial Fluxes
            deallocate(Me%Fluxes%LAIGrowth                                  ) 
            deallocate(Me%Fluxes%LAIDeclineFraction                         ) 
            deallocate(Me%LAISenescence                                     )             
            deallocate(Me%PlantLAIMaxFraction                               ) 
            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                deallocate(Me%LAIBeforeSenescence                           ) 
            endif



        endif

        if (Me%ComputeOptions%Grazing) then
            deallocate(Me%DaysOfGrazing                                     ) 
            deallocate(Me%IsPlantBeingGrazed                                ) 

            deallocate(Me%GrazingFinished                                   )
            deallocate(Me%GrazingOperations                                 )
            
            deallocate(Me%Fluxes%BiomassGrazed                              ) 
            deallocate(Me%Fluxes%BiomassGrazedFraction                      ) 
            deallocate(Me%Fluxes%ToSoil%GrazingBiomassToSoil                ) 

            if (Me%ComputeOptions%ModelNitrogen) then
                deallocate(Me%Fluxes%NitrogenGrazed                         )  
                deallocate(Me%Fluxes%ToSoil%GrazingNitrogenToSoil           ) 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                deallocate(Me%Fluxes%PhosphorusGrazed                       ) 
                deallocate(Me%Fluxes%ToSoil%GrazingPhosphorusToSoil         ) 
            endif

        endif
        
        if (Me%ComputeOptions%Management) then
            deallocate(Me%HarvestOnlyOccurred                               )
            deallocate(Me%HarvestKillOccurred                               )
            deallocate(Me%KillOccurred                                      )
            deallocate(Me%HarvestFinished                                   )
            deallocate(Me%HarvestOperations                                 )
            deallocate(Me%HarvestKillOperations                             ) 
            deallocate(Me%KillOperations                                    ) 
                      
            deallocate(Me%Fluxes%BiomassRemovedInHarvest                    ) 
            deallocate(Me%Fluxes%BiomassHarvestedFraction                   ) 
            deallocate(Me%Fluxes%ToSoil%ManagementBiomassToSoil             ) 
            deallocate(Me%Fluxes%ToSoil%ManagementRootBiomassLeftInSoil     )

            if (Me%ComputeOptions%ModelNitrogen) then            
                deallocate(Me%Fluxes%NitrogenRemovedInHarvest               )  
                deallocate(Me%Fluxes%ToSoil%ManagementNitrogenToSoil        ) 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                deallocate(Me%Fluxes%PhosphorusRemovedInHarvest             ) 
                deallocate(Me%Fluxes%ToSoil%ManagementPhosphorusToSoil      )
            endif

        endif

        if (Me%ComputeOptions%Dormancy) then
            deallocate(Me%IsPlantDormant                                    ) 
            deallocate(Me%PlantGoingDormant                                 )
            
            deallocate(Me%Fluxes%BiomassRemovedInDormancy                   ) 
            deallocate(Me%Fluxes%ToSoil%DormancyBiomassToSoil               ) 
            
            if (Me%ComputeOptions%ModelNitrogen) then
                deallocate(Me%Fluxes%NitrogenRemovedInDormancy              )  
                deallocate(Me%Fluxes%ToSoil%DormancyNitrogenToSoil          ) 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                deallocate(Me%Fluxes%PhosphorusRemovedInDormancy            ) 
                deallocate(Me%Fluxes%ToSoil%DormancyPhosphorusToSoil        ) 
            endif
 
       endif
        
        if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH1 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH2 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH4 )
        endif
        


        if (Me%ComputeOptions%ModelPlantBiomass) then
            deallocate(Me%Fluxes%BiomassGrowth                              ) 
            deallocate(Me%Growth%BiomassGrowthOld                           )  
        
        endif

        if (Me%ComputeOptions%ModelCanopyHeight) then
            deallocate(Me%ChangeCanopyEnabled                               ) 
        endif



    end subroutine DeAllocateVariables  
   
   !----------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Vegetation), pointer          :: AuxObjVegetation
        type (T_Vegetation), pointer          :: PreviousObjVegetation

        !Updates pointers
        if (Me%InstanceID == FirstObjVegetation%InstanceID) then
            FirstObjVegetation => FirstObjVegetation%Next
        else
            PreviousObjVegetation => FirstObjVegetation
            AuxObjVegetation      => FirstObjVegetation%Next
            do while (AuxObjVegetation%InstanceID /= Me%InstanceID)
                PreviousObjVegetation => AuxObjVegetation
                AuxObjVegetation      => AuxObjVegetation%Next
            enddo

            !Now update linked list
            PreviousObjVegetation%Next => AuxObjVegetation%Next

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

    subroutine Ready (ObjVegetation_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjVegetation_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjVegetation_ID > 0) then
            call LocateObjVegetation (ObjVegetation_ID)
            ready_ = VerifyReadLock (mVegetation_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjVegetation (ObjVegetationID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjVegetationID

        !Local-----------------------------------------------------------------

        Me => FirstObjVegetation
        do while (associated (Me))
            if (Me%InstanceID == ObjVegetationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleVegetation - LocateObjVegetation - ERR01'

    end subroutine LocateObjVegetation

    
    !--------------------------------------------------------------------------
    
    subroutine ReadLockExternalVar(ReadAtmosphere)
        
        !Arguments-------------------------------------------------------------
        logical, intent(IN)                     :: ReadAtmosphere       
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        !Begin-----------------------------------------------------------------
        

!        !Gets Time
!        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)              
!        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR010'

        if (ReadAtmosphere) then
            !Gets Air Temperature [ºC]
            call GetAtmosphereProperty  (Me%ObjAtmosphere, Me%ExternalVar%AirTemperature,         &
                                         ID = AirTemperature_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR030'

            call GetAtmosphereProperty  (Me%ObjAtmosphere, Me%ExternalVar%SolarRadiation,         &
                                         ID = SolarRadiation_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR040'

            call GetAtmosphereProperty  (Me%ObjAtmosphere, Me%ExternalVar%RelativeHumidity,       &
                                         ID = RelativeHumidity_, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR041'

        else

            call GetThetaField(Me%ObjPorousMedia, Me%ExternalVar%FieldCapacity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR042'   

            call GetWaterContent(Me%ObjPorousMedia, Me%ExternalVar%SoilWaterContent, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR050'

            call GetHead(Me%ObjPorousMedia, Me%ExternalVar%Head, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR051'

            call GetThetaR(Me%ObjPorousMedia, Me%ExternalVar%ResidualWaterContent, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR052'

            call GetGeometryDistances (Me%ObjGeometry, DWZ = Me%ExternalVar%DWZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR60'

            call GetGridCellArea (Me%ObjHorizontalGrid,                                           & 
                                  GridCellArea = Me%ExternalVar%GridCellArea,                     & 
                                  STAT = STAT_CALL )    
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR080'

            call GetGeometryVolumes(Me%ObjGeometry,                                               &
                                    VolumeZ    = Me%ExternalVar%CellVolume,                       &
                                    STAT       = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                            &
                call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModuleVegetation - ERR0100")

            !Gets a pointer to OpenPoints2D
            call GetOpenPoints2D  (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR110'
            

        endif


    end subroutine ReadLockExternalVar

    !----------------------------------------------------------------------

    subroutine ReadUnlockExternalVar(ReadAtmosphere)
        
        !Arguments-------------------------------------------------------------
        logical, intent(IN)                     :: ReadAtmosphere       
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        if (ReadAtmosphere) then

            call UnGetAtmosphere    (Me%ObjAtmosphere, Me%ExternalVar%AirTemperature, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR010'

            call UnGetAtmosphere    (Me%ObjAtmosphere, Me%ExternalVar%SolarRadiation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR020'

            call UnGetAtmosphere    (Me%ObjAtmosphere, Me%ExternalVar%RelativeHumidity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR021'

        else

            call UnGetPorousMedia    (Me%ObjPorousMedia, Me%ExternalVar%FieldCapacity, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR022'

            call UnGetPorousMedia    (Me%ObjPorousMedia, Me%ExternalVar%SoilWaterContent, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR030'

            call UnGetPorousMedia    (Me%ObjPorousMedia, Me%ExternalVar%Head, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR031'

            call UnGetPorousMedia    (Me%ObjPorousMedia, Me%ExternalVar%ResidualWaterContent, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR032'

            call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%DWZ, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleVegetation - ERR040'

            !GridCellArea
            call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleVegetation - ERR050'

            call UnGetGeometry( Me%ObjGeometry, Me%ExternalVar%CellVolume,   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                                   &
                call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalVar; ModuleVegetation - ERR060")

            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR70'


        endif

    end subroutine ReadUnlockExternalVar

end module ModuleVegetation

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 







