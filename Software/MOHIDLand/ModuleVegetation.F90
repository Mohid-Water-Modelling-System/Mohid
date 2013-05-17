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
! REVISION      : David Brito 2010
! DESCRIPTION   : Module to simulate plant characteristics and growth cycle
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

!Units in vegetation
!   N, P, biomass  : kg/ha
!   Transpiration  : m3/s
!   LAI            : m2leafs/m2area


!DATAFILE
!         [Keyword]                 [Format]  [Units]  [Default]  [Short Description]
! VEGETATION_ID_FILE               : string      -        [-]     !Vegetation distribution grid path
!
! VEGETATION_DT                    : real        s      [86400.]  !Vegetation DT
! INTEGRATION_DT                   : real        s      [ModelDT] !DT to integrate external variables until vegetation is
!                                                                 ! is called (vegetation DT)
! PARAMETERS_FILE                  : string      -         -      !agricultural practices definition
! GROWTH_DATABASE                  : string      -         -      !Growth parameters for each vegetation type - readed in case 
!                                                                  of vegetation growth simulation
! PESTICIDE_DATABASE               : string      -         -      !Readed if growth simulation and PESTICIDE : 1
! FERTILIZER_DATABASE              : string      -         -      !Readed if growth simulation and if FERTILIZATION : 1
! FEDDES_DATABASE                  : string      -         -      !Readed if not usding growth simulation
!
! MODEL_EVOLUTION                  : integer     -         -      !Global properties evolution: 
!                                                                     !1-Read from file
!                                                                     !2-Vegetation growth model
!
! WATER_STRESS                     : 0/1         -         1      !Connects/disconnects water limitation on plant growth
! NITROGEN_STRESS                  : 0/1         -         1      !Connects/disconnects N limitation on plant growth
! PHOSPHORUS_STRESS                : 0/1         -         1      !Connects/disconnects P limitation on plant growth
! TEMPERATURE_STRESS               : 0/1         -         1      !Connects/disconnects temp. limitation on plant growth
! ADJUST_RUE_FOR_CO2               : 0/1         -         1      !Connects/disconnects CO2 limitation on plant growth
! ADJUST_RUE_FOR_VPD               : 0/1         -         1      !Connects/disconnects Vapour Pressure Deficit limitation
!                                                                  plant growth
!
! GRAZING                          : 0/1         -        [0]     !Connects/disconnects grazing
! HARVEST_KILL                     : 0/1         -        [0]     !Connects/disconnects Harvest and/or Kill
! DORMANCY                         : 0/1         -        [0]     !Connects/disconnects dormancy
! FERTILIZATION                    : 0/1         -        [0]     !Connects/disconnects fertilization     
! NUTRIENT_FLUXES_WITH_SOIL        : 0/1         -        [1]     !Connects/disconnects nutrient fluxes with soil
!
! WATER_UPTAKE_METHOD              : integer     -        [1]     !1- according to root profile; 2-SWAT based (exponential
!                                                                  and tresholds)
!   LIMIT_TRANSP_WATER_VEL         : 0/1         -        [0]     !Read if WATER_UPTAKE_METHOD == 1.
!   ROOT_PROFILE                   : integer     -        [1]     !Read if WATER_UPTAKE_METHOD == 1: 
!                                                                   !1-Triangular; 2-Constant; 3-Exponential(SWAT like)
!   WATER_UPTAKE_STRESS_METHOD     : integer     -        [1]     !Read if WATER_UPTAKE_METHOD == 1;1-Feddes;2-VanGenuchten
!   SALINITY_STRESS_METHOD         : integer     -        [1]     !Read if WATER_UPTAKE_METHOD == 1;1-Threshold/Slope;
!                                                                                   2-VanGenuchten (not implemented yet)
!
! NUTRIENT_UPTAKE_METHOD           : integer     -        [2]     !1- uptake is: conc * water uptake; 2- SWAT based 
!                                                                  (independent of water uptake); 3 - NO nutrient uptake
! NUTRIENT_STRESS_METHOD           : integer     -        [2]     !1- effective/optimal; 2- SWAT based
!
! CHANGE_LAI_SENESCENCE            : 0/1         -        [0]     !Changes made to swat code because showed error with 
! CHANGE_CANOPY_HEIGHT             : 0/1         -        [0]       grazing
!
! ATMOSPHERE_OUTPUT                : 0/1         -        [0]     !Output averaged atmosphere properties during dt
! FLUXES_TO_SOIL_OUTPUT            : 0/1         -        [0]     !Output fluxes to soil
!
!
!    
! ATMOSPHERE_CO2                   : real       ppm      [330.]   !Atmosphere CO2 concetrations - should be atmosphere prop
! WATER_UPTAKE_COMPENSATION_FACTOR : real        -        [0.]    !Factor for uptake compensation from lower layers if
!                                                                   !computed layer demand is not met
!                                                                   !If zero there will exist no compensation. If 1. total
!                                                                   !demand no met may come from lower layers
! NITROGEN_DISTRIBUTION_PARAMETER  : real                [-20.]
! PHOSPHORUS_DISTRIBUTION_PARAMETER: real                [-20.]
!
!Potential total HU (yearly HU) -  SUMi=1to12(average monthly temperature in month i * days in month i)
! <begin_TotalPotentialHU>
! INITIALIZATION_METHOD     : CONSTANT
! DEFAULTVALUE              : 5475.
! REMAIN_CONSTANT           : 1
! <end_TotalPotentialHU>

! <beginproperty>
!  See module fillmatrix
! EVOLUTION                        : integer     -         *      !Property evolution: 
!                                                                     !1-Read from file
!                                                                     !2-vegetation growth model
!                                                                     !*Default value is the value set in MODEL_EVOLUTION
!                                                                     !If MODEL_EVOLUTION is set to 1, than property can be 
!                                                                     !only set to 1 also
! <endproperty>

!------------------------------------------------------------------------------------------------------------
! PARAMETERS_FILE - always used
! Arable Land - Trigo
! <beginagriculturalpractice>
! AGRIC_PRACT_ID            : 2                !agriculture practice ID
! NAME                      : Agriculture
! 
! VEGETATION_ID             : 2                !crop ID used in this practice that has correspondence to 
!                                                SWAT crop growth database (see growth database)
! 
! <begintimingparameters>
! !active if growth model used
! PLANTING_JULIANDAY                : -99.      !julian day when planting will occur
! PLANTING_HUBASE                   : 0.15      !Percentage of POTENTIAL YEARLY HU when planting will occur
! MATURITY_HU                       : 1700.     !Total PLANT ACCUMULATED HU when reaching maturity
! <endtimingparameters>
! 
! <beginharvestkillparameters>
! !active if growth model used and in data file HARVEST_KILL : 1
! HARVESTKILL_JULIANDAY             : -99.      !julian day when harvestkill operation occur
! HARVESTKILL_PLANTHU               : 1.2       !Percentage of PLANT ACCUMULATED HU when harvestkill operation occur
! HARVEST_JULIANDAY                 : -99.      !julian day when harvest operation occur
! HARVEST_PLANTHU                   : -99.      !Percentage of PLANT ACCUMULATED HU when harvest operation occur
! HARVEST_EFFICIENCY                : 1.0       !Efficiency for harvest operation (residue if lower than 1)
! KILL_JULIANDAY                    : -99.      !julian day when harvestkill operation occur
! KILL_PLANTHU                      : -99.      !Percentage of PLANT ACCUMULATED HU when kill operation occur
! <endharvestkillparameters>
! 
! <begingrazeparameters>
! !Graze active if growth model used and in data file GRAZING : 1
! GRAZING_START_JULIANDAY           : -99.      !julian day when grazing will occur
! GRAZING_START_PLANTHU             : 0.5       !Percentage of POTENTIAL YEARLY HU when grazing will occur
! GRAZING_DAYS                      : 10        !Days of grazing (continuous)
! MINIMUM_BIOMASS_FOR_GRAZING       : 10.       !minimum biomass (kg/ha) for grazing
! GRAZING_BIOMASS                   : 70.       !grazed biomass (kh/ha.day)
! TRAMPLING_BIOMASS                 : 30.       !biomass not eaten but removed from plant and moved to soil, 
!                                                related to grazing efficiency (kg/ha.day)
! <endgrazeparameters>
! 
! <beginfertilizationparameters>
! !Autofertilization - active if growth model used and in data file FERTILIZATION : 1 and (in data file NITROGEN : 1 
!                                              and NITROGEN_TRESHOLD > 0) or (PHOSPHORUS : 1 and PHOSPHORUS_TRESHOLD > 0)
! FERTILIZER_ID                     : 1      !Fertilizer used in autofertilization (see fertilizer database)
! NITROGEN_TRESHOLD                 : 0.93   !Percentage of stress below which autofertilization starts
! NITROGEN_APPLICATION_MAX          : 50.    !Maximum amount of fertilizer in one application (kg/ha)
! NITROGEN_ANNUAL_MAX               : 300.   !Maximum amount of fertilizer in one year (kg/ha)
! EXPLICIT_PHOSPHORUS               : 1      !1- explicit add phosphorus if needed; 0-add phosphorus if nitrogen needed 
!                                              (SWAT method)
! PHOSPHORUS_TRESHOLD               : 0.93   !only read if EXPLICIT_PHOSPHORUS : 1
! PHOSPHORUS_APPLICATION_MAX        : 10.    !only read if EXPLICIT_PHOSPHORUS : 1
! PHOSPHORUS_ANNUAL_MAX             : 60.    !only read if EXPLICIT_PHOSPHORUS : 1
! 
! !Scheduled Fertilization - not programmed - it should look like pesticide
! FERTILIZATION_JULIANDAY           :-99.
! FERTILIZATION_HU                  :-99.
! <endfertilizationparameters>
! 
! <beginpesticideparameters>
! !Active if growth model used and in data file PESTICIDE : 1
! 
!  <<beginpesticideapp>>
!  PESTICIDE_ID                : 1           !Pesticide used in this application (see pesticide database)
!  PESTICIDE_APPLICATION_JDAY  : -99.        !julian day when pesticide application will occur
!  PESTICIDE_APPLICATION_HU    : 0.10        !Percentage of POTENTIAL YEARLY HU when pesticide application will occur
!  PESTICIDE_APPLICATION_KG_HA : 1.          !Amount of pesticide applied (kg/ha)
!  <<endpesticideapp>>
!
! <endpesticideparameters>
!
! <endagriculturalpractice>


!---------------------------------------------------------------------------------------------------------------
!  GROWTH_DATABASE - used if using growth model
!   <begingrowthdatabase>
!    VEGETATION_ID                     : 1
!    NAME                              : Forest    
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
!    TREE_YEARSTOMATURITY              : -99.
!    TREE_MAXIMUMBIOMASS               : -99.
!    BIOMASS_FRAC_REMOVED_DORMANCY     : 0.30    
!    LAI_MIN_DORMANCY                  : 0.75
!   <endgrowthdatabase>
!      

!--------------------------------------------------------------------------------------------------------------
! PESTICIDE_DATABASE - used if using growth model and PESTICIDE :1
! <beginPesticide>
! PESTICIDE_ID : 1
! PESTICIDE_NAME : generic pesticide 1
! <endPesticide>

!---------------------------------------------------------------------------------------------------------------
! FEDDES_DATABASE - Used if not using growth model
! <beginfeddesdatabase>
! VEGETATION_ID             : 2                !crop ID used in this practice that has correspondence to 
!                                                SWAT crop growth database (see growth database)
! 
! FEDDES_TYPE               : 1 - Normal 2 - With ponint 3 variable consoant the transpiration
! FEDDES_H1                 : -0.1             !higher head for transpiration (saturation and oxygen loss)
! FEDDES_H2                 : -0.25            !1st optimal head for transpiration
! FEDDES_H3                 : -2.0             !2nd optimal head for transpiration
! FEDDES_H3H                : 
! FEDDES_H3L                :
! FEDDES_H4                 : -80.0            !lower head  for transpiration (wilting)
! FEDDES_R2L                :
! FEDDES_R2H                :
! 
!
! USE_SALINITY              : 1 (yes) - 0 (no)
! SALINITY_STRESS_THRESHOLD : 
! SALINITY_STRESS_SLOPE     :  
! STRESS_INTERACTION        : 1 - Additive; 2 - Multiplicative; 3 - Min of all  
! <endfeddesdatabase>            

Module ModuleVegetation

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleHDF5,           only : ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,             &
                                     HDF5WriteData, HDF5FlushMemory, HDF5ReadData, KillHDF5,      &
                                     GetHDF5GroupNumberOfItems
    use ModuleAtmosphere,     only : GetAtmosphereProperty, AtmospherePropertyExists,             &
                                     UnGetAtmosphere
    use ModuleBasinGeometry,  only : GetBasinPoints, UngetBasin
    use ModulePorousMedia,    only : GetWaterContent, GetHead, GetThetaR, GetComputeSoilField,    &
                                     GetThetaField, GetLimitThetaLow, GetUnsatK, UngetPorousMedia
    use ModuleFunctions,      only : SetMatrixValue, ConstructPropertyID, LinearInterpolation,    &
                                     InterpolateValueInTime
    use ModuleHorizontalGrid, only : GetHorizontalGridSize, GetGridCellArea, WriteHorizontalGrid, &
                                     UngetHorizontalGrid, GetGridLatitudeLongitude, GetXYCellZ
!    use ModuleHorizontalMap,  only : GetOpenPoints2D, UngetHorizontalMap
    use ModuleFillMatrix,     only : ConstructFillMatrix, GetDefaultValue, ModifyFillMatrix,      &
                                     KillFillMatrix, GetIfMatrixRemainsConstant
    use ModuleTimeSerie,      only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,               &
                                     StartTimeSerieInput, GetTimeSerieValue,                      &
                                     GetNumberOfTimeSeries, GetTimeSerieLocation,                 &
                                     GetTimeSerieName, TryIgnoreTimeSerie, CorrectsCellsTimeSerie, &
                                     GetTimeSerieDataMatrix, GetTimeSerieDataValues
    use ModuleGridData,       only : ConstructGridData, GetGridData, UngetGridData,               &
                                     KillGridData
    use ModuleGeometry,       only : GetGeometryDistances, GetGeometrySize, GetGeometryKFloor,    &
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
!    private ::              ConstructPropertyID       
    private ::              ConstructPropertyValues
    private ::              ConstructPropertyEvolution   
    private ::              ConstructPropertyOutput  
    private ::          AddProperty
!    private ::      ConstructSoilConcList
!    private ::          ConstructSoilConc
!    private ::              ConstructSoilConcOptions
!    private ::          AddSoilConc
    private ::      CheckOptionsConsistence 
    private ::      ConstructVegetationGrids
    private ::      ConstructVegetationParameters
    private ::          ReadFeddesDatabase
    private ::          ReadTimingParameters
    private ::          ReadGrowthDatabase
    private ::          ReadHarvestKillParameters
    private ::          ReadGrazingParameters
    private ::          ReadFertilizationParameters
    private ::          ReadPesticideParameters
    private ::              GetUniquePesticide
    private ::      ConstructTimeSerie
    private ::      ConstructHDF
    private ::          ConstructHDF5Output 
    private ::      ConstructLog                 

    !Selector
    public  :: GetLeafAreaIndex
    public  :: GetPotLeafAreaIndex
    public  :: GetSpecificLeafStorage
    public  :: GetEVTPCropCoefficient
    public  :: GetVegetationDT
    public  :: GetRootDepth
    public  :: GetRootDepthOld
    public  :: GetCanopyHeight
    public  :: GetTranspiration
    public  :: GetTranspirationBottomLayer
    public  :: GetVegetationOptions
    public  :: GetVegetationGrowing
    public  :: GetNutrientFraction
    public  :: GetCanopyStorageType
!    public  :: GetFeddesH
!    public  :: GetVegetationRootProfile
    public  :: GetVegetationSoilFluxes
    public  :: GetVegetationAerialFluxes
    public  :: SetSoilConcVegetation
    public  :: SetECw
    public  :: UngetVegetation
    public  :: UngetVegetationSoilFluxes
    public  :: UngetVegetationAerialFluxes    

    !Modifier
    public  :: ModifyVegetation
                    !Use model for property evolution
    private ::      CheckPlantState
    private ::      ModifyFluxes
    private ::          PlantRootFluxes
    private ::              WaterUptakeSWAT
    private ::              NitrogenUptakeSWAT
    private ::              NitrogenFixationSWAT
    private ::              PhosphorusUptakeSWAT
    private ::          PlantAerialFluxes
    private ::              BiomassGrowthFromRadiationSWAT
    private ::              LAIGrowthSWAT
    private ::          GrazingFluxes
    private ::          HarvestKillFluxes
    private ::              HarvestOperation
    private ::              HarvestKillOperation
    private ::              KillOperation
    private ::      ModifyModelledProperties
    private ::          UpdateGlobalPlantProperties
    private ::          UpdateRootProperties
    private ::          UpdateStemProperties
    private ::          UpdateLeafProperties
                    !Read property evolution
    private ::      ModifyReadedProperties

    private ::      Modify_OutputHDF
    private ::      Modify_OutputTimeSeries

    !Destructor
    public  :: KillVegetation                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjVegetation 
    
    !Interfaces----------------------------------------------------------------
    private :: UngetVegetation2D
    private :: UngetVegetation3D
    interface  UngetVegetation
        module procedure UngetVegetation2D
        module procedure UngetVegetation2D_I
        module procedure UngetVegetation2D_L
        module procedure UngetVegetation3D
    end interface UngetVegetation

    !Parameters----------------------------------------------------------------
    !Vegetation Points
    integer, parameter                          :: NoVegetationPoint             = 0
    integer, parameter                          :: VegetationPoint               = 1
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
    !Nutrient uptake Method               
    integer, parameter                          :: NutrientUptake_TranspConc     = 1   !nutrient uptake is soil conc * water uptake
    integer, parameter                          :: NutrientUptakeSWAT            = 2
    integer, parameter                          :: NoNutrientUptake              = 3
    !Nutrient stress Method               
    integer, parameter                          :: NutrientStressUptake          = 1   !stress is actual / potential uptake
    integer, parameter                          :: NutrientStressSWAT            = 2   ! nutrient stress is equations from swat
    !Autofertilization N demand method
    integer, parameter                          :: NTarget_                      = 1
    integer, parameter                          :: AnnualMaxApproach_            = 2

    character(StringLength), parameter          :: beginproperty                 = '<beginproperty>'
    character(StringLength), parameter          :: endproperty                   = '<endproperty>'

    character(StringLength), parameter          :: beginsoilconc                 = '<beginsoilconc>'
    character(StringLength), parameter          :: endsoilconc                   = '<endsoilconc>'
    
    !Plant Type
    integer, parameter                          :: NotAPlant                     = 0
    integer, parameter                          :: WarmSeasonAnnualLegume        = 1
    integer, parameter                          :: ColdSeasonAnnualLegume        = 2
    integer, parameter                          :: PerennialLegume               = 3
    integer, parameter                          :: WarmSeasonAnnual              = 4
    integer, parameter                          :: ColdSeasonAnnual              = 5
    integer, parameter                          :: Perennial                     = 6
    integer, parameter                          :: Trees                         = 7
    

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
        integer, dimension(:,:), pointer                :: KFloor
        integer, dimension(:,:), pointer                :: MappingPoints
        integer, dimension(:,:), pointer                :: BasinPoints       
        real, dimension(:,:), pointer                   :: Topography           
        real, dimension(:,:), pointer                   :: AirTemperature                 !ºC
        real, dimension(:,:), pointer                   :: SolarRadiation                 !W/m2
        real, dimension(:,:), pointer                   :: RelativeHumidity               !0-1
        real, dimension(:,:), pointer                   :: PotentialTranspiration         !m/s     
        real, dimension(:,:), pointer                   :: Latitude                       !deg
        real, dimension(:,:,:), pointer                 :: SoilWaterContent               !m3H2O/m3soil
        real, dimension(:,:,:), pointer                 :: Head                           !m
        real, dimension(:,:,:), pointer                 :: ResidualWaterContent           !m3H2O/m3soil                    
        real, dimension(:,:,:), pointer                 :: SoilNitrate                    !ug/m3H2O    
        real, dimension(:,:,:), pointer                 :: SoilPhosphorus                 !ug/m3H2O   
        real, dimension(:,:,:), pointer                 :: FieldCapacity                  !m3H2O/m3soil  
        real, dimension(:,:), pointer                   :: GridCellArea
        real(8), dimension(:,:,:), pointer              :: CellVolume
        type(T_Integration)                             :: Integration
        real, dimension(:,:,:), pointer                 :: ECw        
    end type   T_External


    type       T_Files
        character(len=Pathlength)                       :: ConstructData
        character(len=Pathlength)                       :: Results
        character(len=Pathlength)                       :: InitialFile
        character(len=Pathlength)                       :: FinalFile
!        character(PathLength)                           :: VegetationIDFile

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
        real, dimension(:,:), pointer                   :: Field                => null()
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

!    type T_SoilConc
!        type(T_PropertyID)                              :: ID
!        real, dimension(:,:,:), pointer                 :: Conc                 => null()
!        real                                            :: ConversionFactor     =  1.0
!        logical                                         :: UseForSalinityStress =  .true.    
!        type(T_SoilConc), pointer                       :: Next                 => null()
!        type(T_SoilConc), pointer                       :: Prev                 => null()
!    end type T_SoilConc

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
        real                                            :: LAIShape1
        real                                            :: LAIShape2
        logical                                         :: Evergreen
        real                                            :: FrGrowLAIDecline
        real                                            :: LAIMax
        real                                            :: LAIDecRate
        real                                            :: MaxCanopyHeight
        real                                            :: OptimalHarvestIndex
        real                                            :: MinimumHarvestIndex
        real                                            :: NitrogenFractionInYeld
        real                                            :: PhosphorusFractionInYeld
        real                                            :: TreeYearsToMaturity
        real                                            :: TreeMaximumBiomass
        real                                            :: BiomassFracRemovedInDormancy
        real                                            :: LAIMinDormant
    end type T_GrowthDatabase

    type T_GrazingDatabase
        real                                            :: GrazingStartJulianDay
        real                                            :: GrazingStartPlantHU
        integer                                         :: GrazingDays
        real                                            :: GrazingMinimumBiomass
        real                                            :: GrazingBiomass
        real                                            :: TramplingBiomass
    end type T_GrazingDatabase

    type T_HarvestKillDatabase
        real                                            :: TramplingBiomass
        real                                            :: HarvestKillJulianDay
        real                                            :: HarvestKillPlantHU
        real                                            :: HarvestJulianDay
        real                                            :: HarvestPlantHU
        real                                            :: KillJulianDay
        real                                            :: KillPlantHU
        real                                            :: HarvestEfficiency
    end type T_HarvestKillDatabase
    
    type T_Fertilizer
        integer                                         :: FertilizerID
        real                                            :: FertilizerFracApplyedInSurface
        real                                            :: OrganicFracParticulate
        real                                            :: AmmoniaFracInMineralN
        real                                            :: MineralNFracInFertilizer
        real                                            :: OrganicNFracInFertilizer
        real                                            :: MineralPFracInFertilizer
        real                                            :: OrganicPFracInFertilizer    
    end type T_Fertilizer

    type T_AutoFertilization
        integer                                         :: FertilizerID
        real                                            :: NTreshold
        real                                            :: PTreshold
        logical                                         :: ExplicitPhosphorus
        real                                            :: NitrogenApplicationMax
        real                                            :: NitrogenAnnualMax
        real                                            :: PhosphorusApplicationMax
        real                                            :: PhosphorusAnnualMax
        integer                                         :: NStress
    end type T_AutoFertilization
    
    type T_FertilizerApps
        type(T_Fertilizer)                              :: Fertilizer       
        integer                                         :: FertilizerID
        real                                            :: FertilizerAppJDay
        real                                            :: FertilizerAppHU
        real                                            :: FertilizerAppAmount !kg/ha.day
        character(StringLength)                         :: FertilizerName
        !Continuous application
        logical                                         :: ContFertilizationON
        integer                                         :: ContFertilizationDays !number of days that amount will be repeated
        logical                                         :: ContFertilizationActive  = .false.
        integer                                         :: ContFertilizationAccDays = 0
        !auto fertilization - as one application method
        type(T_AutoFertilization)                       :: Auto       
        
        logical                                         :: FertilizerAppOccurred
        
    end type T_FertilizerApps

    type T_FertilizerDatabase
        integer                                         :: NumberFertilizerApps = 0
        type(T_FertilizerApps), dimension(:), pointer   :: FertilizerApps
    end type T_FertilizerDatabase

    type T_PesticideApps
        integer                                         :: PesticideID
        real                                            :: PesticideAppJDay
        real                                            :: PesticideAppHU
        real                                            :: PesticideAppAmount
        character(StringLength)                         :: PesticideName

        !Continuous application
        logical                                         :: ContPesticideON
        integer                                         :: ContPesticideDays !number of days that amount will be repeated
        logical                                         :: ContPesticideActive  = .false.
        integer                                         :: ContPesticideAccDays = 0
        
        logical                                         :: PesticideAppOccurred
        
    end type T_PesticideApps

    type T_PesticideDatabase
        integer                                         :: NumberPesticideApps = 0
        type(T_PesticideApps), dimension(:), pointer    :: PesticideApps
    end type T_PesticideDatabase

    type       T_TranspirationMOHID
        integer, dimension(:,:), pointer                :: FeddesType
        real, dimension(:,:), pointer                   :: RootFeddesH1
        real, dimension(:,:), pointer                   :: RootFeddesH2
        real, dimension(:,:), pointer                   :: RootFeddesH3
        real, dimension(:,:), pointer                   :: RootFeddesH4
        real, dimension(:,:), pointer                   :: RootFeddesH3L
        real, dimension(:,:), pointer                   :: RootFeddesH3H
        real, dimension(:,:), pointer                   :: RootFeddesR2L
        real, dimension(:,:), pointer                   :: RootFeddesR2H
    end type   T_TranspirationMOHID

    type T_Evolution
        logical                                         :: ReadNeeded        = .false.
        logical                                         :: GrowthModelNeeded = .false.
        logical                                         :: ModelSWAT         = .false.
        logical                                         :: ModelDEB          = .false.
    end type T_Evolution 

    type       T_ComputeOptions
        logical                                         :: Grazing              = .false.
        logical                                         :: HarvestKill          = .false.
        logical                                         :: Dormancy             = .false.
        logical                                         :: Fertilization        = .false.
        logical                                         :: Pesticide            = .false.
        logical                                         :: AutoFertilization    = .false. !fertilization when plant needs
        logical                                         :: FertilizationApp     = .false. !scheduled fertilization  
        logical                                         :: CheckSpinUp    
        integer                                         :: RootProfile          
        integer                                         :: TranspirationMethod
        integer                                         :: NutrientUptakeMethod
        integer                                         :: NutrientStressMethod
        logical                                         :: LimitTPVel
        logical                                         :: WaterUptakeOld  
        real                                            :: VegetationDT
        real                                            :: IntegrationDT 
 !       logical                                         :: IsPlantGrowing
        logical                                         :: ChangeLAISenescence
        logical                                         :: ChangeCanopyHeight
        real                                            :: PlantingHUBase, PlantHUatMaturity
 !       real                                            :: PotentialHUBase
        real                                            :: WaterUptakeCompensationFactor
        logical                                         :: NutrientReduceUptake
        real                                            :: NutrientReduceUptakeFactor
        logical                                         :: NutrientReduceDemand
        real                                            :: NitrogenDistributionParameter
        real                                            :: PhosphorusDistributionParameter
        real                                            :: AtmosphereCO2
        integer                                         :: WaterUptakeStressMethod
        integer                                         :: SalinityStressMethod
        logical                                         :: ModelNitrogen, ModelPhosphorus
        logical                                         :: ModelRootBiomass, ModelPlantBiomass
        logical                                         :: ModelCanopyHeight, ModelTemperatureStress
        logical                                         :: ModelWater
        logical                                         :: WaterStress
        logical                                         :: AtmospherePropertiesOutput
        logical                                         :: FluxesToSoilOutput
        logical                                         :: NutrientFluxesWithSoil

        logical                                         :: AdjustRUEForCO2
        logical                                         :: AdjustRUEForVPD
        
        logical                                         :: UseLAIDecRate

        logical                                         :: Continuous
        logical                                         :: StopOnWrongDate
        
        logical                                         :: VegParticFertilization = .false.

        type(T_TranspirationMOHID)                      :: TranspirationMOHID
        type(T_Evolution)                               :: Evolution
        integer                                         :: GlobalEvolution
    end type   T_ComputeOptions

    type       T_HeatUnits
        real, dimension(:,:), pointer                   :: PotentialHUTotal
        real, dimension(:,:), pointer                   :: PotentialHUBase
        real, dimension(:,:), pointer                   :: PotentialHUBase_Old
        real, dimension(:,:), pointer                   :: PlantHUAccumulated
        real, dimension(:,:), pointer                   :: PlantHUAccumulated_Old
    end type   T_HeatUnits

    type       T_Growth
        real, dimension(:,:), pointer                   :: WaterStress
        real, dimension(:,:), pointer                   :: NitrogenStress
        real, dimension(:,:), pointer                   :: PhosphorusStress
        real, dimension(:,:), pointer                   :: TemperatureStress
        real, dimension(:,:), pointer                   :: GlobalStress
        real, dimension(:,:), pointer                   :: PAR
        real, dimension(:,:), pointer                   :: RUE
        real, dimension(:,:), pointer                   :: PotentialGrowth
        real, dimension(:,:), pointer                   :: PotentialBiomass
        real, dimension(:,:), pointer                   :: BiomassGrowthOld
        integer, dimension(:,:), pointer                :: TreeCurrentYear
        real, dimension(:,:), pointer                   :: TreeFractionToMaturity
        real, dimension(:,:), pointer                   :: TreeMaximumAnnualBiomass
        logical, dimension(:,:), pointer                :: TreeComingFromContinuous
    end type   T_Growth

    type T_SalinityStressParams
        real :: ECt
        real :: EC0
        real :: Slope        
    end type T_SalinityStressParams

    type T_VegetationType
        integer                                         :: ID
        character(StringLength)                         :: Name
        type(T_VegetationType), pointer                 :: Next, Prev
!        type(T_VegetationType), pointer                 :: FirstVegetation
!        type(T_VegetationType), pointer                 :: LastVegetation
        integer                                         :: VegetationsNumber     = 0
        integer                                         :: FeddesType
        real                                            :: RootFeddesH1
        real                                            :: RootFeddesH2
        real                                            :: RootFeddesH3
        real                                            :: RootFeddesH4
        real                                            :: RootFeddesH3L
        real                                            :: RootFeddesH3H
        real                                            :: RootFeddesR2L
        real                                            :: RootFeddesR2H
        integer                                         :: VegetationID
        type (T_TimingDatabase)                         :: TimingDatabase
        type (T_GrowthDatabase)                         :: GrowthDatabase
        type (T_HarvestKillDatabase)                    :: HarvestKillDatabase
        type (T_GrazingDatabase)                        :: GrazingDatabase
        type (T_FertilizerDatabase)                     :: FertilizerDatabase
        type (T_SalinityStressParams)                   :: SalinityStressParams
        type (T_PesticideDatabase)                      :: PesticideDatabase
        logical                                         :: ComputeRoot
        logical                                         :: ComputeStem
        logical                                         :: HasLeaves
        logical                                         :: UseSalinityStress
        integer                                         :: StressInteraction        
    end type T_VegetationType

    type T_FluxesToSoil
        real, dimension(:,:), pointer                   :: GrazingBiomassToSoil
        real, dimension(:,:), pointer                   :: GrazingNitrogenToSoil
        real, dimension(:,:), pointer                   :: GrazingPhosphorusToSoil
        real, dimension(:,:), pointer                   :: HarvestKillBiomassToSoil
        real, dimension(:,:), pointer                   :: HarvestKillNitrogenToSoil
        real, dimension(:,:), pointer                   :: HarvestKillPhosphorusToSoil
        real, dimension(:,:), pointer                   :: KillRootBiomassLeftInSoil  !In Kill operation.
        real, dimension(:,:), pointer                   :: DormancyBiomassToSoil
        real, dimension(:,:), pointer                   :: DormancyNitrogenToSoil
        real, dimension(:,:), pointer                   :: DormancyPhosphorusToSoil
        real, dimension(:,:), pointer                   :: FertilNitrateToSoilSurface
        real, dimension(:,:), pointer                   :: FertilNitrateToSoilSubSurface
        real, dimension(:,:), pointer                   :: FertilAmmoniaToSoilSurface
        real, dimension(:,:), pointer                   :: FertilAmmoniaToSoilSubSurface
        real, dimension(:,:), pointer                   :: FertilOrganicNToSoilSurface
        real, dimension(:,:), pointer                   :: FertilOrganicNParticToFluff
        real, dimension(:,:), pointer                   :: FertilOrganicNToSoilSubSurface
        real, dimension(:,:), pointer                   :: FertilOrganicPToSoilSurface
        real, dimension(:,:), pointer                   :: FertilOrganicPParticToFluff
        real, dimension(:,:), pointer                   :: FertilOrganicPToSoilSubSurface
        real, dimension(:,:), pointer                   :: FertilMineralPToSoilSurface
        real, dimension(:,:), pointer                   :: FertilMineralPToSoilSubSurface
    end type T_FluxesToSoil
    
    type T_FluxesFromSoil
        real, dimension(:,:,:), pointer                 :: WaterUptakeFromSoil
        real, dimension(:,:,:), pointer                 :: NitrogenUptakeFromSoil
        real, dimension(:,:,:), pointer                 :: PhosphorusUptakeFromSoil
    end type T_FluxesFromSoil
    
    !Block needed to indidualize fluxes by pesticide (e.g. all fertilizer put in N and P but each pesticide is a diff molecule).
    type T_PesticideApp
        integer                                        :: PesticideID
        type (T_ID)                                    :: ID
        real, dimension(:,:), pointer                  :: Soil                   !pesticide flux to soil
        real, dimension(:,:), pointer                  :: Vegetation             !pesticide flux to vegetation
!        logical, dimension(:,:), pointer               :: PesticideAppOccurred 
    end type T_PesticideApp    
    
    type T_PesticideFlux
        integer                                        :: UniquePesticides
        type (T_PesticideApp), dimension(:), pointer   :: Application  !for every unique pesticide
    end type T_PesticideFlux

    type T_Fluxes
        real, dimension(:,:  ), pointer                 :: BiomassGrazed
        real, dimension(:,:  ), pointer                 :: NitrogenGrazed
        real, dimension(:,:  ), pointer                 :: PhosphorusGrazed
        real, dimension(:,:  ), pointer                 :: BiomassGrazedFraction
        real, dimension(:,:,:), pointer                 :: WaterUptakeLayer     !m3/s
        real, dimension(:,:  ), pointer                 :: WaterUptake          !m3/s
        real, dimension(:,:,:), pointer                 :: NitrogenUptakeLayer
        real, dimension(:,:  ), pointer                 :: NitrogenUptake
        real, dimension(:,:,:), pointer                 :: PhosphorusUptakeLayer
        real, dimension(:,:  ), pointer                 :: PhosphorusUptake
        real, dimension(:,:  ), pointer                 :: BiomassGrowth
        real, dimension(:,:  ), pointer                 :: LAIChange
        real, dimension(:,:  ), pointer                 :: BiomassRemovedInHarvest
        real, dimension(:,:  ), pointer                 :: NitrogenRemovedInHarvest
        real, dimension(:,:  ), pointer                 :: PhosphorusRemovedInHarvest
        real, dimension(:,:  ), pointer                 :: BiomassHarvestedFraction
        real, dimension(:,:  ), pointer                 :: BiomassRemovedInDormancy
        real, dimension(:,:  ), pointer                 :: NitrogenRemovedInDormancy
        real, dimension(:,:  ), pointer                 :: PhosphorusRemovedInDormancy
        real, dimension(:,:  ), pointer                 :: FertilNitrateInSurface
        real, dimension(:,:  ), pointer                 :: FertilNitrateInSubSurface
        real, dimension(:,:  ), pointer                 :: FertilAmmoniaInSurface
        real, dimension(:,:  ), pointer                 :: FertilAmmoniaInSubSurface
        real, dimension(:,:  ), pointer                 :: FertilOrganicNInSurface
        real, dimension(:,:  ), pointer                 :: FertilOrganicNInSubSurface
        real, dimension(:,:  ), pointer                 :: FertilOrganicNParticInFluff
        real, dimension(:,:  ), pointer                 :: FertilOrganicPInSurface
        real, dimension(:,:  ), pointer                 :: FertilOrganicPInSubSurface
        real, dimension(:,:  ), pointer                 :: FertilOrganicPParticInFluff
        real, dimension(:,:  ), pointer                 :: FertilMineralPInSurface
        real, dimension(:,:  ), pointer                 :: FertilMineralPInSubSurface
        type (T_PesticideFlux)                          :: Pesticides                    

!        real, dimension(:,:  ), pointer                 :: BiomassDormancyFraction
        type(T_FluxesToSoil)                            :: ToSoil
        type(T_FluxesFromSoil)                          :: FromSoil
    end type T_Fluxes

    type T_StateVariables
        real, dimension(:,:  ), pointer                       :: TotalPlantNitrogen
        real, dimension(:,:  ), pointer                       :: TotalPlantPhosphorus
        real, dimension(:,:  ), pointer                       :: TotalPlantBiomass
        real, dimension(:,:  ), pointer                       :: RootBiomass
        real, dimension(:,:  ), pointer                       :: RootDepth
        real, dimension(:,:  ), pointer                       :: LeafAreaIndex
        real, dimension(:,:  ), pointer                       :: PotLeafAreaIndex
        real, dimension(:,:  ), pointer                       :: SpecificLeafStorage
        real, dimension(:,:  ), pointer                       :: EVTPCropCoefficient
        
        real, dimension(:,:  ), pointer                       :: CanopyHeight
    end type T_StateVariables

    private :: T_Vegetation
    type       T_Vegetation
        integer                                                 :: InstanceID
        type(T_Size3D)                                          :: Size, WorkSize
        type(T_Size2D)                                          :: Size2D, WorkSize2D
        type(T_External)                                        :: ExternalVar
        type(T_Files)                                           :: Files
        type(T_OutPut)                                          :: OutPut
        type(T_Time     )                                       :: BeginTime
        type(T_Time     )                                       :: EndTime
        type(T_Time     )                                       :: ActualTime
        type(T_Time     )                                       :: NextCompute
        type(T_Time     )                                       :: NextIntegration
        type(T_Time     )                                       :: LastOutPutHDF5
        type(T_Property    ), pointer                           :: FirstProperty
        type(T_Property    ), pointer                           :: LastProperty
        integer                                                 :: PropertiesNumber     = 0

        type(T_VegetationType    ), pointer                     :: FirstVegetation
        type(T_VegetationType    ), pointer                     :: LastVegetation
        integer                                                 :: VegetationsNumber     = 0

        !DataMatrixes
        integer, dimension(:,:), pointer                        :: AgricPractID !ID's from user - SWAT ID'S
        type (T_PropertyID)                                     :: AgricPract
        
        !ordered ID's to save space (max value = number of different ID's)
        integer, dimension(:,:), pointer                        :: VegetationID    
                                                                                   
        integer                                                 :: ObjEnterData         = 0
        integer                                                 :: ObjTime              = 0
        integer                                                 :: ObjGridData          = 0
        integer                                                 :: ObjHorizontalGrid    = 0
        integer                                                 :: ObjHorizontalMap     = 0
        integer                                                 :: ObjGeometry          = 0
        integer                                                 :: ObjBasinGeometry     = 0
        integer                                                 :: ObjTimeSerie         = 0
        integer                                                 :: ObjTimeSerieToSoil   = 0
        integer                                                 :: ObjTimeSerieAtm      = 0
        integer                                                 :: ObjFillMatrix        = 0
        integer                                                 :: ObjAtmosphere        = 0
        integer                                                 :: ObjPorousMedia       = 0
        integer                                                 :: ObjHDF5              = 0
        type(T_Vegetation), pointer                             :: Next
        
        logical                                                 :: UsePotLAI      = .false.
        logical                                                 :: UseBoundaryLAI = .false.
        
        logical, dimension(:,:), pointer                        :: IsPlantGrowing       
        logical, dimension(:,:), pointer                        :: PlantingOccurred      
        logical, dimension(:,:), pointer                        :: KillOccurred        
        logical, dimension(:,:), pointer                        :: HarvestOnlyOccurred   
        logical, dimension(:,:), pointer                        :: HarvestKillOccurred
        logical, dimension(:,:), pointer                        :: IsPlantDormant       
        logical, dimension(:,:), pointer                        :: PlantGoingDormant    
        logical, dimension(:,:), pointer                        :: IsPlantBeingGrazed
        logical, dimension(:,:), pointer                        :: LAISenescence        
        logical, dimension(:,:), pointer                        :: ChangeCanopyEnabled      
        logical, dimension(:,:), pointer                        :: SoilFluxesActive

        !global variables to be used in several routines
        real, dimension(:,:), pointer                           :: PlantNitrogenFraction
        real, dimension(:,:), pointer                           :: PlantPhosphorusFraction
        real, dimension(:,:), pointer                           :: PlantLAIMaxFraction
        real, dimension(:,:), pointer                           :: LAIBeforeSenescence
        real, dimension(:,:  ), pointer                         :: LAIDeclineFraction
        real, dimension(:,:), pointer                           :: DayLength
        real, dimension(:,:), pointer                           :: MinimumDayLength
        integer, dimension(:,:), pointer                        :: TranspirationBottomLayer
        
        character(len = StringLength)                           :: GrowthDatabase
        character(len = StringLength)                           :: FertilizerDatabase
        character(len = StringLength)                           :: PesticideDatabase
        character(len = StringLength)                           :: FeddesDatabase

        !global counter
        real, dimension(:,:), pointer                           :: DaysOfGrazing                    !counter to days of grazing
        real, dimension(:,:), pointer                           :: AnnualNitrogenFertilized
        real, dimension(:,:), pointer                           :: AnnualPhosphorusFertilized
        real, dimension(:,:), pointer                           :: OptimalTotalPlantNitrogen
        real, dimension(:,:), pointer                           :: OptimalTotalPlantPhosphorus
        integer                                                 :: nIterations                      !counter to atmosphere integration
        
        integer, dimension(:), pointer                          :: PesticideListID
        character(len = StringLength), dimension(:), pointer    :: PesticideListName
       
        type(T_VegetationType), dimension(:), pointer           :: VegetationTypes  => null()
        type(T_ComputeOptions)                                  :: ComputeOptions
        type(T_HeatUnits)                                       :: HeatUnits
        type(T_Growth)                                          :: Growth
        type(T_Fluxes)                                          :: Fluxes
        type(T_StateVariables)                                  :: StateVariables
        
        real, dimension(:,:), pointer                           :: RootDepthOld
        
        !Auto fertilization target
        integer                                                 :: SimulationYear = 1
        real, dimension(:,:), pointer                           :: NitrogenYeldEstimate
        real, dimension(:,:), pointer                           :: NitrogenYeldTarget
        
        real, dimension(:,:,:), pointer                         :: ECw => null()      
        
        logical                                                 :: AllowNegativeLAI   = .false.
        logical                                                 :: CorrectNegativeLAI = .false.
        real                                                    :: ValueInsteadNegativeLAI

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
                                   BasinGeometryID,                 &
                                   CoupledAtmosphere,               &
                                   UsePotLAI,                       &
                                   STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjVegetationID 
        integer                                         :: TimeID         
        integer                                         :: GridDataID     
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: AtmosphereID
        integer                                         :: PorousMediaID
        integer, dimension(:,:), pointer                :: MappingPoints
        integer                                         :: GeometryID
        integer                                         :: BasinGeometryID
        logical                                         :: CoupledAtmosphere
        logical,           intent(OUT)                  :: UsePotLAI
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

            Me%UsePotLAI          = .false.
            Me%UseBoundaryLAI     = .false.
            Me%AllowNegativeLAI   = .false.
            Me%CorrectNegativeLAI = .false.
            
            Me%ValueInsteadNegativeLAI = 0.0
            
            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjAtmosphere     = AssociateInstance (mATMOSPHERE_,     AtmosphereID    )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            
            Me%ExternalVar%MappingPoints => MappingPoints
            Me%ExternalVar%CoupledAtmosphere = CoupledAtmosphere                       

            call ReadVegetationFileNames

            
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetation - ModuleVegetation - ERR01' 
                                                
            call ConstructGlobalVariables   !Read Data file Options
            call AllocateVariables          !Allocate global variables
            call AllocateEvolutionVariables !Allocate evolution method dependent variables
            
            call ConstructPropertyList      !Constructs the list of vegetation properties            
            call CreatePropertiesAlias      !Creates the aliases for the properties

            call ConstructVegetationList
            call ConstructVegetationParameters             
            call CheckOptionsConsistence
            call AllocatePropertyVariables  !Allocate property dependent variables
            
            !Grid operations            
            
            call ConstructVegetationGrids

            !call CheckRootDepth
            
            if (Me%ComputeOptions%Dormancy) then
                call CheckLatitude
            endif
            
            if ((Me%ComputeOptions%ModelWater) .and.                                &
                (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) .and. &
                (Me%ComputeOptions%WaterStress)) then
                call ConstructFeddes
            endif
            
            !Initial conditions read if "real" continuous simulation
            if (Me%ComputeOptions%CheckSpinUp) then
                if (Me%ComputeOptions%Continuous .and. Me%ComputeOptions%StopOnWrongDate    &
                    .and. Me%ComputeOptions%Evolution%ModelSWAT) then
                    call ReadInitialHDF
                else
                    Me%Growth%TreeCurrentYear (:,:) = 0
                    Me%Growth%TreeComingFromContinuous (:,:) = .false.                
                endif
            else if (Me%ComputeOptions%Continuous .and. Me%ComputeOptions%Evolution%ModelSWAT) then
                call ReadInitialHDF
            else
                Me%Growth%TreeCurrentYear (:,:) = 0
                Me%Growth%TreeComingFromContinuous (:,:) = .false.
            endif

            !Output
            call ConstructTimeSerie

            call ConstructHDF


            call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetation - ModuleVegetation - ERR02'


            call ConstructLog

            !First Output
            if (Me%OutPut%HDF_ON)           call Modify_OutputHDF                        

            !Returns ID
            ObjVegetationID  = Me%InstanceID
            UsePotLAI        = Me%UsePotLAI

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

        !Reads the name of the file where to store final data
        call ReadFileName ('VEGETATION_FIN', Me%Files%FinalFile, "Vegetation Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadVegetationFileNames - ModuleVegetation - ERR03'


    end subroutine ReadVegetationFileNames

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag
        real                                    :: auxFactor, Erroraux, DTaux, ModelDT
        real                                    :: PotentialHUTotal
        logical                                 :: Aux
        !Begin------------------------------------------------------------------

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = Me%BeginTime,                 &
                                  EndTime = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR010'
        
        !Actualize the time
        Me%ActualTime  = Me%BeginTime
        Me%NextCompute = Me%ActualTime

        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                                &
                                   Size     = Me%Size2D,                                &
                                   WorkSize = Me%WorkSize2D,                            &
                                   STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalVariables - ModuleVegetation  - ERR020'

        call GetGeometrySize    (Me%ObjGeometry,                                        &    
                                 Size     = Me%Size,                                    &
                                 WorkSize = Me%WorkSize,                                &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation  - ERR030'

        call GetData(Me%ComputeOptions%GlobalEvolution,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'MODEL_EVOLUTION',                                      &                     
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR040'
        if (iflag /= 1) then
            write (*,*) ''
            write (*,*) 'ATTENTION'
            write (*,*) 'The keyword MODEL_EVOLUTION is missing in the vegetation input file.'
            write (*,*) 'The possible values are:'
            write (*,*) '  1 - Vegetation growth will NOT be modeled.'
            write (*,*) '  2 - SWAT method will be used to model vegetation growth.'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR041'
        endif
        select case (Me%ComputeOptions%GlobalEvolution)
        case (1)
            Me%ComputeOptions%Evolution%ReadNeeded        = .true.
            Me%ComputeOptions%Evolution%GrowthModelNeeded = .false.
            Me%ComputeOptions%Evolution%ModelSWAT         = .false.
            Me%ComputeOptions%Evolution%ModelDEB          = .false.        
        case (2)
            Me%ComputeOptions%Evolution%ReadNeeded        = .false.
            Me%ComputeOptions%Evolution%GrowthModelNeeded = .true.
            Me%ComputeOptions%Evolution%ModelSWAT         = .true.
            Me%ComputeOptions%Evolution%ModelDEB          = .false.                
        case default
            write (*,*) ''
            write (*,*) 'ATTENTION'
            write (*,*) 'The keyword GROWTH_METHOD value in the vegetation input file is INVALID.'
            write (*,*) 'The possible values are:'
            write (*,*) '  1 - Vegetation growth will NOT be modeled.'
            write (*,*) '  2 - SWAT method will be used to model vegetation growth.'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR042'        
        end select

        call GetData(Me%ComputeOptions%ModelTemperatureStress,                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TEMPERATURE_STRESS',                             &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR050'

        call GetData(Me%ComputeOptions%ModelNitrogen,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'NITROGEN_STRESS',                                &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR060'

        call GetData(Me%ComputeOptions%ModelPhosphorus,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'PHOSPHORUS_STRESS',                              &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR070'
        
        call GetData(Me%ComputeOptions%ModelWater,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'WATER_UPTAKE',                                   &
                     Default        = .true.,                                           &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR080'        
        
        if (Me%ComputeOptions%ModelWater) then
            call GetData(Me%ComputeOptions%TranspirationMethod,                             &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'WATER_UPTAKE_METHOD',                            &
                         Default        = TranspirationMOHID,                               &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR081'
            
            call GetData(Me%ComputeOptions%WaterStress,                                     &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'WATER_STRESS',                                   &
                         Default        = .true.,                                           &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR082'            
        
            if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then        

                call GetData(Me%ComputeOptions%LimitTPVel,                                  &
                             Me%ObjEnterData, iflag,                                        &
                             Keyword        = 'LIMIT_TRANSP_WATER_VEL',                     &
                             Default        = .false.,                                      &
                             SearchType     = FromFile,                                     &
                             ClientModule   = 'ModuleVegetation',                           &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR100'
            
            
                call GetData(Me%ComputeOptions%RootProfile,                                 &
                             Me%ObjEnterData, iflag,                                        &
                             Keyword        = 'ROOT_PROFILE',                               &
                             Default        = RootTriangular,                               &
                             SearchType     = FromFile,                                     &
                             ClientModule   = 'ModuleVegetation',                           &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR110'

                call GetData(Me%ComputeOptions%WaterUptakeOld,                              &
                             Me%ObjEnterData, iflag,                                        &
                             Keyword        = 'WATER_UPTAKE_OLD',                           &
                             Default        = .false.,                                      &
                             SearchType     = FromFile,                                     &
                             ClientModule   = 'ModuleVegetation',                           &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR120'        

                call GetData(Me%ComputeOptions%WaterUptakeStressMethod,                     &
                             Me%ObjEnterData, iflag,                                        &
                             Keyword        = 'WATER_UPTAKE_STRESS_METHOD',                 &
                             Default        = Feddes,                                       &
                             SearchType     = FromFile,                                     &
                             ClientModule   = 'ModuleVegetation',                           &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR130'
                
                call GetData(Me%ComputeOptions%SalinityStressMethod,                        &
                             Me%ObjEnterData, iflag,                                        &
                             Keyword        = 'SALINITY_STRESS_METHOD',                     &
                             Default        = 1,                                            & !Threshold/Slope
                             SearchType     = FromFile,                                     &
                             ClientModule   = 'ModuleVegetation',                           &
                             STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR140'    
                                
            endif
            
            call GetData(Me%ComputeOptions%WaterUptakeCompensationFactor,                   &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'WATER_UPTAKE_COMPENSATION_FACTOR',               &
                         Default        = 0.,                                               &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR150'
       
        endif
        
        call GetData(Me%ComputeOptions%NutrientUptakeMethod,                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'NUTRIENT_UPTAKE_METHOD',                         &
                     Default        = NutrientUptake_TranspConc,                        &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR160'  
        
        if (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
            call GetData(Me%ComputeOptions%NutrientStressMethod,                            &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'NUTRIENT_STRESS_METHOD',                         &
                         Default        = NutrientStressUptake,                             &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR170'

            call GetData(Me%ComputeOptions%NitrogenDistributionParameter,                   &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'NITROGEN_DISTRIBUTION_PARAMETER',                &
                         Default        = 20.,                                              &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR180'

            call GetData(Me%ComputeOptions%PhosphorusDistributionParameter,                 &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'PHOSPHORUS_DISTRIBUTION_PARAMETER',              &
                         Default        = 20.,                                              &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR190'
            
            call GetData(Me%ComputeOptions%NutrientReduceUptake,                 &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'NUTRIENT_REDUCE_UPTAKE',                         &
                         Default        = .true.,                                           &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR195'

            call GetData(Me%ComputeOptions%NutrientReduceUptakeFactor,                      &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'NUTRIENT_REDUCE_UPTAKE_FACTOR',                  &
                         Default        = 0.5,                                              &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR196'

            call GetData(Me%ComputeOptions%NutrientReduceDemand,                            &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'DEMAND_REDUCE_UPTAKE',                           &
                         Default        = .false.,                                          &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR197'
            
        endif

        call GetData(Me%ComputeOptions%AtmosphereCO2,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ATMOSPHERE_CO2',                                 &
                     Default        = 330.,                                             &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR200'

        call GetData(PotentialHUTotal,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'POTENTIALHUTOTAL',                               &
                     Default        = 0.,                                               &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (iflag /= 0) then
            write(*,*) 
            write(*,*) 'Potential Annual HU  for the simulation is now defined in'
            write(*,*) ' <begin_TotalPotentialHU>/<end_TotalPotentialHU> block'
            write(*,*) 'and not with POTENTIALHUTOTAL keyword'    
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR210'
        endif

        call GetData(Me%ComputeOptions%CheckSpinUp,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'CHECK_SPINUP',                                   &
                     Default        = .true.,                                           &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR211'    

        call GetData(Me%ComputeOptions%Grazing,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'GRAZING',                                        &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR220'

        call GetData(Aux,                                                               &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'MANAGEMENT',                                     &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0170'
        if (iflag /= 0) then
            write(*,*)'Found MANAGEMENT keyword in Vegetation data file'
            write(*,*)'This keyword is obsolete and is now HARVEST_KILL'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR0175'
        endif
        
        call GetData(Me%ComputeOptions%HarvestKill,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'HARVEST_KILL',                                   &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR230'
  
        call GetData(Me%ComputeOptions%Dormancy,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'DORMANCY',                                       &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR240'

        call GetData(Me%ComputeOptions%Fertilization,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'FERTILIZATION',                                  &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR250'

        call GetData(Me%ComputeOptions%Pesticide,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'PESTICIDE',                                      &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR260'
        
        if (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
            call GetData(Me%ComputeOptions%NutrientFluxesWithSoil,                          &
                         Me%ObjEnterData, iflag,                                            &
                         Keyword        = 'NUTRIENT_FLUXES_WITH_SOIL',                      &
                         Default        = .true.,                                           &
                         SearchType     = FromFile,                                         &
                         ClientModule   = 'ModuleVegetation',                               &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR270'
        else
            !PorousMedia properties does not need to compute nutrient fluxes from vegetation
            Me%ComputeOptions%NutrientFluxesWithSoil = .false.
        endif
        
        call GetData(Me%ComputeOptions%ChangeLAISenescence,                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'CHANGE_LAI_SENESCENCE',                          &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR280'

        call GetData(Me%ComputeOptions%ChangeCanopyHeight,                              &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'CHANGE_CANOPY_HEIGHT',                           &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR290'

        call GetData(Me%ComputeOptions%AdjustRUEForCO2,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ADJUST_RUE_FOR_CO2',                             &
                     Default        = .true.,                                           &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR300'

        call GetData(Me%ComputeOptions%AdjustRUEForVPD,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ADJUST_RUE_FOR_VPD',                             &
                     Default        = .true.,                                           &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0310'
        
        call GetData(Me%ComputeOptions%UseLAIDecRate,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'USE_LAI_DEC_RATE',                               &
                     Default        = .false.,                                          &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleVegetation',                               &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR0311'
                
        call GetData(Me%ComputeOptions%Continuous,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'CONTINUOUS',                                         &
                     Default    = .false.,                                              &                                           
                     ClientModule ='ModulePorousMedia',                                 &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR320'
        
        if (Me%ComputeOptions%Continuous) then
            !Reads the name of the file where to read initial data
            call ReadFileName ('VEGETATION_INI', Me%Files%InitialFile, "Vegetation Initial File", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR330'

            call GetData(Me%ComputeOptions%StopOnWrongDate,                             &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType = FromFile,                                         &
                         keyword    = 'STOP_ON_WRONG_DATE',                             &
                         Default    = .true.,                                           &                                           
                         ClientModule ='ModulePorousMedia',                             &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMedia - ERR340'

        endif

        call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR350'

        Me%ExternalVar%DT = ModelDT

        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleVegetation - ERR360'

        call GetData(Me%ComputeOptions%VegetationDT,                                     &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'VEGETATION_DT',                                     &
                     Default      = ModelDT,                                             &
                     SearchType   = FromFile,                                            &
                     ClientModule = 'ModuleVegetation',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR370'
                                   
        call GetData(Me%ComputeOptions%IntegrationDT,                                    &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'INTEGRATION_DT',                                    &
                     Default      = ModelDT,                                             &
                     SearchType   = FromFile,                                            &
                     ClientModule = 'ModuleVegetation',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR380'

        if (Me%ComputeOptions%VegetationDT .lt. Me%ComputeOptions%IntegrationDT) then
            write(*,*) 
            write(*,*) 'Vegetation DT time step is smaller then the external variables integration time step'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR390'
        endif
        
        if (Me%ComputeOptions%IntegrationDT .lt. Me%ExternalVar%DT) then
            write(*,*) 
            write(*,*) 'Integration DT time step is smaller then the model time step'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR400'  
        endif      
        
        if (Me%ComputeOptions%VegetationDT .lt. ModelDT) then
            write(*,*) 
            write(*,*) 'Vegetation DT time step is smaller then model time step'
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR410'

        else

            !Vegetation time step must be a multiple of the model time step
            auxFactor = Me%ComputeOptions%VegetationDT  / ModelDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) 'Vegetation time step must be a multiple of model time step.'
                write(*,*) 'Please review your input data.'
                stop 'ConstructGlobalVariables - ModuleVegetation - ERR420'
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
                stop 'ConstructGlobalVariables - ModuleVegetation - ERR430'
            end if
        endif

        Me%NextCompute = Me%ExternalVar%Now !+ Me%ComputeOptions%VegetationDT

        Me%NextIntegration = Me%ExternalVar%Now !+ Me%ComputeOptions%IntegrationDT
        
        call GetData(Me%AllowNegativeLAI,                                                &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'ALLOW_NEGATIVE_LAI',                                &
                     Default      = .false.,                                             &
                     SearchType   = FromFile,                                            &
                     ClientModule = 'ModuleVegetation',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleVegetation - ERR440'        
        if (.not. Me%AllowNegativeLAI) then 
        
            !This is experimental...      
            call GetData(Me%CorrectNegativeLAI,                                          &
                         Me%ObjEnterData, iflag,                                         &
                         keyword      = 'CORRECT_NEGATIVE_LAI',                          &
                         Default      = .false.,                                         &
                         SearchType   = FromFile,                                        &
                         ClientModule = 'ModuleVegetation',                              &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'ConstructGlobalVariables - ModuleVegetation - ERR450'        
!            if (iflag /= 1) then
!                write (*,*) "ATTENTION"
!                write (*,*) "ALLOW_NEGATIVE_LAI is set to FALSE and"
!                write (*,*) "CORRECT_NEGATIVE_LAI is missing. Check your"
!                write (*,*) "vegetation input file"
!                stop 'ConstructGlobalVariables - ModuleVegetation - ERR460'        
!            endif
            if (Me%CorrectNegativeLAI) then
                call GetData(Me%ValueInsteadNegativeLAI,                                     &
                             Me%ObjEnterData, iflag,                                         &
                             keyword      = 'VALUE_INSTEAD_NEGATIVE_LAI',                    &
                             SearchType   = FromFile,                                        &
                             default      = 0.0,                                             &
                             ClientModule = 'ModuleVegetation',                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                 &
                    stop 'ConstructGlobalVariables - ModuleVegetation - ERR470'  
            endif
        endif

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
                    
                    if (NewProperty%ID%IDNumber .EQ. PotLeafAreaIndex_) then
                        Me%UsePotLAI = .true.
                    endif

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

        call ConstructPropertyEvolution (NewProperty)
        
        call ConstructPropertyValues    (NewProperty)

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
        
        !If "real" continuous cycle and using SWAT model, all properties modelled have to have old active
        if(NewProperty%Evolution == SWAT .and. Me%ComputeOptions%Continuous   & 
           .and. Me%ComputeOptions%StopOnWrongDate .and. .not. NewProperty%Old) then
            write(*,*)'ERROR! Property ', trim(adjustl(NewProperty%ID%name))
            write(*,*)'has not OLD keyword active but global keyword CONTINUOUS is.'
            stop 'ConstructPropertyValues - ModuleVegetation - ERR06'
        endif
        
        if (Me%ComputeOptions%CheckSpinUp) then
            !Do not use old property values if on spin-up periods (final HDF file will have different date) 
            if(Me%ComputeOptions%Continuous .and. .not. Me%ComputeOptions%StopOnWrongDate .and. NewProperty%Old) then
                write(*,*)'ERROR! Property ', trim(adjustl(NewProperty%ID%name))
                write(*,*)'has OLD keyword active but user is running spin-up with "wrongdate" on.'
                stop 'ConstructPropertyValues - ModuleVegetation - ERR06.5'
            endif
        endif
        
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
                                       PointsToFill2D       = Me%ExternalVar%MappingPoints,     &
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

!            if (Me%ComputeOptions%StopOnWrongDate) then

                ! If the property is old then the program is going to try to find a property
                ! with the same name in the Vegetation initial file written in HDF format  
                call ReadOldHDF(NewProperty)
                 
!                call ReadInitialFile(NewProperty)

!            endif

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
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                               
                if (NewProperty%Field(i, j) < 0.) then

                    NewProperty%Field(i, j) = 0.
                
                elseif (NewProperty%Evolution == SWAT .and. NewProperty%Field(i, j) > 0.) then

                    write(*,*) 
                    write(*,*) 'Vegetation growth model is active and'
                    write(*,*) 'Property values are greater than 0.'
                    write(*,*) 'By default there is no plant growing at the beggining of the simulation'
                    stop 'CheckFieldConsistence - ModuleVegetation - ERR010'                                       
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


        inquire (FILE=trim(Me%Files%InitialFile)//"5", EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialFile)//"5",&
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
                     Default      = Me%ComputeOptions%GlobalEvolution,                   &
                     STAT         = STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'ConstructPropertyEvolution - ModuleVegetation - ERR010'

        if (NewProperty%Evolution == ReadValue) then
            Me%ComputeOptions%Evolution%ReadNeeded = .true.
        elseif (Me%ComputeOptions%GlobalEvolution == ReadValue) then
            write(*,*) 
            write(*,*) 'ATTENTION. Property: ', NewProperty%ID%Name 
            write(*,*) 'Global EVOLUTION keyword is set to 1 (no vegetation growth will be modeled)'
            write(*,*) 'But the property EVOLUTION keyword is set to a diferent value.' 
            write(*,*) 'If no vegetation growth will be modeled, all properties'
            write(*,*) 'must be provided (property EVOLUTION = 1)'
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR020'              
        elseif (NewProperty%Evolution /= Me%ComputeOptions%GlobalEvolution) then
            write(*,*) 
            write(*,*) 'ATTENTION. Property: ', NewProperty%ID%Name 
            write(*,*) 'If the property EVOLUTION keyword is set to a DIFFERENT value than 1, '
            write(*,*) 'Than it must be the same than the GLOBAL keyword EVOLUTION.' 
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR030'        
        endif

        if (NewProperty%ID%IDNumber .EQ. BoundaryLeafAreaIndex_) then
            if (NewProperty%Evolution /= ReadValue) then
                write(*,*) 
                write(*,*) "'boundary Leaf area index' property can't be modeled."
                write(*,*) 'EVOLUTION keyword must be set to 1 (= READ)'
                stop 'Construct_PropertyEvolution - ModuleVegetation - ERR040'                
            endif
            
            Me%UseBoundaryLAI = .true.        
        endif    

        if (NewProperty%Evolution == DEB) then
            write(*,*) 
            write(*,*) 'DEB vegetation model is not yet implemented.'
            write(*,*) 'Check EVOLUTION keyword options'
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR050'        
        endif
        
        call GetData(NewProperty%IsConstant,                                             &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'REMAIN_CONSTANT',                                   &
                     ClientModule = 'ModuleVegetation',                                  &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleVegetation - ERR060'

        if (.not. NewProperty%Evolution == ReadValue .and. NewProperty%IsConstant) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Property ', trim(NewProperty%ID%Name) 
            write(*,*    ) 'in Vegetation file has model and constant value options active.'
            write(*,*    ) 'Check EVOLUTION and REMAIN_CONSTANT keywords. '
            stop 'ConstructPropertyEvolution - ModuleVegetation - ERR070'  
        endif

        if (.not. NewProperty%Evolution == ReadValue .and. NewProperty%ID%SolutionFromFile) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Property ', trim(NewProperty%ID%Name)
            write(*,*    ) 'in Vegetation file has model and read value options active.'
            write(*,*    ) 'Check EVOLUTION and FILE_IN_TIME keywords. '
            stop 'ConstructPropertyEvolution - ModuleVegetation - ERR080'  
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
                stop 'Construct_PropertyEvolution - ModuleVegetation - ERR090'
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

!    subroutine ConstructSoilConcList
!
!        !Arguments-------------------------------------------------------------
!
!        !External--------------------------------------------------------------                                  
!        integer :: STAT_CALL
!        integer :: ClientNumber
!        logical :: BlockFound
!                                  
!        !Local-----------------------------------------------------------------
!        type (T_SoilConc), pointer :: NewSoilConc
!
!        !----------------------------------------------------------------------
!        do
!            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                &
!                                        beginsoilconc, endsoilconc, BlockFound,       &
!                                        STAT = STAT_CALL)
!
!            if (STAT_CALL .EQ. SUCCESS_) then
!    
!                if (BlockFound) then                                  
!                
!                    ! Construct a new Soil Concentration Property 
!                    call ConstructSoilConc(NewSoilConc)
!                                             
!                    ! Add new Soil Concentration Property to the Soil Concentration Property List 
!                    call AddSoilConc(NewSoilConc)
!                    
!                else
!                
!                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
!                    if (STAT_CALL .NE. SUCCESS_) &
!                        stop 'ConstructSoilConcList - ModuleVegetation. ERR010.'
!                    exit !No more blocks
!                    
!                end if    
!
!            else
!            
!                write(*,*)  
!                write(*,*) 'Error calling ExtractBlockFromBuffer. '
!                stop 'ConstructSoilConcList - ModuleVegetation - ERR020'
!
!            end if     
!            
!        end do    
!        !----------------------------------------------------------------------
!        
!    end subroutine ConstructSoilConcList
!       
!    !--------------------------------------------------------------------------
!               
!    subroutine ConstructSoilConc (NewSoilConc)
!
!        !Arguments-------------------------------------------------------------
!        type(T_SoilConc), pointer :: NewSoilConc
!
!        !----------------------------------------------------------------------             
!        allocate (NewSoilConc)
!
!        nullify (NewSoilConc%Conc)
!        nullify (NewSoilConc%Prev)
!        nullify (NewSoilConc%Next)        
!
!        call ConstructPropertyID (NewSoilConc%ID, Me%ObjEnterData, FromBlock)
!        call ConstructSoilConcOptions (NewSoilConc)
!        !----------------------------------------------------------------------        
!
!    end subroutine ConstructSoilConc
!
!    !--------------------------------------------------------------------------    
!
!    subroutine ConstructSoilConcOptions (NewSoilConc)
!
!        !Arguments-------------------------------------------------------------
!        type(T_SoilConc), pointer :: NewSoilConc
!
!        !External--------------------------------------------------------------
!        integer                   :: STAT_CALL
!
!        !Local-----------------------------------------------------------------
!        integer                   :: iflag
!
!        !----------------------------------------------------------------------       
!        call GetData(NewSoilConc%UseForSalinityStress,                                  &
!                     Me%ObjEnterData, iflag,                                            &
!                     keyword      = 'USE_FOR_SALINITY_STRESS',                          &
!                     Default      = .true.,                                             &                        
!                     SearchType   = FromBlock,                                          &
!                     ClientModule = 'ModuleVegetation',                                 &
!                     STAT         = STAT_CALL)              
!        if (STAT_CALL .NE. SUCCESS_) &
!            stop 'ConstructSoilConcOptions - ModuleVegetation - ERR010'
!        
!        if (NewSoilConc%UseForSalinityStress) then
!        
!            call GetData(NewSoilConc%ConversionFactor,                                      &
!                         Me%ObjEnterData, iflag,                                            &
!                         keyword      = 'CONVERSION_FACTOR',                                &
!                         SearchType   = FromBlock,                                          &
!                         ClientModule = 'ModuleVegetation',                                 &
!                         STAT         = STAT_CALL)              
!            if (STAT_CALL .NE. SUCCESS_) &
!                stop 'ConstructSoilConcOptions - ModuleVegetation - ERR020' 
!            if (iflag .NE. 1) &
!                stop 'ConstructSoilConcOptions - ModuleVegetation - ERR030' 
!                    
!        else
!        
!            NewSoilConc%ConversionFactor = 0.0
!            
!        endif               
!        !----------------------------------------------------------------------
! 
!    end subroutine ConstructSoilConcOptions
!
!    !--------------------------------------------------------------------------
!     
!    subroutine AddSoilConc (NewSoilConc)
!
!        !Arguments-------------------------------------------------------------
!        type(T_SoilConc), pointer :: NewSoilConc
!
!        !----------------------------------------------------------------------
!        if (.not.associated(Me%FirstSoilConc)) then
!            Me%SoilConcNumber     =  1
!            Me%FirstSoilConc      => NewSoilConc
!            Me%LastSoilConc       => NewSoilConc
!        else
!            NewSoilConc%Prev      => Me%LastSoilConc
!            Me%LastSoilConc%Next  => NewSoilConc
!            Me%LastSoilConc       => NewSoilConc
!            Me%SoilConcNumber     =  Me%SoilConcNumber + 1
!        end if 
!        !----------------------------------------------------------------------
!
!    end subroutine AddSoilConc 

   !--------------------------------------------------------------------------

    ! This subroutine adds a new vegetation type to the vegetation List  
    subroutine AddVegetation(NewVegetation)

        !Arguments-------------------------------------------------------------
        type(T_VegetationType),           pointer     :: NewVegetation

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstVegetation)) then
            Me%VegetationsNumber    = 1
            Me%FirstVegetation      => NewVegetation
            Me%LastVegetation       => NewVegetation
        else
            NewVegetation%Prev      => Me%LastVegetation
            Me%LastVegetation%Next  => NewVegetation
            Me%LastVegetation       => NewVegetation
            Me%VegetationsNumber    = Me%VegetationsNumber + 1
        end if 


    end subroutine AddVegetation 

   !--------------------------------------------------------------------------
    
    subroutine ConstructVegetationList

        !Local----------------------------------------------------------------
        integer                                 :: STAT_CALL !, ObjGD
        integer                                 :: ClientNumber !, i, j
        real, dimension(:,:), pointer           :: AgricPractIDMatrix, DataMatrix
        logical                                 :: BlockFound
        integer                                 :: ObjHDF5, NumberOfInstants, iflag, ObjTimeSerie
        integer                                 :: AgricPractIDScalar, HDFInstant, TimeSerieInstant
        integer                                 :: column
        character(len=StringLength)             :: FileInTime, FileName, field_name
        integer                                 :: HDF5_READ, DataValues

        !Begin----------------------------------------------------------------
    
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationList - ModuleVegetation - ERR00'
        
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = ClientNumber,                         &
                                    block_begin     = '<begin_AgriculturalPractices>',      &
                                    block_end       = '<end_AgriculturalPractices>',        &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. BlockFound) then
                write(*,*)'Missing Block <begin_AgriculturalPractices> / <end_AgriculturalPractices>'
                write(*,*)'This block removes old VEGETATION_ID_FILE and now input is done using'
                write(*,*)'FillMatrix standards: constant, timeserie, ASCII grid or HDF input'
                stop 'ConstructVegetationList - ModuleVegetation - ERR10'
            endif        
            
            allocate (AgricPractIDMatrix(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            AgricPractIDMatrix = FillValueReal
                        
            !fill the first grid (in case of solution from file) or the unique grid
            call ConstructFillMatrix  (PropertyID           = Me%AgricPract,                    &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExternalVar%MappingPoints,     &
                                       Matrix2D             = AgricPractIDMatrix,               &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'ConstructVegetationList - ModuleVegetation - ERR20'
            
!            !Unget Basin Points
!            call UnGetBasin (Me%ObjBasinGeometry, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationList - ModuleVegetation - ERR025'

            if(.not. Me%AgricPract%SolutionFromFile)then

                ! only one grid
                call KillFillMatrix(Me%AgricPract%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructVegetationList - ModuleVegetation - ERR30'
                
                !Check in the only grid the unique vegetations
                call CheckVegetationListMatrix(AgricPractIDMatrix)
            
            else
                
                !several grids - check all values for different agricultural practices
                !Check if HDF or time serie
                call GetData(FileInTime,                                            &
                             Me%ObjEnterData, iflag,                                &
                             SearchType   = FromFile,                               &
                             keyword      = 'FILE_IN_TIME',                         &
                             ClientModule = 'ModuleVegetation',                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ConstructVegetationList - ModuleVegetation - ERR40'

                call GetData(FileName,                                              &
                             Me%ObjEnterData,iflag,                                 &
                             SearchType   = FromFile,                               &
                             keyword      = 'FILENAME',                             &
                             ClientModule = 'ModuleVegetation',                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ConstructVegetationList - ModuleVegetation - ERR50'     
                    
                call GetData(column,                                                &
                             Me%ObjEnterData,iflag,                                 &
                             SearchType   = FromFile,                               &
                             keyword      = 'DATA_COLUMN',                          &
                             default      = 2,                                      &
                             ClientModule = 'ModuleVegetation',                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ConstructVegetationList - ModuleVegetation - ERR51'     
                    
                call GetData(field_name,                                            &
                             Me%ObjEnterData,iflag,                                 &
                             SearchType   = FromFile,                               &
                             keyword      = 'HDF_FIELD_NAME',                       &
                             ClientModule = 'ModuleVegetation',                     &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'ConstructVegetationList - ModuleVegetation - ERR52'                                            

                select case (trim(adjustl(FileInTime)))
                    case ("Hdf",        "HDF",          "hdf")
                        
                        
                        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
                        
                        !Open HDF
                        ObjHDF5 = 0   
                        call ConstructHDF5      (ObjHDF5, trim(FileName), HDF5_READ, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationList - ModuleVegetation - ERR60'                        
                        
                        call GetHDF5GroupNumberOfItems (ObjHDF5, "/Time",     &
                                                        NumberOfInstants,     &
                                                        STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                            & 
                            stop 'ConstructVegetationList - ModuleVegetation - ERR70'

                        ! Sets Limits for the next operations
                        call HDF5SetLimits  (ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                                             Me%WorkSize%JLB, Me%WorkSize%JUB,                       &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                                   &
                            stop 'ConstructVegetationList - ModuleVegetation - ERR75'
                        
                        !Take every grid from hdf and analyse
                        do HDFInstant = 1, NumberOfInstants 
                            
                            !Read every instant and save grid
                            call HDF5ReadData   (ObjHDF5, "/Results/"//field_name,                       &
                                                 field_name,                                             &
                                                 Array2D = AgricPractIDMatrix,                           &
                                                 OutputNumber = HDFInstant,                              &
                                                 STAT    = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_)                                                   &
                                stop 'ConstructVegetationList - ModuleVegetation - ERR80'
                            
                            !and check for any new agricultural practices in the grid 
                            call CheckVegetationListMatrix(AgricPractIDMatrix)
                            
                        enddo
                        
                    case ("Timeserie",  "TIMESERIE",    "timeserie",    "TimeSerie")

                        ObjTimeSerie = 0
                        !Constructs TimeSerie
                        call StartTimeSerieInput(ObjTimeSerie,                              &
                                                 FileName,                                  &
                                                 STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'ConstructVegetationList - ModuleVegetation - ERR90'                        
                        
                        !Get all timeserie
                        call GetTimeSerieDataMatrix(ObjTimeSerie,                           &
                                                    DataMatrix,                             &
                                                    STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'ConstructVegetationList - ModuleVegetation - ERR100' 
                                                       
                        !Get number of timeserie instants (lines)
                        call GetTimeSerieDataValues(ObjTimeSerie,                          &
                                                    DataValues,                            & 
                                                    STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'ConstructVegetationList - ModuleVegetation - ERR110' 
                        
                        !Analyse every instant to check for agricultural practices
                        do TimeSerieInstant = 1, DataValues
                            !Get ID that will be the same for every domain
                            AgricPractIDScalar  = NINT(DataMatrix (TimeSerieInstant, column))
                            
                            !and check if that agricultural practice is already present   
                            call CheckVegetationList(AgricPractIDScalar)
                            
                        enddo
                                                        
                    case default
                        write(*,*) 'Option not found in agricultural practices defined in file'
                        write(*,*) 'Check FILE_IN_TIME keyword in block Agricultural Practices'
                        write(*,*) 'Options allowed are TIMESERIE and HDF'
                        stop 'ConstructVegetationList - ModuleVegetation - ERR120' 
                end select

                
                call KillFillMatrix(Me%AgricPract%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructVegetationList - ModuleVegetation - ERR30'
            
            endif

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationList - ModuleVegetation - ERR70'  
            
            deallocate (AgricPractIDMatrix)
        
        else
            stop 'ConstructVegetationList - ModuleVegetation - ERR130'
        endif            
        
    end subroutine ConstructVegetationList

   !--------------------------------------------------------------------------
   
    subroutine CheckVegetationListMatrix (AgricPractID)
        
        !Arguments-----------------------------------------------------------
        real, dimension(:,:), pointer        :: AgricPractID
        !Local----------------------------------------------------------------
        type (T_VegetationType), pointer        :: VegetationX, VegetationInList
        logical                                 :: FoundVegetation 
        integer                                 :: i, j    
        !Begin----------------------------------------------------------------
        
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                !Found boundary or negative values inside the domain
                if (AgricPractID(i,j) < 0) then
                    write(*,*)'Practice ID negative inside the domain.'
                    write(*,*)'Check in practices grid the ID:', NINT(AgricPractID(i,j))
                    stop 'CheckVegetationListMatrix - ModuleVegetation - ERR010'
                endif
                
                !Search for the agricultural practices ID in the vegetation list
                !if already exists continue, else add to the list.
                FoundVegetation = .false.
                VegetationInList => Me%FirstVegetation
doV:            do while (associated(VegetationInList))

                    if (VegetationInList%ID == NINT(AgricPractID(i,j))) then
                        FoundVegetation = .true.
                        exit doV
                    endif
                    
                    VegetationInList => VegetationInList%Next
                
                enddo doV                
                
                if (.not. FoundVegetation) then
                    
                    allocate (VegetationX)
                    nullify  (VegetationX%Prev,VegetationX%Next)
                    
                    VegetationX%ID           = NINT(AgricPractID(i,j))     !vegetation type ID
                   
                    call AddVegetation(VegetationX)

                endif
                
            endif
        enddo
        enddo

    end subroutine CheckVegetationListMatrix

   !--------------------------------------------------------------------------    


    subroutine CheckVegetationList (AgricPractID)
        
        !Arguments-----------------------------------------------------------
        integer                                 :: AgricPractID
        !Local----------------------------------------------------------------
        type (T_VegetationType), pointer        :: VegetationX, VegetationInList
        logical                                 :: FoundVegetation     
        !Begin----------------------------------------------------------------

        !Found boundary or negative values inside the domain
        if (AgricPractID.lt.0.) then
            write(*,*)'Vegetation ID negative inside the domain.'
            write(*,*)'Check in vegetation grid provided the ID:', AgricPractID
            stop 'CheckVegetationList - ModuleVegetation - ERR10'                
        endif    

        FoundVegetation = .false.

        VegetationInList => Me%FirstVegetation
doV:    do while (associated(VegetationInList))

            if (VegetationInList%ID == AgricPractID) then
                FoundVegetation = .true.
                exit doV
            endif
            
            VegetationInList => VegetationInList%Next

        enddo doV                

        if (.not. FoundVegetation) then
            allocate (VegetationX)
            nullify  (VegetationX%Prev,VegetationX%Next)
            VegetationX%ID           = AgricPractID     !vegetation type ID
!           VegetationX%VegetationID = AuxID(i,j)     !crop ID growing - for now is the same
           
            call AddVegetation(VegetationX)
        
        endif
                

    end subroutine CheckVegetationList

   !--------------------------------------------------------------------------  
    
    subroutine CheckOptionsConsistence

        !Local----------------------------------------------------------------
        integer                                 :: STAT_CALL
        type(T_Property), pointer               :: TotalPlantBiomass, TotalPlantNitrogen
        type(T_Property), pointer               :: TotalPlantPhosphorus, RootBiomass, RootDepth
        type(T_Property), pointer               :: LeafAreaIndex, SpecificLeafStorage, CanopyHeight
        type(T_Property), pointer               :: EVTPCropCoefficient
        logical                                 :: NitrogenFound, PhosphorusFound, RootBiomassFound
        logical                                 :: Grazing, HarvestKill, Dormancy
        !Begin----------------------------------------------------------------
        
        if(Me%ComputeOptions%Evolution%ModelSWAT .and. Me%ComputeOptions%Evolution%ModelDEB) then
            write(*,*    ) 
            write(*,*    ) 'Fatal error ! Properties Evolution can not be modeled' 
            write(*,*    ) 'by SWAT and DEB at the same time. Check EVOLUTION keyword in' 
            write(*,*    ) 'properties block.'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR000'
        endif
        
        !All fluxes in SWAT are computed in daily basis. to change this all fluxes have to be changed
        if ((Me%ComputeOptions%Evolution%ModelSWAT) .and. (Me%ComputeOptions%VegetationDT .ne. 86400.)) then
            write(*,*)
            write(*,*)'Using growth model is obligatory to use'
            write(*,*)'vegetation daily step (86400 s) Check VEGETATION_DT keyword'
            write(*,*)
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR001'
        endif
        
        
        call GetComputeSoilField(Me%ObjPorousMedia, Me%ExternalVar%ComputeSoilField, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsConsistence - ModuleVegetation - ERR05'   
        
        if (Me%ComputeOptions%ModelWater) then
        
            if(Me%ComputeOptions%TranspirationMethod .eq. TranspirationSWAT .and. .not. Me%ExternalVar%ComputeSoilField) then
                write(*,*    ) 
                write(*,*    ) 'Fatal error ! User defined method for water uptake implies soil field capacity' 
                write(*,*    ) 'computation. Check keyword COMPUTE_SOIL_FIELD in porous media file' 
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR010'
            endif
            if(.not. Me%ComputeOptions%Evolution%GrowthModelNeeded .and.                              &
                (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus)) then
                    if (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptakeSWAT) then
                        write(*,*    )
                        write(*,*    ) 'If growth model disconnected cant use that nutrient uptake option'
                        write(*,*    ) 'Check NUTRIENT_UPTAKE_METHOD keyword in vegetation data file'
                        stop 'CheckOptionsConsistence - ModuleVegetation - ERR015'
                    endif
            endif  
            
!            if (Me%ComputeOptions%Evolution%ModelSWAT .and.                                           &
!                (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptake_TranspConc) .and.          &
!                (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressSWAT)) then
!                write(*,*    )
!                write(*,*    ) 'If growth model connected cant use nutrient uptake Q*C option'
!                write(*,*    ) 'and compute stress based on SWAT formulation'
!                write(*,*    ) 'Check NUTRIENT_UPTAKE_METHOD and NUTRIENT_STRESS_METHOD keyword in vegetation data file'
!                stop 'CheckOptionsConsistence - ModuleVegetation - ERR016'        
!            endif
                          
        else
!            if (.not. Me%ComputeOptions%Evolution%GrowthModelNeeded) then
!                write(*,*    ) 
!                write(*,*    ) 'Fatal error ! Water stress and growth model disconnected' 
!                write(*,*    ) 'To use vegetation readed from file can not disconnect water.'
!                write(*,*    ) 'Check WATER_STRESS keyword in vegetation data file' 
!                stop 'CheckOptionsConsistence - ModuleVegetation - ERR020'
!            endif        
             if (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
                if (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptake_TranspConc) then
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Water uptake is disconnected and nutrient uptake method' 
                    write(*,*    ) 'is dependent on water uptake. Check WATER_STRESS and NUTRIENT_UPTAKE_METHOD'
                    write(*,*    ) 'keywords in vegetation data file' 
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR020.1'
                endif  
            endif              
        endif
        
        if (.not. Me%ComputeOptions%Evolution%GrowthModelNeeded) then
            
            if(Me%ComputeOptions%HarvestKill .or. Me%ComputeOptions%Grazing .or. Me%ComputeOptions%Dormancy         &
!               .or. Me%ComputeOptions%Fertilization .or. Me%ComputeOptions%NutrientFluxesWithSoil                   &
               .or. Me%ComputeOptions%Fertilization                                                                  &
               .or. Me%ComputeOptions%Pesticide) then
                
                write(*,*    ) 
                write(*,*    ) 'Fatal error ! Vegetation growth model is disconnected' 
                write(*,*    ) 'but trying to simulate HarvestKill or grazing or dormancy or fertilization' 
                write(*,*    ) 'or fluxes with soil. Check this processes keywords in vegetation data file' 
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR020.5'
            endif
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
            write(*,*    ) 'vegetation'
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

            if(Me%ComputeOptions%ModelNitrogen) then

                !Optional Property
                call SearchProperty(TotalPlantNitrogen, TotalPlantNitrogen_, .false., STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    if(TotalPlantNitrogen%Evolution == ReadValue) then
                        write(*,*    ) 
                        write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                        write(*,*    ) 'total plant nitrogen must be modeled in this conditions.'
                        write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                        stop 'CheckOptionsConsistence - ModuleVegetation - ERR050'
                    endif
                else
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Nitrogen was selected to run but' 
                    write(*,*    ) 'total plant nitrogen property was not found.'
                    write(*,*    ) 'Check NITROGEN_STRESS keyword or property block.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR051'

                endif
            else
                !Optional Property
                call SearchProperty(TotalPlantNitrogen, TotalPlantNitrogen_, .false., STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    write(*,*    ) 
                    write(*,*    ) 'total plant nitrogen property was found but nitrogen is not to' 
                    write(*,*    ) 'to be computed. Check NITROGEN_STRESS or delete the property'
                    write(*,*    ) 'in order to produce a consistent input.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR052'
                endif           
            endif


            if(Me%ComputeOptions%ModelPhosphorus) then
        
                !Optional Property
                call SearchProperty(TotalPlantPhosphorus, TotalPlantPhosphorus_, .false., STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    if(TotalPlantPhosphorus%Evolution == ReadValue) then
                        write(*,*    ) 
                        write(*,*    ) 'Fatal error ! Vegetation growth model was selected to run and' 
                        write(*,*    ) 'total plant phosphorus must be modeled in this conditions.'
                        write(*,*    ) 'Check EVOLUTION keyword consistency in all properties blocks.'
                        stop 'CheckOptionsConsistence - ModuleVegetation - ERR060'
                    endif           
                else
                    write(*,*    ) 
                    write(*,*    ) 'Fatal error ! Phosphorus was selected to run but' 
                    write(*,*    ) 'total plant phosphorus property was not found.'
                    write(*,*    ) 'Check PHOSPHORUS_STRESS keyword or property block.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR061'

                endif
            else
                call SearchProperty(TotalPlantPhosphorus, TotalPlantPhosphorus_, .false., STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    write(*,*    ) 
                    write(*,*    ) 'total plant phosphorus property was found but phosphorus is not to' 
                    write(*,*    ) 'to be computed. Check PHOSPHORUS_STRESS or delete the property'
                    write(*,*    ) 'in order to produce a consistent input.'
                    stop 'CheckOptionsConsistence - ModuleVegetation - ERR062'
                endif
            endif

            !Optional
            HarvestKill       = Me%ComputeOptions%HarvestKill

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
                endif
            elseif (STAT_CALL /= SUCCESS_ .and. HarvestKill) then
                write(*,*    )
                write(*,*    ) 'root biomass is a mandatory Property'
                write(*,*    ) 'if you want to run vegetation model and want to model'
                write(*,*    ) 'HarvestKill practices (e.g. hatvest where aerial biomass must be known'
                stop 'CheckOptionsConsistence - ModuleVegetation - ERR071'

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
                                                .or. HarvestKill .or. Dormancy))  then
                write(*,*    )
                write(*,*    ) 'total plant biomass is a mandatory Property'
                write(*,*    ) 'if you want to run vegetation model and want to model'
                write(*,*    ) 'nitrogen, phosphorus, root biomass, grazing or HarvestKill practices'
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
        
        if(Me%ComputeOptions%Fertilization .and. .not. Me%ComputeOptions%ModelNitrogen                &
           .and. .not. Me%ComputeOptions%ModelPhosphorus) then 
            write(*,*) 
            write(*,*) 'Fertilization is active but plant nitrogen and phosphorus properties not'
            write(*,*) 'Check this options in vegetation data file'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR131'
        endif

        !To do work
        if (Me%ComputeOptions%Grazing .and. Me%ComputeOptions%Continuous .and. Me%ComputeOptions%StopOnWrongDate) then
            write(*,*) 'WARNING'
            write(*,*) 'Grazing operations happening in the end of the last run'
            write(*,*) ' will be descontinued in the next run'
            stop 'CheckOptionsConsistence - ModuleVegetation - ERR140'
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
            
            if(Me%ComputeOptions%ModelWater) then
                write(*,*    ) '---Water Uptake             : ON'
            else
                write(*,*    ) '---Water Uptake             : OFF'
            endif
            
            if(Me%ComputeOptions%ModelWater) then
                if (Me%ComputeOptions%TranspirationMethod .eq. 1) then
                    if (Me%ComputeOptions%WaterStress) then
                        write(*,*    ) '   ---WaterUptake Method    : MOHID based (Feddes ON)'
                    else
                        write(*,*    ) '   ---WaterUptake Method    : MOHID based (Feddes OFF)'
                    endif
                else
                    write(*,*    ) '   ---WaterUptake Method    : SWAT based'
                endif
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
            if(Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
                select case (Me%ComputeOptions%NutrientUptakeMethod)
                    case (1)
                        write(*,*    ) '   ---NutrientUptake Method : Q*Conc'
                    case (2)
                        write(*,*    ) '   ---NutrientUptake Method : SWAT based'
 !                       write(*,*    ) '       WARNING - This option gives mass errors'
                end select
                
                if (Me%ComputeOptions%NutrientStressMethod .eq. 1) then
                    write(*,*    ) '   ---NutrientStress Method : Actual/Pot'
                else
                    write(*,*    ) '   ---NutrientStress Method : SWAT based'
                endif
            endif
            if (Me%ComputeOptions%ModelTemperatureStress) then
                write(*,*    ) '---Temperature Stress       : ON'
            else
                write(*,*    ) '---Temperature Stress       : OFF'
            endif
            if (Me%ComputeOptions%HarvestKill) then
                write(*,*    ) '---HarvestKill Operations   : ON'
            else
                write(*,*    ) '---HarvestKill Operations   : OFF'
            endif              
            if (Me%ComputeOptions%Grazing) then
                write(*,*    ) '---Grazing Operations       : ON'
            else
                write(*,*    ) '---Grazing Operations       : OFF'
            endif                                   
            if (Me%ComputeOptions%Dormancy) then
                write(*,*    ) '---Dormancy                 : ON'
            else
                write(*,*    ) '---Dormancy                 : OFF'
            endif                                   
            if (Me%ComputeOptions%Fertilization) then
                write(*,*    ) '---Fertilization            : ON'
            else
                write(*,*    ) '---Fertilization            : OFF'
            endif                      
            if (Me%ComputeOptions%Pesticide) then
                write(*,*    ) '---Pesticide                : ON'
            else
                write(*,*    ) '---Pesticide                : OFF'
            endif                                  
            if (Me%ComputeOptions%NutrientFluxesWithSoil) then
                write(*,*    ) '---NutrientFluxesWithSoil   : ON'
                write(*,*    ) 
            else
                write(*,*    ) '---NutrientFluxesWithSoil   : OFF'
                write(*,*    ) 
            endif                      
        
        else

            write(*,*    ) 'Vegetation Growth Model not Used'
            write(*,*    ) '---Root readed from file'
            write(*,*    ) '---LAI  readed from file'
            if (Me%ComputeOptions%ModelWater) then
                write(*,*    ) '---Water Uptake             : ON'
                if (Me%ComputeOptions%TranspirationMethod .eq. 1) then
                    if (Me%ComputeOptions%WaterStress) then
                        write(*,*    ) '   ---WaterUptake Method    : MOHID based (Feddes ON)'
                    else
                        write(*,*    ) '   ---WaterUptake Method    : MOHID based (Feddes OFF)'
                    endif                
                else
                    write(*,*    ) '   ---WaterUptakeMethod     : SWAT based'
                endif
            else
                write(*,*    ) '---Water Uptake             : OFF'
            endif
            write(*,*    ) '---Biomass Growth           : OFF'
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
            if(Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
                select case (Me%ComputeOptions%NutrientUptakeMethod)
                    case (1)
                        write(*,*    ) '   ---NutrientUptake Method : 1'
                    case (2)
                        write(*,*    ) '   ---NutrientUptake Method : 2'
                    case (3)
                        write(*,*    ) '   ---NutrientUptake Method : 3'
                end select
            endif
            write(*,*    ) '---Temperature Stress       : OFF'
            write(*,*    ) '---Grazing, HarvestKill'       
            write(*,*    ) '   Dormancy & Fertilization : OFF'
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
    
    subroutine CreatePropertiesAlias
    

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_property), pointer                          :: PropertyX
        integer                                             :: STAT_CALL
        
        !Begin-----------------------------------------------------------------    
        call SearchProperty(PropertyX  , TotalPlantBiomass_    , .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%TotalPlantBiomass => PropertyX%Field        
        call SearchProperty(PropertyX  , TotalPlantNitrogen_   , .false., STAT = STAT_CALL)        
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%TotalPlantNitrogen => PropertyX%Field
        call SearchProperty(PropertyX  , TotalPlantPhosphorus_ , .false., STAT = STAT_CALL)        
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%TotalPlantPhosphorus => PropertyX%Field
        call SearchProperty(PropertyX  , RootBiomass_          , .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%RootBiomass => PropertyX%Field
        call SearchProperty(PropertyX  , RootDepth_            , .false., STAT = STAT_CALL)        
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%RootDepth => PropertyX%Field
        call SearchProperty(PropertyX  , LeafAreaIndex_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%LeafAreaIndex => PropertyX%Field        
        call SearchProperty(PropertyX  , CanopyHeight_         , .false., STAT = STAT_CALL)         
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%CanopyHeight => PropertyX%Field
        call SearchProperty(PropertyX  , SpecificLeafStorage_  , .false., STAT = STAT_CALL)        
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%SpecificLeafStorage => PropertyX%Field
        call SearchProperty(PropertyX  , EVTPCropCoefficient_  , .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)       Me%StateVariables%EVTPCropCoefficient => PropertyX%Field
        call SearchProperty(PropertyX  , PotLeafAreaIndex_     , .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            Me%StateVariables%PotLeafAreaIndex => PropertyX%Field  
        else
            Me%StateVariables%PotLeafAreaIndex => null()
        endif   
        
        !Actualize these state variables using property values
        call SetMatrixValue(Me%RootDepthOld, Me%Size2D,  Me%StateVariables%RootDepth)   
        
        if (Me%ComputeOptions%Evolution%ModelSWAT) then 
            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                call SetMatrixValue(Me%LAIBeforeSenescence,           &
                                    Me%Size2D,                        &
                                    Me%StateVariables%LeafAreaIndex,  &
                                    Me%ExternalVar%MappingPoints)           
            endif    
        endif
    end subroutine CreatePropertiesAlias    
    
    !--------------------------------------------------------------------------
    
    subroutine AllocateVariables
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB
        integer                                             :: KLB, KUB                                  
        
        !Begin-----------------------------------------------------------------    
    
        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB  
        
        allocate(Me%VegetationID (ILB:IUB,JLB:JUB))
        Me%VegetationID (:,:) = FillValueInt

        allocate(Me%AgricPractID (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%AgricPractID (:,:) = FillValueInt
        
        allocate(Me%RootDepthOld (ILB:IUB,JLB:JUB))
        Me%RootDepthOld (:,:) = 0.0

        !Soil Fluxes
        allocate(Me%Fluxes%WaterUptake (ILB:IUB,JLB:JUB))  
        Me%Fluxes%WaterUptake (:,:) = 0.0 
        
        allocate(Me%Fluxes%WaterUptakeLayer (ILB:IUB,JLB:JUB,KLB:KUB))  
        Me%Fluxes%WaterUptakeLayer (:,:,:) = 0.0 
                
        allocate(Me%Growth%WaterStress (ILB:IUB,JLB:JUB)) 
        Me%Growth%WaterStress (:,:) = 1.0
        
        allocate(Me%TranspirationBottomLayer (ILB:IUB,JLB:JUB))        
        if (Me%ComputeOptions%ModelWater) then
            Me%TranspirationBottomLayer = KLB
        else
            Me%TranspirationBottomLayer = FillValueInt
        endif
        
        allocate(Me%Growth%TemperatureStress (ILB:IUB,JLB:JUB)) 
        Me%Growth%TemperatureStress (:,:) = 1.0
        
        allocate(Me%Growth%NitrogenStress (ILB:IUB,JLB:JUB)) 
        Me%Growth%NitrogenStress (:,:) = 1.0
        
        allocate(Me%Growth%PhosphorusStress (ILB:IUB,JLB:JUB))                 
        Me%Growth%PhosphorusStress (:,:) = 1.0
            
        if (Me%ComputeOptions%ModelNitrogen) then
            allocate(Me%Fluxes%NitrogenUptake (ILB:IUB,JLB:JUB))
            Me%Fluxes%NitrogenUptake (:,:) = 0.0 
            
            allocate(Me%Fluxes%NitrogenUptakeLayer (ILB:IUB,JLB:JUB,KLB:KUB))              
            Me%Fluxes%NitrogenUptakeLayer (:,:,:) = 0.0 
        endif

        if (Me%ComputeOptions%ModelPhosphorus) then
            allocate(Me%Fluxes%PhosphorusUptake (ILB:IUB,JLB:JUB))
            Me%Fluxes%PhosphorusUptake (:,:) = 0.0 
            
            allocate(Me%Fluxes%PhosphorusUptakeLayer (ILB:IUB,JLB:JUB,KLB:KUB))              
            Me%Fluxes%PhosphorusUptakeLayer (:,:,:) = 0.0 
        endif   
        
        !Nutrient uptake fluxes in porous media properties
        if (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
            allocate(Me%SoilFluxesActive (ILB:IUB,JLB:JUB))
            Me%SoilFluxesActive (:,:) = .false.
        endif

        allocate(Me%ExternalVar%Integration%SumPotTP (ILB:IUB,JLB:JUB))
        Me%ExternalVar%Integration%SumPotTP (:,:) = 0.0                     
        
        allocate(Me%ExternalVar%Integration%AveragePotTPDuringDT (ILB:IUB,JLB:JUB))
        Me%ExternalVar%Integration%AveragePotTPDuringDT (:,:) = 0.0
    
    end subroutine AllocateVariables
    
    !--------------------------------------------------------------------------

    subroutine AllocateEvolutionVariables

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB
        integer                                             :: KLB, KUB, JulDay                 
        real                                                :: InitialHUBase        
        integer                                             :: Pest
        !Begin-----------------------------------------------------------------

        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
               
        if (Me%ComputeOptions%Evolution%ModelSWAT) then                        
            call JulianDay(Me%ExternalVar%Now, JulDay)
            Me%ExternalVar%JulianDay_Old = JulDay
            
            allocate(Me%PlantingOccurred (ILB:IUB,JLB:JUB))
            Me%PlantingOccurred = .false.            

            allocate(Me%IsPlantGrowing (ILB:IUB,JLB:JUB))
            Me%IsPlantGrowing (:,:) = .false.                                            
            
            allocate(Me%ExternalVar%Integration%SumTemperature (ILB:IUB,JLB:JUB))  
            Me%ExternalVar%Integration%SumTemperature (:,:) = 0.0
            
            allocate(Me%ExternalVar%Integration%SumHumidity (ILB:IUB,JLB:JUB)) 
            Me%ExternalVar%Integration%SumHumidity (:,:) = 0.0  
            
            allocate(Me%ExternalVar%Integration%SumRadiation (ILB:IUB,JLB:JUB))
            Me%ExternalVar%Integration%SumRadiation (:,:) = 0.0
                             
            allocate(Me%HeatUnits%PotentialHUTotal (ILB:IUB,JLB:JUB)) 
            Me%HeatUnits%PotentialHUTotal (:,:) = 0.0           
            
            allocate(Me%HeatUnits%PotentialHUBase_Old (ILB:IUB,JLB:JUB)) 
            Me%HeatUnits%PotentialHUBase_Old(:,:) = 0.0 
            
            allocate(Me%HeatUnits%PlantHUAccumulated (ILB:IUB,JLB:JUB)) 
            Me%HeatUnits%PlantHUAccumulated (:,:) = 0.0
                         
            allocate(Me%HeatUnits%PlantHUAccumulated_Old (ILB:IUB,JLB:JUB)) 
            Me%HeatUnits%PlantHUAccumulated_Old (:,:) = 0.0
            
            allocate(Me%ExternalVar%Integration%AverageAirTempDuringDT (ILB:IUB,JLB:JUB))
            Me%ExternalVar%Integration%AverageAirTempDuringDT (:,:) = 0.0
              
            allocate(Me%ExternalVar%Integration%AverageAirHumidityDuringDT(ILB:IUB,JLB:JUB)) 
            Me%ExternalVar%Integration%AverageAirHumidityDuringDT (:,:) = 0.0
            
            allocate(Me%ExternalVar%Integration%AverageRadiationDuringDT (ILB:IUB,JLB:JUB))                                                
            Me%ExternalVar%Integration%AverageRadiationDuringDT (:,:) = 0.0
            
            !Construct fraction of HU at the beggining of simulation
            !Estimated based on fraction of the year - assumes that HU are equally distributed in the year
            !just to get one idea of where in the year (in terms of HU) the simulation starts
            !HU base is used to see when is time to plant
            allocate(Me%HeatUnits%PotentialHUBase (ILB:IUB,JLB:JUB)) 
            call JulianDay(Me%ExternalVar%Now, JulDay)
            InitialHUBase = JulDay / 365.
            Me%HeatUnits%PotentialHUBase (:,:) = InitialHUBase                                                         

            allocate(Me%Growth%GlobalStress (ILB:IUB,JLB:JUB))
            Me%Growth%GlobalStress (:,:) = 1.0

            allocate(Me%Growth%TreeCurrentYear (ILB:IUB,JLB:JUB))
            Me%Growth%TreeCurrentYear (:,:) = 0
            
            allocate(Me%Growth%TreeComingFromContinuous (ILB:IUB,JLB:JUB))       
            Me%Growth%TreeComingFromContinuous (:,:) = .false.
            
            allocate(Me%Growth%TreeFractionToMaturity (ILB:IUB,JLB:JUB))        
            Me%Growth%TreeFractionToMaturity (:,:) = 0
            
            allocate(Me%Growth%TreeMaximumAnnualBiomass (ILB:IUB,JLB:JUB))
            Me%Growth%TreeMaximumAnnualBiomass (:,:) = 0
          
            !Soil fluxes
            if(Me%ComputeOptions%ModelNitrogen) then
                allocate(Me%PlantNitrogenFraction (ILB:IUB,JLB:JUB)) 
                Me%PlantNitrogenFraction(:,:) = 0.0
                
                allocate(Me%OptimalTotalPlantNitrogen (ILB:IUB,JLB:JUB)) 
                Me%OptimalTotalPlantNitrogen (:,:) = 0.0 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                allocate(Me%PlantPhosphorusFraction (ILB:IUB,JLB:JUB))
                Me%PlantPhosphorusFraction(:,:) = 0.0
                
                allocate(Me%OptimalTotalPlantPhosphorus (ILB:IUB,JLB:JUB))
                Me%OptimalTotalPlantPhosphorus (:,:) = 0.0                 
            endif
            
            !Aerial Fluxes
            allocate(Me%Fluxes%LAIChange (ILB:IUB,JLB:JUB)) 
            Me%Fluxes%LAIChange (:,:) = 0.0
            
            allocate(Me%LAISenescence (ILB:IUB,JLB:JUB)) 
            Me%LAISenescence(:,:) = .false.  
                      
            allocate(Me%PlantLAIMaxFraction (ILB:IUB,JLB:JUB)) 
            Me%PlantLAIMaxFraction (:,:) = 0.0

            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                allocate(Me%LAIBeforeSenescence (ILB:IUB,JLB:JUB)) 
                Me%LAIBeforeSenescence (:,:) = 0.0
            endif


        endif
        
        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

            if (Me%ComputeOptions%Grazing) then
                allocate(Me%DaysOfGrazing (ILB:IUB,JLB:JUB)) 
                Me%DaysOfGrazing (:,:) = 0.0 
                
                allocate(Me%IsPlantBeingGrazed (ILB:IUB,JLB:JUB))                 
                Me%IsPlantBeingGrazed (:,:) = .false.
            
                allocate(Me%Fluxes%BiomassGrazed (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%BiomassGrazed (:,:) = 0.0
                
                allocate(Me%Fluxes%BiomassGrazedFraction (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%BiomassGrazedFraction (:,:) = 0.0
                
                allocate(Me%Fluxes%ToSoil%GrazingBiomassToSoil (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%ToSoil%GrazingBiomassToSoil (:,:) = 0.0
                
                if (Me%ComputeOptions%ModelNitrogen) then
                    allocate(Me%Fluxes%NitrogenGrazed (ILB:IUB,JLB:JUB))  
                    Me%Fluxes%NitrogenGrazed (:,:) = 0.0
                    
                    allocate(Me%Fluxes%ToSoil%GrazingNitrogenToSoil (ILB:IUB,JLB:JUB)) 
                    Me%Fluxes%ToSoil%GrazingNitrogenToSoil (:,:) = 0.0
                endif
            
                if (Me%ComputeOptions%ModelPhosphorus) then
                    allocate(Me%Fluxes%PhosphorusGrazed (ILB:IUB,JLB:JUB)) 
                    Me%Fluxes%PhosphorusGrazed (:,:) = 0.0
                    
                    allocate(Me%Fluxes%ToSoil%GrazingPhosphorusToSoil (ILB:IUB,JLB:JUB)) 
                    Me%Fluxes%ToSoil%GrazingPhosphorusToSoil (:,:) = 0.0
                endif

            endif
        
            if (Me%ComputeOptions%HarvestKill) then
                allocate(Me%HarvestOnlyOccurred (ILB:IUB,JLB:JUB))
                Me%HarvestOnlyOccurred (:,:) = .false.

                allocate(Me%HarvestKillOccurred (ILB:IUB,JLB:JUB))
                Me%HarvestKillOccurred (:,:) = .false.

                allocate(Me%KillOccurred (ILB:IUB,JLB:JUB))
                Me%KillOccurred (:,:) = .false.
                      
                allocate(Me%Fluxes%BiomassRemovedInHarvest (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%BiomassRemovedInHarvest (:,:) = 0.0
                
                allocate(Me%Fluxes%BiomassHarvestedFraction (ILB:IUB,JLB:JUB))
                Me%Fluxes%BiomassHarvestedFraction (:,:) = 0.0
                 
                allocate(Me%Fluxes%ToSoil%HarvestKillBiomassToSoil (ILB:IUB,JLB:JUB))
                Me%Fluxes%ToSoil%HarvestKillBiomassToSoil (:,:) = 0.0
                 
                allocate(Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil (ILB:IUB,JLB:JUB))
                Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil (:,:) = 0.0
            
                if (Me%ComputeOptions%ModelNitrogen) then            
                    allocate(Me%Fluxes%NitrogenRemovedInHarvest (ILB:IUB,JLB:JUB))  
                    Me%Fluxes%NitrogenRemovedInHarvest (:,:) = 0.0
                    
                    allocate(Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil (ILB:IUB,JLB:JUB)) 
                    Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil (:,:) = 0.0
                endif
            
                if (Me%ComputeOptions%ModelPhosphorus) then
                    allocate(Me%Fluxes%PhosphorusRemovedInHarvest (ILB:IUB,JLB:JUB)) 
                    Me%Fluxes%PhosphorusRemovedInHarvest (:,:) = 0.0
                    
                    allocate(Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil (ILB:IUB,JLB:JUB))
                    Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil (:,:) = 0.0
                endif

            endif

            if (Me%ComputeOptions%Dormancy) then
                allocate(Me%DayLength (ILB:IUB,JLB:JUB))
                Me%DayLength (:,:) = 0.0
                
                allocate(Me%MinimumDayLength (ILB:IUB,JLB:JUB))
                Me%MinimumDayLength (:,:) = 0.0
                
                allocate(Me%IsPlantDormant (ILB:IUB,JLB:JUB)) 
                Me%IsPlantDormant (:,:) = .false.
                
                allocate(Me%PlantGoingDormant (ILB:IUB,JLB:JUB))
                Me%PlantGoingDormant (:,:) = .false.
            
                allocate(Me%Fluxes%BiomassRemovedInDormancy (ILB:IUB,JLB:JUB))
                Me%Fluxes%BiomassRemovedInDormancy (:,:) = 0.0
            
                if (Me%ComputeOptions%ModelNitrogen) then
                    allocate(Me%Fluxes%NitrogenRemovedInDormancy (ILB:IUB,JLB:JUB))
                    Me%Fluxes%NitrogenRemovedInDormancy (:,:) = 0.0  
                endif
            
                if (Me%ComputeOptions%ModelPhosphorus) then
                    allocate(Me%Fluxes%PhosphorusRemovedInDormancy (ILB:IUB,JLB:JUB)) 
                    Me%Fluxes%PhosphorusRemovedInDormancy (:,:) = 0.0 
                endif
 
            endif
        
            if (Me%ComputeOptions%Fertilization) then
                if (Me%ComputeOptions%AutoFertilization) then
                    allocate(Me%NitrogenYeldEstimate (ILB:IUB,JLB:JUB))
                    Me%NitrogenYeldEstimate (:,:) = 0.0
                    
                    allocate(Me%NitrogenYeldTarget (ILB:IUB,JLB:JUB))
                    Me%NitrogenYeldTarget (:,:) = 0.0
                endif
                
                if (Me%ComputeOptions%ModelNitrogen) then
                    allocate(Me%AnnualNitrogenFertilized (ILB:IUB,JLB:JUB))
                    Me%AnnualNitrogenFertilized (:,:) = 0.0 
                    
                    allocate(Me%Fluxes%FertilNitrateInSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilNitrateInSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilNitrateInSubSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilNitrateInSubSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilAmmoniaInSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilAmmoniaInSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilAmmoniaInSubSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilAmmoniaInSubSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilOrganicNInSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilOrganicNInSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilOrganicNParticInFluff (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilOrganicNParticInFluff (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilOrganicNInSubSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilOrganicNInSubSurface (:,:) = 0.0
                endif
                if (Me%ComputeOptions%ModelPhosphorus) then
                    allocate(Me%AnnualPhosphorusFertilized (ILB:IUB,JLB:JUB))
                    Me%AnnualPhosphorusFertilized (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilOrganicPInSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilOrganicPInSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilOrganicPInSubSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilOrganicPInSubSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilOrganicPParticInFluff (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilOrganicPParticInFluff (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilMineralPInSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilMineralPInSurface (:,:) = 0.0
                    
                    allocate(Me%Fluxes%FertilMineralPInSubSurface (ILB:IUB,JLB:JUB))
                    Me%Fluxes%FertilMineralPInSubSurface (:,:) = 0.0
                endif

            endif

            if (Me%ComputeOptions%Pesticide) then
                do  Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
                    allocate (Me%Fluxes%Pesticides%Application(Pest)%Soil (ILB:IUB,JLB:JUB))
                    Me%Fluxes%Pesticides%Application(Pest)%Soil (:,:) = 0.0
                    
                    allocate (Me%Fluxes%Pesticides%Application(Pest)%Vegetation (ILB:IUB,JLB:JUB))                        
                    Me%Fluxes%Pesticides%Application(Pest)%Vegetation (:,:) = 0.0
                enddo                            
            endif                              
        endif

        if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then
            allocate(Me%ComputeOptions%TranspirationMOHID%FeddesType (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%FeddesType (:,:) = 0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH1 (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesH1 (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH2 (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesH2 (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3 (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesH3 (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH4 (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesH4 (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3L (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesH3L (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3H (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesH3H (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesR2L (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesR2L (:,:) = 0.0
            
            allocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesR2H (ILB:IUB,JLB:JUB))
            Me%ComputeOptions%TranspirationMOHID%RootFeddesR2H (:,:) = 0.0
        endif        

    end subroutine AllocateEvolutionVariables

    !--------------------------------------------------------------------------
    
    subroutine AllocatePropertyVariables
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                             :: ILB, IUB, JLB, JUB
        integer                                             :: KLB, KUB                                          
        !Begin-----------------------------------------------------------------

        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
            
        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
        
            if (Me%ComputeOptions%ModelPlantBiomass) then
                allocate(Me%Fluxes%BiomassGrowth (ILB:IUB,JLB:JUB)) 
                Me%Fluxes%BiomassGrowth (:,:) = 0.0
                
                allocate(Me%Growth%BiomassGrowthOld (ILB:IUB,JLB:JUB))  
                Me%Growth%BiomassGrowthOld (:,:) = 0.0
                
                allocate(Me%Growth%PAR (ILB:IUB,JLB:JUB))  
                Me%Growth%PAR (:,:) = 0.0
                
                allocate(Me%Growth%RUE (ILB:IUB,JLB:JUB))  
                Me%Growth%RUE (:,:) = 0.0
                
                allocate(Me%Growth%PotentialGrowth (ILB:IUB,JLB:JUB))
                Me%Growth%PotentialGrowth (:,:) = 0.0
                
                allocate(Me%Growth%PotentialBiomass (ILB:IUB,JLB:JUB))  
                Me%Growth%PotentialBiomass (:,:) = 0.0                               
            endif

            if (Me%ComputeOptions%ModelCanopyHeight) then
                allocate(Me%ChangeCanopyEnabled (ILB:IUB,JLB:JUB)) 
                Me%ChangeCanopyEnabled (:,:) = .false.
            endif 
        
        endif
    
    end subroutine AllocatePropertyVariables

    !--------------------------------------------------------------------------

    subroutine ConstructInitialConditions

    end subroutine ConstructInitialConditions

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: nProperties, STAT_CALL
        integer                                             :: nPropertiesToSoil
        integer                                             :: iflag, AddProperties, i
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer                                             :: Pest
        integer                                             :: TimeSerieNumber, dn, Id, Jd
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        character(len=StringLength)                         :: TimeSerieName
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
            if (Me%ComputeOptions%ModelNitrogen) then
                AddProperties = AddProperties + 3
            endif  
            if (Me%ComputeOptions%ModelPhosphorus) then
                AddProperties = AddProperties + 3
            endif 
            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                AddProperties = AddProperties + 3
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

            if (Me%ComputeOptions%ModelNitrogen) then
                nProperties = nProperties + 1
                PropertyList(nProperties) = "Nitrogen Uptake kg/ha"            
                nProperties = nProperties + 1
                PropertyList(nProperties) = "OptimalPlantNitrogen_kg/ha"  
                nProperties = nProperties + 1
                PropertyList(nProperties) = "NitrogenStressFactor"  
            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
                nProperties = nProperties + 1
                PropertyList(nProperties) = "Phosphorus Uptake kg/ha"
                nProperties = nProperties + 1
                PropertyList(nProperties) = "OptimalPlantPhosphorus_kg/ha"  
                nProperties = nProperties + 1
                PropertyList(nProperties) = "PhosphorusStressFactor"               
            endif

            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
            
                nProperties = nProperties + 1 
                PropertyList(nProperties) = "HU Accumulated"
                nProperties = nProperties + 1
                PropertyList(nProperties) = "Potential HU"
                nProperties = nProperties + 1
                PropertyList(nProperties) = "TemperatureStressFactor"
            
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
                                PropertyList, "srvg",                           &
                                WaterPoints2D = Me%ExternalVar%MappingPoints,   &
                                STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_Time_Serie - ModuleVegetation - ERR03'

            !Deallocates PropertyList
            deallocate(PropertyList, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_Time_Serie - ModuleVegetation - ERR04'




            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

              !!!! Time serie with averaged atmosphere properties
                call GetData(Me%ComputeOptions%AtmospherePropertiesOutput,          &
                             Me%ObjEnterData,iflag,                                 &
                             SearchType   = FromFile,                               &
                             keyword      = 'ATMOSPHERE_OUTPUT',                    &
                             ClientModule = 'ModuleVegetation',                     &
                             Default      = .false.,                                &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'Construct_Time_Serie - ModuleVegetation - ERR05' 
            

                if(Me%ComputeOptions%AtmospherePropertiesOutput) then
                    
                    if (Me%ComputeOptions%ModelPlantBiomass) then
                        allocate (PropertyList(7))
                    else
                        allocate (PropertyList(3))
                    endif

                    PropertyList(1) = "Average Air Temperature_ºC"
                    PropertyList(2) = "Average Relative Humidity_-"
                    PropertyList(3) = "Average Solar Radiation_W/m2"

                    if (Me%ComputeOptions%ModelPlantBiomass) then
                        PropertyList(4) = "PAR_(MJ/m2)"
                        PropertyList(5) = "RUE_(kg/ha)/(MJ/m2)"
                        PropertyList(6) = "Potential Growth_(Kg/ha)"
                        PropertyList(7) = "Potential Biomass_(Kg/ha)"
                    endif
            
                    call StartTimeSerie(Me%ObjTimeSerieAtm,                                         &
                                        Me%ObjTime,                                                 &
                                        TimeSerieLocationFile,                                      &
                                        PropertyList, "srvgs",                                      &
                                        WaterPoints2D = Me%ExternalVar%MappingPoints,               &
                                        STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleVegetation - ERR06'

                    deallocate(PropertyList)
                endif



                !Time serie with outputs to soil
                call GetData(Me%ComputeOptions%FluxestoSoilOutput,                  &
                             Me%ObjEnterData,iflag,                                 &
                             SearchType   = FromFile,                               &
                             keyword      = 'FLUXES_TO_SOIL_OUTPUT',                &
                             ClientModule = 'ModuleVegetation',                     &
                             Default      = .false.,                                &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                        &
                    stop 'Construct_Time_Serie - ModuleVegetation - ERR07' 
            
                if (Me%ComputeOptions%FluxestoSoilOutput) then
            
                    nPropertiesToSoil = 0
                    !Count number of properties
                    if (Me%ComputeOptions%Fertilization) then
                        if (Me%ComputeOptions%ModelNitrogen) then
                            nPropertiesToSoil = nPropertiesToSoil + 8
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            nPropertiesToSoil = nPropertiesToSoil + 6
                        endif
                    endif
                    if (Me%ComputeOptions%HarvestKill) then
                        nPropertiesToSoil = nPropertiesToSoil + 2
                        if (Me%ComputeOptions%ModelNitrogen) then
                            nPropertiesToSoil = nPropertiesToSoil + 1
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            nPropertiesToSoil = nPropertiesToSoil + 1
                        endif
                    endif                
                    if (Me%ComputeOptions%Grazing) then
                        nPropertiesToSoil = nPropertiesToSoil + 1
                        if (Me%ComputeOptions%ModelNitrogen) then
                            nPropertiesToSoil = nPropertiesToSoil + 1
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            nPropertiesToSoil = nPropertiesToSoil + 1
                        endif
                    endif   
                    if (Me%ComputeOptions%Dormancy) then
                        nPropertiesToSoil = nPropertiesToSoil + 1
                        if (Me%ComputeOptions%ModelNitrogen) then
                            nPropertiesToSoil = nPropertiesToSoil + 1
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            nPropertiesToSoil = nPropertiesToSoil + 1
                        endif
                    endif 
                    if (Me%ComputeOptions%Pesticide) then
                        do Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
                            nPropertiesToSoil = nPropertiesToSoil + 2
                        enddo
                    endif 

                    allocate(PropertyList(1:nPropertiesToSoil), STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructTimeSerie - ModuleVegetation - ERR08'

                    i = 0
                    !Properties header
                    if (Me%ComputeOptions%Fertilization) then
                        if (Me%ComputeOptions%ModelNitrogen) then
                            i = i + 1
                            PropertyList(i) = "Fert_Nitrate_Surface_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_Nitrate_SubSurface_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_Ammonia_Surface_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_Ammonia_SubSurface_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_OrganicN_Surface_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_OrganicN_SubSurface_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_OrganicN_FluffLayer_KgN/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_YearlyAccumulated_Nitrogen_KgN/ha"
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            i = i + 1
                            PropertyList(i) = "Fert_OrganicP_Surface_KgP/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_OrganicP_SubSurface_KgP/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_MineralP_Surface_KgP/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_MineralP_SubSurface_KgP/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_OrganicP_FluffLayer_KgP/ha"
                            i = i + 1
                            PropertyList(i) = "Fert_YearlyAccumulated_Phosphorus_KgP/ha"                 
                        endif
                    endif
                    if (Me%ComputeOptions%HarvestKill) then
                        i = i + 1
                        PropertyList(i) = "HarvestKill_AerialBiomass_To_Soil_Kg/ha"
                        i = i + 1
                        PropertyList(i) = "HarvestKill_RootBiomass_Left_In_Soil_Kg/ha"                    
                        if (Me%ComputeOptions%ModelNitrogen) then
                            i = i + 1
                            PropertyList(i) = "HarvestKill_Nitrogen_To_Soil_KgN/ha"
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            i = i + 1
                            PropertyList(i) = "HarvestKill_Phosphorus_To_Soil_KgP/ha"
                        endif
                    endif                
                    if (Me%ComputeOptions%Grazing) then
                        i = i + 1
                        PropertyList(i) = "Grazing_Biomass_To_Soil_Kg/ha"
                        if (Me%ComputeOptions%ModelNitrogen) then
                            i = i + 1
                            PropertyList(i) = "Grazing_Nitrogen_To_Soil_KgN/ha"
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            i = i + 1
                            PropertyList(i) = "Grazing_Phosphorus_To_Soil_KgP/ha"
                        endif
                    endif   
                    if (Me%ComputeOptions%Dormancy) then
                        i = i + 1
                        PropertyList(i) = "Dormancy_Biomass_To_Soil_Kg/ha"
                        if (Me%ComputeOptions%ModelNitrogen) then
                            i = i + 1
                            PropertyList(i) = "Dormancy_Nitrogen_To_Soil_KgN/ha"
                        endif
                        if (Me%ComputeOptions%ModelPhosphorus) then
                            i = i + 1
                            PropertyList(i) = "Dormancy_Phosphorus_To_Soil_KgP/ha"
                        endif
                    endif  
                    if (Me%ComputeOptions%Pesticide) then
                        do Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
                            i = i + 1
                            PropertyList(i) = trim(Me%Fluxes%Pesticides%Application(Pest)%ID%Name)//"_Soil_Kg/ha"
                            i = i + 1
                            PropertyList(i) = trim(Me%Fluxes%Pesticides%Application(Pest)%ID%Name)//"_Veg_Kg/ha"  
                        enddo
                    endif 
                    call StartTimeSerie(Me%ObjTimeSerieToSoil,                                      &
                                        Me%ObjTime,                                                 &
                                        TimeSerieLocationFile,                                      &
                                        PropertyList, "srvgf",                                      &
                                        WaterPoints2D = Me%ExternalVar%MappingPoints,               &
                                        STAT           = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleVegetation - ERR06'

                    deallocate(PropertyList)
                
                endif
            
            endif

            ! Time serie with fertilization fluxes to soil
!            if (Me%ComputeOptions%Fertilization) then
!
!                call GetData(Me%ComputeOptions%FertilizationOutput,                 &
!                             Me%ObjEnterData,iflag,                                 &
!                             SearchType   = FromFile,                               &
!                             keyword      = 'FERTILIZATION_OUTPUT',                 &
!                             ClientModule = 'ModuleVegetation',                     &
!                             Default      = .false.,                                &
!                             STAT         = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_)                                        &
!                    stop 'Construct_Time_Serie - ModuleVegetation - ERR07' 
!                
!                if (Me%ComputeOptions%FertilizationOutput) then
!                                
!                    if (Me%ComputeOptions%ModelNitrogen) then
!                        if (Me%ComputeOptions%ModelPhosphorus) then
!                            allocate (PropertyList(12))
!                        else
!                            allocate (PropertyList(7))
!                        endif
!                    else
!                        if (Me%ComputeOptions%ModelPhosphorus) then
!                            allocate (PropertyList(5))
!                        endif
!                    endif
!                    
!                    i = 1                        
!                    if (Me%ComputeOptions%ModelNitrogen) then
!                        PropertyList(i)   = "Fert_Nitrate_Surface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_Nitrate_SubSurface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_Ammonia_Surface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_Ammonia_SubSurface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_OrganicN_Surface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_OrganicN_SubSurface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_YearlyAccumulated_Nitrogen_KgN/ha"
!                        i = i + 1
!                    endif
!
!                    if (Me%ComputeOptions%ModelPhosphorus) then
!                        PropertyList(i) = "Fert_OrganicP_Surface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_OrganicP_SubSurface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_MineralP_Surface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_MineralP_SubSurface_KgN/ha"
!                        i = i + 1
!                        PropertyList(i) = "Fert_YearlyAccumulated_Phosphorus_KgN/ha"
!                        i = i + 1
!                    endif
!        
!                    call StartTimeSerie(Me%ObjTimeSerieFert,                                        &
!                                        Me%ObjTime,                                                 &
!                                        TimeSerieLocationFile,                                      &
!                                        PropertyList, "vgf",                                        &
!                                        WaterPoints2D = Me%ExternalVar%MappingPoints2D,             &
!                                        Instance = 'Vegetation_Fertilization',                      &
!                                        STAT           = STAT_CALL) 
!                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleVegetation - ERR06'
!
!                    deallocate(PropertyList)
!
!                endif
!
!            endif
            
            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR03'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR04'
                
                call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR04'
                
    i1:         if (CoordON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR05'

                    if (Id < 0 .or. Jd < 0) then
                    
                        call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR06'

                        if (IgnoreOK) then
                            write(*,*) 'Time Serie outside the domain - ',trim(TimeSerieName)
                            cycle
                        else
                            stop 'ConstructTimeSerie - PorousMedia - ERR07'
                        endif

                    endif


                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR08'
                    
                    if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then   
                        if(Me%ComputeOptions%AtmospherePropertiesOutput) then
                            call CorrectsCellsTimeSerie(Me%ObjTimeSerieAtm, dn, Id, Jd, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR10'                        
                        endif
                        
                        if(Me%ComputeOptions%FluxesToSoilOutput) then             
                            call CorrectsCellsTimeSerie(Me%ObjTimeSerieToSoil, dn, Id, Jd, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR20'                    
                        endif
                        
                    endif

                endif i1

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          LocalizationI   = Id,                             &
                                          LocalizationJ   = Jd,                             & 
                                          STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR30'

                if (Me%ExternalVar%MappingPoints(Id, Jd) /= VegetationPoint) then
                     write(*,*) 'Time Serie in a outside boundaries cell - ',trim(TimeSerieName)
                endif

            enddo  
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
            Me%OutPut%NextOutPut = 1
            if (STAT_CALL /= SUCCESS_)                                       &
                stop 'ConstructHDF - ModuleVegetation - ERR01' 

            if (Me%OutPut%HDF_ON) then

                call ConstructHDF5Output

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructHDF - ModuleVegetation - ERR02'
            endif 

        endif

    end subroutine ConstructHDF

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

       
        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR01'
        
        !Gets Basin Points
        call GetBasinPoints (Me%ObjBasinGeometry, Me%ExternalVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR01.1'        
        
        
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%Results)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR02'

      
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR02'


        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR05'
        
        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR05'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR07'

        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR08'       


        !Unget Basin Points
        call UnGetBasin (Me%ObjBasinGeometry, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR08.1'

        !UnGets Topography
        call UnGetGridData      (Me%ObjGridData, Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleVegetation - ERR09'


    end subroutine ConstructHDF5Output

   !--------------------------------------------------------------------------

    subroutine ConstructVegetationGrids

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL !, ObjGD
        integer                                 :: i, j, ivt, ClientNumber
        logical                                 :: BlockFound
        type (T_PropertyID)                     :: PotentialHUID !, AgricPract
        real                                    :: NitrogenFractionInYeld, BiomassEnergyRatio
        real                                    :: OptimalHarvestIndex
        real, dimension(:,:), pointer           :: AgricPractIDMatrix
        !Begin------------------------------------------------------------------
        
        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
        
            !! Construct Potential Annual HU
            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR10'

            !Constructs Potential Total HU in one year
            call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                        ClientNumber    = ClientNumber,                         &
                                        block_begin     = '<begin_TotalPotentialHU>',           &
                                        block_end       = '<end_TotalPotentialHU>',             &
                                        BlockFound      = BlockFound,                           &   
                                        STAT            = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                if (.not. BlockFound) then
                    write(*,*)'Missing Block <begin_TotalPotentialHU> / <end_TotalPotentialHU>'
                    stop 'ConstructVegetationGrids - ModuleVegetation - ERR20'
                endif

                !Gets a pointer to OpenPoints2D
                call GetBasinPoints  (Me%ObjHorizontalMap, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR025'
                
                call ConstructFillMatrix  ( PropertyID           = PotentialHUID,                       &
                                            EnterDataID          = Me%ObjEnterData,                     &
                                            TimeID               = Me%ObjTime,                          &
                                            HorizontalGridID     = Me%ObjHorizontalGrid,                &
                                            ExtractType          = FromBlock,                           &
                                            PointsToFill2D       = Me%ExternalVar%BasinPoints,          &
                                            Matrix2D             = Me%HeatUnits%PotentialHUTotal,       &
                                            TypeZUV              = TypeZ_,                              &
                                            STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR30'
                
                call KillFillMatrix       (PotentialHUID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR40'

                call UngetBasin (Me%ObjHorizontalMap, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR50'
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR70'      
            
            else
                stop 'ConstructVegetationGrids - ModuleVegetation - ERR60'
            endif                
            
      

          
                
        endif


!        ObjGD = 0
!        call ConstructGridData  (ObjGD, Me%ObjHorizontalGrid, FileName = Me%Files%VegetationIDFile, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR80'
!        
!        !allocate (AuxID(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
!        call GetGridData (ObjGD, AgricPractIDMatrix, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR90'


        !Get Agricultural Practices ID
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR70'
        
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
                                    ClientNumber    = ClientNumber,                         &
                                    block_begin     = '<begin_AgriculturalPractices>',      &
                                    block_end       = '<end_AgriculturalPractices>',        &
                                    BlockFound      = BlockFound,                           &   
                                    STAT            = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. BlockFound) then
                write(*,*)'Missing Block <begin_AgriculturalPractices> / <end_AgriculturalPractices>'
                write(*,*)'This block removes old VEGETATION_ID_FILE and now input is done using'
                write(*,*)'FillMatrix standards: constant, timeserie, ASCII grid or HDF input'
                stop 'ConstructVegetationGrids - ModuleVegetation - ERR80'
            endif        
            
            allocate (AgricPractIDMatrix(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            AgricPractIDMatrix = FillValueReal
            
            !Gets Basin Points
            call GetBasinPoints (Me%ObjBasinGeometry, Me%ExternalVar%BasinPoints, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR090'  
            
            !fill the first grid  (only real values allowed in fillmatrix)
            call ConstructFillMatrix  (PropertyID           = Me%AgricPract,                    &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExternalVar%BasinPoints,       &
                                       Matrix2D             = AgricPractIDMatrix,               &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'ConstructVegetationGrids - ModuleVegetation - ERR100'


            if(.not. Me%AgricPract%SolutionFromFile)then

                ! only one grid
                call KillFillMatrix(Me%AgricPract%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructVegetationGrids - ModuleVegetation - ERR110'
            
            endif
            
            !Save First ID grid in integer form (FillMatrix only allows reals)
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExternalVar%MappingPoints(i, j) == VegetationPoint) then
                    
                    Me%AgricPractID(i,j) = NINT(AgricPractIDMatrix(i,j))
                    
                endif
            enddo
            enddo            
            
            !Unget Basin Points
            call UnGetBasin (Me%ObjBasinGeometry, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR0120'        


            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationGrids - ModuleVegetation - ERR130'      
        
        else
            stop 'ConstructVegetationGrids - ModuleVegetation - ERR140'
        endif  
        
        !Gets first Vegetation IDs for next actions
        

        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                do ivt = 1, Me%VegetationsNumber
                    !if (Me%VegetationTypes(ivt)%ID == NINT(AgricPractIDMatrix(i, j))) then
                    if (Me%VegetationTypes(ivt)%ID == Me%AgricPractID(i, j)) then
                        !Me%VegetationID(i, j) = NINT(AuxID(i, j))
                        Me%VegetationID(i, j) = ivt
                    endif
!                    if (ivt .eq. Me%VegetationsNumber .and. Me%VegetationID(i,j) < FillValueInt / 2) then
!                        write(*,*)'Vegetation ID not found in vegetation type definition.'
!                  write(*,*)'Check in vegetation grid or in vegetation types in data file, the ID:', NINT(Me%AgricPractID(i,j))
!                        stop 'ConstructVegetationGrid - ModuleVegetation - ERR100'
!                    endif
                enddo

            endif
        enddo
        enddo
        

        !deallocate (AuxID)
        
!        call UnGetGridData (ObjGD, AgricPractIDMatrix, STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR110'
!        
!        call KillGridData (ObjGD, STAT = STAT_CALL)   
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationGrid - ModuleVegetation - ERR120'

        
        if (Me%ComputeOptions%AutoFertilization) then
            !Fill first nitrogen yeld estimate
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                    
                    !Plant type 0 is NO PLANT, so no operations should happen
                    if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                    !kgN/kgYeld 
                    NitrogenFractionInYeld = Me%VegetationTypes(Me%VegetationID(i, j))%GrowthDatabase%NitrogenFractionInYeld
                    !kg/ha/(MJ/m2) = 10-1 g/MJ
                    BiomassEnergyRatio     = Me%VegetationTypes(Me%VegetationID(i, j))%GrowthDatabase%BiomassEnergyRatio
                    OptimalHarvestIndex    = Me%VegetationTypes(Me%VegetationID(i, j))%GrowthDatabase%OptimalHarvestIndex
                    if (OptimalHarvestIndex .lt. 1.) then
                        !KgN/ha
                        Me%NitrogenYeldEstimate(i,j) = 350. * NitrogenFractionInYeld * BiomassEnergyRatio
                    else
                        Me%NitrogenYeldEstimate(i,j) = 1000. * NitrogenFractionInYeld * BiomassEnergyRatio
                    endif
                    
                endif
            enddo
            enddo
        endif

    end subroutine ConstructVegetationGrids

   !--------------------------------------------------------------------------

!    subroutine CheckRootDepth
!
!        !External--------------------------------------------------------------
!        integer                                 :: STAT_CALL
!        integer                                 :: i, j, kTop, kFloor
!        real                                    :: MaxRootDepth, SoilDepth
!        
!        !Begin------------------------------------------------------------------
!
!        call GetGeometryDistances(Me%ObjGeometry,                                       &
!                                  SZZ         = Me%ExternalVar%SZZ,                     &
!                                  STAT        = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR01'
!
!        call GetGeometryKFloor(Me%ObjGeometry,                                          &
!                               Z    = Me%ExternalVar%KFloor,                            &
!                               STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR02'
!
!        do j = Me%Size%JLB, Me%Size%JUB
!        do i = Me%Size%ILB, Me%Size%IUB
!            
!            if (Me%ExternalVar%MappingPoints2D(i, j) == 1) then
!                
!                !Check maximum root depth
!                MaxRootDepth = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaximumRootDepth
!                kTop    = Me%WorkSize%KUB
!                kFloor  = Me%ExternalVar%KFloor(i,j)
!                
!                SoilDepth = (-1. * Me%ExternalVar%SZZ(i,j,kTop)) - (-1.* Me%ExternalVar%SZZ(i,j,kFloor))
!
!                if(MaxRootDepth .gt. SoilDepth) then
!                    
!                    Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaximumRootDepth = SoilDepth
!                    write(*,*)'Maximum Root Depth is greater then soil depth. Maximum root depth was set to soil depth'
!                    write(*,*)'Cell:', i, j
!                
!                endif
!            
!            endif
!        
!        enddo
!        enddo
!
!        call UnGetGeometry( Me%ObjGeometry, Me%ExternalVar%KFloor,       STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR03'
!        
!        call UnGetGeometry( Me%ObjGeometry, Me%ExternalVar%SZZ,          STAT = STAT_CALL )        
!        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR04'
!
!    end subroutine CheckRootDepth

   !--------------------------------------------------------------------------

    subroutine CheckLatitude

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: i, j
        
        !Begin------------------------------------------------------------------

        call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                     &
                                      GridLatitude  = Me%ExternalVar%Latitude,  &
                                      STAT          = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'CheckLatitude - ModuleVegetation - ERR01'


        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                if (Me%ExternalVar%Latitude(i,j) .gt. 90.) then
                    write(*,*) ' Computing vegetation dormancy and Latitude is greater than 90.'
                    write(*,*) ' Please check DTM coordinates and if using grid coordinates (type 5)'
                    write(*,*) ' please check LATITUDE keyword in DTM file'
                    stop 'CheckLatitude - ModuleVegetation ERR20'
                endif
            
            endif
        
        enddo
        enddo

        !Latitude
        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%Latitude, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckLatitude - ModuleVegetation - ERR30'

    end subroutine CheckLatitude

   !--------------------------------------------------------------------------

    subroutine ConstructFeddes

        !Local----------------------------------------------------------------
        integer                                      :: i,j
        integer                                      :: veg_id
        integer                                      :: veg_type
        !Begin----------------------------------------------------------------   
    

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then 
                veg_id = Me%VegetationID(i,j)
                veg_type = Me%VegetationTypes(veg_id)%GrowthDatabase%PlantType
                                      
                !Plant type 0 is NO PLANT, so no operations should happen
                if (veg_type == NotAPlant) cycle
            
                Me%ComputeOptions%TranspirationMOHID%FeddesType(i,j) = &
                    Me%VegetationTypes(veg_id)%FeddesType 

                Me%ComputeOptions%TranspirationMOHID%RootFeddesH1(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesH1 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH2(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesH2 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH3(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesH3 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH4(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesH4 
                
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH3L(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesH3L 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesH3H(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesH3H 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesR2L(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesR2L 
                Me%ComputeOptions%TranspirationMOHID%RootFeddesR2H(i,j) = &
                    Me%VegetationTypes(veg_id)%RootFeddesR2H 
            endif
        enddo
        enddo
    
   
    end subroutine ConstructFeddes

   !--------------------------------------------------------------------------
    
    subroutine ConstructVegetationParameters

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ClientNumber, iflag
        integer                                     :: ivt, VegetationTypeID
        character (Len = StringLength)              :: ParameterFile
        integer                                     :: ParameterObjEnterData
        logical                                     :: VegetationFound, DatabaseFound
        type (T_VegetationType), pointer            :: VegetationType
        integer                                     :: Pest
        !----------------------------------------------------------------------

        !After constructing vegetation list module knows how many vegetations exist
        allocate(Me%VegetationTypes(Me%VegetationsNumber))
        
        !Get Vegetation parameters file 
        call GetData(ParameterFile, Me%ObjEnterData, iflag,                                                       &
                     SearchType     = FromFile,                                                                   &
                     keyword        = 'PARAMETERS_FILE',                                                          &
                     ClientModule   = 'ModuleVegetation',                                                         &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR10'     
        if (iflag /= 1) then
            write(*,*)'Missing PARAMETERS_FILE in vegetation data file'
            stop 'ConstructVegetationParameters - ModuleVegetation - ERR75'            
        endif
        
        ivt = 0
        ParameterObjEnterData = 0
        !Open and save parameter file
        call ConstructEnterData(ParameterObjEnterData, ParameterFile, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR20'        
        
        if (Me%ComputeOptions%Pesticide) then
            !temporary lists used to allocate unique pesticides
            allocate (Me%PesticideListID(100))  
            allocate (Me%PesticideListName(100))   
            Me%Fluxes%Pesticides%UniquePesticides  = 1              
        endif
        
        !Check for every vegetation type (agricultural practice) for info
        VegetationType => Me%FirstVegetation
doV:    do while (associated(VegetationType))
            
            !each property cycle restart from the beggining og the file
            call RewindBuffer(ParameterObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR30'

            
            VegetationFound = .false. 
            
doH1:       do 
                call ExtractBlockFromBuffer(ParameterObjEnterData,                                    &
                                            ClientNumber      = ClientNumber,                         &
                                            block_begin       = '<beginagriculturalpractice>',        &
                                            block_end         = '<endagriculturalpractice>',         &
                                            BlockFound        = DatabaseFound,                        &   
                                            STAT              = STAT_CALL)
HF1:            if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

                    call GetData(VegetationTypeID, ParameterObjEnterData, iflag,                                             &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'AGRIC_PRACT_ID',                                                          &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR40'     
                    if (iflag /= 1) then
                        write(*,*)'Missing AGRIC_PRACT_ID in agricultural practice definition'
                        write(*,*)'Check vegetation definition in',ParameterFile
                        stop 'ConstructVegetationParameters - ModuleVegetation - ERR50'
                    endif
                    
                    if (VegetationTypeID == VegetationType%ID) then
                        
                        VegetationFound = .true.
                        ivt = ivt + 1
                        Me%VegetationTypes(ivt)%ID = VegetationType%ID
                        
                        !Reads Name
                        call GetData(Me%VegetationTypes(ivt)%Name, ParameterObjEnterData,  iflag,             &
                                     SearchType     = FromBlock,                                        &
                                     keyword        = 'NAME',                                           &
                                     ClientModule   = 'ModuleVegetation',                               &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR60'
                        if (iflag /= 1) then
                            write(*,*)'Missing NAME in Vegetation Type definition'
                            stop 'ConstructVegetationParameters - ModuleVegetation - ERR70'
                        endif

                        !Reads Crop ID associated
                        call GetData(Me%VegetationTypes(ivt)%VegetationID, ParameterObjEnterData,  iflag,     &
                                     SearchType     = FromBlock,                                        &
                                     keyword        = 'VEGETATION_ID',                                  &
                                     ClientModule   = 'ModuleVegetation',                               &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR80'
                        if (iflag /= 1) then
                            write(*,*)'Missing VEGETATION_ID in Vegetation Type definition'
                            stop 'ConstructVegetationParameters - ModuleVegetation - ERR90'
                        endif
                        
                        
                        if ((Me%ComputeOptions%ModelWater) .and.                                &
                            (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) .and. &
                            (Me%ComputeOptions%WaterStress)) then

                            !Get Feddes from database
                            call GetData(Me%FeddesDatabase, Me%ObjEnterData, iflag,                                          &
                                         SearchType     = FromFile,                                                          &
                                         keyword        = 'FEDDES_DATABASE',                                                 &
                                         ClientModule   = 'ModuleVegetation',                                                &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR100'  
                            if (iflag /= 1) then
                                write(*,*)'Missing FEDDES_DATABASE in vegetation data file'
                                stop 'ConstructVegetationParameters - ModuleVegetation - ERR110'            
                            endif                        
                            call ReadFeddesDatabase (ivt)
                            
                        endif

                        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

                            call GetData(Me%VegetationTypes(ivt)%HasLeaves, ParameterObjEnterData,  iflag,        &
                                         SearchType     = FromBlock,                                        &
                                         keyword        = 'HAS_LEAVES',                                     &
                                         default        = .true.,                                           &
                                         ClientModule   = 'ModuleVegetation',                               &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR120'
!                            if (iflag /= 1) then
!                                write(*,*)'Missing HAS_LEAVES in Vegetation Type definition'
!                                stop 'ConstructVegetationParameters - ModuleVegetation - ERR32'
!                            endif
                            
                            call ReadTimingParameters(ivt, ClientNumber,ParameterObjEnterData)

                            if (Me%ComputeOptions%HarvestKill) then 
                                call ReadHarvestKillParameters(ivt, ClientNumber,ParameterObjEnterData)
                            endif

                            if (Me%ComputeOptions%Grazing) then 
                                call ReadGrazingParameters(ivt, ClientNumber,ParameterObjEnterData)
                            endif
                            
                            if (Me%ComputeOptions%Fertilization) then 
                            
                                !Get fertilizers from database
                                call GetData(Me%FertilizerDatabase, Me%ObjEnterData, iflag,                                      &
                                             SearchType     = FromFile,                                                          &
                                             keyword        = 'FERTILIZER_DATABASE',                                             &
                                             ClientModule   = 'ModuleVegetation',                                                &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR130'   
                                if (iflag /= 1) then
                                    write(*,*)'Missing FERTILIZER_DATABASE in vegetation data file'
                                    stop 'ConstructVegetationParameters - ModuleVegetation - ERR140'            
                                endif                            
                                call ReadFertilizationParameters(ivt, ClientNumber,ParameterObjEnterData)
                            endif
                            if (Me%ComputeOptions%Pesticide) then 

                                !Get Pesticide from database
                                call GetData(Me%PesticideDatabase, Me%ObjEnterData, iflag,                                       &
                                             SearchType     = FromFile,                                                          &
                                             keyword        = 'PESTICIDE_DATABASE',                                              &
                                             ClientModule   = 'ModuleVegetation',                                                &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR150'     
                                if (iflag /= 1) then
                                    write(*,*)'Missing PESTICIDE_DATABASE in vegetation data file'
                                    stop 'ConstructVegetationParameters - ModuleVegetation - ERR160'            
                                endif                                
                                call ReadPesticideParameters(ivt, ClientNumber,ParameterObjEnterData)
                            endif                    

                            !Get growth database (SWAT growth parameters)
                            call GetData(Me%GrowthDatabase, Me%ObjEnterData, iflag,                                              &
                                         SearchType     = FromFile,                                                              &
                                         keyword        = 'GROWTH_DATABASE',                                                     &
                                         ClientModule   = 'ModuleVegetation',                                                    &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR170'   
                            if (iflag /= 1) then
                                write(*,*)'Missing GROWTH_DATABASE in vegetation data file'
                                stop 'ConstructVegetationParameters - ModuleVegetation - ERR180'            
                            endif                            
                            !Read growth parameters (from SWAT database)
                            call ReadGrowthDatabase(ivt)
                        
                        endif
                        
                        exit doH1
                        
                    endif
                    
                else
                    
                    call Block_Unlock(ParameterObjEnterData, ClientNumber, STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructVegetationParameters - ModuleVegetation - ERR190'
                   
                    exit doH1 
                              
                endif HF1
                
            enddo doH1
             
            if (.not. VegetationFound) then
                write(*,*)
                write(*,*) 'Vegetation ID read in Vegetation Grid'
                write(*,*) 'not found in parameter file.'
                write(*,*) 'Check that ID' , VegetationType%ID
                write(*,*) 'was defined in parameter file', trim(ParameterFile)
                stop 'ConstructVegetationParameters - ModuleVegetation - ERR200' 
            endif 
            
            VegetationType => VegetationType%Next
        
        enddo doV
        
        if (Me%ComputeOptions%Pesticide) then
            !Fill unique global pesticide list
            !allocate one position less because first position was phantom
            allocate(Me%Fluxes%Pesticides%Application(Me%Fluxes%Pesticides%UniquePesticides - 1))
            do Pest = 2, Me%Fluxes%Pesticides%UniquePesticides
                
                Me%Fluxes%Pesticides%Application(Pest - 1)%PesticideID  = Me%PesticideListID(Pest)
                Me%Fluxes%Pesticides%Application(Pest - 1)%ID%Name      = Me%PesticideListName(Pest)
                if (.not. CheckPropertyName (Me%Fluxes%Pesticides%Application(Pest - 1)%ID%Name,     &
                                             Me%Fluxes%Pesticides%Application(Pest - 1)%ID%IDNumber)) then
                    write (*,*)'The property isnt recognized by the model :'
                    write (*,*) trim(Me%Fluxes%Pesticides%Application(Pest - 1)%ID%Name)
                    stop 'ConstructVegetationParameters - ConstructVegetationParameters - ERR210'
                endif
                
            enddo
            Me%Fluxes%Pesticides%UniquePesticides = size(Me%Fluxes%Pesticides%Application)
        endif
       
        call KillEnterData(ParameterObjEnterData, STAT = STAT_CALL)         
        if (STAT_CALL /= SUCCESS_) stop 'ConstructVegetationParameters - ConstructVegetationParameters - ERR220'        
   
    end subroutine ConstructVegetationParameters

   !--------------------------------------------------------------------------

    subroutine ReadFeddesDatabase(ivt)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound, VegetationFound
        integer                                     :: iflag
        integer                                     :: FeddesObjEnterData
        integer                                     :: VegetationID
        integer                                     :: FeddesClientNumber
        !Begin-----------------------------------------------------------------

        FeddesObjEnterData = 0
        !Open and save growth database
        call ConstructEnterData(FeddesObjEnterData, Me%FeddesDatabase, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR010'        
        
           
        VegetationFound = .false. 
            
doH1:   do 
            call ExtractBlockFromBuffer(FeddesObjEnterData,                                      &
                                        ClientNumber      = FeddesClientNumber,                  &
                                        block_begin       = '<beginfeddesdatabase>',             &
                                        block_end         = '<endfeddesdatabase>',               &
                                        BlockFound        = DatabaseFound,                       &   
                                        STAT              = STAT_CALL)
HF1:        if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then
                
                call GetData(VegetationID, FeddesObjEnterData, iflag,                            &
                             SearchType     = FromBlock,                                         &
                             keyword        = 'VEGETATION_ID',                                   &
                             ClientModule   = 'ModuleVegetation',                                &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddedDatabase - ModuleVegetation - ERR020'     
                if (iflag /= 1) then
                    write(*,*)'Missing VEGETATION_ID in Feddes definition'
                    write(*,*)'Check vegetation definition in',Me%FeddesDatabase
                    stop 'ReadFeddesDatabase - ModuleVegetation - ERR030'
                endif
                                    
                if (VegetationID == Me%VegetationTypes(ivt)%VegetationID) then
                    
                    VegetationFound = .true.
                    
                    !Reads Feddes
                    call GetData(Me%VegetationTypes(ivt)%FeddesType, FeddesObjEnterData, iflag,     &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_TYPE',                                    &
                                 default        = 1,                                                &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR040'
                    
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH1, FeddesObjEnterData, iflag,   &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H1',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR050'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H1 in Vegetation Type definition'
                        stop 'ReadFeddesDatabase - ModuleVegetation - ERR50'
                    endif
                    
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH2, FeddesObjEnterData, iflag,   &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H2',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR60'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H2 in Vegetation Type definition'
                        stop 'ReadFeddesDatabase - ModuleVegetation - ERR70'
                    endif
                    
                    select case (Me%VegetationTypes(ivt)%FeddesType) 
                    case (1)
                        call GetData(Me%VegetationTypes(ivt)%RootFeddesH3, FeddesObjEnterData, iflag,  &
                                     SearchType     = FromBlock,                                       &
                                     keyword        = 'FEDDES_H3',                                     &
                                     ClientModule   = 'ModuleVegetation',                              &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR070'
                        if (iflag /= 1) then
                            write(*,*)'Missing FEDDES_H3 in Vegetation Type definition'
                            stop 'ReadFeddesDatabase - ModuleVegetation - ERR080'
                        endif
                    case (2)
                        call GetData(Me%VegetationTypes(ivt)%RootFeddesH3L, FeddesObjEnterData, iflag, &
                                     SearchType     = FromBlock,                                       &
                                     keyword        = 'FEDDES_H3L',                                    &
                                     ClientModule   = 'ModuleVegetation',                              &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR090'
                        if (iflag /= 1) then
                            write(*,*)'Missing FEDDES_H3L in Vegetation Type definition'
                            stop 'ReadFeddesDatabase - ModuleVegetation - ERR100'
                        endif
                        call GetData(Me%VegetationTypes(ivt)%RootFeddesH3H, FeddesObjEnterData, iflag, &
                                     SearchType     = FromBlock,                                       &
                                     keyword        = 'FEDDES_H3H',                                    &
                                     ClientModule   = 'ModuleVegetation',                              &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR110'
                        if (iflag /= 1) then
                            write(*,*)'Missing FEDDES_R2L in Vegetation Type definition'
                            stop 'ReadFeddesDatabase - ModuleVegetation - ERR120'
                        endif
                    
                        call GetData(Me%VegetationTypes(ivt)%RootFeddesR2L, FeddesObjEnterData, iflag, &
                                     SearchType     = FromBlock,                                       &
                                     keyword        = 'FEDDES_R2L',                                    &
                                     ClientModule   = 'ModuleVegetation',                              &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR130'
                        if (iflag /= 1) then
                            write(*,*)'Missing FEDDES_H3H in Vegetation Type definition'
                            stop 'ReadFeddesDatabase - ModuleVegetation - ERR140'
                        endif

                        call GetData(Me%VegetationTypes(ivt)%RootFeddesR2H, FeddesObjEnterData, iflag, &
                                     SearchType     = FromBlock,                                       &
                                     keyword        = 'FEDDES_R2H',                                    &
                                     ClientModule   = 'ModuleVegetation',                              &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR150'
                        if (iflag /= 1) then
                            write(*,*)'Missing FEDDES_R2H in Vegetation Type definition'
                            stop 'ReadFeddesDatabase - ModuleVegetation - ERR160'
                        endif

                    case default
                        write(*,*)'Invalid FEDDES_TYPE in Vegetation Type definition'
                        stop 'ReadFeddesDatabase - ModuleVegetation - ERR170'
                    end select                    
                    
                    call GetData(Me%VegetationTypes(ivt)%RootFeddesH4, FeddesObjEnterData,  iflag,  &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'FEDDES_H4',                                      &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR180'
                    if (iflag /= 1) then
                        write(*,*)'Missing FEDDES_H4 in Vegetation Type definition'
                        stop 'ReadFeddesDatabase - ModuleVegetation - ERR190'
                    endif
                    
                    !Load SALINITY STRESS parameters
                    call GetData(Me%VegetationTypes(ivt)%UseSalinityStress,                         &
                                 FeddesObjEnterData,  iflag,                                        &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'USE_SALINITY_STRESS',                            &
                                 default        = .false.,                                          &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR200'
                    
                    if (Me%VegetationTypes(ivt)%UseSalinityStress) then
                    
                        call GetData(Me%VegetationTypes(ivt)%SalinityStressParams%ECt,                  &
                                     FeddesObjEnterData,  iflag,                                        &
                                     SearchType     = FromBlock,                                        &
                                     keyword        = 'SALINITY_STRESS_THRESHOLD',                      &
                                     ClientModule   = 'ModuleVegetation',                               &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR210'
                        if (iflag .NE. 1) stop 'ReadFeddesDatabase - ModuleVegetation - ERR220'     
                        
                        call GetData(Me%VegetationTypes(ivt)%SalinityStressParams%Slope,                &
                                     FeddesObjEnterData,  iflag,                                        &
                                     SearchType     = FromBlock,                                        &
                                     keyword        = 'SALINITY_STRESS_SLOPE',                          &
                                     ClientModule   = 'ModuleVegetation',                               &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR230'
                        if (iflag .NE. 1) stop 'ReadFeddesDatabase - ModuleVegetation - ERR240' 
                                                               
                        Me%VegetationTypes(ivt)%SalinityStressParams%EC0 = Me%VegetationTypes(ivt)%SalinityStressParams%ECt + &
                                                                      (100 / Me%VegetationTypes(ivt)%SalinityStressParams%Slope)
                    endif

                    call GetData(Me%VegetationTypes(ivt)%StressInteraction,                         &
                                 FeddesObjEnterData,  iflag,                                        &
                                 SearchType     = FromBlock,                                        &
                                 keyword        = 'STRESS_INTERACTION',                             &
                                 default        = 1,                                                &
                                 ClientModule   = 'ModuleVegetation',                               &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddesDatabase - ModuleVegetation - ERR250'
                                        
                    exit doH1
                    
                endif
                
            else
                
                call Block_Unlock(FeddesObjEnterData, FeddesClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFeddedDatabase - ModuleVegetation - ERR260'
               
                exit doH1 
                          
            endif HF1
            
        enddo doH1
             
        if (.not. VegetationFound) then
            write(*,*)
            write(*,*) 'Vegetation ID not found:', Me%VegetationTypes(ivt)%VegetationID
            write(*,*) 'Check that ID is defined in database:',Me%FeddesDatabase
            stop 'ReadFeddesDatabase - ModuleVegetation - ERR270' 
        endif           
    
    
        call KillEnterData(FeddesObjEnterData, STAT = STAT_CALL)         
        if (STAT_CALL /= SUCCESS_) stop 'ReadFeddeshDatabase - ModuleVegetation - ERR280'        

    end subroutine ReadFeddesDatabase

    !--------------------------------------------------------------------------
    
    subroutine ReadTimingParameters(ivt, ClientNumber,ParameterObjEnterData)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        integer                                     :: ParameterObjEnterData
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound
        integer                                     :: iflag, FailRead

        !Begin-----------------------------------------------------------------



        call ExtractBlockFromBlock(ParameterObjEnterData,                                     &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<begintimingparameters>',            &
                                    block_end         = '<endtimingparameters>',              &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

            !Reads timing Parameters

            call GetData(Me%VegetationTypes(ivt)%TimingDatabase%PlantHUAtMaturity, ParameterObjEnterData,  iflag,         &
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
                    
            call GetData(Me%VegetationTypes(ivt)%TimingDatabase%PlantingJulianDay, ParameterObjEnterData,  iflag,         &
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
                        
            call GetData(Me%VegetationTypes(ivt)%TimingDatabase%PlantingHUBase, ParameterObjEnterData,  iflag,            &
                         SearchType     = FromBlockInBlock,                                                               &
                         Default        = -99.,                                                                           &
                         keyword        = 'PLANTING_HUBASE',                                                              &
                         ClientModule   = 'ModuleVegetation',                                                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadTimingDatabase - ModuleVegetation - ERR30'        
            if (iflag /= 1) then
                FailRead = FailRead + 1
                if (FailRead == 2) then
                    write(*,*)'Missing HU or julian day in planting timing'
                    write(*,*)' parameters definition -', trim(Me%VegetationTypes(ivt)%Name)
                    stop 'ReadTimingDatabase - ModuleVegetation - ERR35'
                endif
            endif

        endif HF

    end subroutine ReadTimingParameters

    !--------------------------------------------------------------------------
    
    subroutine ReadGrowthDatabase(ivt)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound, VegetationFound
        integer                                     :: iflag, PlantType
        integer                                     :: GrowthObjEnterData
        integer                                     :: VegetationID
        integer                                     :: GrowthClientNumber
        type(T_GrowthDatabase), pointer             :: gdb !gdb is a pointer to the GrowthDataBase structure
        !Begin-----------------------------------------------------------------

        GrowthObjEnterData = 0
        !Open and save growth database
        call ConstructEnterData(GrowthObjEnterData, Me%GrowthDatabase, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR80'        
        
           
        VegetationFound = .false. 
            
doH1:   do 
            call ExtractBlockFromBuffer(GrowthObjEnterData,                                       &
                                        ClientNumber      = GrowthClientNumber,                   &
                                        block_begin       = '<begingrowthdatabase>',              &
                                        block_end         = '<endgrowthdatabase>',                &
                                        BlockFound        = DatabaseFound,                        &   
                                        STAT              = STAT_CALL)
HF1:        if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then
                
                call GetData(VegetationID, GrowthObjEnterData, iflag,                                                     &
                             SearchType     = FromBlock,                                                                  &
                             keyword        = 'VEGETATION_ID',                                                            &
                             ClientModule   = 'ModuleVegetation',                                                         &
                             STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR90'     
                if (iflag /= 1) then
                    write(*,*)'Missing VEGETATION_ID in Vegetation Type definition'
                    write(*,*)'Check vegetation definition in',Me%GrowthDatabase
                    stop 'ReadGrowthDatabase - ModuleVegetation - ERR30'
                endif
                                    
                if (VegetationID == Me%VegetationTypes(ivt)%VegetationID) then
                    
                    VegetationFound = .true.
                    
                    !Reads Growth Parameters
                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantType, GrowthObjEnterData,  iflag,               &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'PLANT_TYPE',                                                              &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR10'
                    if (iflag /= 1) then
                        write(*,*)'Missing plant type in growth parameters definition -', trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR15'
                    endif
                    
                    !Vegetation TYPE 0 is NO VEGETATION.
                    if (Me%VegetationTypes(ivt)%GrowthDatabase%PlantType == NotAPlant) cycle                    

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionN1, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_NITROGENFRACTION_N1',                                             &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR20'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal nitrogen fraction parameter in growth parameters definition -',           &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR25'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionN2, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_NITROGENFRACTION_N2',                                             &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR30'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal nitrogen fraction parameter in growth parameters definition -',           &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR35'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionN3, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 & 
                                 keyword        = 'OPTIMAL_NITROGENFRACTION_N3',                                             &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR40'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal nitrogen fraction parameter in growth parameters definition -',           &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR45'
                    endif
                    
                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionP1, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_PHOSPHORUSFRACTION_P1',                                           &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR50'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal phosphorus fraction parameter in growth parameters definition -',         &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR55'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionP2, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_PHOSPHORUSFRACTION_P2',                                           &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR60'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal phosphorus fraction parameter in growth parameters definition -',         &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR65'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantFractionP3, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_PHOSPHORUSFRACTION_P3',                                           &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR70'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal phosphorus fraction parameter in growth parameters definition -',         &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR75'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantBaseTemperature, GrowthObjEnterData,  iflag,    &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'BASE_TEMPERATURE',                                                        &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR80'
                    if (iflag /= 1) then
                        write(*,*)'Missing base temperature in growth parameters definition -',                              &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR85'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PlantOptimalTemperature, GrowthObjEnterData,  iflag, &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_TEMPERATURE',                                                     &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR90'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal temperature in growth parameters definition -',                           &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR95'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%ExtinctCoef, GrowthObjEnterData,  iflag,             &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'RADIATION_EXTINCTION_COEF',                                               &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR100'
                    if (iflag /= 1) then
                        write(*,*)'Missing radiation extiction coefficient in growth parameters definition -',               &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR105'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%BiomassEnergyRatio, GrowthObjEnterData,  iflag,      &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'BIOMASS_ENERGY_RATIO',                                                    &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR110'
                    if (iflag /= 1) then
                        write(*,*)'Missing biomass energy ratio in growth parameters definition -',                          &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR115'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%CO2ConcHigh, GrowthObjEnterData,  iflag,             &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'CO2_HIGH',                                                                &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR120'
                    if (iflag /= 1) then
                        write(*,*)'Missing high CO2 concentration in growth parameters definition -',                        &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR125'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%BiomassEnergyRatioHigh, GrowthObjEnterData,  iflag,  &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'BIOMASS_ENERGY_RATIO_HIGH',                                               &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR130'
                    if (iflag /= 1) then
                        write(*,*)'Missing high biomass energy ratio in growth parameters definition -',                     &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR135'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%RUEDeclineRate, GrowthObjEnterData,  iflag,          &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'RUE_DECLINE_RATE',                                                        &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR140'
                    if (iflag /= 1) then
                        write(*,*)'Missing RUE decline rate in growth parameters definition -',                              &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR145'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrLAIMax1, GrowthObjEnterData,  iflag,               &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_LAIMAXFRACTION_1',                                                &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR150'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal fraction of maximum LAI in growth parameters definition -',               &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR155'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrLAIMax2, GrowthObjEnterData,  iflag,               &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_LAIMAXFRACTION_2',                                                &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR160'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal fraction of maximum LAI in growth parameters definition -',               &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR165'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrGrow1, GrowthObjEnterData,  iflag,                 &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'GROWFRACTION_1',                                                          &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR170'
                    if (iflag /= 1) then
                        write(*,*)'Missing grow stage fraction asscoiated with LAI curve in growth parameters definition -', &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR175'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrGrow2, GrowthObjEnterData,  iflag,                 &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'GROWFRACTION_2',                                                          &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR180'
                    if (iflag /= 1) then
                        write(*,*)'Missing grow stage fraction asscoiated with LAI curve in growth parameters definition -', &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR185'
                    endif
                
                    !Computes LAIShape1 and LAIShape2
                    gdb => Me%VegetationTypes(ivt)%GrowthDatabase
                    call ComputeShapeCoefficients (gdb%FrLAIMax1, gdb%FrLAIMax2, gdb%FrGrow1, gdb%FrGrow2, gdb%LAIShape1, gdb%LAIShape2)

                    !LAI decline will not happen to EVERGREEN plants.
                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%Evergreen, GrowthObjEnterData,  iflag,               &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'EVERGREEN',                                                               &
                                 default        = .false.,                                                                   &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR190'
                                        
                    if (.not. Me%VegetationTypes(ivt)%GrowthDatabase%Evergreen) then
                        call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%FrGrowLAIDecline, GrowthObjEnterData,  iflag,        &
                                     SearchType     = FromBlock,                                                                 &
                                     keyword        = 'GROWFRACTION_LAIDECLINE',                                                 &
                                     ClientModule   = 'ModuleVegetation',                                                        &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR192'
                        if (iflag /= 1) then
                            write(*,*)'Missing grow stage fraction asscoiated with LAI curve in growth parameters definition -', &
                                       trim(Me%VegetationTypes(ivt)%Name)
                            stop 'ReadGrowthDatabase - ModuleVegetation - ERR193'
                        endif
                        
                        if (Me%ComputeOptions%UseLAIDecRate) then
                            call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%LAIDecRate, GrowthObjEnterData,  iflag,          &
                                         SearchType     = FromBlock,                                                             &
                                         keyword        = 'LAI_DEC_RATE',                                                        &
                                         ClientModule   = 'ModuleVegetation',                                                    &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR194'
                            if (iflag /= 1) then
                                write(*,*)'Missing LAI Decline rate in growth parameters definition -',                          &
                                           trim(Me%VegetationTypes(ivt)%Name)
                                stop 'ReadGrowthDatabase - ModuleVegetation - ERR195'
                            endif
                        endif                        
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%LAIMax, GrowthObjEnterData,  iflag,                  &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'LAI_MAX',                                                                 &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR200'
                    if (iflag /= 1) then
                        write(*,*)'Missing maximum LAI in growth parameters definition -',                                   &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR205'
                    endif
                    
                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%MaximumRootDepth, GrowthObjEnterData,  iflag,        &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'ROOT_DEPTH_MAX',                                                          &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR210'
                    if (iflag /= 1) then
                        write(*,*)'Missing maximum root depth in growth parameters definition -',                            &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR215'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%MaxCanopyHeight, GrowthObjEnterData,  iflag,         &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'CANOPY_HEIGHT_MAX',                                                       &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR220'
                    if (iflag /= 1) then
                        write(*,*)'Missing maximum canopy height in growth parameters definition -',                         &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR225'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%OptimalHarvestIndex, GrowthObjEnterData,  iflag,     &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'OPTIMAL_HARVEST_INDEX',                                                   &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR230'
                    if (iflag /= 1) then
                        write(*,*)'Missing optimal harvest index in growth parameters definition -',                         &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR235'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%MinimumHarvestIndex, GrowthObjEnterData,  iflag,     &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'MINIMUM_HARVEST_INDEX',                                                   &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR240'
                    if (iflag /= 1) then
                        write(*,*)'Missing minimum harvest index in growth parameters definition -',                         &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR245'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%NitrogenFractionInYeld, GrowthObjEnterData,  iflag,  &
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'YELD_NITROGENFRACTION',                                                   &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR250'
                    if (iflag /= 1) then
                        write(*,*)'Missing nitrogen fraction in yeld in growth parameters definition -',                     &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR255'
                    endif

                    call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%PhosphorusFractionInYeld, GrowthObjEnterData,  iflag,&
                                 SearchType     = FromBlock,                                                                 &
                                 keyword        = 'YELD_PHOSPHORUSFRACTION',                                                 &
                                 ClientModule   = 'ModuleVegetation',                                                        &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR260'
                    if (iflag /= 1) then
                        write(*,*)'Missing phosphorus fraction in yeld in growth parameters definition -',                   &
                                   trim(Me%VegetationTypes(ivt)%Name)
                        stop 'ReadGrowthDatabase - ModuleVegetation - ERR265'
                    endif

                    PlantType = Me%VegetationTypes(ivt)%GrowthDatabase%PlantType
                    if (PlantType == Trees) then
                    
                        call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%TreeYearsToMaturity, GrowthObjEnterData,  iflag, &
                                     SearchType     = FromBlock,                                                             &
                                     keyword        = 'TREE_YEARSTOMATURITY',                                                &
                                     ClientModule   = 'ModuleVegetation',                                                    &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR270'
                        if (iflag /= 1) then
                            write(*,*)'Missing tree years too maturity in growth parameters definition -',                   &
                                       trim(Me%VegetationTypes(ivt)%Name)
                            stop 'ReadGrowthDatabase - ModuleVegetation - ERR275'
                        endif

                        call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%TreeMaximumBiomass, GrowthObjEnterData,  iflag,  &
                                     SearchType     = FromBlock,                                                             &
                                     keyword        = 'TREE_MAXIMUMBIOMASS',                                                 &
                                     ClientModule   = 'ModuleVegetation',                                                    &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR280'
                        if (iflag /= 1) then
                            write(*,*)'Missing tree maximum biomass in growth parameters definition -', &
                                       trim(Me%VegetationTypes(ivt)%Name)
                            stop 'ReadGrowthDatabase - ModuleVegetation - ERR285'
                        endif
                    
                    endif
                    
                    if (Me%ComputeOptions%Dormancy) then
                        call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%BiomassFracRemovedInDormancy, GrowthObjEnterData,   &
                                     iflag,                                                                                     &
                                     SearchType     = FromBlock,                                                                &
                                     keyword        = 'BIOMASS_FRAC_REMOVED_DORMANCY',                                          &
                                     ClientModule   = 'ModuleVegetation',                                                       &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR290'
                        if (iflag /= 1) then
                            write(*,*)'Missing biomass fraction removed in dormancy in growth parameters definition -',     &
                                       trim(Me%VegetationTypes(ivt)%Name)
                            stop 'ReadGrowthDatabase - ModuleVegetation - ERR295'
                        endif
                
                        call GetData(Me%VegetationTypes(ivt)%GrowthDatabase%LAIMinDormant, GrowthObjEnterData,  iflag,         &
                                     SearchType     = FromBlock,                                                     &
                                     keyword        = 'LAI_MIN_DORMANCY',                                                   &
                                     ClientModule   = 'ModuleVegetation',                                                   &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR300'
                        if (iflag /= 1) then
                            write(*,*)'Missing LAI minimum in dormancy in growth parameters definition -',                  &
                                       trim(Me%VegetationTypes(ivt)%Name)
                            stop 'ReadGrowthDatabase - ModuleVegetation - ERR305'
                        endif
                    endif
                    
                    exit doH1
                    
                endif
                
            else
                
                call Block_Unlock(GrowthObjEnterData, GrowthClientNumber, STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR120'
               
                exit doH1 
                          
            endif HF1
            
        enddo doH1
             
        if (.not. VegetationFound) then
            write(*,*)
            write(*,*) 'Vegetation ID not found:', Me%VegetationTypes(ivt)%ID
            write(*,*) 'Check that ID is defined in database:',Me%GrowthDatabase
            stop 'ReadGrowthDatabase - ModuleVegetation - ERR130' 
        endif           
    
    
        call KillEnterData(GrowthObjEnterData, STAT = STAT_CALL)         
        if (STAT_CALL /= SUCCESS_) stop 'ReadGrowthDatabase - ModuleVegetation - ERR140'                

    end subroutine ReadGrowthDatabase

    !--------------------------------------------------------------------------

    subroutine ReadHarvestKillParameters (ivt, ClientNumber,ParameterObjEnterData)

        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        integer                                     :: ParameterObjEnterData
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock(ParameterObjEnterData,                                           &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<beginharvestkillparameters>',  &
                                    block_end         = '<endharvestkillparameters>',    &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

            !Reads Parameters
            
            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%HarvestKillJulianDay, ParameterObjEnterData,  iflag,    &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'HARVESTKILL_JULIANDAY',                                                    &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR100'     

            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%HarvestKillPlantHU, ParameterObjEnterData,  iflag,&
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'HARVESTKILL_PLANTHU',                                                      &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR110'     

            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%HarvestJulianDay, ParameterObjEnterData,  iflag,  &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'HARVEST_JULIANDAY',                                                        &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR120'     


!                allocate (Me%VegetationTypes(ivt)%HarvestKillDatabase%HarvestPlantHU(1))

            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%HarvestPlantHU, ParameterObjEnterData,  iflag, &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'HARVEST_PLANTHU',                                                          &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR130'     

            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%KillJulianDay, ParameterObjEnterData,  iflag,     &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'KILL_JULIANDAY',                                                           &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR140'     

            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%KillPlantHU, ParameterObjEnterData,  iflag,       &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'KILL_PLANTHU',                                                             &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR150'     
            
            call GetData(Me%VegetationTypes(ivt)%HarvestKillDatabase%HarvestEfficiency, ParameterObjEnterData,  iflag, &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = 1.,                                                                         &
                         keyword        = 'HARVEST_EFFICIENCY',                                                       &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadHarvestKillDatabase - ModuleVegetation - ERR160'     
           
        endif HF


    end subroutine ReadHarvestKillParameters

    !--------------------------------------------------------------------------

    subroutine ReadGrazingParameters (ivt, ClientNumber,ParameterObjEnterData)

        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        integer                                     :: ParameterObjEnterData
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock(ParameterObjEnterData,                            &
                                    ClientNumber      = ClientNumber,                &
                                    block_begin       = '<begingrazeparameters>',   &
                                    block_end         = '<endgrazeparameters>',     &
                                    BlockInBlockFound = DatabaseFound,               &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then

            !Reads Parameters
            
            call GetData(Me%VegetationTypes(ivt)%GrazingDatabase%GrazingStartJulianDay, ParameterObjEnterData,  iflag,   &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'GRAZING_START_JULIANDAY',                                                  &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'GrazingDatabase - ModuleVegetation - ERR40'     
            
!            allocate (Me%VegetationTypes(ivt)%GrazingDatabase%GrazingStartPlantHU(1))
            
            call GetData(Me%VegetationTypes(ivt)%GrazingDatabase%GrazingStartPlantHU, ParameterObjEnterData,  iflag,  &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = -99.,                                                                       &
                         keyword        = 'GRAZING_START_PLANTHU',                                                    &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrazingDatabase - ModuleVegetation - ERR50'     

            call GetData(Me%VegetationTypes(ivt)%GrazingDatabase%GrazingDays, ParameterObjEnterData,  iflag,       &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = 0,                                                                          &
                         keyword        = 'GRAZING_DAYS',                                                             &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrazingDatabase - ModuleVegetation - ERR60'     
            
            call GetData(Me%VegetationTypes(ivt)%GrazingDatabase%GrazingMinimumBiomass, ParameterObjEnterData,  iflag,   &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = 0.,                                                                         &
                         keyword        = 'MINIMUM_BIOMASS_FOR_GRAZING',                                              &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrazingDatabase - ModuleVegetation - ERR70'      

            call GetData(Me%VegetationTypes(ivt)%GrazingDatabase%GrazingBiomass, ParameterObjEnterData,  iflag,    &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = 0.,                                                                         &
                         keyword        = 'GRAZING_BIOMASS',                                                          &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrazingDatabase - ModuleVegetation - ERR80'

            call GetData(Me%VegetationTypes(ivt)%GrazingDatabase%TramplingBiomass, ParameterObjEnterData,  iflag,  &
                         SearchType     = FromBlockInBlock,                                                           &
                         Default        = 0.,                                                                         &
                         keyword        = 'TRAMPLING_BIOMASS',                                                        &
                         ClientModule   = 'ModuleVegetation',                                                         &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadGrazingDatabase - ModuleVegetation - ERR90'     


        endif HF


    end subroutine ReadGrazingParameters

    !--------------------------------------------------------------------------
    
    
    subroutine ReadFertilizationParameters (ivt, ClientNumber,ParameterObjEnterData)

        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        integer                                     :: ParameterObjEnterData
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound, AutofertilizationFound
        logical                                     :: FertilizationAppFound
        integer                                     :: iflag
        integer                                     :: FertilizerObjEnterData
        integer                                     :: FertFileClientNumber
        logical                                     :: FertilizerFound
        integer                                     :: FertilizerID
        integer                                     :: NumberOfFertilizerApps, FertApp
        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock(ParameterObjEnterData,                                     &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<beginfertilizationparameters>',     &
                                    block_end         = '<endfertilizationparameters>',       &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF:     if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then


            !Autofertilization block
            call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                        ClientNumber      = ClientNumber,                         &
                                        block_begin       = '<<beginautofertilization>>',         &
                                        block_end         = '<<endautofertilization>>',           &
                                        BlockInBlockInBlockFound = AutofertilizationFound,        &   
                                        STAT              = STAT_CALL)
            if (STAT_CALL.NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR10' 
           
           
            !Fertilization app block
            call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                        ClientNumber      = ClientNumber,                         &
                                        block_begin       = '<<beginfertilizerapp>>',             &
                                        block_end         = '<<endfertilizerapp>>',               &
                                        BlockInBlockInBlockFound = FertilizationAppFound,         &   
                                        STAT              = STAT_CALL)
            if (STAT_CALL.NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR10' 
           
            if (AutofertilizationFound .and. FertilizationAppFound) then
                write(*,*)
                write(*,*)'Error! Autofertilization and Fertilization applications both active'
                write(*,*)'in agricultural practice', Me%VegetationTypes(ivt)%VegetationID
                write(*,*)'Please check that only one block is present.'
                stop 'ReadFertilizationDatabase - ModuleVegetation - ERR20' 
            endif
            
            if (AutofertilizationFound) then
                Me%ComputeOptions%Autofertilization = .true.
                Me%ComputeOptions%FertilizationApp  = .false.
                !one continuous application
                NumberOfFertilizerApps = 1

            else
                Me%ComputeOptions%Autofertilization = .false.
                Me%ComputeOptions%FertilizationApp  = .true.

                call RewindBlock(ParameterObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR25' 

                call ExtractBlockFromBlock(ParameterObjEnterData,                                     &
                                            ClientNumber      = ClientNumber,                         &
                                            block_begin       = '<beginfertilizationparameters>',     &
                                            block_end         = '<endfertilizationparameters>',       &
                                            BlockInBlockFound = DatabaseFound,                        &   
                                            STAT              = STAT_CALL)
                if (STAT_CALL.NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR40' 
            
                !Count number of applications
                NumberOfFertilizerApps = 0
doH:            do 
                    call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                                ClientNumber      = ClientNumber,                         &
                                                block_begin       = '<<beginfertilizerapp>>',             &
                                                block_end         = '<<endfertilizerapp>>',               &
                                                BlockInBlockInBlockFound = FertilizationAppFound,         &   
                                                STAT              = STAT_CALL)
HF2:                if (STAT_CALL == SUCCESS_ .and. FertilizationAppFound) then
                        NumberOfFertilizerApps = NumberOfFertilizerApps + 1
                        
                    else
                        
                        exit doH  
                                  
                    endif HF2
                    
                enddo doH
            
            endif
            
HF6:        if (NumberOfFertilizerApps .gt. 0) then
            
                allocate (Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(NumberOfFertilizerApps))
                Me%VegetationTypes(ivt)%FertilizerDatabase%NumberFertilizerApps = NumberOfFertilizerApps                        


                call RewindBlock(ParameterObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR30' 

                call ExtractBlockFromBlock(ParameterObjEnterData,                                     &
                                            ClientNumber      = ClientNumber,                         &
                                            block_begin       = '<beginfertilizationparameters>',     &
                                            block_end         = '<endfertilizationparameters>',       &
                                            BlockInBlockFound = DatabaseFound,                        &   
                                            STAT              = STAT_CALL)
                if (STAT_CALL.NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR40' 


do2:            do FertApp = 1, NumberOfFertilizerApps
                    
                    !Search fo the block depending on method
                    if (Me%ComputeOptions%Autofertilization) then

                        !Autofertilization block
                        call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                                    ClientNumber      = ClientNumber,                         &
                                                    block_begin       = '<<beginautofertilization>>',         &
                                                    block_end         = '<<endautofertilization>>',           &
                                                    BlockInBlockInBlockFound = AutofertilizationFound,        &   
                                                    STAT              = STAT_CALL)
                        if (STAT_CALL.NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR50' 

                        !AutoFertilization
                        call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%NTreshold,      &
                                     ParameterObjEnterData, iflag,                                                           &
                                     SearchType     = FromBlockInBlockInBlock,                                               &
                                     Default        = -99.,                                                                  &
                                     keyword        = 'NITROGEN_TRESHOLD',                                                   &
                                     ClientModule   = 'ModuleVegetation',                                                    &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR230'     
                    
                        call GetData(                                                                                        &
                             Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%NitrogenApplicationMax, &
                                     ParameterObjEnterData,  iflag,                                                          &
                                     SearchType     = FromBlockInBlockInBlock,                                               &
                                     Default        = -99.,                                                                  &
                                     keyword        = 'NITROGEN_APPLICATION_MAX',                                            &
                                     ClientModule   = 'ModuleVegetation',                                                    &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR240'     

                        call GetData(                                                                                        &
                                  Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%NitrogenAnnualMax, &
                                     ParameterObjEnterData,  iflag,                                                          &
                                     SearchType     = FromBlockInBlockInBlock,                                               &
                                     Default        = -99.,                                                                  &
                                     keyword        = 'NITROGEN_ANNUAL_MAX',                                                 &
                                     ClientModule   = 'ModuleVegetation',                                                    &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR250'     

                        call GetData(                                                                                        &
                                 Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%ExplicitPhosphorus, &
                                     ParameterObjEnterData,  iflag,                                                          &
                                     SearchType     = FromBlockInBlockInBlock,                                               &
                                     Default        = .false.,                                                               &
                                     keyword        = 'EXPLICIT_PHOSPHORUS',                                                 &
                                     ClientModule   = 'ModuleVegetation',                                                    &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR260'     

                        if (Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%ExplicitPhosphorus) then

                            call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%PTreshold, &
                                         ParameterObjEnterData,  iflag,                                                     &
                                         SearchType     = FromBlockInBlockInBlock,                                          &
                                         Default        = -99.,                                                             &
                                         keyword        = 'PHOSPHORUS_TRESHOLD',                                            &
                                         ClientModule   = 'ModuleVegetation',                                               &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR270'     
                    
                            call GetData(                                                                                   &
                          Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%PhosphorusApplicationMax, &
                                         ParameterObjEnterData,                                                             &
                                         iflag,                                                                             &
                                         SearchType     = FromBlockInBlockInBlock,                                          &
                                         Default        = -99.,                                                             &
                                         keyword        = 'PHOSPHORUS_APPLICATION_MAX',                                     &
                                         ClientModule   = 'ModuleVegetation',                                               &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR280'     

                            call GetData(                                                                                         &
                                     Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%PhosphorusAnnualMax, &
                                         ParameterObjEnterData,  iflag,                                                           &
                                         SearchType     = FromBlockInBlockInBlock,                                                &
                                         Default        = -99.,                                                                   &
                                         keyword        = 'PHOSPHORUS_ANNUAL_MAX',                                                &
                                         ClientModule   = 'ModuleVegetation',                                                     &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR290'   
                            
                        endif

                        call GetData(                                                                                         &
                                 Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Auto%NStress,             &
                                     ParameterObjEnterData,  iflag,                                                           &
                                     SearchType     = FromBlockInBlockInBlock,                                                &
                                     Default        = NTarget_,                                                               &
                                     keyword        = 'N_STRESS_TYPE',                                                        &
                                     ClientModule   = 'ModuleVegetation',                                                     &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR300'              
                        
                    elseif (Me%ComputeOptions%FertilizationApp) then

                        !FertilizationApplication block
                        call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                                    ClientNumber      = ClientNumber,                         &
                                                    block_begin       = '<<beginfertilizerapp>>',             &
                                                    block_end         = '<<endfertilizerapp>>',               &
                                                    BlockInBlockInBlockFound = FertilizationAppFound,         &   
                                                    STAT              = STAT_CALL)
                        if (STAT_CALL.NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR100' 
           
                        !day for application on beggining day for application in case of continuous fertilization
                        call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppJDay,        &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = -99.,                                                                       &
                                     keyword        = 'FERTILIZER_APPLICATION_JDAY',                                              &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR305'    
                        
                        !hu for application on beggining hu for application in case of continuous fertilization                    
                        call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppHU,          &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = -99.,                                                                       &
                                     keyword        = 'FERTILIZER_APPLICATION_HU',                                                &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR310'     
                        
                        !Amount applied - kg/ha.day
                        call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppAmount,      &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = -99.,                                                                       &
                                     keyword        = 'FERTILIZER_APPLICATION_KG_HA',                                             &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR320'    

                        !Continuous fertilization?
                        call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationOn,      &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = .false.,                                                                    &
                                     keyword        = 'FERTILIZER_CONT_ON',                                                       &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR340'       
                        
                        
                        if (Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationOn) then
                            !Continuous fertilization during how many days
                            call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationDays,&
                                         ParameterObjEnterData,  iflag,                                                           &
                                         SearchType     = FromBlockInBlockInBlock,                                                &
                                         Default        = 0,                                                                      &
                                         keyword        = 'FERTILIZER_CONT_DAYS',                                                 &
                                         ClientModule   = 'ModuleVegetation',                                                     &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR340'    
                        endif

                    
                    endif

                    !Reads Fertilizer ID in that block already open
                    call GetData(Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerID, &
                                 ParameterObjEnterData,  iflag,                                                   &
                                 SearchType     = FromBlockInBlockInBlock,                                        &
                                 keyword        = 'FERTILIZER_ID',                                                &
                                 ClientModule   = 'ModuleVegetation',                                             &
                                 STAT           = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR110'     
                    if (iflag /= 1) then
                        write(*,*)'Missing FERTILIZER_ID in Vegetation parameter file'
                        stop 'ReadFertilizationDatabase - ModuleVegetation - ERR120'
                    endif    
                 
                    !Read fertilizer database to check that fertilizer exists and read fertilizer data
                    FertilizerObjEnterData = 0
                    !Open and save fertilizer database file
                    call ConstructEnterData(FertilizerObjEnterData, Me%FertilizerDatabase, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR120'  
                    
do3:                do 
                        FertilizerFound = .false.
                        call ExtractBlockFromBuffer(FertilizerObjEnterData,                                   &
                                                    ClientNumber      = FertFileClientNumber,                 &
                                                    block_begin       = '<beginFertilizer>',                  &
                                                    block_end         = '<endFertilizer>',                    &
                                                    BlockFound        = DatabaseFound,                        &   
                                                    STAT              = STAT_CALL)
HF3:                    if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then
                            
                            call GetData(FertilizerID, FertilizerObjEnterData, iflag,                                            &
                                         SearchType     = FromBlock,                                                             &
                                         keyword        = 'FERTILIZER_ID',                                                       &
                                         ClientModule   = 'ModuleVegetation',                                                    &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR130'     
                            if (iflag /= 1) then
                                write(*,*)'Missing FERTILIZER_ID in Fertilizer database:', Me%FertilizerDatabase
                                stop 'ReadFertilizationDatabase - ModuleVegetation - ERR140'
                            endif  
                            
                            Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%FertilizerID =        &
                            FertilizerID
                                              
                            if (FertilizerID == &
                                Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerID) then
                                
                                FertilizerFound = .true.
                                
                                call GetData(                                                                                   &
                        Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%MineralNFracInFertilizer, &
                                             FertilizerObjEnterData,                                                            &
                                             iflag,                                                                             &
                                             SearchType     = FromBlock,                                                        &
                                             Default        = -99.,                                                             &
                                             keyword        = 'MINERAL_N_FRACTION_IN_FERTILIZER',                               &
                                             ClientModule   = 'ModuleVegetation',                                               &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR150'     
                                
                                call GetData(                                                                                   &
                       Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicNFracInFertilizer, &
                                             FertilizerObjEnterData,                                                            &
                                             iflag,                                                                             &
                                             SearchType     = FromBlock,                                                        &
                                             Default        = -99.,                                                             &
                                             keyword        = 'ORGANIC_N_FRACTION_IN_FERTILIZER',                               &
                                             ClientModule   = 'ModuleVegetation',                                               &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR160'     

                                call GetData(                                                                                   &
                          Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%AmmoniaFracInMineralN, &
                                             FertilizerObjEnterData,                                                            &
                                             iflag,                                                                             &
                                             SearchType     = FromBlock,                                                        &
                                             Default        = -99.,                                                             &
                                             keyword        = 'AMMONIA_FRACTION_IN_MINERAL_N',                                  &
                                             ClientModule   = 'ModuleVegetation',                                               &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR170'     
                                
                                call GetData(                                                                                   & 
                       Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%MineralPFracInFertilizer, &
                                             FertilizerObjEnterData,                                                            &
                                             iflag,                                                                             &
                                             SearchType     = FromBlock,                                                        &
                                             Default        = -99.,                                                             &
                                             keyword        = 'MINERAL_P_FRACTION_IN_FERTILIZER',                               &
                                             ClientModule   = 'ModuleVegetation',                                               &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR180'      

                                call GetData(                                                                                    &
                        Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicPFracInFertilizer, &
                                             FertilizerObjEnterData,                                                             &
                                             iflag,                                                                              &
                                             SearchType     = FromBlock,                                                         &
                                             Default        = -99.,                                                              &
                                             keyword        = 'ORGANIC_P_FRACTION_IN_FERTILIZER',                                &
                                             ClientModule   = 'ModuleVegetation',                                                &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR190'  

                                call GetData(                                                                                    &
                  Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%FertilizerFracApplyedInSurface, &
                                             FertilizerObjEnterData,  iflag,                                                     &
                                             SearchType     = FromBlock,                                                         &
                                             Default        = -99.,                                                              &
                                             keyword        = 'FERTILIZER_FRACTION_IN_SURFACE',                                  &
                                             ClientModule   = 'ModuleVegetation',                                                &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR200'     

                                call GetData(                                                                                    &
                       Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicFracParticulate,     &
                                             FertilizerObjEnterData,  iflag,                                                     &
                                             SearchType     = FromBlock,                                                         &
                                             Default        = 0.,                                                              &
                                             keyword        = 'ORGANIC_FRACTION_PARTICULATE',                                    &
                                             ClientModule   = 'ModuleVegetation',                                                &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR210'     
                                
                                !Check if particulate material to runoff properties (fluff layer)                
                                if (                                                                                             &
                Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicFracParticulate .gt. 0.) then
                                    Me%ComputeOptions%VegParticFertilization = .true.
                                endif
                                
                                exit do3
                                
                            endif
                            
                        else
                            
                            call Block_Unlock(FertilizerObjEnterData, FertFileClientNumber, STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadFertilizationDatabase - ModuleVegetation - ERR210'
                           
                            exit do3 
                                      
                        endif HF3
                        
                    enddo do3
                     
                    if (.not. FertilizerFound) then
                        write(*,*)
                        write(*,*) 'Fertilizer ID not found:'
                        write(*,*)  Me%VegetationTypes(ivt)%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerID
                        write(*,*) 'Check fertilizer database in: ', Me%FertilizerDatabase
                        stop 'ReadPesticideDatabase - ModuleVegetation - ERR220' 
                    endif   
                    
                enddo do2
            
            endif HF6
          
        endif HF


    end subroutine ReadFertilizationParameters

    !--------------------------------------------------------------------------        


    subroutine ReadPesticideParameters (ivt, ClientNumber,ParameterObjEnterData)

        !Arguments-------------------------------------------------------------
        integer                                     :: ivt
        integer                                     :: ClientNumber
        integer                                     :: ParameterObjEnterData
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: DatabaseFound, ApplicationFound
        integer                                     :: iflag, NumberOfPesticideApps
        integer, dimension(:), pointer              :: PesticideList
        integer                                     :: NumberOfUniquePesticides, UniquePest
 !       character (Len = StringLength)              :: PesticideDatabase
        integer                                     :: PesticideObjEnterData, PesticideID, PestApp
        logical                                     :: PesticideFound
        integer                                     :: PestAppllication
        integer                                     :: PestFileClientNumber

        !Begin-----------------------------------------------------------------
        
        call ExtractBlockFromBlock(ParameterObjEnterData,                                     &
                                    ClientNumber      = ClientNumber,                         &
                                    block_begin       = '<beginpesticideparameters>',         &
                                    block_end         = '<endpesticideparameters>',           &
                                    BlockInBlockFound = DatabaseFound,                        &   
                                    STAT              = STAT_CALL)
HF4:    if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then        

           
            !Count number of applications
            NumberOfPesticideApps = 0
doH:        do 
                call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                            ClientNumber      = ClientNumber,                         &
                                            block_begin       = '<<beginpesticideapp>>',                &
                                            block_end         = '<<endpesticideapp>>',                  &
                                            BlockInBlockInBlockFound = ApplicationFound,              &   
                                            STAT              = STAT_CALL)
HF:             if (STAT_CALL == SUCCESS_ .and. ApplicationFound) then
                    NumberOfPesticideApps = NumberOfPesticideApps + 1
                    
                else
                    
                    exit doH  
                              
                endif HF
                
            enddo doH
            
        endif HF4
            
        
HF6:    if (NumberOfPesticideApps .gt. 0) then
        
            allocate (Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(NumberOfPesticideApps))
            Me%VegetationTypes(ivt)%PesticideDatabase%NumberPesticideApps = NumberOfPesticideApps
            !Read pesticide application info    
            !Rewind block so it can search again from beggining father block
            call RewindBlock(ParameterObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR30' 

            call ExtractBlockFromBlock(ParameterObjEnterData,                                     &
                                        ClientNumber      = ClientNumber,                         &
                                        block_begin       = '<beginpesticideparameters>',         &
                                        block_end         = '<endpesticideparameters>',           &
                                        BlockInBlockFound = DatabaseFound,                        &   
                                        STAT              = STAT_CALL)
HF5:        if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then        

               
                PestApp = 0
                
doE:            do
                    call ExtractBlockFromBlockFromBlock(ParameterObjEnterData,                            &
                                                ClientNumber      = ClientNumber,                         &
                                                block_begin       = '<<beginpesticideapp>>',                &
                                                block_end         = '<<endpesticideapp>>',                  &
                                                BlockInBlockInBlockFound = ApplicationFound,              &   
                                                STAT              = STAT_CALL)
HF3:                if (STAT_CALL == SUCCESS_ .and. ApplicationFound) then   
                        
                        PestApp = PestApp + 1         
                        
                        call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideID,               &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     keyword        = 'PESTICIDE_ID',                                                             &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR40'    
                        if (iflag /= 1) then
                            write(*,*)'Missing PESTICIDE_ID in Vegetation parameter file'
                            stop 'ReadPesticideDatabase - ModuleVegetation - ERR50'
                        endif  

                        call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideAppJDay,           &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = -99.,                                                                       &
                                     keyword        = 'PESTICIDE_APPLICATION_JDAY',                                               &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR45'    
                                            
                        call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideAppHU,             &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = -99.,                                                                       &
                                     keyword        = 'PESTICIDE_APPLICATION_HU',                                                 &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR50'     

                        call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideAppAmount,         &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = -99.,                                                                       &
                                     keyword        = 'PESTICIDE_APPLICATION_KG_HA',                                              &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR60'     

                        call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%ContPesticideON,            &
                                     ParameterObjEnterData,  iflag,                                                               &
                                     SearchType     = FromBlockInBlockInBlock,                                                    &
                                     Default        = .false.,                                                                    &
                                     keyword        = 'PESTICIDE_CONT_ON',                                                        &
                                     ClientModule   = 'ModuleVegetation',                                                         &
                                     STAT           = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR70'     

                        if (Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%ContPesticideOn) then
                            !Continuous pesticide during how many days
                            call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%ContPesticideDays,      &
                                         ParameterObjEnterData,  iflag,                                                           &
                                         SearchType     = FromBlockInBlockInBlock,                                                &
                                         Default        = 0,                                                                      &
                                         keyword        = 'PESTICIDE_CONT_DAYS',                                                 &
                                         ClientModule   = 'ModuleVegetation',                                                     &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR75'    
                        endif


                    else
                    
                        exit  doE
                              
                    
                    endif HF3
                
                enddo doE
                
               
                !First Count number of different pesticides - only for pesticides inside vegetationtype
                
                allocate (PesticideList(NumberOfPesticideApps))
                NumberOfUniquePesticides = 1
                PesticideList(1) = Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(1)%PesticideID
                
                if (NumberOfPesticideApps .gt. 1) then
                    
                    !Go from second pesticide application to all
do1:                do PestAppllication = 2, NumberOfPesticideApps
                       
                        !go to every unique pesticide list position
do2:                    do UniquePest = 1,  NumberOfUniquePesticides
                            if ((Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestAppllication)%PesticideID)   &
                                 == PesticideList(UniquePest)) then
                                !found equal, go to next in application list    
                                cycle do1
                            endif
                        enddo do2
                        
                        !not found in unique list, so add it
                        NumberOfUniquePesticides = NumberOfUniquePesticides + 1
                        PesticideList(NumberOfUniquePesticides) =                                      &
                        Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestAppllication)%PesticideID
                    enddo do1

                endif                
                
                !Now get the info needed about every pesticide in pesticide list
                
                PesticideObjEnterData = 0
                !Open and save pesticide database
                call ConstructEnterData(PesticideObjEnterData, Me%PesticideDatabase, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR80'        
                
                
                !Check for every application that pesticide is in database and get pesticide info
doP1:           do PestApp = 1, NumberOfPesticideApps  
                    
                    PesticideFound = .false. 
                    
                    !look again from the beggining
                    call RewindBuffer(PesticideObjEnterData, STAT = STAT_CALL)
                    if(STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR85'
                    
doH1:               do 
                        call ExtractBlockFromBuffer(PesticideObjEnterData,                                    &
                                                    ClientNumber      = PestFileClientNumber,                 &
                                                    block_begin       = '<beginPesticide>',                   &
                                                    block_end         = '<endPesticide>',                     &
                                                    BlockFound        = DatabaseFound,                        &   
                                                    STAT              = STAT_CALL)
HF1:                    if (STAT_CALL == SUCCESS_ .and. DatabaseFound) then
                            
                            call GetData(PesticideID, PesticideObjEnterData, iflag,                            &
                                         SearchType     = FromBlock,                                           &
                                         keyword        = 'PESTICIDE_ID',                                      &
                                         ClientModule   = 'ModuleVegetation',                                  &
                                         STAT           = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR90'     
                            if (iflag /= 1) then
                                write(*,*)'Missing PESTICIDE_ID in Pesticide database:', Me%PesticideDatabase
                                stop 'ReadPesticideDatabase - ModuleVegetation - ERR100'
                            endif       
                                               
                            if (PesticideID == Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideID) then
                                
                                PesticideFound = .true.
                                
                                call GetData(Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideName,  &
                                             EnterDataID    = PesticideObjEnterData,                                          &
                                             flag           = iflag,                                                          &
                                             SearchType     = FromBlock,                                                      &
                                             keyword        = 'PESTICIDE_NAME',                                               &
                                             ClientModule   = 'ModuleVegetation',                                             &
                                             STAT           = STAT_CALL)
                                if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR110'

                                exit doH1
                                
                            endif
                            
                        else
                            
                            call Block_Unlock(PesticideObjEnterData, PestFileClientNumber, STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR120'
                           
                            exit doH1 
                                      
                        endif HF1
                        
                    enddo doH1
                     
                     if (.not. PesticideFound) then
                        write(*,*)
                        write(*,*) 'Pesticide ID not found:'
                        write(*,*) Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideID
                        write(*,*) 'Check that ID in pesticide database:', Me%PesticideDatabase
                        stop 'ReadPesticideDatabase - ModuleVegetation - ERR130' 
                    endif           
                
                enddo doP1
                
                call KillEnterData(PesticideObjEnterData, STAT = STAT_CALL)         
                if (STAT_CALL /= SUCCESS_) stop 'ReadPesticideDatabase - ModuleVegetation - ERR140'        
                
                !get global list of unique pesticides (all vegetation types) - for fluxes
                call GetUniquePesticide(ivt, NumberOfPesticideApps)
                      
            endif HF5
        endif  HF6                  
        
    end subroutine ReadPesticideParameters

    !--------------------------------------------------------------------------        

    subroutine GetUniquePesticide(ivt, NumberOfPesticideApps)
    
        !Arguments-------------------------------------------------------------
        integer                                          :: ivt
        integer                                          :: NumberOfPesticideApps
        !Local-----------------------------------------------------------------
        integer                                          :: PestApp, GlobalUniquePest
        !Begin-----------------------------------------------------------------


        !Do the same process for unique pesticide for all the vegetation types 
        !needed for pesticide fluxes
        
        !Go from first pesticide app to all 
do3:    do PestApp = 1, NumberOfPesticideApps
           
            !go to every unique pesticide list position now global list
do4:        do GlobalUniquePest = 1,  Me%Fluxes%Pesticides%UniquePesticides
                if (Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideID      &
                     == Me%PesticideListID(GlobalUniquePest)) then
                    !found equal, go to next in vegetation type unique list    
                    cycle do3
                endif
            enddo do4
            
            !not found in global unique list, so add it
            Me%Fluxes%Pesticides%UniquePesticides = Me%Fluxes%Pesticides%UniquePesticides + 1
            Me%PesticideListID(Me%Fluxes%Pesticides%UniquePesticides) =                           &
             Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideID
            Me%PesticideListName(Me%Fluxes%Pesticides%UniquePesticides) =                           &
             Me%VegetationTypes(ivt)%PesticideDatabase%PesticideApps(PestApp)%PesticideName
        enddo do3
        
    
    end subroutine GetUniquePesticide
    
    !--------------------------------------------------------------------------        

!    subroutine ReadInitialFile
!
!        !Arguments-------------------------------------------------------------
!
!        !Local-----------------------------------------------------------------
!        real                                        :: Year_File, Month_File, Day_File
!        real                                        :: Hour_File, Minute_File, Second_File
!        integer                                     :: InitialFile
!        type (T_Time)                               :: BeginTime, EndTimeFile, EndTime
!        real                                        :: DT_error
!        integer                                     :: STAT_CALL
!
!        !----------------------------------------------------------------------
!
!        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleVegetation - ERR01'
!
!        open(Unit = InitialFile, File = Me%Files%InitialFile, Form = 'UNFORMATTED',     &
!             status = 'OLD', IOSTAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleVegetation - ERR02'
!
!        !Reads Date
!        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
!        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)
!
!        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, EndTime = EndTime, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleVegetation - ERR03'
!        
!        DT_error = EndTimeFile - BeginTime
!
!        !Avoid rounding erros - Frank 08-2001
!        if (abs(DT_error) >= 0.01) then
!            
!            write(*,*) 'The end time of the previous run is different from the start time of this run'
!            write(*,*) 'Date in the file'
!            write(*,*) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
!            write(*,*) 'DT_error', DT_error
!            if (Me%ComputeOptions%StopOnWrongDate) stop 'ReadInitialFile - ModuleVegetation - ERR04'   
!
!        endif
!
!
!        read(InitialFile)Me%IsPlantGrowing
!        read(InitialFile, end=10)Me%IsPlantDormant
!        
!   10   continue
!
!
!        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialFile - ModuleVegetation - ERR05'  
!      
!    end subroutine ReadInitialFile

    !-------------------------------------------------------------------------- 
    
    subroutine ReadInitialHDF

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        integer, dimension(:,:), pointer            :: PlantGrowingInteger, PlantDormantInteger
        character (Len = StringLength)              :: ConvertType 

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


        inquire (FILE=trim(Me%Files%InitialFile)//"5", EXIST = Exist)

cd0:    if (Exist) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialFile)//"5",&
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR01'


            ! Reads from HDF file the Property value and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR02'



            allocate (PlantGrowingInteger(Me%Worksize%ILB:Me%Worksize%IUB, Me%Worksize%JLB:Me%Worksize%JUB))
            
            PlantGrowingInteger (:,:) = 0
            
            call HDF5ReadData   (ObjHDF5, "/Results/"//"IsPlantGrowing",                 &
                                 "IsPlantGrowing",                                       &
                                 Array2D = PlantGrowingInteger,                          &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR03'


            ConvertType = "IntegerToLogical"
            call ConvertLogicalInteger (Me%IsPlantGrowing, PlantGrowingInteger, ConvertType)
            
            deallocate (PlantGrowingInteger)
            
            if (Me%ComputeOptions%Dormancy) then

                allocate (PlantDormantInteger(Me%Worksize%ILB:Me%Worksize%IUB, Me%Worksize%JLB:Me%Worksize%JUB))
                
                PlantDormantInteger (:,:) = 0
                
                call HDF5ReadData   (ObjHDF5, "/Results/"//"IsPlantDormant",                 &
                                     "IsPlantDormant",                                       &
                                     Array2D = PlantDormantInteger,                          &
                                     STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                   &
                    stop 'ReadInitialHDF - ModuleVegetation - ERR04'


                ConvertType = "IntegerToLogical"
                call ConvertLogicalInteger (Me%IsPlantDormant, PlantDormantInteger, ConvertType)
                
                deallocate (PlantDormantInteger)
            
            endif

            call HDF5ReadData   (ObjHDF5, "/Results/"//"HUAccumulated",                  &
                                 "HUAccumulated",                                        &
                                 Array2D = Me%HeatUnits%PlantHUAccumulated,              &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR03'

            !old value = read value
            call SetMatrixValue (Me%HeatUnits%PlantHUAccumulated_Old,                    &
                                 Me%Size2D,                                              &
                                 Me%HeatUnits%PlantHUAccumulated,                        &
                                 Me%ExternalVar%MappingPoints)

            call HDF5ReadData   (ObjHDF5, "/Results/"//"PotentialHUBase",                &
                                 "PotentialHUBase",                                      &
                                 Array2D = Me%HeatUnits%PotentialHUBase,                 &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR04'
            
            !old value = read value
            call SetMatrixValue (Me%HeatUnits%PotentialHUBase_Old,                       &
                                 Me%Size2D,                                              &
                                 Me%HeatUnits%PotentialHUBase,                           &
                                 Me%ExternalVar%MappingPoints)


            call HDF5ReadData   (ObjHDF5, "/Results/"//"TreeCurrentYear",                &
                                 "TreeCurrentYear",                                      &
                                 Array2D = Me%Growth%TreeCurrentYear,                    &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR05'
                            
            Me%Growth%TreeComingFromContinuous (:,:) = .true.                               

            call HDF5ReadData   (ObjHDF5, "/Results/"//"PlantLAIMaxFraction",            &
                                 "PlantLAIMaxFraction",                                  &
                                 Array2D = Me%PlantLAIMaxFraction,                       &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR05.1'

            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

                if (Me%ComputeOptions%Dormancy .or. Me%ComputeOptions%Grazing) then
            
                    if (Me%ComputeOptions%ModelNitrogen) then

                        call HDF5ReadData   (ObjHDF5, "/Results/"//"NitrogenFraction",               &
                                             "NitrogenFraction",                                     &
                                             Array2D = Me%PlantNitrogenFraction,                     &
                                             STAT    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                                   &
                            stop 'ReadInitialHDF - ModuleVegetation - ERR05.2'
                    endif


                    if (Me%ComputeOptions%ModelPhosphorus) then

                        call HDF5ReadData   (ObjHDF5, "/Results/"//"PhosphorusFraction",             &
                                             "PhosphorusFraction",                                   &
                                             Array2D = Me%PlantPhosphorusFraction,                   &
                                             STAT    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                                   &
                            stop 'ReadInitialHDF - ModuleVegetation - ERR05.3'
                    endif
                endif

                if (Me%ComputeOptions%Fertilization) then

                    if (Me%ComputeOptions%ModelNitrogen) then

                        call HDF5ReadData   (ObjHDF5, "/Results/"//"AnnualNitrogenFertilized",       &
                                             "AnnualNitrogenFertilized",                             &
                                             Array2D = Me%AnnualNitrogenFertilized,                  &
                                             STAT    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                                   &
                            stop 'ReadInitialHDF - ModuleVegetation - ERR05.4'
                    endif

                    if (Me%ComputeOptions%ModelPhosphorus) then

                        call HDF5ReadData   (ObjHDF5, "/Results/"//"AnnualPhosphorusFertilized",     &
                                             "AnnualPhosphorusFertilized",                           &
                                             Array2D = Me%AnnualPhosphorusFertilized,                &
                                             STAT    = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                                   &
                            stop 'ReadInitialHDF - ModuleVegetation - ERR05.5'
                    endif
                endif

            endif

            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialHDF - ModuleVegetation - ERR06'

        else
            
            write(*,*)
            stop 'ReadInitialHDF - ModuleVegetation - ERR07'

        end if cd0

    end subroutine ReadInitialHDF

    !--------------------------------------------------------------------------           

    subroutine CheckPlantGrowing

        !Local-----------------------------------------------------------------
        integer                                     :: i,j, PlantType
        !Begin-----------------------------------------------------------------
        
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            
            if (Me%ExternalVar%MappingPoints(i, j) == VegetationPoint) then
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle                
                
                Me%Growth%TreeComingFromContinuous(i,j) = .false.
                PlantType  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                if (PlantType == Trees) then
                    Me%Growth%TreeComingFromContinuous(i,j) = .true.
                endif                
                
                if(Me%StateVariables%TotalPlantBiomass(i,j) .gt. 0.0) then
                    
                    Me%IsPlantGrowing(i,j) = .true.
                    

                else
                    Me%IsPlantGrowing(i,j) = .false.
                endif
            
            endif

        enddo
        enddo


    end subroutine CheckPlantGrowing

    !--------------------------------------------------------------------------           


    subroutine ConvertLogicalInteger(LogicalMatrix, IntegerMatrix, ConversionType)

        !Local-----------------------------------------------------------------
        logical, dimension(:,:), pointer             :: LogicalMatrix
        integer, dimension(:,:), pointer             :: IntegerMatrix
        character (Len = StringLength)               :: ConversionType      
        integer                                     :: i,j
        !Begin-----------------------------------------------------------------
        
        
        if (ConversionType == "LogicalToInteger") then
            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExternalVar%BasinPoints(i, j) == 1) then
                    
                    if(LogicalMatrix(i,j)) then
                        IntegerMatrix(i,j) = 1
                    else
                        IntegerMatrix(i,j) = 0
                    endif

               endif

            enddo
            enddo
        
        elseif (ConversionType == "IntegerToLogical") then

            do j = Me%Size%JLB, Me%Size%JUB
            do i = Me%Size%ILB, Me%Size%IUB
                
                if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                    
                    if(IntegerMatrix(i,j)==1) then
                        LogicalMatrix(i,j) = .true.
                    else
                        LogicalMatrix(i,j) = .false.
                    endif
                
                endif

            enddo
            enddo
            
        else
            write(*,*)'Error in converting logical and integer'
            stop 'ConvertLogicalInteger - Module Vegetation - ERR01'
        
        endif

    end subroutine ConvertLogicalInteger
   
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetVegetationDT(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real                                        :: Scalar
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Scalar  = Me%ComputeOptions%VegetationDT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetVegetationDT
    
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

            Scalar  => Me%StateVariables%LeafAreaIndex

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetLeafAreaIndex
    
    !--------------------------------------------------------------------------

    subroutine GetPotLeafAreaIndex(VegetationID, Scalar, STAT)
                                  
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

            Scalar  => Me%StateVariables%PotLeafAreaIndex

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetPotLeafAreaIndex
    
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

            Scalar  => Me%StateVariables%SpecificLeafStorage

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

            Scalar  => Me%StateVariables%EVTPCropCoefficient

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

            Scalar  => Me%StateVariables%RootDepth

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetRootDepth
    
    !--------------------------------------------------------------------------
    
    !For PMP to compute the sources/sinks with the same rrot depth that gave origin to the fluxes
    subroutine GetRootDepthOld(VegetationID, Scalar, STAT)
                                  
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

            Scalar  => Me%RootDepthOld

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetRootDepthOld
    
    !--------------------------------------------------------------------------
    
        subroutine GetCanopyHeight(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer               :: Scalar
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

            Scalar  => Me%StateVariables%CanopyHeight

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetCanopyHeight
    
    !--------------------------------------------------------------------------
    
    subroutine GetNutrientFraction(VegetationID, NitrogenFraction, PhosphorusFraction, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        real, dimension(:,:), pointer, optional     :: NitrogenFraction
        real, dimension(:,:), pointer, optional     :: PhosphorusFraction
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(NitrogenFraction)) then
                
                call Read_Lock(mVEGETATION_, Me%InstanceID)

                NitrogenFraction  => Me%PlantNitrogenFraction
            
            endif

            if (present(PhosphorusFraction)) then
                
                call Read_Lock(mVEGETATION_, Me%InstanceID)

                PhosphorusFraction  => Me%PlantPhosphorusFraction
            
            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_

    end subroutine GetNutrientFraction
    
    !--------------------------------------------------------------------------
    
    subroutine GetCanopyStorageType(VegetationID, IsConstant, ConstantValue, STAT)

        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        logical, intent(OUT)                        :: IsConstant
        real, optional, intent(OUT)                 :: ConstantValue
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_, STAT_CALL
        type(T_Property), pointer                   :: SpecificLeafStorage        
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mVEGETATION_, Me%InstanceID)

            call SearchProperty(SpecificLeafStorage, SpecificLeafStorage_, .false., STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'GetCanopyStorageType - ModuleVegetation - ERR010'

            IsConstant = SpecificLeafStorage%IsConstant
        
            if (present(ConstantValue)) ConstantValue = SpecificLeafStorage%ConstantValue

            call Read_UnLock(mVegetation_, Me%InstanceID, "UngetCanopyStorageType")
            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if 

        if (present(STAT)) STAT = STAT_
            
    end subroutine GetCanopyStorageType
    
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

    subroutine GetTranspirationBottomLayer(VegetationID, Scalar, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        integer, dimension(:,:), pointer, optional  :: Scalar
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

            Scalar  => Me%TranspirationBottomLayer

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine GetTranspirationBottomLayer

    !--------------------------------------------------------------------------

    subroutine GetVegetationGrowing(VegetationID, Growing, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        logical, dimension(:,:), pointer            :: Growing
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
            
            if (.not. Me%ComputeOptions%Evolution%ModelSWAT) then
                allocate(Me%IsPlantGrowing (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                Me%IsPlantGrowing = .false.
                Growing  => Me%IsPlantGrowing
                !deallocate (Me%IsPlantGrowing)
            endif
            
            Growing  => Me%IsPlantGrowing

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine GetVegetationGrowing

    !--------------------------------------------------------------------------

    subroutine GetVegetationOptions(VegetationID,                              &
                                    ModelWater,                                &
                                    ModelNitrogen,                             &
                                    ModelPhosphorus,                           &
                                    NutrientFluxesWithSoil,                    &
                                    Grazing,                                   &
                                    HarvestKill,                               &
                                    Dormancy,                                  &
                                    Fertilization,                             &
                                    GrowthModel,                               &
                                    ModelCanopyHeight,                         &
                                    NutrientUptakeMethod,                      &
                                    Pesticide,                                 &
                                    NumberOfPesticides,                        &
                                    VegParticFertilization,                    &
                                    STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                           :: VegetationID
        logical, optional                 :: ModelWater
        logical, optional                 :: ModelNitrogen
        logical, optional                 :: ModelPhosphorus
        logical, optional                 :: NutrientFluxesWithSoil
        logical, optional                 :: Grazing
        logical, optional                 :: HarvestKill
        logical, optional                 :: Dormancy
        logical, optional                 :: Fertilization
        logical, optional                 :: GrowthModel
        logical, optional                 :: ModelCanopyHeight
        integer, optional                 :: NutrientUptakeMethod
        logical, optional                 :: Pesticide
        integer, optional                 :: NumberOfPesticides
        logical, optional                 :: VegParticFertilization

        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then

!            call Read_Lock(mVEGETATION_, Me%InstanceID)
            if(present(ModelWater            )) ModelWater             = Me%ComputeOptions%ModelWater
            if(present(ModelNitrogen         )) ModelNitrogen          = Me%ComputeOptions%ModelNitrogen
            if(present(ModelPhosphorus       )) ModelPhosphorus        = Me%ComputeOptions%ModelPhosphorus
            if(present(NutrientFluxesWithSoil)) NutrientFluxesWithSoil = Me%ComputeOptions%NutrientFluxesWithSoil
            if(present(Grazing               )) Grazing                = Me%ComputeOptions%Grazing
            if(present(HarvestKill           )) HarvestKill            = Me%ComputeOptions%HarvestKill
            if(present(Dormancy              )) Dormancy               = Me%ComputeOptions%Dormancy
            if(present(Fertilization         )) Fertilization          = Me%ComputeOptions%Fertilization
            if(present(GrowthModel           )) GrowthModel            = Me%ComputeOptions%Evolution%GrowthModelNeeded
            if(present(ModelCanopyHeight     )) ModelCanopyHeight      = Me%ComputeOptions%ModelCanopyHeight
            if(present(NutrientUptakeMethod  )) NutrientUptakeMethod   = Me%ComputeOptions%NutrientUptakeMethod
            if(present(Pesticide             )) Pesticide              = Me%ComputeOptions%Pesticide
            if(present(NumberOfPesticides    )) NumberOfPesticides     = Me%Fluxes%Pesticides%UniquePesticides
            if(present(VegParticFertilization)) VegParticFertilization = Me%ComputeOptions%VegParticFertilization
            

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine GetVegetationOptions

    !Fluxes to soil and from soil
    !--------------------------------------------------------------------------

    subroutine GetVegetationSoilFluxes(VegetationID,                             &
                                       SoilFluxesActive,                         &
                                       GrazingBiomass,                           &
                                       GrazingNitrogen,                          &
                                       GrazingPhosphorus,                        &
                                       HarvestKillAerialBiomass,                 &
                                       HarvestKillNitrogen,                      &
                                       HarvestKillPhosphorus,                    &
                                       HarvestKillRootBiomass,                   &
                                       DormancyBiomass,                          &
                                       DormancyNitrogen,                         &
                                       DormancyPhosphorus,                       &
                                       FertilNitrateSurface,                     &
                                       FertilNitrateSubSurface,                  &
                                       FertilAmmoniaSurface,                     &
                                       FertilAmmoniaSubSurface,                  &
                                       FertilOrganicNSurface,                    &
                                       FertilOrganicNSubSurface,                 &
                                       FertilOrganicNParticFluff,                &
                                       FertilOrganicPSurface,                    &
                                       FertilOrganicPSubSurface,                 &
                                       FertilOrganicPParticFluff,                &
                                       FertilMineralPSurface,                    &
                                       FertilMineralPSubSurface,                 &
                                       Pest,                                     &
                                       PesticideIDNumber,                        &
                                       PesticideSoil,                            &
                                       NitrogenUptake,                           &
                                       PhosphorusUptake,                         &
                                       STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        logical, dimension(:,:), pointer, optional  :: SoilFluxesActive
        real, dimension(:,:), pointer, optional     :: GrazingBiomass
        real, dimension(:,:), pointer, optional     :: GrazingNitrogen
        real, dimension(:,:), pointer, optional     :: GrazingPhosphorus
        real, dimension(:,:), pointer, optional     :: HarvestKillAerialBiomass
        real, dimension(:,:), pointer, optional     :: HarvestKillNitrogen
        real, dimension(:,:), pointer, optional     :: HarvestKillPhosphorus
        real, dimension(:,:), pointer, optional     :: HarvestKillRootBiomass
        real, dimension(:,:), pointer, optional     :: DormancyBiomass
        real, dimension(:,:), pointer, optional     :: DormancyNitrogen
        real, dimension(:,:), pointer, optional     :: DormancyPhosphorus
        real, dimension(:,:), pointer, optional     :: FertilNitrateSurface
        real, dimension(:,:), pointer, optional     :: FertilNitrateSubSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNParticFluff
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPParticFluff
        real, dimension(:,:), pointer, optional     :: FertilMineralPSurface
        real, dimension(:,:), pointer, optional     :: FertilMineralPSubSurface
        integer, optional                           :: Pest
        integer, optional                           :: PesticideIDNumber
        real, dimension(:,:), pointer, optional     :: PesticideSoil
        real, dimension(:,:,:),pointer,optional     :: NitrogenUptake
        real, dimension(:,:,:),pointer,optional     :: PhosphorusUptake

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
            
            if (present(SoilFluxesActive         )) SoilFluxesActive         => Me%SoilFluxesActive
            if (present(GrazingBiomass           )) GrazingBiomass           => Me%Fluxes%ToSoil%GrazingBiomassToSoil
            if (present(GrazingNitrogen          )) GrazingNitrogen          => Me%Fluxes%ToSoil%GrazingNitrogenToSoil
            if (present(GrazingPhosphorus        )) GrazingPhosphorus        => Me%Fluxes%ToSoil%GrazingPhosphorusToSoil
            if (present(HarvestKillAerialBiomass )) HarvestKillAerialBiomass => Me%Fluxes%ToSoil%HarvestKillBiomassToSoil 
            if (present(HarvestKillNitrogen      )) HarvestKillNitrogen      => Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil
            if (present(HarvestKillPhosphorus    )) HarvestKillPhosphorus    => Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil
            if (present(HarvestKillRootBiomass   )) HarvestKillRootBiomass   => Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil
            if (present(DormancyBiomass          )) DormancyBiomass          => Me%Fluxes%ToSoil%DormancyBiomassToSoil
            if (present(DormancyNitrogen         )) DormancyNitrogen         => Me%Fluxes%ToSoil%DormancyNitrogenToSoil
            if (present(DormancyPhosphorus       )) DormancyPhosphorus       => Me%Fluxes%ToSoil%DormancyPhosphorusToSoil
            if (present(FertilNitrateSurface     )) FertilNitrateSurface     => Me%Fluxes%ToSoil%FertilNitrateToSoilSurface
            if (present(FertilNitrateSubSurface  )) FertilNitrateSubSurface  => Me%Fluxes%ToSoil%FertilNitrateToSoilSubSurface
            if (present(FertilAmmoniaSurface     )) FertilAmmoniaSurface     => Me%Fluxes%ToSoil%FertilAmmoniaToSoilSurface
            if (present(FertilAmmoniaSubSurface  )) FertilAmmoniaSubSurface  => Me%Fluxes%ToSoil%FertilAmmoniaToSoilSubSurface
            if (present(FertilOrganicNSurface    )) FertilOrganicNSurface    => Me%Fluxes%ToSoil%FertilOrganicNToSoilSurface
            if (present(FertilOrganicNSubSurface )) FertilOrganicNSubSurface => Me%Fluxes%ToSoil%FertilOrganicNToSoilSubSurface
            if (present(FertilOrganicPSurface    )) FertilOrganicPSurface    => Me%Fluxes%ToSoil%FertilOrganicPToSoilSurface
            if (present(FertilOrganicPSubSurface )) FertilOrganicPSubSurface => Me%Fluxes%ToSoil%FertilOrganicPToSoilSubSurface
            if (present(FertilMineralPSurface    )) FertilMineralPSurface    => Me%Fluxes%ToSoil%FertilMineralPToSoilSurface
            if (present(FertilMineralPSubSurface )) FertilMineralPSubSurface => Me%Fluxes%ToSoil%FertilMineralPToSoilSubSurface
            if (present(NitrogenUptake           )) NitrogenUptake           => Me%Fluxes%FromSoil%NitrogenUptakeFromSoil
            if (present(PhosphorusUptake         )) PhosphorusUptake         => Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil
            if (present(FertilOrganicNParticFluff)) FertilOrganicNParticFluff=> Me%Fluxes%ToSoil%FertilOrganicNParticToFluff
            if (present(FertilOrganicPParticFluff)) FertilOrganicPParticFluff=> Me%Fluxes%ToSoil%FertilOrganicPParticToFluff
            
            if (present(PesticideSoil            )) then
                if( .not. present(Pest)) stop 'GetVegetationSoilFluxes - ModuleVegetation - ERR10'
                PesticideIDNumber        =  Me%Fluxes%Pesticides%Application(Pest)%ID%IDNumber
                PesticideSoil            => Me%Fluxes%Pesticides%Application(Pest)%Soil
            endif
            
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine GetVegetationSoilFluxes

    !--------------------------------------------------------------------------

    subroutine SetSoilConcVegetation   (ObjVegetationID,                        & 
                                        Nitrate,                                &
                                        InorganicPhosphorus,                    &
                                        STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjVegetationID
        real, dimension(:,:,:), pointer, optional       :: Nitrate
        real, dimension(:,:,:), pointer, optional       :: InorganicPhosphorus
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjVegetationID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            if (present(Nitrate)) Me%ExternalVar%SoilNitrate  => Nitrate

            if (present(InorganicPhosphorus)) Me%ExternalVar%SoilPhosphorus  => InorganicPhosphorus


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetSoilConcVegetation 

    !---------------------------------------------------------------------------

    subroutine SetECw (ObjVegetationID, ECw, STAT)

        !Arguments--------------------------------------------------------------
        integer                         :: ObjVegetationID
        real, dimension(:,:,:), pointer :: ECw
        integer, intent(OUT), optional  :: STAT

        !Local------------------------------------------------------------------
        integer                         :: STAT_, ready_ !, index
        
        !-----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(ObjVegetationID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            Me%ExternalVar%ECw => ECw            
            STAT_ = SUCCESS_

        else
        
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_
        !-----------------------------------------------------------------------
    
    end subroutine SetECw
    
    !---------------------------------------------------------------------------

    !Fluxes to aerial part of vegetation - pesticide in leafs    
    subroutine GetVegetationAerialFluxes(VegetationID,                           &
                                       Pest,                                     &
                                       PesticideIDNumber,                        &
                                       PesticideVegetation,                      &
                                       STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        integer, optional                           :: Pest
        integer, optional                           :: PesticideIDNumber
        real, dimension(:,:), pointer, optional     :: PesticideVegetation

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
            
            if (present(PesticideVegetation      )) then
                if( .not. present(Pest)) stop 'GetVegetationAerialFluxes - ModuleVegetation - ERR10'
                PesticideIDNumber        =  Me%Fluxes%Pesticides%Application(Pest)%ID%IDNumber
                PesticideVegetation      => Me%Fluxes%Pesticides%Application(Pest)%Vegetation
            endif
            
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine GetVegetationAerialFluxes

    !--------------------------------------------------------------------------

    subroutine UnGetVegetationAerialFluxes(VegetationID,                         &
                                           Pest,                                 &
                                           PesticideVegetation,                  &
                                           STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        integer, optional                           :: Pest
        real, dimension(:,:), pointer, optional     :: PesticideVegetation

        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

        
        if ((ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PesticideVegetation      )) then
                if( .not. present(Pest)) stop 'GetVegetationAerialFluxes - ModuleVegetation - ERR10'
                
                nullify(PesticideVegetation)
                
                call Read_UnLock(mVEGETATION_, Me%InstanceID, "UngetVegetationAerialFluxes")
            endif
           
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if 

        if (present(STAT)) STAT = STAT_


    end subroutine UnGetVegetationAerialFluxes

    !--------------------------------------------------------------------------

    subroutine UngetVegetationSoilFluxes(VegetationID,                             &
                                         SoilFluxesActive,                         &
                                         GrazingBiomass,                           &
                                         GrazingNitrogen,                          &
                                         GrazingPhosphorus,                        &
                                         HarvestKillAerialBiomass,                 &
                                         HarvestKillNitrogen,                      &
                                         HarvestKillPhosphorus,                    &
                                         HarvestKillRootBiomass,                   &
                                         DormancyBiomass,                          &
                                         DormancyNitrogen,                         &
                                         DormancyPhosphorus,                       &
                                         FertilNitrateSurface,                     &
                                         FertilNitrateSubSurface,                  &
                                         FertilAmmoniaSurface,                     &
                                         FertilAmmoniaSubSurface,                  &
                                         FertilOrganicNSurface,                    &
                                         FertilOrganicNSubSurface,                 &
                                         FertilOrganicNParticFluff,                &
                                         FertilOrganicPSurface,                    &
                                         FertilOrganicPSubSurface,                 &
                                         FertilOrganicPParticFluff,                &
                                         FertilMineralPSurface,                    &
                                         FertilMineralPSubSurface,                 &
                                         NitrogenUptake,                           &
                                         PhosphorusUptake,                         &
                                         PesticideSoil,                            &
                                         STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: VegetationID
        logical, dimension(:,:), pointer, optional  :: SoilFluxesActive
        real, dimension(:,:), pointer, optional     :: GrazingBiomass
        real, dimension(:,:), pointer, optional     :: GrazingNitrogen
        real, dimension(:,:), pointer, optional     :: GrazingPhosphorus
        real, dimension(:,:), pointer, optional     :: HarvestKillAerialBiomass
        real, dimension(:,:), pointer, optional     :: HarvestKillNitrogen
        real, dimension(:,:), pointer, optional     :: HarvestKillPhosphorus
        real, dimension(:,:), pointer, optional     :: HarvestKillRootBiomass
        real, dimension(:,:), pointer, optional     :: DormancyBiomass
        real, dimension(:,:), pointer, optional     :: DormancyNitrogen
        real, dimension(:,:), pointer, optional     :: DormancyPhosphorus
        real, dimension(:,:), pointer, optional     :: FertilNitrateSurface
        real, dimension(:,:), pointer, optional     :: FertilNitrateSubSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNParticFluff
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPParticFluff
        real, dimension(:,:), pointer, optional     :: FertilMineralPSurface
        real, dimension(:,:), pointer, optional     :: FertilMineralPSubSurface
        real, dimension(:,:,:),pointer,optional     :: NitrogenUptake
        real, dimension(:,:,:),pointer,optional     :: PhosphorusUptake
        real, dimension(:,:), pointer, optional     :: PesticideSoil      

        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(VegetationID, ready_) 

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
 
            if (present(SoilFluxesActive        )) nullify(SoilFluxesActive)
            if (present(GrazingBiomass          )) nullify(GrazingBiomass)
            if (present(GrazingNitrogen         )) nullify(GrazingNitrogen)
            if (present(GrazingPhosphorus       )) nullify(GrazingPhosphorus)
            if (present(HarvestKillAerialBiomass)) nullify(HarvestKillAerialBiomass)
            if (present(HarvestKillNitrogen     )) nullify(HarvestKillNitrogen)
            if (present(HarvestKillPhosphorus   )) nullify(HarvestKillPhosphorus)
            if (present(HarvestKillRootBiomass  )) nullify(HarvestKillRootBiomass)
            if (present(DormancyBiomass         )) nullify(DormancyBiomass)
            if (present(DormancyNitrogen        )) nullify(DormancyNitrogen)
            if (present(DormancyPhosphorus      )) nullify(DormancyPhosphorus)
            if (present(FertilNitrateSurface    )) nullify(FertilNitrateSurface)
            if (present(FertilNitrateSubSurface )) nullify(FertilNitrateSubSurface)
            if (present(FertilAmmoniaSurface    )) nullify(FertilAmmoniaSurface)
            if (present(FertilAmmoniaSubSurface )) nullify(FertilAmmoniaSubSurface)
            if (present(FertilOrganicNSurface   )) nullify(FertilOrganicNSurface)
            if (present(FertilOrganicNSubSurface)) nullify(FertilOrganicNSubSurface)
            if (present(FertilOrganicPSurface   )) nullify(FertilOrganicPSurface)
            if (present(FertilOrganicPSubSurface)) nullify(FertilOrganicPSubSurface)
            if (present(FertilMineralPSurface   )) nullify(FertilMineralPSurface)
            if (present(FertilMineralPSubSurface)) nullify(FertilMineralPSubSurface)
            if (present(NitrogenUptake          )) nullify(NitrogenUptake)
            if (present(PhosphorusUptake        )) nullify(PhosphorusUptake)
            if (present(PesticideSoil           )) nullify(PesticideSoil)
            if (present(FertilOrganicNParticFluff)) nullify(FertilOrganicNParticFluff)
            if (present(FertilOrganicPParticFluff)) nullify(FertilOrganicPParticFluff)

            call Read_UnLock(mVegetation_, Me%InstanceID, "UngetVegetationSoilFluxes")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine UngetVegetationSoilFluxes    

   !----------------------------------------------------------------------


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


    end subroutine UngetVegetation2D    

   !----------------------------------------------------------------------

    subroutine UngetVegetation2D_I(VegetationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: VegetationID
        integer, pointer, dimension(:,:) :: Array
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

            call Read_UnLock(mVegetation_, Me%InstanceID, "UngetVegetation2D_I")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine UngetVegetation2D_I    

   !----------------------------------------------------------------------
   
    subroutine UngetVegetation2D_L(VegetationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: VegetationID
        logical, pointer, dimension(:,:) :: Array
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

            call Read_UnLock(mVegetation_, Me%InstanceID, "UngetVegetation2D_I")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine UngetVegetation2D_L    

   !----------------------------------------------------------------------   

    subroutine UngetVegetation3D(VegetationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: VegetationID
        real, pointer, dimension(:,:,:) :: Array
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

            call Read_UnLock(mVegetation_, Me%InstanceID, "UngetVegetation3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine UngetVegetation3D    

   !----------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyVegetation(ObjVegetationID,                  &
                                MappingPoints,                    &
                                PotentialTranspiration,           &
                                ActualTranspiration,              &
                                ECw,                              &
                                STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                            :: ObjVegetationID         !IN
        integer, dimension(:,:  ), pointer             :: MappingPoints           !IN
        real,    dimension(:,:  ), pointer             :: PotentialTranspiration  !IN
        real,    dimension(:,:,:), pointer             :: ActualTranspiration     !OUT
        real,    dimension(:,:,:), pointer, optional   :: ECw                     !IN
        integer, intent(OUT),      optional            :: STAT                    !OUT

        !Local-----------------------------------------------------------------
        integer                                        :: STAT_, ready_
        integer                                        :: STAT_CALL
        integer                                        :: JulDay
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjVegetationID, ready_)

        if (present(ECw)) then
            Me%ECw => ECw
        else
            Me%ECw => null()
        endif
        
        if (ready_ .EQ. IDLE_ERR_) then

            !Actualize the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyVegetation - ModuleVegetation - ERR01'

            !Sets External Variable
            Me%ExternalVar%MappingPoints => MappingPoints                        
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
                
                else 
                    
                    !Do not use vegetation growth model. Only water uptake is modeled for now
                    !Use readed vegetation properties (LAI, root depth) to uptake water
                    if (Me%ComputeOptions%ModelWater) then
                        call WaterUptake
                    endif
                    
                    !Uptake nitrogen - soil nitrogen concentration with water uptake volume
                    if (Me%ComputeOptions%ModelNitrogen) then
                        call NitrogenUptake_TranspConc
                    endif
                    
                    !!Uptake phosphorus - soil phosphorus concentration with water uptake volume
                    if (Me%ComputeOptions%ModelPhosphorus) then
                        call PhosphorusUptake_TranspConc
                    endif
                    
                endif
            

                !Output
                if (Me%OutPut%HDF_ON)           call Modify_OutputHDF                        
                if (Me%OutPut%TimeSerie_ON)     call Modify_OutPutTimeSeries


                !Get julian day from iteration to keep to next
                call JulianDay(Me%ExternalVar%Now, JulDay)
                Me%ExternalVar%JulianDay_Old = JulDay

                nullify (Me%ExternalVar%MappingPoints)                
                
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
                
        !Begin-----------------------------------------------------------------

        !If dt's the same use values without integration
        if (Me%ComputeOptions%VegetationDT .eq. Me%ComputeOptions%IntegrationDT) then

!            Me%ExternalVar%Integration%AveragePotTPDuringDT => Me%ExternalVar%PotentialTranspiration

            call SetMatrixValue (Me%ExternalVar%Integration%AveragePotTPDuringDT,           &
                                 Me%Size2D,                                                 &
                                 Me%ExternalVar%PotentialTranspiration,                     &
                                 Me%ExternalVar%MappingPoints)

            if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

                call ReadLockExternalVar(ReadAtmosphere = .true.)

!                Me%ExternalVar%Integration%AverageAirTemPDuringDT     => Me%ExternalVar%AirTemperature
                call SetMatrixValue (Me%ExternalVar%Integration%AverageAirTemPDuringDT,      &
                                     Me%Size2D,                                              &
                                     Me%ExternalVar%AirTemperature,                          &
                                     Me%ExternalVar%MappingPoints)
    
!                Me%ExternalVar%Integration%AverageAirHumidityDuringDT => Me%ExternalVar%RelativeHumidity
                call SetMatrixValue (Me%ExternalVar%Integration%AverageAirHumidityDuringDT,  &
                                     Me%Size2D,                                              &
                                     Me%ExternalVar%RelativeHumidity,                        &
                                     Me%ExternalVar%MappingPoints)

!                Me%ExternalVar%Integration%AverageRadiationDuringDT   => Me%ExternalVar%SolarRadiation
                call SetMatrixValue (Me%ExternalVar%Integration%AverageRadiationDuringDT,    &
                                     Me%Size2D,                                              &
                                     Me%ExternalVar%SolarRadiation,                          &
                                     Me%ExternalVar%MappingPoints)    
                
                call ReadUnLockExternalVar(ReadAtmosphere = .true.)

            endif
        

        else
                    
            call ReadLockExternalVar(ReadAtmosphere = .true.)

            Me%nIterations = Me%nIterations + 1
      
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then  
                
                    !Plant type 0 is NO PLANT, so no operations should happen
                    if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

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
                            Me%ExternalVar%Integration%AverageAirTemPDuringDT(i,j)      = &
                                  Me%ExternalVar%Integration%SumTemperature(i,j) / Me%nIterations
                            Me%ExternalVar%Integration%AverageAirHumidityDuringDT(i,j)  = &
                                  Me%ExternalVar%Integration%SumHumidity(i,j)   / Me%nIterations
                            Me%ExternalVar%Integration%AverageRadiationDuringDT(i,j)    = &
                                  Me%ExternalVar%Integration%SumRadiation(i,j)   / Me%nIterations
                            Me%ExternalVar%Integration%SumTemperature(i,j) = 0.0
                            Me%ExternalVar%Integration%SumHumidity   (i,j) = 0.0
                            Me%ExternalVar%Integration%SumRadiation  (i,j) = 0.0 
                        endif
                    endif

                    ! Transpiration always needed
                    Me%ExternalVar%Integration%SumPotTP(i,j) = Me%ExternalVar%Integration%SumPotTP(i,j)                          &
                                                               + Me%ExternalVar%PotentialTranspiration(i,j)
                
                    if (Me%ExternalVar%Now .ge. Me%NextCompute) then
                        Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j)  = Me%ExternalVar%Integration%SumPotTP(i,j)        &
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
        integer                                              :: Pest
        
        !Begin-----------------------------------------------------------------
        if (Me%ComputeOptions%Grazing) then

            Me%Fluxes%BiomassGrazed                                       (:,:) = 0.0 
            Me%Fluxes%BiomassGrazedFraction                               (:,:) = 0.0 
            Me%Fluxes%ToSoil%GrazingBiomassToSoil                         (:,:) = 0.0 
            if (Me%ComputeOptions%ModelNitrogen) then
                Me%Fluxes%NitrogenGrazed                                  (:,:) = 0.0
                Me%Fluxes%ToSoil%GrazingNitrogenToSoil                    (:,:) = 0.0 
            endif                                                   
            if (Me%ComputeOptions%ModelPhosphorus) then             
                Me%Fluxes%PhosphorusGrazed                                (:,:) = 0.0 
                Me%Fluxes%ToSoil%GrazingPhosphorusToSoil                  (:,:) = 0.0 
            endif                                                   

        endif
        
        if (Me%ComputeOptions%HarvestKill) then

            Me%Fluxes%BiomassRemovedInHarvest                             (:,:) = 0.0  
            Me%Fluxes%BiomassHarvestedFraction                            (:,:) = 0.0
            Me%Fluxes%ToSoil%HarvestKillBiomassToSoil                     (:,:) = 0.0 
            Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil                    (:,:) = 0.0
            if (Me%ComputeOptions%ModelNitrogen) then
                Me%Fluxes%NitrogenRemovedInHarvest                        (:,:) = 0.0 
                Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil                (:,:) = 0.0 
            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
                Me%Fluxes%PhosphorusRemovedInHarvest                      (:,:) = 0.0
                Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil              (:,:) = 0.0 
            endif
        endif

        if (Me%ComputeOptions%Dormancy) then
            
            Me%Fluxes%BiomassRemovedInDormancy                          (:,:) = 0.0
!            Me%Fluxes%ToSoil%DormancyBiomassToSoil                      (:,:) = 0.0
            if (Me%ComputeOptions%ModelNitrogen) then
                Me%Fluxes%NitrogenRemovedInDormancy                     (:,:) = 0.0
!                Me%Fluxes%ToSoil%DormancyNitrogenToSoil                 (:,:) = 0.0
            endif                                                
            if (Me%ComputeOptions%ModelPhosphorus) then          
                Me%Fluxes%PhosphorusRemovedInDormancy                   (:,:) = 0.0
!                Me%Fluxes%ToSoil%DormancyPhosphorusToSoil               (:,:) = 0.0
            endif                                                
        endif
        if (Me%ComputeOptions%Fertilization) then
            if (Me%ComputeOptions%ModelNitrogen) then
                Me%Fluxes%FertilNitrateInSurface                 (:,:) = 0.0
                Me%Fluxes%FertilNitrateInSubSurface              (:,:) = 0.0
                Me%Fluxes%FertilAmmoniaInSurface                 (:,:) = 0.0
                Me%Fluxes%FertilAmmoniaInSubSurface              (:,:) = 0.0
                Me%Fluxes%FertilOrganicNInSurface                (:,:) = 0.0
                Me%Fluxes%FertilOrganicNInSubSurface             (:,:) = 0.0
                Me%Fluxes%FertilOrganicNParticInFluff            (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilNitrateToSoilSurface      (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilNitrateToSoilSubSurface   (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilAmmoniaToSoilSurface      (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilAmmoniaToSoilSubSurface   (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilOrganicNToSoilSurface     (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilOrganicNToSoilSubSurface  (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilOrganicNParticToFluff     (:,:) = 0.0

            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
                Me%Fluxes%FertilOrganicPInSurface                (:,:) = 0.0
                Me%Fluxes%FertilOrganicPInSubSurface             (:,:) = 0.0
                Me%Fluxes%FertilMineralPInSurface                (:,:) = 0.0
                Me%Fluxes%FertilMineralPInSubSurface             (:,:) = 0.0
                Me%Fluxes%FertilOrganicPParticInFluff            (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilOrganicPToSoilSurface     (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilOrganicPToSoilSubSurface  (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilMineralPToSoilSurface     (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilMineralPToSoilSubSurface  (:,:) = 0.0
                Me%Fluxes%ToSoil%FertilOrganicPParticToFluff     (:,:) = 0.0

            endif

        endif
        
        if (Me%ComputeOptions%Pesticide) then
            do Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
                Me%Fluxes%Pesticides%Application(Pest)%Soil      (:,:) = 0.0
                Me%Fluxes%Pesticides%Application(Pest)%Vegetation(:,:) = 0.0
            enddo
        endif

        Me%Fluxes%WaterUptake           (:,:  ) = 0.0 
        Me%Fluxes%WaterUptakeLayer      (:,:,:) = 0.0 
        Me%Fluxes%LAIChange             (:,:  ) = 0.0

        if (Me%ComputeOptions%ModelNitrogen) then
            Me%Fluxes%NitrogenUptake        (:,:  ) = 0.0 
            Me%Fluxes%NitrogenUptakeLayer   (:,:,:) = 0.0 
            Me%OptimalTotalPlantNitrogen    (:,:) = 0.0 
        endif
        if (Me%ComputeOptions%ModelPhosphorus) then
            Me%Fluxes%PhosphorusUptake      (:,:  ) = 0.0 
            Me%Fluxes%PhosphorusUptakeLayer (:,:,:) = 0.0 
            Me%OptimalTotalPlantPhosphorus  (:,:) = 0.0 
        endif
        if (Me%ComputeOptions%ModelPlantBiomass) then
            Me%Fluxes%BiomassGrowth (:,:) = 0.0
            Me%Growth%PAR           (:,:) = 0.0
            Me%Growth%RUE           (:,:) = 0.0
            Me%Growth%PotentialGrowth (:,:) = 0.0
            Me%Growth%PotentialBiomass (:,:) = 0.0
        endif

       
    end subroutine NullifyFluxes

    !--------------------------------------------------------------------------

    subroutine CheckPlantState
 
       !Local-----------------------------------------------------------------
        integer                                            :: i,j, PlantType
        logical                                            :: Dormant
        integer                                            :: JulDay, JulDay_Old
        integer                                            :: Pest, FertApp
        
       !Begin-----------------------------------------------------------------
        
        !If dormancy and day changed, compute day lenght
        if (Me%ComputeOptions%Dormancy) then
            call JulianDay(Me%ExternalVar%Now, JulDay)
            JulDay_Old = Me%ExternalVar%JulianDay_Old
        
            if(JulDay .ne. JulDay_Old) then
                call ComputeDayLength
            endif
        endif
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                PlantType  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType                
                
                !Reset global state variables
                if (Me%ComputeOptions%HarvestKill) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        Me%IsPlantGrowing(i,j)      = .false.
                        Me%KillOccurred(i,j)        = .false.        
                        Me%HarvestKillOccurred(i,j) = .false.
                        if (Me%ComputeOptions%Dormancy) then
                            Me%IsPlantDormant(i,j)      = .false.
                        endif
                        
                        if (PlantType == Trees) then                        
                            Me%Growth%TreeCurrentYear(i,j) = 0
                        endif
                    endif
                    if (Me%HarvestOnlyOccurred(i,j)) then
                        Me%HarvestOnlyOccurred(i,j) = .false.
                    endif                
                endif

                
                if (Me%ComputeOptions%Fertilization) then
                    do FertApp = 1, Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%NumberFertilizerApps
                        if (Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppOccurred) then
                            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppOccurred = .false.
                        endif
                    enddo
                endif
                
                if (Me%PlantingOccurred(i,j)) then
                    Me%PlantingOccurred(i,j) = .false.
                endif
                
                if (Me%ComputeOptions%Dormancy) then
                    if (Me%PlantGoingDormant(i,j)) then
                        Me%PlantGoingDormant(i,j) = .false.
                    endif
                endif
                
                if (Me%ComputeOptions%Pesticide) then
                    do Pest = 1, Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%NumberPesticideApps
                        if (Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(Pest)%PesticideAppOccurred) then
                            Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(Pest)%PesticideAppOccurred = .false.
                        endif
                    enddo                
!                    do Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
!                        if (Me%Fluxes%Pesticides%Application(Pest)%PesticideAppOccurred(i,j)) then
!                            Me%Fluxes%Pesticides%Application(Pest)%PesticideAppOccurred(i,j) = .false.
!                        endif
!                    enddo
                endif
                
            endif
        enddo
        enddo
        
        !If agricultural practices defined in timeserie or HDF than update vegetationID grid
        !Need to be after kill variable reset and before any other plant state checks
        if (Me%AgricPract%SolutionFromFile) then
            call CheckVegetationIDChanges        
        endif        

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                              
                !SWAT Base Heat Units counter (for planting schedule)
                if (Me%ComputeOptions%Evolution%ModelSWAT) then                    
                    call ComputePlantAnnualStage(i,j)                
                endif
                
                !Check if planting will occur
                if (.not. Me%IsPlantGrowing(i,j)) then
                    call CheckPlanting(i,j)                
                endif
                
                 !Check if it is pesticide application time
                if (Me%ComputeOptions%Pesticide) then
                    call CheckPlantPesticide(i,j)                
                endif                

                 !Check if it is pesticide application time
                if (Me%ComputeOptions%Fertilization .and. Me%ComputeOptions%FertilizationApp) then                    
                    call CheckPlantFertilization(i,j)                
                endif                


                !Logical to porous media properties interface
                if (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
                    Me%SoilFluxesActive(i,j) = .false.
                    if (Me%IsPlantGrowing(i,j)) then
                        Me%SoilFluxesActive(i,j) = .true.
                    else
                        if (Me%ComputeOptions%Fertilization) then
                            do FertApp = 1, Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%NumberFertilizerApps
                                if (Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppOccurred) then
                                    Me%SoilFluxesActive(i,j) = .true.
                                endif
                            enddo
                       endif
                    endif
                endif

                if (Me%IsPlantGrowing(i,j)) then
                    
                    !Check tree age counter
                    PlantType  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                    if (PlantType == Trees) then                        
                        call CheckTreeAge(i,j)                    
                    endif
                    
                    !Check if it is dormant period
                    if (Me%ComputeOptions%Dormancy) then                        
                        call CheckPlantDormancy(i,j)                    
                    endif
                    
                    !Check if it is grazing period
                    if (Me%ComputeOptions%Grazing) then                        
                        call CheckPlantGrazing(i,j)                    
                    endif
 
                    !Check if HarvestKill will occurr (harvest, kill... etc)
                    if (Me%ComputeOptions%HarvestKill) then                        
                        call CheckPlantHarvestKill(i,j)                    
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
        
        !Update year of simulation for fertilization target - used in autofertilization
        if (Me%ComputeOptions%AutoFertilization) then
            call JulianDay(Me%ExternalVar%Now, JulDay)
            JulDay_Old = Me%ExternalVar%JulianDay_Old
            
            if (JulDay_Old .gt. 364 .and. JulDay .ge. 1) then
                Me%SimulationYear = Me%SimulationYear + 1
            endif
        endif


    end subroutine CheckPlantState

    !--------------------------------------------------------------------------

    subroutine CheckVegetationIDChanges
        
        !Local-----------------------------------------------------------------
        integer                                      :: i, j, ivt, VegetationOld
        integer                                      :: STAT_CALL
        real, dimension(:,:), pointer                :: AgricPractID
        !Begin-----------------------------------------------------------------
        
        allocate (AgricPractID(Me%Worksize%ILB:Me%Worksize%IUB, Me%Worksize%JLB:Me%Worksize%JUB))
        
        !Update agricultural practices MatrixID - from user
        call ModifyFillMatrix (FillMatrixID   = Me%AgricPract%ObjFillMatrix,         &
                               Matrix2D       = AgricPractID,                        &
                               PointsToFill2D = Me%ExternalVar%MappingPoints,        &
                               STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckVegetationIDChanges - ModuleVegetation - ERR01'
        
        !Update VegetationID matrix - ordered ID's
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                Me%AgricPractID(i,j) = NINT(AgricPractID(i, j))
                
                !Serach for the right vegetation ID - ordered ID to save space
                do ivt = 1, Me%VegetationsNumber
                    
                    if (Me%VegetationTypes(ivt)%ID == NINT(AgricPractID(i, j))) then
                        
                        VegetationOld = Me%VegetationID(i, j)
                        Me%VegetationID(i, j) = ivt
                        
                        !changed agricultural practice - kill plant if growing
                        if (Me%IsPlantGrowing(i,j) .and. (VegetationOld .ne. Me%VegetationID(i, j))) then
                            Me%KillOccurred(i,j) = .true.
                        endif
                                                
                    endif
                enddo

            endif
        enddo
        enddo

    
    end subroutine CheckVegetationIDChanges
    
    !--------------------------------------------------------------------------    

    subroutine ModifyFluxes
 
       !Local-----------------------------------------------------------------
       !Begin-----------------------------------------------------------------
 
        call PlantRootFluxes

        call PlantAerialFluxes
 
        !Harvest and kill operations
        if (Me%ComputeOptions%HarvestKill) then
            !HarvestKill fluxes are computed if plant is being managed (kill, harvest)            
            call HarvestKillFluxes
        endif

        if (Me%ComputeOptions%Dormancy) then
            !Dormancy fluxes are computed if plant is going dormant (leaf fall, etc)
            call DormancyFluxes
        endif
    
        if(Me%ComputeOptions%Grazing) then
            !Grazing fluxes are computed if plant is being grazed
            call GrazingFluxes
        endif

        if(Me%ComputeOptions%Fertilization) then
            !Fertilization fluxes are computed (based on plant needs or scheduled)    
            call FertilizationFluxes
        endif
        
        if(Me%ComputeOptions%Pesticide) then
            !Pesticide fluxes are copmputed
            call PesticideFluxes
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
       integer :: i, j
       !Begin-----------------------------------------------------------------
 
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))

!            if (PropertyX%IsConstant) then
!                PropertyX%Field = PropertyX%ConstantValue

            if (PropertyX%ID%SolutionFromFile) then

                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,   &
                                       Matrix2D       = PropertyX%Field,              &
                                       PointsToFill2D = Me%ExternalVar%MappingPoints, &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write (*,*) "ATTENTION"
                    write (*,*) "Was not possible to read property '", trim(PropertyX%ID%Name), "' from file."  
                    stop 'ModifyReadedProperties - ModuleVegetation - ERR010'
                endif
            endif
            
            !This code exists to avoid Leaf Area Index (LAI) negative
            !The best is to preprocess data to avoid negative values (maybe interpolate to fill missing values)
            if (((PropertyX%ID%IDNumber == LeafAreaIndex_) .or. (PropertyX%ID%IDNumber == BoundaryLeafAreaIndex_)) .and. &
                (.not. Me%AllowNegativeLAI)) then
                if (Me%CorrectNegativeLAI) then            
do1:                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:                do i = Me%WorkSize%ILB, Me%WorkSize%IUB 
                        if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then                   
                            if (PropertyX%Field(i, j) < 0.0) then
                                PropertyX%Field(i, j) = Me%ValueInsteadNegativeLAI
                            endif
                        endif
                    enddo do2
                    enddo do1
                else
do3:                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do4:                do i = Me%WorkSize%ILB, Me%WorkSize%IUB   
                        if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then                 
                            if (PropertyX%Field(i, j) < 0.0) then
                                write (*,*) "ATTENTION."
                                write (*,*) "Found negative values for property '", trim(PropertyX%ID%Name), "'"
                                write (*,*) "Check the input file or the options ALLOW_NEGATIVE_LAI and"
                                write (*,*) "CORRECT_NEGATIVE_LAI"
                                Stop "ModifyReadedProperties - ModuleVegetation - ERR020"
                            endif
                        endif
                    enddo do4
                    enddo do3
                endif
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
        !----------------------------------------------------------------------

        Me%HeatUnits%PotentialHUBase_Old(i,j) = Me%HeatUnits%PotentialHUBase(i,j)

        call JulianDay(Me%ExternalVar%Now, JulDay)
        JulDay_Old = Me%ExternalVar%JulianDay_Old

        if (JulDay_Old .gt. 364 .and. JulDay .ge. 1) then
            
            !Base heat units nullified
            Me%HeatUnits%PotentialHUBase(i,j) = 0.0
            
            !update fertilization target
            if (Me%ComputeOptions%AutoFertilization) then
                Me%NitrogenYeldEstimate(i,j) = (Me%NitrogenYeldEstimate(i,j) * Me%SimulationYear + Me%NitrogenYeldTarget(i,j)) &
                                               / (Me%SimulationYear + 1)

            endif
            
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
        !----------------------------------------------------------------------

        call JulianDay(Me%ExternalVar%Now, JulDay)
        JulDay_Old = Me%ExternalVar%JulianDay_Old

        
        !Update if tree was just planted
!        if (Me%Growth%TreeCurrentYear(i,j) .eq. 0) then
!
!            !Tree maximum annual biomass update
!            TreeYearsToMaturity = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeYearsToMaturity
!            
!            if (TreeYearsToMaturity .gt. 0) then
!                
!                Me%Growth%TreeCurrentYear(i,j) = 1           
!            
!                TreeMaximumBiomass = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeMaximumBiomass
!                if (TreeMaximumBiomass .lt. 0.0) then
!                    write(*,*) 'Tree years to maturity was defined but maximum tree biomass did not.'
!                    write(*,*) 'Please check trees growth databases in vegetation file.'
!                    stop 'CheckTreeAge - Module Vegetation - ERR01'
!                endif
!
!                Me%Growth%TreeFractionToMaturity(i,j)   = float(Me%Growth%TreeCurrentYear(i,j)) / float(TreeYearsToMaturity)
!                Me%Growth%TreeMaximumAnnualBiomass(i,j) = Me%Growth%TreeFractionToMaturity(i,j) * TreeMaximumBiomass
!            else
!                Me%Growth%TreeFractionToMaturity(i,j) = 1.
!            endif
!       
!       !Update if tree was just planted, coming from continuos (lack of matrix values). update at beggining of year
        if (Me%PlantingOccurred(i,j) .or. Me%Growth%TreeComingFromContinuous(i,j) .or.                                   &
            (JulDay_Old .gt. 364 .and. JulDay .ge. 1)) then
            
            !Tree maximum annual biomass update
            TreeYearsToMaturity = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeYearsToMaturity
            
            if (TreeYearsToMaturity .gt. 0) then
                
                if (.not. Me%Growth%TreeComingFromContinuous(i,j)) then
                    Me%Growth%TreeCurrentYear(i,j)          = Me%Growth%TreeCurrentYear(i,j) + 1
                else
                    Me%Growth%TreeComingFromContinuous(i,j) = .false.
                endif

                Me%Growth%TreeCurrentYear(i,j)          = min(Me%Growth%TreeCurrentYear(i,j), TreeYearsToMaturity)

                TreeMaximumBiomass                      = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeMaximumBiomass
                if (TreeMaximumBiomass .lt. 0.0) then
                    write(*,*) 'Tree years to maturity was defined but maximum tree biomass did not.'
                    write(*,*) 'Please check trees growth databases in vegetation file.'
                    stop 'CheckTreeAge - Module Vegetation - ERR02'
                endif

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
        character (Len = StringLength)                  :: WarningString
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
            PotentialHUBase_Old = Me%HeatUnits%PotentialHUBase_Old(i,j)
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
            write(*,*    ) 'Warning ! Julian day and Base HU for planting' 
            write(*,*    ) 'inconsistently defined in crop database.'            
            stop 'CheckPlanting - ModuleVegetation - ERR10'  

        endif

        !Nullify Variables for new growth cycle
        if (Me%PlantingOccurred(i,j)) then
            
            WarningString = 'Planting'
            call UpdatePlantGrowingStage (i,j, WarningString)
            
!            if (Me%ComputeOptions%HarvestKill) then
!!                Me%HarvestOperations       = 1
!!                Me%HarvestKillOperations   = 1
!!                Me%KillOperations          = 1
!!                Me%HarvestFinished(i,j)    = .false.     
!            endif           
!            if (Me%ComputeOptions%Grazing) then
!!                Me%GrazingOperations       = 1
!!                Me%GrazingFinished(i,j)    = .false.
!            endif
            
        endif


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
        !Begin-----------------------------------------------------------------

        VegetationID         = Me%VegetationID(i,j)
        PlantHUatMaturity    = Me%VegetationTypes(VegetationID)%TimingDatabase%PlantHUatMaturity
        AverageTempDuringDT  = Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j)
        PlantBaseTemperature = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantBaseTemperature

        !! update accumulated heat units for the plant
        PlantHUVariation = 0.
        if (PlantHUatMaturity .gt. 0.1) then
            AccumulatedTemperature = AverageTempDuringDT * (Me%ComputeOptions%VegetationDT/86400.)
            PlantHUVariation = (AccumulatedTemperature - PlantBaseTemperature) / PlantHUatMaturity
        end if
        if (PlantHUVariation .lt. 0.) then
            PlantHUVariation = 0.
        endif
        Me%HeatUnits%PlantHUAccumulated_Old(i,j) = Me%HeatUnits%PlantHUAccumulated (i,j)
        Me%HeatUnits%PlantHUAccumulated    (i,j) = Me%HeatUnits%PlantHUAccumulated (i,j) + PlantHUVariation
        
        !Debug
        !HUAcc = Me%HeatUnits%PlantHUAccumulated (i,j)

    end subroutine ComputePlantGrowingStage

    !--------------------------------------------------------------------------

    subroutine PlantRootFluxes

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%Evolution%ModelSWAT) then
            
            if (Me%ComputeOptions%ModelWater) then
                
                call WaterUptake
            
            endif
             
         !! if plant is not dormant nor hasn't reached maturity, nutrient assimilation is allowed. 
            if (Me%ComputeOptions%ModelNitrogen) then
                
                call NitrogenUptake

!                call NitrogenFixationSWAT

            endif
            if (Me%ComputeOptions%ModelPhosphorus) then
            
                call PhosphorusUptake

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
            write(*,*) 'Transpiration method not known. Verify vegetation input'
            stop 'WaterUptake - ModuleVegetation - ERR100'
        endif

    end subroutine WaterUptake
    
    !--------------------------------------------------------------------------


    subroutine NitrogenUptake

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptakeSWAT) then
            
            call NitrogenUptakeSWAT
        
        elseif (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptake_TranspConc) then
            
            call NitrogenUptake_TranspConc
        
        else
            write(*,*) 
            write(*,*) 'Nutrient uptake method not known. Verify vegetation input'
            stop 'NitrogenUptake - ModuleVegetation - ERR100'
        endif

    end subroutine NitrogenUptake
    
    !--------------------------------------------------------------------------

    subroutine PhosphorusUptake

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
        
        if (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptakeSWAT) then
            
            call PhosphorusUptakeSWAT
        
        elseif (Me%ComputeOptions%NutrientUptakeMethod .eq. NutrientUptake_TranspConc) then
            
            call PhosphorusUptake_TranspConc
        
        else
            write(*,*) 
            write(*,*) 'Nutrient uptake method not known. Verify vegetation input'
            stop 'PhosphorusUptake - ModuleVegetation - ERR100'
        endif

    end subroutine PhosphorusUptake
    
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

    subroutine WaterUptakeSWAT

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i,j
        real                                            :: TopDepth, BottomDepth
        real                                            :: PotTP
        logical                                         :: FoundRoot
        real                                            :: NormalizationParameter
        real                                            :: RootDepth
        real                                            :: TopUptake, BottomUptake
        real                                            :: PotentialWaterUptake
        real                                            :: EffectiveWaterUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: LayerFieldCapacity
        real                                            :: LayerWaterContent
        real                                            :: SumDemand, SumUptake
        real                                            :: DemandNotMetInUpperLayers
        real                                            :: IncreaseUptake, ReduceUptake
        integer                                         :: k, KUB, KLB
        logical                                         :: UptakeAllowed
        
        !Begin-----------------------------------------------------------------
        
        !Aeration stress from SWAT not implmented yet
        
        !!  compute aeration stress
!       if (sol_sw(j) > sol_sumfc(j)) then
!          satco = (sol_sw(j) - sol_sumfc(j)) / (sol_sumul(j) - 
!     &                                                 sol_sumfc(j))
!          strsa(j) = 1. - (satco / (satco + Exp(.176 - 4.544 *
!     &                                                      satco)))
!        else
!          strsa(j) = 1.
!        end if
                

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                        
                UptakeAllowed = .true.
                if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                    !To by-pass no growing season periods. 
                    if (.not. Me%IsPlantGrowing(i,j)) then
                        UptakeAllowed = .false.
                    endif
                endif
               
!                Me%ExternalVar%PotentialTranspiration(i,j) = 1e-8 !m/s - 1 mm/dia
                ! mm         = m/s * s * 1000mm/m
                PotTP        = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT * 1000.0 
        
cd1:            if (.not. UptakeAllowed .or. PotTP .le. 0.01) then

                    Me%Fluxes%WaterUptake(i,j) = 0.0
                    Me%Growth%WaterStress(i,j) = 1.0

                else
                
                    BottomDepth  = 0.
                    FoundRoot = .false.
                    
!                    if (Me%StateVariables%RootDepth%Field(i,j) .lt. 0.01) then
!                        Me%StateVariables%RootDepth%Field(i,j) = 0.01
!                    endif
                    
                    !                               m
                    RootDepth    = Me%StateVariables%RootDepth(i,j)
!                    RootDepth     = PropRootDepth%Field(i,j)

                    GridCellArea = Me%ExternalVar%GridCellArea(i,j)
                    SumDemand    = 0.
                    SumUptake    = 0.

                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3 :               do k = KUB, KLB, -1
            
                        if (FoundRoot) then
                            Me%TranspirationBottomLayer(i,j) = k - 1
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
                        if (RootDepth .eq. 0.0) then
                            !    mm
                            PotentialWaterUptake = 0.0
!                        elseif (RootDepth .le. 1e-5) then
!                            !       mm
!                            PotentialWaterUptake = PotTP/NormalizationParameter
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
                        Me%Fluxes%WaterUptakeLayer(i,j,k) = EffectiveWaterUptake *  (1./1000.) * GridCellArea    &
                                                             * (1./Me%ComputeOptions%VegetationDT)

                    enddo do3
        
                    !    m3/s
                    Me%Fluxes%WaterUptake(i,j) = SumUptake *  (1./1000.) * GridCellArea * (1./Me%ComputeOptions%VegetationDT)

            !       strsw(j) = strsa (j) * xx / ep_max
                    Me%Growth%WaterStress(i,j) = SumUptake/PotTP
        
                endif cd1

            endif

        enddo do2
        enddo do1

        Me%Fluxes%FromSoil%WaterUptakeFromSoil => Me%Fluxes%WaterUptakeLayer

    end subroutine WaterUptakeSWAT

    !--------------------------------------------------------------------------

    subroutine WaterUptakeMOHID

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        real, dimension(:,:,:), pointer                 :: UnsatK
        integer                                         :: i, j
        real                                            :: SumUptake, RootDepth, BottomDepth
        real                                            :: TopDepth, CenterDepth, TotalCol
        real                                            :: FeddesH1, FeddesH2, FeddesH3, FeddesH4
        real                                            :: Head, Factor, h50, p1, WFactor, SFactor
        real                                            :: Aux, RootHeightInCell, CellBase
        real                                            :: ColToTransp, TranspirationDistribution
        real                                            :: SoilVolume, WaterVolume, TranspVolume
        real                                            :: ActualTranspiration, PotentialTranspiration
        logical                                         :: FoundRoot, UptakeAllowed
        integer                                         :: k, KLB, KUB, RootProfile
        real                                            :: VelocityVolume, NewTheta
        real                                            :: LimitThetaLow
        integer                                         :: StressInteraction
        type(T_SalinityStressParams), pointer           :: StressParams
        logical                                         :: UseSalinityStress
        integer                                         :: STAT_CALL
        real                                            :: P2L, P2H, R2L, R2H, TP        
        
        !Begin-----------------------------------------------------------------
        
        !Fluxes Nullification
        Me%Fluxes%WaterUptakeLayer(:,:,:) = 0.0
        Me%Fluxes%WaterUptake     (:,:  ) = 0.0

        call GetLimitThetaLow(Me%ObjPorousMedia, LimitThetaLow, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WaterUptakeMOHID - ModuleVegetation - ERR010'   
        
        if (Me%ComputeOptions%LimitTPVel) then 

            call GetUnsatK(Me%ObjPorousMedia, UnsatK, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WaterUptakeMOHID - ModuleVegetation - ERR020'   
        
        endif

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                        
                UptakeAllowed = .true.    
                if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                    !To by-pass no growing season periods. 
                    if (.not. Me%IsPlantGrowing(i,j)) then
                        UptakeAllowed = .false.
                    endif
                endif


                !Total Water Column to Remove (Potentialy) 
                !m = m/s * s
                TotalCol    = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j) * Me%ComputeOptions%VegetationDT    
                !  m
                RootDepth    = Me%StateVariables%RootDepth(i,j)
                BottomDepth = 0.0                   
            
cd1:            if (.not. UptakeAllowed .or. TotalCol .eq. 0.0 .or. RootDepth .eq. 0.0) then

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
                            Me%TranspirationBottomLayer(i,j) = k - 1
                            exit do3
                        endif

                        !Layer depths to consider transpiration
                        TopDepth        = BottomDepth
                        BottomDepth     = BottomDepth + Me%ExternalVar%DWZ(i, j, k)

                        !Center of layer
                        CenterDepth    = BottomDepth - Me%ExternalVar%DWZ(i,j,k) / 2.0 


                        if (RootDepth .le. BottomDepth ) then
                            FoundRoot = .true.
                            BottomDepth = RootDepth
                        endif

         
                        !Reduce Factor of EVTP due to water stress
                        if (.not. Me%ComputeOptions%WaterStress) then
                        
                            WFactor = 1.0
                            
                        elseif (Me%ComputeOptions%WaterUptakeStressMethod .eq. Feddes) then   

                            ! stress functions by Feddes
                            FeddesH1 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH1(i,j)
                            FeddesH2 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH2(i,j)
                            FeddesH4 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH4(i,j)
                            Head     = Me%ExternalVar%Head(i, j, k)

                            select case (Me%ComputeOptions%TranspirationMOHID%FeddesType(i,j))
                            case (1)
                                FeddesH3 = Me%ComputeOptions%TranspirationMOHID%RootFeddesH3(i,j)
                            case (2)
                                TP  = Me%ExternalVar%Integration%AveragePotTPDuringDT(i,j)
                                R2L = Me%ComputeOptions%TranspirationMOHID%RootFeddesR2L(i,j)
                                R2H = Me%ComputeOptions%TranspirationMOHID%RootFeddesR2H(i,j)
                                P2L = Me%ComputeOptions%TranspirationMOHID%RootFeddesH3L(i,j)
                                P2H = Me%ComputeOptions%TranspirationMOHID%RootFeddesH3H(i,j)
                                
                                if (TP > R2H) then
                                    FeddesH3 = P2H
                                elseif (TP < R2L) then
                                    FeddesH3 = P2L
                                else
                                    FeddesH3 = P2H + (P2L - P2H) / (R2H - R2L) * (R2H - TP)                                    
                                endif
                                
                            end select
                            
                            if      (Head .gt. FeddesH1) then
                                WFactor = 0.0
                            else if (Head .gt. FeddesH2) then
                                WFactor = LinearInterpolation(FeddesH2, 1.0, FeddesH1, 0.0, Head)
                            else if (Head .gt. FeddesH3) then
                                WFactor = 1.0
                            else if (Head .gt. FeddesH4) then
                                WFactor = LinearInterpolation(FeddesH4, 0.0, FeddesH3, 1.0, Head)
                            else
                                WFactor = 0.0
                            endif
          
                        elseif (Me%ComputeOptions%WaterUptakeStressMethod .eq. Genuchten) then  
        
                            !or stress functions by S-Shape (van Genuchten)
                            h50 =-8.0
                            p1 =3.0
                            WFactor= 1/(1+(Head/h50)**p1)

                        endif

                        call GetVegetationSalinityParameters (I, J, UseSalinityStress, StressParams, StressInteraction)
                
                        if (UseSalinityStress) then
                            select case (Me%ComputeOptions%SalinityStressMethod)
                            case (1) !Threshold/slope                        
                                !call ComputeECw (I, J, K, ECw)                               
                                
                                if (.not. (associated(Me%ExternalVar%ECw))) then
                                    stop 'WaterUptakeMOHID - ModuleVegetation - ERR030'
                                endif
                                
                                if (Me%ExternalVar%ECw(i, j, k) < StressParams%ECt) then
                                    SFactor = 1.0
                                else if (Me%ExternalVar%ECw(i, j, k) > StressParams%EC0) then
                                    SFactor = 0.0
                                else
                                    SFactor = 1.0 - (Me%ExternalVar%ECw(i, j, k) - StressParams%ECt) * 0.01 * StressParams%Slope
                                endif
                            case (2) !S shape
                                write (*,*) 'The S-Shape method for water uptake salinity stress (SALINITY_STRESS_METHOD keyword)'
                                write (*,*) 'isn''t implemented yet. Use 1 (Threshold/slope) instead.'
                                stop
                            case default
                                stop 'WaterUptakeMOHID - ModuleVegetation - ERR040'
                            end select
                    
                            select case (StressInteraction)
                            case (1) !Additive
                                Factor = abs(WFactor - SFactor)
                            case (2) !Multiplicative
                                Factor = Wfactor * SFactor
                            case (3) !Less of
                                Factor = min(WFactor, SFactor)
                            case default
                                stop 'WaterUptakeMOHID - ModuleVegetation - ERR050'
                            end select                       
                        else
                            Factor = WFactor
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

                            !Same method as first but fraction computed with layer numbers - if root exists inside layer 
                            !implicitily occupies all layer. 
                            !Removed because it would need more variables to compute and the geometry is the same as first option.
                !           elseif (RootProfile == 3) then !linear root distribution
                !   
                !                ColToEvap      =  2 * Factor * TotalCol * (1. / layerroot - (prof/layerroot**2))              
                ! 
                            elseif (RootProfile == RootExponential) then  ! exponential (from SWAT)

                !                a =10
                !                ColToEvap = Factor*TotalCol*a*exp(-a*(Me%ExtVar%RootDepth(i,j)-deep2))

                                TranspirationDistribution = (1.0 - exp(-10.0 * BottomDepth/RootDepth))/(1.0 - exp(-10.0))         &
                                                           - (1.0 - exp(-10.0 * TopDepth/RootDepth))/(1.0 - exp(-10.0))

                                ColToTransp = Factor * TotalCol * TranspirationDistribution
                            else
                               write(*,*)'Invalid Root Profile'
                               stop 'WaterUptakeMOHID - ModuleVegetation - ERR060'
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
                        
                            SoilVolume  = (Me%ExternalVar%SoilWaterContent(i,j,k) - Me%ExternalVar%ResidualWaterContent(i,j,k))  &
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
            if (STAT_CALL /= SUCCESS_) stop 'WaterUptakeMOHID - ModuleVegetation - ERR070'   
        
        endif

    end subroutine WaterUptakeMOHID

    !--------------------------------------------------------------------------

!    subroutine ComputeECw (I, J, K, ECw)
!    
!        !Arguments-------------------------------------------------------------
!        integer           :: I, J, K
!        real, intent(OUT) :: ECw
!        
!        !Local-----------------------------------------------------------------
!        type(T_SoilConc), pointer :: SoilConcX
!
!        !----------------------------------------------------------------------
!        ECw = 0.0
!        
!        SoilConcX => Me%FirstSoilConc
!        do while (associated (SoilConcX))
!            ECw = ECw + &
!                  (SoilConcX%Conc(I, J, K) * SoilConcX%ConversionFactor)
!        
!            SoilConcX => SoilConcX%Next
!            
!        end do
!        
!        if (associated(Me%ECw)) then
!            Me%ECw(I, J, K) = ECw
!        endif
!        !----------------------------------------------------------------------
!
!    end subroutine ComputeECw

    !--------------------------------------------------------------------------

    subroutine GetVegetationSalinityParameters (I, J, UseSalinityStress, StressParams, StressInteraction)
    
        !Arguments-------------------------------------------------------------
        integer                               :: I, J
        logical, intent(OUT)                  :: UseSalinityStress
        type(T_SalinityStressParams), pointer :: StressParams
        integer, intent(OUT)                  :: StressInteraction
        
        !Local-----------------------------------------------------------------
        type(T_VegetationType), pointer :: Vegetation_
        
        !----------------------------------------------------------------------
        Vegetation_  => Me%VegetationTypes(Me%VegetationID(I, J))
        
        UseSalinityStress =  Vegetation_%UseSalinityStress
        StressParams      => Vegetation_%SalinityStressParams        
        StressInteraction =  Vegetation_%StressInteraction
        !----------------------------------------------------------------------
        
    end subroutine GetVegetationSalinityParameters
    
    !--------------------------------------------------------------------------

    subroutine NitrogenUptakeSWAT

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i, j
        integer                                         :: VegetationID
        real                                            :: TotalPlantNitrogen, TotalPlantBiomass
        real                                            :: PredictedNitrogenBiomass
        real                                            :: PlantFractionN1, PlantFractionN2
        real                                            :: PlantFractionN3
        real                                            :: b1, b2, b3
        real                                            :: PlantShape1, PlantShape2
        real                                            :: NitrogenDemand, OptimalNContent
        real                                            :: Stress, HUAcc
        real                                            :: TopDepth, BottomDepth
        logical                                         :: FoundRoot
        real                                            :: NormalizationParameter
        real                                            :: RootDepth
        real                                            :: TopUptake, BottomUptake
        real                                            :: PotentialNitrogenUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: DistributionParameter
        real                                            :: LayerNitrogenContent
        real                                            :: SumDemand, SumUptake
        real                                            :: DemandNotMetInUpperLayers
        integer                                         :: k, KUB, KLB, PlantType
        logical                                         :: Dormant
        real                                            :: AmountAvailable
        
       !Begin-----------------------------------------------------------------        

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then

                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                Dormant = .false.
                if (Me%ComputeOptions%Dormancy) then
                    if (Me%IsPlantDormant(i,j)) then
                        Dormant = .true.
                    endif
                endif

                               
                if (Me%IsPlantGrowing(i,j) .and. .not. Dormant             &
                    .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then  
                
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
                    TotalPlantBiomass  = Me%StateVariables%TotalPlantBiomass(i,j)
                    TotalPlantNitrogen = Me%StateVariables%TotalPlantNitrogen(i,j)
                    OptimalNContent = Me%PlantNitrogenFraction(i,j) * TotalPlantBiomass
                    if (OptimalNContent .lt. TotalPlantNitrogen) then
                        OptimalNContent = TotalPlantNitrogen
                    endif
                    
                    Me%OptimalTotalPlantNitrogen(i,j) = OptimalNContent


                    !  kgN/ha
                    NitrogenDemand = OptimalNContent - TotalPlantNitrogen
                    
                    !This is taken from SWAT formulation and usually originates stresses to oscillate between 0 and 1
                    !because if low nutrient is available in one time step stress will be high (much lower than 1), 
                    !plant will not grow and with this formulation in the next step demand will be reduced because 
                    !plant did not grow and stress will be low (closer to 1). And so on so on
                    if (Me%ComputeOptions%NutrientReduceDemand) then
                        NitrogenDemand = min(4. * PlantFractionN3 * Me%Growth%BiomassGrowthOld(i,j), NitrogenDemand)
                    endif

    cd1:            if(NitrogenDemand .lt. 1e-5) then
        
                        Me%Fluxes%NitrogenUptake(i,j) = 0.0
                        Me%Growth%NitrogenStress(i,j) = 1.0
    
                    else
    
                        BottomDepth = 0.0
                        FoundRoot = .false.
                        !                               m
                        RootDepth             = Me%StateVariables%RootDepth(i,j)
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
                            if (RootDepth .eq. 0.0) then
                                PotentialNitrogenUptake = 0.0
                            else
                                !   kgN/ha
                                BottomUptake = NitrogenDemand *(1.0 - exp(-DistributionParameter * BottomDepth/RootDepth))     &
                                                /NormalizationParameter
                                TopUptake = NitrogenDemand *(1.0 - exp(-DistributionParameter * TopDepth/RootDepth))           &
                                                /NormalizationParameter
                                PotentialNitrogenUptake = BottomUptake - TopUptake
                            endif

                            DemandNotMetInUpperLayers = SumDemand - SumUptake
                            !Demand in next iteration
                            SumDemand = SumDemand + PotentialNitrogenUptake
                            
                            Me%Fluxes%NitrogenUptakeLayer(i,j,k) = PotentialNitrogenUptake + DemandNotMetInUpperLayers

                            CellVolume = Me%ExternalVar%CellVolume(i,j,k)
                            
                            !      kgN/ha        = gN/m3H20 * 1E-3kg/g * m3H20/m3cell * m3cell / (m2) * 10000m2/ha 
                            LayerNitrogenContent = Me%ExternalVar%SoilNitrate(i,j,k) * 1E-3                       &
                                                    * Me%ExternalVar%SoilWaterContent(i,j,k)                      &
                                                    * (CellVolume) / (GridCellArea) * 10000

                            
                            if (Me%ComputeOptions%NutrientReduceUptake) then
                               !Avoid taking more than exists because there is also transport that can remove mass
                               !Problem with this SWAT formulation where nutrient is completely disconnected from transpiration
                               !ex. If plant tries to remove per day more than half the amount that exists than nutrients
                               !will run out in two days. at this low amount plant will encounter limitation to uptake
                               
                               AmountAvailable = LayerNitrogenContent * Me%ComputeOptions%NutrientReduceUptakeFactor
                            
                            else
                                AmountAvailable = LayerNitrogenContent
                            endif
                                                           
                            if (Me%Fluxes%NitrogenUptakeLayer(i,j,k) .ge. AmountAvailable) then
                                Me%Fluxes%NitrogenUptakeLayer(i,j,k) = AmountAvailable
                            end if                                     
                            
                            SumUptake = SumUptake + Me%Fluxes%NitrogenUptakeLayer(i,j,k)
                        
                        enddo do3

                        Me%Fluxes%NitrogenUptake(i,j) = SumUptake
                        
                        !Compute Nitrogen Stress
                        PlantType = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                        select case (PlantType)
                            !Legumes
                            case (1,2,3)
                                Me%Growth%NitrogenStress(i,j) = 1.
                            case default
                                
                                if (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressSWAT) then
                                    
                                    PredictedNitrogenBiomass = TotalPlantNitrogen + SumUptake
                                    call ComputeNutrientStress(PredictedNitrogenBiomass,OptimalNContent,Stress)
                            
                                
                                elseif (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressUptake) then
                                    
                                    !stress update
                                    if (NitrogenDemand .gt. 1e-5) then
                                        Stress = SumUptake / NitrogenDemand
                                    else
                                        Stress = 1.0
                                    endif
                                    
                                    Stress = min (Stress, 1.0)
                                    
!                                   Stress = amax1(Stress, NewStress)
!                                   Stress = amin1(Stress, 1.0)
                                else

                                    write(*,*) 'Error in NutrientUptake method. Check vegetation options.'
                                    stop 'NitrogenUptakeSWAT - Module Vegetation ERR01'

                                endif
    
                        end select

                        Me%Growth%NitrogenStress(i,j) = Stress

                    endif cd1

                endif

            endif

        enddo do2
        enddo do1

        Me%Fluxes%FromSoil%NitrogenUptakeFromSoil => Me%Fluxes%NitrogenUptakeLayer

    end subroutine NitrogenUptakeSWAT

    !--------------------------------------------------------------------------

    subroutine NitrogenUptake_TranspConc

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i, j
        integer                                         :: VegetationID
        real                                            :: TotalPlantNitrogen, TotalPlantBiomass
        real                                            :: PredictedNitrogenBiomass
        real                                            :: PlantFractionN1, PlantFractionN2
        real                                            :: PlantFractionN3
        real                                            :: b1, b2, b3
        real                                            :: PlantShape1, PlantShape2
        real                                            :: NitrogenDemand, OptimalNContent
        real                                            :: Stress, HUAcc
        real                                            :: TopDepth, BottomDepth
        logical                                         :: FoundRoot
        real                                            :: RootDepth
        real                                            :: PotentialNitrogenUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: LayerNitrogenContent
        real                                            :: SumUptake
        integer                                         :: k, KUB, KLB, PlantType
        logical                                         :: ComputeCell !Dormant
        
       !Begin-----------------------------------------------------------------        

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then

                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                ComputeCell = .false.
                if (.not. Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                    ComputeCell = .true.
                else
!                    Dormant = .false.
!                    if (Me%ComputeOptions%Dormancy) then
!                        if (Me%IsPlantDormant(i,j)) then
!                            Dormant = .true.
!                        endif
!                    endif
                    
                    !only Plantgrowing can be used because is the condition for water.
                    !if any other condition (as dormant or HUacc then it will generate inconsistency
                    !because this is the model for Q*C
                    if (Me%IsPlantGrowing(i,j)) then
                        ComputeCell = .true.
                    endif
!                    if (Me%IsPlantGrowing(i,j) .and. .not. Dormant             &
!                    .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then
!                        ComputeCell = .true.
!                    endif
                endif 

                if (ComputeCell) then  

                    if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                
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
                        TotalPlantBiomass  = Me%StateVariables%TotalPlantBiomass(i,j)
                        TotalPlantNitrogen = Me%StateVariables%TotalPlantNitrogen(i,j)
                        OptimalNContent = Me%PlantNitrogenFraction(i,j) * TotalPlantBiomass
                        if (OptimalNContent .lt. TotalPlantNitrogen) then
                            OptimalNContent = TotalPlantNitrogen
                        endif

                        Me%OptimalTotalPlantNitrogen(i,j) = OptimalNContent


                        NitrogenDemand = OptimalNContent - TotalPlantNitrogen
                        
                        !This is taken from SWAT formulation and usually originates stresses to oscillate between 0 and 1
                        !because if low nutrient is available in one time step stress will be high (much lower than 1), 
                        !plant will not grow and with this formulation in the next step demand will be reduced because 
                        !plant did not grow and stress will be low (closer to 1). And so on so on
                        if (Me%ComputeOptions%NutrientReduceDemand) then                        
                            NitrogenDemand = min(4. * PlantFractionN3 * Me%Growth%BiomassGrowthOld(i,j), NitrogenDemand)
                        endif
                    endif
    

                    BottomDepth = 0.0
                    FoundRoot = .false.
                    !                               m
                    RootDepth             = Me%StateVariables%RootDepth(i,j)
                    SumUptake             = 0.
                    GridCellArea          = Me%ExternalVar%GridCellArea(i,j)

                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3:                do k = KUB, KLB, -1
    
                        if (FoundRoot) then
                            exit do3
                        endif
    
                        TopDepth = BottomDepth
                        BottomDepth = BottomDepth + Me%ExternalVar%DWZ(i,j,k)
    
                        if (RootDepth .le. BottomDepth ) then
                            BottomDepth = RootDepth
                            FoundRoot = .true.
                        endif
             
                        CellVolume = Me%ExternalVar%CellVolume(i,j,k)
                        
                        !    KgN/ha             =  gN/m3H20 * 1E-3kg/g *  m3/s * s / (m2) * 10000m2/ha 
                        PotentialNitrogenUptake = Me%ExternalVar%SoilNitrate(i,j,k) * 1E-3                    &
                                                  * Me%Fluxes%WaterUptakeLayer(i,j,k)                         &
                                                  * Me%ComputeOptions%VegetationDT                            &
                                                  / (GridCellArea) * 10000
    
                        !      kgN/ha        =  gN/m3H20 * 1E-3kg/g * m3H20/m3cell * m3cell / (m2) * 10000m2/ha                     
                        LayerNitrogenContent = Me%ExternalVar%SoilNitrate(i,j,k) * 1E-3                       &
                                                * Me%ExternalVar%SoilWaterContent(i,j,k)                      &
                                                * (CellVolume) / (GridCellArea) * 10000

                        Me%Fluxes%NitrogenUptakeLayer(i,j,k) = PotentialNitrogenUptake
                        if (LayerNitrogenContent .lt. Me%Fluxes%NitrogenUptakeLayer(i,j,k)) then
                            Me%Fluxes%NitrogenUptakeLayer(i,j,k) = LayerNitrogenContent
                        end if
    
                        SumUptake = SumUptake + Me%Fluxes%NitrogenUptakeLayer(i,j,k)
                    enddo do3

                    Me%Fluxes%NitrogenUptake(i,j) = SumUptake
                        
                    
                    if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                        
                        !Compute Nitrogen Stress
                        PlantType = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                        select case (PlantType)
                            !Legumes
                            case (1,2,3)
                                Me%Growth%NitrogenStress(i,j) = 1.
                            case default
                                
                                if (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressSWAT) then
                                    
                                    PredictedNitrogenBiomass = TotalPlantNitrogen + SumUptake
                                    call ComputeNutrientStress(PredictedNitrogenBiomass,OptimalNContent,Stress)
                            
                                
                                elseif (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressUptake) then
                                    
                                    !SWAT stress update
                                    if (NitrogenDemand .gt. AllmostZero) then
                                        Stress = SumUptake / NitrogenDemand
                                    else
                                        Stress = 1.0
                                    endif
                                    
                                    Stress = min (Stress, 1.0)
!                                   Stress = amax1(Stress, NewStress)
!                                   Stress = amin1(Stress, 1.0)
                                else

                                    write(*,*) 'Error in NutrientUptake method. Check vegetation options.'
                                    stop 'NitrogenUptakeSWAT - Module Vegetation ERR01'

                                endif
    
                        end select

                        Me%Growth%NitrogenStress(i,j) = Stress

                    endif

                endif

            endif

        enddo do2
        enddo do1

        Me%Fluxes%FromSoil%NitrogenUptakeFromSoil => Me%Fluxes%NitrogenUptakeLayer

    end subroutine NitrogenUptake_TranspConc

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
        integer                                         :: i, j
        integer                                         :: VegetationID
        real                                            :: TotalPlantBiomass, TotalPlantPhosphorus
        real                                            :: PredictedPhosphorusBiomass
        real                                            :: PlantFractionP1, PlantFractionP2
        real                                            :: PlantFractionP3
        real                                            :: b1, b2, b3
        real                                            :: PlantShape1, PlantShape2
        real                                            :: PhosphorusDemand, OptimalPContent
        real                                            :: Stress, HUAcc
        real                                            :: TopDepth, BottomDepth
        logical                                         :: FoundRoot
        real                                            :: NormalizationParameter
        real                                            :: RootDepth
        real                                            :: TopUptake, BottomUptake
        real                                            :: PotentialPhosphorusUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: DistributionParameter
        real                                            :: LayerPhosphorusContent
        real                                            :: SumDemand, SumUptake
        real                                            :: DemandNotMetInUpperLayers
        integer                                         :: k, KUB, KLB, PlantType
        logical                                         :: Dormant
        real                                            :: AmountAvailable
        
       !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                        
                Dormant = .false.
                if (Me%ComputeOptions%Dormancy) then
                    if (Me%IsPlantDormant(i,j)) then
                        Dormant = .true.
                    endif
                endif

               
                if (Me%IsPlantGrowing(i,j) .and. .not. Dormant                         &
                    .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then    

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
                    TotalPlantBiomass    = Me%StateVariables%TotalPlantBiomass(i,j)
                    TotalPlantPhosphorus = Me%StateVariables%TotalPlantPhosphorus(i,j)
                    OptimalPContent = Me%PlantPhosphorusFraction(i,j) * TotalPlantBiomass
                    if (OptimalPContent .lt. TotalPlantPhosphorus) then
                        OptimalPContent = TotalPlantPhosphorus
                    endif
        
                    Me%OptimalTotalPlantPhosphorus(i,j) = OptimalPContent


                    !  kgP/ha
                    PhosphorusDemand = OptimalPContent - TotalPlantPhosphorus
                    
                    !This is taken from SWAT formulation and usually originates stresses to oscillate between 0 and 1
                    !because if low nutrient is available in one time step stress will be high (much lower than 1), 
                    !plant will not grow and with this formulation in the next step demand will be reduced because 
                    !plant did not grow and stress will be low (closer to 1). And so on so on
                    if (Me%ComputeOptions%NutrientReduceDemand) then                    
                        PhosphorusDemand = min(4. * PlantFractionP3 * Me%Growth%BiomassGrowthOld(i,j), PhosphorusDemand)
                    endif
                    
                    !Luxury P uptake
                    PhosphorusDemand =  PhosphorusDemand * 1.5
        
    cd1:            if(PhosphorusDemand .lt. 1e-6) then
            
                        Me%Fluxes%PhosphorusUptake(i,j) = 0.0
                        Me%Growth%PhosphorusStress(i,j) = 1.0
        
                    else
        
                        BottomDepth  = 0.0
                        FoundRoot = .false.
                        !                               m
                        RootDepth    = Me%StateVariables%RootDepth(i,j)
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
                            if (RootDepth .eq. 0.0) then
                                PotentialPhosphorusUptake = 0.0
                            else
                                !   kgP/ha
                                BottomUptake = PhosphorusDemand *(1.0 - exp(-DistributionParameter * BottomDepth/RootDepth))    &
                                               /NormalizationParameter
                                TopUptake = PhosphorusDemand *(1.0 - exp(-DistributionParameter * TopDepth/RootDepth))          &
                                                /NormalizationParameter
                                PotentialPhosphorusUptake = BottomUptake - TopUptake
                            endif
  
                            DemandNotMetInUpperLayers = SumDemand - SumUptake
                            !Demand in next iteration
                            SumDemand = SumDemand + PotentialPhosphorusUptake
                            
                            Me%Fluxes%PhosphorusUptakeLayer(i,j,k) = PotentialPhosphorusUptake + DemandNotMetInUpperLayers
                            
                            CellVolume = Me%ExternalVar%CellVolume(i,j,k)
                            
                            !      kgP/ha          = gP/m3H20 * 1E-3kg/g  * m3H20/m3cell * m3cell / (m2) * 10000m2/ha 
                            LayerPhosphorusContent = Me%ExternalVar%SoilPhosphorus(i,j,k)  * 1E-3                   &
                                                     * Me%ExternalVar%SoilWaterContent(i,j,k)                       &
                                                     * (CellVolume) / (GridCellArea) * 10000
                            
                            if (Me%ComputeOptions%NutrientReduceUptake) then
                               !Avoid taking more than exists because there is also transport that can remove mass
                               !Problem with this SWAT formulation where nutrient is completely disconnected from transpiration
                               !ex. If plant tries to remove per day more than half the amount that exists than nutrients
                               !will run out in two days. at this low amount plant will encounter limitation to uptake                                
                                AmountAvailable = LayerPhosphorusContent * Me%ComputeOptions%NutrientReduceUptakeFactor
                            
                            else                            
                                AmountAvailable = LayerPhosphorusContent                           
                            endif
                            
                            if (Me%Fluxes%PhosphorusUptakeLayer(i,j,k) .ge. AmountAvailable) then
                                Me%Fluxes%PhosphorusUptakeLayer(i,j,k) = AmountAvailable
                            end if
                            
                            SumUptake = SumUptake + Me%Fluxes%PhosphorusUptakeLayer(i,j,k)
                        
                        enddo do3

                        Me%Fluxes%PhosphorusUptake(i,j) = SumUptake

                        !Compute Phosphorus Stress
                        PlantType = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                        select case (PlantType)
                            !Legumes
                            case (1,2,3)
                                Me%Growth%PhosphorusStress(i,j) = 1.
                            case default
                                
                                if (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressSWAT) then
                                    
                                    PredictedPhosphorusBiomass = TotalPlantPhosphorus + SumUptake
                                    call ComputeNutrientStress(PredictedPhosphorusBiomass,OptimalPContent,Stress)
                            
                                
                                elseif (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressUptake) then
                                    
                                    !SWAT stress update (adapted from nitrogen)
                                    if (PhosphorusDemand .gt. 1e-6) then
                                        Stress = SumUptake / PhosphorusDemand
                                    else
                                        Stress = 1.0
                                    endif
                                    
                                    Stress = min (Stress, 1.0)
                                    
!                                   Stress = amax1(Stress, NewStress)
!                                   Stress = amin1(Stress, 1.0)
                                else

                                    write(*,*) 'Error in NutrientUptake method. Check vegetation options.'
                                    stop 'NitrogenUptakeSWAT - Module Vegetation ERR01'

                                endif
    
                        end select

                        Me%Growth%PhosphorusStress(i,j) = Stress
                    
                    endif cd1

                endif
            
            endif

        enddo do2
        enddo do1

        Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil => Me%Fluxes%PhosphorusUptakeLayer


    end subroutine PhosphorusUptakeSWAT

    !--------------------------------------------------------------------------

    subroutine PhosphorusUptake_TranspConc

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i, j
        integer                                         :: VegetationID
        real                                            :: TotalPlantBiomass, TotalPlantPhosphorus
        real                                            :: PredictedPhosphorusBiomass
        real                                            :: PlantFractionP1, PlantFractionP2
        real                                            :: PlantFractionP3
        real                                            :: b1, b2, b3
        real                                            :: PlantShape1, PlantShape2
        real                                            :: PhosphorusDemand, OptimalPContent
        real                                            :: Stress, HUAcc
        real                                            :: TopDepth, BottomDepth
        logical                                         :: FoundRoot
        real                                            :: RootDepth
        real                                            :: PotentialPhosphorusUptake
        real                                            :: GridCellArea, CellVolume
        real                                            :: LayerPhosphorusContent
        real                                            :: SumUptake
        integer                                         :: k, KUB, KLB, PlantType
        logical                                         :: ComputeCell !Dormant
         
       !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%MappingPoints (i,j) == VegetationPoint) then
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                        
                ComputeCell = .false.
                if (.not. Me%ComputeOptions%Evolution%GrowthModelNeeded) then
                    ComputeCell = .true.
                else
!                    Dormant = .false.
!                    if (Me%ComputeOptions%Dormancy) then
!                        if (Me%IsPlantDormant(i,j)) then
!                            Dormant = .true.
!                        endif
!                    endif
                    
                    !only Plantgrowing can be used because is the condition for water.
                    !if any other condition (as dormant or HUacc then it will generate inconsistency
                    !because this is the model for Q*C                    
                    if (Me%IsPlantGrowing(i,j)) then
                        ComputeCell = .true.
                    endif
!                    if (Me%IsPlantGrowing(i,j) .and. .not. Dormant             &
!                    .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then
!                        ComputeCell = .true.
!                    endif
                endif 
                
                if (ComputeCell) then    
                    
                    if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

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
                        TotalPlantBiomass    = Me%StateVariables%TotalPlantBiomass(i,j)
                        TotalPlantPhosphorus = Me%StateVariables%TotalPlantPhosphorus(i,j)
                        OptimalPContent = Me%PlantPhosphorusFraction(i,j) * TotalPlantBiomass
                        if (OptimalPContent .lt. TotalPlantPhosphorus) then
                            OptimalPContent = TotalPlantPhosphorus
                        endif

                        Me%OptimalTotalPlantPhosphorus(i,j) = OptimalPContent

                        !  kgP/ha
                        PhosphorusDemand = OptimalPContent - TotalPlantPhosphorus
                            
                        !This is taken from SWAT formulation and usually originates stresses to oscillate between 0 and 1
                        !because if low nutrient is available in one time step stress will be high (much lower than 1), 
                        !plant will not grow and with this formulation in the next step demand will be reduced because 
                        !plant did not grow and stress will be low (closer to 1). And so on so on
                        if (Me%ComputeOptions%NutrientReduceDemand) then                        
                            PhosphorusDemand = min(4. * PlantFractionP3 * Me%Growth%BiomassGrowthOld(i,j), PhosphorusDemand)
                        endif
                        
                        !Luxury P uptake
                        PhosphorusDemand =  PhosphorusDemand * 1.5
                   
                    endif
        
                    BottomDepth  = 0.0
                    FoundRoot = .false.
                    !                               m
                    RootDepth    = Me%StateVariables%RootDepth(i,j)
                    SumUptake    = 0.
                    GridCellArea = Me%ExternalVar%GridCellArea(i,j)
    
                    KUB = Me%WorkSize%KUB
                    KLB = Me%WorkSize%KLB

do3 :               do k = KUB, KLB, -1
        
                        if (FoundRoot) then
                            exit do3
                        endif
        
                        TopDepth = BottomDepth
                        BottomDepth = BottomDepth + Me%ExternalVar%DWZ(i,j,k)
        
                        if (RootDepth .le. BottomDepth ) then
                            BottomDepth = RootDepth
                            FoundRoot = .true.
                        endif
                 
                        CellVolume = Me%ExternalVar%CellVolume(i,j,k)
                        
                        !    KgP/ha               =   gN/m3H20 * 1E-3kg/g * m3/s * s / (m2) * 10000m2/ha 
                        PotentialPhosphorusUptake = Me%ExternalVar%SoilPhosphorus(i,j,k) * 1E-3                    &
                                                    * Me%Fluxes%WaterUptakeLayer(i,j,k)                            &
                                                    * Me%ComputeOptions%VegetationDT                               &
                                                    / (GridCellArea) * 10000
        
                        !      kgP/ha          =  gP/m3H20 * 1E-3kg/g * m3H20/m3cell * m3cell / (m2) * 10000m2/ha 
                        LayerPhosphorusContent = Me%ExternalVar%SoilPhosphorus(i,j,k)  * 1E-3                                &
                                                 * Me%ExternalVar%SoilWaterContent(i,j,k)                                    &
                                                 * (CellVolume) / (GridCellArea) * 10000

                        Me%Fluxes%PhosphorusUptakeLayer(i,j,k) = PotentialPhosphorusUptake
                        
                        if (LayerPhosphorusContent .lt. Me%Fluxes%PhosphorusUptakeLayer(i,j,k)) then
                            Me%Fluxes%PhosphorusUptakeLayer(i,j,k) = LayerPhosphorusContent
                        end if
        
                        SumUptake = SumUptake + Me%Fluxes%PhosphorusUptakeLayer(i,j,k)
                        
                    enddo do3

                    Me%Fluxes%PhosphorusUptake(i,j) = SumUptake

                    if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

                        !Compute Phosphorus Stress
                        PlantType = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                        select case (PlantType)
                            !Legumes
                            case (1,2,3)
                                Me%Growth%PhosphorusStress(i,j) = 1.
                            case default
                                
                                if (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressSWAT) then
                                    
                                    PredictedPhosphorusBiomass = TotalPlantPhosphorus + SumUptake
                                    call ComputeNutrientStress(PredictedPhosphorusBiomass,OptimalPContent,Stress)
                            
                                
                                elseif (Me%ComputeOptions%NutrientStressMethod .eq. NutrientStressUptake) then
                                    
                                    !stress update (adapted from nitrogen)
                                    if (PhosphorusDemand .gt. AllmostZero) then
                                        Stress = SumUptake / PhosphorusDemand
                                    else
                                        Stress = 1.0
                                    endif
                                    
                                    Stress = min (Stress, 1.0)
                                    
!                                   Stress = amax1(Stress, NewStress)
!                                   Stress = amin1(Stress, 1.0)
                                else

                                    write(*,*) 'Error in NutrientUptake method. Check vegetation options.'
                                    stop 'NitrogenUptakeSWAT - Module Vegetation ERR01'

                                endif
    
                        end select

                        Me%Growth%PhosphorusStress(i,j) = Stress

                    endif

                endif
            
            endif

        enddo do2
        enddo do1

        Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil => Me%Fluxes%PhosphorusUptakeLayer


    end subroutine PhosphorusUptake_TranspConc

    !--------------------------------------------------------------------------

    subroutine NitrogenFixationSWAT
!This process is not yet implemented

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
        integer                                           :: i, j
        real                                              :: ExtinctCoef
        integer                                           :: VegetationID
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
        integer                                           :: PlantType
        logical                                           :: Dormant
        
       !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                
                
                !Allow Trees to grow
!                PlantType              =   Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
!                GrowAllowed = .false.
!                select case (PlantType)
!                    case (7)
!                        GrowAllowed = .true.
!                    case default
!                        if (Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then 
!                            GrowAllowed = .true.
!                        endif
!                end select

                !Check if dormant
                Dormant = .false.
                if (Me%ComputeOptions%Dormancy) then
                    if (Me%IsPlantDormant(i,j)) then
                        Dormant = .true.
                    endif
                endif

                if (Me%IsPlantGrowing(i,j) .and. .not. Dormant .and. Me%HeatUnits%PlantHUAccumulated (i,j) .le. 1.) then     
        
                    VegetationID           =   Me%VegetationID(i,j)
                    ExtinctCoef            =   Me%VegetationTypes(VegetationID)%GrowthDatabase%ExtinctCoef
                    BiomassEnergyRatio     =   Me%VegetationTypes(VegetationID)%GrowthDatabase%BiomassEnergyRatio
                    CO2ConcHigh            =   Me%VegetationTypes(VegetationID)%GrowthDatabase%CO2ConcHigh
                    BiomassEnergyRatioHigh =   Me%VegetationTypes(VegetationID)%GrowthDatabase%BiomassEnergyRatioHigh
                    RUEDeclineRate         =   Me%VegetationTypes(VegetationID)%GrowthDatabase%RUEDeclineRate
                    !  m2/m2
                    LAI                    =   Me%StateVariables%LeafAreaIndex(i,j)
    !                LAI                    =   PropLeafAreaIndex%Field(i,j)
                    !  MJ/m2               =   W/m2 * s * 1e-6MJ/J
                    SolarRadiation         =   Me%ExternalVar%Integration%AverageRadiationDuringDT(i,j)  &
                                               * Me%ComputeOptions%VegetationDT * 1e-6
                    PlantType              =   Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
       
                    !! calculate optimal biomass
                    !Photossintetically Active Radiation (MJ/m2)
                    PAR = 0.
                    PAR = .5 * SolarRadiation * (1. - Exp(-ExtinctCoef * (LAI + .05)))
                    
                    RUE = BiomassEnergyRatio

                    if (Me%ComputeOptions%AdjustRUEForCO2) then
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
                    endif

                    if (Me%ComputeOptions%AdjustRUEForVPD) then
                        !! adjust radiation-use efficiency for vapor pressure deficit. 
!                        !Compute VPD (KPa). Equations in ModuleBasin adapted to average atmosphere values
!                        SaturatedVapourPressure  = 0.6108 * exp (17.27 * Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j)  &
!                                                    / (Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) + 237.3))
                        !Compute VPD (Kpa) from SWAT equation
                        SaturatedVapourPressure = exp((16.78 * Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) - 116.9)    &
                                                    / (Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j) + 237.3))
                        
                        ActualVapourPressure  = SaturatedVapourPressure  &
                                               * Me%ExternalVar%Integration%AverageAirHumidityDuringDT(i, j)  
        
                        VapourPressureDeficit = SaturatedVapourPressure - ActualVapourPressure
        
                        !!assumes vapor pressure threshold of 1.0 kPa
                        if (VapourPressureDeficit > 1.0) then
                            Decline = VapourPressureDeficit - 1.0
                            RUE = RUE - RUEDeclineRate * Decline
                            RUE = Max(RUE, 0.27 * BiomassEnergyRatio)
                        end if
                    endif

                    PotentialBiomassGrowth = RUE * PAR
                    if (PotentialBiomassGrowth .lt. 0.) then
                        PotentialBiomassGrowth = 0.
                    endif
                
                    if(Me%ComputeOptions%AtmospherePropertiesOutput) then
                        Me%Growth%PAR(i,j)             = PAR
                        Me%Growth%RUE(i,j)             = RUE
                        Me%Growth%PotentialGrowth(i,j) =  PotentialBiomassGrowth
                        Me%Growth%PotentialBiomass(i,j) = Me%StateVariables%TotalPlantBiomass(i,j) + PotentialBiomassGrowth
                    endif

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
                
                endif

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

        if (Me%ComputeOptions%ModelWater .and. Me%ComputeOptions%WaterStress) then
            WaterStress = Me%Growth%WaterStress(i,j)
        else
            WaterStress = 1.0
        endif
                
        if (Me%ComputeOptions%ModelTemperatureStress) then
            !Temperature stress must be computed now, because it is dependent on the plant characteristics
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

        !Begin-----------------------------------------------------------------        
        BaseTemperature     = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantBaseTemperature
        OptimalTemperature  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantOptimalTemperature
        AverageTempDuringDT = Me%ExternalVar%Integration%AverageAirTempDuringDT(i,j)

        Temperature1 = 0.
        Temperature1 = AverageTempDuringDT - BaseTemperature

        if (Temperature1 .le. 0.) then
            Me%Growth%TemperatureStress(i,j) = 0.
        else
            if (AverageTempDuringDT .gt. OptimalTemperature) then
                Temperature1 = 2. * OptimalTemperature - BaseTemperature - AverageTempDuringDT
            end if

            Temperature2 = 0.
            !Ignored swat code adding 1E-6
            if (Temperature1 .ne. 0.0) then
                Temperature2 = ((OptimalTemperature - AverageTempDuringDT) / (Temperature1)) ** 2
            endif

            if (Temperature2 .le. 200. .and. Temperature1 .gt. 0.) then
                Me%Growth%TemperatureStress(i,j) = Exp(-0.1054 * Temperature2)
            else
                Me%Growth%TemperatureStress(i,j) = 0.
            end if
        endif

    end subroutine ComputeTemperatureStress

    !--------------------------------------------------------------------------

    subroutine LAIGrowthSWAT

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i, j
        integer                                         :: veg_id
        type(T_GrowthDatabase), pointer                 :: gdb !gdb is a pointer to the GrowthDataBase structure        
        real                                            :: HUAcc, HUAcc_Old        
        real                                            :: FractionLAIMax_Old                        
        real                                            :: LAIMax, LAI
        real                                            :: DeltaFractionLAIMax
        real                                            :: LAIDeclineFraction
        real                                            :: GlobalStress
        type (T_Property), pointer                      :: BoundaryLAI
        integer                                         :: STAT_CALL
        logical                                         :: UseBoundaryLAI
        logical                                         :: LAIIsGrowing
                
        !Begin----------------------------------------------------------        

        UseBoundaryLAI = .false.

        if (Me%UseBoundaryLAI) then
            call SearchProperty(BoundaryLAI, BoundaryLeafAreaIndex_, .false., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'LAIGrowthSWAT - ModuleVegetation - ERR001'
            endif
        endif

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then                                         
            
                veg_id = Me%VegetationID(i,j)
                
                gdb => Me%VegetationTypes(veg_id)%GrowthDatabase
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (gdb%PlantType == NotAPlant) cycle                                            

                if (Me%UseBoundaryLAI) then
                    !If LAI is negative than it do not exist for this DT.
                    !In this case, the LAI must be computed.
                    if (BoundaryLAI%Field(i,j) >= 0) then
                        UseBoundaryLAI = .true.
                    else
                        UseBoundaryLAI = .false.
                    endif
                endif

                !====Start checking if growth must be computed====
                !If plant do not has leaves, go to next cell
                if (.not. Me%VegetationTypes(veg_id)%HasLeaves) cycle
                    
                !If plant is dorment, go to next cell
                if (Me%ComputeOptions%Dormancy) then
                    if (Me%IsPlantDormant(i,j)) then                    
                        Me%LAISenescence(i,j) = .false.
                        
                        if (UseBoundaryLAI) then
                            Me%Fluxes%LAIChange(i,j) = BoundaryLAI%Field(i,j) - Me%StateVariables%LeafAreaIndex(i,j)    
                        endif                        
                        
                        cycle
                    endif
                endif                
                                
                !If plant is not growing, go to next cell
                if (.not. Me%IsPlantGrowing(i,j)) then
                    if (UseBoundaryLAI) then
                        Me%Fluxes%LAIChange(i,j) = BoundaryLAI%Field(i,j) - Me%StateVariables%LeafAreaIndex(i,j)    
                    endif
                    
                    cycle
                endif
                
                !If plant already reached maturity, go to next cell
                LAIIsGrowing = .false.
                HUAcc = Me%HeatUnits%PlantHUAccumulated(i,j)                
                if (HUAcc <= 1.0) then
                    LAIIsGrowing = .true.
                else if (HUAcc_Old < 1.0) then
                    LAIIsGrowing = .true.
                else 
                    if (UseBoundaryLAI) then
                        Me%Fluxes%LAIChange(i,j) = BoundaryLAI%Field(i,j) - Me%StateVariables%LeafAreaIndex(i,j)    
                    endif                
                
                    cycle
                endif
                !==================================================                                                                        
                
                HUAcc_Old = Me%HeatUnits%PlantHUAccumulated_Old(i,j)

                !Save old value. It will be used in next step                
                if (Me%PlantingOccurred(i,j)) then
                    FractionLAIMax_Old = 0.0
                else
                    FractionLAIMax_Old = Me%PlantLAIMaxFraction(i,j)
                endif                  
                
                !Update PlantLAIMaxFraction
                Me%PlantLAIMaxFraction(i,j) = HUAcc / (HUAcc + exp(gdb%LAIShape1 - gdb%LAIShape2 * HUAcc))                
                                                
                !Compute new Leaf Area Index                     
                if (LAIIsGrowing) then
                
                    Me%LAISenescence(i,j) = .false.
                
                    if (UseBoundaryLAI) then
                    
                        Me%Fluxes%LAIChange(i,j) = BoundaryLAI%Field(i,j) - Me%StateVariables%LeafAreaIndex(i,j)
                    
                    else                                        
            
                        !Gets the stress calculated before. Only temperature stress is computed at this routine.
                        call ComputeGlobalStress(i,j)
                        GlobalStress = Me%Growth%GlobalStress(i,j)
                        
                        if (gdb%PlantType .eq. 7) then
                            LAIMax = gdb%LAIMax * Me%Growth%TreeFractionToMaturity(i,j)
                        else
                            LAIMax = gdb%LAIMax
                        endif
            
                        DeltaFractionLAIMax = Me%PlantLAIMaxFraction(i,j) - FractionLAIMax_Old                                
                        LAI = Me%StateVariables%LeafAreaIndex(i,j)
                        Me%Fluxes%LAIChange(i,j) = (DeltaFractionLAIMax * LAIMax * (1.0 - exp(5.0 * (LAI - LAIMax)))) * sqrt(GlobalStress)

                        !When coming from dormancy or grazing or harvest LAI growth turns negative because of HU decline
                        !However, LAI decline due to HUAccumulated decline is already computed explicitly 
                        !for grazing and harvesting in leaf are index state variable. In order to avoid double account here
                        !flux can only be positive.
                        if (Me%Fluxes%LAIChange(i,j) .lt. 0.0) then
                            Me%Fluxes%LAIChange(i,j) = 0.00
                        elseif ((LAI + Me%Fluxes%LAIChange(i,j)) > LAIMax) then
                            Me%Fluxes%LAIChange(i,j) = LAIMax - LAI
                        endif
                        
                    endif
    
                else if (.not. gdb%Evergreen) then 
                !if LAI is not growing and is not a "evergreen" plant, then the lai is declining
        
                    Me%LAISenescence(i,j) = .true.
                           
                    if (Me%ComputeOptions%UseLAIDecRate) then                    
                        LAIDeclineFraction = 1.0 - gdb%LAIDecRate                        
                    elseif(.not. Me%ComputeOptions%ChangeLAISenescence) then                    
                        LAIDeclineFraction = (1.0 - HUAcc) / (1.0 - gdb%FrGrowLAIDecline)                        
                    else
                        if (HUAcc .gt. HUAcc_Old) then                        
                            LAIDeclineFraction = (1.0 - HUAcc) / (1.0 - HUAcc_Old)                            
                        else                        
                            LAIDeclineFraction = 1.0                            
                        endif
                    endif
                        
                    if (UseBoundaryLAI) then                        
                        Me%Fluxes%LAIChange(i,j) = BoundaryLAI%Field(i,j) - Me%StateVariables%LeafAreaIndex(i,j)
                    else
                        Me%Fluxes%LAIChange(i,j) = Me%StateVariables%LeafAreaIndex(i,j) * (1.0 - LAIDeclineFraction)
                    endif
                      
                else if (UseBoundaryLAI) then
                !If LAI is not growing and not declining (it's a evergreen plant at maturity), 
                !actualize LAI if using boundary condition
                
                    Me%Fluxes%LAIChange(i,j) = BoundaryLAI%Field(i,j) - Me%StateVariables%LeafAreaIndex(i,j)    
                                        
                endif

            endif

        enddo do2
        enddo do1


    end subroutine LAIGrowthSWAT

    !--------------------------------------------------------------------------

    subroutine CheckPlantHarvestKill(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: HUAcc, HUAcc_Old
        real                                            :: HarvestKillJulianDay 
        real                                            :: HarvestJulianDay, KillJulianDay
        integer                                         :: PlantType, JulDay, JulDay_Old
!        integer                                         :: Op
        real                                            :: HarvestKillPlantHU 
        real                                            :: KillPlantHU
        real                                            :: HarvestHU
        logical                                         :: Dormant
        !Begin-----------------------------------------------------------------


        HUAcc                = Me%HeatUnits%PlantHUAccumulated(i,j)
        HUAcc_Old            = Me%HeatUnits%PlantHUAccumulated_Old(i,j)
        PlantType            = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
        HarvestKillJulianDay = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%HarvestKillJulianDay
        HarvestKillPlantHU   = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%HarvestKillPlantHU
        HarvestJulianDay     = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%HarvestJulianDay
        HarvestHU            = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%HarvestPlantHU
        KillJulianDay        = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%KillJulianDay
        KillPlantHU          = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%KillPlantHU

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
!        if (.not. Me%HarvestFinished(i,j)) then
            
!            Op = Me%HarvestOperations(i,j)
            
            if (HarvestJulianDay .gt. 0.0) then
                if(JulDay .ge. HarvestJulianDay .and. JulDay_Old .lt. HarvestJulianDay) then
            
                    Me%HarvestOnlyOccurred(i,j) = .true.
            
                endif
            elseif (HarvestHU .gt. 0.0) then
                if(HUAcc .ge. HarvestHU .and. HUAcc_Old .lt. HarvestHU) then
            
                    Me%HarvestOnlyOccurred(i,j) = .true.
                
!                    Me%HarvestOperations(i,j) = Me%HarvestOperations(i,j) + 1
!                    if (Me%HarvestOperations(i,j) .gt. size(HarvestHU)) then
!                        Me%HarvestFinished(i,j) = .true. 
!                    endif

                endif
            else
!                Me%HarvestFinished(i,j) = .true. 
            endif
!        endif


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
        
!        deallocate(Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%HarvestPlantHU)
   
    end subroutine CheckPlantHarvestKill

    !--------------------------------------------------------------------------
    
    subroutine HarvestKillFluxes

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i, j
        character (Len = StringLength)                  :: WarningString
        
        !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint  .and. Me%IsPlantGrowing(i,j)) then  

                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                if (Me%HarvestKillOccurred(i,j)) then
                    call HarvestKillOperation(i,j)
                    WarningString = 'Kill'
                    call UpdatePlantGrowingStage(i,j, WarningString)
                endif

                if (Me%HarvestOnlyOccurred(i,j)) then
                    call HarvestOperation(i,j)
                    WarningString = 'Harvest'
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

   
    end subroutine HarvestKillFluxes
    
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

        TotalPlantBiomass = Me%StateVariables%TotalPlantBiomass(i,j)
        AerialBiomass     = TotalPlantBiomass - Me%StateVariables%RootBiomass(i,j)
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
        HarvestEfficiency = Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%HarvestEfficiency
        Clip = 0.0
        Clip = Yeld * (1. - HarvestEfficiency)
        Yeld = Yeld *  HarvestEfficiency
        
        if (Yeld .lt. 0.0) then
            Yeld = 0.0
        endif
        if (Clip .lt. 0.0) then
            Clip = 0.0
        endif

        !!Biomass to soil. The fraction left by clippings
        Me%Fluxes%ToSoil%HarvestKillBiomassToSoil(i,j) = Clip
        
        !StateVariables Fluxes
        Me%Fluxes%BiomassRemovedInHarvest(i,j)    = Yeld + Clip
        
        TotalPlantBiomass        = Me%StateVariables%TotalPlantBiomass(i,j)
        
        !For HU update
        BiomassHarvestedFraction = (Yeld + Clip) / TotalPlantBiomass
        Me%Fluxes%BiomassHarvestedFraction = BiomassHarvestedFraction

        

        if (Me%ComputeOptions%ModelNitrogen) then

            TotalPlantNitrogen       = Me%StateVariables%TotalPlantNitrogen(i,j)
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
            Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil(i,j) = NitrogenClipping
 
            !StateVariables Fluxes
            Me%Fluxes%NitrogenRemovedInHarvest(i,j)   = NitrogenYeld + NitrogenClipping
        
        endif

        if (Me%ComputeOptions%ModelPhosphorus) then
        
            TotalPlantPhosphorus     = Me%StateVariables%TotalPlantPhosphorus(i,j)
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
            Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil(i,j) = PhosphorusClipping

            !StateVariables Fluxes
            Me%Fluxes%PhosphorusRemovedInHarvest(i,j) = PhosphorusYeld + PhosphorusClipping
        
        endif

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
        real                                            :: NitrogenToSoil, PhosphorusToSoil, FertilizationEfficiency
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

        TotalPlantBiomass = Me%StateVariables%TotalPlantBiomass(i,j)
        AerialBiomass     = TotalPlantBiomass - Me%StateVariables%RootBiomass(i,j)
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
        
        !!Biomass to soil. The fraction not removed by yeld because plant will die
        Me%Fluxes%ToSoil%HarvestKillBiomassToSoil(i,j) = Residue
        !The remaining biomass in soil - Not accounted in SWAT because N and P to soil come from all plant
        Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil(i,j) = Me%StateVariables%RootBiomass(i,j)

        if (Me%ComputeOptions%ModelNitrogen) then        
            
            NitrogenYeld   = 0.0
            NitrogenFractionInYeld   = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%NitrogenFractionInYeld            
            TotalPlantNitrogen       = Me%StateVariables%TotalPlantNitrogen(i,j)
            
            NitrogenYeld   = Yeld * NitrogenFractionInYeld
            NitrogenYeld   = min (NitrogenYeld,   0.9 * TotalPlantNitrogen)

            !!Nitrogen to soil. The fraction not removed by yeld because plant will die
            if (TotalPlantNitrogen .gt. NitrogenYeld) then
                NitrogenToSoil = TotalPlantNitrogen - NitrogenYeld
            else
                NitrogenToSoil = 0.0
            endif

            Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil(i,j) = NitrogenToSoil

        endif

        if (Me%ComputeOptions%ModelPhosphorus) then        

            PhosphorusYeld = 0.0
            PhosphorusFractionInYeld = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PhosphorusFractionInYeld
            TotalPlantPhosphorus     = Me%StateVariables%TotalPlantPhosphorus(i,j)

            PhosphorusYeld = Yeld * PhosphorusFractionInYeld
            PhosphorusYeld = min (PhosphorusYeld, 0.9 * TotalPlantPhosphorus)

            !!Phosphorus to soil. The fraction not removed by yeld because plant will die
            if (TotalPlantPhosphorus .gt. PhosphorusYeld) then
                PhosphorusToSoil = TotalPlantPhosphorus - PhosphorusYeld
            else
                PhosphorusToSoil = 0.0
            endif
            
            Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil(i,j) = PhosphorusToSoil

        endif

        !Biomass or nutrients fluxes for plant are not computed because plant dies.
        !Instead a warning variable is constructed and when state variables are updated, 
        !model knows that it needs to kill plant - biomass to zero.

        
        Me%Growth%WaterStress(i,j)       = 1.0
        Me%Growth%TemperatureStress(i,j) = 1.0
        Me%Growth%NitrogenStress(i,j)    = 1.0
        Me%Growth%PhosphorusStress(i,j)  = 1.0
        Me%Growth%TreeCurrentYear(i,j)   = 0
        
        if (Me%ComputeOptions%AutoFertilization) then
            AerialBiomass            = Me%StateVariables%TotalPlantBiomass(i,j) - Me%StateVariables%RootBiomass(i,j)
            NitrogenFractionInYeld   = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%NitrogenFractionInYeld
            FertilizationEfficiency  = 1.0
            
            !kgN/ha
            Me%NitrogenYeldTarget(i,j) = AerialBiomass * NitrogenFractionInYeld * FertilizationEfficiency
        endif
        
    end subroutine HarvestKillOperation

    !--------------------------------------------------------------------------

    subroutine KillOperation(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: AerialBiomass, Residue

        AerialBiomass     = Me%StateVariables%TotalPlantBiomass(i,j) - Me%StateVariables%RootBiomass(i,j)
        Residue           = 0.0

        !Biomass to soil
        Residue = AerialBiomass
        Me%Fluxes%ToSoil%HarvestKillBiomassToSoil(i,j) = Residue
        !The remaining biomass in soil - Not accounted in SWAT because N and P to soil come from all plant
        Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil(i,j) = Me%StateVariables%RootBiomass(i,j)

        if (Me%ComputeOptions%ModelNitrogen) then
            !Nitrogen to soil
            Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil(i,j) = Me%StateVariables%TotalPlantNitrogen(i,j)
        endif
        
        if (Me%ComputeOptions%ModelPhosphorus) then
            !Phosphorus to soil
            Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil(i,j) = Me%StateVariables%TotalPlantPhosphorus(i,j)
        endif

        !Biomass or nutrients fluxes for plant are not computed because plant dies.
        !Instead a warning variable is constructed and when state variables are updated, 
        !model knows that it needs to kill plant.

        Me%Growth%WaterStress(i,j)       = 1.0
        Me%Growth%TemperatureStress(i,j) = 1.0
        Me%Growth%NitrogenStress(i,j)    = 1.0
        Me%Growth%PhosphorusStress(i,j)  = 1.0
        Me%Growth%TreeCurrentYear(i,j)   = 0

    end subroutine KillOperation

    !--------------------------------------------------------------------------

    subroutine ComputeDayLength

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: JulDay
        real                                            :: SunDeclination 
        real                                            :: EarthRotationVelocity, Aux, AuxMin
        real                                            :: LatitudeRad, TresholdForDormancy
        integer                                         :: STAT_CALL, i,j
        
        !Begin-----------------------------------------------------------------
        
        call GetGridLatitudeLongitude(Me%ObjHorizontalGrid,                     &
                                      GridLatitude  = Me%ExternalVar%Latitude,  &
                                      STAT          = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDayLenght - ModuleVegetation - ERR01'


        call JulianDay(Me%ExternalVar%Now, JulDay)
        
        !Sun Declination (from atmosphere module) in radians
        !Possibly taken from Deas, M.L. and Lowney C.L. - "Water Temperature Modelling Review" Sept. 2000.
        SunDeclination = 23.45 * cos((172.0 - JulDay) * 2.0 * PI / 365.0) * PI / 180.0
        
        EarthRotationVelocity = 0.2618 !rad/hr = 15deg/hr

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then

                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                ! Compute day length (h)
                LatitudeRad = Me%ExternalVar%Latitude(i,j) * PI / 180.0
            
                Aux = -tan(SunDeclination) * tan(LatitudeRad)
            
                !Day Lenght equation does not apply for latitudes greater than +/- 66.5
                !When latitude exceeds +/-66.5 in winter there is no sun for 24h (Aux > 1)
                if(Aux .gt. 1.0) then
                
                    Me%DayLength(i,j) = 0.0
            
                ! between 66.5 and -66.5 deg in latitude (Aux between -1 and 1)
                elseif (Aux .ge. -1.0) then 
                
                    !SWAT 2005 equation 1:1.1.6
                    Me%DayLength(i,j) = 2 * acos(Aux) / EarthRotationVelocity
            
                !When latitude exceeds +/-66.5 in summer there is sun for 24h (Aux < -1)
                else
                
                    Me%DayLength(i,j) = 24.0 
            
                endif


                ! Compute minimum day length (h) (when occurs -23.45 deg sun declination (-0.40928rad) in northern 
                !hemisphere and + 23.45 deg sun declination (+0.40928rad) in southern hemisphere)
                if(Me%ExternalVar%Latitude(i,j) .ge. 0.0) then
                    AuxMin = -tan(-0.40928) * tan(LatitudeRad)
                else
                    AuxMin = -tan(0.40928) * tan(LatitudeRad)
                endif
             
                !Minimum Day Lenght equation does not apply for latitudes greater than +/- 66.5 (AuxMin > 1.0)
                !For latitudes greater than +/- 66.5 the minimum day length is 0h (in summer or winter
                !deppending on hemisphere there will be no sunrise)
                if(AuxMin .ge. 1.0) then
                    Me%MinimumDayLength(i,j) = 0.0
                else
                    Me%MinimumDayLength(i,j) = 2 * acos(AuxMin) / EarthRotationVelocity
                endif

                !Compute day lenght treshold from miminum (h) and update mimimum allowed
                if (abs(Me%ExternalVar%Latitude(i,j)) .gt. 40.) then
                    TresholdForDormancy = 1.
                elseif (abs(Me%ExternalVar%Latitude(i,j)) .lt. 20.) then
                    TresholdForDormancy = 0.
                else
                    TresholdForDormancy = (abs(Me%ExternalVar%Latitude(i,j)) - 20.) / 20.
                end if

                Me%MinimumDayLength(i,j) = Me%MinimumDayLength(i,j) + TresholdForDormancy

            endif

        enddo do2
        enddo do1

        !Latitude
        call UngetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%Latitude, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDayLenght - ModuleVegetation - ERR02'


    endsubroutine ComputeDayLength

    !--------------------------------------------------------------------------

    subroutine CheckPlantDormancy(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        integer                                         :: VegetationID, PlantType
        !Begin-----------------------------------------------------------------


        VegetationID = Me%VegetationID(i,j)
        PlantType = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType

        !! check for beginning of dormant season
        if (.not. Me%IsPlantDormant(i,j) .and. (Me%DayLength(i,j) .lt. Me%MinimumDayLength(i,j))) then

          
            select case (PlantType)

                case (3,6,7)
            
                    Me%PlantGoingDormant(i,j) = .true.
                    Me%IsPlantDormant   (i,j) = .true.
                
                !! beginning of cool season annual dormant period
                case (2,5)
                    
                    if (Me%HeatUnits%PlantHUAccumulated (i,j) .lt. 0.75) then
                        Me%PlantGoingDormant(i,j) = .true.
                        Me%IsPlantDormant   (i,j) = .true.
                    endif
                        
            end select
        endif

        !! check for end of dormant season
        if (Me%IsPlantDormant(i,j) .and. (Me%DayLength(i,j) .gt. Me%MinimumDayLength(i,j))) then

            select case (PlantType)

                case (2,3,5,6,7)
            
                    Me%IsPlantDormant(i,j)    = .false.
            
            end select
        
        endif


    end subroutine CheckPlantDormancy

    !--------------------------------------------------------------------------
    
    subroutine DormancyFluxes
       
        !Local-----------------------------------------------------------------        
        integer                                            :: i, j, VegetationID
        integer                                            :: PlantType
        logical                                            :: PlantKilled
        real                                               :: TotalPlantBiomass    
        real                                               :: BiomassRemovedInDormancy
        real                                               :: BiomassFracRemovedInDormancy
        real                                               :: BiomassRemovedInHarvest
        real                                               :: PredictedBiomass
        character (Len = StringLength)                     :: WarningString

        !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                                
                PlantKilled = .false.
                if (Me%ComputeOptions%HarvestKill) then
                    if (Me%HarvestKillOccurred(i,j) .or. Me%KillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif
                
                if (Me%IsPlantGrowing(i,j) .and. Me%PlantGoingDormant(i,j)      &
                    .and. .not. PlantKilled) then     

                    VegetationID = Me%VegetationID(i,j)
                    PlantType    = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
            
                    select case (PlantType)
                        !Trees and perennial
                        case (7, 3, 6)
                    
                            BiomassFracRemovedInDormancy =                                          &
                            Me%VegetationTypes(VegetationID)%GrowthDatabase%BiomassFracRemovedInDormancy
                            TotalPlantBiomass = Me%StateVariables%TotalPlantBiomass(i,j)

                            BiomassRemovedInDormancy    = TotalPlantBiomass * BiomassFracRemovedInDormancy

                            BiomassRemovedInHarvest = 0.0
                            if (Me%ComputeOptions%HarvestKill .and. Me%HarvestOnlyOccurred(i,j)) then
                                BiomassRemovedInHarvest = Me%Fluxes%BiomassRemovedInHarvest(i,j)
                            endif
                        
                            !Test biomass to adapt biomass removed
                            PredictedBiomass = TotalPlantBiomass - BiomassRemovedInHarvest - BiomassRemovedInDormancy 
                        
                            if(PredictedBiomass .lt. 0.0) then
                                BiomassRemovedInDormancy = TotalPlantBiomass - BiomassRemovedInHarvest
                            endif

                            Me%Fluxes%BiomassRemovedInDormancy(i,j)     = BiomassRemovedInDormancy
                            !Duplication for interface

                            if (Me%ComputeOptions%ModelNitrogen) then
                                Me%Fluxes%NitrogenRemovedInDormancy(i,j)    = Me%PlantNitrogenFraction(i,j) *     &
                                                                              Me%Fluxes%BiomassRemovedInDormancy(i,j)

                            endif
                            if (Me%ComputeOptions%ModelPhosphorus) then
                                Me%Fluxes%PhosphorusRemovedInDormancy(i,j)  = Me%PlantPhosphorusFraction(i,j) *   &
                                                                              Me%Fluxes%BiomassRemovedInDormancy(i,j)
                            endif


                            !Flux for HU Accumulated
!                            Me%Fluxes%BiomassDormancyFraction(i,j)      = BiomassRemovedInDormancy / TotalPlantBiomass
                            
                            !update HUAccumulated
                            WarningString = 'Dormancy'
                            call UpdatePlantGrowingStage(i,j, WarningString)
                       
    !                        Me%Growth%WaterStress(i,j)       = 1.0
                            Me%Growth%TemperatureStress(i,j) = 1.0
                            Me%Growth%NitrogenStress(i,j)    = 1.0
                            Me%Growth%PhosphorusStress(i,j)  = 1.0

                        !cool season plants (2,5)
                        case (2,5)
                        
    !                        Me%Growth%WaterStress(i,j)       = 1.0
                            Me%Growth%TemperatureStress(i,j) = 1.0
                            Me%Growth%NitrogenStress(i,j)    = 1.0
                            Me%Growth%PhosphorusStress(i,j)  = 1.0


                    end select

                endif
            endif

        enddo do2
        enddo do1

        !Duplication to compute interfaces with soil
        Me%Fluxes%ToSoil%DormancyBiomassToSoil    => Me%Fluxes%BiomassRemovedInDormancy
        if (Me%ComputeOptions%ModelNitrogen) then
            Me%Fluxes%ToSoil%DormancyNitrogenToSoil   => Me%Fluxes%NitrogenRemovedInDormancy
        endif
        if (Me%ComputeOptions%ModelNitrogen) then
            Me%Fluxes%ToSoil%DormancyPhosphorusToSoil => Me%Fluxes%PhosphorusRemovedInDormancy
        endif
                
   
    end subroutine DormancyFluxes
   
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
        real                                            :: GrazingStartPlantHU
!        integer                                         :: Op
        !Begin-----------------------------------------------------------------
        
!        if (.not. Me%GrazingFinished(i,j)) then
            
!            Op = Me%GrazingOperations(i,j)
            

            !! if cell currently not grazed, check to see if it is time
            !! to initialize grazing
            if (.not.Me%IsPlantBeingGrazed(i,j)) then

                call JulianDay(Me%ExternalVar%Now, JulDay)
                JulDay_Old = Me%ExternalVar%JulianDay_Old

                GrazingStartJulianDate = Me%VegetationTypes(Me%VegetationID(i,j))%GrazingDatabase%GrazingStartJulianDay 
                GrazingStartPlantHU    = Me%VegetationTypes(Me%VegetationID(i,j))%GrazingDatabase%GrazingStartPlantHU
                DT_day                 = Me%ComputeOptions%VegetationDT / (60 * 60 * 24)
                HUAcc                  = Me%HeatUnits%PlantHUAccumulated(i,j)
                HUAcc_Old              = Me%HeatUnits%PlantHUAccumulated_Old(i,j)

                !Check if grazing start is scheduled for day 
                if (GrazingStartJulianDate .gt. 0) then
                
                    if(JulDay .ge. GrazingStartJulianDate .and. JulDay_Old .lt. GrazingStartJulianDate) then
                
                        Me%IsPlantBeingGrazed(i,j) = .true.
                        Me%DaysOfGrazing(i,j) = 1
                
                    endif

                else if (GrazingStartPlantHU .gt. 0.) then
                
                    if(HUAcc .ge. GrazingStartPlantHU .and. HUAcc_Old .lt. GrazingStartPlantHU) then
                
                        Me%IsPlantBeingGrazed(i,j) = .true.
                        Me%DaysOfGrazing(i,j) = 1
                
                    endif
!                else
!                    Me%GrazingFinished(i,j) = .true. 

                endif
                
            else

                !! if not first day of grazing increment total days of grazing by one
                DT_day                 = Me%ComputeOptions%VegetationDT / (60 * 60 * 24)
                Me%DaysOfGrazing(i,j)  = Me%DaysOfGrazing(i,j) + DT_day
                GrazingDays            = Me%VegetationTypes(Me%VegetationID(i,j))%GrazingDatabase%GrazingDays 
                !! check to set if grazing period is over
                if (Me%DaysOfGrazing(i,j) .ge. GrazingDays) then
    
                    Me%IsPlantBeingGrazed(i,j) = .false.
    
                    Me%DaysOfGrazing(i,j) = 0.

!                    Me%GrazingOperations(i,j) = Me%GrazingOperations(i,j) + 1

!                    if (Me%GrazingOperations(i,j) .gt. size(GrazingStartPlantHU)) then
!                        Me%GrazingFinished(i,j) = .true. 
!                    endif

                end if

            end if

!        endif

!        deallocate(Me%VegetationTypes(Me%VegetationID(i,j))%HarvestKillDatabase%GrazingStartPlantHU)
    
    end subroutine CheckPlantGrazing

    !--------------------------------------------------------------------------

    subroutine GrazingFluxes

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                   :: i, j
        real                                      :: TotalPlantBiomass, PredictedBiomass
        real                                      :: BiomassRemovedInDormancy, BiomassRemovedInHarvest
        real                                      :: TotalPlantNitrogen, PredictedNitrogen
        real                                      :: NitrogenRemovedInDormancy, NitrogenRemovedInHarvest
        real                                      :: TotalPlantPhosphorus, PredictedPhosphorus
        real                                      :: PhosphorusRemovedInDormancy, PhosphorusRemovedInHarvest
        real                                      :: GrazingMinimumBiomass
        real                                      :: GrazingBiomass
        real                                      :: TramplingBiomass
        real                                      :: BiomassGrazed
        real                                      :: NitrogenGrazed
        real                                      :: PhosphorusGrazed
        real                                      :: BiomassTrampled
        real                                      :: NitrogenTrampled
        real                                      :: PhosphorusTrampled
        real                                      :: BiomassGrazedFraction
        character (Len = StringLength)            :: WarningString
        logical                                   :: PlantKilled
        
        !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                        
                PlantKilled = .false.
                if (Me%ComputeOptions%HarvestKill) then
                    if (Me%HarvestKillOccurred(i,j) .or. Me%KillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif
                
                if (Me%IsPlantGrowing(i,j) .and. Me%IsPlantBeingGrazed(i,j)      &
                    .and. .not. PlantKilled) then     

                    GrazingMinimumBiomass = Me%VegetationTypes(Me%VegetationID(i,j))%GrazingDatabase%GrazingMinimumBiomass
                    TotalPlantBiomass     = Me%StateVariables%TotalPlantBiomass(i,j)
                
                    BiomassRemovedInDormancy = 0.0
                    if (Me%ComputeOptions%Dormancy) then
                        if(Me%PlantGoingDormant(i,j)) then
                            BiomassRemovedInDormancy = Me%Fluxes%BiomassRemovedInDormancy(i,j)
                        endif
                    endif
                    BiomassRemovedInHarvest = 0.0
                    if (Me%ComputeOptions%HarvestKill) then
                        if(Me%HarvestOnlyOccurred(i,j)) then
                            BiomassRemovedInHarvest = Me%Fluxes%BiomassRemovedInHarvest(i,j)
                        endif
                    endif
                    PredictedBiomass = TotalPlantBiomass - BiomassRemovedInDormancy - BiomassRemovedInHarvest
                 
                    !! graze only if adequate biomass in HRU
                    if (PredictedBiomass .gt. GrazingMinimumBiomass) then

                        !Kg/ha
                        GrazingBiomass        = Me%VegetationTypes(Me%VegetationID(i,j))%GrazingDatabase%GrazingBiomass
                        TramplingBiomass      = Me%VegetationTypes(Me%VegetationID(i,j))%GrazingDatabase%TramplingBiomass
            
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
                        
                            TotalPlantNitrogen    = Me%StateVariables%TotalPlantNitrogen(i,j)

                            NitrogenRemovedInDormancy = 0.0
                            if (Me%ComputeOptions%Dormancy .and. Me%PlantGoingDormant(i,j)) then
                                NitrogenRemovedInDormancy = Me%Fluxes%NitrogenRemovedInDormancy(i,j)
                            endif
                            BiomassRemovedInHarvest = 0.0
                            if (Me%ComputeOptions%HarvestKill .and. Me%HarvestOnlyOccurred(i,j)) then
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

                            TotalPlantPhosphorus  = Me%StateVariables%TotalPlantNitrogen(i,j)

                            PhosphorusRemovedInDormancy = 0.0
                            if (Me%ComputeOptions%Dormancy) then
                                if(Me%PlantGoingDormant(i,j)) then
                                    PhosphorusRemovedInDormancy = Me%Fluxes%PhosphorusRemovedInDormancy(i,j)
                                endif
                            endif
                            PhosphorusRemovedInHarvest = 0.0
                            if (Me%ComputeOptions%HarvestKill) then
                                if(Me%HarvestOnlyOccurred(i,j)) then
                                    PhosphorusRemovedInHarvest = Me%Fluxes%PhosphorusRemovedInHarvest(i,j)
                                endif
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

            endif
        
        enddo do2
        enddo do1


    end subroutine GrazingFluxes

    !--------------------------------------------------------------------------

    subroutine CheckPlantPesticide(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: PotentialHU, PotentialHU_Old
        real                                            :: PesticideAppJulianDay, PesticideAppHU
        integer                                         :: JulDay, JulDay_Old
        integer                                         :: PesticideAppID, PesticideApps
        integer                                         :: PestApp !, UniquePest
        logical                                         :: PesticideApplied, PesticideContinuous
        logical                                         :: PlantBeingContPest
        integer                                         :: AccPesticideDays, TotalPesticideDays
        real                                            :: DT_day
        !Begin-----------------------------------------------------------------

        PotentialHU          = Me%HeatUnits%PotentialHUBase(i,j)
        PotentialHU_Old      = Me%HeatUnits%PotentialHUBase_Old(i,j)
        PesticideApps        = Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%NumberPesticideApps
        
if2:    if (PesticideApps .gt. 0) then
            
            !For every pesticide application
do3:        do PestApp = 1, PesticideApps
                
                PesticideAppID    = Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideID
                PesticideAppJulianDay   = &
                Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideAppJDay
                PesticideAppHU          = &
                 Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideAppHU
                
                
                PesticideContinuous    =  &
                Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%ContPesticideON
                PlantBeingContPest      = &
                Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%ContPesticideActive  
                AccPesticideDays = 0
                
                !if not already continuous pesticide check if it will start a unique or continuous pesticide application
                if (.not. PlantBeingContPest) then
               
                    call JulianDay(Me%ExternalVar%Now, JulDay)
                    JulDay_Old = Me%ExternalVar%JulianDay_Old
                    
                    PesticideApplied    = .false.
                    
                    !! pesticide application check
    if3:            if (PesticideAppJulianDay .gt. 0.0) then
    if4:                if((JulDay .ge. PesticideAppJulianDay) .and. (JulDay_Old .lt. PesticideAppJulianDay)) then
                            
                            PesticideApplied    = .true.

                            if (PesticideContinuous) then
                                PlantBeingContPest   = .true.
                                AccPesticideDays     = AccPesticideDays + 1
                            endif
                            
    !                        !correspondent in global pesticide list
    !do4:                    do UniquePest = 1, Me%Fluxes%Pesticides%UniquePesticides
    !if5:                        if (PesticideAppID == Me%Fluxes%Pesticides%Application(UniquePest)%PesticideID) then
    !                                Me%Fluxes%Pesticides%Application(UniquePest)%PesticideAppOccurred(i,j) = .true.
    !                                exit do4
    !                            endif if5
    !                        enddo do4
                        
                        endif if4
                    elseif (PesticideAppHU .gt. 0.0) then
    if6:                if((PotentialHU .ge. PesticideAppHU) .and. (PotentialHU_Old .lt. PesticideAppHU)) then
                            
                            PesticideApplied    = .true.

                            if (PesticideContinuous) then
                                PlantBeingContPest   = .true.
                                AccPesticideDays     = AccPesticideDays + 1
                            endif
                            
    !                        !correspondent in global pesticide list
    !do5:                    do UniquePest = 1, Me%Fluxes%Pesticides%UniquePesticides
    !if7:                        if (PesticideAppID == Me%Fluxes%Pesticides%Application(UniquePest)%PesticideID) then
    !                                Me%Fluxes%Pesticides%Application(UniquePest)%PesticideAppOccurred(i,j) = .true.
    !                                exit do5
    !                            endif if7
    !                        enddo do5                        
                        
                        endif if6
                    endif if3
                    
                else
                ! continue to apply the daily amount and check if is over

                    AccPesticideDays      = &
                    Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%ContPesticideAccDays
                    DT_day                = Me%ComputeOptions%VegetationDT / (60 * 60 * 24)
                    AccPesticideDays      = AccPesticideDays + DT_day
                    TotalPesticideDays    = &
                    Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%ContPesticideDays 
                    
                    PesticideApplied      = .true.
                    PlantBeingContPest    = .true.
                                   
                    if (AccPesticideDays .gt. TotalPesticideDays) then
                        PesticideApplied   = .false.
                        PlantBeingContPest = .false.
                        AccPesticideDays   = 0
                    endif                
                
                endif
                
                !update variables - too big to use inside do loops
                Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideAppOccurred = &
                 PesticideApplied
                Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%ContPesticideActive  = &
                 PlantBeingContPest
                Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%ContPesticideAccDays = &
                 AccPesticideDays
                
            enddo do3
        endif if2
   
    end subroutine CheckPlantPesticide
    
    !--------------------------------------------------------------------------

    subroutine CheckPlantFertilization(i,j)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        !Local-----------------------------------------------------------------
        real                                            :: PotentialHU, PotentialHU_Old
        real                                            :: FertilizerAppJulianDay, FertilizerAppHU
        integer                                         :: JulDay, JulDay_Old
        integer                                         :: FertilizerApps
        integer                                         :: FertApp
        logical                                         :: FertilizationContinuous, PlantBeingContFert
        logical                                         :: FertilizerApplied
        integer                                         :: AccFertilizationDays, TotalFertilizationDays
        real                                            :: DT_day
        !Begin-----------------------------------------------------------------

        PotentialHU          = Me%HeatUnits%PotentialHUBase(i,j)
        PotentialHU_Old      = Me%HeatUnits%PotentialHUBase_Old(i,j)
        FertilizerApps        = Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%NumberFertilizerApps
        
if2:    if (FertilizerApps .gt. 0) then
            
            !For every pesticide application
do3:        do FertApp = 1, FertilizerApps
                
                FertilizerAppJulianDay   = &
                 Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppJDay
                FertilizerAppHU          = &
                 Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppHU
                FertilizationContinuous  = &
                 Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationON
                
                !if plant is under continuous fertilization
                PlantBeingContFert       = &
                 Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationActive
                
                FertilizerApplied    = .false.
                AccFertilizationDays = 0
                                
                !if not check if it will start a unique or continuous fertilization
                if (.not. PlantBeingContFert) then
                
                    call JulianDay(Me%ExternalVar%Now, JulDay)
                    JulDay_Old = Me%ExternalVar%JulianDay_Old
                    
                    !! pesticide application check
if3:                if (FertilizerAppJulianDay .gt. 0.0) then
if4:                    if(JulDay .ge. FertilizerAppJulianDay .and. JulDay_Old .lt. FertilizerAppJulianDay) then
                            
                            FertilizerApplied = .true.
                            
                            if (FertilizationContinuous) then
                                PlantBeingContFert   = .true.
                                AccFertilizationDays = AccFertilizationDays + 1
                            endif
                        
                        endif if4
                        
                    elseif (FertilizerAppHU .gt. 0.0) then
if6:                    if(PotentialHU .ge. FertilizerAppHU .and. PotentialHU_Old .lt. FertilizerAppHU) then
                            
                            FertilizerApplied = .true.

                            if (FertilizationContinuous) then
                                PlantBeingContFert   = .true.
                                AccFertilizationDays = AccFertilizationDays + 1
                            endif                        
                        
                        endif if6
                    
                    endif if3
               
                ! continue to apply the daily amount and check if is over
                else
                    AccFertilizationDays      = &
                    Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationAccDays
                    DT_day                    = Me%ComputeOptions%VegetationDT / (60 * 60 * 24)
                    AccFertilizationDays      = AccFertilizationDays + DT_day
                    TotalFertilizationDays    = &
                    Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationDays 
                    
                    FertilizerApplied         = .true.
                    PlantBeingContFert        = .true.
                                   
                    if (AccFertilizationDays .gt. TotalFertilizationDays) then
                        FertilizerApplied     = .false.
                        PlantBeingContFert    = .false.
                        AccFertilizationDays  = 0
                    endif
                endif
                
                !update variables - too big to use inside do loops
                Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppOccurred   = &
                 FertilizerApplied
                Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationActive = &
                 PlantBeingContFert
                Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%ContFertilizationAccDays = &
                 AccFertilizationDays  
                
            enddo do3
        endif if2
   
    end subroutine CheckPlantFertilization
    
    !--------------------------------------------------------------------------

    subroutine UpdatePlantGrowingStage (i,j, WarningString )
       
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: i, j
        character (Len = StringLength), intent(IN)      :: WarningString
        !Local-----------------------------------------------------------------
        real                                            :: HUAcc, HarvestedFraction 
        real                                            :: GrazedFraction!, DormancyFraction
        !Begin-----------------------------------------------------------------
        

        if (WarningString .eq. 'Planting') then
            
            Me%HeatUnits%PlantHUAccumulated_Old(i,j) = 0.0            
            Me%HeatUnits%PlantHUAccumulated    (i,j) = 0.0
        
        elseif (WarningString .eq. 'Kill') then
        
            Me%HeatUnits%PlantHUAccumulated_Old(i,j) = Me%HeatUnits%PlantHUAccumulated(i,j)            
            Me%HeatUnits%PlantHUAccumulated    (i,j) = 0.0

        elseif (WarningString .eq. 'Harvest') then

            HUAcc = Me%HeatUnits%PlantHUAccumulated(i,j)

            HarvestedFraction = Me%Fluxes%BiomassHarvestedFraction(i,j)
            HUAcc = HUAcc - (HUAcc * HarvestedFraction)
            Me%HeatUnits%PlantHUAccumulated_Old(i,j)  = Me%HeatUnits%PlantHUAccumulated(i,j) 
            Me%HeatUnits%PlantHUAccumulated(i,j) = max(HUAcc, 0.0)
        
        elseif (WarningString .eq. 'Dormancy') then
            
            !Changed from SWAT. Not zero after dormancy but fraction
            HUAcc = Me%HeatUnits%PlantHUAccumulated(i,j)
            
!            DormancyFraction = Me%Fluxes%BiomassDormancyFraction(i,j)
!            HUAcc = HUAcc - (HUAcc * DormancyFraction)
            Me%HeatUnits%PlantHUAccumulated_Old(i,j) = Me%HeatUnits%PlantHUAccumulated(i,j)            
!            Me%HeatUnits%PlantHUAccumulated    (i,j) = max(HUAcc, 0.0)
            Me%HeatUnits%PlantHUAccumulated    (i,j) = 0.0

        
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
        integer                                         :: i, j, PlantType, TreeYearsToMaturity
        logical                                         :: PlantKilled
        real                                            :: PlantBiomass, BiomassGrowth, BiomassGrazed
        real                                            :: BiomassRemovedInHarvest, BiomassRemovedInDormancy
        real                                            :: PlantNitrogen, NitrogenUptake, NitrogenGrazed
        real                                            :: NitrogenRemovedInHarvest, NitrogenRemovedInDormancy
        real                                            :: PlantPhosphorus, PhosphorusUptake, PhosphorusGrazed
        real                                            :: PhosphorusRemovedInHarvest, PhosphorusRemovedInDormancy
        !Begin------------------------------------------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint .and. Me%IsPlantGrowing(i,j)) then
                     
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                     
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
                if (Me%ComputeOptions%HarvestKill) then
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
                    PlantBiomass      = Me%StateVariables%TotalPlantBiomass(i,j)
                    BiomassGrowth     = Me%Fluxes%BiomassGrowth(i,j)
                endif

                if (Me%ComputeOptions%ModelNitrogen) then
                    PlantNitrogen     = Me%StateVariables%TotalPlantNitrogen(i,j)
                    NitrogenUptake    = Me%Fluxes%NitrogenUptake(i,j)
                endif
                
!                write(*,*), Me%ComputeOptions%ModelPhosphorus
                if (Me%ComputeOptions%ModelPhosphorus) then
                    PlantPhosphorus   = Me%StateVariables%TotalPlantPhosphorus(i,j)
                    PhosphorusUptake  = Me%Fluxes%PhosphorusUptake(i,j)
                endif
            
                PlantKilled = .false.
                if (Me%ComputeOptions%HarvestKill) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif
                
                
                if (PlantKilled .or. Me%PlantingOccurred(i,j)) then
        
                    PlantBiomass    = 0.0
                    PlantNitrogen   = 0.0
                    PlantPhosphorus = 0.0

                else

                    PlantBiomass    = PlantBiomass + BiomassGrowth - BiomassGrazed - BiomassRemovedInHarvest             &
                                      - BiomassRemovedInDormancy
                
                    PlantNitrogen   = PlantNitrogen + NitrogenUptake - NitrogenGrazed - NitrogenRemovedInHarvest         &
                                      - NitrogenRemovedInDormancy

                    PlantPhosphorus = PlantPhosphorus + PhosphorusUptake - PhosphorusGrazed - PhosphorusRemovedInHarvest &
                                      - PhosphorusRemovedInDormancy
                
  
                endif

                PlantType            = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                TreeYearsToMaturity  = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%TreeYearsToMaturity

                !For trees, growth is slower, getting maximum biomass not in a single year. As so growth may be limited.
                if (PlantType == Trees .and. TreeYearsToMaturity .gt. 0 ) then
                
                    PlantBiomass = min (PlantBiomass, Me%Growth%TreeMaximumAnnualBiomass(i,j))

                end if

                if (Me%ComputeOptions%ModelPlantBiomass) then
                    Me%StateVariables%TotalPlantBiomass(i,j)    = PlantBiomass
                endif
                if (Me%ComputeOptions%ModelNitrogen) then
                    Me%StateVariables%TotalPlantNitrogen(i,j)   = PlantNitrogen
                endif
                if (Me%ComputeOptions%ModelPhosphorus) then
                    Me%StateVariables%TotalPlantPhosphorus(i,j) = PlantPhosphorus
                endif

            endif


        enddo do2
        enddo do1

    end subroutine UpdateGlobalPlantProperties

    !--------------------------------------------------------------------------

    subroutine UpdateRootProperties

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------        
        integer                                         :: i, j, kTop, KFloor, STAT_CALL
        logical                                         :: PlantKilled
        real                                            :: MaxRootDepth, RootBiomassFraction, SoilDepth
        integer                                         :: VegetationID, PlantType
        
        !Begin-----------------------------------------------------------------

        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  SZZ         = Me%ExternalVar%SZZ,                     &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR01'

        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z    = Me%ExternalVar%KFloor,                            &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR02'


do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint .and. Me%IsPlantGrowing(i,j)) then
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
            
                !To send to PMP to source/sinks fluxes computed in vegetation
                !Specially needed in the case of kill where it would not compute fluxes in PMP
                Me%RootDepthOld(i,j) = Me%StateVariables%RootDepth(i,j)
                
                PlantKilled = .false.
                if (Me%ComputeOptions%HarvestKill) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif


                if (PlantKilled) then
        
                    if(Me%ComputeOptions%ModelRootBiomass) then
                        Me%StateVariables%RootBiomass(i,j) = 0.0
                    endif
                    Me%StateVariables%RootDepth(i,j)   = 0.0

                else

                    !!Root Biomass
                    if (Me%ComputeOptions%ModelRootBiomass) then
                        
                        RootBiomassFraction = 0.4 - 0.2 * Me%HeatUnits%PlantHUAccumulated(i,j)
    
                        if (RootBiomassFraction .lt. 0.0) then
                            RootBiomassFraction = 0.0
                        endif
                        Me%StateVariables%RootBiomass(i,j) = RootBiomassFraction * Me%StateVariables%TotalPlantBiomass(i,j)

                    endif

                    !!Root Depth
                    VegetationID = Me%VegetationID(i,j)
                    PlantType    = Me%VegetationTypes(VegetationID)%GrowthDatabase%PlantType
                    MaxRootDepth = Me%VegetationTypes(VegetationID)%GrowthDatabase%MaximumRootDepth
                    kTop         = Me%WorkSize%KUB
                    kFloor       = Me%ExternalVar%KFloor(i,j)
                    !SZZ has the altitude (negative) of top face. For the Kfloor - 1 is also defined
                    SoilDepth    = (-1. * Me%ExternalVar%SZZ(i,j,kTop)) - (-1.* Me%ExternalVar%SZZ(i,j,kFloor - 1))
                    if(MaxRootDepth .gt. SoilDepth) then
                        MaxRootDepth = SoilDepth
                    endif
                    
                    select case (PlantType)
                        case (1, 2, 4, 5)
            
                            Me%StateVariables%RootDepth(i,j) = 2.5 * Me%HeatUnits%PlantHUAccumulated (i,j) * MaxRootDepth
            
                            if (Me%StateVariables%RootDepth(i,j) .gt. MaxRootDepth) then
                                Me%StateVariables%RootDepth(i,j) = MaxRootDepth
                            endif

                            if (Me%StateVariables%RootDepth(i,j) .lt. 0.01) then
                                Me%StateVariables%RootDepth(i,j) = 0.01
                            endif
                        case default
            
                            Me%StateVariables%RootDepth(i,j) = MaxRootDepth

                    end select
            
                endif
        
            endif

        enddo do2
        enddo do1 

        call UnGetGeometry( Me%ObjGeometry, Me%ExternalVar%KFloor,       STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR03'
        
        call UnGetGeometry( Me%ObjGeometry, Me%ExternalVar%SZZ,          STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'CheckRootDepth - ModuleVegetation - ERR04'

    end subroutine UpdateRootProperties

    !--------------------------------------------------------------------------

    subroutine UpdateLeafProperties

        !Arguments-------------------------------------------------------------
        
        !Local----------------------------------------------------------        
        integer                                         :: i, j, STAT_CALL
        logical                                         :: PlantKilled, PlantGoingDormant
        real                                            :: LAI, BiomassGrazedFraction, BiomassHarvestedFraction
        real                                            :: LAIChange
        real                                            :: LAIMax, LAIMinDormant
        integer                                         :: PlantType  
        logical                                         :: UseBoundaryLAI      
        type(T_Property), pointer                       :: BoundaryLAI
        integer                                         :: vegID
        type(T_VegetationType)                          :: veg
        !Begin----------------------------------------------------------
                
        UseBoundaryLAI = .false.
                
        if (Me%UseBoundaryLAI) then
            call SearchProperty(BoundaryLAI, BoundaryLeafAreaIndex_, .false., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'UpdateLeafProperties - ModuleVegetation - ERR001'
            endif
        else
            BoundaryLAI => null()
        endif
        
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint .and. Me%IsPlantGrowing(i,j) ) then
                   
                !Plant type 0 is NO PLANT, so no operations should happen
                vegID = Me%VegetationID(i,j)
                veg   = Me%VegetationTypes(vegID)
                 
                if (veg%GrowthDatabase%PlantType == NotAPlant) cycle
                   
                LAIMax = veg%GrowthDatabase%LAIMax 
                if ((Me%UseBoundaryLAI)) then
                    if (BoundaryLAI%Field(i,j) >= 0) then
                        if (BoundaryLAI%Field(i,j) > LAIMax) then
                            LAIMax = BoundaryLAI%Field(i,j)
                        endif
                        
                        UseBoundaryLAI = .true.
                    else
                        UseBoundaryLAI = .false.
                    endif
                endif
                   
                PlantType = veg%GrowthDatabase%PlantType

                if (veg%HasLeaves) then

                    PlantKilled = .false.
                    if (Me%ComputeOptions%HarvestKill) then
                        if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                            PlantKilled = .true.
                        endif
                    endif
                    
                    PlantGoingDormant = .false.
                    if (Me%ComputeOptions%Dormancy) then
                        if (Me%PlantGoingDormant(i,j)) then
                            PlantGoingDormant = .true.
                        endif 
                    endif                   
                
                    if (PlantKilled) then

                        if (UseBoundaryLAI) then
                            Me%StateVariables%LeafAreaIndex(i,j) = Me%StateVariables%LeafAreaIndex(i,j) + LAIChange  
                        else
                            Me%StateVariables%LeafAreaIndex(i,j) = 0.0
                        endif

                    elseif (PlantGoingDormant) then
                    
                        LAIMinDormant = veg%GrowthDatabase%LAIMinDormant
                        
                        if (UseBoundaryLAI) then
                            Me%StateVariables%LeafAreaIndex(i,j) = Me%StateVariables%LeafAreaIndex(i,j) + LAIChange  
                        else                        
                            if (Me%StateVariables%LeafAreaIndex(i,j) .gt. LAIMinDormant) then
                                Me%StateVariables%LeafAreaIndex(i,j) = LAIMinDormant
                            endif
                        endif
                                        
                    else

                        LAI = Me%StateVariables%LeafAreaIndex(i,j)

                        BiomassGrazedFraction = 0.0
                        if (Me%ComputeOptions%Grazing) then
                            if (Me%IsPlantBeingGrazed(i,j)) then
                                BiomassGrazedFraction = Me%Fluxes%BiomassGrazedFraction(i,j)
                            endif
                        endif

                        BiomassHarvestedFraction = 0.0
                        if (Me%ComputeOptions%HarvestKill) then
                            if (Me%HarvestOnlyOccurred(i,j)) then
                                BiomassHarvestedFraction = Me%Fluxes%BiomassHarvestedFraction(i,j)
                            endif
                        endif


                        if (.not. Me%LAISenescence(i,j)) then
            
                            Me%StateVariables%LeafAreaIndex(i,j) = (LAI + Me%Fluxes%LAIChange(i,j)) - &
                                                                   (LAI * BiomassGrazedFraction) -    &
                                                                   (LAI * BiomassHarvestedFraction)
                
                            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                                Me%LAIBeforeSenescence(i,j) = Me%StateVariables%LeafAreaIndex(i,j)                                
                            endif

                        elseif (LAI .ne. 0.0) then                           

                            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                                !LAI Computed from maximum value reached and not from last (SWAT theory). If grazing or 
                                !harvesting occurr during senescence LAI is increased due to the computation formulation.
                                Me%StateVariables%LeafAreaIndex(i,j) = Me%LAIBeforeSenescence(i,j)
                            else
                                !LAI computed from last value. It avoids erorrs when grazing or harvesting occurrs 
                                !during senescence.
                                Me%StateVariables%LeafAreaIndex(i,j) = LAI
                            endif 
                               
                            Me%StateVariables%LeafAreaIndex(i,j) = (Me%StateVariables%LeafAreaIndex(i,j) - LAIChange) - &
                                                                   (LAI * BiomassGrazedFraction) -                      &
                                                                   (LAI * BiomassHarvestedFraction)                                                                      
                            
                        else if (UseBoundaryLAI) then
                        
                                Me%StateVariables%LeafAreaIndex(i,j) = Me%StateVariables%LeafAreaIndex(i,j) + LAIChange                              
                            
                        endif
        
                        if ((PlantType .eq. Trees) .and. (Me%UseBoundaryLAI .EQV. .false.)) then
                        
                            LAIMax = Me%Growth%TreeFractionToMaturity(i,j) * LAIMax
                        
                        endif
        
                    endif                                       

                else
        
                    if (UseBoundaryLAI) then
                        Me%StateVariables%LeafAreaIndex(i,j) = Me%StateVariables%LeafAreaIndex(i,j) + LAIChange  
                    else
                        Me%StateVariables%LeafAreaIndex(i,j) = 0.0
                    endif
                    
                endif
                
                if (Me%StateVariables%LeafAreaIndex(i,j) .gt. LAIMax) then
                    Me%StateVariables%LeafAreaIndex(i,j) = LAIMax
                endif
            
                !If in the future a leaf flux is computed explicitly than LAI can not be
                !taken to zero here. If mass flux .gt. LAI then flux is LAI
                if (Me%StateVariables%LeafAreaIndex(i,j) .lt. 0.0) then   
                    Me%StateVariables%LeafAreaIndex(i,j) = 0.0
                endif                

            endif
            
        enddo do2
        enddo do1

    end subroutine UpdateLeafProperties

    !--------------------------------------------------------------------------

    subroutine UpdateStemProperties

        !Arguments-------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        logical                                         :: PlantKilled
        real                                            :: MaxCanopyHeight, CanopyHeight, NearMaximum
        integer                                         :: PlantType
        real                                            :: BiomassHarvestedFraction
        !Begin-----------------------------------------------------------------
       
do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint .and. Me%IsPlantGrowing(i,j)) then

                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                Me%ChangeCanopyEnabled(i,j) = .false.
                PlantKilled = .false.
                if (Me%ComputeOptions%HarvestKill) then
                    if (Me%KillOccurred(i,j) .or. Me%HarvestKillOccurred(i,j)) then
                        PlantKilled = .true.
                    endif
                endif

                if (Me%ComputeOptions%HarvestKill .and. PlantKilled) then

                    Me%StateVariables%CanopyHeight(i,j) = 0.0

                else

                    PlantType         = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                    MaxCanopyHeight   = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%MaxCanopyHeight

                    !!Compute new canopy height
                    if (PlantType == Trees) then
            
                        Me%StateVariables%CanopyHeight(i,j) = MaxCanopyHeight * Me%Growth%TreeFractionToMaturity(i,j)
        
                    else

                        NearMaximum = MaxCanopyHeight - 0.01
                
                        !Change canopy was enabled because after reaching its maximum, crop canopy height was not affected  
                        !by harvest(or grazing but grazing can be consider only to remove leafs and not reduce stem height)
                        Me%ChangeCanopyEnabled(i,j) = .false.
                        if (Me%StateVariables%CanopyHeight(i,j) .gt. NearMaximum .and. &
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
                            if (Me%ComputeOptions%HarvestKill) then
                                if (Me%HarvestOnlyOccurred(i,j)) then
                                    BiomassHarvestedFraction = Me%Fluxes%BiomassHarvestedFraction(i,j)
                                endif
                            endif
            
                            CanopyHeight = Me%StateVariables%CanopyHeight(i,j)
            
                            Me%StateVariables%CanopyHeight(i,j) = CanopyHeight                                      &
!                                                                    - CanopyHeight * BiomassGrazedFraction         &
                                                                    - CanopyHeight * BiomassHarvestedFraction
                    
                            !If in the future a height reduction is computed explicitly than canopy height can not be
                            !taken to zero here. If mass flux .gt. canopy height then flux is canopy height  
                            if (Me%StateVariables%CanopyHeight(i,j) .lt. 0.0) then
    
                                Me%StateVariables%CanopyHeight(i,j) = 0.0

                            endif    
                
                        else                        
                

                            Me%StateVariables%CanopyHeight(i,j) = MaxCanopyHeight * sqrt(Me%PlantLAIMaxFraction(i,j))
            
                        endif 

                    endif

                endif


            endif

        enddo do2
        enddo do1


    end subroutine UpdateStemProperties

    !--------------------------------------------------------------------------

!    subroutine InterfaceWithSoil
!        !Local-----------------------------------------------------------------
!        !Begin-----------------------------------------------------------------
!        
!        !Fertilization
!        if(Me%ComputeOptions%Fertilization) then
!            call FertilizationFluxes
!        endif
!        
!        !
!        call PlantFluxesToSoil
!
!    end subroutine InterfaceWithSoil

    !--------------------------------------------------------------------------

    subroutine FertilizationFluxes
        !Local-----------------------------------------------------------------
        integer                                            :: PlantType, i,j
        integer                                            :: JulDay, JulDay_Old
        real                                               :: NTreshold, PTreshold
        logical                                            :: ExplicitPhosphorus, NeedNitrogen
        logical                                            :: NeedPhosphorus, ApplyPhosphorus, Dormant
        integer                                            :: NumberOfApplications, FertApp
        !Begin-----------------------------------------------------------------

        !Autofertilization
        if (Me%ComputeOptions%Autofertilization) then
       
do1:        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    
                if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then

                    !Plant type 0 is NO PLANT, so no operations should happen
                    if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle

                    NumberOfApplications     = Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%NumberFertilizerApps
                    if (NumberOfApplications .gt. 0) then
                        
                        PlantType              = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType
                        NTreshold              = &
                        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%NTreshold
                        PTreshold              = &
                        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%PTreshold
                        ExplicitPhosphorus     = &
                        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%ExplicitPhosphorus
                        NeedNitrogen           = .false.
                        NeedPhosphorus         = .false.
                        ApplyPhosphorus        = .false.                    
                        
                        call JulianDay(Me%ExternalVar%Now, JulDay)
                        JulDay_Old = Me%ExternalVar%JulianDay_Old

                        !Year Begining - nullify annual counter
                        if (JulDay_Old .gt. 364 .and. JulDay .ge. 1) then
                            if (Me%ComputeOptions%ModelNitrogen) then
                                Me%AnnualNitrogenFertilized(i,j) = 0.0
                            endif
                            if (Me%ComputeOptions%ModelPhosphorus) then
                                Me%AnnualPhosphorusFertilized(i,j) = 0.0
                            endif
                        endif           
                        
                        Dormant = .false.
                        if (Me%ComputeOptions%Dormancy) then
                            if (Me%IsPlantDormant(i,j)) then
                                Dormant = .true.
                            endif
                        endif
                        
                        if (Me%IsPlantGrowing(i,j) .and. .not. Dormant .and. Me%HeatUnits%PlantHUAccumulated(i,j).le. 1.0) then

                            !Test if plant needs nitrogen to apply in soil
                            if (Me%ComputeOptions%ModelNitrogen .and. (Me%Growth%NitrogenStress(i,j) .lt. NTreshold)) then
                                select case (PlantType)
                                    case (4,5,6,7)
                                    
                                        NeedNitrogen = .true.

                                        !Apply phosphorus contained in the fertilizer (even if not needed by the plant)
                                        !From SWAT - phosphorus is applied when plant goes nitrogen stress
                                        if (Me%ComputeOptions%ModelPhosphorus .and. .not. ExplicitPhosphorus) then
                                            ApplyPhosphorus = .true.
                                        endif
                                        
                                end select
                            endif
                        
                            !New option
                            !Only apply phosphorus if needed by the plant
                            if (Me%ComputeOptions%ModelPhosphorus .and. ExplicitPhosphorus) then
                                if (Me%Growth%PhosphorusStress(i,j) .lt. PTreshold) then
                                    select case (PlantType)
                                        case (4,5,6,7)
                                
                                            NeedPhosphorus = .true.

                                    end select
                                endif
                            endif

                            if (NeedNitrogen .or. ApplyPhosphorus .or. NeedPhosphorus) then
                            
                                call AutoFertilization(i,j, NeedNitrogen, ApplyPhosphorus, NeedPhosphorus)
                        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%FertilizerAppOccurred = .true.

                            endif

                        endif
                    endif
                endif
                
            enddo do2
            enddo do1
        
                
        !Scheduled fertilization with julian day or HU
        elseif(Me%ComputeOptions%FertilizationApp) then

do3:        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do4:        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                   
                    !Plant type 0 is NO PLANT, so no operations should happen
                    if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                   
                    NumberOfApplications = Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%NumberFertilizerApps
                    if (NumberOfApplications .gt. 0) then
                        !check if any application occured
                        do FertApp = 1, NumberOfApplications
                            if (&
                   Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppOccurred) then
                                
                                call ScheduledFertilization(i,j, FertApp)
                        
                            endif
                        enddo
                    endif
                endif
                
            enddo do4
            enddo do3
        
        endif
        
        !Duplication to compute interfaces with soil
        if (Me%ComputeOptions%ModelNitrogen) then
            Me%Fluxes%ToSoil%FertilNitrateToSoilSurface      => Me%Fluxes%FertilNitrateInSurface
            Me%Fluxes%ToSoil%FertilNitrateToSoilSubSurface   => Me%Fluxes%FertilNitrateInSubSurface
            Me%Fluxes%ToSoil%FertilAmmoniaToSoilSurface      => Me%Fluxes%FertilAmmoniaInSurface
            Me%Fluxes%ToSoil%FertilAmmoniaToSoilSubSurface   => Me%Fluxes%FertilAmmoniaInSubSurface
            Me%Fluxes%ToSoil%FertilOrganicNToSoilSurface     => Me%Fluxes%FertilOrganicNInSurface
            Me%Fluxes%ToSoil%FertilOrganicNParticToFluff     => Me%Fluxes%FertilOrganicNParticInFluff
            Me%Fluxes%ToSoil%FertilOrganicNToSoilSubSurface  => Me%Fluxes%FertilOrganicNInSubSurface
        endif

        if (Me%ComputeOptions%ModelPhosphorus) then
            Me%Fluxes%ToSoil%FertilOrganicPToSoilSurface     => Me%Fluxes%FertilOrganicPInSurface
            Me%Fluxes%ToSoil%FertilOrganicPToSoilSubSurface  => Me%Fluxes%FertilOrganicPInSubSurface
            Me%Fluxes%ToSoil%FertilMineralPToSoilSurface     => Me%Fluxes%FertilMineralPInSurface
            Me%Fluxes%ToSoil%FertilMineralPToSoilSubSurface  => Me%Fluxes%FertilMineralPInSubSurface   
            Me%Fluxes%ToSoil%FertilOrganicPParticToFluff     => Me%Fluxes%FertilOrganicPParticInFluff   
        endif

    end subroutine FertilizationFluxes

    !--------------------------------------------------------------------------

    subroutine AutoFertilization(i,j, NeedNitrogen, ApplyPhosphorus, NeedPhosphorus)
        !Argument--------------------------------------------------------------
        integer, intent(IN)                         :: i,j
        logical, intent(IN)                         :: NeedNitrogen, ApplyPhosphorus, NeedPhosphorus
        !Local-----------------------------------------------------------------
        real                                        :: NitrogenApplicationMax, NitrogenAnnualMax
        real                                        :: PhosphorusApplicationMax, PhosphorusAnnualMax
        real                                        :: NitrogenDemand, PhosphorusDemand, FertilizerAmount
        real                                        :: MineralNFracInFertilizer, HUAcc
        real                                        :: FertilizerFracApplyedInSurface, OrganicFracParticulate
        real                                        :: AmmoniaFracInMineralN, OrganicNFracInFertilizer
        real                                        :: MineralPFracInFertilizer, OrganicPFracInFertilizer
        real                                        :: MineralPFrac
        real                                        :: PotentialAnnualNitrogen, PotentialAnnualPhosphorus
        integer                                     :: NStress, KUB, KLB, k
        real                                        :: NitrateSoilAmount, WaterVolume, Area
        real                                        :: PlantNitrogenAmount
        !Begin-----------------------------------------------------------------

        !Compute fertilizer applyance
        if (NeedNitrogen) then
            
            NStress = Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%NStress
            
            if (NStress == NTarget_) then
                
                NitrogenApplicationMax = &
                Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%NitrogenApplicationMax 
                
                KUB = Me%WorkSize%KUB
                KLB = Me%WorkSize%KLB
                
                NitrateSoilAmount = 0.
                !how much nitrate in soil layers - as for previous step
                do k = KLB, KUB
                    
                    if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then
                        
                        !Plant type 0 is NO PLANT, so no operations should happen
                        if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
                        
                        !m3H20 = m3H20/m3cell * m3 cell
                        WaterVolume = Me%ExternalVar%SoilWaterContent(i,j,k) * Me%ExternalVar%CellVolume(i,j,k)
                        Area        = Me%ExternalVar%GridCellArea(i,j)
                        
                        !kgN/ha = gN/m3H20 * 1E-3kg/g * m3H20 / (m2) * 10000m2/ha 
                        NitrateSoilAmount = NitrateSoilAmount + Me%ExternalVar%SoilNitrate(i,j,k) * 1E-3               &
                                             * WaterVolume / (Area) * 10000.
                    endif
                    
                enddo
                PlantNitrogenAmount = Me%StateVariables%TotalPlantNitrogen(i,j)
                
                NitrogenDemand = Me%NitrogenYeldEstimate(i,j) - NitrateSoilAmount - PlantNitrogenAmount

                if (NitrogenDemand .gt. NitrogenApplicationMax) then
                    NitrogenDemand = NitrogenApplicationMax
                endif

                if (NitrogenDemand .lt. 0.) then
                    NitrogenDemand = 0.
                endif
                
            elseif (NStress == AnnualMaxApproach_) then
            
                !KgN/ha
                NitrogenApplicationMax = &
                Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%NitrogenApplicationMax 
                
                HUAcc                  = Me%HeatUnits%PlantHUAccumulated(i,j)
                
                !KgN/ha
                NitrogenDemand = NitrogenApplicationMax * (1. - HUAcc)
                if (NitrogenDemand .gt. NitrogenApplicationMax) then
                    NitrogenDemand = NitrogenApplicationMax
                endif
            
            endif

            NitrogenAnnualMax      = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%NitrogenAnnualMax
        
            !KgN/ha
            PotentialAnnualNitrogen = Me%AnnualNitrogenFertilized(i,j) + NitrogenDemand
            if (PotentialAnnualNitrogen .ge. NitrogenAnnualMax) then
                NitrogenDemand = NitrogenAnnualMax - (PotentialAnnualNitrogen - NitrogenDemand)
            endif

!!         if (targn <= 1.e-6) return
            
            !Kg/ha
            FertilizerAmount = 0.0
            MineralNFracInFertilizer = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%MineralNFracInFertilizer
            
            if (MineralNFracInFertilizer .gt. 0.0001) then
                Me%AnnualNitrogenFertilized(i,j) = Me%AnnualNitrogenFertilized(i,j) + NitrogenDemand
                FertilizerAmount = NitrogenDemand / MineralNFracInFertilizer
            else
                FertilizerAmount = 0.0
            endif
            
            FertilizerFracApplyedInSurface = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%FertilizerFracApplyedInSurface
            OrganicFracParticulate = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%OrganicFracParticulate
            AmmoniaFracInMineralN          = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%AmmoniaFracInMineralN
            OrganicNFracInFertilizer       = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%OrganicNFracInFertilizer
            MineralPFracInFertilizer       = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%MineralPFracInFertilizer
            OrganicPFracInFertilizer       = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%OrganicPFracInFertilizer
            
            !KgN/ha
            Me%Fluxes%FertilNitrateInSurface(i,j)     = FertilizerFracApplyedInSurface * FertilizerAmount                      &
                                                        * MineralNFracInFertilizer * (1. - AmmoniaFracInMineralN)
            Me%Fluxes%FertilNitrateInSubSurface(i,j)  = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount               &
                                                        * MineralNFracInFertilizer * (1. - AmmoniaFracInMineralN)
            Me%Fluxes%FertilAmmoniaInSurface(i,j)     = FertilizerFracApplyedInSurface * FertilizerAmount                      &
                                                        * MineralNFracInFertilizer * AmmoniaFracInMineralN
            Me%Fluxes%FertilAmmoniaInSubSurface(i,j)  = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount               &
                                                        * MineralNFracInFertilizer * AmmoniaFracInMineralN
            Me%Fluxes%FertilOrganicNParticInFluff(i,j) = FertilizerAmount * OrganicNFracInFertilizer * OrganicFracParticulate
                                                                    
            Me%Fluxes%FertilOrganicNInSurface(i,j)    = FertilizerFracApplyedInSurface * FertilizerAmount                      &
                                                        * OrganicNFracInFertilizer * (1. - OrganicFracParticulate)
            Me%Fluxes%FertilOrganicNInSubSurface(i,j) = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount               &
                                                        * OrganicNFracInFertilizer * (1. - OrganicFracParticulate)

            if (ApplyPhosphorus) then
                
                !KgP/ha
                Me%Fluxes%FertilOrganicPParticInFluff(i,j) = FertilizerAmount * OrganicPFracInFertilizer * OrganicFracParticulate
                
                Me%Fluxes%FertilOrganicPInSurface(i,j)    = FertilizerFracApplyedInSurface * FertilizerAmount          &
                                                            * OrganicPFracInFertilizer * (1. - OrganicFracParticulate)
                Me%Fluxes%FertilOrganicPInSubSurface(i,j) = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount   &
                                                            * OrganicPFracInFertilizer * (1. - OrganicFracParticulate)
                !P stress 
                MineralPFrac = 0.0
                if (Me%Growth%PhosphorusStress(i,j) .le. 0.75) then
                    MineralPFrac = MineralNFracInFertilizer / 7.
                else
                    MineralPFrac = MineralPFracInFertilizer
                endif
                
                !KgP/ha
                Me%Fluxes%FertilMineralPInSurface(i,j)    =  FertilizerFracApplyedInSurface * FertilizerAmount          &
                                                             * MineralPFrac 
                Me%Fluxes%FertilMineralPInSubSurface(i,j) =  (1. - FertilizerFracApplyedInSurface) * FertilizerAmount   &
                                                             * MineralPFrac 
            endif

        endif
        
        !only annual max approach
        if (NeedPhosphorus) then

            PhosphorusApplicationMax = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%PhosphorusApplicationMax 
            PhosphorusAnnualMax      = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Auto%PhosphorusAnnualMax
            HUAcc                    = Me%HeatUnits%PlantHUAccumulated(i,j)
            
            PhosphorusDemand = PhosphorusApplicationMax * (1. - HUAcc)
            if (PhosphorusDemand .gt. PhosphorusDemand) then
                PhosphorusDemand = PhosphorusDemand
            endif
                        
            PotentialAnnualPhosphorus = Me%AnnualPhosphorusFertilized(i,j) + PhosphorusDemand
            if (PotentialAnnualPhosphorus .ge. PhosphorusAnnualMax) then
                PhosphorusDemand = PhosphorusAnnualMax - (PotentialAnnualPhosphorus - PhosphorusDemand)
            endif

!!         if (targn <= 1.e-6) return
            
            FertilizerAmount = 0.0
            MineralPFracInFertilizer = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%MineralPFracInFertilizer
            
            if (MineralPFracInFertilizer .gt. 0.0001) then
                Me%AnnualPhosphorusFertilized(i,j) = Me%AnnualPhosphorusFertilized(i,j) + PhosphorusDemand
                FertilizerAmount = PhosphorusDemand / MineralPFracInFertilizer
            else
                FertilizerAmount = 0.0
            endif
            
            FertilizerFracApplyedInSurface = &
         Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%FertilizerFracApplyedInSurface
            OrganicFracParticulate = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%OrganicFracParticulate 
            MineralPFracInFertilizer       = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%MineralPFracInFertilizer
            OrganicPFracInFertilizer       = &
            Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(1)%Fertilizer%OrganicPFracInFertilizer
            
            !KgP/ha
            Me%Fluxes%FertilOrganicPParticInFluff(i,j) = FertilizerAmount * OrganicPFracInFertilizer * OrganicFracParticulate
            
            Me%Fluxes%FertilOrganicPInSurface(i,j)    = FertilizerFracApplyedInSurface * FertilizerAmount          &
                                                        * OrganicPFracInFertilizer * (1. - OrganicFracParticulate)
            Me%Fluxes%FertilOrganicPInSubSurface(i,j) = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount   &
                                                        * OrganicPFracInFertilizer * (1. - OrganicFracParticulate)

            Me%Fluxes%FertilMineralPInSurface(i,j)    = FertilizerFracApplyedInSurface * FertilizerAmount          &
                                                        * MineralPFracInFertilizer 
            Me%Fluxes%FertilMineralPInSubSurface(i,j) = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount   &
                                                        * MineralPFracInFertilizer 

        endif

    end subroutine AutoFertilization

    !--------------------------------------------------------------------------


    subroutine ScheduledFertilization(i,j, FertApp)
        !Argument--------------------------------------------------------------
        integer, intent(IN)                         :: i,j
        integer, intent(IN)                         :: FertApp
        !Local-----------------------------------------------------------------
        real                                        :: FertilizerAmount
        real                                        :: MineralNFracInFertilizer
        real                                        :: FertilizerFracApplyedInSurface, OrganicFracParticulate
        real                                        :: AmmoniaFracInMineralN, OrganicNFracInFertilizer
        real                                        :: MineralPFracInFertilizer, OrganicPFracInFertilizer
        !Begin-----------------------------------------------------------------

        !Compute fertilizer applyance
            !KgN/ha

        MineralNFracInFertilizer = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%MineralNFracInFertilizer
        FertilizerFracApplyedInSurface = &
      Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%FertilizerFracApplyedInSurface
        OrganicFracParticulate = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicFracParticulate
        AmmoniaFracInMineralN          = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%AmmoniaFracInMineralN
        OrganicNFracInFertilizer       = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicNFracInFertilizer
        MineralPFracInFertilizer       = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%MineralPFracInFertilizer
        OrganicPFracInFertilizer       = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%Fertilizer%OrganicPFracInFertilizer
        
        FertilizerAmount               = &
        Me%VegetationTypes(Me%VegetationID(i,j))%FertilizerDatabase%FertilizerApps(FertApp)%FertilizerAppAmount
        
        if (Me%ComputeOptions%ModelNitrogen) then
        
            !KgN/ha
            Me%Fluxes%FertilNitrateInSurface(i,j)     = FertilizerFracApplyedInSurface * FertilizerAmount                      &
                                                        * MineralNFracInFertilizer * (1. - AmmoniaFracInMineralN)
            Me%Fluxes%FertilNitrateInSubSurface(i,j)  = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount               &
                                                        * MineralNFracInFertilizer * (1. - AmmoniaFracInMineralN)
            Me%Fluxes%FertilAmmoniaInSurface(i,j)     = FertilizerFracApplyedInSurface * FertilizerAmount                      &
                                                        * MineralNFracInFertilizer * AmmoniaFracInMineralN
            Me%Fluxes%FertilAmmoniaInSubSurface(i,j)  = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount               &
                                                        * MineralNFracInFertilizer * AmmoniaFracInMineralN
            Me%Fluxes%FertilOrganicNParticInFluff(i,j) = FertilizerAmount * OrganicNFracInFertilizer * OrganicFracParticulate
                                                                    
            Me%Fluxes%FertilOrganicNInSurface(i,j)    = FertilizerFracApplyedInSurface * FertilizerAmount                      &
                                                        * OrganicNFracInFertilizer * (1. - OrganicFracParticulate)
            Me%Fluxes%FertilOrganicNInSubSurface(i,j) = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount               &
                                                        * OrganicNFracInFertilizer * (1. - OrganicFracParticulate)
        endif
        
        if (Me%ComputeOptions%ModelPhosphorus) then
            !KgP/ha
            Me%Fluxes%FertilOrganicPParticInFluff(i,j) = FertilizerAmount * OrganicPFracInFertilizer * OrganicFracParticulate
            
            Me%Fluxes%FertilOrganicPInSurface(i,j)    = FertilizerFracApplyedInSurface * FertilizerAmount          &
                                                        * OrganicPFracInFertilizer * (1. - OrganicFracParticulate)
            Me%Fluxes%FertilOrganicPInSubSurface(i,j) = (1. - FertilizerFracApplyedInSurface) * FertilizerAmount   &
                                                        * OrganicPFracInFertilizer * (1. - OrganicFracParticulate)
           
            !KgP/ha
            Me%Fluxes%FertilMineralPInSurface(i,j)    =  FertilizerFracApplyedInSurface * FertilizerAmount          &
                                                         * MineralPFracInFertilizer 
            Me%Fluxes%FertilMineralPInSubSurface(i,j) =  (1. - FertilizerFracApplyedInSurface) * FertilizerAmount   &
                                                         * MineralPFracInFertilizer 
        endif
        
    end subroutine ScheduledFertilization

    !--------------------------------------------------------------------------

    
    subroutine PesticideFluxes
    
        !Local-----------------------------------------------------------------        
        integer                                            :: i,j
        real                                               :: PesticideAmount
        real                                               :: ExtinctionCoef
        real                                               :: GroundCover
        integer                                            :: NumberPestApp, UniquePest
        integer                                            :: PestApp
        !Begin-----------------------------------------------------------------

do1:    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do2:    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
if1:        if (Me%ExternalVar%MappingPoints(i,j) == VegetationPoint) then    
            
                !Plant type 0 is NO PLANT, so no operations should happen
                if (Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%PlantType == NotAPlant) cycle
            
                NumberPestApp = Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%NumberPesticideApps
if2:            if (NumberPestApp .gt. 0) then
                    
                    !Check unique pesticide list if applications occured
do3:                do PestApp = 1, NumberPestApp
!do3:                do UniquePest = 1, Me%Fluxes%Pesticides%UniquePesticides
if3:                    if (&
                     Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideAppOccurred) then
!if3:                    if (Me%Fluxes%Pesticides%Application(UniquePest)%PesticideAppOccurred(i,j)) then
                            
                            !if application occurred check for that id in agricultural practices and put the 
                            !right amount
do4:                        do UniquePest = 1, Me%Fluxes%Pesticides%UniquePesticides                            
!do4:                        do PestApp = 1, NumberPestApp
if4:                            if (Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideID &
                                     == Me%Fluxes%Pesticides%Application(UniquePest)%PesticideID) then
                                        
                                        !kg/ha.day
                                        PesticideAmount =  &
                        Me%VegetationTypes(Me%VegetationID(i,j))%PesticideDatabase%PesticideApps(PestApp)%PesticideAppAmount
                                        ExtinctionCoef = Me%VegetationTypes(Me%VegetationID(i,j))%GrowthDatabase%ExtinctCoef
                                        GroundCover = 1.0 - exp(-ExtinctionCoef * Me%StateVariables%LeafAreaIndex(i, j)) 
                                        !kg/ha
                                        Me%Fluxes%Pesticides%Application(UniquePest)%Soil(i,j)   = PesticideAmount       &
                                         * (1 - GroundCover)
                                        !kg/ha
                                        Me%Fluxes%Pesticides%Application(UniquePest)%Vegetation(i,j) = PesticideAmount   &
                                         * (GroundCover)  
                                        
                                        exit do4
                                endif if4                                   
                            enddo do4
                           
                        endif if3
                    enddo do3
                
                endif if2

            endif if1
            
        enddo do2
        enddo do1

    endsubroutine PesticideFluxes

    !--------------------------------------------------------------------------
    
    subroutine Modify_OutPutHDF

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        type(T_Time)                                :: Actual
        integer                                     :: OutPutNumber
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePointer
        integer                                     :: ILB, IUB, JLB, JUB                 
        !Begin----------------------------------------------------------------

        Actual   = Me%ExternalVar%Now        
        OutPutNumber = Me%OutPut%NextOutput
        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        if (Actual .GE. Me%OutPut%OutTime(OutPutNumber)) then 
            
            call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            
            TimePointer => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR00'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time",                   &
                                 "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D      = TimePointer,                    &
                                 OutputNumber = OutputNumber,                   &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR01'

            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR02'

!            call GetBasinPoints  (Me%ObjHorizontalMap, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR02'
!
!            !Writes the Open Points
!            call HDF5WriteData   (Me%ObjHDF5, "/Grid/BasinPoints",              &
!                                  "BasinPoints", "-",                            &
!!                                  Array2D = Me%ExternalVar%OpenPoints2D,        &
!                                  Array2D = Me%ExternalVar%BasinPoints,         &
!                                  OutputNumber = OutPutNumber,                  &
!                                  STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR03'
!
!            call UngetBasin (Me%ObjHorizontalMap, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR200'


            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))

                if (PropertyX%OutputHDF) then

                    call HDF5WriteData   (Me%ObjHDF5,                                    &
                                          "/Results/"//PropertyX%ID%Name,                &
                                          PropertyX%ID%Name,                             &
                                          PropertyX%ID%Units,                            &
                                          Array2D = PropertyX%Field,                     &
                                          OutputNumber = OutPutNumber,                   &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR04'

                endif

                PropertyX => PropertyX%Next

            enddo

            if (Me%ComputeOptions%Evolution%ModelSWAT) then
                call HDF5WriteData  (Me%ObjHDF5, "//Results/HUAccumulated",          &
                                     "HUAccumulated", "-",                              &
                                     Array2D      = Me%HeatUnits%PlantHUAccumulated,    &
                                     OutputNumber = OutPutNumber,                       &
                                     STAT         = STAT_CALL)                      
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OutPutHDF - ModuleVegetation - ERR05'                 
            endif

            call HDF5WriteData  (Me%ObjHDF5, "/Results/WaterUptake",&
                                 "WaterUptake", "m3/s",                                  &
                                 Array2D      = Me%Fluxes%WaterUptake,                   &
                                 OutputNumber = OutPutNumber,                            &
                                 STAT         = STAT_CALL)                      
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'OutPutHDF - ModuleVegetation - ERR06'                 

            call HDF5WriteData  (Me%ObjHDF5, "/Results/WaterStress",&
                                 "WaterStress", "-",                                     &
                                 Array2D      = Me%Growth%WaterStress,                   &
                                 OutputNumber = OutPutNumber,                            &
                                 STAT         = STAT_CALL)                      
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'OutPutHDF - ModuleVegetation - ERR07'                 
            
            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleVegetation - ERR08'

            Me%OutPut%NextOutput = OutPutNumber + 1

        
        endif

    end subroutine Modify_OutPutHDF

    !--------------------------------------------------------------------------

    subroutine Modify_OutPutTimeSeries

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        integer                                     :: Pest
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

        if (Me%ComputeOptions%ModelNitrogen) then

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Fluxes%NitrogenUptake, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR04'
            
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%OptimalTotalPlantNitrogen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR05'
            
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%NitrogenStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR06'
            
        endif
        if (Me%ComputeOptions%ModelPhosphorus) then

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Fluxes%PhosphorusUptake, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR07'
            
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%OptimalTotalPlantPhosphorus, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR08'
           
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%PhosphorusStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR09'

        endif


        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
            
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%HeatUnits%PlantHUAccumulated, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR010'

            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%HeatUnits%PotentialHUBase, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR011'
        
            call WriteTimeSerie(Me%ObjTimeSerie, Data2D = Me%Growth%TemperatureStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR012'

            !Atmosphere output file
            !Averaged Atmosphere properties 
            if (Me%ComputeOptions%AtmospherePropertiesOutput) then

                call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%ExternalVar%Integration%AverageAirTemPDuringDT,     &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR020'     

                call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%ExternalVar%Integration%AverageAirHumidityDuringDT, &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR021'     

                call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%ExternalVar%Integration%AverageRadiationDuringDT,   &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR022'     

                if (Me%ComputeOptions%ModelPlantBiomass) then
                    call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%Growth%PAR, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR023'     

                    call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%Growth%RUE, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR024'     

                    call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%Growth%PotentialGrowth, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR025'   
                    
                    call WriteTimeSerie (Me%ObjTimeSerieAtm, Data2D = Me%Growth%PotentialBiomass, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR026' 
                    
                endif  
        
            endif
        
            !Fluxes to Soil output file
            if (Me%ComputeOptions%FluxesToSoilOutput) then

                if (Me%ComputeOptions%Fertilization) then
                    if (Me%ComputeOptions%ModelNitrogen) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilNitrateToSoilSurface,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR031'     

                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilNitrateToSoilSubSurface, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR032'     
 
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilAmmoniaToSoilSurface,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR033'     

                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilAmmoniaToSoilSubSurface, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR034'     

                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilOrganicNToSoilSurface,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR035'     

                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilOrganicNToSoilSubSurface, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR036'   

                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilOrganicNParticToFluff, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR036.1'   
                        
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%AnnualNitrogenFertilized, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR036.5'   
                    endif
                    if (Me%ComputeOptions%ModelPhosphorus) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilOrganicPToSoilSurface,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR037'     
                    
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilOrganicPToSoilSubSurface, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR038'     
    
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilMineralPToSoilSurface,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR039'     
                    
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilMineralPToSoilSubSurface, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR040'     

                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%FertilOrganicPParticToFluff,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR040.1'   
    
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%AnnualPhosphorusFertilized, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR040.5'                  
                    endif
                endif
                if (Me%ComputeOptions%HarvestKill) then
                    call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%HarvestKillBiomassToSoil,        &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR041'     

                    call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil,       &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR042'     

                    if (Me%ComputeOptions%ModelNitrogen) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil,   &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR043'     
                    endif
                    if (Me%ComputeOptions%ModelPhosphorus) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR044'     
                    endif
                endif                
                if (Me%ComputeOptions%Grazing) then
                    call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%GrazingBiomassToSoil,           &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR045'     
                
                    if (Me%ComputeOptions%ModelNitrogen) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%GrazingNitrogenToSoil,      &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR046'     
                    endif
                    if (Me%ComputeOptions%ModelPhosphorus) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%GrazingPhosphorusToSoil,    &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR047'     
                    endif
                endif   
                if (Me%ComputeOptions%Dormancy) then
                    call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%DormancyBiomassToSoil,          &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR048'     
                
                    if (Me%ComputeOptions%ModelNitrogen) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%DormancyNitrogenToSoil,     &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR049'     
                    endif
                    if (Me%ComputeOptions%ModelPhosphorus) then
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%ToSoil%DormancyPhosphorusToSoil,   &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR050'     
                    endif
                endif  
                if (Me%ComputeOptions%Pesticide) then
                    do Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%Pesticides%Application(Pest)%Soil, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR060'
                        call WriteTimeSerie (Me%ObjTimeSerieToSoil, Data2D = Me%Fluxes%Pesticides%Application(Pest)%Vegetation, &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR060'                        
                    enddo                        
                endif

            endif

        endif


!        if (Me%ComputeOptions%Fertilization .and. Me%ComputeOptions%FertilizationOutput) then
!
!            if(Me%ComputeOptions%ModelNitrogen) then
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilNitrateInSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR031'     
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilNitrateInSubSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR032'     
! 
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilAmmoniaInSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR033'     
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilAmmoniaInSubSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR034'     
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilOrganicNInSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR035'     
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilOrganicNInSubSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR036'   
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%AnnualNitrogenFertilized, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR036.5'   
!            
!            endif  
!
!            if(Me%ComputeOptions%ModelPhosphorus) then
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilOrganicPInSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR037'     
!                
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilOrganicPInSubSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR038'     
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilMineralPInSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR039'     
!                
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%Fluxes%FertilMineralPInSubSurface, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR040'     
!
!                call WriteTimeSerie (Me%ObjTimeSerieFert, Data2D = Me%AnnualPhosphorusFertilized, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'OutPutTimeSeries - ModuleVegetation - ERR040.5'   
!            
!            endif
!
!        
!        endif

        

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

                !Write Output for continuous computation
                call Write_FinalVegetation_HDF
 !               call Write_FinalVegetation_File

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

                nUsers = DeassociateInstance (mBasinGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillVegetation - ModuleVegetation - ERR05.1'

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

                if ((Me%ObjTimeSerieToSoil > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerieToSoil, STAT = STAT_)
                    if (STAT_ .NE. SUCCESS_) stop 'KillVegetation - ModuleVegetation - ERR10'
                endif

                if ((Me%ObjTimeSerieAtm > 0)) then
                    call KillTimeSerie(Me%ObjTimeSerieAtm, STAT = STAT_)
                    if (STAT_ .NE. SUCCESS_) stop 'KillVegetation - ModuleVegetation - ERR11'
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

    !   Write the final vegetation results in HDF format  !
  
    subroutine Write_FinalVegetation_HDF

        !Local--------------------------------------------------------------
        type(T_Property),           pointer     :: Property
        integer                                 :: ObjHDF5
        character (Len = StringLength)          :: PropertyName
        integer                                 :: WorkILB, WorkIUB
        integer                                 :: WorkJLB, WorkJUB
        integer                                 :: WorkKLB, WorkKUB
        integer                                 :: STAT_CALL
        integer(4)                              :: HDF5_CREATE
        character (Len = Pathlength)            :: filename
        real, dimension(6), target              :: AuxTime
        real, dimension(:), pointer             :: TimePtr
        type (T_Time)                           :: Actual
        integer, dimension(:,:), pointer        :: PlantGrowingInteger, PlantDormantInteger
        character (Len = StringLength)          :: ConvertType 
        
        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB

        !Gets a pointer to Topography
        call GetGridData        (Me%ObjGridData, Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR01'

        !Gets a pointer to OpenPoints2D
        call GetBasinPoints  (Me%ObjHorizontalMap, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR02'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)


        filename = trim(Me%Files%FinalFile)


        ObjHDF5 = 0
        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,                                                     &
                            trim(filename)//"5",    &
                            HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR10'


        Actual   = Me%ExternalVar%Now
         
        call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),          &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
        !Writes Time
        TimePtr => AuxTime
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR11'

        call HDF5WriteData  (ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                             Array1D = TimePtr, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR12'

       

        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR20'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR25'


        !Writes the Grid
        call HDF5WriteData   (ObjHDF5, "//Grid/Topography", "Topography", "m",           &
                              Array2D = Me%ExternalVar%Topography,                       &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR30'


        call HDF5WriteData  (ObjHDF5, "//Grid/BasinPoints", "BasinPoints",               &
                             "-", Array2D = Me%ExternalVar%BasinPoints,                &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR110'



        Property => Me%FirstProperty

        do while (associated(Property))
      
            PropertyName = trim(adjustl(Property%ID%name))
            
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR120'

            !Final value
            call HDF5WriteData  (ObjHDF5, "/Results/"//Property%ID%Name,                &
                                 Property%ID%Name,                                      &
                                 Property%ID%Units,                                     &
                                 Array2D = Property%Field,                              &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR130'

            
            Property => Property%Next

        enddo

        nullify (Property)

        if (Me%ComputeOptions%Evolution%ModelSWAT) then
            
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                             &
                                 WorkJLB, WorkJUB,                                      &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR140'

            
            allocate (PlantGrowingInteger(Me%Worksize%ILB:Me%Worksize%IUB, Me%Worksize%JLB:Me%Worksize%JUB))
            
            PlantGrowingInteger (:,:) = 0
            
            ConvertType = "LogicalToInteger"
            call ConvertLogicalInteger (Me%IsPlantGrowing, PlantGrowingInteger, ConvertType)
            
            !Final value
            call HDF5WriteData  (ObjHDF5, "/Results/"//"IsPlantGrowing",                &
                                 "IsPlantGrowing",                                       &
                                 "-",                                                  &
                                 Array2D = PlantGrowingInteger,                        &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR141'
            
            deallocate (PlantGrowingInteger)
            
            
            if (Me%ComputeOptions%Dormancy) then            

                allocate (PlantDormantInteger(Me%Worksize%ILB:Me%Worksize%IUB, Me%Worksize%JLB:Me%Worksize%JUB))
                
                PlantDormantInteger (:,:) = 0
                
                ConvertType = "LogicalToInteger"
                call ConvertLogicalInteger (Me%IsPlantDormant, PlantDormantInteger, ConvertType)
                
                !Final value
                call HDF5WriteData  (ObjHDF5, "/Results/"//"IsPlantDormant",                &
                                     "IsPlantDormant",                                       &
                                     "-",                                                  &
                                     Array2D = PlantDormantInteger,                        &
                                     STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR142'
                
                deallocate (PlantDormantInteger)
            
            endif


            !Final value
            call HDF5WriteData  (ObjHDF5, "/Results/"//"HUAccumulated",                 &
                                 "HUAccumulated",                                       &
                                 "HU",                                                  &
                                 Array2D = Me%HeatUnits%PlantHUAccumulated,             &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR150'



            call HDF5WriteData  (ObjHDF5, "/Results/"//"PotentialHUBase",               &
                                 "PotentialHUBase",                                     &
                                 "HU",                                                  &
                                 Array2D = Me%HeatUnits%PotentialHUBase,                &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR170'


            call HDF5WriteData  (ObjHDF5, "/Results/"//"TreeCurrentYear",               &
                                 "TreeCurrentYear",                                     &
                                 "years",                                               &
                                 Array2D = Me%Growth%TreeCurrentYear,                   &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR180'

            call HDF5WriteData  (ObjHDF5, "/Results/"//"PlantLAIMaxFraction",           &
                                 "PlantLAIMaxFraction",                                 &
                                 "-",                                                   &
                                 Array2D = Me%PlantLAIMaxFraction,                      &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR190'
        endif

        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then
            
            if (Me%ComputeOptions%Dormancy .or. Me%ComputeOptions%Grazing) then

                if (Me%ComputeOptions%ModelNitrogen) then
                    !Needed for dormancy and grazing
                    call HDF5WriteData  (ObjHDF5, "/Results/"//"NitrogenFraction",              &
                                         "NitrogenFraction",                                    &
                                         "-",                                                   &
                                         Array2D = Me%PlantNitrogenFraction,                    &
                                         STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                                  &
                        stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR191'
                endif
                if (Me%ComputeOptions%ModelPhosphorus) then
                    call HDF5WriteData  (ObjHDF5, "/Results/"//"PhosphorusFraction",            &
                                         "PhosphorusFraction",                                  &
                                         "-",                                                   &
                                         Array2D = Me%PlantPhosphorusFraction,                  &
                                         STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                                  &
                        stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR192'
                endif            
            endif

            if (Me%ComputeOptions%Fertilization) then

                if (Me%ComputeOptions%ModelNitrogen) then
                    !Needed for dormancy and grazing
                    call HDF5WriteData  (ObjHDF5, "/Results/"//"AnnualNitrogenFertilized",      &
                                         "AnnualNitrogenFertilized",                            &
                                         "kg/ha",                                               &
                                         Array2D = Me%AnnualNitrogenFertilized,                 &
                                         STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                                  &
                        stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR191'
                endif
                if (Me%ComputeOptions%ModelPhosphorus) then
                    call HDF5WriteData  (ObjHDF5, "/Results/"//"AnnualPhosphorusFertilized",    &
                                         "AnnualPhosphorusFertilized",                          &
                                         "kg/ha",                                               &
                                         Array2D = Me%AnnualPhosphorusFertilized,               &
                                         STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                                  &
                        stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR192'
                endif
            
            endif

        endif
  
        call UngetBasin (Me%ObjHorizontalMap, Me%ExternalVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR200'
  
        !UnGets Topography
        call UnGetGridData      (Me%ObjGridData, Me%ExternalVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_FinalVegetation_HDF - ModuleVegetation - ERR205'

   
        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_FinalWaterProperties_HDF - ModuleWaterProperties - ERR210'

        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_FinalWaterProperties_HDF - ModuleWaterProperties - ERR220'

    end subroutine Write_FinalVegetation_HDF

    !--------------------------------------------------------------------------

!    subroutine Write_FinalVegetation_File
!
!        !Arguments-------------------------------------------------------------
!
!        !Local-----------------------------------------------------------------
!        real                                        :: Year_File, Month_File, Day_File
!        real                                        :: Hour_File, Minute_File, Second_File
!        integer                                     :: FinalFile
!        integer                                     :: STAT_CALL
!        type (T_Property), pointer                  :: PropertyX
!
!        !----------------------------------------------------------------------
!
!        call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalVegFile - ModuleVegetation - ERR01'
!
!        open(Unit = FinalFile, File = Me%Files%FinalFile, Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalVegFile - ModuleVegetation - ERR02'
!
!        !Writes Date
!        call ExtractDate(Me%ExternalVar%Now, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)
!        write(FinalFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
!
!        PropertyX => Me%FirstProperty
!        do while (associated(PropertyX))
!
!            write(FinalFile)PropertyX%Field
!                
!            PropertyX => PropertyX%Next
!
!        enddo
!        
!        if (Me%ComputeOptions%Evolution%ModelSWAT) then
!            write(FinalFile)Me%IsPlantGrowing
!            write(FinalFile)Me%IsPlantDormant
!        endif
!
!        call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalVegFile - ModuleVegetation - ERR03'
!        
!
!
!    end subroutine Write_FinalVegetation_File

    !------------------------------------------------------------------------  
   
    subroutine DeAllocateVariables

        !Local------------------------------------------------------------------
        type (T_Property), pointer :: PropertyX
        integer                    :: STAT_CALL
        integer                    :: Pest
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

        if (Me%AgricPract%ObjFillMatrix /= 0) then
            call KillFillMatrix(Me%AgricPract%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'DeAllocateVariables - ModuleVegetation - ERR03'
        endif
        
        !Sets the number of properties equal to the FillValueInt
        Me%PropertiesNumber = FillValueInt

        Nullify   (Me%FirstProperty,Me%LastProperty)
        
        deallocate(Me%VegetationID                                         )
        
        deallocate (Me%AgricPractID                                        )
        
        deallocate(Me%RootDepthOld                                         )
        
        !Soil Fluxes
        deallocate(Me%Fluxes%WaterUptake                                   )  
        deallocate(Me%Fluxes%WaterUptakeLayer                              )  
!        deallocate(Me%Fluxes%FromSoil%WaterUptakeFromSoil                 ) 
        deallocate(Me%TranspirationBottomLayer                             )
         
        if(Me%ComputeOptions%ModelNitrogen) then
            deallocate(Me%Fluxes%NitrogenUptake                            )  
            deallocate(Me%Fluxes%NitrogenUptakeLayer                       )  
        endif
        
        if (Me%ComputeOptions%ModelPhosphorus) then
            deallocate(Me%Fluxes%PhosphorusUptake                          ) 
            deallocate(Me%Fluxes%PhosphorusUptakeLayer                     )  
        endif        
        if (Me%ComputeOptions%ModelNitrogen .or. Me%ComputeOptions%ModelPhosphorus) then
            deallocate(Me%SoilFluxesActive                                         )  
        endif
        deallocate(Me%Growth%WaterStress                                    ) 
        deallocate(Me%ExternalVar%Integration%SumPotTP                      )
        deallocate(Me%ExternalVar%Integration%AveragePotTPDuringDT  )

        if (Me%ComputeOptions%Evolution%ModelSWAT) then
            deallocate(Me%IsPlantGrowing                                    )
            deallocate(Me%PlantingOccurred                                  )
            deallocate(Me%HeatUnits%PotentialHUTotal                        ) 
            deallocate(Me%HeatUnits%PotentialHUBase                         ) 
            deallocate(Me%HeatUnits%PotentialHUBase_Old                     ) 
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
            deallocate(Me%Growth%TreeComingFromContinuous                   )        
            deallocate(Me%Growth%TreeFractionToMaturity                     )        
            deallocate(Me%Growth%TreeMaximumAnnualBiomass                   )        


            !Soil fluxes
            if(Me%ComputeOptions%ModelNitrogen) then
                deallocate(Me%PlantNitrogenFraction                         ) 
                deallocate(Me%OptimalTotalPlantNitrogen                     ) 
!                deallocate(Me%Fluxes%FromSoil%NitrogenUptakeFromSoil        ) 
            endif
            
            if (Me%ComputeOptions%ModelPhosphorus) then
                deallocate(Me%PlantPhosphorusFraction                       )
                deallocate(Me%OptimalTotalPlantPhosphorus                   ) 
!                deallocate(Me%Fluxes%FromSoil%PhosphorusUptakeFromSoil      ) 
            endif
            
            !Aerial Fluxes
            deallocate(Me%Fluxes%LAIChange                                  ) 
             
            deallocate(Me%LAISenescence                                     )             
            deallocate(Me%PlantLAIMaxFraction                               ) 
            if(.not. Me%ComputeOptions%ChangeLAISenescence) then
                deallocate(Me%LAIBeforeSenescence                           ) 
            endif



        endif
        
        if (Me%ComputeOptions%Evolution%GrowthModelNeeded) then

            if (Me%ComputeOptions%Grazing) then
                deallocate(Me%DaysOfGrazing                                     ) 
                deallocate(Me%IsPlantBeingGrazed                                ) 

!                deallocate(Me%GrazingFinished                                   )
!                deallocate(Me%GrazingOperations                                 )
            
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
        
            if (Me%ComputeOptions%HarvestKill) then
                deallocate(Me%HarvestOnlyOccurred                               )
                deallocate(Me%HarvestKillOccurred                               )
                deallocate(Me%KillOccurred                                      )
!                deallocate(Me%HarvestFinished                                   )
!                deallocate(Me%HarvestOperations                                 )
!                deallocate(Me%HarvestKillOperations                             ) 
!                deallocate(Me%KillOperations                                    ) 
                      
                deallocate(Me%Fluxes%BiomassRemovedInHarvest                    ) 
                deallocate(Me%Fluxes%BiomassHarvestedFraction                   ) 
                deallocate(Me%Fluxes%ToSoil%HarvestKillBiomassToSoil            ) 
                deallocate(Me%Fluxes%ToSoil%KillRootBiomassLeftInSoil           )

                if (Me%ComputeOptions%ModelNitrogen) then            
                    deallocate(Me%Fluxes%NitrogenRemovedInHarvest               )  
                    deallocate(Me%Fluxes%ToSoil%HarvestKillNitrogenToSoil       ) 
                endif
            
                if (Me%ComputeOptions%ModelPhosphorus) then
                    deallocate(Me%Fluxes%PhosphorusRemovedInHarvest             ) 
                    deallocate(Me%Fluxes%ToSoil%HarvestKillPhosphorusToSoil     )
                endif

            endif

            if (Me%ComputeOptions%Dormancy) then
                deallocate(Me%DayLength                                         )
                deallocate(Me%MinimumDayLength                                  )
                deallocate(Me%IsPlantDormant                                    ) 
                deallocate(Me%PlantGoingDormant                                 )
            
                deallocate(Me%Fluxes%BiomassRemovedInDormancy                   ) 
    !            deallocate(Me%Fluxes%BiomassDormancyFraction                    ) 

    !            deallocate(Me%Fluxes%ToSoil%DormancyBiomassToSoil               ) 
            
                if (Me%ComputeOptions%ModelNitrogen) then
                    deallocate(Me%Fluxes%NitrogenRemovedInDormancy              )  
    !                deallocate(Me%Fluxes%ToSoil%DormancyNitrogenToSoil          ) 
                endif
            
                if (Me%ComputeOptions%ModelPhosphorus) then
                    deallocate(Me%Fluxes%PhosphorusRemovedInDormancy            ) 
    !                deallocate(Me%Fluxes%ToSoil%DormancyPhosphorusToSoil        ) 
                endif
 
            endif
        
            if (Me%ComputeOptions%Fertilization) then
!                deallocate(Me%FertilizationOccurred                               )
                if (Me%ComputeOptions%ModelNitrogen) then
                    deallocate(Me%Fluxes%FertilNitrateInSurface          )
                    deallocate(Me%Fluxes%FertilNitrateInSubSurface       )
                    deallocate(Me%Fluxes%FertilAmmoniaInSurface          )
                    deallocate(Me%Fluxes%FertilAmmoniaInSubSurface       )
                    deallocate(Me%Fluxes%FertilOrganicNInSurface         )
                    deallocate(Me%Fluxes%FertilOrganicNParticInFluff     )
                    deallocate(Me%Fluxes%FertilOrganicNInSubSurface      )
                endif                                                           
                if (Me%ComputeOptions%ModelPhosphorus) then                     
                    deallocate(Me%Fluxes%FertilOrganicPInSurface         )
                    deallocate(Me%Fluxes%FertilOrganicPParticInFluff     )
                    deallocate(Me%Fluxes%FertilOrganicPInSubSurface      )
                    deallocate(Me%Fluxes%FertilMineralPInSurface         )
                    deallocate(Me%Fluxes%FertilMineralPInSubSurface      )
                endif

            endif        
            
            if (Me%ComputeOptions%Pesticide) then
                do  Pest = 1, Me%Fluxes%Pesticides%UniquePesticides
                    deallocate (Me%Fluxes%Pesticides%Application(Pest)%Soil       )
                    deallocate (Me%Fluxes%Pesticides%Application(Pest)%Vegetation )                        
!                    deallocate (Me%Fluxes%Pesticides%Application(Pest)%PesticideAppOccurred) 
                enddo                            
            endif

            if (Me%ComputeOptions%ModelPlantBiomass) then
                deallocate(Me%Fluxes%BiomassGrowth                              ) 
                deallocate(Me%Growth%BiomassGrowthOld                           )  
                deallocate(Me%Growth%PAR                                        )  
                deallocate(Me%Growth%RUE                                        )  
                deallocate(Me%Growth%PotentialGrowth                            )  
                deallocate(Me%Growth%PotentialBiomass                           )
            endif

            if (Me%ComputeOptions%ModelCanopyHeight) then
                deallocate(Me%ChangeCanopyEnabled                               ) 
            endif        
        
        endif

        if (Me%ComputeOptions%TranspirationMethod == TranspirationMOHID) then
            deallocate(Me%ComputeOptions%TranspirationMOHID%FeddesType   )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH1 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH2 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH4 )
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3L)
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesH3H)
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesR2L)
            deallocate(Me%ComputeOptions%TranspirationMOHID%RootFeddesR2H)
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

!            !Gets a pointer to OpenPoints2D
!            call GetOpenPoints2D  (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleVegetation - ERR110'
            

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

!            call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleVegetation - ERR70'


        endif

    end subroutine ReadUnlockExternalVar

end module ModuleVegetation

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 







