!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : InterfaceWaterAir
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Pedro Chambel Leitao, Luis Fernandes- v4.0
! DESCRIPTION   : Module to compute water-air interface fluxes
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

Module ModuleInterfaceWaterAir

    use ModuleGlobalData
    use ModuleTime
    use ModuleHDF5             
    use ModuleFunctions,            only: SaturatedVaporPressure, ConstructPropertyID,          &
                                          LatentHeat, SensibleHeat, LongWaveDownward,           &
                                          LongWaveUpward, AerationFlux, AerationFlux_CO2,       &
                                          SetMatrixValue, COAREInterfaceMoistureContent,        &
                                          LongWaveUpwardCOARE, LatentHeatOfVaporization,        &
                                          COAREMoistureContentAir, CHUNK_J
                                          
    use ModuleEnterData,            only: ConstructEnterData, GetData, ExtractBlockFromBuffer,  &
                                          Block_Unlock, GetOutPutTime, ReadFileName,            &
                                          KillEnterData   
    use ModuleDrawing                                                    
    use ModuleGridData,             only: GetGridData, UnGetGridData
    use ModuleHorizontalGrid,       only: GetHorizontalGridSize, WriteHorizontalGrid,           &
                                          RotateVectorFieldToGrid, RotateVectorGridToField,     &
                                          GetGridCellArea, UnGetHorizontalGrid, GetXYCellZ,     &
                                          GetDDecompMPI_ID, GetDDecompON,&
                                          GetGridOutBorderPolygon, UngetHorizontalGrid         
    use ModuleHorizontalMap,        only: GetWaterPoints2D, UnGetHorizontalMap          
    use ModuleGeometry,             only: GetGeometrySize, GetGeometryDistances, UnGetGeometry
    use ModuleMap,                  only: GetWaterPoints3D, UnGetMap                   
    use ModuleBoxDif,               only: StartBoxDif, BoxDif, KillBoxDif               
    use ModuleTimeSerie,            only: StartTimeSerie, WriteTimeSerie, KillTimeSerie,        &
                                          GetTimeSerieLocation, CorrectsCellsTimeSerie,         &
                                          GetNumberOfTimeSeries, TryIgnoreTimeSerie, GetTimeSerieName
    use ModuleWaterProperties,      only: GetDensity, GetConcentration, UnGetWaterProperties,   &
                                          GetWaterPropertiesAirOptions, SetSurfaceFlux,         &
                                          GetPropertySurfaceFlux
    use ModuleHydrodynamic,         only: GetHydrodynamicAirOptions,SetAtmosphericPressure,     &
                                          SetSurfaceWaterFlux, SetWindStress,                   &
                                          GetHorizontalVelocity, UngetHydrodynamic

#ifndef _LAGRANGIAN_                                         
#ifdef  _LAGRANGIAN_GLOBAL_                                         
    use ModuleLagrangianGlobal,     only: SetLagrangianShearGlobal, GetLagrangianAirOptionsGlobal,   &
                                          SetLagrangianWindGlobal, SetLagrangianAirTemperature,      &
                                          SetLagSolarRadiationGlobal,  &
                                          SetLagrangianAtmPressureGlobal   
#else
    use ModuleLagrangian,           only: SetLagrangianShear, GetLagrangianAirOptions,          &
                                          SetLagrangianWind, SetLagrangianSolarRadiation,       &
                                          SetLagrangianAtmPressure
#endif    
#endif

    use ModuleTurbGOTM,             only: SetTurbGOTMSurfaceRugosity, SetTurbGOTMWindShearVelocity, &
                                          SetTurbGOTMWaveSurfaceFluxTKE
    use ModuleAtmosphere,           only: GetAtmosphereProperty, AtmospherePropertyExists, UngetAtmosphere, &
                                          GetWindHeight, GetAirMeasurementHeight
    use ModuleFillMatrix,           only: ConstructFillMatrix, ModifyFillMatrix, ModifyFillMatrixVectorial, &
                                          GetDefaultValue, GetIfMatrixRemainsConstant, KillFillMatrix,      &
                                          GetValuesProcessingOptions

#ifndef _WAVES_
    use ModuleWaves,                only: SetWavesWind, GetWaves, UnGetWaves
#endif
    use ModuleStopWatch,           only : StartWatch, StopWatch         


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartInterfaceWaterAir
    private ::      AllocateInstance
    private ::      ReadWaterAirFilesName
    private ::      ConstructGlobalVariables
    private ::          AllocateCOAREVariables
    private ::      ConstructRugosity
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 16/12/2011
    private ::      ConstructWaveFluxTKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    private ::      Construct_PropertyList
    private ::          Construct_Property
    private ::              Construct_PropertyValues
    private ::              Construct_PropertyEvolution
    private ::              Construct_PropertyOutPut
    private ::          Add_Property
    private ::      CheckInternalOptions
    private ::      CheckOptionsWater
    private ::      CheckOptionsAir
    private ::      Construct_Sub_Modules
    private ::          Construct_Time_Serie
    private ::          StartOutputBoxFluxes
    
    private ::      SetSubModulesConstructor

    !Selector
                     
    !Modifier
    public  :: ModifyInterfaceWaterAir
    private ::      ModifyRugosity
    private ::          ComputeWavesRugosity
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    private ::      ModifyWaveFluxTKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    private ::      ModifyLocalAtmVariables
    private ::      ModifyWaterAirFluxes
    private ::          ModifyAlbedo
    private ::          ModifySensibleHeat
    private ::              ComputeSensibleHeat
    private ::          ModifyLatentHeat
    private ::              ComputeLatentHeat
    private ::          ModifyEvaporation
    private ::              ComputeEvaporation
    private ::          ModifyNetLongWaveRadiation
    private ::              ComputeUpLongWaveRad
    private ::              ComputeDownLongWaveRad
    private ::          ComputeCOAREHeatBudget
    private ::              ComputeWarmLayer
    private ::              ComputeCOAREsurfacefluxes
    private ::                  FirstGuess
    private ::                  BulkLoop
    private ::                  Psit
    private ::                  Psiu
    private ::                  AirViscosity
    private ::                  WaterVapourDiffusivity
    private ::                  HeatDiffusivity
    private ::          CheckRadiationOptions
    private ::          CheckLatentSensibleOptions
    private ::          ModifyOxygenFlux
    private ::          ModifyCarbonDioxideFlux
    private ::              ModifyAerationFlux
    private ::              ModifyCO2AerationFlux
    private ::          ModifyNitrateFlux
    private ::          ModifyAmmoniaFlux
    private ::          ModifySurfaceRadiation
    private ::              ComputeSurfaceRadiation
    private ::          ModifyNonSolarFlux
    private ::              ComputeNonSolarFlux
    private ::          ModifyWindShearVelocity
    private ::          ModifyTurbulentKE
    private ::              ComputeTKEWind
    private ::          ModifySurfaceWaterFlux
    private ::      Output_TimeSeries
    private ::      Output_BoxTimeSeries  
    private ::      OutPut_Results_HDF

    !Destructor
    public  :: KillInterfaceWaterAir                                                     
    private ::      DeAllocateInstance


    !Management
    private ::      Ready
    private ::          LocateObjInterfaceWaterAir
    
    private ::              ReadLockExternalGlobal
    private ::              ReadLockExternalWater
    private ::              ReadUnlockExternalGlobal
    private ::              ReadUnlockExternalWater
    
    !Interfaces----------------------------------------------------------------
    
    !Parameters----------------------------------------------------------------
    character(LEN = StringLength), parameter        :: prop_block_begin = '<beginproperty>'
    character(LEN = StringLength), parameter        :: prop_block_end   = '<endproperty>'

    character(LEN = StringLength), parameter        :: rugosity_begin   = '<begin_rugosity>'
    character(LEN = StringLength), parameter        :: rugosity_end     = '<end_rugosity>'

    !Coefficients
    real,    parameter                              :: WaterEmissivity          = 0.96
    real,    parameter                              :: StefanBoltzmann          = 5.669e-08     ![W/m2/K4]
    real,    parameter                              :: BowenCoefficient         = 0.47          ![mmHg/ºC]
    real,    parameter                              :: ReferenceDensity         = 1000.         ![kg/m3]
    real,    parameter                              :: RefLatentHeatOfVaporization = 2.5e6         ![J/kg]

    !Constant for Air-Sea kinetic energy transfer        
    real,    parameter                              :: CDE = 0.63E-06
    
    ! Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Wave height used for surface rugosity parametrization
    integer,    parameter                              :: NoWave                  = 0
    integer,    parameter                              :: SurfaceRugosityFromHS   = 1
    integer,    parameter                              :: SurfaceRugosityFromHSW  = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    
    !Methods to define the wind drag coefficient
    integer, parameter                              :: Constant = 1, WindFunction = 2, ShearVelocity = 3
    
    !Types---------------------------------------------------------------------
    private :: T_Files
    type       T_Files 
         character(len=PathLength)                  :: InputData    = null_str !initialization: Jauch
         character(len=PathLength)                  :: Results      = null_str !initialization: Jauch
         character(len=PathLength)                  :: BoxesFile    = null_str !initialization: Jauch
    end type T_Files

    private :: T_OutPut
    type       T_OutPut
         type (T_Time), pointer, dimension(:)       :: OutTime      => null()
         integer                                    :: NextOutPut   = null_int !initialization: Jauch
         logical                                    :: Yes          =.false.
         integer                                    :: Number       = null_int !initialization: Jauch
    end type T_OutPut
    
    private :: T_Ext_Global
    type       T_Ext_Global
        type(T_Time)                                :: Now
        real,    pointer, dimension(:,:)            :: GridCellArea => null()
        logical                                     :: Backtracking         = .false.
        real,    pointer, dimension(:,:,:)          :: DWZ
    end type T_Ext_Global


    private :: T_Ext_Options
    type       T_Ext_Options
        logical                                     :: HeatFluxYes                  = .false.
        logical                                     :: OxygenFluxYes                = .false.
        logical                                     :: CarbonDioxideFluxYes         = .false.
        logical                                     :: AmmoniaFluxYes               = .false.
        logical                                     :: NitrateFluxYes               = .false.       
        logical                                     :: WQMYes                       = .false.
        logical                                     :: SurfaceWaterFluxYes          = .false.
        logical                                     :: HydrodynamicWindYes          = .false.
        logical                                     :: HydrodynamicAtmPressureYes   = .false.
        logical                                     :: OilYes                       = .false.
        logical                                     :: HNSYes                       = .false.
        logical                                     :: WavesWindYes                 = .false.
        logical                                     :: GOTMWindShearVelocityYes     = .false.
        logical                                     :: Irrigation                   = .false.
        logical                                     :: Precipitation                = .false.
        logical                                     :: LagrangianWindYes            = .false.
        logical                                     :: LagrangianWQMYes             = .false.
        logical                                     :: LagrangianT90Yes             = .false.
        logical                                     :: T90VariableYes               = .false.
        logical                                     :: HydrodynamicMslpYes          = .false.
    end type    T_Ext_Options

    private :: T_Int_Options
    type       T_Int_Options
        logical                                     :: Evaporation                  = .false.
        logical                                     :: WindStress                   = .false.
        logical                                     :: SurfaceWaterFlux             = .false.
        logical                                     :: LatentHeat                   = .false.
        logical                                     :: SensibleHeat                 = .false.
        logical                                     :: NetLongWaveRadiation         = .false.
        logical                                     :: UpwardLongWaveRadiation      = .false.
        logical                                     :: DownwardLongWaveRadiation    = .false.
        logical                                     :: NonSolarFlux                 = .false.
        logical                                     :: OxygenFlux                   = .false.
        logical                                     :: CarbonDioxideFlux            = .false.
        logical                                     :: AmmoniaFlux                  = .false.
        logical                                     :: NitrateFlux                  = .false.
        logical                                     :: SpecificOxygenFlux           = .false.
        logical                                     :: SpecificCarbonDioxideFlux    = .false.        
        logical                                     :: WindShearVelocity            = .false.
        logical                                     :: WindShearVelAllocate         = .false. 
        logical                                     :: TurbulentKineticEnergy       = .false.
        logical                                     :: SurfaceRadiation             = .false.
        logical                                     :: Albedo                       = .false.
    end type   T_Int_Options

    
    private :: T_Ext_Water
    type       T_Ext_Water
        real,    pointer, dimension(:,:,:)          :: WaterTemperature => null()
        real,    pointer, dimension(:,:,:)          :: WaterSalinity    => null()    
        real,    pointer, dimension(:,:,:)          :: WaterVelocity    => null()
        real,    pointer, dimension(:,:,:)          :: WaterDepth       => null()
        real,    pointer, dimension(:,:,:)          :: Density          => null()
        integer, pointer, dimension(:,:  )          :: WaterPoints2D    => null()
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D    => null()
        real,    pointer, dimension(:,:  )          :: Bathymetry       => null()
    end type T_Ext_Water

    private :: T_ExtField
    type       T_ExtField
        real,    pointer, dimension(:,:  )          :: Field    => null()
        logical                                     :: Yes      = .false.      
    end type  T_ExtField


    private :: T_Ext_Atm
    type       T_Ext_Atm
        type(T_ExtField)                            :: AtmosphericPressure
        type(T_ExtField)                            :: Mslp     !Mean Sea Level Pressure
        type(T_ExtField)                            :: SolarRadiation
        type(T_ExtField)                            :: RelativeHumidity     
        type(T_ExtField)                            :: AirTemperature
        type(T_ExtField)                            :: CloudCover
        type(T_ExtField)                            :: WindVelocityU
        type(T_ExtField)                            :: WindVelocityV
        type(T_ExtField)                            :: WindDirection
        type(T_ExtField)                            :: AtmospDeposOxidNO3
        type(T_ExtField)                            :: AtmospDeposReduNH4
    end type T_Ext_Atm

    private :: T_Local_Atm
    type       T_Local_Atm
        type(T_ExtField)                            :: AtmosphericPressure
        type(T_ExtField)                            :: Mslp     !Mean Sea Level Pressure
        type(T_ExtField)                            :: WindVelocityU
        type(T_ExtField)                            :: WindVelocityV
        type(T_ExtField)                            :: WindDirection
    end type   T_Local_Atm

    
    private :: T_Evolution
        type   T_Evolution
        logical                                     :: Variable             = .false.
        real                                        :: DTInterval           = null_Real !initialization: Jauch
        type(T_Time)                                :: LastCompute
        type(T_Time)                                :: NextCompute
    end type T_Evolution
    
    private :: T_Property
    type       T_Property
         type(T_PropertyID)                         :: ID
         real, dimension(:,:), pointer              :: Field                => null() 
        !vectorial field rotated to grid cells - U comp.                
         real, dimension(:,:),  pointer             :: FieldU               => null() 
        !vectorial field rotated to grid cells - V comp.         
         real, dimension(:,:),  pointer             :: FieldV               => null() 
        !vectorial original field - X (zonal component)         
         real, dimension(:,:),  pointer             :: FieldX               => null() 
        !vectorial original field - Y (meridional comp.)                   
         real, dimension(:,:),  pointer             :: FieldY               => null() 
         type(T_Evolution)                          :: Evolution
         integer                                    :: SVPMethod            = 1
         integer                                    :: C1                   = 1
         integer                                    :: C2                   = 1
         logical                                    :: TimeSerie            = .false.
         logical                                    :: BoxTimeSerie         = .false.
         logical                                    :: CEQUALW2             = .false.
         logical                                    :: OutputHDF            = .false.
         logical                                    :: Constant             = .false.
         logical                                    :: FirstActualization   = .true.
         type(T_Property), pointer                  :: Next                 => null()
         type(T_Property), pointer                  :: Prev                 => null()
    end type T_Property

    private :: T_Coupling
    type       T_Coupling
         type(T_Time)                               :: NextCompute
         real                                       :: DT_Compute           = FillValueReal
         logical                                    :: Yes                  = .false.
         integer                                    :: NumberOfProperties   = 0
    end type T_Coupling  

    private :: T_Coupled
    type       T_Coupled
         type(T_Coupling)                           :: TimeSerie
         type(T_Coupling)                           :: BoxTimeSerie
         type(T_Coupling)                           :: CEQUALW2
    end type T_Coupled

    type       T_Rugosity
        type(T_PropertyID)                          :: ID
        real, pointer, dimension (:,:)              :: Field            => null()
        real                                        :: Scalar           = null_real, & !initialization: Jauch
                                                       WavesRelation    = null_real    !initialization: Jauch
!Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! logical                                     :: Constant, ON, WavesFunction

        logical                                     :: Constant, ON      = .false.
        integer                                     :: WavesFunction     = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end type T_Rugosity
    
!Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type       T_WaveFlux
        type(T_PropertyID)                          :: ID
        real, pointer, dimension (:,:)              :: Field            => null()
        logical                                     :: Constant, ON     = .false.
    end type T_WaveFlux
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    private :: T_InterfaceWaterAir
    type       T_InterfaceWaterAir
        integer                                     :: InstanceID   = null_int !initialization: Jauch
        character(PathLength)                       :: ModelName    = null_str !initialization: Jauch
        type(T_Time       )                         :: BeginTime
        type(T_Time       )                         :: EndTime
        type(T_Time       )                         :: ActualTime
        type(T_Time       )                         :: LastTimeInstant
        type(T_Size2D     )                         :: Size2D, WorkSize2D
        type(T_Size3D     )                         :: Size3D, WorkSize3D
        type(T_Files      )                         :: Files
        type(T_OutPut     )                         :: Output
        type(T_Property   ), pointer                :: FirstProperty
        type(T_Property   ), pointer                :: LastProperty
        type(T_Coupled    )                         :: Coupled
        type(T_Ext_Global )                         :: ExternalVar
        type(T_Ext_Water  )                         :: ExtWater
        type(T_Ext_Atm    )                         :: ExtAtm
        type(T_Ext_Options)                         :: ExtOptions
        type(T_Int_Options)                         :: IntOptions
        type(T_Local_Atm  )                         :: LocalAtm
        type(T_Rugosity   )                         :: Rugosity
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        type(T_WaveFlux   )                         :: WaveFluxTKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                         
        logical                                     :: COARE                    = .false.
        real                                        :: ReflectionCoef           = FillValueReal
        real                                        :: CDWIND                   = FillValueReal
        real                                        :: WindHeight               = FillValueReal
        real                                        :: AirMeasurementHeight     = FillValueReal
        real                                        :: WarmingAboveWaterTemp    = 0.
        real                                        :: BoundaryLayerDepth       = FillValueReal
        logical                                     :: DefineCDWIND             = .false.
        integer                                     :: CDWINDMethod             = FillValueInt
        real, pointer, dimension(:,:)               :: SurfaceTemperature       => null()
        real, pointer, dimension(:,:)               :: Fxp                      => null()
        real, pointer, dimension(:,:)               :: WarmLayerThickness       => null()
        real, pointer, dimension(:,:)               :: Tau                      => null()
        real, pointer, dimension(:,:)               :: Tau_ac                   => null()
        real, pointer, dimension(:,:)               :: AccumulatedEnergy        => null()
        real, pointer, dimension(:,:)               :: WarmLayerTempDiff        => null()
        real, pointer, dimension(:,:)               :: Al                       => null()
        real, pointer, dimension(:,:)               :: RainFlux                 => null()
        real, pointer, dimension(:,:)               :: LastLatentHeat           => null()
        real, pointer, dimension(:,:)               :: LastSensibleHeat         => null()
        real, pointer, dimension(:,:)               :: SurfaceAlbedo            => null()
        real(8), pointer, dimension(:,:)            :: Scalar2D                 => null()
        real   , pointer, dimension(:,:)            :: WindShearVelocity        => null()
        integer                                     :: AerationEquation         = FillValueInt
        integer                                     :: CO2AerationEquation      = FillValueInt
        real                                        :: Altitude                 = FillValueReal
        real                                        :: AltitudeCorrection       = FillValueReal
        real                                        :: ReNAtmDep                = FillValueReal
        real                                        :: OxNAtmDep                = FillValueReal 

        integer                                     :: PropertiesNumber         = 0
        logical                                     :: EnergyThreshold          = .false.
        logical                                     :: Jump                     = .true.
        logical                                     :: Jwarm                    = .false.
        real                                        :: Jcool                    = 0.
        
        logical                                     :: EvapReqConv              = .false.
        real                                        :: ConversionFactorEvap     = 1.0

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
                                                                                
        !Instance of ModuleHorizontalGrid                                       
        integer                                     :: ObjHorizontalGrid        = 0
                                                                                
        !Instance of ModuleGridData                                             
        integer                                     :: ObjGridData              = 0
                                                                                
        !Instance of ModuleGeometry                                             
        integer                                     :: ObjGeometry              = 0
        
        !Instance of ModuleHorizontalMap                                    
        integer                                     :: ObjHorizontalMap         = 0
                                                                                
        !Instance of ModuleMap
        integer                                     :: ObjMap                   = 0
        
        !Instance of ModuleHydrodynamic
        integer                                     :: ObjHydrodynamic          = 0

        !Instance of ModuleTurbGOTM
        integer                                     :: ObjTurbGOTM              = 0

        !Instance of ModuleWaves
        integer                                     :: ObjWaves                 = 0
        
        !Instance of ModuleInterfaceWaterAir
        integer                                     :: ObjWaterProperties       = 0

        !Instance of ModuleAtmosphere
        integer                                     :: ObjAtmosphere            = 0
        
        !Instance of ModuleLagrangian 
        integer                                     :: ObjLagrangian            = 0
        
        
        type(T_InterfaceWaterAir), pointer          :: Next => null()
    end type  T_InterfaceWaterAir

    !Global Module Variables
    type (T_InterfaceWaterAir), pointer             :: FirstObjInterfaceWaterAir    => null()
    type (T_InterfaceWaterAir), pointer             :: Me                           => null()

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartInterfaceWaterAir(ModelName,                      &
                                      ObjInterfaceWaterAirID,         &
                                      TimeID,                         &
                                      HorizontalGridID,               &
                                      WaterGridDataID,                &
                                      WaterHorizontalMapID,           &
                                      WaterMapID,                     &
                                      WaterGeometryID,                &
                                      HydrodynamicID,                 &
                                      TurbGOTMID,                     &
                                      WavesID,                        &
                                      WaterPropertiesID,              &
                                      AtmosphereID,                   &
                                      LagrangianID,                   &
                                      STAT)

        !Arguments---------------------------------------------------------------
        character(Len=*)                                :: ModelName
        integer                                         :: ObjInterfaceWaterAirID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: WaterGridDataID
        integer                                         :: WaterHorizontalMapID
        integer                                         :: WaterMapID
        integer                                         :: WaterGeometryID
        integer                                         :: HydrodynamicID        
        integer                                         :: TurbGOTMID
        integer                                         :: WavesID           
        integer                                         :: WaterPropertiesID
        integer                                         :: AtmosphereID
        integer                                         :: LagrangianID
        integer, optional, intent(OUT)                  :: STAT    

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL       

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterfaceWaterAir_)) then
            nullify (FirstObjInterfaceWaterAir)
            call RegisterModule (mInterfaceWaterAir_) 
        endif

        call Ready(ObjInterfaceWaterAirID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            nullify (Me%FirstProperty)
            nullify (Me%LastProperty )

            nullify (Me%LocalAtm%WindVelocityU%Field)
            nullify (Me%LocalAtm%WindVelocityV%Field)
            nullify (Me%LocalAtm%WindDirection%Field)
            nullify (Me%LocalAtm%AtmosphericPressure%Field)
            nullify (Me%LocalAtm%Mslp%Field)    !Mean Sea Level Pressure
            
            !Associates External Instances
            Me%ObjTime                 = AssociateInstance(mTIME_,           TimeID                     )
            Me%ObjHorizontalGrid       = AssociateInstance(mHORIZONTALGRID_, HorizontalGridID           )
            
            !Water column
            Me%ObjGridData             = AssociateInstance(mGRIDDATA_,          WaterGridDataID         )
            Me%ObjHorizontalMap        = AssociateInstance(mHORIZONTALMAP_,     WaterHorizontalMapID    )
            Me%ObjGeometry             = AssociateInstance(mGEOMETRY_,          WaterGeometryID         )
            Me%ObjMap                  = AssociateInstance(mMAP_,               WaterMapID              )
            Me%ObjHydrodynamic         = AssociateInstance(mHYDRODYNAMIC_,      HydrodynamicID          )
            Me%ObjWaterProperties      = AssociateInstance(mWATERPROPERTIES_,   WaterPropertiesID       )
            Me%ObjAtmosphere           = AssociateInstance(mATMOSPHERE_,        AtmosphereID            )

            if(LagrangianID /= 0)then
                Me%ObjLagrangian       = AssociateInstance(mLAGRANGIAN_,        LagrangianID            )
            endif

            if(TurbGOTMID /= 0)then
                Me%ObjTurbGOTM         = AssociateInstance(mTURBGOTM_,          TurbGOTMID              )
            endif

            if(WavesID /= 0)then
                Me%ObjWaves            = AssociateInstance(mWAVES_,             WavesID                 )
            endif



            call ReadLockExternalGlobal
            
            call ReadLockExternalWater

            call ReadWaterAirFilesName

            call ConstructEnterData(Me%ObjEnterData, Me%Files%InPutData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                 &
                stop 'StartInterfaceWaterAir - InterfaceWaterAir - ERR01'

            call ConstructGlobalVariables

            call ConstructRugosity
            
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call ConstructWaveFluxTKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call Construct_PropertyList

            call CheckInternalOptions

            call CheckOptionsWater
            
            call CheckOptionsAir

            call Construct_Sub_Modules

            call ConstructGlobalOutput
            
            if (Me%OutPut%Yes) call Open_HDF5_OutPut_File

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                 &
                stop 'StartInterfaceWaterAir - InterfaceWaterAir - ERR02'

            call ReadUnlockExternalWater

            call ReadUnlockExternalGlobal

            call SetSubModulesConstructor

            call SetSubModulesModifier

            !Returns ID
            ObjInterfaceWaterAirID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleInterfaceWaterAir - StartInterfaceWaterAir - ERR99' 

        endif cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartInterfaceWaterAir
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance
                                                    
        !Local-----------------------------------------------------------------
        type (T_InterfaceWaterAir), pointer    :: NewObjWaterSedInterface
        type (T_InterfaceWaterAir), pointer    :: PreviousObjWaterSedInterface


        !Allocates new instance
        allocate (NewObjWaterSedInterface)
        nullify  (NewObjWaterSedInterface%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjInterfaceWaterAir)) then
            FirstObjInterfaceWaterAir           => NewObjWaterSedInterface
            Me                                  => NewObjWaterSedInterface
        else
            PreviousObjWaterSedInterface        => FirstObjInterfaceWaterAir
            Me                                  => FirstObjInterfaceWaterAir%Next
            do while (associated(Me))
                PreviousObjWaterSedInterface    => Me
                Me                              => Me%Next
            enddo
            Me                                  => NewObjWaterSedInterface
            PreviousObjWaterSedInterface%Next   => NewObjWaterSedInterface
        endif

        Me%InstanceID = RegisterNewInstance (mINTERFACEWATERAIR_)


    end subroutine AllocateInstance
    
    
   
    !--------------------------------------------------------------------------

    
    subroutine ReadWaterAirFilesName

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len = StringLength)       :: Message

        !----------------------------------------------------------------------

        Message ='ASCII file used to construct new water-air properties.'
        Message = trim(Message)

        call ReadFileName('AIRW_DAT', Me%Files%InPutData, Message = Message, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadWaterAirFilesName - ModuleInterfaceWaterAir - ERR01'

        Message   ='Instant fields of water-air properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('AIRW_HDF', Me%Files%Results, Message = Message, TIME_END = Me%EndTime, &
                           Extension = 'arw',                                           &
                           MPI_ID    = GetDDecompMPI_ID(Me%ObjHorizontalGrid),&
                           DD_ON     = GetDDecompON    (Me%ObjHorizontalGrid),&
                           STAT      = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadWaterAirFilesName - ModuleInterfaceWaterAir - ERR02'


    end subroutine ReadWaterAirFilesName
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructGlobalVariables
                                                    
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL, iflag
        logical                             :: isdefined

        !Begin-----------------------------------------------------------------

        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                        &
                                   Size        = Me%Size2D,                     &
                                   WorkSize    = Me%WorkSize2D,                 &
                                   STAT        = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR10'

        call GetGeometrySize(Me%ObjGeometry,                                    &
                             Size       = Me%Size3D,                            &
                             WorkSize   = Me%WorkSize3D,                        &
                             STAT       = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR20'

        call GetComputeTimeLimits(Me%ObjTime,                                   &
                                  BeginTime = Me%BeginTime,                     &
                                  EndTime   = Me%EndTime,                       &
                                  STAT      = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR30'

        !Actualize the time
        Me%ActualTime = Me%BeginTime


        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ObjTime, Me%ExternalVar%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR40'
        

        !Altitude in km
        call GetData(Me%Altitude,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword      ='ALTITUDE',                                          &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleInterfaceWaterAir',                          &
                     Default      = 0.0,                                                &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR50'

        !Mortimer's altitude correction
        Me%AltitudeCorrection = (1. - (Me%Altitude/1000.)/44.3)**5.25

        call GetData(Me%AerationEquation,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword      ='AERATION_METHOD',                                   &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleInterfaceWaterAir',                          &
                     Default      = Gelda_et_al_1996,                                   &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR60'
        
        call GetData(Me%CO2AerationEquation,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword      ='CO2_AERATION_METHOD',                               &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleInterfaceWaterAir',                          &
                     Default      = Borges_et_al_2004,                                  &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR70'
!COARE  checks if the user wants to compute surfaceheatbudget with COARE based algorithm
        call GetData(Me%COARE,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword      ='USE_COARE',                                         &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleInterfaceWaterAir',                          &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR80'
    
i1:    if (Me%COARE) then
        !checks whether the user wants to compute warmlayer 
        call GetData(Me%Jwarm,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword      ='COMPUTE_WARM_LAYER',                                &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleInterfaceWaterAir',                          &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR90'
        !checks whether the user wants to compute cool skin
        call GetData(Me%Jcool,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword      ='COMPUTE_COOL_SKIN',                                 &
                     SearchType   = FromFile,                                           &
                     ClientModule = 'ModuleInterfaceWaterAir',                          &
                     Default      = 0.,                                                 &
                     STAT         = STAT_CALL)            
        if (STAT_CALL  /= SUCCESS_)                                                     &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR100'


        call GetWindHeight(AtmosphereID  = Me%ObjAtmosphere,                            &
                                    Height       = Me%WindHeight,                       &
                                    isdefined    = isdefined,                           &
                                    STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR115'
        if (.not. isdefined) then
            write(*,*) 'WIND_MEASUREMENT_HEIGHT not defined'
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR120'
        endif

        call GetAirMeasurementHeight(AtmosphereID  = Me%ObjAtmosphere,                  &
                                    AirHeight    = Me%AirMeasurementHeight,             &
                                    isdefined    = isdefined,                           &
                                    STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR125'
        if (.not. isdefined) then
            write(*,*) 'AIR_MEASUREMENT_HEIGHT not defined'
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR130'
        endif
        
! --------------------------------Allocation of COARE Variables--------------------------------       
    
            call AllocateCOAREVariables
        
    endif i1
       
    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------
    subroutine AllocateCOAREVariables

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
      allocate (Me%SurfaceTemperature(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR135'
        
      allocate (Me%Fxp(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR140'
      call SetMatrixValue(Me%Fxp, Me%Size2D, 0.5)

      allocate (Me%WarmLayerThickness(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR145'
      call SetMatrixValue(Me%WarmLayerThickness, Me%Size2D, 19.)

      allocate (Me%Tau_ac(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR150'
      call SetMatrixValue(Me%Tau_ac, Me%Size2D, 0.)
      
      allocate (Me%Tau(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR155'
      call SetMatrixValue(Me%Tau, Me%Size2D, 0.)
        
      allocate (Me%AccumulatedEnergy(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR160'
      call SetMatrixValue(Me%AccumulatedEnergy, Me%Size2D, 0.)
        
      allocate (Me%WarmLayerTempDiff(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR165'
      call SetMatrixValue(Me%WarmLayerTempDiff, Me%Size2D, 0.)
    
      allocate (Me%Al(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR170'
      call SetMatrixValue(Me%Al, Me%Size2D, 0.)

      allocate (Me%RainFlux(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR175'
      call SetMatrixValue(Me%RainFlux, Me%Size2D, 0.)
      
    allocate (Me%LastLatentHeat(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR180'
        call SetMatrixValue(Me%LastLatentHeat, Me%Size2D, 0.)
      
    allocate (Me%LastSensibleHeat(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR185'
        call SetMatrixValue(Me%LastSensibleHeat, Me%Size2D, 0.)        
        
        
        
    end subroutine AllocateCOAREVariables
    
    
    subroutine ConstructRugosity

        !Local------------------------------------------------------------------
        integer                        :: ClientNumber, ILB, IUB, JLB, JUB 
        integer                        :: STAT_CALL, iflag
        logical                        :: BlockFound
        
        !Begin------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        Me%Rugosity%ON            = .false.

        call GetData(Me%Rugosity%Scalar,                                &
                     Me%ObjEnterData, iflag,                            &
                     keyword      ='RUGOSITY',                          &
                     SearchType   = FromFile,                           &
                     ClientModule = 'ModuleInterfaceWaterAir',          &
                     Default      = 0.0025,                             &
                     STAT         = STAT_CALL)

        if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR05'

i2:     if (iflag == 1) then

            nullify (Me%Rugosity%Field)
            allocate(Me%Rugosity%Field(ILB:IUB, JLB:JUB))
            Me%Rugosity%Field(:,:) = Me%Rugosity%Scalar
            Me%Rugosity%Constant      = .true.
            Me%Rugosity%ON            = .true.
        
! Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Me%Rugosity%WavesFunction = NoWave
!!!!!!!!!!!!!
        else i2

            !Reads the Atmosphere rugosity
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &
                                        rugosity_begin, rugosity_end, BlockFound,           &
                                        STAT = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR10'

i1:         if(BlockFound)then

                nullify (Me%Rugosity%Field)
                allocate(Me%Rugosity%Field(ILB:IUB, JLB:JUB))
                Me%Rugosity%Field(:,:) = FillValueReal

                Me%Rugosity%ON = .true.

                call ConstructFillMatrix  (PropertyID           = Me%Rugosity%ID,           &
                                           EnterDataID          = Me%ObjEnterData,          &
                                           TimeID               = Me%ObjTime,               &
                                           HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                           ExtractType          = FromBlock,                &
                                           PointsToFill2D       = Me%ExtWater%WaterPoints2D,&
                                           Matrix2D             = Me%Rugosity%Field,        &
                                           TypeZUV              = TypeZ_,                   &
                                           ClientID             = ClientNumber,             &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR20'


                call GetDefaultValue(Me%Rugosity%ID%ObjFillMatrix, Me%Rugosity%Scalar, STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR30'

                call GetIfMatrixRemainsConstant(FillMatrixID    = Me%Rugosity%ID%ObjFillMatrix,&
                                                RemainsConstant = Me%Rugosity%Constant, &
                                                STAT            = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR40'


                call GetData(Me%Rugosity%WavesFunction,                                 &
                             Me%ObjEnterData, iflag,                                    &
                             keyword      ='WAVES_FUNCTION',                            &
                             SearchType   = FromBlock,                                  &
                             ClientModule = 'ModuleInterfaceWaterAir',                  &
! Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             !Default      = .false.,                                &
                             Default      = NoWave,                                 &
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                            &
                    stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR60'

! Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! if (Me%Rugosity%WavesFunction) then
                if (Me%Rugosity%WavesFunction == SurfaceRugosityFromHS  .or.        &
                    Me%Rugosity%WavesFunction == SurfaceRugosityFromHSW   ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    

                    call GetData(Me%Rugosity%WavesRelation,                             &
                                 Me%ObjEnterData, iflag,                                &
                                 keyword      ='WAVES_RELATION',                        &
                                 SearchType   = FromBlock,                              &
                                 ClientModule = 'ModuleInterfaceWaterAir',              &
                                 Default      = 1.,                                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                        &
                        stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR70'

                    if (Me%ObjWaves == 0) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR80'

                    if (Me%Rugosity%Constant) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR90'

                    call ModifyRugosity
            
                endif

                if(.not. Me%Rugosity%ID%SolutionFromFile) then

                    call KillFillMatrix(Me%Rugosity%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR50'
                endif
          

            else i1

                Me%Rugosity%ON = .false.

            endif i1

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructRugosity - ModuleInterfaceWaterAir - ERR100'

        endif i2


    end subroutine ConstructRugosity

    !--------------------------------------------------------------------------
    
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 16/12/2011 

    subroutine ConstructWaveFluxTKE

        !Local------------------------------------------------------------------
        integer                        :: ILB, IUB, JLB, JUB 
        integer                        :: STAT_CALL, iflag
        
        !Begin------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        Me%WaveFluxTKE%ON                = .false.

        ! If WAVEMODEL_TKE_FLUX : 1, the surface TKE flux read in Module Waves will be imported in the 
        ! Module InterfaceWaterAir by subroutine ModifyWaveFluxTKE
        ! If WAVEMODEL_TKE_FLUX : 2, the same + the surface shear velocity is computed in a manner consistent
        ! with the GLM theory. To be used with option WAVE_FORCING_3D: 2 in module Hydrodynamics. 
        call GetData(Me%WaveFluxTKE%ON,                                 &
                     Me%ObjEnterData, iflag,                            &
                     keyword      ='WAVEMODEL_TKE_FLUX',                &
                     SearchType   = FromFile,                           &
                     ClientModule = 'ModuleInterfaceWaterAir',          &
                     Default      = .false.,                            &
                     STAT         = STAT_CALL)

        ! Allocate field
        if (Me%WaveFluxTKE%ON) then
            nullify (Me%WaveFluxTKE%Field)
            allocate(Me%WaveFluxTKE%Field(ILB:IUB, JLB:JUB))
            Me%WaveFluxTKE%Field(:,:) = FillValueReal
        endif

        ! Importation of the field from Module Waves
        if (Me%WaveFluxTKE%ON) then
            call ModifyWaveFluxTKE
        endif

    end subroutine ConstructWaveFluxTKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                                 :: ClientNumber
        integer                                 :: STAT_CALL
        logical                                 :: BlockFound
        type (T_Property), pointer              :: PropertyX        => null()
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
                        stop 'Construct_PropertyList - ModuleInterfaceWaterAir - ERR01'

                    exit do1    !No more blocks
                endif cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_PropertyList - ModuleInterfaceWaterAir - ERR02'
            else cd1
                stop 'Construct_PropertyList - ModuleInterfaceWaterAir - ERR03'
            endif cd1
        end do do1


        !Verifies Wind consistence. Now is done by vectorial prop
        call SearchProperty(PropertyX, WindStressX_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            write(*,*) 'Vectorial Property wind stress is now defined in one single block'
            write(*,*) 'See Documentation on how to implement it'
            stop 'Construct_PropertyList - ModuleInterfaceWaterAir . ERR04'                     
        endif
        call SearchProperty(PropertyX, WindStressY_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            write(*,*) 'Vectorial Property wind stress is now defined in one single block'
            write(*,*) 'See Documentation on how to implement it'
            stop 'Construct_PropertyList - ModuleInterfaceWaterAir . ERR05'                     
        endif        

        !------------------------------------------------------------------------

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------

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
                if (PrintWarning) write (*,*)'Property Not Found in Module InterfaceWaterAir ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchProperty

    
    !----------------------------------------------------------------------    

    subroutine Construct_Property(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty
        integer                         :: ClientNumber

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        integer                         :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB 
        IUB = Me%Size2D%IUB 
        JLB = Me%Size2D%JLB 
        JUB = Me%Size2D%JUB 
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_Property - ModuleInterfaceWaterAir - ERR01' 

        nullify(NewProperty%Prev,NewProperty%Next)

        call ConstructPropertyID        (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call Construct_PropertyEvolution(NewProperty)

        call Construct_PropertyValues   (NewProperty, ClientNumber)

        call Construct_PropertyOutPut   (NewProperty)

    end subroutine Construct_Property

    
    !--------------------------------------------------------------------------


    subroutine Construct_PropertyValues(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property),   pointer                 :: NewProperty
        integer                                     :: ClientNumber

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        integer                                     :: SizeILB, SizeIUB
        integer                                     :: SizeJLB, SizeJUB 
        
        !----------------------------------------------------------------------

        SizeILB = Me%Size2D%ILB
        SizeIUB = Me%Size2D%IUB
        SizeJLB = Me%Size2D%JLB
        SizeJUB = Me%Size2D%JUB

        
!~         if (Check_Vectorial_Property(NewProperty%ID%IDNumber)) then 
        if (NewProperty%ID%IsVectorial) then
                        
            allocate(NewProperty%FieldU(SizeILB:SizeIUB, SizeJLB:SizeJUB))
            NewProperty%FieldU(:,:) = FillValueReal
            
            allocate(NewProperty%FieldV(SizeILB:SizeIUB, SizeJLB:SizeJUB))
            NewProperty%FieldV(:,:) = FillValueReal
            
            allocate(NewProperty%FieldX(SizeILB:SizeIUB, SizeJLB:SizeJUB))
            NewProperty%FieldX(:,:) = FillValueReal
            
            allocate(NewProperty%FieldY(SizeILB:SizeIUB, SizeJLB:SizeJUB))
            NewProperty%FieldY(:,:) = FillValueReal                    
            
            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,           &
                                       EnterDataID          = Me%ObjEnterData,          &
                                       TimeID               = Me%ObjTime,               &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                       ExtractType          = FromBlock,                &
                                       PointsToFill2D       = Me%ExtWater%WaterPoints2D,&
                                       Matrix2DU            = NewProperty%FieldU,       &
                                       Matrix2DV            = NewProperty%FieldV,       &
                                       Matrix2DX            = NewProperty%FieldX,       &
                                       Matrix2DY            = NewProperty%FieldY,       &  
                                       TypeZUV              = TypeZ_,                   &
                                       ClientID             = ClientNumber,             &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR00'
        
        else
        
            allocate (NewProperty%Field (SizeILB:SizeIUB, SizeJLB:SizeJUB))

            NewProperty%Field(:,:) = FillValueReal
                
            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,           &
                                       EnterDataID          = Me%ObjEnterData,          &
                                       TimeID               = Me%ObjTime,               &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,     &
                                       ExtractType          = FromBlock,                &
                                       PointsToFill2D       = Me%ExtWater%WaterPoints2D,&
                                       Matrix2D             = NewProperty%Field,        &
                                       TypeZUV              = TypeZ_,                   &
                                       ClientID             = ClientNumber,             &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR10'

        endif        
        
        call GetIfMatrixRemainsConstant(FillMatrixID    = NewProperty%ID%ObjFillMatrix,&
                                        RemainsConstant = NewProperty%Constant,        &
                                        STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR20'

        if(.not. NewProperty%ID%SolutionFromFile)then

            call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR30'
        endif

        if (NewProperty%ID%IDNumber == WindStress_) then

            !Shear Coefficient at the Atmosphere
            call GetData(Me%DefineCDWIND,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'DEFINE_CDWIND',                                &  
                         default      = .false.,                                        &
                         ClientModule = 'ModuleInterfaceWaterAir',                      &
                         SearchType   = FromBlock,                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR50'
                
            if (iflag == 1) then
            
                if (Me%DefineCDWIND) then
                
                    Me%CDWINDmethod = Constant
                
                else
                
                    Me%CDWINDmethod = WindFunction                
                
                endif
            
            else

                !method to define the wind drag coefficient 
                !1 - constant, 2 - wind function (a*Wind + b), 3 - wind shear velocity 
                call GetData(Me%CDWINDmethod,                                               &
                             Me%ObjEnterData, iflag,                                        &
                             keyword      = 'CDWIND_METHOD',                                &  
                             default      = WindFunction,                                   &
                             ClientModule = 'ModuleInterfaceWaterAir',                      &
                             SearchType   = FromBlock,                                      &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR60'
            
            endif                

            !Shear Coefficient at the Atmosphere
            ! Mellor, Introduction to Physical Oceanography, p52 (1996)
            call GetData(Me%CDWIND,                                                     &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'CDWIND',                                       &  
                         ClientModule = 'ModuleInterfaceWaterAir',                      &
                         default      = 0.0015,                                         &
                         SearchType   = FromBlock,                                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR06'
        
        endif
        
        if (NewProperty%ID%IDNumber == LatentHeat_) then
            !1 - Uses AirTemperature. 2. Calculates SurfaceAirTemperature based on airtemperature and wind velocity
            !Method 2 created to correct high latent heat losses when WaterTemperature >> AirTemperature
            call GetData (NewProperty%SVPMethod,                                             &
                Me%ObjEnterData, iflag,                                            &
                Keyword        = 'SVP_METHOD',                                     &
                Default        = 1,                                                &
                SearchType     = FromBlock,                                        &
                ClientModule   = 'ModuleInterfaceWaterAir',                        &
                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR12'
            
            if (NewProperty%SVPMethod == 2) then
            call GetData (NewProperty%C1,                                             &
                Me%ObjEnterData, iflag,                                            &
                Keyword        = 'KS_LH',                                     &
                Default        = 1,                                                &
                SearchType     = FromBlock,                                        &
                ClientModule   = 'ModuleInterfaceWaterAir',                        &
                STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR15'
            
                call GetData (NewProperty%C2,                                          &
                    Me%ObjEnterData, iflag,                                            &
                    Keyword        = 'POWER_LH',                                       &
                    Default        = 1,                                              &
                    SearchType     = FromBlock,                                        &
                    ClientModule   = 'ModuleInterfaceWaterAir',                        &
                    STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR15'
                
            
            endif
                  
        endif
        
        if (NewProperty%ID%IDNumber == SurfaceRadiation_) then
                      
            call GetData(Me%ReflectionCoef,                                             &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'ALBEDO',                                       &  
                         default      = 0.0,                                            &
                         ClientModule = 'ModuleInterfaceWaterAir',                      &
                         SearchType   = FromBlock,                                      &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR40'
            
            nullify(Me%SurfaceAlbedo)
            allocate (Me%SurfaceAlbedo(Me%Size2D%ILB:Me%Size2D%IUB, Me%Size2D%JLB:Me%Size2D%JUB), STAT=STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'ConstructGlobalVariables - ModuleInterfaceWaterAir - ERR131'
            
        endif
         

    end subroutine Construct_PropertyValues


   !--------------------------------------------------------------------------


    subroutine Construct_PropertyEvolution(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        real                                    :: ModelDT

        !Local-----------------------------------------------------------------
        integer                                 :: iflag
        real                                    :: ErrorAux, auxFactor, DTaux
        logical                                 :: VariableDT

        !----------------------------------------------------------------------


        !Time Step if the property field is variable in time
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR08'

            call GetVariableDT (Me%ObjTime, VariableDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR09'
    
            if (VariableDT) then
                
                NewProperty%Evolution%DTInterval       = ModelDT

            else

            !<BeginKeyword>
                !Keyword          : DT_INTERVAL
                !<BeginDescription>       
                   ! In the future a GET_DT_Hydro to a hydrodynamic module should be done to know the DT value  
                   ! that is been used to computethe the velocity field
                   ! By default the DTinterval is equal to the time step of the hydrodynamic model. In this case 
                   ! is admitted that the hydrodynamic model is computing the velocity field using a ADI method.
                !<EndDescription>
                !Type             : Real 
                !Default          : Time step of hydrodynamic model
                !File keyword     : DISPQUAL
                !Multiple Options : Do not have
                !Search Type      : FromBlock
                !Begin Block      : <beginproperty>
                !End Block        : <endproperty>
            !<EndKeyword>

            call GetData(NewProperty%Evolution%DTInterval,                          &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType     = FromBlock,                                &
                         keyword        = 'DT_INTERVAL',                            &
                         Default        = ModelDT,                                  &
                         ClientModule   = 'ModuleInterfaceWaterAir',                &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
                stop 'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR10'


            if (NewProperty%evolution%DTInterval < (ModelDT)) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop       'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR11'

            elseif (NewProperty%evolution%DTInterval > (ModelDT)) then

                !Property DT  must be a multiple of the ModelDT
                auxFactor = NewProperty%evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) ' Time step error.'
                    stop       'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR12'
                endif

                    ! Run period in seconds
                    DTaux = Me%EndTime - Me%ExternalVar%Now

                    !The run period   must be a multiple of the Property DT
                    auxFactor = DTaux / NewProperty%evolution%DTInterval

                    ErrorAux = auxFactor - int(auxFactor)
                    if (ErrorAux /= 0) then
                        write(*,*) 
                        write(*,*) ' Time step error.'
                        stop       'Construct_PropertyEvolution - ModuleInterfaceWaterAir - ERR13'
                    endif

                endif

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

        !<BeginKeyword>
            !Keyword          : OUTPUT_HDF
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a HDF format file for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%OutputHDF,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'OUTPUT_HDF',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceWaterAir',                        &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceWaterAir - ERR01'
           
        !<BeginKeyword>
            !Keyword          : TIME_SERIE
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a time serie for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%TimeSerie,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'TIME_SERIE',                                     &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceWaterAir',                        &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceWaterAir - ERR02'


        !<BeginKeyword>
            !Keyword          : BOX_TIME_SERIE
            !<BeginDescription>       
                ! Checks out if the user pretends to write a time serie inside each box for this property
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : DISPQUAL
            !Multiple Options : 1 (.true.) , 0 (.false.) 
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%BoxTimeSerie,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'BOX_TIME_SERIE',                                 &
                     Default        = .false.,                                          &
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleInterfaceWaterAir',                        &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleInterfaceWaterAir - ERR03'

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
        endif 

        !----------------------------------------------------------------------

    end subroutine Add_Property 

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: HDF5_CREATE
        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize2D%ILB 
        WorkIUB = Me%WorkSize2D%IUB 
        WorkJLB = Me%WorkSize2D%JLB 
        WorkJUB = Me%WorkSize2D%JUB 

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, trim(Me%Files%Results)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR02'
        
        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR03'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",               &
                              Array2D = Me%ExtWater%Bathymetry,                     &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR04'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",            &
                              Array2D = Me%ExtWater%WaterPoints2D,                  &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR05'

        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
           stop 'Open_HDF5_OutPut_File - ModuleInterfaceWaterAir - ERR06'
       
        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty
        logical                                     :: OutputON
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        nullify(Me%OutPut%OutTime)

        OutputON = OFF

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))
            
            if(CurrentProperty%OutputHDF) OutputON = ON
            CurrentProperty => CurrentProperty%Next

        enddo

        if(OutputON)then

            call GetOutPutTime(Me%ObjEnterData,                                         &
                               CurrentTime      = Me%ExternalVar%Now,                   &
                               EndTime          = Me%EndTime,                           &
                               keyword          = 'OUTPUT_TIME',                        &
                               SearchType       = FromFile,                             &
                               OutPutsTime      = Me%OutPut%OutTime,                    &
                               OutPutsOn        = Me%OutPut%Yes,                        &
                               OutPutsNumber    = Me%OutPut%Number,                     &                               
                               STAT             = STAT_CALL)
                                                
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructGlobalOutput - ModuleInterfaceWaterAir - ERR01' 

            if (Me%OutPut%Yes) then

                Me%OutPut%NextOutPut = 1

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructGlobalOutput - ModuleInterfaceWaterAir - ERR02'
            endif 

        endif

    end subroutine ConstructGlobalOutput

    !------------------------------------------------------------------------
    

    subroutine Construct_Sub_Modules

        !Local-----------------------------------------------------------------
        type (T_Property),           pointer                 :: PropertyX            

        !----------------------------------------------------------------------

        Me%Coupled%TimeSerie%NumberOfProperties            = 0
        Me%Coupled%BoxTimeSerie%NumberOfProperties         = 0

        PropertyX => Me%FirstProperty

do1 :   do while (associated(PropertyX))

            if (PropertyX%TimeSerie) then
                Me%Coupled%TimeSerie%NumberOfProperties             = &
                Me%Coupled%TimeSerie%NumberOfProperties             + 1
                Me%Coupled%TimeSerie%Yes                            = ON
            endif

            if (PropertyX%BoxTimeSerie) then
                
!~                 if (Check_Vectorial_Property(PropertyX%ID%IDNumber)) then
                if (PropertyX%ID%IsVectorial) then
                    write(*,*) 'No box time serie available yet to vectorial properties'
                    write(*,*) 'Property : ', trim(PropertyX%ID%Name)
                    stop 'Construct_Sub_Modules - ModuleInterfaceWaterAir - ERR01'
                endif
                
                Me%Coupled%BoxTimeSerie%NumberOfProperties          = &
                Me%Coupled%BoxTimeSerie%NumberOfProperties          + 1
                Me%Coupled%BoxTimeSerie%Yes                         = ON
            endif

            PropertyX=>PropertyX%Next

        end do do1

        if(Me%Coupled%TimeSerie%Yes)then

            call Construct_Time_Serie
             
        endif
        
        if(Me%Coupled%BoxTimeSerie%Yes)then

            call StartOutputBoxFluxes
               
        endif

    end subroutine Construct_Sub_Modules


    !--------------------------------------------------------------------------


    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
         
        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber  
        type(T_Property), pointer                           :: PropertyX
        integer                                             :: nProperties
        character(len=PathLength)                           :: TimeSerieLocationFile
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        character(len=StringLength)                         :: TimeSerieName
        type (T_Polygon), pointer                           :: ModelDomainLimit

        !----------------------------------------------------------------------

        !First checks out how many properties will have time series
        PropertyX   => Me%FirstProperty
        nProperties =  0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                nProperties = nProperties + 1
!~                 if (Check_Vectorial_Property(PropertyX%ID%IDNumber)) then
                if (PropertyX%ID%IsVectorial) then
                    nProperties = nProperties + 1  !x and y comp
                endif
            endif
            PropertyX=>PropertyX%Next
        enddo

        if (nProperties > 0) then

            !Allocates PropertyList
            allocate(PropertyList(nProperties), STAT = STAT_CALL)
            if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR10'

            !Fills up PropertyList
            PropertyX   => Me%FirstProperty
            nProperties =  0
            do while (associated(PropertyX))
                if (PropertyX%TimeSerie) then
                    
!~                     if (Check_Vectorial_Property(PropertyX%ID%IDNumber)) then
                    if (PropertyX%ID%IsVectorial) then
                        nProperties = nProperties + 1
                        PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name)//" X")
                        nProperties = nProperties + 1
                        PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name)//" Y")     
                    else
                        nProperties = nProperties + 1
                        PropertyList(nProperties) = trim(adjustl(PropertyX%ID%name))
                    endif
                    
                endif
                PropertyX=>PropertyX%Next
            enddo


            call GetData(TimeSerieLocationFile,                                         &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'TIME_SERIE_LOCATION',                          &
                         ClientModule = 'ModuleWaterProperties',                        &
                         Default      = Me%Files%InputData,                             &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR20' 

            call GetGridOutBorderPolygon(HorizontalGridID = Me%ObjHorizontalGrid,       &
                                         Polygon          = ModelDomainLimit,           &
                                         STAT             = STAT_CALL)           
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR25' 

            !Constructs TimeSerie
            call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                            &
                                TimeSerieLocationFile,                                  &
                                PropertyList, "sri",                                    &
                                WaterPoints3D = Me%ExtWater%WaterPoints3D,              &
                                ModelName     = Me%ModelName,                           & 
                                ModelDomain   = ModelDomainLimit,                       &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR30'
            
            call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     Polygon          = ModelDomainLimit,               &
                                     STAT             = STAT_CALL)                          
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleInterfaceWaterAir - ERR35'            

            !Deallocates PropertyList
            deallocate(PropertyList, STAT = STAT_CALL)
            if (STAT_CALL /= 0) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR40'


            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR50'

            do dn = 1, TimeSerieNumber
            
                call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleInterfaceWaterAir - ERR60'
                
                if (IgnoreOK) cycle

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleInterfaceWaterAir - ERR70'
                
                call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleInterfaceWaterAir - ERR80'  
                                                          
                if (CoordON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= OUT_OF_BOUNDS_ERR_) then
                        stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR90'
                    endif                            

                    !if (STAT_CALL == OUT_OF_BOUNDS_ERR_ .or. Id < 0 .or. Jd < 0) then
                
                    !    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    !    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR70'

                    !    if (IgnoreOK) then
                    !        cycle
                    !    else
                    !        stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR80'
                    !    endif

                    !endif


                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleInterfaceWaterAir - ERR90'
                endif

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          LocalizationI   = Id,                             &
                                          LocalizationJ   = Jd,                             & 
                                          STAT     = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleInterfaceWaterAir - ERR120'

                if (Me%ExtWater%WaterPoints3D(Id, Jd, Me%WorkSize3D%KUB) /= WaterPoint) then
                    
                     write(*,*) 'Time Serie in a land cell - ',trim(TimeSerieName),' - ',   &
                                 trim(Me%ModelName), ' ModuleInterfaceWaterAir'
                     write(*,*) 'ConstructTimeSerie - ModuleInterfaceWaterAir - WRN120'
                     
                endif

            enddo

        endif
        
    end subroutine Construct_Time_Serie

    !----------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        integer                                             :: ILB, IUB, JLB, JUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        type(T_Property),                           pointer :: PropertyX
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList
        integer                                             :: nScalars, n

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        !<BeginKeyword>
            !Keyword          : BOXFLUXES
            !<BeginDescription>       
            ! This keyword have two functions if exist fluxes between boxes are compute 
            ! and the value read is the name file where the boxes are defined
            !
            !<EndDescription>
            !Type             : Character 
            !Default          : Do not have
            !File keyword     : SEDPROP
            !Multiple Options : Do not have
            !Search Type      : From File
        !<EndKeyword>
        
        call GetData(Me%Files%BoxesFile,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule = 'ModuleInterfaceWaterAir',                      &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR01'
        if (iflag .EQ. 0)                                                           &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR02'    
        
        inquire(File = Me%Files%BoxesFile, Exist = exist)
        if (exist) then
            inquire(File = Me%Files%BoxesFile, Opened  = Opened)
            if (opened) then
                write(*,*    ) 
                write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Files%BoxesFile))
                write(*,*    ) 'Already opened.'
                stop           'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR03'    
            endif
        else
            write(*,*) 
            write(*,*)     'Could not find the boxes file.'
            write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
            stop           'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR04'    
        endif

        nScalars = Me%Coupled%BoxTimeSerie%NumberOfProperties 

        allocate(ScalarOutputList(nScalars), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR06'

        PropertyX  => Me%FirstProperty
        n = 0
        do while (associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then
                n = n + 1
                ScalarOutputList(n) = trim(PropertyX%ID%name)
            endif 

            PropertyX=>PropertyX%Next
        end do

        call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &
                         TimeID             = Me%ObjTime,                   &
                         HorizontalGridID   = Me%ObjHorizontalGrid,         &
                         BoxesFilePath      = Me%Files%BoxesFile,           &
                         ScalarOutputList   = ScalarOutputList,             &
                         WaterPoints2D      = Me%ExtWater%WaterPoints2D,    &
                         STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR07'

        deallocate(ScalarOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR07'

        allocate(Me%Scalar2D(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleInterfaceWaterAir - ERR12'
        Me%Scalar2D(:,:) = 0.

    end subroutine StartOutputBoxFluxes

    
    !----------------------------------------------------------------------

    
    subroutine CheckInternalOptions

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        type (T_Property), pointer                          :: PropertyY

        !----------------------------------------------------------------------
     
        call Search_Property(PropertyX, LatentHeat_,   STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%LatentHeat             = ON

        call Search_Property(PropertyX, SensibleHeat_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SensibleHeat           = ON

        call Search_Property(PropertyX, Evaporation_,  STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%Evaporation            = ON

        call Search_Property(PropertyX, NetLongWaveRadiation_,         STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%NetLongWaveRadiation     = ON

        call Search_Property(PropertyX, UpwardLongWaveRadiation_,      STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%UpwardLongWaveRadiation   = ON

        call Search_Property(PropertyX, DownwardLongWaveRadiation_,    STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%DownwardLongWaveRadiation = ON

        call Search_Property(PropertyX, NonSolarFlux_,           STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%NonSolarFlux           = ON

        call Search_Property(PropertyX, OxygenFlux_,             STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%OxygenFlux             = ON

        call Search_Property(PropertyX, CarbonDioxideFlux_,      STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%CarbonDioxideFlux      = ON
        
        call Search_Property(PropertyX, SpecificOxygenFlux_,     STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SpecificOxygenFlux     = ON

        call Search_Property(PropertyX, SpecificCarbonDioxideFlux_,  STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SpecificCarbonDioxideFlux  = ON
        
        call Search_Property(PropertyX, AmmoniaFlux_,             STAT = STAT_CALL)  !LLP
        if (STAT_CALL == SUCCESS_) Me%IntOptions%AmmoniaFlux             = ON
        
        call Search_Property(PropertyX, NitrateFlux_,             STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%NitrateFlux             = ON
        
        call Search_Property(PropertyX, TurbulentKineticEnergy_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%TurbulentKineticEnergy = ON

        call Search_Property(PropertyX, SurfaceRadiation_,       STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SurfaceRadiation       = ON
            
        call Search_Property(PropertyY, Albedo_,                 STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%Albedo                 = ON
            
        call Search_Property(PropertyX, SurfaceWaterFlux_,       STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%SurfaceWaterFlux       = ON
        
        call Search_Property(PropertyX, WindStress_,            STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) Me%IntOptions%WindStress        = ON


    end subroutine CheckInternalOptions


    !--------------------------------------------------------------------------
    !Module interface Water Air must provide to all other modules the interface fluxes.
    !If Interface fluxes are not specified in this module but are specified in another
    !module, model stops and writes a warning to the screen
    !--------------------------------------------------------------------------
    subroutine CheckOptionsWater

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX

        !----------------------------------------------------------------------
        
        !Check WaterProperties  
        call GetWaterPropertiesAirOptions(WaterPropertiesID    = Me%ObjWaterProperties,                & 
                                          TemperatureFluxYes   = Me%ExtOptions%HeatFluxYes,            &
                                          OxygenFluxYes        = Me%ExtOptions%OxygenFluxYes,          &
                                          CarbonDioxideFluxYes = Me%ExtOptions%CarbonDioxideFluxYes,   &
                                          AmmoniaFluxYes       = Me%ExtOptions%AmmoniaFluxYes,         &  !LLP
                                          NitrateFluxYes       = Me%ExtOptions%NitrateFluxYes,         &
                                          WQMYes               = Me%ExtOptions%WQMYes,                 &
                                          T90VariableYes       = Me%ExtOptions%T90VariableYes,         &   
                                          STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR10'
 
        !Checks Hydrodynamic 
        call GetHydrodynamicAirOptions(HydrodynamicID      = Me%ObjHydrodynamic,                &
                                       SurfaceWaterFluxYes = Me%ExtOptions%SurfaceWaterFluxYes, &
                                       WindYes             = Me%ExtOptions%HydrodynamicWindYes, &
                                       AtmPressureYes      = Me%ExtOptions%HydrodynamicAtmPressureYes, &
                                       MslpYes             = Me%ExtOptions%HydrodynamicMslpYes, &
                                       STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR20'

#ifndef _LAGRANGIAN_
        if(Me%ObjLagrangian /= 0)then
#ifdef  _LAGRANGIAN_GLOBAL_                                         
        
            !Checks Lagrangian 
            call GetLagrangianAirOptionsGlobal(LagrangianID  = Me%ObjLagrangian,         &
                                         Oil           = Me%ExtOptions%OilYes,           &
                                         HNS           = Me%ExtOptions%HNSYes,           &
                                         Wind          = Me%ExtOptions%LagrangianWindYes,&
                                         WaterQuality  = Me%ExtOptions%LagrangianWQMYes, &
                                         T90Variable   = Me%ExtOptions%LagrangianT90Yes, &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR30'
#else
            call GetLagrangianAirOptions(LagrangianID  = Me%ObjLagrangian,               &
                                         Oil           = Me%ExtOptions%OilYes,           &
                                         Wind          = Me%ExtOptions%LagrangianWindYes,&
                                         WaterQuality  = Me%ExtOptions%LagrangianWQMYes, &
                                         T90Variable   = Me%ExtOptions%LagrangianT90Yes, &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR35'
#endif         
        endif
#endif

#ifndef _WAVES_
        if(Me%ObjWaves    /= 0) then
            if(AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                Me%ExtOptions%WavesWindYes     = ON
            else
                if (AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityX_) .and.    &
                    AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityY_)) then
                    Me%ExtOptions%WavesWindYes = ON
                endif                    
            endif
        endif            
#endif
        if(Me%ObjTurbGOTM /= 0)Me%ExtOptions%GOTMWindShearVelocityYes = ON

        
        if (Me%ExtOptions%HeatFluxYes) then

            call Search_Property(PropertyX, NonSolarFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR40'

            if(.not. PropertyX%ID%SolutionFromFile .and. .not. PropertyX%Constant) then

                call Search_Property(PropertyX, LatentHeat_, .true., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR50'

                call Search_Property(PropertyX, SensibleHeat_, .true., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR60'

                call Search_Property(PropertyX, NetLongWaveRadiation_, .true., STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR70'

                if(.not. PropertyX%ID%SolutionFromFile .and. .not. PropertyX%Constant) then

                    call Search_Property(PropertyX, UpwardLongWaveRadiation_, .true., STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR80'

                    call Search_Property(PropertyX, DownwardLongWaveRadiation_, .true., STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR90'
                endif

            endif

            call Search_Property(PropertyX, SurfaceRadiation_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR100'

            if(.not. PropertyX%ID%SolutionFromFile .and. .not. PropertyX%Constant) then

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_))then
                    write(*,*) 'Missing SolarRadiation in Module Atmosphere '
                    stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR110'
                endif
                

            endif

        endif

        if (Me%ExtOptions%OxygenFluxYes) then

            call Search_Property(PropertyX, OxygenFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR130'
            
        endif

        if (Me%ExtOptions%CarbonDioxideFluxYes) then

            call Search_Property(PropertyX, CarbonDioxideFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR140'         
           
        endif

        if (Me%ExtOptions%AmmoniaFluxYes) then    !LLP

            call Search_Property(PropertyX, AmmoniaFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR150'
            
        endif

        if (Me%ExtOptions%NitrateFluxYes) then

            call Search_Property(PropertyX, NitrateFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR150b'         
           
        endif


        if (Me%ExtOptions%SurfaceWaterFluxYes) then

            call Search_Property(PropertyX, SurfaceWaterFlux_, .true., STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR160'

            Me%ExtOptions%Precipitation = AtmospherePropertyExists(Me%ObjAtmosphere, Precipitation_)
            Me%ExtOptions%Irrigation    = AtmospherePropertyExists(Me%ObjAtmosphere, Irrigation_   )         
            

        endif



        !if (Me%ExtOptions%HydrodynamicWindYes) then
        !
        !    call Search_Property(PropertyX, WindStressX_, .true., STAT = STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR170'
        !
        !
        !    call Search_Property(PropertyY, WindStressY_, .true., STAT = STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_) stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR180'
        !
        !    if ((      PropertyX%ID%SolutionFromFile .and. .not. PropertyY%ID%SolutionFromFile) .or. &
        !        (.not. PropertyX%ID%SolutionFromFile .and.       PropertyY%ID%SolutionFromFile)) then
        !        
        !        write (*,*) 'wind stress X must be given in the same way as wind stress Y'
        !        stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR0190'
        !
        !    endif
        !
        !endif

        !Check if the Atmospheric Pressure property in the Atmosphere module exists
        if (Me%ExtOptions%HydrodynamicAtmPressureYes) then

            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AtmosphericPressure_)) then
                write(*,*) 'Missing AtmosphericPressure in Module Atmosphere'
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR200'
            endif

        endif

        !Check if the Mean Sea Level Pressure property in the Atmosphere module exists
        if (Me%ExtOptions%HydrodynamicMslpYes) then

            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, MeanSeaLevelPressure_)) then
                write(*,*) 'Missing mslp in Module Atmosphere'
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR210'
            endif

        endif

        if (Me%ExtOptions%OilYes .or. Me%ExtOptions%HNSYes) then

            !TO_DO Waves
            !Checks if one of the following Atmospheric properties exist:
            !Atmospheric Pressur or Mslp (Mean Sea Level Pressure)
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, AtmosphericPressure_)          &
                .and. .not. AtmospherePropertyExists(Me%ObjAtmosphere, MeanSeaLevelPressure_)) then

                write(*,*) 'Missing AtmosphericPressure or Mslp in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR220'

            endif

        endif

        if (Me%ExtOptions%LagrangianWindYes) then

            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR230'
            endif

            !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
            !    write(*,*) 'Missing WindVelocity in Module Atmosphere '
            !    stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR240'
            !endif

        endif

        if (Me%ExtOptions%WQMYes           .or. Me%ExtOptions%T90VariableYes .or.       &
            Me%ExtOptions%LagrangianWQMYes .or. Me%ExtOptions%LagrangianT90Yes) then
            
            call Search_Property(PropertyX, SurfaceRadiation_, .true., STAT = STAT_CALL) 
            if (STAT_CALL .ne. SUCCESS_)then
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR249'
            endif
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_))then
                write(*,*) 'Missing SolarRadiation in Module Atmosphere '
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR250'
            endif

        endif
        
        Me%IntOptions%WindShearVelAllocate = .false.
        
        if (Me%ExtOptions%GOTMWindShearVelocityYes .and. Me%IntOptions%WindStress       &
            .or. Me%CDWINDmethod == ShearVelocity) then
            
            Me%IntOptions%WindShearVelocity = .true.

            call Search_Property(PropertyX, WindShearVelocity_, STAT = STAT_CALL) 
            
            if (STAT_CALL == SUCCESS_) then
            
                Me%WindShearVelocity  => PropertyX%Field
            
            else if (STAT_CALL == NOT_FOUND_ERR_) then  
            
                allocate(Me%WindShearVelocity(Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
                Me%IntOptions%WindShearVelAllocate = .true.
                
            else if (STAT_CALL /= SUCCESS_) then
            
                stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR260'
                
            endif
            
        else
        
            Me%IntOptions%WindShearVelocity = .false.
            
        endif


    end subroutine CheckOptionsWater
    
    !--------------------------------------------------------------------------
    !Module interface Water Air needs some atmosphere properties. This subroutine
    !checks if those properties are defined in atmosphere.
    !--------------------------------------------------------------------------
    
    subroutine CheckOptionsAir

        !External--------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB
        
        !if (AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityU_))then
        if (AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
            
            Me%LocalAtm%WindVelocityU%Yes = ON
            nullify (Me%LocalAtm%WindVelocityU%Field)
            allocate(Me%LocalAtm%WindVelocityU%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%WindVelocityU%Field(:,:) = FillValueReal
            
            Me%LocalAtm%WindVelocityV%Yes = ON
            nullify (Me%LocalAtm%WindVelocityV%Field)
            allocate(Me%LocalAtm%WindVelocityV%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%WindVelocityV%Field(:,:) = FillValueReal            

        endif

        !if (AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
        !    
        !    Me%LocalAtm%WindVelocityV%Yes = ON
        !    nullify (Me%LocalAtm%WindVelocityV%Field)
        !    allocate(Me%LocalAtm%WindVelocityV%Field(ILB:IUB, JLB:JUB))
        !    Me%LocalAtm%WindVelocityV%Field(:,:) = FillValueReal
        !endif

        if (AtmospherePropertyExists(Me%ObjAtmosphere, WindDirection_))then
            
            Me%LocalAtm%WindDirection%Yes = ON
            nullify (Me%LocalAtm%WindDirection%Field)
            allocate(Me%LocalAtm%WindDirection%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%WindDirection%Field(:,:) = FillValueReal
        endif
        
        !Allocate field for the Atmospheric Pressure if present
        if (AtmospherePropertyExists(Me%ObjAtmosphere, AtmosphericPressure_)) then
           
            Me%LocalAtm%AtmosphericPressure%Yes = ON
            nullify (Me%LocalAtm%AtmosphericPressure%Field)
            allocate(Me%LocalAtm%AtmosphericPressure%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%AtmosphericPressure%Field(:,:) = FillValueReal
        endif

        !Allocate field for the Mean Sea Level Pressure if present
        if (AtmospherePropertyExists(Me%ObjAtmosphere, MeanSeaLevelPressure_)) then
           
            Me%LocalAtm%Mslp%Yes = ON
            nullify (Me%LocalAtm%Mslp%Field)
            allocate(Me%LocalAtm%Mslp%Field(ILB:IUB, JLB:JUB))
            Me%LocalAtm%Mslp%Field(:,:) = FillValueReal
        endif

        call Search_Property(PropertyX, LatentHeat_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)then
            if((.not. PropertyX%ID%SolutionFromFile) &
                .and. (.not. PropertyX%Constant)) then

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                    write(*,*) 'Missing WindVelocity in Module Atmosphere '
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR10'
                endif

                !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
                !    write(*,*) 'Missing WindVelocity in Module Atmosphere '
                !    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR20'
                !endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                    write(*,*) 'Specify Atmosphere AirTemperature to calculate Interface'
                    write(*,*) 'LatentHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR30'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, RelativeHumidity_)) then
                    write(*,*) 'Specify Atmosphere RelativeHumidity to calculate Interface'
                    write(*,*) 'LatentHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR40'
                endif
            endif

        endif

        call Search_Property(PropertyX, SensibleHeat_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_)then
            if ((.not. PropertyX%ID%SolutionFromFile) &
                 .and. (.not. PropertyX%Constant)) then

                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                    write(*,*) 'Missing WindVelocity X in Module Atmosphere '
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR50'
                endif

                !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
                !    write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
                !    stop 'CheckOptionsWater - ModuleInterfaceWaterAir - ERR60'
                !endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                    write(*,*) 'Specify Atmosphere AirTemperature to calculate Interface'
                    write(*,*) 'SensibleHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR70'
                endif
            endif
        endif

        call Search_Property(PropertyX, DownwardLongWaveRadiation_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if ((.not. PropertyX%ID%SolutionFromFile) &
                 .and. (.not. PropertyX%Constant)) then

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, RelativeHumidity_)) then
                    write(*,*) 'Specify Atmosphere RelativeHumidity to calculate Interface'
                    write(*,*) 'NetLongWaveRadiation'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR80'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, CloudCover_)) then
                    write(*,*) 'Specify Atmosphere CloudCover to calculate Interface'
                    write(*,*) 'NetLongWaveRadiation'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR90'
                endif

                if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, AirTemperature_)) then
                    write(*,*) 'Specify Atmosphere AirTemperature to calculate Interface'
                    write(*,*) 'SensibleHeat'
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR100'
                endif
            endif

        endif


        call Search_Property(PropertyX, OxygenFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR110'
            endif
            !
            !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
            !    write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
            !    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR120'
            !endif

        endif
        
        
        call Search_Property(PropertyX, CarbonDioxideFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR121'
            endif

            !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
            !    write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
            !    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR122'
            !endif

        endif


        call Search_Property(PropertyX, SpecificOxygenFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR123'
            endif

            !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
            !    write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
            !    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR124'
            !endif

        endif
        
        
        call Search_Property(PropertyX, SpecificCarbonDioxideFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                write(*,*) 'Missing WindVelocity in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR125'
            endif

            !if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocityV_))then
            !    write(*,*) 'Missing WindVelocity Y in Module Atmosphere '
            !    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR126'
            !endif

        endif


        call Search_Property(PropertyX, SurfaceRadiation_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not.AtmospherePropertyExists(Me%ObjAtmosphere, SolarRadiation_)) then
                write(*,*) 'Specify Atmosphere SolarRadiation to calculate Interface'
                write(*,*) 'SurfaceRadiation'
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR130'
            endif
        endif
        
        call Search_Property(PropertyX, AmmoniaFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, AtmospDeposReduNH4_))then
                write(*,*) 'Missing Atmospheric Deposition Reduced NH4 in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR140'
            endif

        endif
        
        call Search_Property(PropertyX, NitrateFlux_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            
             if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, AtmospDeposOxidNO3_))then
                write(*,*) 'Missing Atmospheric Deposition Oxidized NO3 in Module Atmosphere '
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR150'
            endif
        endif
        
        if (Me%COARE) then
            if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, PBLHeight_))then
                write(*,*) 'Missing PBLHeight in Module Atmosphere'
                stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR160'
            endif
        endif
        
        !If computing wind stress it is needed wind velocity in atmosphere
        call Search_Property(PropertyX, WindStress_, .false., STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            if (.not. PropertyX%Constant            .and. &
                .not. PropertyX%ID%SolutionFromFile ) then               
                if (.not. AtmospherePropertyExists(Me%ObjAtmosphere, WindVelocity_))then
                    write(*,*) 'Missing WindVelocity in Module Atmosphere '
                    stop 'CheckOptionsAir - ModuleInterfaceWaterAir - ERR170'
                endif

            endif
        endif        
        
    end subroutine CheckOptionsAir


    subroutine SetSubModulesConstructor 

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL 
        
        !Begin-----------------------------------------------------------------

        if(Me%ObjTurbGOTM /= 0)then

            if (.not. Me%Rugosity%ON)                                                   &
                stop 'SetSubModulesConstructor - ModuleInterfaceWaterAir - ERR10'

            call SetTurbGOTMSurfaceRugosity(TurbGOTMID      = Me%ObjTurbGOTM,           &
                                            SurfaceRugosity = Me%Rugosity%Field,        &
                                            STAT            = STAT_CALL)                
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'SetSubModulesConstructor - ModuleInterfaceWaterAir - ERR20'
            
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 16/12/2011
            if (Me%WaveFluxTKE%ON) then

                call SetTurbGOTMWaveSurfaceFluxTKE(TurbGOTMID      = Me%ObjTurbGOTM,                 &
                                                   WaveSurfaceFluxTKE = Me%WaveFluxTKE%Field,        &
                                                   STAT            = STAT_CALL)                
                if (STAT_CALL /= SUCCESS_)                                                           &
                    stop 'SetSubModulesConstructor - ModuleInterfaceWaterAir - ERR20a'
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        endif


    end subroutine SetSubModulesConstructor

    !----------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    subroutine Search_Property(PropertyX, PropertyXID, PrintWarning, STAT)


        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer                   , intent (IN)         :: PropertyXID
        logical,          optional, intent (IN)         :: PrintWarning
        integer,          optional, intent (OUT)        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX)) 
            if (PropertyX%ID%IDNumber == PropertyXID) then
                exit
            else
                PropertyX => PropertyX%Next                 
            endif
        end do

       !A PropertyX was found
       if (associated(PropertyX)) then
            STAT_ = SUCCESS_  
        else
            if (present(PrintWarning)) then
                if (PrintWarning) then
                    write (*,*)'Property Not Found in Module InterfaceWaterAir ',trim(GetPropertyName(PropertyXID))
                endif
            endif
            STAT_  = NOT_FOUND_ERR_  
        endif

        if (present(STAT))STAT = STAT_

    end subroutine Search_Property

    !--------------------------------------------------------------------------
    
    subroutine ReadLockExternalGlobal
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceWaterAir - ERR01'

        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadLockExternalGlobal - ModuleInterfaceWaterAir - ERR02'

    end subroutine ReadLockExternalGlobal

    
    !--------------------------------------------------------------------------

    
    subroutine ReadLockExternalWater
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExtWater%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceWaterAir - ERR01'

        call GetDensity(WaterPropertiesID = Me%ObjWaterProperties,          &
                        Density           = Me%ExtWater%Density,            &
                        CurrentTime       = Me%ExternalVar%Now,             &
                        STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalWater - ModuleInterfaceWaterAir - ERR02'

        call GetGridData(Me%ObjGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterfaceWaterAir - ERR04'

        !WaterPoints3D
        call GetWaterPoints3D(Me%ObjMap, Me%ExtWater%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleInterfaceWaterAir - ERR05'


    end subroutine ReadLockExternalWater


    !----------------------------------------------------------------------

    subroutine ReadUnlockExternalGlobal
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnlockExternalGlobal - ModuleInterfaceWaterAir - ERR01'

        call null_time(Me%ExternalVar%Now)

    end subroutine ReadUnlockExternalGlobal

    
    !----------------------------------------------------------------------


    subroutine ReadUnlockExternalWater
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Ungets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap,  Me%ExtWater%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalGlobal - ModuleInterfaceWaterAir - ERR01'

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%Density, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceWaterAir - ERR02'

        call UnGetGridData(Me%ObjGridData, Me%ExtWater%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceWaterAir - ERR03'

        call UnGetMap(Me%ObjMap, Me%ExtWater%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleInterfaceWaterAir - ERR04'

    end subroutine ReadUnlockExternalWater


    !----------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyInterfaceWaterAir(ObjInterfaceWaterAirID, LagrangianID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceWaterAirID, LagrangianID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterfaceWaterAir", "ModifyInterfaceWaterAir")

        STAT_ = UNKNOWN_

        call Ready(ObjInterfaceWaterAirID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if(Me%ObjLagrangian == 0 .and. LagrangianID /= 0)then
                Me%ObjLagrangian  = AssociateInstance(mLAGRANGIAN_, LagrangianID)
            endif

            call ReadLockExternalGlobal

            call ReadLockExternalWater

            call ModifyLocalAtmVariables

            call ModifyWaterAirFluxes

            call ModifyRugosity

! Modified by Matthias DELPEY - 16/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (Me%WaveFluxTKE%ON) then
                call ModifyWaveFluxTKE
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if(Me%Coupled%TimeSerie%Yes)            &
                call Output_TimeSeries

            if(Me%Coupled%BoxTimeSerie%Yes)         &
                call Output_BoxTimeSeries

            if(Me%OutPut%Yes)                       &
                call OutPut_Results_HDF

            call ReadUnlockExternalWater

            call ReadUnlockExternalGlobal

            call SetSubModulesModifier

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        endif

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceWaterAir", "ModifyInterfaceWaterAir")

    end subroutine ModifyInterfaceWaterAir

    
    !--------------------------------------------------------------------------


    subroutine ModifyLocalAtmVariables

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-------------------------------------------------------------

        !Begin-------------------------------------------------------------
        
        if(Me%LocalAtm%WindVelocityU%Yes .and. Me%LocalAtm%WindVelocityV%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,              &
                                       ScalarU      = Me%ExtAtm%WindVelocityU%Field, &
                                       ScalarV      = Me%ExtAtm%WindVelocityV%Field, &
                                       ID           = WindVelocity_,                 &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR01'   
            
            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%WindVelocityU%Field, Me%Size2D,         &
                                Me%ExtAtm%WindVelocityU%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindVelocityU%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR02'       
            
            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%WindVelocityV%Field, Me%Size2D,         &
                                Me%ExtAtm%WindVelocityV%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindVelocityV%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04'            
            
        endif
        
        !if(Me%LocalAtm%WindVelocityU%Yes)then
        !
        !    call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
        !                               Scalar       = Me%ExtAtm%WindVelocityU%Field,&
        !                               ID           = WindVelocityU_,               &
        !                               STAT         = STAT_CALL)
        !    if (STAT_CALL .ne. SUCCESS_)                                            &
        !        stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR01'
        !    
        !    !Memory duplication
        !    call SetMatrixValue(Me%LocalAtm%WindVelocityU%Field, Me%Size2D,         &
        !                        Me%ExtAtm%WindVelocityU%Field)
        !
        !    call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindVelocityU%Field, STAT = STAT_CALL)
        !    if (STAT_CALL .ne. SUCCESS_)                                            &
        !        stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR02'
        !
        !
        !endif
        !
        !
        !if(Me%LocalAtm%WindVelocityV%Yes)then
        !    call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
        !                               Scalar       = Me%ExtAtm%WindVelocityV%Field,&
        !                               ID           = WindVelocityV_,               &
        !                               STAT         = STAT_CALL)
        !    if (STAT_CALL .ne. SUCCESS_)                                            &
        !        stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR03'
        !
        !    !Memory duplication
        !    call SetMatrixValue(Me%LocalAtm%WindVelocityV%Field, Me%Size2D,         &
        !                        Me%ExtAtm%WindVelocityV%Field)
        !
        !    call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindVelocityV%Field, STAT = STAT_CALL)
        !    if (STAT_CALL .ne. SUCCESS_)                                            &
        !        stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04'
        !
        !endif


        if(Me%LocalAtm%WindDirection%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%WindDirection%Field,&
                                       ID           = WindDirection_,               &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04a'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%WindDirection%Field, Me%Size2D,             &
                                Me%ExtAtm%WindDirection%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%WindDirection%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR04b'

        endif


        !Copy to local variable Atmospheric Pressure
        if(Me%LocalAtm%AtmosphericPressure%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%AtmosphericPressure%Field,&
                                       ID           = AtmosphericPressure_,         &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR05'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%AtmosphericPressure%Field, Me%Size2D,   &
                                Me%ExtAtm%AtmosphericPressure%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AtmosphericPressure%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR06'

        endif

        !Copy to local variable Mean Sea Level Pressure
        if(Me%LocalAtm%Mslp%Yes)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,             &
                                       Scalar       = Me%ExtAtm%Mslp%Field,&
                                       ID           = MeanSeaLevelPressure_,         &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR07'

            !Memory duplication
            call SetMatrixValue(Me%LocalAtm%Mslp%Field, Me%Size2D,   &
                                Me%ExtAtm%Mslp%Field)

            call UngetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%Mslp%Field, STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                            &
                stop 'ModifyLocalAtmVariables - ModuleInterfaceWaterAir - ERR08'

        endif

    end subroutine ModifyLocalAtmVariables

    !--------------------------------------------------------------------------

    subroutine ModifyWaterAirFluxes

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-------------------------------------------------------------
        type(T_Property), pointer                   :: PropertyX, PropertyY, PropertyZ
        type(T_Property), pointer                   :: PropertyW, PropertyA, PropertyB, PropertyC

        !Begin-------------------------------------------------------------

        if(Me%IntOptions%Albedo)then
            call Search_Property(PropertyY, Albedo_,             STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyAlbedo            (PropertyY)
        endif
        
        if(Me%IntOptions%LatentHeat .and. (.not. Me%COARE))then
            call Search_Property(PropertyX, LatentHeat_,             STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyLatentHeat            (PropertyX)
        endif

        if(Me%IntOptions%SensibleHeat .and. (.not. Me%COARE))then
            call Search_Property(PropertyX, SensibleHeat_,             STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifySensibleHeat            (PropertyX)
        endif

        if(Me%IntOptions%Evaporation .and. (.not. Me%COARE))then
            call Search_Property(PropertyX, Evaporation_,            STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyEvaporation           (PropertyX)
        endif

        if(Me%IntOptions%NetLongWaveRadiation .and. (.not. Me%COARE))then
            call Search_Property(PropertyX, NetLongWaveRadiation_,      STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyNetLongWaveRadiation     (PropertyX)
        endif

        if(Me%IntOptions%NonSolarFlux .and. (.not. Me%COARE))then
            call Search_Property(PropertyX, NonSolarFlux_,           STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyNonSolarFlux          (PropertyX)
        endif

        if(Me%IntOptions%OxygenFlux)then
            call Search_Property(PropertyX, OxygenFlux_,             STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyOxygenFlux            (PropertyX)
        endif
           
        if(Me%IntOptions%CarbonDioxideFlux)then
            call Search_Property(PropertyX, CarbonDioxideFlux_,      STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyCarbonDioxideFlux  (PropertyX)
        endif
        
        if(Me%IntOptions%AmmoniaFlux)then
            call Search_Property(PropertyX, AmmoniaFlux_,      STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_)  call ModifyAmmoniaFlux   (PropertyX)
        endif
        
        if(Me%IntOptions%NitrateFlux)then
            call Search_Property(PropertyX, NitrateFlux_,      STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyNitrateFlux  (PropertyX)
        endif

        if(Me%IntOptions%WindShearVelocity)then
            call ModifyWindShearVelocity     
        endif

        if(Me%IntOptions%TurbulentKineticEnergy)then
            call Search_Property(PropertyX, TurbulentKineticEnergy_, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifyTurbulentKE           (PropertyX)
        endif


        if(Me%IntOptions%SurfaceRadiation)then
            call Search_Property(PropertyX, SurfaceRadiation_,       STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifySurfaceRadiation      (PropertyX)
        endif

        if(Me%IntOptions%SurfaceWaterFlux .and. (.not. Me%COARE))then
            call Search_Property(PropertyX, SurfaceWaterFlux_,       STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) call ModifySurfaceWaterFlux      (PropertyX)
        endif
        
        if(Me%IntOptions%WindStress)then
            !call Search_Property(PropertyX, WindStressX_,            STAT = STAT_CALL)
            !if (STAT_CALL == SUCCESS_) then
            !    call Search_Property(PropertyY, WindStressY_,        STAT = STAT_CALL)
            !    if (STAT_CALL == SUCCESS_)  call ModifyWindStress    (PropertyX, PropertyY)
            !endif
            
            call Search_Property(PropertyX, WindStress_,        STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_)  call ModifyWindStress    (PropertyX)            
        endif

        if (Me%COARE)then
        
        call Search_Property(PropertyX, LatentHeat_,                STAT = STAT_CALL)
        call Search_Property(PropertyY, SensibleHeat_,              STAT = STAT_CALL)
        call Search_Property(PropertyZ, NetLongWaveRadiation_,      STAT = STAT_CALL)
        call Search_Property(PropertyW, UpwardLongWaveRadiation_,   STAT = STAT_CALL)
        call Search_Property(PropertyA, DownwardLongWaveRadiation_, STAT = STAT_CALL)
        call Search_Property(PropertyB, SurfaceRadiation_,          STAT = STAT_CALL)
        call Search_Property(PropertyC, NonSolarFlux_,              STAT = STAT_CALL)
        !call Search_Property(PropertyD, Albedo_,                    STAT = STAT_CALL)
        call CheckRadiationOptions(PropertyZ, PropertyW, PropertyA, PropertyB)
        call CheckLatentSensibleOptions (PropertyX, PropertyY)
        
        call ComputeCOAREHeatBudget(PropertyX, PropertyY, PropertyZ, PropertyW, PropertyA, PropertyB, PropertyC)
            
            if(Me%IntOptions%Evaporation)then
                call Search_Property(PropertyX, Evaporation_,            STAT = STAT_CALL)
                if (STAT_CALL == SUCCESS_) call ModifyEvaporation           (PropertyX)
            endif
        
        endif

    end subroutine ModifyWaterAirFluxes

    !--------------------------------------------------------------------------
    subroutine ModifyAlbedo (PropAlbedo)
    
        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropAlbedo

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropAlbedo%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropAlbedo%ID%ObjFillMatrix,  &
                                  Matrix2D          = PropAlbedo%Field,             &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,        &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyAlbedo - ModuleInterfaceWaterAir - ERR01'
            
            call SetMatrixValue(Me%SurfaceAlbedo, Me%Size2D,   &
                                PropAlbedo%Field)
            
        endif    
        
    end subroutine ModifyAlbedo
    
    !----------------------------------------------------------------------------
 
    subroutine ComputeLatentHeat (PropLatentHeat)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat

        !Local-----------------------------------------------------------------
        real,    dimension(:,:  ), pointer          :: UWIND, VWIND
        real,    dimension(:,:,:), pointer          :: WaterTemperature
        real,    dimension(:,:  ), pointer          :: AirTemperature, RelativeHumidity
        real                                        :: WindVelocity
        real                                        :: SurfaceAirTemperature
        real                                        :: SurfaceRelHumidity
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j, STAT_CALL
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR01'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%RelativeHumidity%Field, &
                                   ID           = RelativeHumidity_,                &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR02'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR03'

            
        !Shorten Variables
        UWIND            => Me%LocalAtm%WindVelocityU%Field
        VWIND            => Me%LocalAtm%WindVelocityV%Field
        AirTemperature   => Me%ExtAtm%AirTemperature%Field
        WaterTemperature => Me%ExtWater%WaterTemperature
        RelativeHumidity => Me%ExtAtm%RelativeHumidity%Field

        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Calculates wind^2
                WindVelocity = sqrt(UWIND(i, j)**2. + VWIND(i, j)**2.)

                !We just consider the loss of heat (evaporation), once when water condenses
                !the heat gain is absorved by the air
                !Frank/Ramiro 26-10-2000

                !When the Water temperature > AirTemperature, the flux is negativ
!                PropLatentHeat%Field(i, j) = -1. * (19.0 + 0.95 * WindVelocity2)                  * &
!                                              (SaturatedVaporPressure(WaterTemperature(i, j,KUB)) - &
!                                              RelativeHumidity(i, j)                             * &
!                                              (SaturatedVaporPressure(AirTemperature  (i, j    ))))
                                              
                
!                if (PropLatentHeat%Field(i, j) > 0.) PropLatentHeat%Field(i, j) = 0.

!               When water temperature was higher than that of Air temperature, the latent heat losses were too high.
!               The new code considers that the air temperature near water in this cases is closer to the water temperature
!               value (and never higher than the water temperature) and is dependent on the wind velocity.
!               Variable SurfaceAirTemperature added. The same happens with the relative humidity.
                if ((PropLatentHeat%SVPMethod == 2) .and. (WaterTemperature(i, j, KUB) > AirTemperature(i, j))) then
                !         ºC                       ºC + ºC * (ms-1 / (ms-1 + ms-1)** ()
                    SurfaceAirTemperature = AirTemperature(i, j) + (WaterTemperature(i, j, KUB)- AirTemperature(i, j))*  &
                                            (PropLatentHeat%C1/(PropLatentHeat%C1+WindVelocity))**PropLatentHeat%C2
                    
                    if(RelativeHumidity(i,j) > 0.7) then
                !                              % + % * (ms-1 / (ms-1 + ms-1)** ()
                        SurfaceRelHumidity = RelativeHumidity(i,j) + (1 - RelativeHumidity(i,j))*   &
                                            (PropLatentHeat%C1/(PropLatentHeat%C1+WindVelocity))**PropLatentHeat%C2
                    else
                        SurfaceRelHumidity = RelativeHumidity(i,j)
                    endif                        
                         
                else
                  
                    SurfaceAirTemperature = AirTemperature(i, j)                                                        
                    SurfaceRelHumidity = RelativeHumidity(i,j)
                
                    
                endif
                                             
                PropLatentHeat%Field(i, j) = LatentHeat    (ReferenceDensity, WaterTemperature(i, j, KUB),  &
                                                            SurfaceAirTemperature, SurfaceRelHumidity, &
                                                            WindVelocity)
                   

            endif

        enddo
        enddo

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR04'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%RelativeHumidity%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR05'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeLatentHeat - ModuleInterfaceWaterAir - ERR06'


        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(AirTemperature)
        nullify(WaterTemperature)
        nullify(RelativeHumidity)
    
    end subroutine ComputeLatentHeat


    subroutine ModifyLatentHeat(PropLatentHeat)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropLatentHeat%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropLatentHeat%ID%ObjFillMatrix,  &
                                  Matrix2D          = PropLatentHeat%Field,             &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,        &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyLatentHeat - ModuleInterfaceWaterAir - ERR01'

        elseif(.not. PropLatentHeat%Constant)then
            
            call ComputeLatentHeat(PropLatentHeat)

        endif


    end subroutine ModifyLatentHeat

    
    !------------------------------------------------------------------------------

    
    subroutine ModifyNonSolarFlux(PropNonSolarFlux)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropNonSolarFLux

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropNonSolarFlux%ID%SolutionFromFile) then
            
            call ModifyFillMatrix(FillMatrixID      = PropNonSolarFlux%ID%ObjFillMatrix,    &
                                  Matrix2D          = PropNonSolarFlux%Field,               &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNonSolarFlux - ModuleInterfaceWaterAir - ERR01'

        elseif (.not. PropNonSolarFlux%Constant) then

            call ComputeNonSolarFlux(PropNonSolarFlux)

        endif


    end subroutine ModifyNonSolarFlux
    
    
    !--------------------------------------------------------------------------

    subroutine ComputeNonSolarFlux(PropNonSolarFlux)

        !Arguments-------------------------------------------------------------
        type(T_Property),       pointer         :: PropNonSolarFlux
        
        !Local-----------------------------------------------------------------
        type(T_Property),       pointer         :: PropertyX
        real, dimension(:,:),   pointer         :: LatentHeat
        real, dimension(:,:),   pointer         :: SensibleHeat
        real, dimension(:,:),   pointer         :: NetLongWaveRadiation
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call Search_Property(PropertyX, LatentHeat_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNonSolarFlux - ModuleInterfaceWaterAir - ERR01'
        LatentHeat => PropertyX%Field

        call Search_Property(PropertyX, SensibleHeat_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNonSolarFlux - ModuleInterfaceWaterAir - ERR02'
        SensibleHeat => PropertyX%Field


        call Search_Property(PropertyX, NetLongWaveRadiation_, .true., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeNonSolarFlux - ModuleInterfaceWaterAir - ERR03'
        NetLongWaveRadiation => PropertyX%Field

        PropNonSolarFlux%Field = LatentHeat + SensibleHeat + NetLongWaveRadiation

        nullify(LatentHeat, SensibleHeat, NetLongWaveRadiation)

    end subroutine ComputeNonSolarFlux

    !--------------------------------------------------------------------------    

    
    subroutine ModifySurfaceRadiation (PropSurfaceRadiation)
        
        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropSurfaceRadiation

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropSurfaceRadiation%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropSurfaceRadiation%ID%ObjFillMatrix,  &
                                  Matrix2D          = PropSurfaceRadiation%Field,             &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,        &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifySurfaceRadiation - ModuleInterfaceWaterAir - ERR01a'


        elseif(.not. PropSurfaceRadiation%Constant)then

            call ComputeSurfaceRadiation(PropSurfaceRadiation)

        endif

    end subroutine ModifySurfaceRadiation


    !--------------------------------------------------------------------------    

    
    subroutine ComputeSurfaceRadiation (PropSurfaceRadiation)
        
        !Arguments--------------------------------------------------------------
        type(T_Property), pointer           :: PropSurfaceRadiation
        type(T_Property), pointer           :: PropertyX
        !External--------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, i, j
        real, dimension(:,:), pointer       :: SolarRadiation
        integer                             :: STAT_CALL
        integer                             :: CHUNK

        !Begin-----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB


        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                   Scalar       = SolarRadiation,           &
                                   ID           = SolarRadiation_,          &
                                   STAT         = STAT_CALL)
        if (STAT_CALL .ne. SUCCESS_)                                        &
            stop 'ComputeSurfaceRadiation - ModuleInterfaceWaterAir - ERR02'
        
        call Search_Property(PropertyX, Albedo_,                 STAT = STAT_CALL)
                
        if (STAT_CALL == SUCCESS_) then
                
            call SetMatrixValue(Me%SurfaceAlbedo, Me%Size2D,   &
                                    PropertyX%Field)
        else
            if (Me%ReflectionCoef < 0 .OR. Me%ReflectionCoef > 1) then
                    write(*,*) 'Albedo must be between 0 and 1'
                    stop 'Construct_PropertyValues - ModuleInterfaceWaterAir - ERR40'
            endif
            
            call SetMatrixValue(Me%SurfaceAlbedo, Me%Size2D,   &
                    Me%ReflectionCoef)
        endif
                
                CHUNK = CHUNK_J(JLB, JUB)
                !$OMP PARALLEL PRIVATE(i,j)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j=JLB, JUB
                do i=ILB, IUB
            
                    if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                        PropSurfaceRadiation%Field(i, j) = SolarRadiation(i, j) * (1 - Me%SurfaceAlbedo(i, j))
                                                                
                    endif

                enddo
                enddo
                !$OMP END DO
                !$OMP END PARALLEL

            !if (MonitorPerformance) then
            !    call StopWatch ("ModuleInterfaceWaterAir", "ComputeSurfaceRadiation")
            !endif

        call UnGetAtmosphere(Me%ObjAtmosphere, SolarRadiation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSurfaceRadiation - ModuleInterfaceWaterAir - ERR04'


    end subroutine ComputeSurfaceRadiation

    !--------------------------------------------------------------------------
        
    subroutine ComputeSensibleHeat(PropSensibleHeat)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropSensibleHeat

        !Local-----------------------------------------------------------------
        real,    dimension(:,:  ), pointer          :: UWIND, VWIND
        real,    dimension(:,:,:), pointer          :: WaterTemp
        real,    dimension(:,:  ), pointer          :: AirTemp
        real                                        :: WindVelocity
        integer                                     :: ILB, IUB, JLB, JUB, KUB
        integer                                     :: i, j, STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR01'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR02'


        UWIND            => Me%LocalAtm%WindVelocityU%Field
        VWIND            => Me%LocalAtm%WindVelocityV%Field
        AirTemp          => Me%ExtAtm%AirTemperature%Field
        WaterTemp        => Me%ExtWater%WaterTemperature

        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Calculates the Modulus of the wind
                WindVelocity = sqrt(UWIND(i, j)**2. + VWIND(i, j)**2.)

                PropSensibleHeat%Field(i, j) = SensibleHeat (ReferenceDensity, WaterTemp(i, j, KUB), &
                                                             AirTemp(i, j), WindVelocity)

                !When the Water temperature > AirTemperature, the flux is negativ
!                PropSensibleHeat%Field(i, j) = -1. * BowenCoefficient               *   &
!                                                     (19.0 + 0.95 * WindVelocity2)  *   &
!                                                     (WaterTemp(i, j, KUB) - AirTemp(i, j))
            endif

        enddo
        enddo


        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR03'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeSensibleHeat - ModuleInterfaceWaterAir - ERR04'


        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(AirTemp)
        nullify(WaterTemp)

    end subroutine ComputeSensibleHeat

    
    !--------------------------------------------------------------------------
    
    subroutine ModifySensibleHeat(PropSensibleHeat)

        !Arguments-------------------------------------------------------------

        type(T_Property), pointer :: PropSensibleHeat

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL


        !Begin-----------------------------------------------------------------

        if (PropSensibleHeat%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropSensibleHeat%ID%ObjFillMatrix,    &
                                  Matrix2D          = PropSensibleHeat%Field,               &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifySensibleHeat - ModuleInterfaceWaterAir - ERR01'

        elseif (.not. PropSensibleHeat%Constant) then

            call ComputeSensibleHeat(PropSensibleHeat)

        endif

    end subroutine ModifySensibleHeat
    
    !------------------------------------------------------------------------------

    subroutine ComputeEvaporation(PropEvaporation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropEvaporation

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        integer                                     :: STAT_CALL
        integer                                     :: CHUNK
         
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

        !Looks if the property also calculates latent heat. If so, use the evaporation
        !corresponding to the latent heat, if not check for other sources.
        call Search_Property (PropLatentHeat, LatentHeat_, STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then

            if (MonitorPerformance) then
                call StartWatch ("ModuleInterfaceWaterAir", "ComputeEvaporation")
            endif

            CHUNK = CHUNK_J(JLB, JUB)
            !$OMP PARALLEL PRIVATE(i,j)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                ![m3/s]                    = [J/m2/s] / [J/kg] / [kg/m3] * [m] * [m]
                PropEvaporation%Field(i, j) = PropLatentHeat%Field(i, j)         / &
                                              RefLatentHeatOfVaporization           / &
                                              ReferenceDensity                   * &
                                              Me%ExternalVar%GridCellArea(i, j)
            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
        
            if (MonitorPerformance) then
                call StopWatch ("ModuleInterfaceWaterAir", "ComputeEvaporation")
            endif
            
        endif


    end subroutine ComputeEvaporation

    !--------------------------------------------------------------------------
    
    subroutine ModifyEvaporation(PropEvaporation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropEvaporation

        !External--------------------------------------------------------------
        integer                                     :: i,j, STAT_CALL
        logical                                     :: Accumulate, Interpolate, UseOriginal
        integer                                     :: CHUNK
        !Begin-----------------------------------------------------------------

        if (PropEvaporation%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropEvaporation%ID%ObjFillMatrix,     &
                                  Matrix2D          = PropEvaporation%Field,                &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR01'

        elseif(.not. PropEvaporation%Constant)then

            call ComputeEvaporation(PropEvaporation)

        endif
        
        !if not computed, let the user input it with mm (that is what monitoring stations have).
        !using m3/s user needs to have knowledge of surface area that in reservoirs it changes a lot
        !and also need to know model cells. it does not make sense to only allow input as m3/s
        if (PropEvaporation%ID%SolutionFromFile .or. PropEvaporation%Constant) then
            
            call GetValuesProcessingOptions (PropEvaporation%ID%ObjFillMatrix,    &
                                             Accumulate = Accumulate,               &
                                             Interpolate = Interpolate,             &
                                             UseOriginal = UseOriginal,             &
                                             STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR010'             
            
            
            if (PropEvaporation%FirstActualization) then
        
                if (Interpolate) then
            
                    write(*,*) "Evaporation can't be INTERPOLATED. Check evaporation property options."
                    stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR050'
                
                elseif (Accumulate) then
            
                    if (trim(PropEvaporation%ID%Units) /= 'mm') then
                
                        write(*,*)'Invalid Evaporation Units for accumulated values'
                        write(*,*)'Use mm with ACCUMULATE_VALUES = 1 or '
                        write(*,*)'use a FLUX (ex. mm/hour) with USE_ORIGINAL_VALUES = 1'
                        stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR060'
                    
                    endif
                
                    if (PropEvaporation%Constant) then
                
                        write(*,*)'Invalid Evaporation Units for constant value'
                        write(*,*)'Use a FLUX (ex. mm/hour) when using constant value for evaporation'
                        stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR070'                
                
                    endif
                
                else !use original
                
                    if (trim(PropEvaporation%ID%Units) == 'mm') then
                
                        write(*,*)'Invalid Evaporation Units for original values'
                        write(*,*)'Use mm with ACCUMULATE_VALUES = 1 or '
                        write(*,*)'use a FLUX (ex. mm/hour) with USE_ORIGINAL_VALUES = 1'
                        stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR080'
                    
                    endif                
                
                endif

                if (trim(adjustl(PropEvaporation%ID%Units)) /= 'm3/s') then
            
                    Me%EvapReqConv = .true.

                    select case (PropEvaporation%ID%Units)

                    case ('mm/day')

                        Me%ConversionFactorEvap = 1. / 86400000.        !In m/s

                    case ('mm/hour')

                        Me%ConversionFactorEvap = 1. / 3600000.         !In m/s
                
                    case ('mm')
                
                        if(.not. Accumulate) then
                            write(*,*)'when using "mm" as units for evaporation'
                            write(*,*)'you must use ACCUMULATE_VALUES = 1'
                            stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR090'
                        endif
                    
                        Me%ConversionFactorEvap = 1. / 1000.            !In m/s => Fillmatrix converted from mm to mm/s

                    case default

                        write(*,*)'Invalid Evaporation Units'
                        write(*,*)'Use mm, m3/s, mm/day or mm/hour'
                        stop 'ModifyEvaporation - ModuleInterfaceWaterAir - ERR100'

                    end select 
                
                else
            
                    Me%EvapReqConv = .false.
                
                endif       
            endif

            if (PropEvaporation%FirstActualization .or. .not. PropEvaporation%Constant) then

                if (Me%EvapReqConv) then

                    CHUNK = CHUNK_J(Me%WorkSize2D%JLB, Me%WorkSize2D%JUB)
                    !$OMP PARALLEL PRIVATE(i,j)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                    do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB

                        if (Me%ExtWater%WaterPoints2D(i, j) == 1) then
                            PropEvaporation%Field(i, j) = PropEvaporation%Field(i, j)      * &
                                                            (Me%ExternalVar%GridCellArea(i, j) * &
                                                            Me%ConversionFactorEvap)  !In m3/s                            
                        endif
                    enddo
                    enddo
                    !$OMP END DO
                    !$OMP END PARALLEL             

                endif                  

            endif            
            
        endif
    
    end subroutine ModifyEvaporation

    
    !--------------------------------------------------------------------------

    
    subroutine ModifySurfaceWaterFlux(PropSurfaceWaterFlux)

        !Arguments-------------------------------------------------------------
        type(T_Property),       pointer             :: PropSurfaceWaterFlux

        !Local-----------------------------------------------------------------
        type(T_Property),       pointer             :: PropEvaporation
        real, dimension(:,:),   pointer             :: Precipitation
        real, dimension(:,:),   pointer             :: Irrigation
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call SetMatrixValue(PropSurfaceWaterFlux%Field, Me%Size2D, 0.)

        if(Me%IntOptions%Evaporation)then
            !Gets Evaporation
            call Search_Property(PropEvaporation, Evaporation_, .false., STAT = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_ .and. STAT_CALL .ne. NOT_FOUND_ERR_)then
                stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR01'
            elseif(STAT_CALL .eq. SUCCESS_)then
                PropSurfaceWaterFlux%Field = PropSurfaceWaterFlux%Field + PropEvaporation%Field
            endif
        endif


        if(Me%ExtOptions%Precipitation)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                       Scalar       = Precipitation,            &
                                       ID           = Precipitation_,           &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                        &
                stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR02'
        
            PropSurfaceWaterFlux%Field = PropSurfaceWaterFlux%Field + Precipitation

            call UnGetAtmosphere(Me%ObjAtmosphere, Precipitation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR04'


        endif

        if(Me%ExtOptions%Irrigation)then
            call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                       Scalar       = Irrigation,               &
                                       ID           = Irrigation_,              &
                                       STAT         = STAT_CALL)
            if (STAT_CALL .ne. SUCCESS_)                                        &
                stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR03'
            
        
            PropSurfaceWaterFlux%Field = PropSurfaceWaterFlux%Field + Irrigation
            
            call UnGetAtmosphere(Me%ObjAtmosphere, Irrigation, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifySurfaceWaterFlux - ModuleInterfaceWaterAir - ERR05'

        endif

    end subroutine ModifySurfaceWaterFlux
    
    !--------------------------------------------------------------------------

    subroutine ComputeUpLongWaveRad(PropUpwardLongWaveRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropUpwardLongWaveRadiation

        !Local-----------------------------------------------------------------
        real, dimension   (:,:,:), pointer          :: WaterTemp
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ComputeUpLongWaveRad - ModuleInterfaceWaterAir - ERR10'

       
        WaterTemp        => Me%ExtWater%WaterTemperature


        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                PropUpwardLongWaveRadiation%Field(i, j) = LongWaveUpward (WaterTemp(i, j, KUB))


            endif
        enddo
        enddo

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeUpLongWaveRad - ModuleInterfaceWaterAir - ERR20'

         nullify(WaterTemp)

    end subroutine ComputeUpLongWaveRad

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------


    subroutine ComputeDownLongWaveRad(PropDownwardLongWaveRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropDownwardLongWaveRadiation

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB

       
        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%RelativeHumidity%Field, &
                                   ID           = RelativeHumidity_,                &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR10'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%CloudCover%Field,       &
                                   ID           = CloudCover_,                      &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR20'

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR30'


        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                PropDownwardLongWaveRadiation%Field(i, j) = LongWaveDownward (Me%ExtAtm%CloudCover%Field    (i, j),      &
                                                                              Me%ExtAtm%AirTemperature%Field(i, j)) 

            endif
        enddo
        enddo

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%RelativeHumidity%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR40'
        
        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%CloudCover%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR50'

        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDownLongWaveRad - ModuleInterfaceWaterAir - ERR60'


    end subroutine ComputeDownLongWaveRad

    !--------------------------------------------------------------------------
    
    subroutine ModifyNetLongWaveRadiation(PropNetLongWaveRadiation)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropNetLongWaveRadiation
        
        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: PropUpLongWaveRadiation
        type(T_Property), pointer                   :: PropDownLongWaveRadiation
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

i1:     if (PropNetLongWaveRadiation%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropNetLongWaveRadiation%ID%ObjFillMatrix,   &
                                  Matrix2D          = PropNetLongWaveRadiation%Field,              &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,                   &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNetLongWaveRadiation - ModuleInterfaceWaterAir - ERR10'

        elseif (.not. PropNetLongWaveRadiation%Constant)then i1
            
            call Search_Property(PropUpLongWaveRadiation,       UpwardLongWaveRadiation_,      STAT = STAT_CALL)
            call Search_Property(PropDownLongWaveRadiation,     DownwardLongWaveRadiation_,    STAT = STAT_CALL)
            
i2:         if (PropUpLongWaveRadiation%ID%SolutionFromFile) then

                call ModifyFillMatrix(FillMatrixID      = PropUpLongWaveRadiation%ID%ObjFillMatrix,&
                                      Matrix2D          = PropUpLongWaveRadiation%Field,           &
                                      PointsToFill2D    = Me%ExtWater%WaterPoints2D,               &
                                      STAT              = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNetLongWaveRadiation - ModuleInterfaceWaterAir - ERR20'

            elseif (.not. PropUpLongWaveRadiation%Constant) then i2
                    
                call ComputeUpLongWaveRad(PropUpLongWaveRadiation)

            endif i2


i3:         if (PropDownLongWaveRadiation%ID%SolutionFromFile) then

                call ModifyFillMatrix(FillMatrixID      = PropDownLongWaveRadiation%ID%ObjFillMatrix,&
                                      Matrix2D          = PropDownLongWaveRadiation%Field,           &
                                      PointsToFill2D    = Me%ExtWater%WaterPoints2D,                     &
                                      STAT              = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'ModifyNetLongWaveRadiation - ModuleInterfaceWaterAir - ERR30'

            elseif (.not. PropDownLongWaveRadiation%Constant) then i3

                call ComputeDownLongWaveRad(PropDownLongWaveRadiation)

            endif i3

            PropNetLongWaveRadiation%Field(:,:) = PropUpLongWaveRadiation%Field(:,:) + &
                                                  PropDownLongWaveRadiation%Field(:,:)

        endif i1

    end subroutine ModifyNetLongWaveRadiation
    
    !--------------------------------------------------------------------------

    subroutine ComputeCOAREHeatBudget(PropLatentHeat, PropSensibleHeat, PropNetLongWaveRadiation, PropUpLongWaveRadiation, &
                                      PropDownLongWaveRadiation, PropSurfaceRadiation, PropNonSolarFlux)
    
        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat
        type(T_Property), pointer                   :: PropSensibleHeat
        type(T_Property), pointer                   :: PropNetLongWaveRadiation
        type(T_Property), pointer                   :: PropUpLongWaveRadiation
        type(T_Property), pointer                   :: PropDownLongWaveRadiation
        type(T_Property), pointer                   :: PropSurfaceRadiation
        type(T_Property), pointer                   :: PropNonSolarFlux
        !type(T_Property), pointer                   :: PropAlbedo
    
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
    
        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer          :: WaterTemperature, WaterDensity
        !real,    dimension(:,:  ), pointer          :: SolarRadiation
        real                                        :: CheckNightTime, WaterTemperatureDepth
        real                                        :: Year, Month, Day, Hour, Minute, Second
        Integer, dimension(6    )                   :: AuxTimeInstant
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j, CHUNK
        !Begin-----------------------------------------------------------------
    
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        KUB = Me%WorkSize3D%KUB
    
        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR01'
        
        call GetGeometryDistances (Me%ObjGeometry, DWZ = Me%ExternalVar%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR02'

        !call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
        !                           Scalar       = SolarRadiation,                   &
        !                           ID           = SolarRadiation_,                  &
        !                           STAT         = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR03'
        
        call GetDensity(Me%ObjWaterProperties, Me%ExtWater%Density,                  &
                        Me%ExternalVar%Now,  STAT = STAT_CALL)       
        if (STAT_CALL /= SUCCESS_)stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR04'
        
        
        Me%WarmingAboveWaterTemp = 0.0
        
        call SetDate(Me%LastTimeInstant, 0, 0, 0, 0, 0, 0)
        
!................   Shorten Variables   ..............................................
        WaterDensity        => Me%ExtWater%Density
        WaterTemperature    => Me%ExtWater%WaterTemperature       
!................   fill surfacetemperature with watertemperature values at KUB, and compute surface radiation   ...........
        
        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                Me%SurfaceTemperature (i, j)     = WaterTemperature(i, j, KUB)    
                !PropSurfaceRadiation%Field(i, j) = SolarRadiation(i, j) * (1 - Me%ReflectionCoef)
                !PropSurfaceRadiation%Field(i, j) = SolarRadiation(i, j) * (1 - PropAlbedo%Field(i, j))
            endif
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
!......................   warm layer..................................................
        
        if (Me%Jwarm) Then
        
            if (.not. Me%Jump)then
            
                call ExtractDate(Me%ExternalVar%Now,                         &
                                    Year = Year, Month  = Month,  Day    = Day, &
                                    Hour = Hour, Minute = Minute, Second = Second)

                AuxTimeInstant(1) = int(Year  )
                AuxTimeInstant(2) = int(Month )
                AuxTimeInstant(3) = int(Day   )
                AuxTimeInstant(4) = int(Hour  )
                AuxTimeInstant(5) = int(Minute)
                AuxTimeInstant(6) = int(Second)            
         
                CheckNightTime = AuxTimeInstant(4)*3600.0 + AuxTimeInstant(5)*60.0 + AuxTimeInstant(6)
            
                if ((CheckNightTime .ge. 0.0) .and. (CheckNightTime .le. 21600.0)) Then  !checks if time is between 0-6h
                                                                                         !for warm layer existence and checks if first time through.
                
                        if(Me%ExternalVar%Now .lt. Me%LastTimeInstant) Then !restarts some constants and does not compute warm layer
                                
                            !restart variables
                            call SetMatrixValue(Me%Fxp, Me%Size2D, 0.5)
                            call SetMatrixValue(Me%WarmLayerThickness, Me%Size2D, 19.)
                            call SetMatrixValue(Me%Tau_ac, Me%Size2D, 0.)
                            call SetMatrixValue(Me%AccumulatedEnergy, Me%Size2D, 0.)
                            call SetMatrixValue(Me%WarmLayerTempDiff, Me%Size2D, 0.)
                        
                        else
                
                            call ComputeWarmLayer (PropNetLongWaveRadiation, PropUpLongWaveRadiation,           &
                                                   PropDownLongWaveRadiation, PropSurfaceRadiation, WaterDensity)!, WaterDensity
                    
                        endif

                endif
            
            endif
            
            Me%Jump = .false.
            
        endif  !end warm layer calculation
        
!...................do cycle for surface temperature calculation.........................
        
        !$OMP PARALLEL PRIVATE(i,j,WaterTemperatureDepth)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                
                WaterTemperatureDepth = Me%ExternalVar%DWZ(i, j, KUB) / 2.0
                
                if (Me%WarmLayerThickness(i, j) .LT. WaterTemperatureDepth) then                    
                    Me%SurfaceTemperature(i,j) = Me%SurfaceTemperature(i,j) + Me%WarmLayerTempDiff(i, j) 
                else
                    Me%SurfaceTemperature(i,j) = Me%SurfaceTemperature(i,j) + Me%WarmLayerTempDiff(i, j) * &
                                                 WaterTemperatureDepth / Me%WarmLayerThickness(i, j) 
                endif
    
            endif
            
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
              
!................. Surface fluxes computation .................................        
        ! computes the COARE3.0-based algorithm
        call ComputeCOAREsurfacefluxes(PropLatentHeat, PropSensibleHeat, PropNetLongWaveRadiation, PropUpLongWaveRadiation, &
                                       PropDownLongWaveRadiation, PropSurfaceRadiation, PropNonSolarFlux, WaterDensity)  
    
!................. update last instant             
        Me%LastTimeInstant = Me%ExternalVar%Now
        
!************************************* UNGET *****************************************************
        
        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR05'

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR06'
        
        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%Density, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR07'

        !call UnGetAtmosphere(Me%ObjAtmosphere, SolarRadiation, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREHeatBudget - ModuleInterfaceWaterAir - ERR08'
        
        nullify (WaterTemperature)
        nullify (WaterDensity)
    
    end subroutine ComputeCOAREHeatBudget
    
    subroutine ComputeWarmLayer (PropNetLongWaveRadiation, PropUpLongWaveRadiation,           &
                                 PropDownLongWaveRadiation, PropSurfaceRadiation, WaterDensity)
    
        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropNetLongWaveRadiation
        type(T_Property), pointer                   :: PropUpLongWaveRadiation
        type(T_Property), pointer                   :: PropDownLongWaveRadiation
        type(T_Property), pointer                   :: PropSurfaceRadiation
    
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer              :: WaterDensity
        real                                            :: RichNumb, Ctd1, Ctd2, TotalSurfaceHeatLoss, TotalHeatAbsorbed 
        real                                            :: ModelDT, Cpw, TotalEnergyAbsorbed
        integer                                         :: ILB, IUB, JLB, JUB, KUB, i, j, n, CHUNK
        !Begin-----------------------------------------------------------------
           
        call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                 &
            stop 'Subroutine ComputeWarmLayer - ModuleInterfaceWaterAir - ERR02'
        
        !Size SurfaceTemperature
        ILB   = Me%WorkSize2D%ILB
        JLB   = Me%WorkSize2D%JLB
        IUB   = Me%WorkSize2D%IUB
        JUB   = Me%WorkSize2D%JUB
        KUB   = Me%WorkSize3D%KUB
        Cpw   = 4000.0 ! J/Kg.K
        RichNumb  = 0.65
        
        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j, Ctd1, Ctd2, TotalSurfaceHeatLoss, TotalHeatAbsorbed, TotalEnergyAbsorbed)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)            
        do j = JLB, JUB
        do i = ILB, IUB
        
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
               
                Me%Al(i, j) = 2.1e-5 * (Me%SurfaceTemperature(i,j) + 3.2)**0.79    !Thermal expansion coefficient (1/k)
              
!***************   set warm layer constants  **************
                   

!                                  (kJ/kg.K)   / ((1/k) *    (m/s^2)*    (kg/m^3)        
                Ctd1 = sqrt(2 * RichNumb * Cpw/ (Me%Al(i, j) * 9.8 * WaterDensity(i, j, KUB)))
                !Factor in warm layer thickness equation - IN PAPER: "cool-skin and warm-layer effects on sea surface temperature"
                    
!                                     (1/k) * m/s^2  /                (kg/m^3))               / (kJ/kg.K)
                Ctd2    = sqrt(2 * Me%Al(i, j) * 9.8 / (RichNumb * WaterDensity(i, j, KUB))) / (Cpw**1.5)
                !Factor in warm layer temperature difference equation  - IN PAPER: "cool-skin and warm-layer effects on sea surface temperature"
                    
                PropUpLongWaveRadiation%Field(i, j)  = LongWaveUpwardCOARE(Me%SurfaceTemperature(i, j), Me%Jcool)
                    
                PropNetLongWaveRadiation%Field(i, j) = PropDownLongWaveRadiation%Field(i, j) + PropUpLongWaveRadiation%Field(i, j)
                                 
                TotalSurfaceHeatLoss  = -PropNetLongWaveRadiation%Field(i, j) - Me%LastSensibleHeat(i, j) - &
                                        Me%LastLatentHeat(i, j) - Me%RainFlux(i, j) 
                
                !total cooling at surface (w/m^2) RainFlux, latent and 
                !sensible heat are the ones calculated in the previous timestep

                TotalHeatAbsorbed  = Me%Fxp(i, j) * PropSurfaceRadiation%Field(i, j) - TotalSurfaceHeatLoss
                !total heat abs in warm layer (w/m2)                    
                    
                if (TotalHeatAbsorbed .LT. 50 .and. .not. Me%EnergyThreshold) then !check for threshold       
                    
                else
                    Me%EnergyThreshold = .true.    !indicates threshold crossed  ! dúvida
!                                                               Previous windstress                                                                     
                    Me%Tau_ac(i, j) = Me%Tau_ac(i, j) + max(0.002, Me%Tau(i, j)) * ModelDT   !momentum integral
                    
                    !          J/m2           +         (W/m2)    *    s
                    if ((Me%AccumulatedEnergy(i, j) + TotalHeatAbsorbed * ModelDT) .GT. 0) then  
                        !check threshold for warm layer existence.
                        !AccumulatedEnergy in the amount of energy stored, integrated over time
                        
!..........................................  loop 5 times for fxp  ...............................................
                        
                        do n = 1,5
                        
                        Me%Fxp(i, j) = 1.0 - (0.28 * 0.014 * (1 - exp(-Me%WarmLayerThickness(i, j) / 0.014)) + 0.27 * 0.357 * &
                                       (1 - exp(-Me%WarmLayerThickness(i, j) / 0.357)) + 0.45 * 12.82 *                       &
                                       (1 - exp(-Me%WarmLayerThickness(i, j) / 12.82))) / Me%WarmLayerThickness(i, j)
                        
!                       Where fxp is the average fraction of the solar radiation absorbed in the trapping layer
                        
                       ! J/m2               = (                     W/m2                      ) * s
                        TotalEnergyAbsorbed = (Me%Fxp(i, j) * PropSurfaceRadiation%Field(i, j) - TotalSurfaceHeatLoss) * ModelDT
                        
                                    !J/m2                      + J/m2
                                if (Me%AccumulatedEnergy(i, j) + TotalEnergyAbsorbed .GT. 0) then
                                    Me%WarmLayerThickness(i, j) = min(19.0, Ctd1 * Me%Tau_ac(i, j) /  &
                                                                  sqrt(Me%AccumulatedEnergy(i, j) + TotalEnergyAbsorbed)) 
                                endif
                                
                        enddo
!..........................................  end loop for fxp  ...............................................
                        
                    else    !warm layer wiped out
                        
                            Me%Fxp(i, j)                     = 0.75 
                            Me%WarmLayerThickness(i, j)      = 19.0 
                            TotalEnergyAbsorbed              = (Me%Fxp(i, j) *                                            &
                                                               PropSurfaceRadiation%Field(i, j) - TotalSurfaceHeatLoss) * &
                                                               ModelDT 
                            
                    endif !   end sign check on AccumulatedEnergy
                  
                            Me%AccumulatedEnergy(i, j) = Me%AccumulatedEnergy(i, j) + TotalEnergyAbsorbed  !heat integral
                            
!...........................       compute dt_warm     ..........................................................
                            
                    if (Me%AccumulatedEnergy(i, j) .GT. 0.0) then
                        Me%WarmLayerTempDiff(i, j) = Ctd2 * (Me%AccumulatedEnergy(i, j))**1.5 / Me%Tau_ac(i, j) 
                    else 
                        Me%WarmLayerTempDiff(i, j) = 0.0 
                    endif 
                        
                endif ! end threshold check
              
            endif ! end check for waterpoint
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
    end subroutine ComputeWarmLayer
    
    subroutine ComputeCOAREsurfacefluxes(PropLatentHeat, PropSensibleHeat, PropNetLongWaveRadiation, PropUpLongWaveRadiation, &
                                         PropDownLongWaveRadiation, PropSurfaceRadiation, PropNonSolarFlux, WaterDensity)

        !Arguments-------------------------------------------------------------
    
        type(T_Property), pointer                   :: PropLatentHeat
        type(T_Property), pointer                   :: PropSensibleHeat
        type(T_Property), pointer                   :: PropNetLongWaveRadiation
        type(T_Property), pointer                   :: PropUpLongWaveRadiation
        type(T_Property), pointer                   :: PropDownLongWaveRadiation
        type(T_Property), pointer                   :: PropNonSolarFlux

        type(T_Property), pointer                   :: PropSurfaceRadiation        

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer            :: AirTemperature, RelativeHumidity, Precipitation, PBLHeight
        real,    dimension(:,:,:), pointer          :: WaterDensity
        real,    dimension(:,:  ), pointer          :: UWIND, VWIND
        real                                        :: P, Rgas, Grav, Cpa, Cpw, Visw, Charn, CoolSkinThickness, COAREWindVelocity
        real                                        :: CoolSkinTempDepression, LatHeatVap, DeltaSatVapPressCorr, Bigc, AirDensity
        real                                        :: DeltaTempInterface, WetBulbFactor, WindVelocity
        real                                        :: Visa, ScailingTemperaturePar
        real                                        :: ScailingVelocityPar, CoolSkinHumDepression, DeltaSatVapPress, Tcw
        real                                        :: MoistureContentAir, InterfaceMoistureContent, ScailingMoisturePar
        integer                                     :: ILB, IUB, JLB, JUB, KUB, i, j, Nits
        !Begin-----------------------------------------------------------------

        !Important papers used in this algorithm:
        !Coolskin paper: (Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B. Edson, and G.S. Young, 1996a:
        !The cool skin and the warm layer in bulk flux calculations.)
        !Stability parameter paper : 
        !A. A. GRACHEV* AND C. W. FAIRALLDependence of the Monin–Obukhov Stability Parameter on the Bulk Richardson Number over the Ocean
        
        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                   ID           = AirTemperature_,                  &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR01'
        
        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%RelativeHumidity%Field, &
                                   ID           = RelativeHumidity_,                &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR02'

        call GetAtmosphereProperty(AtmosphereID  = Me%ObjAtmosphere,                &
                                    Scalar       = Precipitation,                   &
                                    ID           = Precipitation_,                  &
                                    STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR03'
        
        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = PBLHeight,                        &
                                   ID           = PBLHeight_,                       &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR04'
        
        
        Rgas    = 287.1    ! J/kg/K
        Grav    = 9.8      ! m/s2
        Cpa     = 1004.67     ! J/Kg.K
        Cpw     = 4000.0        ! J/Kg.K
        Visw    = 1e-6       ! kinematic viscosity of sea-water [m^2/s]
        Tcw     = 0.6         ! thermal conductivity of sea-water [W/m/K]
        P       = 1008.0     ! Atmospheric pressure
        Charn   = 0.011
        
        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB
        KUB = Me%WorkSize3D%KUB

        AirTemperature   => Me%ExtAtm%AirTemperature%Field
        RelativeHumidity => Me%ExtAtm%RelativeHumidity%Field
        UWIND            => Me%LocalAtm%WindVelocityU%Field
        VWIND            => Me%LocalAtm%WindVelocityV%Field
        
do1: do j = JLB, JUB
do2: do i = ILB, IUB
            
i1:     if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
            
            !Calculates wind^2
            WindVelocity                         = sqrt(UWIND(i, j)**2. + VWIND(i, j)**2.)            
            
            LatHeatVap                           = LatentHeatOfVaporization (Me%SurfaceTemperature(i, j))
            
            MoistureContentAir                   = COAREMoistureContentAir (AirTemperature(i, j), P, RelativeHumidity(i, j)) / 1000.
            
            InterfaceMoistureContent             = COAREInterfaceMoistureContent (Me%SurfaceTemperature(i, j), P) / 1000.

            PropUpLongWaveRadiation%Field(i, j)  = LongWaveUpwardCOARE (Me%SurfaceTemperature(i, j), Me%Jcool)
                    
            PropNetLongWaveRadiation%Field(i, j) = PropDownLongWaveRadiation%Field(i, j) + PropUpLongWaveRadiation%Field(i, j) 


        !*************  air constants ************
!            Kg/m3
             AirDensity = P * 100.0 / (Rgas * (AirTemperature(i, j) + 273.16) * (1.0 + 0.61 * MoistureContentAir))
!            m2/s
             Visa = AirViscosity(AirTemperature(i, j))
             
        !************  cool skin constants  *******
             
             Me%Al(i, j)   = 2.1e-5 * (Me%SurfaceTemperature(i,j) + 3.2)**0.79
             !BigConstant present in eq13 of the cool skin paper. Only used to reduce the size of the overall equation in the code.
             Bigc = 16.0 * Grav * Cpw * (WaterDensity(i, j, KUB) * Visw)**3.0 / (Tcw**2.0 * AirDensity**2.0)
             
  !          Correction factor for the specific humidity difference:           
             DeltaSatVapPressCorr = 0.622 * LatHeatVap * InterfaceMoistureContent /    &
             (Rgas * (Me%SurfaceTemperature(i,j) + 273.16)**2) 
     
        !***************   wave parameters  *********
        !     lwave = grav / 2.0 / pi * twave**2   for a future development
             
        !     cwave = grav / 2.0 / pi * twave     for a future development
             
           call FirstGuess (AirTemperature, MoistureContentAir, DeltaSatVapPressCorr, Visa, InterfaceMoistureContent,     &
                            WindVelocity, COAREWindVelocity, CoolSkinHumDepression, DeltaSatVapPress, DeltaTempInterface, &
                            CoolSkinTempDepression, ScailingVelocityPar, ScailingTemperaturePar, ScailingMoisturePar,     &
                            CoolSkinThickness, PBLHeight, Nits, i, j)
             
           if (COAREWindVelocity .GT. 10.0) then
               
              Charn = 0.011 + (COAREWindVelocity - 10.0) / (18.0 - 10.0) * (0.018 - 0.011) 
              
           endif 
           
           if (COAREWindVelocity .GT. 18.0) then
               
              Charn = 0.018
              
           endif 
           
           call BulkLoop(PropLatentHeat, PropSensibleHeat, PropNetLongWaveRadiation, PropSurfaceRadiation,  &
                         AirTemperature, WaterDensity, DeltaTempInterface, ScailingTemperaturePar, Charn,   &
                         ScailingMoisturePar, DeltaSatVapPress, CoolSkinTempDepression, AirDensity, Visa,   &
                         CoolSkinHumDepression, CoolSkinThickness, Bigc, LatHeatVap, ScailingVelocityPar,   &
                         COAREWindVelocity, DeltaSatVapPressCorr, MoistureContentAir, WindVelocity,         &
                         PBLHeight, i, j, KUB, Nits)
           

            Me%Tau(i, j)             = AirDensity * ScailingVelocityPar**2.0 * WindVelocity / COAREWindVelocity  
            !stress
            
            Me%LastSensibleHeat(i, j)= PropSensibleHeat%Field(i, j) 
            
            Me%LastLatentHeat(i, j)  = PropLatentHeat%Field(i, j)
            
            WetBulbFactor            = 1.0 /                                                    &
                                       (1.0 + (DeltaSatVapPressCorr * LatHeatVap *              &
                                       WaterVapourDiffusivity(AirTemperature(i, j))) /          &
                                       (Cpa * HeatDiffusivity(AirTemperature(i, j), AirDensity, Cpa)))
            
            !Precipitation is converted in the beginning to m3/s , 
            !so here it is converted back to mm/h
            Me%RainFlux(i, j)        = -1. * (Precipitation(i, j) * (1000 / Me%ExternalVar%GridCellArea(i, j)) *   &
                                       WetBulbFactor * Cpw * ((Me%SurfaceTemperature(i,j) - AirTemperature(i, j) - &
                                       CoolSkinTempDepression * Me%Jcool) + (InterfaceMoistureContent -            &
                                       MoistureContentAir - CoolSkinHumDepression * Me%Jcool) * LatHeatVap / Cpa))
            
            PropNonSolarFlux%Field(i, j)    = PropLatentHeat%Field(i, j) + PropSensibleHeat%Field(i, j) +          &
                                              PropNetLongWaveRadiation%Field(i, j) + Me%RainFlux(i, j)


        endif i1  !external water point
            
     enddo do2
     enddo do1
    
!************************************  UNGET   ******************************************************
        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%RelativeHumidity%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR05'
        
        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR06'
        
        call UnGetAtmosphere(Me%ObjAtmosphere, Precipitation, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR07'
        
        call UnGetAtmosphere(Me%ObjAtmosphere, PBLHeight, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Subroutine ComputeCOAREsurfacefluxes - ModuleInterfaceWaterAir - ERR08'
        
        nullify (AirTemperature)
        nullify (RelativeHumidity)
        nullify (UWIND)
        nullify (VWIND)
        nullify (PBLHeight)
        
    end subroutine ComputeCOAREsurfacefluxes

    subroutine FirstGuess(AirTemperature, MoistureContentAir, DeltaSatVapPressCorr, Visa, InterfaceMoistureContent,     &
                          WindVelocity, COAREWindVelocity, CoolSkinHumDepression, DeltaSatVapPress, DeltaTempInterface, &
                          CoolSkinTempDepression, ScailingVelocityPar, ScailingTemperaturePar, ScailingMoisturePar,     &
                          CoolSkinThickness, PBLHeight, Nits, i, j)
    
        !Arguments-------------------------------------------------------------     

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer        :: AirTemperature, PBLHeight
        real                                    :: VerticalWindVelocity, COAREWindVelocity10m, WindVelocity, Beta
        real                                    :: Grav, Charn, Von, Fdg, MoistureTransfCoeffNeutral, MOLength
        real                                    :: VelocityRoughnessLgtNeutral, VelocityTransfCoeffNeutral, VelocityTransfCoeff
        real                                    :: TemperatureTransfCoeffNeutral, TemperatureRoughnessLgtNeutral
        real                                    :: TemperatureTransfCoeff, CoeffC, SatRichNumb, RichNumb, StabilityParameter
        real, intent(OUT)                       :: ScailingVelocityPar, ScailingTemperaturePar, ScailingMoisturePar
        real, intent(OUT)                       :: CoolSkinThickness, DeltaSatVapPress, DeltaTempInterface, CoolSkinTempDepression
        real, intent(OUT)                       :: CoolSkinHumDepression, COAREWindVelocity
        integer, intent(IN)                     :: i, j
        real, intent(IN)                        :: DeltaSatVapPressCorr, Visa, InterfaceMoistureContent
        real, intent(IN)                        :: MoistureContentAir
        integer, intent(OUT)                    :: Nits
                                                 
        !Begin-----------------------------------------------------------------
            Beta    = 1.2     ! saline contraction coeff [1/psu]   parameter evaluated from Fairall low windspeed turbulence data
            Fdg     = 1.0     ! Fairall LKB roughness Reynolds number to Von Karman
            Grav    = 9.8     ! m/s2
            Von     = 0.4     ! von karman constant
            Charn   = 0.011 
            
        !****************************First Guess ***************************************************
            
            DeltaTempInterface            = Me%SurfaceTemperature(i,j) - AirTemperature(i, j) - 0.0098 * Me%AirMeasurementHeight 
            DeltaSatVapPress              = InterfaceMoistureContent - MoistureContentAir 
            VerticalWindVelocity          = 0.5 
            CoolSkinTempDepression        = 0.3  
            CoolSkinHumDepression         = DeltaSatVapPressCorr * CoolSkinTempDepression 
            COAREWindVelocity             = sqrt(WindVelocity**2.0 + VerticalWindVelocity**2.0) 
            COAREWindVelocity10m          = COAREWindVelocity * log(10.0 / 1.0e-4)/ log(Me%WindHeight / 1.0e-4) 
            ScailingVelocityPar           = 0.035 * COAREWindVelocity10m 
            VelocityRoughnessLgtNeutral   = Charn * ScailingVelocityPar**2.0 / Grav + 0.11 * Visa / ScailingVelocityPar 
            VelocityTransfCoeffNeutral    = (Von / log(10 / VelocityRoughnessLgtNeutral))**2.0 
            MoistureTransfCoeffNeutral    = 0.00115 
            TemperatureTransfCoeffNeutral = MoistureTransfCoeffNeutral / sqrt(VelocityTransfCoeffNeutral) 
            TemperatureRoughnessLgtNeutral= 10.0 / exp(Von / TemperatureTransfCoeffNeutral) 
            VelocityTransfCoeff           = (Von / log(Me%WindHeight / VelocityRoughnessLgtNeutral))**2.0 
            TemperatureTransfCoeff        = Von / log(Me%AirMeasurementHeight / TemperatureRoughnessLgtNeutral) 
            CoeffC                        = Von * Temperaturetransfcoeff / VelocityTransfCoeff      
            ! CoeffC is a analytical function of the standard bulk exchange coefficient
            SatRichNumb                   = -Me%WindHeight / PBLHeight(i, j) / 0.004 / Beta**3.0 
            RichNumb                      = -Grav * Me%WindHeight / (AirTemperature(i, j) + 273.16)   *  &
                                            ((DeltaTempInterface - CoolSkinTempDepression * Me%Jcool) +  &
                                            0.61 * (AirTemperature(i, j) + 273.16) * DeltaSatVapPress) / &
                                            COAREWindVelocity**2.0 
            Nits = 3
            
            if (RichNumb .LT. 0) then 
                StabilityParameter = CoeffC * RichNumb / (1.0 + RichNumb / SatRichNumb)  
                ! StabilityParameter = ratio of height to Obukhov length
            else 
                StabilityParameter = CoeffC * RichNumb * (1.0 + 27.0 / 9.0 * RichNumb / CoeffC)
            endif 
            
            MOLength = Me%WindHeight / StabilityParameter

            if (StabilityParameter > 50.0) then 
                Nits = 1 
            endif 
            
             ScailingVelocityPar    = COAREWindVelocity * Von /                                                 &
                                      (log(Me%WindHeight / VelocityRoughnessLgtNeutral)-Psiu(StabilityParameter))
             
             ScailingTemperaturePar = -(DeltaTempInterface - CoolSkinTempDepression * Me%Jcool) * Von * Fdg /  &
                                       (log(Me%AirMeasurementHeight / TemperatureRoughnessLgtNeutral) -     &
                                        Psit(Me%AirMeasurementHeight / MOLength))
             
             ScailingMoisturePar    = -(DeltaSatVapPress - DeltaSatVapPressCorr * CoolSkinTempDepression * Me%Jcool) * &
                                       Von * Fdg / (log(Me%AirMeasurementHeight / TemperatureRoughnessLgtNeutral) -    &
                                        Psit(Me%AirMeasurementHeight / MOLength))
             
             CoolSkinThickness      = 0.001
             
        !********************************End First Guess*************************************
    
    end subroutine FirstGuess
                          
    subroutine BulkLoop(PropLatentHeat, PropSensibleHeat, PropNetLongWaveRadiation, PropSurfaceRadiation,  &
                        AirTemperature, WaterDensity, DeltaTempInterface, ScailingTemperaturePar, Charn,   &
                        ScailingMoisturePar, DeltaSatVapPress, CoolSkinTempDepression, AirDensity, Visa,   &
                        CoolSkinHumDepression, CoolSkinThickness, Bigc, LatHeatVap, ScailingVelocityPar,   &
                        COAREWindVelocity, DeltaSatVapPressCorr, MoistureContentAir, WindVelocity,         &
                        PBLHeight, i, j, KUB, Nits)

        !Arguments-------------------------------------------------------------
    
        type(T_Property), pointer               :: PropLatentHeat
        type(T_Property), pointer               :: PropSensibleHeat
        type(T_Property), pointer               :: PropNetLongWaveRadiation
        type(T_Property), pointer               :: PropSurfaceRadiation        

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer        :: AirTemperature, PBLHeight
        real,    dimension(:,:,:), pointer      :: WaterDensity
        real                                    :: Be, Visw, Grav, Von, Cpw, Tcw, Cpa, Fdg, Beta, HeatLossWithoutRainFlux
        real                                    :: FracSolarRadAbs, MoistureRoughnessLgt, MOLength, RoughnessReynolds
        real                                    :: StabilityParameter, SurfaceBuoyancyFlux, TemperatureRoughnessLgt
        real                                    :: TotalCoolingSurface, SaundersCoeff
        real                                    :: VelocityRoughnessLgt, VerticalWindVelocity, VirtualSurfaceCooling
        real, intent(IN)                        :: DeltaTempInterface, Charn, DeltaSatVapPress, MoistureContentAir
        real, intent(IN)                        :: Bigc, Visa, AirDensity, LatHeatVap, DeltaSatVapPressCorr, WindVelocity
        real, intent(INOUT)                     :: ScailingVelocityPar, COAREWindVelocity, CoolSkinTempDepression
        real, intent(INOUT)                     :: CoolSkinHumDepression
        real, intent(INOUT)                     :: CoolSkinThickness, ScailingMoisturePar, ScailingTemperaturePar
        integer                                 :: n
        integer, intent(IN)                     :: i, j, KUB, Nits
        
        !Begin-----------------------------------------------------------------
           Beta = 1.2
           Be   = 0.026     ! salinity expansion coeff
           Visw = 1e-6      ! kinematic viscosity of sea-water [m^2/s]
           Grav = 9.8
           Von  = 0.4       ! von karman constant
           Cpa  = 1004.67   ! J/Kg.K
           Cpw  = 4000.0    ! J/Kg.K
           Visw = 1e-6      ! kinematic viscosity of sea-water [m^2/s]
           Tcw  = 0.6       ! thermal conductivity of sea-water [W/m/K]
           Fdg  = 1.0       ! Fairall LKB roughness Reynolds number to Von Karman
           
      !************************************************** Begin loop *********************************************************     
           do n = 1, Nits 
     
            StabilityParameter = (Von * Grav * Me%WindHeight / (AirTemperature(i, j) + 273.16))    *                  &
                                 (ScailingTemperaturePar * (1.0 + 0.61 * MoistureContentAir) + .61 *                  &
                                 (AirTemperature(i, j) + 273.16)                                   *                  &
                                 ScailingMoisturePar) / (ScailingVelocityPar**2.0) / (1.0 + 0.61 * MoistureContentAir) 

            VelocityRoughnessLgt    = Charn * ScailingVelocityPar**2.0 / Grav + 0.11 * Visa / ScailingVelocityPar  
            RoughnessReynolds       = VelocityRoughnessLgt * ScailingVelocityPar / Visa 
            MOLength                = Me%WindHeight / StabilityParameter 
            MoistureRoughnessLgt    = min(1.15e-4, 5.5e-5 / RoughnessReynolds**0.6) 
            TemperatureRoughnessLgt = MoistureRoughnessLgt
            
            ScailingVelocityPar     = COAREWindVelocity * Von /                      &      ! u*
                                      (log(Me%WindHeight / VelocityRoughnessLgt) -   &
                                      Psiu(Me%WindHeight / MOLength))
            
            ScailingTemperaturePar  = - (DeltaTempInterface - CoolSkinTempDepression * Me%Jcool) * Von * Fdg /     &     ! t*
                                       (log(Me%AirMeasurementHeight / TemperatureRoughnessLgt)               -     &
                                       Psit(Me%AirMeasurementHeight / MOLength))
            !q* 
            ScailingMoisturePar     = - (DeltaSatVapPress - DeltaSatVapPressCorr * CoolSkinTempDepression * Me%Jcool)* &
                                      Von * Fdg / (log(Me%AirMeasurementHeight / MoistureRoughnessLgt)          - &             
                                      Psit(Me%AirMeasurementHeight / MOLength)) 
            
            ! stability parameter paper Eq.3
            SurfaceBuoyancyFlux     = - Grav / (AirTemperature(i, j) + 273.16) * ScailingVelocityPar *        &
                                      (ScailingTemperaturePar + 0.61 * (AirTemperature(i, j) + 273.16) * ScailingMoisturePar)
            
            if (SurfaceBuoyancyFlux .GT. 0.0) then
                
            VerticalWindVelocity = Beta * (SurfaceBuoyancyFlux * PBLHeight(i, j))**0.333  !stability parameter paper Eq.16
            
            else
                
            VerticalWindVelocity = 0.2 
            
            endif
            
            COAREWindVelocity               = sqrt(WindVelocity**2.0 + VerticalWindVelocity**2.0)  
            ! wind velocity with gustiness included
            
            PropSensibleHeat%Field(i, j)    = AirDensity * Cpa * ScailingVelocityPar * ScailingTemperaturePar 
            ! alterado para ponto de vista da água
            
            PropLatentHeat%Field(i, j)      = AirDensity * LatHeatVap * ScailingVelocityPar * ScailingMoisturePar  
            ! alterado para ponto de vista da água
            
            HeatLossWithoutRainFlux         = - PropNetLongWaveRadiation%Field(i, j) - PropSensibleHeat%Field(i, j) - &
                                              PropLatentHeat%Field(i, j)
            
            FracSolarRadAbs                 = PropSurfaceRadiation%Field(i, j) * (0.065 + 11.0 * CoolSkinThickness - 6.6e-5 /  &
                                              CoolSkinThickness * (1.0 - exp(-CoolSkinThickness / 8.0e-4))) ! Eq.16 coolskin
            
            TotalCoolingSurface             = HeatLossWithoutRainFlux - FracSolarRadAbs !coolskin
            
                  ! W/m2.K                       1/k     * w/m2                -              W/m2               * (J/Kg.K) / (J/Kg)
            VirtualSurfaceCooling          = Me%Al(i, j) * TotalCoolingSurface - Be * PropLatentHeat%Field(i, j) * Cpw / LatHeatVap
            ! Eq. 7 coolskin

            if (VirtualSurfaceCooling .GT. 0) then 
                 
            SaundersCoeff         = 6.0 / (1.0 + (Bigc * VirtualSurfaceCooling / ScailingVelocityPar**4.0)**0.75)**0.333
            ! Eq 13 coolskin
            CoolSkinThickness     = SaundersCoeff * Visw / (sqrt(AirDensity / WaterDensity(i, j, KUB)) * ScailingVelocityPar)
            !Eq.11 coolskin

            else
            SaundersCoeff         = 6.0 
            CoolSkinThickness     = min(0.01, SaundersCoeff * Visw /                                  &
                                    (sqrt(AirDensity / WaterDensity(i, j, KUB)) * ScailingVelocityPar))
            !Eq.11 coolskin
            endif 
     
            CoolSkinTempDepression    = TotalCoolingSurface * CoolSkinThickness / Tcw   
            CoolSkinHumDepression     = DeltaSatVapPressCorr * CoolSkinTempDepression 
   
           enddo !bulk iter loop    
    
    end subroutine BulkLoop

    

    subroutine CheckRadiationOptions (PropNetLongWaveRadiation, PropUpLongWaveRadiation, PropDownLongWaveRadiation, &
                                      PropSurfaceRadiation)
    
        !Arguments ------------------------------------------------------------
        type(T_Property), pointer                   :: PropNetLongWaveRadiation
        type(T_Property), pointer                   :: PropUpLongWaveRadiation
        type(T_Property), pointer                   :: PropDownLongWaveRadiation
        type(T_Property), pointer                   :: PropSurfaceRadiation
        
        !Local -----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin -------------------------------------------------------------------
        
        if (PropNetLongWaveRadiation%ID%SolutionFromFile .OR. PropNetLongWaveRadiation%Constant) then
            
                write(*,*) 'Net Long Wave Radiation must be computed when using COARE-based algorithm.'
                stop 'CheckCOARELongWaveOptions - ModuleInterfaceWaterAir - ERR01'
        endif
            
        if (PropUpLongWaveRadiation%ID%SolutionFromFile .OR. PropUpLongWaveRadiation%Constant) then

                write(*,*) 'Upward Long Wave Radiation must be computed when using COARE-based algorithm.'
                stop 'CheckCOARELongWaveOptions - ModuleInterfaceWaterAir - ERR02'

        endif
         
        if (PropDownLongWaveRadiation%ID%SolutionFromFile) then

                call ModifyFillMatrix(FillMatrixID      = PropDownLongWaveRadiation%ID%ObjFillMatrix,&
                                      Matrix2D          = PropDownLongWaveRadiation%Field,           &
                                      PointsToFill2D    = Me%ExtWater%WaterPoints2D,                     &
                                      STAT              = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'CheckCOARELongWaveOptions - ModuleInterfaceWaterAir - ERR03'

        elseif (.not. PropDownLongWaveRadiation%Constant) then

                call ComputeDownLongWaveRad(PropDownLongWaveRadiation)

        endif
        
        if (PropSurfaceRadiation%ID%SolutionFromFile .OR. PropSurfaceRadiation%Constant) then
                write(*,*) 'Surface Radiation must be computed when using COARE-based algorithm.'
                stop 'CheckCOARELongWaveOptions - ModuleInterfaceWaterAir - ERR04'
        endif
        
    end subroutine CheckRadiationOptions
    
    subroutine CheckLatentSensibleOptions(PropLatentHeat, PropSensibleHeat)

        !Arguments ------------------------------------------------------------
        type(T_Property), pointer                   :: PropLatentHeat
        type(T_Property), pointer                   :: PropSensibleHeat
        
        !Local -----------------------------------------------------------------
        
        !Begin -------------------------------------------------------------------
        
        if (PropLatentHeat%ID%SolutionFromFile .OR. PropLatentHeat%Constant) then
            
                write(*,*) 'Latent heat must be computed when using COARE-based algorithm.'
                stop 'CheckLatentSensibleOptions - ModuleInterfaceWaterAir - ERR01'
        endif
            
         if (PropSensibleHeat%ID%SolutionFromFile .OR. PropSensibleHeat%Constant) then

                write(*,*) 'Sensible heat must be computed when using COARE-based algorithm.'
                stop 'CheckLatentSensibleOptions - ModuleInterfaceWaterAir - ERR02'

         endif        
    
    end subroutine CheckLatentSensibleOptions
    
    real function Psit (zet)
    
        !Arguments-------------------------------------------------------------
        real, intent(IN)              :: zet
        !Local -------------------------------------------------
        
        real                          :: F, C, x, Psik, Psic
        !begin -----------------------------------------------------------
        
        if(zet > 0)then
            
            C       = min(50.0 , 0.35 * zet) 
            Psit = - ((1. + 2. / 3. * zet)**1.5 + .6667 * (zet - 14.28) / exp(c) + 8.525)
        else
            
            x       = (1. - (15 * zet))**.5
            Psik    = 2 * log((1 + x) / 2) 
            x       = (1. - (34.15 * zet))**.333
    !       x1      = (1. + 2. * x)  / sqrt(3.)
    !       A       = AIMAG(x1)
    !       B       = REAL(x1)        
            Psic    = 1.5 * log((1. + x + x**2) / 3.) - sqrt(3.) * atan((1. + 2. * x)  / sqrt(3.)) + 4. * atan(1.) / sqrt(3.) 
            F       = zet * zet / (1 + zet**2) 
            Psit    = (1. - f) * psik + f * psic               
            
        endif     
    
    end function Psit
    
    real function Psiu (zet)
    
        !Arguments-------------------------------------------------------------
        real, intent(IN)              :: zet
        !Local ----------------------------------------------------------------
        real                          :: X, Psik, Psic, F, C
        !begin ----------------------------------------------------------------
        
        if(zet>0)then
            
            C     = min(50.0 , 0.35 * zet) 
            Psiu  = -((1. + 1.0 * zet) + .667 * (zet - 14.28) / exp(c) + 8.525)
          
        else
            X      = (1. - 15. * zet)**.25 
            Psik   = 2. * log((1. + x) / 2.) + log((1. + x**2) / 2.) - 2. * atan(x) + 2. * atan(1.) 
            X      = (1. - 10.15 * zet)**.3333 
            Psic   = 1.5 * log((1. + x + x**2) / 3.) - sqrt(3.) * atan((1. + 2. * x) / sqrt(3.)) + 4. * atan(1.) / sqrt(3.) 
            F      = zet**2 / (1 + zet**2) 
            Psiu   = (1 - f) * psik + f * psic    
        endif
    end function Psiu
    
    real function AirViscosity (Temperature)
    
    !Arguments ----------------------------------------------
    Real                         :: Temperature
    
    !Begin --------------------------------------------------
    
     AirViscosity = 1.326e-5 * (1.0 + (6.542e-3 * Temperature) + (8.301e-6 * Temperature**2) - &
                    (4.84e-9 * Temperature**3))
     
    end function AirViscosity
    
    real function WaterVapourDiffusivity(Temperature)
    
    !Arguments ----------------------------------------------
    Real                         :: Temperature
    
    !Begin --------------------------------------------------
    
    WaterVapourDiffusivity = 2.11e-5 * ((Temperature + 273.16) / 273.16)**1.94
    
    end function WaterVapourDiffusivity
    
    real function HeatDiffusivity(Temperature, AirDensity, Cpa)
    
    !Arguments ----------------------------------------------
    Real                         :: Temperature, AirDensity, Cpa
    
    !Begin -------------------------------------------------- 
    
    HeatDiffusivity = (1. + 3.309e-3 * Temperature - 1.44e-6 * Temperature**2) * 0.02411 / (AirDensity * Cpa)
    
    end function HeatDiffusivity
       
    !--------------------------------------------------------------------------
    
    subroutine ModifyWindStress(PropWindStress)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindStress

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterfaceWaterAir", "ModifyWindStress")

        if (PropWindStress%ID%SolutionFromFile) then

            call ModifyFillMatrixVectorial(FillMatrixID      = PropWindStress%ID%ObjFillMatrix,     &
                                  Matrix2DU         = PropWindStress%FieldU,               &
                                  Matrix2DV         = PropWindStress%FieldV,               &
                                  Matrix2DX         = PropWindStress%FieldX,               &
                                  Matrix2DY         = PropWindStress%FieldY,               & 
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,           &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR01'

        endif

        !call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid,      &
        !                             VectorInX         = PropWindStressX%Field,     &
        !                             VectorInY         = PropWindStressY%Field,     &
        !                             VectorOutX        = PropWindStressX%FieldGrid, &
        !                             VectorOutY        = PropWindStressY%FieldGrid, &   
        !                             WaterPoints2D     = Me%ExtWater%WaterPoints2D, &
        !                             RotateX           = .true.,                    &
        !                             RotateY           = .true.,                    &
        !                             STAT              = STAT_CALL)
        !
        !if(STAT_CALL .ne. SUCCESS_) stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR30'


        if (.not. PropWindStress%Constant            .and. &
            .not. PropWindStress%Constant            .and. &
            .not. PropWindStress%ID%SolutionFromFile .and. &
            .not. PropWindStress%ID%SolutionFromFile ) then

            call ComputeTauWind (PropWindStress)

            !need to rotate field computed in ComputeTauWind since we are outside FillMatrix
            call RotateVectorGridToField(HorizontalGridID  = Me%ObjHorizontalGrid,      &
                                         VectorInX         = PropWindStress%FieldU,     &
                                         VectorInY         = PropWindStress%FieldV,     &
                                         VectorOutX        = PropWindStress%FieldX,     &
                                         VectorOutY        = PropWindStress%FieldY,     &   
                                         WaterPoints2D     = Me%ExtWater%WaterPoints2D, &
                                         RotateX           = .true.,                    &
                                         RotateY           = .true.,                    &
                                         STAT              = STAT_CALL)

        endif

        if (MonitorPerformance) call StopWatch ("ModuleInterfaceWaterAir", "ModifyWindStress")

    end subroutine ModifyWindStress

    !--------------------------------------------------------------------------
    
    subroutine ComputeTauWind (PropWindStress)

        !Arguments------------------------------------------------------------
        type(T_Property), pointer                   :: PropWindStress

        !Local----------------------------------------------------------------
        real,    dimension(:,:), pointer            :: UWIND, VWIND
        type(T_Property), pointer                   :: WindShearVelocity
        integer                                     :: IUB, ILB, JUB, JLB, i, j
        real                                        :: VM, CDWIND
        real                                        :: Coef
        integer                                     :: CHUNK, STAT_CALL

        !Begin----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB

        UWIND   => Me%LocalAtm%WindVelocityU%Field
        VWIND   => Me%LocalAtm%WindVelocityV%Field

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ComputeTauWind")
        endif

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,VM,Coef,CDWIND)
cd1:    if(Me%CDWINDmethod == Constant)then
        
            !The next line is done by each thread with OpenMP compilation.
            !This is ok because it is a private variable.
            CDWIND = Me%CDWIND
        
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                    
                    !Compute the square root 
                    VM = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)
                    
                    Coef = CDWIND * VM * Air_Density
                    
                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStress%FieldU(I,J) = Coef * UWIND(I,J)
                    
                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStress%FieldV(I,J) = Coef * VWIND(I,J)

                endif
                
            enddo
            enddo
            !$OMP END DO

        else if (Me%CDWINDmethod == WindFunction) then cd1
        
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                    
                    !Compute the velocity modulus
                    VM = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)  
                    
                    !Compute the wind drag coefficient based on Large & Pond, 1981 
                    !included the modifications for low wind speeds considered in Trenberth et al. (1990) Francisco

                    if (VM <= 1.)then

                        CDWIND = 0.00218
                        
                    elseif (VM >  1.  .and. VM <  3.)then

                        CDWIND = (0.00062 + 0.00156 / VM)                         
                    
                    elseif (VM >=  3.  .and. VM <  10.)then

                        CDWIND = 0.00114

                    elseif (VM >= 10. .and. VM <= 26.)then

                        CDWIND = (0.00049 + 0.000065 * VM)                   

                    elseif (VM > 26.) then
                        
                        !The Large & Pond, 1981 formulation is not valid for wind higher than 26 m/s
                        !A constant value is assumed   
                        CDWIND = 0.00218

                    endif
                    
                    Coef = CDWIND * VM * Air_Density

                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStress%FieldU(I,J) = Coef * UWIND(I,J)

                    !Mellor, Introduction to Physical Oceanography, p52 (1996)---------
                    PropWindStress%FieldV(I,J) = Coef * VWIND(I,J)

                endif

            enddo
            enddo
            !$OMP END DO

        else if (Me%CDWINDmethod == ShearVelocity) then cd1        
        
            call Search_Property(WindShearVelocity, WindShearVelocity_, PrintWarning = .true., STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ComputeTauWind - ModuleInterfaceWaterAir - ERR10'

            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j=JLB, JUB
            do i=ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                    
                    !Compute the velocity modulus
                    VM = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)  
                    
                    ![M/L^2/T] = [L/T]^2/[L/T] * [M/L^3]
                    Coef = WindShearVelocity%Field(i,j)**2 / VM * Air_Density
                    
                    ![M*L/T^2/L^2] = [M/L^2/T]*[L/T]
                    PropWindStress%FieldU(I,J) = Coef * UWIND(I,J)
                    PropWindStress%FieldV(I,J) = Coef * VWIND(I,J)

                endif

            enddo
            enddo
            !$OMP END DO
        endif cd1
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ComputeTauWind")
        endif

        nullify(UWIND)
        nullify(VWIND)

    end subroutine ComputeTauWind        

    !--------------------------------------------------------------------------

    subroutine ModifyTurbulentKE(PropTurbulentKE)

        !Arguments-------------------------------------------------------------
        
        type(T_Property), pointer :: PropTurbulentKE

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (PropTurbulentKE%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = PropTurbulentKE%ID%ObjFillMatrix,   &
                                  Matrix2D          = PropTurbulentKE%Field,              &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,          &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyTurbulentKE - ModuleInterfaceWaterAir - ERR02'

        else if (.not. PropTurbulentKE%Constant) then

            call ComputeTKEWind (PropTurbulentKE) 

        endif

    end subroutine ModifyTurbulentKE

    !--------------------------------------------------------------------------

    subroutine ComputeTKEWind (PropTurbulentKE)

        !Arguments------------------------------------------------------------
        type(T_Property),       pointer             :: PropTurbulentKE

        !Local----------------------------------------------------------------
        real, dimension(:,:),   pointer             :: UWIND, VWIND
        integer                                     :: IUB, ILB, JUB, JLB, i, j
        real                                        :: VM
        integer                                     :: CHUNK

        !Begin----------------------------------------------------------------

        IUB = Me%WorkSize2D%IUB
        ILB = Me%WorkSize2D%ILB
        JUB = Me%WorkSize2D%JUB
        JLB = Me%WorkSize2D%JLB

        UWIND   => Me%LocalAtm%WindVelocityU%Field
        VWIND   => Me%LocalAtm%WindVelocityV%Field

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ComputeTKEWind")
        endif
    
        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j,VM)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do3:    do j=JLB, JUB
do4:    do i=ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Compute the square root 
                VM     = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)       

                ![m/s2 * kg * m / m2 / s]         = [kg/m^3 * (m/s)^3]     
                PropTurbulentKE%Field(I,J) = (Air_Density * CDE * VM ** 3.0) 
            endif

        enddo do4
        enddo do3
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ComputeTKEWind")
        endif

        nullify(UWIND)
        nullify(VWIND)

    end subroutine ComputeTKEWind 
    
  !-------------------------------------------------------------------------- 

    subroutine ModifyOxygenFlux(PropertyO2)

        !Arguments--------------------------------------------------------------
        type (T_Property),        pointer       :: PropertyO2

        !Local------------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Begin------------------------------------------------------------------

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyOxygenFlux - ModuleInterfaceWaterAir - ERR01'

        call ModifyAerationFlux(PropertyO2%Field)

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyOxygenFlux - ModuleInterfaceWaterAir - ERR02'

    end subroutine ModifyOxygenFlux

    !-------------------------------------------------------------------------- 

    subroutine ModifyAerationFlux(Flux)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: Flux

        !Local------------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: UWIND, VWIND
        real,    dimension(:,:,:),    pointer   :: WaterTemperature
        integer                                 :: i,j
        real                                    :: WindVelocity, WaterTemp
        integer                                 :: ILB, IUB, JLB, JUB, KUB
        
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB 
        KUB = Me%WorkSize3D%KUB

        !Surface temperature of the waterbody
        WaterTemperature => Me%ExtWater%WaterTemperature
        UWIND            => Me%LocalAtm%WindVelocityU%Field
        VWIND            => Me%LocalAtm%WindVelocityV%Field

        !The formulation is taken from CE-QUAL-W2 v3.1
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Module of the wind velocity 
                WindVelocity = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)
                
                WaterTemp = WaterTemperature(i, j, KUB)
                
                Flux(i,j) = AerationFlux(Me%AerationEquation, WindVelocity, WaterTemp)
               
            endif

        end do
        end do
         

        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(WaterTemperature)
        
    end subroutine ModifyAerationFlux
    

    !---------------------------------------------------------------------------
    
    subroutine ModifyCarbonDioxideFlux(PropCarbonDioxide)

        !Arguments--------------------------------------------------------------
        type (T_Property),         pointer :: PropCarbonDioxide

        !External---------------------------------------------------------------
        integer                            :: STAT_CALL

        !Local------------------------------------------------------------------
        integer                            :: ILB, IUB, JLB, JUB, KUB
        real,    dimension(:,:,:), pointer :: WaterTemperature
        integer                            :: i, j
        integer                            :: CHUNK
        
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB 
        KUB = Me%WorkSize3D%KUB

        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterTemperature,     &
                              PropertyXIDNumber = Temperature_,                     &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyCarbonDioxideFlux - ModuleInterfaceWaterAir - ERR10'


        If (Me%Coupled%CEQUALW2%Yes) then
        
            !Surface temperature of the waterbody
            WaterTemperature   => Me%ExtWater%WaterTemperature

            call ModifyCO2AerationFlux(PropCarbonDioxide%Field)

            if (MonitorPerformance) then
                call StartWatch ("ModuleInterfaceWaterAir", "ModifyCarbonDioxideFlux")
            endif

            CHUNK = CHUNK_J(JLB, JUB)
            !$OMP PARALLEL PRIVATE(i,j)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                    PropCarbonDioxide%Field(i,j)= PropCarbonDioxide%Field(i,j) * 0.923              * &
                                                  (0.286 * exp(-0.0314*(WaterTemperature(i,j,KUB))) * &
                                                   Me%AltitudeCorrection)

                endif

            enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            if (MonitorPerformance) then
                call StopWatch ("ModuleInterfaceWaterAir", "ModifyCarbonDioxideFlux")
            endif         

            call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModuleInterfaceWaterAir - ERR040'


            !Nullify auxiliar variables
            nullify(WaterTemperature)
        
        else
        
            call ModifyCO2AerationFlux(PropCarbonDioxide%Field)

            call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterTemperature, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModuleInterfaceWaterAir - ERR050'
 
        
        endif
        

    end subroutine ModifyCarbonDioxideFlux

   !-------------------------------------------------------------------------- 
    
    
    subroutine ModifyCO2AerationFlux(Flux)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: Flux

        !Local------------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: UWIND, VWIND
        real,    dimension(:,:,:),    pointer   :: WaterTemperature, WaterSalinity 
        real,    dimension(:,:,:),    pointer   :: Velocity_U, Velocity_V
        integer                                 :: i,j
        real                                    :: WindVelocity, WaterVel, WaterTemp, WaterSal
        integer                                 :: ILB, IUB, JLB, JUB, KUB, STAT_CALL
        
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB 
        IUB = Me%WorkSize2D%IUB 
        JLB = Me%WorkSize2D%JLB 
        JUB = Me%WorkSize2D%JUB 
        KUB = Me%WorkSize3D%KUB

        !Surface temperature of the waterbody
        WaterTemperature => Me%ExtWater%WaterTemperature
        UWIND            => Me%LocalAtm%WindVelocityU%Field
        VWIND            => Me%LocalAtm%WindVelocityV%Field
        
        
        call GetHorizontalVelocity(HydrodynamicID      = Me%ObjHydrodynamic,        &
                                   Velocity_U          = Velocity_U,                &
                                   Velocity_V          = Velocity_V,                &
                                   STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR10'
        
        
        call GetConcentration(WaterPropertiesID = Me%ObjWaterProperties,            &
                              ConcentrationX    = Me%ExtWater%WaterSalinity,        &
                              PropertyXIDNumber = Salinity_,                        &
                              STAT              = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR15'

        WaterSalinity    => Me%ExtWater%WaterSalinity

        
        !The formulation is taken from CE-QUAL-W2 v3.1
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                !Module of the wind velocity 
                WindVelocity = sqrt(UWIND(i,j)**2. + VWIND(i,j)**2.)
                
                WaterTemp = WaterTemperature(i, j, KUB)
                
                WaterSal = WaterSalinity(i, j, KUB)
                
                WaterVel  = SQRT(Velocity_U(i, j, KUB)**2. + Velocity_V(i, j, KUB)**2.)
                
                Flux(i,j) = AerationFlux_CO2(Me%CO2AerationEquation, WindVelocity, WaterVel, WaterTemp, WaterSal)
               
            endif

        end do
        end do
 
             
        call UngetHydrodynamic(Me%ObjHydrodynamic, Velocity_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR020'    
        
        call UngetHydrodynamic(Me%ObjHydrodynamic, Velocity_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR030'     

        call UnGetWaterProperties(Me%ObjWaterProperties, Me%ExtWater%WaterSalinity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyCarbonDioxideFlux - ModifyCO2AerationFlux - ERR040'

        !Nullify auxiliar variables
        nullify(UWIND)
        nullify(VWIND)
        nullify(Velocity_U)
        nullify(Velocity_V)
        nullify(WaterTemperature)
        nullify(WaterSalinity)

        
    end subroutine ModifyCO2AerationFlux
    
    !--------------------------------------------------------------------------  
    !---------------------------------------------------------------------------
    
    subroutine ModifyAmmoniaFlux(PropAmmoniaFlux)

        !Arguments--------------------------------------------------------------
        type (T_Property),         pointer :: PropAmmoniaFlux

        !External---------------------------------------------------------------
        integer                            :: STAT_CALL

        !Local------------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: AtmospDepositionR
        integer                            :: i, j
        integer                            :: ILB, IUB, JLB, JUB
        !Begin------------------------------------------------------------------
       
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

        
       call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AtmospDeposReduNH4%Field, &
                                   ID           = AtmospDeposReduNH4_,                &
                                   STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ModifyAmmoniaFlux - ModuleInterfaceWaterAir - ERR10'
        
        AtmospDepositionR  => Me%ExtAtm%AtmospDeposReduNH4%Field
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                
                PropAmmoniaFlux%Field(i,j) = AtmospDepositionR(i,j)
               
            endif

        end do
        end do  
        
        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AtmospDeposReduNH4%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyAmmoniaFlux - ModuleInterfaceWaterAir - ERR20'

    end subroutine ModifyAmmoniaFlux

   !-------------------------------------------------------------------------- 
  

   
   
   subroutine ModifyNitrateFlux(PropNitrateFlux)

        !Arguments--------------------------------------------------------------
        type (T_Property),         pointer :: PropNitrateFlux

        !External---------------------------------------------------------------
        integer                            :: STAT_CALL

        !Local------------------------------------------------------------------
        real,    dimension(:,:  ),    pointer   :: AtmospDepositionO
        integer                            :: i, j
        integer                            :: ILB, IUB, JLB, JUB
        !Begin------------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

        call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                   Scalar       = Me%ExtAtm%AtmospDeposOxidNO3%Field, &
                                   ID           = AtmospDeposOxidNO3_,                &
                                   STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)stop 'ModifyNitrateFlux - ModuleInterfaceWaterAir - ERR10'
        
        AtmospDepositionO  => Me%ExtAtm%AtmospDeposOxidNO3%Field
       
        
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
                
                PropNitrateFlux%Field(i,j) = AtmospDepositionO(i,j)
               
            endif

        end do
        end do  
        
        call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AtmospDeposOxidNO3%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyNitrateFlux - ModuleInterfaceWaterAir - ERR20'

    end subroutine ModifyNitrateFlux

   !-------------------------------------------------------------------------- 
   

    subroutine ModifyWindShearVelocity

        !Arguments--------------------------------------------------------------


        !External--------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, i, j
        real                                :: WindStressModule
        type(T_Property), pointer           :: WindShearVelocity, WindStress
        integer                             :: STAT_CALL
        integer                             :: CHUNK
        logical                             :: FromFile

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ModifyWindShearVelocity")
        endif
        
        call Search_Property(WindShearVelocity, WindShearVelocity_, STAT = STAT_CALL) 

        if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= NOT_FOUND_ERR_) then               
            stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR10'
        endif
            
        if (STAT_CALL == SUCCESS_) then
            FromFile = WindShearVelocity%ID%SolutionFromFile
        else
            FromFile = .false.
        endif
       
        if (FromFile) then

            call ModifyFillMatrix(FillMatrixID      = WindShearVelocity%ID%ObjFillMatrix,   &
                                  Matrix2D          = Me%WindShearVelocity,                 &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyWindStress - ModuleInterfaceWaterAir - ERR20'
                
        else 

            call Search_Property(WindStress, WindStress_, .true., STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModifyWindShearVelocity - ModuleInterfaceWaterAir - ERR30'

            !call Search_Property(WindStressY, WindStressY_, .true., STAT = STAT_CALL) 
            !if (STAT_CALL /= SUCCESS_) stop 'ModifyWindShearVelocity - ModuleInterfaceWaterAir - ERR40'

            IUB = Me%WorkSize2D%IUB
            ILB = Me%WorkSize2D%ILB
            JUB = Me%WorkSize2D%JUB
            JLB = Me%WorkSize2D%JLB            

            CHUNK = CHUNK_J(JLB, JUB)
            !$OMP PARALLEL PRIVATE(i,j,WindStressModule)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
    do1:    do j=JLB, JUB
    do2:    do i=ILB, IUB
                
                if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then

                    WindStressModule                  = sqrt(WindStress%FieldU(i, j)**2. + &
                                                             WindStress%FieldV(i, j)**2.)

                   Me%WindShearVelocity(i, j) = sqrt(WindStressModule / SigmaDensityReference)
                                                                    
                endif

            enddo do2
            enddo do1
            !$OMP END DO
            !$OMP END PARALLEL            
        
        endif



        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ModifyWindShearVelocity")
        endif

    end subroutine ModifyWindShearVelocity

    !------------------------------------------------------------------------------

    subroutine ComputeWavesRugosity

        !Arguments-------------------------------------------------------------
#ifndef _WAVES_
        !Local-----------------------------------------------------------------
! Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! real, dimension(:,:), pointer               :: WaveHeight
        real, dimension(:,:), pointer               :: WaveHeightForRugosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
        integer                                     :: ILB, IUB, JLB, JUB, i, j
        integer                                     :: STAT_CALL
        integer                                     :: CHUNK

        !Begin-----------------------------------------------------------------

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB

!Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Gets wave height
        !call GetWaves (Me%ObjWaves, WaveHeight = WaveHeight, STAT = STAT_CALL)

        !Gets wave height used for rugosity computation
        if (Me%Rugosity%WavesFunction == SurfaceRugosityFromHS) then
            call GetWaves (Me%ObjWaves, WaveHeight = WaveHeightForRugosity, STAT = STAT_CALL)
        endif

        if (Me%Rugosity%WavesFunction == SurfaceRugosityFromHSW) then
            call GetWaves (Me%ObjWaves, BreakingWaveHeight = WaveHeightForRugosity, STAT = STAT_CALL)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        if (STAT_CALL /= SUCCESS_) stop 'ComputeWavesRugosity - InterfaceWaterAir - ERR10'

        if (MonitorPerformance) then
            call StartWatch ("ModuleInterfaceWaterAir", "ComputeWavesRugosity")
        endif

        CHUNK = CHUNK_J(JLB, JUB)
        !$OMP PARALLEL PRIVATE(i,j)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExtWater%WaterPoints2D(i, j) == WaterPoint) then
! Modified by Matthias DELPEY - 17/08/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 15/12/2011
                ! Me%Rugosity%Field(i, j) = Me%Rugosity%WavesRelation * WaveHeight(i, j)
                if (WaveHeightForRugosity(i,j) > 0.) then
                    Me%Rugosity%Field(i, j) = Me%Rugosity%WavesRelation * WaveHeightForRugosity(i, j)
                else
                    Me%Rugosity%Field(i, j) = Me%Rugosity%WavesRelation * 0.01
                endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            endif
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        if (MonitorPerformance) then
            call StopWatch ("ModuleInterfaceWaterAir", "ComputeWavesRugosity")
        endif
        
        call UnGetWaves (Me%ObjWaves, WaveHeightForRugosity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWavesRugosity - InterfaceWaterAir - ERR20'
    
#endif
    end subroutine ComputeWavesRugosity


    !--------------------------------------------------------------------------
    
    subroutine ModifyRugosity

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        if (Me%Rugosity%ID%SolutionFromFile) then

            call ModifyFillMatrix(FillMatrixID      = Me%Rugosity%ID%ObjFillMatrix,     &
                                  Matrix2D          = Me%Rugosity%Field,                &
                                  PointsToFill2D    = Me%ExtWater%WaterPoints2D,            &
                                  STAT              = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ModifyRugosity - ModuleInterfaceWaterAir - ERR01'

        elseif(.not. Me%Rugosity%Constant)then

 ! Modified by Matthias DELPEY - 15/12/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! if (Me%Rugosity%WavesFunction) then
            if (Me%Rugosity%WavesFunction == SurfaceRugosityFromHS  .or.                &
                Me%Rugosity%WavesFunction == SurfaceRugosityFromHSW   ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                call ComputeWavesRugosity

            endif

        endif
    
    end subroutine ModifyRugosity

    !--------------------------------------------------------------------------
    
! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 28/11/2011 - 16/12/2011
    
    subroutine ModifyWaveFluxTKE

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: DensityReferenceJustHere
        real, dimension(:,:), pointer               :: WaveSurfaceFluxTKE

        !Begin-----------------------------------------------------------------

        DensityReferenceJustHere = 1026.2

        ! Importation of the field from Module Waves
        call GetWaves (Me%ObjWaves, WaveSurfaceFluxTKE = WaveSurfaceFluxTKE, STAT = STAT_CALL)

        
        ! Conversion to right unit to be used in Module TurbGOTM
        Me%WaveFluxTKE%Field(:,:) = WaveSurfaceFluxTKE(:,:) / DensityReferenceJustHere


        call UnGetWaves (Me%ObjWaves,  WaveSurfaceFluxTKE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaveFluxTKE - InterfaceWaterAir - ERR01a'
    
    end subroutine ModifyWaveFluxTKE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX

        !Begin-----------------------------------------------------------------
        
        PropertyX   => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                
                !vectorial property - need to get data in user referential - X and Y
!~                 if (Check_Vectorial_Property(PropertyX%ID%IDNumber)) then 
                if (PropertyX%ID%IsVectorial) then

                    call WriteTimeSerie(Me%ObjTimeSerie, Data2D = PropertyX%FieldX, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleInterfaceWaterAir - ERR010'
            
                    call WriteTimeSerie(Me%ObjTimeSerie, Data2D = PropertyX%FieldY, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_TimeSeries - ModuleInterfaceWaterAir - ERR020'                      
                
                else
                    
                    call WriteTimeSerie(Me%ObjTimeSerie,                            &
                                        Data2D  = PropertyX%Field,                  &
                                        STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                      &
                        stop 'OutPut_TimeSeries - ModuleInterfaceWaterAir - ERR01'
                    
                endif
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
                    call StartWatch ("ModuleInterfaceWaterAir", "OutPut_BoxTimeSeries")
                endif

                CHUNK = CHUNK_J(JLB, JUB)
                !$OMP PARALLEL PRIVATE(I,J)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do J = JLB, JUB
                do I = ILB, IUB

                    if(Me%ExtWater%WaterPoints2D(i,j) == WaterPoint)then

                        Me%Scalar2D(i,j) = PropertyX%Field(i, j) * Me%ExternalVar%GridCellArea(i, j)
                    
                    endif

                end do
                end do
                !$OMP END DO
                !$OMP END PARALLEL
                
                if (MonitorPerformance) then
                    call StopWatch ("ModuleInterfaceWaterAir", "OutPut_BoxTimeSeries")
                endif
                               
                call BoxDif(Me%ObjBoxDif,                           &
                            Me%Scalar2D,                            &
                            "WaterAir "//trim(PropertyX%ID%name),   &
                            Me%ExtWater%WaterPoints2D,              &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  &
                    stop 'OutPut_BoxTimeSeries - ModuleInterfaceWaterAir - ERR01'

                Me%Scalar2D(:,:) = null_real

            endif
            PropertyX=>PropertyX%Next
        enddo

    end subroutine OutPut_BoxTimeSeries

    
    !--------------------------------------------------------------------------


    subroutine OutPut_Results_HDF
        
        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
        real                               :: Year, Month, Day, Hour, Minute, Second
         
        !Local-----------------------------------------------------------------
        type (T_Property), pointer         :: PropertyX
        integer                            :: OutPutNumber
        integer, dimension(6    )          :: TimeAux
        real,    dimension(6    ), target  :: AuxTime
        real,    dimension(:    ), pointer :: TimePtr
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                            :: WorkKLB, WorkKUB
        real(8)                            :: TotalTime, AuxPeriod
        type (T_Time)                      :: Aux
        character(LEN = StringLength)      :: PropertyNameX, PropertyNameY
        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize3D%ILB 
        WorkIUB = Me%WorkSize3D%IUB 
        WorkJLB = Me%WorkSize3D%JLB 
        WorkJUB = Me%WorkSize3D%JUB
        WorkKLB = Me%WorkSize3D%KLB 
        WorkKUB = Me%WorkSize3D%KUB


        OutPutNumber = Me%OutPut%NextOutPut


TOut:   if (Me%ExternalVar%Now >= Me%OutPut%OutTime(OutPutNumber)) then


                if (Me%ExternalVar%BackTracking) then
                    OutPutNumber = Me%OutPut%Number - OutPutNumber + 1 
                endif 
                
                
                if (Me%ExternalVar%BackTracking) then  
                    TotalTime = Me%EndTime         - Me%BeginTime                  
                    AuxPeriod = Me%ExternalVar%Now - Me%BeginTime
                    AuxPeriod = TotalTime          - AuxPeriod
                    
                    Aux = Me%BeginTime + AuxPeriod
                else
                    Aux = Me%ExternalVar%Now
                endif   
                

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
                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR01'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR02'

                PropertyX => Me%FirstProperty

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                        &
                                      Me%WorkSize2D%ILB,                                 &
                                      Me%WorkSize2D%IUB,                                 &
                                      Me%WorkSize2D%JLB,                                 &
                                      Me%WorkSize2D%JUB,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR03'


PropX:          do while (associated(PropertyX))

                    if(PropertyX%OutputHDF)then

                        if (PropertyX%ID%IDNumber /= SpecificOxygenFlux_ ) then
                            
                            !vectorial property - need to get data in user referential - X and Y
!~                             if (Check_Vectorial_Property(PropertyX%ID%IDNumber)) then
                            if (PropertyX%ID%IsVectorial) then
                                
                                !get the correct names of the properties
                                call Get_Vectorial_PropertyNames(PropertyX%ID%IDNumber, PropertyNameX, PropertyNameY)    
                                
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyNameX,     &
                                                     PropertyNameX, PropertyX%ID%Units,          &
                                                     Array2D      = PropertyX%FieldX,            &
                                                     OutputNumber = OutPutNumber,                &
                                                     STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                                       &
                                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR06a'     
                                
                                !just for debbuging
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(PropertyNameX)//"_Grid",     &
                                                     trim(PropertyNameX)//"_Grid", PropertyX%ID%Units,          &
                                                     Array2D      = PropertyX%FieldU,            &
                                                     OutputNumber = OutPutNumber,                &
                                                     STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                                       &
                                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR06b'    
                                
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyNameY,     &
                                                     PropertyNameY, PropertyX%ID%Units,          &
                                                     Array2D      = PropertyX%FieldY,            &
                                                     OutputNumber = OutPutNumber,                &
                                                     STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                                       &
                                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR07a'     
                                
                                !just for debbuging
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(PropertyNameY)//"_Grid",     &
                                                     trim(PropertyNameY)//"_Grid", PropertyX%ID%Units,          &
                                                     Array2D      = PropertyX%FieldV,            &
                                                     OutputNumber = OutPutNumber,                &
                                                     STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                                       &
                                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR07b'                                 
                                
                            else
                            
                                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//PropertyX%ID%Name, &
                                                     PropertyX%ID%Name, PropertyX%ID%Units,      &
                                                     Array2D      = PropertyX%Field,             &
                                                     OutputNumber = OutPutNumber,                &
                                                     STAT         = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                                       &
                                    stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR08a'
                                
                            endif
                           
                         endif
                           
                                                 
                        !Writes everything to disk
                        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'OutPut_Results_HDF - ModuleInterfaceWaterAir - ERR09'

                    endif

                PropertyX => PropertyX%Next

            enddo PropX

            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1


        endif  TOut    

        nullify(PropertyX)


    end subroutine OutPut_Results_HDF


    !--------------------------------------------------------------------------
    
    
    subroutine SetSubModulesModifier

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property), pointer               :: PropertyX
        real, pointer, dimension(:,:)           :: AtmPressure, Precipitation

        !Begin-----------------------------------------------------------------

        
        if(Me%ExtOptions%HeatFluxYes)then

            call Search_Property(PropertyX, NonSolarFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR10'
            
            endif


            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR20'
            


            endif

        endif

        if(Me%ExtOptions%OxygenFluxYes)then

            call Search_Property(PropertyX, OxygenFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR30'
            endif

        endif

        if(Me%ExtOptions%CarbonDioxideFluxYes)then

            call Search_Property(PropertyX, CarbonDioxideFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR40'
            endif

        endif


        if(Me%ExtOptions%AmmoniaFluxYes)then   !LLP

            call Search_Property(PropertyX, AmmoniaFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR45'
            endif

        endif
        
        if(Me%ExtOptions%NitrateFluxYes)then   !LLP

            call Search_Property(PropertyX, NitrateFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR45a'
            endif

        endif

        if(Me%ExtOptions%WQMYes)then

            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR50'
           
            endif

        endif

        if(Me%ExtOptions%T90VariableYes)then

            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceFlux(Me%ObjWaterProperties,                  &
                                    PropertyX%ID%IDNumber,                  &
                                    PropertyX%Field,                        &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR60'
           
            endif

        endif

        if(Me%ExtOptions%SurfaceWaterFluxYes)then

            call Search_Property(PropertyX, SurfaceWaterFlux_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                call SetSurfaceWaterFlux(Me%ObjHydrodynamic,                  &
                                         PropertyX%Field,                     &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR70'
                    
                if (Me%ExtOptions%Precipitation) then
                    
                    call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,         &
                                               Scalar       = Precipitation,            &
                                               ID           = Precipitation_,           &
                                               STAT         = STAT_CALL)
                    if (STAT_CALL .eq. SUCCESS_) then

                        call SetSurfaceFlux(Me%ObjWaterProperties,                      &
                                            Precipitation_,                             &
                                            Precipitation,                              &
                                            STAT = STAT_CALL)

                        if (STAT_CALL /= SUCCESS_) &
                            stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR75'                
                    else

                        stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR77'                                            
                    
                    endif
                    
                   call UnGetAtmosphere(Me%ObjAtmosphere, Precipitation, &
                                        STAT         = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) &
                        stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR79'
                    
                
                endif                    
            endif

        endif

        if(Me%ExtOptions%HydrodynamicWindYes)then

            call Search_Property(PropertyX, WindStress_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then

                !call Search_Property(PropertyY, WindStressY_, STAT = STAT_CALL) 

                if (STAT_CALL == SUCCESS_)then

                    call SetWindStress (Me%ObjHydrodynamic, PropertyX%FieldU, &
                                        PropertyX%FieldV, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR80'

                endif

            endif

        endif


        !Use the Atmospheric Pressure for Sea level pressure reference
        if(Me%ExtOptions%HydrodynamicAtmPressureYes)then

            call SetAtmosphericPressure(Me%ObjHydrodynamic,                         &
                                        Me%LocalAtm%AtmosphericPressure%Field,      &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR90'
        
        !Or Use the Mean Sea Level Pressure for Sea level pressure reference
        elseif(Me%ExtOptions%HydrodynamicMslpYes)then

            call SetAtmosphericPressure(Me%ObjHydrodynamic,                         &
                                        Me%LocalAtm%Mslp%Field,      &
                                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR90'
        
        endif
        
        if (Me%IntOptions%WindShearVelocity)then
            
            call SetTurbGOTMWindShearVelocity(TurbGOTMID        = Me%ObjTurbGOTM,           &
                                              WindShearVelocity = Me%WindShearVelocity,     &
                                              STAT              = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                      &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR100'
            
        endif

#ifndef _WAVES_
        if(Me%ExtOptions%WavesWindYes)then

            call SetWavesWind(WavesID       = Me%ObjWaves,                                      &
                              WindU         = Me%LocalAtm%WindVelocityU%Field,                  &
                              WindV         = Me%LocalAtm%WindVelocityV%Field,                  &
                              STAT          = STAT_CALL)                                    
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR110'

        endif
#endif

i22:    if (Me%ObjLagrangian /= 0) then
#ifndef _LAGRANGIAN_
#ifdef  _LAGRANGIAN_GLOBAL_        
            !Checks Lagrangian 
            call GetLagrangianAirOptionsGlobal(LagrangianID  = Me%ObjLagrangian,                &
                                         Oil           = Me%ExtOptions%OilYes,                  &
                                         HNS           = Me%ExtOptions%HNSYes,                  &
                                         Wind          = Me%ExtOptions%LagrangianWindYes,       &
                                         WaterQuality  = Me%ExtOptions%LagrangianWQMYes,        &
                                         T90Variable   = Me%ExtOptions%LagrangianT90Yes,        &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR120'
#endif
#endif

#ifndef _LAGRANGIAN_
        if(Me%ExtOptions%LagrangianWindYes)then     
#ifdef  _LAGRANGIAN_GLOBAL_        
            !Checks Lagrangian 

            if (Me%ExtOptions%HNSYes) then

                call GetAtmosphereProperty(AtmosphereID = Me%ObjAtmosphere,                 &
                                           Scalar       = Me%ExtAtm%AirTemperature%Field,   &
                                           ID           = AirTemperature_,                  &
                                           STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR121'

                call SetLagrangianAirTemperature(LagrangianID = Me%ObjLagrangian,                   &
                                       ModelName              = Me%ModelName,                       &
                                       AirTemperature         = Me%ExtAtm%AirTemperature%Field,     &
                                       STAT                   = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                          &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR122'
            
                call UnGetAtmosphere(Me%ObjAtmosphere, Me%ExtAtm%AirTemperature%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR123'
            endif
            
            call SetLagrangianWindGlobal(LagrangianID = Me%ObjLagrangian,                       &
                                   ModelName    = Me%ModelName,                                 &
                                   WindX        = Me%LocalAtm%WindVelocityU%Field,              &
                                   WindY        = Me%LocalAtm%WindVelocityV%Field,              &
                                   STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR130'
#else
            call SetLagrangianWind(LagrangianID = Me%ObjLagrangian,                       &
                                   WindX        = Me%LocalAtm%WindVelocityU%Field,              &
                                   WindY        = Me%LocalAtm%WindVelocityV%Field,              &
                                   STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR140'
#endif                
        endif

        if(Me%ExtOptions%LagrangianWQMYes  .or. Me%ExtOptions%LagrangianT90Yes)then
            call Search_Property(PropertyX, SurfaceRadiation_, STAT = STAT_CALL) 
            if (STAT_CALL == SUCCESS_)then        
#ifdef  _LAGRANGIAN_GLOBAL_ 
                call SetLagSolarRadiationGlobal( LagrangianID   = Me%ObjLagrangian,             &
                                                 ModelName      = Me%ModelName,                 &
                                                 SolarRadiation = PropertyX%Field,              &
                                                 STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                      &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR150'
#else
                call SetLagrangianSolarRadiation(LagrangianID   = Me%ObjLagrangian,             &
                                                 SolarRadiation = PropertyX%Field,              &
                                                 STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                      &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR160'
#endif
            endif
        
        endif

        if(Me%ExtOptions%OilYes .or. Me%ExtOptions%HNSYes)then

            !In this situation, if both Atmospheric Pressure and Mean Sea
            !Level Pressure are present, then MSLP will take precedence over
            !Atmospheric Pressure.
            
            nullify(AtmPressure)
                        
            !Checks if Mslp Field is present
            if(Me%LocalAtm%Mslp%Yes) then

                AtmPressure  => Me%LocalAtm%Mslp%Field

            !Checks if Atmospheric Pressure Field is present
            elseif(Me%LocalAtm%AtmosphericPressure%Yes) then

                AtmPressure  => Me%LocalAtm%AtmosphericPressure%Field

            endif

            if(Me%LocalAtm%Mslp%Yes .or. Me%LocalAtm%AtmosphericPressure%Yes) then
               
#ifdef  _LAGRANGIAN_GLOBAL_ 
                call SetLagrangianAtmPressureGlobal(LagrangianID = Me%ObjLagrangian,    &
                                              ModelName    = Me%ModelName,              &
                                              AtmPressure  = AtmPressure,               &
                                              STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR170'
#else
                call SetLagrangianAtmPressure(LagrangianID = Me%ObjLagrangian,          &
                                              AtmPressure  = AtmPressure,               &
                                              STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'SetSubModulesModifier - ModuleInterfaceWaterAir - ERR180'
#endif                        
                nullify(AtmPressure)
            endif

        endif
#endif
        
        endif i22
        


    end subroutine SetSubModulesModifier

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillInterfaceWaterAir(ObjInterfaceWaterAirID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjInterfaceWaterAirID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, STAT_CALL, nUsers              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_           
        type (T_Property),    pointer       :: PropertyX

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjInterfaceWaterAirID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mINTERFACEWATERAIR_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%OutPut%Yes) then

                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR10'
                endif

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR20'
                
                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR30'
                
                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjGridData)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR40'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR50'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR60'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjMap)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR70'

                nUsers = DeassociateInstance(mHYDRODYNAMIC_,    Me%ObjHydrodynamic)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR80'

                nUsers = DeassociateInstance(mATMOSPHERE_,      Me%ObjAtmosphere)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR90'
                
                nUsers = DeassociateInstance(mWATERPROPERTIES_, Me%ObjWaterProperties)
                if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR100'

#ifndef _LAGRANGIAN_
                if(Me%ObjLagrangian /= 0)then
                    nUsers = DeassociateInstance(mLAGRANGIAN_, Me%ObjLagrangian)
                    if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR110'
                endif
#endif
                if(Me%ObjTurbGOTM /= 0)then
                    nUsers = DeassociateInstance(mTURBGOTM_, Me%ObjTurbGOTM)
                    if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR120'
                endif

#ifndef _WAVES_
                if(Me%ObjWaves /= 0)then
                    nUsers = DeassociateInstance(mWAVES_, Me%ObjWaves)
                    if (nUsers == 0) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR130'
                endif
#endif

                if(Me%Coupled%TimeSerie%Yes)then

                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR140'
                endif


                if(associated(Me%LocalAtm%WindVelocityU%Field      )) &
                   deallocate(Me%LocalAtm%WindVelocityU%Field      )
                if(associated(Me%LocalAtm%WindVelocityV%Field      )) &
                   deallocate(Me%LocalAtm%WindVelocityV%Field      )
                if(associated(Me%LocalAtm%WindDirection%Field      )) &
                   deallocate(Me%LocalAtm%WindDirection%Field      )
                if(associated(Me%LocalAtm%AtmosphericPressure%Field)) &
                   deallocate(Me%LocalAtm%AtmosphericPressure%Field)                
                if(associated(Me%LocalAtm%Mslp%Field)) &    
                   deallocate(Me%LocalAtm%Mslp%Field)   
                if(associated(Me%SurfaceTemperature)) &   
                   deallocate(Me%SurfaceTemperature)
                
                if(associated(Me%SurfaceAlbedo))then
                    deallocate (Me%SurfaceAlbedo, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR145'
                    nullify    (Me%SurfaceAlbedo)
                endif
                
                PropertyX => Me%FirstProperty

                do while(associated(PropertyX))
                    
                    
!~                     if (Check_Vectorial_Property(PropertyX%ID%IDNumber)) then
                    if (PropertyX%ID%IsVectorial) then
                        
                        deallocate(PropertyX%FieldU,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR146'    
                        
                        deallocate(PropertyX%FieldV,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR147'    
                        
                        deallocate(PropertyX%FieldX,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR148'    
                        
                        deallocate(PropertyX%FieldY,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_) stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR149'
                        
                    else
                        
                        deallocate(PropertyX%Field,   STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR150'
                        nullify(PropertyX%Field)                        
                    endif


                    if(PropertyX%ID%SolutionFromFile)then

                        call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)&
                            stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR160'
                    endif


                    PropertyX => PropertyX%Next

                end do
!**************************************CoareDeallocate****************************************
                  if(Me%COARE) then
                      
                  deallocate (Me%Fxp, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR180'
                    nullify(Me%Fxp)
                    
                  deallocate (Me%WarmLayerThickness, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR190'
                    nullify(Me%WarmLayerThickness)

                  deallocate (Me%Tau_ac, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR200'
                    nullify(Me%Tau_ac)
      
                  deallocate (Me%Tau, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR210'
                    nullify(Me%Tau)
        
                  deallocate (Me%AccumulatedEnergy, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR220'
                    nullify(Me%AccumulatedEnergy)
        
                  deallocate (Me%WarmLayerTempDiff, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR230'
                    nullify(Me%WarmLayerTempDiff)
    
                  deallocate (Me%Al, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR240'
                    nullify(Me%Al)

                  deallocate (Me%RainFlux, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR250'
                    nullify(Me%RainFlux)
      
                  deallocate (Me%LastLatentHeat, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR260'
                    nullify(Me%LastLatentHeat)
      
                  deallocate (Me%LastSensibleHeat, STAT=STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR270'
                    nullify(Me%LastSensibleHeat)
                    
                  endif
!***********************************************EndCoareDeallocate*****************************        
                if (Me%Rugosity%ON) then

                    deallocate(Me%Rugosity%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR280'
                    nullify(Me%Rugosity%Field)

                    if(Me%Rugosity%ID%SolutionFromFile) then

                        call KillFillMatrix(Me%Rugosity%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)&
                            stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR290'
                    endif

                endif

! Modified by Matthias DELPEY - 18/10/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 16/12/2011
                if (Me%WaveFluxTKE%ON) then

                    deallocate(Me%WaveFluxTKE%Field,   STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)&
                        stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR170a'
                    nullify(Me%WaveFluxTKE%Field)

                    if(Me%WaveFluxTKE%ID%SolutionFromFile) then

                        call KillFillMatrix(Me%WaveFluxTKE%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)&
                            stop 'KillInterfaceWaterAir - ModuleInterfaceWaterAir - ERR180a'
                    endif

                endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                if(Me%Coupled%BoxTimeSerie%Yes)then
                
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillInterfaceWaterAir - ModuleInterfaceSedimentWater - ERR300'

                    deallocate(Me%Scalar2D)
                    nullify   (Me%Scalar2D)
                endif
                
                if (Me%IntOptions%WindShearVelAllocate) then
                    deallocate(Me%WindShearVelocity)
                    nullify   (Me%WindShearVelocity)
                endif

                
                !Deallocates Instance
                call DeallocateInstance

                ObjInterfaceWaterAirID = 0
                STAT_      = SUCCESS_

            endif
        else 
            STAT_ = ready_
        endif cd1

        if (present(STAT)) STAT = STAT_
           

    end subroutine KillInterfaceWaterAir
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_InterfaceWaterAir), pointer          :: AuxObjInterfaceWaterAir
        type (T_InterfaceWaterAir), pointer          :: PreviousObjWaterSedInterface

        !Updates pointers
        if (Me%InstanceID == FirstObjInterfaceWaterAir%InstanceID) then
            FirstObjInterfaceWaterAir           => FirstObjInterfaceWaterAir%Next
        else
            PreviousObjWaterSedInterface        => FirstObjInterfaceWaterAir
            AuxObjInterfaceWaterAir             => FirstObjInterfaceWaterAir%Next
            do while (AuxObjInterfaceWaterAir%InstanceID /= Me%InstanceID)
                PreviousObjWaterSedInterface    => AuxObjInterfaceWaterAir
                AuxObjInterfaceWaterAir         => AuxObjInterfaceWaterAir%Next
            enddo

            !Now update linked list
            PreviousObjWaterSedInterface%Next   => AuxObjInterfaceWaterAir%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjInterfaceWaterAir_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceWaterAir_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjInterfaceWaterAir_ID > 0) then
            call LocateObjInterfaceWaterAir (ObjInterfaceWaterAir_ID)
            ready_ = VerifyReadLock (mINTERFACEWATERAIR_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        endif cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjInterfaceWaterAir (ObjInterfaceWaterAirID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjInterfaceWaterAirID

        !Local-----------------------------------------------------------------

        Me => FirstObjInterfaceWaterAir
        do while (associated (Me))
            if (Me%InstanceID == ObjInterfaceWaterAirID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleInterfaceWaterAir - LocateObjInterfaceWaterAir - ERR01'

    end subroutine LocateObjInterfaceWaterAir

    !--------------------------------------------------------------------------

end module ModuleInterfaceWaterAir

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------






