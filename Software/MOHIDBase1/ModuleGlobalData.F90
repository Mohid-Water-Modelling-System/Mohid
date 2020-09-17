!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Global Data
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which store global parameter
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

Module ModuleGlobalData
    use, intrinsic      :: Iso_C_Binding

    implicit none
    

    public

    !Subroutines---------------------------------------------------------------
    public  ::  CheckPropertyName
    public  ::  GetPropertyName
    public  ::  GetPropertyIDNumber
    private ::      ConstructPropList
    private ::          AddPropList
!~     public  ::  CheckDynamicPropertyName
    public  ::  RegisterDynamicProperty
    public  ::  GetDynamicPropertyIDNumber

    public  ::  Check_Particulate_Property
    public  ::  TranslateTypeZUV
    
    public  ::  UnitsManager
    public  ::  SetError
    public  ::  WriteErrorMessage
    public  ::  LogKeyWord
    public  ::  StartupMohid
    public  ::  ShutdownMohid

    public  ::  GetUsersNumber
    
    !Modules management-----------------------------------------------------------------
    public  ::  RegisterNewInstance
    public  ::  AssociateInstance
    public  ::  DeassociateInstance
    public  ::  Read_Lock
    public  ::  Read_Unlock
    public  ::  VerifyReadLock


    !Interfaces
    private :: SetErrorTime
    private :: SetErrorMessage
    interface  SetError
        module procedure SetErrorTime
        module procedure SetErrorMessage
    end interface SetError
    
    !Parameter-----------------------------------------------------------------
    integer, parameter  :: MaxModules           =  98

#ifdef _INCREASE_MAXINSTANCES
    integer, parameter  :: MaxInstances         = 2000
#else
    integer, parameter  :: MaxInstances         = 500
#endif

    integer, parameter  :: MaxErrorMessages     = 20
    integer             :: NumberOfErrorMessages=0

    integer, parameter  :: StringLength         = 128
    integer, parameter  :: PathLength           = 256
    !max http string length 2018
    integer, parameter  :: LinkLength           = 2000

    !Search in Block, file, etc...
    integer, parameter :: FromFile_                = 1
    integer, parameter :: FromBlock_               = 2
    integer, parameter :: FromBlockInBlock_        = 3
    integer, parameter :: FromBlockInBlockInBlock_ = 4

    integer, parameter :: FromFile                = 1
    integer, parameter :: FromBlock               = 2
    integer, parameter :: FromBlockInBlock        = 3
    integer, parameter :: FromBlockInBlockInBlock = 4

    !Matrix Types (Centered in Z, U, V, W)
    integer, parameter :: TypeZ_                = 1
    integer, parameter :: TypeU_                = 2
    integer, parameter :: TypeV_                = 3
    integer, parameter :: TypeW_                = 4    

#if defined(_SHORT_LINE_LENGTH)
    integer, parameter  :: line_length          = 64
#elif defined(_LONG_LINE_LENGTH)
    integer, parameter  :: line_length          = 1024
#elif defined(_BIG_LINE_LENGTH)
    integer, parameter  :: line_length          = 6144    
#elif defined(_EXTRA_LONG_LINE_LENGTH)
    integer, parameter  :: line_length          = 131072
#elif defined(_EXTRA_SHORT_LINE_LENGTH)
    integer, parameter  :: line_length          = 32
#else
    integer, parameter  :: line_length          = 256
#endif
                                                
    real   , parameter  :: FillValueReal        = -9.9e15
    integer, parameter  :: FillValueInt         = -9999999
    real   , parameter  :: HalfFillValueReal    = -9.9e15/2.0
                                                
    real,    parameter  :: AllmostZeroFraction  = 1.e-5
    real,    parameter  :: AllmostZero          = 1.e-15
    real,    parameter  :: AlmostZero           = 1.e-15

    integer,            parameter :: null_int   = -999999
    real,               parameter :: null_real  = -9.9E15
    character(LEN = 5), parameter :: null_str   = '*****'

    !Time
    real,               parameter :: CyclicTime = -9999

    !Space units
    integer,            parameter :: meters_      = 1
    integer,            parameter :: centimeters_ = 2
    
    !File formats
    integer,            parameter :: HDF5_        = 1
    integer,            parameter :: NetCDF_      = 2
    integer,            parameter :: Binary_      = 3  
    
    !Vertical referentials
    integer,            parameter :: Hydrographic_  = 1
    integer,            parameter :: Topographic_   = 2
    
    !Join matrixes from sub-domains to global
    integer, parameter :: Type1DI = 1, Type1DJ = 2
       
    !character
    character(LEN = 1), parameter :: space      = char(32)   !" "
    character(LEN = 1), parameter :: dot        = char(46)   !"."
    character(LEN = 1), parameter :: delimiter  = char(58)   !":"
    character(LEN = 1), parameter :: semicolumn = char(59)   !";"
    character(LEN = 1), parameter :: tab        = char(9)    !" "
    character(LEN = 1), parameter :: backslash  = char(92)   !"\"
    character(LEN = 1), parameter :: exclamation= char(33)   !"!"

    logical, parameter :: OFF    = .FALSE.
    logical, parameter :: ON     = .TRUE.
    logical, parameter :: IDLE   = .FALSE.
    logical, parameter :: ACTIVE = .TRUE.

    !Errors
    integer, parameter :: UNKNOWN_               =-99
    integer, parameter :: SUCCESS_               = 0
    integer, parameter :: OFF_ERR_               = 1 
    integer, parameter :: ON_ERR_                = 2
    integer, parameter :: FILE_NOT_FOUND_ERR_    = 3 
    integer, parameter :: SIZE_ERR_              = 4
    integer, parameter :: BLOCK_END_ERR_         = 5
    integer, parameter :: KEYWORD_NOT_FOUND_ERR_ = 6
    integer, parameter :: READ_LOCK_ERR_         = 7
    integer, parameter :: IDLE_ERR_              = 9
    integer, parameter :: BLOCK_LOCK_ERR_        = 10
    integer, parameter :: BLOCK_UNLOCK_ERR_      = 11
    integer, parameter :: NBUSERS_ERR_           = 12
    integer, parameter :: CLIENT_NB_ERR_         = 13
    integer, parameter :: UNIT_ERR_              = 14
    integer, parameter :: INPUT_ERR_             = 15
    integer, parameter :: NO_ID_ERR_             = 16
    integer, parameter :: ID_ERR_                = 17
    integer, parameter :: NOT_FOUND_ERR_         = 18
    integer, parameter :: TIME_ERR_              = 19
    integer, parameter :: NOT_ASSOCIATE_         = 20
    integer, parameter :: FILE_EXISTS_ERR_       = 21
    integer, parameter :: OUT_OF_BOUNDS_ERR_     = 22 !Add to use with ModulePhreeqC
    
    !MOHID LAND AND MOHID WATER and MOHID RIVER
    integer, parameter :: MOHIDLAND_      = 1
    integer, parameter :: MOHIDWATER_     = 2
    integer, parameter :: MOHIDRIVER_     = 3
    

    !Types of coordinates
    integer, parameter :: GEOG_             = 1         !Coordenadas Geograficas  
    integer, parameter :: UTM_              = 2         !Coordenadas (UTM)
    integer, parameter :: MIL_PORT_         = 3         !Coordenadas Militares Portuguesas     
    integer, parameter :: SIMPLE_GEOG_      = 4         !Coordenadas Geograficas Simplificadas
                                                        !SIMPLE_GEOG_ - considera a terra como uma esfera
    integer, parameter :: GRID_COORD_       = 5         !Coordenadas relativas à origem da malha
    integer, parameter :: CIRCULAR_         = 6         !Coordenadas circulares (XX- raio, YY - anglo em graus)
    integer, parameter :: NLRD_             = 7         !Coordenadas Netherlands RD
    integer, parameter :: LAMB_CONF_CONIC_  = 8
    integer, parameter :: POLAR_STEREO_     = 9
    integer, parameter :: SPHERE_MERCATOR_  = 10
    integer, parameter :: PAULO_PROJECTION_ = 11
        

    integer, parameter :: PORTUGUESE_UTM_ZONE_ = 29

    !Types of referential
    integer, parameter :: GridCoord_ = 1
    integer, parameter :: Cartesian_ = 2
    integer, parameter :: AlongGrid_ = 3

    !Types of direction referential
    integer, parameter :: NauticalWind_    = 1
    integer, parameter :: NauticalCurrent_ = 2    
    integer, parameter :: CartesianDir_    = 3


    !Types of grid borders 
    integer, parameter :: ComplexPolygon_  = 1
    integer, parameter :: RotatedRectang_  = 2
    integer, parameter :: Rectang_         = 3


    !Mohid Land ground water flow type
    integer, parameter :: GWFlowToChanByCell_               = 1
    integer, parameter :: GWFlowToChanByLayer_              = 2
    
    !Angles Referential
    integer, parameter :: NauticalReferential_              = 1   !0º is from N, 90º is from E
    integer, parameter :: CurrentsReferential_              = 2   !0º is to N, 90º is to E

    !Water Properties 
    integer, parameter :: Density_                          =  0        
    integer, parameter :: Temperature_                      =  1        
    integer, parameter :: Salinity_                         =  2 
    integer, parameter :: Phytoplankton_                    =  3 
    integer, parameter :: Zooplankton_                      =  4 
    integer, parameter :: DOPRefractory_                    =  6
    integer, parameter :: DOPNon_Refractory_                =  7
    integer, parameter :: DONRefractory_                    =  8 
    integer, parameter :: DONNon_Refractory_                =  9  
    integer, parameter :: Inorganic_Phosphorus_             = 10

!____used @ modulelife______________________________________________
    integer, parameter :: POC_                              = 11   
    integer, parameter :: POP_                              = 12   
    integer, parameter :: PON_                              = 13   
    integer, parameter :: PONRefractory_                    = 14       
    integer, parameter :: DOC_                              = 15   
    integer, parameter :: DOP_                              = 16   
    integer, parameter :: DON_                              = 17   
    integer, parameter :: Ammonia_                          = 20 
    integer, parameter :: Nitrate_                          = 21         
    integer, parameter :: Silicate_                         = 22 
    integer, parameter :: BioSilica_                        = 23
    integer, parameter :: CarbonDioxide_                    = 24
    integer, parameter :: Oxygen_                           = 25 
    integer, parameter :: DissolO2PercentSat_               = 26
    integer, parameter :: CO2PartialPressure_               = 27

    !Bivalve species
    integer, parameter :: Bivalve1_                         = 28
    integer, parameter :: Bivalve2_                         = 29
    integer, parameter :: Bivalve3_                         = 30
    integer, parameter :: Bivalve4_                         = 31

    !Bivalve predators
    integer, parameter :: Shrimp_                           = 32
    integer, parameter :: Crab_                             = 33
    integer, parameter :: OysterCatcher_                    = 34
    integer, parameter :: EiderDuck_                        = 35
    integer, parameter :: HerringGull_                      = 36

    integer, parameter :: PhytoChla_                        = 49
    integer, parameter :: OilThickness_                     = 50
    integer, parameter :: Nitrite_                          = 51
    integer, parameter :: BOD_                              = 52
    integer, parameter :: Cohesive_Sediment_                = 53
    integer, parameter :: Oil_                              = 54
    integer, parameter :: Ciliate_                          = 55
    integer, parameter :: Bacteria_                         = 56
    integer, parameter :: ParticulateArsenic_               = 57
    integer, parameter :: DissolvedArsenic_                 = 58
    integer, parameter :: Larvae_                           = 59
    integer, parameter :: Age_                              = 60
    integer, parameter :: Fecal_Coliforms_                  = 61
    integer, parameter :: Fish_                             = 62
    integer, parameter :: FishFood_                         = 63
    integer, parameter :: MacroAlgae_                       = 64
    integer, parameter :: MicroPhytoBenthos_                = 65
    integer, parameter :: AdsorbedAmmonia_                  = 66
    integer, parameter :: RefreactaryOrganicN_              = 67
    integer, parameter :: Ngas_                             = 68
    integer, parameter :: HeterotrophicN_                   = 69
    integer, parameter :: AnaerobicN_                       = 70
    integer, parameter :: AutotrophicN_                     = 71
    integer, parameter :: AnaerobicC_                       = 72
    integer, parameter :: AutotrophicC_                     = 73
    integer, parameter :: HeterotrophicC_                   = 74
    integer, parameter :: LabileOrganicC_                   = 75
    integer, parameter :: RefreactaryOrganicC_              = 76
    integer, parameter :: HeterotrophicP_                   = 7000
    integer, parameter :: AutotrophicP_                     = 7001
    integer, parameter :: AnaerobicP_                       = 7002
    integer, parameter :: LabileOrganicP_                   = 7003
    integer, parameter :: RefreactaryOrganicP_              = 7004
    integer, parameter :: AdsorbedInorganicP_               = 7005
    integer, parameter :: SolubilizingC_                    = 7006
    integer, parameter :: SolubilizingN_                    = 7007
    integer, parameter :: SolubilizingP_                    = 7008
    integer, parameter :: Urea_                             = 7009
    integer, parameter :: AmmoniaGas_                       = 7010
    integer, parameter :: Methane_                          = 7011
    integer, parameter :: SoilDryDensity_                   = 7012
    integer, parameter :: IonicStrength_                    = 7013
    integer, parameter :: PhosphorusAdsortionIndex_         = 7014
    integer, parameter :: AutotrophicPop_                   = 7015
    integer, parameter :: HeterotrophicPop_                 = 7016
    integer, parameter :: AnaerobicPop_                     = 7017
    integer, parameter :: SolPop_                           = 7018
    
    ! Generic metal
    integer, parameter :: ParticulateMetal_                 = 8000
    integer, parameter :: DissolvedMetal_                   = 8001
    
    ! Metals
    integer, parameter :: ParticulateCopper_                = 8002  ! Cu
    integer, parameter :: DissolvedCopper_                  = 8003  
    integer, parameter :: ParticulateZinc_                  = 8004  ! Zn
    integer, parameter :: DissolvedZinc_                    = 8005
    integer, parameter :: ParticulateLead_                  = 8006  ! Pb
    integer, parameter :: DissolvedLead_                    = 8007
    integer, parameter :: ParticulateCadmium_               = 8008  ! Cd
    integer, parameter :: DissolvedCadmium_                 = 8009
    integer, parameter :: ParticulateMercury_               = 8010  !  Hg
    integer, parameter :: DissolvedMercury_                 = 8011
    
!____POM pools (for aquaculture cages)__________________________________________    
    
    integer, parameter :: PON1_                             = 7020
    integer, parameter :: PON2_                             = 7021       
    integer, parameter :: PON3_                             = 7022
    integer, parameter :: PON4_                             = 7023
    integer, parameter :: PON5_                             = 7024

    integer, parameter :: POP1_                             = 7025
    integer, parameter :: POP2_                             = 7026       
    integer, parameter :: POP3_                             = 7027
    integer, parameter :: POP4_                             = 7028
    integer, parameter :: POP5_                             = 7029    
     
        
    
    integer, parameter :: GenericProperty_                  = 77
    integer, parameter :: T90_                              = 78
    integer, parameter :: E_Coli_                           = 79
    integer, parameter :: T90_E_Coli_                       = 799999


    !Production parameters  
    integer, parameter :: GrossProd_                        = 80
    integer, parameter :: NutrientLim_                      = 81
    integer, parameter :: LightLim_                         = 82
    integer, parameter :: TemperatureLim_                   = 83
    integer, parameter :: SalinityLim_                      = 84
    integer, parameter :: NetProd_                          = 85
    integer, parameter :: DiaGrossProd_                     = 800 
    integer, parameter :: DiaNutrientLim_                   = 801 
    integer, parameter :: DiaLightLim_                      = 802 
    integer, parameter :: DiaTemperatureLim_                = 803 
    integer, parameter :: Diatoms_                          = 804 
    integer, parameter :: NLim_                             = 805
    integer, parameter :: PLim_                             = 806
    integer, parameter :: DiaNLim_                          = 807
    integer, parameter :: DiaPLim_                          = 808
    integer, parameter :: DiaSiLim_                         = 809
    integer, parameter :: Excretion_                        = 810
    integer, parameter :: Respiration_                      = 811
    integer, parameter :: NaturalMort_                      = 812
    integer, parameter :: Grazing_                          = 813
    integer, parameter :: MACondition_                      = 814
    integer, parameter :: CarrCapLim_                       = 815
    
    !Drifting macroalgae
    integer, parameter :: DriftingMacroAlgae_               = 850

    ! BenthicEcology -isabella
    integer, parameter :: DepositFeedersC_                  = 851
    integer, parameter :: DepositFeedersN_                  = 852
    integer, parameter :: DepositFeedersP_                  = 853
    integer, parameter :: SuspensionFeedersN_               = 860 
    integer, parameter :: SuspensionFeedersC_               = 861
    integer, parameter :: SuspensionFeedersP_               = 862
    integer, parameter :: MicroPhytoBenthosN_               = 864 
    integer, parameter :: MicroPhytoBenthosC_               = 865
    integer, parameter :: MicroPhytoBenthosP_               = 866
    integer, parameter :: SeagrassesLeaves_                 = 867             
    integer, parameter :: SeagrassesRoots_                  = 868             
    integer, parameter :: SeagrassesN_                      = 869             
    integer, parameter :: SeagrassesP_                      = 870             
    integer, parameter :: LeavesUptakeN_                    = 871  
    integer, parameter :: LeavesUptakeP_                    = 872
    integer, parameter :: LeavesLightFactor_                = 873
    integer, parameter :: RootsUptakeN_                     = 874
    integer, parameter :: NintFactor_                       = 875
    integer, parameter :: PintFactor_                       = 876
    integer, parameter :: NintFactorR_                      = 877
    integer, parameter :: RootsMort_                        = 878  
    integer, parameter :: PintFactorR_                      = 879
    integer, parameter :: RootsUptakeP_                     = 880  
    
    !Hydrodynamic Properties
    integer, parameter :: WaterLevel_                       = 86
    integer, parameter :: VelocityU_                        = 87
    integer, parameter :: VelocityV_                        = 88
    integer, parameter :: WaterFluxX_                       = 89
    integer, parameter :: WaterFluxY_                       = 90
    integer, parameter :: WaterDepth_                       = 900
    integer, parameter :: Volume_                           = 901
    integer, parameter :: PercentageMaxVolume_              = 902
    integer, parameter :: VelocityModulus_                  = 903
    integer, parameter :: Courant_                          = 904
    integer, parameter :: FlowModulus_                      = 905
    integer, parameter :: ShearStress_                      = 906
    integer, parameter :: VelocityDirection_                = 907
    !integer, parameter :: WaterLevelMax_                    = 908
    !integer, parameter :: WaterLevelMin_                    = 909
    integer, parameter :: WaterColumn_                      = 910
    integer, parameter :: ZonalVelocity_                    = 920
    integer, parameter :: MeridionalVelocity_               = 930
    !1 - low tide, 2 - flood, 3 - high tide, 4 - ebb 
    integer, parameter :: TideState_                        = 940
    integer, parameter :: ShearStressX_                     = 950
    integer, parameter :: ShearStressY_                     = 960    
    
    !Assimilation Properties        guillaume nogueira
    integer, parameter :: AltimLevelAnalyzed_               = 4000
    integer, parameter :: AltimTemperatureAnalyzed_         = 4001
    integer, parameter :: AltimSalinityAnalyzed_            = 4002
    integer, parameter :: AltimLevelToAssimilate_           = 4003
    integer, parameter :: VarianceFieldToAssimilate_        = 4004
                                                            
    !Explicit forces                                        
    integer, parameter :: CoriolisX_                        = 91
    integer, parameter :: BaroclinicForceX_                 = 92
    integer, parameter :: HorizontalTransportX_             = 93
                                                            
    integer, parameter :: CoriolisY_                        = 94
    integer, parameter :: BaroclinicForceY_                 = 95
    integer, parameter :: HorizontalTransportY_             = 96

    !MRV. Barotropic velocities for using in Flather condition
    integer, parameter :: BarotropicVelocityU_              = 97
    integer, parameter :: BarotropicVelocityV_              = 98

    !Hydrodynamic Properties
    integer, parameter :: VelocityW_                        = 99

    integer, parameter :: ShearVelocity_                    = 100

    integer, parameter :: ParticulateContaminant_           = 101
    integer, parameter :: DissolvedContaminant_             = 102

    integer, parameter :: Sediment                          = 103

    !Soil DissolvedSodium Adsorption Ratio parameters
    integer, parameter :: DissolvedSodium_                  = 104
    integer, parameter :: DissolvedCalcium_                 = 105
    integer, parameter :: ParticulateSodium_                = 106
    integer, parameter :: ParticulateCalcium_               = 107

    !Hydrodynamic Properties
    integer, parameter :: Vorticity_                        = 108
    integer, parameter :: BaroclinicKE_                     = 109
    integer, parameter :: PerturbationPE_                   = 110
    integer, parameter :: KineticEnergy_                    = 111

    integer, parameter :: BaroclinicVelocityU_              = 113
    integer, parameter :: BaroclinicVelocityV_              = 114
    
    integer, parameter :: ObstacleDragCoef_                 = 200
    
    integer, parameter :: ScraperVelU_                      = 210
    integer, parameter :: ScraperVelV_                      = 211
    integer, parameter :: ScraperVelW_                      = 212
    
    !Geometry properties
    integer, parameter :: VerticalZ_                        = 300
     

    !Interface AirWater Properties
    integer, parameter :: LatentHeat_                       = 500
    integer, parameter :: SensibleHeat_                     = 501     
    integer, parameter :: Evaporation_                      = 502
    integer, parameter :: NetLongWaveRadiation_             = 503
    integer, parameter :: OxygenFlux_                       = 504
    integer, parameter :: WindShearVelocity_                = 505
    integer, parameter :: SurfaceRadiation_                 = 506
    integer, parameter :: WindStressX_                      = 507
    integer, parameter :: WindStressY_                      = 508
    integer, parameter :: SurfaceWaterFlux_                 = 509
    integer, parameter :: NonSolarFlux_                     = 510
    integer, parameter :: TurbulentKineticEnergy_           = 511
    integer, parameter :: CarbonDioxideFlux_                = 512
    integer, parameter :: Albedo_                           = 513
    integer, parameter :: UpwardLongWaveRadiation_          = 514
    integer, parameter :: DownwardLongWaveRadiation_        = 515
    integer, parameter :: ShortWaveSolarRadiation_          = 516
    integer, parameter :: LongWaveSolarRadiation_           = 517
    integer, parameter :: ShortWaveSolarRadiationExtin_     = 518
    integer, parameter :: LongWaveSolarRadiationExtin_      = 519
    integer, parameter :: SpecificOxygenFlux_               = 520
    integer, parameter :: SpecificCarbonDioxideFlux_        = 521
    integer, parameter :: AmmoniaFlux_                      = 522
    integer, parameter :: NitrateFlux_                      = 523
    !Vectorial 
    integer, parameter :: WindStress_                       = 524

    !AtmosProperties
    integer, parameter :: WindVelocityX_                    = 600
    integer, parameter :: WindVelocityY_                    = 601
    integer, parameter :: SolarRadiation_                   = 602
    integer, parameter :: Precipitation_                    = 603
    integer, parameter :: AtmosphericPressure_              = 604
    integer, parameter :: AirTemperature_                   = 605
    integer, parameter :: RelativeHumidity_                 = 606
    integer, parameter :: WindModulos_                      = 607
    integer, parameter :: WindAngle_                        = 608
    integer, parameter :: CloudCover_                       = 609
    integer, parameter :: Irrigation_                       = 610
    integer, parameter :: SunHours_                         = 611
    integer, parameter :: ATMTransmitivity_                 = 612
    integer, parameter :: MeanSeaLevelPressure_             = 613
    integer, parameter :: WindModulus_                      = 614
    integer, parameter :: WindDirection_                    = 615
    integer, parameter :: SpecificHumidity_                 = 616
    integer, parameter :: WindModulusBeaufort_              = 619
    integer, parameter :: WindGust_                         = 624
    integer, parameter :: PBLHeight_                        = 625
    integer, parameter :: Reflectivity_                     = 626
    !vectorial
    integer, parameter :: WindVelocity_                     = 627
    !air quality 
    integer, parameter :: CO2AtmosphericPressure_           = 617
    integer, parameter :: O2AtmosphericPressure_            = 618
    integer, parameter :: HydrogenSulfide_                  = 620
    integer, parameter :: MethylMercaptan_                  = 621
    integer, parameter :: AtmospDeposOxidNO3_               = 622
    integer, parameter :: AtmospDeposReduNH4_               = 623
    integer, parameter :: Visibility_                       = 628
    integer, parameter :: Dust_                             = 629    

    
    !Basin Properties
    integer, parameter :: RefEvapotrans_                    = 708
    integer, parameter :: TotalPlantBiomass_                = 709
    integer, parameter :: TotalPlantNitrogen_               = 710
    integer, parameter :: TotalPlantPhosphorus_             = 711
    integer, parameter :: RootBiomass_                      = 712
    integer, parameter :: RootDepth_                        = 713
    integer, parameter :: LeafAreaIndex_                    = 714
    integer, parameter :: SpecificLeafStorage_              = 715
    integer, parameter :: EVTPCropCoefficient_              = 716
    integer, parameter :: CanopyHeight_                     = 718
    integer, parameter :: PotLeafAreaIndex_                 = 719
    integer, parameter :: BoundaryLeafAreaIndex_            = 720

!____used @ moduleCEQUALW2______________________________________________
    integer, parameter :: RPOM_                             = 2001
    integer, parameter :: LPOM_                             = 2002
    integer, parameter :: LDOM_                             = 2003
    integer, parameter :: RDOM_                             = 2004
    integer, parameter :: PSilica_                          = 2005
    integer, parameter :: DSilica_                          = 2006
    integer, parameter :: ICarbon_                          = 2007
    integer, parameter :: pH_                               = 2008
    integer, parameter :: HCO3_                             = 2009
    integer, parameter :: CO3_                              = 2010      
    integer, parameter :: Algae_1_                          = 2011
    integer, parameter :: Algae_2_                          = 2012
    integer, parameter :: Algae_3_                          = 2013
    integer, parameter :: Algae_4_                          = 2014
    integer, parameter :: Algae_5_                          = 2015
    integer, parameter :: Epiphyton_1_                      = 2016
    integer, parameter :: Epiphyton_2_                      = 2017
    integer, parameter :: Epiphyton_3_                      = 2018
    integer, parameter :: Epiphyton_4_                      = 2019
    integer, parameter :: Epiphyton_5_                      = 2020
    integer, parameter :: Alkalinity_                       = 2040
    integer, parameter :: Detritus_                         = 2041

    ! Rates for output -------------------------------------------

    !Algal Growth Limitations                     
    integer, parameter :: ANLim_                           = 2101
    integer, parameter :: APLim_                           = 2102
    integer, parameter :: ASLim_                           = 2103
    integer, parameter :: ALightLim_                       = 2104
    integer, parameter :: AOverallLim_                     = 2106
    
    integer, parameter :: AGR_                             = 2051
    integer, parameter :: AMR_                             = 2052
    integer, parameter :: AER_                             = 2053
    integer, parameter :: ARR_                             = 2054 
    
    
    !Epiphyte Growth Limitations                         
    integer, parameter :: ENLim_                           = 2107
    integer, parameter :: EPLim_                           = 2108
    integer, parameter :: ESLim_                           = 2109
    integer, parameter :: ELightLim_                       = 2110
    integer, parameter :: EOverallLim_                     = 2112
                                            
    
    !Decay Rates                                
    integer, parameter ::  NH4D_                           = 2113
    integer, parameter ::  NO3D_                           = 2114
    integer, parameter ::  LDOMD_                          = 2115
    integer, parameter ::  RDOMD_                          = 2116
    integer, parameter ::  LPOMD_                          = 2117
    integer, parameter ::  RPOMD_                          = 2118
    integer, parameter ::  LRDOMD_                         = 2119
    integer, parameter ::  LRPOMD_                         = 2120
    integer, parameter ::  CBODD_                          = 2121
    
    !phosphorus
    integer, parameter ::  PO4ER_                          = 2122
    integer, parameter ::  PO4EG_                          = 2123
    integer, parameter ::  PO4AR_                          = 2124
    integer, parameter ::  PO4AG_                          = 2125
    integer, parameter ::  PO4OM_                          = 2126
    integer, parameter ::  PO4BOD_                         = 2127
            
    ! ammonia
    
    integer, parameter ::  NH4ER_                          = 2128 
    integer, parameter ::  NH4EG_                          = 2129
    integer, parameter ::  NH4AR_                          = 2130
    integer, parameter ::  NH4AG_                          = 2131
    integer, parameter ::  NH4OM_                          = 2132
    integer, parameter ::  NH4BOD_                         = 2133
    
    !nitrate
    integer, parameter ::  NO3AG_                          = 2134
    integer, parameter ::  NO3EG_                          = 2135

    !dissolved silica
    integer, parameter ::  DSIAG_                          = 2136
    integer, parameter ::  DSIEG_                          = 2137
    integer, parameter ::  DSID_                           = 2138

    !particulated silica
    integer, parameter ::  PSIAM_                          = 2139
    integer, parameter ::  PSID_                           = 2140 

    !LabDom
    integer, parameter ::  LDOMAP_                         = 2141
    integer, parameter ::  LDOMEP_                         = 2142


    !labPom
    integer, parameter ::  LPOMAP_                         = 2143


    !oxygen
    integer, parameter ::  DOAP_                           = 2144
    integer, parameter ::  DOEP_                           = 2145
    integer, parameter ::  DOAR_                           = 2146
    integer, parameter ::  DOER_                           = 2147
    integer, parameter ::  DOOM_                           = 2148
    integer, parameter ::  DONIT_                          = 2149
    
    !inorganic Carbon
    integer, parameter ::  ICarbonAP_                      = 2150
    integer, parameter ::  ICarbonEP_                      = 2151
    integer, parameter ::  ICarbonBOD_                     = 2152

    !sand transport
    integer, parameter ::  Diameter_                       = 3001
    integer, parameter ::  Percentage_                     = 3002
    integer, parameter ::  D35_                            = 3003
    integer, parameter ::  D50_                            = 3004
    integer, parameter ::  D90_                            = 3005
    integer, parameter ::  BedRock_                        = 3006
    integer, parameter ::  SandTauCritic_                  = 3007
    integer, parameter ::  Sand_                           = 3008
    integer, parameter ::  MappDZ_                         = 3009    

    integer, parameter ::  TransportCapacity_              = 3101 
    integer, parameter ::  TransportCapacityX_             = 3102 
    integer, parameter ::  TransportCapacityY_             = 3103 
    integer, parameter ::  BottomEvolution_                = 3104 
    integer, parameter ::  Newbathymetry_                  = 3105 
    integer, parameter ::  bathymetry_                     = 3106     

    !wave dynamics
    integer, parameter ::  WaveStressX_                    = 3401
    integer, parameter ::  WaveStressY_                    = 3402
    integer, parameter ::  CurrentX_                       = 3403
    integer, parameter ::  CurrentY_                       = 3404
    integer, parameter ::  WaveX_                          = 3405
    integer, parameter ::  WaveY_                          = 3406
    !vectorial
    integer, parameter :: WaveStress_                      = 3407
    
! Modified by Matthias DELPEY - 26/06/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 21/07/2011 - 04/08/2011 - 05/09/2011 - 25/10/2011 - 14/12/2011 
!                             - 16/12/2011 - 02/03/2012
    integer, parameter ::  WaveDriftSpecU_                 = 3408
    integer, parameter ::  WaveDriftSpecV_                 = 3409
    integer, parameter ::  AtmToWaveMomentumU_             = 3410
    integer, parameter ::  AtmToWaveMomentumV_             = 3411
    integer, parameter ::  WaveToOceanMomentumU_           = 3412
    integer, parameter ::  WaveToOceanMomentumV_           = 3413
    integer, parameter ::  StokesDriftU_                   = 3414
    integer, parameter ::  StokesDriftV_                   = 3415
    integer, parameter ::  StokesDriftModulus_             = 3416
    integer, parameter ::  StokesDriftW_                   = 3417
    integer, parameter ::  GlmVelocityU_                   = 3418
    integer, parameter ::  GlmVelocityV_                   = 3419
    integer, parameter ::  GlmVelocityModulus_             = 3420
    integer, parameter ::  GlmVelocityW_                   = 3421
    integer, parameter ::  BreakingWaveHeight_             = 3422
    integer, parameter ::  WaveSurfaceFluxTKE_             = 3423
    
    integer, parameter ::  WaveRad3DX_                     = 3424
    integer, parameter ::  WaveRad3DY_                     = 3425

    integer, parameter ::  WavePressureJ_                  = 3426
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!Monocromatic:
    integer, parameter ::  WaveLength_                     = 3500
    integer, parameter ::  WaveAmplitude_                  = 3501
    integer, parameter ::  WavePeriod_                     = 3502
    integer, parameter ::  WaveDirection_                  = 3503
    !!Statistical wave parametres (WW3,SWAN)
    integer, parameter ::  SignificantWaveHeight_          = 3504
    integer, parameter ::  MeanWaveLength_                 = 3505
    integer, parameter ::  MeanWavePeriod_                 = 3506
    integer, parameter ::  MeanWaveDirection_              = 3507
    integer, parameter ::  MeanDirectionalSpread_          = 3508
    integer, parameter ::  PeakFrequency_                  = 3509
    integer, parameter ::  PeakDirection_                  = 3510
    integer, parameter ::  WindSeaPeakFrequency_           = 3511
    integer, parameter ::  WindSeaPeakDirection_           = 3512
    integer, parameter ::  WaveSwellHeight_                = 3513
    integer, parameter ::  Ubw_                            = 3514
    integer, parameter ::  PeakPeriod_                     = 3520
    integer, parameter ::  PeakDirectionX_                 = 3530    
    integer, parameter ::  PeakDirectionY_                 = 3531        
    integer, parameter ::  SignificantWaveHeightBeaufort_  = 3532
    integer, parameter ::  MeanAbsoluteZeroCrossingPeriod_ = 3533
    integer, parameter ::  MeanAbsoluteWavePeriodEnergy_   = 3534
    integer, parameter ::  WavePower_                      = 3535
    integer, parameter ::  TransportEnergyX_               = 3536
    integer, parameter ::  TransportEnergyY_               = 3537
    integer, parameter ::  DirectionEnergyTransport_       = 3538
    integer, parameter ::  SmoothedPeakPeriod_             = 3539    
    integer, parameter ::  MeanAbsoluteWavePeriod_         = 3540
    
    integer, parameter ::  Swell01_SignificantWaveHeight_  = 3541
    integer, parameter ::  Swell01_WavePeriod_             = 3542
    integer, parameter ::  Swell01_WaveDirection_          = 3543
    
    integer, parameter ::  WindSea_SignificantWaveHeight_  = 3544
    integer, parameter ::  WindSea_WavePeriod_             = 3545
    integer, parameter ::  WindSea_WaveDirection_          = 3546
    
    integer, parameter ::  PeakWaveLength_                 = 3547
!____________________________________________________________________________________
!________________________________________________________exclusive use @ modulelife__

    integer, parameter :: Diatom_C_                         = 5001
    integer, parameter :: Diatom_N_                         = 5002
    integer, parameter :: Diatom_P_                         = 5003
    integer, parameter :: Diatom_Si_                        = 5004
    integer, parameter :: Diatom_Chl_                       = 5005
    integer, parameter :: Mix_Flagellate_C_                 = 5006
    integer, parameter :: Mix_Flagellate_N_                 = 5007
    integer, parameter :: Mix_Flagellate_P_                 = 5008
    integer, parameter :: Mix_Flagellate_Chl_               = 5009
    integer, parameter :: Picoalgae_C_                      = 5010
    integer, parameter :: Picoalgae_N_                      = 5011
    integer, parameter :: Picoalgae_P_                      = 5012
    integer, parameter :: Picoalgae_Chl_                    = 5013
    integer, parameter :: Flagellate_C_                     = 5014
    integer, parameter :: Flagellate_N_                     = 5015
    integer, parameter :: Flagellate_P_                     = 5016
    integer, parameter :: Flagellate_Chl_                   = 5017
    integer, parameter :: Microzooplankton_C_               = 5018
    integer, parameter :: Microzooplankton_N_               = 5019
    integer, parameter :: Microzooplankton_P_               = 5020
    integer, parameter :: Het_Nanoflagellate_C_             = 5021
    integer, parameter :: Het_Nanoflagellate_N_             = 5022
    integer, parameter :: Het_Nanoflagellate_P_             = 5023
    integer, parameter :: Mesozooplankton_C_                = 5024
    integer, parameter :: Mesozooplankton_N_                = 5025
    integer, parameter :: Mesozooplankton_P_                = 5026
    integer, parameter :: Het_Bacteria_C_                   = 5027
    integer, parameter :: Het_Bacteria_N_                   = 5028
    integer, parameter :: Het_Bacteria_P_                   = 5029
    integer, parameter :: DOCsl_                            = 5030  
    integer, parameter :: DOPsl_                            = 5031 
    integer, parameter :: DONsl_                            = 5032


    !Macroalgae rates
    integer, parameter :: MA_GrossFact_                     = 1
    integer, parameter :: MA_TLimFact_                      = 2
    integer, parameter :: MA_LLimFact_                      = 3
    integer, parameter :: MA_NutLimFact_                    = 4
    integer, parameter :: MA_SLimFact_                      = 5
    integer, parameter :: MA_NLimFact_                      = 6
    integer, parameter :: MA_PLimFact_                      = 7
    integer, parameter :: MA_Excretion_                     = 8
    integer, parameter :: MA_Respiration_                   = 9
    integer, parameter :: MA_NaturalMort_                   = 10
    integer, parameter :: MA_Grazing_                       = 11
    integer, parameter :: MA_Condition_                     = 12
    integer, parameter :: MA_CarrCapFact_                   = 13

    integer, parameter :: ConsolidationFlux_                = 9000
    integer, parameter :: Porosity_                         = 9001
    
    !Cohesive sediment fractions - ModuleDrainageNetwork
    integer, parameter ::  TSS_                             = 9101
    integer, parameter ::  COHSED_FINE_                     = 9102
    integer, parameter ::  COHSED_MEDIUM_                   = 9103
    integer, parameter ::  COHSED_COARSE_                   = 9104
    integer, parameter ::  VSS_                             = 9111
    
    !PhreeqC properties ------------------------------------------
    integer, parameter ::  WaterSaturation_                 = 10000
    integer, parameter ::  CellPorosity_                    = 10001
    integer, parameter ::  Pressure_                        = 10002
    integer, parameter ::  SolutionMapping_                 = 10003
    integer, parameter ::  EquilibriumMapping_              = 10004
    integer, parameter ::  ExchangeMapping_                 = 10005
    integer, parameter ::  SurfaceMapping_                  = 10006
    integer, parameter ::  GasPhaseMapping_                 = 10007
    integer, parameter ::  SolidSolutionMapping_            = 10008
    integer, parameter ::  KineticsMapping_                 = 10009
    
    !Solution properties
    integer, parameter :: SolutionMagnesium_               = 10000
    integer, parameter :: SolutionCalcium_                 = 10001
    integer, parameter :: SolutionSodium_                  = 10002
    integer, parameter :: SolutionNitrogenGas_             = 10003
    integer, parameter :: SolutionOxygenGas_               = 10004
    integer, parameter :: SolutionAmmonia_                 = 10005
    integer, parameter :: SolutionNitrate_                 = 10006
    integer, parameter :: SolutionNitrite_                 = 10007
    integer, parameter :: SolutionChlorine_                = 10008
    integer, parameter :: SolutionCarbon_                  = 10009
    integer, parameter :: SolutionPotassium_               = 10010
    integer, parameter :: SolutionAluminium_               = 10011
    integer, parameter :: SolutionSilicium_                = 10012
    
    !rain concentrantion properties
    integer, parameter :: RainMagnesium_                   = 10201
    integer, parameter :: RainCalcium_                     = 10202
    integer, parameter :: RainSodium_                      = 10203 
    integer, parameter :: RainChlorine_                    = 10204
    integer, parameter :: RainAmmonia_                     = 10205
    
    !Exchange properties
    integer, parameter :: eCaX2_                           = 10301
    integer, parameter :: eMgX2_                           = 10302
    integer, parameter :: eNaX_                            = 10303
    integer, parameter :: eNH4X_                           = 10304
    integer, parameter :: eKX_                             = 10305
    
    !Species properties
    integer, parameter :: sCa2_                            = 10501
    integer, parameter :: sCaOH_                           = 10502
    integer, parameter :: sH2_                             = 10503
    integer, parameter :: sMg2_                            = 10504
    integer, parameter :: sMgOH_                           = 10505
    integer, parameter :: sNH4_                            = 10506
    integer, parameter :: sNH3_                            = 10507
    integer, parameter :: sN2_                             = 10508
    integer, parameter :: sNO2_                            = 10509
    integer, parameter :: sNO3_                            = 10510
    integer, parameter :: sNa_                             = 10511
    integer, parameter :: sNaOH_                           = 10512
    integer, parameter :: sO2_                             = 10513
    
    !MasterSpecies Mass Properties
    integer, parameter :: msmSolutionCalcium_              = 10800
    integer, parameter :: msmSolutionMagnesium_            = 10801
    integer, parameter :: msmSolutionSodium_               = 10802
    integer, parameter :: msmSolutionAmmonia_              = 10803     
    
    !Gas phase properties
    integer, parameter :: GasN2_                           = 11000
    integer, parameter :: GasCO2_                          = 11001
    
    !Phases properties
    integer, parameter :: Calcite_                         = 11500 !CaCO3
    integer, parameter :: Dolomite_                        = 11501 !CaMg(CO3)2
    integer, parameter :: Aragonite_                       = 11502 !CaCO3
    integer, parameter :: Halite_                          = 11503 !NaCl
    integer, parameter :: KFeldspar_                       = 11504 !KAlSi3O8

    !Other PhreeqC properties
    integer, parameter :: pE_                              = 12000
    
    !ChainReactions 
    integer, parameter :: SoilVolumetricDensity_           = 13000
    
    !Other properties
    integer, parameter :: SolEC_                            = 14000 !Solution electrical conductivity

    !Pesticides
    integer, parameter :: GenericDissPesticide_1_           = 15001
    integer, parameter :: GenericDissPesticide_2_           = 15002
    integer, parameter :: GenericDissPesticide_3_           = 15003
    integer, parameter :: GenericDissPesticide_4_           = 15004

    integer, parameter :: GenericPartPesticide_1_           = 15005
    integer, parameter :: GenericPartPesticide_2_           = 15006
    integer, parameter :: GenericPartPesticide_3_           = 15007
    integer, parameter :: GenericPartPesticide_4_           = 15008
    
    !number of individuals per cell
    integer, parameter :: IndividualsPerCell_               = 20000
    
    !floating object
    integer, parameter :: FloatingObject_                   = 20001
    integer, parameter :: HumanBody_                        = 20002

    !HNS
    integer, parameter :: HNS_                              = 20003
    
    !WWTPQProperties
    !Oxygen, Nitrate, Ammonia and Alkalinity were already added
    integer, parameter :: SolInertOrgMat_                   = 20500
    integer, parameter :: ReadilyBioSub_                    = 20501            
    integer, parameter :: PartInertOrgMar_                  = 20502
    integer, parameter :: SlowlyBioSub_                     = 20503
    integer, parameter :: HetBio_                           = 20504
    integer, parameter :: AutBio_                           = 20505
    integer, parameter :: PartProd_                         = 20506
    integer, parameter :: SolBioOrgNitrogen_                = 20510
    integer, parameter :: PartBioOrgNitrogen_               = 20511
    !Snow
    integer, parameter :: SnowPack_                         = 30000
    integer, parameter :: DailyAvgTemp_                     = 30001
    integer, parameter :: ForestCoverFraction_              = 30002
    integer, parameter :: SnowSlopeFactor_                  = 30003
    integer, parameter :: SnowPrecipitation_                = 30004
    integer, parameter :: SnowWaterEquivalent_              = 30005
    integer, parameter :: SurfaceDownLatentHeat_            = 30006
    integer, parameter :: SurfaceDownSensibleHeat_          = 30007
    !Irrigation
    integer, parameter :: ApplicationArea_                  = 31000
    integer, parameter :: FixedIrrigation_                  = 31001
    integer, parameter :: AccIrrigation_                    = 31002
    
    
    !Percentage of a cell occupied by particles 
    integer, parameter :: CellPercentContamin_              = 40000

    !Spatial emission discharge
    integer, parameter :: DischPoint_                       = 1
    integer, parameter :: DischLine_                        = 2
    integer, parameter :: DischPolygon_                     = 3
    integer, parameter :: DischXYZPoints_                   = 4
                                                            
    !Flow discharge horizontal distribution                 
    integer, parameter :: DischByCell_                      = 1
    integer, parameter :: DischByWaterColumn_               = 2
    integer, parameter :: DischByVolume_                    = 3

    !Discharge location in the vertical direction
    integer, parameter :: DischBottom_                      = 1
    integer, parameter :: DischSurf_                        = 2
    integer, parameter :: DischDepth_                       = 3
    integer, parameter :: DischLayer_                       = 4
    integer, parameter :: DischUniform_                     = 5
    integer, parameter :: DischProfile_                     = 6


!_______________________________________________________________________________________________

    !Other Properties
    character(StringLength), private, parameter :: Char_SolEC                = 'solution electrical conductivity'
      
    !Name of PhreeqC properties
    character(StringLength), private, parameter :: Char_WaterSaturation      = 'water saturation'
    character(StringLength), private, parameter :: Char_CellPorosity         = 'cell porosity'
    character(StringLength), private, parameter :: Char_Pressure             = 'pressure'
    character(StringLength), private, parameter :: Char_SolutionMapping      = 'solution mapping'
    character(StringLength), private, parameter :: Char_EquilibriumMapping   = 'equilibrium mapping'
    character(StringLength), private, parameter :: Char_ExchangeMapping      = 'exchange mapping'
    character(StringLength), private, parameter :: Char_SurfaceMapping       = 'surface mapping'
    character(StringLength), private, parameter :: Char_GasPhaseMapping      = 'gas phase mapping'
    character(StringLength), private, parameter :: Char_SolidSolutionMapping = 'solid solution mapping'
    character(StringLength), private, parameter :: Char_KineticsMapping      = 'kinetics mapping'
    !character(StringLength), private, parameter :: Char_SolutionMagnesium    = 'solution magnesium'
    !character(StringLength), private, parameter :: Char_SolutionCalcium      = 'solution calcium'
    !character(StringLength), private, parameter :: Char_SolutionSodium       = 'solution sodium'
    !character(StringLength), private, parameter :: Char_SolutionNitrogenGas  = 'solution nitrogen gas'
    !character(StringLength), private, parameter :: Char_SolutionOxygenGas    = 'solution oxygen gas'
    !character(StringLength), private, parameter :: Char_SolutionAmmonia      = 'solution ammonia'
    !character(StringLength), private, parameter :: Char_SolutionNitrate      = 'solution nitrate'
    !character(StringLength), private, parameter :: Char_SolutionNitrite      = 'solution nitrite'
    !character(StringLength), private, parameter :: Char_SolutionChlorine     = 'solution chlorine'
    !character(StringLength), private, parameter :: Char_GasN2                = 'gas n2'
    !character(StringLength), private, parameter :: Char_GasCO2               = 'gas co2'
    !character(StringLength), private, parameter :: Char_pE                   = 'pE'
    !character(StringLength), private, parameter :: Char_eCaX2                = 'CaX2'
    !character(StringLength), private, parameter :: Char_eMgX2                = 'MgX2'
    !character(StringLength), private, parameter :: Char_eNaX                 = 'NaX'
    !character(StringLength), private, parameter :: Char_eKX                  = 'KX'
    !character(StringLength), private, parameter :: Char_eNH4X                = 'AmmoniaX' !'NH4X'
    !character(StringLength), private, parameter :: Char_sCa2                 = 'Ca+2'
    !character(StringLength), private, parameter :: Char_sCaOH                = 'CaOH+'
    !character(StringLength), private, parameter :: Char_sH2                  = 'H2'
    !character(StringLength), private, parameter :: Char_sMg2                 = 'Mg+2'
    !character(StringLength), private, parameter :: Char_sMgOH                = 'MgOH+'
    !character(StringLength), private, parameter :: Char_sNH4                 = 'Ammonia+' !'NH4+'
    !character(StringLength), private, parameter :: Char_sNH3                 = 'NH3'
    !character(StringLength), private, parameter :: Char_sN2                  = 'N2'
    !character(StringLength), private, parameter :: Char_sNO2                 = 'NO2-'
    !character(StringLength), private, parameter :: Char_sNO3                 = 'NO3-'
    !character(StringLength), private, parameter :: Char_sNa                  = 'Na+'
    !character(StringLength), private, parameter :: Char_sNaOH                = 'NaOH+'
    !character(StringLength), private, parameter :: Char_sO2                  = 'O2'
    !character(StringLength), private, parameter :: Char_msmSolutionCalcium   = 'solution calcium mass'
    !character(StringLength), private, parameter :: Char_msmSolutionMagnesium = 'solution magnesium mass'
    !character(StringLength), private, parameter :: Char_msmSolutionSodium    = 'solution sodium mass'
    !character(StringLength), private, parameter :: Char_msmSolutionAmmonia   = 'solution ammonia mass'
    !character(StringLength), private, parameter :: Char_Calcite              = 'calcite'
    !character(StringLength), private, parameter :: Char_Dolomite             = 'dolomite'
    !character(StringLength), private, parameter :: Char_Aragonite            = 'aragonite'
    !character(StringLength), private, parameter :: Char_Halite               = 'halite' 
    !character(StringLength), private, parameter :: Char_KFeldspar            = 'k-feldspar'
    !character(StringLength), private, parameter :: Char_SolutionCarbon       = 'solution carbon' 
    !character(StringLength), private, parameter :: Char_SolutionPotassium    = 'solution potassium'     
    !character(StringLength), private, parameter :: Char_SolutionAluminium    = 'solution aluminium' 
    !character(StringLength), private, parameter :: Char_SolutionSilicium     = 'solution silicium' 
    !character(StringLength), private, parameter :: Char_RainMagnesium        = 'rain magnesium' 
    !character(StringLength), private, parameter :: Char_RainCalcium          = 'rain calcium' 
    !character(StringLength), private, parameter :: Char_RainSodium           = 'rain sodium' 
    !character(StringLength), private, parameter :: Char_RainChlorine         = 'rain chlorine'
    !character(StringLength), private, parameter :: Char_RainAmmonia          = 'rain ammonia'


!_______________________________________________________________________________________________
    
    !Name of MOHID LAND ModuleIrrigation properties
    character(StringLength), private, parameter :: Char_ApplicationArea      = 'application area'
    character(StringLength), private, parameter :: Char_FixedIrrigation      = 'fixed irrigation'
    character(StringLength), private, parameter :: Char_AccIrrigation        = 'acc. irrigation'
      
    !Name of Waterproperties
    character(StringLength), private, parameter :: Char_Density              = 'density'
    character(StringLength), private, parameter :: Char_Temperature          = 'temperature'
    character(StringLength), private, parameter :: Char_Salinity             = 'salinity'
    character(StringLength), private, parameter :: Char_Phytoplankton        = 'phytoplankton'
    character(StringLength), private, parameter :: Char_Zooplankton          = 'zooplankton'
    character(StringLength), private, parameter :: Char_DOPRefractory        = 'dissolved refractory organic phosphorus'
    character(StringLength), private, parameter :: Char_DOPNon_Refractory    = 'dissolved non-refractory organic phosphorus'
    character(StringLength), private, parameter :: Char_DONRefractory        = 'dissolved refractory organic nitrogen'
    character(StringLength), private, parameter :: Char_DONNon_Refractory    = 'dissolved non-refractory organic nitrogen'
    character(StringLength), private, parameter :: Char_Inorganic_Phosphorus = 'inorganic phosphorus'


! Marky Mark 
!_______used @ modulelife__________________________________________________________________________________________________

    character(StringLength), private, parameter :: Char_POC                  = 'particulate organic carbon'
    character(StringLength), private, parameter :: Char_POP                  = 'particulate organic phosphorus'    
    character(StringLength), private, parameter :: Char_PON                  = 'particulate organic nitrogen'
    character(StringLength), private, parameter :: Char_PONRefractory        = 'particulate refractory organic nitrogen'
    character(StringLength), private, parameter :: Char_DOC                  = 'labile dissolved organic carbon'
    character(StringLength), private, parameter :: Char_DOP                  = 'labile dissolved organic phosphorus'    
    character(StringLength), private, parameter :: Char_DON                  = 'labile dissolved organic nitrogen'
    character(StringLength), private, parameter :: Char_DOCsl                = 'semi-labile dissolved organic carbon'
    character(StringLength), private, parameter :: Char_DOPsl                = 'semi-labile dissolved organic phosphorus'    
    character(StringLength), private, parameter :: Char_DONsl                = 'semi-labile dissolved organic nitrogen'
    character(StringLength), private, parameter :: Char_Ammonia              = 'ammonia'
    character(StringLength), private, parameter :: Char_Nitrate              = 'nitrate'
    character(StringLength), private, parameter :: Char_Silicate             = 'silicate acid'
    character(StringLength), private, parameter :: Char_BioSilica            = 'biogenic silica'
    character(StringLength), private, parameter :: Char_CarbonDioxide        = 'carbon dioxide'
    character(StringLength), private, parameter :: Char_Oxygen               = 'oxygen'
    character(StringLength), private, parameter :: Char_DissolO2PercentSat   = 'dissolved oxygen percent saturation'
    character(StringLength), private, parameter :: Char_CO2PartialPressure   = 'dissolved CO2 partial pressure'
    
    character(StringLength), private, parameter :: Char_Diatom_C             = 'diatoms carbon'
    character(StringLength), private, parameter :: Char_Diatom_N             = 'diatoms nitrogen'
    character(StringLength), private, parameter :: Char_Diatom_P             = 'diatoms phosphorus'
    character(StringLength), private, parameter :: Char_Diatom_Si            = 'diatoms silica'
    character(StringLength), private, parameter :: Char_Diatom_Chl           = 'diatoms chlorophyll'

    character(StringLength), private, parameter :: Char_Mix_Flagellate_C     = 'autotrophic flagellates carbon'
    character(StringLength), private, parameter :: Char_Mix_Flagellate_N     = 'autotrophic flagellates nitrogen'
    character(StringLength), private, parameter :: Char_Mix_Flagellate_P     = 'autotrophic flagellates phosphorus'
    character(StringLength), private, parameter :: Char_Mix_Flagellate_Chl   = 'autotrophic flagellates chlorophyll'

    character(StringLength), private, parameter :: Char_Picoalgae_C          = 'picoalgae carbon'
    character(StringLength), private, parameter :: Char_Picoalgae_N          = 'picoalgae nitrogen'
    character(StringLength), private, parameter :: Char_Picoalgae_P          = 'picoalgae phosphorus'
    character(StringLength), private, parameter :: Char_Picoalgae_Chl        = 'picoalgae chlorophyll'

    character(StringLength), private, parameter :: Char_Flagellate_C         = 'flagellates carbon'
    character(StringLength), private, parameter :: Char_Flagellate_N         = 'flagellates nitrogen'
    character(StringLength), private, parameter :: Char_Flagellate_P         = 'flagellates phosphorus'
    character(StringLength), private, parameter :: Char_Flagellate_Chl       = 'flagellates chlorophyll'

    character(StringLength), private, parameter :: Char_Microzooplankton_C   = 'microzooplankton carbon'
    character(StringLength), private, parameter :: Char_Microzooplankton_N   = 'microzooplankton nitrogen'
    character(StringLength), private, parameter :: Char_Microzooplankton_P   = 'microzooplankton phosphorus'

    character(StringLength), private, parameter :: Char_Het_Nanoflagellate_C = 'heterotrophic nanoflagellate carbon'
    character(StringLength), private, parameter :: Char_Het_Nanoflagellate_N = 'heterotrophic nanoflagellate nitrogen'
    character(StringLength), private, parameter :: Char_Het_Nanoflagellate_P = 'heterotrophic nanoflagellate phosphorus'

    character(StringLength), private, parameter :: Char_Mesozooplankton_C    = 'mesozooplankton carbon'
    character(StringLength), private, parameter :: Char_Mesozooplankton_N    = 'mesozooplankton nitrogen'
    character(StringLength), private, parameter :: Char_Mesozooplankton_P    = 'mesozooplankton phosphorus'

    character(StringLength), private, parameter :: Char_Het_Bacteria_C       = 'heterotrophic bacteria carbon'
    character(StringLength), private, parameter :: Char_Het_Bacteria_N       = 'heterotrophic bacteria nitrogen'
    character(StringLength), private, parameter :: Char_Het_Bacteria_P       = 'heterotrophic bacteria phosphorus'
!___________________________________________________________________________________________________________________

    !Name of Bivalve
    character(StringLength), private, parameter :: Char_Bivalve1             = 'bivalve1'
    character(StringLength), private, parameter :: Char_Bivalve2             = 'bivalve2'
    character(StringLength), private, parameter :: Char_Bivalve3             = 'bivalve3'
    character(StringLength), private, parameter :: Char_Bivalve4             = 'bivalve4'

    !Name of Bivalve Predators
    character(StringLength), private, parameter :: Char_Shrimp               = 'shrimp'
    character(StringLength), private, parameter :: Char_Crab                 = 'crab'
    character(StringLength), private, parameter :: Char_OysterCatcher        = 'oystercatcher'
    character(StringLength), private, parameter :: Char_EiderDuck            = 'eider duck'
    character(StringLength), private, parameter :: Char_HerringGull          = 'herring gull'

    character(StringLength), private, parameter :: Char_POM                  = 'particulate organic matter'

    character(StringLength), private, parameter :: Char_Nitrite              = 'nitrite'
    character(StringLength), private, parameter :: Char_BOD                  = 'biochemical oxygen demand'
    character(StringLength), private, parameter :: Char_Cohesive_Sediment    = 'cohesive sediment'
    character(StringLength), private, parameter :: Char_Fecal_Coliforms      = 'fecal coliforms'
    character(StringLength), private, parameter :: Char_E_Coli               = 'escherichia coli'
    character(StringLength), private, parameter :: Char_T90                  = 'T90'
    character(StringLength), private, parameter :: Char_T90_E_Coli           = 'T90 e.coli'

    character(StringLength), private, parameter :: Char_Oil                  = 'oil'
    character(StringLength), private, parameter :: Char_FloatingObject       = 'floating object'
    character(StringLength), private, parameter :: Char_HumanBody            = 'human body'
    character(StringLength), private, parameter :: Char_OilThickness         = 'oil thickness'
    character(StringLength), private, parameter :: Char_HNS                  = 'hns'
    character(StringLength), private, parameter :: Char_Ciliate              = 'ciliate'
    character(StringLength), private, parameter :: Char_Bacteria             = 'bacteria'
    character(StringLength), private, parameter :: Char_ParticulateArsenic   = 'particulate arsenic'
    character(StringLength), private, parameter :: Char_DissolvedArsenic     = 'dissolved arsenic'
    character(StringLength), private, parameter :: Char_Larvae               = 'larvae'
    character(StringLength), private, parameter :: Char_Age                  = 'age'
    character(StringLength), private, parameter :: Char_Fish                 = 'fish'
    character(StringLength), private, parameter :: Char_FishFood             = 'fish food'
    character(StringLength), private, parameter :: Char_MacroAlgae           = 'macroalgae'
    character(StringLength), private, parameter :: Char_DriftingMacroAlgae   = 'drifting macroalgae'
    character(StringLength), private, parameter :: Char_MicroPhytoBenthos    = 'microphytobenthos' !Rosa
    character(StringLength), private, parameter :: Char_GenericProperty      = 'generic property'
    
    ! metals
    character(StringLength), private, parameter :: Char_ParticulateMetal     = 'particulate metal'
    character(StringLength), private, parameter :: Char_DissolvedMetal       = 'dissolved metal'
    character(StringLength), private, parameter :: Char_ParticulateCopper    = 'particulate copper'
    character(StringLength), private, parameter :: Char_DissolvedCopper      = 'dissolved copper'
    character(StringLength), private, parameter :: Char_ParticulateCadmium   = 'particulate cadmium'
    character(StringLength), private, parameter :: Char_DissolvedCadmium     = 'dissolved cadmium'
    character(StringLength), private, parameter :: Char_ParticulateZinc      = 'particulate zinc'
    character(StringLength), private, parameter :: Char_DissolvedZinc        = 'dissolved zinc'
    character(StringLength), private, parameter :: Char_ParticulateMercury   = 'particulate mercury'
    character(StringLength), private, parameter :: Char_DissolvedMercury     = 'dissolved mercury'
    character(StringLength), private, parameter :: Char_ParticulateLead      = 'particulate lead'
    character(StringLength), private, parameter :: Char_DissolvedLead        = 'dissolved lead'
    
    ! benthos
    character(StringLength), private, parameter :: Char_SuspensionFeedersC    = 'suspension feeders carbon'
    character(StringLength), private, parameter :: Char_SuspensionFeedersN    = 'suspension feeders nitrogen'
    character(StringLength), private, parameter :: Char_SuspensionFeedersP    = 'suspension feeders phosphorus'
    character(StringLength), private, parameter :: Char_MicroPhytoBenthosC    = 'microphytobenthos carbon'
    character(StringLength), private, parameter :: Char_MicroPhytoBenthosN    = 'microphytobenthos nitrogen'
    character(StringLength), private, parameter :: Char_MicroPhytoBenthosP    = 'microphytobenthos phosphorus'
    character(StringLength), private, parameter :: Char_DepositFeedersC       = 'deposit feeders carbon'
    character(StringLength), private, parameter :: Char_DepositFeedersN       = 'deposit feeders nitrogen'
    character(StringLength), private, parameter :: Char_DepositFeedersP       = 'deposit feeders phosphorus'
    character(StringLength), private, parameter :: Char_SeagrassesN           = 'seagrasses internal nitrogen'
    character(StringLength), private, parameter :: Char_SeagrassesP           = 'seagrasses internal phosphorus'
    character(StringLength), private, parameter :: Char_SeagrassesLeaves      = 'seagrasses leaves'
    character(StringLength), private, parameter :: Char_SeagrassesRoots       = 'seagrasses roots'
    character(StringLength), private, parameter :: Char_LeavesUptakeN         = 'leavesuptaken'    
    character(StringLength), private, parameter :: Char_LeavesUptakeP         = 'leavesuptakep'   
    character(StringLength), private, parameter :: Char_LeavesLightFactor     = 'leaveslightfactor'
    character(StringLength), private, parameter :: Char_RootsUptakeN          = 'rootsuptaken'
    character(StringLength), private, parameter :: Char_RootsUptakeP          = 'rootsuptakep'
    character(StringLength), private, parameter :: Char_NintFactor            = 'internal nitrogen function leaves'
    character(StringLength), private, parameter :: Char_NintFactorR           = 'internal nitrogen function roots'
    character(StringLength), private, parameter :: Char_PintFactorR           = 'internal phosphorus function roots'
    character(StringLength), private, parameter :: Char_PintFactor            = 'internal phosphorus function leaves'
    character(StringLength), private, parameter :: Char_RootsMort             = 'roots mortality'

    character(StringLength), private, parameter :: Char_GrossProd            = 'grossprod'
    character(StringLength), private, parameter :: Char_NutrientLim          = 'nutrientlim'
    character(StringLength), private, parameter :: Char_NLim                 = 'nitrogenlim'
    character(StringLength), private, parameter :: Char_PLim                 = 'phosphoruslim'
    character(StringLength), private, parameter :: Char_LightLim             = 'lightlim'
    character(StringLength), private, parameter :: Char_TemperatureLim       = 'temperaturelim'
    character(StringLength), private, parameter :: Char_SalinityLim          = 'salinitylim' 
    character(StringLength), private, parameter :: Char_NetProd              = 'netprod'
    character(StringLength), private, parameter :: Char_Excretion            = 'excretion'
    character(StringLength), private, parameter :: Char_Respiration          = 'respiration'
    character(StringLength), private, parameter :: Char_NaturalMort          = 'naturalmort'
    character(StringLength), private, parameter :: Char_Grazing              = 'grazing'
    character(StringLength), private, parameter :: Char_MACondition          = 'macondition'
    character(StringLength), private, parameter :: Char_CarrCapLim           = 'carrcaplim'
    
    character(StringLength), private, parameter :: Char_DiaGrossProd         = 'diagrossprod'      
    character(StringLength), private, parameter :: Char_DiaNutrientLim       = 'dianutrientlim'    
    character(StringLength), private, parameter :: Char_DiaLightLim          = 'dialightlim'       
    character(StringLength), private, parameter :: Char_DiaTemperatureLim    = 'diatemperaturelim' 
    character(StringLength), private, parameter :: Char_Diatoms              = 'diatoms'           
    character(StringLength), private, parameter :: Char_DiaNLim              = 'dianitrogenlim'
    character(StringLength), private, parameter :: Char_DiaPLim              = 'diaphosphoruslim'
    character(StringLength), private, parameter :: Char_DiaSiLim             = 'diasilicalim'

    character(StringLength), private, parameter :: Char_AdsorbedAmmonia      = 'particulated ammonia'
    character(StringLength), private, parameter :: Char_RefreactaryOrganicN  = 'particulated refractory organic nitrogen'
    character(StringLength), private, parameter :: Char_Ngas                 = 'nitrogen gas'
    character(StringLength), private, parameter :: Char_AmmoniaGas           = 'ammonia gas'
    character(StringLength), private, parameter :: Char_Urea                 = 'urea'

    character(StringLength), private, parameter :: Char_Methane              = 'methane'

    character(StringLength), private, parameter :: Char_HeterotrophicC       = 'heterotrophic microorganism carbon'
    character(StringLength), private, parameter :: Char_HeterotrophicN       = 'heterotrophic microorganism nitrogen'
    character(StringLength), private, parameter :: Char_HeterotrophicP       = 'heterotrophic microorganism phosphorus'
    character(StringLength), private, parameter :: Char_AnaerobicC           = 'anaerobic microorganism carbon'
    character(StringLength), private, parameter :: Char_AnaerobicN           = 'anaerobic microorganism nitrogen'
    character(StringLength), private, parameter :: Char_AnaerobicP           = 'anaerobic microorganism phosphorus'
    character(StringLength), private, parameter :: Char_AutotrophicC         = 'autotrophic microorganism carbon'
    character(StringLength), private, parameter :: Char_AutotrophicN         = 'autotrophic microorganism nitrogen'
    character(StringLength), private, parameter :: Char_AutotrophicP         = 'autotrophic microorganism phosphorus'
    character(StringLength), private, parameter :: Char_SolubilizingC        = 'solubilizing microorganism carbon'
    character(StringLength), private, parameter :: Char_SolubilizingN        = 'solubilizing microorganism nitrogen'
    character(StringLength), private, parameter :: Char_SolubilizingP        = 'solubilizing microorganism phosphorus'

    character(StringLength), private, parameter :: Char_LabileOrganicC       = 'particulate labile organic carbon'
    character(StringLength), private, parameter :: Char_RefreactaryOrganicC  = 'particulated refractory organic carbon'

    character(StringLength), private, parameter :: Char_RefreactaryOrganicP  = 'particulated refractory organic phosphorus'
    character(StringLength), private, parameter :: Char_AdsorbedInorganicP   = 'particulated inorganic phosphorus'
    character(StringLength), private, parameter :: Char_SoilDryDensity       = 'soil dry density'
    character(StringLength), private, parameter :: Char_IonicStrength        = 'ionic strength'
    character(StringLength), private, parameter :: Char_PhosphorusAdsortionIndex = 'phosphorus adsortion index'

    character(StringLength), private, parameter :: Char_AutotrophicPop         = 'autotrophic microorganism population'
    character(StringLength), private, parameter :: Char_HeterotrophicPop       = 'heterotrophic microorganism population'
    character(StringLength), private, parameter :: Char_AnaerobicPop           = 'anaerobic microorganism population'
    character(StringLength), private, parameter :: Char_SolPop                 = 'solubilizing microorganism population'

    character(StringLength), private, parameter :: Char_WaterLevel_          = 'water level'
    !character(StringLength), private, parameter :: Char_WaterLevelMax_       = 'water level maximum'   
    !character(StringLength), private, parameter :: Char_WaterLevelMin_       = 'water level minimum'   
    character(StringLength), private, parameter :: Char_VelocityModulus_     = 'velocity modulus'   
    character(StringLength), private, parameter :: Char_VelocityDirection_   = 'velocity direction'   
    character(StringLength), private, parameter :: Char_FlowModulus_         = 'flow modulus'   
    character(StringLength), private, parameter :: Char_VelocityU_           = 'velocity U'
    character(StringLength), private, parameter :: Char_VelocityV_           = 'velocity V'
    character(StringLength), private, parameter :: Char_VelocityW_           = 'velocity W'
    character(StringLength), private, parameter :: Char_ShearStress_         = 'shear stress'
    character(StringLength), private, parameter :: Char_ShearStressX_        = 'shear stress X'
    character(StringLength), private, parameter :: Char_ShearStressY_        = 'shear stress Y'    
    character(StringLength), private, parameter :: Char_WaterColumn_         = 'water column'    
    character(StringLength), private, parameter :: Char_MeridionalVelocity_  = 'meridional velocity'
    character(StringLength), private, parameter :: Char_ZonalVelocity_       = 'zonal velocity'
    character(StringLength), private, parameter :: Char_TideState_           = 'tide state'



!_______used @ moduleWQ for POM pools (aquaculture cages)_____________________________________________________________

    character(StringLength), private, parameter :: Char_PON1                  = 'pon1'
    character(StringLength), private, parameter :: Char_PON2                  = 'pon2'
    character(StringLength), private, parameter :: Char_PON3                  = 'pon3'
    character(StringLength), private, parameter :: Char_PON4                  = 'pon4'
    character(StringLength), private, parameter :: Char_PON5                  = 'pon5'

    character(StringLength), private, parameter :: Char_POP1                  = 'pop1'
    character(StringLength), private, parameter :: Char_POP2                  = 'pop2'
    character(StringLength), private, parameter :: Char_POP3                  = 'pop3'
    character(StringLength), private, parameter :: Char_POP4                  = 'pop4'
    character(StringLength), private, parameter :: Char_POP5                  = 'pop5'

!_______used @ moduleWQ for output
    character(StringLength), private, parameter :: Char_PhytoChla             = 'phytoplankton Chla'


    ! guillaume nogueira
    character(StringLength), private, parameter :: Char_AltimLevelAnalyzed_         = 'water level analyzed for altimetry'
    character(StringLength), private, parameter :: Char_AltimTemperatureAnalyzed_   = 'temperature analyzed for altimetry'
    character(StringLength), private, parameter :: Char_AltimSalinityAnalyzed_      = 'salinity analyzed for altimetry'
    character(StringLength), private, parameter :: Char_AltimLevelToAssimilate_     = 'water level for altimetry assimilation'
    character(StringLength), private, parameter :: Char_VarianceFieldToAssimilate_  = 'variance field for assimilation'

    character(StringLength), private, parameter :: Char_WaterFluxX_          = 'water flux X'
    character(StringLength), private, parameter :: Char_WaterFluxY_          = 'water flux Y'
    character(StringLength), private, parameter :: Char_WaterDepth_          = 'water depth'
    character(StringLength), private, parameter :: Char_CoriolisX_           = 'coriolis X'
    character(StringLength), private, parameter :: Char_BaroclinicForceX_    = 'baroclinic force X'
    character(StringLength), private, parameter :: Char_HorizontalTransportX_= 'horizontal transport X'
    character(StringLength), private, parameter :: Char_CoriolisY_           = 'coriolis Y'
    character(StringLength), private, parameter :: Char_BaroclinicForceY_    = 'baroclinic force Y'
    character(StringLength), private, parameter :: Char_HorizontalTransportY_= 'horizontal transport Y'
    character(StringLength), private, parameter :: Char_BarotropicVelocityU_ = 'barotropic velocity U'
    character(StringLength), private, parameter :: Char_BarotropicVelocityV_ = 'barotropic velocity V'
    character(StringLength), private, parameter :: Char_BaroclinicVelocityU_ = 'baroclinic velocity U'
    character(StringLength), private, parameter :: Char_BaroclinicVelocityV_ = 'baroclinic velocity V'

    character(StringLength), private, parameter :: Char_ObstacleDragCoef     = 'obstacle drag coefficient'
    
    character(StringLength), private, parameter :: Char_ScraperVelU          = 'scraper velocity U'
    character(StringLength), private, parameter :: Char_ScraperVelV          = 'scraper velocity V'    
    character(StringLength), private, parameter :: Char_ScraperVelW          = 'scraper velocity W'    

    character(StringLength), private, parameter :: Char_VerticalZ            = 'vertical z'

    character(StringLength), private, parameter :: Char_ShearVelocity        = 'shear velocity'

    character(StringLength), private, parameter :: Char_Vorticity            = 'vorticity'
    character(StringLength), private, parameter :: Char_BaroclinicKE         = 'baroclinic KE'
    character(StringLength), private, parameter :: Char_PerturbationPE       = 'Perturbation Potential Energy'
    character(StringLength), private, parameter :: Char_KineticEnergy        = 'Kinetic Energy'

    character(StringLength), private, parameter :: Char_ParticulateContaminant='particulate contaminant'
    character(StringLength), private, parameter :: Char_DissolvedContaminant = 'dissolved contaminant'

    character(StringLength), private, parameter :: Char_Sediment             = 'sediment'

    character(StringLength), private, parameter :: Char_DissolvedSodium      = 'dissolved sodium'
    character(StringLength), private, parameter :: Char_DissolvedCalcium     = 'dissolved calcium'

    character(StringLength), private, parameter :: Char_ParticulateSodium    = 'particulate sodium'
    character(StringLength), private, parameter :: Char_ParticulateCalcium   = 'particulate calcium'

!_______used @ moduleCEQUALW2__________________________________________________________________________________________________

    character(StringLength), private, parameter :: Char_RPOM                = 'refractory particulate organic matter'
    character(StringLength), private, parameter :: Char_LPOM                = 'labile particulate organic matter'
    character(StringLength), private, parameter :: Char_LDOM                = 'labile dissolved organic matter'
    character(StringLength), private, parameter :: Char_RDOM                = 'refractory dissolved organic matter'
    character(StringLength), private, parameter :: Char_PSilica             = 'particulate silica'
    character(StringLength), private, parameter :: Char_DSilica             = 'dissolved silica'
    character(StringLength), private, parameter :: Char_ICarbon             = 'inorganic carbon'
    character(StringLength), private, parameter :: Char_pH                  = 'pH'
    character(StringLength), private, parameter :: Char_HCO3                = 'bicarbonate'
    character(StringLength), private, parameter :: Char_CO3                 = 'carbonate'
    character(StringLength), private, parameter :: Char_Algae_1             = 'algae_1'
    character(StringLength), private, parameter :: Char_Algae_2             = 'algae_2'
    character(StringLength), private, parameter :: Char_Algae_3             = 'algae_3'
    character(StringLength), private, parameter :: Char_Algae_4             = 'algae_4'
    character(StringLength), private, parameter :: Char_Algae_5             = 'algae_5'
    character(StringLength), private, parameter :: Char_Epiphyton_1         = 'epiphyton_1'
    character(StringLength), private, parameter :: Char_Epiphyton_2         = 'epiphyton_2'
    character(StringLength), private, parameter :: Char_Epiphyton_3         = 'epiphyton_3'
    character(StringLength), private, parameter :: Char_Epiphyton_4         = 'epiphyton_4'
    character(StringLength), private, parameter :: Char_Epiphyton_5         = 'epiphyton_5'
    character(StringLength), private, parameter :: Char_Alkalinity          = 'alkalinity'
    character(StringLength), private, parameter :: Char_Detritus            = 'detritus'

    ! Rates for output ---------------------------------------------------------------------------------------------------

             
    !Algal Growth Limitations  
                                    
    character(StringLength), private, parameter :: Char_ANLim               = 'ANLIM'
    character(StringLength), private, parameter :: Char_APLim               = 'APLIM'
    character(StringLength), private, parameter :: Char_ASLim               = 'ASLIM'
    character(StringLength), private, parameter :: Char_ALightLim           = 'ALIGHTLIM'
    character(StringLength), private, parameter :: Char_AOverallLim         = 'AOVERALLLIM'
    
    character(StringLength), private, parameter :: Char_AGR                 = 'AGR'
    character(StringLength), private, parameter :: Char_AMR                 = 'AMR'
    character(StringLength), private, parameter :: Char_AER                 = 'AER'
    character(StringLength), private, parameter :: Char_ARR                 = 'ARR'    

    !Epiphyte Growth Limitations                         
    character(StringLength), private, parameter :: Char_ENLim               = 'ENLIM'
    character(StringLength), private, parameter :: Char_EPLim               = 'EPLIM'
    character(StringLength), private, parameter :: Char_ESLim               = 'ESLIM'
    character(StringLength), private, parameter :: Char_ELightLim           = 'ELIGHTLIM'
    character(StringLength), private, parameter :: Char_EOverallLim         = 'EOVERALLLIM'
                                                
       
    !Decay Rates                                
    character(StringLength), private, parameter :: Char_NH4D                = 'NH4D'
    character(StringLength), private, parameter :: Char_NO3D                = 'NO3D'
    character(StringLength), private, parameter :: Char_LDOMD               = 'LDOMD'
    character(StringLength), private, parameter :: Char_RDOMD               = 'RDOMD'
    character(StringLength), private, parameter :: Char_LPOMD               = 'LPOMD'
    character(StringLength), private, parameter :: Char_RPOMD               = 'RPOMD'
    character(StringLength), private, parameter :: Char_LRDOMD              = 'LRDOMD'
    character(StringLength), private, parameter :: Char_LRPOMD              = 'LRPOMD'
    character(StringLength), private, parameter :: Char_CBODD               = 'CBODD'
                                                                            
    !phosphorus                                                             
    character(StringLength), private, parameter :: Char_PO4ER               = 'PO4ER'
    character(StringLength), private, parameter :: Char_PO4EG               = 'PO4EG'
    character(StringLength), private, parameter :: Char_PO4AR               = 'PO4AR'
    character(StringLength), private, parameter :: Char_PO4AG               = 'PO4AG'
    character(StringLength), private, parameter :: Char_PO4OM               = 'PO4OM'
    character(StringLength), private, parameter :: Char_PO4BOD              = 'PO4BOD'
        
    !ammonia
    character(StringLength), private, parameter :: Char_NH4ER               = 'NH4ER'
    character(StringLength), private, parameter :: Char_NH4EG               = 'NH4EG'
    character(StringLength), private, parameter :: Char_NH4AR               = 'NH4AR'
    character(StringLength), private, parameter :: Char_NH4AG               = 'NH4AG'
    character(StringLength), private, parameter :: Char_NH4OM               = 'NH4OM'
    character(StringLength), private, parameter :: Char_NH4BOD              = 'NH4BOD'
                                                                            
    !nitrate                                                                
    character(StringLength), private, parameter :: Char_NO3AG               = 'NO3AG'
    character(StringLength), private, parameter :: Char_NO3EG               = 'NO3EG'
                                                                            
    !dissolved silica                                                       
    character(StringLength), private, parameter :: Char_DSIAG               = 'DSIAG'
    character(StringLength), private, parameter :: Char_DSIEG               = 'DSIEG'
    character(StringLength), private, parameter :: Char_DSID                = 'DSID '
                                                                            
    !particulated silica                                                    
    character(StringLength), private, parameter :: Char_PSIAM               = 'PSIAM'
    character(StringLength), private, parameter :: Char_PSID                = 'PSID '

    !LabDom
    character(StringLength), private, parameter :: Char_LDOMAP              = 'LDOMAP'
    character(StringLength), private, parameter :: Char_LDOMEP              = 'LDOMEP'
                                                                            
                                                                            
    !labPom                                                                 
    character(StringLength), private, parameter :: Char_LPOMAP              = 'LPOMAP'
                                                                            
    !oxygen                                                                 
    character(StringLength), private, parameter :: Char_DOAP                = 'DOAP'
    character(StringLength), private, parameter :: Char_DOEP                = 'DOEP'
    character(StringLength), private, parameter :: Char_DOAR                = 'DOAR'
    character(StringLength), private, parameter :: Char_DOER                = 'DOER'
    character(StringLength), private, parameter :: Char_DOOM                = 'DOOM'
    character(StringLength), private, parameter :: Char_DONIT               = 'DONIT'

    !inorganic Carbon
    character(StringLength), private, parameter :: Char_ICarbonAP           = 'ICARBONAP'
    character(StringLength), private, parameter :: Char_ICarbonEP           = 'ICARBONEP'
    character(StringLength), private, parameter :: Char_ICarbonBOD          = 'ICARBONBOD'
!________________________________________________________________________________________________________________________   
    
    !Zero dimensional modules
    character(StringLength), parameter          :: WaterQualityModel        = 'WaterQuality'
    character(StringLength), parameter          :: SedimentQualityModel     = 'SedimentQualityModel'
    character(StringLength), parameter          :: LifeModel                = 'LifeModel'
    character(StringLength), parameter          :: BFMModel                 = 'BFMModel'
    character(StringLength), parameter          :: CEQUALW2Model            = 'CEQUALW2'
    character(StringLength), parameter          :: BenthicCEQUALW2Model     = 'BenthicCEQUALW2'
    character(StringLength), parameter          :: BenthosModel             = 'Benthos'
    character(StringLength), parameter          :: TempqsimModel            = 'Tempqsim'
    character(StringLength), parameter          :: MacroAlgaeModel          = 'MacroAlgae'
    character(StringLength), parameter          :: PhreeqCModel             = 'PhreeqCModel'
    character(StringLength), parameter          :: BoxDifModel              = 'BoxDifModel'
    character(StringLength), parameter          :: BenthicEcologyModel      = 'BenthicEcology'
    character(StringLength), parameter          :: WWTPQModel               = 'WWTPQ'
    character(StringLength), parameter          :: SeagrassWaterInteractionModel  = 'SeagrassWaterInteraction'
    character(StringLength), parameter          :: SeagrassSedimInteractionModel  = 'SeagrassSedimInteraction'
    character(StringLength), parameter          :: BivalveModel             = 'BivalveModel'

    !Water air interface
    character(StringLength), private, parameter :: Char_LatentHeat               = 'latent heat'
    character(StringLength), private, parameter :: Char_SensibleHeat             = 'sensible heat'
    character(StringLength), private, parameter :: Char_Evaporation              = 'evaporation'
    character(StringLength), private, parameter :: Char_NetLongWaveRadiation     = 'net long wave radiation'
    character(StringLength), private, parameter :: Char_UpwardLongWaveRadiation  = 'upward long wave radiation'
    character(StringLength), private, parameter :: Char_DownwardLongWaveRadiation= 'downward long wave radiation'

    character(StringLength), private, parameter :: Char_ShortWaveSolarRadiation  = 'short wave solar radiation'
    character(StringLength), private, parameter :: Char_LongWaveSolarRadiation   = 'long wave solar radiation'

    character(StringLength), private, parameter :: Char_ShortWaveSolarRadiaExtin = 'short wave solar radiation extinction'
    character(StringLength), private, parameter :: Char_LongWaveSolarRadiaExtin  = 'long wave solar radiation extinction'

    character(StringLength), private, parameter :: Char_OxygenFlux               = 'oxygen flux'
    character(StringLength), private, parameter :: Char_SpecificOxygenFlux       = 'specific oxygen flux'
    character(StringLength), private, parameter :: Char_CarbonDioxideFlux        = 'carbon dioxide flux'
    character(StringLength), private, parameter :: Char_SpecificCarbonDioxideFlux= 'specific carbon dioxide flux'
    character(StringLength), private, parameter :: Char_AmmoniaFlux              = 'ammonia flux'
    character(StringLength), private, parameter :: Char_NitrateFlux              = 'nitrate flux'
    character(StringLength), private, parameter :: Char_WindShearVelocity        = 'wind shear velocity'
    character(StringLength), private, parameter :: Char_SurfaceRadiation         = 'surface radiation'
    character(StringLength), private, parameter :: Char_WindStressX              = 'wind stress X'
    character(StringLength), private, parameter :: Char_WindStressY              = 'wind stress Y'
    character(StringLength), private, parameter :: Char_WindStress               = 'wind stress'
    character(StringLength), private, parameter :: Char_SurfaceWaterFlux         = 'surface water flux'
    character(StringLength), private, parameter :: Char_NonSolarFlux             = 'non solar flux'
    character(StringLength), private, parameter :: Char_TurbulentKineticEnergy   = 'turbulent kinetic energy'
    character(StringLength), private, parameter :: Char_Albedo                   = 'albedo'

    !Atmosphere
    character(StringLength), private, parameter :: Char_WindVelocityX            = 'wind velocity X'
    character(StringLength), private, parameter :: Char_WindVelocityY            = 'wind velocity Y'
    character(StringLength), private, parameter :: Char_WindVelocity             = 'wind velocity' !vectorial
    character(StringLength), private, parameter :: Char_SolarRadiation           = 'solar radiation'
    character(StringLength), private, parameter :: Char_Precipitation            = 'precipitation'
    character(StringLength), private, parameter :: Char_AtmosphericPressure      = 'atmospheric pressure'
    character(StringLength), private, parameter :: Char_AirTemperature           = 'air temperature'
    character(StringLength), private, parameter :: Char_RelativeHumidity         = 'relative humidity'
    character(StringLength), private, parameter :: Char_WindModulos              = 'wind modulos'
    character(StringLength), private, parameter :: Char_WindAngle                = 'wind angle'
    character(StringLength), private, parameter :: Char_CloudCover               = 'cloud cover'
    character(StringLength), private, parameter :: Char_Irrigation               = 'irrigation'
    character(StringLength), private, parameter :: Char_SunHours                 = 'sunshine hours'
    character(StringLength), private, parameter :: Char_ATMTransmitivity         = 'atmospheric transmitivity'
    character(StringLength), private, parameter :: Char_PBLHeight                = 'pbl height'
    character(StringLength), private, parameter :: Char_Reflectivity             = 'reflectivity'    
    

    character(StringLength), private, parameter :: Char_MeanSeaLevelPressure     = 'mean sea level pressure'
    character(StringLength), private, parameter :: Char_WindModulus              = 'wind modulus'
    character(StringLength), private, parameter :: Char_WindModulusBeaufort      = 'wind modulus beaufort'    
    character(StringLength), private, parameter :: Char_WindDirection            = 'wind direction'
    character(StringLength), private, parameter :: Char_SpecificHumidity         = 'specific humidity'    
    character(StringLength), private, parameter :: Char_WindGust                 = 'wind gust'  
    !Air quality    
    character(StringLength), private, parameter :: Char_CO2AtmosphericPressure   = 'CO2 atmospheric pressure'    
    character(StringLength), private, parameter :: Char_O2AtmosphericPressure    = 'O2 atmospheric pressure'    
    character(StringLength), private, parameter :: Char_HydrogenSulfide          = 'hydrogen sulfide'           ! H2S
    character(StringLength), private, parameter :: Char_MethylMercaptan          = 'methyl mercaptan'           ! or methanethiol 
    character(StringLength), private, parameter :: Char_AtmospDeposOxidNO3       = 'atmospheric deposition oxidized NO3' !LLP
    character(StringLength), private, parameter :: Char_AtmospDeposReduNH4       = 'atmospheric deposition reduced NH4'  !LLP
    character(StringLength), private, parameter :: Char_Visibility               = 'visibility'
    character(StringLength), private, parameter :: Char_Dust                     = 'dust'

    

    
    !Sand Transport
    character(StringLength), private, parameter :: Char_Diameter                 = 'sand diameter'
    character(StringLength), private, parameter :: Char_Percentage               = 'sand percentage'
    character(StringLength), private, parameter :: Char_D35                      = 'D35'
    character(StringLength), private, parameter :: Char_D50                      = 'D50'
    character(StringLength), private, parameter :: Char_D90                      = 'D90'
    character(StringLength), private, parameter :: Char_BedRock                  = 'bed rock'
    character(StringLength), private, parameter :: Char_SandTauCritic            = 'sand tau critic'
    character(StringLength), private, parameter :: Char_TransportCapacity        = "transport capacity"
    character(StringLength), private, parameter :: Char_TransportCapacityX       = "transport capacity X"
    character(StringLength), private, parameter :: Char_TransportCapacityY       = "transport capacity Y"
    character(StringLength), private, parameter :: Char_BottomEvolution          = "bottom evolution"
    character(StringLength), private, parameter :: Char_Newbathymetry            = "new bathymetry"
    character(StringLength), private, parameter :: Char_Sand                     = "sand"
    character(StringLength), private, parameter :: Char_bathymetry               = "bathymetry"    
    character(StringLength), private, parameter :: Char_MappDZ                   = "mapping DZ"

    !wave dynamics
    character(StringLength), private, parameter :: Char_WaveStressX              = 'wave stress X'
    character(StringLength), private, parameter :: Char_WaveStressY              = 'wave stress Y'   
    character(StringLength), private, parameter :: Char_WaveStress               = 'wave stress' 
    character(StringLength), private, parameter :: Char_CurrentX                 = 'Current X'
    character(StringLength), private, parameter :: Char_CurrentY                 = 'Current Y'
    character(StringLength), private, parameter :: Char_WaveX                    = 'wave_x'
    character(StringLength), private, parameter :: Char_WaveY                    = 'wave_y'
    
! Modified by Matthias DELPEY - 24/06/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 21/07/2011 - 04/08/2011 - 05/09/2011 - 25/10/2011 - 14/12/2011 
!                             - 16/12/2011 - 02/03/2012
    character(StringLength), private, parameter :: Char_WaveRad3DX               = 'exp radiation stress X'
    character(StringLength), private, parameter :: Char_WaveRad3DY               = 'exp radiation stress Y' 
    
    character(StringLength), private, parameter :: Char_WavePressureJ            = 'wave induced pressure J'
    character(StringLength), private, parameter :: Char_AtmToWaveMomentumU       = 'atmosphere to wave momentum flux X'
    character(StringLength), private, parameter :: Char_AtmToWaveMomentumV       = 'atmosphere to wave momentum flux Y'
    character(StringLength), private, parameter :: Char_WaveToOceanMomentumU     = 'wave to ocean momentum flux X'
    character(StringLength), private, parameter :: Char_WaveToOceanMomentumV     = 'wave to ocean momentum flux Y'
    character(StringLength), private, parameter :: Char_WaveDriftSpecU           = 'stokes drift spectrum X'
    character(StringLength), private, parameter :: Char_WaveDriftSpecV           = 'stokes drift spectrum Y'  
    character(StringLength), private, parameter :: Char_StokesDriftU             = 'velocity Ustokes'
    character(StringLength), private, parameter :: Char_StokesDriftV             = 'velocity Vstokes'
    character(StringLength), private, parameter :: Char_StokesDriftModulus       = 'velocity modulus stokes'
    character(StringLength), private, parameter :: Char_StokesDriftW             = 'velocity Wstokes'  
    character(StringLength), private, parameter :: Char_GlmVelocityU             = 'velocity Uglm'
    character(StringLength), private, parameter :: Char_GlmVelocityV             = 'velocity Vglm'
    character(StringLength), private, parameter :: Char_GlmVelocityModulus       = 'velocity modulus glm'
    character(StringLength), private, parameter :: Char_GlmVelocityW             = 'velocity Wglm'
    character(StringLength), private, parameter :: Char_BreakingWaveHeight       = 'breaking wave height'
    character(StringLength), private, parameter :: Char_WaveSurfaceFluxTKE       = 'wave induced surface TKE flux' 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!Monocromatic:
    character(StringLength), private, parameter :: Char_WaveLength               = 'wave length'
    character(StringLength), private, parameter :: Char_WaveAmplitude            = 'wave amplitude'
    character(StringLength), private, parameter :: Char_WavePeriod               = 'wave period'
    character(StringLength), private, parameter :: Char_WaveDirection            = 'wave direction'
    !!Statistical wave parametres (WW3,SWAN)
    character(StringLength), private, parameter :: Char_SignificantWaveHeight    = 'significant wave height'
    character(StringLength), private, parameter :: Char_SignificantWaveHeightBeaufort    = 'significant wave height beaufort'
    character(StringLength), private, parameter :: Char_MeanWaveLength           = 'mean wave length'
    character(StringLength), private, parameter :: Char_Ubw                      = 'Ubw'
    character(StringLength), private, parameter :: Char_MeanWavePeriod           = 'mean wave period'
    character(StringLength), private, parameter :: Char_MeanWaveDirection        = 'mean wave direction'
    character(StringLength), private, parameter :: Char_MeanDirectionalSpread    = 'mean directional spread'
    character(StringLength), private, parameter :: Char_PeakFrequency            = 'peak frequency'
    character(StringLength), private, parameter :: Char_PeakDirection            = 'peak direction'
    character(StringLength), private, parameter :: Char_PeakDirectionX           = 'peak direction x'    
    character(StringLength), private, parameter :: Char_PeakDirectionY           = 'peak direction y'    
    character(StringLength), private, parameter :: Char_PeakPeriod               = 'peak period'    
    character(StringLength), private, parameter :: Char_WindSeaPeakFrequency     = 'wind sea peak frequency'
    character(StringLength), private, parameter :: Char_WindSeaPeakDirection     = 'wind sea peak direction'
    character(StringLength), private, parameter :: Char_WaveSwellHeight          = 'wave swell height'
    character(StringLength), private, parameter :: Char_MeanAbsoluteZeroCrossingPeriod   = 'mean absolute zero-crossing period'
    character(StringLength), private, parameter :: Char_MeanAbsoluteWavePeriodEnergy     = 'mean absolute wave period energy'
    character(StringLength), private, parameter :: Char_WavePower                = 'wave power'    
    character(StringLength), private, parameter :: Char_TransportEnergyX         = 'transport energy X'
    character(StringLength), private, parameter :: Char_TransportEnergyY         = 'transport energy Y'
    character(StringLength), private, parameter :: Char_DirectionEnergyTransport = 'direction energy transport'
    character(StringLength), private, parameter :: Char_SmoothedPeakPeriod       = 'smoothed peak period'
    character(StringLength), private, parameter :: Char_MeanAbsoluteWavePeriod   = 'mean absolute wave period'
    character(StringLength), private, parameter :: Char_PeakWaveLength           = 'peak wave length'

    character(StringLength), private, parameter :: Char_Swell01_SignificantWaveHeight = 'primary swell significant wave height'
    character(StringLength), private, parameter :: Char_Swell01_WavePeriod            = 'primary swell wave period'
    character(StringLength), private, parameter :: Char_Swell01_WaveDirection         = 'primary swell wave direction'
    character(StringLength), private, parameter :: Char_WindSea_SignificantWaveHeight = 'wind sea significant wave height'
    character(StringLength), private, parameter :: Char_WindSea_WavePeriod            = 'wind sea wave period'
    character(StringLength), private, parameter :: Char_WindSea_WaveDirection         = 'wind sea wave direction'
    !Consolidation
    character(StringLength), private, parameter :: Char_ConsolidationFlux        = 'consolidation flux'
    character(StringLength), private, parameter :: Char_Porosity                 = 'porosity'

    !Basin
    character(StringLength), private, parameter :: Char_RefEvapotrans            = 'reference evapotranspiration'
    character(StringLength), private, parameter :: Char_TotalPlantBiomass        = 'total plant biomass'    
    character(StringLength), private, parameter :: Char_TotalPlantNitrogen       = 'total plant nitrogen'
    character(StringLength), private, parameter :: Char_TotalPlantPhosphorus     = 'total plant phosphorus'
    character(StringLength), private, parameter :: Char_RootBiomass              = 'root biomass'
    character(StringLength), private, parameter :: Char_RootDepth                = 'root depth'
    character(StringLength), private, parameter :: Char_LeafAreaIndex            = 'leaf area index'
    character(StringLength), private, parameter :: Char_SpecificLeafStorage      = 'specific leaf storage'
    character(StringLength), private, parameter :: Char_EVTPCropCoefficient      = 'crop coefficient'
    character(StringLength), private, parameter :: Char_CanopyHeight             = 'canopy height'
    character(StringLength), private, parameter :: Char_PotLeafAreaIndex         = 'potential leaf area index'
    character(StringLength), private, parameter :: Char_BoundaryLeafAreaIndex    = 'boundary leaf area index'

    !Cohesive Fractions - Drainage Network
    character(StringLength), private, parameter :: Char_TSS                      = 'TSS'
    character(StringLength), private, parameter :: Char_Cohsed_fine              = 'cohesive sediment fine'
    character(StringLength), private, parameter :: Char_Cohsed_medium            = 'cohesive sediment medium'
    character(StringLength), private, parameter :: Char_Cohsed_coarse            = 'cohesive sediment coarse'
    character(StringLength), private, parameter :: Char_VSS                      = 'VSS'

    !ChainReactions
    character(StringLength), private, parameter :: Char_SoilVolumetricDensity    = 'soil volumetric density'
    
    !Pesticides
    character(StringLength), private, parameter :: Char_GenericDissPesticide_1       = 'generic dissolved pesticide 1'
    character(StringLength), private, parameter :: Char_GenericDissPesticide_2       = 'generic dissolved pesticide 2'
    character(StringLength), private, parameter :: Char_GenericDissPesticide_3       = 'generic dissolved pesticide 3'
    character(StringLength), private, parameter :: Char_GenericDissPesticide_4       = 'generic dissolved pesticide 4'
    character(StringLength), private, parameter :: Char_GenericPartPesticide_1       = 'generic particulate pesticide 1'
    character(StringLength), private, parameter :: Char_GenericPartPesticide_2       = 'generic particulate pesticide 2'
    character(StringLength), private, parameter :: Char_GenericPartPesticide_3       = 'generic particulate pesticide 3'
    character(StringLength), private, parameter :: Char_GenericPartPesticide_4       = 'generic particulate pesticide 4'
    character(StringLength), private, parameter :: Char_IndividualsPerCell           = 'individuals per cell' 
    character(StringLength), private, parameter :: Char_CellPercentContamin          = 'cell percentage contaminated'
    
    !Snow
    character(StringLength), private, parameter :: Char_SnowPack                     = 'snow pack' 
    character(StringLength), private, parameter :: Char_DailyAvgTemp                 = 'daily average temperature'
    character(StringLength), private, parameter :: Char_ForestCoverFraction          = 'forest cover fraction'
    character(StringLength), private, parameter :: Char_SnowSlopeFactor              = 'snow slope factor'
    character(StringLength), private, parameter :: Char_SnowPrecipitation            = 'snow precipitation'
    character(StringLength), private, parameter :: Char_SnowWaterEquivalent          = 'snow water equivalent'
    character(StringLength), private, parameter :: Char_SurfaceDownLatentHeat        = 'surface downward latent heat'
    character(StringLength), private, parameter :: Char_SurfaceDownSensibleHeat      = 'surface downward sensible heat'    

    !WWTPQProperties
    character(StringLength), private, parameter :: Char_SolInertOrgMat               = 'soluble inert organic matter'
    character(StringLength), private, parameter :: Char_ReadilyBioSub                = 'readily biodegradable substrate'
    character(StringLength), private, parameter :: Char_PartInertOrgMar              = 'particualte inert organic matter'
    character(StringLength), private, parameter :: Char_SlowlyBioSub                 = 'slowly biodegradable substrate'
    character(StringLength), private, parameter :: Char_HetBio                       = 'activate heterotrophic biomass'
    character(StringLength), private, parameter :: Char_AutBio                       = 'activate autotrophic biomass'
    character(StringLength), private, parameter :: Char_PartProd                     = 'particulate products from biomass decay'
    character(StringLength), private, parameter :: Char_SolBioOrgNitrogen            = 'soluble biodegradable organic nitrogen'
    character(StringLength), private, parameter :: Char_PartBioOrgNitrogen           = 'particulate biodegradable organic matter'
    
    !Mapping
    integer, parameter :: Compute         = 1
    integer, parameter :: Covered         = 1
    integer, parameter :: Not_Covered     = 0
    integer, parameter :: Boundary        = 1
    integer, parameter :: Not_Boundary    = 0
    integer, parameter :: Exterior        = 1
    integer, parameter :: WaterPoint      = 1
    integer, parameter :: OpenPoint       = 1
    integer, parameter :: Imposed         = 1

    integer, parameter :: UndefinedPoint   = -999
    integer, parameter :: NoBasinPoint     = 0
    integer, parameter :: BasinPoint       = 1
    integer, parameter :: NoLakePoint      = 0
    integer, parameter :: LakePoint        = 1

    !Open Files
    integer, parameter :: CLOSE_FILE    = 0
    integer, parameter :: OPEN_FILE     = 1

    !Physic constants
    real(8), parameter  :: EarthMass_               = 5.9742e24
    real(8), parameter  :: GravitationalConstant_   = 6.67259e-11
    real(8), parameter  :: mu_                      = 3.9863387178e14 !is the product of the two latter constants
    real(8), parameter  :: EqRadius_                = 6378135.d0 !m
    real(8), parameter  :: MeanRadius_              = 6371000.d0 ! m

    real,    parameter  :: Gravity            = 9.81
    real,    parameter  :: Const_VonKarman    = 0.4
    real,    parameter  :: Air_Density        = 1.2
    real,    parameter  :: WaterDynamicVisc   = 1e-3
    real,    parameter  :: WaterCinematicVisc = 1e-6
    real,    parameter  :: SigmaDensityReference = 1e3
    
    !Reference atmospheric pressure at sea level in Pa
    real,    parameter  :: AtmPressSeaLevelReference = 101325 

    !PV = nRT - thermodynamic (n - number of moles)
    real,    parameter  :: R                  = 8.3144  

    !Avogadro    
    real,    parameter  :: Mole               = 6.02e23    

    !Default specific heat
    real,    parameter  :: SpecificHeatDefault = 4200.  ![joules/kg/Kelvin]

    !Angles                                 
    real(8), parameter  :: Pi               = 3.1415926535897932384626433832795
    ! ARC RADIAN OF 1 DEGREE
    real(8), parameter  :: RAD_DEG          = 0.01745329252
    !Angle units
    integer, parameter  :: Degree_          = 1
    integer, parameter  :: Radian_          = 2
    

    !Zero Degrees Kelvin
    real, parameter     :: AbsoluteZero     = 273.15

    !Error Magnitudes
    integer, parameter :: FATAL_       = 1
    integer, parameter :: WARNING_     = 2

    !Error Types
    integer, parameter :: INTERNAL_    = 1
    integer, parameter :: OUT_OF_MEM_  = 2
    integer, parameter :: KEYWORD_     = 3

    ! Vertical communication between nested models 
    integer, parameter :: FatherSonDifDim    = 0
    integer, parameter :: FatherSonEqualDim  = 1
    integer, parameter :: Father2DSon3D      = 2
    integer, parameter :: Father3DSon2D      = 3

    !Methodologies use to compute sedimentation velocity of fine sediments
    integer, parameter :: WSConstant = 1, SPMFunction = 2, WSSecondaryClarifier = 3, WSPrimaryClarifier = 4, WSFlocs = 6
    
    !Methodologies use to compute sedimentation velocity of sand
    integer, parameter :: WSSand = 5

    !Advection 1D parameters
    integer, parameter :: UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4, CentralDif = 5, LeapFrog = 6
    integer, parameter :: MinMod = 1, VanLeer = 2, Muscl = 3, Superbee = 4, PDM = 5
    real,    parameter :: MinValue = 1.e-16

    !Transport Parameters
    integer, parameter :: velocityX = 1, velocityY = 2, velocityZ = 3, massProperty = 4
    integer, parameter :: NearestNeighbour = 1, centered = 2 !other interpolation methods here 
    
    !Interpolation 2D
    integer, parameter                                      :: Bilinear2D_         = 1
    integer, parameter                                      :: NearestNeighbor2D_  = 2
    
    !Extrapolation parameters
    integer, parameter :: ExtrapolAverage_ = 1, ExtrapolNearstCell_ = 2, ExtrapolConstant_ = 3

    !Filter grid data 
    integer, parameter :: NoFilter = 0, ModifyLax = 1


    !Transport options
    real,    parameter                          :: ExplicitScheme           = 0.
    real,    parameter                          :: ImplicitScheme           = 1.


    !Density compute methods
    integer, parameter                          :: LeendertseState_         = 1
    integer, parameter                          :: UNESCOState_             = 2
    integer, parameter                          :: Linear_                  = 3
    integer, parameter                          :: Mel96State_              = 4
    integer, parameter                          :: JMD95State_              = 5
    integer, parameter                          :: ConstantDensity_         = 6
    integer, parameter                          :: WangState_               = 7

    !Methods to compute baroclinic force 
    integer, parameter                          :: DensityUniform           = 1
    integer, parameter                          :: DensityLinear            = 2
    integer, parameter                          :: Leibniz                  = 3
    integer, parameter                          :: Leibniz2                 = 4
    integer, parameter                          :: MARSALEIX                = 5

    !Datums 
    integer, parameter                          :: CLARKE_1866_DATUM        = 1
    integer, parameter                          :: GRS_80_DATUM             = 2
    integer, parameter                          :: WGS_84_DATUM             = 3
    !Elipoid Hayford 1909 (or International 1909)
    integer, parameter                          ::  ED_50_DATUM             = 4
    integer, parameter                          :: SPHERE_DATUM             = 5

       

    !Methods to compute oxygen
    integer, parameter                          :: Broecker_et_al_1978      = 1
    integer, parameter                          :: Gelda_et_al_1996         = 2
    integer, parameter                          :: Banks_Herrera_1977       = 3
    integer, parameter                          :: Wanninkhof_et_al_1991    = 4
    integer, parameter                          :: Chen_Kanwisher_1963      = 5
    integer, parameter                          :: Cole_Buchak_1993         = 6
    integer, parameter                          :: Banks_1975               = 7
    integer, parameter                          :: Smith_1978               = 8
    integer, parameter                          :: Liss_1973                = 9
    integer, parameter                          :: Downing_Truesdale_1955   = 10
    integer, parameter                          :: Kanwisher_1963           = 11
    integer, parameter                          :: Yu_et_al_1977            = 12
    integer, parameter                          :: Weiler_1974              = 13

    !Methods to compute carbon dioxide fluxes
    integer, parameter                          :: Borges_et_al_2004        = 1
    integer, parameter                          :: Carini_et_al_1996        = 2      
    integer, parameter                          :: Raimond_Cole_2001        = 3
    

    !Number of input instances of the lagrangian module 
    integer, parameter                          :: TotalLagInst_            = 11    


    ! HNS Particle State
    integer, parameter :: Air_Volatilized_            = 1
    integer, parameter :: Air_Evaporated_             = 2
    integer, parameter :: Surface_                    = 3
    integer, parameter :: WaterColumn_Droplet_        = 4
    integer, parameter :: WaterColumn_Dissolved_      = 5
    integer, parameter :: WaterColumn_Sedimented_     = 6
    integer, parameter :: Bottom_Deposited_           = 7
    integer, parameter :: Beached_                    = 8

    !Methods to compute Droplets D50
    integer, parameter  :: UserDefined_             = 1
    integer, parameter  :: Computed_Half_D50_       = 2
    integer, parameter  :: Computed_Classes_Random_ = 3

    !Module IDs
    integer, parameter ::  mGLOBALDATA_             =  1        
    integer, parameter ::  mTIME_                   =  2 
    integer, parameter ::  mENTERDATA_              =  3 
    integer, parameter ::  mFUNCTIONS_              =  4 
    integer, parameter ::  mLUD_                    =  5 
    integer, parameter ::  mWATERQUALITY_           =  6
    integer, parameter ::  mGRIDDATA_               =  7 
    integer, parameter ::  mHDF5_                   =  8            
    integer, parameter ::  mHORIZONTALGRID_         =  9
    integer, parameter ::  mHORIZONTALMAP_          = 10 
    integer, parameter ::  mGEOMETRY_               = 11 
    integer, parameter ::  mSTOPWATCH_              = 12 
    integer, parameter ::  mADVECTIONDIFFUSION_     = 13 
    integer, parameter ::  mMAP_                    = 14 
    integer, parameter ::  mTIMESERIE_              = 15 
    integer, parameter ::  mSOIL_                   = 16
    integer, parameter ::  mBOXDIF_                 = 17
    integer, parameter ::  mASSIMILATION_           = 18
    integer, parameter ::  mSTATISTIC_              = 19
    integer, parameter ::  mDISCHARGES_             = 20
    integer, parameter ::  mATMOSPHERE_             = 21
    integer, parameter ::  mHYDROINTEGRATION_       = 22
    integer, parameter ::  mHYDRODYNAMICFILE_       = 23
    integer, parameter ::  mTOGA_                   = 24
    integer, parameter ::  mGAUGE_                  = 25
    integer, parameter ::  mTRIANGULATION_          = 26
    integer, parameter ::  mGOTM_                   = 27
    integer, parameter ::  mTURBGOTM_               = 28
    integer, parameter ::  mTURBULENCE_             = 29
    integer, parameter ::  mOPENBOUNDARY_           = 30
    integer, parameter ::  mHYDRODYNAMIC_           = 31
    integer, parameter ::  mFREEVERTICALMOVEMENT_   = 32
    integer, parameter ::  mSEDIMENTQUALITY_        = 33
    integer, parameter ::  mINTERFACE_              = 34
    integer, parameter ::  mWATERPROPERTIES_        = 35
    integer, parameter ::  mOIL_                    = 36
    integer, parameter ::  mJET_                    = 37
    integer, parameter ::  mLAGRANGIAN_             = 38
    integer, parameter ::  mMODEL_                  = 39
    integer, parameter ::  mCONSOLIDATION_          = 40
    integer, parameter ::  mSEDIMENTPROPERTIES_     = 41
    integer, parameter ::  mLIFE_                   = 42
    integer, parameter ::  mCEQUALW2_               = 43
    integer, parameter ::  mBASINGEOMETRY_          = 44
    integer, parameter ::  mRUNOFF_                 = 45
    integer, parameter ::  mDRAINAGENETWORK_        = 46
    integer, parameter ::  mINTERFACEWATERAIR_      = 47
    integer, parameter ::  mINTERFACESEDIMENTWATER_ = 48
    integer, parameter ::  mLIGHTEXTINCTION_        = 49
    integer, parameter ::  mRIVERHYDRODYNAMIC_      = 50
    integer, parameter ::  mFILLMATRIX_             = 51
    integer, parameter ::  mWAVES_                  = 52
    integer, parameter ::  mBASIN_                  = 53
    integer, parameter ::  mINFILTRATION_           = 54
    integer, parameter ::  mSOILPROPERTIES_         = 55
    integer, parameter ::  mSOILPLANTAIR_           = 56
    integer, parameter ::  mSOILMACROPORES_         = 57
    integer, parameter ::  mSAND_                   = 58
    integer, parameter ::  mMACROPOREPROPERTIES_    = 59
    integer, parameter ::  mPOROUSMEDIA_            = 60
    integer, parameter ::  mSWAT_                   = 61
    integer, parameter ::  mPROFILE_                = 62
    integer, parameter ::  mBENTHOS_                = 63
    integer, parameter ::  mCLIMATOLOGY_            = 64
    integer, parameter ::  mFIREINDEX_              = 65
    integer, parameter ::  mINTERPOLATION_          = 66
    integer, parameter ::  mVEGETATION_             = 67
    integer, parameter ::  mRESERVOIROPTIMIZATION_  = 68
    integer, parameter ::  mMACROALGAE_             = 69
    integer, parameter ::  mBFM_                    = 70
    integer, parameter ::  mNETCDF_                 = 71
    integer, parameter ::  mSEQUENTIALASSIMILATION_ = 72
    integer, parameter ::  mPOROUSMEDIAPROPERTIES_  = 73
    integer, parameter ::  mPHREEQC_                = 74
    integer, parameter ::  mRUNOFFPROPERTIES_       = 75
    integer, parameter ::  mCHAINREACTIONS_         = 76
    integer, parameter ::  mCUDA_                   = 77
    integer, parameter ::  mField4D_                = 78
    integer, parameter ::  mBENTHICECOLOGY_         = 79
    integer, parameter ::  mWWTPQ_                  = 80
    integer, parameter ::  mSEAGRASSSEDIMINTERAC_   = 81     ! Isabella
    integer, parameter ::  mSEAGRASSWATERINTERAC_   = 82     ! Isabella
    integer, parameter ::  mBivalve_                = 83
    integer, parameter ::  mTimeSeriesAnalyser_     = 84
    integer, parameter ::  mNetworkStatistics_      = 85
    integer, parameter ::  mTimeSeriesOperator_     = 86
    integer, parameter ::  mPressureDifferences_    = 87
    integer, parameter ::  mAnalytical_LDS_         = 88
    integer, parameter ::  mHNS_                    = 89
    integer, parameter ::  mGlueWW3_OBC_            = 90
    integer, parameter ::  mSnow_                   = 91
    integer, parameter ::  mSediment_               = 92
    integer, parameter ::  mReservoirs_             = 93
    integer, parameter ::  mIrrigation_             = 94
    integer, parameter ::  mTURBINE_                = 95
    integer, parameter ::  mLitter_                 = 96
    integer, parameter ::  mTwoWay_                 = 97
    integer, parameter ::  mOutputGrid_             = 98
    
    !Domain decomposition
    integer, parameter :: WestSouth        = 1
    integer, parameter :: EastNorth        = 2
    integer, parameter :: SouthNorth_      = 4
    integer, parameter :: WestEast_        = 5
    
    !Parameter
    integer, parameter :: WaveNameLength = 5
    integer, parameter :: NComponents    = 146
    integer, parameter :: NAdmit         = 19
    
    !PREDICTION method 
    integer, parameter :: Task2000_ = 1, Toga_ = 2        

    type T_Size1D
        integer                 :: ILB            = null_int
        integer                 :: IUB            = null_int
    end type T_Size1D

#ifdef _USE_CUDA
    ! Bind this type to C if CUDA is used
    type, bind(c)   :: T_Size2D
#else
    type T_Size2D
#endif _USE_CUDA
        integer                 :: ILB            = null_int
        integer                 :: IUB            = null_int
        integer                 :: JLB            = null_int
        integer                 :: JUB            = null_int
    end type T_Size2D

#ifdef _USE_CUDA
    ! Bind this type to C if CUDA is used
    type, bind(c)   :: T_Size3D
#else
    type T_Size3D
#endif _USE_CUDA
        integer(C_INT)                 :: ILB            = null_int
        integer(C_INT)                 :: IUB            = null_int
        integer(C_INT)                 :: JLB            = null_int
        integer(C_INT)                 :: JUB            = null_int
        integer(C_INT)                 :: KLB            = null_int
        integer(C_INT)                 :: KUB            = null_int
    end type T_Size3D

    type T_PropertyID
        character(StringLength) :: Name              = null_str
        character(StringLength) :: Units             = null_str
        character(StringLength) :: Description       = null_str
        integer                 :: IDNumber          = null_int    
        integer                 :: ObjFillMatrix     = 0
        logical                 :: SolutionFromFile  = OFF
        logical                 :: IsAngle           = OFF
        logical                 :: IsParticulate     = OFF
        logical                 :: IsVectorial       = OFF
        logical                 :: IsDynamic         = OFF
    end type T_PropertyID

    type T_Instance
        integer                 :: ID         = 0
        integer                 :: Users      = 0
        integer                 :: Readers    = 0
        logical                 :: Read_Lock  = IDLE
    end type T_Instance

    type T_Module
        integer                 :: ID   = null_int  
        character(StringLength) :: Name = null_str 
    end type T_Module

    type (T_Module), dimension(MaxModules), parameter  :: MohidModules =  (/             &
        T_Module(mGLOBALDATA_            , "GlobalData"),            T_Module(mTIME_                   , "Time"),                  &
        T_Module(mENTERDATA_             , "EnterData"),             T_Module(mFUNCTIONS_              , "Functions"),             &
        T_Module(mLUD_                   , "LUD"),                   T_Module(mWATERQUALITY_           , "WaterQuality"),          &
        T_Module(mGRIDDATA_              , "GridData"),              T_Module(mHDF5_                   , "HDF5"),                  &
        T_Module(mHORIZONTALGRID_        , "HorizontalGrid"),        T_Module(mHORIZONTALMAP_          , "HorizontalMap"),         &
        T_Module(mGEOMETRY_              , "Geometry"),              T_Module(mSTOPWATCH_              , "StopWatch"),             &
        T_Module(mADVECTIONDIFFUSION_    , "AdvectionDiffusion"),    T_Module(mMAP_                    , "Map"),                   &
        T_Module(mTIMESERIE_             , "TimeSerie"),             T_Module(mSOIL_                   , "Soil"),                  &
        T_Module(mBOXDIF_                , "BoxDif"),                T_Module(mASSIMILATION_           , "Assimilation"),          &
        T_Module(mSTATISTIC_             , "Statistic"),             T_Module(mDISCHARGES_             , "Discharges"),            &
        T_Module(mATMOSPHERE_            , "Atmosphere"),            T_Module(mHYDROINTEGRATION_       , "HydroIntegration"),      &
        T_Module(mHYDRODYNAMICFILE_      , "HydrodynamicFile"),      T_Module(mTOGA_                   , "Toga"),                  &
        T_Module(mGAUGE_                 , "Gauge"),                 T_Module(mTRIANGULATION_          , "Triangulation"),         &
        T_Module(mGOTM_                  , "GOTM"),                  T_Module(mTURBGOTM_               , "TurbGOTM"),              &
        T_Module(mTURBULENCE_            , "Turbulence"),            T_Module(mOPENBOUNDARY_           , "OpenBoundary"),          &
        T_Module(mHYDRODYNAMIC_          , "Hydrodynamic"),          T_Module(mFREEVERTICALMOVEMENT_   , "FreeVerticalMovement"),  &
        T_Module(mSEDIMENTQUALITY_       , "SedimentQuality"),       T_Module(mINTERFACE_              , "Interface"),             &
        T_Module(mWATERPROPERTIES_       , "WaterProperties"),       T_Module(mOIL_                    , "Oil"),                   &
        T_Module(mJET_                   , "Jet"),                   T_Module(mLAGRANGIAN_             , "Lagrangian"),            &
        T_Module(mMODEL_                 , "Model"),                 T_Module(mCONSOLIDATION_          , "Consolidation"),         &
        T_Module(mSEDIMENTPROPERTIES_    , "SedimentProperties"),    T_Module(mLIFE_                   , "Life"),                  &
        T_Module(mCEQUALW2_              , "CEQUALW2"),              T_Module(mBASINGEOMETRY_          , "BasinGeometry"),         &
        T_Module(mRUNOFF_                , "RunOff"),                T_Module(mDRAINAGENETWORK_        , "DrainageNetwork"),       &
        T_Module(mINTERFACEWATERAIR_     , "InterfaceWaterAir"),     T_Module(mINTERFACESEDIMENTWATER_ , "InterfaceWaterSediment"),&
        T_Module(mLIGHTEXTINCTION_       , "LightExtinction"),       T_Module(mRIVERHYDRODYNAMIC_      , "RiverHydrodynamic"),     &
        T_Module(mFILLMATRIX_            , "FillMatrix"),            T_Module(mWAVES_                  , "Waves"),                 &
        T_Module(mBASIN_                 , "Basin"),                 T_Module(mSOILPROPERTIES_         , "SoilProperties"),        &
        T_Module(mINFILTRATION_          , "Infiltration"),          T_Module(mSOILPLANTAIR_           , "SoilPlantAir"),          &
        T_Module(mSOILMACROPORES_        , "SoilMacropores"),        T_Module(mSAND_                   , "Sand"),                  &
        T_Module(mMACROPOREPROPERTIES_   , "MacroporeProperties"),   T_Module(mPOROUSMEDIA_            , "PorousMedia"),           &
        T_Module(mSWAT_                  , "ModuleSwat"),            T_Module(mPROFILE_                , "Profile"),               &
        T_Module(mBENTHOS_               , "Benthos"),               T_Module(mCLIMATOLOGY_            , "Climatology"),           &
        T_Module(mFIREINDEX_             , "FireIndex"),             T_Module(mINTERPOLATION_          , "Interpolation"),         &
        T_Module(mVEGETATION_            , "Vegetation"),            T_Module(mRESERVOIROPTIMIZATION_  , "ReservoirOptimization"), &
        T_Module(mMACROALGAE_            , "MacroAlgae"),            T_Module(mBFM_                    , "BFM"),                   &
        T_Module(mNETCDF_                , "NETCDF"),                T_Module(mSEQUENTIALASSIMILATION_ , "SequentialAssimilation"),&
        T_Module(mPOROUSMEDIAPROPERTIES_ , "PorousMediaProperties"), T_Module(mPHREEQC_                , "PhreeqC"),               &
        T_Module(mCUDA_                  , "Cuda"),                                                                                &
        T_Module(mRUNOFFPROPERTIES_      , "RunoffProperties"),      T_Module(mCHAINREACTIONS_         , "ChainReactions"),        &
        T_Module(mField4D_               , "Field4D"),               T_Module(mBENTHICECOLOGY_         , "BenthicEcology"),        &
        T_Module(mWWTPQ_                 , "WWTPQ") ,                                                                              &
        T_Module(mSEAGRASSWATERINTERAC_  , "SeagrassWaterInteraction"),                                                            &
        T_Module(mSEAGRASSSEDIMINTERAC_  , "SeagrassSedimInteraction"), T_Module(mBivalve_             , "BivalveModel"),          &
        T_Module(mTimeSeriesAnalyser_    , "TimeSeriesAnalyser"      ), T_Module(mNetworkStatistics_   , "NetworkStatistics"),     &
        T_Module(mTimeSeriesOperator_    , "TimeSeriesOperator")     ,  T_Module(mAnalytical_LDS_      , "Analytical_LDS"),        &
        T_Module(mPressureDifferences_   , "PressureDifferences"),   T_Module(mHNS_                    , "HNS"           ),        &
        T_Module(mGlueWW3_OBC_           , "GlueWW3_OBC"),           T_Module(mSnow_                   , "Snow"          ),        &
        T_Module(mSediment_              , "Sediment"           ),   T_Module(mReservoirs_             , "Reservoirs"    ),        &
        T_Module(mIrrigation_            , "Irrigation"         ),   T_Module(mTURBINE_                , "Turbine"       ),        &
        T_Module(mLitter_                , "Litter"             ),   T_Module(mTwoWay_                 , "TwoWay"        ),        &
        T_Module(mOutputGrid_            , "OuputGrid"          )/)
        

    !Variables
    logical, dimension(MaxModules)                                  :: RegisteredModules  = .false.
    character(LEN=PathLength)                                       :: FilesName          = 'nomfich.dat'
    logical                                                         :: MonitorPerformance = .false.  
    logical                                                         :: MonitorDT          = .false.  
    integer                                                         :: UnitDT             = null_int 
    character(LEN=StringLength), dimension(:), pointer              :: PropNameList          => null()
    character(LEN=StringLength), dimension(:), pointer              :: DynamicPropNameList   => null()
    integer,                     dimension(:), pointer              :: PropNumberList        => null()
    integer,                     dimension(:), pointer              :: DynamicPropNumberList => null()
    integer                                                         :: PropertiesNumber        = 0 
    integer                                                         :: DynamicPropertiesNumber = 0 
    integer                                                         :: HighestPropertyID = 0
    integer, private                                                :: ErrorFileID     = 0
    integer, private                                                :: UsedKeyFileID   = 0
    integer, private                                                :: LogFileID       = 0
    character(LEN=1024)                                             :: OnLineString    = null_str 
    character(StringLength),     dimension(MaxErrorMessages)        :: ErrorMessagesStack
    
    !Used only when OpenMP is ON. 
    !Chunk(x)Factor is used to find the CHUNK size (number of rows/columns/layers)
    !that is passed to each thread. The Chunk(x) is the value of the chunk size and
    !is found by this equation:
    ! Chunk(X) = number_of_(rows or columns or layers) / Chunk(x)Factor
    integer                                                         :: openmp_num_threads = 1
    integer                                                         :: ChunkIFactor = 99999, &
                                                                       ChunkJFactor = 99999, &
                                                                       ChunkKFactor = 99999
    integer                                                         :: ChunkI = 1, &
                                                                       ChunkJ = 1, &
                                                                       ChunkK = 1
      
    type (T_Instance), dimension (MaxModules, MaxInstances), save   :: ObjCollector
    public :: ObjCollector

    contains 

    !--------------------------------------------------------------------------

    logical function ModuleIsRegistered (iModule)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
    
        ModuleIsRegistered = RegisteredModules(iModule)
        
    end function ModuleIsRegistered
    
    !--------------------------------------------------------------------------
    
    subroutine SetInputFullPath (Path)

        !Arguments------------------------------------
        character(LEN=*) :: Path

        !Begin------------------------------------

        if (trim(Path) /= '') then
            FilesName = trim(Path)
        else
            FilesName = 'nomfich.dat'
        endif

    end subroutine SetInputFullPath    
    
    !--------------------------------------------------------------------------
    
    subroutine RegisterModule (iModule)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        logical                                     :: OK

        RegisteredModules(iModule) = .true.

        OK = .true.

#ifdef _ONLINE_ 
        OK = .false.
#endif

#ifdef _GUI_
        OK = .false.
#endif

        !Feedback
        if (OK) write(unit=ErrorFileID, fmt=*)"Registered Module  : ", trim(MohidModules(iModule)%Name)

        
    end subroutine RegisterModule
    
    !--------------------------------------------------------------------------

    integer function RegisterNewInstance (iModule)

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule

        !Local-----------------------------------------------------------------
        integer                                     :: iInstance
        logical                                     :: OK

        !Search first empty ID
        iInstance = 1
        do while (ObjCollector(iModule, iInstance)%Users /= 0)
            iInstance = iInstance + 1
        enddo

        !Set users to 1
        ObjCollector(iModule, iInstance)%Users = 1

        OK = .true.

#ifdef _GUI_ 
        OK =.false.
#endif

#ifdef _ONLINE_ 
        OK =.false.
#endif
        if (OK) then
            !Feedback
            write(unit=ErrorFileID, fmt=*)"Registered Instance: ", trim(MohidModules(iModule)%Name)
            write(unit=ErrorFileID, fmt=*)"Instance ID        : ", iInstance
        endif

        !Returns ID
        RegisterNewInstance = iInstance

    end function RegisterNewInstance

    !--------------------------------------------------------------------------

    integer function AssociateInstance   (iModule, iInstance)

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        integer                                     :: iInstance

        !Local-----------------------------------------------------------------

        !Checks for error
        if (iInstance == 0) then
            write(*,*)'iInstance cannot be zero'
            stop 'AssociateInstance - ModuleGlobalData - ERR01'
        endif                    
       
        !Increase number by one
        ObjCollector(iModule, iInstance)%Users = ObjCollector(iModule, iInstance)%Users + 1

        !Returns ID
        AssociateInstance = iInstance
        
    end function AssociateInstance

    !--------------------------------------------------------------------------

    integer function DeassociateInstance (iModule, iInstance)

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        integer                                     :: iInstance

        !Local-----------------------------------------------------------------

        !Decrease number by one
        ObjCollector(iModule, iInstance)%Users = ObjCollector(iModule, iInstance)%Users - 1

        !Checks for error
        if (ObjCollector(iModule, iInstance)%Users < 0) then
            write(*,*)'Users cannot be negative'
            stop 'DeassociateInstance - ModuleGlobalData - ERR01'
        endif

        !Returns Number of users
        DeassociateInstance = ObjCollector(iModule, iInstance)%Users
        
    end function DeassociateInstance

    !--------------------------------------------------------------------------

    subroutine Read_Lock(iModule, iInstance) 

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        integer                                     :: iInstance

        !----------------------------------------------------------------------

        ObjCollector(iModule, iInstance)%Read_Lock = ACTIVE
        ObjCollector(iModule, iInstance)%readers   = ObjCollector(iModule, iInstance)%readers + 1

        !----------------------------------------------------------------------

    end subroutine Read_Lock

    !--------------------------------------------------------------------------

    subroutine Read_Unlock(iModule, iInstance, RoutineName) 

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        integer                                     :: iInstance
        character(len=*)                            :: RoutineName

        !----------------------------------------------------------------------
              
        !Cannot Read_unlock if Instance is not Read_lock
        if (.not. ObjCollector(iModule, iInstance)%Read_Lock) then
            write(*, *) 'Number of Readers error.'
            write(*, *) 'Routine Name    : ', trim(RoutineName)
            write(*, *) 'Module  Name    : ', trim(MohidModules(iModule)%Name)
            write(*, *) 'Instance  ID    : ', iInstance
            stop        'Read_Unlock - ModuleGlobalData - ERR01'
        end if
            
        !Decreases number of readers
        ObjCollector(iModule, iInstance)%readers = ObjCollector(iModule, iInstance)%readers - 1

        !If number of reades equal to zero set Read_lock to IDLE
        if (ObjCollector(iModule, iInstance)%readers == 0) ObjCollector(iModule, iInstance)%Read_Lock = IDLE

        !if Number of readers is negative, somethink is wrong
        if (ObjCollector(iModule, iInstance)%readers < 0) then
            write(*, *) 'Negative number of readers'
            write(*, *) 'Routine Name    : ', trim(RoutineName)
            write(*, *) 'Module  Name    : ', trim(MohidModules(iModule)%Name)
            write(*, *) 'Instance  ID    : ', iInstance
            stop        'Read_Unlock - ModuleGlobalData - ERR02'
        end if

        !----------------------------------------------------------------------

    end subroutine Read_Unlock

    !--------------------------------------------------------------------------

    integer function VerifyReadLock (iModule, iInstance) 

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        integer                                     :: iInstance

        if  (ObjCollector(iModule, iInstance)%Read_Lock) then 
            VerifyReadLock = READ_LOCK_ERR_
        else 
            VerifyReadLock = IDLE_ERR_
        end if

    end function VerifyReadLock

    !--------------------------------------------------------------------------

    logical function CheckPropertyName (PropertyName, Number, IsDynamic)

        !Arguments-------------------------------------------------------------
        character(len=*), intent (IN)               :: PropertyName
        integer,          intent (OUT), optional    :: Number
        logical,          intent (OUT), optional    :: IsDynamic

        !Local-----------------------------------------------------------------
        integer :: i

        !----------------------------------------------------------------------

        CheckPropertyName = .false.
        if (present(IsDynamic)) IsDynamic = .false.

        call ConstructPropList

        do i=1, PropertiesNumber

            if (PropertyName == PropNameList(i)) then

                if (present(Number)) Number = PropNumberList(i)
                CheckPropertyName = .TRUE.

            endif

        enddo
        
        if (.not. CheckPropertyName) then
            
            if(associated(DynamicPropNameList)) then            

                do i=1, DynamicPropertiesNumber

                    if (PropertyName == DynamicPropNameList(i)) then

                        if (present(Number)) Number = DynamicPropNumberList(i)
                        CheckPropertyName = .TRUE.
                        if (present(IsDynamic)) IsDynamic = .true.

                    endif

                enddo
        
            endif            
        endif

        !----------------------------------------------------------------------

    end function  CheckPropertyName

!~     !--------------------------------------------------------------------------
     
!~     logical function CheckDynamicPropertyName (PropertyName, Number)

!~         !Arguments-------------------------------------------------------------
!~         character(len=*), intent (IN)               :: PropertyName
!~         integer,          intent (OUT), optional    :: Number

!~         !Local-----------------------------------------------------------------
!~         integer :: i

!~         !----------------------------------------------------------------------

!~         CheckDynamicPropertyName = .false.

!~         if(associated(DynamicPropNameList)) then            

!~             do i=1, DynamicPropertiesNumber

!~                 if (PropertyName == DynamicPropNameList(i)) then

!~                     if (present(Number)) Number = DynamicPropNumberList(i)
!~                     CheckDynamicPropertyName = .TRUE.

!~                 endif

!~             enddo
        
!~         endif

!~         !----------------------------------------------------------------------

!~     end function  CheckDynamicPropertyName

    !--------------------------------------------------------------------------
    
    integer function RegisterDynamicProperty(PropertyName) !soffs
      
        !Arguments-------------------------------------------------------------
        character(len=*), intent (IN)                       :: PropertyName
        
        !Local-----------------------------------------------------------------
        integer, save                                       :: UniqueIDNumber
        character(LEN=StringLength), dimension(:), pointer  :: TempPropNameList
        integer,                     dimension(:), pointer  :: TempPropNumberList
        logical                                             :: is_dynamic
        integer                                             :: i
        logical                                             :: property_exists
        
        !----------------------------------------------------------------------
        property_exists = .false.
        
        if (CheckPropertyName(PropertyName, IsDynamic = is_dynamic)) then
            if (is_dynamic) then
do1:            do i=1, DynamicPropertiesNumber
                    if (trim(PropertyName) == trim(DynamicPropNameList(i))) then
                        RegisterDynamicProperty = DynamicPropNumberList(i)
                        property_exists = .true.
                        exit do1
                    endif
                enddo do1
            else
                print *, "Dynamic Property '"//trim(PropertyName)//"' already exists in GLOBAL DATA."
                stop "RegisterDynamicProperty - ModuleGlobalData - ERR010"
            endif
        endif
        
        if(.not. property_exists) then
            if(.not. associated(DynamicPropNameList))then
        
                DynamicPropertiesNumber = 1
                allocate(DynamicPropNameList  (DynamicPropertiesNumber))
                allocate(DynamicPropNumberList(DynamicPropertiesNumber))
                UniqueIDNumber         = HighestPropertyID + 1

                DynamicPropNameList(1)      = trim(PropertyName)
                DynamicPropNumberList(1)    = UniqueIDNumber
            
            else            
                    allocate(TempPropNameList     (DynamicPropertiesNumber))
                    allocate(TempPropNumberList   (DynamicPropertiesNumber))
            
                    TempPropNameList    = DynamicPropNameList
                    TempPropNumberList  = DynamicPropNumberList
            
                    deallocate(DynamicPropNameList, DynamicPropNumberList)
            
                    UniqueIDNumber          = UniqueIDNumber + 1
                    DynamicPropertiesNumber = DynamicPropertiesNumber + 1

                    allocate(DynamicPropNameList  (DynamicPropertiesNumber))
                    allocate(DynamicPropNumberList(DynamicPropertiesNumber))
            
                    DynamicPropNameList  (1:DynamicPropertiesNumber-1) = TempPropNameList  (1:DynamicPropertiesNumber-1)
                    DynamicPropNumberList(1:DynamicPropertiesNumber-1) = TempPropNumberList(1:DynamicPropertiesNumber-1)
           
                    DynamicPropNameList  (DynamicPropertiesNumber) = trim(PropertyName)
                    DynamicPropNumberList(DynamicPropertiesNumber) = UniqueIDNumber
            
                    deallocate(TempPropNameList, TempPropNumberList)

            endif 
            
            RegisterDynamicProperty = UniqueIDNumber 
            
        endif

    end function

    !----------------------------------------------------------------------
    
    integer function GetDynamicPropertyIDNumber (DynamicPropertyName) !soffs

        !Arguments-------------------------------------------------------------
        character(len=*), intent (IN )              :: DynamicPropertyName

        !Local-----------------------------------------------------------------
        integer :: i

        !----------------------------------------------------------------------
        
        GetDynamicPropertyIDNumber = UNKNOWN_

        do i=1, DynamicPropertiesNumber

            if (trim(DynamicPropertyName) == trim(DynamicPropNameList(i))) then

                GetDynamicPropertyIDNumber = DynamicPropNumberList(i)

            endif

        enddo

        if(GetDynamicPropertyIDNumber == UNKNOWN_)then
            write(*,*)'Unknown property'
            stop 'GetDynamicPropertyIDNumber - ModuleGlobalData - ERR01'
        end if

    !----------------------------------------------------------------------

    end function GetDynamicPropertyIDNumber

    !--------------------------------------------------------------------------

    character (len=StringLength) function GetPropertyName (Number)

        !Arguments-------------------------------------------------------------
        integer,          intent (IN )              :: Number

        !Local-----------------------------------------------------------------
        integer :: i
        logical :: found

        !----------------------------------------------------------------------

        call ConstructPropList

        found = .false.
        
do1:    do i=1, PropertiesNumber

            if (Number == PropNumberList(i)) then

                GetPropertyName = PropNameList(i)
                found = .true.
                exit do1

            endif

        enddo do1
        
        if (.not. found) then
            
            if(associated(DynamicPropNameList)) then            

do2:            do i=1, DynamicPropertiesNumber

                    if (Number == DynamicPropNumberList(i)) then

                        GetPropertyName = DynamicPropNameList(i)
                        exit do2

                    endif

                enddo do2
        
            endif            
        endif

        !----------------------------------------------------------------------

    end function GetPropertyName

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    integer function GetPropertyIDNumber (PropertyName, StopActive)

        !Arguments-------------------------------------------------------------
        character(len=*), intent (IN )              :: PropertyName
        logical, optional                           :: StopActive

        !Local-----------------------------------------------------------------
        integer                                     :: i
        logical                                     :: StopActive_

        !----------------------------------------------------------------------
        
        if (present(StopActive)) then
            StopActive_ = StopActive
        else
            StopActive_ = .true. 
        endif
        
        GetPropertyIDNumber = UNKNOWN_

        call ConstructPropList

do1:    do i=1, PropertiesNumber

            if (trim(PropertyName) == trim(PropNameList(i))) then

                GetPropertyIDNumber = PropNumberList(i)
                exit do1

            endif

        enddo do1

        if (GetPropertyIDNumber == UNKNOWN_) then
            
            if(associated(DynamicPropNameList)) then            

do2:            do i=1, DynamicPropertiesNumber

                    if (trim(PropertyName) == trim(DynamicPropNameList(i))) then

                        GetPropertyIDNumber = DynamicPropNumberList(i)
                        exit do2

                    endif

                enddo do2
        
            endif            
        endif
        
        if(GetPropertyIDNumber == UNKNOWN_)then
            write(*,*)'Unknown property: ', PropertyName
            if (StopActive_) then
                stop 'GetPropertyIDNumber - ModuleGlobalData - ERR010'
            endif                
        end if

        !----------------------------------------------------------------------

    end function GetPropertyIDNumber

    !--------------------------------------------------------------------------

    subroutine ConstructPropList

        !Local-----------------------------------------------------------------
        integer       :: ListNumber, i, j

        !Begin-----------------------------------------------------------------
        logical, save :: FirstList = .true.


        
        if (FirstList) then

            FirstList = .false.

            call AddPropList (Temperature_,             Char_Temperature,               ListNumber)
            call AddPropList (Salinity_,                Char_Salinity,                  ListNumber)   
            call AddPropList (Density_,                 Char_Density,                   ListNumber)   
            call AddPropList (Phytoplankton_,           Char_Phytoplankton,             ListNumber)
            call AddPropList (Zooplankton_,             Char_Zooplankton,               ListNumber)
            call AddPropList (DOPRefractory_,           Char_DOPRefractory,             ListNumber)
            call AddPropList (DOPNon_Refractory_,       Char_DOPNon_Refractory,         ListNumber)
            call AddPropList (DONRefractory_,           Char_DONRefractory,             ListNumber)
            call AddPropList (DONNon_Refractory_,       Char_DONNon_Refractory,         ListNumber)
            call AddPropList (Inorganic_Phosphorus_,    Char_Inorganic_Phosphorus,      ListNumber)         
            
            call AddPropList (POC_,                     Char_POC,                       ListNumber)                             
            call AddPropList (POP_,                     Char_POP,                       ListNumber)                      
            call AddPropList (PON_,                     Char_PON,                       ListNumber)                      
            call AddPropList (PONRefractory_,           Char_PONRefractory,             ListNumber)                      
            call AddPropList (DOC_,                     Char_DOC,                       ListNumber)                      
            call AddPropList (DOP_,                     Char_DOP,                       ListNumber)                      
            call AddPropList (DON_,                     Char_DON,                       ListNumber)                      
            call AddPropList (DOCsl_,                   Char_DOCsl,                     ListNumber)                    
            call AddPropList (DOPsl_,                   Char_DOPsl,                     ListNumber)                    
            call AddPropList (DONsl_,                   Char_DONsl,                     ListNumber)                    
            call AddPropList (Ammonia_,                 Char_Ammonia,                   ListNumber)                  
            call AddPropList (Nitrate_,                 Char_Nitrate,                   ListNumber)                  
            call AddPropList (Silicate_,                Char_Silicate,                  ListNumber)                 
            call AddPropList (BioSilica_,               Char_BioSilica,                 ListNumber)               
            call AddPropList (CarbonDioxide_,           Char_CarbonDioxide,             ListNumber)                      
            call AddPropList (Oxygen_,                  Char_Oxygen,                    ListNumber) 
            call AddPropList (DissolO2PercentSat_,      Char_DissolO2PercentSat,        ListNumber)
            call AddPropList (CO2PartialPressure_,      Char_CO2PartialPressure,        ListNumber)
            call AddPropList (PhytoChla_,               Char_PhytoChla,                 ListNumber)
            call AddPropList (Diatom_C_,                Char_Diatom_C,                  ListNumber)     
            call AddPropList (Diatom_N_,                Char_Diatom_N,                  ListNumber)                 
            call AddPropList (Diatom_P_,                Char_Diatom_P,                  ListNumber)                 
            call AddPropList (Diatom_Si_,               Char_Diatom_Si,                 ListNumber)               
            call AddPropList (Diatom_Chl_,              Char_Diatom_Chl,                ListNumber)               
            call AddPropList (Mix_Flagellate_C_,        Char_Mix_Flagellate_C,          ListNumber)               
            call AddPropList (Mix_Flagellate_N_,        Char_Mix_Flagellate_N,          ListNumber)    
            call AddPropList (Mix_Flagellate_P_,        Char_Mix_Flagellate_P,          ListNumber)         
            call AddPropList (Mix_Flagellate_Chl_,      Char_Mix_Flagellate_Chl,        ListNumber)         
            call AddPropList (Picoalgae_C_,             Char_Picoalgae_C,               ListNumber)          
            call AddPropList (Picoalgae_N_,             Char_Picoalgae_N,               ListNumber)          
            call AddPropList (Picoalgae_P_,             Char_Picoalgae_P,               ListNumber)     
            call AddPropList (Picoalgae_Chl_,           Char_Picoalgae_Chl,             ListNumber)      
            call AddPropList (Flagellate_C_,            Char_Flagellate_C,              ListNumber)              
            call AddPropList (Flagellate_N_,            Char_Flagellate_N,              ListNumber)              
            call AddPropList (Flagellate_P_,            Char_Flagellate_P,              ListNumber)              
            call AddPropList (Flagellate_Chl_,          Char_Flagellate_Chl,            ListNumber)              
            call AddPropList (Microzooplankton_C_,      Char_Microzooplankton_C,        ListNumber)     
            call AddPropList (Microzooplankton_N_,      Char_Microzooplankton_N,        ListNumber)             
            call AddPropList (Microzooplankton_P_,      Char_Microzooplankton_P,        ListNumber)             
            call AddPropList (Het_Nanoflagellate_C_,    Char_Het_Nanoflagellate_C,      ListNumber)         
            call AddPropList (Het_Nanoflagellate_N_,    Char_Het_Nanoflagellate_N,      ListNumber) 
            call AddPropList (Het_Nanoflagellate_P_,    Char_Het_Nanoflagellate_P,      ListNumber)  
            call AddPropList (Mesozooplankton_C_,       Char_Mesozooplankton_C,         ListNumber)    
            call AddPropList (Mesozooplankton_N_,       Char_Mesozooplankton_N,         ListNumber)                
            call AddPropList (Mesozooplankton_P_,       Char_Mesozooplankton_P,         ListNumber)               
            call AddPropList (Het_Bacteria_C_,          Char_Het_Bacteria_C,            ListNumber)                 
            call AddPropList (Het_Bacteria_N_,          Char_Het_Bacteria_N,            ListNumber)                 
            call AddPropList (Het_Bacteria_P_,          Char_Het_Bacteria_P,            ListNumber)                          
         
            call AddPropList (Bivalve1_,                Char_Bivalve1,                  ListNumber) !soffs
            call AddPropList (Bivalve2_,                Char_Bivalve2,                  ListNumber)
            call AddPropList (Bivalve3_,                Char_Bivalve3,                  ListNumber)
            call AddPropList (Bivalve4_,                Char_Bivalve4,                  ListNumber)

            call AddPropList (Shrimp_,                  Char_Shrimp,                    ListNumber) !soffs
            call AddPropList (Crab_,                    Char_Crab,                      ListNumber)
            call AddPropList (OysterCatcher_,           Char_OysterCatcher,             ListNumber)
            call AddPropList (EiderDuck_,               Char_EiderDuck,                 ListNumber)
            call AddPropList (HerringGull_,             Char_HerringGull,               ListNumber)

            call AddPropList (Nitrite_,                 Char_Nitrite,                   ListNumber)
            call AddPropList (BOD_,                     Char_BOD,                       ListNumber)
            call AddPropList (Cohesive_Sediment_,       Char_Cohesive_Sediment,         ListNumber)
            call AddPropList (Fecal_Coliforms_,         Char_Fecal_Coliforms,           ListNumber)
            call AddPropList (E_Coli_,                  Char_E_Coli,                    ListNumber)
            call AddPropList (T90_,                     Char_T90,                       ListNumber)
            call AddPropList (T90_E_Coli_,              Char_T90_E_Coli,                ListNumber)
            call AddPropList (Oil_,                     Char_Oil,                       ListNumber)
            call AddPropList (FloatingObject_,          Char_FloatingObject,            ListNumber)
            call AddPropList (HumanBody_,               Char_HumanBody,                 ListNumber)
            call AddPropList (OilThickness_,            Char_OilThickness,              ListNumber)            
            call AddPropList (HNS_,                     Char_HNS,                       ListNumber)
            call AddPropList (Ciliate_,                 Char_Ciliate,                   ListNumber)
            call AddPropList (Bacteria_,                Char_Bacteria,                  ListNumber)
            call AddPropList (ParticulateArsenic_,      Char_ParticulateArsenic,        ListNumber)
            call AddPropList (DissolvedArsenic_,        Char_DissolvedArsenic,          ListNumber)
            call AddPropList (ParticulateZinc_,         Char_ParticulateZinc,           ListNumber)
            call AddPropList (DissolvedZinc_,           Char_DissolvedZinc,             ListNumber)
            call AddPropList (ParticulateCopper_,       Char_ParticulateCopper,         ListNumber)
            call AddPropList (DissolvedCopper_,         Char_DissolvedCopper,           ListNumber)
            call AddPropList (ParticulateCadmium_,      Char_ParticulateCadmium,        ListNumber)
            call AddPropList (DissolvedCadmium_,        Char_DissolvedCadmium,          ListNumber)
            call AddPropList (ParticulateMercury_,      Char_ParticulateMercury,        ListNumber)
            call AddPropList (DissolvedMercury_,        Char_DissolvedMercury,          ListNumber)
            call AddPropList (ParticulateLead_,         Char_ParticulateLead,           ListNumber)
            call AddPropList (DissolvedLead_,           Char_DissolvedLead,             ListNumber)
            call AddPropList (ParticulateMetal_,        Char_ParticulateMetal,          ListNumber)
            call AddPropList (DissolvedMetal_,          Char_DissolvedMetal,            ListNumber)
            
            call AddPropList (Larvae_,                  Char_Larvae,                    ListNumber)
            call AddPropList (Age_,                     Char_Age,                       ListNumber)
            call AddPropList (Fish_,                    Char_Fish,                      ListNumber)
            call AddPropList (FishFood_,                Char_FishFood,                  ListNumber)
            call AddPropList (MacroAlgae_,              Char_MacroAlgae,                ListNumber)
            call AddPropList (DriftingMacroAlgae_,      Char_DriftingMacroAlgae,        ListNumber)
            call AddPropList (MicroPhytoBenthos_,       Char_MicroPhytoBenthos,         ListNumber)
            
            call AddPropList (SuspensionFeedersN_,      Char_SuspensionFeedersN,        ListNumber) ! isab
            call AddPropList (SuspensionFeedersC_,      Char_SuspensionFeedersC,        ListNumber)
            call AddPropList (SuspensionFeedersP_,      Char_SuspensionFeedersP,        ListNumber)
            call AddPropList (MicroPhytoBenthosN_,      Char_MicroPhytoBenthosN,        ListNumber) ! isab
            call AddPropList (MicroPhytoBenthosC_,      Char_MicroPhytoBenthosC,        ListNumber)
            call AddPropList (MicroPhytoBenthosP_,      Char_MicroPhytoBenthosP,        ListNumber)
            call AddPropList (DepositFeedersC_,         Char_DepositFeedersC,           ListNumber)
            call AddPropList (DepositFeedersN_,         Char_DepositFeedersN,           ListNumber)
            call AddPropList (DepositFeedersP_,         Char_DepositFeedersP,           ListNumber)
            call AddPropList (SeagrassesN_,             Char_SeagrassesN,               ListNumber)
            call AddPropList (SeagrassesP_,             Char_SeagrassesP,               ListNumber)
            call AddPropList (SeagrassesRoots_,         Char_SeagrassesRoots,           ListNumber)
            call AddPropList (SeagrassesLeaves_,        Char_SeagrassesLeaves,          ListNumber)

            call AddPropList (AdsorbedAmmonia_,         Char_AdsorbedAmmonia,           ListNumber)
            call AddPropList (RefreactaryOrganicN_,     Char_RefreactaryOrganicN,       ListNumber)
            call AddPropList (Ngas_,                    Char_Ngas,                      ListNumber)
            call AddPropList (AmmoniaGas_,              Char_AmmoniaGas,                ListNumber)
            call AddPropList (Urea_,                    Char_Urea,                      ListNumber)
            call AddPropList (HeterotrophicN_,          Char_HeterotrophicN,            ListNumber)
            call AddPropList (AnaerobicN_,              Char_AnaerobicN,                ListNumber)
            call AddPropList (AutotrophicN_,            Char_AutotrophicN,              ListNumber)
            call AddPropList (SolubilizingN_,           Char_SolubilizingN,             ListNumber)
            call AddPropList (AnaerobicC_,              Char_AnaerobicC,                ListNumber)
            call AddPropList (AutotrophicC_,            Char_AutotrophicC,              ListNumber)
            call AddPropList (HeterotrophicC_,          Char_HeterotrophicC,            ListNumber)
            call AddPropList (SolubilizingC_,           Char_SolubilizingC,             ListNumber)
            call AddPropList (LabileOrganicC_,          Char_LabileOrganicC,            ListNumber)
            call AddPropList (RefreactaryOrganicC_,     Char_RefreactaryOrganicC,       ListNumber)           
            call AddPropList (Methane_,                 Char_Methane,                   ListNumber)           
            call AddPropList (HeterotrophicP_,          Char_HeterotrophicP,            ListNumber)
            call AddPropList (AnaerobicP_,              Char_AnaerobicP,                ListNumber)
            call AddPropList (AutotrophicP_,            Char_AutotrophicP,              ListNumber)
            call AddPropList (SolubilizingP_,           Char_SolubilizingP,             ListNumber)
            call AddPropList (RefreactaryOrganicP_,     Char_RefreactaryOrganicP,       ListNumber)
            call AddPropList (AdsorbedInorganicP_,      Char_AdsorbedInorganicP,        ListNumber)
            call AddPropList (SoilDryDensity_,          Char_SoilDryDensity,            ListNumber)
            call AddPropList (IonicStrength_,           Char_IonicStrength,             ListNumber)
            call AddPropList (PhosphorusAdsortionIndex_,Char_PhosphorusAdsortionIndex,  ListNumber)
            call AddPropList (AutotrophicPop_,          Char_AutotrophicPop,            ListNumber)
            call AddPropList (HeterotrophicPop_,        Char_HeterotrophicPop,          ListNumber)
            call AddPropList (AnaerobicPop_,            Char_AnaerobicPop,              ListNumber)
            call AddPropList (SolPop_,                  Char_SolPop,                    ListNumber)
           
            call AddPropList (GenericProperty_,         Char_GenericProperty,           ListNumber)
            call AddPropList (GrossProd_,               Char_GrossProd,                 ListNumber)
            call AddPropList (NetProd_,                 Char_NetProd,                   ListNumber)
            call AddPropList (Excretion_,               Char_Excretion,                 ListNumber)
            call AddPropList (Respiration_,             Char_Respiration,               ListNumber)
            call AddPropList (NaturalMort_,             Char_NaturalMort,               ListNumber)
            call AddPropList (Grazing_,                 Char_Grazing,                   ListNumber)
            call AddPropList (MACondition_,             Char_MACondition,               ListNumber)
            call AddPropList (NutrientLim_,             Char_NutrientLim,               ListNumber)
            call AddPropList (NLim_,                    Char_NLim,                      ListNumber)
            call AddPropList (PLim_,                    Char_PLim,                      ListNumber)
            call AddPropList (LightLim_ ,               Char_LightLim ,                 ListNumber)
            call AddPropList (TemperatureLim_,          Char_TemperatureLim,            ListNumber)
            call AddPropList (SalinityLim_,             Char_SalinityLim,               ListNumber)
            call AddPropList (CarrCapLim_,              Char_CarrCapLim,                ListNumber)
            call AddPropList (DiaGrossProd_,            Char_DiaGrossProd,              ListNumber)
            call AddPropList (DiaNutrientLim_,          Char_DiaNutrientLim,            ListNumber)
            call AddPropList (DiaNLim_,                 Char_DiaNLim,                   ListNumber)
            call AddPropList (DiaPLim_,                 Char_DiaPLim,                   ListNumber)
            call AddPropList (DiaSiLim_,                Char_DiaSiLim,                  ListNumber)
            call AddPropList (DiaLightLim_ ,            Char_DiaLightLim ,              ListNumber)
            call AddPropList (DiaTemperatureLim_,       Char_DiaTemperatureLim,         ListNumber)
            call AddPropList (Diatoms_,                 Char_Diatoms,                   ListNumber)
            call AddPropList (WaterLevel_,              Char_WaterLevel_,               ListNumber)
            !call AddPropList (WaterLevelMax_,           Char_WaterLevelMax_,            ListNumber)
            !call AddPropList (WaterLevelMin_,           Char_WaterLevelMin_,            ListNumber)
            call AddPropList (VelocityModulus_,         Char_VelocityModulus_,          ListNumber)
            call AddPropList (VelocityDirection_,       Char_VelocityDirection_,        ListNumber)
            call AddPropList (FlowModulus_,             Char_FlowModulus_,              ListNumber)
            call AddPropList (VelocityU_,               Char_VelocityU_,                ListNumber)
            call AddPropList (VelocityV_,               Char_VelocityV_,                ListNumber)
            call AddPropList (VelocityW_,               Char_VelocityW_,                ListNumber)
            call AddPropList (WaterFluxX_,              Char_WaterFluxX_,               ListNumber)
            call AddPropList (WaterFluxY_,              Char_WaterFluxY_,               ListNumber)
            call AddPropList (WaterDepth_,              Char_WaterDepth_,               ListNumber)
            call AddPropList (CoriolisX_,               Char_CoriolisX_,                ListNumber)
            call AddPropList (BaroclinicForceX_,        Char_BaroclinicForceX_,         ListNumber)
            call AddPropList (HorizontalTransportX_,    Char_HorizontalTransportX_,     ListNumber)
            call AddPropList (CoriolisY_,               Char_CoriolisY_,                ListNumber)
            call AddPropList (BaroclinicForceY_  ,      Char_BaroclinicForceY_  ,       ListNumber)
            call AddPropList (HorizontalTransportY_,    Char_HorizontalTransportY_,     ListNumber)
            call AddPropList (BarotropicVelocityU_ ,    Char_BarotropicVelocityU_ ,     ListNumber)
            call AddPropList (BarotropicVelocityV_ ,    Char_BarotropicVelocityV_ ,     ListNumber)
            call AddPropList (WaterColumn_,             Char_WaterColumn_,              ListNumber)
            call AddPropList (MeridionalVelocity_,      Char_MeridionalVelocity_,       ListNumber)
            call AddPropList (ZonalVelocity_,           Char_ZonalVelocity_,            ListNumber)
            call AddPropList (TideState_,               Char_TideState_,                ListNumber)            
           
            !seagrasses rates and limiting functions
            call AddPropList (LeavesUptakeN_ ,          Char_LeavesUptakeN ,            ListNumber)  !Isabella
            call AddPropList (LeavesUptakeP_ ,          Char_LeavesUptakeP ,            ListNumber)
            call AddPropList (LeavesLightFactor_ ,      Char_LeavesLightFactor ,        ListNumber)
            call AddPropList (RootsUptakeN_ ,           Char_RootsUptakeN ,             ListNumber)
            call AddPropList (RootsUptakeP_ ,           Char_RootsUptakeP ,             ListNumber)
            call AddPropList (NintFactor_ ,             Char_NintFactor ,               ListNumber)
            call AddPropList (NintFactorR_ ,            Char_NintFactorR ,              ListNumber)
            call AddPropList (PintFactor_ ,             Char_PintFactor ,               ListNumber)
            call AddPropList (PintFactorR_ ,            Char_PintFactorR ,              ListNumber)
            call AddPropList (RootsMort_ ,              Char_RootsMort ,                ListNumber)

!____POM pools (for aquaculture cages)__________________________________________    
    
            call AddPropList (PON1_,                    Char_PON1,                      ListNumber) 
            call AddPropList (PON2_,                    Char_PON2,                      ListNumber)             
            call AddPropList (PON3_,                    Char_PON3,                      ListNumber)             
            call AddPropList (PON4_,                    Char_PON4,                      ListNumber) 
            call AddPropList (PON5_,                    Char_PON5,                      ListNumber) 
            
            call AddPropList (POP1_,                    Char_POP1,                      ListNumber) 
            call AddPropList (POP2_,                    Char_POP2,                      ListNumber)             
            call AddPropList (POP3_,                    Char_POP3,                      ListNumber)             
            call AddPropList (POP4_,                    Char_POP4,                      ListNumber) 
            call AddPropList (POP5_,                    Char_POP5,                      ListNumber)             
            
            
            
            ! guillaume nogueira
            call AddPropList (AltimLevelAnalyzed_,         Char_AltimLevelAnalyzed_,          ListNumber)
            call AddPropList (AltimTemperatureAnalyzed_,   Char_AltimTemperatureAnalyzed_,    ListNumber)
            call AddPropList (AltimSalinityAnalyzed_,      Char_AltimSalinityAnalyzed_,       ListNumber)
            call AddPropList (AltimLevelToAssimilate_ ,    Char_AltimLevelToAssimilate_ ,     ListNumber)
            call AddPropList (VarianceFieldToAssimilate_ , Char_VarianceFieldToAssimilate_ ,  ListNumber)
            
            call AddPropList (BaroclinicVelocityU_ ,    Char_BaroclinicVelocityU_ ,     ListNumber)
            call AddPropList (BaroclinicVelocityV_ ,    Char_BaroclinicVelocityV_ ,     ListNumber)
            call AddPropList (ObstacleDragCoef_    ,    Char_ObstacleDragCoef     ,     ListNumber)
            
            call AddPropList (ScraperVelU_         ,    Char_ScraperVelU          ,     ListNumber)
            call AddPropList (ScraperVelV_         ,    Char_ScraperVelV          ,     ListNumber)            
            call AddPropList (ScraperVelW_         ,    Char_ScraperVelW          ,     ListNumber)                        
            

            call AddPropList (Vorticity_ ,              Char_Vorticity ,                ListNumber)
            call AddPropList (BaroclinicKE_ ,           Char_BaroclinicKE ,             ListNumber)
            call AddPropList (PerturbationPE_ ,         Char_PerturbationPE ,           ListNumber)
            call AddPropList (KineticEnergy_ ,          Char_KineticEnergy ,            ListNumber)

            call AddPropList (ShearVelocity_,           Char_ShearVelocity,             ListNumber)

            call AddPropList (VerticalZ_,               Char_VerticalZ,                 ListNumber)

            call AddPropList (ParticulateContaminant_,  Char_ParticulateContaminant,    ListNumber)
            call AddPropList (DissolvedContaminant_,    Char_DissolvedContaminant,      ListNumber)

            call AddPropList (Sediment,                 Char_Sediment,                  ListNumber)

            call AddPropList (DissolvedSodium_ ,        Char_DissolvedSodium ,          ListNumber)
            call AddPropList (DissolvedCalcium_,        Char_DissolvedCalcium,          ListNumber)
                                                        
            call AddPropList (ParticulateSodium_ ,      Char_ParticulateSodium ,        ListNumber)
            call AddPropList (ParticulateCalcium_,      Char_ParticulateCalcium,        ListNumber)

            call AddPropList (LatentHeat_,              Char_LatentHeat,                ListNumber)
            call AddPropList (SensibleHeat_,            Char_SensibleHeat,              ListNumber)
            call AddPropList (Evaporation_  ,           Char_Evaporation  ,             ListNumber)
            call AddPropList (NetLongWaveRadiation_,    Char_NetLongWaveRadiation,      ListNumber)
            call AddPropList (UpwardLongWaveRadiation_, Char_UpwardLongWaveRadiation,   ListNumber)
            call AddPropList (DownwardLongWaveRadiation_,Char_DownwardLongWaveRadiation,ListNumber)

            call AddPropList (ShortWaveSolarRadiation_, Char_ShortWaveSolarRadiation,   ListNumber)
            call AddPropList (LongWaveSolarRadiation_,  Char_LongWaveSolarRadiation,    ListNumber)

            call AddPropList (ShortWaveSolarRadiationExtin_, Char_ShortWaveSolarRadiaExtin,   ListNumber)
            call AddPropList (LongWaveSolarRadiationExtin_ , Char_LongWaveSolarRadiaExtin,    ListNumber)
            
            call AddPropList (OxygenFlux_ ,                Char_OxygenFlux ,               ListNumber)
            call AddPropList (SpecificOxygenFlux_ ,        Char_SpecificOxygenFlux ,       ListNumber)         
            call AddPropList (WindShearVelocity_ ,         Char_WindShearVelocity ,        ListNumber)
            call AddPropList (SurfaceRadiation_ ,          Char_SurfaceRadiation ,         ListNumber)
            call AddPropList (WindStressX_,                Char_WindStressX,               ListNumber)
            call AddPropList (WindStressY_,                Char_WindStressY,               ListNumber)
            call AddPropList (WindStress_,                 Char_WindStress,                ListNumber)
            call AddPropList (SurfaceWaterFlux_ ,          Char_SurfaceWaterFlux ,         ListNumber)
            call AddPropList (NonSolarFlux_ ,              Char_NonSolarFlux ,             ListNumber)
            call AddPropList (TurbulentKineticEnergy_,     Char_TurbulentKineticEnergy,    ListNumber)
            call AddPropList (CarbonDioxideFlux_ ,         Char_CarbonDioxideFlux,         ListNumber)
            call AddPropList (SpecificCarbonDioxideFlux_ , Char_SpecificCarbonDioxideFlux, ListNumber)
            call AddPropList (Albedo_ ,                    Char_Albedo,                    ListNumber)
            call AddPropList (AmmoniaFlux_ ,               Char_AmmoniaFlux,               ListNumber)
            call AddPropList (NitrateFlux_ ,               Char_NitrateFlux,               ListNumber)


            call AddPropList (WindVelocityX_,           Char_WindVelocityX          ,      ListNumber)
            call AddPropList (WindVelocityY_,           Char_WindVelocityY          ,      ListNumber)
            call AddPropList (WindVelocity_,            Char_WindVelocity           ,      ListNumber)
            call AddPropList (SolarRadiation_,          Char_SolarRadiation         ,      ListNumber)
            call AddPropList (Precipitation_,           Char_Precipitation          ,      ListNumber)
            call AddPropList (Reflectivity_,            Char_Reflectivity           ,      ListNumber)            
            call AddPropList (AtmosphericPressure_,     Char_AtmosphericPressure    ,      ListNumber)
            call AddPropList (CO2AtmosphericPressure_,  Char_CO2AtmosphericPressure ,      ListNumber)
            call AddPropList (O2AtmosphericPressure_,   Char_O2AtmosphericPressure  ,      ListNumber)
            call AddPropList (AirTemperature_,          Char_AirTemperature         ,      ListNumber)
            call AddPropList (RelativeHumidity_,        Char_RelativeHumidity       ,      ListNumber)
            call AddPropList (SpecificHumidity_,        Char_SpecificHumidity       ,      ListNumber)                       
            call AddPropList (WindModulos_,             Char_WindModulos            ,      ListNumber)
            call AddPropList (WindAngle_,               Char_WindAngle              ,      ListNumber)
            call AddPropList (CloudCover_,              Char_CloudCover             ,      ListNumber)
            call AddPropList (Irrigation_,              Char_Irrigation             ,      ListNumber)
            call AddPropList (SunHours_ ,               Char_SunHours               ,      ListNumber)
            call AddPropList (ATMTransmitivity_ ,       Char_ATMTransmitivity       ,      ListNumber)
            call AddPropList (AtmospDeposOxidNO3_ ,     Char_AtmospDeposOxidNO3     ,      ListNumber)
            call AddPropList (AtmospDeposReduNH4_ ,     Char_AtmospDeposReduNH4     ,      ListNumber)
            call AddPropList (WindGust_ ,               Char_WindGust               ,      ListNumber)
            call AddPropList (PBLHeight_ ,              Char_PBLHeight              ,      ListNumber)
            call AddPropList (MeanSeaLevelPressure_ ,   Char_MeanSeaLevelPressure,      ListNumber)
            call AddPropList (WindModulus_ ,            Char_WindModulus         ,      ListNumber)
            call AddPropList (WindModulusBeaufort_ ,    Char_WindModulusBeaufort ,      ListNumber)            
            call AddPropList (WindDirection_ ,          Char_WindDirection       ,      ListNumber)
            call AddPropList (HydrogenSulfide_ ,        Char_HydrogenSulfide     ,      ListNumber)
            call AddPropList (MethylMercaptan_ ,        Char_MethylMercaptan     ,      ListNumber)            
            call AddPropList (Visibility_ ,             Char_Visibility          ,      ListNumber)
            call AddPropList (Dust_ ,                   Char_Dust                ,      ListNumber)
         
            call AddPropList (RPOM_ ,                   Char_RPOM                ,      ListNumber)
            call AddPropList (LPOM_ ,                   Char_LPOM                ,      ListNumber)
            call AddPropList (LDOM_ ,                   Char_LDOM                ,      ListNumber)
            call AddPropList (RDOM_ ,                   Char_RDOM                ,      ListNumber)
            call AddPropList (PSilica_ ,                Char_PSilica             ,      ListNumber)
            call AddPropList (DSilica_ ,                Char_DSilica             ,      ListNumber)
            call AddPropList (ICarbon_ ,                Char_ICarbon             ,      ListNumber)
            call AddPropList (pH_ ,                     Char_pH                  ,      ListNumber)
            call AddPropList (HCO3_ ,                   Char_HCO3                ,      ListNumber)
            call AddPropList (CO3_ ,                    Char_CO3                 ,      ListNumber)
            call AddPropList (Algae_1_ ,                Char_Algae_1             ,      ListNumber)
            call AddPropList (Algae_2_ ,                Char_Algae_2             ,      ListNumber)
            call AddPropList (Algae_3_ ,                Char_Algae_3             ,      ListNumber)
            call AddPropList (Algae_4_ ,                Char_Algae_4             ,      ListNumber)
            call AddPropList (Algae_5_ ,                Char_Algae_5             ,      ListNumber)
            call AddPropList (Epiphyton_1_ ,            Char_Epiphyton_1         ,      ListNumber)
            call AddPropList (Epiphyton_2_ ,            Char_Epiphyton_2         ,      ListNumber)
            call AddPropList (Epiphyton_3_ ,            Char_Epiphyton_3         ,      ListNumber)
            call AddPropList (Epiphyton_4_ ,            Char_Epiphyton_4         ,      ListNumber)
            call AddPropList (Epiphyton_5_ ,            Char_Epiphyton_5         ,      ListNumber)
            call AddPropList (Alkalinity_ ,             Char_Alkalinity          ,      ListNumber)
            call AddPropList (Detritus_ ,               Char_Detritus            ,      ListNumber)
 
            call AddPropList (ANLim_,                   Char_ANLim,                      ListNumber)       
            call AddPropList (APLim_,                   Char_APLim,                      ListNumber)       
            call AddPropList (ASLim_,                   Char_ASLim,                      ListNumber)       
            call AddPropList (ALightLim_,               Char_ALightLim,                  ListNumber)   
            call AddPropList (AOverallLim_,             Char_AOverallLim,                ListNumber) 
            call AddPropList (AGR_,                     Char_AGR,                        ListNumber)       
            call AddPropList (AMR_,                     Char_AMR,                        ListNumber)       
            call AddPropList (AER_,                     Char_AER,                        ListNumber)       
            call AddPropList (ARR_,                     Char_ARR,                        ListNumber)               
            call AddPropList (ENLim_,                   Char_ENLim,                      ListNumber)               
            call AddPropList (EPLim_,                   Char_EPLim,                      ListNumber)   
            call AddPropList (ESLim_,                   Char_ESLim,                      ListNumber)      
            call AddPropList (ELightLim_,               Char_ELightLim,                  ListNumber)              
            call AddPropList (EOverallLim_,             Char_EOverallLim,                ListNumber)       
            call AddPropList (NH4D_,                    Char_NH4D,                       ListNumber) 
            call AddPropList (NO3D_,                    Char_NO3D,                       ListNumber) 
            call AddPropList (LDOMD_,                   Char_LDOMD,                      ListNumber) 
            call AddPropList (RDOMD_,                   Char_RDOMD,                      ListNumber) 
            call AddPropList (LPOMD_,                   Char_LPOMD,                      ListNumber)            
            call AddPropList (RPOMD_,                   Char_RPOMD,                      ListNumber) 
            call AddPropList (LRDOMD_,                  Char_LRDOMD,                     ListNumber)           
            call AddPropList (LRPOMD_,                  Char_LRPOMD,                     ListNumber)        
            call AddPropList (CBODD_,                   Char_CBODD,                      ListNumber)         
            call AddPropList (PO4ER_,                   Char_PO4ER,                      ListNumber)     
            call AddPropList (PO4EG_,                   Char_PO4EG,                      ListNumber)     
            call AddPropList (PO4AR_,                   Char_PO4AR,                      ListNumber)     
            call AddPropList (PO4AG_,                   Char_PO4AG,                      ListNumber)     
            call AddPropList (PO4OM_,                   Char_PO4OM,                      ListNumber)     
            call AddPropList (PO4BOD_,                  Char_PO4BOD,                     ListNumber)    
            call AddPropList (NH4ER_,                   Char_NH4ER,                      ListNumber)  
            call AddPropList (NH4EG_,                   Char_NH4EG,                      ListNumber)  
            call AddPropList (NH4AR_,                   Char_NH4AR,                      ListNumber)  
            call AddPropList (NH4AG_,                   Char_NH4AG,                      ListNumber)    
            call AddPropList (NH4OM_,                   Char_NH4OM,                      ListNumber)    
            call AddPropList (NH4BOD_,                  Char_NH4BOD,                     ListNumber)   
            call AddPropList (NO3AG_,                   Char_NO3AG,                      ListNumber)   
            call AddPropList (NO3EG_,                   Char_NO3EG,                      ListNumber)   
            call AddPropList (DSIAG_,                   Char_DSIAG,                      ListNumber) 
            call AddPropList (DSIEG_,                   Char_DSIEG,                      ListNumber) 
            call AddPropList (DSID_,                    Char_DSID,                       ListNumber)  
            call AddPropList (PSIAM_,                   Char_PSIAM,                      ListNumber)     
            call AddPropList (PSID_,                    Char_PSID,                       ListNumber)      
            call AddPropList (LDOMAP_,                  Char_LDOMAP,                     ListNumber)    
            call AddPropList (LDOMEP_,                  Char_LDOMEP,                     ListNumber)    
            call AddPropList (LPOMAP_,                  Char_LPOMAP,                     ListNumber)    
            call AddPropList (DOAP_,                    Char_DOAP,                       ListNumber)      
            call AddPropList (DOEP_,                    Char_DOEP,                       ListNumber)      
            call AddPropList (DOAR_,                    Char_DOAR,                       ListNumber)      
            call AddPropList (DOER_,                    Char_DOER,                       ListNumber)      
            call AddPropList (DOOM_,                    Char_DOOM,                       ListNumber)       
            call AddPropList (DONIT_,                   Char_DONIT,                      ListNumber)      
            call AddPropList (ICarbonAP_,               Char_ICarbonAP,                  ListNumber) 
            call AddPropList (ICarbonEP_,               Char_ICarbonEP,                  ListNumber) 
            call AddPropList (ICarbonBOD_,              Char_ICarbonBOD,                 ListNumber)

            call AddPropList (Diameter_,                Char_Diameter,                   ListNumber)
            call AddPropList (Percentage_,              Char_Percentage,                 ListNumber)
            call AddPropList (D35_,                     Char_D35,                        ListNumber)
            call AddPropList (D50_,                     Char_D50,                        ListNumber)
            call AddPropList (D90_,                     Char_D90,                        ListNumber)
            call AddPropList (BedRock_,                 Char_BedRock,                    ListNumber)
            call AddPropList (SandTauCritic_,           Char_SandTauCritic,              ListNumber)
            call AddPropList (Sand_,                    Char_Sand,                       ListNumber)
            call AddPropList (MappDZ_,                  Char_MappDZ,                     ListNumber)            

            call AddPropList (TransportCapacity_,       Char_TransportCapacity,          ListNumber)
            call AddPropList (TransportCapacityX_,      Char_TransportCapacityX,         ListNumber)
            call AddPropList (TransportCapacityY_,      Char_TransportCapacityY,         ListNumber)
            call AddPropList (BottomEvolution_,         Char_BottomEvolution  ,          ListNumber)
            call AddPropList (Newbathymetry_,           Char_Newbathymetry    ,          ListNumber)
            call AddPropList (bathymetry_,              Char_bathymetry    ,             ListNumber)
            !wave dynamics
            call AddPropList (WaveStressX_,             Char_WaveStressX,                ListNumber)
            call AddPropList (WaveStressY_,             Char_WaveStressY,                ListNumber)
            call AddPropList (WaveStress_,              Char_WaveStress,                 ListNumber)
            call AddPropList (CurrentX_,                Char_CurrentX,                   ListNumber)
            call AddPropList (CurrentY_,                Char_CurrentY,                   ListNumber)
            call AddPropList (WaveX_,                   Char_WaveX,                      ListNumber)
            call AddPropList (WaveY_,                   Char_WaveY,                      ListNumber)
            
! Modified by Matthias DELPEY - 24/06/2011 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modified by Matthias DELPEY - 21/07/2011 - 04/08/2011 - 05/09/2011 - 25/10/2011 - 14/12/2011 
!                             - 16/12/2011 - 02/03/2012
            call AddPropList (WaveRad3DX_,              Char_WaveRad3DX,                 ListNumber)
            call AddPropList (WaveRad3DY_,              Char_WaveRad3DY,                 ListNumber)
            
            call AddPropList (WavePressureJ_,           Char_WavePressureJ,              ListNumber)
            call AddPropList (AtmToWaveMomentumU_,      Char_AtmToWaveMomentumU,         ListNumber)
            call AddPropList (AtmToWaveMomentumV_,      Char_AtmToWaveMomentumV,         ListNumber)
            call AddPropList (WaveToOceanMomentumU_,    Char_WaveToOceanMomentumU,       ListNumber)
            call AddPropList (WaveToOceanMomentumV_,    Char_WaveToOceanMomentumV,       ListNumber)
            call AddPropList (WaveDriftSpecU_,          Char_WaveDriftSpecU,             ListNumber)
            call AddPropList (WaveDriftSpecV_,          Char_WaveDriftSpecV,             ListNumber)
            call AddPropList (StokesDriftU_,            Char_StokesDriftU,               ListNumber)
            call AddPropList (StokesDriftV_,            Char_StokesDriftV,               ListNumber)
            call AddPropList (StokesDriftModulus_,      Char_StokesDriftModulus,         ListNumber)
            call AddPropList (StokesDriftW_,            Char_StokesDriftW,               ListNumber)
            call AddPropList (GlmVelocityU_,            Char_GlmVelocityU,               ListNumber)
            call AddPropList (GlmVelocityV_,            Char_GlmVelocityV,               ListNumber)
            call AddPropList (GlmVelocityModulus_,      Char_GlmVelocityModulus,         ListNumber)
            call AddPropList (GlmVelocityW_,            Char_GlmVelocityW,               ListNumber)
            call AddPropList (BreakingWaveHeight_,      Char_BreakingWaveHeight,         ListNumber)
            call AddPropList (WaveSurfaceFluxTKE_,      Char_WaveSurfaceFluxTKE,         ListNumber)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!Monocromatic:
            call AddPropList (WaveLength_,              Char_WaveLength,                 ListNumber)
            call AddPropList (WaveAmplitude_,           Char_WaveAmplitude,              ListNumber)
            call AddPropList (WavePeriod_,              Char_WavePeriod,                 ListNumber)
            call AddPropList (WaveDirection_,           Char_WaveDirection,              ListNumber)
            !!Statistical wave parametres (WW3,SWAN)
            call AddPropList (SignificantWaveHeight_,   Char_SignificantWaveHeight,      ListNumber)
            call AddPropList (SignificantWaveHeightBeaufort_,   Char_SignificantWaveHeightBeaufort,      ListNumber)            
            call AddPropList (MeanWaveLength_,          Char_MeanWaveLength,             ListNumber)
            call AddPropList (MeanWavePeriod_,          Char_MeanWavePeriod,             ListNumber)
            call AddPropList (MeanWaveDirection_,       Char_MeanWaveDirection,          ListNumber)
            call AddPropList (MeanDirectionalSpread_,   Char_MeanDirectionalSpread,      ListNumber)
            call AddPropList (PeakFrequency_,           Char_PeakFrequency,              ListNumber)
            call AddPropList (PeakDirection_,           Char_PeakDirection,              ListNumber)
            call AddPropList (PeakDirectionX_,          Char_PeakDirectionX,             ListNumber)
            call AddPropList (PeakDirectionY_,          Char_PeakDirectionY,             ListNumber)            
            call AddPropList (PeakPeriod_,              Char_PeakPeriod,                 ListNumber)
            call AddPropList (WindSeaPeakFrequency_,    Char_WindSeaPeakFrequency,       ListNumber)
            call AddPropList (WindSeaPeakDirection_,    Char_WindSeaPeakDirection,       ListNumber)
            call AddPropList (WaveSwellHeight_     ,    Char_WaveSwellHeight,            ListNumber)
            call AddPropList (MeanAbsoluteZeroCrossingPeriod_, Char_MeanAbsoluteZeroCrossingPeriod,      ListNumber)
            call AddPropList (MeanAbsoluteWavePeriodEnergy_, Char_MeanAbsoluteWavePeriodEnergy,          ListNumber)
            call AddPropList (WavePower_,               Char_WavePower,                  ListNumber)
            call AddPropList (TransportEnergyX_,        Char_TransportEnergyX,           ListNumber)
            call AddPropList (TransportEnergyY_,        Char_TransportEnergyY,           ListNumber)
            call AddPropList (Ubw_,                     Char_Ubw,                        ListNumber)
            call AddPropList (DirectionEnergyTransport_, Char_DirectionEnergyTransport,  ListNumber)
            call AddPropList (SmoothedPeakPeriod_,      Char_SmoothedPeakPeriod,         ListNumber)
            call AddPropList (MeanAbsoluteWavePeriod_,  Char_MeanAbsoluteWavePeriod,     ListNumber)
            call AddPropList (PeakWaveLength_,          Char_PeakWaveLength,             ListNumber)
            
            call AddPropList (Swell01_SignificantWaveHeight_, Char_Swell01_SignificantWaveHeight, ListNumber)
            call AddPropList (Swell01_WavePeriod_,            Char_Swell01_WavePeriod,            ListNumber)
            call AddPropList (Swell01_WaveDirection_,         Char_Swell01_WaveDirection,         ListNumber)
            
            call AddPropList (WindSea_SignificantWaveHeight_, Char_WindSea_SignificantWaveHeight, ListNumber)
            call AddPropList (WindSea_WavePeriod_,            Char_WindSea_WavePeriod,            ListNumber)
            call AddPropList (WindSea_WaveDirection_,         Char_WindSea_WaveDirection,         ListNumber)            
            
            call AddPropList (ConsolidationFlux_,       Char_ConsolidationFlux,          ListNumber)
            call AddPropList (Porosity_,                Char_Porosity,                   ListNumber)
            
            call AddPropList (ShearStress_,             Char_ShearStress_,               ListNumber)
            call AddPropList (ShearStressX_,            Char_ShearStressX_,              ListNumber)            
            call AddPropList (ShearStressY_,            Char_ShearStressY_,              ListNumber)                        

            call AddPropList (RefEvapotrans_,           Char_RefEvapotrans,              ListNumber)
            call AddPropList (TotalPlantBiomass_,       Char_TotalPlantBiomass,          ListNumber)
            call AddPropList (TotalPlantNitrogen_,      Char_TotalPlantNitrogen,         ListNumber)
            call AddPropList (TotalPlantPhosphorus_,    Char_TotalPlantPhosphorus,       ListNumber)
            call AddPropList (RootBiomass_,             Char_RootBiomass,                ListNumber)
            call AddPropList (RootDepth_,               Char_RootDepth,                  ListNumber)
            call AddPropList (LeafAreaIndex_,           Char_LeafAreaIndex,              ListNumber)
            call AddPropList (SpecificLeafStorage_,     Char_SpecificLeafStorage,        ListNumber)
            call AddPropList (EVTPCropCoefficient_,     Char_EVTPCropCoefficient,        ListNumber)
            call AddPropList (CanopyHeight_,            Char_CanopyHeight,               ListNumber)
            call AddPropList (PotLeafAreaIndex_,        Char_PotLeafAreaIndex,           ListNumber)
            call AddPropList (BoundaryLeafAreaIndex_,   Char_BoundaryLeafAreaIndex,      ListNumber)

            call AddPropList (TSS_,                     Char_TSS,                        ListNumber)
            call AddPropList (COHSED_FINE_,             Char_Cohsed_fine,                ListNumber)
            call AddPropList (COHSED_MEDIUM_,           Char_Cohsed_medium,              ListNumber)
            call AddPropList (COHSED_COARSE_,           Char_Cohsed_coarse,              ListNumber)
            call AddPropList (VSS_,                     Char_VSS,                        ListNumber)

            !PhreeqCRM
            call AddPropList (WaterSaturation_,         Char_WaterSaturation,            ListNumber)
            call AddPropList (Pressure_,                Char_Pressure,                   ListNumber)
            call AddPropList (CellPorosity_,            Char_CellPorosity,               ListNumber)
            !PhreeqC temporary code for tests
            !call AddPropList (SolutionMagnesium_,       Char_SolutionMagnesium,          ListNumber)
            !call AddPropList (SolutionCalcium_,         Char_SolutionCalcium,            ListNumber)
            !call AddPropList (SolutionSodium_,          Char_SolutionSodium,             ListNumber)
            !call AddPropList (SolutionNitrogenGas_,     Char_SolutionNitrogenGas,        ListNumber)
            !call AddPropList (SolutionOxygenGas_,       Char_SolutionOxygenGas,          ListNumber)
            !call AddPropList (SolutionAmmonia_,         Char_SolutionAmmonia,            ListNumber)
            !call AddPropList (SolutionNitrate_,         Char_SolutionNitrate,            ListNumber)
            !call AddPropList (SolutionNitrite_,         Char_SolutionNitrite,            ListNumber)
            !call AddPropList (SolutionChlorine_,        Char_SolutionChlorine,           ListNumber)
            !call AddPropList (GasN2_,                   Char_GasN2,                      ListNumber)
            !call AddPropList (GasCO2_,                  Char_GasCO2,                     ListNumber)
            !call AddPropList (pE_,                      Char_pE,                         ListNumber)
            !call AddPropList (eCaX2_,                   Char_eCaX2,                      ListNumber)
            !call AddPropList (eMgX2_,                   Char_eMgX2,                      ListNumber)
            !call AddPropList (eNaX_,                    Char_eNaX,                       ListNumber)
            !call AddPropList (eNH4X_,                   Char_eNH4X,                      ListNumber)
            !call AddPropList (eKX_,                     Char_eKX,                        ListNumber)
            !call AddPropList (sCa2_,                    Char_sCa2,                       ListNumber)
            !call AddPropList (sCaOH_,                   Char_sCaOH,                      ListNumber)
            !call AddPropList (sH2_,                     Char_sH2,                        ListNumber)
            !call AddPropList (sMg2_,                    Char_sMg2,                       ListNumber)
            !call AddPropList (sMgOH_,                   Char_sMgOH,                      ListNumber)
            !call AddPropList (sNH4_,                    Char_sNH4,                       ListNumber)
            !call AddPropList (sNH3_,                    Char_sNH3,                       ListNumber)
            !call AddPropList (sN2_,                     Char_sN2,                        ListNumber)
            !call AddPropList (sNO2_,                    Char_sNO2,                       ListNumber)
            !call AddPropList (sNO3_,                    Char_sNO3,                       ListNumber)
            !call AddPropList (sNa_,                     Char_sNa,                        ListNumber)
            !call AddPropList (sNaOH_,                   Char_sNaOH,                      ListNumber)
            !call AddPropList (sO2_,                     Char_sO2,                        ListNumber)
            !call AddPropList (msmSolutionCalcium_,      Char_msmSolutionCalcium,         ListNumber)
            !call AddPropList (msmSolutionMagnesium_,    Char_msmSolutionMagnesium,       ListNumber)
            !call AddPropList (msmSolutionSodium_,       Char_msmSolutionSodium,          ListNumber)
            !call AddPropList (msmSolutionAmmonia_,      Char_msmSolutionAmmonia,         ListNumber)
            !call AddPropList (Calcite_,                 Char_Calcite,                    ListNumber)
            !call AddPropList (Dolomite_,                Char_Dolomite,                   ListNumber)
            !call AddPropList (Aragonite_,               Char_Aragonite,                  ListNumber)
            !call AddPropList (Halite_,                  Char_Halite,                     ListNumber)
            !call AddPropList (KFeldspar_,               Char_KFeldspar,                  ListNumber)
            !call AddPropList (SolutionCarbon_,          Char_SolutionCarbon,             ListNumber) 
            !call AddPropList (SolutionPotassium_,       Char_SolutionPotassium,          ListNumber) 
            !call AddPropList (SolutionAluminium_,       Char_SolutionAluminium,          ListNumber) 
            !call AddPropList (SolutionSilicium_,        Char_SolutionSilicium,           ListNumber) 
            !call AddPropList (RainMagnesium_,           Char_RainMagnesium,              ListNumber) 
            !call AddPropList (RainCalcium_,             Char_RainCalcium,                ListNumber) 
            !call AddPropList (RainSodium_,              Char_RainSodium,                 ListNumber)
            !call AddPropList (RainChlorine_,            Char_RainChlorine,               ListNumber)  
            !call AddPropList (RainAmmonia_,             Char_RainAmmonia,                ListNumber)  
            !END of PhreeqC temporary code for tests
            call AddPropList (SoilVolumetricDensity_,   Char_SoilVolumetricDensity,      ListNumber)  
            call AddPropList (SolEC_,                   Char_SolEC,                      ListNumber)
            call AddPropList (GenericDissPesticide_1_,  Char_GenericDissPesticide_1,     ListNumber)
            call AddPropList (GenericDissPesticide_2_,  Char_GenericDissPesticide_2,     ListNumber)
            call AddPropList (GenericDissPesticide_3_,  Char_GenericDissPesticide_3,     ListNumber) 
            call AddPropList (GenericDissPesticide_4_,  Char_GenericDissPesticide_4,     ListNumber) 
            call AddPropList (GenericPartPesticide_1_,  Char_GenericPartPesticide_1,     ListNumber)
            call AddPropList (GenericPartPesticide_2_,  Char_GenericPartPesticide_2,     ListNumber)
            call AddPropList (GenericPartPesticide_3_,  Char_GenericPartPesticide_3,     ListNumber) 
            call AddPropList (GenericPartPesticide_4_,  Char_GenericPartPesticide_4,     ListNumber) 
            call AddPropList (IndividualsPerCell_,      Char_IndividualsPerCell,         ListNumber)             
            !WWTPQProperties
            call AddPropList (SolInertOrgMat_,          Char_SolInertOrgMat,             ListNumber)
            call AddPropList (ReadilyBioSub_,           Char_ReadilyBioSub,              ListNumber)
            call AddPropList (PartInertOrgMar_,         Char_PartInertOrgMar ,           ListNumber)
            call AddPropList (SlowlyBioSub_,            Char_SlowlyBioSub,               ListNumber)
            call AddPropList (HetBio_,                  Char_HetBio,                     ListNumber)
            call AddPropList (AutBio_,                  Char_AutBio,                     ListNumber)
            call AddPropList (PartProd_,                Char_PartProd,                   ListNumber)
            call AddPropList (SolBioOrgNitrogen_,       Char_SolBioOrgNitrogen,          ListNumber)
            call AddPropList (PartBioOrgNitrogen_,      Char_PartBioOrgNitrogen,         ListNumber)
            call AddPropList (SnowPack_,                Char_SnowPack,                   ListNumber)
            call AddPropList (DailyAvgTemp_,            Char_DailyAvgTemp,               ListNumber)
            call AddPropList (ForestCoverFraction_,     Char_ForestCoverFraction,        ListNumber)
            call AddPropList (SnowSlopeFactor_,         Char_SnowSlopeFactor,            ListNumber)
            call AddPropList (SnowPrecipitation_,       Char_SnowPrecipitation,          ListNumber)
            call AddPropList (SnowWaterEquivalent_,     Char_SnowWaterEquivalent,        ListNumber)
            call AddPropList (SurfaceDownLatentHeat_,   Char_SurfaceDownLatentHeat,      ListNumber)
            call AddPropList (SurfaceDownSensibleHeat_, Char_SurfaceDownSensibleHeat,    ListNumber)     
            call AddPropList (CellPercentContamin_,     Char_CellPercentContamin,        ListNumber)                  
            call AddPropList (SolutionMapping_,         Char_SolutionMapping,            ListNumber)                  
            call AddPropList (EquilibriumMapping_,      Char_EquilibriumMapping,         ListNumber)                  
            call AddPropList (ExchangeMapping_,         Char_ExchangeMapping,            ListNumber)                  
            call AddPropList (SurfaceMapping_,          Char_SurfaceMapping,             ListNumber)                  
            call AddPropList (GasPhaseMapping_,         Char_GasPhaseMapping,            ListNumber)                  
            call AddPropList (SolidSolutionMapping_,    Char_SolidSolutionMapping,       ListNumber)                  
            call AddPropList (KineticsMapping_,         Char_KineticsMapping,            ListNumber)  
            call AddPropList (ApplicationArea_,         Char_ApplicationArea,            ListNumber)  
            call AddPropList (FixedIrrigation_,         Char_FixedIrrigation,            ListNumber)  
            call AddPropList (AccIrrigation_,           Char_AccIrrigation,              ListNumber)  
            !Place to add new properties to the names list
        
            !Ends building the property name list
            call AddPropList(FillValueInt, 'No name', ListNumber, EndList = .true.)


            PropertiesNumber = ListNumber

            do i = 1, PropertiesNumber
            do j = 1, PropertiesNumber

                if (i/=j .and. (PropNameList  (i) == PropNameList  (j) .or.           &
                                PropNumberList(i) == PropNumberList(j)))              &
                                stop 'ConstructPropList - GlobalData, ERR01'
            enddo
            enddo
            
        endif    



    end subroutine ConstructPropList
    
    !--------------------------------------------------------------------------

    subroutine AddPropList(PropNumber, PropName, ListNumber, EndList)

        !Arguments-------------------------------------------------------------
        integer,          intent (IN)               :: PropNumber
        character(len=*), intent (IN)               :: PropName
        integer,          intent (INOUT)            :: ListNumber
        logical,          intent (IN),    optional  :: EndList

        !Begin-----------------------------------------------------------------
        logical                                                      :: EndList_
        character(len=StringLength), allocatable, dimension(:), save :: AuxChar
        integer,                     allocatable, dimension(:), save :: AuxInt
        logical                                               , save :: FirstProp = .true.
        
        !----------------------------------------------------------------------


        if (present(EndList)) then 

            EndList_ = EndList

        else

            EndList_ = .false.

        endif

       if (.not. EndList_) then

            if (FirstProp) then

                FirstProp = .false.

                allocate (PropNumberList(1))                
                allocate (PropNameList  (1))
                allocate (AuxInt        (1))
                allocate (AuxChar       (1))
                
                ListNumber = 1
    
            else     

                deallocate   (PropNameList)
                deallocate   (PropNumberList)

                allocate     (PropNameList  (ListNumber + 1))      
                allocate     (PropNumberList(ListNumber + 1))      

                PropNameList  (1 : ListNumber) = AuxChar(1 : ListNumber)
                PropNumberList(1 : ListNumber) = AuxInt (1 : ListNumber)
 
                deallocate   (AuxInt)
                deallocate   (AuxChar)
                allocate     (AuxChar(ListNumber + 1))
                allocate     (AuxInt(ListNumber + 1))

                ListNumber = ListNumber + 1

            endif     

            PropNameList  (ListNumber) = PropName
            PropNumberList(ListNumber) = PropNumber

            AuxChar      (:)          = PropNameList  (:)
            AuxInt       (:)          = PropNumberList(:)

        else       

            deallocate (AuxChar)
            deallocate (AuxInt)

        endif      

        if (PropNumber > HighestPropertyID) HighestPropertyID = PropNumber

    end subroutine AddPropList

    logical function isVSS(Property)
        !Arguments-------------------------------------------------------------
        integer, intent (IN) :: Property

        !----------------------------------------------------------------------

cd1 :   if ((Property == Phytoplankton_         ) .OR.  (Property == Diatoms_               ) .OR.          &
            (Property == Zooplankton_           ) .OR.  (Property == Ciliate_               ) .OR.          &
            (Property == Bacteria_              ) .OR.  (Property == PON_                   ) .OR.          &
            (Property == PONRefractory_         )) then
            
            isVSS = .TRUE.

        else

            isVSS = .FALSE.

        end if cd1

    end function isVSS
    


    !----------------------------------------------------------------------


    logical function Check_Particulate_Property(Property)

        !Arguments-------------------------------------------------------------
        integer, intent (IN) :: Property

        !----------------------------------------------------------------------

cd1 :   if ((Property == POC_                   ) .OR.  (Property == PON_                   ) .OR.          &
            (Property == PONRefractory_         ) .OR.  (Property == Bacteria_              ) .OR.          &
            (Property == POP_                   ) .OR.  (Property == ParticulateArsenic_    ) .OR.          &
            (Property == AdsorbedAmmonia_       ) .OR.  (Property == AnaerobicC_            ) .OR.          &
            (Property == AnaerobicN_            ) .OR.  (Property == AutotrophicC_          ) .OR.          &
            (Property == AutotrophicN_          ) .OR.  (Property == HeterotrophicC_        ) .OR.          &
            (Property == HeterotrophicN_        ) .OR.  (Property == LabileOrganicC_        ) .OR.          &
            (Property == RefreactaryOrganicC_   ) .OR.  (Property == RefreactaryOrganicN_   ) .OR.          &            
            (Property == Cohesive_Sediment_     ) .OR.  (Property == ParticulateContaminant_) .OR.          &  
            (Property == ParticulateSodium_     ) .OR.  (Property == ParticulateCalcium_    ) .OR.          &  
            (Property == Detritus_              ) .OR.  (Property == Phytoplankton_         ) .OR.          &  
            (Property == Diatom_C_              ) .OR.  (Property == Diatom_N_              ) .OR.          &  
            (Property == Diatom_P_              ) .OR.  (Property == Diatom_Si_             ) .OR.          &
            (Property == Diatom_Chl_            ) .OR.  (Property == Mix_Flagellate_C_      ) .OR.          &
            (Property == Mix_Flagellate_N_      ) .OR.  (Property == Mix_Flagellate_P_      ) .OR.          &
            (Property == Mix_Flagellate_Chl_    ) .OR.  (Property == Picoalgae_C_           ) .OR.          &
            (Property == Picoalgae_N_           ) .OR.  (Property == Picoalgae_P_           ) .OR.          &
            (Property == Picoalgae_Chl_         ) .OR.  (Property == Flagellate_C_          ) .OR.          &
            (Property == Flagellate_N_          ) .OR.  (Property == Flagellate_P_          ) .OR.          &
            (Property == Flagellate_Chl_        ) .OR.  (Property == Algae_1_               ) .OR.          &  
            (Property == Algae_2_               ) .OR.  (Property == Algae_3_               ) .OR.          &  
            (Property == Algae_4_               ) .OR.  (Property == Algae_5_               ) .OR.          &
            (Property == Epiphyton_1_           ) .OR.  (Property == Epiphyton_2_           ) .OR.          &
            (Property == Epiphyton_3_           ) .OR.  (Property == Epiphyton_4_           ) .OR.          &  
            (Property == Epiphyton_5_           ) .OR.  (Property == PSilica_               ) .OR.          &  
            (Property == LPOM_                  ) .OR.  (Property == RPOM_                  ) .OR.          &  
            (Property == BioSilica_             ) .OR.  (Property == Diatoms_               ) .OR.          &
            (Property == Larvae_                ) .OR.  (Property == COHSED_FINE_           ) .OR.          &
            (Property == COHSED_MEDIUM_         ) .OR.  (Property == COHSED_COARSE_         ) .OR.          &
            (Property == VSS_                   ) .OR.  (Property == TSS_                   ) .OR.          &
            (Property == MacroAlgae_            ) .OR.  (Property == AnaerobicP_            ) .OR.          &
            (Property == AutotrophicP_          ) .OR.  (Property == HeterotrophicP_        ) .OR.          &
            (Property == RefreactaryOrganicP_   ) .OR.  (Property == AdsorbedInorganicP_    ) .OR.          &
            (Property == SolubilizingN_         ) .OR.  (Property == SolubilizingC_         ) .OR.          &                     
            (Property == PON1_                  ) .OR.  (Property == PON2_                  ) .OR.          &
            (Property == PON3_                  ) .OR.  (Property == PON4_                  ) .OR.          &
            (Property == PON5_                  ) .OR.  (Property == POP1_                  ) .OR.          &
            (Property == POP2_                  ) .OR.  (Property == POP3_                  ) .OR.          &
            (Property == POP4_                  ) .OR.  (Property == POP5_                  ) .OR.          &
            (Property == GenericPartPesticide_1_) .OR.  (Property == GenericPartPesticide_2_) .OR.          &            
            (Property == GenericPartPesticide_3_) .OR.  (Property == GenericPartPesticide_4_) .OR.          &            
            (Property == SolubilizingP_         ) .OR.  (Property == ParticulateMetal_      ) .OR.          &
            (Property == ParticulateCadmium_    ) .OR.  (Property == ParticulateCopper_     ) .OR.          &
            (Property == ParticulateLead_       ) .OR.  (Property == ParticulateZinc_       ) .OR.          &
            (Property == ParticulateMercury_    ) .OR.  (Property == AnaerobicPop_          ) .OR.          &
            (Property == AutotrophicPop_        ) .OR.  (Property == HeterotrophicPop_      ) .OR.          &
            (Property == SolPop_                ) .OR.  (Property == IndividualsPerCell_    ) .OR.          &
            (Property == SuspensionFeedersC_    ) .OR.  (Property == SuspensionFeedersN_    ) .OR.          &   !isab    
            (Property == SuspensionFeedersP_    ) .OR.  (Property == MicroPhytoBenthosN_    ) .OR.          &   !isab    
            (Property == MicroPhytoBenthosP_    ) .OR.  (Property == MicroPhytoBenthosC_    ) .OR.          &
            (Property == SeagrassesN_           ) .OR.  (Property == SeagrassesP_           ) .OR.          &
            (Property == SeagrassesLeaves_      ) .OR.  (Property == SeagrassesRoots_       ) .OR.          &
            (Property == DepositFeedersN_       ) .OR.  (Property == DepositFeedersP_       ) .OR.          &
            (Property == DepositFeedersC_       ) .OR.  (Property == Zooplankton_           ) .OR.          &
            (Property == CellPercentContamin_   )) then

            Check_Particulate_Property = .TRUE. 
        
        else
        
            Check_Particulate_Property = .FALSE.

        end if cd1

    end function Check_Particulate_Property

        
    !--------------------------------------------------------------------------
    
    logical function Check_Angle_Property(Property)

        !Arguments-------------------------------------------------------------
         integer, intent (IN) :: Property

        !----------------------------------------------------------------------

cd1 :   if (Property == WindDirection_           .OR.                                   &
            Property == MeanWaveDirection_       .OR.                                   &
            Property == VelocityDirection_       .OR.                                   &
            Property == WaveDirection_           .OR.                                   &
            Property == MeanDirectionalSpread_   .OR.                                   &
            Property == PeakDirection_           .OR.                                   &
            Property == WindSeaPeakDirection_) then

            Check_Angle_Property = .TRUE. 
                
        else

            Check_Angle_Property = .FALSE.

        end if cd1

    end function Check_Angle_Property   

        
    !--------------------------------------------------------------------------    

    !--------------------------------------------------------------------------
            
    integer function Get_Angle_Referential(PropertyID)

        !Arguments-------------------------------------------------------------
        type(T_PropertyID), pointer :: PropertyID

        !----------------------------------------------------------------------
        if (PropertyID%IsAngle) then
            
cd1 :       if ((PropertyID%IDNumber == WindDirection_) .OR. (PropertyID%IDNumber == MeanWaveDirection_)) then

                Get_Angle_Referential = NauticalReferential_          

            else
            
                stop  'Get_Angle_Referential - ModuleGlobalData - ERR010'

            end if cd1
        else
        
            stop  'Get_Angle_Referential - ModuleGlobalData - ERR020'
        
        endif
    
    end function Get_Angle_Referential

        
    !--------------------------------------------------------------------------    

    !--------------------------------------------------------------------------            
    
    logical function Check_Vectorial_Property(Property)

        !Arguments-------------------------------------------------------------
        integer, intent (IN) :: Property

        !----------------------------------------------------------------------

cd1 :   if ((Property == WindVelocity_) .OR. (Property == WaveStress_) .OR. &
            (Property == WindStress_)) then

            Check_Vectorial_Property = .TRUE. 
         
        else

            Check_Vectorial_Property = .FALSE.

        end if cd1

    end function Check_Vectorial_Property

        
    !--------------------------------------------------------------------------      

    subroutine Get_Vectorial_PropertyNames(PropertyID, PropertyX, PropertyY, PropertyZ)

        !Arguments-------------------------------------------------------------
        integer, intent (IN)                         :: PropertyID
        character(len=*), intent (OUT)               :: PropertyX
        character(len=*), intent (OUT)               :: PropertyY
        character(len=*), intent (OUT), optional     :: PropertyZ
        !integer              :: i

        !----------------------------------------------------------------------

cd1 :   if ((PropertyID == WindVelocity_ )) then

            PropertyX = GetPropertyName(WindVelocityX_)
            PropertyY = GetPropertyName(WindVelocityY_)
            PropertyZ = null_str
        
        else if ((PropertyID == WaveStress_ )) then
            PropertyX = GetPropertyName(WaveStressX_)
            PropertyY = GetPropertyName(WaveStressY_)
            PropertyZ = null_str

        else if ((PropertyID == WindStress_ )) then
            PropertyX = GetPropertyName(WindStressX_)
            PropertyY = GetPropertyName(WindStressY_)
            PropertyZ = null_str            
            
        else
            PropertyX = null_str
            PropertyY = null_str
            PropertyZ = null_str
            
        end if cd1

    end subroutine Get_Vectorial_PropertyNames

        
    !--------------------------------------------------------------------------                
            
    integer function TranslateTypeZUV(Char_TypeZUV)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in)            :: Char_TypeZUV
        
        !----------------------------------------------------------------------

        select case(Char_TypeZUV)

            case("Z", "z")

                TranslateTypeZUV = TypeZ_

            case("U", "u")

                TranslateTypeZUV = TypeU_

            case("V", "v")

                TranslateTypeZUV = TypeV_

            case default

                write(*,*)'Invalid type ZUV'
                stop      'TranslateTypeZUV - ModuleGlobalData - ERR01'

        end select

    end function TranslateTypeZUV


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !  UnitsManager - Routine to index unit files numbers  
    !--------------------------------------------------------------------------
    !   Variables:                                         
    !   OPENCLOSE     - Activate unit (1) or deactivate (0)
    !   UNIT          - Number to the opened unit or to the unit to close   

    subroutine UnitsManager(UNIT, OPENCLOSE, STATUS, STAT)

        !Arguments-------------------------------------------------------------
        integer,                    intent(INOUT)   :: UNIT
        integer,                    intent(IN   )   :: OPENCLOSE
        character(len=*), optional, intent(IN   )   :: STATUS
        integer,          optional, intent(OUT  )   :: STAT

        !External--------------------------------------------------------------
        logical                                     :: opened
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character(len = StringLength)               :: STATUS_        
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

cd1 :   if      (OPENCLOSE .EQ. OPEN_FILE ) then

            Unit   = 8
            Opened = .true.
do1 :       do while (Opened)
                Unit = Unit + 1
                inquire(UNIT = Unit, OPENED = Opened)
            end do do1

            STAT_ = SUCCESS_

        else if (OPENCLOSE .EQ. CLOSE_FILE) then cd1

            inquire(UNIT = Unit, OPENED = opened)

if8 :       if (opened) then

                if (present(STATUS)) then
                    STATUS_ = STATUS
                    if (STATUS_ /= 'KEEP' .and.                                         &
                        STATUS_ /= 'SAVE' .and.                                         &
                        STATUS_ /= 'DELETE') then
                        stop       "UnitsManager - ModuleGlobalData - ERR100"
                    endif
                else
                    STATUS_ = 'KEEP'
                endif

                close(UNIT, STATUS = STATUS_, IOSTAT = STAT_CALL)

if9 :           if (STAT_CALL .NE. SUCCESS_) then
                    write(*,*)  
                    write(*,*) "Error closing file unit ", Unit
                    write(*,*) "UnitsManager - ModuleGlobalData - WRN01"

                    STAT_ = UNIT_ERR_
                else
          
                    STAT_ = SUCCESS_

                end if if9
            else
                write(*,*)  
                write(*,*) "Error closing file unit ", Unit
                write(*,*) "The unit is already closed."
                write(*,*) "UnitsManager - ModuleGlobalData - WRN02"

                STAT_ = UNIT_ERR_
            end if if8

        else
            write(*,*) 
            write(*,*) "OPENCLOSE must be 0 or 1."
            write(*,*) "OPENCLOSE is ", OPENCLOSE
            stop       "UnitsManager - ModuleGlobalData - ERR01"

        end if cd1


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnitsManager

    !--------------------------------------------------------------------------

    subroutine SetErrorMessage(ErrorMagnitude, ErrorType, SmallMessage, Screen)

        !Arguments-------------------------------------------------------------
        integer,            intent(IN)  :: ErrorMagnitude
        integer,            intent(IN)  :: ErrorType
        logical, optional,  intent(IN)  :: Screen
        character(LEN = *), intent(IN)  :: SmallMessage
                                                                                   
        !Local-----------------------------------------------------------------
        character(10)                   :: StrErrorMagnitude = "          "
        character(10)                   :: StrErrorType      = "          "
        logical                         :: Screen_
                     
        !----------------------------------------------------------------------

cd5 :   if (present(Screen)) then
            Screen_ = Screen
        else cd5
            Screen_ = ON
        end if cd5


cd1 :   if      (ErrorMagnitude == FATAL_  ) then
            StrErrorMagnitude = "FATAL"
        else if (ErrorMagnitude == WARNING_) then cd1
            StrErrorMagnitude = "WARNING"
        else
            StrErrorMagnitude = "UNKNOWN"
        end if cd1


cd2 :   if      (ErrorType == INTERNAL_  ) then
            StrErrorType  = "INTERNAL"
        else if (ErrorType == OUT_OF_MEM_) then cd2
            StrErrorType  = "OUT_OF_MEM"
        else if (ErrorType == KEYWORD_   ) then cd2
            StrErrorType  = "KEYWORD"
        else
            StrErrorType  = "UNKNOWN"
        end if cd2


        !Writes a message to the error file
        !Ricardo - dont put here any trim or adjustl or blank lines, ok?
        write(ErrorFileID, "(A)") StrErrorMagnitude//space//StrErrorType//space//trim(adjustl(SmallMessage))//dot
  
cd7 :   if (Screen_) then
            write(*,                *     ) 
            write(*,                 "(A)") trim(adjustl(StrErrorMagnitude))//semicolumn//   &
                                            space//trim(adjustl(StrErrorType))//semicolumn// &
                                            space//trim(adjustl(SmallMessage))//dot
            write(*,                *     ) 
        end if cd7

        !Places Message on the Message Stack
        call PlaceErrorMessageOnStack(trim(adjustl(SmallMessage)))

        !If the error message isn"t a warning, stop the execution of this process
cd3 :   if (ErrorMagnitude == FATAL_) then
            close (unit=ErrorFileID  )
            close (unit=UsedKeyFileID)
           
            stop "Program stop anormaly. See Log file."
        end if cd3

    end subroutine SetErrorMessage
        
    !--------------------------------------------------------------------------

    subroutine SetErrorTime(ErrorMagnitude, ErrorType, SmallMessage,        &
                            Year, Month, Day, Hour, Minute, Second, Screen)

        !Arguments-------------------------------------------------------------
        integer,            intent(IN)  :: ErrorMagnitude
        integer,            intent(IN)  :: ErrorType
        character(LEN = *), intent(IN)  :: SmallMessage
        real,            intent(IN)     :: Year, Month, Day, Hour, Minute, Second
        logical, optional,  intent(IN)  :: Screen

        !Local-----------------------------------------------------------------
        character(10)                   :: StrErrorMagnitude = "          "
        character(10)                   :: StrErrorType      = "          "
        logical                         :: Screen_
                     
        !----------------------------------------------------------------------

cd5 :   if (present(Screen)) then
            Screen_ = Screen
        else cd5
            Screen_ = ON
        end if cd5

cd1 :   if      (ErrorMagnitude == FATAL_  ) then
            StrErrorMagnitude = "FATAL"
        else if (ErrorMagnitude == WARNING_) then cd1
            StrErrorMagnitude = "WARNING"
        else
            StrErrorMagnitude = "UNKNOWN"
        end if cd1

cd2 :   if      (ErrorType == INTERNAL_  ) then
            StrErrorType  = "INTERNAL"
        else if (ErrorType == OUT_OF_MEM_) then cd2
            StrErrorType  = "OUT_OF_MEM"
        else if (ErrorType == KEYWORD_   ) then cd2
            StrErrorType  = "KEYWORD"
        else
            StrErrorType  = "UNKNOWN"
        end if cd2

        !Writes a message to the error file
        write(ErrorFileID, 555) StrErrorMagnitude//space//StrErrorType//space// &
                                trim(adjustl(SmallMessage))//space,             &
                                int(Year), int(Month),  int(Day),               &
                                int(Hour), int(Minute), int(Second)
  
        !Places Message on the Message Stack
        call PlaceErrorMessageOnStack(trim(adjustl(SmallMessage)))

cd7 :   if (Screen_) then
            write(*,*) 
            write(*,555) trim(adjustl(StrErrorMagnitude))//semicolumn//     &
                  space//trim(adjustl(StrErrorType))//semicolumn//          &
                  space//trim(adjustl(SmallMessage)),                       &
                  int(Year), int(Month), int(Day), int(Hour), int(Minute), int(Second)
            write(*,*) 
        end if cd7

        555 format(1x, A60,(i4,":"),4(i2, ":"), i2)

        !If the error message isn"t a warning, stop the execution of this process
cd3 :   if (ErrorMagnitude == FATAL_) then
            close (unit=ErrorFileID  )
            close (unit=UsedKeyFileID)
            stop "Program stop anormaly. See Log file."
        end if cd3

    end subroutine SetErrorTime

    !----------------------------------------------------------------------
    
    !Messages are stored in a stack: Last written = Highest number.
    !First written is discarded if the stack size is exceeded.
    subroutine PlaceErrorMessageOnStack(ErrorMessage)
    
        !Arguments-------------------------------------------------------------
        character(len=*)                            :: ErrorMessage
        
        !Local-----------------------------------------------------------------
        integer                                     :: i

        if (NumberOfErrorMessages == 0) then
            do i=1, MaxErrorMessages
                ErrorMessagesStack(i)=' '
            enddo
        endif

        NumberOfErrorMessages =NumberOfErrorMessages + 1

        if(NumberOfErrorMessages.gt.MaxErrorMessages)then
            NumberOfErrorMessages=MaxErrorMessages
            do i=1, NumberOfErrorMessages-1
                ErrorMessagesStack(i)=ErrorMessagesStack(i+1)
            enddo
        endif

        ErrorMessagesStack(NumberOfErrorMessages)=ErrorMessage
       
    end subroutine PlaceErrorMessageOnStack

    !----------------------------------------------------------------------

    subroutine WriteDTLog (ModelName, iter, DT, i, j, k, PropertyName)

        !Arguments---------------------------------------------------------
        character(len=*)                :: ModelName
        integer                         :: iter
        real                            :: DT
        integer, optional               :: i, j, k
        character(len=*), optional      :: PropertyName
        
        !Local-----------------------------------------------------------------
        

        if (MonitorDT) then
            if (present(k)) then
                write (UnitDT, 211) iter, DT, i, j, k, trim(ModelName), trim(PropertyName)
            elseif (present(j)) then
                write (UnitDT, 212) iter, DT, i, j, trim(ModelName), trim(PropertyName)
            elseif (present(i)) then
                write (UnitDT, 213) iter, DT, i, trim(ModelName), trim(PropertyName)
            else
                write (UnitDT, 214) trim(ModelName), iter, DT
                !write (UnitDT, '(A25, I10,F12.4)') ModuleName, iter, DT
            end if
        end if

211 format(i8, ', DT: ', f9.3, ', I: ', i4, ', J: ', i4, ', K: ', i4, ', Model - Property: ', a, ' - ', a)
212 format(i8, ', DT: ', f9.3, ', I: ', i4, ', J: ', i4, ', Model - Property: ', a, ' - ', a)
213 format(i8, ', DT: ', f9.3, ', I: ', i4, ', Model - Property: ', a, ' - ', a)
214 format(i8, ', DT: ', f9.3, a)

    end subroutine WriteDTLog

    !--------------------------------------------------------------------------

    subroutine WriteDTLog_ML (ModuleName, iter, DT, Message)

        !Arguments---------------------------------------------------------
        character(len=*)                :: ModuleName
        integer                         :: iter
        real                            :: DT
        character(len=*), optional      :: Message
        
        !Local-----------------------------------------------------------------
       
        if (MonitorDT) then
            if (present(Message)) then
                write (UnitDT, '(A25, I10,F12.4,A132)') ModuleName, iter, DT, Message
            else
                write (UnitDT, '(A25, I10,F12.4)') ModuleName, iter, DT 
            end if
        end if

    end subroutine WriteDTLog_ML

    !--------------------------------------------------------------------------

    subroutine StartupMohid(ModelName)

        !Arguments-------------------------------------------------------------
        character(len=*)                    :: ModelName

        !Local-----------------------------------------------------------------
#ifndef _ONLINE_ 
        integer                             :: STAT_CALL
        integer                             :: Counter
        character(LEN=4)                    :: Number
#endif        
        integer                             :: N

        integer, dimension(:), allocatable  :: Seed

        !Begin-----------------------------------------------------------------

        !Initializes Random number generator
        !The use of the "Seed" guarantees the always the same sequence of random
        !numbers is generated. This is important for debug proposes - Frank 03/05
        call RANDOM_SEED(SIZE = N)
        allocate        (Seed(N))
        Seed = 1
        call RANDOM_SEED(PUT  = Seed)

#ifndef _OUTPUT_OFF_
        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"      AUTHOR   : IST/MARETEC, Marine Modelling Group      "
        write(*, *)"      WWW      : http://www.mohid.com"
        write(*, *)                    
        write(*, *)                    
        write(*, *)"Copyright (C) 1985, 1998, 2002, 2006."
        write(*, *)"Instituto Superior Tecnico, Technical University of Lisbon"
        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Constructing "//ModelName
        write(*, *)"Please Wait..."
#endif

#ifndef _ONLINE_ 

        !Gets a unit for the Error File    
        call UnitsManager(ErrorFileID, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'StartupMohid - GlobalData - ERR01'

         
        Counter  = 1
do1:    do
            Number = '    '
            write(Number, fmt='(i4)')Counter

            open(UNIT   = ErrorFileID,                                    &
                 FILE   = 'Error_and_Messages_'//trim(adjustl(Number))//'.log',  &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1


        !Gets a unit for the UsedKeyWordFile
        call UnitsManager(UsedKeyFileID, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'StartupMohid - GlobalData - ERR02'



        Counter  = 1
do2:    do
            Number = '    '
            write(Number, fmt='(i4)')Counter

            open(UNIT   = UsedKeyFileID,                                  &
                 FILE   = 'UsedKeyWords_'//trim(adjustl(Number))//'.dat',        &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do2
            else
                Counter = Counter + 1
            end if

        enddo do2

#endif
        !----------------------------------------------------------------------

    end subroutine StartupMohid

    !--------------------------------------------------------------------------

    subroutine ShutdownMohid (ModelName, ElapsedSeconds, TotalCPUTime,            &
                              WorkCycleElapsed, WorkCycleCPUTime)

        !Arguments-------------------------------------------------------------
        character(len=*)                :: ModelName
        real                            :: ElapsedSeconds, TotalCPUTime
        real, optional                  :: WorkCycleElapsed, WorkCycleCPUTime
        
        !Local---------------------------------------------------------------
        integer                         :: ElapsedHours, ElapsedMinutes, ElapsedSecremain

        !----------------------------------------------------------------------
#ifndef _ONLINE_ 

        integer                         :: STAT_CALL

        !Write a commentb in the error and messages to indicate that model run went ok
        write(unit=ErrorFileID, fmt=*)"SUCCESSFUL  : 1"
        
        !Close Errors and messages
        call UnitsManager(ErrorFileID, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ShutdownMohid - GlobalData - ERR01'

        call UnitsManager(UsedKeyFileID, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ShutdownMohid - GlobalData - ERR02'

#endif

        ElapsedHours = INT(Elapsedseconds/3600)
        ElapsedMinutes = INT((ElapsedSeconds-ElapsedHours*3600)/60)
        ElapsedSecremain = INT((ElapsedSeconds-ElapsedMinutes*60-ElapsedHours*3600))
#ifndef _OUTPUT_OFF_
        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Program "//ModelName//" successfully terminated"
        write(*, *)                    
        write(*, *)
        write(*, 110)ElapsedSeconds, ElapsedHours, ElapsedMinutes,ElapsedSecremain
        write(*, 120)TotalCPUTime
        write(*, 130)100.*TotalCPUTime/ElapsedSeconds
        if (present(WorkCycleElapsed) .and. present(WorkCycleCPUTime)) then
            write(*, 140)WorkCycleElapsed
            write(*, 150)WorkCycleCPUTime
            if(WorkCycleElapsed/=.0) then
                write(*, 160)100.*WorkCycleCPUTime/WorkCycleElapsed
            endif
        endif
        write(*, *)
        write(*, *)"----------------------------------------------------------"

#endif

#ifndef _OPENMI_        
        !This Stop Statement has been removed because it makes OpenMI Destructor work improberly
        stop
#endif

    110 format(1x, "Total Elapsed Time     : ",f14.2," ",i3,"h ",i2,"min ",i2,"s",/)
    120 format(1x, "Total CPU time         : ",f14.2,/)
    130 format(1x, "CPU usage (%)          : ",f14.2,/)
    140 format(1x, "Workcycle Elapsed Time : ",f14.2,/)
    150 format(1x, "Workcycle CPU time     : ",f14.2,/)
    160 format(1x, "Workcycle CPU usage (%): ",f14.2,/)

        !----------------------------------------------------------------------

    end subroutine ShutdownMohid

    !--------------------------------------------------------------------------

    subroutine SaveRunInfo (ModelName, ElapsedSeconds, TotalCPUTime, OutputFile, &
                            WorkCycleElapsed, WorkCycleCPUTime)
                                         
        !Arguments-------------------------------------------------------------
        character(len=*)                :: ModelName
        character(len=*)                :: OutputFile
        real                            :: ElapsedSeconds, TotalCPUTime
        real, optional                  :: WorkCycleElapsed, WorkCycleCPUTime
        integer                         :: ElapsedHours, ElapsedMinutes, ElapsedSecremain
        integer                         :: UnitOutput
        integer                         :: STAT_CALL
        !----------------------------------------------------------------------

        ElapsedHours = INT(Elapsedseconds/3600)
        ElapsedMinutes = INT((ElapsedSeconds-ElapsedHours*3600)/60)
        ElapsedSecremain = INT((ElapsedSeconds-ElapsedMinutes*60-ElapsedHours*3600))

        call UnitsManager (UnitOutput, OPEN_FILE)      
        open(UNIT = UnitOutput, FILE = OutputFile, STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'SaveRunInfo - GlobalData - ERR010'

        write(UnitOutput, *)"-------------------------- MOHID -------------------------"
        write(UnitOutput, *)
        write(UnitOutput, *)"Program "//ModelName//" successfully terminated"
        write(UnitOutput, *)                    
        write(UnitOutput, *)
        write(UnitOutput, 110)ElapsedSeconds, ElapsedHours, ElapsedMinutes,ElapsedSecremain
        write(UnitOutput, 120)TotalCPUTime
        write(UnitOutput, 130)100.*TotalCPUTime/ElapsedSeconds
        if (present(WorkCycleElapsed) .and. present(WorkCycleCPUTime)) then
            write(UnitOutput, 140)WorkCycleElapsed
            write(UnitOutput, 150)WorkCycleCPUTime
            if(WorkCycleElapsed/=.0) then
                write(UnitOutput, 160)100.*WorkCycleCPUTime/WorkCycleElapsed
            endif
        endif
        write(*, *)
        write(*, *)"----------------------------------------------------------"
                
        call UnitsManager (UnitOutput, CLOSE_FILE)

    110 format(1x, "Total Elapsed Time     : ",f14.2," ",i3,"h ",i2,"min ",i2,"s",/)
    120 format(1x, "Total CPU time         : ",f14.2,/)
    130 format(1x, "CPU usage (%)          : ",f14.2,/)
    140 format(1x, "Workcycle Elapsed Time : ",f14.2,/)
    150 format(1x, "Workcycle CPU time     : ",f14.2,/)
    160 format(1x, "Workcycle CPU usage (%): ",f14.2,/)

        !----------------------------------------------------------------------
            
    end subroutine SaveRunInfo

    !--------------------------------------------------------------------------

    integer function GetUsersNumber (iModule, iInstance)

        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
        integer                                     :: iInstance

        !Local-----------------------------------------------------------------

        !Instance users number
        GetUsersNumber = ObjCollector(iModule, iInstance)%Users

    end function GetUsersNumber
    
    !--------------------------------------------------------------------------

    subroutine WriteErrorMessage(keyword, text, STAT)

        !Arguments-------------------------------------------------------------
        character(LEN = *), intent(IN )             :: keyword               
        character(LEN = *), intent(IN )             :: text
        integer, optional,  intent(OUT)             :: STAT    

        !Local-----------------------------------------------------------------
        integer                                     :: LengthText, LengthKeyword
        logical                                     :: OK

        !----------------------------------------------------------------------
        LengthText    = LEN_TRIM(text   )
        LengthKeyword = LEN_TRIM(keyword)

        OK = .true.

#ifdef _GUI_ 
        OK =.false.
#endif

#ifdef _ONLINE_ 
        OK =.false.
#endif
        if (OK) then
            write (ErrorFileID,*    )             
            write (ErrorFileID,*    )             
            write (ErrorFileID,"(A)") text
            write (ErrorFileID, 100 ) keyword
    100     format("Keyword: ", A, " not found"/)
        endif    

        if (present(STAT)) STAT = SUCCESS_
        !----------------------------------------------------------------------
    
    end subroutine WriteErrorMessage

    !--------------------------------------------------------------------------

    subroutine LogKeyWord (keyword, SearchType, ClientModule, Default, CaseSensitive, Value)
    
        !Arguments-------------------------------------------------------------
        character (len=*)                           :: keyword
        integer                                     :: SearchType
        character (len=*)                           :: ClientModule
        character (len=*)                           :: Default
        logical                                     :: CaseSensitive
        character (len=*)                           :: Value

        !Local-----------------------------------------------------------------
        character (len=line_length)                 :: lkeyword
        character (len=line_length)                 :: lSearchType
        character (len=line_length)                 :: lClientModule
        character (len=line_length)                 :: lDefault
        character (len=line_length)                 :: lValue
        character (len=3          )                 :: lCaseSensitive

        logical                                     :: OK

        lkeyword      = adjustl (keyword     )
        lClientModule = adjustl (ClientModule)
        lDefault      = adjustl (Default     )
        lValue        = adjustl (Value       )

        select case (SearchType)
        case (FromFile_)
            lSearchType = "From File"
        case (FromBlock_)
            lSearchType = "From Block"
        case (FromBlockInBlock_)
            lSearchType = "From Block In Block"
        case (FromBlockInBlockInBlock_)
            lSearchType = "From Block In Block In Block"            
        end select

        if (CaseSensitive) then
            lCaseSensitive = 'YES'
        else
            lCaseSensitive = 'NO '
        endif

        OK = .true.

#ifdef _GUI_ 
        OK =.false.
#endif

#ifdef _ONLINE_ 
        OK =.false.
#endif
        if (OK) then
            write  (unit = UsedKeyFileID, fmt = 10)lkeyword, lClientModule, lDefault, lSearchType, lValue, lCaseSensitive
       10   format (1x, a25, 1x, a25, 1x, a25, 1x, a25, 1x, a60, 1x, a3)
        endif

    end subroutine LogKeyWord 


#ifdef _OPENMI_

    !DEC$ IFDEFINED (VF66)
    !dec$ attributes dllexport::GetPropertyNameByID
    !DEC$ ELSE
    !dec$ attributes dllexport,alias:"_GETPROPERTYNAMEBYID"::GetPropertyNameByID
    !DEC$ ENDIF
    logical function GetPropertyNameByID(PropertyID, PropertyName)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: PropertyID
        character(len=*)                            :: PropertyName        
        
        !Local-----------------------------------------------------------------

        PropertyName = GetPropertyName(PropertyID)
        GetPropertyNameByID = .true.

        return

    end function GetPropertyNameByID

    

#endif


end module ModuleGlobalData

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 

