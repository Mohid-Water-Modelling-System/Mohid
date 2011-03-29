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

    implicit none

    public

    !Subroutines---------------------------------------------------------------
    public  ::  CheckPropertyName
    public  ::  GetPropertyName
    public  ::  GetPropertyIDNumber
    private ::      ConstructPropList
    private ::          AddPropList
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
    integer, parameter  :: MaxModules           =  76

#ifdef _INCREASE_MAXINSTANCES
    integer, parameter  :: MaxInstances         = 1000
#else
    integer, parameter  :: MaxInstances         = 500
#endif

    integer, parameter  :: MaxErrorMessages     = 20
    integer             :: NumberOfErrorMessages=0

    integer, parameter  :: StringLength         = 128
    integer, parameter  :: PathLength           = 256


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
   
    !character
    character(LEN = 1), parameter :: space      = char(32)   !" "
    character(LEN = 1), parameter :: dot        = char(46)   !"."
    character(LEN = 1), parameter :: delimiter  = char(58)   !":"
    character(LEN = 1), parameter :: semicolumn = char(59)   !";"
    character(LEN = 1), parameter :: tab        = char(9)    !" "
    character(LEN = 1), parameter :: backslash  = char(92)   !"\"

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


    !Types of grid borders 
    integer, parameter :: ComplexPolygon_  = 1
    integer, parameter :: RotatedRectang_  = 2
    integer, parameter :: Rectang_         = 3


    !Mohid Land ground water flow type
    integer, parameter :: GWFlowToChanByCell_               = 1
    integer, parameter :: GWFlowToChanByLayer_              = 2

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

    !Drifting macroalgae
    integer, parameter :: DriftingMacroAlgae_               = 850


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
    integer, parameter :: CO2AtmosphericPressure_           = 617
    integer, parameter :: O2AtmosphericPressure_            = 618


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

    integer, parameter ::  TransportCapacity_              = 3101 
    integer, parameter ::  TransportCapacityX_             = 3102 
    integer, parameter ::  TransportCapacityY_             = 3103 
    integer, parameter ::  BottomEvolution_                = 3104 
    integer, parameter ::  Newbathymetry_                  = 3105 

    !wave dynamics
    integer, parameter ::  WaveStressX_                    = 3401
    integer, parameter ::  WaveStressY_                    = 3402
    integer, parameter ::  CurrentX_                       = 3403
    integer, parameter ::  CurrentY_                       = 3404
    !!Monocromatic:
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


    integer, parameter :: ConsolidationFlux_                = 9000
    integer, parameter :: Porosity_                         = 9001
    
    !Cohesive sediment fractions - Module Drainage Network
    integer, parameter ::  TSS_                             = 9101
    integer, parameter ::  COHSED_FINE_                     = 9102
    integer, parameter ::  COHSED_MEDIUM_                   = 9103
    integer, parameter ::  COHSED_COARSE_                   = 9104
    integer, parameter ::  VSS_                             = 9111
    
    !PhreeqC properties ------------------------------------------
    
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
    
    !Pesticides
    integer, parameter :: GenericPesticide_1_              = 15001
    integer, parameter :: GenericPesticide_2_              = 15002
    integer, parameter :: GenericPesticide_3_              = 15003

    !Spatial emission discharge
    integer, parameter :: DischPoint_                       = 1
    integer, parameter :: DischLine_                        = 2
    integer, parameter :: DischPolygon_                     = 3
                                                            
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

    !Other properties
    integer, parameter :: SolEC_                            = 14000 !Solution electrical conductivity

!_______________________________________________________________________________________________

    !Other Properties
    character(StringLength), private, parameter :: Char_SolEC                = 'solution electrical conductivity'
      
    !Name of PhreeqC properties
    character(StringLength), private, parameter :: Char_SolutionMagnesium    = 'solution magnesium'
    character(StringLength), private, parameter :: Char_SolutionCalcium      = 'solution calcium'
    character(StringLength), private, parameter :: Char_SolutionSodium       = 'solution sodium'
    character(StringLength), private, parameter :: Char_SolutionNitrogenGas  = 'solution nitrogen gas'
    character(StringLength), private, parameter :: Char_SolutionOxygenGas    = 'solution oxygen gas'
    character(StringLength), private, parameter :: Char_SolutionAmmonia      = 'solution ammonia'
    character(StringLength), private, parameter :: Char_SolutionNitrate      = 'solution nitrate'
    character(StringLength), private, parameter :: Char_SolutionNitrite      = 'solution nitrite'
    character(StringLength), private, parameter :: Char_SolutionChlorine     = 'solution chlorine'
    character(StringLength), private, parameter :: Char_GasN2                = 'gas n2'
    character(StringLength), private, parameter :: Char_GasCO2               = 'gas co2'
    character(StringLength), private, parameter :: Char_pE                   = 'pE'
    character(StringLength), private, parameter :: Char_eCaX2                = 'CaX2'
    character(StringLength), private, parameter :: Char_eMgX2                = 'MgX2'
    character(StringLength), private, parameter :: Char_eNaX                 = 'NaX'
    character(StringLength), private, parameter :: Char_eKX                  = 'KX'
    character(StringLength), private, parameter :: Char_eNH4X                = 'AmmoniaX' !'NH4X'
    character(StringLength), private, parameter :: Char_sCa2                 = 'Ca+2'
    character(StringLength), private, parameter :: Char_sCaOH                = 'CaOH+'
    character(StringLength), private, parameter :: Char_sH2                  = 'H2'
    character(StringLength), private, parameter :: Char_sMg2                 = 'Mg+2'
    character(StringLength), private, parameter :: Char_sMgOH                = 'MgOH+'
    character(StringLength), private, parameter :: Char_sNH4                 = 'Ammonia+' !'NH4+'
    character(StringLength), private, parameter :: Char_sNH3                 = 'NH3'
    character(StringLength), private, parameter :: Char_sN2                  = 'N2'
    character(StringLength), private, parameter :: Char_sNO2                 = 'NO2-'
    character(StringLength), private, parameter :: Char_sNO3                 = 'NO3-'
    character(StringLength), private, parameter :: Char_sNa                  = 'Na+'
    character(StringLength), private, parameter :: Char_sNaOH                = 'NaOH+'
    character(StringLength), private, parameter :: Char_sO2                  = 'O2'
    character(StringLength), private, parameter :: Char_msmSolutionCalcium   = 'solution calcium mass'
    character(StringLength), private, parameter :: Char_msmSolutionMagnesium = 'solution magnesium mass'
    character(StringLength), private, parameter :: Char_msmSolutionSodium    = 'solution sodium mass'
    character(StringLength), private, parameter :: Char_msmSolutionAmmonia   = 'solution ammonia mass'
    character(StringLength), private, parameter :: Char_Calcite              = 'calcite'
    character(StringLength), private, parameter :: Char_Dolomite             = 'dolomite'
    character(StringLength), private, parameter :: Char_Aragonite            = 'aragonite'
    character(StringLength), private, parameter :: Char_Halite               = 'halite' 
    character(StringLength), private, parameter :: Char_KFeldspar            = 'k-feldspar'
    character(StringLength), private, parameter :: Char_SolutionCarbon       = 'solution carbon' 
    character(StringLength), private, parameter :: Char_SolutionPotassium    = 'solution potassium'     
    character(StringLength), private, parameter :: Char_SolutionAluminium    = 'solution aluminium' 
    character(StringLength), private, parameter :: Char_SolutionSilicium     = 'solution silicium' 
    character(StringLength), private, parameter :: Char_RainMagnesium        = 'rain magnesium' 
    character(StringLength), private, parameter :: Char_RainCalcium          = 'rain calcium' 
    character(StringLength), private, parameter :: Char_RainSodium           = 'rain sodium' 
    character(StringLength), private, parameter :: Char_RainChlorine         = 'rain chlorine'
    character(StringLength), private, parameter :: Char_RainAmmonia          = 'rain ammonia'


!_______________________________________________________________________________________________

      
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


    character(StringLength), private, parameter :: Char_Nitrite              = 'nitrite'
    character(StringLength), private, parameter :: Char_BOD                  = 'biochemical oxygen demand'
    character(StringLength), private, parameter :: Char_Cohesive_Sediment    = 'cohesive sediment'
    character(StringLength), private, parameter :: Char_Fecal_Coliforms      = 'fecal coliforms'
    character(StringLength), private, parameter :: Char_E_Coli               = 'escherichia coli'
    character(StringLength), private, parameter :: Char_T90                  = 'T90'
    character(StringLength), private, parameter :: Char_T90_E_Coli           = 'T90 e.coli'

    character(StringLength), private, parameter :: Char_Oil                  = 'oil'
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
    

    character(StringLength), private, parameter :: Char_GrossProd            = 'grossprod'
    character(StringLength), private, parameter :: Char_NutrientLim          = 'nutrientlim'
    character(StringLength), private, parameter :: Char_NLim                 = 'nitrogenlim'
    character(StringLength), private, parameter :: Char_PLim                 = 'phosphoruslim'
    character(StringLength), private, parameter :: Char_LightLim             = 'lightlim'
    character(StringLength), private, parameter :: Char_TemperatureLim       = 'temperaturelim'
    character(StringLength), private, parameter :: Char_SalinityLim          = 'salinitylim' 
    character(StringLength), private, parameter :: Char_NetProd              = 'netprod'

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


    ! guillaume nogueira
    character(StringLength), private, parameter :: Char_AltimLevelAnalyzed_         = 'water level analyzed for altimetry'
    character(StringLength), private, parameter :: Char_AltimTemperatureAnalyzed_   = 'temperature analyzed for altimetry'
    character(StringLength), private, parameter :: Char_AltimSalinityAnalyzed_      = 'salinity analyzed for altimetry'
    character(StringLength), private, parameter :: Char_AltimLevelToAssimilate_     = 'water level for altimetry assimilation'
    character(StringLength), private, parameter :: Char_VarianceFieldToAssimilate_  = 'variance field for assimilation'

    character(StringLength), private, parameter :: Char_WaterFluxX_          = 'water flux X'
    character(StringLength), private, parameter :: Char_WaterFluxY_          = 'water flux Y'
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
    character(StringLength), private, parameter :: Char_WindShearVelocity        = 'wind shear velocity'
    character(StringLength), private, parameter :: Char_SurfaceRadiation         = 'surface radiation'
    character(StringLength), private, parameter :: Char_WindStressX              = 'wind stress X'
    character(StringLength), private, parameter :: Char_WindStressY              = 'wind stress Y'
    character(StringLength), private, parameter :: Char_SurfaceWaterFlux         = 'surface water flux'
    character(StringLength), private, parameter :: Char_NonSolarFlux             = 'non solar flux'
    character(StringLength), private, parameter :: Char_TurbulentKineticEnergy   = 'turbulent kinetic energy'
    character(StringLength), private, parameter :: Char_Albedo                   = 'albedo'

    !Atmosphere
    character(StringLength), private, parameter :: Char_WindVelocityX            = 'wind velocity X'
    character(StringLength), private, parameter :: Char_WindVelocityY            = 'wind velocity Y'
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

    character(StringLength), private, parameter :: Char_MeanSeaLevelPressure     = 'mean sea level pressure'
    character(StringLength), private, parameter :: Char_WindModulus              = 'wind modulus'
    character(StringLength), private, parameter :: Char_WindDirection            = 'wind direction'
    character(StringLength), private, parameter :: Char_SpecificHumidity         = 'specific humidity'

    character(StringLength), private, parameter :: Char_CO2AtmosphericPressure   = 'CO2 atmospheric pressure'
    character(StringLength), private, parameter :: Char_O2AtmosphericPressure    = 'O2 atmospheric pressure'

    
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

    !wave dynamics
    character(StringLength), private, parameter :: Char_WaveStressX              = 'wave stress X'
    character(StringLength), private, parameter :: Char_WaveStressY              = 'wave stress Y'   
    character(StringLength), private, parameter :: Char_CurrentX                 = 'Current X'
    character(StringLength), private, parameter :: Char_CurrentY                 = 'Current Y'
    !!Monocromatic:
    character(StringLength), private, parameter :: Char_WaveAmplitude            = 'wave amplitude'
    character(StringLength), private, parameter :: Char_WavePeriod               = 'wave period'
    character(StringLength), private, parameter :: Char_WaveDirection            = 'wave direction'
    !!Statistical wave parametres (WW3,SWAN)
    character(StringLength), private, parameter :: Char_SignificantWaveHeight    = 'significant wave height'
    character(StringLength), private, parameter :: Char_MeanWaveLength           = 'mean wave length'
    character(StringLength), private, parameter :: Char_MeanWavePeriod           = 'mean wave period'
    character(StringLength), private, parameter :: Char_MeanWaveDirection        = 'mean wave direction'
    character(StringLength), private, parameter :: Char_MeanDirectionalSpread    = 'mean directional spread'
    character(StringLength), private, parameter :: Char_PeakFrequency            = 'peak frequency'
    character(StringLength), private, parameter :: Char_PeakDirection            = 'peak direction'
    character(StringLength), private, parameter :: Char_WindSeaPeakFrequency     = 'wind sea peak frequency'
    character(StringLength), private, parameter :: Char_WindSeaPeakDirection     = 'wind sea peak direction'
    character(StringLength), private, parameter :: Char_WaveSwellHeight          = 'wave swell height'

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

    !Cohesive Fractions - Drainage Network
    character(StringLength), private, parameter :: Char_TSS                      = 'TSS'
    character(StringLength), private, parameter :: Char_Cohsed_fine              = 'cohesive sediment fine'
    character(StringLength), private, parameter :: Char_Cohsed_medium            = 'cohesive sediment medium'
    character(StringLength), private, parameter :: Char_Cohsed_coarse            = 'cohesive sediment coarse'
    character(StringLength), private, parameter :: Char_VSS                      = 'VSS'

    !ChainReactions
    character(StringLength), private, parameter :: Char_SoilVolumetricDensity    = 'soil volumetric density'
    
    !Pesticides
    character(StringLength), private, parameter :: Char_GenericPesticide_1       = 'generic pesticide 1'
    character(StringLength), private, parameter :: Char_GenericPesticide_2       = 'generic pesticide 2'
    character(StringLength), private, parameter :: Char_GenericPesticide_3       = 'generic pesticide 3'

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
    real(8), parameter  :: Pi = 3.1415926535897932384626433832795

    !Zero Degrees Kelvin
    real, parameter     :: AbsoluteZero = 273.15

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
    integer, parameter :: WSConstant = 1, SPMFunction = 2

    !Advection 1D parameters
    integer, parameter :: UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4, CentralDif = 5, LeapFrog = 6
    integer, parameter :: MinMod = 1, VanLeer = 2, Muscl = 3, Superbee = 4, PDM = 5
    real,    parameter :: MinValue = 1.e-16

    !Transport Parameters
    integer, parameter :: velocityX = 1, velocityY = 2, velocityZ = 3, massProperty = 4
    integer, parameter :: NearestNeighbour = 1, centered = 2 !other interpolation methods here 

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

    type T_Size1D
        integer                 :: ILB            = null_int
        integer                 :: IUB            = null_int
    end type T_Size1D

    type T_Size2D
        integer                 :: ILB            = null_int
        integer                 :: IUB            = null_int
        integer                 :: JLB            = null_int
        integer                 :: JUB            = null_int
    end type T_Size2D

    type T_Size3D
        integer                 :: ILB            = null_int
        integer                 :: IUB            = null_int
        integer                 :: JLB            = null_int
        integer                 :: JUB            = null_int
        integer                 :: KLB            = null_int
        integer                 :: KUB            = null_int
    end type T_Size3D

    type T_PropertyID
        character(StringLength) :: Name              = null_str
        character(StringLength) :: Units             = null_str
        character(StringLength) :: Description       = null_str
        integer                 :: IDNumber          = null_int    
        integer                 :: ObjFillMatrix     = 0
        logical                 :: SolutionFromFile  = OFF
    end type T_PropertyID

    type T_Instance
        integer                 :: ID         = 0
        integer                 :: Users      = 0
        integer                 :: Readers    = 0
        logical                 :: Read_Lock  = IDLE
    end type T_Instance

    type T_Module
        integer                 :: ID
        character(StringLength) :: Name
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
        T_Module(mRUNOFFPROPERTIES_      , "RunoffProperties"),      T_Module(mCHAINREACTIONS_         , "ChainReactions") /)

    !Variables
    logical, dimension(MaxModules)                                  :: RegisteredModules = .false.
    character(LEN=PathLength)                                       :: FilesName = 'nomfich.dat'
    logical                                                         :: MonitorPerformance
    logical                                                         :: MonitorDT
    integer                                                         :: UnitDT
    character(LEN=StringLength), dimension(:), pointer              :: PropNameList
    integer,                     dimension(:), pointer              :: PropNumberList
    integer                                                         :: PropertiesNumber
    integer, private                                                :: ErrorFileID     = 0
    integer, private                                                :: UsedKeyFileID   = 0
    integer, private                                                :: LogFileID       = 0
    character(LEN=1024)                                             :: OnLineString
    character(StringLength),     dimension(MaxErrorMessages)        :: ErrorMessagesStack
      
    type (T_Instance), dimension (MaxModules, MaxInstances), save   :: ObjCollector
    private :: ObjCollector

    contains 

    !--------------------------------------------------------------------------

    logical function ModuleIsRegistered (iModule)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iModule
    
        ModuleIsRegistered = RegisteredModules(iModule)
        
    end function ModuleIsRegistered
    
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

    logical function CheckPropertyName (PropertyName, Number)

        !Arguments-------------------------------------------------------------
        character(len=*), intent (IN)               :: PropertyName
        integer,          intent (OUT), optional    :: Number

        !Local-----------------------------------------------------------------
        integer :: i

        !----------------------------------------------------------------------

        CheckPropertyName = .false.

        call ConstructPropList

        do i=1, PropertiesNumber

            if (PropertyName == PropNameList(i)) then

                if (present(Number)) Number = PropNumberList(i)
                CheckPropertyName = .TRUE.

            endif

        enddo

        !----------------------------------------------------------------------

    end function  CheckPropertyName

    !--------------------------------------------------------------------------

    character (len=StringLength) function GetPropertyName (Number)

        !Arguments-------------------------------------------------------------
        integer,          intent (IN )              :: Number

        !Local-----------------------------------------------------------------
        integer :: i

        !----------------------------------------------------------------------

        call ConstructPropList

        do i=1, PropertiesNumber

            if (Number == PropNumberList(i)) then

                GetPropertyName = PropNameList(i)

            endif

        enddo

        !----------------------------------------------------------------------

    end function GetPropertyName

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    integer function GetPropertyIDNumber (PropertyName)

        !Arguments-------------------------------------------------------------
        character(len=*), intent (IN )              :: PropertyName

        !Local-----------------------------------------------------------------
        integer :: i

        !----------------------------------------------------------------------
        
        GetPropertyIDNumber = UNKNOWN_

        call ConstructPropList

        do i=1, PropertiesNumber

            if (trim(PropertyName) == trim(PropNameList(i))) then

                GetPropertyIDNumber = PropNumberList(i)

            endif

        enddo

        if(GetPropertyIDNumber == UNKNOWN_)then
            write(*,*)'Unknown property: ', PropertyName
            stop 'GetPropertyIDNumber - ModuleGlobalData - ERR010'
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
         
            call AddPropList (Nitrite_,                 Char_Nitrite,                   ListNumber)
            call AddPropList (BOD_,                     Char_BOD,                       ListNumber)
            call AddPropList (Cohesive_Sediment_,       Char_Cohesive_Sediment,         ListNumber)
            call AddPropList (Fecal_Coliforms_,         Char_Fecal_Coliforms,           ListNumber)
            call AddPropList (E_Coli_,                  Char_E_Coli,                    ListNumber)
            call AddPropList (T90_,                     Char_T90,                       ListNumber)
            call AddPropList (T90_E_Coli_,              Char_T90_E_Coli,                ListNumber)
            call AddPropList (Oil_,                     Char_Oil,                       ListNumber)
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
            
            call AddPropList (GenericProperty_,         Char_GenericProperty,           ListNumber)
            call AddPropList (GrossProd_,               Char_GrossProd,                 ListNumber)
            call AddPropList (NetProd_,                 Char_NetProd,                   ListNumber)
            call AddPropList (NutrientLim_,             Char_NutrientLim,               ListNumber)
            call AddPropList (NLim_,                    Char_NLim,                      ListNumber)
            call AddPropList (PLim_,                    Char_PLim,                      ListNumber)
            call AddPropList (LightLim_ ,               Char_LightLim ,                 ListNumber)
            call AddPropList (TemperatureLim_,          Char_TemperatureLim,            ListNumber)
            call AddPropList (SalinityLim_,             Char_SalinityLim,               ListNumber)
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
            call AddPropList (CoriolisX_,               Char_CoriolisX_,                ListNumber)
            call AddPropList (BaroclinicForceX_,        Char_BaroclinicForceX_,         ListNumber)
            call AddPropList (HorizontalTransportX_,    Char_HorizontalTransportX_,     ListNumber)
            call AddPropList (CoriolisY_,               Char_CoriolisY_,                ListNumber)
            call AddPropList (BaroclinicForceY_  ,      Char_BaroclinicForceY_  ,       ListNumber)
            call AddPropList (HorizontalTransportY_,    Char_HorizontalTransportY_,     ListNumber)
            call AddPropList (BarotropicVelocityU_ ,    Char_BarotropicVelocityU_ ,     ListNumber)
            call AddPropList (BarotropicVelocityV_ ,    Char_BarotropicVelocityV_ ,     ListNumber)


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
            call AddPropList (ObstacleDragCoef_ ,       Char_ObstacleDragCoef ,         ListNumber)

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
            call AddPropList (SurfaceWaterFlux_ ,          Char_SurfaceWaterFlux ,         ListNumber)
            call AddPropList (NonSolarFlux_ ,              Char_NonSolarFlux ,             ListNumber)
            call AddPropList (TurbulentKineticEnergy_,     Char_TurbulentKineticEnergy,    ListNumber)
            call AddPropList (CarbonDioxideFlux_ ,         Char_CarbonDioxideFlux,         ListNumber)
            call AddPropList (SpecificCarbonDioxideFlux_ , Char_SpecificCarbonDioxideFlux, ListNumber)
            call AddPropList (Albedo_ ,                    Char_Albedo,                    ListNumber)


            call AddPropList (WindVelocityX_,           Char_WindVelocityX          ,      ListNumber)
            call AddPropList (WindVelocityY_,           Char_WindVelocityY          ,      ListNumber)
            call AddPropList (SolarRadiation_,          Char_SolarRadiation         ,      ListNumber)
            call AddPropList (Precipitation_,           Char_Precipitation          ,      ListNumber)
            call AddPropList (AtmosphericPressure_,     Char_AtmosphericPressure    ,      ListNumber)
            call AddPropList (CO2AtmosphericPressure_,  Char_CO2AtmosphericPressure ,      ListNumber)
            call AddPropList (O2AtmosphericPressure_,   Char_O2AtmosphericPressure  ,      ListNumber)
            call AddPropList (AirTemperature_,          Char_AirTemperature         ,      ListNumber)
            call AddPropList (RelativeHumidity_,        Char_RelativeHumidity       ,      ListNumber)
            call AddPropList (WindModulos_,             Char_WindModulos            ,      ListNumber)
            call AddPropList (WindAngle_,               Char_WindAngle              ,      ListNumber)
            call AddPropList (CloudCover_,              Char_CloudCover             ,      ListNumber)
            call AddPropList (Irrigation_,              Char_Irrigation             ,      ListNumber)
            call AddPropList (SunHours_ ,               Char_SunHours               ,      ListNumber)
            call AddPropList (ATMTransmitivity_ ,       Char_ATMTransmitivity       ,      ListNumber)

            call AddPropList (MeanSeaLevelPressure_ ,   Char_MeanSeaLevelPressure,      ListNumber)
            call AddPropList (WindModulus_ ,            Char_WindModulus         ,      ListNumber)
            call AddPropList (WindDirection_ ,          Char_WindDirection       ,      ListNumber)

         
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

            call AddPropList (TransportCapacity_,       Char_TransportCapacity,          ListNumber)
            call AddPropList (TransportCapacityX_,      Char_TransportCapacityX,         ListNumber)
            call AddPropList (TransportCapacityY_,      Char_TransportCapacityY,         ListNumber)
            call AddPropList (BottomEvolution_,         Char_BottomEvolution  ,          ListNumber)
            call AddPropList (Newbathymetry_,           Char_Newbathymetry    ,          ListNumber)

            !wave dynamics
            call AddPropList (WaveStressX_,             Char_WaveStressX,                ListNumber)
            call AddPropList (WaveStressY_,             Char_WaveStressY,                ListNumber)
            call AddPropList (CurrentX_,                Char_CurrentX,                   ListNumber)
            call AddPropList (CurrentY_,                Char_CurrentY,                   ListNumber)

            !!Monocromatic:
            call AddPropList (WaveAmplitude_,           Char_WaveAmplitude,              ListNumber)
            call AddPropList (WavePeriod_,              Char_WavePeriod,                 ListNumber)
            call AddPropList (WaveDirection_,           Char_WaveDirection,              ListNumber)
            !!Statistical wave parametres (WW3,SWAN)
            call AddPropList (SignificantWaveHeight_,   Char_SignificantWaveHeight,      ListNumber)
            call AddPropList (MeanWaveLength_,          Char_MeanWaveLength,             ListNumber)
            call AddPropList (MeanWavePeriod_,          Char_MeanWavePeriod,             ListNumber)
            call AddPropList (MeanWaveDirection_,       Char_MeanWaveDirection,          ListNumber)
            call AddPropList (MeanDirectionalSpread_,   Char_MeanDirectionalSpread,      ListNumber)
            call AddPropList (PeakFrequency_,           Char_PeakFrequency,              ListNumber)
            call AddPropList (PeakDirection_,           Char_PeakDirection,              ListNumber)
            call AddPropList (WindSeaPeakFrequency_,    Char_WindSeaPeakFrequency,       ListNumber)
            call AddPropList (WindSeaPeakDirection_,    Char_WindSeaPeakDirection,       ListNumber)
            call AddPropList (WaveSwellHeight_     ,    Char_WaveSwellHeight,            ListNumber)
            
            call AddPropList (ConsolidationFlux_,       Char_ConsolidationFlux,          ListNumber)
            call AddPropList (Porosity_,                Char_Porosity,                   ListNumber)
            
            call AddPropList (ShearStress_,             Char_ShearStress_,               ListNumber)

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

            call AddPropList (TSS_,                     Char_TSS,                        ListNumber)
            call AddPropList (COHSED_FINE_,             Char_Cohsed_fine,                ListNumber)
            call AddPropList (COHSED_MEDIUM_,           Char_Cohsed_medium,              ListNumber)
            call AddPropList (COHSED_COARSE_,           Char_Cohsed_coarse,              ListNumber)
            call AddPropList (VSS_,                     Char_VSS,                        ListNumber)

            !PhreeqC temporary code for tests
            call AddPropList (SolutionMagnesium_,       Char_SolutionMagnesium,          ListNumber)
            call AddPropList (SolutionCalcium_,         Char_SolutionCalcium,            ListNumber)
            call AddPropList (SolutionSodium_,          Char_SolutionSodium,             ListNumber)
            call AddPropList (SolutionNitrogenGas_,     Char_SolutionNitrogenGas,        ListNumber)
            call AddPropList (SolutionOxygenGas_,       Char_SolutionOxygenGas,          ListNumber)
            call AddPropList (SolutionAmmonia_,         Char_SolutionAmmonia,            ListNumber)
            call AddPropList (SolutionNitrate_,         Char_SolutionNitrate,            ListNumber)
            call AddPropList (SolutionNitrite_,         Char_SolutionNitrite,            ListNumber)
            call AddPropList (SolutionChlorine_,        Char_SolutionChlorine,           ListNumber)
            call AddPropList (GasN2_,                   Char_GasN2,                      ListNumber)
            call AddPropList (GasCO2_,                  Char_GasCO2,                     ListNumber)
            call AddPropList (pE_,                      Char_pE,                         ListNumber)
            call AddPropList (eCaX2_,                   Char_eCaX2,                      ListNumber)
            call AddPropList (eMgX2_,                   Char_eMgX2,                      ListNumber)
            call AddPropList (eNaX_,                    Char_eNaX,                       ListNumber)
            call AddPropList (eNH4X_,                   Char_eNH4X,                      ListNumber)
            call AddPropList (eKX_,                     Char_eKX,                        ListNumber)
            call AddPropList (sCa2_,                    Char_sCa2,                       ListNumber)
            call AddPropList (sCaOH_,                   Char_sCaOH,                      ListNumber)
            call AddPropList (sH2_,                     Char_sH2,                        ListNumber)
            call AddPropList (sMg2_,                    Char_sMg2,                       ListNumber)
            call AddPropList (sMgOH_,                   Char_sMgOH,                      ListNumber)
            call AddPropList (sNH4_,                    Char_sNH4,                       ListNumber)
            call AddPropList (sNH3_,                    Char_sNH3,                       ListNumber)
            call AddPropList (sN2_,                     Char_sN2,                        ListNumber)
            call AddPropList (sNO2_,                    Char_sNO2,                       ListNumber)
            call AddPropList (sNO3_,                    Char_sNO3,                       ListNumber)
            call AddPropList (sNa_,                     Char_sNa,                        ListNumber)
            call AddPropList (sNaOH_,                   Char_sNaOH,                      ListNumber)
            call AddPropList (sO2_,                     Char_sO2,                        ListNumber)
            call AddPropList (msmSolutionCalcium_,      Char_msmSolutionCalcium,         ListNumber)
            call AddPropList (msmSolutionMagnesium_,    Char_msmSolutionMagnesium,       ListNumber)
            call AddPropList (msmSolutionSodium_,       Char_msmSolutionSodium,          ListNumber)
            call AddPropList (msmSolutionAmmonia_,      Char_msmSolutionAmmonia,         ListNumber)
            call AddPropList (Calcite_,                 Char_Calcite,                    ListNumber)
            call AddPropList (Dolomite_,                Char_Dolomite,                   ListNumber)
            call AddPropList (Aragonite_,               Char_Aragonite,                  ListNumber)
            call AddPropList (Halite_,                  Char_Halite,                     ListNumber)
            call AddPropList (KFeldspar_,               Char_KFeldspar,                  ListNumber)
            call AddPropList (SolutionCarbon_,          Char_SolutionCarbon,             ListNumber) 
            call AddPropList (SolutionPotassium_,       Char_SolutionPotassium,          ListNumber) 
            call AddPropList (SolutionAluminium_,       Char_SolutionAluminium,          ListNumber) 
            call AddPropList (SolutionSilicium_,        Char_SolutionSilicium,           ListNumber) 
            call AddPropList (RainMagnesium_,           Char_RainMagnesium,              ListNumber) 
            call AddPropList (RainCalcium_,             Char_RainCalcium,                ListNumber) 
            call AddPropList (RainSodium_,              Char_RainSodium,                 ListNumber)
            call AddPropList (RainChlorine_,            Char_RainChlorine,               ListNumber)  
            call AddPropList (RainAmmonia_,             Char_RainAmmonia,                ListNumber)  
            !END of PhreeqC temporary code for tests
            call AddPropList (SoilVolumetricDensity_,   Char_SoilVolumetricDensity,      ListNumber)  
            call AddPropList (SolEC_,                   Char_SolEC,                      ListNumber)
            call AddPropList (GenericPesticide_1_,      Char_GenericPesticide_1,         ListNumber)
            call AddPropList (GenericPesticide_2_,      Char_GenericPesticide_2,         ListNumber)
            call AddPropList (GenericPesticide_3_,      Char_GenericPesticide_3,         ListNumber) 
            
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
 
                deallocate   (AuxInt, AuxChar)
                allocate     (AuxChar(ListNumber + 1), AuxInt(ListNumber + 1))       

                ListNumber = ListNumber + 1

            endif     

            PropNameList  (ListNumber) = PropName
            PropNumberList(ListNumber) = PropNumber

            AuxChar      (:)          = PropNameList  (:)
            AuxInt       (:)          = PropNumberList(:)

        else       

            deallocate (AuxChar, AuxInt)

        endif      


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
            (Property == SolubilizingP_         ) .OR.  (Property == ParticulateMetal_      ) .OR.          &
            (Property == ParticulateCadmium_    ) .OR.  (Property == ParticulateCopper_     ) .OR.          &
            (Property == ParticulateLead_       ) .OR.  (Property == ParticulateZinc_       ) .OR.          &
            (Property == ParticulateMercury_    )    ) then

            Check_Particulate_Property = .TRUE.    
        
        else

            Check_Particulate_Property = .FALSE.

        end if cd1

    end function Check_Particulate_Property

        
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

    subroutine UnitsManager(UNIT, OPENCLOSE, STAT)

        !Arguments-------------------------------------------------------------
        integer,           intent(INOUT)    :: UNIT
        integer,           intent(IN   )    :: OPENCLOSE
        integer, optional, intent(OUT  )    :: STAT

        !External--------------------------------------------------------------
        logical                             :: opened
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

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
                close(UNIT, IOSTAT = STAT_CALL)

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

    subroutine WriteDTLog (ModuleName, iter, DT)

        !Arguments-------------------------------------------------------------
        character(len=*)                :: ModuleName
        integer                         :: iter
        real                            :: DT

        if (MonitorDT) then
            write (UnitDT, '(A25, I10,F12.4)') ModuleName, iter, DT 
        end if

    end subroutine WriteDTLog

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

#ifndef _OPENMI_
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
#ifndef _OPENMI_
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

