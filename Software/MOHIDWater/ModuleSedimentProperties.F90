!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : SedimentProperties
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute sediment properties processes
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
!DataFile
!
!   OUTPUT_TIME                 : sec. sec. sec.    -           !Output Time
!   TIME_SERIE_LOCATION         : char              -           !Path to time serie location file
!   BOXFLUXES                   : char              -           !If specified computes box integration
!                                                               !based on boxes file defined by this keyword
!<beginproperty>
!   NAME                        : char              -           !Property name
!   UNITS                       : char              -           !Property units
!   DESCRIPTION                 : char              -           !Small description of the property
!   See module FillMatrix       : -                 -           !Initialization of concentration values
!                                                               !NOTE: Dissolved concentration   = mass/porewater_volume
!                                                               !      Particulate concentration = mass/drysediment_volume
!   PARTICULATE                 : 0/1               [0]         !Property physical state: 0 - Dissolved ; 1 - Particulate
!   OLD                         : 0/1               [0]         !Initialization from previous run (overrides FillMatrix)
!   MIN_VALUE                   : real              -           !Minimum allowed value of property concentration 
!   IS_COEF                     : real              1e-3        !Conversion factor to I.S. units (1e-3 = mg/l)
!   ADVECTION_DIFFUSION         : 0/1               [0]         !Compute advection-diffusion
!       DIFFUSION_METHOD        : 1/2               [1]         !Method to compute diffusion coefficeient correction
!                                                               !for the sediments. 1 - Berner, 1980 ; 2 - Soetaert, 1996
!       MOLECULAR_DIFF_COEF     : real              [0 m2/s]    !Infinite dilution molecular diffusion coefficient
!   BIOTURBATION                : 0/1               [0]         !Compute bioturbation
!       BIOTURBATION_COEF       : real              [0 m2/s]    !Bioturbation diffusion coefficient 
!       BIOTURBATION_DEPTH      : real              [0.1 m]     !Depth till which bioturbation diffusion is constant (m)
!       BIOTURBATION_DECAY_COEF : real              [0.01]      !Decay factor to compute decay of bioturbation effect
!
!   PARTITION                   : 0/1               [0]         !Compute partition between dissolved-particulate phases
!       PARTITION_COUPLE        : char              []          !Name of the property (oposite phase) to compute partition
!       <<begin_partition_fraction>>
!       See module FillMatrix   : -                 []          !Initialization of partition fractions (values between 0-1)
!       <<end_partition_fraction>>
!
!       <<begin_partition_rate>>
!       See module FillMatrix   : -                 [s-1]       !Initialization of partition rates values
!       <<end_partition_rate>>'

!   SEDIMENT_QUALITY            : 0/1               [0]         !Compute sediment quality processes
!   SURFACE_FLUXES              : 0/1               [0]         !Compute fluxes from sediment surface
!   DT_INTERVAL                 : real              [ModelDT]   !Property evolution time step
!   TIME_SERIE                  : 0/1               [0]         !Ouputs results in time series  
!   BOX_TIME_SERIE              : 0/1               [0]         !Ouputs results in box time series
!   OUTPUT_HDF                  : 0/1               [0]         !Ouputs results in HDF5 format
!<endproperty>
!
!<begin_sed_dry_density>
!   See module FillMatrix       : -                 -           !Initialization of sediment dry density values
!<end_sed_dry_density>
!     
!<begin_vert_diff_coef>
!   See module FillMatrix       : -                 -           !Initialization of sediment vertical diffusivity
!<end_vert_diff_coef>
    
!<begin_horiz_diff_coef>
!   See module FillMatrix       : -                 -           !Initialization of horizontal vertical diffusivity
!<end_horiz_diff_coef>

!<beginrate>
!   NAME                        : char              -           !Name of the rate to perform output
!   DESCRIPTION                 : char              -           !Description of the rate to perform output
!   FIRSTPROP                   : char              -           !Name of the first property involved in the rate
!   SECONDPROP                  : char              -           !Name of the second property involved in the rate
!<endrate>


!review this

Module ModuleSedimentProperties

    use ModuleGlobalData
    use ModuleTime
    use ModuleFunctions,            only: ConstructPropertyID , SetMatrixValue
    use ModuleHDF5,                 only: ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,          &
                                          HDF5WriteData, HDF5FlushMemory, HDF5ReadData, KillHDF5
    use ModuleEnterData,            only: ReadFileName, ConstructEnterData, GetData,                &
                                          ExtractBlockFromBuffer, Block_Unlock, RewindBuffer,       &
                                          GetOutPutTime, KillEnterData, ExtractBlockFromBlock,      &
                                          RewindBlock           
    use ModuleTimeSerie,            only: StartTimeSerie, WriteTimeSerie, KillTimeSerie
    use ModuleGridData,             only: GetGridData, UngetGridData
    use ModuleHorizontalMap,        only: GetOpenPoints2D, GetWaterPoints2D, GetBoundaries,         &
                                          UnGetHorizontalMap       
    use ModuleHorizontalGrid,       only: GetHorizontalGrid, WriteHorizontalGrid, GetGridCellArea,  &
                                          UnGetHorizontalGrid
    use ModuleGeometry,             only: GetGeometrySize, UnGetGeometry, GetGeometryVolumes,       &
                                          GetGeometryDistances, GetGeometryKtop
    use ModuleMap,                  only: GetWaterPoints3D, GetOpenPoints3D, GetComputeFaces3D,     &
                                          GetLandPoints3D, UngetMap             
    use ModuleBoxDif,               only: StartBoxDif, GetBoxes, GetNumberOfBoxes, UngetBoxDif,     &
                                          BoxDif, KillBoxDif        
    use ModuleAdvectionDiffusion,   only: GetBoundaryConditionList
    use ModuleConsolidation,        only: UngetConsolidation, GetConsolidationWaterPercentage,      &
                                          GetConsolidationDrySedVolume, GetConsolidationTortuosity, &
                                          GetConsolidationWaterVolume, GetConsolidationWaterFluxes, &
                                          GetConsolidationVelocityZ, GetConsolidationVelocityXY,    &
                                          GetConsolidationDepth, GetConsolidationPorosity,          &
                                          SetSedimentDryDensity, GetConsolidationKTopState,         &
                                          GetConsolidationMinThickness
    use ModuleInterface,            only: ConstructInterface, Modify_Interface, GetRateFlux, KillInterface
    use ModuleFillMatrix,           only: ConstructFillMatrix, GetDefaultValue, KillFillMatrix
    use ModuleLUD,                  only: StartLUD, LUD, KillLUD
    use ModuleStopWatch,            only: StartWatch, StopWatch

    implicit none 

    private

    !subroutines---------------------------------------------------------------

    !Constructor
    public  :: Construct_SedimentProperties
    private ::      AllocateInstance
    private ::      Construct_GlobalVariables
    private ::          ReadSedimentPropertiesFilesName
    private ::      Construct_PropertyList
    private ::          Construct_Property  
    private ::              Construct_PropertyValues
    private ::                  ReadOldConcBoundariesHDF
    private ::              Construct_PropertyEvolution
    private ::                  Read_Advec_Difus_Parameters
    private ::                  Read_Bioturbation_Parameters
    private ::                  Read_Partition_Parameters
    private ::                  Construct_Property_Diffusivity
    private ::              Construct_PropertyState
    private ::              Construct_PropertyOutPut
    private ::          Add_Property
    private ::      Construct_SedimentRateList
    private ::          Construct_SedimentRate 
    private ::              Construct_SedimentRateID
    private ::              Construct_SedimentRateValues
    private ::          Add_SedimentRate
    private ::      Construct_Sub_Modules
    private ::          CoupleAdvectionDiffusion
    private ::          ConstructPartition
    private ::          CoupleSedimentQuality
    private ::          StartOutputBoxFluxes
    private ::          Construct_Time_Serie
    private ::          ConstructGlobalOutput
    private ::          Open_HDF5_OutPut_File
    private ::     ConstructLog
    private ::          CoupleSeagrassesRoots   !Isabella
    private ::             RootsOccupation      !Isabella
    private ::             DistributeRoots      !Isabella
    private ::     Read_Old_Properties_2D       !Isabella

    !Selector
    public  :: SedimentPropertyExists
    private ::      Search_Property
    public  :: GetSedimentConcentration
    public  :: GetSedimentDryDensity
    public  :: GetSedimentPropertyOptions
    public  :: UngetSedimentProperties
    public  :: SetFluxToSedimentProperties
    public  :: SetSedimentWaterFlux
    public  :: GetSeagrassesRootsRates
    public  :: GetRootsArray2D

    !Modifier
    public  :: SedimentProperties_Evolution
    private ::      Modify_SurfaceBoundaryFluxes
    private ::      Bioturbation_Processes
    private ::          Compute_Bioturbation
    private ::      Advection_Diffusion_Processes
    private ::          ComputeDiffusivity
    private ::      Partition_Processes
    private ::      SedimentQuality_Processes
    private ::          ComputeDissolvedToParticulate3D
    private ::      OutPut_Results_HDF
    private ::      OutPut_TimeSeries
    private ::      OutPut_BoxTimeSeries
    private ::      TimeStepActualization
    private ::      Actualize_Time_Evolution
    private ::      SeagrassesRoots_Processes

    !Destructor
    public  :: KillSedimentProperties
    private ::      DeallocateInstance
    private ::      Write_Final_SedProperties_HDF 
   
    !Management
    private ::      Ready
    private ::          LocateObjSedimentProperties

    private ::              ReadLockExternalVar
    private ::              ReadUnlockExternalVar

    !Interfaces----------------------------------------------------------------

    private :: UngetSedimentProperties3D
    interface  UngetSedimentProperties
        module procedure UngetSedimentProperties3D
        module procedure UngetSedimentProperties2Dreal8
        module procedure UngetSedimentProperties2Dreal4
    end interface UngetSedimentProperties

    !Parameter-----------------------------------------------------------------
    integer, parameter                          :: DirectionX           = 1
    integer, parameter                          :: DirectionY           = 2
    
    character(LEN = StringLength), parameter    :: prop_block_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter    :: prop_block_end       = '<endproperty>'
    character(LEN = StringLength), parameter    :: rate_block_begin     = '<beginrate>'
    character(LEN = StringLength), parameter    :: rate_block_end       = '<endrate>'

    !Correction to compute mol. dif. coefs.
    integer, parameter                          :: Berner_1980          = 1
    integer, parameter                          :: Soetaert_1996        = 2

                                                                        
    !Types---------------------------------------------------------------------
    type       T_ID
        integer                                 :: IDNumber     = null_int !inicialization: Carina
        character(LEN = StringLength)           :: Name         = null_str  !inicialization: Carina
        character(LEN = StringLength)           :: Description  = null_str  !inicialization: Carina
        character(LEN = StringLength)           :: Units        = null_str  !inicialization: Carina
    end type T_ID

    type       T_Property_3D
         type(T_PropertyID)                     :: ID
         real                                   :: Scalar  = null_real !inicialization: Carina
         real, pointer, dimension (:,:,:)       :: Field   => null() !inicialization: Carina
    end type   T_Property_3D

    type       T_AdvectionDiffusion_Parameters
        integer                                 :: BoundaryCondition   = null_int !inicialization: Carina
        integer                                 :: Diffusion_Method    = null_int !inicialization: Carina
        real                                    :: Molecular_Diff_Coef = null_real !inicialization: Carina
    end type T_AdvectionDiffusion_Parameters

    type       T_Partition                      
        type(T_Property_3D)                     :: Rate
        type(T_Property_3D)                     :: Fraction 
        character(LEN = StringLength)           :: Couple     = null_str  !inicialization: Carina
        integer                                 :: Couple_ID  = null_int
    end type T_Partition

    type       T_Bioturbation 
        real                                    :: DefaultCoef = null_real !inicialization: Carina
        real                                    :: BioDepth    = null_real !inicialization: Carina
        real                                    :: DecayCoef   = null_real !inicialization: Carina
        real, pointer, dimension(:,:,:)         :: Coef        => null() !inicialization: Carina
    end type   T_Bioturbation
    
    type       T_Evolution
        logical                                 :: Variable             = .false.
        logical                                 :: SedimentQuality      = .false.
        logical                                 :: Partitioning         = .false.
        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: ComputeBioturbation  = .false.
        logical                                 :: SurfaceFluxes        = .false.
        logical                                 :: MinConcentration     = .false.
        logical                                 :: SeagrassesRoots      = .false.
        type(T_AdvectionDiffusion_Parameters)   :: Advec_Difus_Parameters
        type(T_Partition)                       :: Partition
        type(T_Bioturbation)                    :: Bioturbation
        type(T_Time)                            :: LastCompute
        type(T_Time)                            :: NextCompute
        real                                    :: DTInterval           = FillValueReal
    end type T_Evolution

    type       T_Property
        type(T_PropertyID)                      :: ID
        type(T_Evolution)                       :: Evolution
        real                                    :: Scalar               = FillValueReal
        real                                    :: IScoefficient        = null_real
        logical                                 :: Particulate          = .false.
        logical                                 :: Old                  = .false.
        logical                                 :: TimeSerie            = .false.
        logical                                 :: BoxTimeSerie         = .false.
        logical                                 :: OutputHDF            = .false.
        real, pointer, dimension(:,:,:)         :: Concentration        => null() !inicialization: Carina
        real, pointer, dimension(:,:,:)         :: HorizontalDiffusivity => null() !inicialization: Carina
        real, pointer, dimension(:,:,:)         :: VerticalDiffusivity  => null() !inicialization: Carina
        real, pointer, dimension(:,:,:)         :: Mass_Created         => null() !inicialization: Carina
        real, pointer, dimension(:,:  )         :: BoundaryFlux         => null() !inicialization: Carina
        real                                    :: MinValue             = FillValueReal
        type(T_Property), pointer               :: Next                 => null() !inicialization: Carina
        type(T_Property), pointer               :: Prev                 => null() !inicialization: Carina
    end type T_Property

    type       T_SedimentRate
        type (T_ID)                             :: ID
        type (T_ID)                             :: FirstProp
        type (T_ID)                             :: SecondProp
        real, pointer, dimension(:,:,:)         :: Field         => null() !inicialization: Carina
        real, pointer, dimension(:,:,:)         :: Field2        => null() !inicialization: Carina
        type(T_SedimentRate), pointer           :: Next          => null() !inicialization: Carina
        type(T_SedimentRate), pointer           :: Prev          => null() !inicialization: Carina
    end type T_SedimentRate

    type  T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime        => null() !inicialization: Carina   
         integer                                :: NextOutPut     = null_int !inicialization: Carina
         logical                                :: Yes            = .false. !inicialization: Carina
    end type   T_OutPut


    type       T_Files
         character(len=StringLength)            :: Initial     = null_str !inicialization: Carina
         character(len=StringLength)            :: Final       = null_str !inicialization: Carina
         character(len=StringLength)            :: HDFResults  = null_str !inicialization: Carina
         character(len=StringLength)            :: InputData   = null_str !inicialization: Carina
         character(len=StringLength)            :: BoxesFile   = null_str !inicialization: Carina
    end type T_Files

    type       T_Coupling
         type(T_time)                           :: NextCompute
         real                                   :: DT_Compute           = FillValueReal
         logical                                :: Yes                  = .false.
         integer                                :: NumberOfProperties   = 0
    end type T_Coupling  

    type       T_Coupled
         type(T_Coupling)                       :: AdvectionDiffusion
         type(T_Coupling)                       :: Bioturbation
         type(T_Coupling)                       :: SedimentQuality      
         type(T_Coupling)                       :: Partition
         type(T_Coupling)                       :: SurfaceFluxes
         type(T_Coupling)                       :: MinimumConcentration
         type(T_Coupling)                       :: TimeSerie
         type(T_Coupling)                       :: BoxTimeSerie
         type(T_Coupling)                       :: OutputHDF
         type(T_Coupling)                       :: SeagrassesRoots
    end type T_Coupled

    type       T_External
        type(T_Time)                            :: Now
        real,    pointer, dimension(:,:  )      :: GridCellArea      => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:  )      :: WaterFlux         => null() !inicialization: Carina   
        real,    pointer, dimension(:,:,:)      :: DWZ               => null() !inicialization: Carina   
        real,    pointer, dimension(:,:,:)      :: DZZ               => null() !inicialization: Carina   
        real,    pointer, dimension(:,:,:)      :: SZZ               => null() !inicialization: Carina   
        real,    pointer, dimension(:,:,:)      :: WaterPercentage   => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: VolumeZ           => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: WaterVolume       => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: WaterVolumeOld    => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: DrySedimentVolume => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: DrySedimentVolumeOld  => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: WaterFluxX        => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: WaterFluxY        => null() !inicialization: Carina   
        real(8), pointer, dimension(:,:,:)      :: WaterFluxZ        => null() !inicialization: Carina   
        real   , pointer, dimension(:,:,:)      :: Velocity_U        => null() !inicialization: Carina   
        real   , pointer, dimension(:,:,:)      :: Velocity_V        => null() !inicialization: Carina   
        real   , pointer, dimension(:,:,:)      :: Velocity_W        => null() !inicialization: Carina   
        real   , pointer, dimension(:,:,:)      :: Tortuosity        => null() !inicialization: Carina   
        integer, pointer, dimension(:,:,:)      :: ComputeFacesU3D   => null() !inicialization: Carina   
        integer, pointer, dimension(:,:,:)      :: ComputeFacesV3D   => null() !inicialization: Carina   
        integer, pointer, dimension(:,:,:)      :: ComputeFacesW3D   => null() !inicialization: Carina   
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D      => null() !inicialization: Carina   
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D     => null() !inicialization: Carina   
        integer, pointer, dimension(:,:  )      :: WaterPoints2D     => null() !inicialization: Carina   
        integer, pointer, dimension(:,:,:)      :: LandPoints3D      => null() !inicialization: Carina   
        integer, pointer, dimension(:,:  )      :: BoundaryPoints2D  => null() !inicialization: Carina   
        real,    pointer, dimension(:,:  )      :: Bathymetry        => null() !inicialization: Carina   
        real,    pointer, dimension(:,:  )      :: XX_IE             => null() !inicialization: Carina   
        real,    pointer, dimension(:,:  )      :: YY_IE             => null() !inicialization: Carina
        real   , pointer, dimension(:,:,:)      :: Depth             => null() !inicialization: Carina   
        real   , pointer, dimension(:,:,:)      :: Porosity          => null() !inicialization: Carina   
        integer, pointer, dimension(:,:  )      :: KTop              => null() !inicialization: Carina   
        integer, pointer, dimension(:,:  )      :: KTopState         => null() !inicialization: Carina
        real                                    :: MinLayerThickness = null_real !inicialization: Carina   
        real                                    :: DT                = null_real !inicialization: Carina   
    end type T_External

   type       T_SeagrassesRoots
        type(T_PropertyID)                      :: ID
        real,    pointer, dimension(:,:  )      :: Biomass       => null() !inicialization: Carina    !gdw/m2
        real,    pointer, dimension(:,:  )      :: Length        => null() !inicialization: Carina
        real,    pointer, dimension(:,:,:)      :: Occupation    => null() !inicialization: Carina
        real,    pointer, dimension(:,:,:)      :: NintFactor3DR => null() !inicialization: Carina
        real,    pointer, dimension(:,:  )      :: NintFactor2DR => null() !inicialization: Carina
        real,    pointer, dimension(:,:,:)      :: PintFactor3DR => null() !inicialization: Carina
        real,    pointer, dimension(:,:  )      :: PintFactor2DR => null() !inicialization: Carina
        real,    pointer, dimension(:,:,:)      :: RootsMort3DR  => null() !inicialization: Carina
        real,    pointer, dimension(:,:  )      :: RootsMort2DR  => null() !inicialization: Carina
        real,    pointer, dimension(:,:,:)      :: UptakeNH4s3D  => null() !inicialization: Carina 
        real,    pointer, dimension(:,:,:)      :: UptakePO4s3D  => null() !inicialization: Carina
        real                                    :: DefaultValue, LBRatio = null_real !inicialization: Carina
        real(8), pointer, dimension(:,:,:)      :: Volume        => null() !inicialization: Carina
    end type   T_SeagrassesRoots


    type      T_SedimentProperties 
        integer                                 :: InstanceID
        type(T_Size3D      )                    :: Size
        type(T_Size3D      )                    :: WorkSize
        type(T_Files       )                    :: Files
        type(T_Coupled     )                    :: Coupled
        type(T_Time        )                    :: BeginTime
        type(T_Time        )                    :: EndTime
        type(T_Time        )                    :: NextCompute 
        type(T_External    )                    :: ExternalVar
        type(T_OutPut      )                    :: OutPut
        type(T_Property    ), pointer           :: FirstProperty
        type(T_Property    ), pointer           :: LastProperty
        type(T_SedimentRate), pointer           :: FirstSedimentRate
        type(T_SedimentRate), pointer           :: LastSedimentRate
        type(T_Property_3D)                     :: Sediment_DryDensity
        type(T_Property_3D)                     :: HorizontalDiffCoef
        type(T_Property_3D)                     :: VerticalDiffCoef
        type(T_SeagrassesRoots)                 :: SeagrassesRoots
        integer                                 :: PropertiesNumber     = 0
        integer                                 :: SedimentRatesNumber  = 0
        real(8), pointer, dimension(:,:  )      :: PartitionMatrix      => null() !inicialization: Carina
        real   , pointer, dimension(:    )      :: IndependentTerm      => null() !inicialization: Carina
        real   , pointer, dimension(:,:,:)      :: DissolvedToParticulate3D => null() !inicialization: Carina
        real                                    :: ResidualTime         = 0.
        real(8), pointer, dimension(:,:,:)      :: MassFluxesX          => null() !inicialization: Carina
        real(8), pointer, dimension(:,:,:)      :: MassFluxesY          => null() !inicialization: Carina
        real(8), pointer, dimension(:,:,:)      :: MassFluxesZ          => null() !inicialization: Carina
        real(8), pointer, dimension(:,:,:)      :: CellMass             => null() !inicialization: Carina
        real(8), pointer, dimension(:    )      :: Ti_Coef              => null() !inicialization: Carina
		real(8), pointer, dimension(:    )      :: D_Coef               => null() !inicialization: Carina
		real(8), pointer, dimension(:    )      :: E_Coef               => null() !inicialization: Carina
		real(8), pointer, dimension(:    )      :: F_Coef               => null() !inicialization: Carina
		real(8), pointer, dimension(:    )      :: VolZOld, VolZ        => null() !inicialization: Carina
        real,    pointer, dimension(:    )      :: DZZ, DiffCoef, Concentration, DWZ  => null() !inicialization: Carina

        !Instance of ModuleHDF5        
        integer                                 :: ObjHDF5              = 0

        !Instance of ModuleTimeSerie
        integer                                 :: ObjTimeSerie         = 0

        !Instance of Module_EnterData
        integer                                 :: ObjEnterData         = 0
       
        !Instance of ModuleGridData
        integer                                 :: ObjGridData          = 0
      
        !Instance of ModuleHorizontalGrid
        integer                                 :: ObjHorizontalGrid    = 0

        !Instance of ModuleGeometry
        integer                                 :: ObjGeometry          = 0

        !Instance of ModuleHorizontalMap
        integer                                 :: ObjHorizontalMap     = 0

        !Instance of ModuleMap
        integer                                 :: ObjMap               = 0
      
        !Instance of ModuleTime
        integer                                 :: ObjTime              = 0

        !Instance of ModuleBoxDif
        integer                                 :: ObjBoxDif            = 0

        !Instance of ModuleSedimentHydrodynamic
        integer                                 :: ObjConsolidation     = 0
        
        !Instance of ModuleInterface
        integer                                 :: ObjInterface         = 0
        
        !Instance of ModuleInterface
        integer                                 :: ObjLUD               = 0  
        
        
        !Instance of ModuleSeagrassesRoots
        integer                                 :: ObjSeagrassSedimInteraction     = 0  

        !Collection of instances
        type(T_SedimentProperties), pointer     :: Next => null() !inicialization: Carina

    end type T_SedimentProperties
    
    !Global Module Variables
    type (T_SedimentProperties),    pointer     :: FirstObjSedimentProperties   => null() !inicialization: Carina
    type (T_SedimentProperties),    pointer     :: Me      => null() !inicialization: Carina

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine Construct_SedimentProperties(SedimentPropertiesID,           &
                                            TimeID,                         &
                                            HorizontalGridID,               &
                                            GridDataID,                     &
                                            HorizontalMapID,                &
                                            MapID,                          &
                                            GeometryID,                     &
                                            ConsolidationID,                &
                                            STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: SedimentPropertiesID
        integer                                         :: TimeID
        integer                                         :: HorizontalGridID
        integer                                         :: GridDataID
        integer                                         :: HorizontalMapID
        integer                                         :: MapID
        integer                                         :: GeometryID
        integer                                         :: ConsolidationID
        integer, optional, intent(OUT)                  :: STAT   
          
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: ready_         
                       
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_

        !Begin-----------------------------------------------------------------


        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSedimentProperties_)) then
            nullify (FirstObjSedimentProperties)
            call RegisterModule (mSedimentProperties_) 
        endif

        call Ready(SedimentPropertiesID, ready_)

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            !Allocates a new Instance
            call AllocateInstance

            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )
            Me%ObjConsolidation  = AssociateInstance (mCONSOLIDATION_,  ConsolidationID )

            call ReadLockExternalVar

            !Initialize the water properties number   
            Me%PropertiesNumber = 0

            !Initialize the water properties list   
            nullify (Me%FirstProperty)
            nullify (Me%LastProperty )

            call GetComputeTimeLimits(Me%ObjTime,               &
                                      EndTime   = Me%EndTime,   &
                                      BeginTime = Me%BeginTime, &
                                      STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'Construct_SedimentProperties - ModuleSedimentProperties. ERR01'

            !Actualize the time
            Me%ExternalVar%Now = Me%BeginTime
            Me%NextCompute     = Me%BeginTime

            call ReadSedimentPropertiesFilesName

            call ConstructEnterData(Me%ObjEnterData, Me%Files%InPutData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)    &
                stop 'Construct_SedimentProperties - ModuleSedimentProperties - ERR01'

            
            call Construct_GlobalVariables

            call Construct_PropertyList

            call Construct_SedimentRateList

            call Construct_Sub_Modules

            call ConstructLog

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)    &
                stop 'Construct_SedimentProperties - ModuleSedimentProperties - ERR02'

            call ReadUnlockExternalVar

            call SetSubModulesConstructor

            SedimentPropertiesID = Me%InstanceID

            STAT_ = SUCCESS_
        
        else 
            
            stop 'Construct_SedimentProperties - ModuleSedimentProperties - ERR99' 

        end if cd0

        if (present(STAT)) STAT = STAT_

    end subroutine Construct_SedimentProperties

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_SedimentProperties), pointer           :: NewSedimentProperties
        type (T_SedimentProperties), pointer           :: PreviousSedimentProperties

        !Allocates new instance
        allocate (NewSedimentProperties)
        nullify  (NewSedimentProperties%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSedimentProperties)) then
            FirstObjSedimentProperties      => NewSedimentProperties
            Me                              => NewSedimentProperties
        else
            PreviousSedimentProperties      => FirstObjSedimentProperties
            Me                              => FirstObjSedimentProperties%Next
            do while (associated(Me))
                PreviousSedimentProperties  => Me
                Me                          => Me%Next
            enddo
            Me                              => NewSedimentProperties
            PreviousSedimentProperties%Next => NewSedimentProperties
        endif

        Me%InstanceID = RegisterNewInstance (mSEDIMENTPROPERTIES_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        integer                                             :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        type(T_Property    ),                       pointer :: PropertyX
        type(T_SedimentRate),                       pointer :: SedimentRateX
        character(len=StringLength), dimension(:),  pointer :: ScalarOutputList
        character(len=StringLength), dimension(:),  pointer :: FluxesOutputList
        integer                                             :: nScalars, nFluxes, n

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        ! This keyword have two functions if exist fluxes between boxes are compute 
        ! and the value read is the name file where the boxes are defined
        call GetData(Me%Files%BoxesFile,                                            &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule = 'ModuleSedimentProperties',                     &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR01'
        if (iflag .EQ. 0)                                                           &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR02'    
        
        inquire(File = Me%Files%BoxesFile, Exist = exist)
        if (exist) then
            inquire(File = Me%Files%BoxesFile, Opened  = Opened)
            if (opened) then
                write(*,*    ) 
                write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Files%BoxesFile))
                write(*,*    ) 'Already opened.'
                stop           'StartOutputBoxFluxes - ModuleSedimentProperties - ERR03'    
            end if
        else
            write(*,*) 
            write(*,*)     'Could not find the boxes file.'
            write(*,'(A)') 'BoxFileName = ', Me%Files%BoxesFile
            stop           'StartOutputBoxFluxes - ModuleSedimentProperties - ERR04'    
        end if

        nScalars = Me%Coupled%BoxTimeSerie%NumberOfProperties + Me%SedimentRatesNumber
        nFluxes  = Me%Coupled%BoxTimeSerie%NumberOfProperties

        allocate(ScalarOutputList(nScalars), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR05'

        allocate(FluxesOutputList(nFluxes), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR06'

        PropertyX   => Me%FirstProperty
        n = 0
        do while (associated(PropertyX))

            if (PropertyX%BoxTimeSerie) then
                n = n + 1
                ScalarOutputList(n) ="Sed_"//trim(PropertyX%ID%name)
                FluxesOutputList(n) ="Sed_"//trim(PropertyX%ID%name)
            end if 

            PropertyX=>PropertyX%Next
        end do

        SedimentRateX => Me%FirstSedimentRate
        do while(associated(SedimentRateX))
            n = n + 1
            ScalarOutputList(n) ="Sed_"//trim(SedimentRateX%ID%name)
            SedimentRateX => SedimentRateX%Next
        end do


        call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &
                         TimeID             = Me%ObjTime,                   &
                         HorizontalGridID   = Me%ObjHorizontalGrid,         &
                         BoxesFilePath      = Me%Files%BoxesFile,           &
                         FluxesOutputList   = FluxesOutputList,             &
                         ScalarOutputList   = ScalarOutputList,             &
                         WaterPoints3D      = Me%ExternalVar%WaterPoints3D, &
                         Size3D             = Me%Size,                      &
                         WorkSize3D         = Me%WorkSize,                  &
                         STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR07'

        deallocate(FluxesOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR08'

        deallocate(ScalarOutputList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR09'

        allocate(Me%MassFluxesX(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ER10'
        Me%MassFluxesX(:,:,:) = 0.

        allocate(Me%MassFluxesY(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR11'
        Me%MassFluxesY(:,:,:) = 0.

        allocate(Me%MassFluxesZ(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR12'
        Me%MassFluxesZ(:,:,:) = 0.

        allocate(Me%CellMass(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'StartOutputBoxFluxes - ModuleSedimentProperties - ERR13'
        Me%CellMass(:,:,:) = 0.

    end subroutine StartOutputBoxFluxes

    
    !--------------------------------------------------------------------------

    
    subroutine Construct_GlobalVariables

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_Property_3D), pointer            :: Scalar3D
        
        !----------------------------------------------------------------------

        call GetGeometrySize(Me%ObjGeometry,                                    &
                             Size     = Me%Size,                                &
                             WorkSize = Me%WorkSize,                            &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'Construct_GlobalVariables - ModuleSedimentProperties - ERR01'

        Scalar3D => Me%Sediment_DryDensity
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                               block_begin = '<begin_sed_dry_density>',         &
                               block_end   = '<end_sed_dry_density>')

        Scalar3D => Me%VerticalDiffCoef
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                               block_begin = '<begin_vert_diff_coef>',          &
                               block_end   = '<end_vert_diff_coef>')
       
        Scalar3D => Me%HorizontalDiffCoef
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                              block_begin = '<begin_horiz_diff_coef>',          &
                              block_end   = '<end_horiz_diff_coef>')

        call GetConsolidationMinThickness(Me%ObjConsolidation, Me%ExternalVar%MinLayerThickness, &
                                          STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'Construct_GlobalVariables - ModuleSedimentProperties - ERR02'



    end subroutine Construct_GlobalVariables

    
    !--------------------------------------------------------------------------


    subroutine ConstructScalar3D(Scalar3D, ExtractType, ClientNumber, block_begin, block_end)

        !Arguments-------------------------------------------------------------
        type(T_Property_3D), pointer        :: Scalar3D
        integer, intent(in)                 :: ExtractType
        integer, intent(in), optional       :: ClientNumber
        character(len=*)                    :: block_begin, block_end

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: BlockClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        
        !----------------------------------------------------------------------
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB


        select case(ExtractType)

            case(FromBlock)

                call ExtractBlockFromBuffer(Me%ObjEnterData, BlockClientNumber, block_begin, block_end,  &
                                            BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR01'

            case(FromBlockInBlock)

                if(.not. present(ClientNumber))then
                    stop 'ConstructScalar3D - ModuleSedimentProperties - ERR02'
                end if
                
                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber, block_begin, block_end,  &
                                           BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR03'

        end select

        if(BlockFound)then

            allocate(Scalar3D%Field(ILB:IUB, JLB:JUB, KLB:KUB))

            call ConstructFillMatrix  (PropertyID           = Scalar3D%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill3D       = Me%ExternalVar%WaterPoints3D,     &
                                       Matrix3D             = Scalar3D%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR04'


            call GetDefaultValue(Scalar3D%ID%ObjFillMatrix, Scalar3D%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR05'

            call KillFillMatrix(Scalar3D%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR06'
            
            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, BlockClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR07'
            end if
        
        else

            stop 'ConstructScalar3D - ModuleSedimentProperties - ERR08'

        end if
        
        select case(ExtractType)

            case(FromBlock)

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR09'

            case(FromBlockInBlock)

                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleSedimentProperties - ERR10'

        end select

   
    end subroutine ConstructScalar3D

    !--------------------------------------------------------------------------
    
    subroutine Construct_PropertyList

        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_Property), pointer          :: NewProperty

        !------------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
                                        ClientNumber    = ClientNumber,     &
                                        block_begin     = prop_block_begin, &
                                        block_end       = prop_block_end,   &
                                        BlockFound      = BlockFound,       &
                                        STAT            = STAT_CALL)
cd1 :       if (STAT_CALL .EQ. SUCCESS_) then    
cd2 :           if (BlockFound) then                                                  
                    
                    Call Construct_Property(NewProperty, ClientNumber)

                    Call Add_Property(NewProperty)

                else cd2

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_PropertyList - ModuleSedimentProperties - ERR01'
                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_PropertyList - ModuleSedimentProperties - ERR02'
            else cd1
                stop       'Construct_PropertyList - ModuleSedimentProperties - ERR03'
            end if cd1
        end do do1

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------

    subroutine Construct_SedimentRateList

        !External----------------------------------------------------------------
        integer                                :: ClientNumber
        integer                                :: STAT_CALL
        logical                                :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_SedimentRate),    pointer      :: NewSedimentRate

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
                    
                    !Construct a New Sediment Rate
                    Call Construct_SedimentRate(NewSedimentRate)

                    !Add new Rate to the Sediment Rates List 
                    Call Add_SedimentRate(NewSedimentRate)

                else
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'Construct_SedimentRateList - ModuleSedimentProperties - ERR01'

                    exit do1    !No more blocks
                end if


            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_SedimentRateList - ModuleSedimentProperties - ERR02'
            else
                stop       'Construct_SedimentRateList - ModuleSedimentProperties - ERR03'
            end if
        end do do1
         
    end subroutine Construct_SedimentRateList

    !----------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5,                                              &
                           trim(Me%Files%HDFResults)//"5",                          &
                           HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR10'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR20'

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",               &
                              Array2D = Me%ExternalVar%Bathymetry,                  &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR30'


        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",            &
                              Array3D = Me%ExternalVar%WaterPoints3D,               &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR40'

        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB+1, WorkJLB,              &
                              WorkJUB+1, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR50'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX", "m",              &
                              Array2D = Me%ExternalVar%XX_IE,                       &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR60'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY", "m",              &
                              Array2D = Me%ExternalVar%YY_IE,                       &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR70'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'Open_HDF5_OutPut_File - ModuleSedimentProperties - ERR80'

    end subroutine Open_HDF5_OutPut_File

   
    !----------------------------------------------------------------------
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

        ILB = Me%WorkSize%ILB 
        JLB = Me%WorkSize%JLB
        IUB = Me%WorkSize%IUB 
        JUB = Me%WorkSize%JUB

        ObjHDF5 = 0


        inquire(File = trim(Me%Files%Initial)//"5", Exist = exist)
        
        if(.not. exist)then
            write(*,*) 
            write(*,*)     'Could not find the final InterfaceSedimentWater file.'
            write(*,'(A)') 'BoxFileName = ', trim(Me%Files%Initial)//"5"
            stop           'Read_Old_Properties_2D - ModuleSedimentProperty - ERR00'    
        end if

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5, trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'read_Old_Properties_2D - ModuleSedimentProperty - ERR01'
            
        !Reads from HDF file the Property concentration and open boundary values
        call HDF5SetLimits  (ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'read_Old_Properties_2D - ModuleSedimentProperty - ERR02'
            
        call HDF5ReadData(ObjHDF5, "/Results/"//PropertyName, &
                          PropertyName, Array2D = Scalar_2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'read_Old_Properties_2D - ModuleSedimentProperty - ERR03'            
        
        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'read_Old_Properties_2D - ModuleSedimentProperty - ERR04'            

    end subroutine Read_Old_Properties_2D
    
    
    !----------------------------------------------------------------------
    
    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: nProperties, STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !----------------------------------------------------------------------

        !Allocates PropertyList
        allocate(PropertyList(1:Me%Coupled%TimeSerie%NumberOfProperties), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Time_Serie - ModuleSedimentProperties - ERR01'

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
        
        call GetData(TimeSerieLocationFile,                                 &
                     Me%ObjEnterData,iflag,                                 &
                     SearchType   = FromFile,                               &
                     keyword      = 'TIME_SERIE_LOCATION',                  &
                     ClientModule = 'ModuleSedimentProperties',             &
                     Default      = Me%Files%InputData,                     &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                        &
            stop 'Construct_Time_Serie - ModuleSedimentProperties - ERR02' 


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie,                                &
                            Me%ObjTime,                                     &
                            TimeSerieLocationFile,                          &
                            PropertyList, "sre",                            &
                            WaterPoints3D = Me%ExternalVar%WaterPoints3D,   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Time_Serie - ModuleSedimentProperties - ERR03'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_Time_Serie - ModuleSedimentProperties - ERR04'

        
    end subroutine Construct_Time_Serie

    !------------------------------------------------------------------------

    subroutine Construct_Sub_Modules

        !Local-----------------------------------------------------------------
        type (T_Property),           pointer                 :: PropertyX

        !----------------------------------------------------------------------
        
        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%Evolution%AdvectionDiffusion) then
                Me%Coupled%AdvectionDiffusion%NumberOfProperties    = &
                Me%Coupled%AdvectionDiffusion%NumberOfProperties    + 1
                Me%Coupled%AdvectionDiffusion%Yes                   = ON
            endif

            if (PropertyX%Evolution%SedimentQuality) then
                Me%Coupled%SedimentQuality%NumberOfProperties       = &
                Me%Coupled%SedimentQuality%NumberOfProperties       + 1
                Me%Coupled%SedimentQuality%Yes                      = ON
            endif

            if (PropertyX%Evolution%Partitioning) then
                Me%Coupled%Partition%NumberOfProperties             = &
                Me%Coupled%Partition%NumberOfProperties             + 1
                Me%Coupled%Partition%Yes                            = ON
            endif

            if (PropertyX%Evolution%ComputeBioturbation) then
                Me%Coupled%Bioturbation%NumberOfProperties          = &
                Me%Coupled%Bioturbation%NumberOfProperties          + 1
                Me%Coupled%Bioturbation%Yes                         = ON
            endif


            if (PropertyX%Evolution%SurfaceFluxes) then
                Me%Coupled%SurfaceFluxes%NumberOfProperties         = &
                Me%Coupled%SurfaceFluxes%NumberOfProperties         + 1
                Me%Coupled%SurfaceFluxes%Yes                        = ON
            endif

            if (PropertyX%Evolution%MinConcentration) then
                Me%Coupled%MinimumConcentration%NumberOfProperties  = &
                Me%Coupled%MinimumConcentration%NumberOfProperties  + 1
                Me%Coupled%MinimumConcentration%Yes                 = ON
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
            
            if (PropertyX%Evolution%SeagrassesRoots) then                       ! Isabella
                Me%Coupled%SeagrassesRoots%NumberOfProperties       = &
                Me%Coupled%SeagrassesRoots%NumberOfProperties       + 1
                Me%Coupled%SeagrassesRoots%Yes                      = ON
            endif


            PropertyX=>PropertyX%Next

        end do

        if (Me%Coupled%AdvectionDiffusion%Yes .or. Me%Coupled%Bioturbation%Yes) then
            call CoupleAdvectionDiffusion
        endif

        if (Me%Coupled%SedimentQuality%Yes)then
            call CoupleSedimentQuality
        end if

        if (Me%Coupled%Partition%Yes)then
            call ConstructPartition
        end if

        if (Me%Coupled%TimeSerie%Yes)then
            call Construct_Time_Serie
        end if

        if (Me%Coupled%BoxTimeSerie%Yes)then
            call StartOutputBoxFluxes
        end if

        if (Me%Coupled%OutputHDF%Yes)then
            call ConstructGlobalOutput
        end if
        
        if (Me%Coupled%SeagrassesRoots%Yes)then    ! isabella
            call CoupleSeagrassesRoots
        end if
        
        !Opens HDF format output file
        if (Me%OutPut%Yes) call Open_HDF5_Output_File

    end subroutine Construct_Sub_Modules
    
    
    !--------------------------------------------------------------------------

    subroutine ConstructLog


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrentProperty


        write(*, *)"-------------------SEDIMENTPROPERTIES --------------------"
        write(*, *)
        write(*, *)"Num of Properties : ", Me%PropertiesNumber
        write(*, *)

        CurrentProperty => Me%FirstProperty
        do while (associated(CurrentProperty))

            write(*, *)"Property          : ", trim(CurrentProperty%ID%Name)
            write(*, *)"---Sed. Quality   : ", CurrentProperty%Evolution%SedimentQuality
            write(*, *)"---Partitioning   : ", CurrentProperty%Evolution%Partitioning
            write(*, *)"---Adv. Diff.     : ", CurrentProperty%Evolution%AdvectionDiffusion
            write(*, *)"---Bioturbation   : ", CurrentProperty%Evolution%ComputeBioturbation
            write(*, *)"---Surface Fluxes : ", CurrentProperty%Evolution%SurfaceFluxes
            write(*, *)"---Seagrassesroots: ", CurrentProperty%Evolution%SeagrassesRoots
            write(*, *)

            CurrentProperty=>CurrentProperty%Next
        enddo

    end subroutine ConstructLog

    
    !--------------------------------------------------------------------------


    subroutine CoupleAdvectionDiffusion

        !Local-----------------------------------------------------------------
        integer                                             :: KLB, KUB
        integer                                             :: STAT_CALL

        !----------------------------------------------------------------------
        
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        allocate(Me%DiffCoef(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR01'
        
        allocate(Me%Ti_Coef(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR02'
        
        allocate(Me%D_Coef(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR03'
        
        allocate(Me%E_Coef(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR04'
        
        allocate(Me%F_Coef(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR05'
    
        allocate(Me%Concentration(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR06'

        allocate(Me%VolZ(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR07'
        
        allocate(Me%VolZOld(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR08'
        
        allocate(Me%DZZ(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR09'

        allocate(Me%DWZ(KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'CoupleAdvectionDiffusion - ModuleSedimentProperties - ERR10'

        Me%Coupled%AdvectionDiffusion%NextCompute = Me%ExternalVar%Now

    
    end subroutine CoupleAdvectionDiffusion

    !--------------------------------------------------------------------------
    
    subroutine CoupleSedimentQuality

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer                                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer, pointer, dimension(:)                      :: SedimentQualityPropertyList
        integer                                             :: STAT_CALL
        real                                                :: SedimentQualityDT
        integer                                             :: Index = 0
        !----------------------------------------------------------------------
        
        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        Index = 0

        nullify (SedimentQualityPropertyList)
        allocate(SedimentQualityPropertyList(1:Me%Coupled%SedimentQuality%NumberOfProperties))

        PropertyX => Me%FirstProperty
            
        do while(associated(PropertyX))
            
            if(PropertyX%Evolution%SedimentQuality)then
                Index = Index + 1
                SedimentQualityPropertyList(Index)  = PropertyX%ID%IDNumber
            end if

            PropertyX => PropertyX%Next
          
        enddo
        
        nullify(PropertyX)


        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = SedimentQualityModel,          &
                                DT                  = SedimentQualityDT,             &
                                PropertiesList      = SedimentQualityPropertyList,   &
                                WaterPoints3D       = Me%ExternalVar%WaterPoints3D,  &
                                Size3D              = Me%WorkSize,                   &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleSedimentQuality - ModuleSedimentProperties - ERR01'
        
        Me%Coupled%SedimentQuality%DT_Compute  = SedimentQualityDT 
        Me%Coupled%SedimentQuality%NextCompute = Me%ExternalVar%Now
        
        deallocate(SedimentQualityPropertyList)
        nullify   (SedimentQualityPropertyList)
          
        nullify (Me%DissolvedToParticulate3D)
        allocate(Me%DissolvedToParticulate3D(ILB:IUB, JLB:JUB, KLB:KUB))
        Me%DissolvedToParticulate3D(:,:,:) = null_real

        Me%ResidualTime = 0.

    end subroutine CoupleSedimentQuality
    
    
    !--------------------------------------------------------------------------
    
     subroutine CoupleSeagrassesRoots
        
        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        
        integer, pointer, dimension(:)                      :: SeagrassesRootsPropertyList
        integer                                             :: STAT_CALL, iflag, ClientNumber
        real                                                :: SeagrassesRootsDT
        integer                                             :: Index = 0
        integer                                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                             :: i, j
        logical                                             :: BlockFound

        !----------------------------------------------------------------------
        
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        !gC/m2
        call GetData(Me%SeagrassesRoots%DefaultValue,                       &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ROOTS_MASS',                                &
                     Default        = 0.001,                                            &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleSedimentProperties',                          &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSeagrassesRoots - ModuleSedimentProperties- ERR02'
        if (iflag == 0)            stop 'CoupleSeagrassesRoots - ModuleSedimentProperties- ERR03'
        
        !m
        call GetData(Me%SeagrassesRoots%LBratio,                           &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'ROOTS_LBRATIO',                           &
                     Default        = 0.003,                                            &
                     SearchType     = FromFile,                                         &
                     ClientModule   = 'ModuleSedimentProperties',                          &
                     STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSeagrassesRoots - ModuleSedimentProperties- ERR04'
        if (iflag == 0)            stop 'CoupleSeagrassesRoots - ModuleSedimentProperties- ERR05'
        
              
        

        ! Roots
        allocate(Me%SeagrassesRoots%Biomass  (ILB:IUB, JLB:JUB)) !gdw/m2
        Me%SeagrassesRoots%Biomass   (:,:) = Me%SeagrassesRoots%DefaultValue
        
        
        
        
        
             call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = "<begin_roots_biomass>",        &
                                        block_end       = "<end_roots_biomass>",          &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
cd11 :       if (STAT_CALL .EQ. SUCCESS_     ) then    
cd12 :       if (BlockFound) then                                                  


                call ConstructFillMatrix  (PropertyID           = Me%SeagrassesRoots%ID,           &
                                           EnterDataID          = Me%ObjEnterData,                  &
                                           TimeID               = Me%ObjTime,                       &
                                           HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                           ExtractType          = FromBlock,                        &
                                           PointsToFill2D       = Me%ExternalVar%WaterPoints2D,     &
                                           Matrix2D             = Me%SeagrassesRoots%Biomass,      &
                                           TypeZUV              = TypeZ_,                           &
                                           STAT                 = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                          &
                    stop 'CoupleSeagrassesRoots - ModuleSedimentProperties - ERR07'


                call KillFillMatrix(Me%SeagrassesRoots%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'CoupleSeagrassesRoots - ModuleSedimentProperties - ERR08'

            endif cd12
            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'CoupleSeagrassesRoots - ModuleSedimentProperties - ERR08a'
            
            endif cd11
        
        
        
        
        allocate(Me%SeagrassesRoots%Length  (ILB:IUB, JLB:JUB)) 
        do j = JLB, JUB
        do i = ILB, IUB
        
           Me%SeagrassesRoots%Length   (i,j) = Me%SeagrassesRoots%Biomass(i,j) * &
                                                            Me%SeagrassesRoots%LBRatio

        enddo
        enddo
        allocate(Me%SeagrassesRoots%Occupation    (ILB:IUB, JLB:JUB, KLB:KUB)) 
        Me%SeagrassesRoots%Occupation(:,:,:) = 0.
        

        allocate (Me%SeagrassesRoots%NintFactor3DR(ILB:IUB,JLB:JUB,KLB:KUB))
        Me%SeagrassesRoots%NintFactor3DR(:,:,:) =0.
        
        allocate (Me%SeagrassesRoots%NintFactor2DR(ILB:IUB,JLB:JUB))
        Me%SeagrassesRoots%NintFactor2DR(:,:) =0.
        
        allocate (Me%SeagrassesRoots%PintFactor3DR(ILB:IUB,JLB:JUB,KLB:KUB))
        Me%SeagrassesRoots%PintFactor3DR(:,:,:) =0.
        
        allocate (Me%SeagrassesRoots%PintFactor2DR(ILB:IUB,JLB:JUB))
        Me%SeagrassesRoots%PintFactor2DR(:,:) =0.

        
        allocate (Me%SeagrassesRoots%UptakeNH4s3D(ILB:IUB,JLB:JUB,KLB:KUB))
        Me%SeagrassesRoots%UptakeNH4s3D(:,:,:) =0.
        
        allocate (Me%SeagrassesRoots%UptakePO4s3D(ILB:IUB,JLB:JUB,KLB:KUB))
        Me%SeagrassesRoots%UptakePO4s3D(:,:,:) =0.
        
        allocate(Me%SeagrassesRoots%Volume(ILB:IUB,JLB:JUB,KLB:KUB))
        Me%SeagrassesRoots%Volume(:,:,:) =0.
        
        allocate (Me%SeagrassesRoots%RootsMort2DR(ILB:IUB,JLB:JUB))
        Me%SeagrassesRoots%RootsMort2DR(:,:) =0.
        
        allocate (Me%SeagrassesRoots%RootsMort3DR(ILB:IUB,JLB:JUB,KLB:KUB))
        Me%SeagrassesRoots%RootsMort3DR(:,:,:) =0.
       
        ! Forcings


        Index = 0

        nullify (SeagrassesRootsPropertyList)
        allocate(SeagrassesRootsPropertyList(1:Me%Coupled%SeagrassesRoots%NumberOfProperties))

        PropertyX => Me%FirstProperty
            
 DoProd1:       do while(associated(PropertyX))

            if(PropertyX%Evolution%SeagrassesRoots)then

                Index                               = Index + 1
                SeagrassesRootsPropertyList(Index) = PropertyX%ID%IDNumber
                
                
                if(PropertyX%ID%IDNumber == SeagrassesRoots_)then
                
                
                                
                  if(PropertyX%Evolution%AdvectionDiffusion)then
                    
                    write(*,*)
                    write(*,*)'Seagrasses Roots can not have ADVECTION_DIFFUSION'
                    write(*,*)'Property : ', trim(adjustl(PropertyX%ID%Name))
                    write(*,*)'ADVECTION_DIFFUSION will be switched off'
                    
                    PropertyX%Evolution%AdvectionDiffusion = OFF

                end if

                    if(PropertyX%Old)then
                        call Read_Old_Properties_2D(Me%SeagrassesRoots%Biomass, "seagrasses roots biomass" )
                    else
                        call RootsOccupation(Me%SeagrassesRoots)
                        call DistributeRoots(PropertyX, Me%SeagrassesRoots)
                    end if

                end if

               
                
                


            end if
            
            PropertyX => PropertyX%Next
        enddo DoProd1
        
        nullify(PropertyX)

        call ConstructInterface(InterfaceID         = Me%ObjSeagrassSedimInteraction,   &
                                TimeID              = Me%ObjTime,                        &
                                SinksSourcesModel   = SeagrassSedimInteractionModel,                   &
                                DT                  = SeagrassesRootsDT,                &
                                PropertiesList      = SeagrassesRootsPropertyList,      &
                                WaterPoints3D       = Me%ExternalVar%WaterPoints3D,      &
                                Size3D              = Me%WorkSize,                       &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'CoupleSeagrassesRoots - ModuleSedimentProperties - ERR01'
        
        Me%Coupled%SeagrassesRoots%DT_Compute   = SeagrassesRootsDT 
        Me%Coupled%SeagrassesRoots%NextCompute  = Me%ExternalVar%Now
        
        deallocate(SeagrassesRootsPropertyList)
        nullify   (SeagrassesRootsPropertyList)


    end subroutine CoupleSeagrassesRoots
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    subroutine  RootsOccupation(SeagrassesRoots)  
    
  
          
        !Arguments------------------------------------------------------------- 
        type (T_SeagrassesRoots)                          :: SeagrassesRoots
        
        !Local----------------------------------------------------------------- 
        real, allocatable, dimension(:,:  )     :: SedimentColumnZ
        integer                                 :: i, j, k, KTop
        integer                                 :: ILB, IUB, JLB, JUB, KUB, KLB
        integer                                 :: STAT_CALL
        real                                    :: Remaining_root_Length
        
        !Begin----------------------------------------------------------------- 

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KUB = Me%WorkSize%KUB  
        KLB = Me%WorkSize%KLB

        !SedimentColumnZ
        
         allocate(SedimentColumnZ(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
         SedimentColumnZ(ILB:IUB, JLB:JUB) =0.
        

            do j = JLB, JUB
            do i = ILB, IUB
 
                     do k = KLB, Me%Externalvar%KTop(i,j)
            
                       SedimentcolumnZ(i,j)=SedimentcolumnZ(i,j)+Me%Externalvar%DWZ(i,j,k)
                     
                     enddo
         
            enddo
            enddo

        call SetMatrixValue(SeagrassesRoots%Occupation, Me%WorkSize, 0.)

        !if running in 3D
        if(Me%WorkSize%KUB > 1)then

            do j = JLB, JUB
            do i = ILB, IUB
               ! i have a doubt if it should be :
               !                 ktop= Me%ExternalVar%KTop(i, j)
               !                 if (Me%ExternalVar%OpenPoints3D(i, j, ktop) == OpenPoint) then
               
               
                if (Me%ExternalVar%OpenPoints3D(i, j, KUB) == OpenPoint) then

                    !KTop = Me%ExternalVar%KTop(i, j)
                    SeagrassesRoots%Length(i,j)=SeagrassesRoots%Biomass(i,j)*SeagrassesRoots%LBratio
                    if(SeagrassesRoots%Length(i,j) .ge. SedimentColumnZ(i,j))then

                        SeagrassesRoots%Occupation(i,j,KLB:Me%ExternalVar%KTop(i, j)) = 1

                    else

                        k = Me%ExternalVar%KTop(i, j)    ! index for the top sediment layer

                        Remaining_root_Length = SeagrassesRoots%Length(i,j)

                        do while(k .GE. KLB)

                            
                            if  (Remaining_root_Length >= Me%ExternalVar%DWZ(i,j,k))then

                                SeagrassesRoots%Occupation(i,j,k)   = 1.
                                Remaining_root_Length               = Remaining_root_Length - Me%ExternalVar%DWZ(i,j,k)
                                k                                   = k - 1

                            else
                                SeagrassesRoots%Occupation(i,j,k)   = Remaining_root_Length / Me%ExternalVar%DWZ(i,j,k)
                                k                                   = k - 1
                            end if

                        end do

                    end if

                endif

            enddo
            enddo
        
        else !if running in 2D (this way is faster)

            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%OpenPoints3D(i, j, KUB) == OpenPoint) then

                    KTop                 = Me%ExternalVar%KTop(i, j)
                    SeagrassesRoots%Length(i,j)=SeagrassesRoots%Biomass(i,j)*SeagrassesRoots%LBratio

                    if (SeagrassesRoots%Length(i,j) .ge. SedimentColumnZ(i,j))then
                        
                        SeagrassesRoots%Occupation(i,j,KLB:KTop) = 1.

                    else

                        SeagrassesRoots%Occupation(i,j,KLB:KTop) = SeagrassesRoots%Length(i,j) / &
                                                                      SedimentColumnZ(i,j)
                    end if

                endif

            enddo
            enddo

        end if

          deallocate(SedimentColumnZ)
         
    end  subroutine RootsOccupation  
         
         
!-------------------------------------------------------------------------------


     !--------------------------------------------------------------------------

    subroutine DistributeRoots(PropertyX, SeagrassesRoots)
        
        !Arguments------------------------------------------------------------- 
        type (T_Property),      pointer         :: PropertyX
        type (T_SeagrassesRoots)                :: SeagrassesRoots
       
        !Local----------------------------------------------------------------- 
        real, dimension(:,:  ), pointer         :: SedimentColumnZ
        integer                                 :: i, j, k, ktop  
        integer                                 :: ILB, IUB, JLB, JUB, KLB
        integer                                 :: STAT_CALL
        
        !Begin----------------------------------------------------------------- 

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB

          !SedimentColumnZ
        
         allocate(SedimentColumnZ(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
         SedimentColumnZ(ILB:IUB, JLB:JUB) =0.
        

            do j = JLB, JUB
            do i = ILB, IUB
 
                     do k = KLB, Me%Externalvar%KTop(i,j)
            
                       SedimentcolumnZ(i,j)=SedimentcolumnZ(i,j)+Me%Externalvar%DWZ(i,j,k)
                     
                     enddo
         
            enddo
            enddo

        if(Me%WorkSize%KUB > 1)then

            do j = JLB, JUB
            do i = ILB, IUB
               
               ktop = Me%ExternalVar%KTop(i, j)
                
                if (Me%ExternalVar%OpenPoints3D(i, j, ktop) == OpenPoint) then

                    
                    
                    SeagrassesRoots%Length(i,j)=SeagrassesRoots%Biomass(i,j)*SeagrassesRoots%LBratio

                    do k = KLB, ktop
                        !gC/m3 = gC/m2 * m2 / m3 * m / m
                        PropertyX%Concentration(i,j,k) = SeagrassesRoots%Occupation(i,j,k) * &
                                                         SeagrassesRoots%Biomass(i,j)      * &
                                                         Me%ExternalVar%GridCellArea(i,j)  / &
                                                         Me%ExternalVar%VolumeZ(i,j,k)     * &
                                                        (Me%ExternalVar%DWZ(i,j,k)         / &
                                                         SeagrassesRoots%Length(i,j))

                    enddo

                endif

            enddo
            enddo

        else

            do j = JLB, JUB
            do i = ILB, IUB


                  ktop = Me%ExternalVar%KTop(i, j)
                
                
                if (Me%ExternalVar%OpenPoints3D(i, j, ktop) == OpenPoint) then

        
                    do k = KLB, ktop

                        !gC/m3 = gC/m2 * m2 / m3 * m / m
                        PropertyX%Concentration(i,j,k) = SeagrassesRoots%Biomass(i,j) * &
                                                         Me%ExternalVar%GridCellArea(i,j)  / &
                                                         Me%ExternalVar%VolumeZ(i,j,k)

                    enddo



                endif

            enddo
            enddo

        endif

   deallocate(SedimentColumnZ)
   
   
    end subroutine DistributeRoots
    
    !---------------------------------------------------------------------------------------------


    
    subroutine ConstructPartition

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                   :: Property
        type(T_Property), pointer                   :: DissolvedProperty
        type(T_Property), pointer                   :: ParticulateProperty
        integer                                     :: STAT_CALL
        character(len=StringLength)                 :: PartPropName
        real                                        :: TotalPartition
        real                                        :: Error
        integer                                     :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB
        integer                                     :: i, j, k, Couple_ID

        !----------------------------------------------------------------------
        
        WILB = Me%WorkSize%ILB 
        WIUB = Me%WorkSize%IUB 
        WJLB = Me%WorkSize%JLB 
        WJUB = Me%WorkSize%JUB 
        WKLB = Me%WorkSize%KLB
        WKUB = Me%WorkSize%KUB

        Property => Me%FirstProperty
            
do1:    do while(associated(Property))

            if (Property%Evolution%Partitioning .and. .not. Property%Particulate) then

                DissolvedProperty => Property 

                PartPropName = trim(DissolvedProperty%Evolution%Partition%Couple)
              
                if (.not. CheckPropertyName(PartPropName, Couple_ID)) then

                    write(*,*)
                    write(*,*) 'The property name is not recognised by the model'
                    stop       'ConstructPartition - ModuleSedimentProperties - ERR01'
                     
                else

                    DissolvedProperty%Evolution%Partition%Couple_ID = Couple_ID

                end if

                call Search_Property(ParticulateProperty, PropertyXID = Couple_ID, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructPartition - ModuleSedimentProperties - ERR02'

                if(ParticulateProperty%Evolution%DTInterval .NE. &
                   DissolvedProperty%Evolution%DTInterval)       &
                    stop 'ConstructPartition - ModuleSedimentProperties - ERR03'
                
                do k = WKLB, WKUB
                do j = WJLB, WJUB
                do i = WILB, WIUB
                    if(DissolvedProperty%Evolution%Partition%Rate%Field(i,j,k) .ne.    &
                       ParticulateProperty%Evolution%Partition%Rate%Field(i,j,k))then

                        stop 'ConstructPartition - ModuleSedimentProperties - ERR04'

                    end if
                enddo
                enddo
                enddo

                TotalPartition = DissolvedProperty%Evolution%Partition%Fraction%Scalar + &
                                 ParticulateProperty%Evolution%Partition%Fraction%Scalar

                Error = abs(1. - TotalPartition)

                if (Error > 0.001)&
                    stop 'ConstructPartition - ModuleSedimentProperties - ERR05'

            endif

            nullify(DissolvedProperty, ParticulateProperty)

            Property => Property%Next
          
        enddo do1
          
        nullify(Property)

        !initialize system solver
        call StartLUD(Me%ObjLUD, SizeLB = 0, SizeUB = 3, &
                      WorkSizeLB = 1, WorkSizeUB = 2,    & 
                      STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPartition - ModuleSedimentProperties - ERR06'
            
        
        allocate(Me%PartitionMatrix(0:3, 0:3), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPartition - ModuleSedimentProperties - ERR07'
            

        allocate(Me%IndependentTerm(0:3     ), STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPartition - ModuleSedimentProperties - ERR08'

    end subroutine ConstructPartition


    !--------------------------------------------------------------------------


    subroutine Construct_Property(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer           :: NewProperty
        integer, intent(in)                 :: ClientNumber

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_)stop 'Construct_Property - ModuleSedimentProperties - ERR00'
        
        nullify(NewProperty%Prev,NewProperty%Next)
        nullify(NewProperty%Concentration        )
        nullify(NewProperty%BoundaryFlux         )
        nullify(NewProperty%HorizontalDiffusivity)
        nullify(NewProperty%VerticalDiffusivity  )


        call ConstructPropertyID            (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call Construct_PropertyState        (NewProperty)

        call Construct_PropertyValues       (NewProperty)

        call Construct_PropertyEvolution    (NewProperty, ClientNumber)

        call Construct_PropertyOutPut       (NewProperty)


    end subroutine Construct_Property

    !--------------------------------------------------------------------------

    subroutine Construct_SedimentRate(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_SedimentRate), pointer       :: NewSedimentRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewSedimentRate, STAT = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                          &
            stop 'Construct_SedimentRate - ModuleSedimentProperties - ERR01' 

        nullify(NewSedimentRate%Field, NewSedimentRate%Prev, NewSedimentRate%Next)

        call Construct_SedimentRateID       (NewSedimentRate)

        call Construct_SedimentRateValues   (NewSedimentRate)


    end subroutine Construct_SedimentRate

     !--------------------------------------------------------------------------
    
    !This subroutine reads all the information needed to construct the property ID          
    subroutine Construct_SedimentRateID(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_SedimentRate), pointer       :: NewSedimentRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: iflag, PropNumber
        logical                             :: CheckName
        type (T_Property), pointer          :: PropertyX
      
        !----------------------------------------------------------------------
        
        !First Property defined in a rate relation
        call GetData(NewSedimentRate%FirstProp%name,                                    &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'FIRSTPROP',                                        &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     SearchType   = FromBlock,                                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR01' 
        if (iflag==0)                                                                   &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR02' 

        !Check if the property name is valid
        CheckName = CheckPropertyName(NewSedimentRate%FirstProp%Name, Number = PropNumber)
        if (CheckName) then
            NewSedimentRate%FirstProp%IDnumber = PropNumber
        else
            write(*,*)
            write(*,*) 'The first property name is not recognised by the model.'
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR03' 
        end if 

        !call Search_Property(PropertyX, PropertyXID = PropNumber, STAT = STAT_CALL)                                  
        !if (STAT_CALL /= SUCCESS_)                                                     &
        !    stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR04' 
        
        !second Property defined in a rate relation
        call GetData(NewSedimentRate%SecondProp%name,                                    &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'SECONDPROP',                                        &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR05' 
        if (iflag==0)                                                                    &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR06' 
      
        ! Check if the property name is valid OR not
        CheckName = CheckPropertyName(NewSedimentRate%SecondProp%name, Number = PropNumber)
        if (CheckName) then
            NewSedimentRate%SecondProp%IDnumber = PropNumber
        else
            write(*,*)
            write(*,*) 'The Second property name is not recognised by the model.'
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR07' 
        end if
        
        call Search_Property(PropertyX,PropertyXID = PropNumber, STAT = STAT_CALL)                                  
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR08' 
  
  
        !Rate description ex: zooplankton grazing over phytoplankton
        call GetData(NewSedimentRate%ID%Description,                                     &
                     Me%ObjEnterData, iflag,                                             &
                     keyword      = 'DESCRIPTION',                                       &
                     Default      = 'No description was given.',                         &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     SearchType   = FromBlock,                                           &
                     STAT         = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR09' 

        call GetData(NewSedimentRate%ID%Name,                                              &
                     Me%ObjEnterData, iflag,                                               &
                     keyword      = 'NAME',                                                &
                     ClientModule = 'ModuleSedimentProperties',                            &
                     SearchType   = FromBlock,                                             &
                     Default      = 'No name was given to sediment rate.',                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                         &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR10' 
        if (iflag==0)                                                                      &
            stop 'Construct_SedimentRateID - ModuleSedimentProperties - ERR11' 

    end subroutine Construct_SedimentRateID

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct the property values       
    ! in the domain and in the boundaries            

    subroutine Construct_PropertyValues(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),              pointer      :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        integer                                     :: ILB,IUB
        integer                                     :: JLB,JUB
        integer                                     :: KLB,KUB
        integer                                     :: WorkSizeILB, WorkSizeIUB
        integer                                     :: WorkSizeJLB, WorkSizeJUB
        integer                                     :: WorkSizeKLB, WorkSizeKUB
        
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB
        WorkSizeKLB = Me%WorkSize%KLB
        WorkSizeKUB = Me%WorkSize%KUB

        allocate(NewProperty%Concentration(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR03'
             
        NewProperty%Concentration(:,:,:) = FillValueReal
        
        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR05'


        !This variable is to avoid negative concentration for properties which
        !have concentrations close to zero.
        call GetData(NewProperty%MinValue,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MIN_VALUE',                                        &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR06'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR07'
            NewProperty%Mass_Created(:,:,:) = 0.
        end if
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not.NewProperty%Old) then
            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill3D       = Me%ExternalVar%WaterPoints3D,     &
                                       Matrix3D             = NewProperty%Concentration,        &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR07'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR08'
            end if

            call CheckFieldConsistence (NewProperty)

        else

            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Water properties initial file written in HDF format  
            call ReadOldConcBoundariesHDF(NewProperty)

        end if

        !This coeficient can only be used when the 
        !relation between the IS units and units that 
        !the user wants to use is linear
        call GetData(NewProperty%IScoefficient,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword        = 'IS_COEF',                                        &
                     Default        = 1.e-3,                                            &      
                     SearchType     = FromBlock,                                        &
                     ClientModule   = 'ModuleSedimentProperties',                       &
                     STAT           = STAT_CALL)            
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_PropertyValues - ModuleSedimentProperties - ERR09' 

    end subroutine Construct_PropertyValues

    
    !--------------------------------------------------------------------------
    
    
    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer               :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: i,j,k
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        !Verification if the values read are lower than zero in water points
        do I = ILB, IUB
        do J = JLB, JUB
        do K = KLB, KUB
            
            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then
                               
                if (Me%ExternalVar%VolumeZ(i, j, k) == 0.) then

                    NewProperty%Concentration(i, j, k) = 0.
                    
                end if

            else

                NewProperty%Concentration(i, j, k) = FillValueReal

            endif

        enddo
        enddo
        enddo

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------

    
    subroutine Construct_SedimentRateValues(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_Sedimentrate), pointer       :: NewSedimentRate

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        allocate(NewSedimentRate%Field(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_SedimentRateValues - ModuleSedimentProperties - ERR01' 
        NewSedimentRate%Field(:,:,:) = FillValueReal

    end subroutine Construct_SedimentRateValues
    
    !--------------------------------------------------------------------------
    
    !This subroutine reads all the information needed to construct the property                          
    ! evolution parameters             
    subroutine Construct_PropertyEvolution(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer                   :: NewProperty
        integer, intent(in)                         :: ClientNumber

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        real                                        :: ErrorAux, AuxFactor, DTAux
        real                                        :: ModelDT

        !----------------------------------------------------------------------

        !By default the transport due to advection and Diffusion 
        !computed for all properties
        call GetData(NewProperty%Evolution%AdvectionDiffusion,                           &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'ADVECTION_DIFFUSION',                               &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR10'

        if(NewProperty%Evolution%AdvectionDiffusion)NewProperty%Evolution%Variable = .true.
        
        !if(NewProperty%Particulate .and. NewProperty%Evolution%AdvectionDiffusion)then
        !    write(*,*)
        !    write(*,*)'Particulate properties cannot have option ADVECTION_DIFFUSION'
        !    write(*,*)'activated in SedimentProperties. Please review your options in'
        !    write(*,*)'property '//trim(NewProperty%ID%Name)
        !    write(*,*)
        !    stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR20'
        !end if

        !Property afected by bioturbation
        call GetData(NewProperty%Evolution%ComputeBioturbation,                          &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'BIOTURBATION',                                      &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR30'

        if(NewProperty%Evolution%ComputeBioturbation)NewProperty%Evolution%Variable = .true.

        if(NewProperty%Evolution%ComputeBioturbation)then
            
            call Read_Bioturbation_Parameters   (NewProperty)

        end if

        if (NewProperty%Evolution%AdvectionDiffusion .or. NewProperty%Evolution%ComputeBioturbation)then

            call Read_Advec_Difus_Parameters    (NewProperty)

            call Construct_Property_Diffusivity (NewProperty)
        
        end if
        
        
          call GetData(NewProperty%Evolution%SeagrassesRoots,                          &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SEAGROOT',                                      &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR35'
            
        if(NewProperty%Evolution%SeagrassesRoots)NewProperty%Evolution%Variable = .true.

        !Property has partition as sink and source
        call GetData(NewProperty%Evolution%Partitioning,                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'PARTITION',                                         &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR40'

        if (NewProperty%Evolution%Partitioning)NewProperty%Evolution%Variable = .true.

        !Partition parameters  
        if (NewProperty%Evolution%Partitioning)                                          &
            call Read_Partition_Parameters (NewProperty, ClientNumber)

        ! Property has the sediment quality model as sink and source
        call GetData(NewProperty%Evolution%SedimentQuality,                              &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SEDIMENT_QUALITY',                                  &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR50'

        if (NewProperty%Evolution%SedimentQuality)NewProperty%Evolution%Variable = .true.
        

        !This property has fluxes at the sediment water interface? no - 0;  yes - 1
        call GetData(NewProperty%Evolution%SurfaceFluxes,                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SURFACE_FLUXES',                                    &
                     ClientModule = 'ModuleSedimentProperties',                          &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR70'
        
        if (NewProperty%Evolution%SurfaceFluxes)NewProperty%Evolution%Variable = .true.
            

        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR80'


            Me%ExternalVar%DT = ModelDT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DT_INTERVAL',                                   &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModuleSedimentProperties',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR90'
                                       
            
            if (NewProperty%Evolution%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR100'

            elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR110'
                endif

                !Run period in seconds
                DTaux = Me%EndTime - Me%ExternalVar%Now

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%Evolution%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'Construct_PropertyEvolution - ModuleSedimentProperties - ERR120'
                end if
            endif

            NewProperty%Evolution%NextCompute = Me%ExternalVar%Now + NewProperty%Evolution%DTInterval

        else

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%Evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution     


    !--------------------------------------------------------------------------


    subroutine Read_Advec_Difus_Parameters(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: NullFlux           
        integer                                 :: ImposedValue              
        integer                                 :: NullGradient                        

        !Local-----------------------------------------------------------------
        integer                                 :: iflag, BoundaryCondition

        !----------------------------------------------------------------------
        
        ! molecular diffusion coefficient corrected for sediments compute method
        call GetData(NewProperty%Evolution%Advec_Difus_Parameters%Diffusion_Method,     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DIFFUSION_METHOD',                                 &
                     Default      = Berner_1980,                                        &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Advec_Difus_Parameters - ModuleSedimentProperties - ERR01'


        !molecular diffusion coefficient
        call GetData(NewProperty%Evolution%Advec_Difus_Parameters%Molecular_Diff_Coef,  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MOLECULAR_DIFF_COEF',                              &
                     Default      = 0.0,                                                &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Advec_Difus_Parameters - ModuleSedimentProperties - ERR02'

        call GetBoundaryConditionList(MassConservation           = NullFlux,            &
                                      ImposedValue               = ImposedValue,        &
                                      NullGradient               = NullGradient)

        !Multiple Options : 1-NullFlux,2-ImposedValue,4-NullGradient
        call GetData(BoundaryCondition,                                                 &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BOUNDARY_CONDITION',                               &
                     Default      = NullFlux,                                           &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Advec_Difus_Parameters - ModuleSedimentProperties - ERR12'

        ! By default it's imposed a value dependent only from the exterior
        ! value and of the decay time. However this method doesn't conserve mass 
        ! when the water fluxes near the frontier are dominant

        if (BoundaryCondition /= NullFlux                      .and.                    &
            BoundaryCondition /= ImposedValue                  .and.                    &
            BoundaryCondition /= NullGradient)                                          &
            stop 'Read_Advec_Difus_Parameters - ModuleSedimentProperties - ERR13'

        NewProperty%Evolution%Advec_Difus_Parameters%BoundaryCondition = BoundaryCondition


    end subroutine Read_Advec_Difus_Parameters


    !--------------------------------------------------------------------------


    subroutine Read_Bioturbation_Parameters(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: iflag

        !Local-------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        
        !Begin-------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        !Standard bioturbation diffusion coefficient
        call GetData(NewProperty%Evolution%Bioturbation%DefaultCoef,                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BIOTURBATION_COEF',                                &
                     Default      = 0.,                                                 &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Bioturbation_Parameters - ModuleSedimentProperties - ERR10'


        !Depth were bioturbation is important
        call GetData(NewProperty%Evolution%Bioturbation%BioDepth,                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BIOTURBATION_DEPTH',                               &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Bioturbation_Parameters - ModuleSedimentProperties - ERR20'

        !Coefficient to compute bioturbation effect decay in depth
        call GetData(NewProperty%Evolution%Bioturbation%DecayCoef,                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BIOTURBATION_DECAY_COEF',                          &
                     Default      = 0.01,                                               &
                     ClientModule = 'ModuleSedimentProperties',                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Read_Bioturbation_Parameters - ModuleSedimentProperties - ERR30'
        if(NewProperty%Evolution%Bioturbation%DecayCoef .le. 0)then
            write(*,*)'BIOTURBATION_DECAY_COEF cannot be ZERO or NEGATIVE'
            stop 'Read_Bioturbation_Parameters - ModuleSedimentProperties - ERR40'
        end if
        
        
        allocate(NewProperty%Evolution%Bioturbation%Coef(ILB:IUB, JLB:JUB, KLB:KUB),STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Read_Bioturbation_Parameters - ModuleSedimentProperties - ERR50' 
        NewProperty%Evolution%Bioturbation%Coef(:,:,:) = null_real
        
    end subroutine Read_Bioturbation_Parameters
    
    !--------------------------------------------------------------------------

    
    subroutine Construct_Property_Diffusivity (NewProperty)

        !Arguments---------------------------------------------------------
        type(T_Property),    pointer                :: NewProperty

        !Local-------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: STAT_CALL
        
        !Begin-------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        allocate(NewProperty%HorizontalDiffusivity  (ILB:IUB, JLB:JUB, KLB:KUB),STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_Property_Diffusivity - ModuleSedimentProperties - ERR01' 
        NewProperty%HorizontalDiffusivity         (:,:,:) = 0.

        allocate(NewProperty%VerticalDiffusivity    (ILB:IUB, JLB:JUB, KLB:KUB),STAT = STAT_CALL)  
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Construct_Property_Diffusivity - ModuleSedimentProperties - ERR02' 
        NewProperty%VerticalDiffusivity           (:,:,:) = 0.

    end subroutine Construct_Property_Diffusivity     

    
    !--------------------------------------------------------------------------
    
    
    subroutine Read_Partition_Parameters(NewProperty, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_property),           pointer     :: NewProperty
        integer                                 :: ClientNumber

        !External--------------------------------------------------------------
        integer                                 :: iflag
        integer                                 :: STAT_CALL
        integer                                 :: ILB, IUB
        integer                                 :: JLB, JUB
        integer                                 :: KLB, KUB
        
        !Local-----------------------------------------------------------------
        type(T_Property_3D), pointer            :: Scalar3D

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        Scalar3D => NewProperty%Evolution%Partition%Fraction
        call ConstructScalar3D(Scalar3D,                                            &
                               ExtractType = FromBlockInBlock,                      &
                               ClientNumber= ClientNumber,                          &
                               block_begin = '<<begin_partition_fraction>>',        &
                               block_end   = '<<end_partition_fraction>>')

        Scalar3D => NewProperty%Evolution%Partition%Rate
        call ConstructScalar3D(Scalar3D,                                            &
                               ExtractType  = FromBlockInBlock,                     &
                               ClientNumber = ClientNumber,                         &
                               block_begin  = '<<begin_partition_rate>>',           &
                               block_end    = '<<end_partition_rate>>')
        
        ! Name of property (dissolved/particulated) to couple  
        call GetData(NewProperty%Evolution%Partition%Couple,                        &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'PARTITION_COUPLE',                             & 
                     ClientModule ='ModuleSedimentProperties',                      &
                     SearchType   = FromBlock,                                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
            stop 'Read_Partition_Parameters - ModuleSedimentProperties - ERR10'
        if (iflag .NE. 1)                                                           &
            stop 'Read_Partition_Parameters - ModuleSedimentProperties - ERR20'

    end subroutine Read_Partition_Parameters
    
    !------------------------------------------------------------------------------

    subroutine Construct_PropertyState(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL, iflag

        !----------------------------------------------------------------------
        
        call GetData(NewProperty%Particulate,                                       &
                     Me%ObjEnterData,  iflag,                                       &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'PARTICULATE',                                  &
                     ClientModule = 'ModuleSedimentProperties',                     &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModuleSedimentProperties - ERR01'
        if(iflag == 0)              stop 'Construct_PropertyState - ModuleSedimentProperties - ERR02'

        if (NewProperty%Particulate)then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'Construct_PropertyState - ModuleSedimentProperties - ERR03'
            end if
        endif

    end subroutine Construct_PropertyState

    !--------------------------------------------------------------------------

    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),    pointer        :: NewProperty

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        !Checks out if the user pretends to write a time serie for this property
        call GetData(NewProperty%TimeSerie,                         &
                     Me%ObjEnterData, iflag,                        &
                     Keyword      = 'TIME_SERIE',                   &
                     ClientModule = 'ModuleSedimentProperties',     &
                     Default      = .false.,                        &
                     SearchType   = FromBlock,                      &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSedimentProperties - ERR01'
        
        ! Checks out if the user pretends to write a 
        ! time serie inside each box for this property
        call GetData(NewProperty%BoxTimeSerie,                      &
                     Me%ObjEnterData, iflag,                        &
                     Keyword      = 'BOX_TIME_SERIE',               &
                     Default      = .false.,                        &
                     SearchType   = FromBlock,                      &
                     ClientModule = 'ModuleSedimentProperties',     &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSedimentProperties - ERR02'

        !Checks out if the user pretends to write a outputs 
        !in HDF for this property
        call GetData(NewProperty%OutputHDF,                         &
                     Me%ObjEnterData, iflag,                        &
                     Keyword      = 'OUTPUT_HDF',                   &
                     ClientModule = 'ModuleSedimentProperties',     &
                     Default      = .false.,                        &
                     SearchType   = FromBlock,                      &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSedimentProperties - ERR03'
        
    end subroutine Construct_PropertyOutPut

    !--------------------------------------------------------------------------
    

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

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
            stop 'ConstructGlobalOutput - ModuleSedimentProperties - ERR01' 

        if (Me%OutPut%Yes) then

            Me%OutPut%NextOutPut = 1

        else
            write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
            write(*,*)'one property has HDF format outputs.'
            stop 'ConstructGlobalOutput - ModuleSedimentProperties - ERR02'
        endif 


    end subroutine ConstructGlobalOutput


    !--------------------------------------------------------------------------

    
    subroutine ReadSedimentPropertiesFilesName

        !External--------------------------------------------------------------
        integer                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character(len=StringLength) :: Message
        
        !Begin-----------------------------------------------------------------
        
        
        call ReadFileName('SED_DAT', Me%Files%InPutData, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) &
            stop 'ReadSedimentPropertiesFilesName - ModuleSedimentProperties - ERR01'


        Message   = 'Instant fields of sediment properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SED_HDF', Me%Files%HDFResults,               &
                          Message = Message, Time_end = Me%EndTime,     &
                          Extension = 'sed', STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) &
            stop 'ReadSedimentPropertiesFilesName - ModuleSedimentProperties - ERR02'


        Message   ='Sediment properties final values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SED_FIN', Me%Files%Final,                    &
                          Message = Message, Time_end = Me%EndTime,     &
                          Extension = 'sef', STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) &
            stop 'ReadSedimentPropertiesFilesName - ModuleSedimentProperties - ERR03'

        !Sediment properties initial values in HDF format
        Message   ='sediment properties initial values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SED_INI', Me%Files%Initial,                  &
                           Message = Message, Time_end = Me%EndTime,    &
                           STAT = STAT_CALL)
cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_   ) then 

            call SetError(FATAL_, INTERNAL_,                            & 
            'Initial file not found - ReadSedimentPropertiesFilesName - ModuleSedimentProperties - ERR04')

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then cd1

            call SetError(WARNING_, KEYWORD_,                           &
            'Keyword SED_INI not found - ReadSedimentPropertiesFilesName - ModuleSedimentProperties - ERR05', &
            Screen = .false.)

        else if (STAT_CALL .EQ. SUCCESS_              ) then cd1

            continue
        
        else cd1

            call SetError(FATAL_, INTERNAL_, &
                'Calling ReadFileName - ReadSedimentPropertiesFilesName; ModuleSedimentProperties. ERR06.')

        end if cd1  
                                                                             
    end subroutine ReadSedimentPropertiesFilesName

    !--------------------------------------------------------------------------
    !If the user want's to use the values of a previous   
    ! run the read the property values form the final      
    ! results file of a previous run. By default this      
    ! file is in HDF format                                

    subroutine ReadOldConcBoundariesHDF(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character (Len=StringLength)                :: PropertyName
        logical                                     :: Exist
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        !----------------------------------------------------------------------

        inquire (FILE=trim(Me%Files%Initial)//"5", EXIST = Exist)

cd0:    if (Exist) then

            ObjHDF5 = 0

            !Gets File Access Code
            call GetHDF5FileAccess(HDF5_READ = HDF5_READ)

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5, trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOldConcBoundariesHDF - ModuleSedimentProperties - ERR01'

            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Concentration(:,:,:) = FillValueReal

            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                                 &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOldConcBoundariesHDF - ModuleSedimentProperties - ERR02'
                

            call HDF5ReadData   (ObjHDF5, "/Concentration/"//PropertyName,                  &
                                 PropertyName,                                              &
                                 Array3D = NewProperty%Concentration,                       &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOldConcBoundariesHDF - ModuleSedimentProperties - ERR03'
                
            if (Me%OutPut%Yes) then
                call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadOldConcBoundariesHDF - ModuleSedimentProperties - ERR04'
            endif
                

        else
            
            stop 'ReadOldConcBoundariesHDF - ModuleSedimentProperties - ERR05'

        end if cd0

    end subroutine ReadOldConcBoundariesHDF


    !--------------------------------------------------------------------------
    
    ! This subroutine adds a new property to the Water Property List  
    subroutine Add_Property(NewProperty)

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


    end subroutine Add_Property 


    !--------------------------------------------------------------------------

    
    subroutine Add_SedimentRate(NewSedimentRate)

        !Arguments-------------------------------------------------------------
        type(T_SedimentRate), pointer       :: NewSedimentRate

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstSedimentRate)) then
            Me%SedimentRatesNumber      = 1
            Me%FirstSedimentRate        => NewSedimentRate
            Me%LastSedimentRate         => NewSedimentRate
        else
            NewSedimentRate%Prev        => Me%LastSedimentRate
            Me%LastSedimentRate%Next    => NewSedimentRate
            Me%LastSedimentRate         => NewSedimentRate
            Me%SedimentRatesNumber      = Me%SedimentRatesNumber + 1
        end if 

    end subroutine Add_SedimentRate 

    
    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine SedimentProperties_Evolution (SedimentPropertiesID, STAT)
                                              
        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        integer, optional, intent(OUT)              :: STAT
        
        !External-----------------------------------------------------------------
        integer                                     :: ready_
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_)

cd0 :   if (ready_ .EQ. IDLE_ERR_) then

            call TimeStepActualization

            call ReadLockExternalVar

            if (Me%Coupled%Bioturbation%Yes)            &
                call Bioturbation_Processes

            if (Me%Coupled%SurfaceFluxes%Yes)           &
                call Modify_SurfaceBoundaryFluxes

            if (Me%Coupled%AdvectionDiffusion%Yes)      &
                call Advection_Diffusion_Processes
                
            if (Me%Coupled%SeagrassesRoots%Yes)         &     ! Isabella
                call SeagrassesRoots_Processes


            !if (Me%Coupled%SedimentQuality%Yes)         &
            !    call SedimentQuality_Processes

            if (Me%Coupled%Partition%Yes)               &
                call Partition_Processes
            
            if (Me%Coupled%MinimumConcentration%Yes)    &
                call SetMinimumConcentration

            if (Me%Coupled%TimeSerie%Yes)               &
                call OutPut_TimeSeries

            if (Me%Coupled%BoxTimeSerie%Yes)            &
                call OutPut_BoxTimeSeries
            
            if (Me%Coupled%OutputHDF%Yes)               &
                call OutPut_Results_HDF
            
            call Actualize_Time_Evolution

            call ReadUnlockExternalVar

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd0

        if (present(STAT)) STAT = STAT_        

    end subroutine SedimentProperties_Evolution

    !--------------------------------------------------------------------------


    subroutine SedimentQuality_Processes
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property),          pointer     :: PropertyX
        type (T_SedimentRate),      pointer     :: SedimentRateX
        integer                                 :: i, j, k
        
        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "SedimentQuality_Processes")


        !Computes residual conversion factor
        call ComputeDissolvedToParticulate3D

        if (Me%ExternalVar%Now .GE. Me%Coupled%SedimentQuality%NextCompute) then
            
            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))

                call Modify_Interface(InterfaceID               = Me%ObjInterface,                  &
                                      PropertyID                = PropertyX%ID%IDNumber,            &
                                      Concentration             = PropertyX%Concentration,          &
                                      WaterPoints3D             = Me%ExternalVar%WaterPoints3D,     &
                                      OpenPoints3D              = Me%ExternalVar%OpenPoints3D,      &
                                      WaterPercentage           = Me%ExternalVar%WaterPercentage,   &
                                      DissolvedToParticulate3D  = Me%DissolvedToParticulate3D,      &
                                      STAT                      = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                        &
                    stop 'SedimentQuality_Processes - ModuleSedimentProperties - ERR01'

                PropertyX => PropertyX%Next

            end do
            
            Me%Coupled%SedimentQuality%NextCompute = Me%Coupled%SedimentQuality%NextCompute +       &
                                                     Me%Coupled%SedimentQuality%DT_Compute

        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if (PropertyX%Evolution%SedimentQuality) then

                if (Me%ExternalVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%Concentration,          &
                                          WaterPoints3D = Me%ExternalVar%WaterPoints3D,     &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'SedimentQuality_Processes - ModuleSedimentProperties - ERR02'

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        SedimentRateX => Me%FirstSedimentRate

        do while (associated(SedimentRateX))

            call GetRateFlux(InterfaceID    = Me%ObjInterface,                          &
                             FirstProp      = SedimentRateX%FirstProp%IDNumber,         &
                             SecondProp     = SedimentRateX%SecondProp%IDNumber,        &
                             RateFlux3D     = SedimentRateX%Field,                      &
                             WaterPoints3D  = Me%ExternalVar%WaterPoints3D,             &
                             STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'SedimentQuality_Processes - ModuleSedimentProperties - ERR03'

            call Search_Property(PropertyX,                                             &
                                 PropertyXID = SedimentRateX%FirstProp%IDnumber,        &
                                 STAT = STAT_CALL)                      
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'SedimentQuality_Processes - ModuleSedimentProperties - ERR04'

            if(PropertyX%Particulate)then

                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if(Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                        SedimentRateX%Field(i,j,k) = SedimentRateX%Field(i,j,k)  * &
                                       Me%ExternalVar%DrySedimentVolume (i,j,k)  * &
                                       Me%Sediment_DryDensity%Field(i,j,k)       / &
                                       Me%Coupled%SedimentQuality%DT_Compute
                    end if
                end do
                end do
                end do

            else

                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                    if(Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                        SedimentRateX%Field(i,j,k) = SedimentRateX%Field(i,j,k)        * &
                                                     Me%ExternalVar%WaterVolume (i,j,k)/ &
                                                     Me%Coupled%SedimentQuality%DT_Compute

                    end if
                end do
                end do
                end do

            end if

            call BoxDif(Me%ObjBoxDif,                                       &
                        SedimentRateX%Field,                                &
                        SedimentRateX%ID%Name,                              &
                        Me%ExternalVar%WaterPoints3D,                       &
                        STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                    &
                stop 'SedimentQuality_Processes - ModuleSedimentProperties - ERR05'

            SedimentRateX => SedimentRateX%Next
        enddo

        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "SedimentQuality_Processes")

    end subroutine SedimentQuality_Processes

    
    !--------------------------------------------------------------------------
       !--------------------------------------------------------------------------
    subroutine SeagrassesRoots_Processes
        !Local----------------------------------------------------------------- 
        type (T_Property),       pointer         :: PropertyX
        type (T_SedimentRate  ), pointer         :: SedimentRateX
        integer                                  :: i, j, k, KTop
        integer                                  :: ILB, IUB, JLB, JUB,KLB, KUB
        integer                                  :: STAT_CALL


        

        !Begin----------------------------------------------------------------- 

            ILB = Me%WorkSize%ILB 
            IUB = Me%WorkSize%IUB 
            JLB = Me%WorkSize%JLB 
            JUB = Me%WorkSize%JUB 
            KUB = Me%WorkSize%KUB
            KLB = Me%WorkSize%KLB
        
        
        

        
        PropertyX => Me%FirstProperty

       do while (associated(PropertyX))

                if (PropertyX%Evolution%SeagrassesRoots) then
                
                    if(PropertyX%ID%IDNumber == SeagrassesRoots_)then
                                    
                        call RootsOccupation(Me%SeagrassesRoots)
                        call DistributeRoots(PropertyX, Me%SeagrassesRoots)
                       
                    end if
                
               
                
                end if 
                
                PropertyX=>PropertyX%Next
       end do
       
       
            !CHUNK = CHUNK_J(JLB, JUB)
            
            ! total occupation

            if (Me%WorkSize%KUB >1) then
            
                        do j=JLB, JUB
                        do i=ILB, IUB
                        do k=KLB, KUB
                        
                        
                             if (Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
                           
                                !factors assumed to be the constant over the sediment column 
                                Me%SeagrassesRoots%NintFactor3DR(i,j,k)=Me%SeagrassesRoots%NintFactor2DR(i,j)
                                Me%SeagrassesRoots%PintFactor3DR(i,j,k)=Me%SeagrassesRoots%PintFactor2DR(i,j)
                                ! Roots mortality was calculated in ModuleBenthicEcology
                                ! the roots mortality is a flux that will be added to the bottom cell of the sediment column
                                ! in module SeagrassSedimentInteraction
                                if (k==KLB) then
                                Me%SeagrassesRoots%RootsMort3DR(i,j,k)=Me%SeagrassesRoots%RootsMort2DR(i,j)  !Kg DW/day
                                else                                      
                                Me%SeagrassesRoots%RootsMort3DR(i,j,k)= 0.   !Kg DW/day
                                endif     
                            
                            endif   
                        
                        enddo
                        enddo
                        enddo
            
            else
            
            ! 2D
                        do j=JLB, JUB
                        do i=ILB, IUB

                        KTop                 = Me%ExternalVar%KTop(i, j)
                        
                        if (Me%ExternalVar%OpenPoints3D(i, j, KTop) == OpenPoint) then
                       
                        !factors assumed to be the constant over the sediment column 
                        Me%SeagrassesRoots%NintFactor3DR(i,j,KTop)=Me%SeagrassesRoots%NintFactor2DR(i,j)
                        Me%SeagrassesRoots%PintFactor3DR(i,j,KTop)=Me%SeagrassesRoots%PintFactor2DR(i,j)
                        !kgrams of dry weight per day
                        ! if 2d, it can only be added to the top layer of the sediment
                        Me%SeagrassesRoots%RootsMort3DR(i,j,KTop)=Me%SeagrassesRoots%RootsMort2DR(i,j) 
                                                               
                        
                        endif        
                        
                        enddo
                        enddo

            
            endif

            
            PropertyX => Me%FirstProperty

        if (Me%ExternalVar%Now .GE. Me%Coupled%SeagrassesRoots%NextCompute) then

      
   


          
            do while(associated(PropertyX))
            
                if(PropertyX%ID%IDNumber == SeagrassesRoots_)then
                
                Me%SeagrassesRoots%Volume => Me%ExternalVar%VolumeZ
            
                   call Modify_Interface(InterfaceID    = Me%ObjSeagrassSedimInteraction,        &
                                      PropertyID        = PropertyX%ID%IDNumber,                 &
                                      Concentration     = PropertyX%Concentration,               &
                                      WaterPoints3D     = Me%ExternalVar%WaterPoints3D,          &
                                      OpenPoints3D      = Me%ExternalVar%OpenPoints3D,           &
                                      DWZ               = Me%ExternalVar%DWZ,                    &
                                      NintFac3DR        = Me%SeagrassesRoots%NintFactor3DR,      & 
                                      PintFac3DR        = Me%SeagrassesRoots%PintFactor3DR,      & 
                                      SedimCellVol3D    = Me%SeagrassesRoots%Volume,             &
                                      RootsMort         = Me%SeagrassesRoots%RootsMort3DR,      &
                                      STAT              = STAT_CALL)
                     if (STAT_CALL .NE. SUCCESS_)                                                &
                     stop 'SeagrassesRoots_Processes - ModuleSedimentProperties - ERR02'
           
              else
           
                     call Modify_Interface(InterfaceID  = Me%ObjSeagrassSedimInteraction,   &
                                      PropertyID        = PropertyX%ID%IDNumber,                 &
                                      Concentration     = PropertyX%Concentration,               &
                                      WaterPoints3D     = Me%ExternalVar%WaterPoints3D,          &
                                      OpenPoints3D      = Me%ExternalVar%OpenPoints3D,           & 
                                      DWZ               = Me%ExternalVar%DWZ,               &
                                      STAT              = STAT_CALL)
                         if (STAT_CALL .NE. SUCCESS_)                                                &
                         stop 'SeagrassesRoots_Processes - ModuleSedimentProperties - ERR02'
           
             endif
           
                PropertyX => PropertyX%Next

            end do

            Me%Coupled%SeagrassesRoots%NextCompute = Me%Coupled%SeagrassesRoots%NextCompute + &
                                                      Me%Coupled%SeagrassesRoots%DT_Compute
            
        end if 




        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%Evolution%SeagrassesRoots) then
    
                if (Me%ExternalVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    call Modify_Interface(InterfaceID   = Me%ObjSeagrassSedimInteraction,  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%Concentration,          &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          WaterPoints3D = Me%ExternalVar%WaterPoints3D,     &
                                          OpenPoints3D  = Me%ExternalVar%OpenPoints3D,      &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'SeagrassesRoots_Processes - ModuleSedimentProperties - ERR03'

                        if(PropertyX%ID%IDNumber == SeagrassesRoots_)then
                        !Integrate roots distribution in the sediment column in to kgC/m2
                        !call IntegrateRoots(PropertyX, Me%SeagrassesRoots)
                    
                    end if
                endif

            endif
            
            PropertyX=>PropertyX%Next

        enddo


        SedimentRateX => Me%FirstSedimentRate
        
        do while (associated(SedimentRateX))

            !if(SedimentRateX%Model == 'SeagrassSedimInteraction')then

                call GetRateFlux(InterfaceID    = Me%ObjSeagrassSedimInteraction,          &
                                 FirstProp      = SedimentRateX%FirstProp%IDNumber,               &
                                 SecondProp     = SedimentRateX%SecondProp%IDNumber,              &
                                 RateFlux3D     = SedimentRateX%Field,                            &
                                 WaterPoints3D  = Me%ExternalVar%WaterPoints3D,             &
                                 STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                &
                    stop 'SeagrassesRoots_Processes - ModuleSedimentProperties - ERR04'

          SedimentRateX%Field2 =>SedimentRateX%Field

          if (SedimentRateX%FirstProp%IDNumber==RootsUptakeN_) then
             
             
          Me%SeagrassesRoots%UptakeNH4s3D =SedimentRateX%Field  !gN/day 
                         
          elseif (SedimentRateX%FirstProp%IDNumber==RootsUptakeP_) then
          
                  Me%SeagrassesRoots%UptakePO4s3D =SedimentRateX%Field  !gN/day 
          else
                    
                where (Me%ExternalVar%WaterPoints3D == WaterPoint) &
                     SedimentRateX%Field =  SedimentRateX%Field * Me%ExternalVar%VolumeZ / &
                                    Me%Coupled%SeagrassesRoots%DT_Compute

                call BoxDif(Me%ObjBoxDif,                                                   &
                             SedimentRateX%Field,                                                  &
                            trim(SedimentRateX%ID%Name),                                          &
                            Me%ExternalVar%WaterPoints3D,                                   &
                            STAT  = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                &
                    stop 'SeagrassesRoots_Processes - ModuleSedimentProperties - ERR05'

            end if  
            
            
            
            
        !endif
             SedimentRateX=> SedimentRateX%Next

        enddo

        nullify( SedimentRateX)
        nullify(PropertyX)
     
     
     end   subroutine SeagrassesRoots_Processes
 
    
    !---------------------------------------------------------------------------------------------

    
    subroutine ComputeDissolvedToParticulate3D

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
         
        !Local----------------------------------------------------------------- 
        integer                                 :: i, j, k
        real                                    :: DT, InstantValue, ResidualValue

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "ComputeDissolvedToParticulate3D")


        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'ComputeDissolvedToParticulate3D - ModuleSedimentProperties - ERR01'


        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                InstantValue                       = Me%ExternalVar%WaterVolume(i,j,k)   / &
                                                    (Me%Sediment_DryDensity%Field(i,j,k) * &
                                                     Me%ExternalVar%DrySedimentVolume(i,j,k))

                ResidualValue                      = Me%DissolvedToParticulate3D(i,j,k)

                Me%DissolvedToParticulate3D(i,j,k) = (ResidualValue * Me%ResidualTime +     &
                                                      InstantValue * DT) / (Me%ResidualTime + DT)
                                                       
            end if
        end do
        end do
        end do

        Me%ResidualTime = Me%ResidualTime + DT
        
        
        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "ComputeDissolvedToParticulate3D")


    end subroutine ComputeDissolvedToParticulate3D

    !--------------------------------------------------------------------------
    ! This subroutine is responsable for computing the  
    ! Advection Diffusion processes of the dissolved properties                      

    subroutine Advection_Diffusion_Processes

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer           :: Property

        !External--------------------------------------------------------------
        !integer                             :: STAT_CALL    

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: WIUB,WJUB,WILB,WJLB, WKUB,WKLB
        integer                             :: i, j, k
        !real,    dimension(:    ), pointer  :: FluxZ
        !real                                :: Area
        real                                ::DT, BottomFace, TopFace
        !real(8), dimension(:,:,:), pointer  :: AdvFluxZ
        !real(8), dimension(:,:,:), pointer  :: DifFluxZ
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "Advection_Diffusion_Processes")

        WIUB = Me%WorkSize%IUB; ILB = Me%Size%ILB
        WJUB = Me%WorkSize%JUB; IUB = Me%Size%IUB
        WILB = Me%WorkSize%ILB; JLB = Me%Size%JLB 
        WJLB = Me%WorkSize%JLB; JUB = Me%Size%JUB 
        WKUB = Me%WorkSize%KUB; KLB = Me%Size%KLB
        WKLB = Me%WorkSize%KLB; KUB = Me%Size%KUB
        
        Property => Me%FirstProperty

        do while(associated(Property))

            if (Property%Evolution%AdvectionDiffusion) then

                if(Me%ExternalVar%Now .ge. Property%Evolution%NextCompute)then

                    call ComputeDiffusivity(Property)

                    DT = Property%Evolution%DTInterval

                    do j = WJLB, WJUB
                    do i = WILB, WIUB

                        WKUB = Me%ExternalVar%KTop(i, j)
                        KUB  = WKUB + 1

                        !Area                      = Me%ExternalVar%GridCellArea     (i, j         )
                        Me%Concentration(KLB:KUB) = Property%Concentration          (i, j, KLB:KUB)
                        Me%DiffCoef     (KLB:KUB) = Property%VerticalDiffusivity    (i, j, KLB:KUB)

                        Me%DWZ           = Me%ExternalVar%DWZ             (i, j, KLB:KUB)
                        Me%DZZ           = Me%ExternalVar%DZZ             (i, j, KLB:KUB)
                        !Me%VolZ          = Me%ExternalVar%WaterVolume     (i, j, KLB:KUB)
                        !Me%VolZOld       = Me%ExternalVar%WaterVolumeOld  (i, j, KLB:KUB)

                        if(Me%ExternalVar%ComputeFacesW3D(i, j, WKUB) == 1)then

                            !null gradient
                            Me%Concentration(WKLB - 1) = Property%Concentration(i, j, WKLB)
                            Me%Concentration(WKUB + 1) = Property%Concentration(i, j, WKUB)

                            do k = WKLB, WKUB

                                !Me%Ti_Coef(k)  = Me%Concentration(k) * Me%VolZOld(k) / Me%VolZ(k)
                                Me%Ti_Coef(k)  = Me%Concentration(k)

                            enddo

                            Me%D_Coef (WKLB) = 0.
                            Me%E_Coef (WKLB) = 1.
                            Me%F_Coef (WKLB) = 0.

                            if(Property%Particulate)then

                                do k = WKLB+1, WKUB

                                    BottomFace = Me%DiffCoef(k) * Me%DWZ(k) * (1.- Me%ExternalVar%Porosity(i,j,k)) * &
                                                 DT / (Me%DZZ(k-1))
                                    TopFace    = Me%DiffCoef(k) * Me%DWZ(k) * (1.- Me%ExternalVar%Porosity(i,j,k)) * &
                                                 DT / (Me%DZZ(k  ))

                                    Me%D_Coef (k) = - BottomFace
                                    Me%E_Coef (k) =   BottomFace + TopFace + 1
                                    Me%F_Coef (k) = - TopFace

                                enddo


                            else
                                
                                do k = WKLB+1, WKUB

                                    BottomFace = Me%DiffCoef(k) * Me%DWZ(k) * Me%ExternalVar%Porosity(i,j,k) * &
                                                 DT / (Me%DZZ(k-1))
                                    TopFace    = Me%DiffCoef(k) * Me%DWZ(k) * Me%ExternalVar%Porosity(i,j,k) * &
                                                 DT / (Me%DZZ(k  ))

                                    Me%D_Coef (k) = - BottomFace
                                    Me%E_Coef (k) =   BottomFace + TopFace + 1
                                    Me%F_Coef (k) = - TopFace

                                enddo



                            end if

                            Me%Ti_Coef(WKLB) = Me%Ti_Coef(WKLB) - Me%D_Coef(WKLB) * Me%Concentration(WKLB)
                            Me%Ti_Coef(WKUB) = Me%Ti_Coef(WKUB) - Me%F_Coef(WKUB) * Me%Concentration(WKUB)

                            call Tridiagonal(Me%Concentration, WKLB, WKUB)

                            do k = WKLB, WKUB

                                Property%Concentration(i, j, k) = Me%Concentration(k)

                            enddo

                        endif

                    enddo
                    enddo

                    !Wflux_Z    = Me%ExternalVar%WaterFluxZ
                    !VolumeZOld = Me%ExternalVar%WaterVolumeOld
                    !VolumeZ    = Me%ExternalVar%WaterVolume

                    if (Property%BoxTimeSerie) then

                        !do K = KLB, KUB
                        !do J = JLB, JUB
                        !do I = ILB, IUB
                            
                        !    if(Me%ExternalVar%WaterPoints3D(i,j,k) == 1)then

                                !Me%MassFluxesZ(i,j,k) = AdvFluxZ(i,j,k) + DifFluxZ (i,j,k)

                         !   endif

                        !end do
                        !end do
                        !end do

                       !call BoxDif(Me%ObjBoxDif,                                               &
                        !            Me%MassFluxesX,                                             &
                        !            Me%MassFluxesY,                                             &
                        !            Me%MassFluxesZ,                                             &
                        !            "Sed_"//trim(adjustl(Property%ID%name)),                    &
                        !            Me%ExternalVar%WaterPoints3D,                               &
                        !            STAT = STAT_CALL)
                        !if (STAT_CALL /= SUCCESS_)  &
                         !   stop 'Advection_Diffusion_Processes - ModuleSedimentProperties - ERR04'

                    end if
            
                end if

            end if

            Property => Property%Next

        end do
        
        
        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "Advection_Diffusion_Processes")

    end subroutine Advection_Diffusion_Processes


    !--------------------------------------------------------------------------

    subroutine TriDiagonal(C, KLB, KUB)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:), pointer         :: C
        integer, intent(in)                 :: KLB, KUB

        !Local-----------------------------------------------------------------
        real                                :: x
        integer                             :: i, j

        !Begin-----------------------------------------------------------------

        do i = KLB+1, KUB
            x             = Me%D_Coef(i)  / Me%E_Coef (i - 1)
            Me%E_Coef(i)  = Me%E_Coef(i)  - Me%F_Coef (i - 1) * x
            Me%Ti_Coef(i) = Me%Ti_Coef(i) - Me%Ti_Coef(i - 1) * x
        enddo

        C(KUB) = Me%Ti_Coef(KUB) / Me%E_Coef(KUB)

        do i = KLB, KUB-1
            j    = KUB - i
            C(j) = (Me%Ti_Coef(j) - Me%F_Coef(j) * C(j + 1)) / Me%E_Coef(j)
        enddo

    end subroutine TriDiagonal



    !--------------------------------------------------------------------------
    ! This subroutine is responsable for computing the  
    ! Diffusion processes of the particulate properties
    ! depending of bioturbation effects

    subroutine Bioturbation_Processes

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer           :: Property

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------


        !if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "Bioturbation_Processes")

        Property => Me%FirstProperty

        do while(associated(Property))

            if (Property%Evolution%ComputeBioturbation) then

                if(Me%ExternalVar%Now .ge. Property%Evolution%NextCompute)then

                    call Compute_Bioturbation(Property)

                end if

            end if

            Property => Property%Next

        end do

        !if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "Bioturbation_Processes")


    end subroutine Bioturbation_Processes


    !--------------------------------------------------------------------------
    
    
    subroutine Modify_SurfaceBoundaryFluxes

        !Local-----------------------------------------------------------------
        type (T_Property),  pointer             :: Property
        integer                                 :: WKUB
        integer                                 :: i, j
        real                                    :: OldMass, NewMass
        real                                    :: WaterColumnVolume

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "Modify_SurfaceBoundaryFluxes")

        Property => Me%FirstProperty

        do while(associated(Property))
            
            if(Property%Evolution%SurfaceFluxes)then

                if(Property%Particulate)then
                
                    if (Me%ExternalVar%Now .ge. Me%NextCompute)then
                    
                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                            OldMass = 0.
                            NewMass = 0.

                            WKUB = Me%ExternalVar%KTop(i,j)

                            if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint ) then

                                if(Me%ExternalVar%KTopState(i, j) == 1)then
                                    
                                    OldMass = Property%Concentration(i, j, WKUB-1)      * &
                                              Me%ExternalVar%MinLayerThickness          * &
                                              Me%ExternalVar%GridCellArea(i,j)          * &
                                              Me%Sediment_DryDensity%Field(i,j,WKUB)    * &
                                              (1.- Me%ExternalVar%Porosity(i,j,WKUB))   * &
                                              Property%ISCoefficient

                                    NewMass = Property%BoundaryFlux(i,j)                * &
                                              Me%ExternalVar%GridCellArea(i,j)          * &
                                              Me%ExternalVar%DT

                                    Property%Concentration(i, j, WKUB) =                  &
                                                                    (OldMass + NewMass) / &
                                             Me%ExternalVar%DrySedimentVolume(i,j,WKUB) / &
                                             Me%Sediment_DryDensity%Field(i,j,WKUB)     / &
                                             Property%ISCoefficient

                                elseif(Me%ExternalVar%KTopState(i, j) == -1)then

                                    OldMass = Property%Concentration(i, j, WKUB+1)          * &
                                              Me%ExternalVar%MinLayerThickness              * &
                                              Me%ExternalVar%GridCellArea(i,j)              * &
                                              Me%Sediment_DryDensity%Field(i,j,WKUB+1)      * &
                                              (1.- Me%ExternalVar%Porosity(i,j,WKUB+1))     * &
                                              Property%ISCoefficient

                                    NewMass = Property%Concentration(i, j, WKUB)            * &
                                              Me%ExternalVar%DrySedimentVolumeOld(i,j,WKUB) * &
                                              Me%Sediment_DryDensity%Field(i,j,WKUB)        * &
                                              Property%ISCoefficient

                                    Property%Concentration(i, j, WKUB) =                      &
                                                                    (OldMass + NewMass)     / &
                                             Me%ExternalVar%DrySedimentVolume(i,j,WKUB)     / &
                                             Me%Sediment_DryDensity%Field(i,j,WKUB)         / &
                                             Property%ISCoefficient

                                    Property%Concentration(i, j, WKUB+1) = 0.

                                else

                                    OldMass = Property%Concentration(i, j, WKUB)            * &
                                              Me%ExternalVar%DrySedimentVolumeOld(i,j,WKUB) * &
                                              Me%Sediment_DryDensity%Field(i,j,WKUB)        * &
                                              Property%ISCoefficient      
                                    
                                    NewMass = Property%BoundaryFlux(i,j)                    * &
                                              Me%ExternalVar%GridCellArea(i,j)              * &
                                              Me%ExternalVar%DT

                                    
                                    Property%Concentration(i, j, WKUB) = (OldMass + NewMass)/ &
                                             Me%ExternalVar%DrySedimentVolume(i,j,WKUB)     / &
                                             Me%Sediment_DryDensity%Field(i,j,WKUB)         / &
                                             Property%ISCoefficient
                                    
                               end if

                            endif

                        enddo
                        enddo

                    end if

                else

                    if (Me%ExternalVar%Now .ge. Me%NextCompute)then

                        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                            
                            OldMass = 0.
                            NewMass = 0.

                            WKUB = Me%ExternalVar%KTop(i,j)

                            if (Me%ExternalVar%WaterPoints3D(i,j,WKUB) == WaterPoint ) then

                                if(Me%ExternalVar%KTopState(i, j) == 1)then
                                    
                                    !if layer is new then concentration is equal to the concentration 
                                    !layer below 
                                    Property%Concentration(i,j, WKUB) = &
                                            Property%Concentration(i,j, WKUB-1)
                                    

                                    WaterColumnVolume = -1 * Me%ExternalVar%WaterFlux(i,j)  * &
                                                        Me%ExternalVar%DT


                                    !review this - approximation is made for water volume for 
                                    !newly formed layer
                                    OldMass  = Property%Concentration(i,j, WKUB)            * &
                                               Property%ISCoefficient                       * &
                                               (Me%ExternalVar%WaterVolume(i,j,WKUB) - WaterColumnVolume)
                                else

                                    OldMass  = Property%Concentration(i,j, WKUB)            * &
                                               Me%ExternalVar%WaterVolumeOld(i,j,WKUB)      * &
                                               Property%ISCoefficient
                                end if

                                NewMass = Property%BoundaryFlux(i,j)                        * &
                                          Me%ExternalVar%GridCellArea(i,j)                  * &
                                          Me%ExternalVar%DT

                                if(Me%ExternalVar%KTopState(i, j) == -1)then

                                    !when erosion occurs remove top layer
                                    NewMass = NewMass + Property%Concentration(i,j,WKUB+1)  * &
                                              Me%ExternalVar%WaterVolumeOld(i,j,WKUB+1)     * &
                                              Property%ISCoefficient

                                    Property%Concentration(i,j,WKUB+1) = 0.

                                end if

                                Property%Concentration(i, j, WKUB) = (OldMass + NewMass)    / &
                                                        Me%ExternalVar%WaterVolume(i,j,WKUB)/ &
                                                        Property%ISCoefficient
                                
                                if(Property%Concentration(i, j, WKUB)<0.)then
                                    Property%Concentration(i, j, WKUB) = 0.
                                    write(*,*) 'negative dissolved concentration surface boundary fluxes'
                                    write(*,*)i,j,WKUB
                                endif 
                                
                            endif
                        enddo
                        enddo

                    end if

                end if

            end if

            Property => Property%Next

        end do
          
        nullify (Property)

        Me%NextCompute = Me%NextCompute + Me%ExternalVar%DT

        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "Modify_SurfaceBoundaryFluxes")

    end subroutine Modify_SurfaceBoundaryFluxes
 
    !--------------------------------------------------------------------------
    
    subroutine Partition_Processes


        !Local-----------------------------------------------------------------
        type (T_Property),  pointer             :: Property
        type (T_Property),  pointer             :: ParticulateProperty
        type (T_Property),  pointer             :: DissolvedProperty
        real                                    :: DissolvedFraction
        real                                    :: ParticulateFraction
        real                                    :: dt, Rate
        integer                                 :: STAT_CALL, PartPropNumber
        integer                                 :: i, j, k
        real, dimension(:), pointer             :: SystemResult

        !Begin-----------------------------------------------------------------


        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "Partition_Processes")

        Property => Me%FirstProperty

        do while(associated(Property))
            
            if(.not. Property%Particulate .and. Property%Evolution%Partitioning)then
                
                DissolvedProperty => Property

                if (Me%ExternalVar%Now .ge. DissolvedProperty%Evolution%NextCompute)then

                    PartPropNumber = DissolvedProperty%Evolution%Partition%Couple_ID

                    call Search_Property(ParticulateProperty, PropertyXID=PartPropNumber,STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                           &
                        stop 'Partition_Processes - ModuleSedimentProperties - ERR01'

                    !converts property concentration in whole cell
                    
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Me%ExternalVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == WaterPoint ) then

                            do k = Me%WorkSize%KLB, Me%ExternalVar%Ktop(i,j)

                                !kg of property per m3 of interstitial water
                                !                 [kg/m3cell]      =  [kg/m3water]/[m3cell]*[m3water]
                                DissolvedProperty%Concentration(i,j,k)=                               &
                                                            DissolvedProperty%Concentration(i,j,k)  / & 
                                                            Me%ExternalVar%VolumeZ(i,j,k)           * &
                                                            Me%ExternalVar%WaterVolume(i,j,k)       * &
                                                            DissolvedProperty%ISCoefficient

                                !kg of property per kg of dry sediment
                                !              [kg/m3cell]     = [kg/kgsed]/[m3cell]*[m3sed]*[kgsed/m3sed]
                                ParticulateProperty%Concentration(i,j,k) = &
                                                            ParticulateProperty%Concentration(i,j,k)/ &
                                                            Me%ExternalVar%VolumeZ(i,j,k)           * &
                                                            Me%ExternalVar%DrySedimentVolume(i,j,k) * &
                                                            Me%Sediment_DryDensity%Field(i,j,k)     * &
                                                            ParticulateProperty%ISCoefficient
                            enddo
                        endif
                    enddo
                    enddo
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Me%ExternalVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == WaterPoint)then

                            do k = Me%WorkSize%KLB, Me%ExternalVar%KTop(i,j)


                                DissolvedFraction   = DissolvedProperty%Evolution%Partition%Fraction%Field(i,j,k)

                                ParticulateFraction = ParticulateProperty%Evolution%Partition%Fraction%Field(i,j,k)

                                dt                  = DissolvedProperty%Evolution%DTInterval

                                Rate                = DissolvedProperty%Evolution%Partition%Rate%Field(i,j,k)

                                !     D(t+dt) - D(t)
                                !   ----------------- = k (rd * P - rp * D)
                                !           dt
                                !           
                                !     P(t+dt) - P(t)
                                !   ----------------- = k (rp * D - rd * P)
                                !           dt
                                !
                                !     D(t+dt) [ 1 + rp * k * dt ] + P(t+dt) [  - rd * k * dt]     = D(t)
                                !
                                !     D(t+dt) [   - rp * k * dt ] + P(t+dt) [1 + rd * k * dt]     = P(t)
                                !     __                                       __   _       _     _    _
                                !    |                                           | |         | = |      |
                                !    | 1 + rp * k * dt            - rd * k * dt  | | D(t+dt) |   | D(t) |
                                !    |                                           | |         |   |      |
                                !    |   - rp * k * dt          1 + rd * k * dt  | | P(t+dt) |   | P(t) |
                                !    |__                                       __| |_       _|   |_    _|
                                !

                                !Fill matrix
                                Me%PartitionMatrix (1, 1) = 1. + ParticulateFraction * Rate * dt
                                Me%PartitionMatrix (1, 2) =    - DissolvedFraction   * Rate * dt
                                Me%PartitionMatrix (2, 1) =    - ParticulateFraction * Rate * dt
                                Me%PartitionMatrix (2, 2) = 1. + DissolvedFraction   * Rate * dt

                                Me%IndependentTerm(1)     = DissolvedProperty%Concentration(i,j,k)
                                Me%IndependentTerm(2)     = ParticulateProperty%Concentration(i,j,k)

                                call LUD(LUD_ID         = Me%ObjLUD,             &
                                         Coefficients   = Me%PartitionMatrix,    &
                                         IndTerm        = Me%IndependentTerm,    &
                                         x              = SystemResult,          &
                                         STAT           = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_)                       &
                                    stop 'Partition_Processes - ModuleSedimentProperties - ERR02'

                                DissolvedProperty%Concentration(i,j,k)      = SystemResult(1)
                                ParticulateProperty%Concentration(i,j,k)    = SystemResult(2)
                            enddo


                        endif
                    enddo
                    enddo
                    !converts property concentration into standard units
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                        if (Me%ExternalVar%WaterPoints3D(i,j,Me%WorkSize%KUB) == WaterPoint)then

                            do k = Me%WorkSize%KLB, Me%ExternalVar%KTop(i,j)

                                !kg/m3water = kg/m3cell * m3cell / m3water
                                DissolvedProperty%Concentration(i,j,k) =                                  &
                                                             DissolvedProperty%Concentration(i, j, k)   * &
                                                             Me%ExternalVar%VolumeZ(i,j,k)              / &
                                                             Me%ExternalVar%WaterVolume(i,j,k)          / &
                                                             DissolvedProperty%ISCoefficient

                                !kg/kgsed = kg/m3cell * m3cell / m3sed / (kgsed/m3sed)
                                ParticulateProperty%Concentration(i,j,k) =                                &
                                                             ParticulateProperty%Concentration(i,j,k)   * &
                                                             Me%ExternalVar%VolumeZ(i,j,k)              / &
                                                             Me%ExternalVar%DrySedimentVolume(i,j,k)    / &
                                                             Me%Sediment_DryDensity%Field(i,j,k)        / &
                                                             ParticulateProperty%ISCoefficient
                            enddo
                        endif
                    enddo
                    enddo

                end if

                nullify (DissolvedProperty, ParticulateProperty)
            
            end if

            Property => Property%Next

        end do
          
        nullify (Property)

        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "Partition_Processes")

    end subroutine Partition_Processes

    
    !--------------------------------------------------------------------------


    subroutine SetMinimumConcentration
       
        !Local--------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: i, j, k

        !Begin----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "SetMinimumConcentration")

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))
            
            if (PropertyX%Evolution%MinConcentration) then
            
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then
    
                        if (PropertyX%Concentration(i, j, k) < PropertyX%MinValue) then

                            PropertyX%Mass_Created(i, j, k)     = &
                                        PropertyX%Mass_Created(i, j, k)                         + &
                                        (PropertyX%MinValue - PropertyX%Concentration(i, j, k)) * &
                                        Me%ExternalVar%VolumeZ (i, j, k)

                            PropertyX%Concentration(i, j, k)    = PropertyX%MinValue
                        
                        endif
               
                    endif

                enddo
                enddo
                enddo

            end if

            PropertyX => PropertyX%Next
        end do

        nullify(PropertyX)
        
        
        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "SetMinimumConcentration")


    end subroutine SetMinimumConcentration

    
    !--------------------------------------------------------------------------

    
    subroutine ComputeDiffusivity (PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, DiffusionMethod
        real,    dimension(:,:,:), pointer          :: VerticalDiffCoef
        real,    dimension(:,:,:), pointer          :: HorizontalDiffCoef
        real                                        :: BioturbationDiffCoef
        real                                        :: TurbulentVertDiffCoef
        real                                        :: TurbulentHorzDiffCoef
        real                                        :: MolecularDiffCoef
        real                                        :: CorrMolecularDiffCoef

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "ComputeDiffusivity")

        VerticalDiffCoef     => Me%VerticalDiffCoef%Field
        HorizontalDiffCoef   => Me%HorizontalDiffCoef%Field
        MolecularDiffCoef    =  PropertyX%Evolution%Advec_Difus_Parameters%Molecular_Diff_Coef
        DiffusionMethod      =  PropertyX%Evolution%Advec_Difus_Parameters%Diffusion_Method

        do k=Me%WorkSize%KLB, Me%WorkSize%KUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then

                TurbulentVertDiffCoef = VerticalDiffCoef(i,j,k)                     * &
                                        abs(Me%ExternalVar%Velocity_W(i,j,k)        * &
                                        Me%ExternalVar%ComputeFacesW3D(I,J,K))      + &
                                        HorizontalDiffCoef(i,j,k)                   * &
                                        abs((Me%ExternalVar%Velocity_V(i,j,k)       * &
                                        Me%ExternalVar%ComputeFacesV3D(I,J,K))      + &
                                        (Me%ExternalVar%Velocity_U(i,j,k)           * &
                                        Me%ExternalVar%ComputeFacesU3D(I,J,K))/2) 

                TurbulentHorzDiffCoef =  VerticalDiffCoef(i,j,k)                    * &
                                         abs((Me%ExternalVar%Velocity_V(i,j,k)      * &
                                         Me%ExternalVar%ComputeFacesV3D(I,J,K))     + &
                                         (Me%ExternalVar%Velocity_U(i,j,k)          * &
                                         Me%ExternalVar%ComputeFacesU3D(I,J,K))/2)  + &
                                         HorizontalDiffCoef(i,j,k)                  * &
                                         abs(Me%ExternalVar%Velocity_W(i,j,k)       * &
                                         Me%ExternalVar%ComputeFacesW3D(I,J,K)) 
                
                if(DiffusionMethod == Berner_1980)then
                    
                    !Berner, 1980
                    CorrMolecularDiffCoef = MolecularDiffCoef / (Me%ExternalVar%Tortuosity(i,j,k)**2)
                
                elseif(DiffusionMethod == Soetaert_1996)then
                    
                    !Soetaert, 1996
                    CorrMolecularDiffCoef = MolecularDiffCoef * (Me%ExternalVar%Porosity(i,j,k)**2)

                end if



                if(PropertyX%Evolution%ComputeBioturbation)then
                    BioturbationDiffCoef = PropertyX%Evolution%Bioturbation%Coef(i,j,k)
                else
                    BioturbationDiffCoef = 0.
                end if

                PropertyX%HorizontalDiffusivity(i,j,k) = CorrMolecularDiffCoef + TurbulentHorzDiffCoef
                                                          

                PropertyX%VerticalDiffusivity  (i,j,k) = CorrMolecularDiffCoef + TurbulentVertDiffCoef + &
                                                         BioturbationDiffCoef

            endif

        enddo
        enddo
        enddo

        nullify(VerticalDiffCoef, HorizontalDiffCoef)

        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "ComputeDiffusivity")

    end subroutine ComputeDiffusivity


    !--------------------------------------------------------------------------

    
    subroutine Compute_Bioturbation (PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
        real, dimension(:,:,:), pointer             :: BioturbationDiffCoef
        real, dimension(:,:,:), pointer             :: Depth
        real                                        :: BioDepth
        real                                        :: DecayCoef
        real                                        :: DefaultCoef

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "Compute_Bioturbation")

        Depth                => Me%ExternalVar%Depth    !review this
        BioturbationDiffCoef => PropertyX%Evolution%Bioturbation%Coef
        BioDepth             =  PropertyX%Evolution%Bioturbation%BioDepth
        DefaultCoef          =  PropertyX%Evolution%Bioturbation%DefaultCoef
        DecayCoef            =  PropertyX%Evolution%Bioturbation%DecayCoef

        do k=Me%WorkSize%KLB, Me%WorkSize%KUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then

                if(Depth(i,j,k) .le. BioDepth)then

                    BioturbationDiffCoef(i,j,k) = DefaultCoef

                else

                    BioturbationDiffCoef(i,j,k) = DefaultCoef * exp(-(Depth(i,j,k)-BioDepth)/DecayCoef)

                endif

                if(PropertyX%Particulate)then

                    PropertyX%HorizontalDiffusivity(i,j,k) = 0.

                    PropertyX%VerticalDiffusivity  (i,j,k) = BioturbationDiffCoef(i,j,k)

                end if

            endif

        enddo
        enddo
        enddo

        nullify(BioturbationDiffCoef, Depth)


        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "Compute_Bioturbation")


    end subroutine Compute_Bioturbation

    !--------------------------------------------------------------------------

    ! This subroutine is responsable for defining       
    ! the next time to actualize the value of each      
    ! property                                          

    subroutine Actualize_Time_Evolution

        !Local-----------------------------------------------------------------
        type (T_Property), pointer          :: Property
        type (T_Time    )                   :: Actual

        !----------------------------------------------------------------------

        Actual = Me%ExternalVar%Now

        Property => Me%FirstProperty

        do while(associated(Property))

            if (Property%Evolution%Variable) then
                if (Actual.GE.Property%Evolution%NextCompute) then
                    Property%Evolution%LastCompute = Property%Evolution%NextCompute
                    Property%Evolution%NextCompute = Property%Evolution%NextCompute &
                                                   + Property%Evolution%DTInterval
                end if
            end if

            Property => Property%Next

        end do

    end subroutine Actualize_Time_Evolution

    !--------------------------------------------------------------------------

    subroutine TimeStepActualization

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type (T_Property), pointer          :: PropertyX
        logical                             :: VariableDT
        real                                :: NewDT

        !Begin-----------------------------------------------------------------

        call GetVariableDT (Me%ObjTime, VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'TimeStepActualization - ModuleSedimentProperties - ERR01'


        call GetComputeTimeStep(Me%ObjTime, NewDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'TimeStepActualization - ModuleSedimentProperties - ERR02'


        if (VariableDT) then

            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))

                PropertyX%Evolution%DTInterval = NewDT
                
                PropertyX => PropertyX%Next
            end do

        endif

    end subroutine TimeStepActualization
    
    !--------------------------------------------------------------------------

    subroutine OutPut_Results_HDF

        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
        real                               :: Year, Month, Day, Hour, Minute, Second

        !Local-----------------------------------------------------------------
        type (T_Property), pointer         :: PropertyX
        logical                            :: FirstTime
        integer                            :: OutPutNumber
        integer, dimension(6)              :: TimeAux
        real,    dimension(6), target      :: AuxTime
        real,    dimension(:),  pointer    :: TimePtr
        integer                            :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                            :: WorkKLB, WorkKUB

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "OutPut_Results_HDF")


        ILB = Me%Size%ILB;  WorkILB = Me%WorkSize%ILB
        IUB = Me%Size%IUB;  WorkIUB = Me%WorkSize%IUB
        JLB = Me%Size%JLB;  WorkJLB = Me%WorkSize%JLB
        JUB = Me%Size%JUB;  WorkJUB = Me%WorkSize%JUB
        KLB = Me%Size%KLB;  WorkKLB = Me%WorkSize%KLB
        KUB = Me%Size%KUB;  WorkKUB = Me%WorkSize%KUB

        FirstTime = .true.

        OutPutNumber = Me%OutPut%NextOutPut  

        if (Me%ExternalVar%Now >= Me%OutPut%OutTime(OutPutNumber)) then 

            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 

                if (PropertyX%OutPutHDF) then

                    call ExtractDate(Me%ExternalVar%Now, &
                                     Year = Year, Month  = Month,  Day    = Day, &
                                     Hour = Hour, Minute = Minute, Second = Second)

                    TimeAux(1) = int(Year  )
                    TimeAux(2) = int(Month )
                    TimeAux(3) = int(Day   )
                    TimeAux(4) = int(Hour  )
                    TimeAux(5) = int(Minute)
                    TimeAux(6) = int(Second)

                    if (FirstTime) then 

                        !Writes current time
                        call ExtractDate   (Me%ExternalVar%Now,                         &
                                            AuxTime(1), AuxTime(2), AuxTime(3),         &
                                            AuxTime(4), AuxTime(5), AuxTime(6))
                        TimePtr => AuxTime
                        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR01'

                        call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                             "Time", "YYYY/MM/DD HH:MM:SS",             &
                                             Array1D = TimePtr,                         &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR02'

                        !Writes SZZ
                        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,     &
                                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR03'


                        call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical", &
                                             "m", Array3D = Me%ExternalVar%SZZ,         &
                                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR04'

                        !Writes OpenPoints
                        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,              &
                                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,        &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR05'

                        call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",&
                                             "-", Array3D = Me%ExternalVar%OpenPoints3D,&
                                             OutputNumber = OutPutNumber,  STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR06'

                        FirstTime = .false.

                    endif

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//PropertyX%ID%Name,                  &
                                       PropertyX%ID%Name,                               &
                                       PropertyX%ID%Units,                              &
                                       Array3D      = PropertyX%Concentration,          &
                                       OutputNumber = OutPutNumber,                     &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR07'

                    if (PropertyX%Evolution%MinConcentration)then

                        call HDF5WriteData  (Me%ObjHDF5,                                &
                                             "/Mass_created/"//PropertyX%ID%Name,       &
                                             PropertyX%ID%Name,                         &
                                             PropertyX%ID%Units,                        &
                                             Array3D = PropertyX%Mass_created,          &
                                             OutputNumber = OutPutNumber,               &
                                             STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                      &
                            stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR08'

                    endif
                    
                                        
                     if (PropertyX%ID%IDNumber==SeagrassesRoots_) then
                
                    call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                   &
                                         WorkJLB, WorkJUB, WorkKLB, WorkKUB,             &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR09'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/"//"seagrasses roots biomass",  &
                                         "seagrasses roots biomass", "gdw/m2",                  &
                                         Array2D = Me%SeagrassesRoots%Biomass, &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR10'
                endif

                    !Writes everything to disk
                    call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'OutPut_Results_HDF - ModuleSedimentProperties - ERR09'


                endif
      
                PropertyX => PropertyX%Next

            enddo

            Me%OutPut%NextOutPut = OutPutNumber + 1
        
        endif

        
        
        nullify(PropertyX)

        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "OutPut_Results_HDF")

    end subroutine OutPut_Results_HDF

    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "OutPut_TimeSeries")


        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data3D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModuleSedimentProperties - ERR01'

            endif
            PropertyX=>PropertyX%Next
        enddo
        
        
        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "OutPut_TimeSeries")

    end subroutine OutPut_TimeSeries
    
    !--------------------------------------------------------------------------

    subroutine OutPut_BoxTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSedimentProperties", "OutPut_BoxTimeSeries")

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        KLB = Me%Size%KLB 
        KUB = Me%Size%KUB 

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%BoxTimeSerie) then

                Me%CellMass(:,:,:) = 0.

                if(PropertyX%Particulate)then

                    do K = KLB, KUB
                    do J = JLB, JUB
                    do I = ILB, IUB
                        
                        Me%CellMass(i,j,k) = PropertyX%Concentration(i, j, k)          * &
                                             Me%ExternalVar%DrySedimentVolume(i, j, k) * &
                                             Me%Sediment_DryDensity%Field(i,j,k)
                    end do
                    end do
                    end do


                else

                    do K = KLB, KUB
                    do J = JLB, JUB
                    do I = ILB, IUB
                        
                        Me%CellMass(i,j,k) = PropertyX%Concentration(i, j, k) * &
                                             Me%ExternalVar%WaterVolume(i, j, k)
                    end do
                    end do
                    end do

                end if

                
                call BoxDif(Me%ObjBoxDif, Me%CellMass,                      &
                            "Sed_"//trim(adjustl(PropertyX%ID%name)),       &
                            Me%ExternalVar%WaterPoints3D,                   &
                            STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)  &
                    stop 'OutPut_BoxTimeSeries - ModuleSedimentProperties - ERR01'

                Me%CellMass(:,:,:) = null_real

            endif

            PropertyX=>PropertyX%Next

        enddo

        if (MonitorPerformance) call StopWatch ("ModuleSedimentProperties", "OutPut_BoxTimeSeries")


    end subroutine OutPut_BoxTimeSeries
    
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine UngetSedimentProperties3D(SedimentPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: SedimentPropertiesID  

        real, pointer, dimension(:,:,:) :: Array

        integer, optional, intent(OUT)  :: STAT
     

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_    
                
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mSEDIMENTPROPERTIES_, Me%InstanceID, "UngetSedimentProperties3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
            

        !----------------------------------------------------------------------

    end subroutine UngetSedimentProperties3D
    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine UngetSedimentProperties2Dreal8(SedimentPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                          :: SedimentPropertiesID  

        real(8), pointer, dimension(:,:) :: Array

        integer, optional, intent(OUT)   :: STAT
     

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mSEDIMENTPROPERTIES_, Me%InstanceID, "UngetSedimentProperties2Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetSedimentProperties2Dreal8

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine UngetSedimentProperties2Dreal4(SedimentPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                          :: SedimentPropertiesID  

        real(4), pointer, dimension(:,:) :: Array

        integer, optional, intent(OUT)   :: STAT
     

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mSEDIMENTPROPERTIES_, Me%InstanceID, "UngetSedimentProperties2Dreal4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetSedimentProperties2Dreal4

    !--------------------------------------------------------------------------
    

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
    
    
    subroutine GetSedimentPropertyOptions(SedimentPropertiesID, PropIDNumber, DTInterval, &
                                          SurfaceFluxes, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        integer,           intent(IN )              :: PropIDNumber
        real,              intent(OUT)              :: DTInterval
        logical,           intent(OUT)              :: SurfaceFluxes
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_, STAT_CALL
        type(T_Property), pointer                   :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropIDNumber, STAT = STAT_CALL)

            if (STAT_CALL == SUCCESS_) then
                
                DTInterval    = PropertyX%Evolution%DTInterval

                SurfaceFluxes = PropertyX%Evolution%SurfaceFluxes

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else cd1
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
            

        !----------------------------------------------------------------------
     
    end subroutine GetSedimentPropertyOptions
   
   
    !----------------------------------------------------------------------

    
    logical function SedimentPropertyExists(SedimentPropertiesID, PropertyXID, STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        integer,                    intent (IN)     :: PropertyXID
        integer         , optional, intent (OUT)    :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_ 
        integer                                     :: STAT_CALL
        type(T_Property), pointer                   :: PropertyX

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            SedimentPropertyExists = .false.

            call Search_Property(PropertyX, PropertyXID = PropertyXID, STAT = STAT_CALL)

            if (STAT_CALL == SUCCESS_) SedimentPropertyExists = .true.



            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
            

        !----------------------------------------------------------------------
     
    end function SedimentPropertyExists

    !----------------------------------------------------------------------


    subroutine GetSedimentConcentration(SedimentPropertiesID, ConcentrationX, PropertyXIDNumber, &
                                PropertyXUnits, PropertyXISCoef, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        real, pointer, dimension(:,:,:)             :: ConcentrationX
        character(LEN = *), optional, intent(OUT)   :: PropertyXUnits
        integer,                      intent(IN )   :: PropertyXIDNumber
        real,               optional, intent(OUT)   :: PropertyXISCoef
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_CALL              
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: UnitsSize
        integer                                     :: STAT_    

        !------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call Read_Lock(mSEDIMENTPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)

            if (STAT_CALL == SUCCESS_) then
                
                ConcentrationX => PropertyX%Concentration

                if (present(PropertyXUnits)) then 
                   UnitsSize      = LEN (PropertyXUnits)
                   PropertyXUnits = PropertyX%ID%Units(1:UnitsSize)
                end if

                if (present(PropertyXISCoef)) then 
                   PropertyXISCoef = PropertyX%IScoefficient
                end if

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetSedimentConcentration


    !--------------------------------------------------------------------------
    
       subroutine GetRootsArray2D(SedimentPropertiesID, Array, ArrayID, STAT)  
    !

        !Arguments---------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        real, pointer, dimension(:,:)               :: Array
        integer,            optional, intent(OUT)   :: STAT
        integer,            optional, intent(IN)    :: ArrayID

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mSEDIMENTPROPERTIES_, Me%InstanceID) 
           
           select case(ArrayID)
           
           case(SeagrassesRoots_)

            Array => Me%SeagrassesRoots%Biomass
           
          case(NintFactorR_)

            Array => Me%SeagrassesRoots%NintFactor2DR
            
          case(RootsMort_)

            Array => Me%SeagrassesRoots%RootsMort2DR

          
          end select

            STAT_ = SUCCESS_
            
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetRootsArray2D
    !--------------------------------------------------------------------------
    
   

    subroutine GetSeagrassesRootsRates(SedimentPropertiesID, RateID, Rateflux, STAT)  
   
        !Arguments---------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        real, pointer, dimension(:,:,:)             :: Rateflux
        integer,                      intent(IN)    :: RateID
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mSEDIMENTPROPERTIES_, Me%InstanceID) 
            
            select case(RateID)

            
            case(RootsUptakeN_)
           
            Rateflux => Me%SeagrassesRoots%UptakeNH4s3D ! gN/day  
         
           case(RootsUptakeP_)
           
            Rateflux => Me%SeagrassesRoots%UptakePO4s3D ! gN/day  
           
            end select

            STAT_ = SUCCESS_
            
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetSeagrassesRootsRates
     !----------------------------------------------------------------------------

    subroutine GetSedimentDryDensity(SedimentPropertiesID, SedimentDryDensity,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        real, dimension(:,:,:),          pointer    :: SedimentDryDensity
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mSEDIMENTPROPERTIES_, Me%InstanceID)

            SedimentDryDensity => Me%Sediment_DryDensity%Field


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_


    end subroutine GetSedimentDryDensity

    !--------------------------------------------------------------------------

    subroutine SetFluxToSedimentProperties(SedimentPropertiesID, PropertyID, Flux, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        integer,                      intent(IN)    :: PropertyID
        real, pointer, dimension(:,:)               :: Flux
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_CALL
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: STAT_    

        !------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then
            
            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyID, STAT = STAT_CALL)
                                 
            if (STAT_CALL == SUCCESS_) then
                
                PropertyX%BoundaryFlux => Flux

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetFluxToSedimentProperties

    !--------------------------------------------------------------------------

    subroutine SetSedimentWaterFlux(SedimentPropertiesID, WaterFlux, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        real(8), pointer, dimension(:,:)            :: WaterFlux
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then
                
            Me%ExternalVar%WaterFlux => WaterFlux

            STAT_ = SUCCESS_

        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_


    end subroutine SetSedimentWaterFlux


    !--------------------------------------------------------------------------
    
    subroutine SetSubModulesConstructor
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        call SetSedimentDryDensity(Me%ObjConsolidation, Me%Sediment_DryDensity%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'SetSubModulesConstructor - ModuleSedimentProperties - ERR01'
    
    end subroutine SetSubModulesConstructor
    
    
    !--------------------------------------------------------------------------

    subroutine ReadLockExternalVar
        
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR01'

        !WaterPoints3D
        call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR02'

        !OpenPoints3D
        call GetOpenPoints3D(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR03'

        !XX_IE and YY_IE
        call GetHorizontalGrid (Me%ObjHorizontalGrid,                                   &
                                XX_IE = Me%ExternalVar%XX_IE,                           &
                                YY_IE = Me%ExternalVar%YY_IE,                           &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR04'

        !BoundaryPoints2D
        call GetBoundaries(Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR05'

        !WaterPercentage
        call GetConsolidationWaterPercentage (Me%ObjConsolidation, Me%ExternalVar%WaterPercentage, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR06'

        !DrySedimentVolume      
        call GetConsolidationDrySedVolume(Me%ObjConsolidation,                                          &
                                          DrySedimentVolume     = Me%ExternalVar%DrySedimentVolume,     &
                                          DrySedimentVolumeOld  = Me%ExternalVar%DrySedimentVolumeOld,  &
                                          STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR07'
                    
        !WaterVolume
        call GetConsolidationWaterVolume(Me%ObjConsolidation,                           &
                                         WaterVolume    = Me%ExternalVar%WaterVolume,   &
                                         WaterVolumeOld = Me%ExternalVar%WaterVolumeOld,&
                                         STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR08'

        !ComputeFaces3D
        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesU3D = Me%ExternalVar%ComputeFacesU3D,        &
                               ComputeFacesV3D = Me%ExternalVar%ComputeFacesV3D,        &
                               ComputeFacesW3D = Me%ExternalVar%ComputeFacesW3D,        &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR09'

        !LandPoints3D
        call GetLandPoints3D(Me%ObjMap, Me%ExternalVar%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR11'

        !WaterFluxes
        call GetConsolidationWaterFluxes(Me%ObjConsolidation,                           & 
                                         WaterFluxX = Me%ExternalVar%WaterFluxX,        & 
                                         WaterFluxY = Me%ExternalVar%WaterFluxY,        & 
                                         WaterFluxZ = Me%ExternalVar%WaterFluxZ,        &
                                         STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR12'

        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR13'

        call GetGeometryVolumes(Me%ObjGeometry, VolumeZ = Me%ExternalVar%VolumeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR14'
                                    
        call GetConsolidationTortuosity(Me%ObjConsolidation, Me%ExternalVar%Tortuosity, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR15'


        call GetConsolidationVelocityZ(Me%ObjConsolidation,                             &
                                       Velocity_W = Me%ExternalVar%Velocity_W,          &
                                       STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR16'

        call GetConsolidationVelocityXY(Me%ObjConsolidation,                            &
                                        Velocity_U = Me%ExternalVar%Velocity_U,         &
                                        Velocity_V = Me%ExternalVar%Velocity_V,         &
                                        STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR17'


        call GetGeometryDistances (Me%ObjGeometry, DWZ = Me%ExternalVar%DWZ,            &
                                   DZZ = Me%ExternalVar%DZZ, SZZ = Me%ExternalVar%SZZ,  &
                                   STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR18'
        
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR19'

        
        call GetGridData(Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT_CALL)     
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR20'

        call GetConsolidationDepth(Me%ObjConsolidation,                                 &
                                   CellCenterDepth = Me%ExternalVar%Depth,              &
                                   STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR21'
        
        call GetConsolidationPorosity(Me%ObjConsolidation,                              &
                                      Porosity = Me%ExternalVar%Porosity,               &
                                      STAT     =  STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR22'

        call GetGeometryKTop(Me%ObjGeometry,                                            &
                             KTopZ  = Me%ExternalVar%KTop,                              &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR23'
        
        
        call GetConsolidationKTopState(Me%ObjConsolidation, Me%ExternalVar%KTopState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSedimentProperties - ERR24'

    end subroutine ReadLockExternalVar


    !--------------------------------------------------------------------------


    subroutine ReadUnlockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------


        call UnGetMap(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR01'
        
        call UnGetMap(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR02'

        call UnGetMap(Me%ObjMap, Me%ExternalVar%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR03'

        call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR04'
        
        call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR05'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR06'
        
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR07'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%WaterPercentage, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR08'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%WaterVolume, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR09'
        
        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%WaterVolumeOld, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR10'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%DrySedimentVolume, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR11'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%DrySedimentVolumeOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR12'
 
        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%WaterFluxX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR13'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%WaterFluxY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR14'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%WaterFluxZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR15'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR16'
        
        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%VolumeZ, STAT = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR18'

        call UngetConsolidation (Me%ObjConsolidation, Me%ExternalVar%Tortuosity, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR19'


        call UngetConsolidation (Me%ObjConsolidation, Me%ExternalVar%Velocity_W, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR20'

        call UngetConsolidation (Me%ObjConsolidation, Me%ExternalVar%Velocity_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR21'

        call UngetConsolidation (Me%ObjConsolidation, Me%ExternalVar%Velocity_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR22'

        call UngetHorizontalMap (Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR23'
                               
        call UnGetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesU3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR24'

        call UnGetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesV3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR25'

        call UnGetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesW3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR26'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR27'

        call UnGetGridData(Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT_CALL)     
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR28'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%Depth, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR29'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%Porosity, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR30'
       
        call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%DZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR31'

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%KTop, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR32'

        call UngetConsolidation(Me%ObjConsolidation, Me%ExternalVar%KTopState, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleSedimentProperties - ERR33'

        call null_time(Me%ExternalVar%Now)

    end subroutine ReadUnlockExternalVar


    !--------------------------------------------------------------------------




    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillSedimentProperties(SedimentPropertiesID, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SedimentPropertiesID
        type (T_Time)                       :: EndTime, BeginTime
        integer, optional, intent(OUT)      :: STAT

        integer :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                             :: STAT_,ready_, nUsers
        type (T_Property),  pointer         :: PropertyX

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(SedimentPropertiesID, ready_)    

cd1:    if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSEDIMENTPROPERTIES_,  Me%InstanceID)

            if (nUsers == 0) then

                call GetComputeTimeLimits(Me%ObjTime, EndTime = EndTime, &
                                          BeginTime = BeginTime, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillSedimentProperties - ModuleSedimentProperties - ERR00'
                               
                ! Actualized the time
                Me%ExternalVar%Now = EndTime

                call ReadLockExternalVar

                !Writes the final results in HDF format
                call Write_Final_SedProperties_HDF

                call ReadUnlockExternalVar 

                if(Me%OutPut%Yes)then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR01'
                end if

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR02'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjGridData)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR03'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR04'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR05'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR06'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjMap)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR07'

                nUsers = DeassociateInstance(mCONSOLIDATION_,   Me%ObjConsolidation)
                if (nUsers == 0) stop 'KillSedimentProperties - ModuleSedimentProperties - ERR08'


                if (Me%Coupled%BoxTimeSerie%Yes) then
                    
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR09'

                    deallocate(Me%MassFluxesX, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR10'
                    nullify(Me%MassFluxesX)

                    deallocate(Me%MassFluxesY, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR11'
                    nullify(Me%MassFluxesY)

                    deallocate(Me%MassFluxesZ, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR12'
                    nullify(Me%MassFluxesZ)

                    deallocate(Me%CellMass, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR13'
                    nullify(Me%CellMass)

                    Me%Coupled%BoxTimeSerie%Yes = .false.

                end if
                    
                if (associated(Me%Sediment_DryDensity%Field))then
                    deallocate(Me%Sediment_DryDensity%Field, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR14'
                    nullify   (Me%Sediment_DryDensity%Field)
                end if

                if (associated(Me%VerticalDiffCoef%Field))then
                    deallocate(Me%VerticalDiffCoef%Field, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR15'
                    nullify   (Me%VerticalDiffCoef%Field)
                end if

                if (associated(Me%HorizontalDiffCoef%Field))then
                    deallocate(Me%HorizontalDiffCoef%Field, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR16'
                    nullify   (Me%HorizontalDiffCoef%Field)
                end if


                if(Me%Coupled%AdvectionDiffusion%Yes)then

                    deallocate(Me%DiffCoef, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR001'
        
                    deallocate(Me%Ti_Coef, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR002'
        
                    deallocate(Me%D_Coef, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR003'
        
                    deallocate(Me%E_Coef, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR004'
        
                    deallocate(Me%F_Coef, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR005'
    
                    deallocate(Me%Concentration, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR006'

                    deallocate(Me%VolZ, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR007'
        
                    deallocate(Me%VolZOld, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR008'
        
                    deallocate(Me%DZZ, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR009'

                    deallocate(Me%DWZ, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_)stop 'KillSedimentProperties - ModuleSedimentProperties - ERR010'

                end if
 
 
           if(Me%Coupled%SeagrassesRoots%Yes)then
              
                 call KillInterface(Me%ObjSeagrassSedimInteraction, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR011'
               
                 deallocate(Me%SeagrassesRoots%Biomass                    ) 
                 deallocate(Me%SeagrassesRoots%Length                     )
                 deallocate(Me%SeagrassesRoots%Occupation                 )
                 deallocate(Me%SeagrassesRoots%NintFactor3DR              )
                 deallocate(Me%SeagrassesRoots%PintFactor3DR              )
                 deallocate(Me%SeagrassesRoots%NintFactor2DR              )
                 deallocate(Me%SeagrassesRoots%UptakeNH4s3D               )
                 deallocate(Me%SeagrassesRoots%UptakePO4s3D               )
                 deallocate(Me%SeagrassesRoots%RootsMort3DR               )
                 deallocate(Me%SeagrassesRoots%RootsMort2DR               )
                 !deallocate(Me%SeagrassesRoots%Volume                     )

            
                 nullify(Me%SeagrassesRoots%Biomass                       ) 
                 nullify(Me%SeagrassesRoots%Length                        )
                 nullify(Me%SeagrassesRoots%Occupation                    )
                 nullify(Me%SeagrassesRoots%NintFactor3DR                 )
                 nullify(Me%SeagrassesRoots%PintFactor3DR                 )
                 nullify(Me%SeagrassesRoots%NintFactor2DR                 )
                 nullify(Me%SeagrassesRoots%UptakeNH4s3D                  )
                 nullify(Me%SeagrassesRoots%UptakePO4s3D                  )
                 nullify(Me%SeagrassesRoots%RootsMort3DR                  )
                 nullify(Me%SeagrassesRoots%RootsMort2DR                  )
                ! nullify(Me%SeagrassesRoots%Volume                        )
 
                ! Roots


              endif

                !Kills the TimeSerie
                if (Me%Coupled%TimeSerie%Yes) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR17'
                    Me%Coupled%TimeSerie%Yes = .false.
                endif

                if (associated(Me%OutPut%OutTime)) then
                    deallocate(Me%OutPut%OutTime, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR18'
                    nullify   (Me%OutPut%OutTime)
                end if


                Me%PropertiesNumber = FillValueInt

                !Deallocates all the water properties 
                PropertyX => Me%FirstProperty

do1 :           do while(associated(PropertyX))  

                    deallocate(PropertyX%Concentration, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                               &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR19'
                    nullify   (PropertyX%Concentration)

                    if (associated(PropertyX%HorizontalDiffusivity)) then
                        deallocate(PropertyX%HorizontalDiffusivity, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                               &
                            stop 'KillSedimentProperties - ModuleSedimentProperties - ERR20'
                        nullify   (PropertyX%HorizontalDiffusivity)
                    end if 

                    if (associated(PropertyX%VerticalDiffusivity)) then
                        deallocate(PropertyX%VerticalDiffusivity, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                               &
                            stop 'KillSedimentProperties - ModuleSedimentProperties - ERR21'
                        nullify   (PropertyX%VerticalDiffusivity)
                    end if 

                    if (associated(PropertyX%Evolution%Partition%Fraction%Field))then
                        deallocate(PropertyX%Evolution%Partition%Fraction%Field, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                               &
                            stop 'KillSedimentProperties - ModuleSedimentProperties - ERR25'
                        nullify   (PropertyX%Evolution%Partition%Fraction%Field)
                    end if

                    if (associated(PropertyX%Mass_created)) then
                        deallocate(PropertyX%Mass_created, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                               &
                            stop 'KillSedimentProperties - ModuleSedimentProperties - ERR29'
                        nullify   (PropertyX%Mass_created)
                    endif

                    if (associated(PropertyX%Evolution%Bioturbation%Coef))then
                        deallocate(PropertyX%Evolution%Bioturbation%Coef, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                               &
                            stop 'KillSedimentProperties - ModuleSedimentProperties - ERR30'
                        nullify   (PropertyX%Evolution%Bioturbation%Coef)
                    end if


                    PropertyX => PropertyX%Next

                end do do1

                nullify(Me%FirstProperty,Me%LastProperty)


                if (Me%Coupled%Partition%Yes) then

                    deallocate(Me%PartitionMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR31'
                    nullify   (Me%PartitionMatrix)
                    
                    deallocate(Me%IndependentTerm, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR32'
                    nullify   (Me%IndependentTerm)
                    

                    call KillLUD(Me%ObjLUD, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR33'
                end if

                if (Me%Coupled%SedimentQuality%Yes) then
                    
                    deallocate(Me%DissolvedToParticulate3D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR34'
                    nullify   (Me%DissolvedToParticulate3D)

                
                    call KillInterface(Me%ObjInterface, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                  &
                        stop 'KillSedimentProperties - ModuleSedimentProperties - ERR35'
                    
                end if

                call DeallocateInstance

                SedimentPropertiesID  = 0
                STAT_                 = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


    end subroutine KillSedimentProperties

    !--------------------------------------------------------------------------
  
    subroutine Write_Final_SedProperties_HDF

        !Local--------------------------------------------------------------
        type(T_Property),           pointer     :: Property
        integer                                 :: ObjHDF5
        character (Len = StringLength)          :: PropertyName
        integer                                 :: WorkILB, WorkIUB
        integer                                 :: WorkJLB, WorkJUB
        integer                                 :: WorkKLB, WorkKUB
        integer                                 :: STAT_CALL
        integer(4)                              :: HDF5_CREATE
        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        ObjHDF5 = 0
        !Opens HDF5 File
        call ConstructHDF5 (ObjHDF5,                                                     &
                            trim(Me%Files%Final)//"5",                                   &
                            HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR01'


        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                        &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR02'

        call HDF5WriteData   (ObjHDF5, "/Grid", "WaterPoints3D", "-",                    &
                              Array3D = Me%ExternalVar%WaterPoints3D,                    &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR03'

        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB+1, WorkJLB,                      &
                              WorkJUB+1, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR04'

        call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionX", "m",                      &
                              Array2D = Me%ExternalVar%XX_IE,                            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR05'

        call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionY", "m",                      &
                              Array2D = Me%ExternalVar%YY_IE,                            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR06'

        !Writes SZZ
        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                         &
                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR07'

        call HDF5WriteData  (ObjHDF5, "/Grid", "VerticalZ",                              &
                             "m", Array3D = Me%ExternalVar%SZZ,                          &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR08'

        !Writes OpenPoints
        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                                  &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,                         &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR09'

        call HDF5WriteData  (ObjHDF5, "/Grid", "OpenPoints",                             &
                             "-", Array3D = Me%ExternalVar%OpenPoints3D,                 &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR10'



        Property => Me%FirstProperty

        do while (associated(Property))
      
            PropertyName = trim(adjustl(Property%ID%name))
            
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB,                     &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR11'

            !Final concentration
            call HDF5WriteData  (ObjHDF5, "/Concentration/"//Property%ID%Name,           &
                                 Property%ID%Name,                                       &
                                 Property%ID%Units,                                      &
                                 Array3D = Property%Concentration,                       &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR12'


            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB,                     &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR13'


            if (associated(Property%Mass_created))then

                call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                          &
                                     WorkJLB, WorkJUB, WorkKLB, WorkKUB,                 &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                               &
                    stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR15'

                call HDF5WriteData  (ObjHDF5,                                            &
                                     "/Mass_created/"//Property%ID%Name,                 &
                                     Property%ID%Name,                                   &
                                     Property%ID%Units,                                  &
                                     Array3D = Property%Mass_created,                    &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                               &
                    stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR16'

            endif
            
             if (Property%ID%IDNumber==SeagrassesRoots_) then
                
                    call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                             &
                                         WorkJLB, WorkJUB, WorkKLB, WorkKUB,                    &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR17'

                    call HDF5WriteData  (ObjHDF5, "/Results/"//"seagrasses roots biomass",       &
                                         "seagrasses roots biomass", "gdw/m2",                    &
                                         Array2D = Me%SeagrassesRoots%Biomass)
                    if (STAT_CALL /= SUCCESS_) stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR18'
                    

               
               endif 

            Property => Property%Next

        enddo

        nullify (Property)
   
        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR17'

        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'Write_Final_SedProperties_HDF - ModuleSedimentProperties - ERR18'

    end subroutine Write_Final_SedProperties_HDF

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Local-----------------------------------------------------------------
        type (T_SedimentProperties), pointer    :: AuxSedimentProperties
        type (T_SedimentProperties), pointer    :: PreviousSedimentProperties

        !Updates pointers
        if (Me%InstanceID == FirstObjSedimentProperties%InstanceID) then
            FirstObjSedimentProperties      => FirstObjSedimentProperties%Next
        else
            PreviousSedimentProperties      => FirstObjSedimentProperties
            AuxSedimentProperties           => FirstObjSedimentProperties%Next
            do while (AuxSedimentProperties%InstanceID /= Me%InstanceID)
                PreviousSedimentProperties  => AuxSedimentProperties
                AuxSedimentProperties       => AuxSedimentProperties%Next
            enddo

            !Now update linked list
            PreviousSedimentProperties%Next => AuxSedimentProperties%Next

        endif
        
        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (SedimentPropertiesID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentPropertiesID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (SedimentPropertiesID > 0) then
            call LocateObjSedimentProperties (SedimentPropertiesID)
            ready_ = VerifyReadLock (mSEDIMENTPROPERTIES_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready


    !--------------------------------------------------------------------------

    subroutine LocateObjSedimentProperties (SedimentPropertiesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: SedimentPropertiesID

        !Local-----------------------------------------------------------------

        Me => FirstObjSedimentProperties
        do while (associated (Me))
            if (Me%InstanceID == SedimentPropertiesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleSedimentProperties - LocateObjSedimentProperties - ERR01'

    end subroutine LocateObjSedimentProperties

    !--------------------------------------------------------------------------   
  
end module ModuleSedimentProperties

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
