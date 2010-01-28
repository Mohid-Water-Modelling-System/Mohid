!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : PorousMediaProperties
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as PorousMediaProperties to create new modules
!
!------------------------------------------------------------------------------

!
!Units in porous media properties
!   Transported properties (soluble)  : g/m3 (or mg/l)  (needs to convert concentrations to SedimentQuality and PREEQC at entrance and exit)
!   Adsorbed properties (non soluble) : ug/kgsoil       (needs to convert concentrations to PREEQC at entrance and exit)
!
Module ModulePorousMediaProperties

    use ModuleGlobalData
    use ModuleStopWatch
    use ModuleFunctions
    use ModuleTime
    use ModuleHDF5
    use ModuleEnterData
    use ModuleProfile,          only : StartProfile, WriteProfile, KillProfile
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData,    &
                                       KillGridData
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie         
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, GetGridCellArea,               &
                                       WriteHorizontalGrid, UnGetHorizontalGrid
    use ModuleBasinGeometry,    only : GetBasinPoints, GetRiverPoints,  UnGetBasin 
                                       
    use ModuleFillMatrix,         only : ConstructFillMatrix, GetDefaultValue,             &
                                         KillFillMatrix, ModifyFillMatrix
    use ModuleGeometry
    use ModuleMap
    use ModulePorousMedia,        only : GetOldWaterContent, GetWaterContent, GetFluxU,    &
                                         GetFluxV, GetFluxW, GetUnsatW, GetUnsatV,         &
                                         GetUnsatU, UnGetPorousMedia,                      &
                                         GetThetaS, GetGWFlowToChannels,   &
                                         GetGWLayer, GetPotentialInfiltration
    use ModuleInterface,          only : ConstructInterface, Modify_Interface
    use ModuleAdvectionDiffusion, only : StartAdvectionDiffusion, AdvectionDiffusion,      &
                                         GetAdvFlux, GetDifFlux, GetBoundaryConditionList, &
                                         UngetAdvectionDiffusion, KillAdvectionDiffusion

   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !PorousMediaProperties
    public  :: ConstructPorousMediaProperties
    private ::      AllocateInstance
    private ::      ReadFileNames
    private ::      ReadGlobalOptions
    private ::      Construct_PropertyList
    private ::      ConstructHDF
    private ::      ConstructTimeSerie
    private ::      CoupleSoilQuality
#ifdef _PHREEQC_       
    private ::      CoupleSoilChemistry
#endif
    private ::      StartAdvectionDiffusion

    !Selector
    public  :: GetPMPnProperties
    public  :: GetPMPPropertiesIDByIdx    
    public  :: GetPMPConcentration
    public  :: GetPMPCoupled
    public  :: SetDNConcPMP              !PMP gets conc from Drainage network
    public  :: SetInfColConcPMP          !PMP gets infcol conc from basin
    public  :: SetVegetationPMProperties !PMP gets conc from vegetation
    public  :: SetWindVelocity                 
    public  :: UnGetPorousMediaProperties
    
    !Modifier
    public  :: ModifyPorousMediaProperties
    private ::      ActualizePropertiesFromFile
    private ::      InterfaceFluxes
    private ::      AdvectionDiffusionProcesses_PMP !Explicit in porous media properties module
    private ::          ModifyAdvectionDiffusion_W  !Vertical direction
    private ::          ModifyAdvectionDiffusion_3D 
    private ::      AdvectionDiffusionProcesses_AD  !Using Module Advection Diffusion
    private ::          ComputeVolumes  
    private ::          AdvectionDiffusion  
    private ::      SoilQualityProcesses
#ifdef _PHREEQC_    
    private ::      SoilChemistryProcesses
#endif        
    private ::      OutPut_TimeSeries
    private ::      Output_HDF
    
    !Destructor
    public  :: KillPorousMediaProperties                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjPorousMediaProperties 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetPorousMediaProperties3D_I
    private :: UnGetPorousMediaProperties3D_R8
    private :: UnGetPorousMediaProperties3D_R4
    interface  UnGetPorousMediaProperties
        module procedure UnGetPorousMediaProperties3D_I
        module procedure UnGetPorousMediaProperties3D_R8
        module procedure UnGetPorousMediaProperties3D_R4
    end interface  UnGetPorousMediaProperties


    !Parameters-----------------------------------------------------------------

    real,    parameter :: WaterReferenceDensity = 1000. ![kg/m3]
    
    integer, parameter :: DirectionX = 1
    integer, parameter :: DirectionY = 2

    character(LEN = StringLength), parameter :: prop_block_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter :: prop_block_end       = '<endproperty>'
    integer, parameter                :: AdvDif_ModulePMP_  = 1
    integer, parameter                :: AdvDif_ModuleAD_   = 2
    integer, parameter                :: AdvDif_Upwind_     = 1
    integer, parameter                :: AdvDif_CentralDif_ = 2
    integer, parameter                :: AdvDif_Diff_Jury_  = 1
    integer, parameter                :: AdvDif_Diff_Old_   = 2    
    !Types---------------------------------------------------------------------
    
    private :: T_PorousMediaProperties

    type T_RelatedID
        integer                       :: IDNumber = -1
        character(LEN = StringLength) :: name     = ''   
    end type T_RelatedID

    type T_ID
        integer                       :: IDNumber
        character(LEN = StringLength) :: name
        character(LEN = StringLength) :: description
        character(LEN = StringLength) :: units
    end type T_ID

    type T_Property_3D
        type(T_PropertyID)               :: ID
        real, pointer, dimension (:,:,:) :: Field
        real                             :: Scalar
    end type T_Property_3D


    type T_ExtVar
        !Map
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D
        integer, pointer, dimension(:,:,:)      :: LandPoints3D
        integer, dimension(:,:), pointer        :: BasinPoints
        integer, dimension(:,:), pointer        :: RiverPoints
        integer, pointer, dimension(:,:,:)      :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:)      :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:)      :: ComputeFacesW3D         !from basin
        real                                        :: PorousMediapropDT
        type(T_Time)                                :: Now
        type(T_Time)                                :: BeginTime
        type(T_Time)                                :: EndTime
   
        ! from porousMedia
        real,    dimension(:,:,:), pointer          :: UnSatW
        real,    dimension(:,:,:), pointer          :: UnSatV
        real,    dimension(:,:,:), pointer          :: UnSatU
        real(8), pointer, dimension(:,:,:)          :: WaterContent
        real(8), pointer, dimension(:,:,:)          :: WaterContentOld
!        real(8), pointer, dimension(:,:)            :: WaterColumn
        real(8), pointer, dimension(:,:)            :: InfiltrationColumn
        real(8), pointer, dimension(:,:,:)          :: CellVolume
        real(8), pointer, dimension(:,:,:)          :: CellWaterMass
        real(8), dimension(:,:,:), pointer          :: FluxU
        real(8), dimension(:,:,:), pointer          :: FluxV
        real(8), dimension(:,:,:), pointer          :: FluxW
        real,    pointer, dimension(:,:,:)          :: ThetaS        
        real   , pointer, dimension(:,:  )          :: Area
        real                                        :: DT
        real   , pointer, dimension(:,:,:)          :: DWZ
        real   , pointer, dimension(:,:,:)          :: DZZ
!        real   , pointer, dimension(:,:,:)          :: DZI
!        real   , pointer, dimension(:,:,:)          :: DZE
!        real   , pointer, dimension(:,:  )          :: DVX
        real   , pointer, dimension(:,:  )          :: DZX
!        real   , pointer, dimension(:,:  )          :: DUY
        real   , pointer, dimension(:,:  )          :: DZY
        real   , pointer, dimension(:,:  )          :: DXX
        real   , pointer, dimension(:,:  )          :: DYY        
        real   , pointer, dimension(:,:  )          :: Topography  
        real ,   pointer, dimension(:,:,:)          :: SZZ
        real(8), pointer, dimension(:,:)            :: FlowToChannels
        integer, pointer, dimension(:,:)            :: GWLayer
        
        !from vegetation
        logical                                     :: ComputeVegInterfaceFluxes 
        logical, dimension(:,:  ), pointer          :: SoilFluxesActive
        real,    dimension(:,:  ), pointer          :: GrazingBiomass
        real,    dimension(:,:  ), pointer          :: GrazingNitrogen
        real,    dimension(:,:  ), pointer          :: GrazingPhosphorus
        real,    dimension(:,:  ), pointer          :: ManagementAerialBiomass
        real,    dimension(:,:  ), pointer          :: ManagementNitrogen
        real,    dimension(:,:  ), pointer          :: ManagementPhosphorus
        real,    dimension(:,:  ), pointer          :: ManagementRootBiomass
        real,    dimension(:,:  ), pointer          :: DormancyBiomass
        real,    dimension(:,:  ), pointer          :: DormancyNitrogen
        real,    dimension(:,:  ), pointer          :: DormancyPhosphorus
        real,    dimension(:,:  ), pointer          :: FertilNitrateSurface
        real,    dimension(:,:  ), pointer          :: FertilNitrateSubSurface
        real,    dimension(:,:  ), pointer          :: FertilAmmoniaSurface
        real,    dimension(:,:  ), pointer          :: FertilAmmoniaSubSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicNSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicNSubSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicPSurface
        real,    dimension(:,:  ), pointer          :: FertilOrganicPSubSurface
        real,    dimension(:,:  ), pointer          :: FertilMineralPSurface
        real,    dimension(:,:  ), pointer          :: FertilMineralPSubSurface
        real,    dimension(:,:  ), pointer          :: NitrogenFraction
        real,    dimension(:,:  ), pointer          :: PhosphorusFraction
        real,    dimension(:,:,:), pointer          :: NitrogenUptake
        real,    dimension(:,:,:), pointer          :: PhosphorusUptake
        real,    dimension(:,:  ), pointer          :: RootDepth
        logical                                     :: Grazing
        logical                                     :: Management
        logical                                     :: Dormancy
        logical                                     :: Fertilization
        logical                                     :: ModelNitrogen
        logical                                     :: ModelPhosphorus
        logical                                     :: GrowthModel
        logical                                     :: CoupledVegetation
        real                                        :: VegetationDT
        
        logical                                     :: CoupledDN  = .false.
        !from basin
        real,    dimension(:,:  ), pointer          :: WindVelocity2D  !m/s
        real,    dimension(:,:,:), pointer          :: WindVelocity3D  !km/day
     end type T_ExtVar

    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime
        integer                                 :: NextOutPut
        integer                                 :: Number
        logical                                 :: Yes = .false.
        logical                                 :: TimeSerie_ON
        logical                                 :: HDF_ON
        logical                                 :: Profile_ON
    end type T_OutPut

    type T_AdvectionDiffusion   
        !--For AdvectionDiffusion module use
        integer                                :: BoundaryCondition
        real                                   :: SchmidtNumberH
        real                                   :: SchmidtCoefV
        real                                   :: SchmidtBackgroundV
        real                                   :: DiffusionH_imp_exp
        real                                   :: ImplicitH_direction
        logical                                :: Nulldif          = .false.
        logical                                :: NumericStability = .false.
        real                                   :: VolumeRelMax
        integer                                :: AdvMethodH, TVDLimitationH
        integer                                :: AdvMethodV, TVDLimitationV
        logical                                :: Upwind2H, Upwind2V  
        logical                                :: Adv_Dif_Explicit                      
        !--For both models use
        real                                   :: Molecular_Diff_Coef 
                           
    end type T_AdvectionDiffusion

    type T_Partition
        logical                                 :: NonLinear
        character(LEN = StringLength)           :: NonLinear_ks_Units
        type(T_Property_3D)                     :: Nu            
        type(T_Property_3D)                     :: Be          
        type(T_Property_3D)                     :: ks
        type(T_Property_3D)                     :: PartitionRate
        type(T_Property_3D)                     :: Fraction 
        character (LEN = StringLength)          :: Partition_Couple
    end type T_Partition

    type T_Evolution
        logical                                 :: Variable = .false.
        real                                    :: DTInterval
        type(T_Time)                            :: LastCompute
        type(T_Time)                            :: NextCompute
        logical                                 :: SoilQuality
        logical                                 :: SoilChemistry
        logical                                 :: Partitioning
        logical                                 :: CationExchangeProcess
        logical                                 :: ChemEquilibriumProcess
        logical                                 :: AdvectionDiffusion
        logical                                 :: SoilWaterFluxes
        logical                                 :: Macropores
        logical                                 :: MinConcentration
        type (T_AdvectionDiffusion)             :: AdvDiff
        type (T_Partition                    )  :: Partition
    end type T_Evolution
    



    type T_Files
        character(PathLength)                   :: InitialFile
        character(PathLength)                   :: DataFile
        character(PathLength)                   :: FinalFile
        character(PathLength)                   :: TransientHDF
        character(PathLength)                   :: DataSedimentQualityFile
        integer                                 :: AsciiUnit        
    end type T_Files    

    type T_Property
        type (T_PropertyID)                     :: ID
        real(8), dimension(:,:,:), pointer      :: Concentration            => null()
        real(8), dimension(:,:,:), pointer      :: ConcentrationOld         => null()
        real(8), dimension(:,:), pointer        :: ConcentrationOnInfColumn      => null()
        real(8), dimension(:,:), pointer        :: ConcentrationDN               => null()
        real, pointer, dimension(:,:,:)         :: Mass_Created
        real(8),    pointer, dimension(:,:,:)   :: ViscosityU
        real(8),    pointer, dimension(:,:,:)   :: ViscosityV
!        real,    pointer, dimension(:,:,:)      :: DiffusivityW          
        type (T_Property), pointer              :: Next, Prev                     => null()
        logical                                 :: Particulate
        type (T_Evolution)                      :: Evolution
        real(8), pointer, dimension(:,:,:)      :: Diffusivity
        real(8), pointer, dimension(:,:,:)      :: Diff_Turbulence_H
        real(8), pointer, dimension(:,:,:)      :: Diff_Turbulence_V
        real(8), pointer, dimension(:,:,:)      :: Viscosity

        logical                                 :: Old     = .false.
        real                                    :: MinValue        = FillValueReal
        logical                                 :: TimeSerie        = .false.
        logical                                 :: BoxTimeSerie     = .false.
        logical                                 :: BoxTimeSerie2D   = .false.
        logical                                 :: OutputHDF        = .false.
!        type (T_RelatedID)                      :: RelatedID
        
    end type T_Property

    type T_Coupled
        logical                                 :: SoilQuality          = .false. !Sediment source/sink model (Sediment Quality)
        real                                    :: SoilQuality_DT
        type (T_Time)                           :: SoilQuality_NextCompute

#ifdef _PHREEQC_        
        logical                                 :: SoilChemistry        = .false.  !Chemical reactions model (PhreeqC)
        real                                    :: SoilChemistry_DT
        type (T_Time)                           :: SoilChemistry_NextCompute
#endif        
        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: MinConcentration     = .false.
    end type T_Coupled

    type T_PorousMediaProperties
        integer                                     :: ObjTime                   = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjAdvectionDiffusion     = 0
        integer                                     :: ObjBasinGeometry     = 0
        integer                                     :: ObjPorousMedia       = 0
        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjMap               = 0
        integer                                     :: ObjGridData          = 0
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjtimeSerie         = 0
        integer                                     :: ObjSedimentQuality   = 0
        integer                                     :: ObjHDF5              = 0
        integer                                     :: ObjBottomTopography  = 0
        integer                                     :: ObjProfile           = 0
        integer                                     :: ObjInterface         = 0
#ifdef _PHREEQC_        
        integer                                     :: ObjPhreeqC                = 0
        integer                                     :: ObjInterfaceSoilChemistry = 0 
#endif        
        type (T_ExtVar)                             :: ExtVar
        logical                                     :: CheckGlobalMass      
        type (T_Files)                              :: Files
        type (T_OutPut)                             :: OutPut
        type (T_Property), pointer                  :: FirstProperty    => null() !Lúcia
        type (T_Property), pointer                  :: LastProperty        
        type (T_PorousMediaProperties), pointer     :: Next             => null() !Lúcia
        type (T_Coupled)                            :: Coupled
        type (T_Time)                               :: LastOutputHDF5

        logical                                     :: PorousMediaProperties
        real,    pointer, dimension(:,:,:)          :: Volume   
        integer                                     :: PropertiesNumber    = 0
        real   , pointer, dimension(:,:,:)          :: DissolvedToParticulate3D
        real                                        :: ResidualTime
        
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D

        real(8), dimension(:, :, :),  pointer       :: Matrix  
        
          
        type(T_Property_3D)                         :: Disper_Trans
        type(T_Property_3D)                         :: Disper_Longi
               
        real(8), pointer, dimension(:,:,:)          :: MassFluxesX
        real(8), pointer, dimension(:,:,:)          :: MassFluxesY
        real(8), pointer, dimension(:,:,:)          :: MassFluxesZ

        integer                                     :: AdvDiff_Module        ! 1 - PorousMediaProperties, 2 - AdvectionDiffusion
        integer                                     :: AdvDiff_Test
        integer                                     :: AdvDiff_SpatialMethod ! 1 - Upwind; 2-Central Differences
        logical                                     :: AdvDiff_Explicit      !
        logical                                     :: AdvDiff_CheckCoefs    !
        logical                                     :: AdvDiff_ComputeTransport3D
        integer                                     :: AdvDiff_DiffMethod    !  1 - Jury based, 2 - AdvectionDiffusion    
        !--For PorousMediaProperties Advection-Diffusion Method
        real,    pointer, dimension(:,:,:)          :: DifusionNumber
        real,    pointer, dimension(:,:,:)          :: ReynoldsMNumber    
        
       
        real(8), pointer, dimension(:,:,:)          :: WaterVolume
        real(8), pointer, dimension(:,:,:)          :: WaterVolumeOld                
        real(8), pointer, dimension(:,:,:)          :: WaterVolumeCorr
!        real(8), pointer, dimension(:,:,:)          :: WaterVolumeOldCorr   
        real(8), pointer, dimension(:,:,:)          :: FluxWCorr             

    end type  T_PorousMediaProperties

    !Global Module Variables
    type (T_PorousMediaProperties), pointer                         :: FirstObjPorousMediaProperties
    type (T_PorousMediaProperties), pointer                         :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructPorousMediaProperties(ObjPorousMediaPropertiesID,                 &
                                              ComputeTimeID,                              &
                                              HorizontalGridID,                           &
                                              HorizontalMapID,                            &
                                              BasinGeometryID,                            &
                                              PorousMediaID,                              &
                                              GeometryID,                                 &
                                              MapID,                                      &
                                              CoupledDN,                                  &
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: BasinGeometryID
        integer                                         :: PorousMediaID
        integer                                         :: GeometryID
        integer                                         :: MapID
        logical, optional                               :: CoupledDN
        integer, optional, intent(OUT)                  :: STAT 
        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_,STAT_CALL
        !------------------------------------------------------------------------
                                    

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mPorousMediaProperties_)) then
            nullify (FirstObjPorousMediaProperties)
            call RegisterModule (mPorousMediaProperties_) 
        endif

        call Ready(ObjPorousMediaPropertiesID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then


            call AllocateInstance
            
            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ComputeTimeID   )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMap_,            MapID           )
        

            if (present(CoupledDN)) then
                Me%ExtVar%CoupledDN = CoupledDN
            endif
            
            
            call ReadFileNames


            !Constructs the DataFile
            call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR01'                
           
            call ReadGlobalOptions

            call AllocateVariables

            call Construct_PropertyList
            
!            if (Me%Coupled%SoilQuality) then
!                call Construct_InitialFields
!            endif
        
            call ConstructHDF    
    
            call ConstructTimeSerie
            
            if (Me%Coupled%AdvectionDiffusion .and. (Me%AdvDiff_Module == AdvDif_ModulePMP_).and. Me%AdvDiff_CheckCoefs) then
                call ConstructAsciiFile
            endif
            
    !       call ConstructProfileOutput   em teste
            
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR010'



            !Couple nutrient, carbon and oxygen sources and sinks model
            if (Me%Coupled%SoilQuality) then
                call CoupleSoilQuality
            endif
            
#ifdef _PHREEQC_            
            !Couple soil chemical model
            if (Me%Coupled%SoilChemistry) then
                call CoupleSoilChemistry
            endif
#endif
            
            if (Me%Coupled%AdvectionDiffusion .and. Me%AdvDiff_Module == AdvDif_ModuleAD_) then !Uses AdvectionDiffusion module
            
                call StartAdvectionDiffusion(Me%ObjAdvectionDiffusion, &
                                             Me%ObjGeometry,           &
                                             Me%ObjHorizontalMap,      &
                                             Me%ObjHorizontalGrid,     &
                                             Me%ObjTime,               &
                                             STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR50'

            endif

            !Returns ID
            ObjPorousMediaPropertiesID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR060' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructPorousMediaProperties
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_PorousMediaProperties), pointer                         :: NewObjPorousMediaProperties
        type (T_PorousMediaProperties), pointer                         :: PreviousObjPorousMediaProp


        !Allocates new instance
        allocate (NewObjPorousMediaProperties)
        nullify  (NewObjPorousMediaProperties%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPorousMediaProperties)) then
            FirstObjPorousMediaProperties         => NewObjPorousMediaProperties
            Me                    => NewObjPorousMediaProperties
        else
            PreviousObjPorousMediaProp      => FirstObjPorousMediaProperties
            Me                    => FirstObjPorousMediaProperties%Next
            do while (associated(Me))
                PreviousObjPorousMediaProp  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjPorousMediaProperties
            PreviousObjPorousMediaProp%Next => NewObjPorousMediaProperties
        endif

        Me%InstanceID = RegisterNewInstance (mPorousMediaProperties_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ReadFileNames

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
!        integer                                     :: iflag

        !Reads the name of the data file from nomfich
        call ReadFileName ('POROUS_PROP_DATA', Me%Files%DataFile, "PorousMedia Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('POROUS_PROP_HDF', Me%Files%TransientHDF, "PorousMedia HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01b'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUS_PROP_FIN', Me%Files%FinalFile, "PorousMedia Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01c'
   

    end subroutine ReadFileNames
    
    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        type(T_Property_3D), pointer            :: Scalar3D
        !Begin-----------------------------------------------------------------

        !Geometry Size
        call GetGeometrySize    (Me%ObjGeometry,             &    
                                 Size     = Me%Size,         &
                                 WorkSize = Me%WorkSize,     &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR10'

        Me%Size2D%ILB = Me%Size%ILB
        Me%Size2D%IUB = Me%Size%IUB
        Me%Size2D%JLB = Me%Size%JLB
        Me%Size2D%JUB = Me%Size%JUB

        call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR20'

        call GetComputeTimeLimits(Me%ObjTime,                      &
                                  EndTime   = Me%ExtVar%EndTime,   &
                                  BeginTime = Me%ExtVar%BeginTime, &
                                  STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR30'

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

        !Needs keyword description and value options definition
        call GetData(Me%AdvDiff_Module,                            &   
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'ADVDIFF_MODULE',              &
                     Default      = AdvDif_ModulePMP_,             &
                     ClientModule = 'ModulePorousMediaProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR40'


		if (Me%AdvDiff_Module == AdvDif_ModulePMP_) then !Uses ModulePorousMediaProperties for Advection-Diffusion calculation

	        call GetData(Me%AdvDiff_SpatialMethod,                     &   !Lúcia
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromFile,                      &
	                     keyword      = 'ADVDIFF_SPATIAL_METHOD',      &
	                     Default      = AdvDif_CentralDif_,            &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR50'

            call GetData(Me%AdvDiff_ComputeTransport3D,                                      &   
                         Me%ObjEnterData, iflag,                                              &
                         SearchType = FromFile,                                               &
                         keyword    = 'ADVDIFF_3D',                                           &
                         Default    = .false.,                                                &                                           
                         ClientModule ='ModulePorousMediaProperties',                         &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR60'
            
            if (Me%AdvDiff_ComputeTransport3D) then
                call GetData(Me%AdvDiff_DiffMethod,                                               &   
                             Me%ObjEnterData, iflag,                                              &
                             SearchType = FromFile,                                               &
                             keyword    = 'ADVDIFF_DIFF_METHOD',                                  &
                             Default    = AdvDif_Diff_Jury_,                                                      &                                           
                             ClientModule ='ModulePorousMediaProperties',                         &
                             STAT       = STAT_CALL)            
                if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR070'            
	        endif
            
            call GetData(Me%AdvDiff_CheckCoefs,                        &   !Eduardo
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromFile,                      &
                         keyword      = 'ADVDIFF_CHECK_COEFS',         &
                         Default      = .false.,                       &
                         ClientModule = 'ModulePorousMediaProperties', &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR80'

		elseif (Me%AdvDiff_Module == AdvDif_ModuleAD_) then ! Uses ModuleAdvectionDiffusion for Advection-Diffusion calculation

	        call GetData(Me%AdvDiff_Explicit,                          &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromFile,                      &
	                     keyword      = 'ADVDIFF_EXPLICIT',            &
	                     Default      = .true.,                	       &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_)  &
	            stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR90'
	            
	    else
	        write(*,*)'Advection diffusion module to be used unrecognized,'
	        write(*,*)'Please check ADVDIFF_MODULE keyword'
	        stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR100'

		endif

        Scalar3D => Me%Disper_Longi
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                              block_begin = '<begin_dispersion_long>',          &
                              block_end   = '<end_dispersion_long>')
        
        Scalar3D => Me%Disper_Trans
        call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                              block_begin = '<begin_dispersion_trans>',         &
                              block_end   = '<end_dispersion_trans>')            

    end subroutine ReadGlobalOptions        

    !--------------------------------------------------------------------------

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB 
        integer                                         :: KLB, KUB 
        integer                                         :: STAT_CALL
        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KUB = Me%Size%KUB
        KLB = Me%Size%KLB
        
               
        !Water Content---------------------------------------------------------
        allocate (Me%Volume                  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%DifusionNumber          (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%ReynoldsMNumber         (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%ExtVar%WindVelocity3D   (ILB:IUB,JLB:JUB,KLB:KUB))
        
#ifdef _PHREEQC_
        allocate (Me%ExtVar%CellWaterMass    (ILB:IUB,JLB:JUB,KLB:KUB))
#endif        

        allocate (Me%WaterVolume          (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%WaterVolumeOld       (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%WaterVolumeCorr          (ILB:IUB,JLB:JUB,KLB:KUB))
    !    allocate (Me%WaterVolumeOldCorr       (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%FluxWCorr                (ILB:IUB,JLB:JUB,KLB:KUB))
        
        Me%WaterVolume          = 0.
        Me%WaterVolumeOld       = 0.
        Me%WaterVolumeCorr      = 0.
        Me%FluxWCorr            = 0.
        

    end subroutine AllocateVariables


    !--------------------------------------------------------------------------
    
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
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR01'


            case(FromBlockInBlock)

                if(.not. present(ClientNumber))then
                    stop 'ConstructScalar3D - ModuleSoilProperties - ERR02'
                end if
                
                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber, block_begin, block_end,  &
                                           BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR03'

        end select

        if(BlockFound)then

            allocate(Scalar3D%Field(ILB:IUB, JLB:JUB, KLB:KUB))

            call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR03.1'

            call ConstructFillMatrix  (PropertyID           = Scalar3D%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill3D       = Me%ExtVar%WaterPoints3D,          &
                                       Matrix3D             = Scalar3D%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR04'


            call GetDefaultValue(Scalar3D%ID%ObjFillMatrix, Scalar3D%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR05'

            call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR5.1'

            call KillFillMatrix(Scalar3D%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR06'

            if(ExtractType == FromBlockInBlock)then
                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR07'
            end if


            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, BlockClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR08'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR09'
            end if
        
        else
            write(*,*) 'Block not present:', block_begin, block_end
            stop 'ConstructScalar3D - ModulePorousMediaProperties - ERR10'
        end if

   
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
                    
                    !Construct a New Property 
                    Call Construct_Property(NewProperty)

                    !Add new Property to the SoilProperties List 
                    Call Add_Property(NewProperty)

                else cd2

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_PropertyList - ModulePorousMediaProeprties - ERR01'
                    exit do1    !No more blocks
                
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_PropertyList - ModulePorousMediaProeprties - ERR02'
            
            else cd1
                
                stop 'Construct_PropertyList - ModulePorousMediaProeprties - ERR03'
            
            end if cd1
        
        end do do1

    end subroutine Construct_PropertyList

    !----------------------------------------------------------------------------    
    
        subroutine Construct_Property(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer           :: NewProperty

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_)stop 'Construct_Property - ModulePorousMediaProeprties - ERR00'
        
        nullify(NewProperty%Prev, NewProperty%Next)
        nullify(NewProperty%Concentration)
        nullify(NewProperty%ConcentrationOnInfColumn)

        call ConstructPropertyID            (NewProperty%ID, Me%ObjEnterData, FromBlock)

        call Construct_PropertyState        (NewProperty)

        call Construct_PropertyValues       (NewProperty)

        call Construct_PropertyEvolution    (NewProperty)

        call Construct_PropertyOutPut       (NewProperty)

    end subroutine Construct_Property
    
    !-------------------------------------------------------------------------------    
    
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

    subroutine Construct_PropertyState(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer       :: NewProperty

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL, iflag

        !----------------------------------------------------------------------
        

        !<BeginKeyword>
            !Keyword          : PARTICULATE
            !<BeginDescription>
            !<EndDescription>
            !Type             : logical   
            !Default          : Dissolved
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : From Block
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Particulate,                                            &
                     Me%ObjEnterData,  iflag,                                            &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'PARTICULATE',                                       &
                     ClientModule = 'ModulePorousMediaProeprties',                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR01'
        if(iflag == 0)              stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR02'

        !Not used the function so this text was commented
!        if (NewProperty%Particulate)then
!            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
!                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
!                write(*,*) 'recognised as PARTICULATE'
!                stop 'Construct_PropertyState - ModulePorousMediaProeprties - ERR03'
!            end if
!        endif

    end subroutine Construct_PropertyState

    !--------------------------------------------------------------------------

    subroutine Construct_PropertyEvolution(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        real                                        :: ErrorAux, AuxFactor, DTAux
        real                                        :: ModelDT
        !----------------------------------------------------------------------

        call GetData(NewProperty%Evolution%AdvectionDiffusion,                           &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'ADVECTION_DIFFUSION',                               &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR10'

        if (NewProperty%Evolution%AdvectionDiffusion) then
            Me%Coupled%AdvectionDiffusion = .true.
            NewProperty%Evolution%Variable = .true.
        endif
        
        if (NewProperty%Evolution%AdvectionDiffusion) then

            call ReadAdvectionDiffusionParameters (NewProperty)
        
        end if

        call ConstructPropertyDiffusivity (NewProperty)


        !<BeginKeyword>
            !Keyword          : SOIL_QUALITY
            !<BeginDescription>
               ! Property has the Soil quality model (sediment quality) as sink and source
            !<EndDescription>
            !Type             : Boolean
            !Default          : .false.
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%SoilQuality,                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SOIL_QUALITY',                                      &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR20'

        if (NewProperty%Evolution%SoilQuality) then
            Me%Coupled%SoilQuality     = .true.
            NewProperty%Evolution%Variable = .true.
        endif


#ifdef _PHREEQC_
        !<BeginKeyword>
            !Keyword          : SOIL_CHEMISTRY
            !<BeginDescription>
               ! Property has the Soil chemistry model (PHREEQC)
            !<EndDescription>
            !Type             : Boolean
            !Default          : .false.
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%SoilChemistry,                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SOIL_CHEMISTRY',                                    &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR30'

        if (NewProperty%Evolution%SoilChemistry) then
            Me%Coupled%SoilChemistry       = .true.
            NewProperty%Evolution%Variable = .true.
        endif
#endif

        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR40'

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModulePorousMediaProperties',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR050'
                                       
            
            if (NewProperty%Evolution%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR60'

            elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR70'
                endif

                !Run period in seconds
                DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%Evolution%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR80'
                end if
            endif

            NewProperty%Evolution%NextCompute = Me%ExtVar%Now + NewProperty%Evolution%DTInterval

        else

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%Evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution     

    !--------------------------------------------------------------------------
    
    subroutine ReadAdvectionDiffusionParameters (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External--------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: MassConservation
        integer                   :: ImposedValue
        integer                   :: NullGradient, CyclicBoundary
        integer                   :: Orlanski, MassConservNullGrad

        !Local-----------------------------------------------------------------
        integer                   :: iflag, BoundaryCondition

        !----------------------------------------------------------------------

        call GetData(NewProperty%Evolution%AdvDiff%Molecular_Diff_Coef, &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'ADVDIFF_MOLECULAR_DIFF_COEF',      &
                     Default      = 0.0,                                &
                     ClientModule = 'ModulePorousMediaProperties',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR10'

cd1:    if (Me%AdvDiff_Module == AdvDif_ModuleAD_) then

	        call GetData(NewProperty%Evolution%AdvDiff%NumericStability, &
	                     Me%ObjEnterData, iflag,                         &
	                     SearchType   = FromBlock,                       &
	                     keyword      = 'ADVDIFF_NUM_STABILITY',         &
	                     Default      = .FALSE.,                         &
	                     ClientModule = 'ModulePorousMediaProperties',   &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR20'

	        call GetData(NewProperty%Evolution%AdvDiff%SchmidtNumberH, &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                     &
	                     keyword      = 'ADVDIFF_SCHMIDT_NUMBER_H',    &
	                     Default      = 1.0,                           &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR30'

	        call GetData(NewProperty%Evolution%AdvDiff%SchmidtCoefV,   &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                     &
	                     keyword      = 'ADVDIFF_SCHMIDT_COEF_V',      &
	                     Default      = 1.0,                           &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR40'

	        call GetData(NewProperty%Evolution%AdvDiff%SchmidtBackgroundV, &
	                     Me%ObjEnterData, iflag,                           &
	                     SearchType   = FromBlock,                         &
	                     keyword      = 'ADVDIFF_SCHMIDT_BACKGROUND_V',    &
	                     Default      = 0.,                                &
	                     ClientModule = 'ModulePorousMediaProperties',     &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR50'

	        call GetData(NewProperty%Evolution%AdvDiff%NullDif,        &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                     &
	                     keyword      = 'ADVDIFF_NULLDIF',             &
	                     Default      = .false.,                       &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR60'

	        call GetBoundaryConditionList(MassConservation    = MassConservation,    &
	                                      ImposedValue        = ImposedValue,        &
	                                      NullGradient        = NullGradient,        &
	                                      Orlanski            = Orlanski,            &
	                                      MassConservNullGrad = MassConservNullGrad, &
	                                      CyclicBoundary      = CyclicBoundary)

	        call GetData(BoundaryCondition,                            &
	                     Me%ObjEnterData,  iflag,                      &
	                     SearchType   = FromBlock,                     &
	                     keyword      = 'ADVDIFF_BOUNDARY_CONDITION',  &
	                     Default      = MassConservation,              &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR70'

	        ! By default it's imposed a value dependent only from the exterior
	        ! value and of the decay time. However this method doesn't conserve mass
	        ! when the water fluxes near the frontier are dominant

	        if (BoundaryCondition /= MassConservation     .and. &
	            BoundaryCondition /= ImposedValue         .and. &
	            BoundaryCondition /= NullGradient         .and. &
	            BoundaryCondition /= CyclicBoundary       .and. &
	            BoundaryCondition /= Orlanski             .and. &
	            BoundaryCondition /= MassConservNullGrad) &
	            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR80'

	        NewProperty%Evolution%AdvDiff%BoundaryCondition = BoundaryCondition

	        !By default the horizontal Diffusion discretization is explicit
	        NewProperty%Evolution%AdvDiff%DiffusionH_imp_exp  = ExplicitScheme

	        NewProperty%Evolution%AdvDiff%ImplicitH_Direction = DirectionX

	        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodH,     &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                      &
	                     keyword      = 'ADVDIFF_METHOD_H',            &
	                     Default      = UpwindOrder1,                  &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)

	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR90'

	        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationH, &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                      &
	                     keyword      = 'ADVDIFF_TVD_LIMIT_H',         &
	                     Default      = Superbee,                      &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)

	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR100'

	        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodV,     &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                      &
	                     keyword      = 'ADVDIFF_METHOD_V',            &
	                     Default      = UpwindOrder1,                  &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)

	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR110'

	        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationV, &
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromBlock,                      &
	                     keyword      = 'ADVDIFF_TVD_LIMIT_V',         &
	                     Default      = Superbee,                      &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)

	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR120'


	        call GetData(NewProperty%Evolution%AdvDiff%VolumeRelMax,   &
	                     Me%ObjEnterData, iflag,                       &
	                     Keyword      = 'ADVDIFF_VOLUME_RELATION_MAX', &
	                     Default      = 5.,                            &
	                     SearchType   = FromBlock,                      &
	                     ClientModule = 'ModulePorousMediaProperties', &
	                     STAT         = STAT_CALL)

	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR130'


	        if (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
	            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3 .or.&
	            NewProperty%Evolution%AdvDiff%AdvMethodH == P2_TVD) then
	            NewProperty%Evolution%AdvDiff%Upwind2H = .true.
	        else
	            NewProperty%Evolution%AdvDiff%Upwind2H = .false.
	        endif

	        if (NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder2 .or.&
	            NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder3 .or.&
	            NewProperty%Evolution%AdvDiff%AdvMethodV == P2_TVD) then
	            NewProperty%Evolution%AdvDiff%Upwind2V = .true.
	        else
	            NewProperty%Evolution%AdvDiff%Upwind2V = .false.
	        endif


	        if (.not. Me%AdvDiff_Explicit .and.&
	           (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
	            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3)) then

	            write(*,*) 'If the advection of mass in the horizontal is implicit'
	            write(*,*) 'the advection method can not be a second or third order upwind'
	            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR140.'

	        endif

	        if (.not. Me%AdvDiff_Explicit .and.&
	           (NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder2 .or.&
	            NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder3)) then

	            write(*,*) 'If the advection of mass in the vertical is implicit'
	            write(*,*) 'the advection method can not be a second or third order upwind'
	            stop 'ReadAdvectionDiffusionParameters - ModulePorousMediaProperties - ERR150.'

	        endif

		endif cd1

    end subroutine ReadAdvectionDiffusionParameters

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyDiffusivity (NewProperty)

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
        allocate (NewProperty%Diffusivity (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR10'
        NewProperty%Diffusivity       = 0. 
        
        if (Me%AdvDiff_Module == AdvDif_ModuleAD_) then !Module advection-Diffusion is used
           
            allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR20'
            NewProperty%Viscosity         = 0.
            
            allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR30'
            NewProperty%Diff_Turbulence_H = 0.

            allocate (NewProperty%Diff_Turbulence_V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR40'
            NewProperty%Diff_Turbulence_V = 0.                   
            
       elseif (Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) then !explicit calculation inside module porous media properties
            allocate (NewProperty%ViscosityU (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR60'
            
            NewProperty%ViscosityU  = 0.0
            
            allocate (NewProperty%ViscosityV (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR70'
        
            NewProperty%ViscosityV  = 0.0!     
            
            if (Me%AdvDiff_DiffMethod == AdvDif_Diff_Old_) then 

                allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR80'
                NewProperty%Viscosity         = 0.
                
                allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR90'
                NewProperty%Diff_Turbulence_H = 0.

                allocate (NewProperty%Diff_Turbulence_V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModulePorousMediaProperties - ERR100'
                NewProperty%Diff_Turbulence_V = 0.                                   
           
               
            endif
            
        endif

    end subroutine ConstructPropertyDiffusivity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
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
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR10'
        NewProperty%Concentration(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR20'
        NewProperty%ConcentrationOld(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOnInfColumn(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR40'
        NewProperty%ConcentrationOnInfColumn(:,:) = 0.

        
        if (Me%ExtVar%CoupledDN) then
            allocate(NewProperty%ConcentrationDN(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR50'
            NewProperty%ConcentrationDN(:,:) = FillValueReal
        endif

         !Eduardo Jauch 12nov2009
!        call GetData(NewProperty%RelatedID%name,                                            &
!                     Me%ObjEnterData, iflag,                                                &
!                     SearchType   = FromBlock,                                              &
!                     keyword      = 'RELATED',                                              &
!                     Default      = '',                                                     &                        
!                     ClientModule = 'ModuleSoilProperties',                                 &
!                     STAT         = STAT_CALL)              
!        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR80'
!        if ((iflag .NE. 0) .AND. (NewProperty%RelatedID%name .NE. '')) then
!            NewProperty%RelatedID%IDNumber = GetPropertyIDNumber(NewProperty%RelatedID%name)    
!        endif
!        !ToDo: After load all the properties, must check if "related properties" exist...        
!        !end

        call GetData(NewProperty%MinValue,                                                  &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'MIN_VALUE',                                            &
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR100'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
            Me%Coupled%MinConcentration = .true.
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR110'
            NewProperty%Mass_Created(:,:,:) = 0.
        endif

        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleVegetation',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR120'
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not. NewProperty%Old) then

            !Get water points
            call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR130'

            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill3D       = Me%ExtVar%WaterPoints3D,     &
                                       Matrix3D             = NewProperty%Concentration,        &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR140'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR0150'
            end if

!Lúcia
            call SetMatrixValue(NewProperty%ConcentrationOld, Me%Size, NewProperty%Concentration,Me%ExtVar%WaterPoints3D)

            call CheckFieldConsistence (NewProperty)

            call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR160'

        else

            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Water properties initial file written in HDF format  
            call ReadOldConcBoundariesHDF(NewProperty)

        end if   

    end subroutine Construct_PropertyValues

      !--------------------------------------------------------------------------

    subroutine ConstructProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL, iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
        integer                                             :: nProperties
        integer                                             :: n
        type (T_Property), pointer                          :: PropertyX

        nProperties = Me%PropertiesNumber 

        !Allocates PropertyList
        allocate(PropertyList(nProperties, 2))
       
        n=1
        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))

            !Fills up PropertyList
            PropertyList(n, 1) = trim(PropertyX%ID%Name)
            PropertyList(n, 2) = "m3/m3"

            n=n+1

            PropertyX=>PropertyX%Next

        enddo

        !----------------------------------------------------------------------

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMediaProperties - ERR02' 
        
        !Starts Profile for Theta / ThetaF
        call StartProfile  (ProfileID       = Me%ObjProfile,                            &
                            ObjTime         = Me%ObjTime,                               &
                            ProfileDataFile = trim(TimeSerieLocationFile),              &
                            WaterPoints2D   = Me%ExtVar%BasinPoints,                    &
                            nProperties     = Me%PropertiesNumber ,                                        &
                            PropertyList    = PropertyList,                             &
                            KUB             = Me%WorkSize%KUB,                          &
                            ClientName      = "PorousMedia",                            &
                            STAT            = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMediaProperties - ERR03' 


        deallocate (PropertyList)
        
    end subroutine ConstructProfileOutput

    !---------------------------------------------------------------------------

    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),    pointer        :: NewProperty

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------


        call GetData(NewProperty%TimeSerie,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'TIME_SERIE',                                        &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR01'
        

        call GetData(NewProperty%BoxTimeSerie,                                           &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'BOX_TIME_SERIE',                                    &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR02'


        call GetData(NewProperty%BoxTimeSerie2D,                                           &
                     Me%ObjEnterData, iflag,                                               &
                     Keyword      = 'BOX_TIME_SERIE2D',                                    &
                     Default      = .false.,                                               &
                     SearchType   = FromBlock,                                             &
                     ClientModule = 'ModulePorousMediaProperties',                         &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR03'



        call GetData(NewProperty%OutputHDF,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'OUTPUT_HDF',                                        &
                     ClientModule = 'ModulePorousMediaProperties',                       &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModulePorousMediaProperties - ERR04'
        
    end subroutine Construct_PropertyOutPut
   
   !---------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        integer                                             :: nProperties
        integer                                             :: STAT_CALL
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: n
        !Begin------------------------------------------------------------------
        
        !Counts the number of Properties which has timeserie option set to true (2x for inf col concentration)
        PropertyX => Me%FirstProperty
        nProperties = 0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                nProperties = nProperties + 2
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates PropertyList
        allocate(PropertyList(nProperties))
        
        !Property names
        n=1
        PropertyX  => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                PropertyList(n)  = trim(PropertyX%ID%Name)
                n=n+1
            endif
            PropertyX=>PropertyX%Next
        enddo

        !Property names for infil column
        PropertyX  => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                PropertyList(n)  = trim(PropertyX%ID%Name) //'_in_InfilColumn'
                n=n+1
            endif
            PropertyX=>PropertyX%Next
        enddo

        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMediaProperties',                      &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMediaProperties - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerie_ON = .true.
        else
            Me%OutPut%TimeSerie_ON = .false.
        endif
        
        !Get water points
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR03.7'


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srp",                                        &
                            WaterPoints3D = Me%ExtVar%WaterPoints3D,                    &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMediaProperties - ERR02' 

        !Unget
        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR085'

        !Deallocates PropertyList
        deallocate(PropertyList)
       
    end subroutine ConstructTimeSerie

    !--------------------------------------------------------------------------
    
    subroutine ConstructAsciiFile            

        !Local-----------------------------------------------------------------
        integer               :: status
        integer               :: STAT_CALL                
        integer               :: Counter
        character(LEN=4)      :: Number

        call UnitsManager(Me%Files%AsciiUnit, OPEN_FILE, STAT = status) 
        if (status /= SUCCESS_) stop "ConstructAsciiOutPut - ModulePorousMediaProperties - ERR01"

        Counter  = 1
do1:     do
            Number = '    '
            write(Number, fmt='(i4)')Counter
            open(UNIT   = Me%Files%AsciiUnit,                                      &
                 FILE   = '..\res\PMP_ADCoefs_'//trim(adjustl(Number))//'.log', &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1
        
        if (.not. Me%AdvDiff_ComputeTransport3D) then
            write (Me%Files%AsciiUnit, FMT=*) 'YY     MM   DD   HH   MM     SS     i   j   k    Coefa   Coefb  Coefc  CoefD '
        else
            write (Me%Files%AsciiUnit, FMT=*) 'YY     MM   DD   HH   MM     SS     i   j   k    cofA_W cofB cofC_W cofD_W cofA_U cofC_U cofA_V cofC_V '
        endif
    
    end subroutine ConstructAsciiFile
    
    !--------------------------------------------------------------------------
    

!    subroutine Construct_InitialFields
!
!        !External-----------------------------------------------------------------
!        integer                                             :: STAT_CALL
!        integer                                             :: ClientNumber
!        logical                                             :: BlockFound
!        type (T_PropertyID)                                 :: SoilDryDensityID
!        type (T_PropertyID)                                 :: SalinityID
!        type (T_PropertyID)                                 :: pHID
!        type (T_PropertyID)                                 :: IonicStrengthID
!        type (T_PropertyID)                                 :: PhosphorusAdsortionIndexID
!
!        
!        !Begin-----------------------------------------------------------------
!
!
!        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR00'
!
!        !Constructs Soil Dry Density
!        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
!                                    ClientNumber    = ClientNumber,                         &
!                                    block_begin     = '<beginsoildrydensity>',              &
!                                    block_end       = '<endsoildrydensity>',                &
!                                    BlockFound      = BlockFound,                           &   
!                                    STAT            = STAT_CALL)
!        if (STAT_CALL == SUCCESS_) then
!            if (.not. BlockFound) then
!                write(*,*)'Missing Block <beginsoildrydensity> / <endsoildrydensity>'
!                stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR10'
!            endif
!            
!            call ConstructFillMatrix  ( PropertyID           = SoilDryDensityID,            &
!                                        EnterDataID          = Me%ObjEnterData,             &
!                                        TimeID               = Me%ObjTime,                  &
!                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
!                                        ExtractType          = FromBlock,                   &
!                                        PointsToFill3D       = Me%ExtVar%OpenPoints3D,      &
!                                        Matrix2D             = Me%SoilDryDensity,           &
!                                        TypeZUV              = TypeZ_,                      &
!                                        STAT                 = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR20'
!            
!            call KillFillMatrix       (SoilDryDensityID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR30'
!
!        else
!            stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR40'
!        endif
!
!
!        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR50'
!
!        !Constructs Salinity
!        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
!                                    ClientNumber    = ClientNumber,                         &
!                                    block_begin     = '<beginsalinity>',                    &
!                                    block_end       = '<endsalinity>',                      &
!                                    BlockFound      = BlockFound,                           &   
!                                    STAT            = STAT_CALL)
!        if (STAT_CALL == SUCCESS_) then
!            if (.not. BlockFound) then
!                write(*,*)'Missing Block <beginsalinity> / <endsalinity>'
!                stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR60'
!            endif
!            
!            call ConstructFillMatrix  ( PropertyID           = SalinityID,                  &
!                                        EnterDataID          = Me%ObjEnterData,             &
!                                        TimeID               = Me%ObjTime,                  &
!                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
!                                        ExtractType          = FromBlock,                   &
!                                        PointsToFill3D       = Me%ExtVar%OpenPoints3D,      &
!                                        Matrix2D             = Me%Salinity,                 &
!                                        TypeZUV              = TypeZ_,                      &
!                                        STAT                 = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR70'
!            
!            call KillFillMatrix       (SalinityID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR80'
!
!        else
!            stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR90'
!        endif
!
!
!        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR100'
!
!        !Constructs pH
!        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
!                                    ClientNumber    = ClientNumber,                         &
!                                    block_begin     = '<beginph>',                          &
!                                    block_end       = '<endph>',                            &
!                                    BlockFound      = BlockFound,                           &   
!                                    STAT            = STAT_CALL)
!        if (STAT_CALL == SUCCESS_) then
!            if (.not. BlockFound) then
!                write(*,*)'Missing Block <beginph> / <endph>'
!                stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR110'
!            endif
!            
!            call ConstructFillMatrix  ( PropertyID           = pHID,                        &
!                                        EnterDataID          = Me%ObjEnterData,             &
!                                        TimeID               = Me%ObjTime,                  &
!                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
!                                        ExtractType          = FromBlock,                   &
!                                        PointsToFill3D       = Me%ExtVar%OpenPoints3D,      &
!                                        Matrix2D             = Me%pH,                       &
!                                        TypeZUV              = TypeZ_,                      &
!                                        STAT                 = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR120'
!            
!            call KillFillMatrix       (pHID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR130'
!
!        else
!            stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR140'
!        endif
!
!
!
!        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR100'
!
!        !Constructs Ionic Strength
!        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
!                                    ClientNumber    = ClientNumber,                         &
!                                    block_begin     = '<beginionicstrength>',               &
!                                    block_end       = '<endionicstrength>',                 &
!                                    BlockFound      = BlockFound,                           &   
!                                    STAT            = STAT_CALL)
!        if (STAT_CALL == SUCCESS_) then
!            if (.not. BlockFound) then
!                write(*,*)'Missing Block <beginionicstrength> / <endionicstrength>'
!                stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR110'
!            endif
!            
!            call ConstructFillMatrix  ( PropertyID           = IonicStrengthID,             &
!                                        EnterDataID          = Me%ObjEnterData,             &
!                                        TimeID               = Me%ObjTime,                  &
!                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
!                                        ExtractType          = FromBlock,                   &
!                                        PointsToFill3D       = Me%ExtVar%OpenPoints3D,      &
!                                        Matrix2D             = Me%IonicStrength,            &
!                                        TypeZUV              = TypeZ_,                      &
!                                        STAT                 = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR120'
!            
!            call KillFillMatrix       (IonicStrengthID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR130'
!
!        else
!            stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR140'
!        endif
!
!
!        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'InitialFields - ModulePorousMedia - ERR150'
!
!        !Constructs Phosphorus Adsortion Index
!        call ExtractBlockFromBuffer(Me%ObjEnterData,                                        &
!                                    ClientNumber    = ClientNumber,                         &
!                                    block_begin     = '<beginphosphorusadsortionindex>',    &
!                                    block_end       = '<endphosphorusadsortionindex>',      &
!                                    BlockFound      = BlockFound,                           &   
!                                    STAT            = STAT_CALL)
!        if (STAT_CALL == SUCCESS_) then
!            if (.not. BlockFound) then
!                write(*,*)'Missing Block <beginphosphorusadsortionindex> / <endphosphorusadsortionindex>'
!                stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR160'
!            endif
!            
!            call ConstructFillMatrix  ( PropertyID           = PhosphorusAdsortionIndexID,  &
!                                        EnterDataID          = Me%ObjEnterData,             &
!                                        TimeID               = Me%ObjTime,                  &
!                                        HorizontalGridID     = Me%ObjHorizontalGrid,        &
!                                        ExtractType          = FromBlock,                   &
!                                        PointsToFill3D       = Me%ExtVar%OpenPoints3D,      &
!                                        Matrix2D             = Me%PhosphorusAdsortionIndex, &
!                                        TypeZUV              = TypeZ_,                      &
!                                        STAT                 = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR170'
!            
!            call KillFillMatrix       (PhosphorusAdsortionIndexID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR180'
!
!        else
!            stop 'Construct_InitialFields - ModulePorousMediaProperties - ERR190'
!        endif
!
!
!    end subroutine Construct_InitialFields

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
                               CurrentTime = Me%ExtVar%BeginTime,            &
                               EndTime     = Me%ExtVar%EndTime,              &
                               keyword     = 'OUTPUT_TIME',                  &
                               SearchType  = FromFile,                       &
                               OutPutsTime = Me%OutPut%OutTime,              &
                               OutPutsOn   = Me%OutPut%HDF_ON,               &
                               OutPutsNumber = Me%OutPut%Number,             &
                               STAT        = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                       &
                stop 'ConstructHDF - ModulePorousMediaProperties - ERR01' 

            if (Me%OutPut%HDF_ON) then

                Me%OutPut%NextOutPut = 1

                call Open_HDF5_OutPut_File

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructHDF - ModulePorousMediaProperties - ERR02'
            endif 

        endif

    end subroutine ConstructHDF

  


    !--------------------------------------------------------------------------

     subroutine Open_HDF5_OutPut_File        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB,KLB,KUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        !Begin-----------------------------------------------------------------

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMediaProperties - ERR01'

      
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMediaProperties - ERR010'


        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMediaProperties - ERR020'
        
        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR030'
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR040'  

        
        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMediaProperties - ERR050'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMediaProperties - ERR060'


        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMediaProperties - ERR070'       


        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR80'  

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR90'


    end subroutine Open_HDF5_OutPut_File   
   
    !--------------------------------------------------------------------------

   
    subroutine ReadOldConcBoundariesHDF(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property), pointer                   :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        character (Len=StringLength)                :: PropertyName
        logical                                     :: EXIST
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: ObjHDF5
        integer                                     :: HDF5_READ
        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        KLB = Me%Size%KLB 
        KUB = Me%Size%KUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

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
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR01'


            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Concentration(:,:,:) = FillValueReal


            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB,                     &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Concentration/"//NewProperty%ID%Name,        &
                                 NewProperty%ID%Name,                                    &
                                 Array3D = NewProperty%Concentration,                    &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR03'


            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR06'

        else
            
            write(*,*)
            stop 'ReadOldConcBoundariesHDF - ModulePorousMediaProperties - ERR07'

        end if cd0

    end subroutine ReadOldConcBoundariesHDF


    !--------------------------------------------------------------------------


    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer               :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: Counter
        integer                                 :: i,j,k
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        logical                                 :: StopSubroutine = .false.
        integer                                 :: UnitAux, kaux
        real                                    :: Aux, Sum
        integer                                 :: STAT_CALL

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
            
            if (Me%ExtVar%WaterPoints3D(i, j, k) == WaterPoint) then
                               
                if (NewProperty%Concentration(i, j, k) < 0.) then

                    StopSubroutine = .true.
                    Aux            = -1
                    kaux           = k + 1

                    do while (Aux < 0)

                        Aux = NewProperty%Concentration(i, j, kaux) 

                        kaux = kaux + 1

                        if (kaux > KUB)  then 

                            Counter = 0
                            Sum     = 0
                            
                            if (NewProperty%Concentration(i-1, j, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i-1, j, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i-1, j, k)
                            
                            endif

                        
                            if (NewProperty%Concentration(i+1, j, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i+1, j, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i+1, j, k)
                            
                            endif

                            if (NewProperty%Concentration(i, j-1, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i, j-1, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i, j-1, k)
                            
                            endif

                            if (NewProperty%Concentration(i, j+1, k) >= 0 .and.          &
                                Me%ExtVar%WaterPoints3D(i, j+1, k) == WaterPoint) then

                                Counter = Counter + 1
                                Sum     = Sum + NewProperty%Concentration(i, j+1, k)
                            
                            endif
                                  

                            if (Counter > 0) then                                        

                                Aux = Sum / real(Counter)
                                exit

                            else

                                stop 'Subroutine CheckFieldConsistence; PorousMediaProperties. ERR01.'

                            endif

                        endif

                    enddo

                    NewProperty%Concentration(i, j, k) = Aux

                endif

            else

                NewProperty%Concentration(i, j, k) = FillValueReal

            endif

        enddo
        enddo
        enddo

        if (StopSubroutine) then                                                   
            
            call UnitsManager(UnitAux, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - PorousMediaProperties - ERR02' 

            open(UnitAux, FILE = trim(NewProperty%ID%name)//'.new',                 &
                 FORM = 'FORMATTED', STATUS = 'UNKNOWN', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - PorousMediaProperties - ERR03' 

            write(UnitAux,*) '<ConcentrationBegin>'
           
            do I = ILB, IUB
            do J = JLB, JUB
            do K = KLB, KUB
            
                write(UnitAux,*) NewProperty%Concentration(i, j, k)

            enddo
            enddo
            enddo

            write(UnitAux,*) '<ConcentrationEnd>'

            call UnitsManager(UnitAux, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - PorousMediaProperties - ERR04' 

            write(*,*) 'A new concentration file was created for property: ', trim(NewProperty%ID%Name)
            write(*,*) 'Run again with this new file ', trim(NewProperty%ID%name)//'.new'
            stop 'CheckFieldConsistence - PorousMediaProperties - ERR05'  

        endif

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------


    subroutine CoupleSoilQuality        

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer, pointer, dimension(:)                      :: SoilQualityPropertyList
        integer                                             :: STAT_CALL
        real                                                :: SoilQualityDT
        integer                                             :: nProp = 0 

        !Begin------------------------------------------------------------------

        !Counts the number of Properties which has WaterQuality option set to true
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilQuality) then
                nProp = nProp + 1
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates Array to hold IDs
        allocate (SoilQualityPropertyList(1:nProp))

        !Fills Array
        PropertyX => Me%FirstProperty
        nProp = 0
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilQuality) then
                nProp = nProp + 1
                SoilQualityPropertyList(nProp) = PropertyX%ID%IDNumber
            endif
            PropertyX => PropertyX%Next
        enddo

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModulePorousMediaProperties - ERR01'

        !Start Interface
        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = SedimentQualityModel,          &
                                DT                  = SoilQualityDT,                 &
                                PropertiesList      = SoilQualityPropertyList,       &
                                WaterPoints3D       = Me%ExtVar%WaterPoints3D,       &
                                Size3D              = Me%WorkSize,                   &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleSoilQuality - ModulePorousMediaProperties - ERR02'


        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModulePorousMediaProperties - ERR03'


        deallocate (SoilQualityPropertyList)

        Me%Coupled%SoilQuality_DT          = SoilQualityDT 
        Me%Coupled%SoilQuality_NextCompute = Me%ExtVar%Now    

        nullify (Me%DissolvedToParticulate3D)
        allocate(Me%DissolvedToParticulate3D(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB, &
                 Me%WorkSize%KLB:Me%WorkSize%KUB))
        Me%DissolvedToParticulate3D(:,:,:) = null_real

        Me%ResidualTime = 0.
    
    end subroutine CoupleSoilQuality

    !--------------------------------------------------------------------------


#ifdef _PHREEQC_
    !--------------------------------------------------------------------------
    subroutine CoupleSoilChemistry        

        !Local-----------------------------------------------------------------
        type(T_Property), pointer                           :: PropertyX
        integer, pointer, dimension(:)                      :: SoilChemistryPropertyList
        integer                                             :: STAT_CALL
        real                                                :: SoilChemistryDT
        integer                                             :: nProp = 0 

        !Begin------------------------------------------------------------------

        !Counts the number of Properties which has SoilChemistry option set to true
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilChemistry) then
                nProp = nProp + 1
            endif
            PropertyX => PropertyX%Next
        enddo

        !Allocates Array to hold IDs
        allocate (SoilChemistryPropertyList(1:nProp))

        !Fills Array
        PropertyX => Me%FirstProperty
        nProp = 0
        do while (associated(PropertyX))
            if (PropertyX%Evolution%SoilChemistry) then
                nProp = nProp + 1
                SoilChemistryPropertyList(nProp) = PropertyX%ID%IDNumber
            endif
            PropertyX => PropertyX%Next
        enddo

        !Question: What does this function? Is it necessary to SoilChemistry process or Interface?
        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR01'

        !Start Interface
        call ConstructInterface(InterfaceID         = Me%ObjInterfaceSoilChemistry,  &
                                TimeID              = Me%ObjTime,                    &
                                SinksSourcesModel   = PhreeqCModel,                  &
                                DT                  = SoilChemistryDT,               &
                                PropertiesList      = SoilChemistryPropertyList,     &
                                WaterPoints3D       = Me%ExtVar%WaterPoints3D,       &
                                Size3D              = Me%WorkSize,                   &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                   &
            stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR02'

        !Question: What does this function? 
        call UnGetMap (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModulePorousMediaProperties - ERR03'

        deallocate (SoilChemistryPropertyList)

        Me%Coupled%SoilChemistry_DT          = SoilChemistryDT 
        Me%Coupled%SoilChemistry_NextCompute = Me%ExtVar%Now    

            
    end subroutine CoupleSoilChemistry
    !--------------------------------------------------------------------------
#endif


    !--------------------------------------------------------------------------
    subroutine Search_Property(PropertyX, PropertyXID, STAT)

        !Arguments-------------------------------------------------------------
        type(T_Property),           pointer             :: PropertyX
        integer,                    intent (IN)         :: PropertyXID
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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
!
    !---------------------------------------------------------------------------

    subroutine GetPMPCoupled(PorousMediaPropertiesID, &
                             SoilQuality,             &
#ifdef _PHREEQC_
                             SoilChemistry,           &
#endif                             
                             STAT) 

        !Arguments-------------------------------------------------------------
        integer                        :: PorousMediaPropertiesID
        integer, optional, intent(OUT) :: STAT
        logical, optional, intent(OUT) :: SoilQuality        
#ifdef _PHREEQC_
        logical, optional, intent(OUT) :: SoilChemistry
#endif

        !External--------------------------------------------------------------

        integer :: ready_              

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(SoilQuality   )) SoilQuality    = Me%Coupled%SoilQuality
            
#ifdef _PHREEQC_            
            if (present(SoilChemistry )) SoilChemistry  = Me%Coupled%SoilChemistry
#endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetPMPCoupled
    !--------------------------------------------------------------------------
    
    subroutine GetPMPConcentration(PorousMediaPropertiesID, ConcentrationX, PropertyXIDNumber, &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        real, pointer, dimension(:,:,:)             :: ConcentrationX
        character(LEN = *), optional, intent(OUT)   :: PropertyXUnits
        integer,                      intent(IN )   :: PropertyXIDNumber
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_CALL              
        type(T_Property), pointer                   :: PropertyX
        integer                                     :: UnitsSize
        integer                                     :: STAT_    

        !------------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mPOROUSMEDIAPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                ConcentrationX => PropertyX%concentration

                if (present(PropertyXUnits)) then 
                   UnitsSize      = LEN (PropertyXUnits)
                   PropertyXUnits = PropertyX%ID%Units(1:UnitsSize)
                end if

                STAT_ = SUCCESS_
            else
                STAT_ = STAT_CALL
            end if
        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine GetPMPConcentration

    !--------------------------------------------------------------------------------

    subroutine SetWindVelocity (PorousMediaPropertiesID, WindModulus, STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        real, dimension(:,:), pointer               :: WindModulus

        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            

            Me%ExtVar%WindVelocity2D   => WindModulus


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetWindVelocity 

    !---------------------------------------------------------------------------

 
    !--------------------------------------------------------------------------    

    subroutine SetVegetationPMProperties(PorousMediaPropertiesID,                  &
                                         SoilFluxesActive,                         &
                                         GrazingBiomass,                           &
                                         GrazingNitrogen,                          &
                                         GrazingPhosphorus,                        &
                                         ManagementAerialBiomass,                  &
                                         ManagementNitrogen,                       &
                                         ManagementPhosphorus,                     &
                                         ManagementRootBiomass,                    &
                                         DormancyBiomass,                          &
                                         DormancyNitrogen,                         &
                                         DormancyPhosphorus,                       &
                                         FertilNitrateSurface,                     &
                                         FertilNitrateSubSurface,                  &
                                         FertilAmmoniaSurface,                     &
                                         FertilAmmoniaSubSurface,                  &
                                         FertilOrganicNSurface,                    &
                                         FertilOrganicNSubSurface,                 &
                                         FertilOrganicPSurface,                    &
                                         FertilOrganicPSubSurface,                 &
                                         FertilMineralPSurface,                    &
                                         FertilMineralPSubSurface,                 &
                                         NitrogenUptake,                           &
                                         PhosphorusUptake,                         &
                                         Grazing,                                  &
                                         Management,                               &
                                         Dormancy,                                 &
                                         Fertilization,                            &
                                         NutrientFluxesWithSoil,                   &
                                         RootDepth,                                &
                                         ModelNitrogen,                            &
                                         ModelPhosphorus,                          &
                                         GrowthModel,                              &
                                         CoupledVegetation,                        &
                                         NitrogenFraction,                         &
                                         PhosphorusFraction,                       &
                                         VegetationDT,                             &
                                         STAT)
                                  
        !Arguments--------------------------------------------------------------
        integer                                     :: PorousMediaPropertiesID
        logical, dimension(:,:), pointer, optional  :: SoilFluxesActive
        real, dimension(:,:), pointer, optional     :: GrazingBiomass
        real, dimension(:,:), pointer, optional     :: GrazingNitrogen
        real, dimension(:,:), pointer, optional     :: GrazingPhosphorus
        real, dimension(:,:), pointer, optional     :: ManagementAerialBiomass
        real, dimension(:,:), pointer, optional     :: ManagementNitrogen
        real, dimension(:,:), pointer, optional     :: ManagementPhosphorus
        real, dimension(:,:), pointer, optional     :: ManagementRootBiomass
        real, dimension(:,:), pointer, optional     :: DormancyBiomass
        real, dimension(:,:), pointer, optional     :: DormancyNitrogen
        real, dimension(:,:), pointer, optional     :: DormancyPhosphorus
        real, dimension(:,:), pointer, optional     :: FertilNitrateSurface
        real, dimension(:,:), pointer, optional     :: FertilNitrateSubSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSurface
        real, dimension(:,:), pointer, optional     :: FertilAmmoniaSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicNSubSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSurface
        real, dimension(:,:), pointer, optional     :: FertilOrganicPSubSurface
        real, dimension(:,:), pointer, optional     :: FertilMineralPSurface
        real, dimension(:,:), pointer, optional     :: FertilMineralPSubSurface
        real, dimension(:,:), pointer, optional     :: NitrogenFraction
        real, dimension(:,:), pointer, optional     :: PhosphorusFraction
        real, dimension(:,:), pointer, optional     :: RootDepth
        real, dimension(:,:,:),pointer,optional     :: NitrogenUptake
        real, dimension(:,:,:),pointer,optional     :: PhosphorusUptake
        logical,  optional                          :: Grazing
        logical,  optional                          :: Management
        logical,  optional                          :: Dormancy
        logical,  optional                          :: Fertilization
        logical,  optional                          :: NutrientFluxesWithSoil
        logical,  optional                          :: ModelNitrogen
        logical,  optional                          :: ModelPhosphorus
        logical,  optional                          :: GrowthModel
        logical,  optional                          :: CoupledVegetation
        real,     optional                          :: VegetationDT
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_        
        integer                                     :: STAT_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
            
            if (present(NutrientFluxesWithSoil  )) then
                Me%ExtVar%ComputeVegInterfaceFluxes   = NutrientFluxesWithSoil
            endif

            if (present(SoilFluxesActive        )) then
                Me%ExtVar%SoilFluxesActive         => SoilFluxesActive
            endif
            if (present(Grazing)) then
                Me%ExtVar%Grazing                  = Grazing
                if (present(GrazingBiomass          )) Me%ExtVar%GrazingBiomass           => GrazingBiomass
                if (present(GrazingNitrogen         )) Me%ExtVar%GrazingNitrogen          => GrazingNitrogen
                if (present(GrazingPhosphorus       )) Me%ExtVar%GrazingPhosphorus        => GrazingPhosphorus
            endif
            if (present(Management)) then
                Me%ExtVar%Management               = Management
                if (present(ManagementAerialBiomass )) Me%ExtVar%ManagementAerialBiomass  => ManagementAerialBiomass
                if (present(ManagementNitrogen      )) Me%ExtVar%ManagementNitrogen       => ManagementNitrogen
                if (present(ManagementPhosphorus    )) Me%ExtVar%ManagementPhosphorus     => ManagementPhosphorus
                if (present(ManagementRootBiomass   )) Me%ExtVar%ManagementRootBiomass    => ManagementRootBiomass
            endif
            if (present(Dormancy)) then
                Me%ExtVar%Dormancy                 = Dormancy
                if (present(DormancyBiomass         )) Me%ExtVar%DormancyBiomass          => DormancyBiomass
                if (present(DormancyNitrogen        )) Me%ExtVar%DormancyNitrogen         => DormancyNitrogen
                if (present(DormancyPhosphorus      )) Me%ExtVar%DormancyPhosphorus       => DormancyPhosphorus
            endif
            if (present(Fertilization)) then
                Me%ExtVar%Fertilization            = Fertilization
                if (present(FertilNitrateSurface    )) Me%ExtVar%FertilNitrateSurface     => FertilNitrateSurface
                if (present(FertilNitrateSubSurface )) Me%ExtVar%FertilNitrateSubSurface  => FertilNitrateSubSurface
                if (present(FertilAmmoniaSurface    )) Me%ExtVar%FertilAmmoniaSurface     => FertilAmmoniaSurface
                if (present(FertilAmmoniaSubSurface )) Me%ExtVar%FertilAmmoniaSubSurface  => FertilAmmoniaSubSurface
                if (present(FertilOrganicNSurface   )) Me%ExtVar%FertilOrganicNSurface    => FertilOrganicNSurface
                if (present(FertilOrganicNSubSurface)) Me%ExtVar%FertilOrganicNSubSurface => FertilOrganicNSubSurface
                if (present(FertilOrganicPSurface   )) Me%ExtVar%FertilOrganicPSurface    => FertilOrganicPSurface
                if (present(FertilOrganicPSubSurface)) Me%ExtVar%FertilOrganicPSubSurface => FertilOrganicPSubSurface
                if (present(FertilMineralPSurface   )) Me%ExtVar%FertilMineralPSurface    => FertilMineralPSurface
                if (present(FertilMineralPSubSurface)) Me%ExtVar%FertilMineralPSubSurface => FertilMineralPSubSurface
            endif

            if (present(NitrogenUptake          )) Me%ExtVar%NitrogenUptake           => NitrogenUptake
            if (present(PhosphorusUptake        )) Me%ExtVar%PhosphorusUptake         => PhosphorusUptake
            if (present(RootDepth               )) Me%ExtVar%RootDepth                => RootDepth

            if (present(ModelNitrogen           )) Me%ExtVar%ModelNitrogen            =  ModelNitrogen
            if (present(ModelPhosphorus         )) Me%ExtVar%ModelPhosphorus          =  ModelPhosphorus

            if (present(GrowthModel             )) Me%ExtVar%GrowthModel              =  GrowthModel
            if (present(CoupledVegetation       )) Me%ExtVar%CoupledVegetation        =  CoupledVegetation

            if (present(VegetationDT            )) Me%ExtVar%VegetationDT             =  VegetationDT

            if (present(NitrogenFraction        )) Me%ExtVar%NitrogenFraction         => NitrogenFraction
            if (present(PhosphorusFraction      )) Me%ExtVar%PhosphorusFraction       => PhosphorusFraction


            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if

        if (present(STAT))STAT = STAT_

    end subroutine SetVegetationPMProperties 

    !---------------------------------------------------------------------------
    
    subroutine SetDNConcPMP (PorousMediaPropertiesID, PropertyID, DNConcentration, ChannelsID, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: PropertyID   
        real, dimension (:), pointer                    :: DNConcentration
        integer, dimension(:, :), pointer               :: ChannelsID
        integer                                         :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
        
            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcPMP - ModulePorousMediaProperties - ERR01'
            
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyID, STAT = STAT_)
            if (STAT_ == SUCCESS_) then
            
		        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
		        do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
                        PropertyX%ConcentrationDN (i,j) = DNConcentration (ChannelsID(i,j))
                    endif
                enddo
                enddo                

            else
                write(*,*) 'Looking for Drainage Network Property in Porous Media Properties', GetPropertyName(PropertyID)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetPMPConcDrainageNetwork - ModuleDrainageNetwork - ERR010'
            end if

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR040'               


        else
            STAT_ = ready_
        end if

        STAT = STAT_        
                     
    end subroutine SetDNConcPMP

    !---------------------------------------------------------------------------
    

    subroutine SetInfColConcPMP (PorousMediaPropertiesID, PropertyXIDNumber, ConcentrationX, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: PropertyXIDNumber   
        real, dimension (:,:), pointer                  :: ConcentrationX
        integer                                         :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetInfColConcPMP - ModulePorousMediaProperties - ERR01'
        
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_)
            if (STAT_ == SUCCESS_) then
            
		        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
		        do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                        PropertyX%ConcentrationOnInfColumn (i,j) = ConcentrationX (i,j)
                    endif
                enddo
                enddo                

            else
                write(*,*) 'Looking for Atmosphere Property in Porous Media Properties', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetInfColConcPMP - ModulePorousMediaProperties - ERR010'
            end if
           
            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetInfColConcPMP - ModulePorousMediaProperties - ERR040'               


        else
            STAT_ = ready_
        end if

        STAT = STAT_        
                     
    end subroutine SetInfColConcPMP

    !--------------------------------------------------------------------------- 
    
    subroutine GetPMPnProperties (PorousMediaPropertiesID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            nProperties       = Me%PropertiesNumber
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetPMPnProperties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetPMPPropertiesIDByIdx (PorousMediaPropertiesID, Idx, ID, PropAdvDiff, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: PorousMediaPropertiesID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT)                            :: PropAdvDiff
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(PorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstProperty
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID          = CurrProp%ID%IDNumber
            PropAdvDiff = CurrProp%Evolution%AdvectionDiffusion
            
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetPMPPropertiesIDByIdx

    !---------------------------------------------------------------------------

    subroutine UnGetPorousMediaProperties3D_I(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID, "UnGetPorousMediaProperties3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMediaProperties3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetPorousMediaProperties3D_R4(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(4), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID,  "UnGetPorousMediaProperties3D_R4")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMediaProperties3D_R4


    !--------------------------------------------------------------------------

    subroutine UnGetPorousMediaProperties3D_R8(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID,  "UnGetPorousMediaProperties3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPorousMediaProperties3D_R8


    !--------------------------------------------------------------------------

!    subroutine UnGetPorousMediaProperties3D_R8i(ObjPorousMediaPropertiesID, Array, STAT)
!
!        !Arguments-------------------------------------------------------------
!        integer                                         :: ObjPorousMediaPropertiesID
!        real(8), dimension(:, :), pointer               :: Array
!        integer, intent(OUT), optional                  :: STAT
!
!        !Local-----------------------------------------------------------------
!        integer                                         :: STAT_, ready_
!
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(ObjPorousMediaPropertiesID, ready_)
!
!        if (ready_ .EQ. READ_LOCK_ERR_) then
!
!            nullify(Array)
!            call Read_Unlock(mPorousMediaProperties_, Me%InstanceID,  "UnGetPorousMediaProperties3D_R8")
!
!
!            STAT_ = SUCCESS_
!        else               
!            STAT_ = ready_
!        end if
!
!        if (present(STAT)) STAT = STAT_
!
!    end subroutine UnGetPorousMediaProperties3D_R8i


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyPorousMediaProperties(ObjPorousMediaPropertiesID,          &
                                           STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaPropertiesID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_,STAT_CALL
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then


            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMediaProperties - ModulePorousMediaProperties - ERR02'
            
            !Actualize the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyPorousMediaProperties - ModulePorousMediaProperties - ERR03'
                      
            call ReadLockExternalVar
            
            !Eduardo Jauch
            !Actualize properties if evolution from file
            call ActualizePropertiesFromFile

            if (Me%Coupled%AdvectionDiffusion) then

                !Nutrient sources and sinks from vegetation
                call InterfaceFluxes
                
				if (Me%AdvDiff_Module == AdvDif_ModulePMP_) then !Advection and Diffusion from ModulePorousMediaProperties
                    
                    call AdvectionDiffusionProcesses_PMP
 
               	elseif (Me%AdvDiff_Module == AdvDif_ModuleAD_) then !Advection and Diffusion from ModuleAdvectionDiffusion
            
                    call AdvectionDiffusionProcesses_AD

    	        endif
            
            endif

            if (Me%Coupled%SoilQuality) then
                call SoilQualityProcesses
            endif

#ifdef _PHREEQC_
            if (Me%Coupled%SoilChemistry) then
                call SoilChemistryProcesses
            endif
#endif            

            if (Me%Coupled%MinConcentration) then
                call SetLimitsConcentration 
            endif

            if (Me%Output%Timeserie_ON) then
                call OutPut_TimeSeries
            endif

            if (Me%Output%HDF_ON) then
                call OutPut_HDF
            endif

    !       call ProfileOutput    em teste no construct e no kill

            call Actualize_Time_Evolution
        
            call ReadUnlockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyPorousMediaProperties

    !-----------------------------------------------------------------------------

   
    !This method implies that there is no diffusivity in top face (interaction with water column)
    subroutine ModifyDiffusivity_Old(PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer           :: PropertyX

        !External--------------------------------------------------------------
        integer                              :: i, j, k, CHUNK
        real, pointer, dimension(:,:,:)      :: ThetaOld, Porosity
        real                                 :: WaterContent_Face, Porosity_Face

        !Begin----------------------------------------------------------------------

        call ModifyTurbulence(PropertyX)

        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K,WaterContent_Face,Porosity_Face)

!        PropertyX => Me%FirstProperty
        ThetaOld  => Me%ExtVar%WaterContentOld
        Porosity  => Me%ExtVar%ThetaS
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k=Me%WorkSize%KLB, Me%WorkSize%KUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%OpenPoints3D(i, j, k) == OpenPoint) then
                
!                        !Compute diffusivity in faces                   
!                        Theta_Face    = min (ThetaOld(i,j,k) , ThetaOld(i,j,k-1))
!                        Porosity_Face = min (Porosity(i,j,k) , Porosity(i,j,k-1))
                WaterContent_Face    = ThetaOld(i,j,k) 
                Porosity_Face = Porosity(i,j,k) 
                
                !Vertical Diffusivity
                !m2/s = m2/s * [-]  + m2/s 
                PropertyX%Diffusivity(i,j,k) = Me%ExtVar%ComputeFacesW3D(I,J,K) *                &
                                               (PropertyX%Evolution%AdvDiff%Molecular_Diff_Coef * & 
                                               Tortuosity (WaterContent_Face,Porosity_Face) +           &
                                               PropertyX%Diff_Turbulence_V (i,j,k))
                
                !Horizontal Diffusivity is called viscosity to maintain the format of when is called in Module Advection Diffusion
                !m2/s = m2/s * [-]  + m2/s
                PropertyX%Viscosity(i,j,k)   = PropertyX%Evolution%AdvDiff%Molecular_Diff_Coef * &
                                               Tortuosity (WaterContent_Face,Porosity_Face) +           &
                                               PropertyX%Diff_Turbulence_H (i,j,k)
                
                if (Me%AdvDiff_Module == AdvDif_ModulePMP_) then ! need ViscosityU and ViscosityV
                    PropertyX%ViscosityU(i,j,k) = PropertyX%Viscosity(i,j,k) * Me%ExtVar%ComputeFacesU3D(I,J,K)
                    PropertyX%ViscosityV(i,j,k) = PropertyX%Viscosity(i,j,k) * Me%ExtVar%ComputeFacesV3D(I,J,K)
                    if (K == Me%WorkSize%KUB) then !need also diffusivity in face KUB +1
                        if (Me%ExtVar%InfiltrationColumn(i,j) .gt. 0.) then
                            PropertyX%Diffusivity(i,j,k+1) = PropertyX%Evolution%AdvDiff%Molecular_Diff_Coef * & 
                                                           Tortuosity (WaterContent_Face,Porosity_Face) +           &
                                                           PropertyX%Diff_Turbulence_V (i,j,k)   
                        else
                            PropertyX%Diffusivity(i,j,k+1) = 0.0
                        endif
                    endif                             
                    
                endif
                                                                       

            endif

        enddo
        enddo
        enddo

        !$OMP END DO
        !$OMP END PARALLEL


    end subroutine ModifyDiffusivity_Old

    !--------------------------------------------------------------------------

    subroutine ModifyTurbulence (PropertyX)

        !Arguments------------------------------------------------------------------
        type (T_Property), pointer :: PropertyX

        !Local----------------------------------------------------------------------
        integer                    :: i, j, k, CHUNK
        real                       :: VelMedU, VelMedV, VelMedW

        !Begin----------------------------------------------------------------------

        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k=Me%WorkSize%KLB, Me%WorkSize%KUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%OpenPoints3D(i, j, k) == OpenPoint) then
                
!                VelMedU = 0.5 * (Me%ExtVar%UnsatU(i,j,k) * Me%ExtVar%ComputeFacesU3D(I,J,K) +  Me%ExtVar%UnsatU(i,j+1,k) * Me%ExtVar%ComputeFacesU3D(I,J+1,K))
!                VelMedV = 0.5 * (Me%ExtVar%UnsatV(i,j,k) * Me%ExtVar%ComputeFacesV3D(I,J,K) +  Me%ExtVar%UnsatV(i+1,j,k) * Me%ExtVar%ComputeFacesV3D(I+1,J+,K))
!                VelMedW = 0.5 * (Me%ExtVar%UnsatW(i,j,k) * Me%ExtVar%ComputeFacesW3D(I,J,K) +  Me%ExtVar%UnsatW(i,j,k+1) * Me%ExtVar%ComputeFacesW3D(I,J+,K+1))
! 
!                !m2/s = m * m/s
!                PropertyX%Diff_Turbulence_V(i,j,k) = Me%Disper_Longi%Field(i,j,k) * abs(VelMedW)                                  &
!                                                    + Me%Disper_Trans%Field(i,j,k) * abs(0.5 * (VelMedU + VelMedV))
!                !m2/s = m * m/s
!                PropertyX%Diff_Turbulence_H(i,j,k) = Me%Disper_Longi%Field(i,j,k) *  abs(0.5 * (VelMedU + VelMedV)                &
!                                                    + Me%Disper_Trans%Field(i,j,k) * abs(VelMedW)
 
                
                !m2/s = m * m/s
                PropertyX%Diff_Turbulence_V(i,j,k) = Me%Disper_Longi%Field(i,j,k) *                                               &
		                                                 abs(Me%ExtVar%UnsatW(i,j,k) * Me%ExtVar%ComputeFacesW3D(I,J,K)) +        &
		                                             Me%Disper_Trans%Field(i,j,k) *                                               &
		                                                 abs((Me%ExtVar%UnsatV(i,j,k) * Me%ExtVar%ComputeFacesV3D(I,J,K)) +       &
		                                                     (Me%ExtVar%UnsatU(i,j,k) * Me%ExtVar%ComputeFacesU3D(I,J,K)) / 2.)
                !m2/s = m * m/s
                PropertyX%Diff_Turbulence_H(i,j,k) = Me%Disper_Longi%Field(i,j,k) *                                               &
		                                                 abs((Me%ExtVar%UnsatV(i,j,k) * Me%ExtVar%ComputeFacesV3D(I,J,K)) + &
		                                                     (Me%ExtVar%UnsatU(i,j,k) * Me%ExtVar%ComputeFacesU3D(I,J,K)) / 2.) +  &
		                                             Me%Disper_Trans%Field(i,j,k) *                                               &
		                                                 abs(Me%ExtVar%UnsatW(i,j,k) * Me%ExtVar%ComputeFacesW3D(I,J,K))
               
            endif

        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
        
        !---------------------------------------------------------------------------

    end subroutine ModifyTurbulence

    !--------------------------------------------------------------------------
    
	subroutine AdvectionDiffusionTopBoundary

        !Local--------------------------------------------------------------
        type (T_Property), pointer         :: PropertyX
        integer                            :: i, j, k, CHUNK
        real(8), dimension(:,:,:), pointer :: WaterVolume
        real(8)                            :: InfVolume, MassOnFluxW

        !Begin----------------------------------------------------------------------

        !CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !!$OMP PARALLEL PRIVATE(I,J,K,InfVolume)

        PropertyX => Me%FirstProperty
        k = Me%WorkSize%KUB

        if (Me%AdvDiff_Explicit) then
            WaterVolume => Me%WaterVolumeCorr
        else
            WaterVolume => Me%WaterVolumeOld
        endif

do1:    do while (associated(PropertyX))

			if (PropertyX%Evolution%AdvectionDiffusion) then

                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
		        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
		        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

		        	if (Me%ExtVar%OpenPoints3D(i, j, k) == OpenPoint) then

                         InfVolume = abs(Me%ExtVar%FluxW(i, j, k+1)) * Me%ExtVar%DT
    					
    					if (Me%ExtVar%FluxW(i, j, k+1) < 0.) then						
                            
                            !g = g/m3 * m3
!                            MassOnFluxW = PropertyX%ConcentrationOnWaterColumn(i, j) * InfVolume
                            MassOnFluxW = PropertyX%ConcentrationOninfColumn(i, j) * InfVolume
                            
                            !g/m3  = (g/m3 * m3 + g) /((+ m3 + m3))
		        			PropertyX%Concentration(i, j, k) = ((PropertyX%Concentration(i, j, k) * WaterVolume(i, j, k))            &
		        			                                    + MassOnFluxW) / ((WaterVolume(i, j, k) + InfVolume))
                        else
                            !g = g/m3 * m3
                            MassOnFluxW = PropertyX%Concentration(i, j, Me%WorkSize%KUB) * InfVolume
                            
                            !g/m3  = (g/m3 * m3 + g) /((+ m3 + m3))
		        			PropertyX%Concentration(i, j, k) = ((PropertyX%Concentration(i, j, k) * WaterVolume(i, j, k))            &
		        			                                    - MassOnFluxW) / ((WaterVolume(i, j, k) - InfVolume))
                        		        			                                    
						endif

		        	endif

		        enddo
		        enddo
                !!$OMP END DO
                
			endif

			PropertyX => PropertyX%Next

	    enddo do1
        
        !!$OMP END PARALLEL

	end subroutine AdvectionDiffusionTopBoundary

    !--------------------------------------------------------------------------

    subroutine ComputeVolumes

        !Local-----------------------------------------------------------------
        integer :: i, j, k, CHUNK        

        !----------------------------------------------------------------------
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        !$OMP PARALLEL PRIVATE(I,J,K)
        
        !Compute volumes and correct top volume taking FluxW(KUB+1) because it would be interpreted by module advection diffusion
        !as an additional water flux with the conc of C(i,j,k)
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(i,j,k) == WaterPoint) then             
                 
                Me%WaterVolumeOld(i,j,k)     = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)
 !               Me%WaterVolumeOldCorr(i,j,k) = Me%WaterVolumeCorr(i,j,k)
                Me%WaterVolume(i,j,k)        = Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)
                Me%WaterVolumeCorr(i,j,k)    = Me%WaterVolume(i,j,k)
                if (k == Me%WorkSize%KUB) then
                    !m3 = m3 - m3/s * s
                    Me%WaterVolumeCorr(i,j,k) = Me%WaterVolumeCorr(i,j,k) + Me%ExtVar%FluxW(i,j,k+1)  * Me%ExtVar%DT
                endif

            endif
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL 
               
        !Correct fluxw - take FluxW(KUB+1) because it would be interpreted by module advection diffusion
        !as an additional water flux with the conc of C(i,j,k)
        k = Me%WorkSize%KUB
        call SetMatrixValue (Me%FluxWCorr, Me%Size, Me%ExtVar%FluxW)
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%OpenPoints3D(i,j,k) == OpenPoint) then  
                Me%FluxWCorr(i,j,k+1) = 0.0
            endif
        enddo
        enddo        

   
    end subroutine ComputeVolumes
    
   !----------------------------------------------------------------------


    subroutine AdvectionDiffusionProcesses_AD
    
        !External--------------------------------------------------------------
        integer                             :: STAT_CALL    
        real(8), pointer, dimension(:,:,:)  :: AdvFluxX
        real(8), pointer, dimension(:,:,:)  :: AdvFluxY
        real(8), pointer, dimension(:,:,:)  :: AdvFluxZ
        real(8), pointer, dimension(:,:,:)  :: DifFluxX
        real(8), pointer, dimension(:,:,:)  :: DifFluxY
        real(8), pointer, dimension(:,:,:)  :: DifFluxZ

        !Local-----------------------------------------------------------------
        type(T_Property), pointer           :: Property
        type (T_Time)                       :: Actual
        real                                :: ImpExp_AdvXX, ImpExp_AdvYY           
        integer                             :: i, j, k
        real                                :: AdvectionV_imp_exp  
        real                                :: DiffusionV_imp_exp  
        real                                :: AdvectionH_imp_exp  
        real(8)                             :: f
        real(8), dimension(:,:,:), pointer  :: FluxW

        !----------------------------------------------------------------------      
        Actual = Me%ExtVar%Now

        call ComputeVolumes
    
        Property => Me%FirstProperty

do1 :   do while (associated(Property))

cd1 :       if (Property%Evolution%AdvectionDiffusion) then

                call ModifyDiffusivity_Old(Property)

                if (Me%AdvDiff_Explicit) then

                    AdvectionV_imp_exp = ExplicitScheme
                    DiffusionV_imp_exp = ExplicitScheme
                    AdvectionH_imp_exp = ExplicitScheme
                
                else

                    AdvectionV_imp_exp = ImplicitScheme
                    DiffusionV_imp_exp = ImplicitScheme
                    AdvectionH_imp_exp = ImplicitScheme

                endif
                    
                if(AdvectionH_imp_exp == ImplicitScheme) then

                    if(Property%Evolution%AdvDiff%ImplicitH_Direction == DirectionX)then
                                                   
                        !Direction X implicit
                        ImpExp_AdvXX = ImplicitScheme 
                        ImpExp_AdvYY = ExplicitScheme 

                        Property%Evolution%AdvDiff%ImplicitH_Direction = DirectionY

                    else 
                    
                        !Direction Y implicit
                        ImpExp_AdvXX = ExplicitScheme 
                        ImpExp_AdvYY = ImplicitScheme 

                        Property%Evolution%AdvDiff%ImplicitH_Direction = DirectionX

                    endif 
            
                else ! Horizontal Advection Explicit

                    ImpExp_AdvXX = ExplicitScheme 
                    ImpExp_AdvYY = ExplicitScheme 

                endif

				if (.not. Me%AdvDiff_Explicit) then

	            	call AdvectionDiffusionTopBoundary  
                
                endif

                call AdvectionDiffusion(Me%ObjAdvectionDiffusion,                                            &
				                        PROP                = Property%Concentration,                        &
				                        schmidt_H           = Property%Evolution%AdvDiff%SchmidtNumberH,     &
				                        SchmidtCoef_V       = Property%Evolution%AdvDiff%SchmidtCoefV,       &
				                        SchmidtBackground_V = Property%Evolution%AdvDiff%SchmidtBackgroundV, &
				                        AdvMethodH          = Property%Evolution%AdvDiff%AdvMethodH,         &
				                        TVDLimitationH      = Property%Evolution%AdvDiff%TVDLimitationH,     &
				                        AdvMethodV          = Property%Evolution%AdvDiff%AdvMethodV,         &
				                        TVDLimitationV      = Property%Evolution%AdvDiff%TVDLimitationV,     &
				                        Upwind2H            = Property%Evolution%AdvDiff%Upwind2H,           &
				                        Upwind2V            = Property%Evolution%AdvDiff%Upwind2V,           &
				                        VolumeRelMax        = Property%Evolution%AdvDiff%VolumeRelMax,       &
!				                        DTProp              = Property%Evolution%DTInterval,                 &
                                        DTProp              = Me%ExtVar%DT,                                  &
				                        ImpExp_AdvV         = AdvectionV_imp_exp,                            &
				                        ImpExp_DifV         = DiffusionV_imp_exp,                            &
				                        ImpExp_AdvXX        = ImpExp_AdvXX,                                  &
				                        ImpExp_AdvYY        = ImpExp_AdvYY,                                  &
				                        ImpExp_DifH         = Property%Evolution%AdvDiff%DiffusionH_imp_exp, &
				                        NullDif             = Property%Evolution%AdvDiff%NullDif,            &
				                        Wflux_X             = Me%ExtVar%FluxU,                               &
				                        Wflux_Y             = Me%ExtVar%FluxV,                               &
!				                        Wflux_Z             = Me%ExtVar%FluxW,                               &
				                        Wflux_Z             = Me%FluxWCorr,                                  &
				                        VolumeZOld          = Me%WaterVolumeOld,                             &
!				                        VolumeZ             = Me%WaterVolume,                                &
!				                        VolumeZOld          = Me%WaterVolumeOldCorr,                         &
				                        VolumeZ             = Me%WaterVolumeCorr,                            &
				                        OpenPoints3D        = Me%ExtVar%OpenPoints3D,                        &
				                        LandPoints3D        = Me%ExtVar%LandPoints3D,                        &
				                        ComputeFacesU3D     = Me%ExtVar%ComputeFacesU3D,                     &
				                        ComputeFacesV3D     = Me%ExtVar%ComputeFacesV3D,                     &
				                        ComputeFacesW3D     = Me%ExtVar%ComputeFacesW3D,                     &
				                        Visc_H              = Property%Viscosity,                            &
				                        Diff_V              = Property%Diffusivity,                          &
				                        CellFluxes          = .true.,                                        &
				                        BoundaryCondition   = Property%Evolution%AdvDiff%BoundaryCondition,  &
				                        NumericStability    = Property%Evolution%AdvDiff%NumericStability,   &
				                        STAT                = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'AdvectionDiffusionProcesses_B - ModulePorousMediaProperties - ERR10'
           

				if (Me%AdvDiff_Explicit) then

	            	call AdvectionDiffusionTopBoundary  

                endif
           
            end if cd1

            Property => Property%Next

        end do do1
        nullify(Property)

        !-------------------------------------------------------------------------    
        
    end subroutine AdvectionDiffusionProcesses_AD
        
    !--------------------------------------------------------------------------
    
    subroutine ActualizePropertiesFromFile
    
        !Local--------------------------------------------------------------------        
        type (T_Property), pointer :: PropertyX 
        integer                    :: STAT_CALL   
        
        !-------------------------------------------------------------------------    

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%ID%SolutionFromFile) then
            
                call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                       Matrix3D       = PropertyX%Concentration,    &
                                       PointsToFill3D = Me%ExtVar%WaterPoints3D,    &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModulePorousMediaProperties - ERR01'
            
            endif
            
            PropertyX => PropertyX%Next
            
        enddo
    
        !-------------------------------------------------------------------------    
    
    end subroutine ActualizePropertiesFromFile
    
	!-------------------------------------------------------------------------    

    subroutine InterfaceFluxes
        !Local--------------------------------------------------------------------
        !Begin--------------------------------------------------------------------

        if (Me%ExtVar%CoupledVegetation) then
            if (Me%ExtVar%ComputeVegInterfaceFluxes) then
                call VegetationInterfaceFluxes
            endif
        endif

!        if (Me%ExtVar%CoupledDN) then
!            call DrainageNetworkInterfaceFluxes
!        endif

!        if (Me%ExtVar%CoupledRunoff) then
!            if (Me%ExtVar%ComputeRunoffInterfaceFluxes) then
!                call RunoffInterfaceFluxes
!            endif
!        endif


    end subroutine InterfaceFluxes
 
    !-----------------------------------------------------------------------------
    ! This routine solves mass sources ans sinks due to vegetation and implicitly take or add mass. 
    
    subroutine VegetationInterfaceFluxes 

        !Local--------------------------------------------------------------------
        integer                                     :: i, j, k!, CHUNK
        real                                        :: Area, RootDepth
        logical                                     :: FoundEnd
        real                                        :: BottomDepth, TopDepth
        real                                        :: GrazingNotCarbon, GrazingNitrogen, GrazingPhosphorus
        real                                        :: GrazingBiomass, GrazingCarbon
        real                                        :: DormancyNotCarbon, DormancyNitrogen, DormancyPhosphorus
        real                                        :: DormancyBiomass, DormancyCarbon
        real                                        :: ManagementNotCarbon, ManagementNitrogen, ManagementPhosphorus
        real                                        :: ManagementAerialBiomass, ManagementCarbon
        real                                        :: ManagementRootNotCarbon, ManagementRootNitrogen, ManagementRootPhosphorus
        real                                        :: ManagementRootBiomass, ManagementRootCarbon
        real                                        :: NitrogenFraction, PhosphorusFraction, RootDistribution
        real                                        :: FertilizationAmmonia, FertilizationNitrate
        real                                        :: FertilizationOrganicN, FertilizationOrganicP
        real                                        :: FertilizationMineralP
        real                                        :: NitrogenUptake, PhosphorusUptake
        real                                        :: ModelDT, VegDT, CellWaterVolume, CellSoilMass
        type (T_Property), pointer                  :: Property
        integer                                     :: STAT_CALL

        !Begin--------------------------------------------------------------------

                
        !!CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR01'  


        !!!$OMP PARALLEL PRIVATE(I,J,K)
        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
        if (Me%ExtVar%BasinPoints(i,j) == BasinPoint .and. Me%ExtVar%SoilFluxesActive(i,j)) then
            
            Area         = Me%ExtVar%Area(i,j)
            RootDepth    = Me%ExtVar%RootDepth(i,j)
            FoundEnd     = .false.
            BottomDepth  = 0.
      
do3:        do K = Me%WorkSize%KUB, Me%WorkSize%KLB, -1                
            
                
                if (FoundEnd) then
                    exit do3
                endif
                
                TopDepth    = BottomDepth
                BottomDepth = BottomDepth + Me%ExtVar%DWZ(i,j,k)
                !If found root, let compute, will exit next iteration
                if (BottomDepth .ge. RootDepth) then
                    FoundEnd = .true.
                    BottomDepth = RootDepth
                endif
                    
                if (Me%ExtVar%GrowthModel) then

                    !Fluxes only occuring in surface (aerial biomass residue from grazing, dormancy and management; 
                    !surface fertilization)
                    if (k == Me%WorkSize%KUB) then
                
                        !Grazing
                        GrazingNotCarbon  = 0.0
                        GrazingNitrogen   = 0.0
                        GrazingPhosphorus = 0.0
                        if (Me%ExtVar%Grazing) then
                
                            if (Me%ExtVar%ModelNitrogen) then
                            
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                GrazingNitrogen  = Me%ExtVar%GrazingNitrogen(i,j) * 1e9 * Area / 10000.
                            
                                GrazingNotCarbon = GrazingNotCarbon + GrazingNitrogen
                
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                GrazingPhosphorus = Me%ExtVar%GrazingPhosphorus(i,j) * 1e9 * Area / 10000.

                                GrazingNotCarbon  = GrazingNotCarbon + GrazingPhosphorus
                
                            endif                          
                
                            !      ug       = Kg/ha * 1E9mg/kg * (m2) * 1ha/10000m2                     
                            GrazingBiomass = Me%ExtVar%GrazingBiomass(i,j) * 1e9 * Area / 10000.
                
                            GrazingCarbon  = GrazingBiomass - GrazingNotCarbon

                        endif

        !                !Dormancy
                        DormancyNotCarbon  = 0.0
                        DormancyNitrogen   = 0.0
                        DormancyPhosphorus = 0.0
                        if (Me%ExtVar%Dormancy) then
                
                            if (Me%ExtVar%ModelNitrogen) then

                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                DormancyNitrogen  = Me%ExtVar%DormancyNitrogen(i,j) * 1e9 * Area / 10000.

                                DormancyNotCarbon = DormancyNotCarbon + DormancyNitrogen
                
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                DormancyPhosphorus = Me%ExtVar%DormancyPhosphorus(i,j) * 1e9 * Area / 10000.

                                DormancyNotCarbon  = DormancyNotCarbon + DormancyPhosphorus
                
                            endif                          
                
                            !      ug       = Kg/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                            DormancyBiomass = Me%ExtVar%DormancyBiomass(i,j) * 1e9 * Area / 10000.
                
                            DormancyCarbon  = DormancyBiomass - DormancyNotCarbon

                        endif


        !                !Management
                        ManagementNotCarbon  = 0.0
                        ManagementNitrogen   = 0.0
                        ManagementPhosphorus = 0.0
                        if (Me%ExtVar%Management) then
                
                            if (Me%ExtVar%ModelNitrogen) then

                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                ManagementNitrogen  = Me%ExtVar%ManagementNitrogen(i,j) * 1e9 * Area / 10000.

                                ManagementNotCarbon = ManagementNotCarbon + ManagementNitrogen
                
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                ManagementPhosphorus = Me%ExtVar%ManagementPhosphorus(i,j) * 1e9 * Area / 10000.

                                ManagementNotCarbon  = ManagementNotCarbon + ManagementPhosphorus
                
                            endif                          
                
                            ManagementAerialBiomass = Me%ExtVar%ManagementAerialBiomass(i,j) * 1e9 * Area / 10000.
                
                            !      ug       = Kg/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                            ManagementCarbon  = ManagementAerialBiomass - ManagementNotCarbon

                        endif

        !                !Fertilization in Surface
                        FertilizationAmmonia    = 0.0
                        FertilizationNitrate    = 0.0
                        FertilizationOrganicN   = 0.0
                        FertilizationOrganicP   = 0.0
                        FertilizationMineralP   = 0.0
                        if (Me%ExtVar%Fertilization) then
                

                            if (Me%ExtVar%ModelNitrogen) then
                    
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationNitrate  = Me%ExtVar%FertilNitrateSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationAmmonia  = Me%ExtVar%FertilAmmoniaSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationOrganicN = Me%ExtVar%FertilOrganicNSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationOrganicP = Me%ExtVar%FertilOrganicPSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationMineralP = Me%ExtVar%FertilMineralPSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                          
                
                        endif
                
                    !Fluxes only occuring in subsurface (fertilization in sub surface)
                    elseif (k == Me%WorkSize%KUB - 1) then

        !                !Fertilization in SubSurface
                        FertilizationAmmonia    = 0.0
                        FertilizationNitrate    = 0.0
                        FertilizationOrganicN   = 0.0
                        FertilizationOrganicP   = 0.0
                        FertilizationMineralP   = 0.0
                        if (Me%ExtVar%Fertilization) then
                

                            if (Me%ExtVar%ModelNitrogen) then
                    
                                !      gN       = KgN/ha * 1E6g/kg * (m2) * 1ha/10000m2                     
                                FertilizationNitrate  = Me%ExtVar%FertilNitrateSubSurface(i,j) * 1e6 * Area / 10000.
                                FertilizationAmmonia  = Me%ExtVar%FertilAmmoniaSubSurface(i,j) * 1e6 * Area / 10000.
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2 
                                FertilizationOrganicN = Me%ExtVar%FertilOrganicNSubSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                    
                            if (Me%ExtVar%ModelPhosphorus) then

                                !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationOrganicP = Me%ExtVar%FertilOrganicPSubSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationMineralP = Me%ExtVar%FertilMineralPSubSurface(i,j) * 1e9 * Area / 10000.
                    
                            endif                          
                
                        endif

                    endif

                    !Root death to soil (occurrs in all layers until root end )
                    if (RootDepth .gt. 0.0) then
                              !      ug       = Kg/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                        ManagementRootBiomass = Me%ExtVar%ManagementRootBiomass(i,j) * 1e9 * Area / 10000.
                        !Root distribution (based on SWAT formulation)
                        RootDistribution = (1.0 - exp(-10.0 * BottomDepth/RootDepth))/(1.0 - exp(-10.0))                   &
                                           - (1.0 - exp(-10.0 * TopDepth/RootDepth))/(1.0 - exp(-10.0))

                        ManagementRootNotCarbon  = 0.0
                        ManagementRootNitrogen   = 0.0
                        ManagementRootPhosphorus = 0.0
                        if (Me%ExtVar%ModelNitrogen) then
                            NitrogenFraction        = Me%ExtVar%NitrogenFraction (i,j)
                            ManagementRootNitrogen  = ManagementRootBiomass * NitrogenFraction * RootDistribution
                            ManagementRootNotCarbon = ManagementRootNotCarbon + ManagementRootNitrogen
                        endif

                        if (Me%ExtVar%ModelPhosphorus) then
                            PhosphorusFraction        = Me%ExtVar%PhosphorusFraction (i,j)
                            ManagementRootPhosphorus  = ManagementRootBiomass * PhosphorusFraction * RootDistribution
                            ManagementRootNotCarbon   = ManagementRootNotCarbon + ManagementRootPhosphorus
                        endif
                
                        !     mg
                        ManagementRootCarbon = ManagementRootBiomass - ManagementRootNotCarbon

                    endif
                endif


                !Plant Uptake (occurrs in all layers until root end )
                if (Me%ExtVar%ModelNitrogen) then
                    
                    !      gN       = KgN/ha * 1E6g/kg * (m2) * 1ha/10000m2                     
                    NitrogenUptake = Me%ExtVar%NitrogenUptake(i,j,k) * 1e6 * Area / 10000.
                    
                endif

                if (Me%ExtVar%ModelPhosphorus) then
                    
                    !      gP       = KgP/ha * 1E6g/kg * (m2) * 1ha/10000m2                     
                    PhosphorusUptake = Me%ExtVar%PhosphorusUptake(i,j,k) * 1e6 * Area / 10000.

                endif


                !s
                ModelDT         = Me%ExtVar%DT
                VegDT           = Me%ExtVar%VegetationDT
                !m3             = m3H20/m3cell * m3cell
                CellWaterVolume = Me%ExtVar%WaterContentOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 

                !Soil mass to compute organic and microorganisms pools
                call SearchProperty(Property, SoilDryDensity_        , .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
                !kgsoil         = kg/m3  * m3cell
                CellSoilMass    = Property%ConcentrationOld(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                
                if (Me%ExtVar%GrowthModel) then

                    ! Property Calculation
                    !!Carbon
                    call SearchProperty (Property, RefreactaryOrganicC_        , .false., STAT = STAT_CALL)    
                    if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR80'
                
                    !         ug/kgsoil            = ug/kgsoil + ug / kgsoil
                    Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingCarbon + DormancyCarbon       &
                                                        + ManagementCarbon + ManagementRootCarbon) * ModelDT / VegDT)            &
                                                        / CellSoilMass)

    !                call SearchProperty (Property, LabileOrganicC_        , .false., STAT = STAT_CALL)    
    !                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
    !
    !                Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicC)               &
    !                                                * ModelDT / VegDT) / CellSoilMass)
                endif

                !!Nitrogen
                if (Me%ExtVar%ModelNitrogen) then
                    
                    if (Me%ExtVar%GrowthModel) then
                        call SearchProperty (Property, RefreactaryOrganicN_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR100'
                    
                        !         ug/kgsoil 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingNitrogen + DormancyNitrogen       &
                                                           + ManagementNitrogen + ManagementRootNitrogen) * ModelDT / VegDT)             &
                                                           / CellSoilMass)


                        call SearchProperty (Property, PON_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR110'
                    
                        !         ug/kgsoil 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicN)                   &
                                                           * ModelDT / VegDT) / CellSoilMass)


                        call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'

                        !         g/m3                = g/m3 + g / m3H20
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationNitrate - NitrogenUptake)    &
                                                           * ModelDT / VegDT) / CellWaterVolume)


                        call SearchProperty (Property, Ammonia_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR130'

                        !         g/m3                = g/m3 + g / m3H20 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationAmmonia)                    &
                                                           * ModelDT / VegDT) / CellWaterVolume)
                    else

                        call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'

                        !         g/m3                = g/m3 + g / m3H20 
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) - (((NitrogenUptake)                          &
                                                           * ModelDT / VegDT) / CellWaterVolume)
                    endif

                endif
                
                !Phosphorus
                if (Me%ExtVar%ModelPhosphorus) then

                    if (Me%ExtVar%GrowthModel) then

                        call SearchProperty (Property, RefreactaryOrganicP_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR140'

                        !         ug/kgsoil  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((GrazingPhosphorus + DormancyPhosphorus       &
                                                           + ManagementPhosphorus + ManagementRootNitrogen) * ModelDT / VegDT)               &
                                                           / CellSoilMass)


                        call SearchProperty (Property, POP_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR150'

                        !         ug/kgsoil  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationOrganicP)                       &
                                                           * ModelDT / VegDT) / CellSoilMass)


                        call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'

                        !         g/m3                    = g/m3 + g / m3H20  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) + (((FertilizationMineralP - PhosphorusUptake)   &
                                                           * ModelDT / VegDT) / CellWaterVolume)   
                    else
                        
                        call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'

                        !         g/m3                    = g/m3 + g / m3H20  
                        Property%ConcentrationOld (i,j,k) = Property%ConcentrationOld (i,j,k) - (((PhosphorusUptake)                         &
                                                           * ModelDT / VegDT) / CellWaterVolume)                        
                     
                    endif             

                endif

            enddo do3

        endif
        
        enddo
        enddo

        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR1070'  



    end subroutine VegetationInterfaceFluxes

    !-----------------------------------------------------------------------------
    
    ! This routine solves mass sources ans sinks due to drainage network and implicitly take or add mass. 

!    subroutine DrainageNetworkInterfaceFluxes 
!        !Local--------------------------------------------------------------------
!        integer                                 :: STAT_CALL, i, j, k        
!        real, dimension(:, :), pointer          :: FlowToChannels
!        real, dimension(:, :, :), pointer       :: ThetaOld
!        integer, dimension(:, :), pointer       :: GWLayer!, ChannelsID
!        type (T_Property), pointer              :: PropertyX
!!        real, dimension (:), pointer            :: DNConcentration
!        real                                    :: WaterVolumeOld, MassFlow, OldMass, NewMass
!
!        !begin--------------------------------------------------------------------    
!
!        call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkInterfaceFluxes - ModulePorousMediaProperties - ERR01'
!
!        call GetGWFlowToChannels   (Me%ObjPorousMedia, FlowToChannels, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkInterfaceFluxes - ModulePorousMediaProperties - ERR10'
!
!        call GetGWLayer   (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkInterfaceFluxes - ModulePorousMediaProperties - ERR20'
!
!        PropertyX => Me%FirstProperty
!        
!
!        ThetaOld => Me%ExtVar%WaterContentOld
!
!        do while (associated(PropertyX))        
!        
!            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!                if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
!                    k = GWLayer(i,j)
!                    if (FlowToChannels(i,j) .gt. 0.0) then !transport to channel - looses mass
!                        
!                        !g = m3/s * s * g/m3 
!                        MassFlow = FlowToChannels(i,j) * Me%ExtVar%DT * PropertyX%ConcentrationOld(i,j,k)
!                        
!                    
!                    else ! transport to soil - gains mass
!                        
!                        !g = m3/s * s * g/m3 
!                        MassFlow = FlowToChannels(i,j) * Me%ExtVar%DT * PropertyX%ConcentrationDN(i,j)
!                        
!                    endif
!                    !m3 = - * m3
!                    WaterVolumeOld = ThetaOld(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)
!                    !g = g/m3 * m3
!                    OldMass        = PropertyX%ConcentrationOld(i,j,k) * WaterVolumeOld
!                    NewMass        = OldMass - MassFlow
!                    !g/m3 = (g) / (m3Old + (m3/s Flux * s))
!                    PropertyX%ConcentrationOld(i,j,k) = NewMass / (WaterVolumeOld + (FlowToChannels(i,j) * Me%ExtVar%DT))
!
!                endif
!            enddo
!            enddo        
!
!            PropertyX => PropertyX%Next
!
!        enddo
!
!        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR040'               
!
!        call UnGetPorousMedia (Me%ObjPorousMedia, FlowToChannels, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR050'               
!
!        call UnGetPorousMedia (Me%ObjPorousMedia, GWlayer, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR060'               
!
!        
!    end subroutine DrainageNetworkInterfaceFluxes
    
    !-----------------------------------------------------------------------------

    subroutine AdvectionDiffusionProcesses_PMP

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        
        !begin--------------------------------------------------------------------

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%Evolution%AdvectionDiffusion) then

                if (Me%AdvDiff_ComputeTransport3D) then
                    
                    call ModifyAdvectionDiffusion_3D(PropertyX)
                
                else
                    !1D vertical
                    call ModifyAdvectionDiffusion_W(PropertyX)
                
                endif

            endif


            PropertyX => PropertyX%Next

        enddo

    end subroutine AdvectionDiffusionProcesses_PMP
    
    !-----------------------------------------------------------------------------


    subroutine ModifyAdvectionDiffusion_W (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k!, CHUNK
        real(8)                                     :: Area, AdvTermB_Top, AdvTermC
        real(8)                                     :: AdvTermD, AdvTermB_Bottom
        real(8)                                     :: AdvTermA
        real(8)                                     :: DifTerm_Top, Difterm_Bottom
        real(8)                                     :: aux 
        real(8)                                     :: cofA,cofB,cofC, cofD, cofInterfaceDN
        real(8)                                     :: ConcTop, ConcInInterfaceDN
        real(8), pointer, dimension(:,:,:)          :: FluxW, FluxU, FluxV
        real   , pointer, dimension(:,:,:)          :: DWZ, DZZ, Theta, ThetaOld, Porosity
        logical                                     :: ComputeCofC, ComputeCofD
        !Begin-----------------------------------------------------------------
   
        !!CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

        CurrProperty => PropertyX
        
        call ModifyDiffusivity_New(CurrProperty)
       
        FluxW     => Me%ExtVar%FluxW
    
        DWZ       => Me%ExtVar%DWZ
        DZZ       => Me%ExtVar%DZZ
        Theta     => Me%ExtVar%WaterContent
        ThetaOld  => Me%ExtVar%WaterContentOld
        

        !!!$OMP PARALLEL PRIVATE(I,J,K)
        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                Area      = Me%ExtVar%Area(i, j)
                
                Me%Volume(i,j,k)= Theta(i,j,k)*Me%ExtVar%Cellvolume(i,j,k)

               !!!FLUXES WITH Drainage Network
               CofInterfaceDN    = 0.0
               ConcInInterfaceDN = 0.0
                if (Me%ExtVar%CoupledDN) then
                    if(K == Me%ExtVar%GWLayer(i,j)) then
                       !Flux between river and runoff
                        if (Me%ExtVar%RiverPoints(I,J) == BasinPoint) then                        
                            
                            ! Positive flow -> looses mass
                            cofInterfaceDN = - aux * Me%ExtVar%FlowToChannels(i,j)
                            
                            ! mass going to channel -> conc from soil
                            if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
                                ConcInInterfaceDN =  CurrProperty%ConcentrationOld(i,j,k)
                            
                            !mass coming from channel -> conc from DN
                            elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                                ConcInInterfaceDN = CurrProperty%ConcentrationDN(i,j)
                            endif
                        
                        endif
                    endif
                endif

                !!TRANSPIRATION FLUXES IN POROUSMEDIA (Done in Routine VegetationInterfaceFluxes)
                    !Nitrate or Ortophosphate may be disconnected from transpiration flow (if SWAT orginal equations choosed)
                    !so TranspirationFlow * aux * Conc(i,j,k) may not be the uptake term
                    !in the future compute in ModuleVegetation fluxes as g/s (in routine get) instead of kg/ha
                    !and add here VegetationUptakeTerm = -[Nitrate or Ortophosphate flux] * aux  [g/s * s/m3 = g/m3]
                    !Mass sources from vegetation can still remain in Routine VegetationInterfaceFluxes
                
                !!EVAPORATION FLUXES IN POROUSMEDIA (Not needed to be accounted because it is not a mass flux)
                
                !!FLUXES IN X AND Y DIRECTION                
                if (Me%AdvDiff_SpatialMethod==AdvDif_CentralDif_) then ! diferenças centrais

                    !s/m3
                    aux      = (Me%ExtVar%DT/Me%Volume(i,j,k))

                    if (K == Me%WorkSize%KUB) then
                        ComputeCofD = .true.
                        ComputeCofC = .false.
                        cofC             = 0.0
                        DZZ(i,j,k)       = (DWZ(i,j,k)/2.) + (Me%ExtVar%InfiltrationColumn(i,j) / 2.)
!                        ConcTop          = CurrProperty%ConcentrationOnWaterColumn(i,j)
                        ConcTop          = CurrProperty%ConcentrationOnInfColumn(i,j)
                        
                        if (FluxW(i,j,k+1) .lt. 0.0) then !Negative - downwards, entering the soil 
                            AdvTermB_Top = 0.0
                            AdvTermD     = (aux * FluxW(i,j,k+1))
                            ! g/m3 = g / (m3/s * s)
!                            ConcTop  = CurrProperty%MassOnFluxW(i, j) / (abs(FluxW(i,j,k+1)) * Me%ExtVar%DT)
                            
                        else !Positive (upwards, exiting the soil) or zero
                            AdvTermB_Top = (aux * FluxW(i,j,k+1))
                            AdvTermD     = 0.0
                            !g/m3
!                            ConcTop  = PropertyX%ConcentrationOnWaterColumn(i, j)
                            
                        endif                      
                    else 
                        ComputeCofD = .false.
                        ComputeCofC = .true.
                        cofD             = 0.0

                        AdvTermB_Top     = (aux * FluxW(i,j,k+1) / 2.)
                        AdvTermC         = (aux * FluxW(i,j,k+1) / 2.)
                    endif
                    
                    AdvTermA = (aux * FluxW(i,j,k) / 2.)
                   !AdvTermB_Top    = already defined
                    AdvTermB_Bottom = (aux * FluxW(i,j,k  ) / 2.)
                   !AdvTermC        = already defined
                   !AdvTermD        = already defined
                       

                elseif (Me%AdvDiff_SpatialMethod==AdvDif_Upwind_) then ! upwind

                    !s/m3
                    aux      = (Me%ExtVar%DT/Me%Volume(i,j,k))

                    if (K == Me%WorkSize%KUB) then
                        ComputeCofD = .true.
                        ComputeCofC = .false.
                        cofC             = 0.0
                        DZZ(i,j,k)       = (DWZ(i,j,k)/2.) + (Me%ExtVar%InfiltrationColumn(i,j) / 2.)
!                        ConcTop  = CurrProperty%ConcentrationOnWaterColumn(i, j)
                        ConcTop          = CurrProperty%ConcentrationOnInfColumn(i,j)
                        
                    else 
                        ComputeCofD = .false.
                        ComputeCofC = .true.
                        cofD             = 0.0
                    endif
                    
                    if (FluxW(i,j,k) .lt. 0.0) then !Bottom face, Negative - downwards
                        AdvTermA        = 0.0
                        AdvTermB_Bottom = aux * FluxW(i,j,k)
                    else !Positive - upwards or zero.
                        AdvTermA        = aux * FluxW(i,j,k)
                        AdvTermB_Bottom = 0.0
                    endif

                    if (FluxW(i,j,k+1) .lt. 0.0) then !Top face, Negative - downwards
                        AdvTermC        = aux * FluxW(i,j,k+1)
                        AdvTermD        = aux * FluxW(i,j,k+1)
                        AdvTermB_Top    = 0.0
                    else !Positive - upwards or zero.
                        AdvTermC        = 0.0
                        AdvTermD        = 0.0
                        AdvTermB_Top    = aux * FluxW(i,j,k+1)
                    endif
                        
                endif
                
                DifTerm_Top    = CurrProperty%Diffusivity(i,j,k+1) * Area * aux / DZZ(i,j,k  )
                DifTerm_Bottom = CurrProperty%Diffusivity(i,j,k  ) * Area * aux / DZZ(i,j,k-1)
                
                cofA =  AdvTermA                                                                  &
                        + DifTerm_Bottom 
            
                cofB = ThetaOld(i,j,k)/Theta(i,j,k)                                               &
                       - AdvTermB_Top                                                             &
                       + AdvTermB_Bottom                                                          &
                       - DifTerm_Bottom                                                           &
                       - DifTerm_Top
                
                if (ComputeCofC) then    
                    cofC = - AdvTermC                                                             &
                           + DifTerm_Top
                endif
                
                if (ComputeCofD) then    
                    cofD = - AdvTermD                                                             &
                           + DifTerm_Top
                endif
                
                !Concentration computation based on the coefs (upwind or central differences)
                CurrProperty%Concentration(i,j,k)= cofA  * CurrProperty%ConcentrationOld(i,j,k-1) * Me%ExtVar%ComputeFacesW3D(i,j,k  ) &
                                                  + cofB * CurrProperty%ConcentrationOld(i,j,k  )                                      &
                                                  + cofC * CurrProperty%ConcentrationOld(i,j,k+1) * Me%ExtVar%ComputeFacesW3D(i,j,k+1) &
                                                  + cofD * ConcTop                                                                     &
                                                  + CofInterfaceDN * ConcInInterfaceDN                       
                
    
                Me%DifusionNumber(i,j,k) = cofB

                Me%ReynoldsMNumber(i,j,k)= cofA
                

                !Check if any of the coeffs get a negative value. If true, stop program
                if ((Me%AdvDiff_CheckCoefs) .AND. ((cofA < 0.0) .OR. (cofB < 0.0) .OR. (cofC < 0.0) .OR. (cofD < 0.0))) then
                    
                    call LogCoefs(i,j,k,cofA,cofB,cofC,cofD)
    
                endif
              
            endif
        enddo
        enddo
        enddo
        !!!$OMP END DO
        !!!$OMP END PARALLEL
                        
        call SetMatrixValue (CurrProperty%ConcentrationOld,      Me%Size,   CurrProperty%Concentration,          Me%ExtVar%WaterPoints3D)


    end subroutine ModifyAdvectionDiffusion_W
    
    !---------------------------------------------------------------------------

    subroutine ModifyAdvectionDiffusion_3D (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k!, CHUNK
        real(8)                                     :: Area_Top_W, Area_Bottom_W, Area_Top_U, Area_Bottom_U
        real(8)                                     :: Area_Top_V, Area_Bottom_V
        real(8)                                     :: AdvTermA_W, AdvTermA_U, AdvTermA_V
        real(8)                                     :: AdvTermB_Top_W, AdvTermB_Top_U, AdvTermB_Top_V
        real(8)                                     :: AdvTermB_Bottom_W, AdvTermB_Bottom_U, AdvTermB_Bottom_V
        real(8)                                     :: AdvTermC_W, AdvTermC_U, AdvTermC_V        
        real(8)                                     :: AdvTermD_W, AdvTermD_U, AdvTermD_V
        real(8)                                     :: DifTerm_Top_W, DifTerm_Bottom_W, DifTerm_Top_U
        real(8)                                     :: DifTerm_Bottom_U, DifTerm_Top_V, DifTerm_Bottom_V
        real(8)                                     :: aux 
        real(8)                                     :: cofA_W,cofB_W,cofC_W, cofD_W
        real(8)                                     :: cofA_U,cofB_U,cofC_U
        real(8)                                     :: cofA_V,cofB_V,cofC_V
        real(8)                                     :: cofB, cofInterfaceDN
        real(8)                                     :: ConcTop, ConcInInterfaceDN
        real(8), pointer, dimension(:,:,:)          :: FluxW, FluxU, FluxV
        real   , pointer, dimension(:,:,:)          :: DWZ, DZZ, Theta, ThetaOld, Porosity !, DZE, DZI
        real   , pointer, dimension(:,:  )          :: DZX, DZY, DXX, DYY !, DVX, DUY
        logical                                     :: ComputeCofC_W, ComputeCofD_W
        !Begin-----------------------------------------------------------------
   
        !!CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

        CurrProperty => PropertyX


        if (Me%AdvDiff_DiffMethod == AdvDif_Diff_Jury_) then ! new formulation based on Jury
            call ModifyDiffusivity_New(CurrProperty)
        elseif (Me%AdvDiff_DiffMethod == AdvDif_Diff_Old_) then !old formulation to couple module advection diffusion
            call ModifyDiffusivity_Old(CurrProperty)
        else
	        write(*,*)'Diffusion method to be used unrecognized,'
	        write(*,*)'Please check ADVDIFF_DIFF_METHOD keyword'
	        stop 'ModifyAdvectionDiffusion3D - ModulePorousMediaProperties - ERR100'        
        endif
        
        FluxW     => Me%ExtVar%FluxW
        FluxU     => Me%ExtVar%FluxU
        FluxV     => Me%ExtVar%FluxV        
        DWZ       => Me%ExtVar%DWZ
        DZZ       => Me%ExtVar%DZZ
!        DVX       => Me%ExtVar%DVX
!        DUY       => Me%ExtVar%DUY
        DZX       => Me%ExtVar%DZX
        DZY       => Me%ExtVar%DZY
        DXX       => Me%ExtVar%DXX
        DYY       => Me%ExtVar%DYY
!        DZE       => Me%ExtVar%DZE
!        DZI       => Me%ExtVar%DZI                
        Theta     => Me%ExtVar%WaterContent
        ThetaOld  => Me%ExtVar%WaterContentOld
        

        !!!$OMP PARALLEL PRIVATE(I,J,K)
        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                Area_Top_W         = Me%ExtVar%Area(i, j)
                Area_Bottom_W      = Me%ExtVar%Area(i, j)
                Area_Top_U         = DYY(i,j+1) * DWZ(i,j,k)
                Area_Bottom_U      = DYY(i,j  ) * DWZ(i,j,k)
                Area_Top_V         = DXX(i+1,j) * DWZ(i,j,k)
                Area_Bottom_V      = DXX(i,j  ) * DWZ(i,j,k)                
!                Area_Top_U         = DYY(i,j+1) * DZE(i,j+1,k)
!                Area_Bottom_U      = DYY(i,j  ) * DZE(i,j,k  )
!                Area_Top_V         = DXX(i+1,j) * DZI(i+1,j,k)
!                Area_Bottom_V      = DXX(i,j  ) * DZI(i,j,k  )
                
                Me%Volume(i,j,k)= Theta(i,j,k) * Me%ExtVar%Cellvolume(i,j,k)

               !!!FLUXES WITH Drainage Network
               CofInterfaceDN    = 0.0
               ConcInInterfaceDN = 0.0
                if (Me%ExtVar%CoupledDN) then 
                    if(K == Me%ExtVar%GWLayer(i,j)) then
                       !Flux between river and runoff
                        if (Me%ExtVar%RiverPoints(I,J) == BasinPoint) then                        
                            
                            ! Positive flow -> looses mass
                            cofInterfaceDN = - aux * Me%ExtVar%FlowToChannels(i,j)
                            
                            ! mass going to channel -> conc from soil
                            if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
                                ConcInInterfaceDN =  CurrProperty%ConcentrationOld(i,j,k)
                            
                            !mass coming from channel -> conc from DN
                            elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                                ConcInInterfaceDN = CurrProperty%ConcentrationDN(i,j)
                            endif
                        
                        endif
                    endif
                endif

                !!TRANSPIRATION FLUXES IN POROUSMEDIA (Done in Routine VegetationInterfaceFluxes)
                    !Nitrate or Ortophosphate may be disconnected from transpiration flow (if SWAT orginal equations choosed)
                    !so TranspirationFlow * aux * Conc(i,j,k) may not be the uptake term
                    !in the future compute in ModuleVegetation fluxes as g/s (in routine get) instead of kg/ha
                    !and add here VegetationUptakeTerm = -[Nitrate or Ortophosphate flux] * aux  [g/s * s/m3 = g/m3]
                    !Mass sources from vegetation can still remain in Routine VegetationInterfaceFluxes
                                    
                !!EVAPORATION FLUXES IN POROUSMEDIA (Not needed to be accounted because it is not a mass flux)
                
                !!FLUXES IN X AND Y DIRECTION 
                if (Me%AdvDiff_SpatialMethod==AdvDif_CentralDif_) then ! diferenças centrais

                    !s/m3
                    aux      = (Me%ExtVar%DT/Me%Volume(i,j,k))

                    if (K == Me%WorkSize%KUB) then
                        ComputeCofD_W = .true.
                        ComputeCofC_W = .false.
                        cofC_W           = 0.0
                        DZZ(i,j,k)       = DWZ(i,j,k)/2. + Me%ExtVar%InfiltrationColumn(i,j) / 2.
!                        ConcTop  = CurrProperty%ConcentrationOnWaterColumn(i, j)
                        ConcTop          = CurrProperty%ConcentrationOnInfColumn(i,j)
                        
                        if (FluxW(i,j,k+1) .lt. 0.0) then !Negative - downwards, entering the soil 
                            AdvTermB_Top_W = 0.0
                            AdvTermD_W     = (aux * FluxW(i,j,k+1))
                            ! g/m3 = g / (m3/s * s)
!                            ConcTop          = CurrProperty%MassOnFluxW(i, j) / (abs(FluxW(i,j,k+1)) * Me%ExtVar%DT)
                            
                        else !Positive (upwards, exiting the soil) or zero
                            AdvTermB_Top_W = (aux * FluxW(i,j,k+1))
                            AdvTermD_W     = 0.0
                            !g/m3
!                            ConcTop  = PropertyX%ConcentrationOnWaterColumn(i, j)

                        endif                      
                    else 
                        ComputeCofD_W = .false.
                        ComputeCofC_W = .true.
                        cofD_W             = 0.0
                        ConcTop            = 0.0
                        
                        AdvTermB_Top_W     = (aux * FluxW(i,j,k+1) / 2.)
                        AdvTermC_W         = (aux * FluxW(i,j,k+1) / 2.)
                    endif
                    
                    AdvTermA_W        = (aux * FluxW(i,j,k  ) / 2.)
                    AdvTermA_U        = (aux * FluxU(i,j,k  ) / 2.)
                    AdvTermA_V        = (aux * FluxV(i,j,k  ) / 2.)
                   !AdvTermB_Top_W    = already defined
                    AdvTermB_Bottom_W = (aux * FluxW(i,j,k  ) / 2.)
                    AdvTermB_Top_U    = (aux * FluxU(i,j+1,k) / 2.)
                    AdvTermB_Bottom_U = (aux * FluxU(i,j  ,k) / 2.)
                    AdvTermB_Top_V    = (aux * FluxV(i+1,j,k) / 2.)
                    AdvTermB_Bottom_V = (aux * FluxV(i  ,j,k) / 2.)
                   !AdvTermC_W        = already defined
                    AdvTermC_U        = (aux * FluxU(i,j+1,k) / 2.)
                    AdvTermC_V        = (aux * FluxV(i+1,j,k) / 2.)
                    
                        
                elseif (Me%AdvDiff_SpatialMethod==AdvDif_Upwind_) then ! upwind

                    !s/m3
                    aux      = (Me%ExtVar%DT/Me%Volume(i,j,k))

                    if (K == Me%WorkSize%KUB) then
                        ComputeCofD_W = .true.
                        ComputeCofC_W = .false.
                        cofC_W           = 0.0
                        DZZ(i,j,k)       = DWZ(i,j,k)/2. + Me%ExtVar%InfiltrationColumn(i,j) / 2.
!                       ConcTop  = CurrProperty%ConcentrationOnWaterColumn(i, j)
                        ConcTop          = CurrProperty%ConcentrationOnInfColumn(i,j)
                        

                    else 
                        ComputeCofD_W = .false.
                        ComputeCofC_W = .true.
                        cofD_W             = 0.0
                        ConcTop            = 0.0
                    endif
                    
                    !DirecW face k
                    if (FluxW(i,j,k) .lt. 0.0) then !Bottom face, Negative - exiting
                        AdvTermA_W        = 0.0
                        AdvTermB_Bottom_W = aux * FluxW(i,j,k)
                    else !Positive - entering or zero.
                        AdvTermA_W        = aux * FluxW(i,j,k)
                        AdvTermB_Bottom_W = 0.0
                    endif

                    !DirecU face j
                    if (FluxU(i,j,k) .lt. 0.0) then !Left face, Negative - exiting
                        AdvTermA_U        = 0.0
                        AdvTermB_Bottom_U = aux * FluxU(i,j,k)
                    else !Positive - entering or zero.
                        AdvTermA_U        = aux * FluxU(i,j,k)
                        AdvTermB_Bottom_U = 0.0
                    endif

                    !DirecV face i
                    if (FluxV(i,j,k) .lt. 0.0) then !Left face, Negative - exiting
                        AdvTermA_V        = 0.0
                        AdvTermB_Bottom_V = aux * FluxV(i,j,k)
                    else !Positive - entering or zero.
                        AdvTermA_V        = aux * FluxV(i,j,k)
                        AdvTermB_Bottom_V = 0.0
                    endif
                    
                    !DirecW face k+1
                    if (FluxW(i,j,k+1) .lt. 0.0) then !Top face, Negative - entering
                        AdvTermC_W        = aux * FluxW(i,j,k+1)
                        AdvTermD_W        = aux * FluxW(i,j,k+1)
                        AdvTermB_Top_W    = 0.0
                    else !Positive - exiting or zero.
                        AdvTermC_W        = 0.0
                        AdvTermD_W        = 0.0
                        AdvTermB_Top_W    = aux * FluxW(i,j,k+1)
                    endif

                    !DirecU face j+1
                    if (FluxU(i,j+1,k) .lt. 0.0) then !Right face, Negative - entering
                        AdvTermC_U        = aux * FluxU(i,j+1,k)
                        AdvTermD_U        = aux * FluxU(i,j+1,k)
                        AdvTermB_Top_U    = 0.0
                    else !Positive - exiting or zero.
                        AdvTermC_U        = 0.0
                        AdvTermD_U        = 0.0
                        AdvTermB_Top_U    = aux * FluxU(i,j+1,k)
                    endif                    

                    !DirecV face i+1
                    if (FluxV(i+1,j,k) .lt. 0.0) then !Right face, Negative - entering
                        AdvTermC_V        = aux * FluxV(i+1,j,k)
                        AdvTermD_V        = aux * FluxV(i+1,j,k)
                        AdvTermB_Top_V    = 0.0
                    else !Positive - exiting or zero.
                        AdvTermC_V        = 0.0
                        AdvTermD_V        = 0.0
                        AdvTermB_Top_V    = aux * FluxV(i+1,j,k)
                    endif                    
                        

                endif
                
                DifTerm_Top_W    = CurrProperty%Diffusivity(i,j,k+1) * Area_Top_W    * aux / DZZ(i,j,k  )
                DifTerm_Bottom_W = CurrProperty%Diffusivity(i,j,k  ) * Area_Bottom_W * aux / DZZ(i,j,k-1)
                DifTerm_Top_U    = CurrProperty%ViscosityU(i,j+1,k)  * Area_Top_U    * aux / DZX(i,j  )
                DifTerm_Bottom_U = CurrProperty%ViscosityU(i,j,k  )  * Area_Bottom_U * aux / DZX(i,j-1)
                DifTerm_Top_V    = CurrProperty%ViscosityV(i+1,j,k)  * Area_Top_V    * aux / DZY(i  ,j)
                DifTerm_Bottom_V = CurrProperty%ViscosityV(i,j,k  )  * Area_Bottom_V * aux / DZY(i-1,j) 
                
                cofA_W =  AdvTermA_W                                                          &
                          + DifTerm_Bottom_W 

                cofA_U = AdvTermA_U                                                           &
                          + DifTerm_Bottom_U 
                
                cofA_V = AdvTermA_V                                                           &
                          + DifTerm_Bottom_V 
            
                cofB_W = - AdvTermB_Top_W                                                     &
                         + AdvTermB_Bottom_W                                                  &
                         - DifTerm_Bottom_W                                                   &
                         - DifTerm_Top_W

                cofB_U = - AdvTermB_Top_U                                                     &
                         + AdvTermB_Bottom_U                                                  &
                         - DifTerm_Bottom_U                                                   &
                         - DifTerm_Top_U
                
                cofB_V = - AdvTermB_Top_V                                                     &
                         + AdvTermB_Bottom_V                                                  &
                         - DifTerm_Bottom_V                                                   &
                         - DifTerm_Top_V         
                
                if (ComputeCofC_W) then    
                    cofC_W = - AdvTermC_W                                                     &
                             + DifTerm_Top_W
                endif
                 
                cofC_U = - AdvTermC_U                                                         &
                         + DifTerm_Top_U                   

                cofC_V = - AdvTermC_V                                                         &
                         + DifTerm_Top_V    
                
                if (ComputeCofD_W) then    
                    cofD_W = - AdvTermD_W                                                     &
                             + DifTerm_Top_W 
                endif


                CofB = ThetaOld(i,j,k)/Theta(i,j,k) + cofB_W + cofB_U + cofB_V                                  
                
                CurrProperty%Concentration(i,j,k)=  cofA_W * CurrProperty%ConcentrationOld(i,j,k-1)* Me%ExtVar%ComputeFacesW3D(i,j,k  ) &
                                                  + cofC_W * CurrProperty%ConcentrationOld(i,j,k+1)* Me%ExtVar%ComputeFacesW3D(i,j,k+1) &
                                                  + cofD_W * ConcTop                                                                    &
                                                  + cofA_U * CurrProperty%ConcentrationOld(i,j-1,k)* Me%ExtVar%ComputeFacesU3D(i,j  ,k) &
                                                  + cofC_U * CurrProperty%ConcentrationOld(i,j+1,k)* Me%ExtVar%ComputeFacesU3D(i,j+1,k) &
                                                  + cofA_V * CurrProperty%ConcentrationOld(i-1,j,k)* Me%ExtVar%ComputeFacesV3D(i  ,j,k) &
                                                  + cofC_V * CurrProperty%ConcentrationOld(i+1,j,k)* Me%ExtVar%ComputeFacesV3D(i+1,j,k) &
                                                  + cofB   * CurrProperty%ConcentrationOld(i,j,k  )                                     &
                                                  + cofInterfaceDN * ConcInInterfaceDN              
                                                  
          
                
                !Check if any of the coeffs get a negative value. If true, stop program
                if ((Me%AdvDiff_CheckCoefs) .AND. ((cofA_W < 0.0) .OR. (cofC_W < 0.0) .OR. (cofD_W < 0.0) .OR.     &
                                                   (cofA_U < 0.0) .OR. (cofC_U < 0.0) .OR.                         &
                                                   (cofA_V < 0.0) .OR. (cofC_V < 0.0) .OR. (cofB < 0.0) )) then
                    
                    call LogCoefs (i,j,k,cofA_W,cofB,cofC_W,cofD_W,cofA_U,cofC_U,cofA_V,cofC_V)                                                  
               
                endif
              
            endif
        enddo
        enddo
        enddo
        !!!$OMP END DO
        !!!$OMP END PARALLEL
                        
        call SetMatrixValue (CurrProperty%ConcentrationOld,      Me%Size,   CurrProperty%Concentration,          Me%ExtVar%WaterPoints3D)


    end subroutine ModifyAdvectionDiffusion_3D
    
    !---------------------------------------------------------------------------
    
    subroutine ModifyDiffusivity_New(CurrProperty)
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        !Local-----------------------------------------------------------------
        integer                                     :: I,J,K, CHUNK
        real   , pointer, dimension(:,:  )          :: WaterCol
        real   , pointer, dimension(:,:,:)          :: ThetaOld, Porosity
        real(8), pointer, dimension(:,:,:)          :: UnsatW, FluxW, UnsatU, UnsatV
        real                                        :: WaterContent_Face, Porosity_Face
        real                                        :: DiffCoef        
        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)
        
        !$OMP PARALLEL PRIVATE(I,J,K,DiffCoef,WaterContent_Face,Porosity_Face)

        ThetaOld  => Me%ExtVar%WaterContentOld
        Porosity  => Me%ExtVar%ThetaS
        WaterCol  => Me%ExtVar%InfiltrationColumn
        UnsatW    => Me%ExtVar%UnsatW
        DiffCoef  = CurrProperty%Evolution%AdvDiff%Molecular_Diff_Coef
        
        if (Me%AdvDiff_ComputeTransport3D) then 
            UnsatU    => Me%ExtVar%UnsatU
            UnsatV    => Me%ExtVar%UnsatV
        endif
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                !condition so that does not produce NAN, in WaterContetFace in the boundary
                !in boundary, diffusivity is zero (computefaces) but in release version model when evaluates disp/watercontent_face
                !produces NAN. In debug version this does not appen and the result is the total evaluation (zero).
                if (Me%ExtVar%WaterPoints3D(I,J,K-1) /= WaterPoint) then
                    WaterContent_Face = ThetaOld(i,j,k)
                else
                    WaterContent_Face = min(ThetaOld(i,j,k),ThetaOld(i,j,k-1))
                endif
                Porosity_Face     = min(Porosity(i,j,k),Porosity(i,j,k-1))
                
                !Only compute diffusivity in compute faces W (faces tath are not boundaries), else is zero
                CurrProperty%Diffusivity(i,j,k)  = Me%ExtVar%ComputeFacesW3D(i,j,k  ) *                                             &
                                                    (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))                        &
                                                     +(abs(UnsatW(i,j,k  )) * Me%Disper_Longi%Field(i,j,k) / WaterContent_Face)
                
                if (Me%AdvDiff_ComputeTransport3D) then 
                    if (Me%ExtVar%WaterPoints3D(I,J-1,K) /= WaterPoint) then
                        WaterContent_Face = ThetaOld(i,j,k)
                    else                                                                                       
                        WaterContent_Face = min(ThetaOld(i,j,k),ThetaOld(i,j-1,k))
                    endif
                    Porosity_Face     = min(Porosity(i,j,k),Porosity(i,j-1,k))
                    
                    !Only compute diffusivity in compute faces U (faces tath are not boundaries), else is zero
                    CurrProperty%ViscosityU(i,j,k)  = Me%ExtVar%ComputeFacesU3D(i,j,k  ) *                                         &
                                                        (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))                    &
                                                         +(abs(UnsatU(i,j,k  )) * Me%Disper_Trans%Field(i,j,k) / WaterContent_Face)
                    if (Me%ExtVar%WaterPoints3D(I-1,J,K) /= WaterPoint) then
                        WaterContent_Face = ThetaOld(i,j,k)
                    else                                                           
                        WaterContent_Face = min(ThetaOld(i,j,k),ThetaOld(i-1,j,k))
                    endif
                    Porosity_Face     = min(Porosity(i,j,k),Porosity(i-1,j,k))
                    
                    !Only compute diffusivity in compute faces V (faces tath are not boundaries), else is zero
                    CurrProperty%ViscosityV(i,j,k)  = Me%ExtVar%ComputeFacesV3D(i,j,k  ) *                                         &
                                                        (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))                    &
                                                         +(abs(UnsatV(i,j,k  )) * Me%Disper_Trans%Field(i,j,k) / WaterContent_Face)
                endif
                                                         

            endif                                                                        
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            K = Me%WorkSize%KUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                !UpperFace computation
                WaterContent_Face = ThetaOld(i,j,k)
                Porosity_Face     = Porosity(i,j,k)
                !Only compute diffusivity in compute faces W and if there is transport in the interface
                if (WaterCol(i,j) .gt. 0.) then
                    CurrProperty%Diffusivity(i,j,k+1)  =  (DiffCoef * Tortuosity(WaterContent_Face, Porosity_Face))             &
                                                          +(abs(UnsatW(i,j,k  )) * Me%Disper_Longi%Field(i,j,k) / WaterContent_Face)
                else
                    CurrProperty%Diffusivity(i,j,k+1)  = 0.0
                endif
            endif                                                                        
        enddo
        enddo

    
    end subroutine ModifyDiffusivity_New

    !---------------------------------------------------------------------------

    subroutine SoilQualityProcesses

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property),          pointer     :: PropertyX
        type (T_Property),          pointer     :: SoilDryDensity, Salinity, pH
        type (T_Property),          pointer     :: IonicStrength, PhosphorusAdsortionIndex

!        type (T_SoilRate),      pointer     :: SoilRateX
!        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB 
!        integer                                 :: i, j, k
        
        !Begin-----------------------------------------------------------------
        
!        WIUB = Me%WorkSize%IUB
!        WJUB = Me%WorkSize%JUB
!        WILB = Me%WorkSize%ILB
!        WJLB = Me%WorkSize%JLB
!        WKUB = Me%WorkSize%KUB
!        WKLB = Me%WorkSize%KLB
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "SoilQualityProcesses")

        call ComputeDissolvedToParticulate3D
        
        !Properties not modified by sediment quality (not state variables) but needed in argument
        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property soil dry density not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR00'
        endif

        call SearchProperty(Salinity, Salinity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property salinity not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR10'
        endif

        call SearchProperty(pH, pH_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property pH not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR20'
        endif

        call SearchProperty(IonicStrength, IonicStrength_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property ionic strength not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR30'

        endif

        call SearchProperty(PhosphorusAdsortionIndex, PhosphorusAdsortionIndex_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property phosphorus sdsortion index not found in porous media properties'
            stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR40'
        endif

        call ComputeWindVelocity
 
        
        if (Me%ExtVar%Now .GE. Me%Coupled%SoilQuality_NextCompute) then
            
            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))
                

                call Modify_Interface(InterfaceID               = Me%ObjInterface,                         &
                                      PropertyID                = PropertyX%ID%IDNumber,                   &
                                      Concentration             = PropertyX%Concentration,                 &
                                      WaterPoints3D             = Me%ExtVar%WaterPoints3D,                 &
                                      OpenPoints3D              = Me%ExtVar%OpenPoints3D,                  &
                                      WaterPercentage           = Me%ExtVar%WaterContent,                  &
                                      DissolvedToParticulate3D  = Me%DissolvedToParticulate3D,             &
                                      SoilDryDensity            = SoilDryDensity%Concentration,            &
                                      Salinity                  = Salinity%Concentration,                  &
                                      pH                        = pH%Concentration,                        &
                                      IonicStrength             = IonicStrength%Concentration,             &
                                      PhosphorusAdsortionIndex  = PhosphorusAdsortionIndex%Concentration,  &
                                      WindVelocity              = Me%ExtVar%WindVelocity3D,                &
                                      STAT                      = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                                               &
                    stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR050'
                

                PropertyX => PropertyX%Next
                

            end do
            
            Me%Coupled%SoilQuality_NextCompute = Me%Coupled%SoilQuality_NextCompute +       &
                                                     Me%Coupled%SoilQuality_DT

        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))


            if (PropertyX%Evolution%SoilQuality) then

                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%Concentration,          &
                                          WaterPoints3D = Me%ExtVar%WaterPoints3D,          &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'SoilQuality_Processes - ModulePorousMediaProperties - ERR060'

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SoilQualityProcesses")

    end subroutine SoilQualityProcesses
    
    !-----------------------------------------------------------------------------    

    subroutine ComputeWindVelocity

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                      :: i, j, k
        
        !Begin-----------------------------------------------------------------

        call SetMatrixValue(Me%ExtVar%WindVelocity3D, Me%Size, 0.0, Me%ExtVar%WaterPoints3D)

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            k = Me%WorkSize%KUB
        
            if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
                
                ! km/day                        =  m/s  * 1E-3km/m * 86400s/day 
                Me%ExtVar%WindVelocity3D(i,j,k) = Me%ExtVar%WindVelocity2D(i,j) * 1E-3 * 86400
            endif

        enddo
        enddo

    end subroutine ComputeWindVelocity

    !-----------------------------------------------------------------------------    
    
    
    subroutine SetLimitsConcentration

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k!, CHUNK
        
        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "SetLimitsConcentration")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB 


        Property => Me%FirstProperty  

do1 :   do while (associated(Property))
cd1 :       if (Property%Evolution%MinConcentration) then
                
!                CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)
                
!                !$OMP PARALLEL SHARED(CHUNK, Property) PRIVATE(I,J,K)
!                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)

                do k=Me%WorkSize%KLB, Me%WorkSize%KUB
                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%WaterPoints3D(i, j, k) == 1) then
    
                        if (Property%Concentration(i, j, k) < Property%MinValue) then
                            
                            ! mass created
                            Property%Mass_created(i, j, k) = Property%Mass_Created(i, j, k)   +  &
                                                   (Property%MinValue                -  &
                                                    Property%Concentration(i, j, k)) *  (Me%ExtVar%WaterContent(i,j,k) * &
                                                    Me%ExtVar%CellVolume (i, j, k))

                            Property%Concentration(i, j, k) = Property%MinValue
                            
                        endif

                    endif

                enddo
                enddo
                enddo
                
!                !$OMP END DO NOWAIT
!                !$OMP END PARALLEL
                
            endif cd1
                
        Property => Property%Next
        end do do1

        nullify(Property)

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SetLimitsConcentration")


    end subroutine SetLimitsConcentration

    !--------------------------------------------------------------------------
    

    subroutine ComputeDissolvedToParticulate3D

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL        
         
        !Local----------------------------------------------------------------- 
        integer                                 :: i, j, k
        real                                    :: DT, InstantValue, ResidualValue
        type(T_Property), pointer               :: SoilDryDensity!, DrySedimentVolume
        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "ComputeDissolvedToParticulate3D")


        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR01'

        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModulePorousMediaProperties - ERR02'


        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%ExtVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                ! [m3water/kgsed] = [m3water/m3cell]*[m3cell] / ([kgsed/m3cell] * [m3cell]) 
                InstantValue                       = Me%ExtVar%WaterContent(i,j,k) *  Me%ExtVar%CellVolume(i,j,k) / &
                                                    (SoilDryDensity%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k))

                ResidualValue                      = Me%DissolvedToParticulate3D(i,j,k)

                Me%DissolvedToParticulate3D(i,j,k) = (ResidualValue * Me%ResidualTime +     &
                                                      InstantValue * DT) / (Me%ResidualTime + DT)
                                                       
            end if
        end do
        end do
        end do

        Me%ResidualTime = Me%ResidualTime + DT
        
        
        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "ComputeDissolvedToParticulate3D")


    end subroutine ComputeDissolvedToParticulate3D
    !--------------------------------------------------------------------------

#ifdef _PHREEQC_
    !--------------------------------------------------------------------------
    subroutine SoilChemistryProcesses

        !External--------------------------------------------------------------
        integer :: STAT_CALL        
        
        !Local----------------------------------------------------------------- 
        type (T_Property), pointer                   :: PropertyX, pH
        integer                                      :: I, J, K
        real             , pointer, dimension(:,:,:) :: Theta
        real(8)          , pointer, dimension(:,:,:) :: Volume
        
        !Begin-----------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModulePorousMediaProperties", "SoilChemistryProcesses")

        Theta  => Me%ExtVar%WaterContent
        Volume => Me%ExtVar%Cellvolume
        
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                        
                Me%ExtVar%CellWaterMass(I, J, K) = Theta(I, J, K) * Volume(I, J, K) * WaterReferenceDensity    
                        
            end if
                    
        end do
        end do
        end do
        

        call SearchProperty(pH, pH_ , .false., STAT = STAT_CALL)    
            
        if (STAT_CALL /= SUCCESS_) then
            write(*,*) 'property pH not found in porous media properties'
            stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR01'
        endif

 
        !Question: Why here is used SoilChemistry_NextCompute and after is used NextCompute?
        !          This is related to the possibility of different DT's beetwen 0D model and MOHID general DT?
        !          
        if (Me%ExtVar%Now .GE. Me%Coupled%SoilChemistry_NextCompute) then
            
            PropertyX => Me%FirstProperty

            do while(associated(PropertyX))
                
                call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry , &
                                      PropertyID    = PropertyX%ID%IDNumber        , &
                                      WaterMass     = Me%ExtVar%CellWaterMass      , &                                                                       
                                      Concentration = PropertyX%Concentration      , &
                                      WaterPoints3D = Me%ExtVar%WaterPoints3D      , &
                                      pH            = pH%Concentration             , &
                                      OpenPoints3D  = Me%ExtVar%OpenPoints3D       , &
                                      STAT          = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR02'
                
                PropertyX => PropertyX%Next
                
            end do
            
            Me%Coupled%SoilChemistry_NextCompute = Me%Coupled%SoilChemistry_NextCompute + &
                                                   Me%Coupled%SoilChemistry_DT

        end if

        PropertyX => Me%FirstProperty

        do while(associated(PropertyX))

            if (PropertyX%Evolution%SoilChemistry) then

                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then

                    call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry,                  &
                                          PropertyID    = PropertyX%ID%IDNumber,            &
                                          Concentration = PropertyX%Concentration,          &
                                          WaterPoints3D = Me%ExtVar%WaterPoints3D,          &
                                          DTProp        = PropertyX%Evolution%DTInterval,   &
                                          STAT          = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                            &
                        stop 'SoilChemistryProcesses - ModulePorousMediaProperties - ERR03'

                end if

            end if

            PropertyX => PropertyX%Next
            
        end do

        if (MonitorPerformance) call StopWatch ("ModulePorousMediaProperties", "SoilChemistryProcesses")
        
        !End-------------------------------------------------------------------
        
    end subroutine SoilChemistryProcesses   
    !-----------------------------------------------------------------------------    
#endif

    !--------------------------------------------------------------------------
    ! This subroutine is responsable for defining       
    ! the next time to actualize the value of each      
    ! property                                          
    subroutine Actualize_Time_Evolution

        !Local--------------------------------------------------------------
        type (T_Property), pointer :: Property
        type (T_Time    )          :: Actual

        !----------------------------------------------------------------------

        Property => Me%FirstProperty  

        Actual = Me%ExtVar%Now

do1 :   do while (associated(Property))

cd1:        if (Property%Evolution%Variable) then
cd2 :       if (Actual.GE.Property%Evolution%NextCompute) then
                    Property%Evolution%LastCompute = Property%Evolution%NextCompute
                    Property%Evolution%NextCompute = Property%Evolution%NextCompute &
                                                   + Property%Evolution%DTInterval
            end if cd2
            end if cd1


            Property => Property%Next
        end do do1   

        nullify(Property)


    end subroutine Actualize_Time_Evolution

    
    !--------------------------------------------------------------------------


    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX

        !----------------------------------------------------------------------

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data3D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR01'

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data2D = PropertyX%ConcentrationOnInfColumn,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModulePorousMediaProperties - ERR02'

            

            endif
            PropertyX=>PropertyX%Next
        enddo

    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------


    subroutine OutPut_HDF

        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        integer                                     :: STAT_CALL
        type(T_Time)                                :: Actual, LastTime, EndTime
        integer                                     :: OutPutNumber
        real, dimension(6), target                  :: AuxTime
        real, dimension(:), pointer                 :: TimePtr
             
        !Begin----------------------------------------------------------------

        Actual   = Me%ExtVar%Now
        EndTime  = Me%ExtVar%EndTime
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
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR00'

                    call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                         Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR01'
           
                    Me%LastOutPutHDF5 = Actual
       
                endif First

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                      Me%WorkSize%ILB,                           &
                                      Me%WorkSize%IUB,                           &
                                      Me%WorkSize%JLB,                           &
                                      Me%WorkSize%JUB,                           &
                                      Me%WorkSize%KLB,                           &
                                      Me%WorkSize%KUB,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR02'

                !Writes the Open Points
                call HDF5WriteData   (Me%ObjHDF5, "//Grid/OpenPoints",              &
                                      "OpenPoints", "-",                            &
                                      Array3D = Me%ExtVar%OpenPoints3D,             &
                                      OutputNumber = OutPutNumber,                  &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR03'


                PropertyX => Me%FirstProperty
                do while (associated(PropertyX))

                    if (PropertyX%OutputHDF) then
 
                        call HDF5WriteData   (Me%ObjHDF5,                                    &
                                              "/Results/"//trim(PropertyX%ID%Name),          &
                                              trim(PropertyX%ID%Name),                       &
                                              trim(PropertyX%ID%Units),                      &
                                              Array3D = PropertyX%Concentration,             &
                                              OutputNumber = OutPutNumber,                   &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR04'

                    endif

                    PropertyX => PropertyX%Next

                enddo

                Me%OutPut%NextOutput = OutPutNumber + 1

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModulePorousMediaproperties - ERR06'
            
            endif  TOut
        endif  TNum

    end subroutine OutPut_HDF

    !----------------------------------------------------------------------------

    real function Tortuosity(WC, Porosity)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: WC
        real, intent(IN)                            :: Porosity
        !Begin----------------------------------------------------------------

        !tortuosity = (WC**(10/3))/(porosity)
        Tortuosity = (WC**(7/3))/(Porosity**2)
        

    end function Tortuosity   


    !----------------------------------------------------------------------------

    subroutine ProfileOutput

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                             :: STAT_CALL
        type (T_Property), pointer                          :: PropertyX
        integer                                             :: nProperties
        integer                                             :: n                                    

        nProperties = Me%PropertiesNumber 


     
        n=1
        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))

        call WriteProfile(Me%ObjProfile,                                        &
                          Data3D = PropertyX%Concentration,                              &
                          SZZ    = Me%ExtVar%SZZ,                               &
                          STAT   = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ProfileOutput - ModulePorousMedia - ERR01'

        PropertyX=>PropertyX%Next

        endDo
    end subroutine ProfileOutput


    !------------------------------------------------------------------------------ 

    subroutine LogCoefs(i,j,k,Coef1, Coef2, Coef3, Coef4, Coef5, Coef6, Coef7, Coef8)

        !Arguments-------------------------------------------------------------
        integer                                     :: i,j,k
        real(8)                                     :: Coef1, Coef2, Coef3
        real(8), optional                           :: Coef4, Coef5, Coef6, Coef7, Coef8
 
        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second

        call ExtractDate(Me%ExtVar%Now, Year, Month, Day, Hour, Minute, Second)
        
        if (.not. present(Coef5)) then
            write (Me%Files%AsciiUnit, fmt=1000) Year, Month, Day, Hour, Minute, Second,    &
                                                 i,j,k, Coef1, Coef2, Coef3, Coef4
        else
            write (Me%Files%AsciiUnit, fmt=1000) Year, Month, Day, Hour, Minute, Second,    &
                                                 i,j,k, Coef1, Coef2, Coef3, Coef4, Coef5,  & 
                                                 Coef6, Coef7, Coef8
        endif

        1000 format(f5.0, f5.0, f5.0, f5.0, f5.0, f12.5, i3, i3, i3, 8f13.8)

    end subroutine LogCoefs

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillPorousMediaProperties(ObjPorousMediaPropertiesID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjPorousMediaPropertiesID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers,STAT_CALL  
        type(T_property), pointer           :: PropertyX
        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPorousMediaProperties_,  Me%InstanceID)


            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 
                if(PropertyX%ID%SolutionFromFile)then

                    call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR00'
                end if
                
                PropertyX => PropertyX%Next
            end do 
            
            if (associated(Me%Disper_Longi%Field))then
                deallocate(Me%Disper_Longi%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR01'
                nullify   (Me%Disper_Longi%Field)
            end if

            if (associated(Me%Disper_Trans%Field))then
                deallocate(Me%Disper_Trans%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillPorousMediaProperties - ModulePorousMediaProperties - ERR02'
                nullify   (Me%Disper_Trans%Field)
            end if
            
            if (nUsers == 0) then

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - PorousmediaProperties - ERR05'
                endif

               
                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - PorousmediaProperties  - ERR08'
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR08'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR11'
                
                nUsers = DeassociateInstance (mPOROUSMEDIA_,  Me%ObjPorousMedia)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR12'

                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR13'

                nUsers = DeassociateInstance (mMAP_,  Me%ObjMap)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR14'

                if ((Me%Coupled%AdvectionDiffusion) .and. (Me%AdvDiff_Module == AdvDif_ModuleAD_)) then
                    call KillAdvectionDiffusion(Me%ObjAdvectionDiffusion, STAT = STAT_CALL)
                    
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillPorousMedia - PorousmediaProperties - ERR15'
                end if

                
                call DeallocateVariables

                !Deallocates Instance
                call DeallocateInstance ()

                ObjPorousMediaPropertiesID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillPorousMediaProperties
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_PorousMediaProperties), pointer          :: AuxObjPorousMediaProperties
        type (T_PorousMediaProperties), pointer          :: PreviousObjPorousMediaProp

        !Updates pointers
        if (Me%InstanceID == FirstObjPorousMediaProperties%InstanceID) then
            FirstObjPorousMediaProperties => FirstObjPorousMediaProperties%Next
        else
            PreviousObjPorousMediaProp => FirstObjPorousMediaProperties
            AuxObjPorousMediaProperties      => FirstObjPorousMediaProperties%Next
            do while (AuxObjPorousMediaProperties%InstanceID /= Me%InstanceID)
                PreviousObjPorousMediaProp => AuxObjPorousMediaProperties
                AuxObjPorousMediaProperties      => AuxObjPorousMediaProperties%Next
            enddo

            !Now update linked list
            PreviousObjPorousMediaProp%Next => AuxObjPorousMediaProperties%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance
    
    !--------------------------------------------------------------------------

    subroutine DeallocateVariables

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !Water Content---------------------------------------------------------
        deallocate (Me%Volume                  )
        deallocate (Me%DifusionNumber          )
        deallocate (Me%ReynoldsMNumber         )
        deallocate (Me%ExtVar%WindVelocity3D   )
        
#ifdef _PHREEQC_
        deallocate (Me%ExtVar%CellWaterMass)
#endif        

        deallocate (Me%WaterVolume)
        deallocate (Me%WaterVolumeOld)
        deallocate (Me%WaterVolumeCorr)
 !       deallocate (Me%WaterVolumeOldCorr)
        deallocate (Me%FluxWCorr)
      
    
    end subroutine DeallocateVariables 

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjPorousMediaProperties_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaProperties_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjPorousMediaProperties_ID > 0) then
            call LocateObjPorousMediaProperties (ObjPorousMediaProperties_ID)
            ready_ = VerifyReadLock (mPorousMediaProperties_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjPorousMediaProperties (ObjPorousMediaPropertiesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaPropertiesID

        !Local-----------------------------------------------------------------

        Me => FirstObjPorousMediaProperties
        do while (associated (Me))
            if (Me%InstanceID == ObjPorousMediaPropertiesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModulePorousMediaProperties - LocateObjPorousMediaProperties - ERR01'

    end subroutine LocateObjPorousMediaProperties

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
                if (PrintWarning) write (*,*)'Property Not Found in Module PorousMediaProperties ', &
                                              trim(GetPropertyName(PropertyXIDNumber))
            endif
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchProperty

    !--------------------------------------------------------------------------


    subroutine ReadLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        call GetOldWaterContent (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR010'

        call GetWaterContent    (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR020'

        if ((Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) .or. Me%AdvDiff_Module == AdvDif_ModuleAD_) then 
            call GetFluxU           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR030'

            call GetFluxV           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR040'
            
            call GetUnsatV          (Me%ObjPorousMedia, Me%ExtVar%UnsatV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR061'

            call GetUnsatU          (Me%ObjPorousMedia, Me%ExtVar%UnsatU, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR062'   
        endif  
        
        call GetFluxW           (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR050'

        call GetUnsatW          (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'
        
        call GetThetaS          (Me%ObjPorousMedia, Me%ExtVar%ThetaS, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR65'

        call GetPotentialInfiltration (Me%ObjPorousMedia, Me%ExtVar%InfiltrationColumn, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - Module ModulePorousMediaProperties. ERR10.'

!        call GetWaterColumn     (Me%ObjPorousMedia, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR64'


        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR070'

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR080'

        !LandPoints3D
        call GetLandPoints3D    (Me%ObjMap, Me%ExtVar%LandPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSoilProperties - ERR11'

        !OpenPoints3D
        call GetOpenPoints3D    (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR085'

        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ    = Me%ExtVar%CellVolume,                      &
                                STAT       = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR090'
                       
        call GetGeometryDistances (Me%ObjGeometry,                                      &
                                  SZZ         = Me%ExtVar%SZZ,                          &
                                  DWZ         = Me%ExtVar%DWZ,                          &
                                  DZZ         = Me%ExtVar%DZZ,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR100'

        if (Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) then

!            call GetGeometryDistances (Me%ObjGeometry,                                      &
!                                      DZE         = Me%ExtVar%DZE,                          &
!                                      DZI         = Me%ExtVar%DZI,                          &
!                                      STAT        = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR100.1'
        
            call GetHorizontalGrid(Me%ObjHorizontalGrid,                                    &
!                                      DUY         = Me%ExtVar%DUY,                          &
                                      DZY         = Me%ExtVar%DZY,                          &
!                                      DVX         = Me%ExtVar%DVX,                          &
                                      DZX         = Me%ExtVar%DZX,                          &
                                      DXX         = Me%ExtVar%DXX,                          &
                                      DYY         = Me%ExtVar%DYY,                          &
                                      STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR101'
        endif

        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR110'


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesW3D = Me%ExtVar%ComputeFacesW3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR120'

        if ((Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) .or. Me%AdvDiff_Module == AdvDif_ModuleAD_) then
            call GetComputeFaces3D(Me%ObjMap,                                               &
                                   ComputeFacesU3D = Me%ExtVar%ComputeFacesU3D,             &
                                   ComputeFacesV3D = Me%ExtVar%ComputeFacesV3D,             &
                                   STAT            = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR130'
        endif
        
        if (Me%ExtVar%CoupledDN) then
            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR0140'

            call GetGWFlowToChannels   (Me%ObjPorousMedia, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR150'

            call GetGWLayer   (Me%ObjPorousMedia, Me%ExtVar%GWlayer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR160'
        endif
        
    end subroutine ReadLockExternalVar

    !-----------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR010'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR020'

        if ((Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) .or. Me%AdvDiff_Module == AdvDif_ModuleAD_) then 
            call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR030'

            call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR040'
            
            call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%UnsatV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR061'
            
            call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%UnsatU, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR064'
        endif
        
        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR050'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR060'

        call UnGetPorousMedia           (Me%ObjPorousMedia,Me%ExtVar%ThetaS, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR061'
       
        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR070'
        
!        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterColumn, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR063'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%InfiltrationColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR063'


        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR080'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%LandPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR081'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR085'
        
        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%CellVolume,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR090'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%SZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR100'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DWZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR110'
        
        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR111'


        if (Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) then

!            call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DZI,  STAT = STAT_CALL )        
!            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR111.1'
!            
!            call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DZE,  STAT = STAT_CALL )        
!            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR111.2'
            
!            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVX, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR112'
 
             call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR113'

!            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUY, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR114'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR115'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DXX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR116'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DYY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR117'            

        endif

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR120'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR130'
        
        if ((Me%AdvDiff_Module == AdvDif_ModulePMP_ .and. Me%AdvDiff_ComputeTransport3D) .or. Me%AdvDiff_Module == AdvDif_ModuleAD_) then
            call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesU3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR140'

            call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesV3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR150'
        endif       

        if (Me%ExtVar%CoupledDN) then
            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0160'               

            call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0170'               

            call UnGetPorousMedia (Me%ObjPorousMedia, Me%ExtVar%GWlayer, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModulePorousMediaProperties - ERR0180'      
        
        endif
        
    endsubroutine ReadUnlockExternalVar


end module ModulePorousMediaProperties

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








