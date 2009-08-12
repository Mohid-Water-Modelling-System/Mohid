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
    use ModuleBasinGeometry,    only : GetBasinPoints,  UnGetBasin 
                                       
    use ModuleFillMatrix,       only : ConstructFillMatrix, KillFillMatrix
    use ModuleGeometry
    use ModuleMap
    use ModulePorousMedia,      only : GetOldWaterContent, GetWaterContent, GetFluxU,    &
                                       GetFluxV, GetFluxWOld, GetUnsatWOld, UnGetPorousMedia
    use ModuleInterface,        only : ConstructInterface, Modify_Interface

   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !PorousMediaProperties
    public  :: ConstructPorousMediaProperties
    private ::      AllocateInstance

    !Selector
    public  :: GetConcentration
    public  :: SetVegetationPMProperties
!    public  :: GetNextPorousMediaPropDT
    public  :: GetPMPCoupled
    public  :: SetWindVelocity                 
    public  :: UnGetPorousMediaProperties
    
    !Modifier
    public  :: ModifyPorousMediaProperties
    private ::      AdvectionDiffusionProcesses
    private ::          ModifyAdvectionDiffusion
    private ::      SoilQualityProcesses
    private ::      SoilChemistryProcesses
    
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
!        module procedure UnGetPorousMediaProperties3D_R8i
    end interface  UnGetPorousMediaProperties

    !Types---------------------------------------------------------------------
    
    private :: T_PorousMediaProperties
    

    !PropI
!    character(LEN = StringLength), parameter :: char_nitrite        = trim(adjustl('nitrite'            ))

    !PropII
!    character(LEN = StringLength), parameter :: char_nitrate      = trim(adjustl('nitrate'))                        
    character(LEN = StringLength), parameter    :: prop_block_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter    :: prop_block_end       = '<endproperty>'

    type       T_ID
        integer                                 :: IDNumber
        character(LEN = StringLength)           :: name
        character(LEN = StringLength)           :: description
        character(LEN = StringLength)           :: units
    end type T_ID

    type       T_Property_3D
         type(T_PropertyID)                     :: ID
         real, pointer, dimension (:,:,:)       :: Field
         real                                   :: Scalar
    end type   T_Property_3D


    type T_ExtVar
        !from basin
!        real   , dimension(:,:), pointer            :: CropEvapotrans       => null()
!        real   , dimension(:,:), pointer            :: Transpiration        => null() ! Lúcia nova var
!        real   , dimension(:,:), pointer            :: Evaporation          => null() ! Lúcia nova var
!        real(8), dimension(:,:), pointer            :: PotentialInfCol      => null()
!        integer                                     :: EvapoTranspirationMethod
        real                                        :: PorousMediapropDT
        type(T_Time)                                :: Now
        type(T_Time)                                :: BeginTime
        type(T_Time)                                :: EndTime
   
        ! from porousMedia
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        integer, pointer, dimension(:,:,:)          :: OpenPoints3D
        integer, dimension(:,:), pointer            :: BasinPoints
        real,    pointer, dimension(:,:,:)          :: WaterContentOld
        real,    dimension(:,:,:), pointer          :: UnSatW
        real,    pointer, dimension(:,:,:)          :: WaterContent
        real(8), pointer, dimension(:,:,:)          :: CellVolume
        real(8), dimension(:,:,:), pointer          :: FluxU
        real(8), dimension(:,:,:), pointer          :: FluxV
        real(8), dimension(:,:,:), pointer          :: FluxW
        real   , pointer, dimension(:,:  )          :: Area
        real                                        :: DT
        real   , pointer, dimension(:,:,:)          :: DWZ
        integer, pointer, dimension(:,:,:)          :: ComputeFacesW3D
        real   , pointer, dimension(:,:  )          :: Topography  
        real ,   pointer, dimension(:,:,:)          :: SZZ     

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
    integer                                         :: SpatialMethod
    real,    pointer, dimension(:,:,:)              :: DifusionNumber
    real,    pointer, dimension(:,:,:)              :: ReynoldsMNumber                   
    end type T_AdvectionDiffusion

    type       T_Partition                      
        logical                                 :: NonLinear
        character(LEN = StringLength)           :: NonLinear_ks_Units
        type(T_Property_3D)                     :: Nu            
        type(T_Property_3D)                     :: Be          
        type(T_Property_3D)                     :: ks
        type(T_Property_3D)                     :: PartitionRate
        type(T_Property_3D)                     :: Fraction 
        character (LEN = StringLength)          :: Partition_Couple
    end type T_Partition

    type       T_Evolution
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
    end type T_Files    

    type T_Property
        type (T_PropertyID)                     :: ID
        real, dimension(:,:,:), pointer         :: Concentration            => null()
        real, dimension(:,:,:), pointer         :: ConcentrationIni         => null()
        real, dimension(:,:,:), pointer         :: ConcentrationOld         => null()
        real, dimension(:,:  ), pointer         :: UpperConcentration       => null()
        type (T_Property), pointer              :: Next, Prev                     => null()
        logical                                 :: Particulate
        type (T_Evolution)                      :: Evolution
        real, pointer, dimension(:,:)           :: SedPropAdvFlux
        logical                                 :: Old     = .false.
        real                                    :: MinValue        = FillValueReal
        real, pointer, dimension(:,:,:)         :: Mass_Created
        logical                                 :: TimeSerie        = .false.
        logical                                 :: BoxTimeSerie     = .false.
        logical                                 :: BoxTimeSerie2D   = .false.
        logical                                 :: OutputHDF        = .false.
        real                                    :: DifCoef
        real                                    :: Dispersivity
        real                                    :: RainConc
        real                                    :: BottomConc
    end type T_Property

    type T_Coupled
        logical                                 :: SoilQuality          = .false. !Sediment source/sink model (Sediment Quality)
        real                                    :: SoilQuality_DT
        type (T_Time)                           :: SoilQuality_NextCompute
        logical                                 :: SoilChemistry        = .false.  !Sediment reactions model (PHREEQC)
        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: MinConcentration     = .false.
    end type T_Coupled

    type       T_PorousMediaProperties
        integer                                     :: ObjTime              = 0
!        integer                                     :: ObjTopography        = 0  
        integer                                     :: ObjHorizontalGrid    = 0
!        integer                                     :: ObjHorizontalMap     = 0
!        integer                                     :: ObjDrainageNetwork   = 0
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
        type (T_ExtVar)                             :: ExtVar
        logical                                     :: CheckGlobalMass      
        type (T_Files)                              :: Files
        type (T_AdvectionDiffusion)                 :: AdvDiff
        type (T_OutPut)                             :: OutPut
        type (T_Property), pointer                  :: FirstProperty    => null() !Lúcia
        type (T_Property), pointer                  :: LastProperty        
        type (T_PorousMediaProperties), pointer     :: Next             => null() !Lúcia
        type (T_Coupled)                            :: Coupled
        type (T_Time)                               :: LastOutputHDF5

        real,    pointer, dimension(:,:,:)          :: PropI
        real,    pointer, dimension(:,:,:)          :: PropII
        real,    pointer, dimension(:,:,:)          :: PropInew 
        logical                                     :: PorousMediaProperties
        real,    pointer, dimension(:,:,:)          :: Volume   
        integer                                     :: PropertiesNumber    = 0
        real   , pointer, dimension(:,:,:)          :: DissolvedToParticulate3D
        real                                        :: ResidualTime
!        real   , pointer, dimension(:,:,:)          :: SoilDryDensity
!        real   , pointer, dimension(:,:,:)          :: Salinity
!        real   , pointer, dimension(:,:,:)          :: pH
!        real   , pointer, dimension(:,:,:)          :: IonicStrength
!        real   , pointer, dimension(:,:,:)          :: PhosphorusAdsortionIndex
        
        integer                                     :: InstanceID
        type (T_Size3D)                             :: Size, WorkSize
        type (T_Size2D)                             :: Size2D

        real(8), dimension(:, :, :),  pointer       :: Matrix

    end type  T_PorousMediaProperties

    !Global Module Variables
    type (T_PorousMediaProperties), pointer                         :: FirstObjPorousMediaProperties
    type (T_PorousMediaProperties), pointer                         :: Me

!    integer                                         :: mPorousMediaProperties_ = 0 !just to compile

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
!                                              HorizontalMapID,                            &
!                                              TopographyID,                               &
                                              BasinGeometryID,                            &
!                                              DrainageNetworkID,                          &
                                              PorousMediaID,                              &
                                              GeometryID,                                 &
                                              MapID,                                      &
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
!        integer                                         :: HorizontalMapID
!        integer                                         :: TopographyID
        integer                                         :: BasinGeometryID
!        integer                                         :: DrainageNetworkID
        integer                                         :: PorousMediaID
        integer                                         :: GeometryID
        integer                                         :: MapID
        integer, optional, intent(OUT)                  :: STAT 
        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_,STAT_CALL
!        integer, pointer, dimension(:,:,:)              :: WaterPoints
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
!            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       TopographyID    ) 
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
!            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMap_,            MapID           )
        
!            if (DrainageNetworkID /= 0)                                                     &
!                Me%ObjDrainageNetwork  = AssociateInstance (MDRAINAGENETWORK_,  DrainageNetworkID )
            

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

    !       call ConstructProfileOutput   em teste
            
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR010'



            !Couple nutrient, carbon and oxygen sources and sinks model
            if (Me%Coupled%SoilQuality) then
                call CoupleSoilQuality
            endif
            
            !Couple soil chemical model
            if (Me%Coupled%SoilChemistry) then
                call CoupleSoilChemistry
            endif

            

            !Returns ID
            ObjPorousMediaPropertiesID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR040' 

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
        !Begin-----------------------------------------------------------------

        !Geometry Size
        call GetGeometrySize    (Me%ObjGeometry,             &    
                                 Size     = Me%Size,         &
                                 WorkSize = Me%WorkSize,     &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR03'

        Me%Size2D%ILB = Me%Size%ILB
        Me%Size2D%IUB = Me%Size%IUB
        Me%Size2D%JLB = Me%Size%JLB
        Me%Size2D%JUB = Me%Size%JUB

        call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR004'

        call GetComputeTimeLimits(Me%ObjTime,                      &
                                  EndTime   = Me%ExtVar%EndTime,   &
                                  BeginTime = Me%ExtVar%BeginTime, &
                                  STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ReadGlobalOptions - ModulePorousMediaProperties - ERR005'

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

        !Needs keyword description and value options definition
        call GetData(Me%AdvDiff%SpatialMethod,                                            &   !Lúcia
                     Me%ObjEnterData, iflag,                                              &
                     SearchType = FromFile,                                               &
                     keyword    = 'SPATIAL_METHOD',                                       &
                     Default    = 0,                                                      &                                           
                     ClientModule ='ModulePorousMediaProperties',                         &
                     STAT       = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModulePorousMediaProeprties - ERR06'

 
    
    end subroutine ReadGlobalOptions        

    !--------------------------------------------------------------------------

    subroutine AllocateVariables        
        
        !Local-----------------------------------------------------------------        
        integer                                         :: ILB, IUB, JLB,  JUB 
        integer                                         :: KLB, KUB 

        !Bounds
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KUB = Me%Size%KUB
        KLB = Me%Size%KLB
        
               
        !Water Content---------------------------------------------------------
        allocate (Me%PropI                   (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%PropInew                (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%PropII                  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%Volume                  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%AdvDiff%DifusionNumber  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%AdvDiff%ReynoldsMNumber (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (Me%ExtVar%WindVelocity3D   (ILB:IUB,JLB:JUB,KLB:KUB))

!        if (Me%ExtVar%VegetationCoupled) then
!            allocate (Me%ExtVar%SoilFluxesActive               (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%GrazingBiomass                 (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%GrazingNitrogen                (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%GrazingPhosphorus              (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%ManagementAerialBiomass        (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%ManagementNitrogen             (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%ManagementPhosphorus           (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%ManagementRootBiomass          (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%DormancyBiomass                (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%DormancyNitrogen               (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%DormancyPhosphorus             (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilNitrateSurface           (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilNitrateSubSurface        (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilAmmoniaSurface           (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilAmmoniaSubSurface        (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilOrganicNSurface          (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilOrganicNSubSurface       (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilOrganicPSurface          (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilOrganicPSubSurface       (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilMineralPSurface          (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%FertilMineralPSubSurface       (ILB:IUB,JLB:JUB))
!            allocate (Me%ExtVar%NitrogenUptake         (ILB:IUB,JLB:JUB,KLB:KUB))
!            allocate (Me%ExtVar%PhosphorusUptake       (ILB:IUB,JLB:JUB,KLB:KUB))
!        endif

    endsubroutine AllocateVariables

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
        
        nullify(NewProperty%Prev,NewProperty%Next)
        nullify(NewProperty%Concentration        )

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
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR01'

        if (NewProperty%Evolution%AdvectionDiffusion) then
            Me%Coupled%AdvectionDiffusion = .true.
            NewProperty%Evolution%Variable = .true.
        endif
        
        !The next keywords are not used so will be commented for now. Possibly to remove                    
!        call GetData(NewProperty%Evolution%Partitioning,                                 &
!                     Me%ObjEnterData, iflag,                                             &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'PARTITION',                                         &
!                     ClientModule = 'ModulePorousMediaProperties',                       &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR010'
!
!        if (NewProperty%Evolution%Partitioning) then
!            NewProperty%Evolution%Variable = .true.
!        endif
!
!        !Partition parameters  
!        if (NewProperty%Evolution%Partitioning)                                          &
!          call Read_Partition_Parameters (NewProperty, ClientNumber)
!
!        !<BeginKeyword>
!            !Keyword          : CATION_EXCHANGE
!            !<BeginDescription>       
!               ! Property has cation exchange as sink and source
!            !<EndDescription>
!            !Type             : Boolean 
!            !Default          : .false.
!            !File keyword     : SEDPROP
!            !Multiple Options : 1 (.true.), 0 (.false.)
!            !Search Type      : FromBlock
!            !Begin Block      : <beginproperty>
!            !End Block        : <endproperty>
!        !<EndKeyword>
!        call GetData(NewProperty%Evolution%CationExchangeProcess,                        &
!                     Me%ObjEnterData, iflag,                                             &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'CATION_EXCHANGE',                                   &
!                     ClientModule = 'ModulePorousMediaProperties',                       &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR020'
!
!        if (NewProperty%Evolution%CationExchangeProcess) then
!            NewProperty%Evolution%Variable = .true.
!        endif
!
!        !<BeginKeyword>
!            !Keyword          : CHEMICAL_EQUILIBRIUM
!            !<BeginDescription>       
!               ! Property has chemical equilibrium as sink and source
!            !<EndDescription>
!            !Type             : Boolean 
!            !Default          : .false.
!            !File keyword     : SEDPROP
!            !Multiple Options : 1 (.true.), 0 (.false.)
!            !Search Type      : FromBlock
!            !Begin Block      : <beginproperty>
!            !End Block        : <endproperty>
!        !<EndKeyword>
!        call GetData(NewProperty%Evolution%ChemEquilibriumProcess,                       &
!                     Me%ObjEnterData, iflag,                                             &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'CHEMICAL_EQUILIBRIUM',                              &
!                     ClientModule = 'ModulePorousMediaProperties',                       &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR030'
!
!        if (NewProperty%Evolution%ChemEquilibriumProcess) then 
!            NewProperty%Evolution%Variable = .true.
!        endif


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
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR040'

        if (NewProperty%Evolution%SoilQuality) then
            Me%Coupled%SoilQuality     = .true.
            NewProperty%Evolution%Variable = .true.
        endif


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
            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR050'

        if (NewProperty%Evolution%SoilChemistry) then
            Me%Coupled%SoilChemistry       = .true.
            NewProperty%Evolution%Variable = .true.
        endif


        !The next keywords are not used so will be commented for now. Possibly to remove                    
!        !<BeginKeyword>
!            !Keyword          : Soil_WATER_FLUXES
!            !<BeginDescription>       
!               !  This property has fluxes at the Soil water interface? no - 0;  yes - 1
!            !<EndDescription>
!            !Type             : Boolean 
!            !Default          : .false.
!            !File keyword     : SEDPROP
!            !Multiple Options : 1 (.true.), 0 (.false.)
!            !Search Type      : FromBlock
!            !Begin Block      : <beginproperty>
!            !End Block        : <endproperty>
!        !<EndKeyword>
!
!        call GetData(NewProperty%Evolution%SoilWaterFluxes,                              &
!                     Me%ObjEnterData,iflag,                                              &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'SOIL_WATER_FLUXES',                                 &
!                     ClientModule = 'ModulePorousMediaProperties',                       &
!                     Default      = .false.,                                             &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR060'
!        
!        if (NewProperty%Evolution%SoilWaterFluxes) then
!            NewProperty%Evolution%Variable = .true.
!        endif
!
!        !<BeginKeyword>
!            !Keyword          : MACROPORES
!            !<BeginDescription>       
!               !  This property has fluxes with macropores? no - 0;  yes - 1
!            !<EndDescription>
!            !Type             : Boolean 
!            !Default          : .false.
!            !File keyword     : SEDPROP
!            !Multiple Options : 1 (.true.), 0 (.false.)
!            !Search Type      : FromBlock
!            !Begin Block      : <beginproperty>
!            !End Block        : <endproperty>
!        !<EndKeyword>
!
!        call GetData(NewProperty%Evolution%Macropores,                                   &
!                     Me%ObjEnterData,iflag,                                              &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'MACROPORES',                                        &
!                     ClientModule = 'ModulePorousMediaProperties',                       &
!                     Default      = .false.,                                             &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR070'
!        
!            ! Quanto muito, constroi-se uma função que leia os parametros dentyro das mesmas características da propriedade 
!        
! !       if (NewProperty%Evolution%AdvectionDiffusion)then   
!
!  !          call Read_Advec_Difus_Parameters    (NewProperty)
!
!   !         call Construct_Property_Diffusivity (NewProperty)
!        
!    !    end if
!
        
        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR0100'

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModulePorousMediaProperties',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR0110'
                                       
            
            if (NewProperty%Evolution%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR0120'

            elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR0130'
                endif

                !Run period in seconds
                DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%Evolution%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'Construct_PropertyEvolution - ModulePorousMediaProperties - ERR140'
                end if
            endif

            NewProperty%Evolution%NextCompute = Me%ExtVar%Now + NewProperty%Evolution%DTInterval

        else

            call null_time(NewProperty%Evolution%NextCompute)

            NewProperty%Evolution%DTInterval = FillValueReal

        endif

    end subroutine Construct_PropertyEvolution     


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
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR01'
        NewProperty%Concentration(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR010'
        NewProperty%Concentration(:,:,:) = FillValueReal


        allocate(NewProperty%SedPropAdvFlux(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR020'
        NewProperty%SedPropAdvFlux(:,:) = 0.0
        


        call GetData(NewProperty%DifCoef,                                                   &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'DIFFUSION_COEF',                                       &
                     Default      = .1,                                                     &                        
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR030'

        call GetData(NewProperty%Dispersivity,                                              &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'DISPERSIVITY',                                         &
                     Default      = .1,                                                     &                        
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR040'

        call GetData(NewProperty%RainConc,                                                  &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'RAIN_CONC',                                            &
                     Default      = .1,                                                     &                        
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR050'
        
        call GetData(NewProperty%BottomConc,                                                &
                     Me%ObjEnterData, iflag,                                                &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'BOTTOM_CONC',                                          &
                     Default      = .1,                                                     &                        
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR060'

        call GetData(NewProperty%MinValue,                                                  &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'MIN_VALUE',                                            &
                     ClientModule = 'ModuleSoilProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR070'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
            Me%Coupled%MinConcentration = .true.
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR080'
            NewProperty%Mass_Created(:,:,:) = 0.
        end if

        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleVegetation',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR85'
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not.NewProperty%Old) then

            !Get water points
            call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR86'

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
                stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR090'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR0100'
            end if


!Lúcia
            call SetMatrixValue(NewProperty%ConcentrationOld, Me%Size, NewProperty%Concentration,Me%ExtVar%WaterPoints3D)

            call CheckFieldConsistence (NewProperty)

            call UnGetMap(Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModulePorousMediaProperties - ERR110'

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
        
        !Counts the number of Properties which has timeserie option set to true
        PropertyX => Me%FirstProperty
        nProperties = 0
        do while (associated(PropertyX))
            if (PropertyX%TimeSerie) then
                nProperties = nProperties + 1
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
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR03'

        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR03.5'
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR03.6'  

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR03.7'



        !Writes the Grid
        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Topography", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR04'


        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR07'

        !Water Points
        call HDF5WriteData   ( Me%ObjHDF5,  "/Grid", "WaterPoints3D", "-",                  &
                               Array3D = Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR08'
              
                
        !Flushes All pending HDF5 commands
        call HDF5FlushMemory    (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR10'



        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR085'

        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR12'  

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModulePorousMediaProperties - ERR13'


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

    subroutine CoupleSoilChemistry        

        !Local-----------------------------------------------------------------
        !Begin-----------------------------------------------------------------
    
    end subroutine CoupleSoilChemistry

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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
!    subroutine GetNextPorousMediaPropDT (ObjPorousMediaPropertiesID, PorousMediaPropDT, STAT)
!
!        !Arguments--------------------------------------------------------------
!        integer                                         :: ObjPorousMediaID
!        real, intent(OUT)                               :: PorousMEdiaPropDT
!        integer, intent(OUT), optional                  :: STAT
!
!        !Local------------------------------------------------------------------
!        integer                                         :: STAT_CALL, ready_
!
!        !-----------------------------------------------------------------------
!
!        STAT_CALL = UNKNOWN_
!
!        call Ready(ObjPorousMediaID, ready_)
!
!        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
!
!            PorousMediaPropDT        = Me%ExtVar%PorousMediaPropDT
!
!            STAT_CALL = SUCCESS_
!        else 
!            STAT_CALL = ready_
!        end if
!
!        if (present(STAT)) STAT = STAT_CALL
!
!    end subroutine GetNextPorousMediapropDT

    !---------------------------------------------------------------------------

    subroutine GetPMPCoupled(PorousMediaPropertiesID,           &
                             SoilQuality,                       &
                             SoilChemistry,                     &
                             STAT) 

        !Arguments-------------------------------------------------------------
        integer                        :: PorousMediaPropertiesID
        integer, optional, intent(OUT) :: STAT
        logical, optional, intent(OUT) :: SoilQuality, SoilChemistry

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
            if (present(SoilChemistry )) SoilChemistry  = Me%Coupled%SoilChemistry


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetPMPCoupled
    !--------------------------------------------------------------------------
    
    subroutine GetConcentration(PorousMediaPropertiesID, ConcentrationX, PropertyXIDNumber, &
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
            
    end subroutine GetConcentration

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
            
            !Nullification in the case that not defined (not present)
!            Me%ExtVar%SoilFluxesActive (:,:)   =    .false. 
!            call SetMatrixValue(Me%ExtVar%GrazingBiomass          , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%GrazingNitrogen         , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%GrazingPhosphorus       , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%ManagementAerialBiomass , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%ManagementNitrogen      , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%ManagementPhosphorus    , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%ManagementRootBiomass   , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%DormancyBiomass         , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%DormancyNitrogen        , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%DormancyPhosphorus      , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilNitrateSurface    , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilNitrateSubSurface , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilAmmoniaSurface    , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilAmmoniaSubSurface , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilOrganicNSurface   , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilOrganicNSubSurface, Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilOrganicPSurface   , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilOrganicPSubSurface, Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilMineralPSurface   , Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%FertilMineralPSubSurface, Me%Size2D,   0.0 , Me%ExtVar%BasinPoints)
!            call SetMatrixValue(Me%ExtVar%NitrogenUptake          , Me%Size  ,   0.0 , Me%ExtVar%WaterPoints3D)
!            call SetMatrixValue(Me%ExtVar%PhosphorusUptake        , Me%Size  ,   0.0 , Me%ExtVar%WaterPoints3D)


!            Me%ExtVar%GrazingBiomass              (:,:) = 0.0 
!            Me%ExtVar%GrazingNitrogen             (:,:) = 0.0
!            Me%ExtVar%GrazingPhosphorus           (:,:) = 0.0
!            Me%ExtVar%ManagementAerialBiomass     (:,:) = 0.0
!            Me%ExtVar%ManagementNitrogen          (:,:) = 0.0
!            Me%ExtVar%ManagementPhosphorus        (:,:) = 0.0
!            Me%ExtVar%ManagementRootBiomass       (:,:) = 0.0
!            Me%ExtVar%DormancyBiomass             (:,:) = 0.0
!            Me%ExtVar%DormancyNitrogen            (:,:) = 0.0
!            Me%ExtVar%DormancyPhosphorus          (:,:) = 0.0
!            Me%ExtVar%FertilNitrateSurface        (:,:) = 0.0
!            Me%ExtVar%FertilNitrateSubSurface     (:,:) = 0.0
!            Me%ExtVar%FertilAmmoniaSurface        (:,:) = 0.0
!            Me%ExtVar%FertilAmmoniaSubSurface     (:,:) = 0.0
!            Me%ExtVar%FertilOrganicNSurface       (:,:) = 0.0
!            Me%ExtVar%FertilOrganicNSubSurface    (:,:) = 0.0
!            Me%ExtVar%FertilOrganicPSurface       (:,:) = 0.0
!            Me%ExtVar%FertilOrganicPSubSurface    (:,:) = 0.0
!            Me%ExtVar%FertilMineralPSurface       (:,:) = 0.0
!            Me%ExtVar%FertilMineralPSubSurface    (:,:) = 0.0
!
!            Me%ExtVar%NitrogenUptake            (:,:,:) = 0.0
!            Me%ExtVar%PhosphorusUptake          (:,:,:) = 0.0

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
!                                           InfiltrationColumn,                  &
!                                          EvapotranspirationMethod,            &
!                                           PotentialEVTP,                       &
!                                           Evaporation,                         &
!                                           Transpiration,                       &
                                           STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjPorousMediaPropertiesID
        integer, intent(OUT), optional              :: STAT
!        real(8), dimension (:,:), pointer           :: InfiltrationColumn
!        real, dimension(:, :), pointer              :: PotentialEVTP, Evaporation, Transpiration 
!        integer                                     :: EvapotranspirationMethod

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_,STAT_CALL
!        real                                        :: PorousMediaDT
!        real(8), dimension(:, :), pointer           :: Infiltration
!        real(8), dimension(:, :), pointer           :: EfectiveEVTP,EfectiveEVTP2,plantwaterstress


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
            
            !Nutrient sources and sinks from vegetation
            call InterfaceFluxes

            if (Me%Coupled%AdvectionDiffusion) then
                call AdvectionDiffusionProcesses
            endif

            if (Me%Coupled%SoilQuality) then
                call SoilQualityProcesses
            endif

            if (Me%Coupled%SoilChemistry) then
                call SoilChemistryProcesses
            endif

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

    subroutine InterfaceFluxes
        !Local--------------------------------------------------------------------
        !Begin--------------------------------------------------------------------

        if (Me%ExtVar%CoupledVegetation) then
            if (Me%ExtVar%ComputeVegInterfaceFluxes) then
                call VegetationInterfaceFluxes
            endif
        endif

!        if (Me%ExtVar%CoupledRunoff) then
!            if (Me%ExtVar%ComputeRunoffInterfaceFluxes) then
!                call RunoffInterfaceFluxes
!            endif
!        endif


    end subroutine InterfaceFluxes
 
    !-----------------------------------------------------------------------------

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
                    
                                !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                                FertilizationNitrate  = Me%ExtVar%FertilNitrateSubSurface(i,j) * 1e9 * Area / 10000.
                                FertilizationAmmonia  = Me%ExtVar%FertilAmmoniaSubSurface(i,j) * 1e9 * Area / 10000.
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
                    
                    !      ugN       = KgN/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                    NitrogenUptake = Me%ExtVar%NitrogenUptake(i,j,k) * 1e9 * Area / 10000.

                endif

                if (Me%ExtVar%ModelPhosphorus) then
                    
                    !      ugP       = KgP/ha * 1E9ug/kg * (m2) * 1ha/10000m2                     
                    PhosphorusUptake = Me%ExtVar%PhosphorusUptake(i,j,k) * 1e9 * Area / 10000.

                endif


                !s
                ModelDT         = Me%ExtVar%DT
                VegDT           = Me%ExtVar%VegetationDT
                !m3             = m3H20/m3cell * m3cell
                CellWaterVolume = Me%ExtVar%WaterContent(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 

                !Soil mass to compute organic and microorganisms pools
                call SearchProperty(Property, SoilDryDensity_        , .false., STAT = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
                !kgsoil         = kg/m3  * m3cell
                CellSoilMass    = Property%Concentration(i,j,k) * Me%ExtVar%CellVolume(i,j,k) 
                
                if (Me%ExtVar%GrowthModel) then

                    ! Property Calculation
                    !!Carbon
                    call SearchProperty (Property, RefreactaryOrganicC_        , .false., STAT = STAT_CALL)    
                    if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR80'
                
                    !         ug/kgsoil            = ug/kgsoil + ug / kgsoil
                    Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((GrazingCarbon + DormancyCarbon       &
                                                     + ManagementCarbon + ManagementRootCarbon) * ModelDT / VegDT)            &
                                                     / CellSoilMass)

    !                call SearchProperty (Property, LabileOrganicC_        , .false., STAT = STAT_CALL)    
    !                if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR90'
    !
    !                Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((FertilizationOrganicC)               &
    !                                                * ModelDT / VegDT) / CellSoilMass)
                endif

                !!Nitrogen
                if (Me%ExtVar%ModelNitrogen) then
                    
                    if (Me%ExtVar%GrowthModel) then
                        call SearchProperty (Property, RefreactaryOrganicN_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR100'
                    
                        !         ug/kgsoil 
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((GrazingNitrogen + DormancyNitrogen       &
                                                        + ManagementNitrogen + ManagementRootNitrogen) * ModelDT / VegDT)             &
                                                        / CellSoilMass)


                        call SearchProperty (Property, PON_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR110'
                    
                        !         ug/kgsoil 
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((FertilizationOrganicN)                   &
                                                        * ModelDT / VegDT) / CellSoilMass)


                        call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'

                        !         ug/m3                = ug/m3 + ug / m3H20
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((FertilizationNitrate - NitrogenUptake)    &
                                                        * ModelDT / VegDT) / CellWaterVolume)


                        call SearchProperty (Property, Ammonia_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR130'

                        !         ug/m3 
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((FertilizationAmmonia)                    &
                                                        * ModelDT / VegDT) / CellWaterVolume)
                    else

                        call SearchProperty (Property, Nitrate_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR120'

                        !         ug/m3 
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) - (((NitrogenUptake)                          &
                                                        * ModelDT / VegDT) / CellWaterVolume)
                    endif

                endif
                
                !Phosphorus
                if (Me%ExtVar%ModelPhosphorus) then

                    if (Me%ExtVar%GrowthModel) then

                        call SearchProperty (Property, RefreactaryOrganicP_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR140'

                        !         ug/kgsoil  
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((GrazingPhosphorus + DormancyPhosphorus       &
                                                        + ManagementPhosphorus + ManagementRootNitrogen) * ModelDT / VegDT)               &
                                                        / CellSoilMass)


                        call SearchProperty (Property, POP_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR150'

                        !         ug/kgsoil  
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((FertilizationOrganicP)                       &
                                                        * ModelDT / VegDT) / CellSoilMass)


                        call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'

                        !         ug/m3 
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) + (((FertilizationMineralP - PhosphorusUptake)   &
                                                        * ModelDT / VegDT) / CellWaterVolume)   
                    else
                        
                        call SearchProperty (Property, Inorganic_Phosphorus_        , .false., STAT = STAT_CALL)    
                        if (STAT_CALL /= SUCCESS_) stop 'InterfaceFluxes - ModulePorousMediaProperties - ERR160'

                        !         ug/m3 
                        Property%Concentration (i,j,k) = Property%Concentration (i,j,k) - (((PhosphorusUptake)                         &
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


    subroutine AdvectionDiffusionProcesses

        !Local--------------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX
        
        !begin--------------------------------------------------------------------

        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%Evolution%AdvectionDiffusion) then

                 call ModifyAdvectionDiffusion(PropertyX)

            endif


            PropertyX => PropertyX%Next

        enddo

    end subroutine AdvectionDiffusionProcesses
    
    !-----------------------------------------------------------------------------

    subroutine ModifyAdvectionDiffusion (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k!, CHUNK
!        real                                        :: OldMass, NewMass
!        real                                        :: QbeforeI,QafterI   
!        real                                        :: QbeforeJ,QafterJ   
!        real                                        :: QbeforeK,QafterK
!        real                                        :: CbeforeI,CafterI
!        real                                        :: CbeforeJ,CafterJ
!        real                                        :: CbeforeK,CafterK        
        real                                        :: Area!, dif 
        real                                        :: aux, cofA,cofB,cofC
        real                                        :: difbefore,difafter
        real                                        :: WaterContentBefore,WaterContentAfter
        real                                        :: CO,CKmax
        real(8), pointer, dimension(:,:,:)          :: FluxW
        real   , pointer, dimension(:,:,:)          :: DWZ, Theta, ThetaOld
        !Begin-----------------------------------------------------------------
   
        !!CHUNK = CHUNK_K(Me%WorkSize%KLB, Me%WorkSize%KUB)

        CurrProperty => PropertyX

        !!!$OMP PARALLEL PRIVATE(I,J,K)
        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%WaterPoints3D(I,J,K) == WaterPoint) then
                
                Area      = Me%ExtVar%Area(i, j)
                
                CO        = CurrProperty%BottomConc
                CKmax     = CurrProperty%RainConc
                
                FluxW     => Me%ExtVar%FluxW
                DWZ       => Me%ExtVar%DWZ
                Theta     => Me%ExtVar%WaterContent
                ThetaOld  => Me%ExtVar%WaterContentOld

                Me%Volume(i,j,k)= Theta(i,j,k)*Me%ExtVar%Cellvolume(i,j,k)
            

                if (Me%AdvDiff%SpatialMethod==2) then ! diferenças centrais


                    aux      = (Me%extvar%DT/((Me%extvar%cellvolume(i,j,k))*Theta(i,j,k)))


                    if (K==Me%WorkSize%KUB) then

                        WaterContentBefore = min(ThetaOld(i,j,k),ThetaOld(i,j,k-1))
                        WaterContentAfter  = min(ThetaOld(i,j,k),ThetaOld(i,j,k))

                        difbefore  = CurrProperty%DifCoef*Tortuosity(WaterContentBefore)                                               &
                                     +(abs(FluxW(i,j,k)/Area)*CurrProperty%Dispersivity)/WaterContentBefore
                        difafter   = CurrProperty%DifCoef*Tortuosity(WaterContentAfter)                                                &
                                     +(abs(FluxW(i,j,k)/Area)*CurrProperty%Dispersivity)/WaterContentAfter

                        cofB= ThetaOld(i,j,k)/Theta(i,j,k)                                                                             &
                              + (aux*FluxW(i,j,k)*DWZ(i,j,k-1)/(DWZ(i,j,k-1)+DWZ(i,j,k)))                                              &
                              - difbefore*Area*aux/(0.5*(DWZ(i,j,k)+DWZ(i,j,k-1)))                                                 

                        cofA= (aux*FluxW(i,j,k)*DWZ(i,j,k)/(DWZ(i,j,k)+DWZ(i,j,k-1)))                                                  &
                              +((aux*difbefore*Area)/(0.5*(DWZ(i,j,k)+DWZ(i,j,k-1)))) 

                        cofC=(+aux*Me%ExtVar%UnsatW(i,j,k+1)) ! Velocity signal

                        CurrProperty%Concentration(i,j,k)= cofA*CurrProperty%ConcentrationOld(i,j,k-1)                                 &
                                                           +cofB*CurrProperty%ConcentrationOld(i,j,k)                                  &
                                                           +cofC*CKmax

                        Me%AdvDiff%DifusionNumber(i,j,k)= cofB

                        Me%AdvDiff%ReynoldsMNumber(i,j,k)= cofA

                    elseif (K==1)       then
            
                        WaterContentBefore = min(ThetaOld(i,j,k),ThetaOld(i,j,k))
                        WaterContentAfter  = min(ThetaOld(i,j,k),ThetaOld(i,j,k+1))

                        difbefore  = CurrProperty%DifCoef*Tortuosity(WaterContentBefore)                                               &
                                     +(abs(FluxW(i,j,k))*CurrProperty%Dispersivity)/WaterContentBefore
                        difafter   = CurrProperty%DifCoef*Tortuosity(WaterContentAfter)                                                &
                                     +(abs(FluxW(i,j,k+1))*CurrProperty%Dispersivity)/WaterContentAfter
                
                        cofB= ThetaOld(i,j,k)/Theta(i,j,k)                                                                             &
                              - (aux*FluxW(i,j,k+1)*DWZ(i,j,k+1)/(DWZ(i,j,k+1)+DWZ(i,j,k)) )                                           &
                              - difafter*Area*aux/(0.5*(DWZ(i,j,k)+DWZ(i,j,k+1)))                                                      &                     
                              + aux*FluxW(i,j,k)                                                 

    !                   cofA = aux*FluxW(i,j,k)
                        cofA=0     

                        cofC =  -(aux*FluxW(i,j,k+1)*DWZ(i,j,k)/(DWZ(i,j,k)+DWZ(i,j,k+1)))                                             &
                                +(difafter*Area*aux)/(0.5*(DWZ(i,j,k)+DWZ(i,j,k+1)))


                        CurrProperty%Concentration(i,j,k)= cofA*CurrProperty%ConcentrationOld(i,j,k-1)                                 &
                                                           +cofB*CurrProperty%ConcentrationOld(i,j,k)                                  &
                                                           +cofC*CurrProperty%ConcentrationOld(i,j,k+1)
            
                        Me%AdvDiff%DifusionNumber(i,j,k)= cofB

                        Me%AdvDiff%ReynoldsMNumber(i,j,k)= cofA
      
                    else

                        WaterContentBefore = min(ThetaOld(i,j,k),ThetaOld(i,j,k-1))
                        WaterContentAfter  = min(ThetaOld(i,j,k),ThetaOld(i,j,k+1))

                        difbefore  = CurrProperty%DifCoef*Tortuosity(WaterContentBefore)                                               &
                                     +(abs(FluxW(i,j,k))*CurrProperty%Dispersivity)/WaterContentBefore
                        difafter   = CurrProperty%DifCoef*Tortuosity(WaterContentAfter)                                                &
                                     +(abs(FluxW(i,j,k+1))*CurrProperty%Dispersivity)/WaterContentAfter

            

                        cofB= ThetaOld(i,j,k)/Theta(i,j,k)                                                                             &
                             - (aux*FluxW(i,j,k+1)*DWZ(i,j,k+1)/(DWZ(i,j,k+1)+DWZ(i,j,k)))                                             &
                             + (aux*FluxW(i,j,k)*DWZ(i,j,k-1)/(DWZ(i,j,k-1)+DWZ(i,j,k)))                                               &
                             - difbefore*Area*aux/(0.5*(DWZ(i,j,k)+DWZ(i,j,k-1)))                                                      &
                             - difafter*Area*aux/(0.5*(DWZ(i,j,k)+DWZ(i,j,k+1)))

                        cofA = (aux*FluxW(i,j,k)*DWZ(i,j,k)/(DWZ(i,j,k-1)+DWZ(i,j,k)))                                                 &
                              +((aux*difbefore*Area)/(0.5*(DWZ(i,j,k)+DWZ(i,j,k-1)))) 
                            

                        cofC =  -(aux*FluxW(i,j,k+1)*DWZ(i,j,k)/(DWZ(i,j,k)+DWZ(i,j,k+1)))                                             &
                                +(difafter*Area*aux)/(0.5*(DWZ(i,j,k)+DWZ(i,j,k+1)))


                        CurrProperty%Concentration(i,j,k)= cofA*CurrProperty%ConcentrationOld(i,j,k-1)                                 &
                                                           +cofB*CurrProperty%ConcentrationOld(i,j,k)                                  &
                                                           +cofC*CurrProperty%ConcentrationOld(i,j,k+1)
            
                        Me%AdvDiff%DifusionNumber(i,j,k)= cofB

                        Me%AdvDiff%ReynoldsMNumber(i,j,k)= cofA
                    
!                   If  (Me%propInew(i,j,k)<0) stop 'The is no numerical solution'
        
                            
                    endif

                endif
              
            endif
        enddo
        enddo
        enddo
        !!!$OMP END DO
        !!!$OMP END PARALLEL
                        
        call SetMatrixValue (CurrProperty%ConcentrationOld,      Me%Size,   CurrProperty%Concentration,          Me%ExtVar%WaterPoints3D)

        ! a concentração do tempo t tem de ser agora actualizada para o proximo ciclo
        !call SetMatrixValue (Me%PropI,      Me%Size,   Me%PropInew,          Me%ExtVar%WaterPoints3D)



    end subroutine ModifyAdvectionDiffusion
    
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
                    stop 'SoilQualityProcesses - ModulePorousMediaProperties - ERR01'
                

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
                        stop 'SoilQuality_Processes - ModulePorousMediaProperties - ERR02'

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

    subroutine SoilChemistryProcesses

        !Local--------------------------------------------------------------------
        !begin--------------------------------------------------------------------

    end subroutine SoilChemistryProcesses
    
    !-----------------------------------------------------------------------------    


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
                                                   + Property%evolution%DTInterval
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

    real function Tortuosity(WC)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: WC

        !local
        real                                        :: porosity


        porosity = 0.43

        !tortuosity = (WC**(10/3))/(porosity)
        tortuosity = (WC**(7/3))/(porosity**2)
        
        !Local-------------------------------------------------------------------


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

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPorousMediaProperties_,  Me%InstanceID)

            if (nUsers == 0) then

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - PorousmediaProperties - ERR05'
                endif

                !Deassociates External Instances
!                if (Me%ObjDrainageNetwork /= 0) then
!                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
!                    if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR06'
!                endif                
                
                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - PorousmediaProperties  - ERR08'
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR08'

!                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjTopography)
!                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR09'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR10'

!                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
!                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR11'
                
                nUsers = DeassociateInstance (mPOROUSMEDIA_,  Me%ObjPorousMedia)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR12'

                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR13'

                nUsers = DeassociateInstance (mMAP_,  Me%ObjMap)
                if (nUsers == 0) stop 'KillPorousMedia - PorousmediaProperties - ERR14'


!                call KillPorousMedia (Me%ObjPorousMedia, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR15'

                
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
        deallocate (Me%PropI                   )
        deallocate (Me%PropInew                )
        deallocate (Me%PropII                  )
        deallocate (Me%Volume                  )
        deallocate (Me%AdvDiff%DifusionNumber  )
        deallocate (Me%AdvDiff%ReynoldsMNumber )
        deallocate (Me%ExtVar%WindVelocity3D )

!        if (Me%ExtVar%VegetationCoupled) then
!            deallocate (Me%ExtVar%SoilFluxesActive          )
!            deallocate (Me%ExtVar%GrazingBiomass            )
!            deallocate (Me%ExtVar%GrazingNitrogen           )
!            deallocate (Me%ExtVar%GrazingPhosphorus         )
!            deallocate (Me%ExtVar%ManagementAerialBiomass   )
!            deallocate (Me%ExtVar%ManagementNitrogen        )
!            deallocate (Me%ExtVar%ManagementPhosphorus      )
!            deallocate (Me%ExtVar%ManagementRootBiomass     )
!            deallocate (Me%ExtVar%DormancyBiomass           )
!            deallocate (Me%ExtVar%DormancyNitrogen          )
!            deallocate (Me%ExtVar%DormancyPhosphorus        )
!           deallocate (Me%ExtVar%FertilNitrateSurface      )
!            deallocate (Me%ExtVar%FertilNitrateSubSurface   )
!            deallocate (Me%ExtVar%FertilAmmoniaSurface      )
!            deallocate (Me%ExtVar%FertilAmmoniaSubSurface   )
!            deallocate (Me%ExtVar%FertilOrganicNSurface     )
!            deallocate (Me%ExtVar%FertilOrganicNSubSurface  )
!            deallocate (Me%ExtVar%FertilOrganicPSurface     )
!            deallocate (Me%ExtVar%FertilOrganicPSubSurface  )
!            deallocate (Me%ExtVar%FertilMineralPSurface     )
!           deallocate (Me%ExtVar%FertilMineralPSubSurface  )
!            deallocate (Me%ExtVar%NitrogenUptake            )
!           deallocate (Me%ExtVar%PhosphorusUptake          )
!        endif

    
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

        call GetFluxU           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR030'

        call GetFluxV           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR040'

        call GetFluxWOld        (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR050'

        call GetUnsatWOld       (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'

        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR070'

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR080'

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
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR100'


        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR110'


        call GetComputeFaces3D(Me%ObjMap,                                               &
                               ComputeFacesW3D = Me%ExtVar%ComputeFacesW3D,             &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR120'


    end subroutine ReadLockExternalVar

    !-----------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR010'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR020'
        
        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR030'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR040'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR050'

        call UnGetPorousMedia           (Me%ObjPorousMedia,Me%ExtVar%UnsatW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'
        
        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR070'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR080'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%OpenPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModulePorousMediaProperties - ERR085'
        
        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%CellVolume,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR090'


        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%SZZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR100'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%DWZ,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR110'

        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR120'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR130'


    endsubroutine ReadUnlockExternalVar


end module ModulePorousMediaProperties

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








