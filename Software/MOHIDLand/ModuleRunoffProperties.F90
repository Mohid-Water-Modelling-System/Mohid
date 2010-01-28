!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : RunoffProperties
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as Runoff Properties 
!
!------------------------------------------------------------------------------

!
!Units in runoff properties
!   Transported properties (soluble)  : g/m3 (or mg/l)  
!   Adsorbed properties (non soluble) : ug/kgsoil       
!
Module ModuleRunoffProperties

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
    use ModuleHorizontalGrid,   only : GetHorizontalGridSize, GetHorizontalGrid,         &
                                       GetGridCellArea,                                  &
                                       WriteHorizontalGrid, UnGetHorizontalGrid
    use ModuleBasinGeometry,    only : GetBasinPoints, GetRiverPoints,  UnGetBasin 
                                       
    use ModuleFillMatrix,         only : ConstructFillMatrix, GetDefaultValue,             &
                                         KillFillMatrix, ModifyFillMatrix
    use ModuleGeometry
    use ModuleHorizontalMap,      only : GetComputeFaces2D, UngetHorizontalMap
    use ModuleRunoff,             only : GetOverLandFlow, UnGetRunoff, GetRunoffWaterColumn,  &
                                         GetFlowToChannels, GetRunoffCenterVelocity,          &
                                         GetRunoffWaterColumnOld, GetRunoffWaterColumn       
!    use ModuleInterface,          only : ConstructInterface, Modify_Interface
!    use ModuleAdvectionDiffusion, only : StartAdvectionDiffusion, AdvectionDiffusion,      &
!                                         GetAdvFlux, GetDifFlux, GetBoundaryConditionList, &
!                                         UngetAdvectionDiffusion, KillAdvectionDiffusion

   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !RunoffProperties
    public  :: ConstructRunoffProperties
    private ::      AllocateInstance
    private ::      ReadFileNames
    private ::      ReadGlobalOptions
    private ::      Construct_PropertyList
    private ::      ConstructHDF
    private ::      ConstructTimeSerie
!    private ::      CoupleSoilQuality
!#ifdef _PHREEQC_       
!    private ::      CoupleSoilChemistry
!#endif    
!    private ::      StartAdvectionDiffusion

    !Selector
    public  :: GetRPnProperties
    public  :: GetRPPropertiesIDByIdx    
    public  :: GetRPConcentration
    public  :: SetDNConcRP       !RP gets DN conc 
    public  :: SetBasinConcRP    !RP gets Basin WC conc (updated each time a module changes water column)
!    public  :: SetWindVelocity
    public  :: UnGetRunoffProperties
    
    !Modifier
    public  :: ModifyRunoffProperties
!    private ::      InterfaceFluxes                 !Gets Conc from DN
    private ::      ActualizePropertiesFromFile
    private ::      AdvectionDiffusionProcesses_RP !Explicit in runoff properties module
    private ::          ModifyAdvectionDiffusion 
!    private ::      AdvectionDiffusionProcesses_AD  !Using Module Advection Diffusion
!    private ::          ComputeVolumes  
!    private ::          AdvectionDiffusion  
!    private ::      SoilQualityProcesses
!#ifdef _PHREEQC_    
!    private ::      SoilChemistryProcesses
!#endif        
    private ::      OutPut_TimeSeries
    private ::      Output_HDF
    
    !Destructor
    public  :: KillRunoffProperties                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjRunoffProperties 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetRunoffProperties2D_I
    private :: UnGetRunoffProperties2D_R8
    private :: UnGetRunoffProperties2D_R4
    interface  UnGetRunoffProperties
        module procedure UnGetRunoffProperties2D_I
        module procedure UnGetRunoffProperties2D_R8
        module procedure UnGetRunoffProperties2D_R4
    end interface  UnGetRunoffProperties


    !Parameters----------------------------------------------------------------

    real,    parameter                :: WaterReferenceDensity = 1000. ![kg/m3]
    
    integer, parameter                :: DirectionX = 1
    integer, parameter                :: DirectionY = 2

    character(LEN = StringLength), parameter :: prop_block_begin     = '<beginproperty>'
    character(LEN = StringLength), parameter :: prop_block_end       = '<endproperty>'
    integer, parameter                :: AdvDif_ModuleRP_  = 1
    integer, parameter                :: AdvDif_ModuleAD_   = 2
    integer, parameter                :: AdvDif_Upwind_     = 1
    integer, parameter                :: AdvDif_CentralDif_ = 2    
    integer, parameter                :: AdvDif_Diff_Jury_  = 1
    integer, parameter                :: AdvDif_Diff_Old_   = 2        
    !Types---------------------------------------------------------------------
    
    private :: T_RunoffProperties

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

    type T_Property_2D
        type(T_PropertyID)               :: ID
        real, pointer, dimension (:,:)   :: Field
        real                             :: Scalar
    end type T_Property_2D


    type T_ExtVar
        !Map
        integer, pointer, dimension(:,:,:)      :: LandPoints3D
        integer, dimension(:,:), pointer        :: BasinPoints
        integer, dimension(:,:), pointer        :: RiverPoints
        real                                        :: RunoffpropDT
        type(T_Time)                                :: Now
        type(T_Time)                                :: BeginTime
        type(T_Time)                                :: EndTime
   
        ! from Runoff
        real,    dimension(:,:), pointer           :: CenterVelV
        real,    dimension(:,:), pointer           :: CenterVelU
        real(8), dimension(:,:), pointer           :: FlowToChannels
        real(8), pointer, dimension(:,:)            :: WaterColumn
        real(8), pointer, dimension(:,:)            :: WaterColumnOld
        real(8), pointer, dimension(:,:)            :: CellVolume
        real(8), pointer, dimension(:,:)            :: CellWaterMass
        real(8), dimension(:,:), pointer            :: FluxU
        real(8), dimension(:,:), pointer            :: FluxV
        real(8), pointer, dimension(:,:  )          :: Area
        real                                        :: DT
        real(8), pointer, dimension(:,:  )          :: DZY
        real(8), pointer, dimension(:,:  )          :: DZX
        real(8), pointer, dimension(:,:  )          :: DXX
        real(8), pointer, dimension(:,:  )          :: DYY        
!        real(8), pointer, dimension(:,:  )          :: DUY
!        real(8), pointer, dimension(:,:  )          :: DVY
!        real   , pointer, dimension(:,:  )          :: Topography  
        
       
        real(8), pointer, dimension(:,:)            :: InfiltrationFlux

        logical                                     :: CoupledDN  = .false.

        !from basin
        real,    dimension(:,:  ), pointer          :: WindVelocity2D  !m/s
        real,    dimension(:,:  ), pointer          :: WindVelocity  !km/day
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
!        !--For AdvectionDiffusion module use
!        integer                                :: BoundaryCondition
!        real                                   :: SchmidtNumberH
!        real                                   :: SchmidtCoefV
!        real                                   :: SchmidtBackgroundV
!        real                                   :: DiffusionH_imp_exp
!        real                                   :: ImplicitH_direction
!        logical                                :: Nulldif          = .false.
!        logical                                :: NumericStability = .false.
!        real                                   :: VolumeRelMax
!        integer                                :: AdvMethodH, TVDLimitationH
!        integer                                :: AdvMethodV, TVDLimitationV
!        logical                                :: Upwind2H, Upwind2V  
!        logical                                :: Adv_Dif_Explicit                      
        !--For both models use
        real                                   :: Molecular_Diff_Coef 
                           
    end type T_AdvectionDiffusion

    type T_Partition
        logical                                 :: NonLinear
        character(LEN = StringLength)           :: NonLinear_ks_Units
        type(T_Property_2D)                     :: Nu            
        type(T_Property_2D)                     :: Be          
        type(T_Property_2D)                     :: ks
        type(T_Property_2D)                     :: PartitionRate
        type(T_Property_2D)                     :: Fraction 
        character (LEN = StringLength)          :: Partition_Couple
    end type T_Partition

    type T_Evolution
        logical                                 :: Variable = .false.
        real                                    :: DTInterval
        type(T_Time)                            :: LastCompute
        type(T_Time)                            :: NextCompute
!        logical                                 :: SoilQuality
!        logical                                 :: SoilChemistry
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
        real(8), dimension(:,:), pointer        :: Concentration            => null()
        real(8), dimension(:,:), pointer        :: ConcentrationOld         => null()
        real(8), dimension(:,:), pointer        :: ConcentrationDN               => null()

        real(8), dimension(:,:), pointer        :: MassOnWaterColumn    => null()

        real, pointer, dimension(:,:)           :: Mass_Created
        real(8),    pointer, dimension(:,:)     :: ViscosityU
        real(8),    pointer, dimension(:,:)     :: ViscosityV
        type (T_Property), pointer              :: Next, Prev                     => null()
        logical                                 :: Particulate
        type (T_Evolution)                      :: Evolution
        real(8), pointer, dimension(:,:)        :: Diff_Turbulence_H
        real(8), pointer, dimension(:,:)        :: Viscosity
        real(8), pointer, dimension(:,:)        :: Diffusivity

        logical                                 :: Old     = .false.
        real                                    :: MinValue        = FillValueReal
        logical                                 :: TimeSerie        = .false.
        logical                                 :: BoxTimeSerie     = .false.
        logical                                 :: BoxTimeSerie2D   = .false.
        logical                                 :: OutputHDF        = .false.
!        type (T_RelatedID)                      :: RelatedID
        
    end type T_Property

    type T_Coupled
!        logical                                 :: SoilQuality          = .false. !Sediment source/sink model (Sediment Quality)
!        real                                    :: SoilQuality_DT
!        type (T_Time)                           :: SoilQuality_NextCompute

!#ifdef _PHREEQC_        
!        logical                                 :: SoilChemistry        = .false.  !Chemical reactions model (PhreeqC)
!        real                                    :: SoilChemistry_DT
!        type (T_Time)                           :: SoilChemistry_NextCompute
!#endif        
        logical                                 :: AdvectionDiffusion   = .false.
        logical                                 :: MinConcentration     = .false.
    end type T_Coupled

    type T_RunoffProperties
        integer                                     :: ObjTime                   = 0
        integer                                     :: ObjHorizontalGrid    = 0
!        integer                                     :: ObjAdvectionDiffusion     = 0
        integer                                     :: ObjBasinGeometry     = 0
        integer                                     :: ObjRunoff            = 0
!        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjGridData          = 0
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjtimeSerie         = 0
!        integer                                     :: ObjSedimentQuality   = 0
        integer                                     :: ObjHDF5              = 0
!        integer                                     :: ObjBottomTopography  = 0
        integer                                     :: ObjProfile           = 0
        integer                                     :: ObjInterface         = 0
!#ifdef _PHREEQC_        
!        integer                                     :: ObjPhreeqC                = 0
!        integer                                     :: ObjInterfaceSoilChemistry = 0 
!#endif        
        type (T_ExtVar)                             :: ExtVar
        logical                                     :: CheckGlobalMass      
        type (T_Files)                              :: Files
        type (T_OutPut)                             :: OutPut
        type (T_Property), pointer                  :: FirstProperty    => null() !Lúcia
        type (T_Property), pointer                  :: LastProperty        
        type (T_RunoffProperties), pointer          :: Next             => null() !Lúcia
        type (T_Coupled)                            :: Coupled
        type (T_Time)                               :: LastOutputHDF5

        logical                                     :: RunoffProperties
        real,    pointer, dimension(:,:)            :: Volume   
        integer                                     :: PropertiesNumber    = 0
        real   , pointer, dimension(:,:)            :: DissolvedToParticulate2D
        real                                        :: ResidualTime
        
        integer                                     :: InstanceID
        type (T_Size2D)                             :: Size, WorkSize

        type(T_Property_2D)                         :: Disper_Trans
        type(T_Property_2D)                         :: Disper_Longi
               
        integer                                     :: AdvDiff_Module        ! 1 - RunoffProperties, 2 - AdvectionDiffusion
        integer                                     :: AdvDiff_SpatialMethod !
        logical                                     :: AdvDiff_Explicit      !
        logical                                     :: AdvDiff_CheckCoefs    !
        integer                                     :: AdvDiff_DiffMethod    !  1 - Jury based, 2 - AdvectionDiffusion    
        !--For RunoffProperties Advection-Diffusion Method
!        real,    pointer, dimension(:,:)          :: DifusionNumber
!        real,    pointer, dimension(:,:)          :: ReynoldsMNumber    
        
      
!        real(8), pointer, dimension(:,:)          :: WaterVolume
!        real(8), pointer, dimension(:,:)          :: WaterVolumeOld                
!        real(8), pointer, dimension(:,:)          :: WaterVolumeCorr

    end type  T_RunoffProperties

    !Global Module Variables
    type (T_RunoffProperties), pointer                         :: FirstObjRunoffProperties
    type (T_RunoffProperties), pointer                         :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTORCONSTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructRunoffProperties(ObjRunoffPropertiesID,                           &
                                              ComputeTimeID,                              &
                                              HorizontalGridID,                           &
                                              HorizontalMapID,                            &
                                              BasinGeometryID,                            &
                                              RunoffID,                                   &
!                                              GeometryID,                                 &
                                              CoupledDN,                                  &
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: BasinGeometryID
        integer                                         :: RunoffID
!        integer                                         :: GeometryID
        logical, optional                               :: CoupledDN 
        integer, optional, intent(OUT)                  :: STAT 
        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_,STAT_CALL
        !------------------------------------------------------------------------
                                    

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mRunoffProperties_)) then
            nullify (FirstObjRunoffProperties)
            call RegisterModule (mRunoffProperties_) 
        endif

        call Ready(ObjRunoffPropertiesID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            

            call AllocateInstance

            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ComputeTimeID   )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjRunoff         = AssociateInstance (mRUNOFF_,              RunoffID   )
!            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
        

            if (present(CoupledDN)) then
                Me%ExtVar%CoupledDN  = CoupledDN
            endif            
            
            

            call ReadFileNames


            !Constructs the DataFile
            call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR01'                
           
            call ReadGlobalOptions

            call AllocateVariables

            call Construct_PropertyList
            
!            if (Me%Coupled%SoilQuality) then
!                call Construct_InitialFields
!            endif
        
            call ConstructHDF    
    
            call ConstructTimeSerie
            
            if (Me%Coupled%AdvectionDiffusion .and. Me%AdvDiff_CheckCoefs) then
                call ConstructAsciiFile
            endif

            
            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR010'



            !Couple nutrient, carbon and oxygen sources and sinks model
!            if (Me%Coupled%SoilQuality) then
!                call CoupleSoilQuality
!            endif
            
!#ifdef _PHREEQC_            
!            !Couple soil chemical model
!            if (Me%Coupled%SoilChemistry) then
!                call CoupleSoilChemistry
!            endif
!#endif
            
!            if (Me%Coupled%AdvectionDiffusion .and. Me%AdvDiff_Module == AdvDif_ModuleAD_) then !Uses AdvectionDiffusion module
!            
!                call StartAdvectionDiffusion(Me%ObjAdvectionDiffusion, &
!                                             Me%ObjGeometry,           &
!                                             Me%ObjHorizontalMap,      &
!                                             Me%ObjHorizontalGrid,     &
!                                             Me%ObjTime,               &
!                                             STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR50'
!
!            endif

            !Returns ID
            ObjRunoffPropertiesID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructRunoffProperties - ModuleRunoffProperties - ERR040' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructRunoffProperties
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_RunoffProperties), pointer                         :: NewObjRunoffProperties
        type (T_RunoffProperties), pointer                         :: PreviousObjRunoffProp


        !Allocates new instance
        allocate (NewObjRunoffProperties)
        nullify  (NewObjRunoffProperties%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjRunoffProperties)) then
            FirstObjRunoffProperties         => NewObjRunoffProperties
            Me                    => NewObjRunoffProperties
        else
            PreviousObjRunoffProp      => FirstObjRunoffProperties
            Me                    => FirstObjRunoffProperties%Next
            do while (associated(Me))
                PreviousObjRunoffProp  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjRunoffProperties
            PreviousObjRunoffProp%Next => NewObjRunoffProperties
        endif

        Me%InstanceID = RegisterNewInstance (mRunoffProperties_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    subroutine ReadFileNames

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
!        integer                                     :: iflag

        !Reads the name of the data file from nomfich
        call ReadFileName ('RUNOFF_PROP_DATA', Me%Files%DataFile, "Runoff Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoffProperties - ERR01'

        !Reads the name of the transient HDF file from nomfich
        call ReadFileName ('RUNOFF_PROP_HDF', Me%Files%TransientHDF, "Runoff HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoffProperties - ERR01b'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('RUNOFF_PROP_FIN', Me%Files%FinalFile, "Runoff Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModuleRunoffProperties - ERR01c'
   

    end subroutine ReadFileNames
    
    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        type(T_Property_2D), pointer                :: Scalar2D
        !Begin-----------------------------------------------------------------

        !Geometry Size
        call GetHorizontalGridSize (Me%ObjHorizontalGrid,                            &
                                    Size     = Me%Size,                              &
                                    WorkSize = Me%WorkSize,                          &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR01'


        call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR20'

        call GetComputeTimeLimits(Me%ObjTime,                      &
                                  EndTime   = Me%ExtVar%EndTime,   &
                                  BeginTime = Me%ExtVar%BeginTime, &
                                  STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)    &
                stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR30'

        ! Sets the last output equal to zero 
        call SetDate(Me%LastOutPutHDF5, 0, 0, 0, 0, 0, 0)

        !Needs keyword description and value options definition
        call GetData(Me%AdvDiff_Module,                            &   
                     Me%ObjEnterData, iflag,                       &
                     SearchType   = FromFile,                      &
                     keyword      = 'ADVDIFF_MODULE',              &
                     Default      = AdvDif_ModuleRP_,             &
                     ClientModule = 'ModuleRunoffProperties', &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR40'


		if (Me%AdvDiff_Module == AdvDif_ModuleRP_) then !Uses ModuleRunoffProperties for Advection-Diffusion calculation

	        call GetData(Me%AdvDiff_SpatialMethod,                     &   !Lúcia
	                     Me%ObjEnterData, iflag,                       &
	                     SearchType   = FromFile,                      &
	                     keyword      = 'ADVDIFF_SPATIAL_METHOD',      &
	                     Default      = AdvDif_CentralDif_,            &
	                     ClientModule = 'ModuleRunoffProperties', &
	                     STAT         = STAT_CALL)
	        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR50'

            call GetData(Me%AdvDiff_CheckCoefs,                        &   !Eduardo
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromFile,                      &
                         keyword      = 'ADVDIFF_CHECK_COEFS',         &
                         Default      = .false.,                       &
                         ClientModule = 'ModuleRunoffProperties', &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR80'

            call GetData(Me%AdvDiff_DiffMethod,                                               &   
                         Me%ObjEnterData, iflag,                                              &
                         SearchType = FromFile,                                               &
                         keyword    = 'ADVDIFF_DIFF_METHOD',                                  &
                         Default    = AdvDif_Diff_Jury_,                                      &                                           
                         ClientModule ='ModuleRunoffProperties',                              &
                         STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleRunoffProeprties - ERR085'            

		elseif (Me%AdvDiff_Module == AdvDif_ModuleAD_) then  ! Uses ModuleAdvectionDiffusion for Advection-Diffusion calculation

!	        call GetData(Me%AdvDiff_Explicit,                          &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromFile,                      &
!	                     keyword      = 'ADVDIFF_EXPLICIT',            &
!	                     Default      = .true.,                	       &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_)  &
	            stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR90'
        
        else
	        write(*,*)'Advection diffusion module to be used unrecognized,'
	        write(*,*)'Please check ADVDIFF_MODULE keyword'
	        stop 'ReadGlobalOptions - ModuleRunoffProperties - ERR100'        
        
		endif

!        Scalar2D => Me%Disper_Longi
!        call ConstructScalar2D(Scalar2D, ExtractType = FromBlock,               &
!                              block_begin = '<begin_dispersion_long>',          &
!                              block_end   = '<end_dispersion_long>')
!        
        Scalar2D => Me%Disper_Trans
        call ConstructScalar2D(Scalar2D, ExtractType = FromBlock,               &
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

               
        !Water Content---------------------------------------------------------
!       allocate (Me%DifusionNumber          (ILB:IUB,JLB:JUB))
!        allocate (Me%ReynoldsMNumber         (ILB:IUB,JLB:JUB))
!        allocate (Me%ExtVar%WindVelocity   (ILB:IUB,JLB:JUB))
        
!#ifdef _PHREEQC_
!        allocate (Me%ExtVar%CellWaterMass    (ILB:IUB,JLB:JUB,))
!#endif        

 !       allocate (Me%WaterVolume          (ILB:IUB,JLB:JUB))
!        allocate (Me%WaterVolumeOld       (ILB:IUB,JLB:JUB))
!        allocate (Me%WaterVolumeCorr          (ILB:IUB,JLB:JUB))
    !    allocate (Me%WaterVolumeOldCorr       (ILB:IUB,JLB:JUB,KLB:KUB))
        
       
 !       Me%WaterVolume          = 0.
 !       Me%WaterVolumeOld       = 0.
!        Me%WaterVolumeCorr      = 0.
        

    end subroutine AllocateVariables


    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
  
    subroutine ConstructScalar2D(Scalar2D, ExtractType, ClientNumber, block_begin, block_end)

        !Arguments-------------------------------------------------------------
        type(T_Property_2D), pointer        :: Scalar2D
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

        select case(ExtractType)

            case(FromBlock)

                call ExtractBlockFromBuffer(Me%ObjEnterData, BlockClientNumber, block_begin, block_end,  &
                                            BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleRunoffProperties - ERR01'


            case(FromBlockInBlock)

                if(.not. present(ClientNumber))then
                    stop 'ConstructScalar2D - ModuleSoilProperties - ERR02'
                end if
                
                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber, block_begin, block_end,  &
                                           BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR20'

        end select

        if(BlockFound)then

            allocate(Scalar2D%Field(ILB:IUB, JLB:JUB))

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR030'

            call ConstructFillMatrix  (PropertyID           = Scalar2D%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill2D       = Me%ExtVar%BasinPoints,            &
                                       Matrix2D             = Scalar2D%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR040'


            call GetDefaultValue(Scalar2D%ID%ObjFillMatrix, Scalar2D%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR050'

            call UnGetBasin(Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR60'

            call KillFillMatrix(Scalar2D%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR065'

            if(ExtractType == FromBlockInBlock)then
                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR070'
            end if


            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, BlockClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR080'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar2D - ModuleRunoffProperties - ERR090'
            end if
        
        else
            write(*,*) 'Block not present:', block_begin, block_end
            stop 'ConstructScalar2D - ModuleRunoffProperties - ERR100'
        end if

   
    end subroutine ConstructScalar2D

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
                        stop 'Construct_PropertyList - ModuleRunoffProperties - ERR01'
                    exit do1    !No more blocks
                
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_PropertyList - ModuleRunoffProperties - ERR02'
            
            else cd1
                
                stop 'Construct_PropertyList - ModuleRunoffProperties - ERR03'
            
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
        if(STAT_CALL .NE. SUCCESS_)stop 'Construct_Property - ModuleRunoffProperties - ERR00'
        
        nullify(NewProperty%Prev, NewProperty%Next)
        nullify(NewProperty%Concentration)
        nullify(NewProperty%MassOnWaterColumn)

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
                     ClientModule = 'ModuleRunoffProperties',                       &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModuleRunoffProperties - ERR01'
        if(iflag == 0)              stop 'Construct_PropertyState - ModuleRunoffProperties - ERR02'

        !Not used the function so this text was commented
!        if (NewProperty%Particulate)then
!            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
!                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
!                write(*,*) 'recognised as PARTICULATE'
!                stop 'Construct_PropertyState - ModuleRunoffProeprties - ERR03'
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
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR10'

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

!        call GetData(NewProperty%Evolution%SoilQuality,                                  &
!                     Me%ObjEnterData,iflag,                                              &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'SOIL_QUALITY',                                      &
!                     ClientModule = 'ModuleRunoffProperties',                            &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR20'
!
!        if (NewProperty%Evolution%SoilQuality) then
!            Me%Coupled%SoilQuality     = .true.
!            NewProperty%Evolution%Variable = .true.
!        endif


!#ifdef _PHREEQC_
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

!        call GetData(NewProperty%Evolution%SoilChemistry,                                &
!                     Me%ObjEnterData,iflag,                                              &
!                     SearchType   = FromBlock,                                           &
!                     keyword      = 'SOIL_CHEMISTRY',                                    &
!                     ClientModule = 'ModuleRunoffProperties',                            &
!                     default      = OFF,                                                 &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_)                                                     &
!            stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR30'
!
!        if (NewProperty%Evolution%SoilChemistry) then
!            Me%Coupled%SoilChemistry       = .true.
!            NewProperty%Evolution%Variable = .true.
!        endif
!#endif

        !Property time step
        if (NewProperty%Evolution%Variable) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR40'

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModuleRunoffProperties',                        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR050'
                                       
            
            if (NewProperty%Evolution%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR60'

            elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR70'
                endif

                !Run period in seconds
                DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%Evolution%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'Construct_PropertyEvolution - ModuleRunoffProperties - ERR80'
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
                     ClientModule = 'ModuleRunoffProperties',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR10'

cd1:    if (Me%AdvDiff_Module == AdvDif_ModuleAD_) then

!	        call GetData(NewProperty%Evolution%AdvDiff%NumericStability, &
!	                     Me%ObjEnterData, iflag,                         &
!	                     SearchType   = FromBlock,                       &
!	                     keyword      = 'ADVDIFF_NUM_STABILITY',         &
!	                     Default      = .FALSE.,                         &
!	                     ClientModule = 'ModuleRunoffProperties',   &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR20'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%SchmidtNumberH, &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                     &
!	                     keyword      = 'ADVDIFF_SCHMIDT_NUMBER_H',    &
!	                     Default      = 1.0,                           &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR30'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%SchmidtCoefV,   &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                     &
!	                     keyword      = 'ADVDIFF_SCHMIDT_COEF_V',      &
!	                     Default      = 1.0,                           &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR40'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%SchmidtBackgroundV, &
!	                     Me%ObjEnterData, iflag,                           &
!	                     SearchType   = FromBlock,                         &
!	                     keyword      = 'ADVDIFF_SCHMIDT_BACKGROUND_V',    &
!	                     Default      = 0.,                                &
!	                     ClientModule = 'ModuleRunoffProperties',     &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR50'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%NullDif,        &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                     &
!	                     keyword      = 'ADVDIFF_NULLDIF',             &
!	                     Default      = .false.,                       &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR60'
!
!	        call GetBoundaryConditionList(MassConservation    = MassConservation,    &
!	                                      ImposedValue        = ImposedValue,        &
!	                                      NullGradient        = NullGradient,        &
!	                                      Orlanski            = Orlanski,            &
!	                                      MassConservNullGrad = MassConservNullGrad, &
!	                                      CyclicBoundary      = CyclicBoundary)
!
!	        call GetData(BoundaryCondition,                            &
!	                     Me%ObjEnterData,  iflag,                      &
!	                     SearchType   = FromBlock,                     &
!	                     keyword      = 'ADVDIFF_BOUNDARY_CONDITION',  &
!	                     Default      = MassConservation,              &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!	        if (STAT_CALL .NE. SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR70'
!
!	        ! By default it's imposed a value dependent only from the exterior
!	        ! value and of the decay time. However this method doesn't conserve mass
!	        ! when the water fluxes near the frontier are dominant
!
!	        if (BoundaryCondition /= MassConservation     .and. &
!	            BoundaryCondition /= ImposedValue         .and. &
!	            BoundaryCondition /= NullGradient         .and. &
!	            BoundaryCondition /= CyclicBoundary       .and. &
!	            BoundaryCondition /= Orlanski             .and. &
!	            BoundaryCondition /= MassConservNullGrad) &
!	            stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR80'
!
!	        NewProperty%Evolution%AdvDiff%BoundaryCondition = BoundaryCondition
!
!	        !By default the horizontal Diffusion discretization is explicit
!	        NewProperty%Evolution%AdvDiff%DiffusionH_imp_exp  = ExplicitScheme
!
!	        NewProperty%Evolution%AdvDiff%ImplicitH_Direction = DirectionX
!
!	        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodH,     &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                      &
!	                     keyword      = 'ADVDIFF_METHOD_H',            &
!	                     Default      = UpwindOrder1,                  &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!
!	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR90'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationH, &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                      &
!	                     keyword      = 'ADVDIFF_TVD_LIMIT_H',         &
!	                     Default      = Superbee,                      &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!
!	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR100'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%AdvMethodV,     &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                      &
!	                     keyword      = 'ADVDIFF_METHOD_V',            &
!	                     Default      = UpwindOrder1,                  &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!
!	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR110'
!
!	        call GetData(NewProperty%Evolution%AdvDiff%TVDLimitationV, &
!	                     Me%ObjEnterData, iflag,                       &
!	                     SearchType   = FromBlock,                      &
!	                     keyword      = 'ADVDIFF_TVD_LIMIT_V',         &
!	                     Default      = Superbee,                      &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!
!	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR120'
!
!
!	        call GetData(NewProperty%Evolution%AdvDiff%VolumeRelMax,   &
!	                     Me%ObjEnterData, iflag,                       &
!	                     Keyword      = 'ADVDIFF_VOLUME_RELATION_MAX', &
!	                     Default      = 5.,                            &
!	                     SearchType   = FromBlock,                      &
!	                     ClientModule = 'ModuleRunoffProperties', &
!	                     STAT         = STAT_CALL)
!
!	        if (STAT_CALL /= SUCCESS_) stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR130'
!
!
!	        if (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
!	            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3 .or.&
!	            NewProperty%Evolution%AdvDiff%AdvMethodH == P2_TVD) then
!	            NewProperty%Evolution%AdvDiff%Upwind2H = .true.
!	        else
!	            NewProperty%Evolution%AdvDiff%Upwind2H = .false.
!	        endif
!
!	        if (NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder2 .or.&
!	            NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder3 .or.&
!	            NewProperty%Evolution%AdvDiff%AdvMethodV == P2_TVD) then
!	            NewProperty%Evolution%AdvDiff%Upwind2V = .true.
!	        else
!	            NewProperty%Evolution%AdvDiff%Upwind2V = .false.
!	        endif
!
!
!	        if (.not. Me%AdvDiff_Explicit .and.&
!	           (NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder2 .or.&
!	            NewProperty%Evolution%AdvDiff%AdvMethodH == UpwindOrder3)) then
!
!	            write(*,*) 'If the advection of mass in the horizontal is implicit'
!	            write(*,*) 'the advection method can not be a second or third order upwind'
!	            stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR140.'
!
!	        endif
!
!	        if (.not. Me%AdvDiff_Explicit .and.&
!	           (NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder2 .or.&
!	            NewProperty%Evolution%AdvDiff%AdvMethodV == UpwindOrder3)) then
!
!	            write(*,*) 'If the advection of mass in the vertical is implicit'
!	            write(*,*) 'the advection method can not be a second or third order upwind'
!	            stop 'ReadAdvectionDiffusionParameters - ModuleRunoffProperties - ERR150.'
!
!	        endif

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
        allocate (NewProperty%Diffusivity (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR10'
        NewProperty%Diffusivity       = 0. 
        
        if (Me%AdvDiff_Module == AdvDif_ModuleAD_) then !Module advection-Diffusion is used
           
            allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR20'
            NewProperty%Viscosity         = 0.
            
            allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR30'
            NewProperty%Diff_Turbulence_H = 0.

!            allocate (NewProperty%Diff_Turbulence_V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR40'
!            NewProperty%Diff_Turbulence_V = 0.                   
            
       elseif (Me%AdvDiff_Module == AdvDif_ModuleRP_) then !explicit calculation inside module porous media properties
            allocate (NewProperty%ViscosityU (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR60'
            
            NewProperty%ViscosityU  = 0.0
            
            allocate (NewProperty%ViscosityV (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR70'
        
            NewProperty%ViscosityV  = 0.0!     
            
            if (Me%AdvDiff_DiffMethod == AdvDif_Diff_Old_) then 

                allocate (NewProperty%Viscosity (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR20'
                NewProperty%Viscosity         = 0.
                
                allocate (NewProperty%Diff_Turbulence_H (ILB:IUB, JLB:JUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR30'
                NewProperty%Diff_Turbulence_H = 0.

!                allocate (NewProperty%Diff_Turbulence_V (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyDiffusivity - ModuleRunoffProperties - ERR40'
!                NewProperty%Diff_Turbulence_V = 0.                                   
           
               
            endif
            
        endif


    end subroutine ConstructPropertyDiffusivity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine Construct_PropertyValues(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property),              pointer      :: NewProperty

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j

        !Local-----------------------------------------------------------------
        integer                                     :: iflag
        integer                                     :: ILB,IUB
        integer                                     :: JLB,JUB
        integer                                     :: WorkSizeILB, WorkSizeIUB
        integer                                     :: WorkSizeJLB, WorkSizeJUB
        
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB

        allocate(NewProperty%Concentration(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR10'
        NewProperty%Concentration(:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR20'
        NewProperty%ConcentrationOld(:,:) = FillValueReal

      
        if (Me%ExtVar%CoupledDN) then
            allocate(NewProperty%ConcentrationDN(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModulePorousMediaProperties - ERR40'
            NewProperty%ConcentrationDN(:,:) = FillValueReal
        endif


        allocate(NewProperty%MassOnWaterColumn(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR30'
        NewProperty%MassOnWaterColumn(:,:) = 0.
      
!        !Eduardo Jauch 12nov2009
!        call GetData(NewProperty%RelatedID%name,                                            &
!                     Me%ObjEnterData, iflag,                                                &
!                     SearchType   = FromBlock,                                              &
!                     keyword      = 'RELATED',                                              &
!                     Default      = '',                                                     &                        
!                     ClientModule = 'ModuleRunoffProperties',                                 &
!                     STAT         = STAT_CALL)              
!        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR50'
!        if ((iflag .NE. 0) .AND. (NewProperty%RelatedID%name .NE. '')) then
!            NewProperty%RelatedID%IDNumber = GetPropertyIDNumber(NewProperty%RelatedID%name)    
!        endif
!        !ToDo: After load all the properties, must check if "related properties" exist...        
!        !end


        call GetData(NewProperty%MinValue,                                                  &
                     Me%ObjEnterData,iflag,                                                 &
                     SearchType   = FromBlock,                                              &
                     keyword      = 'MIN_VALUE',                                            &
                     ClientModule = 'ModuleRunoffProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR100'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
            Me%Coupled%MinConcentration = .true.
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR110'
            NewProperty%Mass_Created(:,:) = 0.
        endif

        !This variable is a logic one is true if the property is old
        !and the user wants to continue the run with results of a previous run.
        call GetData(NewProperty%Old,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'OLD',                                              &
                     Default      = .false.,                                            &                        
                     SearchType   = FromBlock,                                          &
                     ClientModule = 'ModuleRunoffProperties',                                 &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR120'
          
        ! if the property is not 'OLD' the property values in the domain and 
        ! in the boundaries are initialized
        ! if it's true ('OLD') this same values are read from the final file of the
        ! previous run
        if (.not. NewProperty%Old) then

            !Get water points
            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR130'

            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,                   &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExtVar%BasinPoints,            &
                                       Matrix2D             = NewProperty%Concentration,        &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR140'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR0150'
            end if

            call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR151'            
            
            !initial concentration based on initial water column
	        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
	        do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then 
                    if (Me%ExtVar%WaterColumn(i,j) .gt. 0.0) then
                        NewProperty%ConcentrationOld(i,j) = NewProperty%Concentration(i,j)           
                    else
                        NewProperty%Concentration(i,j) = 0.0
                        NewProperty%ConcentrationOld(i,j) = NewProperty%Concentration(i,j)
                   endif
                endif
            enddo
            enddo
            
!            call SetMatrixValue(NewProperty%ConcentrationOld, Me%Size, NewProperty%Concentration,Me%ExtVar%BasinPoints)

            call UnGetRunoff     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Construct_PropertyValues - ModuleRunoffProperties - ERR152'            


            call CheckFieldConsistence (NewProperty)

            call UnGetBasin(Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPropertyValues - ModuleRunoffProperties - ERR160'

        else

            ! If the property is old then the program is going to try to find a property
            ! with the same name in the Water properties initial file written in HDF format  
            call ReadOldConcBoundariesHDF(NewProperty)

        end if   

    end subroutine Construct_PropertyValues

      !--------------------------------------------------------------------------


    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),    pointer        :: NewProperty

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------


        call GetData(NewProperty%TimeSerie,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'TIME_SERIE',                                        &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR01'
        

        call GetData(NewProperty%BoxTimeSerie,                                           &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'BOX_TIME_SERIE',                                    &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR02'


        call GetData(NewProperty%BoxTimeSerie2D,                                           &
                     Me%ObjEnterData, iflag,                                               &
                     Keyword      = 'BOX_TIME_SERIE2D',                                    &
                     Default      = .false.,                                               &
                     SearchType   = FromBlock,                                             &
                     ClientModule = 'ModuleRunoffProperties',                              &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR03'



        call GetData(NewProperty%OutputHDF,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'OUTPUT_HDF',                                        &
                     ClientModule = 'ModuleRunoffProperties',                            &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleRunoffProperties - ERR04'
        
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
                     ClientModule = 'ModuleRunoffProperties',                           &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerie_ON = .true.
        else
            Me%OutPut%TimeSerie_ON = .false.
        endif
        
        !Get water points
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR020'


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srr",                                        &
                            WaterPoints2D = Me%ExtVar%BasinPoints,                      &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - ModuleRunoffProperties - ERR030' 

        !Unget
        call UnGetBasin                   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR040'

        !Deallocates PropertyList
        deallocate(PropertyList)
       
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
                               CurrentTime = Me%ExtVar%BeginTime,            &
                               EndTime     = Me%ExtVar%EndTime,              &
                               keyword     = 'OUTPUT_TIME',                  &
                               SearchType  = FromFile,                       &
                               OutPutsTime = Me%OutPut%OutTime,              &
                               OutPutsOn   = Me%OutPut%HDF_ON,               &
                               OutPutsNumber = Me%OutPut%Number,             &
                               STAT        = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                       &
                stop 'ConstructHDF - ModuleRunoffProperties - ERR01' 

            if (Me%OutPut%HDF_ON) then

                Me%OutPut%NextOutPut = 1

                call Open_HDF5_OutPut_File

            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructHDF - ModuleRunoffProperties - ERR02'
            endif 

        endif

    end subroutine ConstructHDF

  
    !--------------------------------------------------------------------------

     subroutine Open_HDF5_OutPut_File        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        !Begin-----------------------------------------------------------------

        !Bounds
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5, trim(Me%Files%TransientHDF)//"5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunoffProperties - ERR02'

      
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunoffProperties - ERR020'


        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunoffProperties - ERR030'
        
!        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR040'
        
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR050'  

       
        !Writes the Grid
!        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",           &
!                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunoffProperties - ERR060'

        !WriteBasinPoints
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "BasinPoints", "-",          &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunoffProperties - ERR070'



        !Flushes All pending HDF5 commands
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModuleRunoffProperties - ERR080'       




        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR90'  

!        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
!        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleRunoffProperties - ERR100'


    end subroutine Open_HDF5_OutPut_File   
   
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
                 FILE   = '..\res\RP_ADCoefs_'//trim(adjustl(Number))//'.log', &
                 STATUS = "REPLACE",                                      &
                 IOSTAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                exit do1
            else
                Counter = Counter + 1
            end if
        enddo do1
        
        write (Me%Files%AsciiUnit, FMT=*) 'YY     MM   DD   HH   MM     SS     i   j    cofA_U cofB cofC_U cofA_V cofC_V'
    
    end subroutine ConstructAsciiFile
    
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
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR01'


            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Concentration(:,:) = FillValueReal


            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB,                                       &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Concentration/"//NewProperty%ID%Name,        &
                                 NewProperty%ID%Name,                                    &
                                 Array2D = NewProperty%Concentration,                    &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR03'


            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR06'

        else
            
            write(*,*)
            stop 'ReadOldConcBoundariesHDF - ModuleRunoffProperties - ERR07'

        end if cd0

    end subroutine ReadOldConcBoundariesHDF


    !--------------------------------------------------------------------------


    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_property), pointer               :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: Counter
        integer                                 :: i,j
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


            
        !Verification if the values read are lower than zero in water points
        do I = ILB, IUB
        do J = JLB, JUB
            
            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                               
                if (NewProperty%Concentration(i, j) < 0.) then
                    
                    StopSubroutine = .true.
                    NewProperty%Concentration(i, j) = 0.0

                endif

            else

                NewProperty%Concentration(i, j) = FillValueReal

            endif

        enddo
        enddo

        if (StopSubroutine) then                                                   
            
            call UnitsManager(UnitAux, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR02' 

            open(UnitAux, FILE = trim(NewProperty%ID%name)//'.new',                 &
                 FORM = 'FORMATTED', STATUS = 'UNKNOWN', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR03' 

            write(UnitAux,*) '<ConcentrationBegin>'
           
            do I = ILB, IUB
            do J = JLB, JUB
            
                write(UnitAux,*) NewProperty%Concentration(i, j)

            enddo
            enddo

            write(UnitAux,*) '<ConcentrationEnd>'

            call UnitsManager(UnitAux, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR04' 

            write(*,*) 'A new concentration file was created for property: ', trim(NewProperty%ID%Name)
            write(*,*) 'Run again with this new file ', trim(NewProperty%ID%name)//'.new'
            stop 'CheckFieldConsistence - ModuleRunoffProperties - ERR05'  

        endif

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------


!    subroutine CoupleSoilQuality        
!
!        !Local-----------------------------------------------------------------
!        type(T_Property), pointer                           :: PropertyX
!        integer, pointer, dimension(:)                      :: SoilQualityPropertyList
!        integer                                             :: STAT_CALL
!        real                                                :: SoilQualityDT
!        integer                                             :: nProp = 0 
!
!        !Begin------------------------------------------------------------------
!
!        !Counts the number of Properties which has WaterQuality option set to true
!        PropertyX => Me%FirstProperty
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilQuality) then
!                nProp = nProp + 1
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        !Allocates Array to hold IDs
!        allocate (SoilQualityPropertyList(1:nProp))
!
!        !Fills Array
!        PropertyX => Me%FirstProperty
!        nProp = 0
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilQuality) then
!                nProp = nProp + 1
!                SoilQualityPropertyList(nProp) = PropertyX%ID%IDNumber
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModuleRunoffProperties - ERR01'
!
!        !Start Interface
!        call ConstructInterface(InterfaceID         = Me%ObjInterface,               &
!                                TimeID              = Me%ObjTime,                    &
!                                SinksSourcesModel   = SedimentQualityModel,          &
!                                DT                  = SoilQualityDT,                 &
!                                PropertiesList      = SoilQualityPropertyList,       &
!                                WaterPoints2D       = Me%ExtVar%BasinPoints,         &
!                                Size2D              = Me%WorkSize,                   &
!                                STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                   &
!            stop 'CoupleSoilQuality - ModuleRunoffProperties - ERR02'
!
!
!        call UnGetBasin                   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilQuality - ModuleRunoffProperties - ERR03'
!
!
!        deallocate (SoilQualityPropertyList)
!
!        Me%Coupled%SoilQuality_DT          = SoilQualityDT 
!        Me%Coupled%SoilQuality_NextCompute = Me%ExtVar%Now    
!
!        nullify (Me%DissolvedToParticulate2D)
!        allocate(Me%DissolvedToParticulate2D(Me%WorkSize%ILB:Me%WorkSize%IUB, Me%WorkSize%JLB:Me%WorkSize%JUB))
!        Me%DissolvedToParticulate2D(:,:) = null_real
!
!        Me%ResidualTime = 0.
!    
!    end subroutine CoupleSoilQuality
!
!    !--------------------------------------------------------------------------


!#ifdef _PHREEQC_
!    !--------------------------------------------------------------------------
!    subroutine CoupleSoilChemistry        
!
!        !Local-----------------------------------------------------------------
!        type(T_Property), pointer                           :: PropertyX
!        integer, pointer, dimension(:)                      :: SoilChemistryPropertyList
!        integer                                             :: STAT_CALL
!        real                                                :: SoilChemistryDT
!        integer                                             :: nProp = 0 
!
!        !Begin------------------------------------------------------------------
!
!        !Counts the number of Properties which has SoilChemistry option set to true
!        PropertyX => Me%FirstProperty
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilChemistry) then
!                nProp = nProp + 1
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        !Allocates Array to hold IDs
!        allocate (SoilChemistryPropertyList(1:nProp))
!
!        !Fills Array
!        PropertyX => Me%FirstProperty
!        nProp = 0
!        do while (associated(PropertyX))
!            if (PropertyX%Evolution%SoilChemistry) then
!                nProp = nProp + 1
!                SoilChemistryPropertyList(nProp) = PropertyX%ID%IDNumber
!            endif
!            PropertyX => PropertyX%Next
!        enddo
!
!        !Question: What does this function? Is it necessary to SoilChemistry process or Interface?
!        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModuleRunoffProperties - ERR01'
!
!        !Start Interface
!        call ConstructInterface(InterfaceID         = Me%ObjInterfaceSoilChemistry,  &
!                                TimeID              = Me%ObjTime,                    &
!                                SinksSourcesModel   = PhreeqCModel,                  &
!                                DT                  = SoilChemistryDT,               &
!                                PropertiesList      = SoilChemistryPropertyList,     &
!                                WaterPoints2D       = Me%ExtVar%BasinPoints,         &
!                                Size3D              = Me%WorkSize,                   &
!                                STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)                                                   &
!            stop 'CoupleSoilChemistry - ModuleRunoffProperties - ERR02'
!
!        !Question: What does this function? 
!        call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'CoupleSoilChemistry - ModuleRunoffProperties - ERR03'
!
!        deallocate (SoilChemistryPropertyList)
!
!        Me%Coupled%SoilChemistry_DT          = SoilChemistryDT 
!        Me%Coupled%SoilChemistry_NextCompute = Me%ExtVar%Now    
!
!            
!    end subroutine CoupleSoilChemistry
!    !--------------------------------------------------------------------------
!#endif


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

!    subroutine GetRPCoupled(RunoffPropertiesID, &
!                             SoilQuality,             &
!#ifdef _PHREEQC_
!                             SoilChemistry,           &
!#endif                             
!                             STAT) 
!
!        !Arguments-------------------------------------------------------------
!        integer                        :: RunoffPropertiesID
!        integer, optional, intent(OUT) :: STAT
!        logical, optional, intent(OUT) :: SoilQuality        
!#ifdef _PHREEQC_
!        logical, optional, intent(OUT) :: SoilChemistry
!#endif
!
!        !External--------------------------------------------------------------
!
!        integer :: ready_              
!
!        !Local-----------------------------------------------------------------
!
!        integer :: STAT_              !Auxiliar local variable
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(RunoffPropertiesID, ready_)
!        
!cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.     &
!            (ready_ .EQ. READ_LOCK_ERR_)) then
!
!            if (present(SoilQuality   )) SoilQuality    = Me%Coupled%SoilQuality
!            
!#ifdef _PHREEQC_            
!            if (present(SoilChemistry )) SoilChemistry  = Me%Coupled%SoilChemistry
!#endif
!
!            STAT_ = SUCCESS_
!        else 
!            STAT_ = ready_
!        end if cd1
!
!        if (present(STAT)) STAT = STAT_
!
!        !----------------------------------------------------------------------
!
!    end subroutine GetRPCoupled
!    !--------------------------------------------------------------------------
    
    subroutine GetRPConcentration(RunoffPropertiesID, ConcentrationX, PropertyXIDNumber, &
                                PropertyXUnits, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: RunoffPropertiesID
        real, pointer, dimension(:,:)             :: ConcentrationX
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

        call Ready(RunoffPropertiesID, ready_) 
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                   &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            call Read_Lock(mRunoffPROPERTIES_, Me%InstanceID) 

            nullify(PropertyX)

            call Search_Property(PropertyX, PropertyXID = PropertyXIDNumber, STAT = STAT_CALL)
            if (STAT_CALL == SUCCESS_) then
                
                ConcentrationX => PropertyX%Concentration

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
            
    end subroutine GetRPConcentration

    !--------------------------------------------------------------------------------

!    subroutine SetWindVelocity (RunoffPropertiesID, WindModulus, STAT)
!                                  
!        !Arguments--------------------------------------------------------------
!        integer                                     :: RunoffPropertiesID
!        real, dimension(:,:), pointer               :: WindModulus
!
!        integer, optional, intent(OUT)              :: STAT
!
!        !Local-----------------------------------------------------------------
!        integer                                     :: ready_        
!        integer                                     :: STAT_
!        
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(RunoffPropertiesID, ready_)
!
!        if (ready_ .EQ. IDLE_ERR_)then
!            
!
!            Me%ExtVar%WindVelocity2D   => WindModulus
!
!
!            STAT_ = SUCCESS_
!        else
!            STAT_ = ready_
!        end if
!
!        if (present(STAT))STAT = STAT_
!
!    end subroutine SetWindVelocity 

    !---------------------------------------------------------------------------

   
    subroutine SetDNConcRP (RunoffPropertiesID, PropertyID, DNConcentration, ChannelsID, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer                                         :: PropertyID   
        real, dimension (:), pointer                    :: DNConcentration
        integer, dimension(:, :), pointer               :: ChannelsID
        integer                                         :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_, ready_, i, j
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then
        
            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcRP - ModuleRunoffProperties - ERR01'
            
           
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
                write(*,*) 'Looking for Drainage Network Property in Runoff Properties', GetPropertyName(PropertyID)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetDNConcRP - ModuleRunoffProperties - ERR010'
            end if

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetDNConcRP - ModuleRunoffProperties - ERR040'               


        else
            STAT_ = ready_
        end if

        STAT = STAT_        
                     
    end subroutine SetDNConcRP

    !---------------------------------------------------------------------------

    subroutine SetBasinConcRP   (RunoffPropertiesID, BasinConcentration,      &
                                                PropertyXIDNumber, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        real(8), dimension(:, :), pointer               :: BasinConcentration
        integer                                         :: PropertyXIDNumber
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: j ,i
        integer                                         :: STAT_, ready_
        type(T_Property), pointer                       :: PropertyX

        !-----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_)then

            call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR01'            
           
            nullify(PropertyX)

            call SearchProperty(PropertyX, PropertyXIDNumber = PropertyXIDNumber, STAT = STAT_)

            if (STAT_ == SUCCESS_) then
            
		        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
		        do i=Me%WorkSize%ILB, Me%WorkSize%IUB                    
                    if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then
                        PropertyX%ConcentrationOld (i,j) = BasinConcentration (i,j)
                    endif
                enddo
                enddo                
            else
                write(*,*) 'Looking for Runoff Property in Runoff Property ???', GetPropertyName(PropertyXIDNumber)
                write(*,*) 'but not found. Link between WQ in modules can not be done.'
                stop 'SetBasinConcRP - ModuleDrainageNetwork - ERR010'
            end if

            call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_)
            if (STAT_ /= SUCCESS_) stop 'SetBasinConcRP - ModuleRunoffProperties - ERR010'

        else
            STAT_ = ready_
        end if

        STAT = STAT_

    end subroutine SetBasinConcRP 

    !---------------------------------------------------------------------------

    
    subroutine GetRPnProperties (RunoffPropertiesID, nProperties, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer                                         :: nProperties
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_
        !-----------------------------------------------------------------------


        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            nProperties       = Me%PropertiesNumber
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetRPnProperties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    subroutine GetRPPropertiesIDByIdx (RunoffPropertiesID, Idx, ID,PropAdvDiff, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: RunoffPropertiesID
        integer, intent(IN)                             :: Idx
        integer, intent(OUT)                            :: ID
        logical, intent(OUT)                            :: PropAdvDiff
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_, i
        type (T_Property), pointer                      :: CurrProp

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(RunoffPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            CurrProp => Me%FirstProperty
            do i = 1, idx - 1
                CurrProp => CurrProp%Next
            enddo

            ID        = CurrProp%ID%IDNumber
            PropAdvDiff = CurrProp%Evolution%AdvectionDiffusion
            
            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetRPPropertiesIDByIdx

    !---------------------------------------------------------------------------

    subroutine UnGetRunoffProperties2D_I(ObjRunoffPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID
        integer, dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRunoffProperties_, Me%InstanceID, "UnGetRunoffProperties2D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunoffProperties2D_I

    !--------------------------------------------------------------------------

    subroutine UnGetRunoffProperties2D_R4(ObjRunoffPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRunoffProperties_, Me%InstanceID,  "UnGetRunoffProperties2D_R4")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunoffProperties2D_R4


    !--------------------------------------------------------------------------

    subroutine UnGetRunoffProperties2D_R8(ObjRunoffPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjRunoffPropertiesID
        real(8), dimension(:, :), pointer                :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mRunoffProperties_, Me%InstanceID,  "UnGetRunoffProperties2D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetRunoffProperties2D_R8


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyRunoffProperties(ObjRunoffPropertiesID,          &
                                           STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunoffPropertiesID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_,STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunoffProperties - ModuleRunoffProperties - ERR02'
            
            !Actualize the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyRunoffProperties - ModuleRunoffProperties - ERR03'
                      
            call ReadLockExternalVar
            
            !Eduardo Jauch
            !Actualize properties if evolution from file
            call ActualizePropertiesFromFile
            
            if (Me%Coupled%AdvectionDiffusion) then

                !sources and sinks from porous media and Drainage network
!                call InterfaceFluxes
                
				if (Me%AdvDiff_Module == AdvDif_ModuleRP_) then !Advection and Diffusion from ModuleRunoffProperties
                    
                    call AdvectionDiffusionProcesses_RP
 
!               elseif (Me%AdvDiff_Module == AdvDif_ModuleAD_) then !Advection and Diffusion from ModuleAdvectionDiffusion
!            
!                   call AdvectionDiffusionProcesses_AD

    	        endif
            
            endif

!            if (Me%Coupled%SoilQuality) then
!                call SoilQualityProcesses
!            endif
!
!#ifdef _PHREEQC_
!            if (Me%Coupled%SoilChemistry) then
!                call SoilChemistryProcesses
!            endif
!#endif            

            if (Me%Coupled%MinConcentration) then
                call SetLimitsConcentration 
            endif

            if (Me%Output%Timeserie_ON) then
                call OutPut_TimeSeries
            endif

            if (Me%Output%HDF_ON) then
                call OutPut_HDF
            endif


            call Actualize_Time_Evolution
        
            call ReadUnlockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyRunoffProperties

    !-----------------------------------------------------------------------------
    
    !This method implies that there is no diffusivity in top face (interaction with water column)
    subroutine ModifyDiffusivity_Old(PropertyX)

        !Arguments-------------------------------------------------------------
        type (T_Property), pointer           :: PropertyX

        !External--------------------------------------------------------------
        integer                              :: i, j,  CHUNK

        !Begin----------------------------------------------------------------------

        call ModifyTurbulence(PropertyX)

        
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J)

               
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
                !Vertical Diffusivity
                !m2/s = m2/s * [-]  + m2/s 
                PropertyX%Diffusivity(i,j) = 0.0
                
                !Horizontal Diffusivity is called viscosity to maintain the format of when is called in Module Advection Diffusion
                !m2/s = m2/s * [-]  + m2/s
                PropertyX%Viscosity(i,j)   = PropertyX%Evolution%AdvDiff%Molecular_Diff_Coef &
                                             + PropertyX%Diff_Turbulence_H (i,j)
                
                if (Me%AdvDiff_Module == AdvDif_ModuleRP_) then ! need ViscosityU and ViscosityV
                    PropertyX%ViscosityU(i,j) = PropertyX%Viscosity(i,j) 
                    PropertyX%ViscosityV(i,j) = PropertyX%Viscosity(i,j)
                    
                endif
                                                                       

            endif

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
        integer                    :: i, j, CHUNK
        real                       :: VelMedU, VelMedV, VelMedW, VelU, VelV

        !Begin----------------------------------------------------------------------

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J)
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
                
!                VelMedU = 0.5 * (Me%ExtVar%UnsatU(i,j,k) * Me%ExtVar%ComputeFacesU3D(I,J,K) +  Me%ExtVar%UnsatU(i,j+1,k) * Me%ExtVar%ComputeFacesU3D(I,J+1,K))
!                VelMedV = 0.5 * (Me%ExtVar%UnsatV(i,j,k) * Me%ExtVar%ComputeFacesV3D(I,J,K) +  Me%ExtVar%UnsatV(i+1,j,k) * Me%ExtVar%ComputeFacesV3D(I+1,J+,K))
! 

!                !m2/s = m * m/s
!                PropertyX%Diff_Turbulence_H(i,j,k) = Me%Disper_Longi%Field(i,j,k) *  abs(0.5 * (VelMedU + VelMedV)                &
!                                                    + Me%Disper_Trans%Field(i,j,k) * abs(VelMedW)
                VelU = 0.5 * (Me%ExtVar%CenterVelU(i,j)+Me%ExtVar%CenterVelU(i,j-1))
                VelV = 0.5 * (Me%ExtVar%CenterVelV(i,j)+Me%ExtVar%CenterVelV(i-1,j))
                !m2/s = m * m/s
                PropertyX%Diff_Turbulence_H(i,j) = Me%Disper_Longi%Field(i,j) *                              &
		                                                 abs(VelV + VelU  / 2.)
               
            endif

        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        
        
        !---------------------------------------------------------------------------

    end subroutine ModifyTurbulence

    !--------------------------------------------------------------------------
    
!	subroutine AdvectionDiffusionTopBoundary
!
!        !Local--------------------------------------------------------------
!        type (T_Property), pointer         :: PropertyX
!        integer                            :: i, j, k, CHUNK
!        real(8), dimension(:,:), pointer   :: WaterVolume
!        real(8)                            :: InfVolume
!
!        !Begin----------------------------------------------------------------------
!
!        !CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
!        
!        !!$OMP PARALLEL PRIVATE(I,J,K,InfVolume)
!
!        PropertyX => Me%FirstProperty
!
!        if (Me%AdvDiff_Explicit) then
!            WaterVolume => Me%WaterVolumeCorr
!        else
!            WaterVolume => Me%WaterVolumeOld
!        endif
!
!do1:    do while (associated(PropertyX))
!
!			if (PropertyX%Evolution%AdvectionDiffusion) then
!
!                !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!		        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
!		        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
!
!		        	if (Me%ExtVar%BasinPoints(i, j) == BasinPoint) then
!
!!    					call ModifyWaterColumnConcentration (i, j, PropertyX)
!!                        if (Me%AdvDiff_Explicit) then
!!                            InfVolume = 0.0
!!                        else
!!!!!!!!!!!!!!!!!                            InfVolume = abs(Me%ExtVar%FluxW(i, j, k+1)) * Me%ExtVar%DT
!!                        endif
!    					
!!!!    					if (Me%ExtVar%FluxW(i, j, k+1) < 0) then						
!!!!                            
!!!!                            !g/m3  = (g/m3 * m3 + g) /((+ m3 + m3))
!!!!		        			PropertyX%Concentration(i, j) = ((PropertyX%Concentration(i, j) * WaterVolume(i, j, k))            &
!!!!		        			                                    + PropertyX%MassOnFluxW(i, j)) / ((WaterVolume(i, j, k) + InfVolume))
!!!!                        else
!!!!                            !g/m3  = (g/m3 * m3 + g) /((+ m3 + m3))
!!!!		        			PropertyX%Concentration(i, j) = ((PropertyX%Concentration(i, j) * WaterVolume(i, j, k))            &
!!!!		        			                                    - PropertyX%MassOnFluxW(i, j)) / ((WaterVolume(i, j, k) - InfVolume))
!!!!                        		        			                                    
!!!!						endif
!
!		        	endif
!
!		        enddo
!		        enddo
!                !!$OMP END DO
!                
!			endif
!
!			PropertyX => PropertyX%Next
!
!	    enddo do1
!        
!        !!$OMP END PARALLEL
!
!	end subroutine AdvectionDiffusionTopBoundary

    !--------------------------------------------------------------------------

!    subroutine ComputeVolumes
!
!        !Local-----------------------------------------------------------------
!        integer :: i, j, CHUNK        
!
!        !----------------------------------------------------------------------
!        
!        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
!        !$OMP PARALLEL PRIVATE(I,J)
!        
!        !Compute volumes and correct top volume taking FluxW(KUB+1) because it would be interpreted by module advection diffusion
!        !as an additional water flux with the conc of C(i,j,k)
!        
!        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
!        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if (Me%ExtVar%BasinPoints(i,j) == BasinPoint) then             
!                 
! !               Me%WaterVolumeOld(i,j)     = Me%ExtVar%WaterContentOld(i,j) * Me%ExtVar%Cellvolume(i,j)
! !               Me%WaterVolume(i,j)        = Me%ExtVar%WaterContent(i,j) * Me%ExtVar%Cellvolume(i,j)
! !               Me%WaterVolumeCorr(i,j)    = Me%WaterVolume(i,j)
! !               if (k == Me%WorkSize%KUB) then
! !                   !m3 = m3 - m3/s * s
! !                   Me%WaterVolumeCorr(i,j) = Me%WaterVolumeCorr(i,j) + Me%ExtVar%FluxW(i,j)  * Me%ExtVar%DT
! !               endif
!
!            endif
!        enddo
!        enddo
!        !$OMP END DO
!        !$OMP END PARALLEL 
!               
!        !Correct fluxw - take FluxW(KUB+1) because it would be interpreted by module advection diffusion
!        !as an additional water flux with the conc of C(i,j,k)
!
!
!        
! !       k = Me%WorkSize%KUB
! !       
! !       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
! !       do i = Me%WorkSize%ILB, Me%WorkSize%IUB
! !           if (Me%ExtVar%OpenPoints3D(i,j,k) == OpenPoint) then             
! !           
! !               Me%WaterColumnVolumeOld(i,j) = Me%WaterColumnVolume(i, j)
! !               Me%WaterColumnVolume(i, j) = Me%ExtVar%WaterColumn(i, j) * Me%ExtVar%Area(i, j)
! !
! !           endif
! !       enddo
! !       enddo
!        
!           
!    
!    end subroutine ComputeVolumes
!    
   !----------------------------------------------------------------------

!    subroutine AdvectionDiffusionProcesses_AD
!    
!        !External--------------------------------------------------------------
!        integer                             :: STAT_CALL    
!        real(8), pointer, dimension(:,:,:)  :: AdvFluxX
!        real(8), pointer, dimension(:,:,:)  :: AdvFluxY
!        real(8), pointer, dimension(:,:,:)  :: AdvFluxZ
!        real(8), pointer, dimension(:,:,:)  :: DifFluxX
!        real(8), pointer, dimension(:,:,:)  :: DifFluxY
!        real(8), pointer, dimension(:,:,:)  :: DifFluxZ
!
!        !Local-----------------------------------------------------------------
!        type(T_Property), pointer           :: Property
!        type (T_Time)                       :: Actual
!        real                                :: ImpExp_AdvXX, ImpExp_AdvYY           
!        integer                             :: i, j, k
!        real                                :: AdvectionV_imp_exp  
!        real                                :: DiffusionV_imp_exp  
!        real                                :: AdvectionH_imp_exp  
!        real(8)                             :: f
!        real(8), dimension(:,:,:), pointer  :: FluxW
!
!        !----------------------------------------------------------------------      
!        Actual = Me%ExtVar%Now
!
!        call ComputeVolumes
!!        call ModifyDiffusivity
!    
!        Property => Me%FirstProperty
!!        call ModifyTopBoundaryCondition (Property) 
!
!do1 :   do while (associated(Property))
!
!cd1 :       if (Property%Evolution%AdvectionDiffusion) then
!
!                call ModifyDiffusivity_Old(Property)
!
!                if (Me%AdvDiff_Explicit) then
!
!                    AdvectionV_imp_exp = ExplicitScheme
!                    DiffusionV_imp_exp = ExplicitScheme
!                    AdvectionH_imp_exp = ExplicitScheme
!                
!                else
!
!                    AdvectionV_imp_exp = ImplicitScheme
!                    DiffusionV_imp_exp = ImplicitScheme
!                    AdvectionH_imp_exp = ImplicitScheme
!
!                endif
!                    
!                if(AdvectionH_imp_exp == ImplicitScheme) then
!
!                    if(Property%Evolution%AdvDiff%ImplicitH_Direction == DirectionX)then
!                                                   
!                        !Direction X implicit
!                        ImpExp_AdvXX = ImplicitScheme 
!                        ImpExp_AdvYY = ExplicitScheme 
!
!                        Property%Evolution%AdvDiff%ImplicitH_Direction = DirectionY
!
!                    else 
!                    
!                        !Direction Y implicit
!                        ImpExp_AdvXX = ExplicitScheme 
!                        ImpExp_AdvYY = ImplicitScheme 
!
!                        Property%Evolution%AdvDiff%ImplicitH_Direction = DirectionX
!
!                    endif 
!            
!                else ! Horizontal Advection Explicit
!
!                    ImpExp_AdvXX = ExplicitScheme 
!                    ImpExp_AdvYY = ExplicitScheme 
!
!                endif
!
!                !for debug purposes, delete before any other tests
!!                Property%Viscosity = 0.
!!                Property%Diffusivity = 0.
!!                Property%Evolution%AdvDiff%BoundaryCondition = 0
!!                f = Me%ExtVar%FluxW(2, 2, 11)
!!                Me%ExtVar%FluxW(2, 2, 11) = 0 
!
!				if (.not. Me%AdvDiff_Explicit) then
!
!	            	call AdvectionDiffusionTopBoundary  
!                
!                endif
!
!!                call AdvectionDiffusion(Me%ObjAdvectionDiffusion,                                            &
!!				                        PROP                = Property%Concentration,                        &
!!				                        schmidt_H           = Property%Evolution%AdvDiff%SchmidtNumberH,     &
!!				                        SchmidtCoef_V       = Property%Evolution%AdvDiff%SchmidtCoefV,       &
!!				                        SchmidtBackground_V = Property%Evolution%AdvDiff%SchmidtBackgroundV, &
!!				                        AdvMethodH          = Property%Evolution%AdvDiff%AdvMethodH,         &
!!				                        TVDLimitationH      = Property%Evolution%AdvDiff%TVDLimitationH,     &
!!				                        AdvMethodV          = Property%Evolution%AdvDiff%AdvMethodV,         &
!!				                        TVDLimitationV      = Property%Evolution%AdvDiff%TVDLimitationV,     &
!!				                        Upwind2H            = Property%Evolution%AdvDiff%Upwind2H,           &
!!				                        Upwind2V            = Property%Evolution%AdvDiff%Upwind2V,           &
!!				                        VolumeRelMax        = Property%Evolution%AdvDiff%VolumeRelMax,       &
!!!				                        DTProp              = Property%Evolution%DTInterval,                 &
!!                                        DTProp              = Me%ExtVar%DT,                                  &
!!				                        ImpExp_AdvV         = AdvectionV_imp_exp,                            &
!!				                        ImpExp_DifV         = DiffusionV_imp_exp,                            &
!!				                        ImpExp_AdvXX        = ImpExp_AdvXX,                                  &
!!				                        ImpExp_AdvYY        = ImpExp_AdvYY,                                  &
!!				                        ImpExp_DifH         = Property%Evolution%AdvDiff%DiffusionH_imp_exp, &
!!				                        NullDif             = Property%Evolution%AdvDiff%NullDif,            &
!!				                        Wflux_X             = Me%ExtVar%FluxU,                               &
!!				                        Wflux_Y             = Me%ExtVar%FluxV,                               &
!!!				                        Wflux_Z             = Me%ExtVar%FluxW,                               &
!!!				                        Wflux_Z             = Me%FluxWCorr,                                  &
!!				                        VolumeZOld          = Me%WaterVolumeOld,                             &
!!!				                        VolumeZ             = Me%WaterVolume,                                &
!!!				                        VolumeZOld          = Me%WaterVolumeOldCorr,                         &
!!				                        VolumeZ             = Me%WaterVolumeCorr,                            &
!!!!!				                        OpenPoints2D        = Me%ExtVar%BasinPoints,                         &
!!!!!				                        LandPoints3D        = Me%ExtVar%LandPoints3D,                        &
!!!!!				                        ComputeFacesU2D     = Me%ExtVar%ComputeFacesU2D,                     &
!!!!!				                        ComputeFacesV2D     = Me%ExtVar%ComputeFacesV2D,                     &
!!!!!				                        ComputeFacesW2D     = Me%ExtVar%ComputeFacesW2D,                     &
!!				                        Visc_H              = Property%Viscosity,                            &
!!				                        Diff_V              = Property%Diffusivity,                          &
!!				                        CellFluxes          = .true.,                                        &
!!				                        BoundaryCondition   = Property%Evolution%AdvDiff%BoundaryCondition,  &
!!				                        NumericStability    = Property%Evolution%AdvDiff%NumericStability,   &
!!				                        STAT                = STAT_CALL)
!!
! !               if (STAT_CALL /= SUCCESS_) stop 'AdvectionDiffusionProcesses - ModuleRunoffProperties - ERR10'
!           
!!                Me%ExtVar%FluxW(2, 2, 11) = f
!
!				if (Me%AdvDiff_Explicit) then
!
!	            	call AdvectionDiffusionTopBoundary  
!
!                endif
!           
!            end if cd1
!
!            Property => Property%Next
!
!        end do do1
!        nullify(Property)
!
!        !-------------------------------------------------------------------------    
!        
!    end subroutine AdvectionDiffusionProcesses_AD
        
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
                                       Matrix2D       = PropertyX%Concentration,    &
                                       PointsToFill2D = Me%ExtVar%BasinPoints,    &
                                       STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleRunoffProperties - ERR01'
            
            endif
            
            PropertyX => PropertyX%Next
            
        enddo
    
        !-------------------------------------------------------------------------    
    
    end subroutine ActualizePropertiesFromFile
    
	!-------------------------------------------------------------------------    

 !   subroutine InterfaceFluxes
        !Local--------------------------------------------------------------------
        !Begin--------------------------------------------------------------------

!        if (Me%ExtVar%CoupledVegetation) then
!            if (Me%ExtVar%ComputeVegInterfaceFluxes) then
!                call VegetationInterfaceFluxes
!            endif
!        endif

!        if (Me%ExtVar%CoupledDN) then
!            call DrainageNetworkInterfaceFluxes
!        endif

!        if (Me%ExtVar%CoupledPMP) then
!           call PorousMediaPropertiesInterfaceFluxes
!        endif


!    end subroutine InterfaceFluxes
 
     !-----------------------------------------------------------------------------
    
    ! This routine solves mass sources ans sinks due to drainage network and implicitly take or add mass. 

!    subroutine DrainageNetworkInterfaceFluxes 
!        !Local--------------------------------------------------------------------
!        integer                                 :: STAT_CALL, i, j, k        
!        real, dimension(:, :), pointer          :: FlowToChannels, WaterColumnOld
!        integer, dimension(:, :), pointer       :: GWLayer!, ChannelsID
!        type (T_Property), pointer              :: PropertyX
!!        real, dimension (:), pointer            :: DNConcentration
!        real                                    :: WaterVolumeOld, MassFlow, OldMass, NewMass
!
!        !begin--------------------------------------------------------------------    
!
!        call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkInterfaceFluxes - ModuleRunoffProperties - ERR01'
!
!        call GetFlowToChannels   (Me%ObjRunoff, FlowToChannels, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkInterfaceFluxes - ModuleRunoffProperties - ERR10'
!
!        call GetRunoffWaterColumnOld   (Me%ObjRunoff, WaterColumnOld, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'DrainageNetworkInterfaceFluxes - ModuleRunoffProperties - ERR20'
!
!
!
!        PropertyX => Me%FirstProperty
!        
!        do while (associated(PropertyX))        
!        
!            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!                if (Me%ExtVar%RiverPoints(i,j) == WaterPoint) then
!                    
!                    if (FlowToChannels(i,j) .gt. 0.0) then !transport to channel - looses mass
!                        
!                        !g = m3/s * s * g/m3 
!                        MassFlow = FlowToChannels(i,j) * Me%ExtVar%DT * PropertyX%ConcentrationOld(i,j)
!                        
!                    
!                    else ! transport to soil - gains mass
!                        
!                        !g = m3/s * s * g/m3 
!                        MassFlow = FlowToChannels(i,j) * Me%ExtVar%DT * PropertyX%ConcentrationDN(i,j)
!                        
!                    endif
!                    !m3 = m * m2
!                    WaterVolumeOld = WaterColumnOld(i,j) * Me%ExtVar%Area(i, j)
!                    !g = g/m3 * m3
!                    OldMass        = PropertyX%ConcentrationOld(i,j) * WaterVolumeOld
!                    NewMass        = OldMass - MassFlow
!                    
!                    !g/m3 =  g / (m3Old - (m3/s flux * s)
!                    PropertyX%ConcentrationOld(i,j) = NewMass / (WaterVolumeOld - (FlowToChannels(i,j) * Me%ExtVar%DT))
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
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR040'               
!
!        call UnGetRunoff (Me%ObjRunoff, FlowToChannels, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR050'               
!
!      
!        call UnGetRunoff (Me%ObjRunoff, WaterColumnOld, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR060'               
!
!           
!        
!    end subroutine DrainageNetworkInterfaceFluxes
    
    !-----------------------------------------------------------------------------


    subroutine AdvectionDiffusionProcesses_RP

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

    end subroutine AdvectionDiffusionProcesses_RP
    
    !-----------------------------------------------------------------------------

    subroutine ModifyAdvectionDiffusion (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k!, CHUNK
        real(8)                                     :: Area_Vertical, Area_Top_U, Area_Top_V
        real(8)                                     :: Area_Bottom_U, Area_Bottom_V
        real(8)                                     :: AdvTermB_Top, AdvTermC
        real(8)                                     :: AdvTermA_U, AdvTermA_V
        real(8)                                     :: AdvTermB_Top_U, AdvTermB_Top_V
        real(8)                                     :: AdvTermB_Bottom_U, AdvTermB_Bottom_V
        real(8)                                     :: AdvTermC_U, AdvTermC_V
        real(8)                                     :: DifTerm_Top_U, DifTerm_Top_V
        real(8)                                     :: DifTerm_Bottom_U, DifTerm_Bottom_V         
        real(8)                                     :: aux 
        real(8)                                     :: cofA_U,cofB_U,cofC_U
        real(8)                                     :: cofA_V,cofB_V,cofC_V, CofB
        real(8)                                     :: cofInterfaceDN, ConcInInterfaceDN
        real                                        :: ConcTop
        real(8), pointer, dimension(:,:)            :: FluxU, FluxV
        real(8), pointer, dimension(:,:  )          :: DZX, DZY, DXX, DYY
        real(8), pointer, dimension(:,:)            :: WaterColumn, WaterColumnOld
        integer                                     :: STAT_
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
	        stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR10'                            
        endif
        
        
        FluxU          => Me%ExtVar%FluxU
        FluxV          => Me%ExtVar%FluxV        
        DZX            => Me%ExtVar%DZX
        DZY            => Me%ExtVar%DZY
        DXX            => Me%ExtVar%DXX
        DYY            => Me%ExtVar%DYY
        WaterColumn    => Me%ExtVar%WaterColumn
        WaterColumnOld => Me%ExtVar%WaterColumnOld

        !!!$OMP PARALLEL PRIVATE(I,J,K)
        !!!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
                
                if (WaterColumn(i,j) .gt. AlmostZero) then
                    Area_Vertical      = Me%ExtVar%Area(i, j)
                    !m2 = WC m * face m
                    Area_Top_U         = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i,j+1))) * DYY(i,j+1)
                    Area_Top_V         = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i+1,j))) * DXX(i+1,j)
                    Area_Bottom_U      = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i,j-1))) * DYY(i,j  )
                    Area_Bottom_V      = (0.5 * (WaterColumnOld(i,j) + WaterColumnOld(i-1,j))) * DXX(i,j  )
                    
                    !s/m3 = s / (m * m2)
                    aux      = (Me%ExtVar%DT/(WaterColumn(i,j)* Area_Vertical))
                   
                   !!!FLUXES WITH Drainage Network
                   CofInterfaceDN    = 0.0
                   ConcInInterfaceDN = 0.0
                    if (Me%ExtVar%CoupledDN) then
                       !Flux between river and runoff
                        if (Me%ExtVar%RiverPoints(I,J) == BasinPoint) then                        
                            
                            ! Positive flow -> looses mass
                            cofInterfaceDN = - aux * Me%ExtVar%FlowToChannels(i,j)
                            
                            ! mass going to channel -> conc from runoff
                            if (Me%ExtVar%FlowToChannels(i,j) .gt. 0.0) then
                                ConcInInterfaceDN =  CurrProperty%ConcentrationOld(i,j)
                            
                            !mass coming from channel -> conc from DN
                            elseif (Me%ExtVar%FlowToChannels(i,j) .lt. 0.0) then
                                ConcInInterfaceDN = CurrProperty%ConcentrationDN(i,j)
                            endif
                        
                        endif
                    endif
                    
                    !!DISCHARGE FLUXES IN RUNOFF (NOT YET DONE)
                    
                    !!BOUNDARY FLUXES IN RUNOFF (NOT YET DONE)
                    
                    !!FLUXES IN X AND Y DIRECTION                        
                    if (Me%AdvDiff_SpatialMethod==AdvDif_CentralDif_) then ! diferenças centrais

                        
                        AdvTermA_U        = (aux * FluxU(i,j  ) / 2.) 
                        AdvTermA_V        = (aux * FluxV(i,j  ) / 2.)
                        AdvTermB_Top_U    = (aux * FluxU(i,j+1) / 2.) 
                        AdvTermB_Bottom_U = (aux * FluxU(i,j  ) / 2.)
                        AdvTermB_Top_V    = (aux * FluxV(i+1,j) / 2.) 
                        AdvTermB_Bottom_V = (aux * FluxV(i  ,j) / 2.)
                        AdvTermC_U        = (aux * FluxU(i,j+1) / 2.) 
                        AdvTermC_V        = (aux * FluxV(i+1,j) / 2.)
                        

                    elseif (Me%AdvDiff_SpatialMethod==AdvDif_Upwind_) then ! upwind

                        !DirecU face j
                        if (FluxU(i,j) .lt. 0.0) then !Left face, Negative - exiting
                            AdvTermA_U        = 0.0
                            AdvTermB_Bottom_U = aux * FluxU(i,j)
                        else !Positive - entering or zero.
                            AdvTermA_U        = aux * FluxU(i,j)
                            AdvTermB_Bottom_U = 0.0
                        endif

                        !DirecV face i
                        if (FluxV(i,j) .lt. 0.0) then !Left face, Negative - exiting
                            AdvTermA_V        = 0.0
                            AdvTermB_Bottom_V = aux * FluxV(i,j)
                        else !Positive - entering or zero.
                            AdvTermA_V        = aux * FluxV(i,j)
                            AdvTermB_Bottom_V = 0.0
                        endif
                        
                        !DirecU face j+1
                        if (FluxU(i,j+1) .lt. 0.0) then !Right face, Negative - entering
                            AdvTermC_U        = aux * FluxU(i,j+1)
                            AdvTermB_Top_U    = 0.0
                        else !Positive - exiting or zero.
                            AdvTermC_U        = 0.0
                            AdvTermB_Top_U    = aux * FluxU(i,j+1)
                        endif                    

                        !DirecV face i+1
                        if (FluxV(i+1,j) .lt. 0.0) then !Right face, Negative - entering
                            AdvTermC_V        = aux * FluxV(i+1,j)
                            AdvTermB_Top_V    = 0.0
                        else !Positive - exiting or zero.
                            AdvTermC_V        = 0.0
                            AdvTermB_Top_V    = aux * FluxV(i+1,j)
                        endif                    
                           

                    endif
                    
                    DifTerm_Top_U    = CurrProperty%ViscosityU(i,j+1) * Area_Top_U    * aux / DZX(i,j  )
                    DifTerm_Bottom_U = CurrProperty%ViscosityU(i,j  ) * Area_Bottom_U * aux / DZX(i,j-1)
                    DifTerm_Top_V    = CurrProperty%ViscosityV(i+1,j) * Area_Top_V    * aux / DZY(i  ,j)
                    DifTerm_Bottom_V = CurrProperty%ViscosityV(i  ,j) * Area_Bottom_V * aux / DZY(i-1,j)
                    
                    cofA_U = AdvTermA_U                                                          &
                              + DifTerm_Bottom_U 
                    
                    cofA_V = AdvTermA_V                                                          &
                              + DifTerm_Bottom_V 
                
                    cofB_U = - AdvTermB_Top_U                                                    &
                             + AdvTermB_Bottom_U                                                 &
                             - DifTerm_Bottom_U                                                  &
                             - DifTerm_Top_U
                    
                    cofB_V = - AdvTermB_Top_V                                                    &
                             + AdvTermB_Bottom_V                                                 &
                             - DifTerm_Bottom_V                                                  &
                             - DifTerm_Top_V        
                     
                    cofC_U = - AdvTermC_U                                                        &
                             + DifTerm_Top_U          

                    cofC_V = - AdvTermC_V                                                        &
                             + DifTerm_Top_V    
                    
                    CofB = ((WaterColumnOld(i,j)*Area_Vertical) / (WaterColumn(i,j)*Area_Vertical)) + cofB_U + cofB_V
                    
                    CurrProperty%Concentration(i,j)=  cofA_U * CurrProperty%ConcentrationOld(i,j-1)  &
                                                     + cofC_U * CurrProperty%ConcentrationOld(i,j+1) &
                                                     + cofA_V * CurrProperty%ConcentrationOld(i-1,j) &
                                                     + cofC_V * CurrProperty%ConcentrationOld(i+1,j) &
                                                     + cofB   * CurrProperty%ConcentrationOld(i,j)   &
                                                     + CofInterfaceDN * ConcInInterfaceDN  
                else !No Volume
                    CurrProperty%Concentration(i,j) = 0.0                                                         
                endif
                
                !Check if any of the coeffs get a negative value. If true, stop program
                if ((Me%AdvDiff_CheckCoefs) .AND. ((cofA_U < 0.0) .OR. (cofB < 0.0) .OR. (cofC_U < 0.0) .OR.    &
                                                   (cofA_V < 0.0) .OR. (cofC_V < 0.0) )) then
                                                   
                    call LogCoefs(i,j,cofA_U,cofB,cofC_U,cofA_V,cofC_V)
                
                endif
              
            endif
        enddo
        enddo
        !!!$OMP END DO
        !!!$OMP END PARALLEL

                        
        call SetMatrixValue (CurrProperty%ConcentrationOld,      Me%Size,   CurrProperty%Concentration,          Me%ExtVar%BasinPoints)



    end subroutine ModifyAdvectionDiffusion
    
    !---------------------------------------------------------------------------
    
    subroutine ModifyDiffusivity_New(CurrProperty)
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        !Local-----------------------------------------------------------------
        integer                                     :: I,J, CHUNK
        real                                        :: DiffCoef, velU, VelV        
        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(I,J,DiffCoef,VelU,VelV)

        DiffCoef  = CurrProperty%Evolution%AdvDiff%Molecular_Diff_Coef
        
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
                
                if ((Me%ExtVar%WaterColumnOld(i,j) .gt. AlmostZero) .and. (Me%ExtVar%WaterColumnOld(i,j-1) .gt. AlmostZero)) then
                    VelU = 0.5 * (Me%ExtVar%CenterVelU(i,j)+Me%ExtVar%CenterVelU(i,j-1))
                
                    CurrProperty%ViscosityU(i,j)  = (DiffCoef +(abs(VelU) * Me%Disper_Trans%Field(i,j)))
                else
                    CurrProperty%ViscosityU(i,j)  = 0.0
                endif                                                   
                
                if ((Me%ExtVar%WaterColumnOld(i,j) .gt. AlmostZero) .and. (Me%ExtVar%WaterColumnOld(i-1,j) .gt. AlmostZero)) then
                    VelV = 0.5 * (Me%ExtVar%CenterVelV(i,j)+Me%ExtVar%CenterVelV(i-1,j))
                    
                    CurrProperty%ViscosityV(i,j)  = (DiffCoef +(abs(VelV) * Me%Disper_Trans%Field(i,j)))
                else
                    CurrProperty%ViscosityV(i,j)  = 0.0
                endif                                                    

            endif                                                                        
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL


    
    end subroutine ModifyDiffusivity_New

    !---------------------------------------------------------------------------
    
!    subroutine SoilQualityProcesses
!
!        !External--------------------------------------------------------------
!        integer                                 :: STAT_CALL        
!        
!        !Local----------------------------------------------------------------- 
!        type (T_Property),          pointer     :: PropertyX
!        type (T_Property),          pointer     :: SoilDryDensity, Salinity, pH
!        type (T_Property),          pointer     :: IonicStrength, PhosphorusAdsortionIndex
!
!!        type (T_SoilRate),      pointer     :: SoilRateX
!!        integer                                 :: WILB, WIUB, WJLB, WJUB, WKLB, WKUB 
!!        integer                                 :: i, j, k
!        
!        !Begin-----------------------------------------------------------------
!        
!!        WIUB = Me%WorkSize%IUB
!!        WJUB = Me%WorkSize%JUB
!!        WILB = Me%WorkSize%ILB
!!        WJLB = Me%WorkSize%JLB
!!        WKUB = Me%WorkSize%KUB
!!        WKLB = Me%WorkSize%KLB
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "SoilQualityProcesses")
!
!        call ComputeDissolvedToParticulate2D
!        
!        !Properties not modified by sediment quality (not state variables) but needed in argument
!        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property soil dry density not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR00'
!        endif
!
!        call SearchProperty(Salinity, Salinity_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property salinity not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR10'
!        endif
!
!        call SearchProperty(pH, pH_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property pH not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR20'
!        endif
!
!        call SearchProperty(IonicStrength, IonicStrength_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property ionic strength not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR30'
!
!        endif
!
!        call SearchProperty(PhosphorusAdsortionIndex, PhosphorusAdsortionIndex_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property phosphorus sdsortion index not found in porous media properties'
!            stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR40'
!        endif
!
!        call ComputeWindVelocity
! 
!        
!        if (Me%ExtVar%Now .GE. Me%Coupled%SoilQuality_NextCompute) then
!            
!            PropertyX => Me%FirstProperty
!
!            do while(associated(PropertyX))
!                
!
!!                call Modify_Interface(InterfaceID               = Me%ObjInterface,                         &
!!                                      PropertyID                = PropertyX%ID%IDNumber,                   &
!!                                      Concentration             = PropertyX%Concentration,                 &
!!                                      WaterPoints2D             = Me%ExtVar%BasinPoints,                   &
!!                                      OpenPoints2D              = Me%ExtVar%BasinPoints,                   &
!!                                      WaterPercentage           = Me%ExtVar%WaterContent,                  &
!!                                      DissolvedToParticulate3D  = Me%DissolvedToParticulate2D,             &
!!                                      SoilDryDensity            = SoilDryDensity%Concentration,            &
!!                                      Salinity                  = Salinity%Concentration,                  &
!!                                      pH                        = pH%Concentration,                        &
!!                                      IonicStrength             = IonicStrength%Concentration,             &
!!                                      PhosphorusAdsortionIndex  = PhosphorusAdsortionIndex%Concentration,  &
!!                                      WindVelocity              = Me%ExtVar%WindVelocity,                  &
!!                                      STAT                      = STAT_CALL)
!!                if (STAT_CALL .NE. SUCCESS_)                                                               &
!!                    stop 'SoilQualityProcesses - ModuleRunoffProperties - ERR01'
!                
!
!                PropertyX => PropertyX%Next
!                
!
!            end do
!            
!            Me%Coupled%SoilQuality_NextCompute = Me%Coupled%SoilQuality_NextCompute +       &
!                                                     Me%Coupled%SoilQuality_DT
!
!        end if
!
!        PropertyX => Me%FirstProperty
!
!        do while(associated(PropertyX))
!
!
!            if (PropertyX%Evolution%SoilQuality) then
!
!                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then
!
!                    call Modify_Interface(InterfaceID   = Me%ObjInterface,                  &
!                                          PropertyID    = PropertyX%ID%IDNumber,            &
!                                          Concentration = PropertyX%Concentration,          &
!                                          WaterPoints2D = Me%ExtVar%BasinPoints,          &
!                                          DTProp        = PropertyX%Evolution%DTInterval,   &
!                                          STAT          = STAT_CALL)
!                    if (STAT_CALL .NE. SUCCESS_)                                            &
!                        stop 'SoilQuality_Processes - ModuleRunoffProperties - ERR02'
!
!                end if
!
!            end if
!
!            PropertyX => PropertyX%Next
!            
!        end do
!
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "SoilQualityProcesses")
!
!    end subroutine SoilQualityProcesses
!    
!    !-----------------------------------------------------------------------------    
!
!    subroutine ComputeWindVelocity
!
!        !Arguments-------------------------------------------------------------
!
!        !External--------------------------------------------------------------
!        integer                                      :: i, j, k
!        
!        !Begin-----------------------------------------------------------------
!
!        call SetMatrixValue(Me%ExtVar%WindVelocity, Me%Size, 0.0, Me%ExtVar%BasinPoints)
!
!        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
!
!            if (Me%ExtVar%BasinPoints(i, j) == 1) then
!                
!                ! km/day                        =  m/s  * 1E-3km/m * 86400s/day 
!                Me%ExtVar%WindVelocity(i,j) = Me%ExtVar%WindVelocity2D(i,j) * 1E-3 * 86400
!            endif
!
!        enddo
!        enddo
!
!    end subroutine ComputeWindVelocity
!
!    !----------------------------------------------------------------------------- 

    
    subroutine SetLimitsConcentration

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        type (T_Property), pointer                  :: Property
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: i, j, k!, CHUNK
        
        !Begin----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "SetLimitsConcentration")

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 

        Property => Me%FirstProperty  

do1 :   do while (associated(Property))
cd1 :       if (Property%Evolution%MinConcentration) then
                
!                CHUNK = CHUNK_J(Me%Size%JLB, Me%Size%JUB)
                
!                !$OMP PARALLEL SHARED(CHUNK, Property) PRIVATE(I,J)
!                !$OMP DO SCHEDULE(DYNAMIC, CHUNK)

                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
                do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExtVar%BasinPoints(i, j) == 1) then
    
                        if (Property%Concentration(i, j) < Property%MinValue) then
                            
                            ! mass created (g) = g + (g/m3)* (m * m2)
                            Property%Mass_created(i, j) = Property%Mass_Created(i, j)   +  &
                                                   (Property%MinValue                -  &
                                                    Property%Concentration(i, j)) *  (Me%ExtVar%WaterColumn(i,j) * &
                                                    Me%ExtVar%Area (i, j))

                            Property%Concentration(i, j) = Property%MinValue
                            
                        endif

                    endif

                enddo
                enddo
                
!                !$OMP END DO NOWAIT
!                !$OMP END PARALLEL
                
            endif cd1
                
        Property => Property%Next
        end do do1

        nullify(Property)

        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "SetLimitsConcentration")


    end subroutine SetLimitsConcentration

    !--------------------------------------------------------------------------
    
!    subroutine ComputeDissolvedToParticulate2D
!
!        !External--------------------------------------------------------------
!        integer                                 :: STAT_CALL        
!         
!        !Local----------------------------------------------------------------- 
!        integer                                 :: i, j, k
!        real                                    :: DT, InstantValue, ResidualValue
!        type(T_Property), pointer               :: SoilDryDensity!, DrySedimentVolume
!        !Begin-----------------------------------------------------------------
!        
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "ComputeDissolvedToParticulate3D")
!
!
!        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_)  &
!            stop 'ComputeDissolvedToParticulate3D - ModuleRunoffProperties - ERR01'
!
!        call SearchProperty(SoilDryDensity, SoilDryDensity_        , .false., STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) stop 'ComputeDissolvedToParticulate3D - ModuleRunoffProperties - ERR02'
!
!
!        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
!            if(Me%ExtVar%BasinPoints(i,j) .eq. BasinPoint)then
!
!                ! [m3water/kgsed] = [m3water/m3cell]*[m3cell] / ([kgsed/m3cell] * [m3cell]) 
!                InstantValue                       = Me%ExtVar%WaterContent(i,j) *  Me%ExtVar%CellVolume(i,j) / &
!                                                    (SoilDryDensity%Concentration(i,j) * Me%ExtVar%CellVolume(i,j))
!
!                ResidualValue                      = Me%DissolvedToParticulate2D(i,j)
!
!                Me%DissolvedToParticulate2D(i,j) = (ResidualValue * Me%ResidualTime +     &
!                                                      InstantValue * DT) / (Me%ResidualTime + DT)
!                                                       
!            end if
!        end do
!        end do
!
!        Me%ResidualTime = Me%ResidualTime + DT
!        
!        
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "ComputeDissolvedToParticulate3D")
!
!
!    end subroutine ComputeDissolvedToParticulate2D
!    !--------------------------------------------------------------------------

!#ifdef _PHREEQC_
!    !--------------------------------------------------------------------------
!    subroutine SoilChemistryProcesses
!
!        !External--------------------------------------------------------------
!        integer :: STAT_CALL        
!        
!        !Local----------------------------------------------------------------- 
!        type (T_Property), pointer                   :: PropertyX, pH
!        integer                                      :: I, J, K
!        real             , pointer, dimension(:,:,:) :: Theta
!        real(8)          , pointer, dimension(:,:,:) :: Volume
!        
!        !Begin-----------------------------------------------------------------
!        if (MonitorPerformance) call StartWatch ("ModuleRunoffProperties", "SoilChemistryProcesses")
!
!        Theta  => Me%ExtVar%WaterContent
!        Volume => Me%ExtVar%Cellvolume
!        
!        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
!        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!        
!            if (Me%ExtVar%BasinPoints(I,J) == BasinPoint) then
!                        
!                Me%ExtVar%CellWaterMass(I, J, K) = Theta(I, J, K) * Volume(I, J, K) * WaterReferenceDensity    
!                        
!            end if
!                    
!        end do
!        end do
!        end do
!        
!
!        call SearchProperty(pH, pH_ , .false., STAT = STAT_CALL)    
!            
!        if (STAT_CALL /= SUCCESS_) then
!            write(*,*) 'property pH not found in porous media properties'
!            stop 'SoilChemistryProcesses - ModuleRunoffProperties - ERR01'
!        endif
!
!        !Question: Is this necessary to any other process than the SoilQualityProcesses?
!        !call ComputeDissolvedToParticulate3D
!        
!        !Question: Is this necessary to any other process than the SoilQualityProcess?
!        !call ComputeWindVelocity
! 
! 
!        !Question: Why here is used SoilChemistry_NextCompute and after is used NextCompute?
!        !          This is related to the possibility of different DT's beetwen 0D model and MOHID general DT?
!        !          
!        if (Me%ExtVar%Now .GE. Me%Coupled%SoilChemistry_NextCompute) then
!            
!            PropertyX => Me%FirstProperty
!
!            do while(associated(PropertyX))
!                
!                call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry , &
!                                      PropertyID    = PropertyX%ID%IDNumber        , &
!                                      WaterMass     = Me%ExtVar%CellWaterMass      , &                                                                       
!                                      Concentration = PropertyX%Concentration      , &
!                                      WaterPoints2D = Me%ExtVar%BasinPoints      , &
!                                      pH            = pH%Concentration             , &
!                                      OpenPoints2D  = Me%ExtVar%BasinPoints       , &
!                                      STAT          = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_) &
!                    stop 'SoilChemistryProcesses - ModuleRunoffProperties - ERR02'
!                
!                PropertyX => PropertyX%Next
!                
!            end do
!            
!            Me%Coupled%SoilChemistry_NextCompute = Me%Coupled%SoilChemistry_NextCompute + &
!                                                   Me%Coupled%SoilChemistry_DT
!
!        end if
!
!        PropertyX => Me%FirstProperty
!
!        do while(associated(PropertyX))
!
!            if (PropertyX%Evolution%SoilChemistry) then
!
!                if (Me%ExtVar%Now .GE. PropertyX%Evolution%NextCompute) then
!
!                    call Modify_Interface(InterfaceID   = Me%ObjInterfaceSoilChemistry,                  &
!                                          PropertyID    = PropertyX%ID%IDNumber,            &
!                                          Concentration = PropertyX%Concentration,          &
!                                          WaterPoints2D = Me%ExtVar%BasinPoints,            &
!                                          DTProp        = PropertyX%Evolution%DTInterval,   &
!                                          STAT          = STAT_CALL)
!                    if (STAT_CALL .NE. SUCCESS_)                                            &
!                        stop 'SoilChemistryProcesses - ModuleRunoffProperties - ERR03'
!
!                end if
!
!            end if
!
!            PropertyX => PropertyX%Next
!            
!        end do
!
!        if (MonitorPerformance) call StopWatch ("ModuleRunoffProperties", "SoilChemistryProcesses")
!        
!        !End-------------------------------------------------------------------
!        
!    end subroutine SoilChemistryProcesses   
!    !-----------------------------------------------------------------------------    
!#endif

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
                                    Data2D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModuleRunoffProperties - ERR01'

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
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR00'

                    call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",      &
                                         Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR01'
           
                    Me%LastOutPutHDF5 = Actual
       
                endif First

                !Sets limits for next write operations
                call HDF5SetLimits   (Me%ObjHDF5,                                &
                                      Me%WorkSize%ILB,                           &
                                      Me%WorkSize%IUB,                           &
                                      Me%WorkSize%JLB,                           &
                                      Me%WorkSize%JUB,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR02'

                !Writes the Open Points
                call HDF5WriteData   (Me%ObjHDF5, "//Grid/BasinPoints",              &
                                      "OpenPoints", "-",                            &
                                      Array2D = Me%ExtVar%BasinPoints,             &
                                      OutputNumber = OutPutNumber,                  &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR03'


                PropertyX => Me%FirstProperty
                do while (associated(PropertyX))

                    if (PropertyX%OutputHDF) then
 
                        call HDF5WriteData   (Me%ObjHDF5,                                    &
                                              "/Results/"//trim(PropertyX%ID%Name),          &
                                              trim(PropertyX%ID%Name),                       &
                                              trim(PropertyX%ID%Units),                      &
                                              Array2D = PropertyX%Concentration,             &
                                              OutputNumber = OutPutNumber,                   &
                                              STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR04'

                    endif

                    PropertyX => PropertyX%Next

                enddo

                Me%OutPut%NextOutput = OutPutNumber + 1

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutHDF - ModuleRunoffproperties - ERR06'
            
            endif  TOut
        endif  TNum

    end subroutine OutPut_HDF

    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------

    subroutine LogCoefs(i,j,Coef1, Coef2, Coef3, Coef4, Coef5)

        !Arguments-------------------------------------------------------------
        integer                                     :: i,j
        real(8)                                     :: Coef1, Coef2, Coef3, Coef4, Coef5

        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day
        real                                        :: Hour, Minute, Second

        call ExtractDate(Me%ExtVar%Now, Year, Month, Day, Hour, Minute, Second)
        
        write (Me%Files%AsciiUnit, fmt=1000) Year, Month, Day, Hour, Minute, Second,    &
                                             i,j, Coef1, Coef2, Coef3, Coef4, Coef5

        1000 format(f5.0, f5.0, f5.0, f5.0, f5.0, f12.5, i3, i3, 5f13.8)

    end subroutine LogCoefs

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillRunoffProperties(ObjRunoffPropertiesID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjRunoffPropertiesID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers,STAT_CALL  
        type(T_property), pointer           :: PropertyX
        

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjRunoffPropertiesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mRunoffProperties_,  Me%InstanceID)


            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 
                if(PropertyX%ID%SolutionFromFile)then

                    call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillRunoffProperties - ModuleRunoffProperties - ERR00'
                end if
                
                PropertyX => PropertyX%Next
            end do 
            
            if (associated(Me%Disper_Longi%Field))then
                deallocate(Me%Disper_Longi%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillRunoffProperties - ModuleRunoffProperties - ERR01'
                nullify   (Me%Disper_Longi%Field)
            end if

            if (associated(Me%Disper_Trans%Field))then
                deallocate(Me%Disper_Trans%Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                               &
                    stop 'KillRunoffProperties - ModuleRunoffProperties - ERR02'
                nullify   (Me%Disper_Trans%Field)
            end if
            
            if (nUsers == 0) then

                !Kills the TimeSerie
                if (Me%ObjTimeSerie /= 0) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillRunoff - ModuleRunoffProperties - ERR05'
                endif

                
                if (Me%OutPut%HDF_ON) then                    
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
                    if (STAT_ /= SUCCESS_) stop 'KillVegetation - ModuleRunoffProperties  - ERR08'
                endif
                
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR08'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR11'
                
                nUsers = DeassociateInstance (mRunoff_,  Me%ObjRunoff)
                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR12'

!                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
!                if (nUsers == 0) stop 'KillRunoff - ModuleRunoffProperties - ERR13'


!                if ((Me%Coupled%AdvectionDiffusion) .and. (Me%AdvDiff_Module == AdvDif_ModuleAD_)) then
!                    call KillAdvectionDiffusion(Me%ObjAdvectionDiffusion, STAT = STAT_CALL)
!                    
!                    if (STAT_CALL /= SUCCESS_) &
!                        stop 'KillRunoff - ModuleRunoffProperties - ERR15'
!                end if

                
                call DeallocateVariables

                !Deallocates Instance
                call DeallocateInstance ()

                ObjRunoffPropertiesID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillRunoffProperties
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_RunoffProperties), pointer          :: AuxObjRunoffProperties
        type (T_RunoffProperties), pointer          :: PreviousObjRunoffProp

        !Updates pointers
        if (Me%InstanceID == FirstObjRunoffProperties%InstanceID) then
            FirstObjRunoffProperties => FirstObjRunoffProperties%Next
        else
            PreviousObjRunoffProp => FirstObjRunoffProperties
            AuxObjRunoffProperties      => FirstObjRunoffProperties%Next
            do while (AuxObjRunoffProperties%InstanceID /= Me%InstanceID)
                PreviousObjRunoffProp => AuxObjRunoffProperties
                AuxObjRunoffProperties      => AuxObjRunoffProperties%Next
            enddo

            !Now update linked list
            PreviousObjRunoffProp%Next => AuxObjRunoffProperties%Next

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
!        deallocate (Me%DifusionNumber          )
!        deallocate (Me%ReynoldsMNumber         )
!        deallocate (Me%ExtVar%WindVelocity   )
        
!#ifdef _PHREEQC_
!        deallocate (Me%ExtVar%CellWaterMass)
!#endif        

!        deallocate (Me%WaterVolume)
!        deallocate (Me%WaterVolumeOld)
!        deallocate (Me%WaterVolumeCorr)
 !       deallocate (Me%WaterVolumeOldCorr)
      
   
    end subroutine DeallocateVariables 

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjRunoffProperties_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunoffProperties_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjRunoffProperties_ID > 0) then
            call LocateObjRunoffProperties (ObjRunoffProperties_ID)
            ready_ = VerifyReadLock (mRunoffProperties_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjRunoffProperties (ObjRunoffPropertiesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjRunoffPropertiesID

        !Local-----------------------------------------------------------------

        Me => FirstObjRunoffProperties
        do while (associated (Me))
            if (Me%InstanceID == ObjRunoffPropertiesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleRunoffProperties - LocateObjRunoffProperties - ERR01'

    end subroutine LocateObjRunoffProperties

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
                if (PrintWarning) write (*,*)'Property Not Found in Module RunoffProperties ', &
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

        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR01'        
        
        call GetOverLandFlow      (Me%ObjRunoff, Me%ExtVar%FluxU, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR030'

        call GetRunoffCenterVelocity (Me%ObjRunoff, Me%ExtVar%CenterVelU, Me%ExtVar%CenterVelV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR040'


        call GetRunoffWaterColumn     (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR64'

        call GetRunoffWaterColumnOld   (Me%ObjRunoff, Me%ExtVar%WaterColumnOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR65'
        

        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR070'


        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                    &
!                                  DVY         = Me%ExtVar%DVY,                          &
!                                  DUY         = Me%ExtVar%DUY,                          &
                                  DXX         = Me%ExtVar%DXX,                          &
                                  DYY         = Me%ExtVar%DYY,                          &
                                  DZY         = Me%ExtVar%DZY,                          &
                                  DZX         = Me%ExtVar%DZX,                          &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR101'

!        call GetGridData  (Me%ObjGeometry, Me%ExtVar%Topography, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleRunoffProperties - ERR110'

        if (Me%ExtVar%CoupledDN) then

            call GetRiverPoints   (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR01'

            call GetFlowToChannels   (Me%ObjRunoff, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR10'            
            
        endif    


    end subroutine ReadLockExternalVar

    !-----------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        call UngetBasin           (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR01'        

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%FluxU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR030'

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%FluxV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR040'
        
        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%CenterVelU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR050'    
        
        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%CenterVelV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR060'             
        
        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR063'

        call UnGetRunoff           (Me%ObjRunoff, Me%ExtVar%WaterColumnOld, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR064'

        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalVar - ModuleRunoffProperties - ERR070'
        


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR112'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR113'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR114'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR115'

!        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DVY, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR116'

!        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExtVar%DUY, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR117'


!        call UnGetGeometry              (Me%ObjGeometry, Me%ExtVar%Topography,  STAT = STAT_CALL )        
!        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleRunoffProperties - ERR120'

        if (Me%ExtVar%CoupledDN) then

            call UnGetBasin (Me%ObjBasinGeometry, Me%ExtVar%RiverPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR140'               
            
            call UnGetRunoff (Me%ObjRunoff, Me%ExtVar%FlowToChannels, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyAdvectionDiffusion - ModuleRunoffProperties - ERR150'     
  
            
        endif


    endsubroutine ReadUnlockExternalVar


end module ModuleRunoffProperties

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 








