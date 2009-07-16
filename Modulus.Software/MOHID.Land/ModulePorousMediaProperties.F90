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
   implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !PorousMediaProperties
    public  :: ConstructPorousMediaProperties
    private ::      AllocateInstance

    !Selector

    public  :: GetPropInfiltration
    public  :: GetPropEfectiveEVTP 
    public  :: GetNextPorousMediaPropDT
    public  :: GetPorousMediaPropertiesPointer
    public  :: GetPorousMediaPropertiesInteger
    public  :: UnGetPorousMediaProperties
                     
    
    !Modifier
    public  :: ModifyPorousMediaProperties
    private :: CalculateAdvectionDiffusion

    !Destructor
    public  :: KillPorousMediaProperties                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjPorousMediaProperties 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetPorousMediaProperties3D_I
    private :: UnGetPorousMediaProperties3D_R8
    interface  UnGetPorousMediaProperties
        module procedure UnGetPorousMediaProperties3D_I
        module procedure UnGetPorousMediaProperties3D_R8
        module procedure UnGetPorousMediaProperties3D_R8i
    end interface  UnGetPorousMediaProperties

    !Types---------------------------------------------------------------------
    
    private :: T_PorousMediaProperties
    

    !PropI
    character(LEN = StringLength), parameter :: char_nitrite        = trim(adjustl('nitrite'            ))

    !PropII
    character(LEN = StringLength), parameter :: char_nitrate      = trim(adjustl('nitrate'))                        
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
        real   , dimension(:,:), pointer            :: CropEvapotrans       => null()
        real   , dimension(:,:), pointer            :: Transpiration        => null() ! Lúcia nova var
        real   , dimension(:,:), pointer            :: Evaporation          => null() ! Lúcia nova var
        real(8), dimension(:,:), pointer            :: PotentialInfCol      => null()
        integer                                     :: EvapoTranspirationMethod
        real                                        :: PorousMediapropDT
        type(T_Time)                                :: Now
        type(T_Time)                                :: BeginTime
        type(T_Time)                                :: EndTime
   
        ! from porousMedia
        real(8), dimension(:, :), pointer           :: Infiltration
        real(8), dimension(:, :), pointer           :: EfectiveEVTP
        real(8), dimension(:, :), pointer           :: EfectiveEVTP2
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        integer, dimension(:,:), pointer        :: BasinPoints
        real,    pointer, dimension(:,:,:)          :: WaterContentold
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
     end type T_ExtVar

    type T_OutPut
        type (T_Time), pointer, dimension(:)    :: OutTime
        integer                                 :: NextOutPut
        logical                                 :: Yes = .false.
        logical                                 :: TimeSerieON
        logical                                 :: ProfileON
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
        real                                    :: rainconc
        real                                    :: bottomconc
    end type T_Property


    type       T_PorousMediaProperties
        integer                                     :: ObjTime
        integer                                     :: ObjTopography
        integer                                     :: ObjHorizontalGrid
        integer                                     :: ObjHorizontalMap
        integer                                     :: ObjDrainageNetwork
        integer                                     :: ObjBasinGeometry
        integer                                     :: ObjPorousMedia
        integer                                     :: ObjGeometry
        integer                                     :: ObjMap
        integer                                     :: ObjGridData
        integer                                     :: ObjEnterData
        integer                                     :: ObjtimeSerie
        integer                                     :: ObjSedimentQuality
        integer                                     :: Objhdf5
        integer                                     :: ObjBottomTopography
        integer                                     :: ObjProfile
        type (T_ExtVar)                             :: ExtVar
        logical                                     :: CheckGlobalMass      
        type (T_Files)                              :: Files
        type (T_AdvectionDiffusion)                 :: AdvDiff
        type (T_OutPut)                             :: OutPut
        type (T_Property), pointer                  :: FirstProperty    => null() !Lúcia
        type(T_Property    ), pointer               :: LastProperty        
        type (T_PorousMediaProperties), pointer     :: Next             => null() !Lúcia

        real,    pointer, dimension(:,:,:)          :: PropI
        real,    pointer, dimension(:,:,:)          :: PropII
        real,    pointer, dimension(:,:,:)          :: PropInew 
        logical                                     :: PorousMediaProperties
        real,    pointer, dimension(:,:,:)          :: Volume   
        integer                                     :: PropertiesNumber    = 0

        
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
                                              HorizontalMapID,                            &
                                              TopographyID,                               &
                                              BasinGeometryID,                            &
                                              DrainageNetworkID,                          &
                                              PorousMediaID,                              &
                                              GeometryID,                                 &
                                              MapID,                                      &
                                              STAT)
     
        !Arguments---------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID 
        integer, optional, intent(OUT)                  :: STAT 
        integer                                         :: ComputeTimeID
        integer                                         :: HorizontalGridID
        integer                                         :: HorizontalMapID
        integer                                         :: TopographyID
        integer                                         :: BasinGeometryID
        integer                                         :: DrainageNetworkID
        integer                                         :: PorousMediaID
        integer                                         :: GeometryID
        integer                                         :: MapID
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
            Me%ObjTopography     = AssociateInstance (mGRIDDATA_,       TopographyID    ) 
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjPorousMedia    = AssociateInstance (mPOROUSMEDIA_,    PorousMediaID   )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMap_,            MapID           )
        
            if (DrainageNetworkID /= 0)                                                     &
                Me%ObjDrainageNetwork  = AssociateInstance (MDRAINAGENETWORK_,  DrainageNetworkID )
            
            !Geometry Size
            call GetGeometrySize    (Me%ObjGeometry,             &    
                                     Size     = Me%Size,         &
                                     WorkSize = Me%WorkSize,     &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR05'

            Me%Size2D%ILB = Me%Size%ILB
            Me%Size2D%IUB = Me%Size%IUB
            Me%Size2D%JLB = Me%Size%JLB
            Me%Size2D%JUB = Me%Size%JUB

       
            call GetWaterPoints3D(Me%ObjMap,                                                &
                                  Me%ExtVar%WaterPoints3D,                                  &
                                  STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                                      &
                stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR06'

            call GetComputeCurrentTime(Me%ObjTime, Me%ExtVar%Now, STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR07'

            call GetComputeTimeLimits(Me%ObjTime,                      &
                                      EndTime   = Me%ExtVar%EndTime,   &
                                      BeginTime = Me%ExtVar%BeginTime, &
                                      STAT      = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)    &
                    stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR08'
            

            call ReadFileNames

            !Constructs the DataFile
            call ConstructEnterData (Me%ObjEnterData, Me%Files%DataFile, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR09'                
           
            call ReadGlobalOptions

            call AllocateVariables

            call Construct_PropertyList
        
    !       call ConstructHDF5Output    em teste no modifier e no kill
    
            call ConstructTimeSerie

    !       call ConstructProfileOutput   em teste


            call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMedia - ModulePorousMedia - ERR02'

            call UnGetMap (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructPorousMediaProperties - ModulePorousMediaProperties - ERR10'

            !Returns ID
            ObjPorousMediaPropertiesID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModulePorousMediaProperties - ConstructPorousMediaProperties - ERR011' 

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
        call ReadFileName ('POROUSPROP_HDF', Me%Files%TransientHDF, "PorousMedia HDF File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01b'
                
        !Reads the name of the file where to store final data
        call ReadFileName ('POROUSPROPERTIES_FIN', Me%Files%FinalFile, "PorousMedia Final File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadDataFile - ModulePorousMediaProperties - ERR01c'
   


    !   Reads the name of the sedimentQuality file 


    end subroutine ReadFileNames
    
    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------        
        integer                                     :: STAT_CALL
        integer                                     :: iflag


        call GetData(Me%AdvDiff%SpatialMethod,                                            &   !Lúcia
                     Me%ObjEnterData, iflag,                                            &
                     SearchType = FromFile,                                             &
                     keyword    = 'SPATIAL_METHOD',                                       &
                     Default    = 0,                                               &                                           
                     ClientModule ='ModulePorousMediaProperties',                                 &
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
        allocate (ME%Volume                  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (ME%AdvDiff%DifusionNumber  (ILB:IUB,JLB:JUB,KLB:KUB))
        allocate (ME%AdvDiff%ReynoldsMNumber (ILB:IUB,JLB:JUB,KLB:KUB))

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
                        stop 'Construct_PropertyList - ModuleSoilProperties - ERR01'
                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_PropertyList - ModuleSoilProperties - ERR02'
            else cd1
                stop       'Construct_PropertyList - ModuleSoilProperties - ERR03'
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
        if(STAT_CALL .NE. SUCCESS_)stop 'Construct_Property - ModuleSoilProperties - ERR00'
        
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
                     ClientModule = 'ModuleSoilProperties',                          &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'Construct_PropertyState - ModuleSoilProperties - ERR01'
        if(iflag == 0)              stop 'Construct_PropertyState - ModuleSoilProperties - ERR02'

        if (NewProperty%Particulate)then
            if(.not. Check_Particulate_Property(NewProperty%ID%IDNumber)) then 
                write(*,*) 'Property '//trim(NewProperty%ID%Name)// 'is not'
                write(*,*) 'recognised as PARTICULATE'
                stop 'Construct_PropertyState - ModuleSoilProperties - ERR03'
            end if
        endif

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
                     ClientModule = 'ModuleSoilProperties',                          &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_)                                                      &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR01'

        if(NewProperty%Evolution%AdvectionDiffusion)NewProperty%Evolution%Variable = .true.
        
        call GetData(NewProperty%Evolution%Partitioning,                                 &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'PARTITION',                                         &
                     ClientModule = 'ModuleSoilProperties',                          &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR03'

        if (NewProperty%Evolution%Partitioning) NewProperty%Evolution%Variable = .true.

        !Partition parameters  
        if (NewProperty%Evolution%Partitioning)                                          &
  !          call Read_Partition_Parameters (NewProperty, ClientNumber)

        !<BeginKeyword>
            !Keyword          : CATION_EXCHANGE
            !<BeginDescription>       
               ! Property has cation exchange as sink and source
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%Evolution%CationExchangeProcess,                        &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'CATION_EXCHANGE',                                   &
                     ClientModule = 'ModuleSoilProperties',                              &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR03'

        if (NewProperty%Evolution%CationExchangeProcess)NewProperty%Evolution%Variable = .true.

        !<BeginKeyword>
            !Keyword          : CHEMICAL_EQUILIBRIUM
            !<BeginDescription>       
               ! Property has chemical equilibrium as sink and source
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%Evolution%ChemEquilibriumProcess,                        &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'CHEMICAL_EQUILIBRIUM',                                   &
                     ClientModule = 'ModuleSoilProperties',                              &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR03'

        if (NewProperty%Evolution%ChemEquilibriumProcess)NewProperty%Evolution%Variable = .true.
        !<BeginKeyword>
            !Keyword          : SOIL_QUALITY
            !<BeginDescription>       
               ! Property has the Soil quality model as sink and source
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%SoilQuality,                              &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SOIL_QUALITY',                                  &
                     ClientModule = 'ModuleSoilProperties',                          &
                     default      = OFF,                                                 &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR04'

        if (NewProperty%Evolution%SoilQuality) NewProperty%Evolution%Variable = .true.
        
        !<BeginKeyword>
            !Keyword          : Soil_WATER_FLUXES
            !<BeginDescription>       
               !  This property has fluxes at the Soil water interface? no - 0;  yes - 1
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%SoilWaterFluxes,                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'SOIL_WATER_FLUXES',                             &
                     ClientModule = 'ModuleSoilProperties',                          &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR05'
        
        if (NewProperty%Evolution%SoilWaterFluxes)NewProperty%Evolution%Variable = .true.

        !<BeginKeyword>
            !Keyword          : MACROPORES
            !<BeginDescription>       
               !  This property has fluxes with macropores? no - 0;  yes - 1
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.), 0 (.false.)
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%Evolution%Macropores,                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromBlock,                                           &
                     keyword      = 'MACROPORES',                             &
                     ClientModule = 'ModuleSoilProperties',                          &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR05'
        
            ! Quanto muito, constroi-se uma função que leia os parametros dentyro das mesmas características da propriedade 
        
 !       if (NewProperty%Evolution%AdvectionDiffusion)then   

  !          call Read_Advec_Difus_Parameters    (NewProperty)

   !         call Construct_Property_Diffusivity (NewProperty)
        
    !    end if

        !Property time step
        if (NewProperty%Evolution%Variable) then

            ModelDT = Me%ExtVar%DT

            call GetData(NewProperty%Evolution%DTInterval,                               &
                         Me%ObjEnterData, iflag,                                         &
                         SearchType   = FromBlock,                                       &
                         keyword      = 'DTINTERVAL',                                    &
                         Default      = ModelDT,                                         &
                         ClientModule = 'ModuleSoilProperties',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                                 &
                stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR08'
                                       
            
            if (NewProperty%Evolution%DTInterval < ModelDT) then
                write(*,*) 
                write(*,*) 'Property time step is smaller then model time step'
                stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR09'

            elseif (NewProperty%Evolution%DTInterval > ModelDT) then 

                !Property time step must be a multiple of the model time step
                auxFactor = NewProperty%Evolution%DTInterval  / ModelDT

                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) 'Property time step must be a multiple of model time step.'
                    write(*,*) 'Please review your input data.'
                    stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR10'
                endif

                !Run period in seconds
                DTaux = Me%ExtVar%EndTime - Me%ExtVar%Now

                !The run period   must be a multiple of the Property DT
                auxFactor = DTaux / NewProperty%Evolution%DTInterval

                ErrorAux = auxFactor - int(auxFactor)
                if (ErrorAux /= 0) then

                    write(*,*) 
                    write(*,*) 'Property time step is not a multiple of model time step.'
                    stop 'Construct_PropertyEvolution - ModuleSoilProperties - ERR11'
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
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR03'
        NewProperty%Concentration(:,:,:) = FillValueReal

        allocate(NewProperty%ConcentrationOld(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR03'
        NewProperty%Concentration(:,:,:) = FillValueReal


        allocate(NewProperty%SedPropAdvFlux(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR03'
        NewProperty%SedPropAdvFlux(:,:) = 0.0
        


        call GetData(NewProperty%DifCoef,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DIFFUSION_COEF',                                              &
                     Default      = .1,                                            &                        
                     ClientModule = 'ModuleSoilProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR05'

        call GetData(NewProperty%Dispersivity,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DISPERSIVITY',                                              &
                     Default      = .1,                                            &                        
                     ClientModule = 'ModuleSoilProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR05'

        call GetData(NewProperty%rainconc,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'RAIN_CONC',                                              &
                     Default      = .1,                                            &                        
                     ClientModule = 'ModuleSoilProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR05'
        
        call GetData(NewProperty%bottomconc,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BOTTOM_CONC',                                              &
                     Default      = .1,                                            &                        
                     ClientModule = 'ModuleSoilProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR05'


        call GetData(NewProperty%MinValue,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MIN_VALUE',                                        &
                     ClientModule = 'ModuleSoilProperties',                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_PropertyValues - ModuleSoilProperties - ERR06'
        if (iflag==1)  then
            NewProperty%Evolution%MinConcentration = ON
        else
            NewProperty%Evolution%MinConcentration = OFF
        endif

        if(NewProperty%Evolution%MinConcentration)then
            allocate(NewProperty%Mass_Created(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)&
                stop 'Construct_PropertyValues - ModuleSoilProperties - ERR07'
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
                                       PointsToFill3D       = Me%ExtVar%WaterPoints3D,     &
                                       Matrix3D             = NewProperty%Concentration,        &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                          &
                stop 'Construct_PropertyValues - ModuleSoilProperties - ERR07'

            if(.not. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_PropertyValues - ModuleSoilProperties - ERR08'
            end if


!Lúcia
            call SetMatrixValue(NewProperty%ConcentrationOld, Me%Size, NewProperty%Concentration,Me%ExtVar%WaterPoints3D)



            call CheckFieldConsistence (NewProperty)

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

        endDo

        !----------------------------------------------------------------------

            call GetData(TimeSerieLocationFile,                                             &
                         Me%ObjEnterData,iflag,                                             &
                         SearchType   = FromFile,                                           &
                         keyword      = 'TIME_SERIE_LOCATION',                              &
                         ClientModule = 'ModulePorousMedia',                                &
                         Default      = Me%Files%DataFile,                                  &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR02' 
            
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
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProfileOutput - ModulePorousMedia - ERR03' 


        deallocate (PropertyList)
        
    end subroutine ConstructProfileOutput

    !---------------------------------------------------------------------------

    subroutine Construct_PropertyOutPut(NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_Property),    pointer        :: NewProperty

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        !<BeginKeyword>
            !Keyword          : TIME_SERIE
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a time serie for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>
        call GetData(NewProperty%TimeSerie,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'TIME_SERIE',                                        &
                     ClientModule = 'ModuleSoilProperties',                          &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSoilProperties - ERR01'
        
        !<BeginKeyword>
            !Keyword          : BOX_TIME_SERIE
            !<BeginDescription>       
                ! Checks out if the user pretends to write a time serie inside each box for this property
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.) , 0 (.false.) 
            !Search Type      : FromBlock
            !Begin Block     : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%BoxTimeSerie,                                           &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'BOX_TIME_SERIE',                                    &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModuleSoilProperties',                          &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSoilProperties - ERR02'

        !<BeginKeyword>
            !Keyword          : BOX_TIME_SERIE2D
            !<BeginDescription>       
                ! Checks out if the user pretends to write a time serie for the cumulative 
                ! flux of this property in the boundary
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : 1 (.true.) , 0 (.false.) 
            !Search Type      : FromBlock
            !Begin Block     : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%BoxTimeSerie2D,                                           &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'BOX_TIME_SERIE2D',                                    &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     ClientModule = 'ModuleSoilProperties',                          &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSoilProperties - ERR02'

        !<BeginKeyword>
            !Keyword          : OUTPUT_HDF
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a outputs in HDF for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SEDPROP
            !Multiple Options : Do not have
            !Search Type      : FromBlock
            !Begin Block      : <beginproperty>
            !End Block        : <endproperty>
        !<EndKeyword>

        call GetData(NewProperty%OutputHDF,                                              &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword      = 'OUTPUT_HDF',                                        &
                     ClientModule = 'ModuleSoilProperties',                          &
                     Default      = .false.,                                             &
                     SearchType   = FromBlock,                                           &
                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'Construct_PropertyOutPut - ModuleSoilProperties - ERR03'
        
    end subroutine Construct_PropertyOutPut
   
   !---------------------------------------------------------------------------
   
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
                stop 'ReadOldConcBoundariesHDF - ModuleSoilProperties - ERR01'


            PropertyName = trim(adjustl(NewProperty%ID%name))

            NewProperty%Concentration(:,:,:) = FillValueReal


            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                              &
                                 WorkJLB, WorkJUB, WorkKLB, WorkKUB,                     &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleSoilProperties - ERR02'

            call HDF5ReadData   (ObjHDF5, "/Concentration/"//NewProperty%ID%Name,        &
                                 NewProperty%ID%Name,                                    &
                                 Array3D = NewProperty%Concentration,                    &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleSoilProperties - ERR03'


            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadOldConcBoundariesHDF - ModuleSoilProperties - ERR06'

        else
            
            write(*,*)
            stop 'ReadOldConcBoundariesHDF - ModuleSoilProperties - ERR07'

        end if cd0

    end subroutine ReadOldConcBoundariesHDF


    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        type (T_Property), pointer              :: PropertyX

        !----------------------------------------------------------------------

        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))
  !          if (PropertyX%TimeSerie) then

                call WriteTimeSerie(Me%ObjTimeSerie,                    &
                                    Data3D = PropertyX%Concentration,   &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                              &
                    stop 'OutPut_TimeSeries - ModuleSoilProperties - ERR01'

  !          endif
            PropertyX=>PropertyX%Next
        enddo

    end subroutine OutPut_TimeSeries




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

        
        nProperties = Me%PropertiesNumber 

        !Allocates PropertyList
        allocate(PropertyList(nProperties))

        
        n=1
        PropertyX  => Me%FirstProperty

        do while (associated(PropertyX))

        PropertyList(n)  = trim(PropertyX%ID%Name)
        n=n+1

        PropertyX=>PropertyX%Next
        End Do


        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModulePorousMedia',                                &
                     Default      = Me%Files%DataFile,                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR01' 

        if (iflag == 1) then
            Me%OutPut%TimeSerieON = .true.
        else
            Me%OutPut%TimeSerieON = .false.
        endif

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "srp",                                        &
                            WaterPoints3D = Me%ExtVar%WaterPoints3D,                    &
                            STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructTimeSerie - PorousMedia - ERR02' 

        !Deallocates PropertyList
        deallocate(PropertyList)
       
    end subroutine ConstructTimeSerie

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

                                stop 'Subroutine CheckFieldConsistence; ModuleSoilProperties. ERR01.'

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
                stop 'CheckFieldConsistence - ModuleSoilProperties - ERR01' 

            open(UnitAux, FILE = trim(NewProperty%ID%name)//'.new',                 &
                 FORM = 'FORMATTED', STATUS = 'UNKNOWN', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'CheckFieldConsistence - ModuleSoilProperties - ERR02' 

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
                stop 'CheckFieldConsistence - ModuleSoilProperties - ERR03' 

            write(*,*) 'A new concentration file was created for property: ', trim(NewProperty%ID%Name)
            write(*,*) 'Run again with this new file ', trim(NewProperty%ID%name)//'.new'
            stop 'CheckFieldConsistence - ModuleSoilProperties - ERR04'  

        endif

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------


    
     !-------------------------------------------------------------------------

     subroutine ConstructHDF5Output        

        !Local-----------------------------------------------------------------
        integer                                             :: ILB,IUB,JLB,JUB,KLB,KUB    
        integer                                             :: STAT_CALL
        integer                                             :: HDF5_CREATE
        real, dimension(:, :), pointer                      :: BottomData

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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR01'

        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits      (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR03'

        !Writes the Grid
        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Topography", "m",                    &
                              Array2D = Me%ExtVar%Topography, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR04'

        call GetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR05'

        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR06'

        call HDF5WriteData      (Me%ObjHDF5, "/Grid", "BasinPoints", "-",                   &
                              Array2D = Me%ExtVar%BasinPoints, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR07'

        !Water Points
        call HDF5WriteData   ( Me%ObjHDF5,  "/Grid", "WaterPoints3D", "-",                  &
                               Array3D = Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR08'
              
        call GetBasinPoints   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR09'  

                
        !Flushes All pending HDF5 commands
        call HDF5FlushMemory    (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR10'

        call UnGetGridData(Me%ObjBottomTopography, BottomData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5Output - ModulePorousMedia - ERR11'

        call UnGetBasin   (Me%ObjBasinGeometry, Me%ExtVar%BasinPoints, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR12'  


    end subroutine ConstructHDF5Output

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     subroutine GetPropInfiltration (ObjPorousMediaPropertiesID, Infiltration, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(8), dimension(:, :), pointer               :: Infiltration
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaPropertiesID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMediaProperties_, Me%InstanceID)
            
            Infiltration => Me%ExtVar%Infiltration

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPropInfiltration

    !--------------------------------------------------------------------------

    subroutine GetPropEfectiveEVTP (ObjPorousMediaPropertiesID, EfectiveEVTP, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(8), dimension(:, :), pointer               :: EfectiveEVTP
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        
        call Ready(ObjPorousMediaPropertiesID, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMediaProperties_, Me%InstanceID)
            
            EfectiveEVTP => Me%ExtVar%EfectiveEVTP

!            EfectiveEVTP2 => Me%ExtVar%EfectiveEVTP2

!           PlantWaterStress => Me%PlantWaterStress

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPropEfectiveEVTP

    !--------------------------------------------------------------------------
   
    
    subroutine GetNextPorousMediaPropDT (ObjPorousMediaID, PorousMediaPropDT, STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: ObjPorousMediaID
        real, intent(OUT)                               :: PorousMEdiaPropDT
        integer, intent(OUT), optional                  :: STAT

        !Local------------------------------------------------------------------
        integer                                         :: STAT_CALL, ready_

        !-----------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ObjPorousMediaID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            PorousMediaPropDT        = Me%ExtVar%PorousMediaPropDT

            STAT_CALL = SUCCESS_
        else 
            STAT_CALL = ready_
        end if

        if (present(STAT)) STAT = STAT_CALL

    end subroutine GetNextPorousMediapropDT
    
    
    !--------------------------------------------------------------------------
    subroutine GetPorousMediaPropertiesPointer (ObjPorousMediaPropertiesID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPorousMediaProperties_, Me%InstanceID)

            Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPorousMediaPropertiesPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetPorousMediaPropertiesInteger (ObjPorousMediaPropertiesID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetPorousMediaPropertiesInteger

    !--------------------------------------------------------------------------

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

    subroutine UnGetPorousMediaProperties3D_R8i(ObjPorousMediaPropertiesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjPorousMediaPropertiesID
        real(8), dimension(:, :), pointer               :: Array
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

    end subroutine UnGetPorousMediaProperties3D_R8i


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
        type (T_Property), pointer                  :: PropertyX


        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjPorousMediaPropertiesID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeTimeStep     (Me%ObjTime, Me%ExtVar%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMedia - ERR02'

            call ReadLockExternalVar

            PropertyX => Me%FirstProperty

            do while (associated(PropertyX))

                if (PropertyX%Evolution%AdvectionDiffusion) then

                     call CalculateAdvectionDiffusion(PropertyX)

                endif

                PropertyX => PropertyX%Next

            enddo


            call OutPut_TimeSeries

    !       call ProfileOutput    em teste no construct e no kill
        
            call ReadUnlockExternalVar


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyPorousMediaProperties

    !-----------------------------------------------------------------------------
    
    subroutine CalculateAdvectionDiffusion (PropertyX)
    
        !Arguments-------------------------------------------------------------
        type (T_Property), pointer                  :: PropertyX


        !Local-----------------------------------------------------------------
        type (T_Property), pointer                  :: CurrProperty
        integer                                     :: i, j, k, CHUNK
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

                Me%Volume(i,j,k)= Theta(i,j,k)*ME%extvar%Cellvolume(i,j,k)
            

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



    end subroutine CalculateAdvectionDiffusion
    
    !---------------------------------------------------------------------------
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
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR05'
                endif

                !Deassociates External Instances
                if (Me%ObjDrainageNetwork /= 0) then
                    nUsers = DeassociateInstance (mDRAINAGENETWORK_, Me%ObjDrainageNetwork)
                    if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR06'
                endif                

       !         call KillHDF5 (Me%ObjHDF5, STAT = STAT_)
        !         if (STAT_ /= SUCCESS_) stop 'KillPorousMedia - Porousmedia - ERR05'  ! em teste no contruct e no modifier

                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR07'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR08'

                nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjTopography)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR09'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR10'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR11'
                
                nUsers = DeassociateInstance (mPOROUSMEDIA_,  Me%ObjPorousMedia)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR12'

                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR13'

                nUsers = DeassociateInstance (mMAP_,  Me%ObjMap)
                if (nUsers == 0) stop 'KillPorousMedia - Porousmedia - ERR14'


!                call KillPorousMedia (Me%ObjPorousMedia, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'KillBasin - ModuleBasin - ERR15'


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


    subroutine ReadLockExternalVar                

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: STAT
        !Begin-----------------------------------------------------------------

        call GetOldWaterContent (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR010'

        call GetWaterContent    (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR020'

        call GetFluxU           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR030'

        call GetFluxV           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR040'

        call GetFluxWOld        (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR050'

        call GetUnsatWOld       (Me%ObjPorousMedia, Me%ExtVar%UnsatW, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'

        call GetGridCellArea    (Me%ObjHorizontalGrid,                                     & 
                                 GridCellArea = Me%ExtVar%Area,                            & 
                                 STAT = STAT_CALL )    
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR070'

        call GetWaterPoints3D   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR080'

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
        integer                                     :: STAT
        !Begin-----------------------------------------------------------------

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContentOld, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR010'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%WaterContent, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR020'
        
        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxU, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR030'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxV, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR040'

        call UnGetPorousMedia           (Me%ObjPorousMedia, Me%ExtVar%FluxW, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR050'

        call UnGetPorousMedia           (Me%ObjPorousMedia,Me%ExtVar%UnsatW, STAT)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR060'
        
        call UnGetHorizontalGrid        (Me%ObjHorizontalGrid,Me%ExtVar%Area,STAT)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR070'

        call UnGetMap                   (Me%ObjMap, Me%ExtVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModulePorousMediaProperties - ERR080'
        
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








