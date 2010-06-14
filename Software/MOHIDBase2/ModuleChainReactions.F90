!----------------------------------------------------------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!----------------------------------------------------------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ChainReactions
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2010
! REVISION      : Eduardo Jauch - v4.0
! DESCRIPTION   : Zero-dimensional model for simple chain reactions
!                 
!----------------------------------------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------------------------------------

Module ModuleChainReactions

    use ModuleFunctions            
    use ModuleEnterData
    use ModuleGlobalData 
    use ModuleGeometry  
    use ModuleFillMatrix
    use ModuleBasinGeometry   
    use ModuleMap       

    implicit none

    private 

    !Subroutines-------------------------------------------------------------------------------------------------------------------
    !Constructor
    public  :: StartChainReactions
    private ::      AllocateInstance
    private ::      InitializeInstance
!    private ::      ReadFileNames
    private ::      ReadInputFile
    private ::          ReadChainReactionsOptions
    private ::              ConstructSinkList
    private ::                  ConstructSinkProperty
    private ::                      ConstructSinkOptions
    private ::                      ConstructSinkValues
    private ::                          ConstructSinkValues2D
    private ::                          ConstructSinkValues3D
    private ::                  AddSinkProperty 
    private ::                  LinkProducts
    private ::                      ConstructProduct
    private ::              CheckSinkProperties 
    private ::              ConstructPartitionList
    private ::                  ConstructPartitionProperty
    private ::                      ConstructPartitionOptions
    private ::                      ConstructPartitionValues
    private ::                          ConstructPartitionValues2D
    private ::                          ConstructPartitionValues3D    
    private ::                  AddPartition
    private ::      ResetSinkProperty
    private ::      CreatePropertiesList  
    private ::      SearchSinkProperty                      
                
    !Selector
    public  :: GetCRPropertiesList
    public  :: SetCRPropertyConcentration
    private ::      SetPropertyConcentration2D
    private ::      SetPropertyConcentration3D
    public  :: UngetChainReactions           
                        
    !Modifier
    public  :: ModifyChainReactions 
    private ::      ModifyChainReactions2D
    private ::      ModifyChainReactions3D
    
    private ::          ActualizePropertiesFromFile
    private ::          GetSinkrate
    private ::          GetPartitionCoef
        
    !Destructor
    public  :: KillChainReactions                                                     
    private ::      DeallocateInstance

    !Management
    private :: Ready
    private ::      LocateObjChainReactions
        
    !Interfaces
    interface SetCRPropertyConcentration
        module procedure SetPropertyConcentration2D
        module procedure SetPropertyConcentration3D
    end interface SetCRPropertyConcentration
    
    interface ModifyChainReactions
        module procedure ModifyChainReactions2D
        module procedure ModifyChainReactions3D
    end interface ModifyChainReactions
        
    !Constants---------------------------------------------------------------------------------------------------------------------
    integer,                       parameter :: MaxStrLength       = StringLength      !StringLength is defined in ModuleGlobalData     
    
    integer,                       parameter :: UnknownGeometry    = 0
    integer,                       parameter :: Geometry2D         = 1
    integer,                       parameter :: Geometry3D         = 2
    
    integer,                       parameter :: NoReaction         = 0
    integer,                       parameter :: ZeroOrder          = 1
    integer,                       parameter :: FirstOrder         = 2
    
    integer,                       parameter :: RateUnknown        = 0
    integer,                       parameter :: RateDefault        = 1
    integer,                       parameter :: RateAlternative    = 2
    
    integer,                       parameter :: NoSink             = 0
    integer,                       parameter :: ConstantSink       = 1 !Single value
    integer,                       parameter :: VariableSink       = 2
    
    integer,                       parameter :: NoPartition        = 0
    integer,                       parameter :: ConstantPartition  = 1 !Single value
    integer,                       parameter :: VariablePartition  = 2
    
    integer,                       parameter :: PartCoefSimple     = 1
    integer,                       parameter :: PartCoefAdsorption = 2 
    
    integer,                       parameter :: WaterPhase         = 1
    integer,                       parameter :: SolidPhase         = 2
    
    character(LEN = StringLength), parameter :: prop_block_begin   = '<beginproperty>'
    character(LEN = StringLength), parameter :: prop_block_end     = '<endproperty>'
       
    character(LEN = StringLength), parameter :: sink_block_begin   = '<beginsink>'
    character(LEN = StringLength), parameter :: sink_block_end     = '<endsink>'

    character(LEN = StringLength), parameter :: part_block_begin   = '<beginpartition>'
    character(LEN = StringLength), parameter :: part_block_end     = '<endpartition>'

    !Types-------------------------------------------------------------------------------------------------------------------------            
       
    type T_Property2D
        real, pointer, dimension(:,:)   :: SinkRate       !ZERO order: mg/s - FIRST order: 1/s
        real, pointer, dimension(:,:)   :: PartitionCoef  
        real, pointer, dimension(:,:)   :: Concentration  !In mg/L
    end type T_Property2D
    
    type T_Property3D
        real, pointer, dimension(:,:,:) :: SinkRate       !ZERO order: mg/s - FIRST order: 1/s
        real, pointer, dimension(:,:,:) :: PartitionCoef  
        real, pointer, dimension(:,:,:) :: Concentration  !In mg/L
    end type T_Property3D
    
    type T_Property
        type(T_PropertyID)        :: ID
        integer                   :: ProductID 
        character(MaxStrLength)   :: ProductName
        integer                   :: RateOrder           
        integer                   :: RateMethod
        real                      :: NewMass
        real                      :: Sink
        real                      :: Source
        real                      :: Mass
        logical                   :: IsOnlyProduct
        integer                   :: SinkEvolution         !0 - no Sink, 1 - constant, 2 - variable
        real                      :: SinkRate              !SinkRate if SinkEvolution is constant
        integer                   :: Phase
        logical                   :: NewConcCalculated    
        type(T_Property2D)        :: G2D
        type(T_Property3D)        :: G3D        
        type(T_Property), pointer :: Next 
        type(T_Property), pointer :: Prev         
        type(T_Property), pointer :: Product        
    end type T_Property
    
    type T_Partition
        type(T_PropertyID)             :: ID
        type(T_Property), pointer      :: Pair1
        type(T_Property), pointer      :: Pair2         
        real                           :: PartitionCoef
        integer                        :: PartitionEvolution
        type(T_Property2D)             :: G2D
        type(T_Property3D)             :: G3D         
        type(T_Partition), pointer :: Next
        type(T_Partition), pointer :: Prev  
    end type T_Partition
    
    type T_External
        integer, pointer, dimension(:,:,:) :: WaterPoints3D
        integer, pointer, dimension(:,:)   :: WaterPoints2D
        real,    pointer, dimension(:,:,:) :: Theta3D
        real,    pointer, dimension(:,:)   :: Theta2D
        real,    pointer, dimension(:,:,:) :: SoilDensity3D
        real,    pointer, dimension(:,:)   :: SoilDensity2D
        integer, pointer, dimension(:)     :: Properties
    end type T_External
          
    type T_ChainReactionsOptions
        integer :: GeometryType  !0 -> Unknown; 1 -> 2D; 2 -> 3D
        logical :: Compensate
        integer :: RateMethod
        integer :: RateOrder
    end type 
       
    type T_Files
        character(PathLength) :: DataFile
    end type T_Files           
       
    type T_ChainReactions
        private
        integer                             :: InstanceID                    !ID of the ModuleChainReactions instance 
        type(T_ChainReactions), pointer     :: Next                          !Collection of instances of ModuleChainReactions
        
        character(MaxStrLength)             :: CallingModule

        integer                             :: ObjEnterData         = 0
        integer                             :: ObjTime              = 0
        integer                             :: ObjGeometry          = 0
        integer                             :: ObjHorizontalGrid    = 0
        integer                             :: ObjHorizontalMap     = 0
        integer                             :: ObjBasinGeometry     = 0
        integer                             :: ObjMap               = 0

        type(T_ChainReactionsOptions)       :: Options
        type(T_Files)                       :: Files
                
        type(T_External)                    :: Ext                           !Pointers to Water Mass, Properties Values 
                                                                             !and other required data 
        
        type(T_Property), pointer           :: FirstProperty               
        type(T_Property), pointer           :: LastProperty
        type(T_Partition), pointer          :: FirstPartition
        type(T_Partition), pointer          :: LastPartition
        integer, pointer, dimension(:)      :: PropertiesList                !List of properties ID's
        integer                             :: PropertiesCount               !Number of properties
        
        integer                             :: ClientNumber
        
        type (T_Size3D)                     :: Size
        type (T_Size3D)                     :: WorkSize      
        
    end type T_ChainReactions


    !Global Module Variables-------------------------------------------------------------------------------------------------------
    type (T_ChainReactions), pointer :: FirstObjChainReactions
    type (T_ChainReactions), pointer :: Me      
       
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartChainReactions(ChainReactionsID,    & 
                                   Properties,          &  
                                   TimeID,              &
                                   HorizontalGridID,    &
                                   HorizontalMapID,     &
                                   BasinGeometryID,     &
                                   GeometryID,          &
                                   MapID,               &                                    
                                   GeometryType,        &
                                   CallingModule,       &
                                   FileName,            &
                                   PropertiesCount,     &
                                   STAT)
        
        !Arguments-----------------------------------------------------------------------------------------------------------------        
        integer                        :: ChainReactionsID
        integer, dimension(:), pointer :: Properties
        integer                        :: TimeID
        integer                        :: HorizontalGridID
        integer                        :: HorizontalMapID
        integer                        :: BasinGeometryID
        integer                        :: GeometryID
        integer                        :: MapID
        integer                        :: GeometryType
        character(LEN=*),  intent(IN)  :: CallingModule
        character(LEN=*),  intent(IN)  :: FileName
        integer,           intent(OUT) :: PropertiesCount
        integer, optional, intent(OUT) :: STAT     

        !Local---------------------------------------------------------------------------------------------------------------------
        integer            :: STAT_CALL, STAT_
        integer            :: ready_         

        !--------------------------------------------------------------------------------------------------------------------------
        
        write(*,*)
        write(*,*)'------------------CHAIN REACTIONS MODULE------------------'
        write(*,*)'Starting module, wait...'       

        STAT_ = UNKNOWN_        

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mCHAINREACTIONS_)) then
        
            nullify (FirstObjChainReactions)
            call RegisterModule (mCHAINREACTIONS_) 
            
        endif
        
        call Ready(ChainReactionsID, ready_)

        if (ready_ .EQ. OFF_ERR_) then
            
            !Allocate a new instance of ModuleChainReactions (this module...)
            call AllocateInstance 
            
            !Initialize this instance variables
            call InitializeInstance
            
            Me%CallingModule = CallingModule
            
            Me%Ext%Properties => Properties
                                                     
            !Associate External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID)
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjBasinGeometry  = AssociateInstance (mBASINGEOMETRY_,  BasinGeometryID )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMap_,            MapID           )                   
                    
            !Get geometry size
            call GetGeometrySize    (Me%ObjGeometry,          &    
                                     Size     = Me%Size,      &
                                     WorkSize = Me%WorkSize,  &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'StartChainReactions - ModuleChainReactions - ERR010'

            if (.NOT. ((GeometryType .EQ. Geometry2D) .OR. (GeometryType .EQ. Geometry3D))) then
                write (*,*) 
                write (*,*) 'Unknown geometry passed to ModuleChainReactions '
                write (*,*) 'Must be set to 1 (2D) or 2 (3D) '
                stop 'StartChainReactions - ModuleChainReactions. ERR020'
            endif
        
            Me%Options%GeometryType = GeometryType
            
!            !Read files names
!            call ReadFileNames
            
            !Read ModuleChainReactions Input File  
            call ReadInputFile (FileName)
                       
            !Creates the property list to be passed to calling modules
            call CreatePropertiesList
            
            !Returns ID
            ChainReactionsID = Me%InstanceID            

            !Returns Number of Properties
            PropertiesCount = Me%PropertiesCount

            STAT_ = SUCCESS_
            
        else     
            
            stop 'Subroutine StartChainReactions; ModuleChainReactions. ERR070'

        end if    

        if (present(STAT)) STAT = STAT_

        write(*,*)
        write(*,*)'Done.'
        write(*,*)

        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine StartChainReactions
    !------------------------------------------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine AllocateInstance

        !Local---------------------------------------------------------------------------------------------------------------------
        type (T_ChainReactions), pointer :: NewObjChainReactions
        type (T_ChainReactions), pointer :: PreviousObjChainReactions

        !--------------------------------------------------------------------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewObjChainReactions)
        nullify  (NewObjChainReactions%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjChainReactions)) then
        
            FirstObjChainReactions => NewObjChainReactions
            Me                     => NewObjChainReactions
            
        else
        
            PreviousObjChainReactions => FirstObjChainReactions
            Me                        => FirstObjChainReactions%Next
            
            do while (associated(Me))
            
                PreviousObjChainReactions => Me
                Me                        => Me%Next
                
            enddo
            
            Me                             => NewObjChainReactions
            PreviousObjChainReactions%Next => NewObjChainReactions
            
        endif

        Me%InstanceID = RegisterNewInstance (mCHAINREACTIONS_)
        
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine AllocateInstance    
    !------------------------------------------------------------------------------------------------------------------------------
       
   
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine InitializeInstance
    
        !--------------------------------------------------------------------------------------------------------------------------
        
        Me%PropertiesCount = 0 
        
        nullify(Me%FirstProperty)
        nullify(Me%LastProperty)
        nullify(Me%FirstPartition)
        nullify(Me%LastPartition) 
                     
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine InitializeInstance    
    !------------------------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ReadInputFile (FileName)
     
        !Arguments-----------------------------------------------------------------------------------------------------------------
        character(LEN = *) :: FileName      
        
        !Local---------------------------------------------------------------------------------------------------------------------
        integer            :: STAT_CALL
        
        !Begin---------------------------------------------------------------------------------------------------------------------
        !Associate EnterData
        call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Subroutine ReadInputFile; module ModuleChainReactions. ERR010.'      
                
        call ReadChainReactionsOptions                
        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Subroutine ReadInputFile; ModuleSedimentQuality. ERR020.'
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ReadInputFile
    !------------------------------------------------------------------------------------------------------------------------------
    
   
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ReadChainReactionsOptions

        !Local-----------------------------------------------------------------
        integer                 :: FromFile
        integer                 :: STAT_CALL
        integer                 :: flag
        character(MaxStrLength) :: Text
        
        !Begin-----------------------------------------------------------------               
        call GetExtractType (FromFile = FromFile)

        call GetData(Me%Options%RateOrder,                  &
                     Me%ObjEnterData, flag,                 &
                     SearchType   = FromFile,               &
                     keyword      = 'RATE_ORDER',           &
                     default      = NoReaction,             & 
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)
                     
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'Subroutine ReadChainReactionsOptions; Module ModuleChainReactions. ERR010.'                

        if (Me%Options%RateOrder .NE. NoReaction) then
            if (Me%Options%RateOrder .EQ. ZeroOrder) then
                Text = 'ZERO Order Rate type'
            else
                Text = 'FIRST Order Rate type'
            endif
            
            write(*,*)
            write(*,*) 'Warning: All properties will be set to be of ', trim(Text)
        endif

        if (Me%Options%RateOrder .EQ. FirstOrder) then
            
            call GetData(Me%Options%RateMethod,                 &
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'RATE_METHOD',          &
                         default      = RateUnknown,            & 
                         ClientModule = 'ModuleChainReactions', &
                         STAT         = STAT_CALL)
                         
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'Subroutine ReadChainReactionsOptions; Module ModuleChainReactions. ERR010.'                
                
            if (Me%Options%RateMethod .NE. RateUnknown) then
                if (Me%Options%RateMethod .EQ. RateDefault) then
                    Text = 'Default First Order equation'
                else
                    Text = 'Alternative First Order equation'
                endif
                
                write(*,*)
                write(*,*) 'Warning: All properties will be set to use the ', trim(Text)
            endif

        else
        
            Me%Options%RateMethod = RateUnknown    
        
        endif
                        
        call ConstructSinkList  
        call LinkProducts
        call CheckSinkProperties 
        
        call ConstructPartitionList
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ReadChainReactionsOptions         
    !------------------------------------------------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine ConstructSinkList 

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                    :: STAT_CALL
        logical                    :: BlockFound
        type (T_Property), pointer :: NewProperty

        !Begin---------------------------------------------------------------------------------------------------------------------
        do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
                                        ClientNumber    = Me%ClientNumber,  &
                                        block_begin     = prop_block_begin, &
                                        block_end       = prop_block_end,   &
                                        BlockFound      = BlockFound,       &
                                        STAT            = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then    

                if (BlockFound) then                                                  
                    
                    !Construct a New Property 
                    Call ConstructSinkProperty(NewProperty)

                    !Add new Property to the Properties List 
                    Call AddSinkProperty(NewProperty)                                        
                    
                else

                    call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'ConstructPropertyList - ModuleChainReactions - ERR010'
                    exit !No more blocks
                
                end if    

            elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'ConstructPropertyList - ModuleChainReactions - ERR020'
            
            else    
                
                stop 'ConstructPropertyList - ModuleChainReactions - ERR030'
            
            end if    
        
        end do            
        !--------------------------------------------------------------------------------------------------------------------------
                
    end subroutine ConstructSinkList
    !------------------------------------------------------------------------------------------------------------------------------   


    !------------------------------------------------------------------------------------------------------------------------------   
    subroutine ConstructSinkProperty(NewProperty)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL

        !Begin---------------------------------------------------------------------------------------------------------------------             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructSinkProperty - ModuleChainReactions - ERR010'
        
        call ResetSinkProperty (NewProperty)                        
        call ConstructPropertyID (NewProperty%ID, Me%ObjEnterData, FromBlock)
        call ConstructSinkOptions (NewProperty)
        call ConstructSinkValues (NewProperty)
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ConstructSinkProperty    
    !------------------------------------------------------------------------------------------------------------------------------    

    
    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructSinkOptions (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL, iflag

        !Begin---------------------------------------------------------------------------------------------------------------------                 
        call GetData(NewProperty%SinkEvolution,             &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'SINK_EVOLUTION',       &
                     Default      = ConstantSink,           &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)
                     
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkOptions - ModuleChainReactions - ERR040'            
        
        if (NewProperty%SinkEvolution .NE. NoSink) then
        
            NewProperty%IsOnlyProduct = .false.
            
            if (Me%Options%RateOrder .NE. NoReaction) then
                NewProperty%RateOrder = Me%Options%RateOrder
            else
                call GetData(NewProperty%RateOrder,                 &   
                             Me%ObjEnterData, iflag,                &
                             SearchType   = FromBlock,              &
                             keyword      = 'RATE_ORDER',           &
                             Default      = ZeroOrder,              &
                             ClientModule = 'ModuleChainReactions', &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructSinkOptions - ModuleChainReactions - ERR010'
            endif
        
            if (Me%Options%RateMethod .NE. RateUnknown) then
                NewProperty%RateMethod = Me%Options%RateMethod
            else
                call GetData(NewProperty%RateMethod,                &   
                             Me%ObjEnterData, iflag,                &
                             SearchType   = FromBlock,              &
                             keyword      = 'RATE_METHOD',          &
                             Default      = RateDefault,            &
                             ClientModule = 'ModuleChainReactions', &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructSinkOptions - ModuleChainReactions - ERR020'
            endif
                       
            call GetData(NewProperty%ProductName,               &   
                         Me%ObjEnterData, iflag,                &
                         SearchType   = FromBlock,              &
                         keyword      = 'PRODUCT_NAME',         &
                         ClientModule = 'ModuleChainReactions', &
                         STAT         = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructSinkOptions - ModuleChainReactions - ERR020'
                    
            if (iFlag .EQ. 0) then
                write(*,*)
                write(*,*) 'You MUST provide a PRODUCT property name '
                write(*,*) '(keyword PRODUCT_NAME) for property ', trim(NewProperty%ID%Name)
                stop 'ConstructSinkOptions - ModuleChainReactions - ERR030'
            endif

            NewProperty%ProductID = GetPropertyIDNumber(NewProperty%ProductName)
        
        else
        
            NewProperty%IsOnlyProduct = .true.
        
        endif
                
        call GetData(NewProperty%Phase,                     &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'PHASE',                &
                     Default      = WaterPhase,             &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)
                     
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkOptions - ModuleChainReactions - ERR050'  
        
        !--------------------------------------------------------------------------------------------------------------------------
        
    end subroutine ConstructSinkOptions
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructSinkValues (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer           :: NewProperty
        
        !Local---------------------------------------------------------------------------------------------------------------------
        logical                             :: BlockFound
        integer                             :: iflag
        integer                             :: STAT_CALL

        !Begin---------------------------------------------------------------------------------------------------------------------             
    
        if (NewProperty%SinkEvolution .EQ. ConstantSink) then

            call GetData(NewProperty%SinkRate,                  &   
                         Me%ObjEnterData, iflag,                &
                         SearchType   = FromBlock,              &
                         keyword      = 'SINK_RATE',            &
                         ClientModule = 'ModuleChainReactions', &
                         STAT         = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructSinkValues - ModuleChainReactions - ERR010'
            
            if (iFlag .EQ. 0) then
                write(*,*)
                write(*,*) 'You MUST provide a value for keyword SINK_RATE for '
                write(*,*) 'property ', trim(NewProperty%ID%Name)
                stop 'ConstructSinkValues - ModuleChainReactions - ERR020'
            endif
                    
        elseif (NewProperty%SinkEvolution .EQ. VariableSink) then
        
            call ExtractBlockFromBlock(Me%ObjEnterData,     & 
                                       Me%ClientNumber,     &
                                       sink_block_begin,    &
                                       sink_block_end,      &
                                       BlockFound,          &
                                       STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructSinkValues - ConstructSinkValues - ERR030'
                
            if (.NOT. BlockFound) then
                write(*,*)
                write(*,*) 'You MUST provide the block for Sink Rate for '
                write(*,*) 'property ', trim(NewProperty%ID%Name)
                write(*,*) 'Block: <beginsink> <endsink>'
                stop 'ConstructSinkValues - ModuleChainReactions - ERR040'
            endif
        
            if (Me%Options%GeometryType .EQ. Geometry2D) then
                call ConstructSinkValues2D (NewProperty)
            else
                call ConstructSinkValues3D (NewProperty)
            endif    
            
            call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructSinkValues - ConstructSinkValues - ERR040'
                    
        endif            
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructSinkValues
    !------------------------------------------------------------------------------------------------------------------------------            



    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructSinkValues2D (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer           :: NewProperty

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: ILB,IUB
        integer                             :: JLB,JUB
        integer                             :: WorkSizeILB, WorkSizeIUB
        integer                             :: WorkSizeJLB, WorkSizeJUB       

        !Begin---------------------------------------------------------------------------------------------------------------------             
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB

        allocate(NewProperty%G2D%SinkRate(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR010'
        NewProperty%G2D%SinkRate(:,:) = 0.0
    
        !Get water points
        call GetBasinPoints (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR020'
    
        call ConstructFillMatrix (PropertyID       = NewProperty%ID,            &
                                  EnterDataID      = Me%ObjEnterData,           &
                                  TimeID           = Me%ObjTime,                &
                                  HorizontalGridID = Me%ObjHorizontalGrid,      &
                                  ExtractType      = FromBlockInBlock,          &
                                  PointsToFill2D   = Me%Ext%WaterPoints2D,      &
                                  Matrix2D         = NewProperty%G2D%SinkRate,  &
                                  TypeZUV          = TypeZ_,                    &
                                  STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR030'

        if(.NOT. NewProperty%ID%SolutionFromFile)then

            call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'ConstructSinkValues2D - ModuleChainReactions - ERR040'
                
        end if  
        
        call UnGetBasin (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR050'                              
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructSinkValues2D
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructSinkValues3D (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: ILB,IUB
        integer                   :: JLB,JUB
        integer                   :: KLB,KUB        
        integer                   :: WorkSizeILB, WorkSizeIUB
        integer                   :: WorkSizeJLB, WorkSizeJUB       
        integer                   :: WorkSizeKLB, WorkSizeKUB

        !Begin---------------------------------------------------------------------------------------------------------------------             
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

        allocate(NewProperty%G3D%SinkRate(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructSinkValues3D - ModuleChainReactions - ERR010'
        NewProperty%G3D%SinkRate(:,:,:) = 0.0
    
    
        !Get water points
        call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkValues3D - ModuleChainReactions - ERR020'
                    
        call ConstructFillMatrix  (PropertyID       = NewProperty%ID,           &
                                   EnterDataID      = Me%ObjEnterData,          &
                                   TimeID           = Me%ObjTime,               &
                                   HorizontalGridID = Me%ObjHorizontalGrid,     &
                                   GeometryID       = Me%ObjGeometry,           &
                                   ExtractType      = FromBlockInBlock,         &
                                   PointsToFill3D   = Me%Ext%WaterPoints3D,     &
                                   Matrix3D         = NewProperty%G3D%SinkRate, &
                                   TypeZUV          = TypeZ_,                   &
                                   STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkValues3D - ModuleChainReactions - ERR030'

        if(.NOT. NewProperty%ID%SolutionFromFile)then

            call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'ConstructSinkValues3D - ModuleChainReactions - ERR040'
                                          
        end if    

        call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructSinkValues3D - ModuleChainReactions - ERR050'                
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructSinkValues3D
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine AddSinkProperty(NewProperty)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: NewProperty

        !Begin---------------------------------------------------------------------------------------------------------------------

        ! Add to the Property List a new property
        if (.NOT. associated(Me%FirstProperty)) then           
            Me%FirstProperty => NewProperty
            Me%LastProperty  => NewProperty
            
            Me%PropertiesCount = 1
        else
            NewProperty%Prev     => Me%LastProperty
            Me%LastProperty%Next => NewProperty
            Me%LastProperty      => NewProperty
            
            Me%PropertiesCount = Me%PropertiesCount + 1
        end if 
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine AddSinkProperty 
    !------------------------------------------------------------------------------------------------------------------------------ 


    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine LinkProducts
    
        !Local---------------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: PropertyX
        type(T_Property), pointer :: ProductX
        type(T_Property), pointer :: NewProperty
        integer                   :: STAT

        !Begin---------------------------------------------------------------------------------------------------------------------

        PropertyX => Me%FirstProperty
        
        do while (associated(PropertyX))
            
            if (.NOT. PropertyX%IsOnlyProduct) then

                call SearchSinkProperty(PropertyX%ProductID, ProductX, STAT)
                
                if (STAT .NE. SUCCESS_) then
                
                    nullify (NewProperty)
                    
                    call ConstructProduct (NewProperty,  PropertyX%ProductName,  PropertyX%ProductID)
                    call AddSinkProperty (NewProperty)  
                                      
                    PropertyX%Product => NewProperty                    
                    
                else
                
                    PropertyX%Product => ProductX
                    
                endif
                
            endif
            
            PropertyX => PropertyX%Next
        
        enddo

        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine LinkProducts
    !------------------------------------------------------------------------------------------------------------------------------ 
    
    
    !------------------------------------------------------------------------------------------------------------------------------   
    subroutine ConstructProduct(NewProperty, Name, IDNumber)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty
        character(LEN=*)          :: Name
        integer                   :: IDNumber

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL

        !Begin---------------------------------------------------------------------------------------------------------------------             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructProduct - ModuleChainReactions - ERR010'
        
        call ResetSinkProperty(NewProperty)
        
        NewProperty%IsOnlyProduct = .true.        
        NewProperty%ID%Name       = Name
        NewProperty%ID%IDNumber   = IDNumber
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ConstructProduct    
    !------------------------------------------------------------------------------------------------------------------------------    

    
    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine CheckSinkProperties 
    
        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: Index
        integer                   :: UB, LB
        type(T_Property), pointer :: PropertyX
        logical                   :: Found

        !Begin---------------------------------------------------------------------------------------------------------------------        
        LB = lbound(Me%Ext%Properties, 1)
        UB = ubound(Me%Ext%Properties, 1)
        
        PropertyX => Me%FirstProperty
        
        do while (associated(PropertyX))
            Found        = .false.
               
            do Index = LB, UB
                if (Me%Ext%Properties(Index) .EQ. PropertyX%ID%IDNumber) then
                    Found = .true.
                    exit
                endif                
            end do
            
            if (.NOT. Found) then 
                write(*,*)
                write(*,*) 'Property ', trim(PropertyX%ID%Name), ' was not found in the ', trim(Me%CallingModule)
                write(*,*) 'list of properties '
                stop 'CheckSinkProperties - ModuleChainReactions - ERR010' 
            endif
            
            PropertyX => PropertyX%Next
        end do
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine CheckSinkProperties
    !------------------------------------------------------------------------------------------------------------------------------ 

    
    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine ConstructPartitionList 

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                     :: STAT_CALL
        logical                     :: BlockFound
        type(T_Partition), pointer  :: NewPartition

        !Begin---------------------------------------------------------------------------------------------------------------------
        do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
                                        ClientNumber    = Me%ClientNumber,  &
                                        block_begin     = part_block_begin, &
                                        block_end       = part_block_end,   &
                                        BlockFound      = BlockFound,       &
                                        STAT            = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then    

                if (BlockFound) then                                                  
                    
                    !Construct a New Property 
                    Call ConstructPartitionProperty(NewPartition)

                    !Add new Property to the Properties List 
                    Call AddPartition(NewPartition)                                        
                    
                else

                    call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'ConstructPartitionList - ModuleChainReactions - ERR010'
                    exit !No more blocks
                
                end if    

            elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'ConstructPartitionList - ModuleChainReactions - ERR020'
            
            else    
                
                stop 'ConstructPartitionList - ModuleChainReactions - ERR030'
            
            end if    
        
        end do            
        !--------------------------------------------------------------------------------------------------------------------------
                
    end subroutine ConstructPartitionList
    !------------------------------------------------------------------------------------------------------------------------------   


    !------------------------------------------------------------------------------------------------------------------------------   
    subroutine ConstructPartitionProperty (NewPartition)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer :: NewPartition

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL

        !Begin---------------------------------------------------------------------------------------------------------------------             
        allocate (NewPartition, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructPartitionProperty - ModuleChainReactions - ERR010'
        
        nullify(NewPartition%Prev)
        nullify(NewPartition%Next)
        nullify(NewPartition%Pair1)
        nullify(NewPartition%Pair2)
       
        NewPartition%PartitionEvolution = NoPartition
        NewPartition%PartitionCoef      = 0.0
        
        call ConstructPartitionID (NewPartition%ID, Me%ObjEnterData, FromBlock)
        call ConstructPartitionOptions (NewPartition)
        call ConstructPartitionValues (NewPartition)
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ConstructPartitionProperty    
    !------------------------------------------------------------------------------------------------------------------------------    


    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine ConstructPartitionID (PartitionID, ObjEnterData, ExtractType)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type (T_PropertyID)                         :: PartitionID
        integer                                     :: ObjEnterData
        integer                                     :: ExtractType

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: STAT_CALL

        !Partition Name
        call GetData(PartitionID%Name, ObjEnterData, flag,  &
                     SearchType   = ExtractType,            &
                     keyword      = 'NAME',                 &
                     Default      = '***',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionID - ModuleChainReactions - ERR010'

        !Units
        call GetData(PartitionID%Units, ObjEnterData, flag, &
                     SearchType   = ExtractType,            &
                     keyword      = 'UNITS',                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionID - ModuleChainReactions - ERR030'

        !Description
        call GetData(PartitionID%Description, ObjEnterData, flag,   &
                     SearchType   = ExtractType,                    &
                     keyword      = 'DESCRIPTION',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPartitionID - ModuleChainReactions - ERR040'
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ConstructPartitionID
    !------------------------------------------------------------------------------------------------------------------------------    

    
    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPartitionOptions (NewPartition)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer :: NewPartition

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL, iflag
        character(MaxStrLength)   :: PairName
        integer                   :: Phase

        !Begin---------------------------------------------------------------------------------------------------------------------                             
        call GetData(NewPartition%PartitionEvolution,       &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'EVOLUTION',            &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)                     
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR010'              
        if (iFlag .EQ. 0) then
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR020'  
        endif
            
        call GetData(PairName,                              &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'PAIR1',                &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)                   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR030'              
        if (iFlag .EQ. 0) then
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR040'  
        endif
        
        call GetData(Phase,                                 &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'PAIR1_PHASE',          &
                     Default      = WaterPhase,             &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)                   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR050'              
        
        call ConstructPairProperty(NewPartition%Pair1, PairName, Phase, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR060'  
            
        call GetData(PairName,                              &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'PAIR2',                &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)                   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR070'              
        if (iFlag .EQ. 0) then
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR080'  
        endif
        
        call GetData(Phase,                                 &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'PAIR2_PHASE',          &
                     Default      = WaterPhase,             &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)                   
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR090'              
        
        call ConstructPairProperty(NewPartition%Pair2, PairName, Phase, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR100'   
        !--------------------------------------------------------------------------------------------------------------------------
        
    end subroutine ConstructPartitionOptions
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPairProperty (PropertyX, PropertyName, Phase, STAT)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: PropertyX
        character(MaxStrLength)   :: PropertyName
        integer                   :: Phase
        integer, optional         :: STAT
        
        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL, STAT_
        integer                   :: Index
        integer                   :: UB, LB
        logical                   :: Found
        integer                   :: PropertyID
        type(T_property), pointer :: NewProperty
        
        !Begin---------------------------------------------------------------------------------------------------------------------        
        STAT_ = UNKNOWN_
        
        LB = lbound(Me%Ext%Properties, 1)
        UB = ubound(Me%Ext%Properties, 1)
       
        PropertyID = GetPropertyIDNumber (PropertyName) 
        call SearchSinkProperty (PropertyID, PropertyX, STAT_CALL)
        
        if  (STAT_CALL .NE. SUCCESS_) then
        
            Found = .false.        
            do Index = LB, UB
                if (Me%Ext%Properties(Index) .EQ. PropertyID) then
                    Found = .true.
                    exit
                endif                
            end do
       
            if (.NOT. Found) then 
                write(*,*)
                write(*,*) 'Property ', trim(PropertyName), ' was not found in the ', trim(Me%CallingModule)
                write(*,*) 'list of properties '
                stop 'ConstructPairProperty - ModuleChainReactions - ERR010' 
            endif
            
            allocate (NewProperty)
            call ResetSinkProperty (NewProperty)
            
            NewProperty%ID%Name     = PropertyName
            NewProperty%ID%IDNumber = PropertyID 
            NewProperty%Phase       = Phase
            
            call AddSinkProperty (NewProperty)
            
            STAT_ = SUCCESS_
            
        else
        
            STAT_ = SUCCESS_
            
        endif
        
        if (present(STAT)) STAT = STAT_                
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructPairProperty
    !------------------------------------------------------------------------------------------------------------------------------        
    

    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPartitionValues (NewPartition)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer          :: NewPartition
        
        !Local---------------------------------------------------------------------------------------------------------------------
        integer                             :: iflag
        integer                             :: STAT_CALL

        !Begin---------------------------------------------------------------------------------------------------------------------             
        if (NewPartition%PartitionEvolution .EQ. ConstantPartition) then
        
            call GetData(NewPartition%PartitionCoef,            &   
                         Me%ObjEnterData, iflag,                &
                         SearchType   = FromBlock,              &
                         keyword      = 'PARTITION_COEF',       &
                         ClientModule = 'ModuleChainReactions', &
                         STAT         = STAT_CALL)
                         
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructPartitionValues - ModuleChainReactions - ERR010'
            
            if (iFlag .EQ. 0) then
                write(*,*)
                write(*,*) 'You MUST provide a value for keyword PARTITION_COEF for '
                write(*,*) 'property ', trim(NewPartition%ID%Name)
                stop 'ConstructPartitionValues - ModuleChainReactions - ERR020'
            endif

        else
               
            if (Me%Options%GeometryType .EQ. Geometry2D) then
                call ConstructPartitionValues2D (NewPartition)
            else
                call ConstructPartitionValues3D (NewPartition)
            endif    
                    
        endif    
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructPartitionValues
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPartitionValues2D (NewPartition)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer          :: NewPartition

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: ILB,IUB
        integer                             :: JLB,JUB
        integer                             :: WorkSizeILB, WorkSizeIUB
        integer                             :: WorkSizeJLB, WorkSizeJUB       

        !Begin---------------------------------------------------------------------------------------------------------------------             
        !Boundaries
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        WorkSizeILB = Me%WorkSize%ILB
        WorkSizeIUB = Me%WorkSize%IUB
        WorkSizeJLB = Me%WorkSize%JLB
        WorkSizeJUB = Me%WorkSize%JUB

        allocate(NewPartition%G2D%PartitionCoef(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR010'
        NewPartition%G2D%PartitionCoef(:,:) = 0.0
    
        !Get water points
        call GetBasinPoints (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR020'
    
        call ConstructFillMatrix (PropertyID       = NewPartition%ID,                   &
                                  EnterDataID      = Me%ObjEnterData,                   &
                                  TimeID           = Me%ObjTime,                        &
                                  HorizontalGridID = Me%ObjHorizontalGrid,              &
                                  ExtractType      = FromBlock,                         &
                                  PointsToFill2D   = Me%Ext%WaterPoints2D,              &
                                  Matrix2D         = NewPartition%G2D%PartitionCoef,    &
                                  TypeZUV          = TypeZ_,                            &
                                  STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR030'

        if(.NOT. NewPartition%ID%SolutionFromFile)then

            call KillFillMatrix(NewPartition%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR040'
                
        end if  
        
        call UnGetBasin (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR050'                              
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructPartitionValues2D
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPartitionValues3D (NewPartition)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer :: NewPartition

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: ILB,IUB
        integer                   :: JLB,JUB
        integer                   :: KLB,KUB        
        integer                   :: WorkSizeILB, WorkSizeIUB
        integer                   :: WorkSizeJLB, WorkSizeJUB       
        integer                   :: WorkSizeKLB, WorkSizeKUB

        !Begin---------------------------------------------------------------------------------------------------------------------             
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

        allocate(NewPartition%G3D%PartitionCoef(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR010'
        NewPartition%G3D%PartitionCoef(:,:,:) = 0.0
    
    
        !Get water points
        call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR020'
                    
        call ConstructFillMatrix (PropertyID       = NewPartition%ID,                   &
                                  EnterDataID      = Me%ObjEnterData,                   &       
                                  TimeID           = Me%ObjTime,                        &
                                  HorizontalGridID = Me%ObjHorizontalGrid,              &       
                                  GeometryID       = Me%ObjGeometry,                    &
                                  ExtractType      = FromBlockInBlock,                  &
                                  PointsToFill3D   = Me%Ext%WaterPoints3D,              &
                                  Matrix3D         = NewPartition%G3D%PartitionCoef,    &
                                  TypeZUV          = TypeZ_,                            &
                                  STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR030'

        if(.NOT. NewPartition%ID%SolutionFromFile)then

            call KillFillMatrix(NewPartition%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR040'
                                          
        end if    

        call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR050'                
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructPartitionValues3D
    !------------------------------------------------------------------------------------------------------------------------------            
        
    
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine AddPartition (NewPartition)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer :: NewPartition
        
        !Begin---------------------------------------------------------------------------------------------------------------------             
        if (.NOT. associated(Me%FirstPartition)) then           
            Me%FirstPartition => NewPartition
            Me%LastPartition  => NewPartition
        else
            NewPartition%Prev     => Me%LastPartition
            Me%LastPartition%Next => NewPartition
            Me%LastPartition      => NewPartition
        end if 
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine AddPartition
    !------------------------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine ResetSinkProperty (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: NewProperty
        
        !Begin---------------------------------------------------------------------------------------------------------------------    
        nullify(NewProperty%Prev)
        nullify(NewProperty%Next)
        nullify(NewProperty%Product)        

        NewProperty%IsOnlyProduct = .true.
        NewProperty%RateOrder     = NoReaction        
        NewProperty%RateMethod    = RateUnknown
        
        NewProperty%SinkEvolution = NoSink
        NewProperty%SinkRate      = 0.0
        
        call ResetMass (NewProperty)      
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ResetSinkProperty
    !------------------------------------------------------------------------------------------------------------------------------ 
    
    
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ResetMass (PropertyX)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: PropertyX        
        
        !Begin---------------------------------------------------------------------------------------------------------------------
        PropertyX%Sink    = 0.0
        PropertyX%Mass    = 0.0
        PropertyX%NewMass = 0.0
        PropertyX%Source  = 0.0    
        
        PropertyX%NewConcCalculated = .false. 
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ResetMass
    !------------------------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine CreatePropertiesList
    
        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: Index
        type(T_Property), pointer :: PropertyX

        !Begin---------------------------------------------------------------------------------------------------------------------        
        allocate (Me%PropertiesList(Me%PropertiesCount))
        
        PropertyX => Me%FirstProperty
        
        Index = 1
        do while (associated(PropertyX))

            Me%PropertiesList(Index) = PropertyX%ID%IDNumber
            PropertyX => PropertyX%Next
            
            Index = Index + 1

        enddo                
        !--------------------------------------------------------------------------------------------------------------------------        
    
    end subroutine CreatePropertiesList
    !------------------------------------------------------------------------------------------------------------------------------ 


    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine SearchSinkProperty (PropertyID, PropertyX, STAT)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                     :: PropertyID
        type(T_Property), pointer   :: PropertyX
        integer, optional           :: STAT
        
        !Local---------------------------------------------------------------------------------------------------------------------
        integer                     :: STAT_CALL
        
        !Begin---------------------------------------------------------------------------------------------------------------------        
        
        PropertyX => Me%FirstProperty

        do while (associated(PropertyX) )
            if (PropertyX%ID%IDNumber .EQ. PropertyID) then
                exit
            endif

            PropertyX => PropertyX%Next
        enddo        
        
        if (associated(PropertyX)) then
            STAT_CALL = SUCCESS_  
        else
            STAT_CALL = NOT_FOUND_ERR_  
        end if 

        if (present(STAT)) STAT = STAT_CALL                
        !--------------------------------------------------------------------------------------------------------------------------        
        
    end subroutine SearchSinkProperty

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !------------------------------------------------------------------------------------------------------------------------------   
    subroutine GetCRPropertiesList (ChainReactionsID, PropertiesList, PropertiesCount, STAT)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                        :: ChainReactionsID
        integer, pointer, dimension(:) :: PropertiesList
        integer, intent(OUT)           :: PropertiesCount                       
        integer, optional, intent(OUT) :: STAT
 
        !Local---------------------------------------------------------------------------------------------------------------------
        integer :: ready_                     
        integer :: STAT_CALL              !Auxiliar local variable
        
        !Begin---------------------------------------------------------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)    
        
cd1:    if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mCHAINREACTIONS_, Me%InstanceID)

            PropertiesList  => Me%PropertiesList
            PropertiesCount =  Me%PropertiesCount 

            STAT_CALL = SUCCESS_

        else 
        
            STAT_CALL = ready_
            
        end if cd1

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine GetCRPropertiesList
    !------------------------------------------------------------------------------------------------------------------------------    
    
    
    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine SetPropertyConcentration2D (ChainReactionsID, PropertyID, Concentration, STAT)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                        :: ChainReactionsID
        integer                        :: PropertyID
        real, dimension (:,:), pointer :: Concentration
        integer, intent(OUT), optional :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: ready_       
        type(T_Property), pointer :: PropertyX

        !Begin---------------------------------------------------------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchSinkProperty(PropertyID, PropertyX, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'SetPropertyConcentrations2D - ModuleChainReactions - ERR010'
            
            PropertyX%G2D%Concentration => Concentration

            STAT_CALL = SUCCESS_
            
        else
                       
            STAT_CALL = ready_
            
        end if

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine SetPropertyConcentration2D
    !------------------------------------------------------------------------------------------------------------------------------    
    
    
    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine SetPropertyConcentration3D (ChainReactionsID, PropertyID, Concentration, STAT)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                          :: ChainReactionsID
        integer                          :: PropertyID
        real, dimension (:,:,:), pointer :: Concentration
        integer, intent(OUT), optional   :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: ready_       
        type(T_Property), pointer :: PropertyX

        !Begin---------------------------------------------------------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchSinkProperty(PropertyID, PropertyX, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'SetPropertyConcentrations3D - ModuleChainReactions - ERR010'
            
            PropertyX%G3D%Concentration => Concentration

            STAT_CALL = SUCCESS_
            
        else
                       
            STAT_CALL = ready_
            
        end if

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine SetPropertyConcentration3D
    !------------------------------------------------------------------------------------------------------------------------------    
        
    
    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine UnGetChainReactions(ChainReactionsID, Array, STAT)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                        :: ChainReactionsID
        integer, dimension(:), pointer :: Array
        integer, intent(OUT), optional :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer :: STAT_CALL, ready_

        !Begin---------------------------------------------------------------------------------------------------------------------

        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mCHAINREACTIONS_, Me%InstanceID, "UnGetChainReactions")

            STAT_CALL = SUCCESS_
            
        else
                       
            STAT_CALL = ready_
            
        end if

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine UnGetChainReactions
    !------------------------------------------------------------------------------------------------------------------------------
    
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyChainReactions2D (ChainReactionsID,    &
                                       Theta,               &
                                       SoilDensity,         &
                                       DT,                  &
                                       STAT)  

        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                                 :: ChainReactionsID
        real, pointer, dimension(:,:)           :: Theta
        real, pointer, dimension(:,:), optional :: SoilDensity
        real                                    :: DT
        integer, optional,  intent(OUT)         :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: I, J  
        integer                         :: ILB,IUB
        integer                         :: JLB,JUB
        integer                         :: ready_
        type(T_Property), pointer       :: PropertyX, ProductX
        type(T_Partition), pointer      :: PairX
        real, pointer, dimension(:,:)   :: k1, k2
        real                            :: ChangeOfMass
        real                            :: Kd
        real                            :: SinkRate

        !Begin---------------------------------------------------------------------------------------------------------------------             
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB

            call ActualizePropertiesFromFile

            Me%Ext%Theta2D => Theta
            if (present(SoilDensity)) then
                Me%Ext%SoilDensity2D => SoilDensity
            else
                Me%Ext%SoilDensity2D => Theta
            endif
            
            !Get water points
            call GetBasinPoints (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) &
                stop 'ModifyChainReactions2D - ModuleChainReactions - ERR010'
            
            do I = ILB, IUB
            do J = JLB, JUB
                                
                if (Me%Ext%WaterPoints2D(I, J) .EQ. WaterPoint) then
                            
                    PropertyX => Me%FirstProperty
                    
                    do while (associated(PropertyX))
                                                
                        if (PropertyX%Phase .EQ. WaterPhase) then
                            PropertyX%Mass = PropertyX%G2D%Concentration(I, J) * Me%Ext%Theta2D(i, j)
                        else
                            PropertyX%Mass = PropertyX%G2D%Concentration(I, J) * Me%Ext%SoilDensity2D(I, J)
                        endif
                                                 
                        if (.NOT. PropertyX%IsOnlyProduct) then
                                               
                            ProductX => PropertyX%Product
                            
                            call GetSinkRate(PropertyX, SinkRate, I, J)
                                                    
                            if (PropertyX%RateOrder .EQ. ZeroOrder) then
                                                       
                                PropertyX%Sink = min(SinkRate * Me%Ext%Theta2D(i, j) * DT / 86400, PropertyX%Mass)
                                
                            elseif (PropertyX%RateOrder .EQ. FirstOrder) then
                                select case (PropertyX%RateMethod)
                                case (RateDefault)
                                    PropertyX%Sink = min(PropertyX%Mass - PropertyX%Mass * &
                                                         exp(-SinkRate * DT / 86400), PropertyX%Mass)
                                case (RateAlternative)
                                    PropertyX%Sink = min(PropertyX%Mass - PropertyX%Mass * &
                                                        ((1 - SinkRate) ** (DT / 86400)), PropertyX%Mass)
                                end select
                            endif
                                                    
                            !     mg        =       mg                 mg      
                            ProductX%Source = ProductX%Source + PropertyX%Sink
                                                
                        endif                                                
                                                
                        PropertyX => PropertyX%Next
                        
                    enddo
                                            
                    PairX => Me%FirstPartition
                    
                    do while (associated(PairX))
                    
                        call GetPartitionCoef (PairX, Kd, I, J)
                    
                        ChangeOfMass = (PairX%Pair1%Sink + PairX%Pair2%Sink) - (PairX%Pair1%Source + PairX%Pair2%Source) 
                        
                        if (PairX%Pair1%Phase .EQ. WaterPhase) then
                            k1 => Me%Ext%Theta2D
                        else
                            k1 => Me%Ext%SoilDensity2D
                        endif
                        
                        if (PairX%Pair2%Phase .EQ. WaterPhase) then
                            k2 => Me%Ext%Theta2D
                        else
                            k2 => Me%Ext%SoilDensity2D
                        endif

                        PairX%Pair1%G2D%Concentration(I, J) = - ChangeOfMass / (k1(I, J) + Kd * k2(I, J)) + &
                                                                PairX%Pair1%G2D%Concentration(I, J)
                        PairX%Pair2%G2D%Concentration(I, J) = Kd * PairX%Pair1%G2D%Concentration(I, J) 
                    
                        PairX => PairX%Next
                    
                    enddo


                    PropertyX => Me%FirstProperty
                    
                    do while (associated(PropertyX))
                        
                        if (.NOT. PropertyX%NewConcCalculated) then
                        
                            PropertyX%NewMass = PropertyX%Mass + PropertyX%Source - PropertyX%Sink
                            
                            if (PropertyX%Phase .EQ. WaterPhase) then
                                PropertyX%G2D%Concentration(I, J) = PropertyX%NewMass / Me%Ext%Theta2D(I, J)
                            else
                                PropertyX%G2D%Concentration(I, J) = PropertyX%NewMass / Me%Ext%SoilDensity2D(I, J)
                            endif
                        
                        endif
                        
                        call ResetMass(PropertyX)
                                                
                        PropertyX => PropertyX%Next
                        
                    enddo

                endif 
            
            enddo
            enddo
                                      
            call UnGetBasin (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ModifyChainReactions2D - ModuleChainReactions - ERR020'                    
                                                         
            STAT_CALL = SUCCESS_
            
        else              
         
            STAT_CALL = ready_
            
        end if 

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ModifyChainReactions2D
    !------------------------------------------------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ModifyChainReactions3D (ChainReactionsID,    &
                                       Theta,               &
                                       SoilDensity,         &
                                       DT,                  &
                                       STAT)  

        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                                     :: ChainReactionsID
        real, pointer, dimension(:,:,:)             :: Theta
        real, pointer, dimension(:,:,:), optional   :: SoilDensity
        real                                        :: DT
        integer, optional,  intent(OUT)             :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: I, J, K  
        integer                         :: ILB, IUB
        integer                         :: JLB, JUB
        integer                         :: KLB, KUB
        integer                         :: ready_
        type(T_Property), pointer       :: PropertyX, ProductX
        type(T_Partition), pointer      :: PairX
        real, pointer, dimension(:,:,:) :: k1, k2
        real                            :: ChangeOfMass
        real                            :: Kd
        real                            :: SinkRate

        !Begin---------------------------------------------------------------------------------------------------------------------             
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            ILB = Me%WorkSize%ILB
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB
            JUB = Me%WorkSize%JUB
            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB

            call ActualizePropertiesFromFile

            Me%Ext%Theta3D          => Theta
            if (present(SoilDensity)) then
                Me%Ext%SoilDensity3D => SoilDensity
            else
                Me%Ext%SoilDensity3D => Theta
            endif
                        
            !Get water points
            call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) &
                stop 'ModifyChainReactions2D - ModuleChainReactions - ERR010'
            
            do I = ILB, IUB
            do J = JLB, JUB
            do K = KLB, KUB
                                
                if (Me%Ext%WaterPoints3D(I, J, K) .EQ. WaterPoint) then
                            
                    PropertyX => Me%FirstProperty
                    
                    do while (associated(PropertyX))
                                                
                        if (PropertyX%Phase .EQ. WaterPhase) then
                            PropertyX%Mass = PropertyX%G3D%Concentration(I, J, K) * Me%Ext%Theta3D(i, j, K)
                        else
                            PropertyX%Mass = PropertyX%G3D%Concentration(I, J, K) * Me%Ext%SoilDensity3D(I, J, K)
                        endif
                                                 
                        if (.NOT. PropertyX%IsOnlyProduct) then
                                               
                            ProductX => PropertyX%Product
                            
                            call GetSinkRate(PropertyX, SinkRate, I, J, K)
                                                    
                            if (PropertyX%RateOrder .EQ. ZeroOrder) then
                                                       
                                PropertyX%Sink = min(SinkRate * Me%Ext%Theta3D(I, J, K) * DT / 86400, PropertyX%Mass)
                                
                            elseif (PropertyX%RateOrder .EQ. FirstOrder) then
                                select case (PropertyX%RateMethod)
                                case (RateDefault)
                                    PropertyX%Sink = min(PropertyX%Mass - PropertyX%Mass * &
                                                         exp(-SinkRate * DT / 86400), PropertyX%Mass)
                                case (RateAlternative)
                                    PropertyX%Sink = min(PropertyX%Mass - PropertyX%Mass * &
                                                        ((1 - SinkRate) ** (DT / 86400)), PropertyX%Mass)
                                end select
                            endif
                                                    
                            !     mg        =       mg                 mg      
                            ProductX%Source = ProductX%Source + PropertyX%Sink
                                                
                        endif                                                
                                                
                        PropertyX => PropertyX%Next
                        
                    enddo
                                            
                    PairX => Me%FirstPartition
                    
                    do while (associated(PairX))
                    
                        call GetPartitionCoef (PairX, Kd, I, J, K)
                    
                        ChangeOfMass = (PairX%Pair1%Sink + PairX%Pair2%Sink) - (PairX%Pair1%Source + PairX%Pair2%Source) 
                        
                        if (PairX%Pair1%Phase .EQ. WaterPhase) then
                            k1 => Me%Ext%Theta3D
                        else
                            k1 => Me%Ext%SoilDensity3D
                        endif
                        
                        if (PairX%Pair2%Phase .EQ. WaterPhase) then
                            k2 => Me%Ext%Theta3D
                        else
                            k2 => Me%Ext%SoilDensity3D
                        endif

                        PairX%Pair1%G3D%Concentration(I, J, K) = - ChangeOfMass / (k1(I, J, K) + Kd * k2(I, J, K)) + &
                                                                   PairX%Pair1%G3D%Concentration(I, J, K)
                        PairX%Pair2%G3D%Concentration(I, J, K) = Kd * PairX%Pair1%G3D%Concentration(I, J, K) 
                    
                        PairX => PairX%Next
                    
                    enddo


                    PropertyX => Me%FirstProperty
                    
                    do while (associated(PropertyX))
                        
                        if (.NOT. PropertyX%NewConcCalculated) then
                        
                            PropertyX%NewMass = PropertyX%Mass + PropertyX%Source - PropertyX%Sink
                            
                            if (PropertyX%Phase .EQ. WaterPhase) then
                                PropertyX%G3D%Concentration(I, J, K) = PropertyX%NewMass / Me%Ext%Theta3D(I, J, K)
                            else
                                PropertyX%G3D%Concentration(I, J, K) = PropertyX%NewMass / Me%Ext%SoilDensity3D(I, J, K)
                            endif
                        
                        endif
                        
                        call ResetMass(PropertyX)
                                                
                        PropertyX => PropertyX%Next
                        
                    enddo

                endif 
            
            enddo
            enddo
            enddo
                                      
            call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR020'                    
                                                         
            STAT_CALL = SUCCESS_
            
        else              
         
            STAT_CALL = ready_
            
        end if 

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ModifyChainReactions3D
    !------------------------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ActualizePropertiesFromFile
    
        !Local---------------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer   :: PropertyX 
        type(T_Partition), pointer  :: PartitionX
        integer                     :: STAT_CALL   
        
        !Begin---------------------------------------------------------------------------------------------------------------------    
        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if ((PropertyX%SinkEvolution .EQ. VariableSink) .AND. (PropertyX%ID%SolutionFromFile)) then
            
                if (Me%Options%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                           Matrix3D       = PropertyX%G3D%SinkRate,     &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,       &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR010'
                else
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                           Matrix2D       = PropertyX%G2D%SinkRate,     &       
                                           PointsToFill2D = Me%Ext%WaterPoints2D,       &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR020'
                endif
                            
            endif
            
            PropertyX => PropertyX%Next
            
        enddo    
        
        PartitionX => Me%FirstPartition

        do while (associated(PartitionX))

            if ((PartitionX%PartitionEvolution .EQ. VariablePartition) .AND. (PropertyX%ID%SolutionFromFile)) then
            
                if (Me%Options%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PartitionX%ID%ObjFillMatrix,    &
                                           Matrix3D       = PartitionX%G3D%PartitionCoef,   &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,           &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR030'
                else
                    call ModifyFillMatrix (FillMatrixID   = PartitionX%ID%ObjFillMatrix,    &
                                           Matrix2D       = PartitionX%G2D%PartitionCoef,   &
                                           PointsToFill2D = Me%Ext%WaterPoints2D,           &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR040'
                endif
                            
            endif
            
            PartitionX => PartitionX%Next
            
        enddo            
        !--------------------------------------------------------------------------------------------------------------------------    
    
    end subroutine ActualizePropertiesFromFile    
    !------------------------------------------------------------------------------------------------------------------------------    


    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine GetSinkRate (PropertyX, SinkRate, I, J, K)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer   :: PropertyX
        real                        :: SinkRate
        integer                     :: I, J
        integer, optional           :: K
    
        !Begin---------------------------------------------------------------------------------------------------------------------
        if (PropertyX%SinkEvolution .EQ. ConstantSink) then
        
            SinkRate = PropertyX%SinkRate
            
        else
        
            if (present(K)) then
            
                SinkRate = PropertyX%G3D%SinkRate(I, J, K)
            
            else
            
                SinkRate = PropertyX%G2D%SinkRate(I, J)
            
            endif
        
        endif
        !--------------------------------------------------------------------------------------------------------------------------    
        
    end subroutine GetSinkRate
    !------------------------------------------------------------------------------------------------------------------------------    

        
    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine GetPartitionCoef (PartitionX, PartitionCoef, I, J, K)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Partition), pointer  :: PartitionX
        real                        :: PartitionCoef
        integer                     :: I, J
        integer, optional           :: K
    
        !Begin---------------------------------------------------------------------------------------------------------------------
        if (PartitionX%PartitionEvolution .EQ. ConstantPartition) then
        
            PartitionCoef = PartitionX%PartitionCoef
            
        else
        
            if (present(K)) then
            
                PartitionCoef = PartitionX%G3D%PartitionCoef(I, J, K)
            
            else
            
                PartitionCoef = PartitionX%G2D%PartitionCoef(I, J)
            
            endif
        
        endif
        !--------------------------------------------------------------------------------------------------------------------------    
        
    end subroutine GetPartitionCoef
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine KillChainReactions(ChainReactionsID, STAT)

        !Arguments---------------------------------------------------------------
        integer                         :: ChainReactionsID
        integer, optional, intent(OUT)  :: STAT

        !External----------------------------------------------------------------
        integer :: ready_ 

        !Local-------------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: nUsers 
        type(T_Property), pointer       :: PropertyX
        type(T_Partition), pointer      :: PartitionX

        !------------------------------------------------------------------------                      

        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mCHAINREACTIONS_,  Me%InstanceID)

            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 
            
                if((PropertyX%SinkEvolution .EQ. VariableSink) .AND. PropertyX%ID%SolutionFromFile) then

                    call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillChainReactions - ModuleChainReactions - ERR010'
                end if
                
                PropertyX => PropertyX%Next
                
            end do 

            PartitionX => Me%FirstPartition
            
            do while (associated(PartitionX)) 
                 
                if((PartitionX%PartitionEvolution .EQ. VariablePartition) .AND. PartitionX%ID%SolutionFromFile) then

                    call KillFillMatrix(PartitionX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillChainReactions - ModuleChainReactions - ERR020'
                end if

                PartitionX => PartitionX%Next
                
            end do 

cd2 :       if (nUsers == 0) then
                              
                nUsers = DeassociateInstance (mTIME_, Me%ObjTime)
                if (nUsers == 0) stop 'KillChainReactions - ModuleChainReactions - ERR030'

                nUsers = DeassociateInstance (mBASINGEOMETRY_, Me%ObjBasinGeometry)
                if (nUsers == 0) stop 'KillChainReactions - ModuleChainReactions - ERR040'

                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillChainReactions - ModuleChainReactions - ERR050'

                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillChainReactions - ModuleChainReactions - ERR060'
                
                nUsers = DeassociateInstance (mGEOMETRY_,  Me%ObjGeometry)
                if (nUsers == 0) stop 'KillChainReactions - ModuleChainReactions - ERR070'

                nUsers = DeassociateInstance (mMAP_,  Me%ObjMap)
                if (nUsers == 0) stop 'KillChainReactions - ModuleChainReactions - ERR080'
                                                                            
                call DeallocateInstance 

                ChainReactionsID = 0
                STAT_CALL = SUCCESS_

            end if cd2
        
        else cd1
        
            STAT_CALL = ready_
        
        end if cd1


        if (present(STAT)) STAT = STAT_CALL
    !------------------------------------------------------------------------

    end subroutine KillChainReactions
    !------------------------------------------------------------------------


    !------------------------------------------------------------------------
    subroutine DeallocateInstance 

        !Local-----------------------------------------------------------------
        type (T_ChainReactions), pointer :: AuxObjChainReactions
        type (T_ChainReactions), pointer :: PreviousObjChainReactions

        !Updates pointers
        if (Me%InstanceID == FirstObjChainReactions%InstanceID) then
        
            FirstObjChainReactions => FirstObjChainReactions%Next
            
        else
        
            PreviousObjChainReactions => FirstObjChainReactions
            AuxObjChainReactions      => FirstObjChainReactions%Next
            
            do while (AuxObjChainReactions%InstanceID /= Me%InstanceID)
            
                PreviousObjChainReactions => AuxObjChainReactions
                AuxObjChainReactions      => AuxObjChainReactions%Next
                
            enddo

            !Now update linked list
            PreviousObjChainReactions%Next => AuxObjChainReactions%Next

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


    !----------------------------------------------------------------------------
    subroutine Ready (ChainReactionsID, ready_) 

        !Arguments-------------------------------------------------------------
        integer :: ChainReactionsID
        integer :: ready_
        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ChainReactionsID > 0) then
            
            call LocateObjChainReactions (ChainReactionsID)
            ready_ = VerifyReadLock (mCHAINREACTIONS_, Me%InstanceID)

        else

            ready_ = OFF_ERR_

        end if cd1
        !----------------------------------------------------------------------

    end subroutine Ready
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine LocateObjChainReactions (ChainReactionsID)

    !Arguments-------------------------------------------------------------
        integer :: ChainReactionsID
    !--------------------------------------------------------------------------

        Me => FirstObjChainReactions
        do while (associated (Me))
            if (Me%InstanceID == ChainReactionsID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))   &
            stop 'Subroutine LocateObjChainReactions; Module ModuleChainReactions. ERR001.'
    !--------------------------------------------------------------------------
    
    end subroutine LocateObjChainReactions
    !--------------------------------------------------------------------------


end module ModuleChainReactions