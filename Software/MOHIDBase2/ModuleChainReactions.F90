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
!
! UNITS:
!     ZERO order rate  : mg/s 
!     FIRST order rate : 1/s
!

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
    private ::              ConstructPropertiesList
    private ::                  ConstructProperty
    private ::                      ConstructPropertyKd
    private ::                      ConstructPropertySink
    private ::                          ConstructValues
    private ::                  AddPropertyToList
    private ::                  MakeConnections
    private ::              CheckProperties 
    private ::      ResetProperty
    private ::          ResetCalculatedValues
    private ::      CreatePropertiesList  
    private ::      SearchProperty 
    public  :: InitCRSoilPhase                     
                
    !Selector
    public  :: GetCRPropertiesList
    public  :: SetCRPropertyConcentration
    private ::      SetPropertyConcentration2D
    private ::      SetPropertyConcentration3D
    public  :: UngetChainReactions           
                        
    !Modifier
    public  :: ModifyChainReactions 
    !private ::      ModifyChainReactions2D
    private ::      ModifyChainReactions3D
    
    private ::          ActualizePropertiesFromFile
!    private ::          GetSinkrate
        
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
!        module procedure ModifyChainReactions2D
        module procedure ModifyChainReactions3D
    end interface ModifyChainReactions
        
    !Constants---------------------------------------------------------------------------------------------------------------------
    integer,                       parameter :: MaxStrLength       = StringLength      !StringLength is defined in ModuleGlobalData     
    
    integer,                       parameter :: UnknownGeometry    = 0
    integer,                       parameter :: Geometry2D         = 1
    integer,                       parameter :: Geometry3D         = 2    
    
    integer,                       parameter :: NoEvolution        = 0
    integer,                       parameter :: ConstantEvolution  = 1 !Single value
    integer,                       parameter :: VariableEvolution  = 2
       
    character(LEN = StringLength), parameter :: prop_block_begin   = '<begin_property>'
    character(LEN = StringLength), parameter :: prop_block_end     = '<end_property>'
       
    character(LEN = StringLength), parameter :: sink_block_begin   = '<begin_sink>'
    character(LEN = StringLength), parameter :: sink_block_end     = '<end_sink>'

    character(LEN = StringLength), parameter :: kd_block_begin     = '<begin_kd>'
    character(LEN = StringLength), parameter :: kd_block_end       = '<end_kd>'

    !Types-------------------------------------------------------------------------------------------------------------------------
                
    type T_PropertyAttribute
        integer                         :: Evolution    = NoEvolution
        real                            :: DefaultValue = 0.0 
        real, pointer, dimension(:,:)   :: Values2D       
        real, pointer, dimension(:,:,:) :: Values3D       
    end type T_PropertyAttribute             
                 
    type T_SPConnection
        type(T_Property), pointer     :: Source      => null()
        integer                       :: ProductID 
        character(MaxStrLength)       :: ProductName  
        type(T_SPConnection), pointer :: Next        => null()
        type(T_SPConnection), pointer :: Prev        => null()
    end type T_SPConnection
    
    type T_Sink
        type(T_PropertyAttribute) :: Soil_0  !ZERO order sink from soil phase
        type(T_PropertyAttribute) :: Soil_1  !FIRST order sink from soil phase
        type(T_PropertyAttribute) :: Water_0 !ZERO order sink from water phase
        type(T_PropertyAttribute) :: Water_1 !FIRST order sink from the water phase
    end type T_Sink
        
    type T_Calc
        real :: OldSoilMass  = 0.0
        real :: SoilMassLost = 0.0
        real :: OldMass      = 0.0
        real :: MassLost     = 0.0
        real :: MassGained   = 0.0
    end type T_Calc
        
    type T_SoilPhase
        type(T_PropertyAttribute) :: Kd   
        type(T_PropertyAttribute) :: Mass 
    end type T_SoilPhase
        
    type T_Property
        type(T_PropertyID)        :: ID
        type(T_Sink)              :: Sink
        type(T_SoilPhase)         :: Soil
        type(T_PropertyAttribute) :: Concentration 
        type(T_Property), pointer :: Product => null()
        type(T_Calc)              :: Calc
        type(T_Property), pointer :: Next => null()
        type(T_Property), pointer :: Prev => null()      
    end type T_Property
       
    type T_External
        integer, pointer, dimension(:,:,:) :: WaterPoints3D
        integer, pointer, dimension(:,:)   :: WaterPoints2D
        real,    pointer, dimension(:,:,:) :: Theta3D
        real,    pointer, dimension(:,:)   :: Theta2D
        real,    pointer, dimension(:,:,:) :: SoilDensity3D
        real,    pointer, dimension(:,:)   :: SoilDensity2D
        integer, pointer, dimension(:)     :: Properties
    end type T_External 
       
    type T_ChainReactions
        private
        integer                             :: InstanceID                    !ID of the ModuleChainReactions instance 
        type(T_ChainReactions), pointer     :: Next => null()                !Collection of instances of ModuleChainReactions
        
        character(MaxStrLength)             :: CallingModule

        integer                             :: ObjEnterData         = 0
        integer                             :: ObjTime              = 0
        integer                             :: ObjGeometry          = 0
        integer                             :: ObjHorizontalGrid    = 0
        integer                             :: ObjHorizontalMap     = 0
        integer                             :: ObjBasinGeometry     = 0
        integer                             :: ObjMap               = 0

        integer                             :: GeometryType         = 0      !0 -> Unknown; 1 -> 2D; 2 -> 3D        
        character(PathLength)               :: DataFile   
                     
        type(T_External)                    :: Ext                           !Pointers to Water Mass, Properties Values 
                                                                             !and other required data         
        type(T_Property), pointer           :: FirstProperty => null()            
        type(T_Property), pointer           :: LastProperty  => null()
        
        type(T_SPConnection), pointer       :: FirstConnection => null()
        type(T_SPConnection), pointer       :: LastConnection  => null()

        integer, pointer, dimension(:)      :: PropertiesList  => null()     !List of properties ID's
        integer                             :: PropertiesCount = 0           !Number of properties
        
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
        
            Me%GeometryType = GeometryType
                     
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
            
            stop 'Subroutine StartChainReactions; ModuleChainReactions. ERR030'

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
        
        nullify(Me%FirstConnection)
        nullify(Me%LastConnection)                             
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
!        integer                 :: STAT_CALL
!        integer                 :: flag
!        character(MaxStrLength) :: Text
        
        !Begin-----------------------------------------------------------------               
        call GetExtractType (FromFile = FromFile)

!        call GetData(Me%Options%RateOrder,                  &
!                     Me%ObjEnterData, flag,                 &
!                     SearchType   = FromFile,               &
!                     keyword      = 'RATE_ORDER',           &
!                     default      = NoReaction,             & 
!                     ClientModule = 'ModuleChainReactions', &
!                     STAT         = STAT_CALL)
!                     
!        if (STAT_CALL .NE. SUCCESS_) &
!            stop 'Subroutine ReadChainReactionsOptions; Module ModuleChainReactions. ERR010.'                
!
!        if (Me%Options%RateOrder .NE. NoReaction) then
!            if (Me%Options%RateOrder .EQ. ZeroOrder) then
!                Text = 'ZERO Order Rate type'
!            else
!                Text = 'FIRST Order Rate type'
!            endif
!            
!            write(*,*)
!            write(*,*) 'Warning: All properties will be set to be of ', trim(Text)
!        endif
!
!        if (Me%Options%RateOrder .EQ. FirstOrder) then
!            
!            call GetData(Me%Options%RateMethod,                 &
!                         Me%ObjEnterData, flag,                 &
!                         SearchType   = FromFile,               &
!                         keyword      = 'RATE_METHOD',          &
!                         default      = RateUnknown,            & 
!                         ClientModule = 'ModuleChainReactions', &
!                         STAT         = STAT_CALL)
!                         
!            if (STAT_CALL .NE. SUCCESS_) &
!                stop 'Subroutine ReadChainReactionsOptions; Module ModuleChainReactions. ERR010.'                
!                
!            if (Me%Options%RateMethod .NE. RateUnknown) then
!                if (Me%Options%RateMethod .EQ. RateDefault) then
!                    Text = 'Default First Order equation'
!                else
!                    Text = 'Alternative First Order equation'
!                endif
!                
!                write(*,*)
!                write(*,*) 'Warning: All properties will be set to use the ', trim(Text)
!            endif
!
!        else
!        
!            Me%Options%RateMethod = RateUnknown    
!        
!        endif
                        
        call ConstructPropertiesList  
        call MakeConnections
        call CheckProperties 
        
!        call ConstructPartitionList
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ReadChainReactionsOptions         
    !------------------------------------------------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine ConstructPropertiesList 

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                       :: STAT_CALL
        logical                       :: BlockFound
        type(T_Property), pointer     :: NewProperty

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
                    Call ConstructProperty(NewProperty)

                    !Add new Property to the Properties List 
                    Call AddPropertyToList(NewProperty)                                        
                    
                else

                    call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_) &
                        stop 'ConstructPropertyList - ModuleChainReactions - ERR010'
                    exit !No more blocks
                
                end if    

            elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'ConstructPropertyList - ModuleChainReactions - ERR030'
            
            else    
                
                stop 'ConstructPropertyList - ModuleChainReactions - ERR040'
            
            end if    
        
        end do            
        !--------------------------------------------------------------------------------------------------------------------------
                
    end subroutine ConstructPropertiesList
    !------------------------------------------------------------------------------------------------------------------------------   


    !------------------------------------------------------------------------------------------------------------------------------   
    subroutine ConstructProperty(NewProperty)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL, iflag
        character(StringLength)   :: product_name

        !Begin---------------------------------------------------------------------------------------------------------------------             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructProperty - ModuleChainReactions - ERR010'
        
        call ResetProperty (NewProperty)                        
        call ConstructPropertyID (NewProperty%ID, Me%ObjEnterData, FromBlock)
        call ConstructPropertyKd (newProperty)
        call ConstructPropertySink (NewProperty)
        
        call GetData(product_name,                          &   
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromBlock,              &
                     keyword      = 'PRODUCT_NAME',         &
                     Default      = '',                     &
                     ClientModule = 'ModuleChainReactions', &
                     STAT         = STAT_CALL)                     
        if (STAT_CALL /= SUCCESS_) &
            stop 'ConstructProperty - ModuleChainReactions - ERR020'
                
        if (iFlag .NE. 0) then
            call AddConnectionToList (NewProperty, product_name)
        endif
              
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine ConstructProperty    
    !------------------------------------------------------------------------------------------------------------------------------    


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPropertyKd (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External------------------------------------------------------------------------------------------------------------------
        integer                   :: STAT_CALL, iflag
        logical                   :: BlockInBlockFound

        !Begin---------------------------------------------------------------------------------------------------------------------                 
        call ExtractBlockFromBlock (Me%ObjEnterData,                       &
                                    ClientNumber      = Me%ClientNumber,   &
                                    block_begin       = kd_block_begin,    &
                                    block_end         = kd_block_end,      &
                                    BlockInBlockFound = BlockInBlockFound, &
                                    STAT              = STAT_CALL)
        if (STAT_CALL .EQ. SUCCESS_) then    

            if (BlockInBlockFound) then                                                            
                
                call GetData(NewProperty%Soil%Kd%Evolution,         &   
                             Me%ObjEnterData, iflag,                &
                             SearchType   = FromBlockInBlock,       &
                             keyword      = 'EVOLUTION',            &
                             Default      = ConstantEvolution,      &
                             ClientModule = 'ModuleChainReactions', &
                             STAT         = STAT_CALL)                             
                if (STAT_CALL /= SUCCESS_) &
                    stop 'ConstructPropertyKd - ModuleChainReactions - ERR010'            
                
                call GetData(NewProperty%Soil%Kd%DefaultValue,      &   
                             Me%ObjEnterData, iflag,                &
                             SearchType   = FromBlockInBlock,       &
                             keyword      = 'DEFAULT',              &
                             Default      = 0.0,                    &
                             ClientModule = 'ModuleChainReactions', &
                             STAT         = STAT_CALL)                             
                if (STAT_CALL .NE. SUCCESS_) &
                    stop 'ConstructPropertyKd - ModuleChainReactions - ERR020' 
                                            
                if (Me%GeometryType .EQ. Geometry2D) then
                    call ConstructValues (NewProperty,                                 &
                                          NewProperty%Soil%Kd,                         &
                                          Matrix2D     = NewProperty%Soil%Kd%Values2D, &
                                          DefaultValue = NewProperty%Soil%Kd%DefaultValue)
                                     
                    NewProperty%Soil%Mass%Evolution = ConstantEvolution
                    call ConstructValues (NewProperty,                                            &
                                          NewProperty%Soil%Mass,                         &
                                          Matrix2D     = NewProperty%Soil%Mass%Values2D, &
                                          DefaultValue = 0.0)                                          
                else
                    call ConstructValues (NewProperty,                                 &
                                          NewProperty%Soil%Kd,                         &
                                          Matrix3D     = NewProperty%Soil%Kd%Values3D, &
                                          DefaultValue = NewProperty%Soil%Kd%DefaultValue)
                                          
                    NewProperty%Soil%Mass%Evolution = ConstantEvolution
                    call ConstructValues (NewProperty,                                            &
                                          NewProperty%Soil%Mass,                         &
                                          Matrix3D     = NewProperty%Soil%Mass%Values3D, &
                                          DefaultValue = 0.0)                                            
                endif                 
                
!                select case (NewProperty%kd%Evolution)                
!                    case (VariablePartition)                
!                        if (Me%GeometryType .EQ. Geometry2D) then
!                            call ConstructValues (NewProperty, &
!                                                  NewProperty%Kd%Evolution,              &
!                                                  Matrix2D = NewProperty%Kd%Values2D)
!                        else
!                            call ConstructValues (NewProperty, Matrix3D = NewProperty%Kd%Values3D)
!                        endif    
!                        
!                    case (ConstantPartition)
!                        call GetData(NewProperty%Kd%DefaultValue,           &   
!                                     Me%ObjEnterData, iflag,                &
!                                     SearchType   = FromBlockInBlock,       &
!                                     keyword      = 'DEFAULT',              &
!                                     ClientModule = 'ModuleChainReactions', &
!                                     STAT         = STAT_CALL)
!                                     
!                        if ((STAT_CALL .NE. SUCCESS_) .OR. (iFlag .EQ. 0)) &
!                            stop 'ConstructPropertyKd - ModuleChainReactions - ERR020'            
!                    
!                    case (NoPartition)
!                        !Nothing to do
!                        
!                    case default
!                        stop 'ConstructPropertyKd - ModuleChainReactions - ERR030'
!                end select                
            else

                NewProperty%Soil%Kd%Evolution    = NoEvolution
                NewProperty%Soil%Kd%DefaultValue = 0.0
                
                if (Me%GeometryType .EQ. Geometry2D) then
                    call ConstructValues (NewProperty,                                 &
                                          NewProperty%Soil%Kd,                         &
                                          Matrix2D     = NewProperty%Soil%Kd%Values2D, &
                                          DefaultValue = NewProperty%Soil%Kd%DefaultValue)
                else
                    call ConstructValues (NewProperty,                                 &
                                          NewProperty%Soil%Kd,                         &
                                          Matrix3D     = NewProperty%Soil%Kd%Values3D, &
                                          DefaultValue = NewProperty%Soil%Kd%DefaultValue)
                endif                 
            
            end if  
            
            call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructPropertyKd - ConstructSinkValues - ERR040'              
       
        else    
            
            stop 'ConstructPropertyKd - ModuleChainReactions - ERR050'
        
        end if                   
        !--------------------------------------------------------------------------------------------------------------------------
        
    end subroutine ConstructPropertyKd
    !------------------------------------------------------------------------------------------------------------------------------            

   
    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructPropertySink (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer :: NewProperty

        !External------------------------------------------------------------------------------------------------------------------
        integer                            :: STAT_CALL, iflag
        integer                            :: sink_type
        Type(T_PropertyAttribute), pointer :: sink
        logical                            :: BlockInBlockFound

        !Begin---------------------------------------------------------------------------------------------------------------------                 
        do
            call ExtractBlockFromBlock (Me%ObjEnterData,                       &
                                        ClientNumber      = Me%ClientNumber,   &
                                        block_begin       = sink_block_begin,  &
                                        block_end         = sink_block_end,    &
                                        BlockInBlockFound = BlockInBlockFound, &
                                        STAT              = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then    

                if (BlockInBlockFound) then  
                         
                    call GetData(sink_type,                             &   
                                 Me%ObjEnterData, iflag,                &
                                 SearchType   = FromBlockInBlock,       &
                                 keyword      = 'TYPE',                 &
                                 ClientModule = 'ModuleChainReactions', &
                                 STAT         = STAT_CALL)                             
                    if ((STAT_CALL /= SUCCESS_) .OR. (iflag .EQ. 0)) &
                        stop 'ConstructPropertySink - ModuleChainReactions - ERR010'            
                            
                    select case (sink_type)
                        case (0)
                            sink => NewProperty%Sink%Soil_0
                        case (1)
                            sink => NewProperty%Sink%Soil_1
                        case (2)
                            sink => NewProperty%Sink%Water_0
                        case (3)
                            sink => NewProperty%Sink%Water_1
                        case default
                            stop 'ConstructPropertySink - ModuleChainReactions - ERR020'
                    end select
                            
                    call GetData(sink%Evolution,                        &   
                                 Me%ObjEnterData, iflag,                &
                                 SearchType   = FromBlockInBlock,       &
                                 keyword      = 'EVOLUTION',            &
                                 Default      = ConstantEvolution,      &
                                 ClientModule = 'ModuleChainReactions', &
                                 STAT         = STAT_CALL)                             
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructPropertySink - ModuleChainReactions - ERR030'            
                
                    call GetData(sink%DefaultValue,                     &   
                                 Me%ObjEnterData, iflag,                &
                                 SearchType   = FromBlockInBlock,       &
                                 keyword      = 'DEFAULT',              &
                                 Default      = 0.0,                    &
                                 ClientModule = 'ModuleChainReactions', &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructPropertySink - ModuleChainReactions - ERR040'  
                                      
                    if (Me%GeometryType .EQ. Geometry2D) then
                        call ConstructValues (NewProperty,                  &
                                              sink,                         &
                                              Matrix2D     = sink%Values2D, &
                                              DefaultValue = sink%DefaultValue)
                    else
                        call ConstructValues (NewProperty,                  &
                                              sink,                         &
                                              Matrix3D     = sink%Values3D, &
                                              DefaultValue = sink%DefaultValue)
                    endif    
!                                                                         
!                    select case (sink%Evolution)                
!                        case (VariableSink)                
!                            if (Me%GeometryType .EQ. Geometry2D) then
!                                call ConstructValues (NewProperty, Matrix2D = sink%Values2D)
!                            else
!                                call ConstructValues (NewProperty, Matrix3D = sink%Values3D)
!                            endif    
!                            
!                        case (ConstantPartition)
!                            call GetData(sink%DefaultValue,                     &   
!                                         Me%ObjEnterData, iflag,                &
!                                         SearchType   = FromBlockInBlock,       &
!                                         keyword      = 'DEFAULT',              &
!                                         Default      = 0.0,                    &
!                                         ClientModule = 'ModuleChainReactions', &
!                                         STAT         = STAT_CALL)
!                                         
!                            if ((STAT_CALL /= SUCCESS_) .OR. (iFlag .EQ. 0)) &
!                                stop 'ConstructPropertySink - ModuleChainReactions - ERR040'            
!                        
!                        case (NoPartition)
!                            !Nothing to do
!                            
!                        case default
!                            stop 'ConstructPropertySink - ModuleChainReactions - ERR050'
!                    end select                                                      

!                    call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
!                    if (STAT_CALL /= SUCCESS_) &
!                        stop 'ConstructPropertySink - ConstructSinkValues - ERR060' 
                                                                
                else

                    call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'ConstructPropertySink - ConstructSinkValues - ERR070' 

                    exit !No more blocks
                
                end if    
            
            end if
            
        end do
        !--------------------------------------------------------------------------------------------------------------------------
        
    end subroutine ConstructPropertySink
    !------------------------------------------------------------------------------------------------------------------------------            


!    !------------------------------------------------------------------------------------------------------------------------------        
!    subroutine ConstructValues2D (NewProperty, Matrix)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_property), pointer           :: NewProperty
!        real, dimension(:,:), pointer       :: Matrix
!
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                             :: STAT_CALL
!        integer                             :: ILB,IUB
!        integer                             :: JLB,JUB
!        integer                             :: WorkSizeILB, WorkSizeIUB
!        integer                             :: WorkSizeJLB, WorkSizeJUB       
!
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        !Boundaries
!        ILB = Me%Size%ILB
!        IUB = Me%Size%IUB
!        JLB = Me%Size%JLB
!        JUB = Me%Size%JUB
!
!        WorkSizeILB = Me%WorkSize%ILB
!        WorkSizeIUB = Me%WorkSize%IUB
!        WorkSizeJLB = Me%WorkSize%JLB
!        WorkSizeJUB = Me%WorkSize%JUB
!
!        allocate(NewProperty%G2D%SinkRate(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) &
!            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR010'
!        Matrix(:,:) = 0.0
!    
!        !Get water points
!        call GetBasinPoints (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR020'
!    
!        call ConstructFillMatrix (PropertyID       = NewProperty%ID,            &
!                                  EnterDataID      = Me%ObjEnterData,           &
!                                  TimeID           = Me%ObjTime,                &
!                                  HorizontalGridID = Me%ObjHorizontalGrid,      &
!                                  ExtractType      = FromBlockInBlock,          &
!                                  PointsToFill2D   = Me%Ext%WaterPoints2D,      &
!                                  Matrix2D         = Matrix,                    &
!                                  TypeZUV          = TypeZ_,                    &
!                                  STAT             = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR030'
!
!        if(.NOT. NewProperty%ID%SolutionFromFile)then
!
!            call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_)&
!                stop 'ConstructSinkValues2D - ModuleChainReactions - ERR040'
!                
!        end if  
!        
!        call UnGetBasin (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructSinkValues2D - ModuleChainReactions - ERR050'                              
!        !--------------------------------------------------------------------------------------------------------------------------
!    
!    end subroutine ConstructValues2D
!    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------        
    subroutine ConstructValues (NewProperty, Attribute, Matrix2D, Matrix3D, DefaultValue)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_property), pointer                 :: NewProperty
        type(T_PropertyAttribute)                 :: Attribute
        real, dimension(:,:), pointer, optional   :: Matrix2D
        real, dimension(:,:,:), pointer, optional :: Matrix3D
        real, optional                            :: DefaultValue

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

        if (present(Matrix2D) .AND. present(Matrix3D)) &
            stop 'ConstructValues - ModuleChainReactions - ERR010'
        
        if ((.NOT. present(Matrix2D)) .AND. (.NOT. present(Matrix3D))) &
            stop 'ConstructValues - ModuleChainReactions - ERR020'
        
        if (present (Matrix2D)) then
            allocate(Matrix2D(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'ConstructValues - ModuleChainReactions - ERR030'
            Matrix2D(:,:) = 0.0        
        else        
            allocate(Matrix3D(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'ConstructValues - ModuleChainReactions - ERR040'
            Matrix3D(:,:,:) = 0.0
        endif
    
        if (Attribute%Evolution .EQ. VariableEvolution) then
            !Get water points
            call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructValues - ModuleChainReactions - ERR050'
                        
            if (present (Matrix2D)) then
                call ConstructFillMatrix  (PropertyID       = NewProperty%ID,           &
                                           EnterDataID      = Me%ObjEnterData,          &
                                           TimeID           = Me%ObjTime,               &
                                           HorizontalGridID = Me%ObjHorizontalGrid,     &
                                           ExtractType      = FromBlockInBlock,         &
                                           PointsToFill2D   = Me%Ext%WaterPoints2D,     &
                                           Matrix2D         = Matrix2D,                 &
                                           TypeZUV          = TypeZ_,                   &
                                           STAT             = STAT_CALL)
            else
                call ConstructFillMatrix  (PropertyID       = NewProperty%ID,           &
                                           EnterDataID      = Me%ObjEnterData,          &
                                           TimeID           = Me%ObjTime,               &
                                           HorizontalGridID = Me%ObjHorizontalGrid,     &
                                           GeometryID       = Me%ObjGeometry,           &
                                           ExtractType      = FromBlockInBlock,         &
                                           PointsToFill3D   = Me%Ext%WaterPoints3D,     &
                                           Matrix3D         = Matrix3D,                 &
                                           TypeZUV          = TypeZ_,                   &
                                           STAT             = STAT_CALL)
            endif
                                                   
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructValues - ModuleChainReactions - ERR060'

            if(.NOT. NewProperty%ID%SolutionFromFile)then

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructValues - ModuleChainReactions - ERR070'
                
                If (present(DefaultValue) .and. (DefaultValue > 0.0)) then
                
                    Attribute%Evolution = ConstantEvolution
                    
                    if (present (Matrix2D)) then
                        Matrix2D = DefaultValue
                    else
                        Matrix3D = DefaultValue
                    endif        
                
                else
                
                    Attribute%Evolution = NoEvolution
                    
                    if (present (Matrix2D)) then
                        Matrix2D = 0.0
                    else
                        Matrix3D = 0.0
                    endif        
                
                endif
                                              
            end if    

            call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ConstructValues - ModuleChainReactions - ERR080' 
                
        elseif (Attribute%Evolution .EQ. ConstantEvolution) then
            if (.NOT. present(DefaultValue)) &
                stop 'ConstructValues - ModuleChainReactions - ERR090' 
                
            if (present (Matrix2D)) then
                Matrix2D = DefaultValue
            else
                Matrix3D = DefaultValue
            endif
        else
            if (present (Matrix2D)) then
                Matrix2D = 0.0
            else
                Matrix3D = 0.0
            endif        
        endif              
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ConstructValues
    !------------------------------------------------------------------------------------------------------------------------------            


    !------------------------------------------------------------------------------------------------------------------------------    
    subroutine AddPropertyToList(NewProperty)

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

    end subroutine AddPropertyToList
    !------------------------------------------------------------------------------------------------------------------------------ 


    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine AddConnectionToList(NewProperty, ProductName)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: NewProperty
        character(LEN = *)        :: ProductName
        
        !Local---------------------------------------------------------------------------------------------------------------------
        type(T_SPConnection), pointer :: NewConnection
        integer                       :: STAT

        !Begin---------------------------------------------------------------------------------------------------------------------

        allocate (NewConnection, STAT = STAT)
        if (STAT /= 0) &
            stop 'AddConnectionToList - ModuleChainReactions - ERR010'

        NewConnection%Source => NewProperty
        NewConnection%ProductName = ProductName
        NewConnection%ProductID = GetPropertyIDNumber(ProductName)
        
        ! Add to the Property List a new property
        if (.NOT. associated(Me%FirstConnection)) then           
            Me%FirstConnection => NewConnection
            Me%LastConnection  => NewConnection

        else
            NewConnection%Prev     => Me%LastConnection
            Me%LastConnection%Next => NewConnection
            Me%LastConnection      => NewConnection
        end if 
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine AddConnectionToList
    !------------------------------------------------------------------------------------------------------------------------------ 
    
    
    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine MakeConnections
    
        !Local---------------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer     :: PropertyX
        type(T_SPConnection), pointer :: ConnectionX
        integer                       :: STAT

        !Begin---------------------------------------------------------------------------------------------------------------------

        ConnectionX => Me%FirstConnection
        
        do while (associated(ConnectionX))            
            call SearchProperty(ConnectionX%ProductID, PropertyX, STAT)
            if (STAT .NE. SUCCESS_) &
                stop 'MakeConnections - ModuleChainReactions - ERR010'

            ConnectionX%Source%Product => PropertyX           
            ConnectionX => ConnectionX%Next      
        enddo

        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine MakeConnections
    !------------------------------------------------------------------------------------------------------------------------------ 
    
    
!    !------------------------------------------------------------------------------------------------------------------------------   
!    subroutine ConstructProduct(NewProperty, Name, IDNumber)
!
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_property), pointer :: NewProperty
!        character(LEN=*)          :: Name
!        integer                   :: IDNumber
!
!        !External------------------------------------------------------------------------------------------------------------------
!        integer                   :: STAT_CALL
!
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        allocate (NewProperty, STAT = STAT_CALL)            
!        if(STAT_CALL .NE. SUCCESS_) &
!            stop 'ConstructProduct - ModuleChainReactions - ERR010'
!        
!        call ResetProperty(NewProperty)
!        
!        NewProperty%IsOnlyProduct = .true.        
!        NewProperty%ID%Name       = Name
!        NewProperty%ID%IDNumber   = IDNumber
!        !--------------------------------------------------------------------------------------------------------------------------
!
!    end subroutine ConstructProduct    
!    !------------------------------------------------------------------------------------------------------------------------------    

    
    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine CheckProperties 
    
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
            Found = .false.
               
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
                stop 'CheckProperties - ModuleChainReactions - ERR010' 
            endif
            
            PropertyX => PropertyX%Next
        end do
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine CheckProperties
    !------------------------------------------------------------------------------------------------------------------------------ 

    
!    !------------------------------------------------------------------------------------------------------------------------------    
!    subroutine ConstructPartitionList 
!
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                     :: STAT_CALL
!        logical                     :: BlockFound
!        type(T_Partition), pointer  :: NewPartition
!
!        !Begin---------------------------------------------------------------------------------------------------------------------
!        do
!            call ExtractBlockFromBuffer(Me%ObjEnterData,                    &
!                                        ClientNumber    = Me%ClientNumber,  &
!                                        block_begin     = part_block_begin, &
!                                        block_end       = part_block_end,   &
!                                        BlockFound      = BlockFound,       &
!                                        STAT            = STAT_CALL)
!            if (STAT_CALL .EQ. SUCCESS_) then    
!
!                if (BlockFound) then                                                  
!                    
!                    !Construct a New Property 
!                    Call ConstructPartitionProperty(NewPartition)
!
!                    !Add new Property to the Properties List 
!                    Call AddPartition(NewPartition)                                        
!                    
!                else
!
!                    call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 
!
!                    if (STAT_CALL .NE. SUCCESS_) &
!                        stop 'ConstructPartitionList - ModuleChainReactions - ERR010'
!                    exit !No more blocks
!                
!                end if    
!
!            elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
!                
!                write(*,*)  
!                write(*,*) 'Error calling ExtractBlockFromBuffer. '
!                stop       'ConstructPartitionList - ModuleChainReactions - ERR020'
!            
!            else    
!                
!                stop 'ConstructPartitionList - ModuleChainReactions - ERR030'
!            
!            end if    
!        
!        end do            
!        !--------------------------------------------------------------------------------------------------------------------------
!                
!    end subroutine ConstructPartitionList
!    !------------------------------------------------------------------------------------------------------------------------------   


!    !------------------------------------------------------------------------------------------------------------------------------   
!    subroutine ConstructPartitionProperty (NewPartition)
!
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Partition), pointer :: NewPartition
!
!        !External------------------------------------------------------------------------------------------------------------------
!        integer                   :: STAT_CALL
!
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        allocate (NewPartition, STAT = STAT_CALL)            
!        if(STAT_CALL .NE. SUCCESS_) &
!            stop 'ConstructPartitionProperty - ModuleChainReactions - ERR010'
!        
!        nullify(NewPartition%Prev)
!        nullify(NewPartition%Next)
!        nullify(NewPartition%Pair1)
!        nullify(NewPartition%Pair2)
!       
!        NewPartition%PartitionEvolution = NoPartition
!        NewPartition%PartitionCoef      = 0.0
!        
!        call ConstructPartitionID (NewPartition%ID, Me%ObjEnterData, FromBlock)
!        call ConstructPartitionOptions (NewPartition)
!        call ConstructPartitionValues (NewPartition)
!        !--------------------------------------------------------------------------------------------------------------------------
!
!    end subroutine ConstructPartitionProperty    
!    !------------------------------------------------------------------------------------------------------------------------------    
!
!
!    !------------------------------------------------------------------------------------------------------------------------------    
!    subroutine ConstructPartitionID (PartitionID, ObjEnterData, ExtractType)
!
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type (T_PropertyID)                         :: PartitionID
!        integer                                     :: ObjEnterData
!        integer                                     :: ExtractType
!
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                                     :: flag
!        integer                                     :: STAT_CALL
!
!        !Partition Name
!        call GetData(PartitionID%Name, ObjEnterData, flag,  &
!                     SearchType   = ExtractType,            &
!                     keyword      = 'NAME',                 &
!                     Default      = '***',                  &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionID - ModuleChainReactions - ERR010'
!
!        !Units
!        call GetData(PartitionID%Units, ObjEnterData, flag, &
!                     SearchType   = ExtractType,            &
!                     keyword      = 'UNITS',                &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionID - ModuleChainReactions - ERR030'
!
!        !Description
!        call GetData(PartitionID%Description, ObjEnterData, flag,   &
!                     SearchType   = ExtractType,                    &
!                     keyword      = 'DESCRIPTION',                  &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'ConstructPartitionID - ModuleChainReactions - ERR040'
!        !--------------------------------------------------------------------------------------------------------------------------
!
!    end subroutine ConstructPartitionID
!    !------------------------------------------------------------------------------------------------------------------------------    
!
!    
!    !------------------------------------------------------------------------------------------------------------------------------        
!    subroutine ConstructPartitionOptions (NewPartition)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Partition), pointer :: NewPartition
!
!        !External------------------------------------------------------------------------------------------------------------------
!        integer                   :: STAT_CALL, iflag
!        character(MaxStrLength)   :: PairName
!        integer                   :: Phase
!
!        !Begin---------------------------------------------------------------------------------------------------------------------                             
!        call GetData(NewPartition%PartitionEvolution,       &   
!                     Me%ObjEnterData, iflag,                &
!                     SearchType   = FromBlock,              &
!                     keyword      = 'EVOLUTION',            &
!                     ClientModule = 'ModuleChainReactions', &
!                     STAT         = STAT_CALL)                     
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR010'              
!        if (iFlag .EQ. 0) then
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR020'  
!        endif
!            
!        call GetData(PairName,                              &   
!                     Me%ObjEnterData, iflag,                &
!                     SearchType   = FromBlock,              &
!                     keyword      = 'PAIR1',                &
!                     ClientModule = 'ModuleChainReactions', &
!                     STAT         = STAT_CALL)                   
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR030'              
!        if (iFlag .EQ. 0) then
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR040'  
!        endif
!        
!        call GetData(Phase,                                 &   
!                     Me%ObjEnterData, iflag,                &
!                     SearchType   = FromBlock,              &
!                     keyword      = 'PAIR1_PHASE',          &
!                     Default      = WaterPhase,             &
!                     ClientModule = 'ModuleChainReactions', &
!                     STAT         = STAT_CALL)                   
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR050'              
!        
!        call ConstructPairProperty(NewPartition%Pair1, PairName, Phase, STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR060'  
!            
!        call GetData(PairName,                              &   
!                     Me%ObjEnterData, iflag,                &
!                     SearchType   = FromBlock,              &
!                     keyword      = 'PAIR2',                &
!                     ClientModule = 'ModuleChainReactions', &
!                     STAT         = STAT_CALL)                   
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR070'              
!        if (iFlag .EQ. 0) then
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR080'  
!        endif
!        
!        call GetData(Phase,                                 &   
!                     Me%ObjEnterData, iflag,                &
!                     SearchType   = FromBlock,              &
!                     keyword      = 'PAIR2_PHASE',          &
!                     Default      = WaterPhase,             &
!                     ClientModule = 'ModuleChainReactions', &
!                     STAT         = STAT_CALL)                   
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR090'              
!        
!        call ConstructPairProperty(NewPartition%Pair2, PairName, Phase, STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionOptions - ModuleChainReactions - ERR100'   
!        !--------------------------------------------------------------------------------------------------------------------------
!        
!    end subroutine ConstructPartitionOptions
!    !------------------------------------------------------------------------------------------------------------------------------            
!
!
!    !------------------------------------------------------------------------------------------------------------------------------        
!    subroutine ConstructPairProperty (PropertyX, PropertyName, Phase, STAT)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Property), pointer :: PropertyX
!        character(MaxStrLength)   :: PropertyName
!        integer                   :: Phase
!        integer, optional         :: STAT
!        
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                   :: STAT_CALL, STAT_
!        integer                   :: Index
!        integer                   :: UB, LB
!        logical                   :: Found
!        integer                   :: PropertyID
!        type(T_property), pointer :: NewProperty
!        
!        !Begin---------------------------------------------------------------------------------------------------------------------        
!        STAT_ = UNKNOWN_
!        
!        LB = lbound(Me%Ext%Properties, 1)
!        UB = ubound(Me%Ext%Properties, 1)
!       
!        PropertyID = GetPropertyIDNumber (PropertyName) 
!        call SearchSinkProperty (PropertyID, PropertyX, STAT_CALL)
!        
!        if  (STAT_CALL .NE. SUCCESS_) then
!        
!            Found = .false.        
!            do Index = LB, UB
!                if (Me%Ext%Properties(Index) .EQ. PropertyID) then
!                    Found = .true.
!                    exit
!                endif                
!            end do
!       
!            if (.NOT. Found) then 
!                write(*,*)
!                write(*,*) 'Property ', trim(PropertyName), ' was not found in the ', trim(Me%CallingModule)
!                write(*,*) 'list of properties '
!                stop 'ConstructPairProperty - ModuleChainReactions - ERR010' 
!            endif
!            
!            allocate (NewProperty)
!            call ResetSinkProperty (NewProperty)
!            
!            NewProperty%ID%Name     = PropertyName
!            NewProperty%ID%IDNumber = PropertyID 
!            NewProperty%Phase       = Phase
!            
!            call AddSinkProperty (NewProperty)
!            
!            STAT_ = SUCCESS_
!            
!        else
!        
!            STAT_ = SUCCESS_
!            
!        endif
!        
!        if (present(STAT)) STAT = STAT_                
!        !--------------------------------------------------------------------------------------------------------------------------
!    
!    end subroutine ConstructPairProperty
!    !------------------------------------------------------------------------------------------------------------------------------        
!    
!
!    !------------------------------------------------------------------------------------------------------------------------------        
!    subroutine ConstructPartitionValues (NewPartition)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Partition), pointer          :: NewPartition
!        
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                             :: iflag
!        integer                             :: STAT_CALL
!
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        if (NewPartition%PartitionEvolution .EQ. ConstantPartition) then
!        
!            call GetData(NewPartition%PartitionCoef,            &   
!                         Me%ObjEnterData, iflag,                &
!                         SearchType   = FromBlock,              &
!                         keyword      = 'PARTITION_COEF',       &
!                         ClientModule = 'ModuleChainReactions', &
!                         STAT         = STAT_CALL)
!                         
!            if (STAT_CALL /= SUCCESS_) &
!                stop 'ConstructPartitionValues - ModuleChainReactions - ERR010'
!            
!            if (iFlag .EQ. 0) then
!                write(*,*)
!                write(*,*) 'You MUST provide a value for keyword PARTITION_COEF for '
!                write(*,*) 'property ', trim(NewPartition%ID%Name)
!                stop 'ConstructPartitionValues - ModuleChainReactions - ERR020'
!            endif
!
!        else
!               
!            if (Me%Options%GeometryType .EQ. Geometry2D) then
!                call ConstructPartitionValues2D (NewPartition)
!            else
!                call ConstructPartitionValues3D (NewPartition)
!            endif    
!                    
!        endif    
!        !--------------------------------------------------------------------------------------------------------------------------
!    
!    end subroutine ConstructPartitionValues
!    !------------------------------------------------------------------------------------------------------------------------------            
!
!
!    !------------------------------------------------------------------------------------------------------------------------------        
!    subroutine ConstructPartitionValues2D (NewPartition)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Partition), pointer          :: NewPartition
!
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                             :: STAT_CALL
!        integer                             :: ILB,IUB
!        integer                             :: JLB,JUB
!        integer                             :: WorkSizeILB, WorkSizeIUB
!        integer                             :: WorkSizeJLB, WorkSizeJUB       
!
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        !Boundaries
!        ILB = Me%Size%ILB
!        IUB = Me%Size%IUB
!        JLB = Me%Size%JLB
!        JUB = Me%Size%JUB
!
!        WorkSizeILB = Me%WorkSize%ILB
!        WorkSizeIUB = Me%WorkSize%IUB
!        WorkSizeJLB = Me%WorkSize%JLB
!        WorkSizeJUB = Me%WorkSize%JUB
!
!        allocate(NewPartition%G2D%PartitionCoef(ILB:IUB, JLB:JUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) &
!            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR010'
!        NewPartition%G2D%PartitionCoef(:,:) = 0.0
!    
!        !Get water points
!        call GetBasinPoints (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR020'
!    
!        call ConstructFillMatrix (PropertyID       = NewPartition%ID,                   &
!                                  EnterDataID      = Me%ObjEnterData,                   &
!                                  TimeID           = Me%ObjTime,                        &
!                                  HorizontalGridID = Me%ObjHorizontalGrid,              &
!                                  ExtractType      = FromBlock,                         &
!                                  PointsToFill2D   = Me%Ext%WaterPoints2D,              &
!                                  Matrix2D         = NewPartition%G2D%PartitionCoef,    &
!                                  TypeZUV          = TypeZ_,                            &
!                                  STAT             = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR030'
!
!        if(.NOT. NewPartition%ID%SolutionFromFile)then
!
!            call KillFillMatrix(NewPartition%ID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_)&
!                stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR040'
!                
!        end if  
!        
!        call UnGetBasin (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionValues2D - ModuleChainReactions - ERR050'                              
!        !--------------------------------------------------------------------------------------------------------------------------
!    
!    end subroutine ConstructPartitionValues2D
!    !------------------------------------------------------------------------------------------------------------------------------            
!
!
!    !------------------------------------------------------------------------------------------------------------------------------        
!    subroutine ConstructPartitionValues3D (NewPartition)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Partition), pointer :: NewPartition
!
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                   :: STAT_CALL
!        integer                   :: ILB,IUB
!        integer                   :: JLB,JUB
!        integer                   :: KLB,KUB        
!        integer                   :: WorkSizeILB, WorkSizeIUB
!        integer                   :: WorkSizeJLB, WorkSizeJUB       
!        integer                   :: WorkSizeKLB, WorkSizeKUB
!
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        !Boundaries
!        ILB = Me%Size%ILB
!        IUB = Me%Size%IUB
!        JLB = Me%Size%JLB
!        JUB = Me%Size%JUB
!        KLB = Me%Size%KLB
!        KUB = Me%Size%KUB
!
!        WorkSizeILB = Me%WorkSize%ILB
!        WorkSizeIUB = Me%WorkSize%IUB
!        WorkSizeJLB = Me%WorkSize%JLB
!        WorkSizeJUB = Me%WorkSize%JUB
!        WorkSizeKLB = Me%WorkSize%KLB
!        WorkSizeKUB = Me%WorkSize%KUB
!
!        allocate(NewPartition%G3D%PartitionCoef(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) &
!            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR010'
!        NewPartition%G3D%PartitionCoef(:,:,:) = 0.0
!    
!    
!        !Get water points
!        call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR020'
!                    
!        call ConstructFillMatrix (PropertyID       = NewPartition%ID,                   &
!                                  EnterDataID      = Me%ObjEnterData,                   &       
!                                  TimeID           = Me%ObjTime,                        &
!                                  HorizontalGridID = Me%ObjHorizontalGrid,              &       
!                                  GeometryID       = Me%ObjGeometry,                    &
!                                  ExtractType      = FromBlockInBlock,                  &
!                                  PointsToFill3D   = Me%Ext%WaterPoints3D,              &
!                                  Matrix3D         = NewPartition%G3D%PartitionCoef,    &
!                                  TypeZUV          = TypeZ_,                            &
!                                  STAT             = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR030'
!
!        if(.NOT. NewPartition%ID%SolutionFromFile)then
!
!            call KillFillMatrix(NewPartition%ID%ObjFillMatrix, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_)&
!                stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR040'
!                                          
!        end if    
!
!        call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) &
!            stop 'ConstructPartitionValues3D - ModuleChainReactions - ERR050'                
!        !--------------------------------------------------------------------------------------------------------------------------
!    
!    end subroutine ConstructPartitionValues3D
!    !------------------------------------------------------------------------------------------------------------------------------            
!        
!    
!    !------------------------------------------------------------------------------------------------------------------------------
!    subroutine AddPartition (NewPartition)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Partition), pointer :: NewPartition
!        
!        !Begin---------------------------------------------------------------------------------------------------------------------             
!        if (.NOT. associated(Me%FirstPartition)) then           
!            Me%FirstPartition => NewPartition
!            Me%LastPartition  => NewPartition
!        else
!            NewPartition%Prev     => Me%LastPartition
!            Me%LastPartition%Next => NewPartition
!            Me%LastPartition      => NewPartition
!        end if 
!        !--------------------------------------------------------------------------------------------------------------------------
!    
!    end subroutine AddPartition
!    !------------------------------------------------------------------------------------------------------------------------------
    
    
    !------------------------------------------------------------------------------------------------------------------------------ 
    subroutine ResetProperty (NewProperty)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: NewProperty
        
        !Begin---------------------------------------------------------------------------------------------------------------------    
        nullify(NewProperty%Prev)
        nullify(NewProperty%Next)
        nullify(NewProperty%Product)  
        
        NewProperty%Sink%Soil_0%Evolution = NoEvolution      
        NewProperty%Sink%Soil_1%Evolution = NoEvolution
        NewProperty%Sink%Water_0%Evolution = NoEvolution
        NewProperty%Sink%Water_1%Evolution = NoEvolution

!        NewProperty%IsOnlyProduct = .true.
!        NewProperty%RateOrder     = NoReaction        
!        NewProperty%RateMethod    = RateUnknown
        
!        NewProperty%SinkEvolution = NoSink
!        NewProperty%SinkRate      = 0.0
        
        call ResetCalculatedValues (NewProperty)      
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ResetProperty
    !------------------------------------------------------------------------------------------------------------------------------ 
    
    
    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ResetCalculatedValues (PropertyX)
    
        !Arguments-----------------------------------------------------------------------------------------------------------------
        type(T_Property), pointer :: PropertyX        
        
        !Begin---------------------------------------------------------------------------------------------------------------------
        PropertyX%Calc%OldMass      = 0.0
        PropertyX%Calc%MassLost     = 0.0
        PropertyX%Calc%MassGained   = 0.0
        PropertyX%Calc%OldSoilMass  = 0.0
        PropertyX%Calc%SoilMassLost = 0.0
        !--------------------------------------------------------------------------------------------------------------------------
    
    end subroutine ResetCalculatedValues
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
    subroutine SearchProperty (PropertyID, PropertyX, STAT)

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
        
    end subroutine SearchProperty

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
    subroutine InitCRSoilPhase (ChainReactionsID, SoilDensity2D, SoilDensity3D, STAT)

        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                                   :: ChainReactionsID
        real, dimension(:,:), pointer, optional   :: SoilDensity2D
        real, dimension(:,:,:), pointer, optional :: SoilDensity3D
        integer, optional, intent(OUT) :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer :: ready_                     
        integer :: STAT_CALL              !Auxiliar local variable
        integer :: I, J, K
        type(T_Property), pointer :: PropertyX
        
        !Begin---------------------------------------------------------------------------------------------------------------------
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)    
        
cd1:    if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            if ((Me%GeometryType == 1) .AND. (.NOT. present(SoilDensity2D))) then
                stop 'InitCRSoilPhase - ModuleChainReactions - ERR010'
            elseif ((Me%GeometryType == 2) .AND. (.NOT. present(SoilDensity3D))) then
                stop 'InitCRSoilPhase - ModuleChainReactions - ERR020'  
            endif

            call Read_Lock(mCHAINREACTIONS_, Me%InstanceID)

            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))
            
                if (PropertyX%Soil%Kd%Evolution /= NoEvolution) then                
                    if (Me%GeometryType == 1) then
                    
                        !Get water points
                        call GetBasinPoints (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'InitCRSoilPhase - ModuleChainReactions - ERR030'
                            
                        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                        
                            if (Me%Ext%WaterPoints2D(I, J) .EQ. WaterPoint) then
                                PropertyX%Soil%Mass%Values2D(I, J) = PropertyX%Concentration%Values2D(I, J) * &
                                                                     SoilDensity2D(I, J) * 0.001 * &
                                                                     PropertyX%Soil%Kd%Values2D(I, J) 
                            endif                
                        enddo
                        enddo
                        
                        call UnGetBasin (Me%ObjBasinGeometry, Me%Ext%WaterPoints2D, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'InitCRSoilPhase - ModuleChainReactions - ERR040'                    
                        
                    
                    else
                    
                        !Get water points
                        call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'InitCRSoilPhase - ModuleChainReactions - ERR050'
                    
                        do I = Me%WorkSize%ILB, Me%WorkSize%IUB
                        do J = Me%WorkSize%JLB, Me%WorkSize%JUB
                        do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                        
                            if (Me%Ext%WaterPoints3D(I, J, K) .EQ. WaterPoint) then
                                PropertyX%Soil%Mass%Values3D(I, J, K) = PropertyX%Concentration%Values3D(I, J, K) * &
                                                                        SoilDensity3D(I, J, K) * 0.001 * &
                                                                        PropertyX%Soil%Kd%Values3D(I, J, K) 
                            endif
                                                                                                  
                        enddo
                        enddo
                        enddo

                        call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) &
                            stop 'InitCRSoilPhase - ModuleChainReactions - ERR060'                    

                    endif
                endif
                
                PropertyX => PropertyX%Next
            enddo
            
            call Read_Unlock(mCHAINREACTIONS_, Me%InstanceID, "InitCRSoilPhase")            
            STAT_CALL = SUCCESS_
            
        else 
        
            STAT_CALL = ready_
            
        end if cd1

        if (present(STAT)) STAT = STAT_CALL
        !--------------------------------------------------------------------------------------------------------------------------

    end subroutine InitCRSoilPhase
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

            call SearchProperty(PropertyID, PropertyX, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'SetPropertyConcentrations2D - ModuleChainReactions - ERR010'
            
            PropertyX%Concentration%Values2D => Concentration

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

            call SearchProperty(PropertyID, PropertyX, STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) &
                stop 'SetPropertyConcentrations3D - ModuleChainReactions - ERR010'
            
            PropertyX%Concentration%Values3D => Concentration

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

!    subroutine ModifyChainReactions2D (ChainReactionsID, &
!                                       WaterVolume,      &
!                                       DT,               &
!                                       SoilMass,         &                                       
!                                       STAT)  
!
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        integer                                 :: ChainReactionsID
!        real, pointer, dimension(:,:)           :: WaterVolume      !L
!        real                                    :: DT
!        real, pointer, dimension(:,:), optional :: SoilMass         !kg
!        integer, optional,  intent(OUT)         :: STAT
!
!        !Local---------------------------------------------------------------------------------------------------------------------
!        integer                         :: STAT_CALL
!        integer                         :: I, J  
!        integer                         :: ILB,IUB
!        integer                         :: JLB,JUB
!        integer                         :: ready_
!        type(T_Property), pointer       :: PropertyX
!        type(T_Property), pointer       :: Product
!!        real                            :: ChangeOfMass
!        real                            :: MassLost_w0, MassLost_w1
!        real                            :: MassLost_s0, MassLost_s1
!        real                            :: SoilMassLost
!        real                            :: SolidMass
!        real                            :: NewMass
!
!        !Begin---------------------------------------------------------------------------------------------------------------------
!        STAT_CALL = UNKNOWN_
!
!        call Ready(ChainReactionsID, ready_)
!
!        if (ready_ .EQ. IDLE_ERR_) then
!
!            call ActualizePropertiesFromFile
!
!            Me%Ext%Theta2D => WaterContent
!            if (present(SoilDensity)) Me%Ext%SoilDensity2D => SoilDensity
!                        
!            !Get water points
!            call GetBasinPoints (Me%ObjMap, Me%Ext%WaterPoints2D, STAT = STAT_CALL)        
!            if (STAT_CALL /= SUCCESS_) &
!                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR010'
!            
!            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
!            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
!                                
!                if (Me%Ext%WaterPoints2D(I, J) .EQ. WaterPoint) then
!
!                    theta = WaterContent (I, J) 
!                    
!                    PropertyX => Me%FirstProperty                    
!                    do while (associated(PropertyX))                        
!                        
!                        if (PropertyX%Soil%Kd%Evolution /= NoEvolution) then
!                            mass = PropertyX%Concentration%Values3D(I, J) * theta +          &
!                                   PropertyX%Soil%Mass%Values3D(I, J, K)
!                        
!                            PropertyX%Concentration%Values3D(I, J, K) =                         &
!                                    mass /                                                      &      
!                                    (theta + SoilDensity(I, J, K) * 0.001 *                     &
!                                    PropertyX%Soil%Kd%Values3D(I, J, K))                                    
!
!                            PropertyX%Soil%Mass%Values3D(I, J, K)     =                         &
!                                    PropertyX%Concentration%Values3D(I, J, K) *                 &
!                                    PropertyX%Soil%Kd%Values3D(I, J, K) *                       & 
!                                    SoilDensity(I, J, K) * 0.001                                                            
!                        endif
!                                                
!                        PropertyX => PropertyX%Next                        
!                    enddo
!
!                    PropertyX => Me%FirstProperty                                                                              
!                    do while (associated(PropertyX))
!                    
!                        PropertyX%Calc%OldMass = PropertyX%Concentration%Values3D(I, J, K) * theta
!                                               
!                        if (PropertyX%Soil%Kd%Evolution .NE. NoEvolution) then
!                            if (.NOT. present(SoilDensity)) &                            
!                                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR020' 
!                                                           
!                            PropertyX%Calc%OldSoilMass = PropertyX%Soil%Mass%Values3D (I, J, K)
!                        else
!                            PropertyX%Calc%OldSoilMass = 0.0
!                        endif
!                                                
!                        if (associated(PropertyX%Product)) then
!                                                                                                                                                     
!                            Product => PropertyX%Product
!                                                        
!                            if (PropertyX%Sink%Water_1%Evolution /= NoEvolution) then
!                                MassLost_w1 = PropertyX%Calc%OldMass - PropertyX%Calc%OldMass * &
!                                              exp(-PropertyX%Sink%Water_1%Values3D(I, J, K) * DT / 86400)
!                            else
!                                MassLost_w1 = 0.0
!                            endif
!                                
!
!                            if (PropertyX%Sink%Soil_1%Evolution /= NoEvolution) then
!                                MassLost_s1 = PropertyX%Calc%OldSoilMass - PropertyX%Calc%OldSoilMass * &
!                                              exp(-PropertyX%Sink%Soil_1%Values3D(I, J, K) * DT / 86400)
!                            else
!                                MassLost_s1 = 0.0
!                            endif
!                                          
!                            PropertyX%Calc%SoilMassLost = min(MassLost_s1, PropertyX%Calc%OldSoilMass)                                                                                                                
!                            PropertyX%Calc%MassLost     = min(MassLost_w1, PropertyX%Calc%OldMass)
!                            
!                            Product%Calc%MassGained = Product%Calc%MassGained + PropertyX%Calc%MassLost + &
!                                                      PropertyX%Calc%SoilMassLost
!                        endif                                                
!                                                
!                        PropertyX => PropertyX%Next                                 
!                        
!                    enddo
!                        
!                    PropertyX => Me%FirstProperty
!                    
!                    do while (associated(PropertyX))                        
!                        
!                        if (PropertyX%Soil%Kd%Evolution /= NoEvolution) then
!                            PropertyX%Concentration%Values3D(I, J, K) =                         &
!                                    PropertyX%Concentration%Values3D(I, J, K) +                 &
!                                    (PropertyX%Calc%MassGained - (PropertyX%Calc%MassLost +     &
!                                    PropertyX%Calc%SoilMassLost)) /                             &      
!                                    (theta + SoilDensity(I, J, K) * 0.001 *                     &
!                                    PropertyX%Soil%Kd%Values3D(I, J, K))                                    
!
!                            PropertyX%Soil%Mass%Values3D(I, J, K)     =                         &
!                                    PropertyX%Concentration%Values3D(I, J, K) *                 &
!                                    PropertyX%Soil%Kd%Values3D(I, J, K) *                       & 
!                                    SoilDensity(I, J, K) * 0.001                                    
!                        else            
!                            PropertyX%Concentration%Values3D(I, J, K) =                         &
!                                    (PropertyX%Calc%OldMass + PropertyX%Calc%MassGained -       &
!                                    PropertyX%Calc%MassLost) / theta                          
!                        endif
!                                                       
!                        call ResetCalculatedValues(PropertyX)
!                                                
!                        PropertyX => PropertyX%Next                        
!                    enddo
!
!                endif 
!                                      
!            enddo
!            enddo
!            enddo
!                                      
!            call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
!            if (STAT_CALL /= SUCCESS_) &
!                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR030'                    
!                                                         
!            STAT_CALL = SUCCESS_
!            
!        else              
!         
!            STAT_CALL = ready_
!            
!        end if 
!
!        if (present(STAT)) STAT = STAT_CALL
!        !--------------------------------------------------------------------------------------------------------------------------
!
!    end subroutine ModifyChainReactions2D
!    !------------------------------------------------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------------------------------------------------
    subroutine ModifyChainReactions3D (ChainReactionsID, &
                                       WaterContent,     &
                                       DT,               &
                                       SoilDensity,      &
                                       STAT)  

        !Arguments-----------------------------------------------------------------------------------------------------------------
        integer                                     :: ChainReactionsID
        real, pointer, dimension(:,:,:)             :: WaterContent
        real                                        :: DT
        real, pointer, dimension(:,:,:), optional   :: SoilDensity
        integer, optional,  intent(OUT)             :: STAT

        !Local---------------------------------------------------------------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: I, J, K 
        integer                         :: ready_
        type(T_Property), pointer       :: PropertyX
        type(T_Property), pointer       :: Product
        
        real                            :: theta
!        real                            :: soil_density
!        real                            :: soil_kd
        
        real                            :: MassLost_w1
        real                            :: MassLost_s1
        
!        real                            :: soil_conc
        
        real                            :: mass

        !Begin---------------------------------------------------------------------------------------------------------------------             
        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call ActualizePropertiesFromFile

            Me%Ext%Theta3D => WaterContent
            if (present(SoilDensity)) Me%Ext%SoilDensity3D => SoilDensity
                        
            !Get water points
            call GetWaterPoints3D (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) &
                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR010'
            
            do I = Me%WorkSize%ILB, Me%WorkSize%IUB
            do J = Me%WorkSize%JLB, Me%WorkSize%JUB
            do K = Me%WorkSize%KLB, Me%WorkSize%KUB
                                
                if (Me%Ext%WaterPoints3D(I, J, K) .EQ. WaterPoint) then

                    theta = WaterContent (I, J, K) 
                    
                    PropertyX => Me%FirstProperty                    
                    do while (associated(PropertyX))                        
                        
                        if (PropertyX%Soil%Kd%Evolution /= NoEvolution) then
                            mass = PropertyX%Concentration%Values3D(I, J, K) * theta +          &
                                   PropertyX%Soil%Mass%Values3D(I, J, K)
                        
                            PropertyX%Concentration%Values3D(I, J, K) =                         &
                                    mass /                                                      &      
                                    (theta + SoilDensity(I, J, K) * 0.001 *                     &
                                    PropertyX%Soil%Kd%Values3D(I, J, K))                                    

                            PropertyX%Soil%Mass%Values3D(I, J, K)     =                         &
                                    PropertyX%Concentration%Values3D(I, J, K) *                 &
                                    PropertyX%Soil%Kd%Values3D(I, J, K) *                       & 
                                    SoilDensity(I, J, K) * 0.001                                                            
                        endif
                                                
                        PropertyX => PropertyX%Next                        
                    enddo

                    PropertyX => Me%FirstProperty                                                                              
                    do while (associated(PropertyX))
                    
                        PropertyX%Calc%OldMass = PropertyX%Concentration%Values3D(I, J, K) * theta
                                               
                        if (PropertyX%Soil%Kd%Evolution .NE. NoEvolution) then
                            if (.NOT. present(SoilDensity)) &                            
                                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR020' 
                                                           
                            PropertyX%Calc%OldSoilMass = PropertyX%Soil%Mass%Values3D (I, J, K)
                        else
                            PropertyX%Calc%OldSoilMass = 0.0
                        endif
                                                
                        if (associated(PropertyX%Product)) then

                            Product => PropertyX%Product
                                                        
                            if (PropertyX%Sink%Water_1%Evolution /= NoEvolution) then
                                MassLost_w1 = PropertyX%Calc%OldMass - PropertyX%Calc%OldMass * &
                                              exp(-PropertyX%Sink%Water_1%Values3D(I, J, K) * DT / 86400)
                            else
                                MassLost_w1 = 0.0
                            endif
                                

                            if (PropertyX%Sink%Soil_1%Evolution /= NoEvolution) then
                                MassLost_s1 = PropertyX%Calc%OldSoilMass - PropertyX%Calc%OldSoilMass * &
                                              exp(-PropertyX%Sink%Soil_1%Values3D(I, J, K) * DT / 86400)
                            else
                                MassLost_s1 = 0.0
                            endif
                                          
                            PropertyX%Calc%SoilMassLost = min(MassLost_s1, PropertyX%Calc%OldSoilMass)
                            PropertyX%Calc%MassLost     = min(MassLost_w1, PropertyX%Calc%OldMass)
                            
                            Product%Calc%MassGained = Product%Calc%MassGained + PropertyX%Calc%MassLost + &
                                                      PropertyX%Calc%SoilMassLost
                        endif                                                
                                                
                        PropertyX => PropertyX%Next                                 
                        
                    enddo
                        
                    PropertyX => Me%FirstProperty
                    
                    do while (associated(PropertyX))                        
                        
                        if (PropertyX%Soil%Kd%Evolution /= NoEvolution) then
                            PropertyX%Concentration%Values3D(I, J, K) =                         &
                                    PropertyX%Concentration%Values3D(I, J, K) +                 &
                                    (PropertyX%Calc%MassGained - (PropertyX%Calc%MassLost +     &
                                    PropertyX%Calc%SoilMassLost)) /                             &      
                                    (theta + SoilDensity(I, J, K) * 0.001 *                     &
                                    PropertyX%Soil%Kd%Values3D(I, J, K))                                    

                            PropertyX%Soil%Mass%Values3D(I, J, K)     =                         &
                                    PropertyX%Concentration%Values3D(I, J, K) *                 &
                                    PropertyX%Soil%Kd%Values3D(I, J, K) *                       & 
                                    SoilDensity(I, J, K) * 0.001                                    
                        else            
                            PropertyX%Concentration%Values3D(I, J, K) =                         &
                                    (PropertyX%Calc%OldMass + PropertyX%Calc%MassGained -       &
                                    PropertyX%Calc%MassLost) / theta                          
                        endif
                                                       
                        call ResetCalculatedValues(PropertyX)
                                                
                        PropertyX => PropertyX%Next                        
                    enddo

                endif 
                                      
            enddo
            enddo
            enddo
                                      
            call UnGetMap (Me%ObjMap, Me%Ext%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'ModifyChainReactions3D - ModuleChainReactions - ERR030'                    
                                                         
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
!        type(T_Partition), pointer  :: PartitionX
        integer                     :: STAT_CALL   
        
        !Begin---------------------------------------------------------------------------------------------------------------------    
        PropertyX => Me%FirstProperty

        do while (associated(PropertyX))

            if (PropertyX%Sink%Water_0%Evolution .EQ. VariableEvolution) then            
                if (Me%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,      &
                                           Matrix3D       = PropertyX%Sink%Water_0%Values3D, &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,            &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR010'
                else
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,      &
                                           Matrix2D       = PropertyX%Sink%Water_0%Values2D, &  
                                           PointsToFill2D = Me%Ext%WaterPoints2D,            &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR020'
                endif                            
            endif
            
            if (PropertyX%Sink%Water_1%Evolution .EQ. VariableEvolution) then            
                if (Me%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,      &
                                           Matrix3D       = PropertyX%Sink%Water_1%Values3D, &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,            &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR030'
                else
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,      &
                                           Matrix2D       = PropertyX%Sink%Water_1%Values2D, &  
                                           PointsToFill2D = Me%Ext%WaterPoints2D,            &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR040'
                endif                            
            endif
                        
            if (PropertyX%Sink%Soil_0%Evolution .EQ. VariableEvolution) then            
                if (Me%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,     &
                                           Matrix3D       = PropertyX%Sink%Soil_0%Values3D, &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,           &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR010'
                else
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,     &
                                           Matrix2D       = PropertyX%Sink%Soil_0%Values2D, &  
                                           PointsToFill2D = Me%Ext%WaterPoints2D,           &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR020'
                endif  
            endif          
                
            if (PropertyX%Sink%Soil_1%Evolution .EQ. VariableEvolution) then            
                if (Me%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,     &
                                           Matrix3D       = PropertyX%Sink%Soil_1%Values3D, &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,           &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR010'
                else
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix,     &
                                           Matrix2D       = PropertyX%Sink%Soil_1%Values2D, &  
                                           PointsToFill2D = Me%Ext%WaterPoints2D,           &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR020'
                endif                                   
            endif             
            
            if (PropertyX%Soil%Kd%Evolution .EQ. VariableEvolution) then            
                if (Me%GeometryType .EQ. Geometry3D) then
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                           Matrix3D       = PropertyX%Soil%Kd%Values3D, &
                                           PointsToFill3D = Me%Ext%WaterPoints3D,       &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR010'
                else
                    call ModifyFillMatrix (FillMatrixID   = PropertyX%ID%ObjFillMatrix, &
                                           Matrix2D       = PropertyX%Soil%Kd%Values2D, &  
                                           PointsToFill2D = Me%Ext%WaterPoints2D,       &
                                           STAT           = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ActualizePropertiesFromFile - ModuleChainReactions - ERR020'
                endif                                   
            endif        
                                   
            PropertyX => PropertyX%Next
            
        enddo             
        !--------------------------------------------------------------------------------------------------------------------------    
    
    end subroutine ActualizePropertiesFromFile    
    !------------------------------------------------------------------------------------------------------------------------------    


!    !------------------------------------------------------------------------------------------------------------------------------    
!    subroutine GetSinkRate (PropertyX, SinkRate, I, J, K)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!        type(T_Property), pointer   :: PropertyX
!        real                        :: SinkRate
!        integer                     :: I, J
!        integer, optional           :: K
!    
!        !Begin---------------------------------------------------------------------------------------------------------------------
!        if (PropertyX%SinkEvolution .EQ. ConstantSink) then
!        
!            SinkRate = PropertyX%SinkRate
!            
!        else
!        
!            if (present(K)) then
!            
!                SinkRate = PropertyX%G3D%SinkRate(I, J, K)
!            
!            else
!            
!                SinkRate = PropertyX%G2D%SinkRate(I, J)
!            
!            endif
!        
!        endif
!        !--------------------------------------------------------------------------------------------------------------------------    
!        
!    end subroutine GetSinkRate
!    !------------------------------------------------------------------------------------------------------------------------------    

        
!    !------------------------------------------------------------------------------------------------------------------------------    
!    subroutine GetPartitionCoef (PartitionX, PartitionCoef, I, J, K)
!    
!        !Arguments-----------------------------------------------------------------------------------------------------------------
!!        type(T_Partition), pointer  :: PartitionX
!        real                        :: PartitionCoef
!        integer                     :: I, J
!        integer, optional           :: K
!    
!        !Begin---------------------------------------------------------------------------------------------------------------------
!        if (PartitionX%PartitionEvolution .EQ. ConstantPartition) then
!        
!            PartitionCoef = PartitionX%PartitionCoef
!            
!        else
!        
!            if (present(K)) then
!            
!                PartitionCoef = PartitionX%G3D%PartitionCoef(I, J, K)
!            
!            else
!            
!                PartitionCoef = PartitionX%G2D%PartitionCoef(I, J)
!            
!            endif
!        
!        endif
!        !--------------------------------------------------------------------------------------------------------------------------    
!        
!    end subroutine GetPartitionCoef
    
    
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
!        type(T_Partition), pointer      :: PartitionX

        !------------------------------------------------------------------------                      

        STAT_CALL = UNKNOWN_

        call Ready(ChainReactionsID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mCHAINREACTIONS_,  Me%InstanceID)

            PropertyX => Me%FirstProperty
            
            do while (associated(PropertyX)) 
            
                if(PropertyX%ID%SolutionFromFile) then

                    call KillFillMatrix(PropertyX%ID%ObjFillMatrix, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)&
                        stop 'KillChainReactions - ModuleChainReactions - ERR010'
                end if
                
                PropertyX => PropertyX%Next
                
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