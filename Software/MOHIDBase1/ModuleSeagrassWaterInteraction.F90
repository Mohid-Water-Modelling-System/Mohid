!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : SeagrassWaterInteraction
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : September 2010
! REVISION      : Isabella Ascione Kenov - v1.0
! ---------------------------------------------------------------------------
! DESCRIPTION
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
!DT                : 120.

!NITROGEN            : 1          
!PHOSPHORUS          : 1

!<begin_SeagrassesRates>
!VMAXNH4W            : 0.0017   !maximum ammonia uptake rate by leaves(gN/gdw/day)
!VMAXNO3W            : 0.0017   !maximum nitrate uptake rate by leaves(gN/gdw/day)
!VMAXPO4W            : 0.21e-3  !maximum phosphate uptake rate by leaves(gP/gdw/day)
!KNO3W               : 0.12     !half saturation constant for nitrate uptake by leaves  (gN/m3)
!KNH4W               : 0.13     !half saturation constant for ammonia uptake by leaves(gN/m3)
!KPO4W               : 0.017    !half saturation constant for phosphate uptake by leaves  (gP/m3)
!<end_SeagrassesRates>

Module ModuleSeagrassWaterInteraction

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions
    use ModuleTime

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartSeagrassWaterInteraction
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ReadWaterInteractionParameters
    private ::      PropertyIndexNumber
    private ::      ConstructPropertyList

    !Selector
    public  :: GetSeagrassWaterInteractionPropertyList
    public  :: GetDTSeagrassWaterInteraction
    public  :: GetSeagrassWaterInteractionSize
    public  :: GetSeagrassWaterInteractionPropIndex
    public  :: GetSeagrassWaterInteractionRateFlux
    public  :: UnGetSeagrassWaterInteractionRateFlux
    public  :: UnGetSeagrassWaterInteraction

    !Modifier
    public  :: ModifySeagrassWaterInteraction
    private ::      ComputeWaterInteraction


    !Destructor
    public  :: KillSeagrassWaterInteraction                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjSeagrassWaterInteraction 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    type     T_External
    type(T_Time     )                       :: Now


        real,       pointer, dimension(:,:)         :: Mass
        integer,    pointer, dimension(:  )         :: OpenPoints
        real,       pointer, dimension(:)           :: Temperature
        real,       pointer, dimension(:)           :: Leaves
        real,       pointer, dimension(:)           :: NintFactor
        real,       pointer, dimension(:)           :: PintFactor
        real(8),       pointer, dimension(:)           :: WaterCellVol
        real,       pointer, dimension(:  )         :: SWRadiation
        real,       pointer, dimension(:  )         :: Thickness ! Thickness of the water column
        real,       pointer, dimension(:  )         :: Occupation
        real,       pointer, dimension(:  )         :: SWLightExctintionCoef

    end type T_External

    type      T_PropIndex
        integer                                     :: Ammonia              = null_int        
        integer                                     :: Phosphate            = null_int
        integer                                     :: Oxygen               = null_int
        integer                                     :: Leaves               = null_int
        integer                                     :: Nitrate              = null_int
    end type T_PropIndex

    type      T_RateIndex  
        integer                                     :: Phosphate            = null_int
        integer                                     :: Nitrate              = null_int
    end type T_RateIndex

    type T_StoredIndex
         integer                                    :: UptakeNH4NO3w        = 1
         integer                                    :: UptakePO4w           = 2
         integer                                    :: LightFactor            = 3
    end Type T_StoredIndex
    
    
    type     T_ComputeOptions
        logical                                     :: Nitrogen             = .false.
        logical                                     :: Phosphorus           = .false.
    end type T_ComputeOptions


   type     T_Parameters
        real                                        :: GrowthMaxRate        = null_real 
        real                                        :: PhotosyntRatioO2C    = null_real
        real                                        :: MinimumConcentration = null_real
        real                                        :: MinimumOxygen        = null_real
        real                                        :: VmaxNH4              = null_real
        real                                        :: VmaxPO4              = null_real
        real                                        :: VmaxNO3              = null_real
        real                                        :: KNH4                 = null_real
        real                                        :: KNO3                 = null_real
        real                                        :: KPO4                 = null_real
        real                                        :: Photoinhibition      = null_real
  end type     T_Parameters

    type       T_Seagrasses       
        integer                                     :: ObjTime              = 0 !Instance of ModuleTime
        integer                                     :: InstanceID
        type (T_Size1D)                             :: Prop
        type (T_Size1D)                             :: Size
        integer                                     :: ObjEnterData         = 0
        real                                        :: DT, DTDay
        character(len=StringLength)                 :: SedimentModel
        integer, dimension(:    ), pointer          :: PropertyList
        real,    dimension(:,:,:), pointer          :: Matrix
        real,    dimension(:,:)  , pointer          :: Rates
        type(T_PropIndex     )                      :: PropIndex
        type(T_External      )                      :: ExternalVar
        type(T_ComputeOptions)                      :: ComputeOptions
        type(T_Seagrasses    ), pointer             :: Next
        type(T_Parameters    )                      :: Parameters
        type(T_StoredIndex    )                      :: StoredIndex
    end type  T_Seagrasses

    !Global Module Variables
    type (T_Seagrasses), pointer                    :: FirstObjSeagrassWaterInteraction
    type (T_Seagrasses), pointer                    :: Me
    

    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartSeagrassWaterInteraction(ObjSeagrassWaterInteractionID, FileName, ILB, IUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjSeagrassWaterInteractionID 
        character(len=*)                                :: FileName
        integer          , intent(IN )                  :: ILB, IUB
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSeagrassWaterInterac_)) then
            nullify (FirstObjSeagrassWaterInteraction)
            call RegisterModule (mSeagrassWaterInterac_) 
        endif

        call Ready(ObjSeagrassWaterInteractionID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%Size%ILB = ILB
            Me%Size%IUB = IUB

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR01'

            call ReadData

           

            call PropertyIndexNumber
        
            call ConstructPropertyList

            call ConstructRates
 
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR02'


            !Returns ID
            ObjSeagrassWaterInteractionID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'StartSeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartSeagrassWaterInteraction
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Seagrasses), pointer                         :: NewObjSeagrassWaterInteraction
        type (T_Seagrasses), pointer                         :: PreviousObjSeagrassWaterInteraction


        !Allocates new instance
        allocate (NewObjSeagrassWaterInteraction)
        nullify  (NewObjSeagrassWaterInteraction%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSeagrassWaterInteraction)) then
            FirstObjSeagrassWaterInteraction         => NewObjSeagrassWaterInteraction
            Me                      => NewObjSeagrassWaterInteraction
        else
            PreviousObjSeagrassWaterInteraction      => FirstObjSeagrassWaterInteraction
            Me                      => FirstObjSeagrassWaterInteraction%Next
            do while (associated(Me))
                PreviousObjSeagrassWaterInteraction  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjSeagrassWaterInteraction
            PreviousObjSeagrassWaterInteraction%Next => NewObjSeagrassWaterInteraction
        endif

        Me%InstanceID = RegisterNewInstance (mSeagrassWaterInterac_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------
    
    
    subroutine ReadData
       
        !Local-----------------------------------------------------------------
        integer                                         :: ClientNumber, STAT_CALL
        logical                                         :: BlockFound
        
        !Begin-----------------------------------------------------------------

        call ConstructGlobalVariables
        

        
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                    ClientNumber    = ClientNumber,                 &
                                    block_begin     = '<begin_SeagrassesRates>',         &
                                    block_end       = '<end_SeagrassesRates>',           &
                                    BlockFound      = BlockFound,                   &
                                    STAT            = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then
            
            if (BlockFound) then

                call ReadWaterInteractionParameters

            else
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                write(*,*) 'Could not find <begin_SeagrassesRates>...<end_SeagrassesRates> block'
                stop 'ReadData - ModuleSeagrassWaterInteraction - ERR01'

            end if

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            stop       'ReadData - ModuleSeagrassWaterInteraction - ERR02'
        
        else
            stop       'ReadData - ModuleSeagrassWaterInteraction - ERR03'
        end if

        
    end subroutine ReadData

    !--------------------------------------------------------------------------


    subroutine ConstructGlobalVariables
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%DT,                                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DT',                                               &
                     Default      = 3600.,                                              &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassWaterInteraction - ERR01'

        Me%DTDay = Me%DT / (3600. * 24.)

   
                 
        call GetData(Me%ComputeOptions%Nitrogen,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NITROGEN',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassWaterInteraction - ERR03'


        call GetData(Me%ComputeOptions%Phosphorus,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHOSPHORUS',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassWaterInteraction - ERR04'

    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------

    subroutine ReadWaterInteractionParameters

        !Arguments-------------------------------------------------------------
        !type(T_Parameters)                              :: Parameters

        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

  

      
         call GetData(Me%Parameters%PhotosyntRatioO2C,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PHOTOSOC',                                         &
                     Default      = 32.0/12.0,                                          &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR35'
        
           
       !Seagrasses  maximum NH4 uptake rate  
        call GetData(Me%Parameters%VmaxNH4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'VMAXNH4W',                                &
                     Default      = 1.7e-3,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR60'
       
  !Seagrasses  maximum NO3 uptake rate  
        call GetData(Me%Parameters%VmaxNO3,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'VMAXNO3W',                                &
                     Default      = 1.07e-3,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR60'
       !Seagrasses  maximum PO4 uptake rate  
        call GetData(Me%Parameters%VmaxPO4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'VMAXPO4W',                                &
                     Default      = 2.1e-4,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR65'
     
        
     !Seagrasses  Half saturation const. for NH4 uptake   
        call GetData(Me%Parameters%KNH4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'KNH4W',                                &
                     Default      = 0.13,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR70'
     
      !Seagrasses  Half saturation const. for NO3 uptake   
        call GetData(Me%Parameters%KNO3,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'KNO3W',                                &
                     Default      = 0.12,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR70'
     
    !Seagrasses  Half saturation const. for PO4 uptake   
        call GetData(Me%Parameters%KPO4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'KPO4W',                                &
                     Default      = 0.017,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR75'

        call GetData(Me%Parameters%MinimumConcentration,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'LEAVES_MINCONC',                               &
                     Default      = 1.e-5,                                         &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR80'

        !Minimum oxygen concentration allowed in mg/l
        call GetData(Me%Parameters%MinimumOxygen,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MIN_OXYGEN',                                       &
                     Default      = 1.0e-5,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR85'


     call GetData(Me%Parameters%Photoinhibition,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PHOTOIN',                                          &
                     Default      = 100.0,                                               &
                     ClientModule = 'ModuleSeagrassWaterInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadWaterInteractionParameters - ModuleSeagrassWaterInteraction - ERR90'
   
   

    end subroutine ReadWaterInteractionParameters

    !--------------------------------------------------------------------------
   

    subroutine PropertyIndexNumber
        
        !Begin-----------------------------------------------------------------
        
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        !Oxygen
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB
        
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Leaves             = Me%Prop%IUB

        if(Me%ComputeOptions%Nitrogen)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB
        
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate        = Me%Prop%IUB
     

        end if

        if(Me%ComputeOptions%Phosphorus)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

      

        end if


    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyList

        !Begin-----------------------------------------------------------------

        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
       
        Me%PropertyList(Me%PropIndex%Oxygen)                = Oxygen_
        
        Me%PropertyList(Me%PropIndex%Leaves)                = SeagrassesLeaves_
        
        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%Ammonia)           = Ammonia_
            Me%PropertyList(Me%PropIndex%Nitrate)           = Nitrate_
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%Phosphate)         = Inorganic_Phosphorus_
        end if

    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------
    
    subroutine ConstructRates
     

        !Begin-----------------------------------------------------------------
        allocate(Me%Matrix  (Me%Size%ILB:Me%Size%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))

     !RUptakeNH4NO3w =1
     !RUptakePO4w =2
     !LightFactor   =3
     
     allocate(Me%Rates(Me%Size%ILB:Me%Size%IUB, 1:3))  !
    

    end subroutine ConstructRates


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine GetSeagrassWaterInteractionPropertyList(Life_ID, PropertyList, STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: Life_ID
        integer, dimension(:), pointer                          :: PropertyList
        integer, optional, intent(OUT)                          :: STAT

        !External--------------------------------------------------------------
        integer                                                 :: ready_              

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Life_ID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSeagrassWaterInterac_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetSeagrassWaterInteractionPropertyList

    !--------------------------------------------------------------------------
    
    subroutine GetDTSeagrassWaterInteraction(SeagrassWaterInteraction_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassWaterInteraction_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassWaterInteraction_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTSeagrassWaterInteraction
    
    !--------------------------------------------------------------------------

    !subroutine GetSeagrassesLeavesOptions(SeagrassesLeaves_ID, Salinity, STAT)

        !Arguments-------------------------------------------------------------
     !   integer                             :: SeagrassesLeaves_ID
     !   logical, optional, intent(OUT)      :: Salinity
      !  integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
     !   integer                             :: ready_              

        !Local-----------------------------------------------------------------
     !   integer                             :: STAT_
       
        !----------------------------------------------------------------------

     !   STAT_ = UNKNOWN_

     !   call Ready(SeagrassesLeaves_ID, ready_)    
        
!cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
!            (ready_ .EQ. READ_LOCK_ERR_)) then

 !           if (present(Salinity)) then
 !               if(Me%SeagrassesLeaves%SalinityEffect)then
 !                   Salinity = .true.
 !               else
 !                   Salinity = .false.
 !               end if
 !           end if

 !           STAT_ = SUCCESS_
 !       else 
 !           STAT_ = ready_
 !       end if cd1

 !       if (present(STAT))STAT = STAT_

 !   end subroutine GetSeagrassesLeavesOptions

    !--------------------------------------------------------------------------

    subroutine GetSeagrassWaterInteractionSize(SeagrassWaterInteraction_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassWaterInteraction_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassWaterInteraction_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetSeagrassWaterInteractionSize
    
    !--------------------------------------------------------------------------

    subroutine GetSeagrassWaterInteractionPropIndex (SeagrassWaterInteraction_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassWaterInteraction_ID
        integer,           intent(IN )      :: PropertyIDNumber
        integer,           intent(OUT)      :: PropertyIndex
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, CurrentIndex
        logical                             :: found
               
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassWaterInteraction_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            found = .false.
            do CurrentIndex = Me%Prop%ILB,Me%Prop%IUB

                if (PropertyIDNumber.eq. Me%PropertyList(CurrentIndex))then
                    PropertyIndex = CurrentIndex
                    found = .true.
                    exit
                end if
            
            end do

            if(.not. found)then
                STAT_ = NOT_FOUND_ERR_
            else
                STAT_ = SUCCESS_
            endif

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetSeagrassWaterInteractionPropIndex

    !--------------------------------------------------------------------------
    
   subroutine GetSeagrassWaterInteractionRateFlux(SeagrassWaterInteractionID, FirstProp, SecondProp,RateFlux, STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassWaterInteractionID
        integer,           intent(IN )      :: FirstProp
        integer,           intent(IN )      :: SecondProp
        real,    dimension(:), pointer      :: RateFlux
        integer, optional, intent(OUT)      :: STAT


        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

       call Ready(SeagrassWaterInteractionID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSeagrassWaterInterac_, Me%InstanceID)
      
            select case(FirstProp)
 
                case (LeavesUptakeN_)

                    RateFlux => Me%Rates    (:, Me%StoredIndex%UptakeNH4NO3w )
               
               case (LeavesUptakeP_)

                    RateFlux => Me%Rates    (:, Me%StoredIndex%UptakePO4w )
               
              case (LeavesLightFactor_)

                    RateFlux => Me%Rates    (:, Me%StoredIndex%LightFactor )  

       
                case default

                    RateFlux => Me%Matrix    (:, FirstProp, SecondProp      )

            end select
             

    

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetSeagrassWaterInteractionRateFlux

    !--------------------------------------------------------------------------
    
  
!-----------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------
    subroutine UnGetSeagrassWaterInteraction(SeagrassWaterInteractionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: SeagrassWaterInteractionID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassWaterInteractionID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSeagrassWaterInterac_, Me%InstanceID, "UnGetSeagrassWaterInteraction")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSeagrassWaterInteraction
    
    !--------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------------------

    subroutine UnGetSeagrassWaterInteractionRateFlux(SeagrassWaterInteractionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: SeagrassWaterInteractionID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassWaterInteractionID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSeagrassWaterInterac_, Me%InstanceID, "UnGetSeagrassWaterInteractionRateFlux")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSeagrassWaterInteractionRateFlux
    
    !-------------------------------------------------------------------------------
   


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifySeagrassWaterInteraction(ObjSeagrassWaterInteractionID, &
                                OpenPoints,           &
                                Mass,                 &
                                Temperature,          &
                                NintFactor,           &
                                PintFactor,           &
                                SWRadiation,          &
                                SWLightExctintionCoef,&
                                Thickness,            &
                                Occupation,           &
                                WaterCellVol,         &
                                STAT)
                                
 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSeagrassWaterInteractionID
  
        integer, dimension(:  ), pointer, optional  :: OpenPoints
        real,    dimension(:,:), pointer            :: Mass
        real,    dimension(:), pointer              :: NintFactor
        real,    dimension(:), pointer              :: PintFactor
        real(8),    dimension(:), pointer              :: WaterCellVol
        real,    dimension(:), pointer              :: Temperature
        real,    dimension(:  ), pointer, optional  :: SWRadiation
        real,    dimension(:  ), pointer, optional  :: SWLightExctintionCoef
        real,    dimension(:  ), pointer, optional  :: Thickness
        real,    dimension(:  ), pointer, optional  :: Occupation

        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: Index

        !Begin-----------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(ObjSeagrassWaterInteractionID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
                     
             Me%ExternalVar%Temperature  => Temperature
            if (.not. associated(Me%ExternalVar%Temperature))               &
                stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR00'
            
            end if
            
            Me%ExternalVar%Mass         => Mass
            if (.not. associated(Me%ExternalVar%Mass))                      &
                stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR01'
                
            if(present(OpenPoints))then
                Me%ExternalVar%OpenPoints  => OpenPoints
                if (.not. associated(Me%ExternalVar%OpenPoints)) &
                    stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR02'
           
            
             
             Me%ExternalVar%NintFactor         => NintFactor
            if (.not. associated(Me%ExternalVar%NintFactor))                      &
                stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR04'
           
           Me%ExternalVar%PintFactor         => PintFactor
            if (.not. associated(Me%ExternalVar%PintFactor))                      &
                stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR05'
                
       
                
            Me%ExternalVar%WaterCellVol         => WaterCellVol
            if (.not. associated(Me%ExternalVar%WaterCellVol))                      &
                stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR06'
            
          
          if(present(SWRadiation))then
                Me%ExternalVar%SWRadiation          => SWRadiation
                if (.not. associated(Me%ExternalVar%SWRadiation))           &
                    stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR07'
            end if

            if(present(SWLightExctintionCoef))then
                Me%ExternalVar%SWLightExctintionCoef=> SWLightExctintionCoef
                if (.not. associated(Me%ExternalVar%SWLightExctintionCoef)) &
                    stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR08'
            end if

            if(present(Thickness))then
                Me%ExternalVar%Thickness  => Thickness
                if (.not. associated(Me%ExternalVar%Thickness)) &
                    stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR09'
            end if
            
            
            if(present(Occupation))then
                Me%ExternalVar%Occupation  => Occupation
                if (.not. associated(Me%ExternalVar%Occupation)) &
                    stop 'ModifySeagrassWaterInteraction - ModuleSeagrassWaterInteraction - ERR09'
            end if            


            do Index = Me%Size%ILB, Me%Size%IUB

                !Reset rates matrix
                Me%Matrix(Index,:,:) = 0.

                call ComputeWaterInteraction(Index)
                
               
            enddo

 
            nullify(Me%ExternalVar%Mass                 )
            nullify(Me%ExternalVar%Temperature          )
            nullify(Me%ExternalVar%SWRadiation          )
            nullify(Me%ExternalVar%SWLightExctintionCoef)
            nullify(Me%ExternalVar%Thickness            )
            nullify(Me%ExternalVar%Occupation     )
            nullify(Me%ExternalVar%NintFactor           )
            nullify(Me%ExternalVar%PintFactor           )
            nullify(Me%ExternalVar%WaterCellVol         )

            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySeagrassWaterInteraction

    !--------------------------------------------------------------------------
    
    subroutine ComputeWaterInteraction(Index)

        !Arguments-------------------------------------------------------------
       
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------

        integer                                     :: NH4
        integer                                     :: PO4
        integer                                     :: NO3 
        integer                                     :: L 
        integer                                     :: O2
        real                                        :: UptakeNH4w
        real                                        :: UptakeNO3w
        real                                        :: UptakePO4w
        real                                        :: NintFactor
        real                                        :: PintFactor
        real, parameter                             :: minthickness = 0.001
        real                                        :: LightLimitingFactor        
        real                                        :: DZ1, DZ2, radiation_at_top_canopy

        !Description-------------------------------------------------------------
        ! This subroutine calculates the uptake of ammonia, nitrate and phosphate
        ! in the water column by seagrasses leaves. It also calculates the 
        ! light limitation  for seagrasses leaves. 
        ! After calculating the uptake, the concentrations of ammonia, 
        ! nitrate and phosphate in the water column are updated. 
        ! The ammonia uptake, the nitrate uptake 
        !  and light limitation  are 
        ! stored in the array Me%Rates. 
        !Begin-----------------------------------------------------------------

        NH4     = Me%PropIndex%Ammonia   ! Index of the property Ammonia
        PO4     = Me%PropIndex%Phosphate ! Index of the property Phosphate
        NO3     = Me%PropIndex%Nitrate   ! Index of the property Nitrate
        L       = Me%PropIndex%Leaves    ! Index of the property Seagrasses Leaves
        O2      = Me%PropIndex%Oxygen    ! Index of the property Oxygen

        !Light Limitation Factor (dimensionless)
         
        ! Thickness is the thickness of the water layer k 
        ! TopRadiation is the solar radiation on the top of the water layer k
        ! Thickness(i,j,k) = DWZ(i,j,k)= SZZ(i,j,k-1)-SZZ(i,j,k) (see module geometry)
        ! 
         
        ! DZ1 is the distance (m) between the top of the cell and the top of the canopy
        ! (minthickness is used to avoid division by 0 if DZ1 is 0)
        DZ1= max(minthickness, (1. - Me%ExternalVar%Occupation(index))*Me%ExternalVar%Thickness(index)) 
           
        ! (minthickness is used to avoid division by 0 if DZ2 is 0)
        ! DZ2 is seagrass height in the cell
        DZ2= max(minthickness,Me%ExternalVar%Occupation(index)*Me%ExternalVar%Thickness(index))
        
        if (DZ1 == minthickness) then
        ! the height of canopy reaches the top of the cell, so the radiation at top of cell is used
            radiation_at_top_canopy = Me%ExternalVar%SWRadiation(index)  
        else
            radiation_at_top_canopy = Me%ExternalVar%SWRadiation(index)*exp(-DZ1*Me%ExternalVar%SWLightExctintionCoef(index))
        end if
         
        !It is assumed that the light extinction coefficient is uniform in the cell (it is rough approximation)
         
        LightLimitingFactor =                                                                                      &
                          PhytoLightLimitationFactor(Thickness       = DZ2,                                         &
                                                     TopRadiation    = radiation_at_top_canopy,                     &
                                                     PExt            = Me%ExternalVar%SWLightExctintionCoef(index), &
                                                     Photoinhibition = Me%Parameters%Photoinhibition)
         
        
              
        
    
        ! NintFactor and PintFactor have been calculated 
        ! in the module BenthicEcology , 
        ! and redistributed over the water column as 3d arrays
        ! in the module WaterProperties
        
        ! NintFactor and PintFactor are dimensionless      
        NintFactor=max(0., Me%ExternalVar%NintFactor(Index))                          
        PintFactor=max(0., Me%ExternalVar%PintFactor(Index))   
        !
        ! The uptake of ammonia follows the Michaelis-Menten kinetics
        ! and depends on the internal nitrogen content of the plant(NintFactor).

        !
        !gN/m3/day
        UptakeNH4w = Me%Parameters%VmaxNH4                               *  &  ! gN/gdw/day *
                     NintFactor                                          *  &  ! []        *
                     Me%ExternalVar%Mass(L, Index)                       *  &  ! gdw/m3     *
                     Me%ExternalVar%Mass(NH4, Index)                     /  &  ! gN/m3     /
                     (Me%ExternalVar%Mass(NH4, Index)                    +  &  ! (gN/m3    +
                     Me%Parameters%KNH4)                                       ! +gN/m3)  

        !The uptake of nitrate follows the Michaelis-Menten kinetics
        ! and depends on the internal nitrogen content of the plant(NintFactor).
        !
        !gN/m3/day
        UptakeNO3w = Me%Parameters%VmaxNO3                               *  &  ! gN/gdw/day *
                     NintFactor                                          *  &  ! []        *
                     Me%ExternalVar%Mass(L, Index)                       *  &  ! gdw/m3     *
                     Me%ExternalVar%Mass(NO3, Index)                     /  &  ! gN/m3     /
                     (Me%ExternalVar%Mass(NO3, Index)                    +  &  ! (gN/m3    +
                     Me%Parameters%KNO3)                                       ! +gN/m3)   

        !The uptake of phosphate follows the Michaelis-Menten kinetics
        ! and depends on the internal nitrogen content of the plant(PintFactor).

    
       !gP/m3/day
        UptakePO4w = Me%Parameters%VmaxPO4                               *  &  ! gP/gdw/day *
                     PintFactor                                          *  &  ! []        *
                     Me%ExternalVar%Mass(L, Index)                       *  &  ! gdw/m3     *
                     Me%ExternalVar%Mass(PO4, Index)                     /  &  ! gP/m3     /
                     (Me%ExternalVar%Mass(PO4, Index)                    +  &  ! (gP/m3    +
                     Me%Parameters%KPO4)                                       ! +gP/m3)   

      !write(*,*), 'NO3 before  ', Me%ExternalVar%Mass(NO3, Index)
      ! write(*,*), 'NintFactor  ', NintFactor
      !  write(*,*), 'Leaves  ', Me%ExternalVar%Mass(L, Index) 

      !Ammonia is consumed by Seagrasses Leaves  
      !gN/m3
      Me%ExternalVar%Mass(NH4, Index)  = Me%ExternalVar%Mass(NH4, Index) -  &  ! gN/m3     -
                                         UptakeNH4w                      *  &  ! gN/m3/day *
                                         Me%DTDay                              ! day       *                        

      ! Nitrate is consumed by Seagrases Leaves 
      !gN/m3 
      Me%ExternalVar%Mass(NO3, Index)  = Me%ExternalVar%Mass(NO3, Index) -  &  ! gN/m3     -
                                         UptakeNO3w                      *  &  ! gN/m3/day *
                                         Me%DTDay                              ! day
     
      !Phosphateis consumed by Seagrasses Leaves 
      !gP/m3       
      Me%ExternalVar%Mass(PO4, Index)  = Me%ExternalVar%Mass(PO4, Index) -  &  ! gP/m3     -
                                         UptakePO4w                      *  &  ! gP/m3/day *
                                         Me%DTDay                              ! day
     
    
     !write(*,*), 'NO3 after  ', Me%ExternalVar%Mass(NO3, Index)
     
     ! Uptake rates and light  limitation function are stored in the array Me%Rates
     ! gN/day                                         =  gN/m3/day             * m3
     Me%Rates(Index, Me%StoredIndex%UptakeNH4NO3w   ) = (UptakeNH4w+UptakeNO3w)* Me%ExternalVar%WaterCellVol(index) 
     ! gP/day                                         =  gP/m3/day             * m3
     Me%Rates(Index, Me%StoredIndex%UptakePO4w      ) =  UptakePO4w            * Me%ExternalVar%WaterCellVol(index)  
     ! dimensionless
     Me%Rates(Index, Me%StoredIndex%LightFactor      ) = LightLimitingFactor     ! dimensionless
     
    ! write(*,*), Me%Rates(Index, Me%StoredIndex%UptakeNH4NO3w   )
   
    end subroutine ComputeWaterInteraction
    
   
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillSeagrassWaterInteraction(ObjSeagrassWaterInteractionID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjSeagrassWaterInteractionID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSeagrassWaterInteractionID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSeagrassWaterInterac_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%PropertyList)
                deallocate(Me%Matrix      )
                deallocate(Me%Rates       )

                !Deallocates Instance
                call DeallocateInstance ()


                ObjSeagrassWaterInteractionID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine KillSeagrassWaterInteraction
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Seagrasses), pointer          :: AuxObjSeagrassWaterInteraction
        type (T_Seagrasses), pointer          :: PreviousObjSeagrassWaterInteraction

        !Updates pointers
        if (Me%InstanceID == FirstObjSeagrassWaterInteraction%InstanceID) then
            FirstObjSeagrassWaterInteraction => FirstObjSeagrassWaterInteraction%Next
        else
            PreviousObjSeagrassWaterInteraction => FirstObjSeagrassWaterInteraction
            AuxObjSeagrassWaterInteraction      => FirstObjSeagrassWaterInteraction%Next
            do while (AuxObjSeagrassWaterInteraction%InstanceID /= Me%InstanceID)
                PreviousObjSeagrassWaterInteraction => AuxObjSeagrassWaterInteraction
                AuxObjSeagrassWaterInteraction      => AuxObjSeagrassWaterInteraction%Next
            enddo

            !Now update linked list
            PreviousObjSeagrassWaterInteraction%Next => AuxObjSeagrassWaterInteraction%Next

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

    subroutine Ready (ObjSeagrassWaterInteraction_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSeagrassWaterInteraction_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjSeagrassWaterInteraction_ID > 0) then
            call LocateObjSeagrassWaterInteraction (ObjSeagrassWaterInteraction_ID)
            ready_ = VerifyReadLock (mSeagrassWaterInterac_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjSeagrassWaterInteraction (ObjSeagrassWaterInteractionID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSeagrassWaterInteractionID

        !Local-----------------------------------------------------------------

        Me => FirstObjSeagrassWaterInteraction
        do while (associated (Me))
            if (Me%InstanceID == ObjSeagrassWaterInteractionID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleSeagrassWaterInteraction - LocateObjSeagrassWaterInteraction - ERR01'

    end subroutine LocateObjSeagrassWaterInteraction

    !--------------------------------------------------------------------------

end module ModuleSeagrassWaterInteraction

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------







