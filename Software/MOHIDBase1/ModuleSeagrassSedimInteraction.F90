!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : SeagrassSedimInteraction
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : September 2010
! REVISION      : Isabella Ascione Kenov - v1.0
!------------------------------------------------------------------------------
! DESCRIPTION   : Module to compute interactions
! between seagrasses roots and sediment.
! The biomass of roots is not calculated here, but in the module ModuleBenthicEcology
! Presently this module also includes subroutines to calculate the mineralization
! in the sediment since the module ModuleSedimentQuality (which calculates 
! sediment diagenesis)is not active in Mohid.
! further developments include the separation of the mineralization processes 
! from this module
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

!MIN_OXY             : 1.e-5        ! Minimum oxygen concentration (mg/l)
!PONDECAY            : 0.02          ! PON mineralization rate (1/day)
!POPDECAY            : 0.03          ! POP mineralization rate (1/day)
!NCRATIO             : 0.18          ! N/C ratio  (gN/gC)
!PCRATIO             : 0.024         ! P/C ratio  (gP/gC)
!GNKGDW              : 16.           ! gN to kg DW ratio in seagrasses
!GPKGDW              : 1.8           ! gP to kg DW ratio in seagrasses

!<begin_SeagrassesRates>
!VMAXNH4S            : 0.0017    !maximum ammonia uptake rate by roots(gN/gdw/day)
!VMAXPO4S            : 0.21e-3   !maximum phosphate uptake rate by roots (gP/gdw/day)
!KNH4S               : 0.9       !half saturation constant for ammonia uptake by roots (gN/m3)
!KPO4S               : 0.017     !half saturation constant for phosphate uptake by roots  (gP/m3)
!<end_SeagrassesRates>

Module ModuleSeagrassSedimInteraction

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions
    use ModuleTime

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartSeagrassSedimInteraction
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ReadSeagrassSedimInteractionParameters
    private ::          ReadMineralisationParameters
    private ::      PropertyIndexNumber
    private ::      ConstructPropertyList

    !Selector
    public  :: GetSeagrassSedimInteractionPropertyList
    public  :: GetDTSeagrassSedimInteraction
    public  :: GetSeagrassSedimInteractionSize
    public  :: GetSeagrassSedimInteractionPropIndex
    public  :: GetSeagrassSedimInteractionRateFlux
    public  :: UnGetSeagrassSedimInteraction
    public  :: UnGetSeagrassSedimInteractionRateFlux


    !Modifier
    public  :: ModifySeagrassSedimInteraction
    private ::      ComputeSeagrassSedimInteraction
    private ::      ComputePONMineralisation
    private ::      ComputePOPMineralisation

    !Destructor
    public  :: KillSeagrassSedimInteraction                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjSeagrassSedimInteraction 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    type     T_External
    type(T_Time     )                       :: Now


        real,       pointer, dimension(:,:)         :: Mass
        integer,    pointer, dimension(:  )         :: OpenPoints
        real,       pointer, dimension(:)           :: Temperature
        real,       pointer, dimension(:)           :: NintFactorR
        real,       pointer, dimension(:)           :: PintFactorR
        real,       pointer, dimension(:)           :: RootsMort
        real(8),       pointer, dimension(:)           :: SedimCellVol
    end type T_External

    type      T_PropIndex
        integer                                     :: Ammonia              = null_int        
        integer                                     :: Phosphate            = null_int
        integer                                     :: PON                  = null_int         
        integer                                     :: POP                  = null_int
        integer                                     :: Oxygen               = null_int
        integer                                     :: Nitrate              = null_int
        integer                                     :: Roots                = null_int
    end type T_PropIndex

    type      T_RateIndex     
        integer                                     :: Phosphate            = null_int
        integer                                     :: PON                  = null_int         
        integer                                     :: POP                  = null_int
        integer                                     :: Nitrate              = null_int
    end type T_RateIndex

    type T_StoredIndex
         integer                                    :: UptakeNH4s           = 1
         integer                                    :: UptakePO4s           = 2
    end Type T_StoredIndex
    
    
    type     T_ComputeOptions
        logical                                     :: Nitrogen             = .false.
        logical                                     :: Phosphorus           = .false.
    end type T_ComputeOptions


   type     T_Parameters
        real                                        :: MinimumOxygen        = null_real
        real                                        :: PONDecayRate         = null_real
        real                                        :: NC_Ratio             = null_real
        real                                        :: POPDecayRate         = null_real
        real                                        :: PC_Ratio             = null_real
        real                                        :: VmaxNH4              = null_real
        real                                        :: VmaxNO3              = null_real
        real                                        :: VmaxPO4              = null_real
        real                                        :: KNH4                 = null_real
        real                                        :: KPO4                 = null_real
        real                                        :: gNKgdw               = null_real
        real                                        :: gPKgdw               = null_real
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
        type(T_StoredIndex    )                     :: StoredIndex
    end type  T_Seagrasses

    !Global Module Variables
    type (T_Seagrasses), pointer                    :: FirstObjSeagrassSedimInteraction
    type (T_Seagrasses), pointer                    :: Me
    

    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartSeagrassSedimInteraction(ObjSeagrassSedimInteractionID, FileName, ILB, IUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjSeagrassSedimInteractionID 
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
        if (.not. ModuleIsRegistered(mSeagrassSedimInterac_)) then
            nullify (FirstObjSeagrassSedimInteraction)
            call RegisterModule (mSeagrassSedimInterac_) 
        endif

        call Ready(ObjSeagrassSedimInteractionID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%Size%ILB = ILB
            Me%Size%IUB = IUB

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR01'

            call ReadData

           

            call PropertyIndexNumber
        
            call ConstructPropertyList

            call ConstructRates
 
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR02'


            !Returns ID
            ObjSeagrassSedimInteractionID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'StartSeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartSeagrassSedimInteraction
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Seagrasses), pointer                         :: NewObjSeagrassSedimInteraction
        type (T_Seagrasses), pointer                         :: PreviousObjSeagrassSedimInteraction


        !Allocates new instance
        allocate (NewObjSeagrassSedimInteraction)
        nullify  (NewObjSeagrassSedimInteraction%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSeagrassSedimInteraction)) then
            FirstObjSeagrassSedimInteraction         => NewObjSeagrassSedimInteraction
            Me                      => NewObjSeagrassSedimInteraction
        else
            PreviousObjSeagrassSedimInteraction      => FirstObjSeagrassSedimInteraction
            Me                      => FirstObjSeagrassSedimInteraction%Next
            do while (associated(Me))
                PreviousObjSeagrassSedimInteraction  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjSeagrassSedimInteraction
            PreviousObjSeagrassSedimInteraction%Next => NewObjSeagrassSedimInteraction
        endif

        Me%InstanceID = RegisterNewInstance (mSeagrassSedimInterac_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------
    
    
    subroutine ReadData
       
        !Local-----------------------------------------------------------------
        integer                                         :: ClientNumber, STAT_CALL
        logical                                         :: BlockFound
        
        !Begin-----------------------------------------------------------------

        call ConstructGlobalVariables
        
        call ReadMineralisationParameters
        
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                    ClientNumber    = ClientNumber,                 &
                                    block_begin     = '<begin_SeagrassesRates>',         &
                                    block_end       = '<end_SeagrassesRates>',           &
                                    BlockFound      = BlockFound,                   &
                                    STAT            = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then
            
            if (BlockFound) then

                call ReadSeagrassSedimInteractionParameters

            else
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                write(*,*) 'Could not find <begin_SeagrassesRates>...<end_SeagrassesRates> block'
                stop 'ReadData - ModuleSeagrassSedimInteraction - ERR01'

            end if

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            stop       'ReadData - ModuleSeagrassSedimInteraction - ERR02'
        
        else
            stop       'ReadData - ModuleSeagrassSedimInteraction - ERR03'
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
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassSedimInteraction - ERR01'

        Me%DTDay = Me%DT / (3600. * 24.)

       ! call GetData(Me%SedimentModel,                                                   &
       !              Me%ObjEnterData, iflag,                                            &
       !              SearchType   = FromFile,                                           &
       !              keyword      = 'SEDIMENT_MODEL',                                    &
       !              ClientModule = 'SedimentQualityModel',                                 &
       !              STAT         = STAT_CALL)
       ! if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassSedimInteraction - ERR02'
       ! if(iflag==0)then
       !     write(*,*)'Please define the sediment model to couple with ModuleSeagrassSedimInteraction'
       !     stop 'ConstructGlobalVariables - ModuleSeagrassSedimInteraction - ERR20'
       ! end if

       ! if((Me%SedimentModel .ne. SedimentQualityModel))then
       !     write(*,*)'sediment model to couple with ModuleSeagrassSedimInteraction must be:'
       !     write(*,*)trim(SedimentQualityModel)
       !     stop 'ConstructGlobalVariables - ModuleSeagrassSedimInteraction - ERR30'
       ! endif
                 
        call GetData(Me%ComputeOptions%Nitrogen,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NITROGEN',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassSedimInteraction - ERR03'


        call GetData(Me%ComputeOptions%Phosphorus,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHOSPHORUS',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleSeagrassSedimInteraction - ERR04'

    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------

    subroutine ReadSeagrassSedimInteractionParameters

        !Arguments-------------------------------------------------------------
        !type(T_Parameters)                              :: Parameters

        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        
           
       !Seagrasses  maximum NH4 uptake rate  gN/gdw/day
        call GetData(Me%Parameters%VmaxNH4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'VMAXNH4S',                                &
                     Default      = 1.7e-3,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR60'
       

        
     !Seagrasses  Half saturation const. for NH4 uptake   gN/m3
        call GetData(Me%Parameters%KNH4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'KNH4S',                                &
                     Default      = 0.9,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR70'
        
        
          !Seagrasses  maximum PO4 uptake rate  gP/gdw/day
        call GetData(Me%Parameters%VmaxPO4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'VMAXPO4S',                                &
                     Default      = 0.21e-3,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR72'
       

        
     !Seagrasses  Half saturation const. for PO4 uptake   gP/m3
        call GetData(Me%Parameters%KPO4,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'KPO4S',                                &
                     Default      = 0.017,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR74'
   

        !Minimum oxygen concentration allowed in mg/l (or g/m3, it is the same)
        call GetData(Me%Parameters%MinimumOxygen,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MIN_OXYGEN',                                       &
                     Default      = 1e-5,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR85'
        
        !Ratio gN to KgDW
        call GetData(Me%Parameters%gNKgdw,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'GNKGDW',                                &
                     Default      = 16.,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR89'
   
  !Ratio gP to KgDW
        call GetData(Me%Parameters%gPKgdw,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'GPKGDW',                                &
                     Default      = 1.8,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSeagrassSedimInteractionParameters - ModuleSeagrassSedimInteraction - ERR98'

    end subroutine ReadSeagrassSedimInteractionParameters

    !--------------------------------------------------------------------------
    
    
    subroutine ReadMineralisationParameters
    
       !Arguments-------------------------------------------------------------
        !type(T_Parameters)                              :: Parameters

        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
    
    
    call GetData(Me%Parameters%MinimumOxygen,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'MIN_OXY',                                       &
                     Default      = 1.e-5,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
     if(STAT_CALL .NE. SUCCESS_) stop 'ReadMineralisationParameters - ModuleSeagrassSedimInteraction - ERR01'
            !PON decayRate
    
    
     call GetData(Me%Parameters%PONDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'PONDECAY',                                       &
                     Default      = 0.01,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
    if(STAT_CALL .NE. SUCCESS_) stop 'ReadMineralisationParameters - ModuleSeagrassSedimInteraction - ERR02'
    
    
     !POP decayRate
     call GetData(Me%Parameters%POPDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'POPDECAY',                                       &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
    if(STAT_CALL .NE. SUCCESS_) stop 'ReadMineralisationParameters - ModuleSeagrassSedimInteraction - ERR03'   
    
                !NC Ratio
     call GetData(Me%Parameters%NC_Ratio,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'NCRATIO',                                       &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
    if(STAT_CALL .NE. SUCCESS_) stop 'ReadMineralisationParameters - ModuleSeagrassSedimInteraction - ERR04'   
        
     call GetData(Me%Parameters%PC_Ratio,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                          &
                     keyword      = 'PCRATIO',                                       &
                     Default      = 0.024,                                               &
                     ClientModule = 'ModuleSeagrassSedimInteraction',                                 &
                     STAT         = STAT_CALL)
    if(STAT_CALL .NE. SUCCESS_) stop 'ReadMineralisationParameters - ModuleSeagrassSedimInteraction - ERR05'   
           
        
        
    end subroutine ReadMineralisationParameters

    subroutine PropertyIndexNumber
        
        !Begin-----------------------------------------------------------------
        
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        !Oxygen
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB
        
        ! Roots
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Roots              = Me%Prop%IUB

        if(Me%ComputeOptions%Nitrogen)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB
        
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate        = Me%Prop%IUB
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%PON            = Me%Prop%IUB

        end if

        if(Me%ComputeOptions%Phosphorus)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POP            = Me%Prop%IUB

        end if


    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyList

        !Begin-----------------------------------------------------------------

        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
       
        Me%PropertyList(Me%PropIndex%Oxygen)                = Oxygen_
        Me%PropertyList(Me%PropIndex%Roots)                = SeagrassesRoots_

        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%Ammonia)           = Ammonia_
            Me%PropertyList(Me%PropIndex%Nitrate)           = Nitrate_
            Me%PropertyList(Me%PropIndex%PON)               = PON_
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%Phosphate)         = Inorganic_Phosphorus_
            Me%PropertyList(Me%PropIndex%POP)               = POP_
        end if

    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------
    
    subroutine ConstructRates
     

        !Begin-----------------------------------------------------------------
        allocate(Me%Matrix  (Me%Size%ILB:Me%Size%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))

     !UptakeNH4s     : 1
     !UptakePO4s     : 2

  
        allocate(Me%Rates(Me%Size%ILB:Me%Size%IUB, 1:2))  ! ATTENZIONE
    

    end subroutine ConstructRates


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine GetSeagrassSedimInteractionPropertyList(Life_ID, PropertyList, STAT)

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

            call Read_Lock(mSeagrassSedimInterac_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetSeagrassSedimInteractionPropertyList

    !--------------------------------------------------------------------------
    
    subroutine GetDTSeagrassSedimInteraction(SeagrassSedimInteraction_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassSedimInteraction_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassSedimInteraction_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTSeagrassSedimInteraction
    
  
    !--------------------------------------------------------------------------

    subroutine GetSeagrassSedimInteractionSize(SeagrassSedimInteraction_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassSedimInteraction_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassSedimInteraction_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetSeagrassSedimInteractionSize
    
    !--------------------------------------------------------------------------

    subroutine GetSeagrassSedimInteractionPropIndex (SeagrassSedimInteraction_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassSedimInteraction_ID
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

        call Ready(SeagrassSedimInteraction_ID, ready_)    
        
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

    end subroutine GetSeagrassSedimInteractionPropIndex

    !--------------------------------------------------------------------------
    
   subroutine GetSeagrassSedimInteractionRateFlux(SeagrassSedimInteractionID, FirstProp, SecondProp, RateFlux, STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: SeagrassSedimInteractionID
        integer,           intent(IN )      :: FirstProp, SecondProp
        real,    dimension(:), pointer      :: RateFlux
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

       call Ready(SeagrassSedimInteractionID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSeagrassSedimInterac_, Me%InstanceID)

            select case(FirstProp)
            
            case(RootsUptakeN_)
                
                RateFlux => Me%Rates    (:, Me%StoredIndex%UptakeNH4s ) 
            
            case(RootsUptakeP_)
                
                RateFlux => Me%Rates    (:, Me%StoredIndex%UptakePO4s ) 
                
     !           case(SGGrossProd_)

    !                RateFlux => Me%Parameters(:, SecondProp, SG_GrossFact_  )

     !           case(SGTemperatureLim_)
                   
    !                RateFlux => Me%Parameters(:, SecondProp, SG_TLimFact_   )

   !             case(SGNutrientLim_)

    !                RateFlux => Me%Parameters(:, SecondProp, SG_NutLimFact_ )


    !            case(SGSalinityLim_)

    !                RateFlux => Me%Parameters(:, SecondProp, SG_SLimFact_   )

    !            case(SGNLim_)

    !                RateFlux => Me%Parameters(:, SecondProp, SG_NLimFact_   )

    !            case(SGPLim_)

   !                 RateFlux => Me%Parameters(:, SecondProp, SG_PLimFact_   )


               
                case default

                    RateFlux => Me%Matrix    (:, FirstProp, SecondProp      )

            end select

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetSeagrassSedimInteractionRateFlux

    !--------------------------------------------------------------------------
    
  

    !--------------------------------------------------------------------------------


    subroutine UnGetSeagrassSedimInteraction(SeagrassSedimInteractionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: SeagrassSedimInteractionID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassSedimInteractionID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSeagrassSedimInterac_, Me%InstanceID, "UnGetSeagrassSedimInteraction")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSeagrassSedimInteraction
    
    !--------------------------------------------------------------------------

    subroutine UnGetSeagrassSedimInteractionRateFlux(SeagrassSedimInteractionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: SeagrassSedimInteractionID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(SeagrassSedimInteractionID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSeagrassSedimInterac_, Me%InstanceID, "UnGetSeagrassSedimInteractionRateFlux")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSeagrassSedimInteractionRateFlux

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifySeagrassSedimInteraction(ObjSeagrassSedimInteractionID, &
                                OpenPoints,  &
                                Mass,        &
                                NintFactorR, &
                                RootsMort,   &
                                PintFactorR, &
                                SedimCellVol,STAT)
                                
 
        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSeagrassSedimInteractionID
  
        integer, dimension(:  ), pointer, optional  :: OpenPoints
        real,    dimension(:,:), pointer            :: Mass
        real,    dimension(:), pointer              :: NintFactorR
        real,    dimension(:), pointer              :: PintFactorR
        real(8),    dimension(:), pointer              :: SedimCellVol
        real,    dimension(:), pointer              :: RootsMort
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: Index

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSeagrassSedimInteractionID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
                     
            Me%ExternalVar%Mass         => Mass
            if (.not. associated(Me%ExternalVar%Mass))                      &
                stop 'ModifySeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR01'
                
            if(present(OpenPoints))then
                Me%ExternalVar%OpenPoints  => OpenPoints
                if (.not. associated(Me%ExternalVar%OpenPoints)) &
                    stop 'ModifySeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR02'
            
             
             Me%ExternalVar%NintFactorR         => NintFactorR
            if (.not. associated(Me%ExternalVar%NintFactorR))                      &
                stop 'ModifySeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR03'
                
                
           Me%ExternalVar%PintFactorR         => PintFactorR
            if (.not. associated(Me%ExternalVar%PintFactorR))                      &
                stop 'ModifySeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR04'
           
           
           Me%ExternalVar%RootsMort         => RootsMort
            if (.not. associated(Me%ExternalVar%RootsMort))                      &
                stop 'ModifySeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR05'
                
            Me%ExternalVar%SedimCellVol         => SedimCellVol
            if (.not. associated(Me%ExternalVar%SedimCellVol))                      &
                stop 'ModifySeagrassSedimInteraction - ModuleSeagrassSedimInteraction - ERR06'
            
            
         
         
            end if

            do Index = Me%Size%ILB, Me%Size%IUB

                !Reset rates matrix
                Me%Matrix(Index,:,:) = 0.

                call ComputeSeagrassSedimInteraction(Index)
                
                call ComputePONMineralisation(Index)
               if (Me%ComputeOptions%Phosphorus) then
                call ComputePOPMineralisation(Index)
              endif
            enddo

 
           ! nullify(Me%ExternalVar%Mass                 )


            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySeagrassSedimInteraction

    !--------------------------------------------------------------------------
    
    subroutine ComputeSeagrassSedimInteraction(Index)

        !Arguments-------------------------------------------------------------
       ! type(T_Parameters)                          :: Parameters
        
        
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: NH4, PO4, PON, POP, R
        real                                        :: UptakeNH4s
        real                                        :: UptakePO4s
        real                                        :: NintFactorR, PintFactorR
        real                                        :: MortalityN, MortalityP

        ! Description
        ! This subroutine calculates the uptake of ammonia in the sediment 
        ! by seagrasses roots. After the uptake, the ammonia 
        ! concentration is updated in the sediment column. 
        ! Dead roots re-integrate the existing pool of organic matter 
        ! in the sediment.
        ! The uptake of ammonia is stored in the array Me%Rates for further use
        ! by other modules. This module is not calculating the sources-sinks terms
        ! for roots. The module ModuleBenthicEcology calculates the sources-sinks terms
        ! for roots.
        !Begin-----------------------------------------------------------------

        NH4     = Me%PropIndex%Ammonia  ! Index of the property Ammonia
        PO4     = Me%PropIndex%Phosphate  ! Index of the property Phosphate
        PON     = Me%PropIndex%PON      ! Index of the property Particulate Organic Nitrogen
        POP     = Me%PropIndex%POP      ! Index of the property Particulate Organic Phosphorus
        R       = Me%PropIndex%Roots    ! Index of the property Roots
 
            
             
        ! uptake limitation due to internal nitrogen content
        ! (dimensionless)   
        NintFactorR=Me%ExternalVar%NintFactorR(index) 
        
        PintFactorR=Me%ExternalVar%PintFactorR(index) 
        !  gN/m3/day =        KgDW/day             * gN/kgDW                   / m3
        MortalityN=Me%ExternalVar%RootsMort(index) * (Me%Parameters%gNKgdw) /Me%ExternalVar%SedimCellVol(index) 
        
        !  gP/m3/day =        KgDW/day             * gP/kgDW                   / m3
        MortalityP=Me%ExternalVar%RootsMort(index) * (Me%Parameters%gPKgdw) /Me%ExternalVar%SedimCellVol(index) 
        
        !
        !  uptake limitation due to 
        ! internal nitrogen content (Me%ExternalVar%NintFactor) was 
        ! calculated in the module ModuleBenthicEcology
        ! 
      
        !gN/m3/day
        UptakeNH4s = Me%Parameters%VmaxNH4                               *  &  ! gN/gdw/day *
                     NintFactorR                                         *  &  ! []        *
                     Me%ExternalVar%Mass(R, Index)                       *  &  ! gdw/m3    *
                     Me%ExternalVar%Mass(NH4, Index)                     /  &  ! gN/m3    /
                     (Me%ExternalVar%Mass(NH4, Index)                    +  &  ! (gN/m3   +
                     Me%Parameters%KNH4)                                       ! +gN/m3)  *

   ! 
      ! Ammonia is consumed by roots uptake:
      ! gN/m3   
      Me%ExternalVar%Mass(NH4, Index)  = Me%ExternalVar%Mass(NH4, Index) -  & ! gN/m3      - 
                                         UptakeNH4s                      *  & ! gN/m3/day  *
                                         Me%DTDay                             ! day



     ! The uptake of ammonia is stored in the array Me%Rates for further use
     ! by the module ModuleBenthicEcology
     ! Me%ExternalVar%SedimCellVol(index) is the volume of the sediment cell
     ! gN/day                                      = gN/m3/day * m3
     Me%Rates(Index, Me%StoredIndex%UptakeNH4s   ) = UptakeNH4s*Me%ExternalVar%SedimCellVol(index)  !
     
     
     
     !gP/m3/day
        UptakePO4s = Me%Parameters%VmaxPO4                               *  &  ! gP/gdw/day *
                     PintFactorR                                         *  &  ! []        *
                     Me%ExternalVar%Mass(R, Index)                       *  &  ! gdw/m3    *
                     Me%ExternalVar%Mass(PO4, Index)                     /  &  ! gP/m3    /
                     (Me%ExternalVar%Mass(PO4, Index)                    +  &  ! (gP/m3   +
                     Me%Parameters%KPO4)                                       ! +gP/m3)  *

   ! 
      ! Phosphate is consumed by roots uptake:
      ! gN/m3   
      Me%ExternalVar%Mass(PO4, Index)  = Me%ExternalVar%Mass(PO4, Index) -  & ! gN/m3      - 
                                         UptakePO4s                      *  & ! gN/m3/day  *
                                         Me%DTDay                             ! day



     ! The uptake of ammonia is stored in the array Me%Rates for further use
     ! by the module ModuleBenthicEcology
     ! Me%ExternalVar%SedimCellVol(index) is the volume of the sediment cell
     ! gN/day                                      = gN/m3/day * m3
     Me%Rates(Index, Me%StoredIndex%UptakeNH4s   ) = UptakeNH4s*Me%ExternalVar%SedimCellVol(index)  !
    
     ! The uptake of phosphate is stored in the array Me%Rates for further use
     ! by the module ModuleBenthicEcology
     ! Me%ExternalVar%SedimCellVol(index) is the volume of the sediment cell
     ! gP/day                                      = gP/m3/day * m3
     Me%Rates(Index, Me%StoredIndex%UptakePO4s   ) = UptakePO4s*Me%ExternalVar%SedimCellVol(index)  !
     
     
    ! Mortality of roots is added to PON and POP in the sediment
     
     Me%ExternalVar%Mass(PON, Index)=Me%ExternalVar%Mass(PON, Index) + MortalityN * Me%DTDay 
     Me%ExternalVar%Mass(POP, Index)=Me%ExternalVar%Mass(POP, Index) + MortalityP * Me%DTDay 
   
    end subroutine ComputeSeagrassSedimInteraction
    
    
    
    subroutine ComputePONMineralisation(index)
    
    !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: AM, PON, O2
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation
        real                                        :: AnaerobicMineralisation
        real                                        :: AerobicMineralisation

        !Begin-----------------------------------------------------------------

        AM  = Me%PropIndex%Ammonia
        PON = Me%PropIndex%PON
        O2  = Me%PropIndex%Oxygen
        
        
        
        !oxygen units are given in mg/l
        OxygenLimitation = max(Me%ExternalVar%Mass(O2, Index), Me%Parameters%MinimumOxygen)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Parameters%PONDecayRate
                            


       if (Me%ExternalVar%Mass(O2, Index) <=Me%Parameters%MinimumOxygen) then
       ! Anaerobic mineralisation occurs at a rate that is the 30% of the normal rate (CAEDYM model)
       AnaerobicMineralisation = 0.3 * MineralizationRate*Me%ExternalVar%Mass(PON, Index)* Me%DTDay
       
       else
       AnaerobicMineralisation = 0.
       end if
       
       AerobicMineralisation = MineralizationRate * OxygenLimitation*Me%ExternalVar%Mass(PON, Index)* Me%DTDay
       
        !gN * day * day-1 (what passes from PON to ammonia)
        Me%Matrix(Index, PON, AM) =  AerobicMineralisation+ AnaerobicMineralisation 
                                    
         
       Me%ExternalVar%Mass(AM,  Index) = Me%ExternalVar%Mass(AM , Index) + Me%Matrix(Index, PON, AM) 

       Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) - Me%Matrix(Index, PON, AM) 

            !what is consumed of oxygen due to mineralization of PON
      OxygenConsumption               = AerobicMineralisation * 1. / Me%Parameters%NC_Ratio * &
                                        32. / 12.

      Me%ExternalVar%Mass(O2, Index ) = Me%ExternalVar%Mass(O2, Index ) - OxygenConsumption   
                                    
    end subroutine ComputePONMineralisation
    
    
    subroutine ComputePOPMineralisation(index)
    
    
        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: IP, POP, O2
        real                                        :: MineralizationRate
        real                                        :: AnaerobicMineralisation
        real                                        :: AerobicMineralisation
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------
        
        IP  = Me%PropIndex%Phosphate
        POP = Me%PropIndex%POP
        O2  = Me%PropIndex%Oxygen
        
            !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max(Me%ExternalVar%Mass(O2,Index), Me%Parameters%MinimumOxygen)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Parameters%POPDecayRate  
                         

             if (Me%ExternalVar%Mass(O2, Index) <=Me%Parameters%MinimumOxygen) then
       ! Anaerobic mineralisation occurs at a rate that is the 30% of the normal rate (CAEDYM model)
       AnaerobicMineralisation = 0.3 * MineralizationRate*Me%ExternalVar%Mass(POP, Index)* Me%DTDay
       
       else
       AnaerobicMineralisation = 0.
       end if
       
       AerobicMineralisation = MineralizationRate * OxygenLimitation*Me%ExternalVar%Mass(POP, Index)* Me%DTDay
        !kgP * day * day-1 (what passes from POP to inorganic phosphorus)
        
        
        Me%Matrix(Index, POP, IP) = AnaerobicMineralisation + AerobicMineralisation
        
                         
        Me%ExternalVar%Mass(IP,  Index) = Me%ExternalVar%Mass(IP , Index) + Me%Matrix(Index, POP, IP)

        Me%ExternalVar%Mass(POP, Index) = Me%ExternalVar%Mass(POP, Index) - Me%Matrix(Index, POP, IP)

        OxygenConsumption               = AerobicMineralisation * 1. / Me%Parameters%PC_Ratio * &
                                              32. / 12.
        
       Me%ExternalVar%Mass(O2, Index ) = Me%ExternalVar%Mass(O2, Index ) - OxygenConsumption
    
    end subroutine ComputePOPMineralisation
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillSeagrassSedimInteraction(ObjSeagrassSedimInteractionID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjSeagrassSedimInteractionID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSeagrassSedimInteractionID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSeagrassSedimInterac_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%PropertyList)
                deallocate(Me%Matrix      )
                deallocate(Me%Rates       )

                !Deallocates Instance
                call DeallocateInstance ()


                ObjSeagrassSedimInteractionID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine KillSeagrassSedimInteraction
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Seagrasses), pointer          :: AuxObjSeagrassSedimInteraction
        type (T_Seagrasses), pointer          :: PreviousObjSeagrassSedimInteraction

        !Updates pointers
        if (Me%InstanceID == FirstObjSeagrassSedimInteraction%InstanceID) then
            FirstObjSeagrassSedimInteraction => FirstObjSeagrassSedimInteraction%Next
        else
            PreviousObjSeagrassSedimInteraction => FirstObjSeagrassSedimInteraction
            AuxObjSeagrassSedimInteraction      => FirstObjSeagrassSedimInteraction%Next
            do while (AuxObjSeagrassSedimInteraction%InstanceID /= Me%InstanceID)
                PreviousObjSeagrassSedimInteraction => AuxObjSeagrassSedimInteraction
                AuxObjSeagrassSedimInteraction      => AuxObjSeagrassSedimInteraction%Next
            enddo

            !Now update linked list
            PreviousObjSeagrassSedimInteraction%Next => AuxObjSeagrassSedimInteraction%Next

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

    subroutine Ready (ObjSeagrassSedimInteraction_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSeagrassSedimInteraction_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjSeagrassSedimInteraction_ID > 0) then
            call LocateObjSeagrassSedimInteraction (ObjSeagrassSedimInteraction_ID)
            ready_ = VerifyReadLock (mSeagrassSedimInterac_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjSeagrassSedimInteraction (ObjSeagrassSedimInteractionID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSeagrassSedimInteractionID

        !Local-----------------------------------------------------------------

        Me => FirstObjSeagrassSedimInteraction
        do while (associated (Me))
            if (Me%InstanceID == ObjSeagrassSedimInteractionID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleSeagrassSedimInteraction - LocateObjSeagrassSedimInteraction - ERR01'

    end subroutine LocateObjSeagrassSedimInteraction

    !--------------------------------------------------------------------------

end module ModuleSeagrassSedimInteraction

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------







