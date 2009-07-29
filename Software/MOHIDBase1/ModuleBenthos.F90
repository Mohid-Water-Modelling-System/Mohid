!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Benthos
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2005
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute simple benthic biogeochemical processes
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
!DataFile
!   DT                          : real           [3600]         !Time step compute biogeochemical processes



Module ModuleBenthos

    use ModuleGlobalData
    use ModuleEnterData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartBenthos
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ReadOrganicMatterParameters
    private ::          ReadNitrogenParameters
    private ::          ReadPhosphorusParameters
    private ::          ReadSilicaParameters
    private ::          ReadPhytoParameters
    private ::          ReadBactParameters
    private ::          ReadDiatomsParameters
    private ::      PropertyIndexNumber
    private ::      ConstructPropertyList

    !Selector
    public  :: GetBenthosPropertyList
    public  :: GetDTBenthos
    public  :: GetBenthosSize
    public  :: GetBenthosPropIndex
    public  :: GetBenthosRateFlux
    public  :: UnGetBenthos
    public  :: UnGetBenthosRateFlux

    !Modifier
    public  :: ModifyBenthos
    private ::      ComputeBenthicPhyto
    private ::      ComputeBenthicBacteria
    private ::      ComputeBenthicDiatoms
    private ::      ComputeBenthicNitrogen
    private ::      ComputeBenthicPhosphorus
    private ::      ComputeBenthicSilica

    !Destructor
    public  :: KillBenthos                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjBenthos 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    type      T_External
        real,       pointer, dimension(:  )         :: Oxygen
        real,       pointer, dimension(:  )         :: Temperature
        real,       pointer, dimension(:,:)         :: Mass      
    end type T_External

    type       T_PropIndex
        integer                                     :: Ammonia          = null_int        
        integer                                     :: Phosphate        = null_int
        integer                                     :: DissolvedSilica  = null_int
        integer                                     :: BioSilica        = null_int
        integer                                     :: PON              = null_int         
        integer                                     :: POP              = null_int
        integer                                     :: Oxygen           = null_int
        integer                                     :: Phyto            = null_int
        integer                                     :: Bacteria         = null_int
        integer                                     :: PONr             = null_int
        integer                                     :: DONnr            = null_int
        integer                                     :: Diatoms          = null_int
    end type T_PropIndex

    type     T_OrganicMatter
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
    end type T_OrganicMatter

    type     T_Nitrogen
        real                                        :: PONDecayRate     = null_real 
        real                                        :: PONDecayTFactor  = null_real
    end type T_Nitrogen
    
    type     T_Phosphorus
        real                                        :: POPDecayRate     = null_real 
        real                                        :: POPDecayTFactor  = null_real
    end type T_Phosphorus
    
    type     T_Oxygen
        real                                        :: Minimum          = null_real 
    end type T_Oxygen

    type     T_Silica
        real                                        :: BioSiDecayRate   = null_real 
    end type T_Silica
    
    type     T_Phyto
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
        real                                        :: MortalityRate    = null_real 
    end type T_Phyto

    type     T_Bact
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
        real                                        :: MortalityRate    = null_real
        real                                        :: BacMaxUptake     = null_real
        real                                        :: NSatConstBac     = null_real
        real                                        :: BotReactorDepth  = null_real
        real                                        :: BacMinSub        = null_real
        real                                        :: TOptBacteriaMin  = null_real
        real                                        :: TOptBacteriaMax  = null_real
        real                                        :: TBacteriaMin     = null_real
        real                                        :: TBacteriaMax     = null_real
        real                                        :: BK1              = null_real
        real                                        :: BK2              = null_real
        real                                        :: BK3              = null_real
        real                                        :: BK4              = null_real
    end type T_Bact

    type     T_Diatoms
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
        real                                        :: SiC_Ratio        = null_real
        real                                        :: MortalityRate    = null_real 
    end type T_Diatoms

    type     T_ComputeOptions
        logical                                     :: Nitrogen         = .false.
        logical                                     :: Phosphorus       = .false.
        logical                                     :: Silica           = .false.
        logical                                     :: Diatoms          = .false.
        logical                                     :: Phyto            = .false.
        logical                                     :: Bacteria         = .false.
    end type T_ComputeOptions

    type       T_Benthos
        integer                                     :: InstanceID
        type (T_Size1D)                             :: Prop
        type (T_Size1D)                             :: Size
        integer                                     :: ObjEnterData     = 0
        real                                        :: DT, DTDay
        character(len=StringLength)                 :: PelagicModel
        integer, dimension(:    ),  pointer         :: PropertyList
        real,    dimension(:,:,:),  pointer         :: Matrix
        type(T_PropIndex     )                      :: PropIndex
        type(T_Nitrogen      )                      :: Nitrogen
        type(T_Phosphorus    )                      :: Phosphorus
        type(T_Silica        )                      :: Silica
        type(T_Diatoms       )                      :: Diatoms
        type(T_Phyto         )                      :: Phyto
        type(T_Bact          )                      :: Bact
        type(T_Oxygen        )                      :: Oxygen
        type(T_OrganicMatter )                      :: OrganicMatter
        type(T_External      )                      :: ExternalVar
        type(T_ComputeOptions)                      :: ComputeOptions
        type(T_Benthos       ),     pointer         :: Next
    end type  T_Benthos

    !Global Module Variables
    type (T_Benthos), pointer                       :: FirstObjBenthos
    type (T_Benthos), pointer                       :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartBenthos(ObjBenthosID, FileName, ILB, IUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjBenthosID 
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
        if (.not. ModuleIsRegistered(mBenthos_)) then
            nullify (FirstObjBenthos)
            call RegisterModule (mBenthos_) 
        endif

        call Ready(ObjBenthosID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%Size%ILB = ILB
            Me%Size%IUB = IUB

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartBenthos - ModuleBenthos - ERR01'

            call ReadData

            call PropertyIndexNumber
        
            call ConstructPropertyList

            call ConstructRates
 
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartBenthos - ModuleBenthos - ERR02'


            !Returns ID
            ObjBenthosID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'StartBenthos - ModuleBenthos - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartBenthos
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Benthos), pointer                         :: NewObjBenthos
        type (T_Benthos), pointer                         :: PreviousObjBenthos


        !Allocates new instance
        allocate (NewObjBenthos)
        nullify  (NewObjBenthos%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjBenthos)) then
            FirstObjBenthos         => NewObjBenthos
            Me                      => NewObjBenthos
        else
            PreviousObjBenthos      => FirstObjBenthos
            Me                      => FirstObjBenthos%Next
            do while (associated(Me))
                PreviousObjBenthos  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjBenthos
            PreviousObjBenthos%Next => NewObjBenthos
        endif

        Me%InstanceID = RegisterNewInstance (mBenthos_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------
    
    
    subroutine ReadData

        call ConstructGlobalVariables

        call ReadOrganicMatterParameters

        call ReadNitrogenParameters

        call ReadPhosphorusParameters
        
        call ReadSilicaParameters

        call ReadPhytoParameters
        
        call ReadBactParameters
       
        call ReadDiatomsParameters

        call ReadOxygenParameters

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
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR01'

        Me%DTDay = Me%DT / (3600. * 24.)

        call GetData(Me%PelagicModel,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PELAGIC_MODEL',                                    &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR10'
        if(iflag==0)then
            write(*,*)'Please define the pelagic model to couple with ModuleBenthos'
            stop 'ConstructGlobalVariables - ModuleBenthos - ERR20'
        end if

        if((Me%PelagicModel .ne. WaterQualityModel .and. Me%PelagicModel .ne. LifeModel))then
            write(*,*)'Pelagic model to couple with ModuleBenthos must be one of the following:'
            write(*,*)trim(WaterQualityModel)
            write(*,*)trim(LifeModel)
            stop 'ConstructGlobalVariables - ModuleBenthos - ERR30'
        endif
                 
        call GetData(Me%ComputeOptions%Nitrogen,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NITROGEN',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR40'


        call GetData(Me%ComputeOptions%Phosphorus,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHOSPHORUS',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR50'


        call GetData(Me%ComputeOptions%Silica,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SILICA',                                           &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR60'
      
       
        call GetData(Me%ComputeOptions%Diatoms,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DIATOMS',                                          &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR80'
       
        call GetData(Me%ComputeOptions%Phyto,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO',                                            &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR90'

        call GetData(Me%ComputeOptions%Bacteria,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BACT',                                             &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthos - ERR100'

    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------

    subroutine ReadOrganicMatterParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%OrganicMatter%NC_Ratio,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NC_RATIO',                                         &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOrganicMatterParameters - ModuleBenthos - ERR01'

        call GetData(Me%OrganicMatter%PC_Ratio,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PC_RATIO',                                         &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOrganicMatterParameters - ModuleBenthos - ERR10'

    end subroutine ReadOrganicMatterParameters

    !--------------------------------------------------------------------------
   
    subroutine ReadNitrogenParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Nitrogen%PONDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PON_DECAY_RATE',                                   &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadNitrogenParameters - ModuleBenthos - ERR01'
        
        call GetData(Me%Nitrogen%PONDecayTFactor,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PON_DECAY_TFACTOR',                                &
                     Default      = 1.02,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadNitrogenParameters - ModuleBenthos - ERR02'

    end subroutine ReadNitrogenParameters

    !--------------------------------------------------------------------------

    subroutine ReadPhosphorusParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Phosphorus%POPDecayRate,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POP_DECAY_RATE',                                   &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhosphorusParameters - ModuleBenthos - ERR01'

        call GetData(Me%Phosphorus%POPDecayTFactor,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POP_DECAY_TFACTOR',                                &
                     Default      = 1.08,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhosphorusParameters - ModuleBenthos - ERR02'


    end subroutine ReadPhosphorusParameters

    !--------------------------------------------------------------------------

    subroutine ReadSilicaParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Silica%BioSiDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BIOSI_DECAY_RATE',                                 &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSilicaParameters - ModuleBenthos - ERR01'


    end subroutine ReadSilicaParameters


    !--------------------------------------------------------------------------

    subroutine ReadPhytoParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetData(Me%Phyto%MortalityRate,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO_MORTALITY',                                  &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhytoParameters - ModuleBenthos - ERR01'

        call GetData(Me%Phyto%NC_Ratio,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO_NC_RATIO',                                   &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhytoParameters - ModuleBenthos - ERR10'


        call GetData(Me%Phyto%PC_Ratio,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO_PC_RATIO',                                   &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhytoParameters - ModuleBenthos - ERR20'


    end subroutine ReadPhytoParameters

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadBactParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetData(Me%Bact%MortalityRate,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BACT_MORTALITY',                                   &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR01'

        call GetData(Me%Bact%NC_Ratio,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BACT_NC_RATIO',                                    &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR10'


        call GetData(Me%Bact%PC_Ratio,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BACT_PC_RATIO',                                    &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR20'

        call GetData(Me%Bact%BacMaxUptake,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BMAXUPTA',                                         &
                     Default      = 0.20,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR30'

        call GetData(Me%Bact%NSatConstBac,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BACNCONS',                                         &
                     Default      = 0.08,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR40'

        call GetData(Me%Bact%BotReactorDepth,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BOTTOMRDEPTH',                                     &
                     Default      = 0.3,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR50'

        call GetData(Me%Bact%BacMinSub,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BACMINSUB',                                        &
                     Default      = 0.01,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR60'
 
        !TOptBacteriaMin, minimum temperature of the optimal interval for the Bacteria growth, 
        call GetData(Me%Bact%TOptBacteriaMin,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword = 'TOPTBMIN',                                              &
                     default      = 24.8,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR70' 

        !TOptBacteriaMax, maximum temperature of the optimal interval for the Bacteria growth, 
        call GetData(Me%Bact%TOptBacteriaMax,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TOPTBMAX',                                         &
                     default      = 25.1,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR80' 

        !TBacteriaMin, minimum tolerable temperature of the  interval for the Bacteria growth, 
        call GetData(Me%Bact%TBacteriaMin,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TBMIN',                                            &
                     default      = 5.0,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR90' 

        !TBacteriaMax, maximum tolerable temperature of the  interval for the Bacteria growth, 
        call GetData(Me%Bact%TBacteriaMax,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TBMAX',                                            &
                     default      = 35.0,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR100' 


        !BK1, constant to control temperature response curve shape
        call GetData(Me%Bact%BK1,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TBCONST1',                                         &
                     default      = 0.05,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR110' 

        !BK2, constant to control temperature response curve shape
        call GetData(Me%Bact%BK2,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TBCONST2',                                         &
                     default      = 0.98,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR120' 


        !BK3, constant to control temperature response curve shape
        call GetData(Me%Bact%BK3,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TBCONST3',                                         &
                     default      = 0.98,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR130' 


        !BK4, constant to control temperature response curve shape
        call GetData(Me%Bact%BK4,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TBCONST4',                                         &
                     default      = 0.02,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadBactParameters - ModuleBenthos - ERR140' 
    end subroutine ReadBactParameters

    !--------------------------------------------------------------------------

    subroutine ReadDiatomsParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetData(Me%Diatoms%MortalityRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DIATOMS_MORTALITY',                                &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadDiatomsParameters - ModuleBenthos - ERR01'


        call GetData(Me%Diatoms%SiC_Ratio,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DIATOMS_SIC_RATIO',                                &
                     Default      = 0.07,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadDiatomsParameters - ModuleBenthos - ERR10'


        call GetData(Me%Diatoms%NC_Ratio,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DIATOMS_NC_RATIO',                                 &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadDiatomsParameters - ModuleBenthos - ERR20'


        call GetData(Me%Diatoms%PC_Ratio,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DIATOMS_PC_RATIO',                                 &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadDiatomsParameters - ModuleBenthos - ERR30'

    end subroutine ReadDiatomsParameters

    !--------------------------------------------------------------------------

    subroutine ReadOxygenParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Oxygen%Minimum,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MIN_OXYGEN',                                       &
                     Default      = 1e-5,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOxygenParameters - ModuleBenthos - ERR01'


    end subroutine ReadOxygenParameters


    !--------------------------------------------------------------------------

    subroutine PropertyIndexNumber
        
        !Begin-----------------------------------------------------------------
        
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0
        
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB

        if(Me%ComputeOptions%Nitrogen)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB
        
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%PON            = Me%Prop%IUB

        end if

        if(Me%ComputeOptions%Phosphorus)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POP            = Me%Prop%IUB

        end if

        if(Me%ComputeOptions%Silica)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DissolvedSilica= Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%BioSilica      = Me%Prop%IUB

        end if

        if(Me%ComputeOptions%Phyto)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phyto          = Me%Prop%IUB

        end if

        if(Me%ComputeOptions%Bacteria)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Bacteria       = Me%Prop%IUB

            if(Me%ComputeOptions%Nitrogen) then            
            
                Me%Prop%IUB             = Me%Prop%IUB + 1
                Me%PropIndex%PONr       = Me%Prop%IUB
                
                Me%Prop%IUB             = Me%Prop%IUB + 1
                Me%PropIndex%DONnr      = Me%Prop%IUB

            end if
        end if

        if(Me%ComputeOptions%Diatoms)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Diatoms        = Me%Prop%IUB

        end if


    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyList

        !Begin-----------------------------------------------------------------

        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
        Me%PropertyList(Me%PropIndex%Oxygen)                = Oxygen_

        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%Ammonia)           = Ammonia_
            Me%PropertyList(Me%PropIndex%PON)               = PON_
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%Phosphate)         = Inorganic_Phosphorus_
            Me%PropertyList(Me%PropIndex%POP)               = POP_
        end if

        if(Me%ComputeOptions%Silica)then

            if(Me%PelagicModel .eq. WaterQualityModel) then
                Me%PropertyList(Me%PropIndex%DissolvedSilica)   = DSilica_
            endif

            if(Me%PelagicModel .eq. LifeModel) then            
                Me%PropertyList(Me%PropIndex%DissolvedSilica)   = Silicate_
            endif

            Me%PropertyList(Me%PropIndex%BioSilica)         = BioSilica_
        
        end if

        if(Me%ComputeOptions%Phyto)then
            Me%PropertyList(Me%PropIndex%Phyto)             = Phytoplankton_
        end if

        if(Me%ComputeOptions%Bacteria)then
            
            Me%PropertyList(Me%PropIndex%Bacteria)          = Bacteria_
            
            if(Me%ComputeOptions%Nitrogen)then
                Me%PropertyList(Me%PropIndex%PONr)          = PONRefractory_
                Me%PropertyList(Me%PropIndex%DONnr)         = DONNon_Refractory_
            end if

        end if

        if(Me%ComputeOptions%Diatoms)then
            Me%PropertyList(Me%PropIndex%Diatoms)           = Diatoms_
        end if


    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------
    
    subroutine ConstructRates

        !Begin-----------------------------------------------------------------

        allocate(Me%Matrix(Me%Size%ILB:Me%Size%IUB,                           &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))
    
    end subroutine ConstructRates


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine GetBenthosPropertyList(Life_ID, PropertyList, STAT)

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

            call Read_Lock(mBenthos_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthosPropertyList

    !--------------------------------------------------------------------------
    
    subroutine GetDTBenthos(Benthos_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Benthos_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Benthos_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTBenthos
    
    !--------------------------------------------------------------------------

    subroutine GetBenthosSize(Benthos_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Benthos_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Benthos_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthosSize
    
    !--------------------------------------------------------------------------

    subroutine GetBenthosPropIndex (Benthos_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: Benthos_ID
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

        call Ready(Benthos_ID, ready_)    
        
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

    end subroutine GetBenthosPropIndex

    !--------------------------------------------------------------------------

    subroutine GetBenthosRateFlux(BenthosID, FirstProp, SecondProp, RateFlux, STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: BenthosID
        integer,           intent(IN )      :: FirstProp, SecondProp
        real,    dimension(:), pointer      :: RateFlux
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthosID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBenthos_, Me%InstanceID)

            RateFlux => Me%Matrix(:, FirstProp, SecondProp)

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBenthosRateFlux

    !--------------------------------------------------------------------------


    subroutine UnGetBenthos(BenthosID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BenthosID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthosID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBenthos_, Me%InstanceID, "UnGetBenthos3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBenthos

    !--------------------------------------------------------------------------

    subroutine UnGetBenthosRateFlux(BenthosID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BenthosID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthosID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBenthos_, Me%InstanceID, "UnGetBenthosRateFlux")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBenthosRateFlux


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyBenthos(ObjBenthosID, Temperature, Oxygen, OpenPoints, Mass, WaterVolume, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthosID
        real,    dimension(:  ), pointer            :: Temperature, Oxygen
        integer, dimension(:  ), pointer, optional  :: OpenPoints
        real,    dimension(:,:), pointer            :: Mass
        real(8), dimension(:  ), pointer, optional  :: WaterVolume
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: Index
        logical                                     :: Compute
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBenthosID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Temperature  => Temperature
            if (.not. associated(Me%ExternalVar%Temperature))       &
                stop 'ModifyBenthos - ModuleBenthos - ERR01'

            Me%ExternalVar%Oxygen       => Oxygen
            if (.not. associated(Me%ExternalVar%Oxygen))            &
                stop 'ModifyBenthos - ModuleBenthos - ERR02'
              
            Me%ExternalVar%Mass         => Mass
            if (.not. associated(Me%ExternalVar%Mass))              &
                stop 'ModifyBenthos - ModuleBenthos - ERR03'


            do Index = Me%Size%ILB, Me%Size%IUB

                if (present(OpenPoints)) then
                    if (OpenPoints(Index) == OpenPoint) then
                        Compute = .true.
                    else
                        Compute = .false.
                    endif
                else
                    Compute = .true.
                endif

                Me%Matrix(Index,:,:) = 0.

                if(Compute)then
                    
                    if(Me%ComputeOptions%Phyto     ) call ComputeBenthicPhyto       (Index)
                    
                    if(Me%ComputeOptions%Bacteria  ) call ComputeBenthicBacteria    (Index, WaterVolume)
                  
                    if(Me%ComputeOptions%Diatoms   ) call ComputeBenthicDiatoms     (Index)
                    
                    if(Me%ComputeOptions%Nitrogen  ) call ComputeBenthicNitrogen    (Index)
                    
                    if(Me%ComputeOptions%Phosphorus) call ComputeBenthicPhosphorus  (Index)
                    
                    if(Me%ComputeOptions%Silica    ) call ComputeBenthicSilica      (Index)

                endif

            enddo

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyBenthos

    !--------------------------------------------------------------------------
    
    subroutine ComputeBenthicNitrogen(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: AM, PON, O2
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------

        AM  = Me%PropIndex%Ammonia
        PON = Me%PropIndex%PON
        O2  = Me%PropIndex%Oxygen

        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max(Me%ExternalVar%Oxygen(Index)*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Nitrogen%PONDecayRate       *  &
                             Me%Nitrogen%PONDecayTFactor    ** &
                            (Me%ExternalVar%Temperature(Index) - 20.0)

       
        !kgN * day * day-1 (what passes from PON to ammonia)
        Me%Matrix(Index, PON, AM) = Me%ExternalVar%Mass(PON, Index) * Me%DTDay * &
                                    MineralizationRate * OxygenLimitation

        Me%ExternalVar%Mass(AM,  Index) = Me%ExternalVar%Mass(AM , Index) + Me%Matrix(Index, PON, AM)

        Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) - Me%Matrix(Index, PON, AM)

        !what is consumed of oxygen due to mineralization of PON
        OxygenConsumption               = Me%Matrix(Index, PON, AM) * 1. / Me%OrganicMatter%NC_Ratio * &
                                          32. / 12.

        Me%ExternalVar%Mass(O2, Index ) = Me%ExternalVar%Mass(O2, Index ) - OxygenConsumption

    
    end subroutine ComputeBenthicNitrogen
    
    !--------------------------------------------------------------------------


    subroutine ComputeBenthicPhosphorus(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: IP, POP, O2
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------
        
        IP  = Me%PropIndex%Phosphate
        POP = Me%PropIndex%POP
        O2  = Me%PropIndex%Oxygen
        
        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max(Me%ExternalVar%Oxygen(Index)*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Phosphorus%POPDecayRate       *  &
                             Me%Phosphorus%POPDecayTFactor    ** &
                            (Me%ExternalVar%Temperature(Index) - 20.0)


        !kgP * day * day-1 (what passes from POP to inorganic phosphorus)
        Me%Matrix(Index, POP, IP) = Me%ExternalVar%Mass(POP, Index) * Me%DTDay * &
                         MineralizationRate * OxygenLimitation


        Me%ExternalVar%Mass(IP,  Index) = Me%ExternalVar%Mass(IP , Index) + Me%Matrix(Index, POP, IP)

        Me%ExternalVar%Mass(POP, Index) = Me%ExternalVar%Mass(POP, Index) - Me%Matrix(Index, POP, IP)

        OxygenConsumption               = Me%Matrix(Index, POP, IP) * 1. / Me%OrganicMatter%PC_Ratio * &
                                          32. / 12.

        Me%ExternalVar%Mass(O2, Index ) = Me%ExternalVar%Mass(O2, Index ) - OxygenConsumption

                         
    end subroutine ComputeBenthicPhosphorus
    
    !--------------------------------------------------------------------------


    subroutine ComputeBenthicSilica(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: BioSi, Sil

        !Begin-----------------------------------------------------------------
        
        Sil   = Me%PropIndex%DissolvedSilica
        BioSi = Me%PropIndex%BioSilica

        !kg * day * day-1 (what passes from biogenic silica to inorganic dissolved silica)
        Me%Matrix(Index, BioSi, Sil) = Me%ExternalVar%Mass(BioSi, Index) * Me%DTDay * Me%Silica%BioSiDecayRate

        Me%ExternalVar%Mass(Sil, Index)   = Me%ExternalVar%Mass(Sil,   Index) + Me%Matrix(Index, BioSi, Sil)

        Me%ExternalVar%Mass(BioSi, Index) = Me%ExternalVar%Mass(BioSi, Index) - Me%Matrix(Index, BioSi, Sil)

    
    end subroutine ComputeBenthicSilica
    
    !--------------------------------------------------------------------------

    subroutine ComputeBenthicPhyto(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: PON, POP, Phyto
        real                                        :: Mortality
        
        !Begin-----------------------------------------------------------------
        
        PON   = Me%PropIndex%PON
        POP   = Me%PropIndex%POP
        Phyto = Me%PropIndex%Phyto

        !kg * day * day-1
        Mortality = Me%ExternalVar%Mass(Phyto, Index) * Me%DTDay * Me%Phyto%MortalityRate

        Me%ExternalVar%Mass(Phyto, Index) = Me%ExternalVar%Mass(Phyto, Index) - Mortality

        if(Me%ComputeOptions%Nitrogen)then
            
            !what passes from Phyto to PON
            Me%Matrix(Index, Phyto, PON)      = Mortality * Me%Phyto%NC_Ratio

            Me%ExternalVar%Mass(PON,   Index) = Me%ExternalVar%Mass(PON,   Index) + &
                                                Me%Matrix(Index, Phyto, PON)

        end if

        if(Me%ComputeOptions%Phosphorus)then

            !what passes from Phyto to POP
            Me%Matrix(Index, Phyto, POP)      = Mortality * Me%Phyto%PC_Ratio

            Me%ExternalVar%Mass(POP,   Index) = Me%ExternalVar%Mass(POP,   Index) + &
                                                Me%Matrix(Index, Phyto, POP)

        end if


    end subroutine ComputeBenthicPhyto
    
    !--------------------------------------------------------------------------

    subroutine ComputeBenthicBacteria(Index, WaterVolume)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        real(8), dimension(:), pointer              :: WaterVolume

        !Local-----------------------------------------------------------------
        integer                                     :: PON
        integer                                     :: BAC
        integer                                     :: DONnr
        integer                                     :: PONr
!        integer                                     :: O
        real                                        :: s1, s2, xa, xb, ya, yb
        real                                        :: Mortality                        = null_real
        real                                        :: BacteriaPONUptake                = null_real
        real                                        :: BacteriaDONUptake                = null_real      
        real                                        :: BacteriatotalUptakeM             = null_real
        real                                        :: BacteriaPONUptakeM               = null_real
        real                                        :: BacteriaDONUptakeM               = null_real
        !Temporary-------------------------------------------------------------
        real                                        :: TBacteriaLimitationFactor        = null_real                          
        real                                        :: BottomPONConc                    = null_real
        real                                        :: BottomDONConc                    = null_real
        
        !Begin-----------------------------------------------------------------    

        PON   = Me%PropIndex%PON
        PONr  = Me%PropIndex%PONr
        BAC   = Me%PropIndex%Bacteria
        DONnr = Me%PropIndex%DONnr

        !Concentration of bottom_PON in terms of the bottom reactor 
        ! WaterVolume is the VolumeNew of the Node and the reactor Volume is taken as a percentage
        !of this Volume
        ! *1000 is the conversion from kg/m3=g/l to mg/l

        BottomPONConc       = Me%ExternalVar%Mass(PON, index)                                    &
                            /(WaterVolume(index)*Me%Bact%BotReactorDepth)*1000.
        BottomDONConc       = Me%ExternalVar%Mass(DONnr, index)                                  &
                            /(WaterVolume(index))*1000.

       !TBacteriaLimitationFactor : limitation by temperature
        s1 = (1. / (Me%Bact%TOptBacteriaMin - Me%Bact%TBacteriaMin)) &
        * log((Me%Bact%BK2 * (1.0 - Me%Bact%BK1))                 &
           / (Me%Bact%BK1 * (1.0 - Me%Bact%BK2)))

        s2 = (1. / (Me%Bact%TBacteriaMax - Me%Bact%TOptBacteriaMax)) &
        * log((Me%Bact%BK3 * (1.0 - Me%Bact%BK4))                 &
           / (Me%Bact%BK4 * (1.0 - Me%Bact%BK3)))

        ya = exp(s1 * (Me%ExternalVar%Temperature(index) - Me%Bact%TBacteriaMin))
        yb = exp(s2 * (Me%Bact%TBacteriaMax - Me%ExternalVar%Temperature(index)))

        xa = (Me%Bact%BK1 * ya) / (1.0 + Me%Bact%BK1 * (ya - 1.0))
        xb = (Me%Bact%BK4 * yb) / (1.0 + Me%Bact%BK4 * (yb - 1.0))

        TBacteriaLimitationFactor = xa * xb

        if ((BottomPONConc > Me%Bact%BacMinSub) .OR.                                                 &
            (BottomDONConc > Me%Bact%BacminSub)) then         

            if (BottomPONConc > Me%Bact%BacMinSub) then         
                ! Bacteria PON uptake (1/d)
                BacteriaPONUptake   = TBacteriaLimitationFactor                                      &
                                     * Me%Bact%BacMaxUptake                                          &
                                     * BottomPONConc                                                 &
                                     / (Me%Bact%NSatConstBac                                         &
                                       + BottomPONConc) 
                ! Bacteria PON uptake (kgC/d)                                    
                BacteriaPONUptakeM  = BacteriaPONUptake                                              &
                                    * Me%ExternalVar%Mass(BAC, Index) 
                                                                   
                if (BacteriaPONUptakeM > Me%ExternalVar%Mass(PON, index)) then
                    BacteriaPONUptakeM = Me%ExternalVar%Mass(PON, index)   
                endif
                ! Bacteria PON uptake (kgN/d)                                    
                BacteriaPONUptakeM  = BacteriaPONUptakeM                                             &
                                    * Me%Bact%NC_Ratio

            else
                BacteriaPONUptake     = 0.0
                BacteriaPONUptakeM    = 0.0
            endif

            if (BottomDONConc > Me%Bact%BacMinSub) then                                                                                                          
                ! Bacteria DON uptake (1/d)
                BacteriaDONUptake   =  TBacteriaLimitationFactor                                     &
                                    * Me%Bact%BacMaxUptake                                           &
                                    * BottomDONConc                                                  &
                                    / (Me%Bact%NSatConstBac                                          &
                                      + BottomDONConc)
                ! Bacteria DON uptake (kgC/d)   
                BacteriaDONUptakeM  = BacteriaDONUptake                                              &
                                    * Me%ExternalVar%Mass(BAC, Index)                                

                if (BacteriaDONUptakeM > Me%ExternalVar%Mass(DONnr, index)) then
                    BacteriaDONUptakeM = Me%ExternalVar%Mass(DONnr, index)   
                endif
                ! Bacteria PON uptake (kgN/d)                                    
                BacteriaDONUptakeM  = BacteriaDONUptakeM                                             &
                                    * Me%Bact%NC_Ratio

            else
                 BacteriaDONUptake     = 0.0
                 BacteriaDONUptakeM    = 0.0
            endif
        else
            BacteriaPONUptake     = 0.0
            BacteriaPONUptakeM    = 0.0
            BacteriaDONUptake     = 0.0
            BacteriaDONUptakeM    = 0.0            
        endif

        !BacteriaTotalUptake, uptake in kgN/d
        BacteriaTotalUptakeM = (BacteriaDONUptakeM                                                 &
                              + BacteriaPONUptakeM)                                           
                                          

        ! Bacteria total uptake in kgC/d
        BacteriatotalUptakeM= BacteriatotalUptakeM                                                 &
                            / Me%Bact%NC_Ratio                                


        !Organic nitrogen from dead bacteria in kgC/d
        Mortality           = Me%ExternalVar%Mass(BAC, Index)                                    &
                            * Me%Bact%MortalityRate                                          
                                         


        Me%ExternalVar%Mass(BAC, Index)    = Me%ExternalVar%Mass(BAC, Index)                   &
                                           +(BacteriatotalUptakeM                              &
                                            -Mortality)                                        &
                                           * Me%DTDay 


        Me%ExternalVar%Mass(PONr, Index)   = Me%ExternalVar%Mass(PONr, Index)                  &
                                           + Mortality *  Me%Bact%NC_Ratio                     &
                                           * Me%DTDay 


        Me%ExternalVar%Mass(PON, Index)   = Me%ExternalVar%Mass(PON, Index)                    &
                                          - BacteriaPONUptakeM                                 &
                                          * Me%DTDay 

 
        Me%ExternalVar%Mass(DONnr, Index) = Me%ExternalVar%Mass(DONnr, Index)                  &
                                          - BacteriaDONUptakeM                                 &
                                          * Me%DTDay                                    
!Debug
!if(Index == 198 .AND. Me%ExternalVar%Mass(DONnr, Index) .LE. 0.0) then
!write(*,*) "BacteriaDONUptake ", BacteriaDONUptake, "BottomDONConc ", BottomDONConc
!write(*,*) "BacteriaTotalUptake ", BacteriaTotalUptake, "Mortality ", Mortality
!write(*,*) "Me%ExternalVar%Mass(BAC, Index)   ", Me%ExternalVar%Mass(BAC, Index)
!write(*,*) "Me%ExternalVar%Mass(PON, Index)   ", Me%ExternalVar%Mass(PON, Index)
!write(*,*) "Me%ExternalVar%Mass(DONnr, Index) ", Me%ExternalVar%Mass(DONnr, Index)
!write(*,*)
!endif                                        

    end subroutine ComputeBenthicBacteria
    
    !--------------------------------------------------------------------------

    subroutine ComputeBenthicDiatoms(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: PON, POP, BioSi, Diatoms
        real                                        :: Mortality
        
        !Begin-----------------------------------------------------------------
        
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        BioSi   = Me%PropIndex%BioSilica
        Diatoms = Me%PropIndex%Diatoms

        !kg * day * day-1
        Mortality = Me%ExternalVar%Mass(Diatoms, Index) * Me%DTDay * Me%Diatoms%MortalityRate

        Me%ExternalVar%Mass(Diatoms, Index) = Me%ExternalVar%Mass(Diatoms, Index) - Mortality

        if(Me%ComputeOptions%Nitrogen)then
            
            !what passes from Diatoms to PON
            Me%Matrix(Index, Diatoms, PON)      = Mortality * Me%Diatoms%NC_Ratio

            Me%ExternalVar%Mass(PON,     Index) = Me%ExternalVar%Mass(PON,     Index) + &
                                                  Me%Matrix(Index, Diatoms, PON)

        end if

        if(Me%ComputeOptions%Phosphorus)then

            !what passes from Diatoms to POP
            Me%Matrix(Index, Diatoms, POP)      = Mortality * Me%Diatoms%PC_Ratio

            Me%ExternalVar%Mass(POP,     Index) = Me%ExternalVar%Mass(POP,     Index) + &
                                                  Me%Matrix(Index, Diatoms, POP)

        end if

        if(Me%ComputeOptions%Silica)then
            
            !what passes from Diatoms to Biogenic Silica
            Me%Matrix(Index, Diatoms, BioSi)    = Mortality * Me%Diatoms%SiC_Ratio

            Me%ExternalVar%Mass(BioSi,   Index) = Me%ExternalVar%Mass(BioSi,   Index) + &
                                                  Me%Matrix(Index, Diatoms, BioSi)

        end if


    end subroutine ComputeBenthicDiatoms

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillBenthos(ObjBenthosID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjBenthosID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBenthosID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBenthos_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%PropertyList)

                !Deallocates Instance
                call DeallocateInstance ()


                ObjBenthosID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine KillBenthos
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Benthos), pointer          :: AuxObjBenthos
        type (T_Benthos), pointer          :: PreviousObjBenthos

        !Updates pointers
        if (Me%InstanceID == FirstObjBenthos%InstanceID) then
            FirstObjBenthos => FirstObjBenthos%Next
        else
            PreviousObjBenthos => FirstObjBenthos
            AuxObjBenthos      => FirstObjBenthos%Next
            do while (AuxObjBenthos%InstanceID /= Me%InstanceID)
                PreviousObjBenthos => AuxObjBenthos
                AuxObjBenthos      => AuxObjBenthos%Next
            enddo

            !Now update linked list
            PreviousObjBenthos%Next => AuxObjBenthos%Next

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

    subroutine Ready (ObjBenthos_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthos_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBenthos_ID > 0) then
            call LocateObjBenthos (ObjBenthos_ID)
            ready_ = VerifyReadLock (mBenthos_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjBenthos (ObjBenthosID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthosID

        !Local-----------------------------------------------------------------

        Me => FirstObjBenthos
        do while (associated (Me))
            if (Me%InstanceID == ObjBenthosID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleBenthos - LocateObjBenthos - ERR01'

    end subroutine LocateObjBenthos

    !--------------------------------------------------------------------------

end module ModuleBenthos

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------







