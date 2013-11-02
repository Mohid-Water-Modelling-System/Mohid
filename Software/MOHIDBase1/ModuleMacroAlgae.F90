!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : MacroAlgae
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
!
!DataFile
!   DT                  : real          [3600]          !Time step compute biogeochemical processes
!   PELAGIC_MODEL       : string    [WaterQuality]      !Pelagic biogeochemical module coupled

!   NITROGEN            : 0/1           [0]
!   PHOSPHORUS          : 0/1

!<begin_macroalgae>
!   GROWMAX             : real          [0.4]           !macroalgae maximum growth rate
!   TOPTMIN             : real          [20.]           !macroalgae optimum minimum temperature for growth
!   TOPTMAX             : real          [25.]           !macroalgae optimum maximum temperature for growth
!   TMIN                : real          [5]             !macroalgae minimum temperature for growth
!   TMAX                : real          [40.]           !macroalgae maximum temperature for growth
!   TCONST1             : real          [0.05]          !constant to control temperature response curve shape
!   TCONST2             : real          [0.98]          !constant to control temperature response curve shape
!   TCONST3             : real          [0.98]          !constant to control temperature response curve shape
!   TCONST4             : real          [0.02 ]         !constant to control temperature response curve shape
!   PHOTOIN             : real          [90.]           !macroalgae optimum radiation value
!   ENDREPC             : real          [0.009]         !macroalgae endogenous respiration rate 
!   PHOTORES            : real          [0.018]         !macroalgae photorespiration rate
!   EXCRCONS            : real          [0.008]         !macroalgae excretion rate
!   MORTMAX             : real          [0.003]         !macroalgae natural mortality rate
!   MORTCON             : real          [0.03]          !macroalgae mortality half saturation constant
!   GRAZCONS            : real          [0.00008]       !grazing rate over macroalgae
!   SOLEXCR             : real          [0.25]          !fraction of soluble inorganic material excreted by macroalgae
!   DISSDON             : real          [0.25]          !fraction of dissolved organic material excreted by macroalgae
!   NSATCONS            : real          [0.065]         !nitrogen half-saturation constant for macroalgae
!   PSATCONS            : real          [0.001]         !phosphorus half-saturation constant for macroalgae
!   RATIONC             : real          [0.18]          !macroalgae nitrogen/carbon ratio
!   RATIOPC             : real          [0.024]         !macroalgae phosphorus/carbon ratio
!   MACROALGAE_MINCONC  : real          [1e-16]         !minimum residual value for macroalgae abundance
!   MIN_OXYGEN          : real          [1e-8]          !minimum oxygen concentration for macroalgae growth
!   DEPLIM              : real          [5e-6]          !maximum SPM deposition flux for macroalgae to grow (kg m-2 s-1)
!   EROCRITSS           : real          [0.144]         !critical shear stress for macroalgae dettachment to occur (in Pa)
!   SALT_EFFECT         : real          [0]             !include salinity limitation on macroalgae growth
!   SALTOPT             : real          [20]            !macroalgae optimum salinity for growth
!   SALTCRIT            : real          [5.]            !macroalgae critical salinity limit growth
!   SALTMIN             : real          [0. ]           !macroalgae minimum salinity for growth
!   SALTMAX             : real          [45.]           !macroalgae maximum salinity for growth
!   BEACHED_MORT_RATE   : real          [0.01]          !beached drifting macroalgae mortality rate
!<end_macroalgae>

Module ModuleMacroAlgae

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleFunctions

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartMacroAlgae
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ReadMacroAlgaeParameters
    private ::          ConsistentOptions
    private ::      PropertyIndexNumber
    private ::      ConstructPropertyList

    !Selector
    public  :: GetMacroAlgaePropertyList
    public  :: GetDTMacroAlgae
    public  :: GetMacroAlgaeSize
    public  :: GetMacroAlgaeOptions
    public  :: GetMacroAlgaePropIndex
    public  :: GetMacroAlgaeRateFlux
    public  :: UnGetMacroAlgae
    public  :: UnGetMacroAlgaeRateFlux

    !Modifier
    public  :: ModifyMacroAlgae
    private ::      ComputeMacroAlgae

    !Destructor
    public  :: KillMacroAlgae                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjMacroAlgae 
    
    !Interfaces----------------------------------------------------------------

    !Types---------------------------------------------------------------------
    type     T_External
        real,       pointer, dimension(:  )         :: Temperature
        real,       pointer, dimension(:  )         :: Salinity
        real,       pointer, dimension(:  )         :: ShearStress
        real,       pointer, dimension(:  )         :: SPMDepositionFlux
        real,       pointer, dimension(:  )         :: SWRadiation
        real,       pointer, dimension(:  )         :: Thickness
        real,       pointer, dimension(:  )         :: SWLightExctintionCoef
        real,       pointer, dimension(:,:)         :: Mass
        integer,    pointer, dimension(:  )         :: OpenPoints
        real,       pointer, dimension(:  )         :: Occupation
    end type T_External

    type      T_PropIndex
        integer                                     :: Ammonia              = null_int        
        integer                                     :: Phosphate            = null_int
        integer                                     :: DissolvedSilica      = null_int
        integer                                     :: BioSilica            = null_int
        integer                                     :: PON                  = null_int         
        integer                                     :: POP                  = null_int
        integer                                     :: Oxygen               = null_int
        integer                                     :: Phyto                = null_int
        integer                                     :: Diatoms              = null_int
        integer                                     :: MacroAlgae           = null_int
        integer                                     :: DriftingMacroAlgae   = null_int
        integer                                     :: Nitrate              = null_int
        integer                                     :: DONnr                = null_int
        integer                                     :: DOPnr                = null_int
    end type T_PropIndex

    type      T_RateIndex
        integer                                     :: MA                   = null_int        
        integer                                     :: Phosphate            = null_int
        integer                                     :: DissolvedSilica      = null_int
        integer                                     :: BioSilica            = null_int
        integer                                     :: PON                  = null_int         
        integer                                     :: POP                  = null_int
        integer                                     :: Oxygen               = null_int
        integer                                     :: Phyto                = null_int
        integer                                     :: Diatoms              = null_int
        integer                                     :: MacroAlgae           = null_int
        integer                                     :: DriftingMacroAlgae   = null_int
        integer                                     :: Nitrate              = null_int
        integer                                     :: DONnr                = null_int
        integer                                     :: DOPnr                = null_int
    end type T_RateIndex

    type     T_Parameters
        real                                        :: GrowthMaxRate        = null_real 
        real                                        :: TOptMin              = null_real
        real                                        :: TOptMax              = null_real
        real                                        :: TMin                 = null_real
        real                                        :: TMax                 = null_real
        real                                        :: K1                   = null_real
        real                                        :: K2                   = null_real
        real                                        :: K3                   = null_real
        real                                        :: K4                   = null_real
        real                                        :: Photoinhibition      = null_real
        real                                        :: EndogRespCons        = null_real
        real                                        :: PhotoRespFactor      = null_real
        real                                        :: ExcretionCons        = null_real
        real                                        :: MortMaxRate          = null_real
        real                                        :: MortSatCons          = null_real
        real                                        :: GrazCons             = null_real
        real                                        :: InorgExcrFraction    = null_real
        real                                        :: DissOrgExcrFraction  = null_real
        real                                        :: NSatCons             = null_real
        real                                        :: PSatCons             = null_real
        real                                        :: RatioNC              = null_real
        real                                        :: RatioPC              = null_real
        real                                        :: ErosCritShear        = null_real
        real                                        :: SPMDepositionLimit   = null_real
        real                                        :: NitrateRatioO2N      = null_real
        real                                        :: IPRatioO2P           = null_real
        real                                        :: PhotosyntRatioO2C    = null_real
        real                                        :: MaximumSalt          = null_real
        real                                        :: MinimumSalt          = null_real
        real                                        :: CriticalSalt         = null_real
        real                                        :: OptimumSalt          = null_real
        real                                        :: MinimumConcentration = null_real
        real                                        :: MinimumOxygen        = null_real
        real                                        :: BeachedMortalityRate = null_real
        logical                                     :: SalinityEffect       = .false.
    end type T_Parameters

    
    type     T_ComputeOptions
        logical                                     :: Nitrogen             = .false.
        logical                                     :: Phosphorus           = .false.
    end type T_ComputeOptions

    type       T_MacroAlgae
        integer                                     :: InstanceID
        type (T_Size1D)                             :: Prop
        type (T_Size1D)                             :: Size
        integer                                     :: ObjEnterData         = 0
        real                                        :: DT, DTDay
        character(len=StringLength)                 :: PelagicModel
        integer, dimension(:    ), pointer          :: PropertyList
        real,    dimension(:,:,:), pointer          :: Matrix
        real,    dimension(:,:,:), pointer          :: Parameters
        type(T_PropIndex     )                      :: PropIndex
        type(T_Parameters    )                      :: Attached
        type(T_Parameters    )                      :: Drifting
        type(T_External      )                      :: ExternalVar
        type(T_ComputeOptions)                      :: ComputeOptions
        type(T_MacroAlgae    ), pointer             :: Next
    end type  T_MacroAlgae

    !Global Module Variables
    type (T_MacroAlgae), pointer                    :: FirstObjMacroAlgae
    type (T_MacroAlgae), pointer                    :: Me

    !--------------------------------------------------------------------------

    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartMacroAlgae(ObjMacroAlgaeID, FileName, ILB, IUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjMacroAlgaeID 
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
        if (.not. ModuleIsRegistered(mMacroAlgae_)) then
            nullify (FirstObjMacroAlgae)
            call RegisterModule (mMacroAlgae_) 
        endif

        call Ready(ObjMacroAlgaeID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%Size%ILB = ILB
            Me%Size%IUB = IUB

            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartMacroAlgae - ModuleMacroAlgae - ERR01'

            call ReadData

            call ConsistentOptions

            call PropertyIndexNumber
        
            call ConstructPropertyList

            call ConstructRates
 
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'StartMacroAlgae - ModuleMacroAlgae - ERR02'


            !Returns ID
            ObjMacroAlgaeID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'StartMacroAlgae - ModuleMacroAlgae - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartMacroAlgae
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_MacroAlgae), pointer                         :: NewObjMacroAlgae
        type (T_MacroAlgae), pointer                         :: PreviousObjMacroAlgae


        !Allocates new instance
        allocate (NewObjMacroAlgae)
        nullify  (NewObjMacroAlgae%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjMacroAlgae)) then
            FirstObjMacroAlgae         => NewObjMacroAlgae
            Me                      => NewObjMacroAlgae
        else
            PreviousObjMacroAlgae      => FirstObjMacroAlgae
            Me                      => FirstObjMacroAlgae%Next
            do while (associated(Me))
                PreviousObjMacroAlgae  => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjMacroAlgae
            PreviousObjMacroAlgae%Next => NewObjMacroAlgae
        endif

        Me%InstanceID = RegisterNewInstance (mMacroAlgae_)


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
                                    block_begin     = '<begin_macroalgae>',         &
                                    block_end       = '<end_macroalgae>',           &
                                    BlockFound      = BlockFound,                   &
                                    STAT            = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then
            
            if (BlockFound) then

                call ReadMacroAlgaeParameters(Me%Attached)

            else
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                write(*,*) 'Could not find <begin_macroalgae>...<end_macroalgae> block'
                stop 'ReadData - ModuleMacroAlgae - ERR01'

            end if

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            stop       'ReadData - ModuleMacroAlgae - ERR02'
        
        else
            stop       'ReadData - ModuleMacroAlgae - ERR03'
        end if

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                    ClientNumber    = ClientNumber,                 &
                                    block_begin     = '<begin_driftingmacroalgae>', &
                                    block_end       = '<end_driftingmacroalgae>',   &
                                    BlockFound      = BlockFound,                   &
                                    STAT            = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then
            
            if (BlockFound) then

                call ReadMacroAlgaeParameters(Me%Drifting)

            else
                
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)
                write(*,*) 'Could not find <begin_driftingmacroalgae>...<end_driftingalgae> block'
                stop 'ReadData - ModuleMacroAlgae - ERR04'

            end if

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            stop       'ReadData - ModuleMacroAlgae - ERR05'
        else
            stop       'ReadData - ModuleMacroAlgae - ERR04'
        end if

    end subroutine ReadData

    !--------------------------------------------------------------------------

    subroutine ConsistentOptions
        
        !Begin-----------------------------------------------------------------

        Me%Drifting%MinimumConcentration = null_real
        Me%Drifting%SPMDepositionLimit   = null_real
        Me%Drifting%ErosCritShear        = null_real
        Me%Attached%BeachedMortalityRate = null_real

    
    end subroutine ConsistentOptions

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
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleMacroAlgae - ERR01'

        Me%DTDay = Me%DT / (3600. * 24.)

        call GetData(Me%PelagicModel,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PELAGIC_MODEL',                                    &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleMacroAlgae - ERR02'
        if(iflag==0)then
            write(*,*)'Please define the pelagic model to couple with ModuleMacroAlgae'
            stop 'ConstructGlobalVariables - ModuleMacroAlgae - ERR20'
        end if

        if((Me%PelagicModel .ne. WaterQualityModel))then
            write(*,*)'Pelagic model to couple with ModuleMacroAlgae must be:'
            write(*,*)trim(WaterQualityModel)
            stop 'ConstructGlobalVariables - ModuleMacroAlgae - ERR30'
        endif
                 
        call GetData(Me%ComputeOptions%Nitrogen,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NITROGEN',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleMacroAlgae - ERR03'


        call GetData(Me%ComputeOptions%Phosphorus,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHOSPHORUS',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleMacroAlgae - ERR04'

    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------

    subroutine ReadMacroAlgaeParameters(Parameters)

        !Arguments-------------------------------------------------------------
        type(T_Parameters)                              :: Parameters

        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Parameters%GrowthMaxRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'GROWMAX',                                          &
                     Default      = 0.4,                                                &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR01'

        call GetData(Parameters%MaximumSalt,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SALTMAX',                                          &
                     Default      = 45.0,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR02'

        call GetData(Parameters%MinimumSalt,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SALTMIN',                                          &
                     Default      = 0.0,                                                &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR03'
        
        
        call GetData(Parameters%CriticalSalt,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SALTCRIT',                                         &
                     Default      = 5.0,                                                &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR04'

        call GetData(Parameters%OptimumSalt,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SALTOPT',                                          &
                     Default      = 20.0,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR05'

        call GetData(Parameters%TOptMin,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TOPTMIN',                                          &
                     Default      = 20.0,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR06'

        call GetData(Parameters%TOptMax,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TOPTMAX',                                          &
                     Default      = 25.0,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR07'

        call GetData(Parameters%TMin,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TMIN',                                             &
                     Default      = 5.0,                                                &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR08'

        call GetData(Parameters%TMax,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TMAX',                                             &
                     Default      = 40.0,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR09'


        call GetData(Parameters%K1,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TCONST1',                                          &
                     Default      = 0.05,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR10'

        call GetData(Parameters%K2,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TCONST2',                                          &
                     Default      = 0.98,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR11'
        
        call GetData(Parameters%K3,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TCONST3',                                          &
                     Default      = 0.98,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR12'
        
        call GetData(Parameters%K4,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TCONST4',                                          &
                     Default      = 0.02,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR13'


        call GetData(Parameters%Photoinhibition,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PHOTOIN',                                          &
                     Default      = 90.0,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR14'


        call GetData(Parameters%EndogRespCons,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'ENDREPC',                                          &
                     Default      = 0.009,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR15'

        call GetData(Parameters%PhotoRespFactor,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PHOTORES',                                         &
                     Default      = 0.018,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR16'


        call GetData(Parameters%ExcretionCons,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'EXCRCONS',                                         &
                     Default      = 0.008,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR17'

        call GetData(Parameters%MortMaxRate,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MORTMAX',                                          &
                     Default      = 0.003,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR18'
       
       
        call GetData(Parameters%MortSatCons,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MORTCON',                                          &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR19'

        call GetData(Parameters%GrazCons,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'GRAZCONS',                                         &
                     Default      = 0.00008,                                            &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR20'
                     
        call GetData(Parameters%InorgExcrFraction,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SOLEXCR',                                          &
                     Default      = 0.25,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR21'
               
        call GetData(Parameters%DissOrgExcrFraction,                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DISSDON',                                          &
                     Default      = 0.25,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR22'

        call GetData(Parameters%NSatCons,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NSATCONS',                                         &
                     Default      = 0.065,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR23'

        call GetData(Parameters%PSatCons,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PSATCONS',                                         &
                     Default      = 0.001,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR24'
                        
        call GetData(Parameters%RatioNC,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'RATIONC',                                          &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR25'
            
        call GetData(Parameters%RatioPC,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'RATIOPC',                                          &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR26'
         
        !Critical shear stress for macroalgae dettachment to occur (in Pa)
        call GetData(Parameters%ErosCritShear,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'EROCRITSS',                                        &
                     Default      = 0.144,                                              &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR27'

        !Maximum SPM deposition flux with which macroalgae is able to grow (kg m-2 s-1)
        call GetData(Parameters%SPMDepositionLimit,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DEPLIM',                                           &
                     Default      = 5e-6,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR28'


        call GetData(Parameters%NitrateRatioO2N,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NITONRAT',                                         &
                     Default      = 48.0 / 14.0,                                        &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR29'

        call GetData(Parameters%IPRatioO2P,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PHOSOPRAT',                                        &
                     Default      = 64.0 / 31.0,                                        &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR30'
                   
        call GetData(Parameters%PhotosyntRatioO2C,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PHOTOSOC',                                         &
                     Default      = 32.0/12.0,                                          &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR31'
        
        call GetData(Parameters%SalinityEffect,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'SALT_EFFECT',                                      &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR32'

        call GetData(Parameters%MinimumConcentration,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MACROALGAE_MINCONC',                               &
                     Default      = AlmostZero,                                         &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR33'

        !Minimum oxygen concentration allowed in mg/l
        call GetData(Parameters%MinimumOxygen,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'MIN_OXYGEN',                                       &
                     Default      = 1e-5,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR34'

        !Mortality rate of beached drifting macroalgae 
        call GetData(Parameters%BeachedMortalityRate,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BEACHED_MORT_RATE',                                &
                     Default      = 0.01,                                               &
                     ClientModule = 'ModuleMacroAlgae',                                 &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadMacroAlgaeParameters - ModuleMacroAlgae - ERR35'



    end subroutine ReadMacroAlgaeParameters

    !--------------------------------------------------------------------------

    subroutine PropertyIndexNumber
        
        !Begin-----------------------------------------------------------------
        
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        !MacroAlgae
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%MacroAlgae         = Me%Prop%IUB

        !DriftingMacroAlgae
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%DriftingMacroAlgae = Me%Prop%IUB


        !Oxygen
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB

        if(Me%ComputeOptions%Nitrogen)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB
        
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Nitrate        = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DONnr          = Me%Prop%IUB
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%PON            = Me%Prop%IUB

        end if

        if(Me%ComputeOptions%Phosphorus)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DOPnr          = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POP            = Me%Prop%IUB

        end if


    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyList

        !Begin-----------------------------------------------------------------

        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
        Me%PropertyList(Me%PropIndex%MacroAlgae)            = MacroAlgae_
        Me%PropertyList(Me%PropIndex%DriftingMacroAlgae)    = DriftingMacroAlgae_
        Me%PropertyList(Me%PropIndex%Oxygen)                = Oxygen_

        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%Ammonia)           = Ammonia_
            Me%PropertyList(Me%PropIndex%Nitrate)           = Nitrate_
            Me%PropertyList(Me%PropIndex%DONnr)             = DONNon_Refractory_
            Me%PropertyList(Me%PropIndex%PON)               = PON_
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%Phosphate)         = Inorganic_Phosphorus_
            Me%PropertyList(Me%PropIndex%DOPnr)             = DOPNon_Refractory_
            Me%PropertyList(Me%PropIndex%POP)               = POP_
        end if

    end subroutine ConstructPropertyList

    !--------------------------------------------------------------------------
    
    subroutine ConstructRates

        !Local-----------------------------------------------------------------
        integer                                             :: nParameters
        integer                                             :: nMacroAlgaeTypes

        !Begin-----------------------------------------------------------------
        allocate(Me%Matrix  (Me%Size%ILB:Me%Size%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))

        !Attached and drifting macroalgae
        nMacroAlgaeTypes = 2

        !Gross production              for each macro algae type    : +1
        !Light       limiting factor for each macro algae type    : +1
        !Nutrients   limiting factor for each macro algae type    : +1
        !Temperature limiting factor for each macro algae type    : +1
        !Salinity    limiting factor for each macro algae type    : +1
        !Nitrogen    limiting factor for each macro algae type    : +1
        !Phosphorus  limiting factor for each macro algae type    : +1
        nParameters = 7

        !Allocates with 2nd dimension size = 2
        !if J = 1 Attached Macroalgae
        !if J = 2 Drifting Macroalgae
        allocate(Me%Parameters(Me%Size%ILB:Me%Size%IUB, 1:nMacroAlgaeTypes, 1:nParameters))
        Me%Parameters = 0.

    end subroutine ConstructRates


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine GetMacroAlgaePropertyList(Life_ID, PropertyList, STAT)

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

            call Read_Lock(mMacroAlgae_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetMacroAlgaePropertyList

    !--------------------------------------------------------------------------
    
    subroutine GetDTMacroAlgae(MacroAlgae_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: MacroAlgae_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(MacroAlgae_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTMacroAlgae
    
    !--------------------------------------------------------------------------

    subroutine GetMacroAlgaeOptions(MacroAlgae_ID, Salinity, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: MacroAlgae_ID
        logical, optional, intent(OUT)      :: Salinity
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(MacroAlgae_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Salinity)) then
                if(Me%Attached%SalinityEffect .or. Me%Drifting%SalinityEffect)then
                    Salinity = .true.
                else
                    Salinity = .false.
                end if
            end if

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetMacroAlgaeOptions

    !--------------------------------------------------------------------------

    subroutine GetMacroAlgaeSize(MacroAlgae_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: MacroAlgae_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(MacroAlgae_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetMacroAlgaeSize
    
    !--------------------------------------------------------------------------

    subroutine GetMacroAlgaePropIndex (MacroAlgae_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: MacroAlgae_ID
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

        call Ready(MacroAlgae_ID, ready_)    
        
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

    end subroutine GetMacroAlgaePropIndex

    !--------------------------------------------------------------------------
    
    subroutine GetMacroAlgaeRateFlux(MacroAlgaeID, FirstProp, SecondProp, RateFlux, STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: MacroAlgaeID
        integer,           intent(IN )      :: FirstProp, SecondProp
        real,    dimension(:), pointer      :: RateFlux
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(MacroAlgaeID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mMacroAlgae_, Me%InstanceID)

            select case(FirstProp)
                
                case(GrossProd_)

                    RateFlux => Me%Parameters(:, SecondProp, MA_GrossFact_  )

                case(TemperatureLim_)
                    
                    RateFlux => Me%Parameters(:, SecondProp, MA_TLimFact_   )

                case(NutrientLim_)

                    RateFlux => Me%Parameters(:, SecondProp, MA_NutLimFact_ )

                case(LightLim_)

                    RateFlux => Me%Parameters(:, SecondProp, MA_LLimFact_   )

                case(SalinityLim_)

                    RateFlux => Me%Parameters(:, SecondProp, MA_SLimFact_   )

                case(NLim_)

                    RateFlux => Me%Parameters(:, SecondProp, MA_NLimFact_   )

                case(PLim_)

                    RateFlux => Me%Parameters(:, SecondProp, MA_PLimFact_   )

                case default

                    RateFlux => Me%Matrix    (:, FirstProp, SecondProp      )

            end select

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetMacroAlgaeRateFlux

    !--------------------------------------------------------------------------

    subroutine UnGetMacroAlgae(MacroAlgaeID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: MacroAlgaeID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(MacroAlgaeID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mMacroAlgae_, Me%InstanceID, "UnGetMacroAlgae")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetMacroAlgae
    
    !--------------------------------------------------------------------------

    subroutine UnGetMacroAlgaeRateFlux(MacroAlgaeID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: MacroAlgaeID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(MacroAlgaeID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mMacroAlgae_, Me%InstanceID, "UnGetMacroAlgaeRateFlux")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetMacroAlgaeRateFlux

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyMacroAlgae(ObjMacroAlgaeID, Temperature, Salinity, OpenPoints, &
                                ShearStress, SPMDepositionFlux, SWRadiation,        &
                                SWLightExctintionCoef, Thickness, Occupation, Mass, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjMacroAlgaeID
        real,    dimension(:  ), pointer            :: Temperature
        real,    dimension(:  ), pointer, optional  :: Salinity
        real,    dimension(:  ), pointer, optional  :: ShearStress
        real,    dimension(:  ), pointer, optional  :: SPMDepositionFlux
        real,    dimension(:  ), pointer, optional  :: SWRadiation
        real,    dimension(:  ), pointer, optional  :: SWLightExctintionCoef
        real,    dimension(:  ), pointer, optional  :: Thickness
        real,    dimension(:  ), pointer, optional  :: Occupation
        integer, dimension(:  ), pointer, optional  :: OpenPoints
        real,    dimension(:,:), pointer            :: Mass
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: Index

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjMacroAlgaeID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Temperature  => Temperature
            if (.not. associated(Me%ExternalVar%Temperature))               &
                stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR01'
              
            Me%ExternalVar%Mass         => Mass
            if (.not. associated(Me%ExternalVar%Mass))                      &
                stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR02'

            if(present(Salinity))then
                Me%ExternalVar%Salinity             => Salinity
                if (.not. associated(Me%ExternalVar%Salinity))              &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR03'
            end if

            if(present(ShearStress))then
                Me%ExternalVar%ShearStress          => ShearStress
                if (.not. associated(Me%ExternalVar%ShearStress))           &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR04'
            end if

            if(present(SPMDepositionFlux))then
                Me%ExternalVar%SPMDepositionFlux    => SPMDepositionFlux
                if (.not. associated(Me%ExternalVar%SPMDepositionFlux))     &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR05'
            end if

            if(present(SWRadiation))then
                Me%ExternalVar%SWRadiation          => SWRadiation
                if (.not. associated(Me%ExternalVar%SWRadiation))           &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR06'
            end if

            if(present(SWLightExctintionCoef))then
                Me%ExternalVar%SWLightExctintionCoef=> SWLightExctintionCoef
                if (.not. associated(Me%ExternalVar%SWLightExctintionCoef)) &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR07'
            end if

            if(present(Thickness))then
                Me%ExternalVar%Thickness  => Thickness
                if (.not. associated(Me%ExternalVar%Thickness)) &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR08'
            end if
            
           if(present(Occupation))then
                Me%ExternalVar%Occupation  => Occupation
                if (.not. associated(Me%ExternalVar%Occupation)) &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR08'
            end if
            
            if(present(OpenPoints))then
                Me%ExternalVar%OpenPoints  => OpenPoints
                if (.not. associated(Me%ExternalVar%OpenPoints)) &
                    stop 'ModifyMacroAlgae - ModuleMacroAlgae - ERR09'
            end if

            do Index = Me%Size%ILB, Me%Size%IUB

                !Reset rates matrix
                Me%Matrix(Index,:,:) = 0.

                call ComputeMacroAlgae(Me%Attached, OFF, Index)

                call ComputeMacroAlgae(Me%Drifting, ON,  Index)

            enddo

            nullify(Me%ExternalVar%Temperature          )
            nullify(Me%ExternalVar%Mass                 )
            nullify(Me%ExternalVar%Salinity             )
            nullify(Me%ExternalVar%ShearStress          )
            nullify(Me%ExternalVar%SPMDepositionFlux    )
            nullify(Me%ExternalVar%SWRadiation          )
            nullify(Me%ExternalVar%SWLightExctintionCoef)
            nullify(Me%ExternalVar%Thickness            )
            nullify(Me%ExternalVar%Occupation           )

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyMacroAlgae

    !--------------------------------------------------------------------------
    
    subroutine ComputeMacroAlgae(Parameters, Drifting, Index)

        !Arguments-------------------------------------------------------------
        type(T_Parameters)                          :: Parameters
        logical, intent(IN)                         :: Drifting
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: PON, DONnr
        integer                                     :: POP, DOPnr
        integer                                     :: NH4, NO3, IP
        integer                                     :: O2, MA, DMA
        real                                        :: TempLimitingFactor 
        real                                        :: LightLimitingFactor
        real                                        :: NutLimitingFactor  
        real                                        :: NitrogenLimitFactor  
        real                                        :: PhosphorusLimitFactor
        real                                        :: SaltLimitingFactor 
        real                                        :: MacAlgGrossGrowRate        
        real                                        :: s1,s2, ya, yb, xa, xb
        real                                        :: sx, exp_salt
        real                                        :: MacAlgNonGrazingMortalityRate
        real                                        :: MacAlgExcretionRate
        real                                        :: MacAlgEndResp, MacAlgPhotoResp
        real                                        :: MacAlgRespirationRate
        real                                        :: MacAlgNitrogenExcretionRate
        real                                        :: x1, x2, x3, x4
        real                                        :: MacAlgAmmoniaPreferenceFactor
        real                                        :: MacAlgPhosphorusExcretionRate
        real                                        :: NitrateOxygen, IPOxygen
        real                                        :: OxygenRespRate, PhotoOxygenProduction
        integer                                     :: Zone
        integer, parameter                          :: NoLimitation = 1
        integer, parameter                          :: Erosion      = 2
        integer, parameter                          :: NoCompute    = 3
        integer, parameter                          :: Beached      = 4
        real                                        :: MacroAlgaeMassOld
        real                                        :: DeadMass
        real, parameter                             :: minthickness = 0.001
        real                                        :: DZ1, DZ2, radiation_at_top_canopy
        
        !Begin-----------------------------------------------------------------
        
        PON     = Me%PropIndex%PON
        DONnr   = Me%PropIndex%DONnr
        POP     = Me%PropIndex%POP
        DOPnr   = Me%PropIndex%DOPnr
        O2      = Me%PropIndex%Oxygen
        NH4     = Me%PropIndex%Ammonia
        NO3     = Me%PropIndex%Nitrate
        IP      = Me%PropIndex%Phosphate

        if(Drifting)then
            MA      = Me%PropIndex%DriftingMacroAlgae
        else
            MA      = Me%PropIndex%MacroAlgae
            DMA     = Me%PropIndex%DriftingMacroAlgae
        endif

        if(Me%ExternalVar%OpenPoints(Index) == OpenPoint)then

            if(Drifting)then
                Zone = NoLimitation
            else
                !Selects type of limitation due to SPM deposition or shear stress
                if    (Me%ExternalVar%SPMDepositionFlux(Index) > Parameters%SPMDepositionLimit)then
                    Zone = NoCompute
                elseif(Me%ExternalVar%ShearStress(Index)       > Parameters%ErosCritShear     )then
                    Zone = Erosion
                else
                    Zone = NoLimitation
                end if

            end if

        else

            if(Drifting)then
                Zone = Beached
            else
                Zone = NoCompute
            end if

        end if

        select case(Zone)

            case(NoLimitation)

                !Temperature limiting factor
                s1 = (1.0 / (Parameters%TOptMin - Parameters%TMin)) * &
                      log((Parameters%K2 * (1.0 - Parameters%K1))   / &
                          (Parameters%K1 * (1.0 - Parameters%K2)))

                s2 = (1.0 / (Parameters%TMax - Parameters%TOptMax))  * &
                      log((Parameters%K3 * (1.0 - Parameters%K4))   / &
                          (Parameters%K4 * (1.0 - Parameters%K3)))

                ya = exp(s1 * (Me%ExternalVar%Temperature(Index) - Parameters%TMin))
                yb = exp(s2 * (Parameters%TMax - Me%ExternalVar%Temperature(Index)))

                xa = (Parameters%K1 * ya) / (1.0 + Parameters%K1 * (ya - 1.0))
                xb = (Parameters%K4 * yb) / (1.0 + Parameters%K4 * (yb - 1.0))

                TempLimitingFactor = xa * xb


                !Nitrogen limiting factor
                if(Me%ComputeOptions%Nitrogen)then
                    
                    NitrogenLimitFactor   = (Me%ExternalVar%Mass(NH4, Index)  +  &
                                             Me%ExternalVar%Mass(NO3, Index)) /  &
                                            (Parameters%NSatCons              +  &
                                             Me%ExternalVar%Mass(NH4, Index)  +  &
                                             Me%ExternalVar%Mass(NO3, Index))
                else
                    NitrogenLimitFactor    = 1.0   
                end if
                
                !Phosphorus limiting factor
                if(Me%ComputeOptions%Phosphorus)then
                    
                    PhosphorusLimitFactor   = Me%ExternalVar%Mass(IP, Index) / &
                                             (Me%ExternalVar%Mass(IP, Index) + &
                                              Parameters%PSatCons)
                else
                    PhosphorusLimitFactor   = 1.0   
                end if


       
                !Nutrients limiting factor
                NutLimitingFactor = min(NitrogenLimitFactor, PhosphorusLimitFactor) 
                
                ! DZ1 is the distance (m) between the top of the cell and the top of the canopy
                ! (minthickness is used to avoid division by 0 if DZ1 is 0)
                DZ1= max(minthickness, (1. - Me%ExternalVar%Occupation(index))*Me%ExternalVar%Thickness(index)) 
                ! (minthickness is used to avoid division by 0 if DZ2 is 0)
                ! DZ2 is macroalgae height in the cell
                DZ2= max(minthickness,Me%ExternalVar%Occupation(index)*Me%ExternalVar%Thickness(index))
        
                if (DZ1 == minthickness) then
                ! the height of canopy reaches the top of the cell, so the radiation at top of cell is used
                    radiation_at_top_canopy = Me%ExternalVar%SWRadiation(index)  
                else
                    radiation_at_top_canopy = Me%ExternalVar%SWRadiation(index) * &
                                              exp(-DZ1*Me%ExternalVar%SWLightExctintionCoef(index))
                 end if
         
               !It is assumed that the light extinction coefficient is uniform in the cell (it is rough approximation)
         
                LightLimitingFactor =                                                                                      &
                          PhytoLightLimitationFactor(Thickness       = DZ2,                                         &
                                                     TopRadiation    = radiation_at_top_canopy,                     &
                                                     PExt            = Me%ExternalVar%SWLightExctintionCoef(index), &
                                                     Photoinhibition = Parameters%Photoinhibition)


                !Salt LimitingFactor
                if (Parameters%SalinityEffect) then
           
                    if (Me%ExternalVar%Salinity(Index) .GE. Parameters%CriticalSalt) then
           
                        if (Me%ExternalVar%Salinity(Index)  .LT. Parameters%OptimumSalt) then
                           Sx = Parameters%MinimumSalt
                           exp_salt = 2.5
                        else
                           Sx = Parameters%MaximumSalt
                           exp_salt = 2.0
                        end if

                        SaltLimitingFactor = 1. - ((Me%ExternalVar%Salinity(Index) - &
                                                    Parameters%OptimumSalt)        / &
                                                   (Sx - Parameters%OptimumSalt))** exp_salt
                    else
                        SaltLimitingFactor =      (Me%ExternalVar%Salinity(Index)  - &
                                                   Parameters%MinimumSalt)         / &
                                                  (Parameters%OptimumSalt          - &
                                                   Parameters%MinimumSalt) 
                    end if

                    if(SaltLimitingFactor < 0.) SaltLimitingFactor = 0.0

                else 

                    SaltLimitingFactor = 1.0

                end if
 
                !Compute gross growth rate 
                MacAlgGrossGrowRate = Parameters%GrowthMaxRate  * &
                                      TempLimitingFactor        * &
                                      NutLimitingFactor         * &
                                      LightLimitingFactor       * &
                                      SaltLimitingFactor


                if    (MacAlgGrossGrowRate  < 0.0)then

                    write(*,*)'Negative macroalgae gross growth rate.'
                    write(*,*)"TempLimitingFactor", TempLimitingFactor 
                    write(*,*)"NutLimitingFactor",  NutLimitingFactor  
                    write(*,*)"LightLimitingFactor",LightLimitingFactor
                    write(*,*)"SaltLimitingFactor", SaltLimitingFactor

                    MacAlgGrossGrowRate = -1.0 / null_real  !Avoid division by zero below

                elseif(MacAlgGrossGrowRate == 0.0)then

                    MacAlgGrossGrowRate = -1.0 / null_real  !Avoid division by zero below

                end if

                !Minimum oxygen conditions must occur to macroalgae to grow
                if(Me%ExternalVar%Mass(O2, Index) .lt. Parameters%MinimumOxygen)then
                    MacAlgGrossGrowRate = -1.0 / null_real  !Avoid division by zero below
                end if

                !Macroalgae natural mortality rate
                MacAlgNonGrazingMortalityRate =  Parameters%MortMaxRate                 * &
                                               ((Me%ExternalVar%Mass(MA, Index)         / &
                                                 MacAlgGrossGrowRate)                   / &
                                                (Parameters%MortSatCons                 + &
                                                 Me%ExternalVar%Mass(MA, Index)         / & 
                                                 MacAlgGrossGrowRate))

                !Macroalgae excretion rate 
                MacAlgExcretionRate = Parameters%ExcretionCons * MacAlgGrossGrowRate    * &
                                     (1.0 - LightLimitingFactor)


                !Macroalgae endogenous respiration
                MacAlgEndResp = Parameters%EndogRespCons * &
                                exp(0.069 * Me%ExternalVar%Temperature(Index))
                        

                if(Me%ExternalVar%Mass(O2, Index) .lt. Parameters%MinimumOxygen)then
                    MacAlgEndResp = -1.0 / null_real !Avoid division by zero below
                endif

                !Macroalgae photo respiration rate
                MacAlgPhotoResp = Parameters%PhotoRespFactor * MacAlgGrossGrowRate                      
        
                !Macroalgae respiration rate
                MacAlgRespirationRate = MacAlgEndResp + MacAlgPhotoResp


                !Nitrogen excretions and NH4/NO3 preference factor
                if(Me%ComputeOptions%Nitrogen)then

                    !1/day * mgN/mgC
                    MacAlgNitrogenExcretionRate = (MacAlgExcretionRate  + MacAlgRespirationRate + &
                                                   Parameters%GrazCons) * Parameters%RatioNC

                    !MacAlgAmmoniaPreferenceFactor
                    x1 = Me%ExternalVar%Mass(NH4, Index) * Me%ExternalVar%Mass(NO3, Index)

                    x2 = (Parameters%NSatCons + Me%ExternalVar%Mass(NH4, Index)) * &
                         (Parameters%NSatCons + Me%ExternalVar%Mass(NO3, Index))

                    x3 = Parameters%NSatCons * Me%ExternalVar%Mass(NH4, Index)

                    x4 = (Me%ExternalVar%Mass(NH4, Index) + Me%ExternalVar%Mass(NO3, Index)) * &
                         (Parameters%NSatCons             + Me%ExternalVar%Mass(NO3, Index))

                    if ((x1 == 0.) .AND. (x3 == 0.)) then
                        MacAlgAmmoniaPreferenceFactor = 0.                 
                    else 
                        MacAlgAmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
                    end if

                end if

                !Phosphorus excretions
                if(Me%ComputeOptions%Phosphorus)then

                    MacAlgPhosphorusExcretionRate = (MacAlgExcretionRate  + MacAlgRespirationRate + &
                                                     Parameters%GrazCons) * Parameters%RatioPC

                end if

                !Oxygen production rate due to nitrate uptake
                if(Me%ComputeOptions%Nitrogen)then
                    
                    NitrateOxygen = MacAlgGrossGrowRate * Parameters%NitrateRatioO2N * &
                                    Parameters%RatioNC * (1 - MacAlgAmmoniaPreferenceFactor)
                else
                    NitrateOxygen = 0.       
                end if        

                !Oxygen production rate due to inorganic phosphorus uptake
                if (Me%ComputeOptions%Phosphorus) then
                    IPOxygen      = MacAlgGrossGrowRate * Parameters%IPRatioO2P * &
                                    Parameters%RatioPC 
                else
                    IPOxygen = 0.0
                end if

                !Oxygen photosynthetic production
                PhotoOxygenProduction = MacAlgGrossGrowRate * Parameters%PhotosyntRatioO2C

                !Oxygen loss by respiration 
                OxygenRespRate = MacAlgRespirationRate * Parameters%PhotosyntRatioO2C 

                !Store macroalgae old concentration
                MacroAlgaeMassOld = Me%ExternalVar%Mass(MA, Index)

                !Compute gross production
                Me%Parameters(Index, MA, MA_GrossFact_  ) = MacAlgGrossGrowRate     * &
                                                            MacroAlgaeMassOld       * Me%DT
                !Compute limiting factors                                            
                Me%Parameters(Index, MA, MA_TLimFact_   ) = TempLimitingFactor      * Me%DT
                Me%Parameters(Index, MA, MA_NLimFact_   ) = NitrogenLimitFactor     * Me%DT
                Me%Parameters(Index, MA, MA_PLimFact_   ) = PhosphorusLimitFactor   * Me%DT
                Me%Parameters(Index, MA, MA_NutLimFact_ ) = NutLimitingFactor       * Me%DT
                Me%Parameters(Index, MA, MA_LLimFact_   ) = LightLimitingFactor     * Me%DT
                Me%Parameters(Index, MA, MA_SLimFact_   ) = SaltLimitingFactor      * Me%DT


                !New macroalgae mass (kg = kg + day-1 * kg * day)
                Me%ExternalVar%Mass(MA, Index)      = Me%ExternalVar%Mass(MA, Index)        + &
                                                     (MacAlgGrossGrowRate                   - &
                                                      MacAlgRespirationRate                 - &
                                                      MacAlgExcretionRate                   - &
                                                      MacAlgNonGrazingMortalityRate         - &
                                                      Parameters%GrazCons)                  * &
                                                      MacroAlgaeMassOld                     * &
                                                      Me%DTDay

                if (Me%ComputeOptions%Nitrogen)then

                    !what passes from ammonia to macroalgae
                    Me%Matrix(Index, NH4, MA) = MacAlgGrossGrowRate                         * &
                                                MacAlgAmmoniaPreferenceFactor               * &
                                                Parameters%RatioNC                          * &
                                                MacroAlgaeMassOld                           * &
                                                Me%DTDay

                    !what passes from macroalgae to ammonia 
                    Me%Matrix(Index, MA, NH4) = Parameters%InorgExcrFraction                * &
                                                MacAlgNitrogenExcretionRate                 * &
                                                MacroAlgaeMassOld                           * &
                                                Me%DTDay

                    !ammonia mass balance
                    Me%ExternalVar%Mass(NH4, Index) = Me%ExternalVar%Mass(NH4, Index)       + &
                                                      Me%Matrix(Index, MA, NH4)             - &
                                                      Me%Matrix(Index, NH4, MA)


                    !what passes from nitrate to macroalgae
                    Me%Matrix(Index, NO3, MA) = MacAlgGrossGrowRate                         * &
                                               (1.-MacAlgAmmoniaPreferenceFactor)           * &
                                                Parameters%RatioNC                          * &
                                                MacroAlgaeMassOld                           * &
                                                Me%DTDay

                    !nitrate mass balance
                    Me%ExternalVar%Mass(NO3, Index) = Me%ExternalVar%Mass(NO3, Index)       - &
                                                      Me%Matrix(Index, NO3, MA)


                   
                    !what passes from macroalgae to DONnr 
                    Me%Matrix(Index, MA, DONnr)     = (Parameters%DissOrgExcrFraction       * &
                                                      (1.- Parameters%InorgExcrFraction)    * &
                                                      MacAlgNitrogenExcretionRate)          * &
                                                      MacroAlgaeMassOld * Me%DTDay
                   
                    !DONnr mass balance
                    Me%ExternalVar%Mass(DONnr,Index)= Me%ExternalVar%Mass(DONnr, Index)     + &
                                                      Me%Matrix(Index, MA, DONnr)

                    !what passes from macroalgae to PON 
                    Me%Matrix(Index, MA, PON)       = ((1.- Parameters%DissOrgExcrFraction) * &
                                                       (1.- Parameters%InorgExcrFraction)   * &
                                                        MacAlgNitrogenExcretionRate         + &
                                                        MacAlgNonGrazingMortalityRate       * &
                                                        Parameters%RatioNC)                 * &
                                                        MacroAlgaeMassOld * Me%DTDay

                    Me%ExternalVar%Mass(PON,Index)  = Me%ExternalVar%Mass(PON, Index)       + &
                                                      Me%Matrix(Index, MA, PON)
                endif


                if (Me%ComputeOptions%Phosphorus) then
                    
                    !what passes from inorganic phosphorus to macroalgae 
                    Me%Matrix(Index, IP, MA)        = Parameters%RatioPC                    * &
                                                      MacAlgGrossGrowRate                   * &
                                                      MacroAlgaeMassOld * Me%DTDay

                    !what passes from macroalgae to inorganic phosphorus 
                    Me%Matrix(Index, MA, IP)        = Parameters%InorgExcrFraction          * &
                                                      MacAlgPhosphorusExcretionRate         * &
                                                      MacroAlgaeMassOld * Me%DTDay
                   
                    !inorganic phosphorus mass balance
                    Me%ExternalVar%Mass(IP, Index)  = Me%ExternalVar%Mass(IP, Index)        + &
                                                      Me%Matrix(Index, MA, IP)              - &
                                                      Me%Matrix(Index, IP, MA)

                    !what passes from macroalgae to DOPnr
                    Me%Matrix(Index, MA, DOPnr)     = Parameters%DissOrgExcrFraction        * &
                                                      (1.- Parameters%InorgExcrFraction)    * &
                                                      MacAlgPhosphorusExcretionRate         * &
                                                      MacroAlgaeMassOld * Me%DTDay
                    !DOPnr mass balance
                    Me%ExternalVar%Mass(DOPnr,Index)= Me%ExternalVar%Mass(DOPnr, Index)     + &
                                                      Me%Matrix(Index, MA, DOPnr)

                    !what passes from macroalgae to POP
                    Me%Matrix(Index, MA, POP)       = ((1.- Parameters%DissOrgExcrFraction) * &
                                                     (1.- Parameters%InorgExcrFraction)     * &
                                                      MacAlgPhosphorusExcretionRate         + &
                                                      MacAlgNonGrazingMortalityRate         * &
                                                      Parameters%RatioPC)                   * &
                                                      MacroAlgaeMassOld * Me%DTDay

                    !POP mass balance
                    Me%ExternalVar%Mass(POP,Index)  = Me%ExternalVar%Mass(POP, Index)       + &
                                                      Me%Matrix(Index, MA, POP)
                end if

                !what passes from oxygen to macroalgae
                Me%Matrix(Index, O2, MA)            = OxygenRespRate * MacroAlgaeMassOld * Me%DTDay


                !what passes from macroalgae to oxygen
                Me%Matrix(Index, MA, O2)            = (PhotoOxygenProduction + IPOxygen + NitrateOxygen) * &
                                                       MacroAlgaeMassOld * Me%DTDay
                !oxygen mass balance
                Me%ExternalVar%Mass(O2,Index)       = Me%ExternalVar%Mass(O2, Index)        + &
                                                      Me%Matrix(Index, MA, O2)              - &
                                                      Me%Matrix(Index, O2, MA)

            case(NoCompute)

                !NOTHING HAPPENS

            case(Erosion)

                !Not all macroalgae mass is dettached to make sure 
                !there's always enough to grow back again (minimum concentration)
                if(Me%ExternalVar%Mass(MA, Index) > Parameters%MinimumConcentration)then
                    
                    !what passes from attached macroalgae to drifting macroalgae
                    Me%Matrix(Index, MA, DMA)   = Me%ExternalVar%Mass(MA, Index) - &
                                                  Parameters%MinimumConcentration
                else
                    !what passes from attached macroalgae to drifting macroalgae
                    Me%Matrix(Index, MA, DMA)   = 0.

                end if
                                
                Me%ExternalVar%Mass(MA, Index)  = Parameters%MinimumConcentration

                Me%ExternalVar%Mass(DMA, Index) = Me%ExternalVar%Mass(DMA, Index) + &
                                                  Me%Matrix(Index, MA, DMA)
            case(Beached)

                !implicit algorithm
                
                !Store macroalgae old concentration
                MacroAlgaeMassOld = Me%ExternalVar%Mass(MA, Index)                 

                Me%ExternalVar%Mass(MA, Index)  = MacroAlgaeMassOld  / &
                                                  (1 + Parameters%BeachedMortalityRate * Me%DTDay)

                DeadMass = Me%ExternalVar%Mass(MA, Index) - MacroAlgaeMassOld

                if (Me%ComputeOptions%Nitrogen) then
                    
                    !what passes from drifting macroalgae to PON
                    Me%Matrix(Index, MA, PON)       = DeadMass * Parameters%RatioNC

                    Me%ExternalVar%Mass(PON,Index)  = Me%ExternalVar%Mass(PON, Index)       + &
                                                      Me%Matrix(Index, MA, PON)

                end if

                if (Me%ComputeOptions%Phosphorus) then
                    
                    !what passes from drifting macroalgae to POP
                    Me%Matrix(Index, MA, POP)       = DeadMass * Parameters%RatioPC

                    Me%ExternalVar%Mass(POP,Index)  = Me%ExternalVar%Mass(POP, Index)       + &
                                                      Me%Matrix(Index, MA, POP)

                end if

        end select


    end subroutine ComputeMacroAlgae

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillMacroAlgae(ObjMacroAlgaeID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjMacroAlgaeID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjMacroAlgaeID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mMacroAlgae_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%PropertyList)
                deallocate(Me%Matrix      )
                deallocate(Me%Parameters  )

                !Deallocates Instance
                call DeallocateInstance ()


                ObjMacroAlgaeID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine KillMacroAlgae
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_MacroAlgae), pointer          :: AuxObjMacroAlgae
        type (T_MacroAlgae), pointer          :: PreviousObjMacroAlgae

        !Updates pointers
        if (Me%InstanceID == FirstObjMacroAlgae%InstanceID) then
            FirstObjMacroAlgae => FirstObjMacroAlgae%Next
        else
            PreviousObjMacroAlgae => FirstObjMacroAlgae
            AuxObjMacroAlgae      => FirstObjMacroAlgae%Next
            do while (AuxObjMacroAlgae%InstanceID /= Me%InstanceID)
                PreviousObjMacroAlgae => AuxObjMacroAlgae
                AuxObjMacroAlgae      => AuxObjMacroAlgae%Next
            enddo

            !Now update linked list
            PreviousObjMacroAlgae%Next => AuxObjMacroAlgae%Next

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

    subroutine Ready (ObjMacroAlgae_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjMacroAlgae_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjMacroAlgae_ID > 0) then
            call LocateObjMacroAlgae (ObjMacroAlgae_ID)
            ready_ = VerifyReadLock (mMacroAlgae_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjMacroAlgae (ObjMacroAlgaeID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjMacroAlgaeID

        !Local-----------------------------------------------------------------

        Me => FirstObjMacroAlgae
        do while (associated (Me))
            if (Me%InstanceID == ObjMacroAlgaeID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleMacroAlgae - LocateObjMacroAlgae - ERR01'

    end subroutine LocateObjMacroAlgae

    !--------------------------------------------------------------------------

end module ModuleMacroAlgae

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------







