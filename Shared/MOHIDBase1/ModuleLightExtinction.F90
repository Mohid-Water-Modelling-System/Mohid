!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : LightExtinction
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute light extinction cofficients in the water column
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


Module ModuleLightExtinction

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleTime
    use ModuleTimeSerie
    use ModuleFunctions

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructLightExtinction
    private ::      AllocateInstance
    private ::      ReadOptions
    private ::          ConstructShortWaveParameters
    private ::          ConstructLongWaveParameters
    private ::          ConstructPlanktonParameters
    private ::          VerifyOptions
    private ::      AllocateVariables


    !Selector
    public  :: GetShortWaveExtinctionField
    public  :: GetLongWaveExtinctionCoef
    public  :: GetRadiationPercentages
    public  :: GetLightExtinctionOptions
    public  :: UnGetLightExtinction
                     
    
    !Modifier
    public  :: ModifyLightExtinctionField
    private ::      ComputeShortWaveExtinction
    private ::          CheckResetField
    private ::          ComputeUnitsCoef
    private ::          FromASCIIFile
    private ::          Compute_Parsons_Ocean
    private ::          Compute_Portela_Tagus
    private ::          Compute_Combined_ParsonsPortela
    private ::          Compute_Multiparameters
    private ::      ComputeLongWaveExtinction

    !Destructor
    public  :: KillLightExtinction                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjLightExtinction 
    
    !Interfaces----------------------------------------------------------------
    private :: ConstructLightExtinction1D
    private :: ConstructLightExtinction3D
    interface  ConstructLightExtinction
        module procedure ConstructLightExtinction1D
        module procedure ConstructLightExtinction3D
    end interface  ConstructLightExtinction

    private :: ModifyLightExtinctionField1D
    private :: ModifyLightExtinctionField3D
    interface ModifyLightExtinctionField
        module procedure ModifyLightExtinctionField1D
        module procedure ModifyLightExtinctionField3D
    end interface ModifyLightExtinctionField


    private :: GetShortWaveExtinctionField1D
    private :: GetShortWaveExtinctionField3D
    interface  GetShortWaveExtinctionField
        module procedure GetShortWaveExtinctionField1D
        module procedure GetShortWaveExtinctionField3D
    end interface  GetShortWaveExtinctionField

    private :: UnGetLightExtinction3D
    interface  UnGetLightExtinction
        module procedure UnGetLightExtinction1D
        module procedure UnGetLightExtinction3D
    end interface  UnGetLightExtinction
    
    !Parameter-----------------------------------------------------------------

    !Light extinction compute methods
    integer, parameter                              :: Constant                = 1
    integer, parameter                              :: Parsons_Ocean           = 2
    integer, parameter                              :: Portela_Tagus           = 3
    integer, parameter                              :: Combined_ParsonsPortela = 4
    integer, parameter                              :: ASCIIFile               = 5
    integer, parameter                              :: Multiparameters         = 6

    !Types---------------------------------------------------------------------    
    type       T_External
        type(T_Time)                                :: Now
        integer, pointer, dimension(:    )          :: RiverPoints1D
        integer, pointer, dimension(:,:,:)          :: WaterPoints3D
        real   , pointer, dimension(:    )          :: Concentration1D
        real   , pointer, dimension(:    )          :: PhytoConcentration1D
        real   , pointer, dimension(:    )          :: SPMConcentration1D
        real   , pointer, dimension(:,:,:)          :: Concentration3D
        real   , pointer, dimension(:,:,:)          :: PhytoConcentration3D
        real   , pointer, dimension(:,:,:)          :: SPMConcentration3D
        real   , pointer, dimension(:,:,:)          :: MacroAlgaeOccupation
        integer                                     :: PropertyID              
        real                                        :: UnitsCoef
        real                                        :: ExtinctionParameter
        real                                        :: Phyto_UnitsCoef, SPM_UnitsCoef
        logical                                     :: Phyto_OK = OFF, SPM_OK = OFF
    end type T_External
    
    type       T_Wave
        real, pointer, dimension(:)                 :: ExtinctionCoefField1D
        real, pointer, dimension(:,:,:)             :: ExtinctionCoefField3D
        integer                                     :: ExtinctionType       = FillValueInt
        real                                        :: ExtinctionCoef       = FillValueReal
        real                                        :: ExtCoefWater         = FillValueReal
        real                                        :: Percentage           = FillValueReal               
    end type T_Wave

    type       T_PlanktonParameters
        real                                        :: RatioC_Chla          = null_real
    end type   T_PlanktonParameters


    type       T_Needs
        logical                                     :: Phyto                = .false.
        logical                                     :: SPM                  = .false.
        logical                                     :: Parameters           = .false.
        logical                                     :: Concentrations       = .false.
    end type   T_Needs

    type       T_LightExtinction
        integer                                     :: InstanceID
        type (T_Size1D)                             :: Size1D
        type (T_Size3D)                             :: Size, WorkSize
        character(LEN = StringLength)               :: FileName
        type (T_Wave)                               :: ShortWave
        type (T_Wave)                               :: LongWave
        integer                                     :: SWColumn
        integer                                     :: LWColumn
        real                                        :: UnitsCoef
        logical                                     :: Is3D
        logical                                     :: MacroAlgaeExtinction = .false.
        type(T_External)                            :: ExternalVar
        type(T_Time)                                :: LastCompute
        type(T_Needs)                               :: Needs
        type(T_PlanktonParameters)                  :: PlanktonParameters

        !Instance of ModuleTime
        integer                                     :: ObjTime              = 0

        !Instance of ModuleEnterData
        integer                                     :: ObjEnterData         = 0
        
        !Instance of ModuleTimeSerie
        integer                                     :: ObjTimeSerie         = 0

        type(T_LightExtinction), pointer            :: Next
    end type  T_LightExtinction

    !Global Module Variables
    type (T_LightExtinction), pointer               :: FirstObjLightExtinction
    type (T_LightExtinction), pointer               :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructLightExtinction1D(LightExtinctionID, TimeID, EnterDataID,        &
                                          Size1D, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: LightExtinctionID
        integer                                         :: TimeID
        integer                                         :: EnterDataID
        type(T_Size1D)                                  :: Size1D
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_, nUsers       

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLightExtinction_)) then
            nullify (FirstObjLightExtinction)
            call RegisterModule (mLightExtinction_) 
        endif

        call Ready(LightExtinctionID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Is3D = .false.

            Me%ObjTime      = AssociateInstance (mTIME_,      TimeID)
            Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)
            Me%Size1D       = Size1D

            call ReadOptions

            call AllocateVariables

            call SetDate(Me%LastCompute, 0, 0, 0, 0, 0, 0)

            nUsers = DeassociateInstance(mTIME_     , Me%ObjTime     )
            if (nUsers == 0) stop 'ModuleLightExtinction - ConstructLightExtinction - ERR00' 

            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ModuleLightExtinction - ConstructLightExtinction - ERR01' 
            
            !Returns ID
            LightExtinctionID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleLightExtinction - ConstructLightExtinction3D - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructLightExtinction1D
 
    !--------------------------------------------------------------------------

    subroutine ConstructLightExtinction3D(LightExtinctionID, TimeID, EnterDataID,        &
                                          Size, WorkSize, MacroAlgaeON, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: LightExtinctionID
        integer                                         :: TimeID
        integer                                         :: EnterDataID
        type(T_Size3D)                                  :: Size, WorkSize
        logical, optional                               :: MacroAlgaeON
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_, nUsers       

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_


        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mLightExtinction_)) then
            nullify (FirstObjLightExtinction)
            call RegisterModule (mLightExtinction_) 
        endif


        call Ready(LightExtinctionID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Is3D = .true.

            Me%ObjTime      = AssociateInstance (mTIME_,      TimeID)
            Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)
            Me%Size         = Size
            Me%WorkSize     = WorkSize
            
            if(present(MacroAlgaeON))then
                Me%MacroAlgaeExtinction = MacroAlgaeON
            else
                Me%MacroAlgaeExtinction = OFF
            endif

            call ReadOptions

            call AllocateVariables

            call SetDate(Me%LastCompute, 0, 0, 0, 0, 0, 0)

            nUsers = DeassociateInstance(mTIME_     , Me%ObjTime     )
            if (nUsers == 0) stop 'ModuleLightExtinction - ConstructLightExtinction - ERR00' 

            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ModuleLightExtinction - ConstructLightExtinction - ERR01' 
            
            !Returns ID
            LightExtinctionID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleLightExtinction - ConstructLightExtinction3D - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructLightExtinction3D
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_LightExtinction), pointer                         :: NewObjLightExtinction
        type (T_LightExtinction), pointer                         :: PreviousObjLightExtinction


        !Allocates new instance
        allocate (NewObjLightExtinction)
        nullify  (NewObjLightExtinction%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjLightExtinction)) then
            FirstObjLightExtinction         => NewObjLightExtinction
            Me                              => NewObjLightExtinction
        else
            PreviousObjLightExtinction      => FirstObjLightExtinction
            Me                              => FirstObjLightExtinction%Next
            do while (associated(Me))
                PreviousObjLightExtinction  => Me
                Me                          => Me%Next
            enddo
            Me                              => NewObjLightExtinction
            PreviousObjLightExtinction%Next => NewObjLightExtinction
        endif

        Me%InstanceID = RegisterNewInstance (mLIGHTEXTINCTION_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    subroutine AllocateVariables
        
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------


        if(Me%Is3D)then
            allocate (Me%ShortWave%ExtinctionCoefField3D(Me%Size%ILB:Me%Size%IUB,       &
                      Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
            Me%ShortWave%ExtinctionCoefField3D(:,:,:) = Me%ShortWave%ExtinctionCoef
        else
            allocate (Me%ShortWave%ExtinctionCoefField1D(Me%Size1D%ILB:Me%Size1D%IUB))
            Me%ShortWave%ExtinctionCoefField1D(:)     = Me%ShortWave%ExtinctionCoef
        endif



    end subroutine AllocateVariables

    !----------------------------------------------------------------------

    subroutine ReadOptions

        !Begin-----------------------------------------------------------------

        call ConstructShortWaveParameters

        call ConstructLongWaveParameters

        call ConstructPlanktonParameters

        call VerifyOptions

    end subroutine ReadOptions


    !----------------------------------------------------------------------
    
    
    subroutine ConstructShortWaveParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%ShortWave%ExtinctionType,                                   &
                     Me%ObjEnterData,  iflag,                                       &
                     SearchType     = FromFile,                                     &
                     keyword        = 'SW_EXTINCTION_TYPE',                         &
                     default        = Constant,                                     &
                     ClientModule   = 'ModuleLightExtinction',                      &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructShortWaveParameters - ModuleLightExtinction - ERR01'

        if (Me%ShortWave%ExtinctionType == Multiparameters) then

            call GetData(Me%ShortWave%ExtCoefWater,                                 &
                         Me%ObjEnterData,  iflag,                                   &
                         SearchType     = FromFile,                                 &
                         keyword        = 'SW_KW',                                  &
                         Default        = 0.08,                                     &
                         ClientModule   = 'ModuleLightExtinction',                  &
                         STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructShortWaveParameters - ModuleLightExtinction - ERR03'

        end if

        call GetData(Me%ShortWave%ExtinctionCoef,                                   &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType     = FromFile,                                     &
                     keyword        = 'SW_EXTINCTION_COEF',                         &
                     Default        = 1/20.,                                        &
                     ClientModule   = 'ModuleLightExtinction',                      &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructShortWaveParameters - ModuleLightExtinction - ERR04'


        call GetData(Me%ShortWave%Percentage,                                       &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType     = FromFile,                                     &
                     keyword        = 'SW_PERCENTAGE',                              &
                     Default        = 0.6,                                          &
                     ClientModule   = 'ModuleLightExtinction',                      &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructShortWaveParameters - ModuleLightExtinction - ERR05'

        if (Me%ShortWave%ExtinctionType == ASCIIFile) then

            call GetData(Me%SWColumn,                                               &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType     = FromFile,                                 &
                         keyword        = 'SW_EXTINCTION_COLUMN',                   &
                         ClientModule   = 'ModuleLightExtinction',                  &
                         STAT           = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructShortWaveParameters - ModuleLightExtinction - ERR06'
            if (iflag /= 1)              stop 'ConstructShortWaveParameters - ModuleLightExtinction - ERR07'

        endif
                                                                             
    end subroutine ConstructShortWaveParameters


    !----------------------------------------------------------------------


    subroutine ConstructLongWaveParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%LongWave%ExtinctionType,                                    &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType     = FromFile,                                     &
                     keyword        = 'LW_EXTINCTION_TYPE',                         &
                     Default        = Constant,                                     &
                     ClientModule   = 'ModuleLightExtinction',                      &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLongWaveParameters - ModuleLightExtinction - ERR01'

        call GetData(Me%LongWave%ExtinctionCoef,                                    &
                     Me%ObjEnterData,                                               &
                     iflag,                                                         &
                     SearchType     = FromFile,                                     &
                     keyword        = 'LW_EXTINCTION_COEF',                         &
                     Default        = 1/3.,                                         &
                     ClientModule   = 'ModuleLightExtinction',                      &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLongWaveParameters - ModuleLightExtinction - ERR02'

        call GetData(Me%LongWave%Percentage,                                        &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType = FromFile,                                         &
                     keyword    = 'LW_PERCENTAGE',                                  &
                     Default    = 0.4,                                              &
                     ClientModule = 'ModuleLightExtinction',                        &
                     STAT       = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLongWaveParameters - ModuleLightExtinction - ERR03'

        if (Me%LongWave%ExtinctionType == ASCIIFile) then

            call GetData(Me%LWColumn,                                               &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType     = FromFile,                                 &
                         keyword        = 'LW_EXTINCTION_COLUMN',                   &
                         ClientModule   = 'ModuleLightExtinction',                  &
                         STAT           = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructLongWaveParameters - ModuleLightExtinction - ERR04'
            if (iflag /= 1)              stop 'ConstructLongWaveParameters - ModuleLightExtinction - ERR05'

        endif

                                                                             
    end subroutine ConstructLongWaveParameters

    !----------------------------------------------------------------------

    subroutine ConstructPlanktonParameters

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%PlanktonParameters%RatioC_Chla,                             &
                     Me%ObjEnterData,  iflag,                                       &
                     SearchType     = FromFile,                                     &
                     keyword        = 'RATIO_C_CHLA',                               &
                     Default        = 60.,                                          &
                     ClientModule   = 'ModuleLightExtinction',                      &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            stop 'ConstructPlanktonParameters - ModuleLightExtinction - ERR01'

    end subroutine ConstructPlanktonParameters

    !----------------------------------------------------------------------

    subroutine VerifyOptions

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        character(StringLength)                 :: SWLWExtFile
        real                                    :: SWPercentage, LWPercentage
        real                                    :: TotalPercentage

        !Begin-----------------------------------------------------------------


        SWPercentage    = 100. * Me%ShortWave%Percentage  
        LWPercentage    = 100. * Me%LongWave%Percentage 

        TotalPercentage = SWPercentage + LWPercentage
                          
        if (abs(TotalPercentage - 100) > 1.e-3) stop 'VerifyOptions - ModuleLightExtinction - ERR01'

        if (Me%ShortWave%ExtinctionType == ASCIIFile .or.                       &
            Me%LongWave%ExtinctionType  == ASCIIFile) then

            call GetData(SWLWExtFile,                                           &
                         Me%ObjEnterData, iflag,                                &
                         SearchType     = FromFile,                             &
                         keyword        = 'SW_LW_EXTINCTION_FILE',              &
                         ClientModule   = 'ModuleLightExtinction',              &
                         STAT           = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - VerifyOptions - ERR02'
            if (iflag /= 1)              stop 'ReadOptions - VerifyOptions - ERR03'
                
            !Starts TimeSerieInput
            call StartTimeSerieInput(Me%ObjTimeSerie,                           &
                                     SWLWExtFile,                               &
                                     Me%ObjTime,                                &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleLightExtinction - ERR04'
                

        endif

        if (Me%ShortWave%ExtinctionType == Parsons_Ocean          ) then
            Me%Needs%Phyto          = .true.
            Me%Needs%Concentrations = .true.
        end if

        if (Me%ShortWave%ExtinctionType == Combined_ParsonsPortela) then
            Me%Needs%Phyto          = .true.
            Me%Needs%SPM            = .true.
            Me%Needs%Concentrations = .true.
        end if


        if (Me%ShortWave%ExtinctionType == Portela_Tagus          ) then
            Me%Needs%SPM            = .true.
            Me%Needs%Concentrations = .true.
        end if

        if (Me%ShortWave%ExtinctionType == Multiparameters        ) then
            Me%Needs%Parameters     = .true.
            Me%Needs%Concentrations = .true.
        end if

        if(Me%MacroAlgaeExtinction .and. Me%ShortWave%ExtinctionType .ne. Multiparameters)then
            write(*,*)'When running with macroalgae, the short wave radiation '
            write(*,*)'extinction type must be MULTIPARAMETERS. Please choose '
            write(*,*)'option SW_EXTINCTION_TYPE : 6'
            stop 'VerifyOptions - ModuleLightExtinction - ERR05'

        end if

    end subroutine VerifyOptions

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    subroutine GetShortWaveExtinctionField1D (LightExtinctionID, ShortWaveLightExtinctionField, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LightExtinctionID
        real, dimension(:),  pointer                :: ShortWaveLightExtinctionField
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mLIGHTEXTINCTION_, Me%InstanceID)

            ShortWaveLightExtinctionField => Me%ShortWave%ExtinctionCoefField1D

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetShortWaveExtinctionField1D
    
    !--------------------------------------------------------------------------

    subroutine GetShortWaveExtinctionField3D (LightExtinctionID, ShortWaveLightExtinctionField, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: LightExtinctionID
        real, dimension(:, :, :),  pointer          :: ShortWaveLightExtinctionField
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mLIGHTEXTINCTION_, Me%InstanceID)

            ShortWaveLightExtinctionField => Me%ShortWave%ExtinctionCoefField3D

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetShortWaveExtinctionField3D

    !--------------------------------------------------------------------------
   
    subroutine GetLongWaveExtinctionCoef (LightExtinctionID, LongWaveExtinctionCoef, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: LightExtinctionID
        real                                            :: LongWaveExtinctionCoef
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            LongWaveExtinctionCoef = Me%LongWave%ExtinctionCoef

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetLongWaveExtinctionCoef
    
    !--------------------------------------------------------------------------
   
    subroutine GetRadiationPercentages (LightExtinctionID, SWPercentage, LWPercentage, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: LightExtinctionID
        real,    intent(OUT), optional                  :: SWPercentage, LWPercentage
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(LWPercentage)) LWPercentage = Me%LongWave%Percentage
            if (present(SWPercentage)) SWPercentage = Me%ShortWave%Percentage

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetRadiationPercentages

    !--------------------------------------------------------------------------
    
    subroutine GetLightExtinctionOptions (LightExtinctionID, NeedsPhyto, NeedsSPM, &
                                          NeedsParameters, NeedsConcentrations, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: LightExtinctionID
        logical             , optional                  :: NeedsPhyto
        logical             , optional                  :: NeedsSPM
        logical             , optional                  :: NeedsParameters
        logical             , optional                  :: NeedsConcentrations
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if(present(NeedsPhyto         )) NeedsPhyto          = Me%Needs%Phyto
            if(present(NeedsSPM           )) NeedsSPM            = Me%Needs%SPM
            if(present(NeedsParameters    )) NeedsParameters     = Me%Needs%Parameters
            if(present(NeedsConcentrations)) NeedsConcentrations = Me%Needs%Concentrations

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetLightExtinctionOptions

    !--------------------------------------------------------------------------

    subroutine UnGetLightExtinction1D(ObjLightExtinctionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjLightExtinctionID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLightExtinctionID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mLIGHTEXTINCTION_, Me%InstanceID,  "UnGetLightExtinction1D")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetLightExtinction1D

    !--------------------------------------------------------------------------

    subroutine UnGetLightExtinction3D(ObjLightExtinctionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjLightExtinctionID
        real, dimension(:, :, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLightExtinctionID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mLIGHTEXTINCTION_, Me%InstanceID,  "UnGetLightExtinction3D")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetLightExtinction3D


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyLightExtinctionField1D(LightExtinctionID, RiverPoints1D,         &
                                            CurrentTime, PropertyID, Concentration,   &
                                            UnitsCoef, ExtinctionParameter, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: LightExtinctionID
        type(T_Time)                                    :: CurrentTime
        integer,    dimension(:    ), pointer           :: RiverPoints1D
        integer,                               optional :: PropertyID
        real,       dimension(:)    , pointer, optional :: Concentration
        real,                                  optional :: UnitsCoef
        real,                                  optional :: ExtinctionParameter
        integer, intent(OUT),                  optional :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%RiverPoints1D    => RiverPoints1D
            Me%ExternalVar%Now              =  CurrentTime

            if(present(PropertyID         )) Me%ExternalVar%PropertyID          =  PropertyID
            if(present(Concentration      )) Me%ExternalVar%Concentration1D     => Concentration
            if(present(UnitsCoef          )) Me%ExternalVar%UnitsCoef           =  UnitsCoef
            if(present(ExtinctionParameter)) Me%ExternalVar%ExtinctionParameter =  ExtinctionParameter

            if (Me%ShortWave%ExtinctionType /= Constant) &
                call ComputeShortWaveExtinction

            call ComputeLongWaveExtinction

            if(present(Concentration      )) nullify(Me%ExternalVar%Concentration1D)

            Me%ExternalVar%PropertyID          =  null_int
            Me%ExternalVar%UnitsCoef           =  null_real
            Me%ExternalVar%ExtinctionParameter =  null_real

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLightExtinctionField1D

    !----------------------------------------------------------------------

    subroutine ModifyLightExtinctionField3D(LightExtinctionID, WaterPoints3D,         &
                                            CurrentTime, PropertyID, Concentration,   &
                                            UnitsCoef, ExtinctionParameter,           &
                                            MacroAlgaeOccupation, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: LightExtinctionID
        type(T_Time)                                    :: CurrentTime
        integer,    dimension(:,:,:), pointer           :: WaterPoints3D
        integer,                               optional :: PropertyID
        real,       dimension(:,:,:), pointer, optional :: Concentration, MacroAlgaeOccupation
        real,                                  optional :: UnitsCoef
        real,                                  optional :: ExtinctionParameter
        integer, intent(OUT),                  optional :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(LightExtinctionID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%WaterPoints3D    => WaterPoints3D
            Me%ExternalVar%Now              =  CurrentTime

            if(present(PropertyID           )) Me%ExternalVar%PropertyID          =  PropertyID
            if(present(Concentration        )) Me%ExternalVar%Concentration3D     => Concentration
            if(present(UnitsCoef            )) Me%ExternalVar%UnitsCoef           =  UnitsCoef
            if(present(ExtinctionParameter  )) Me%ExternalVar%ExtinctionParameter =  ExtinctionParameter
            if(present(MacroAlgaeOccupation )) then
                if (associated(MacroAlgaeOccupation)) Me%ExternalVar%MacroAlgaeOccupation=> MacroAlgaeOccupation
            endif
            

            if (Me%ShortWave%ExtinctionType /= Constant) &
                call ComputeShortWaveExtinction

            call ComputeLongWaveExtinction

            if(present(Concentration        )) nullify(Me%ExternalVar%Concentration3D     )
            if(present(MacroAlgaeOccupation )) nullify(Me%ExternalVar%MacroAlgaeOccupation)

            Me%ExternalVar%PropertyID          =  null_int
            Me%ExternalVar%UnitsCoef           =  null_real
            Me%ExternalVar%ExtinctionParameter =  null_real

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyLightExtinctionField3D

    !----------------------------------------------------------------------

    subroutine ComputeShortWaveExtinction


        call CheckResetField

        call ComputeUnitsCoef

        select case(Me%ShortWave%ExtinctionType)

            case(ASCIIFile)

                call FromASCIIFile

            case(Parsons_Ocean)

                call Compute_Parsons_Ocean

            case(Portela_Tagus)

                call Compute_Portela_Tagus

            case(Combined_ParsonsPortela)

                call Compute_Combined_ParsonsPortela

            case(Multiparameters)

                call Compute_Multiparameters

            case default 

                write (*,*) 
                write (*,*) 'Unknown type of extinction coefficient.'
                stop        'ComputeShortWaveExtinction - ModuleLightExtinction - ERR02'

        end select


    end subroutine ComputeShortWaveExtinction

    !--------------------------------------------------------------------------

    subroutine CheckResetField

        !Begin-----------------------------------------------------------------

        if(Me%ExternalVar%Now .gt. Me%LastCompute) then

            if(Me%ShortWave%ExtinctionType == Multiparameters) then

                if (Me%Is3D) then
                    call SetMatrixValue(Me%ShortWave%ExtinctionCoefField3D, Me%WorkSize, &
                                        Me%ShortWave%ExtCoefWater)
                else
                    Me%ShortWave%ExtinctionCoefField1D = Me%ShortWave%ExtCoefWater
                endif

            else 
                
                if (Me%Is3D) then
                    call SetMatrixValue(Me%ShortWave%ExtinctionCoefField3D, Me%WorkSize, 0.)
                else
                    Me%ShortWave%ExtinctionCoefField1D = 0.
                endif

            end if

            Me%LastCompute = Me%ExternalVar%Now

        end if


    end subroutine CheckResetField

    !--------------------------------------------------------------------------

    subroutine ComputeUnitsCoef

        !Begin-----------------------------------------------------------------

        select case(Me%ExternalVar%PropertyID)

            case(Phytoplankton_, Diatoms_)

                !Phytoplankton and Diatoms units must be in mg Chl/m3 so conversion factor 
                !must be multiplied by 1e6 and by the conversion between Chla and carbon (1/60. defaul)
                Me%UnitsCoef    = Me%ExternalVar%UnitsCoef * 1e6 * 1/Me%PlanktonParameters%RatioC_Chla

                !Chlorophyll units must be in mg/m3 so conversion factor 
                !must be multiplied by 1e6.
            case(Diatom_Chl_, Mix_Flagellate_Chl_, Picoalgae_Chl_, Flagellate_Chl_)
                
                Me%UnitsCoef    = Me%ExternalVar%UnitsCoef * 1e6

            case default

                !Other units must be in mg/l so conversion factor 
                !must be multiplied by 1000.
                Me%UnitsCoef    = Me%ExternalVar%UnitsCoef * 1000.

        end select


    end subroutine ComputeUnitsCoef

    !--------------------------------------------------------------------------

    subroutine FromASCIIFile

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        logical                                 :: TimeCycle
        
        !Local-----------------------------------------------------------------
        type(T_Time)                            :: Time1,Time2
        real                                    :: Value1, Value2, NewValue
        integer                                 :: i

        !Begin-----------------------------------------------------------------

        !Gets Value arround current instant    
        call GetTimeSerieValue(Me%ObjTimeSerie,                             &
                               Me%ExternalVar%Now,                          &
                               Me%SWColumn,                                 &
                               Time1, Value1, Time2, Value2, TimeCycle,     &
                               STAT = STAT_CALL)
        if (STAT_CALL/= SUCCESS_)                                           &
            stop 'FromASCIIFile - ModuleLightExtinction - ERR01'

        if (TimeCycle) then
            NewValue = Value1
        else
            !Interpolates Value for current instant
            call InterpolateValueInTime(Me%ExternalVar%Now, Time1,          &
                                        Value1, Time2, Value2,              &
                                        Me%ShortWave%ExtinctionCoef)
        endif

        if (Me%Is3D) then

            call SetMatrixValue(Me%ShortWave%ExtinctionCoefField3D,         &
                                Me%WorkSize, Me%ShortWave%ExtinctionCoef,   &
                                Me%ExternalVar%WaterPoints3D)
        else

            do i = Me%Size1D%ILB, Me%Size1D%IUB
                if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then
                    Me%ShortWave%ExtinctionCoefField1D(i)       = Me%ShortWave%ExtinctionCoef
                endif
            enddo
            
        endif


    end subroutine FromASCIIFile

    !--------------------------------------------------------------------------

    subroutine Compute_Parsons_Ocean

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k
        real,   dimension(:)    , pointer       :: Concentration1D
        real,   dimension(:,:,:), pointer       :: Concentration3D

        !Begin-----------------------------------------------------------------


        if(Me%ExternalVar%PropertyID == Phytoplankton_ .or. Me%ExternalVar%PropertyID == Diatoms_)then

            if (Me%Is3D) then

                Concentration3D => Me%ExternalVar%Concentration3D

                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then 
                
                        !factor de extincao para alto mar, funcao da concentracao de fitoplancton
                        !Eq. 59 pp 95 Parsons et al(Biol.Ocea.Proc.) C[mg Chl a/m3],tem que haver 
                        !conversão de unidades pq o modelo calcula o fito em mgC/L.

                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) =                                    &
                                        0.04 + 0.0088 * (Concentration3D(i, j, k) * Me%UnitsCoef) +    &
                                        0.54 * ((Concentration3D(i, j, k) * Me%UnitsCoef) ** (2.0 / 3.0))

                    end if

                enddo
                enddo
                enddo

            else

                Concentration1D => Me%ExternalVar%Concentration1D

                do i = Me%Size1D%ILB, Me%Size1D%IUB

                    if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then 
                
                        !factor de extincao para alto mar, funcao da concentracao de fitoplancton
                        !Eq. 59 pp 95 Parsons et al(Biol.Ocea.Proc.) C[mg Chl a/m3],tem que haver 
                        !conversão de unidades pq o modelo calcula o fito em mgC/L.

                        Me%ShortWave%ExtinctionCoefField1D(i) =                                        &
                                        0.04 + 0.0088 * (Concentration1D(i) * Me%UnitsCoef) +          &
                                        0.54 * ((Concentration1D(i) * Me%UnitsCoef) ** (2.0 / 3.0))

                    end if

                enddo
               
            endif

        else
            
            write (*,*) 
            write (*,*) 'Property '//trim(GetPropertyName(Me%ExternalVar%PropertyID))
            write (*,*) 'cannot be used to compute Parsons_Ocean light extinction type.'
            stop        'Compute_Parsons_Ocean - ModuleLightExtinction - ERR01'

        end if



    end subroutine Compute_Parsons_Ocean


    !--------------------------------------------------------------------------

    subroutine Compute_Combined_ParsonsPortela

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k
        real,   dimension(:)    , pointer       :: PhytoConcentration1D, SPMConcentration1D
        real,   dimension(:,:,:), pointer       :: PhytoConcentration3D, SPMConcentration3D
        real                                    :: Phyto_UnitsCoef, SPM_UnitsCoef

        !Begin------------------------------------------------------------------

       
        if(Me%ExternalVar%PropertyID == Phytoplankton_ .or. Me%ExternalVar%PropertyID == Diatoms_)then

            Me%ExternalVar%Phyto_UnitsCoef    =  Me%UnitsCoef
            Me%ExternalVar%Phyto_OK           =  .true.
            if (Me%Is3D) then
                Me%ExternalVar%PhytoConcentration3D => Me%ExternalVar%Concentration3D
            else
                Me%ExternalVar%PhytoConcentration1D => Me%ExternalVar%Concentration1D
            endif    

        elseif(Me%ExternalVar%PropertyID == Cohesive_Sediment_)then

            Me%ExternalVar%SPM_UnitsCoef      =  Me%UnitsCoef
            Me%ExternalVar%SPM_OK             =  .true.
            if (Me%Is3D) then
                Me%ExternalVar%SPMConcentration3D   => Me%ExternalVar%Concentration3D
            else
                Me%ExternalVar%SPMConcentration1D   => Me%ExternalVar%Concentration1D
            endif

        else
            
            write (*,*) 
            write (*,*) 'Property '//trim(GetPropertyName(Me%ExternalVar%PropertyID))
            write (*,*) 'cannot be used to compute Combined_ParsonsPortela light extinction type.'
            stop        'Compute_Combined_ParsonsPortela - ModuleLightExtinction - ERR01'

        end if

        
        if(Me%ExternalVar%Phyto_OK .and. Me%ExternalVar%SPM_OK)then

            if (Me%Is3D) then 
                PhytoConcentration3D => Me%ExternalVar%PhytoConcentration3D
                SPMConcentration3D   => Me%ExternalVar%SPMConcentration3D
                Phyto_UnitsCoef      =  Me%ExternalVar%Phyto_UnitsCoef
                SPM_UnitsCoef        =  Me%ExternalVar%SPM_UnitsCoef


                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then 

                        !factor de extincao para alto mar, funcao da concentracao de fitoplancton
                        !Eq. 59 pp 95 Parsons et al(Biol.Ocea.Proc.) C[mg Chl a/m3],tem que haver
                        !conversão de unidades pq o modelo calcula o fito em mgC/L.    

                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) =                                         &
                                        0.04 + 0.0088 * (PhytoConcentration3D(i, j, k) * Phyto_UnitsCoef) + &
                                        0.54 * (PhytoConcentration3D(i, j, k) * Phyto_UnitsCoef) ** (2.0 / 3.0) 
                                    
                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) =                                         &
                                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) * 0.7 +                   &
                                        0.036 * 0.5 * (SPMConcentration3D(i, j, k) * SPM_UnitsCoef)
                    end if

                enddo
                enddo
                enddo

            else
                
                PhytoConcentration1D => Me%ExternalVar%PhytoConcentration1D
                SPMConcentration1D   => Me%ExternalVar%SPMConcentration1D
                Phyto_UnitsCoef      =  Me%ExternalVar%Phyto_UnitsCoef
                SPM_UnitsCoef        =  Me%ExternalVar%SPM_UnitsCoef

                do i = Me%Size1D%ILB, Me%Size1D%IUB

                    if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then 

                        !factor de extincao para alto mar, funcao da concentracao de fitoplancton
                        !Eq. 59 pp 95 Parsons et al(Biol.Ocea.Proc.) C[mg Chl a/m3],tem que haver
                        !conversão de unidades pq o modelo calcula o fito em mgC/L.    

                        Me%ShortWave%ExtinctionCoefField1D(i) =                                             &
                                        0.04 + 0.0088 * (PhytoConcentration1D(i) * Phyto_UnitsCoef) +       &
                                        0.54 * (PhytoConcentration1D(i) * Phyto_UnitsCoef) ** (2.0 / 3.0) 
                                    
                        Me%ShortWave%ExtinctionCoefField1D(i) =                                             &
                                        Me%ShortWave%ExtinctionCoefField1D(i) * 0.7 +                       &
                                        0.036 * 0.5 * (SPMConcentration1D(i) * SPM_UnitsCoef)
                    end if

                enddo
    
            endif

            Me%ExternalVar%Phyto_OK = .false.; Me%ExternalVar%SPM_OK = .false.
            
            nullify(PhytoConcentration1D, SPMConcentration1D) 
            nullify(PhytoConcentration3D, SPMConcentration3D) 

        end if


    end subroutine Compute_Combined_ParsonsPortela

    !--------------------------------------------------------------------------

    subroutine Compute_Portela_Tagus
        
        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k

        !Begin-----------------------------------------------------------------


        if(Me%ExternalVar%PropertyID == Cohesive_Sediment_)then

            if (Me%Is3D) then
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then 
                
                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) = 1.24 + 0.036 * &
                                        (Me%ExternalVar%Concentration3D(i, j, k) * Me%UnitsCoef)

                    end if

                enddo
                enddo
                enddo
            else
                do i = Me%Size1D%ILB, Me%Size1D%IUB

                    if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then 
                
                        Me%ShortWave%ExtinctionCoefField1D(i) = 1.24 + 0.036 * &
                                        (Me%ExternalVar%Concentration1D(i) * Me%UnitsCoef)

                    end if

                enddo
                
            endif

        else

            write (*,*) 
            write (*,*) 'Property '//trim(GetPropertyName(Me%ExternalVar%PropertyID))
            write (*,*) 'cannot be used to compute Portela_Tagus light extinction type.'
            stop        'Compute_Portela_Tagus - ModuleLightExtinction - ERR01'

        end if


    end subroutine Compute_Portela_Tagus


    !--------------------------------------------------------------------------

    subroutine Compute_Multiparameters

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k

        !Begin-----------------------------------------------------------------

if3D:   if (Me%Is3D) then

            if(Me%ExternalVar%PropertyID .ne. MacroAlgae_)then

                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then 

                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) =                     &
                                            Me%ShortWave%ExtinctionCoefField3D(i,j,k) + &
                                            Me%ExternalVar%ExtinctionParameter        * &
                                           (Me%ExternalVar%Concentration3D(i, j, k)   * Me%UnitsCoef)
                    end if

                enddo
                enddo
                enddo

            else 
                
                do k = Me%WorkSize%KLB, Me%WorkSize%KUB
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then 

                        Me%ShortWave%ExtinctionCoefField3D(i,j,k) =                           &
                                            Me%ShortWave%ExtinctionCoefField3D(i,j,k)       + &
                                            Me%ExternalVar%ExtinctionParameter              * &
                                           (Me%ExternalVar%Concentration3D(i, j, k)         * &
                                            Me%ExternalVar%MacroAlgaeOccupation(i, j, k)    * Me%UnitsCoef)
                    end if

                enddo
                enddo
                enddo

            end if 

        else if3D

            do i = Me%Size1D%ILB, Me%Size1D%IUB

                if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then 
                
                    Me%ShortWave%ExtinctionCoefField1D(i) =                         &
                                        Me%ShortWave%ExtinctionCoefField1D(i)     + &
                                        Me%ExternalVar%ExtinctionParameter        * &
                                        (Me%ExternalVar%Concentration1D(i)        * Me%UnitsCoef)
                end if

            enddo

        endif if3D

    end subroutine Compute_Multiparameters

    !--------------------------------------------------------------------------

    subroutine ComputeLongWaveExtinction

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        logical                                 :: TimeCycle
        
        !Local-----------------------------------------------------------------
        type(T_Time)                            :: Time1,Time2
        real                                    :: Value1, Value2, NewValue
        
        !Begin-----------------------------------------------------------------

        if (Me%LongWave%ExtinctionType == ASCIIFile) then

            !Gets Value arround current instant    
            call GetTimeSerieValue(Me%ObjTimeSerie,                         &
                                   Me%ExternalVar%Now,                      &
                                   Me%LWColumn,                             &
                                   Time1, Value1, Time2, Value2, TimeCycle, &
                                   STAT = STAT_CALL)
                if (STAT_CALL/= SUCCESS_)                                   &
                    stop 'ComputeLongWaveExtinction - ModuleLightExtinction - ERR01'

            if (TimeCycle) then
                NewValue = Value1
            else
                !Interpolates Value for current instant
                call InterpolateValueInTime(Me%ExternalVar%Now, Time1,   &
                                            Value1, Time2, Value2, NewValue)
            endif

            Me%LongWave%ExtinctionCoef = NewValue

        else if (Me%LongWave%ExtinctionType == Constant) then
            
                !Do not compute any time and spatial evolution                

        else
            write (*,*) 
            write (*,*) 'Unknown type of extinction coefficient.'
            stop 'ComputeLongWaveExtinction - ModuleLightExtinction - ERR02'
        endif

    end subroutine ComputeLongWaveExtinction


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillLightExtinction(ObjLightExtinctionID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjLightExtinctionID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjLightExtinctionID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mLIGHTEXTINCTION_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !Kills Input Time Serie 
                if (Me%ShortWave%ExtinctionType == ASCIIFile .or.                       &
                    Me%LongWave%ExtinctionType  == ASCIIFile) then
                    call KillTimeSerie (Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModuleLightExtinction - KillLightExtinction - ERR01' 
                endif

                if (Me%Is3D) then
                    deallocate(Me%ShortWave%ExtinctionCoefField3D)         
                    nullify   (Me%ShortWave%ExtinctionCoefField3D) 
                else
                    deallocate(Me%ShortWave%ExtinctionCoefField1D)         
                    nullify   (Me%ShortWave%ExtinctionCoefField1D) 
                endif

                !Deallocates Instance
                call DeallocateInstance

                ObjLightExtinctionID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillLightExtinction
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_LightExtinction), pointer          :: AuxObjLightExtinction
        type (T_LightExtinction), pointer          :: PreviousObjLightExtinction

        !Updates pointers
        if (Me%InstanceID == FirstObjLightExtinction%InstanceID) then
            FirstObjLightExtinction => FirstObjLightExtinction%Next
        else
            PreviousObjLightExtinction => FirstObjLightExtinction
            AuxObjLightExtinction      => FirstObjLightExtinction%Next
            do while (AuxObjLightExtinction%InstanceID /= Me%InstanceID)
                PreviousObjLightExtinction => AuxObjLightExtinction
                AuxObjLightExtinction      => AuxObjLightExtinction%Next
            enddo

            !Now update linked list
            PreviousObjLightExtinction%Next => AuxObjLightExtinction%Next

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

    subroutine Ready (ObjLightExtinction_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjLightExtinction_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjLightExtinction_ID > 0) then
            call LocateObjLightExtinction (ObjLightExtinction_ID)
            ready_ = VerifyReadLock (mLIGHTEXTINCTION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjLightExtinction (ObjLightExtinctionID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjLightExtinctionID

        !Local-----------------------------------------------------------------

        Me => FirstObjLightExtinction
        do while (associated (Me))
            if (Me%InstanceID == ObjLightExtinctionID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleLightExtinction - LocateObjLightExtinction - ERR01'

    end subroutine LocateObjLightExtinction

    !--------------------------------------------------------------------------

end module ModuleLightExtinction

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------






