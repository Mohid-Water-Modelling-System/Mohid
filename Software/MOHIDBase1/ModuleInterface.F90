!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : Interface 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as link to zero-dimensional sinks and sources modules
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

Module ModuleInterface

    use ModuleGlobalData
    use ModuleTime
    use ModuleWaterQuality
    use ModuleSedimentQuality
    use ModuleCEQUALW2
    use ModuleLife
    use ModuleBenthos
    use ModuleMacroAlgae
    use ModuleEnterData, only: ReadFileName

#ifdef _PHREEQC_ 
    use ModulePhreeqCData   
    use ModulePhreeqC
#endif   

#ifdef _BFM_    
    use ModuleBFM
#endif

    implicit none

    private

    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  :: ConstructInterface
    private ::      AllocateInstance
    private ::      ReadInterfaceFilesName
    private ::      StartSinksSourcesModel
    private ::      AllocateVariables                                                
    private ::      Check_Options
    private ::          FindProperty

    !Modifier 
    public  :: Modify_Interface
    private ::      FillMassTempSalinity
    private ::      InputData
    private ::      PropertyIndexNumber
    private ::      UnfoldMatrix
    public  :: SetSOD

    !Selector
    public  :: GetRateFlux

    public  :: GetWQRatio
    
#ifdef _PHREEQC_       
    public  :: GetPhreeqCID
#endif

    !Destructor 
    public  :: KillInterface
    private ::      DeallocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjInterface

    !Interfaces
    private :: ConstructInterface1D
    private :: ConstructInterface2D
    private :: ConstructInterface3D
    interface  ConstructInterface
        module procedure ConstructInterface1D
        module procedure ConstructInterface2D
        module procedure ConstructInterface3D
    end interface  ConstructInterface


    private :: Modify_Interface1D
    private :: Modify_Interface2D
    private :: Modify_Interface3D
    interface  Modify_Interface
        module procedure Modify_Interface1D
        module procedure Modify_Interface2D
        module procedure Modify_Interface3D
    end interface  Modify_Interface

    
    private :: FillMassTempSalinity1D
    private :: FillMassTempSalinity2D
    private :: FillMassTempSalinity3D
    interface  FillMassTempSalinity
        module procedure FillMassTempSalinity1D
        module procedure FillMassTempSalinity2D
        module procedure FillMassTempSalinity3D
    end interface  FillMassTempSalinity

    private :: InputData1D
    private :: InputData2D
    private :: InputData3D
    interface  InputData
        module procedure InputData1D
        module procedure InputData2D
        module procedure InputData3D
    end interface  InputData

    private :: GetRateFlux2D
    private :: GetRateFlux3D
    interface  GetRateFlux
        module procedure GetRateFlux2D
        module procedure GetRateFlux3D
    end interface  GetRateFlux
    
    private :: UnfoldMatrix1D_I
    private :: UnfoldMatrix1D_R
    private :: UnfoldMatrix2D_I
    private :: UnfoldMatrix2D_R
    private :: UnfoldMatrix3D_I
    private :: UnfoldMatrix3D_R
    interface  UnfoldMatrix
        module procedure UnfoldMatrix1D_I
        module procedure UnfoldMatrix1D_R
        module procedure UnfoldMatrix2D_I
        module procedure UnfoldMatrix2D_R
        module procedure UnfoldMatrix3D_I
        module procedure UnfoldMatrix3D_R
    end interface  UnfoldMatrix


    !Type----------------------------------------------------------------------
    type       T_External
        type(T_Time     )                       :: Now
        integer, dimension(:      ), pointer    :: RiverPoints1D
        integer, dimension(:, :   ), pointer    :: WaterPoints2D
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        integer, dimension(:      ), pointer    :: OpenPoints1D
        integer, dimension(:, :   ), pointer    :: OpenPoints2D
        integer, dimension(:, :, :), pointer    :: OpenPoints3D
        real,    dimension(:      ), pointer    :: DWZ1D
        real,    dimension(:, :, :), pointer    :: DWZ, ShearStress, SPMFlux
        logical                                 :: Vertical1D = .false.
    end type T_External

    type      T_Interface
        private
        integer                                 :: InstanceID
        character(PathLength)                   :: FileName
        character(StringLength)                 :: SinksSourcesModel
        type(T_Size3D)                          :: Size3D
        type(T_Size2D)                          :: Size2D
        type(T_Size1D)                          :: Size1D
        type(T_Size1D)                          :: Array
        type(T_Size1D)                          :: Prop
        type(T_External)                        :: ExternalVar
        real,    pointer, dimension(:,:  )      :: Mass
        real,    pointer, dimension(:,:  )      :: ConcentrationIncrement
        real,    pointer, dimension(:    )      :: Salinity
        real,    pointer, dimension(:    )      :: Alkalinity
        real,    pointer, dimension(:    )      :: Temperature
        real,    pointer, dimension(:    )      :: Oxygen
        real,    pointer, dimension(:    )      :: SOD
        logical                                 :: UseSOD = .false.
        real,    pointer, dimension(:    )      :: ShortWaveTop
        real,    pointer, dimension(:    )      :: ShortWaveAverage
        real,    pointer, dimension(:    )      :: LightExtCoefField
        real,    pointer, dimension(:    )      :: Thickness
        real,    pointer, dimension(:    )      :: SPMFlux
        real,    pointer, dimension(:    )      :: ShearStress
        real,    pointer, dimension(:    )      :: FishFood
        real,    pointer, dimension(:    )      :: WaterPercentage
        real,    pointer, dimension(:    )      :: DissolvedToParticulate
        real,    pointer, dimension(:    )      :: SoilDryDensity
        real,    pointer, dimension(:    )      :: pH
#ifdef _PHREEQC_        
        real,    pointer, dimension(:    )      :: pE
        real,    pointer, dimension(:    )      :: SolutionVolume  
        real,    pointer, dimension(:    )      :: SolidMass   
#endif        
        real,    pointer, dimension(:    )      :: IonicStrength
        real,    pointer, dimension(:    )      :: PhosphorusAdsortionIndex
        real,    pointer, dimension(:    )      :: WindVelocity
        integer, pointer, dimension(:    )      :: OpenPoints
        logical, pointer, dimension(:    )      :: AddedProperties

        !Instance of ModuleTime
        integer                                 :: ObjTime              = 0

        !Instance of ModuleWaterQuality         
        integer                                 :: ObjWaterQuality      = 0

        !Instance of ModuleSedimentQuality      
        integer                                 :: ObjSedimentQuality   = 0

        !Instance of ModuleCEQUALW2
        integer                                 :: ObjCEQUALW2          = 0

        !Instance of ModuleBenthos
        integer                                 :: ObjBenthos           = 0

        !Instance of ModuleLife
        integer                                 :: ObjLife              = 0
#ifdef _BFM_    
        !Instance of ModuleBFM
        integer                                 :: ObjBFM               = 0
#endif    
        !Instance of ModuleMacroAlgae
        integer                                 :: ObjMacroAlgae        = 0

#ifdef _PHREEQC_
        !Instance of ModulePhreeqC
        integer                                 :: ObjPhreeqC           = 0
#endif   

        !Collection of instances                
        type(T_Interface),          pointer     :: Next        

    end type T_Interface

        
    !Global Module Variables
    type (T_Interface), pointer                 :: FirstObjInterface
    type (T_Interface), pointer                 :: Me

    !--------------------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructInterface3D(InterfaceID,                           &
                                    TimeID,                                &
                                    SinksSourcesModel,                     &
                                    DT,PropertiesList,                     &
                                    WaterPoints3D,                         &
                                    Size3D,                                &
                                    Vertical1D,                            &
                                    STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: InterfaceID
        integer                                                 :: TimeID
        character(len=StringLength)                             :: SinksSourcesModel
        integer, dimension(:    ), pointer                      :: PropertiesList
        real,   intent (OUT)                                    :: DT
        integer,                     dimension(:,:,:), pointer  :: WaterPoints3D
        type(T_Size3D)                                          :: Size3D
        logical,intent (IN),  optional                          :: Vertical1D
        integer,intent (OUT), optional                          :: STAT     

        !External--------------------------------------------------------------
        integer                                                 :: ready_

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterface_)) then
            nullify (FirstObjInterface)
            call RegisterModule (mInterface_) 
        endif

        call Ready(InterfaceID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            !Allocates a new Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTime                      =  AssociateInstance (mTIME_, TimeID)
            Me%SinksSourcesModel            =  SinksSourcesModel
            Me%ExternalVar%WaterPoints3D    => WaterPoints3D
            Me%Size3D                       =  Size3D 
            Me%Array%ILB                    = 1
            Me%Array%IUB                    = sum(Me%ExternalVar%WaterPoints3D)

            if (present(Vertical1D)) then
                Me%ExternalVar%Vertical1D = Vertical1D
            else
                Me%ExternalVar%Vertical1D = .false.
            endif
            

            !Path to data file
            call ReadInterfaceFilesName

            !Start sinks and sources model
            call StartSinksSourcesModel(DT)

            !Verify compute options
            call Check_Options(PropertiesList)

            !Allocate variables global to whole module
            call AllocateVariables

            STAT_ = SUCCESS_
            
            !Returns ID
            InterfaceID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleInterface - ConstructInterface3D - ERR01' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructInterface3D

    !--------------------------------------------------------------------------

    subroutine ConstructInterface2D(InterfaceID,                           &
                                    TimeID,                                &
                                    SinksSourcesModel,                     &
                                    DT,PropertiesList,                     &
                                    WaterPoints2D,                         &
                                    Size2D,                                &
                                    STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: InterfaceID
        integer                                                 :: TimeID
        character(len=StringLength)                             :: SinksSourcesModel
        integer, dimension(:    ), pointer                      :: PropertiesList
        real,   intent (OUT)                                    :: DT
        integer,                     dimension(:,:  ), pointer  :: WaterPoints2D
        type(T_Size2D)                                          :: Size2D
        integer,intent (OUT), optional                          :: STAT     

        !External--------------------------------------------------------------
        integer                                                 :: ready_

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterface_)) then
            nullify (FirstObjInterface)
            call RegisterModule (mInterface_) 
        endif


        call Ready(InterfaceID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            !Allocates a new Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTime                      =  AssociateInstance (mTIME_, TimeID)
            Me%SinksSourcesModel            =  SinksSourcesModel
            Me%ExternalVar%WaterPoints2D    => WaterPoints2D
            Me%Size2D                       =  Size2D
            Me%Array%ILB                    = 1
            Me%Array%IUB                    = sum(Me%ExternalVar%WaterPoints2D)

            !Path to data file
            call ReadInterfaceFilesName

            !Start sinks and sources model
            call StartSinksSourcesModel(DT)
 
            !Verify compute options
            call Check_Options(PropertiesList)

            !Allocate variables global to whole module
            call AllocateVariables

            STAT_ = SUCCESS_
            
            !Returns ID
            InterfaceID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleInterface - ConstructInterface2D - ERR01' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructInterface2D

    !--------------------------------------------------------------------------

    subroutine ConstructInterface1D(InterfaceID,                           &
                                    TimeID,                                &
                                    SinksSourcesModel,                     &
                                    DT,PropertiesList,                     &
                                    RiverPoints1D,                         &
                                    Size1D,                                &
                                    STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: InterfaceID
        integer                                                 :: TimeID
        character(len=StringLength)                             :: SinksSourcesModel
        integer, dimension(:    ), pointer                      :: PropertiesList
        real,   intent (OUT)                                    :: DT
        integer, dimension(:), pointer                          :: RiverPoints1D
        type(T_Size1D)                                          :: Size1D
        integer,intent (OUT), optional                          :: STAT     

        !External--------------------------------------------------------------
        integer                                                 :: ready_

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mInterface_)) then
            nullify (FirstObjInterface)
            call RegisterModule (mInterface_) 
        endif

        call Ready(InterfaceID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            !Allocates a new Instance
            call AllocateInstance

            !Associates External Instances
            Me%ObjTime                      =  AssociateInstance (mTIME_, TimeID)
            Me%SinksSourcesModel            =  SinksSourcesModel
            Me%ExternalVar%RiverPoints1D    => RiverPoints1D
            Me%Size1D                       =  Size1D
            Me%Array%ILB                    = 1
            Me%Array%IUB                    = sum(Me%ExternalVar%RiverPoints1D)

            !Path to data file
            call ReadInterfaceFilesName

            !Start sinks and sources model
            call StartSinksSourcesModel(DT)

            !Verify compute options
            call Check_Options(PropertiesList)

            !Allocate variables global to whole module
            call AllocateVariables

            STAT_ = SUCCESS_
            
            !Returns ID
            InterfaceID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleInterface - ConstructInterface2D - ERR01' 

        end if cd0

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructInterface1D

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_Interface), pointer           :: NewInterface
        type (T_Interface), pointer           :: PreviousInterface

        !Allocates new instance
        allocate (NewInterface)
        nullify  (NewInterface%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjInterface)) then
            FirstObjInterface      => NewInterface
            Me                     => NewInterface
        else
            PreviousInterface      => FirstObjInterface
            Me                     => FirstObjInterface%Next
            do while (associated(Me))
                PreviousInterface  => Me
                Me                 => Me%Next
            enddo
            Me                     => NewInterface
            PreviousInterface%Next => NewInterface
        endif

        Me%InstanceID = RegisterNewInstance (mINTERFACE_)

    end subroutine AllocateInstance


    !--------------------------------------------------------------------------

    
    logical function  FindProperty(PropertiesList, PropertyID)

        !Arguments-------------------------------------------------------------
        integer, pointer, dimension(:)              :: PropertiesList
        integer, intent (IN)                        :: PropertyID

        !Local-----------------------------------------------------------------
        integer                                     :: i
        
        !Begin-----------------------------------------------------------------
        
        FindProperty = .false.

        do i = 1, size(PropertiesList)  
            if (PropertyID == PropertiesList(i))then
                FindProperty =.true.
                exit
            end if
        enddo

        !----------------------------------------------------------------------

    end function FindProperty


    !--------------------------------------------------------------------------

    
    subroutine AllocateVariables

        !External--------------------------------------------------------------
        integer                       :: STAT_CALL

        !Local-----------------------------------------------------------------
        integer                       :: ArrayLB, ArrayUB
        integer                       :: PropLB, PropUB

        !----------------------------------------------------------------------

        PropLB  = Me%Prop%ILB
        PropUB  = Me%Prop%IUB

        ArrayLB = Me%Array%ILB
        ArrayUB = Me%Array%IUB

        
        !Allocate global variables--------------------------------------------- 
        allocate(Me%Mass(PropLB:PropUB, ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR01'

        allocate(Me%ConcentrationIncrement(PropLB:PropUB, ArrayLB:ArrayUB), &
                                                     STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR02'

        allocate(Me%Temperature(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR03'

        allocate(Me%OpenPoints(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR04'
        
        allocate(Me%AddedProperties(PropLB:PropUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR05'
        
        
        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR07'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR08'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR09'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR10'

                allocate(Me%FishFood(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR11'
                       
                Me%Salinity           = FillValueReal
                Me%FishFood           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal


            case (SedimentQualityModel)
                
                allocate (Me%WaterPercentage(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR12'
                
                Me%WaterPercentage = FillValueReal

                allocate (Me%DissolvedToParticulate(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR13'
                
                Me%DissolvedToParticulate = FillValueReal
         
                allocate (Me%SoilDryDensity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR14'
                
                Me%SoilDryDensity = FillValueReal

                allocate (Me%Salinity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR15'
                
                Me%Salinity = FillValueReal

                allocate (Me%pH(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR15'
                
                Me%pH = FillValueReal

                allocate (Me%IonicStrength(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR15'
                
                Me%IonicStrength = FillValueReal

                allocate (Me%PhosphorusAdsortionIndex(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR15'
                
                Me%PhosphorusAdsortionIndex = FillValueReal

                allocate (Me%WindVelocity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR15'
                
                Me%WindVelocity = FillValueReal

            case (CEQUALW2Model)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR14'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR15'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR16'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR17'

                allocate(Me%Alkalinity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR18'
                       
                Me%Salinity           = FillValueReal
                Me%Alkalinity         = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal

            case (LifeModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR19'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR20'

                allocate(Me%ShortWaveAverage(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR20a'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR21'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR22'

                allocate(Me%Alkalinity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR23'
                       
                Me%Salinity           = FillValueReal
                Me%Alkalinity         = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveAverage   = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
                
#ifdef _BFM_    
            case (BFMModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR19a'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR20a'
        
                allocate(Me%ShortWaveRadiation(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR21a'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR22a'

                Me%Salinity           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
#endif

            case (BenthicCEQUALW2Model)                
                allocate(Me%Oxygen(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR24'
                
                allocate(Me%SOD (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR25'
            
            case (BenthosModel )
            
                allocate(Me%Oxygen(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR24'

            case (MacroAlgaeModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR30'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR40'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR50'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR60'
              
                allocate(Me%ShearStress(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR70'
                
                allocate(Me%SPMFlux(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR80'

                Me%Salinity           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
                Me%ShearStress        = FillValueReal
                Me%SPMFlux           = FillValueReal


#ifdef _PHREEQC_
            case (PhreeqCModel)
                                               
                allocate (Me%SolutionVolume(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR26'
                
                allocate (Me%pH(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR27'
                
                allocate (Me%pE(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR28'

                allocate (Me%Temperature(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR29B'
                
                allocate (Me%SolidMass(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR29C'                

                Me%SolutionVolume  = FillValueReal
                Me%pH              = FillValueReal
                Me%pE              = FillValueReal
                Me%Temperature     = FillValueReal
                Me%SolidMass       = FillValueReal
#endif

            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'AllocateVariables - ModuleInterface - ERR90'
        end select

        Me%Mass                       = FillValueReal
        Me%ConcentrationIncrement     = FillValueReal
        Me%Temperature                = FillValueReal
        Me%OpenPoints                 = FillValueInt
        Me%AddedProperties            = .false.

        !----------------------------------------------------------------------

    end subroutine AllocateVariables   
    

    !--------------------------------------------------------------------------


    subroutine ReadInterfaceFilesName 

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        character(Len = StringLength)               :: Message

        !Begin-----------------------------------------------------------------
        

        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)

                Message  =trim('Water Quality Data File')

                call ReadFileName('WQDATA',Me%FileName, Message = Message, STAT = STAT_CALL)

cd1 :           if(STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then

                    call ReadFileName('DISPQUAL',                    &
                                      Me%FileName,                   &
                                      Message = Message,             &
                                      STAT    = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - ModuleInterface - ERR01' 

                else if (STAT_CALL .EQ. SUCCESS_) then cd1
                    continue
                else
                    stop 'ReadInterfaceFilesName - ModuleInterface - ERR02' 
                end if cd1


            case (SedimentQualityModel)

                Message  = trim('SedimentQuality Data File')
                call ReadFileName('SQ_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR03' 
            

            case(CEQUALW2Model)
                
                Message  = trim('CEQUALW2 Data File')
                call ReadFileName('CEQUALW2_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR04' 
            

            case(LifeModel)
                
                Message  = trim('Life Data File')
                call ReadFileName('LIFE_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR05' 
#ifdef _BFM_    
            case(BFMModel)
                
                Message  = trim('BFM Data File')
                call ReadFileName('BFM_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR05a' 
#endif   

            case(BenthicCEQUALW2Model)

                Message  = trim('Benthic CEQUALW2 Data File')
                call ReadFileName('BENTHIC_CEQUALW2_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR06' 

            case(BenthosModel)

                Message  = trim('Benthos Data File')
                call ReadFileName('BENTHOS_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR07' 
            
            case(MacroAlgaeModel)

                Message  = trim('MacroAlgae Data File')
                call ReadFileName('MACROALGAE_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR08' 

#ifdef _PHREEQC_
                
            case (PhreeqCModel)
            
                Message = trim('PhreeqC Data File')
                call ReadFileName('PHREEQC_DATA', Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR10' 
            
#endif
            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'ReadInterfaceFilesName - ModuleInterface - ERR09' 
        end select


        Me%FileName = trim(Me%FileName)

        !----------------------------------------------------------------------

    end subroutine ReadInterfaceFilesName 

    !--------------------------------------------------------------------------
    
    subroutine StartSinksSourcesModel (DT)

        !Arguments-------------------------------------------------------------
        real, intent(OUT)                       :: DT
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: PropLB, PropUB

        !Begin-----------------------------------------------------------------

        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)

                !Construct WaterQuality Model
                call StartWaterQuality(Me%ObjWaterQuality,                     &
                                       Me%FileName,                            &
                                       STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR01'
 
                !Construct mass fluxes between properties
                call Construct_WQRateFlux(Me%ObjWaterQuality,                  &
                                          Me%Array%ILB,                        &
                                          Me%Array%IUB,                        &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR02'
               
                !Number of properties involved
                call GetWaterQualitySize(Me%ObjWaterQuality,                   &
                                         PropLB = PropLB,                      &
                                         PropUB = PropUB,                      &
                                         STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR03'

                !Number of properties involved
                Me%Prop%ILB     = PropLB
                Me%Prop%IUB     = PropUB

                !Get DT from WaterQuality model to exit as argument to WaterProperties
                call GetDTWQM(Me%ObjWaterQuality, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR04'


            case (SedimentQualityModel)

                !Construct SedimentQuality Model
                call StartSedimentQuality(Me%ObjSedimentQuality,               &
                                          Me%FileName,                         &
                                          STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR05'
                
                !Construct mass fluxes between properties
                call Construct_SQRateFlux(Me%ObjSedimentQuality,               &
                                          Me%Array%ILB,                        &
                                          Me%Array%IUB,                        &
                                         STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR06'
                
                !Number of properties involved
                call GetSedimentQualitySize(Me%ObjSedimentQuality,             &
                                            PropLB = PropLB,                   &
                                            PropUB = PropUB,                   &
                                            STAT   = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR07'

                
                !Number of properties involved
                Me%Prop%ILB     = PropLB
                Me%Prop%IUB     = PropUB

                !Get DT from SedimentQuality model to exit as argument to SedimentProperties
                call GetDTSedimentQuality(Me%ObjSedimentQuality, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR08'


            case(LifeModel)

                !Construct Life Model
                call ConstructLife(Me%ObjLife,                                  &
                                   FileName = Me%FileName,                      &
                                   STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR09'
 
                !Number of properties involved
                call GetLifeSize(Me%ObjLife,                                    &
                                     PropLB = PropLB,                           &
                                     PropUB = PropUB,                           &
                                     STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR010'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from Life model to exit as argument to WaterProperties
                call GetDTLife(Me%ObjLife, DT = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR011'
#ifdef _BFM_    
            case(BFMModel)

                !Construct BFM Model
                call ConstructBFM(Me%ObjBFM,                                    &
                                  FileName  = Me%FileName,                      &
                                  ILB       = Me%Array%ILB,                     &
                                  IUB       = Me%Array%IUB,                     &
                                  STAT      = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR09a'
 
                !Number of properties involved
                call GetBFMSize(Me%ObjBFM,                                      &
                                PropLB = PropLB,                                &
                                PropUB = PropUB,                                &
                                STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR010a'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from BFM model to exit as argument to WaterProperties
                call GetDTBFM(Me%ObjBFM, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR011a'
#endif   
            case(CEQUALW2Model)

                !Construct CEQUALW2 Model
                call StartCEQUALW2(Me%ObjCEQUALW2,                              &
                                   Me%Array%ILB,                                &
                                   Me%Array%IUB,                                &
                                   FileName = Me%FileName,                      &
                                   Model = CEQUALW2Model,                       & 
                                   STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR09'
 
                !Number of properties involved
                call GetCEQUALW2Size(Me%ObjCEQUALW2,                            &
                                     PropLB = PropLB,                           &
                                     PropUB = PropUB,                           &
                                     STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR010'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from CEQUALW2 model to exit as argument to WaterProperties
                call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR011'


            case(BenthicCEQUALW2Model)
                
                !Construct CEQUALW2 Model
                call StartCEQUALW2(Me%ObjCEQUALW2,                              &
                                   Me%Array%ILB,                                &
                                   Me%Array%IUB,                                &
                                   FileName = Me%FileName,                      &
                                   Model = BenthicCEQUALW2Model,                & 
                                   STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR012'
               
                !Number of properties involved
                call GetCEQUALW2Size(Me%ObjCEQUALW2,                            &
                                     PropLB = PropLB,                           &
                                     PropUB = PropUB,                           &
                                     STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR012'
                 
                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from CEQUALW2 model to exit as argument to WaterProperties
                call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR014'
            
            case(BenthosModel)


                !Construct Benthos Model
                call StartBenthos (Me%ObjBenthos,                               &
                                   FileName = Me%FileName,                      &
                                   ILB      = Me%Array%ILB,                     &
                                   IUB      = Me%Array%IUB,                     &
                                   STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR015'
 
                !Number of properties involved
                call GetBenthosSize (Me%ObjBenthos,                             &
                                     PropLB = PropLB,                           &
                                     PropUB = PropUB,                           &
                                     STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR016'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from CEQUALW2 model to exit as argument to WaterProperties
                call GetDTBenthos(Me%ObjBenthos, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR017'
            
            case(MacroAlgaeModel)

                !Construct MacroAlgae Model
                call StartMacroAlgae (Me%ObjMacroAlgae,                                 &
                                      FileName = Me%FileName,                           &
                                      ILB      = Me%Array%ILB,                          &
                                      IUB      = Me%Array%IUB,                          &
                                      STAT     = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR018'
 
                !Number of properties involved
                call GetMacroAlgaeSize (Me%ObjMacroAlgae,                               &
                                        PropLB = PropLB,                                &
                                        PropUB = PropUB,                                &
                                        STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR019'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from MacroAlgae model to exit as argument to WaterProperties
                call GetDTMacroAlgae(Me%ObjMacroAlgae, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR020'

#ifdef _PHREEQC_

            case (PhreeqCModel)
            
                !Construct PhreeqC Model
                call StartPhreeqC (Me%ObjPhreeqC, FileName = Me%FileName, STAT = STAT_CALL)
                    
                !ToDo: Because I don't know the number of properties, maybe this is inappropriated               
                !Get number of properties involved
!                call GetPhreeqCSize(Me%ObjPhreeqC,   &
!                                    PropLB = PropLB, &
!                                    PropUB = PropUB, &
!                                    STAT   = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR21'
!                
!                !Store number of properties involved
!                Me%Prop%ILB = PropLB
!                Me%Prop%IUB = PropUB

                call GetPhreeqCDT(Me%ObjPhreeqC, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR22'
                
#endif
            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'StartSinksSourcesModel - ModuleInterface - ERR00' 
        end select

        !----------------------------------------------------------------------


    end subroutine StartSinksSourcesModel    
    
    
    !--------------------------------------------------------------------------
     

    subroutine Check_Options(PropertiesList)

        !Arguments-------------------------------------------------------------
        integer, dimension(:), pointer                      :: PropertiesList   

        !External--------------------------------------------------------------
        type(T_Time)                                         :: EndTime
        logical                                              :: Zoo, Phyto
        logical                                              :: Diatoms 
        logical                                              :: Nitrogen, Phosphorus
        logical                                              :: Silica  
        logical                                              :: Oxygen, BOD
        logical                                              :: Carbon, Sol_Bacteria
        logical                                              :: Bacteria, Ciliate 
        logical                                              :: Larvae
        logical                                              :: Pompools
        real                                                 :: DT
        integer                                              :: STAT_CALL

        !Local-----------------------------------------------------------------
        real                                                 :: ErrorAux, auxFactor 
        real                                                 :: RunPeriod, Dtlag
        integer, dimension(:), pointer                       :: CEQUALW2List
        integer, dimension(:), pointer                       :: MacroAlgaeList
        integer                                              :: i,PropLB, PropUB
        integer, dimension(:), pointer                       :: BenthosList, LifeList
#ifdef _BFM_  
        integer, dimension(:), pointer                       :: BFMList
#endif
        !----------------------------------------------------------------------
        
        call GetComputeCurrentTime(Me%ObjTime,                  &
                                   Me%ExternalVar%Now,          &
                                   STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR00' 

        !Get end time
        call GetComputeTimeLimits(Me%ObjTime, EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR01' 
            
        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)
                !Get water quality model time step
                call GetDTWQM(Me%ObjWaterQuality, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Check_Options - ModuleInterface - ERR02' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the WQ DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'Water Quality wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN01.'
                endif 

                !Sinks and sources compute options 
                call GetWQOptions(Me%ObjWaterQuality,   Phyto            = Phyto,           &
                                                        Nitrogen         = Nitrogen,        &
                                                        Diatoms          = Diatoms,         &
                                                        Silica           = Silica,          & 
                                                        Phosphorus       = Phosphorus,      &
                                                        Zoo              = Zoo,             &
                                                        BOD              = BOD,             &
                                                        Oxygen           = Oxygen,          &
                                                        Bacteria         = Bacteria,        &
                                                        Ciliate          = Ciliate,         &
                                                        Larvae           = Larvae,          &
                                                        Pompools         = Pompools,        &
                                                        STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR03' 


cd1 :           if (Phyto) then 
                    if (.not.FindProperty(PropertiesList, Phytoplankton_))          &
                        stop 'WQM needs property phytoplankton - Check_Options'
                end if cd1

cd2 :           if (Zoo) then
                    if (.not.FindProperty(PropertiesList, Zooplankton_))            &
                        stop 'WQM needs property zooplankton - Check_Options'
                end if cd2

cd3 :           if (Larvae) then
                    if (.not.FindProperty(PropertiesList, Larvae_))                 &
                        stop 'WQM needs property Larvae - Check_Options'
                end if cd3

cd5 :           if (Phosphorus) then

                    if (.not.FindProperty(PropertiesList, POP_))                    &
                        stop 'WQM needs property particulate organic phosphorus - Check_Options'

                    if (.not.FindProperty(PropertiesList, DOPRefractory_))          &
                        stop 'WQM needs property dissolved refractory organic phosphorus - Check_Options'

                    if (.not. FindProperty(PropertiesList, DOPNon_Refractory_))     &
                        stop 'WQM needs property dissolved non-refractory organic phosphorus - Check_Options'

                    if (.not.FindProperty(PropertiesList, Inorganic_Phosphorus_))   &
                        stop 'WQM needs property inorganic phosphorus - Check_Options'
                end if cd5

cd6 :           if (Nitrogen) then
                    if (.not.FindProperty(PropertiesList, PON_))                    &
                        stop 'WQM needs property particulate organic nitrogen - Check_Options'
                        
                    if (Bacteria .and. .not.FindProperty(PropertiesList, PONRefractory_))          &
                        stop 'WQM needs property particulate refractory organic nitrogen - Check_Options'

                    if (.not.FindProperty(PropertiesList, DONRefractory_))          &
                        stop 'WQM needs property dissolved refractory organic nitrogen - Check_Options'

                    if (.not. FindProperty(PropertiesList, DONNon_Refractory_))     &
                        stop 'WQM needs property dissolved non-refractory organic nitrogen - Check_Options'

                    if (.not.FindProperty(PropertiesList, Ammonia_))                &
                        stop 'WQM needs property ammonia - Check_Options'

                    if (.not.FindProperty(PropertiesList, Nitrate_))                &
                        stop 'WQM needs property nitrate - Check_Options'

                    if (.not.FindProperty(PropertiesList, Nitrite_))                &
                        stop 'WQM needs property nitrite - Check_Options'
                end if cd6

cd7 :           if (BOD) then
                    if (.not.FindProperty(PropertiesList, BOD_))                    &
                        stop 'WQM needs property biochemical oxygen demand - Check_Options'
                end if cd7


cd8 :           if (Bacteria) then
                    if (.not.FindProperty(PropertiesList, Bacteria_))               &
                        stop 'WQM needs property bacteria - Check_Options'
                end if cd8


cd9 :           if (Ciliate) then
                    if (.not.FindProperty(PropertiesList, Ciliate_))                &
                        stop 'WQM needs property ciliate - Check_Options'
                end if cd9


cd12 :          if (Diatoms) then
                    if (.not.FindProperty(PropertiesList, Diatoms_))                &
                        stop 'WQM needs property diatoms - Check_Options'
                end if cd12

cd13 :          if (Silica) then
                    if (.not.FindProperty(PropertiesList, DSilica_))                &
                        stop 'WQM needs property diatoms - Check_Options'

                    if (.not.FindProperty(PropertiesList, BioSilica_))               &
                        stop 'WQM needs property diatoms - Check_Options'
                end if cd13

                !If one wants to run just the age property oxygen is not needed - Frank 09-2003
                if (Phosphorus .or. BOD .or. Ciliate .or. Nitrogen .or. Bacteria) then
                    if (.not.FindProperty(PropertiesList, Oxygen_))                     &
                        stop 'WQM needs property oxygen - Check_Options'
                endif

cd50 :         if(.NOT.(Pompools)) then

                if (FindProperty(PropertiesList, PON1_) .OR. FindProperty(PropertiesList, PON2_) .OR. &
                    FindProperty(PropertiesList, PON3_) .OR. FindProperty(PropertiesList, PON4_) .OR. &
                    FindProperty(PropertiesList, PON5_) .OR. FindProperty(PropertiesList, POP1_) .OR. &
                    FindProperty(PropertiesList, POP2_) .OR. FindProperty(PropertiesList, POP3_) .OR. &
                    FindProperty(PropertiesList, POP4_) .OR. FindProperty(PropertiesList, POP5_))     &
                        stop 'WQM needs the POM pools option activated - Check_Options'

               end if cd50

cd60 :         if (Pompools) then
                 
                 if (Nitrogen) then
                    if (.not.FindProperty(PropertiesList, PON1_))                    &
                        stop 'WQM with POM pools needs property PON1 - Check_Options'
                    
                    if (.not.FindProperty(PropertiesList, PON2_))                    &
                        stop 'WQM with POM pools needs property PON2 - Check_Options'    
                    
                    if (.not.FindProperty(PropertiesList, PON3_))                    &
                        stop 'WQM with POM pools needs property PON3 - Check_Options'
                    
                    if (.not.FindProperty(PropertiesList, PON4_))                    &
                        stop 'WQM with POM pools needs property PON4 - Check_Options'
                        
                    if (.not.FindProperty(PropertiesList, PON5_))                    &
                        stop 'WQM with POM pools needs property PON5 - Check_Options'
                            
                  end if
                  
                  if (Phosphorus) then
                     
                     if (.not.FindProperty(PropertiesList, POP1_))                    &
                        stop 'WQM with POM pools needs property POP1 - Check_Options'
                    
                    if (.not.FindProperty(PropertiesList, POP2_))                    &
                        stop 'WQM with POM pools needs property POP2 - Check_Options'    
                    
                    if (.not.FindProperty(PropertiesList, POP3_))                    &
                        stop 'WQM with POM pools needs property POP3 - Check_Options'
                    
                    if (.not.FindProperty(PropertiesList, POP4_))                    &
                        stop 'WQM with POM pools needs property POP4 - Check_Options'
                        
                    if (.not.FindProperty(PropertiesList, POP5_))                    &
                        stop 'WQM with POM pools needs property POP5 - Check_Options'
                        
                  end if
                  
                  if (.NOT. Zoo) then
                    stop 'WQM with POM pools needs property Zoo - Check_Options'
                  end if 
                  
                end if cd60
            
            case (SedimentQualityModel)

                !Get SedimentQuality model time step
                call GetDTSedimentQuality(Me%ObjSedimentQuality, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR04' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'SedimentQuality wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN02.'
                endif
                
              !Sinks and sources compute options 
                call GetSQOptions(Me%ObjSedimentQuality,  Nitrogen    = Nitrogen  ,   &
                                                          Carbon      = Carbon    ,   &
                                                          Phosphorus  = Phosphorus,   &
                                                          Sol_Bacteria= Sol_Bacteria, &
                                                          Oxygen      = Oxygen    ,   &
                                                          STAT        = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR05' 

cd10 :           if (Nitrogen) then
                   
                    if (.not.FindProperty(PropertiesList, Ammonia_)              )   &
                        stop 'SQM needs property "ammonia" - Check_Options'

                    if (.not.FindProperty(PropertiesList, Nitrate_)              )   &
                        stop 'SQM needs property "nitrate" - Check_Options'

                    !if (.not.FindProperty(PropertiesList, AdsorbedAmmonia_)     )   &
                     !   STOP 'WQM needs property Adsorbed Ammonia - Check_Options'
                    
                    if (.not.FindProperty(PropertiesList, RefreactaryOrganicN_)  )   &
                        stop 'SQM needs property "particulated refractory organic nitrogen" - Check_Options'

                    if (.not.FindProperty(PropertiesList, PON_           )       )   &
                        stop 'SQM needs property "particulate organic nitrogen" - Check_Options'

                    if (.not.FindProperty(PropertiesList, Ngas_)                 )   &
                        stop 'SQM needs property "nitrogen gas" - Check_Options'

                    if (.not.FindProperty(PropertiesList, HeterotrophicN_)       )   &
                        stop 'SQM needs property "heterotrophic microorganism nitrogen" - Check_Options'

                    if (.not.FindProperty(PropertiesList, AnaerobicN_)           )   &
                        stop 'SQM needs property "anaerobic microorganism nitrogen" - Check_Options'

                    if (.not.FindProperty(PropertiesList, AutotrophicN_)         )   &
                        STOP 'SQM needs property "autotrophic microorganism nitrogen" - Check_Options'

                    if (.not.FindProperty(PropertiesList, Urea_)                 )   &
                        STOP 'SQM needs property "urea" - Check_Options'

!                    if (.not.FindProperty(PropertiesList, AmmoniaGas_)          )   &
!                        STOP 'SQM needs property "ammonia gas" - Check_Options'
                    if (Sol_Bacteria) then
                        if (.not.FindProperty(PropertiesList, SolubilizingN_)    )   &
                            STOP 'SQM needs property "solubilizing microorganism nitrogen" - Check_Options'
                    end if

                end if cd10

cd11 :          if (Carbon) then
                   
                    if (.not.FindProperty(PropertiesList, LabileOrganicC_)       )   &
                        stop 'SQM needs property "particulate labile organic carbon" - Check_Options'

                    if (.not.FindProperty(PropertiesList, RefreactaryOrganicC_)  )   &
                        stop 'SQM needs property "particulated refractory organic carbon" - Check_Options'

                    if (.not.FindProperty(PropertiesList, HeterotrophicC_)       )   &
                        stop 'SQM needs property "heterotrophic microorganism carbon" - Check_Options'

                    if (.not.FindProperty(PropertiesList, AnaerobicC_)           )   &
                        stop 'SQM needs property "anaerobic microorganism carbon" - Check_Options'
                    
                    if (.not.FindProperty(PropertiesList, AutotrophicC_)         )   &
                        stop 'SQM needs property "autotrophic microorganism carbon" - Check_Options'
                    
!                    if (.not.FindProperty(PropertiesList, Methane_)         )   &
!                        stop 'SQM needs property "methane" - Check_Options'

                    if (.not.FindProperty(PropertiesList, CarbonDioxide_)         )   &
                        stop 'SQM needs property "carbon dioxide" - Check_Options'
                    
                    if (Sol_Bacteria) then
                        if (.not.FindProperty(PropertiesList, SolubilizingC_)    )   &
                            STOP 'SQM needs property "solubilizing microorganism carbon" - Check_Options'
                    end if


                end if cd11

cd14 :          if (Phosphorus) then
                   
                    if (.not.FindProperty(PropertiesList, Inorganic_Phosphorus_)              )   &
                        stop 'SQM needs property "inorganic phosphorus" - Check_Options'

                    if (.not.FindProperty(PropertiesList, AdsorbedInorganicP_)              )   &
                        stop 'SQM needs property "particulated inorganic phosphorus" - Check_Options'

                    if (.not.FindProperty(PropertiesList, RefreactaryOrganicP_)  )   &
                        stop 'SQM needs property "particulated refractory organic phosphorus" - Check_Options'

                    if (.not.FindProperty(PropertiesList, POP_           )       )   &
                        stop 'SQM needs property "particulate organic phosphorus" - Check_Options'

                    if (.not.FindProperty(PropertiesList, HeterotrophicP_)       )   &
                        stop 'SQM needs property "heterotrophic microorganism phosphorus" - Check_Options'

                    if (.not.FindProperty(PropertiesList, AnaerobicP_)           )   &
                        stop 'SQM needs property "anaerobic microorganism phosphorus" - Check_Options'

                    if (.not.FindProperty(PropertiesList, AutotrophicP_)         )   &
                        STOP 'SQM needs property "autotrophic microorganism phosphorus" - Check_Options'

                    if (Sol_Bacteria) then
                        if (.not.FindProperty(PropertiesList, SolubilizingP_)    )   &
                            STOP 'SQM needs property "solubilizing microorganism phosphorus" - Check_Options'
                    end if

                end if cd14

                if (.not.FindProperty(PropertiesList, Oxygen_))                     &
                    stop 'SQM needs property oxygen - Check_Options'         
            
             
            case(CEQUALW2Model, BenthicCEQUALW2Model)
                
                !Get CEQUALW2 model time step
                call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR05' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'CEQUALW2 wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN03.'
                endif
                
                call GetCEQUALW2PropertyList (Me%ObjCEQUALW2, CEQUALW2List, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR06'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, CEQUALW2List(i))) then
                        write(*,*) 'Property ',GetPropertyName(CEQUALW2List(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR07'    
                    end if
                end do

                call UngetCEQUALW2 (Me%ObjCEQUALW2, CEQUALW2List, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR08'


            case(LifeModel)

                !Get Life model time step
                call GetDTLife(Me%ObjLife, DT = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR09' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'Life wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN03.'
                endif
                
                call GetLifePropertyList (Me%ObjLife, LifeList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR10'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, LifeList(i))) then
                        write(*,*) 'Property ',GetPropertyName(LifeList(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR11'    
                    end if
                end do

                call UngetLife (Me%ObjLife, LifeList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR12'
#ifdef _BFM_  
            case(BFMModel)

                !Get BFM model time step
                call GetDTBFM(Me%ObjBFM, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR09a' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'Life wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN03a.'
                endif
                
                call GetBFMPropertyList (Me%ObjBFM, BFMList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR10a'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, BFMList(i))) then
                        write(*,*) 'Property ',GetPropertyName(BFMList(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options - ModuleInterface- ERR11a'    
                    end if
                end do

                call UngetBFM (Me%ObjBFM, BFMList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR12a'
#endif
            case(BenthosModel)

                !Get Benthos model time step
                call GetDTBenthos(Me%ObjBenthos, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR13' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'Benthos wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN04.'
                endif
                
                call GetBenthosPropertyList (Me%ObjBenthos, BenthosList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR14'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, BenthosList(i))) then
                        write(*,*) 'Property ',GetPropertyName(BenthosList(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR15'    
                    end if
                end do

                call UngetBenthos (Me%ObjBenthos, BenthosList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR16'

            case(MacroAlgaeModel)

                !Get MacroAlgae model time step
                call GetDTMacroAlgae(Me%ObjMacroAlgae, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR17' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'MacroAlgae wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN05.'
                endif
                
                call GetMacroAlgaePropertyList (Me%ObjMacroAlgae, MacroAlgaeList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR18'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, MacroAlgaeList(i))) then
                        write(*,*) 'Property ',GetPropertyName(MacroAlgaeList(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR19'    
                    end if
                end do

                call UngetMacroAlgae (Me%ObjMacroAlgae, MacroAlgaeList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR20'


#ifdef _PHREEQC_

            case (PhreeqCModel)

                !Get PhreeqC model time step
                call GetPhreeqCDT(Me%ObjPhreeqC, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR21' 
        
                !Run period in seconds
                RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the SedimentQuality DT
                auxFactor = RunPeriod / DT

                ErrorAux  = auxFactor - int(auxFactor)
                
                if (ErrorAux /= 0) then
                    Dtlag = int(ErrorAux * DT)
                    write(*,*) 
                    write(*,*) 'DTSECONDS is not multiple of the run period.'
                    write(*,*) 'PhreeqC wont be computed in the last', Dtlag, ' seconds.'
                    write(*,*) 'Check_Options - ModuleInterface - WRN06.'
                endif
                
                if (.NOT. FindProperty(PropertiesList, Temperature_)) &
                    stop 'PhreeqC needs property "temperature" - Check_Options'

                if (.NOT. FindProperty(PropertiesList, pH_)) &
                    stop 'PhreeqC needs property "ph" - Check_Options'

                if (.NOT. FindProperty(PropertiesList, pE_)) &
                    stop 'PhreeqC needs property "pe" - Check_Options'
                
                !Store number of properties involved
                Me%Prop%ILB = 1
                Me%Prop%IUB = SIZE(PropertiesList) - 3                       

#endif
            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR50' 
        end select


        call null_time   (Me%ExternalVar%Now)

        !----------------------------------------------------------------------

    end subroutine Check_Options

 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    subroutine GetWQRatio(interfaceID, PropertyID, Ratio, STAT)
    
        !Arguments-----------------------------------------------------------------
        integer                                         :: InterfaceID
        integer                                         :: PropertyID
        real                                            :: Ratio
        integer, optional,  intent(OUT)                 :: STAT
        
        !External--------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nProperty

        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_
        
        call Ready(InterfaceID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            !Number indexed to each property 
            nProperty = PropertyIndexNumber(PropertyID  )

            call GetNCRatio(Me%ObjWaterQuality,                          &
                            nProperty,                                   &
                            Ratio, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetMyRatio - ModuleInterface - ERR01'
        
        end if
        
        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetWQRatio
    !--------------------------------------------------------------------------
   
#ifdef _PHREEQC_    
    !--------------------------------------------------------------------------
    subroutine GetPhreeqCID (InterfaceID,  PhreeqCID, STAT)
    
        !Arguments-----------------------------------------------------------------
        integer                                         :: InterfaceID
        integer,            intent(OUT)                 :: PhreeqCID
        integer, optional,  intent(OUT)                 :: STAT
        
        !External--------------------------------------------------------------
        integer                                         :: ready_

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_

        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_
        
        call Ready(InterfaceID, ready_)
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
            PhreeqCID = Me%ObjPhreeqC
            STAT_ = SUCCESS_                    
        end if
        
        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------
            
    end subroutine GetPhreeqCID 
    !--------------------------------------------------------------------------
#endif
   
    !--------------------------------------------------------------------------
    subroutine GetRateFlux3D(InterfaceID,                            & 
                             FirstProp,                              &
                             SecondProp,                             &
                             RateIndex,                              &
                             RateFlux3D,                             &
                             WaterPoints3D,                          &
                             STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: InterfaceID
        real,               dimension(:,:,:),   pointer :: RateFlux3D
        integer,            dimension(:,:,:),   pointer :: WaterPoints3D
        integer, optional,  intent(IN)                  :: FirstProp
        integer, optional,  intent(IN)                  :: SecondProp
        integer, optional,  intent(IN)                  :: RateIndex
        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL          

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nFirstProp, nSecondProp          
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: Index
        integer                                         :: i, j, k
        real,    dimension(:), pointer                  :: RateFlux

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)   
         
        ILB     = Me%Size3D%ILB     
        IUB     = Me%Size3D%IUB     
        JLB     = Me%Size3D%JLB    
        JUB     = Me%Size3D%JUB    
        KLB     = Me%Size3D%KLB   
        KUB     = Me%Size3D%KUB   

        nullify (RateFlux)

cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(FirstProp)) then

                !Number indexed to each property 
                nFirstProp  = PropertyIndexNumber(FirstProp )
                nSecondProp = PropertyIndexNumber(SecondProp)

                !WaterPoints3D
                Me%ExternalVar%WaterPoints3D => WaterPoints3D
                
                select case (Me%SinksSourcesModel)

                    case (WaterQualityModel)

                        call GetWQPropRateFlux( Me%ObjWaterQuality,                     &
                                                nFirstProp, nSecondProp,                &
                                                RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR01'

                    case (SedimentQualityModel)
                    
                        call GetPropRateFlux(Me%ObjSedimentQuality,                     &
                                             nFirstProp, nSecondProp,                   &
                                             RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR02' 

                    case(LifeModel)
#ifdef _BFM_  
                    case(BFMModel)
#endif  

                    case (MacroAlgaeModel)
                    
                        call GetMacroAlgaeRateFlux(Me%ObjMacroAlgae,                          &
                                                   nFirstProp, nSecondProp,                   &
                                                   RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR03' 



                end select

            elseif (present(RateIndex)) then


                call GetCEQUALW2RateFlux( Me%ObjCeQualW2,                               &
                                          RateIndex, RateFlux,                          &
                                          STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR04'


            endif

            
            !Number indexed to 3D cell in the vector 
            Index = 0 

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then

                    Index = Index + 1
                    RateFlux3D(i, j, k) = RateFlux (Index)

                end if
            end do
            end do
            end do


            select case (Me%SinksSourcesModel)

                case (WaterQualityModel)

                    call UnGetWQPropRateFlux(Me%ObjWaterQuality, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR05' 
                
                case (SedimentQualityModel)

                    call UngetPropRateFlux(Me%ObjSedimentQuality, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR06' 

                case(CEQUALW2Model, BenthicCEQUALW2Model)
               
                    call UnGetCEQUALW2RateFlux(Me%ObjCeQualW2, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR07'
                
                case (MacroAlgaeModel)

                    call UnGetMacroAlgaeRateFlux(Me%ObjMacroAlgae, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR08' 

            end select

            STAT_ = SUCCESS_
        else
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetRateFlux3D

    !--------------------------------------------------------------------------


    subroutine GetRateFlux2D(InterfaceID,                            & 
                             FirstProp,                              &
                             SecondProp,                             &
                             RateFlux2D,                             &
                             WaterPoints2D,                          &
                             STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: InterfaceID
        real,               dimension(:,:  ),   pointer :: RateFlux2D
        integer,            dimension(:,:  ),   pointer :: WaterPoints2D
        integer,            intent(IN)                  :: FirstProp
        integer,            intent(IN)                  :: SecondProp
        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        integer                                         :: ready_         

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nFirstProp, nSecondProp          
        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: Index
        integer                                         :: i, j, STAT_CALL
        real,    dimension(:), pointer                  :: RateFlux

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)   

cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            ILB     = Me%Size2D%ILB     
            IUB     = Me%Size2D%IUB     
            JLB     = Me%Size2D%JLB    
            JUB     = Me%Size2D%JUB    

            nullify (RateFlux  )
            nullify (RateFlux2D)

            !Number indexed to each property 
            nFirstProp  = PropertyIndexNumber(FirstProp )
            nSecondProp = PropertyIndexNumber(SecondProp)

            !WaterPoints2D
            Me%ExternalVar%WaterPoints2D => WaterPoints2D
                
            select case (Me%SinksSourcesModel)

                case (BenthosModel)

                    call GetBenthosRateFlux(Me%ObjBenthos,                              &
                                            nFirstProp, nSecondProp,                    &
                                            RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux2D - ModuleInterface - ERR01'
                
            end select

           !Number indexed to 2D cell in the vector 
           Index = 0 

           do j = JLB, JUB
           do i = ILB, IUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then

                    Index = Index + 1
                    RateFlux2D(i, j) = RateFlux (Index)

                end if
            end do
            end do

            select case (Me%SinksSourcesModel)

                case (BenthosModel)

                    call UnGetBenthosRateFlux(Me%ObjBenthos, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux2D - ModuleInterface - ERR02'
                
            end select

            STAT_ = SUCCESS_
        else

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetRateFlux2D

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Modify_Interface3D(InterfaceID, PropertyID, Concentration,               &
                                  WaterPoints3D, OpenPoints3D, DWZ,                     &
                                  ShearStress, SPMFlux,                                 &
                                  ShortWaveAverage, ShortWaveTop,                       &
                                  LightExtCoefField, WaterPercentage,                   &
                                  DissolvedToParticulate3D, SoilDryDensity, Salinity,   &
                                  pH, IonicStrength, PhosphorusAdsortionIndex,          &
#ifdef _PHREEQC_                                  
                                  SolutionVolume, ConversionSelector, SolidMass,        &
#endif                                   
                                  WindVelocity,  DTProp, STAT)
                                 
        !Arguments-------------------------------------------------------------
        integer                                         :: InterfaceID
        integer, intent(IN)                             :: PropertyID
        real,              dimension(:,:,:), pointer    :: Concentration
        real,    optional, dimension(:,:,:), pointer    :: ShortWaveAverage
        real,    optional, dimension(:,:,:), pointer    :: ShortWaveTop
        real,    optional, dimension(:,:,:), pointer    :: LightExtCoefField
        real   , optional, dimension(:,:,:), pointer    :: WaterPercentage
        real,    optional, dimension(:,:,:), pointer    :: DissolvedToParticulate3D
        real,    optional, dimension(:,:,:), pointer    :: SoilDryDensity
        real,    optional, dimension(:,:,:), pointer    :: Salinity
        real,    optional, dimension(:,:,:), pointer    :: pH
#ifdef _PHREEQC_                                         
        real,    optional, dimension(:,:,:), pointer    :: SolutionVolume
        real,    optional, dimension(:,:,:), pointer    :: SolidMass
        integer, optional                               :: ConversionSelector
!        real,    optional, dimension(:,:,:), pointer    :: pE
#endif         
        real,    optional, dimension(:,:,:), pointer    :: IonicStrength
        real,    optional, dimension(:,:,:), pointer    :: PhosphorusAdsortionIndex
        real,    optional, dimension(:,:,:), pointer    :: WindVelocity
        integer,           dimension(:,:,:), pointer    :: WaterPoints3D
        integer, optional, dimension(:,:,:), pointer    :: OpenPoints3D
        real,    optional, dimension(:,:,:), pointer    :: DWZ, ShearStress, SPMFlux
        real,    optional,  intent(IN)                  :: DTProp
        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        real                                            :: DT
        integer                                         :: STAT_CALL
        integer                                         :: ready_
        logical                                         :: ReadyToCompute

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nProperty
        integer                                         :: Index 
        integer                                         :: i, j, k
        integer                                         :: prop, JulDay
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB     
        integer                                         :: PropLB, PropUB, ArrayLB, ArrayUB 
        real                                            :: DTProp_
        logical                                         :: Increment

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB     = Me%Size3D%ILB
            IUB     = Me%Size3D%IUB
            JLB     = Me%Size3D%JLB
            JUB     = Me%Size3D%JUB
            KLB     = Me%Size3D%KLB
            KUB     = Me%Size3D%KUB

            PropLB  = Me%Prop%ILB
            PropUB  = Me%Prop%IUB

            ArrayLB = Me%Array%ILB
            ArrayUB = Me%Array%IUB

            if(present(OpenPoints3D ))Me%ExternalVar%OpenPoints3D     => OpenPoints3D
            if(present(DWZ          ))Me%ExternalVar%DWZ              => DWZ
            if(present(ShearStress  ))Me%ExternalVar%ShearStress      => ShearStress
            if(present(SPMFlux      ))Me%ExternalVar%SPMFlux          => SPMFlux
            
            Me%ExternalVar%WaterPoints3D => WaterPoints3D
            
            ReadyToCompute = .false.

            Increment = OFF

cd7 :       if (present(DTProp))then
                DTProp_   = DTProp
                Increment = ON  !If Increment is ON  then the concentrations is will be changed
            else
                DTProp_   = FillValueReal
                Increment = OFF !If Increment is OFF then the sinks and sources module will be computed
            end if cd7


cd5 :       if (.not. Increment) then 

                call FillMassTempSalinity(PropertyID, Concentration, ReadyToCompute)

cd4 :           if (ReadyToCompute) then

                    call UnfoldMatrix(Me%ExternalVar%OpenPoints3D, Me%OpenPoints, Vertical1D = Me%ExternalVar%Vertical1D)

                    !Stores the concentration before changing them
                    Me%ConcentrationIncrement = Me%Mass

                    select case (Me%SinksSourcesModel)

                        case(WaterQualityModel)

                            call UnfoldMatrix(ShortWaveTop,Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField, Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ,Me%Thickness)


                            call WaterQuality(Me%ObjWaterQuality,                   &
                                              Me%Salinity,                          &
                                              Me%Temperature,                       &
                                              Me%ShortWaveTop,                      &
                                              Me%LightExtCoefField,                 &
                                              Me%Thickness,                         &
                                              Me%Mass,                              &
                                              Me%Array%ILB,                         &
                                              Me%Array%IUB,                         &
                                              Me%OpenPoints,                        &
                                              FishFood = Me%FishFood,               &
                                              STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR05'
                        
                        case(SedimentQualityModel)

                            call UnfoldMatrix(WaterPercentage, Me%WaterPercentage)
                            call UnfoldMatrix(DissolvedToParticulate3D, Me%DissolvedToParticulate)

                            call UnfoldMatrix(SoilDryDensity, Me%SoilDryDensity)
                            call UnfoldMatrix(Salinity, Me%Salinity)
                            call UnfoldMatrix(pH, Me%pH)
                            call UnfoldMatrix(IonicStrength, Me%IonicStrength)
                            call UnfoldMatrix(PhosphorusAdsortionIndex, Me%PhosphorusAdsortionIndex)
                            call UnfoldMatrix(WindVelocity, Me%WindVelocity)


                            call SedimentQuality(Me%ObjSedimentQuality,            &
                                                 Me%Temperature,                   &
                                                 Me%Mass,                          &
                                                 Me%WaterPercentage,               &
                                                 Me%DissolvedToParticulate,        &
                                                 Me%Array%ILB,                     &
                                                 Me%Array%IUB,                     &
                                                 Me%SoilDryDensity,                &
                                                 Me%Salinity,                      &
                                                 Me%pH,                            &
                                                 Me%IonicStrength,                 &
                                                 Me%PhosphorusAdsortionIndex,      &
                                                 Me%WindVelocity,                  &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR06'

                        case(LifeModel)

                            call UnfoldMatrix(ShortWaveAverage,  Me%ShortWaveAverage)
                            call UnfoldMatrix(ShortWaveTop,      Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField, Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ,Me%Thickness)
                            
                            call GetComputeCurrentTime(Me%ObjTime,                  &
                                                       Me%ExternalVar%Now,          &
                                                       STAT = STAT_CALL)                    
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERRXX' 


                            call JulianDay(Me%ExternalVar%Now, JulDay)

                            call ModifyLife(Me%ObjLife,                     &
                                      Me%Salinity,                          &
                                      Me%Temperature,                       &
                                      Me%ShortWaveAverage,                  &
                                      Me%ShortWaveTop,                      &
                                      Me%LightExtCoefField,                 &
                                      Me%Thickness,                         &
                                      Me%Mass,                              &
                                      Me%OpenPoints,                        &
                                      Me%Array,                             &
                                      JulDay,                               &
                                    STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR07'
#ifdef _BFM_  
                        case(BFMModel)
                           
                            call UnfoldMatrix(ShortWaveTop,Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField, Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ,Me%Thickness)


                            call ModifyBFM(Me%ObjBFM,                       &
                                           Me%Salinity,                     &
                                           Me%Temperature,                  &
                                           !Me%ShortWaveTop,                &
                                           !Me%LightExtCoefField,            &
                                           !Me%Thickness,                    &
                                           Me%OpenPoints,                   &
                                           Me%Mass,                         &
                                           STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR07'
#endif  

                        case(CEQUALW2Model)

                            call UnfoldMatrix(ShortWaveTop,      Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField, Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ,Me%Thickness)
                            

                            call CEQUALW2(Me%ObjCEQUALW2,                       &
                                          Me%Salinity,                          &
                                          Me%Temperature,                       &
                                          Me%Alkalinity,                        &
                                          Me%ShortWaveTop,                      &
                                          Me%LightExtCoefField,                 &
                                          Me%Thickness,                         &
                                          Me%Mass,                              &
                                          Me%OpenPoints,                        &
                                          STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR08'
                        
                         case(MacroAlgaeModel)

                            call UnfoldMatrix(ShortWaveTop,               Me%ShortWaveTop  )
                            call UnfoldMatrix(LightExtCoefField,          Me%LightExtCoefField   )
                            call UnfoldMatrix(Me%ExternalVar%DWZ,         Me%Thickness  )
                            call UnfoldMatrix(Me%ExternalVar%ShearStress, Me%ShearStress)
                            call UnfoldMatrix(Me%ExternalVar%SPMFlux,     Me%SPMFlux    ) 

                            call ModifyMacroAlgae(ObjMacroAlgaeID       = Me%ObjMacroAlgae,       &
                                                  Temperature           = Me%Temperature,         &
                                                  Salinity              = Me%Salinity,            &
                                                  OpenPoints            = Me%OpenPoints,          &
                                                  ShearStress           = Me%ShearStress,         &
                                                  SPMDepositionFlux     = Me%SPMFlux,             &
                                                  SWRadiation           = Me%ShortWaveTop,        &
                                                  SWLightExctintionCoef = Me%LightExtCoefField,   &
                                                  Thickness             = Me%Thickness,           &
                                                  Mass                  = Me%Mass,                &
                                                  STAT                  = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR08'

#ifdef _PHREEQC_
                        case(PhreeqcModel)
                         
                            call UnfoldMatrix(SolutionVolume, Me%SolutionVolume)
                            
                            if (present(SolidMass)) then
                            
                                call UnfoldMatrix(SolidMass, Me%SolidMass)
                                call ModifyPhreeqC(PhreeqCID = Me%ObjPhreeqC,               &
                                                   PropertiesValues = Me%Mass,              & 
                                                   SolutionVolume = Me%SolutionVolume,      &
                                                   SolutionTemperature = Me%Temperature,    &
                                                   SolutionpH = Me%pH,                      &
                                                   SolutionpE = Me%pE,                      &
                                                   SolidMass = Me%SolidMass,                &
                                                   CellsArrayLB = Me%Array%ILB,             &
                                                   CellsArrayUB = Me%Array%IUB,             &
                                                   OpenPoints = Me%OpenPoints,              & 
                                                   ConversionSelector = ConversionSelector, &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR14'

                            else

                                call ModifyPhreeqC(PhreeqCID = Me%ObjPhreeqC,               &
                                                   PropertiesValues = Me%Mass,              &   
                                                   SolutionVolume = Me%SolutionVolume,      &
                                                   SolutionTemperature = Me%Temperature,    &
                                                   SolutionpH = Me%pH,                      &
                                                   SolutionpE = Me%pE,                      &
                                                   CellsArrayLB = Me%Array%ILB,             &
                                                   CellsArrayUB = Me%Array%IUB,             &
                                                   OpenPoints = Me%OpenPoints,              & 
                                                   ConversionSelector = ConversionSelector, &
                                                   STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR14'
                            endif

#endif
                    end select

do7 :               do prop  = PropLB,  PropUB
do6 :               do index = ArrayLB, ArrayUB
                        Me%ConcentrationIncrement(prop, index) = Me%Mass(prop, index) - &
                                                                 Me%ConcentrationIncrement(prop, index)
                    end do do6
                    end do do7
                
                end if cd4

            elseif (Increment) then cd5

                Index = 0

                nProperty = PropertyIndexNumber(PropertyID)

                if (nProperty < 0) &
                    stop 'Modify_Interface3D - ModuleInterface - ERR80'

                select case (Me%SinksSourcesModel)

                    case(WaterQualityModel)
                        
                        call GetDTWQM(Me%ObjWaterQuality, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR09'

do1 :                   do k = KLB, KUB
do2 :                   do j = JLB, JUB
do3 :                   do i = ILB, IUB
cd2 :                       if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then
                                    Concentration(i, j, k) = Concentration( i, j, k)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd2
                        end do do3
                        end do do2
                        end do do1

                    case(SedimentQualityModel)
                        
                        call GetDTSedimentQuality(Me%ObjSedimentQuality, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR10'                

do8 :                   do k = KLB, KUB
do9 :                   do j = JLB, JUB
do10:                   do i = ILB, IUB
cd8 :                       if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
                                Index = Index + 1
                                Concentration(i, j, k) = Concentration( i, j, k)                     + &
                                                         Me%ConcentrationIncrement(nProperty, Index) * &
                                                         DTProp / DT 
                            end if cd8
                        end do do10
                        end do do9
                        end do do8



                    case(CEQUALW2Model)

                        call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR11'

do11 :                   do k = KLB, KUB
do12 :                   do j = JLB, JUB
do13 :                   do i = ILB, IUB
cd9 :                       if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then
                                    Concentration(i, j, k) = Concentration( i, j, k)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd9
                        end do do13
                        end do do12
                        end do do11


                   
                    case(LifeModel)

                        call GetDTLife(Me%ObjLife, DT = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR12'

do15 :                   do k = KLB, KUB
do16 :                   do j = JLB, JUB
do17 :                   do i = ILB, IUB
cd14 :                       if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then
                                    Concentration(i, j, k) = Concentration( i, j, k)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd14
                        end do do17
                        end do do16
                        end do do15
#ifdef _BFM_  
                    case(BFMModel)

                        call GetDTBFM(Me%ObjBFM, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR13a'

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
                                Index = Index + 1

                                !Concentrations are in all waterpoints
                                Concentration(i, j, k) = Concentration( i, j, k)      + &
                                Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 

                            end if 
                        end do 
                        end do 
                        end do 
#endif

                    case(MacroAlgaeModel)


                        call GetDTMacroAlgae(Me%ObjMacroAlgae, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR13'

                        do k = KLB, KUB
                        do j = JLB, JUB
                        do i = ILB, IUB
                            if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
                                Index = Index + 1

                                !Concentrations are in all waterpoints
                                Concentration(i, j, k) = Concentration( i, j, k)      + &
                                Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 

                            end if 
                        end do 
                        end do 
                        end do 
                        
#ifdef _PHREEQC_

                    case(PhreeqcModel)
                         
                        call GetPhreeqCDT(Me%ObjPhreeqC, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR15'                

do18:                   do k = KLB, KUB
do19:                   do j = JLB, JUB
do20:                   do i = ILB, IUB

cd15:                       if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then

                                Index = Index + 1
                                Concentration(i, j, k) = Concentration( i, j, k) +                     &
                                                         Me%ConcentrationIncrement(nProperty, Index) * &
                                                         DTProp / DT 
                                                         
                            end if cd15
                            
                        end do do20
                        end do do19
                        end do do18

#endif                   
                end select

            end if cd5

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_
            
        !----------------------------------------------------------------------

    end subroutine Modify_Interface3D

    !--------------------------------------------------------------------------

    subroutine Modify_Interface2D(InterfaceID, PropertyID, Concentration,               &
                                  WaterPoints2D, OpenPoints2D, Oxygen2D,                &
                                  DTProp, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: InterfaceID
        integer,            intent(IN)                  :: PropertyID
        real, optional,     intent(IN)                  :: DTProp
        real,              dimension(:,:  ), pointer    :: Concentration         
        integer,           dimension(:,:  ), pointer    :: WaterPoints2D
        integer, optional, dimension(:,:  ), pointer    :: OpenPoints2D
        real   , optional, dimension(:,:  ), pointer    :: Oxygen2D

        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL
        logical                                         :: ReadyToCompute

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nProperty
        integer                                         :: Index 
        integer                                         :: prop
        integer                                         :: ILB, IUB, JLB, JUB, i, j  
        integer                                         :: PropLB, PropUB, ArrayLB, ArrayUB 
        real                                            :: DTProp_, DT
        logical                                         :: Increment

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB     = Me%Size2D%ILB
            IUB     = Me%Size2D%IUB
            JLB     = Me%Size2D%JLB
            JUB     = Me%Size2D%JUB

            PropLB  = Me%Prop%ILB
            PropUB  = Me%Prop%IUB

            ArrayLB = Me%Array%ILB
            ArrayUB = Me%Array%IUB

            Me%ExternalVar%WaterPoints2D        => WaterPoints2D
            
            if(present(OpenPoints2D))then
                Me%ExternalVar%OpenPoints2D     => OpenPoints2D
            end if

            if(present(Oxygen2D))then
                call UnfoldMatrix(Oxygen2D, Me%Oxygen)
            end if                        
            
            Increment = OFF
            ReadyToCompute = .false.

cd7 :       if (present(DTProp))then
                DTProp_   = DTProp
                Increment = ON  !If Increment is ON  then the concentrations is will be changed
            else
                DTProp_   = FillValueReal
                Increment = OFF !If Increment is OFF then the sinks and sources module will be computed
            end if cd7


cd5 :       if (.not. Increment) then 

                call FillMassTempSalinity(PropertyID, Concentration, ReadyToCompute)

cd4 :           if (ReadyToCompute) then

                    call UnfoldMatrix(Me%ExternalVar%OpenPoints2D, Me%OpenPoints)

                    !Stores the concentration before changing them
                    Me%ConcentrationIncrement = Me%Mass

                    select case (Me%SinksSourcesModel)

                        case(BenthicCEQUALW2Model)
                            
                            if ( Me%UseSOD ) then
                            
                                call CEQUALW2Benthic(Me%ObjCEQUALW2,                      &
                                                    Me%Temperature,                       &
                                                    Me%Oxygen,                            &
                                                    Me%Mass,                              &                                                    
                                                    Me%OpenPoints,                        &
                                                    SODRate  = Me%SOD,                    &
                                                    STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERR01'
                            
                            else
                                    call CEQUALW2Benthic(Me%ObjCEQUALW2,                  &
                                                    Me%Temperature,                       &
                                                    Me%Oxygen,                            &
                                                    Me%Mass,                              &
                                                    Me%OpenPoints,                        &
                                                    STAT     = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERR01a'
                            endif
                        
                        case(BenthosModel)

                            call ModifyBenthos  (Me%ObjBenthos,                        &
                                                 Me%Temperature,                       &
                                                 Me%Oxygen,                            &
                                                 Me%OpenPoints,                        &
                                                 Me%Mass,                              &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERR02'

                    end select

do7 :               do prop  = PropLB,  PropUB
do6 :               do index = ArrayLB, ArrayUB
                        Me%ConcentrationIncrement(prop, index) = Me%Mass(prop, index) - &
                                                                 Me%ConcentrationIncrement(prop, index)
                    end do do6
                    end do do7
                
                end if cd4

            elseif (Increment) then cd5

                Index = 0

                nProperty = PropertyIndexNumber(PropertyID)

                if (nProperty < 0) &
                    stop 'Modify_Interface2D - ModuleInterface - ERR80'


                select case (Me%SinksSourcesModel)

                    case(BenthicCEQUALW2Model)

                        call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERR03'


                    case(BenthosModel)

                        call GetDTBenthos(Me%ObjBenthos, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERR04'

                end select

                do j = JLB, JUB
                do i = ILB, IUB
                    if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then
                        Index = Index + 1
                        
                        if (Me%ExternalVar%OpenPoints2D(i, j) == 1) then
                            Concentration(i, j) = Concentration(i, j)      + &
                            Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                        end if

                    end if
                end do
                end do


            end if cd5

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_
            
        !----------------------------------------------------------------------

    end subroutine Modify_Interface2D
    
    !--------------------------------------------------------------------------

    subroutine Modify_Interface1D(InterfaceID, PropertyID, Concentration,         &
                                  RiverPoints1D, OpenPoints1D, DWZ,               &
                                  ShortWaveAverage, ShortWaveTop, LightExtCoefField, &
                                  WaterPercentage, DissolvedToParticulate1D,      &
                                  SoilDryDensity, Salinity, pH, IonicStrength,    &
                                  PhosphorusAdsortionIndex, WindVelocity,         &
                                  Oxygen1D, WaterVolume, DTProp, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: InterfaceID
        integer, intent(IN)                             :: PropertyID
        real,              dimension(:), pointer        :: Concentration
        real,    optional, dimension(:), pointer        :: ShortWaveAverage
        real,    optional, dimension(:), pointer        :: ShortWaveTop
        real,    optional, dimension(:), pointer        :: LightExtCoefField
        real   , optional, dimension(:), pointer        :: WaterPercentage
        real,    optional, dimension(:), pointer        :: DissolvedToParticulate1D
        real,    optional, dimension(:,:,:), pointer    :: SoilDryDensity
        real,    optional, dimension(:,:,:), pointer    :: Salinity
        real,    optional, dimension(:,:,:), pointer    :: pH
        real,    optional, dimension(:,:,:), pointer    :: IonicStrength
        real,    optional, dimension(:,:,:), pointer    :: PhosphorusAdsortionIndex
        real,    optional, dimension(:,:,:), pointer    :: WindVelocity
        integer,           dimension(:), pointer        :: RiverPoints1D
        integer, optional, dimension(:), pointer        :: OpenPoints1D
        real,    optional, dimension(:), pointer        :: DWZ
        real,    optional, dimension(:), pointer        :: Oxygen1D
        real(8), optional, dimension(:), pointer        :: WaterVolume
        real,    optional,  intent(IN)                  :: DTProp
        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        real                                            :: DT
        integer                                         :: STAT_CALL
        integer                                         :: ready_
        logical                                         :: ReadyToCompute

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nProperty
        integer                                         :: Index 
        integer                                         :: i
        integer                                         :: prop, JulDay
        integer                                         :: ILB, IUB
        integer                                         :: PropLB, PropUB, ArrayLB, ArrayUB 
        real                                            :: DTProp_
        logical                                         :: Increment

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ILB     = Me%Size1D%ILB
            IUB     = Me%Size1D%IUB

            PropLB  = Me%Prop%ILB
            PropUB  = Me%Prop%IUB

            ArrayLB = Me%Array%ILB
            ArrayUB = Me%Array%IUB

            if(present(OpenPoints1D ))Me%ExternalVar%OpenPoints1D     => OpenPoints1D
            if(present(DWZ          ))Me%ExternalVar%DWZ1D            => DWZ
            
            if(present(Oxygen1D))then
                call UnfoldMatrix(Oxygen1D, Me%Oxygen)
            end if

            Me%ExternalVar%RiverPoints1D => RiverPoints1D

            Increment = OFF
            ReadyToCompute = .false.

cd7 :       if (present(DTProp))then
                DTProp_   = DTProp
                Increment = ON  !If Increment is ON  then the concentrations is will be changed
            else
                DTProp_   = FillValueReal
                Increment = OFF !If Increment is OFF then the sinks and sources module will be computed
            end if cd7


cd5 :       if (.not. Increment) then 

                call FillMassTempSalinity(PropertyID, Concentration, ReadyToCompute)

cd4 :           if (ReadyToCompute) then

                    call UnfoldMatrix(Me%ExternalVar%OpenPoints1D, Me%OpenPoints)

                    !Stores the concentration before changing them
                    Me%ConcentrationIncrement = Me%Mass

                    select case (Me%SinksSourcesModel)

                        case(WaterQualityModel)

                            call UnfoldMatrix(ShortWaveTop,         Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField,    Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ1D, Me%Thickness)


                            call WaterQuality(Me%ObjWaterQuality,                   &
                                              Me%Salinity,                          &
                                              Me%Temperature,                       &
                                              Me%ShortWaveTop,                      &
                                              Me%LightExtCoefField,                 &
                                              Me%Thickness,                         &
                                              Me%Mass,                              &
                                              Me%Array%ILB,                         &
                                              Me%Array%IUB,                         &
                                              Me%OpenPoints,                        &
                                              FishFood = Me%FishFood,               &
                                              STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR05'
                        
                        case(SedimentQualityModel)

                            call UnfoldMatrix(WaterPercentage, Me%WaterPercentage)
                            call UnfoldMatrix(DissolvedToParticulate1D, Me%DissolvedToParticulate)

                            call UnfoldMatrix(SoilDryDensity, Me%SoilDryDensity)
                            call UnfoldMatrix(Salinity, Me%Salinity)
                            call UnfoldMatrix(pH, Me%pH)
                            call UnfoldMatrix(IonicStrength, Me%IonicStrength)
                            call UnfoldMatrix(PhosphorusAdsortionIndex, Me%PhosphorusAdsortionIndex)
                            call UnfoldMatrix(WindVelocity, Me%WindVelocity)

                            call SedimentQuality(Me%ObjSedimentQuality,             &
                                                 Me%Temperature,                    &
                                                 Me%Mass,                           &
                                                 Me%WaterPercentage,                &
                                                 Me%DissolvedToParticulate,         &
                                                 Me%Array%ILB,                      &
                                                 Me%Array%IUB,                      &
                                                 Me%SoilDryDensity,                 &
                                                 Me%Salinity,                       &
                                                 Me%pH,                             &
                                                 Me%IonicStrength,                  &
                                                 Me%PhosphorusAdsortionIndex,       &
                                                 Me%WindVelocity,                   &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR06'

                        case(LifeModel)

                            call UnfoldMatrix(ShortWaveAverage,     Me%ShortWaveAverage)
                            call UnfoldMatrix(ShortWaveTop,         Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField,    Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ1D, Me%Thickness)
                            
                            call GetComputeCurrentTime(Me%ObjTime,                  &
                                                       Me%ExternalVar%Now,          &
                                                       STAT = STAT_CALL)                    
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERRXX' 


                            call JulianDay(Me%ExternalVar%Now, JulDay)

                            call ModifyLife(Me%ObjLife,                     &
                                      Me%Salinity,                          &
                                      Me%Temperature,                       &
                                      Me%ShortWaveAverage,                  &
                                      Me%ShortWaveTop,                      &
                                      Me%LightExtCoefField,                 &
                                      Me%Thickness,                         &
                                      Me%Mass,                              &
                                      Me%OpenPoints,                        &
                                      Me%Array,                             &
                                      JulDay,                               &
                                    STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR07'
#ifdef _BFM_  
                        case(BFMModel)
                        
                            call UnfoldMatrix(ShortWaveTop,         Me%ShortWaveRadiation)
                            call UnfoldMatrix(LightExtCoefField,    Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ1D, Me%Thickness)

                            call ModifyBFM(Me%ObjBFM,                       &
                                           Me%Salinity,                     &
                                           Me%Temperature,                  &
                                           !Me%ShortWaveTop,                &
                                           !Me%LightExtCoefField,            &
                                           !Me%Thickness,                    &
                                           Me%OpenPoints,                   &
                                           Me%Mass,                         &
                                           STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR07a'
#endif  
                        case(CEQUALW2Model)

                            call UnfoldMatrix(ShortWaveTop,          Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField,     Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ1D,  Me%Thickness)
                            

                            call CEQUALW2(Me%ObjCEQUALW2,                       &
                                          Me%Salinity,                          &
                                          Me%Temperature,                       &
                                          Me%Alkalinity,                        &
                                          Me%ShortWaveTop,                      &
                                          Me%LightExtCoefField,                 &
                                          Me%Thickness,                         &
                                          Me%Mass,                              &
                                          Me%OpenPoints,                        &
                                          STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR08'
                            
                        case(BenthosModel)

                            call ModifyBenthos  (Me%ObjBenthos,                        &
                                                 Me%Temperature,                       &
                                                 Me%Oxygen,                            &
                                                 Me%OpenPoints,                        &
                                                 Me%Mass,                              &
                                                 WaterVolume = WaterVolume,            &
                                                 STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR09'


                    end select

do7 :               do prop  = PropLB,  PropUB
do6 :               do index = ArrayLB, ArrayUB
                        Me%ConcentrationIncrement(prop, index) = Me%Mass(prop, index) - &
                                                                 Me%ConcentrationIncrement(prop, index)
                    end do do6
                    end do do7
                
                end if cd4

            elseif (Increment) then cd5

                Index = 0

                nProperty = PropertyIndexNumber(PropertyID)

                if (nProperty < 0) &
                    stop 'Modify_Interface1D - ModuleInterface - ERR80'


                select case (Me%SinksSourcesModel)

                    case(WaterQualityModel)
                        
                        call GetDTWQM(Me%ObjWaterQuality, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR09'

do3 :                   do i = ILB, IUB
cd2 :                       if (Me%ExternalVar%RiverPoints1D(i) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints1D(i) == 1) then
                                    Concentration(i) = Concentration( i)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd2
                        end do do3

                    case(SedimentQualityModel)
                        
                        call GetDTSedimentQuality(Me%ObjSedimentQuality, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR10'                

do10:                   do i = ILB, IUB
cd8 :                       if (Me%ExternalVar%RiverPoints1D(i) == 1) then
                                Index = Index + 1
                                Concentration(i) = Concentration( i)                     + &
                                                         Me%ConcentrationIncrement(nProperty, Index) * &
                                                         DTProp / DT 
                            end if cd8
                        end do do10



                    case(CEQUALW2Model)

                        call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR11'

do13 :                   do i = ILB, IUB
cd9 :                       if (Me%ExternalVar%RiverPoints1D(i) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints1D(i) == 1) then
                                    Concentration(i) = Concentration( i)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd9
                        end do do13


                   
                    case(LifeModel)

                        call GetDTLife(Me%ObjLife, DT = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR12'

do17 :                   do i = ILB, IUB
cd14 :                       if (Me%ExternalVar%RiverPoints1D(i) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints1D(i) == 1) then
                                    Concentration(i) = Concentration(i)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd14
                        end do do17
#ifdef _BFM_  
                    case(BFMModel)
                        call GetDTBFM(Me%ObjBFM, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR12a'

                        do i = ILB, IUB
                            
                            if (Me%ExternalVar%RiverPoints1D(i) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints1D(i) == 1) then
                                    Concentration(i) = Concentration(i)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if
                        end do
#endif  
                    case(BenthosModel)

                        call GetDTBenthos(Me%ObjBenthos, DTSecond = DT, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR13'

do18 :                  do i = ILB, IUB
cd10 :                      if (Me%ExternalVar%RiverPoints1D(i) == 1) then
                                Index = Index + 1
                                !Concentrations are only actualized in OpenPoints because of instability
                                !in waterpoints that are not openpoints
                                if (Me%ExternalVar%OpenPoints1D(i) == 1) then
                                    Concentration(i) = Concentration( i)      + &
                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
                                end if
                            end if cd10
                        end do do18

                end select

            end if cd5

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_
            
        !----------------------------------------------------------------------

    end subroutine Modify_Interface1D

    !--------------------------------------------------------------------------
    
        
    subroutine SetSOD (SOD, OpenPoints2D, WaterPoints2D )
    
        !Arguments-------------------------------------------------------------        
        real   , dimension(:,:  ), pointer      :: SOD        
        integer, dimension(:,:  ), pointer      :: OpenPoints2D
        integer, dimension(:,:  ), pointer      :: WaterPoints2D
        
        !External--------------------------------------------------------------        

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        
        Me%UseSOD = .true.
        
        Me%ExternalVar%OpenPoints2D        => OpenPoints2D
        Me%ExternalVar%WaterPoints2D       => WaterPoints2D       
        
        call UnfoldMatrix(SOD, Me%SOD)
    
        !nullify (Me%ExternalVar%OpenPoints2D)
        !nullify (Me%ExternalVar%WaterPoints2D)
    end subroutine SetSOD
    
    !--------------------------------------------------------------------------
    
    subroutine FillMassTempSalinity3D(PropertyID, Concentration, Ready) 

        !Arguments-------------------------------------------------------------
        integer,                intent(IN )     :: PropertyID
        real, dimension(:,:,:), pointer         :: Concentration
        logical,                intent(OUT)     :: Ready

        !External--------------------------------------------------------------
        integer                                 :: numZoo, numPhyto, numDiatoms 
        integer                                 :: numSiBio, numSiDiss 
        integer                                 :: numAmmonia, numNitrate, numNitrite
        integer                                 :: numDONRefractory
        integer                                 :: numDONNonRefractory
        integer                                 :: numPartOrganicNitrogen, numPartOrganicNitrogenRef
        integer                                 :: numPONitrogen1, numPONitrogen2, numPONitrogen3 
        integer                                 :: numPONitrogen4, numPONitrogen5
        integer                                 :: numPOPhosphorus1, numPOPhosphorus2, numPOPhosphorus3         
        integer                                 :: numPOPhosphorus4, numPOPhosphorus5
        integer                                 :: numOxygen, numBOD
        integer                                 :: numDOPRefractory, numDOPNonRefractory
        integer                                 :: numPartOrganicPhosphorus, numInorganicPhosphorus 
        integer                                 :: numBacteria, numCiliate
        integer                                 :: numLarvae
        integer                                 :: numHeterotrophicN, numHeterotrophicC
        integer                                 :: numAutotrophicP, numHeterotrophicP
        integer                                 :: numAutotrophicN, numAutotrophicC
        integer                                 :: numAnaerobicN, numAnaerobicC, numAnaerobicP        
        integer                                 :: numLabil_OM_N, numLabil_OM_C, numLabil_OM_P
        integer                                 :: numRefractOM_N, numRefractOM_C, numRefractOM_P
        integer                                 :: numNgas, numInorganicP_soluble, numInorganicP_fix
        integer                                 :: numSol_C, numSol_N, numSol_P, numCO2, numUrea
!        integer                                 :: numAmmoniaGas, numMethane
        integer                                 :: STAT_CALL
        logical                                 :: lZoo, lPhyto, lDiatoms 
        logical                                 :: lNitrogen, lPhosphorus, lSilica
        logical                                 :: lPompools 
        logical                                 :: lOxygen, lSalinity, lBOD
        logical                                 :: lBacteria, lCiliate
        logical                                 :: lLarvae
        logical                                 :: lCarbon, lSol_Bacteria

        !Local-----------------------------------------------------------------
        integer                                 :: i, PropLB, PropUB, IndexNumber 
        logical                                 :: TemperatureAdded = .false.
        logical                                 :: SalinityAdded = .false.
        logical                                 :: FishFoodAdded = .false.
        logical                                 :: AlkalinityAdded = .false.
#ifdef _PHREEQC_
        logical                                 :: pHAdded = .false.
        logical                                 :: pEAdded = .false.
!        logical                                 :: SoilDryDensityAdded = .false.
#endif        
        !----------------------------------------------------------------------

        PropLB  = Me%Prop%ILB
        PropUB  = Me%Prop%IUB

        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)

                call GetWQPropIndex(Me%ObjWaterQuality,                                         &
                                     Zoo                             = numZoo,                  &
                                     Larvae                          = numLarvae,               &
                                     Phyto                           = numPhyto,                &
                                     Diatoms                         = numDiatoms,              &
                                     Ammonia                         = numAmmonia,              &
                                     Nitrate                         = numNitrate,              &
                                     Nitrite                         = numNitrite,              &
                                     BiogenicSilica                  = numSiBio,                & 
                                     DissolvedSilica                 = numSiDiss,               & 
                                     DissOrganicNitrogenRefractory   = numDONRefractory,        &
                                     DONNonRefractory                = numDONNonRefractory,     &
                                     PartOrganicNitrogen             = numPartOrganicNitrogen,  &
                                     PartOrganicNitrogenRefractory   = numPartOrganicNitrogenRef,  &
                                     PONitrogen1                     = numPONitrogen1,          &
                                     PONitrogen2                     = numPONitrogen2,          &
                                     PONitrogen3                     = numPONitrogen3,          &
                                     PONitrogen4                     = numPONitrogen4,          &
                                     PONitrogen5                     = numPONitrogen5,          &                                
                                     Oxygen                          = numOxygen,               &
                                     BOD                             = numBOD,                  &
                                     DissOrganicPhosphorusRefractory = numDOPRefractory,        &
                                     DOPNonRefractory                = numDOPNonRefractory,     &
                                     PartOrganicPhosphorus           = numPartOrganicPhosphorus,&
                                     POPhosphorus1                   = numPOPhosphorus1,        &
                                     POPhosphorus2                   = numPOPhosphorus2,        &
                                     POPhosphorus3                   = numPOPhosphorus3,        &
                                     POPhosphorus4                   = numPOPhosphorus4,        &
                                     POPhosphorus5                   = numPOPhosphorus5,        &
                                     InorganicPhosphorus             = numInorganicPhosphorus,  &
                                     Bacteria                        = numBacteria,             &
                                     Ciliate                         = numCiliate,              &
                                     STAT                            = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity3D - ModuleInterface - ERR01' 

                call GetWQOptions(Me%ObjWaterQuality,                                           &
                                   Zoo            = lZoo,                                       &
                                   Larvae         = lLarvae,                                    &
                                   Phyto          = lPhyto,                                     &
                                   Diatoms        = lDiatoms,                                   & 
                                   Silica         = lSilica,                                    & 
                                   Nitrogen       = lNitrogen,                                  &
                                   Phosphorus     = lPhosphorus,                                &
                                   Oxygen         = lOxygen,                                    &
                                   Salinity       = lSalinity,                                  &
                                   BOD            = lBOD,                                       & 
                                   Bacteria       = lBacteria,                                  &
                                   Ciliate        = lCiliate,                                   &
                                   Pompools       = lPompools,                                  &
                                   STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity3D - ModuleInterface - ERR02'


cd16 :          if (lPhyto) then 
cd17 :          if (PropertyID == Phytoplankton_) then
                    call InputData(Concentration,numPhyto)
                    Me%AddedProperties(numPhyto)    = .TRUE.
                end if cd17
                end if cd16

cd18 :          if (lZoo) then
cd19 :          if (PropertyID == Zooplankton_) then
                    call InputData(Concentration,numZoo)
                    Me%AddedProperties(numZoo)      = .TRUE.
                end if cd19
                end if cd18

cd118 :         if (lLarvae) then
cd119 :         if (PropertyID == Larvae_) then
                    call InputData(Concentration,numLarvae)
                    Me%AddedProperties(numLarvae)   = .TRUE.
                end if cd119
                end if cd118

cd20 :          if (lPhosphorus) then

                    if (PropertyID == POP_) then
                        call InputData(Concentration, numPartOrganicPhosphorus)
                        Me%AddedProperties(numPartOrganicPhosphorus)    = .TRUE.
                    end if

                    if (PropertyID== DOPRefractory_) then
                        call InputData(Concentration, numDOPRefractory)
                        Me%AddedProperties(numDOPRefractory)            = .TRUE.
                    end if

                    if (PropertyID== DOPNon_Refractory_) then
                        call InputData(Concentration, numDOPNonRefractory)
                        Me%AddedProperties(numDOPNonRefractory)         = .TRUE.
                    end if

cd22 :              if (PropertyID== Inorganic_Phosphorus_) then
                        call InputData(Concentration, numInorganicPhosphorus)
                        Me%AddedProperties(numInorganicPhosphorus)      = .TRUE.
                    end if cd22
                    
                end if cd20


cd9 :           if (lNitrogen) then
cd10 :              if (PropertyID== PON_) then
                        call InputData(Concentration, numPartOrganicNitrogen)
                        Me%AddedProperties(numPartOrganicNitrogen)      = .TRUE.
                    end if cd10

cd11 :              if (PropertyID== PONRefractory_) then
                        call InputData(Concentration, numPartOrganicNitrogenRef)
                        Me%AddedProperties(numPartOrganicNitrogenRef)   = .TRUE.
                    end if cd11

cd12 :              if (PropertyID== DONRefractory_) then
                        call InputData(Concentration, numDONRefractory)
                        Me%AddedProperties(numDONRefractory)            = .TRUE.
                    end if cd12

cd13 :              if (PropertyID== DONNon_Refractory_) then
                        call InputData(Concentration,numDONNonRefractory)
                        Me%AddedProperties(numDONNonRefractory)         = .TRUE.
                    end if cd13

cd14 :              if (PropertyID== Ammonia_) then
                        call InputData(Concentration,numAmmonia)
                        Me%AddedProperties(numAmmonia)                  = .TRUE.
                    end if cd14

cd15 :              if (PropertyID== Nitrate_) then
                        call InputData(Concentration,numNitrate)
                        Me%AddedProperties(numNitrate)                  = .TRUE.
                    end if cd15                                                              

cd1126 :            if (PropertyID== Nitrite_) then
                        call InputData(Concentration,numNitrite)
                        Me%AddedProperties(numNitrite)                  = .TRUE.
                    end if cd1126
                end if cd9


cd1110 :       if (lPompools) then

cd1111 :         if (lNitrogen) then
                    
cd1116 :            if (PropertyID== PON1_) then
                        call InputData(Concentration, numPONitrogen1)
                        Me%AddedProperties(numPONitrogen1)              = .TRUE.
                    end if cd1116                    

cd1117 :             if (PropertyID== PON2_) then
                        call InputData(Concentration, numPONitrogen2)
                        Me%AddedProperties(numPONitrogen2)              = .TRUE.
                    end if cd1117 
                    
cd1118 :            if (PropertyID== PON3_) then
                        call InputData(Concentration, numPONitrogen3)
                        Me%AddedProperties(numPONitrogen3)              = .TRUE.
                    end if cd1118 
                    
cd1119 :             if (PropertyID== PON4_) then
                        call InputData(Concentration, numPONitrogen4)
                        Me%AddedProperties(numPONitrogen4)              = .TRUE.
                    end if cd1119                                                             

cd1120 :             if (PropertyID== PON5_) then
                        call InputData(Concentration, numPONitrogen5)
                        Me%AddedProperties(numPONitrogen5)              = .TRUE.
                    end if cd1120   
                 
                 end if cd1111   

cd1112 :         if (lPhosphorus) then
                 
cd1121 :             if (PropertyID== POP1_) then
                        call InputData(Concentration, numPOPhosphorus1)
                        Me%AddedProperties(numPOPhosphorus1)            = .TRUE.
                    end if cd1121                    

cd1122 :             if (PropertyID== POP2_) then
                        call InputData(Concentration, numPOPhosphorus2)
                        Me%AddedProperties(numPOPhosphorus2)            = .TRUE.
                    end if cd1122 
                    
cd1123 :             if (PropertyID== POP3_) then
                        call InputData(Concentration, numPOPhosphorus3)
                        Me%AddedProperties(numPOPhosphorus3)            = .TRUE.
                    end if cd1123 
                    
cd1124 :             if (PropertyID== POP4_) then
                        call InputData(Concentration, numPOPhosphorus4)
                        Me%AddedProperties(numPOPhosphorus4)            = .TRUE.
                    end if cd1124                                                             

cd1125 :             if (PropertyID== POP5_) then
                        call InputData(Concentration, numPOPhosphorus5)
                        Me%AddedProperties(numPOPhosphorus5)            = .TRUE.
                    end if cd1125                 
                 
                 end if cd1112
                                                           
                end if cd1110


cd7 :           if (lBOD) then
cd8 :           if (PropertyID== BOD_) then
                    call InputData(Concentration,numBOD)
                    Me%AddedProperties(numBOD)      = .TRUE.
                end if cd8
                end if cd7



cd23 :          if (lBacteria) then
cd24 :          if (PropertyID== Bacteria_) then
                    call InputData(Concentration,numBacteria)
                    Me%AddedProperties(numBacteria) = .TRUE.
                end if cd24
                end if cd23



cd25 :          if (lCiliate) then
cd26 :          if (PropertyID== Ciliate_) then
                    call InputData(Concentration,numCiliate)
                    Me%AddedProperties(numCiliate)  = .TRUE.
                end if cd26
                end if cd25


                !The oxygen is always compute. If it's not defined in the waterproperties 
                ! module then the saturation value is compute
cd6 :           if (PropertyID== Oxygen_) then
                    call InputData(Concentration,numOxygen)
                    Me%AddedProperties(numOxygen)   = .TRUE.
                end if cd6

                if (lSalinity) then
                
cd4 :           if (PropertyID== Salinity_) then
                        call UnfoldMatrix(Concentration, Me%Salinity)
                        SalinityAdded       =.TRUE.
                    end if cd4
                end if

cd3 :           if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if cd3

cd99:           if (PropertyID== FishFood_) then
                    call UnfoldMatrix(Concentration, Me%FishFood)
                    FishFoodAdded       =.TRUE.
                end if cd99

cd116 :         if (lDiatoms) then 
cd117 :         if (PropertyID == Diatoms_) then
                    call InputData(Concentration,numDiatoms)
                    Me%AddedProperties(numDiatoms)    = .TRUE.
                end if cd117
                end if cd116


cd150 :          if (lSilica) then

                    if (PropertyID == DSilica_) then
                        call InputData(Concentration, numSiDiss)
                        Me%AddedProperties(numSiDiss)    = .TRUE.
                    end if

                    if (PropertyID== BioSilica_) then
                        call InputData(Concentration, numSiBio)
                        Me%AddedProperties(numSiBio)            = .TRUE.
                    end if

                end if cd150




cd2 :           if ((SalinityAdded .or. (.not. lSalinity)) .AND. TemperatureAdded) then
                    Ready = .TRUE.

do1 :               do i = PropLB, PropUB
cd1 :               if (.NOT. Me%AddedProperties(i)) then
                        Ready = .FALSE.
                        exit do1
                    end if cd1
                    end do do1

                    if (Ready) Me%AddedProperties = .FALSE.
                end if cd2        
        
        case (SedimentQualityModel)

            call GetPropIndex(  Me%ObjSedimentQuality,                                      & 
                                HeterotrophicN                  = numHeterotrophicN,         &
                                HeterotrophicC                  = numHeterotrophicC,         &
                                AutotrophicN                    = numAutotrophicN,          &
                                AutotrophicC                    = numAutotrophicC,          &
                                AnaerobicN                      = numAnaerobicN,            &
                                AnaerobicC                      = numAnaerobicC,            &
                                Labil_OM_C                      = numLabil_OM_C,            &
                                Labil_OM_N                      = numLabil_OM_N,            &
                                RefractOM_C                     = numRefractOM_C,           &
                                RefractOM_N                     = numRefractOM_N,           &
                                Ammonia                         = numAmmonia,               &
                                Nitrate                         = numNitrate,               &
                                Ngas                            = numNgas,                  &
                                Oxygen                          = numOxygen,                &
                                HeterotrophicP                  = numHeterotrophicP,         &  
                                AutotrophicP                    = numAutotrophicP,          &
                                AnaerobicP                      = numAnaerobicP,            &
                                Labil_OM_P                      = numLabil_OM_P,            &
                                RefractOM_P                     = numRefractOM_P,           &
                                Inorganic_P_soluble             = numInorganicP_soluble,    &
                                Inorganic_P_fix                 = numInorganicP_fix,        &
                                SolC                            = numSol_C,                 &
                                SolN                            = numSol_N,                 &
                                SolP                            = numSol_P,                 &
                                CO2                             = numCO2,                   &
                                Urea                            = numUrea,                  &
!                                AmmoniaGas                      = numAmmoniaGas,            &
!                                Methane                         = numMethane,               &
                                STAT                            = STAT_CALL )          
             if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity3D - ModuleInterface - ERR03' 


             call GetSQOptions( Me%ObjSedimentQuality,                                      &
                                Nitrogen                        = lNitrogen,                &
                                Carbon                          = lCarbon,                  &
                                Phosphorus                      = lPhosphorus,              &
                                Sol_Bacteria                    = lSol_Bacteria,            &
                                Oxygen                          = lOxygen,                  &
                                STAT                            = STAT_CALL  )

             if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity - ModuleInterface - ERR04' 



cd27 :          if (lNitrogen) then

cd28 :              if (PropertyID== PON_                  ) then
                        call InputData(Concentration, numLabil_OM_N)
                        Me%AddedProperties(numLabil_OM_N)         = .TRUE.
                    end if cd28

cd29 :              if (PropertyID== RefreactaryOrganicN_  ) then
                        call InputData(Concentration, numRefractOM_N)
                        Me%AddedProperties(numRefractOM_N)        = .TRUE.
                    end if cd29

cd30 :              if (PropertyID== Ammonia_              ) then
                        call InputData(Concentration,numAmmonia)
                        Me%AddedProperties(numAmmonia)            = .TRUE.
                    end if cd30

cd31 :              if (PropertyID== Nitrate_              ) then
                        call InputData(Concentration,numNitrate)
                        Me%AddedProperties(numNitrate)            = .TRUE.
                    end if cd31

cd32 :              if (PropertyID== Ngas_                 ) then
                        call InputData(Concentration,numNgas)
                        Me%AddedProperties(numNgas)               = .TRUE.
                    end if cd32

cd33 :              if (PropertyID== HeterotrophicN_       ) then
                        call InputData(Concentration,numHeterotrophicN)
                        Me%AddedProperties(numHeterotrophicN)      = .TRUE.
                    end if cd33

cd34 :              if (PropertyID== AutotrophicN_         ) then
                        call InputData(Concentration,numAutotrophicN)
                        Me%AddedProperties(numAutotrophicN)    = .TRUE.
                    end if cd34

cd35 :              if (PropertyID== AnaerobicN_           ) then
                        call InputData(Concentration,numAnaerobicN)
                        Me%AddedProperties(numAnaerobicN)         = .TRUE.
                    end if cd35

cd350 :             if (PropertyID== Urea_                 ) then
                        call InputData(Concentration,numUrea)
                        Me%AddedProperties(numUrea)               = .TRUE.
                    end if cd350

!cd351 :             if (PropertyID== AmmoniaGas_           ) then
!                        call InputData(Concentration,numAmmoniaGas)
!                        Me%AddedProperties(numAmmoniaGas)         = .TRUE.
!                    end if cd351

cd352 :             if (lSol_Bacteria) then

cd353 :                 if (PropertyID== SolubilizingN_   ) then
                            call InputData(Concentration,numSol_N)
                            Me%AddedProperties(numSol_N)          = .TRUE.
                        end if cd353                    
                    
                    end if cd352
                    
                
                end if cd27


cd36 :          if (lCarbon) then

cd37 :              if (PropertyID== LabileOrganicC_       ) then
                        call InputData(Concentration, numLabil_OM_C)
                        Me%AddedProperties(numLabil_OM_C)         = .TRUE.
                    end if cd37

cd38 :              if (PropertyID== RefreactaryOrganicC_  ) then
                        call InputData(Concentration, numRefractOM_C)
                        Me%AddedProperties(numRefractOM_C)        = .TRUE.
                    end if cd38

cd39 :              if (PropertyID== HeterotrophicC_       ) then
                        call InputData(Concentration,numHeterotrophicC)
                        Me%AddedProperties(numHeterotrophicC)      = .TRUE.
                    end if cd39

cd40 :              if (PropertyID== AnaerobicC_           ) then
                        call InputData(Concentration,numAnaerobicC)
                        Me%AddedProperties(numAnaerobicC)         = .TRUE.
                    end if cd40

cd41 :              if (PropertyID== AutotrophicC_         ) then
                        call InputData(Concentration,numAutotrophicC)
                        Me%AddedProperties(numAutotrophicC)    = .TRUE.
                    end if cd41

cd410 :             if (PropertyID== CarbonDioxide_        ) then
                        call InputData(Concentration,numCO2)
                        Me%AddedProperties(numCO2)             = .TRUE.
                    end if cd410

!cd411 :             if (PropertyID== Methane_              ) then
!                        call InputData(Concentration,numMethane)
!                        Me%AddedProperties(numMethane)         = .TRUE.
!                    end if cd411

cd412 :             if (lSol_Bacteria) then

cd413 :                 if (PropertyID== SolubilizingC_   ) then
                            call InputData(Concentration,numSol_C)
                            Me%AddedProperties(numSol_C)          = .TRUE.
                        end if cd413                    
                    
                    end if cd412

                end if cd36


cd512 :         if (lPhosphorus) then

cd513:              if (PropertyID== POP_                  ) then
                        call InputData(Concentration, numLabil_OM_P)
                        Me%AddedProperties(numLabil_OM_P)         = .TRUE.
                    end if cd513

cd414 :             if (PropertyID== RefreactaryOrganicP_  ) then
                        call InputData(Concentration, numRefractOM_P)
                        Me%AddedProperties(numRefractOM_P)        = .TRUE.
                    end if cd414

cd415 :              if (PropertyID== Inorganic_Phosphorus_) then
                        call InputData(Concentration,numInorganicP_soluble)
                        Me%AddedProperties(numInorganicP_soluble) = .TRUE.
                    end if cd415

cd416 :              if (PropertyID== AdsorbedInorganicP_  ) then
                        call InputData(Concentration,numInorganicP_fix)
                        Me%AddedProperties(numInorganicP_fix)     = .TRUE.
                    end if cd416

cd417 :             if (PropertyID== AutotrophicP_         ) then
                        call InputData(Concentration,numAutotrophicP)
                        Me%AddedProperties(numAutotrophicP)       = .TRUE.
                    end if cd417                

cd418 :             if (PropertyID== HeterotrophicP_       ) then
                        call InputData(Concentration,numHeterotrophicP)
                        Me%AddedProperties(numHeterotrophicP)     = .TRUE.
                    end if cd418                

cd419 :             if (PropertyID== AnaerobicP_           ) then
                        call InputData(Concentration,numAnaerobicP)
                        Me%AddedProperties(numAnaerobicP)         = .TRUE.
                    end if cd419                 

cd420 :             if (lSol_Bacteria) then

cd421 :                 if (PropertyID== SolubilizingP_    ) then
                            call InputData(Concentration,numSol_P)
                            Me%AddedProperties(numSol_P)          = .TRUE.
                        end if cd421                    
                    
                    end if cd420
                
                end if cd512

                !The oxygen is always compute. If it's not defined in the waterproperties 
                ! module then the saturation value is compute
cd42 :          if (PropertyID== Oxygen_)then
                    call InputData(Concentration,numOxygen)
                    Me%AddedProperties(numOxygen)                 = .TRUE.
                end if cd42


cd43 :          if (PropertyID== Temperature_)then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded =.TRUE.
                end if cd43


cd44 :          if (TemperatureAdded) then
                    Ready = .TRUE.

do2  :              do i = PropLB, PropUB
cd45 :                  if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit do2
                    end if cd45
                    end do do2
                    if (Ready) Me%AddedProperties = .FALSE.
                
                end if cd44

            case(LifeModel)
                
                call GetLifePropIndex(Me%ObjLife,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR05' 
                end if
                     
                               
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

       
                if (SalinityAdded .AND. TemperatureAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if

#ifdef _BFM_  
            case(BFMModel)
                
                call GetBFMPropIndex(Me%ObjBFM,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR05a' 
                end if
                     
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

       
                if (SalinityAdded .AND. TemperatureAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if
#endif
 
            case(CEQUALW2Model) 
                
                call GetCEQUALW2PropIndex(Me%ObjCEQUALW2,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR05' 
                end if
                     
                               
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

                if (PropertyID== Alkalinity_) then
                    call UnfoldMatrix(Concentration, Me%Alkalinity)
                    AlkalinityAdded       =.TRUE.
                end if 

                if (SalinityAdded .AND. TemperatureAdded .AND. AlkalinityAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if

            case(MacroAlgaeModel)
                
                call GetMacroAlgaePropIndex(Me%ObjMacroAlgae,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR06' 
                end if
                     
                               
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

       
                if (SalinityAdded .AND. TemperatureAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if

#ifdef _PHREEQC_
            !ToDo: Need changes
            case (PhreeqCModel)

                select case (PropertyID)
                    case (Temperature_)
                        call UnfoldMatrix(Concentration, Me%Temperature)
                        TemperatureAdded = .TRUE.                    
                    case (pH_)
                        call UnfoldMatrix(Concentration, Me%pH)
                        pHAdded = .TRUE.                                        
                    case (pE_)
                        call UnfoldMatrix(Concentration, Me%pE)
                        pEAdded = .TRUE. 
!                    case (SoilDryDensity_)
!                        !SoilDryDensity will be used only in these cases:
!                        if ((PhreeqCSimOptions%Exchanger .EQ. 1) .OR. (PhreeqCSimOptions%SolidPhase)) then
!                            call UnfoldMatrix(Concentration, Me%SoilDryDensity)
!                        endif
!                        SoilDryDensityAdded = .TRUE.                                                            
                    case default
                        call GetPhreeqCPropIndex(Me%ObjPhreeqC, PropertyID, IndexNumber, STAT=STAT_CALL)
                        if (STAT_CALL == SUCCESS_) then
                        
                            call InputData(Concentration, IndexNumber)
                            Me%AddedProperties(IndexNumber) = .true.
                            
                        else if ((STAT_CALL .NE. SUCCESS_) .AND. (STAT_CALL .NE. NOT_FOUND_ERR_)) then
                            stop 'FillMassTempSalinity3D - ModuleInterface - ERR07' 
                        end if
                end select            
                
                if (TemperatureAdded .AND. pHAdded .AND. pEAdded) then ! .AND. SoilDryDensityAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB                       
                        if (.NOT. Me%AddedProperties(i)) then                        
                            Ready = .false.
                            exit                             
                        end if                         
                    end do 
                
                    if (Ready) Me%AddedProperties = .false.
                    
                end if
#endif
        end select

        !----------------------------------------------------------------------

    end Subroutine FillMassTempSalinity3D

    
    !--------------------------------------------------------------------------

    subroutine FillMassTempSalinity2D(PropertyID, Concentration, Ready) 

        !Arguments-------------------------------------------------------------
        integer,                intent(IN )     :: PropertyID
        real, dimension(:,:),   pointer         :: Concentration
        logical,                intent(OUT)     :: Ready

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: PropLB, PropUB, IndexNumber, STAT_CALL
        logical                                 :: TemperatureAdded = .false.
        integer                                 :: i
        !----------------------------------------------------------------------

        PropLB  = Me%Prop%ILB
        PropUB  = Me%Prop%IUB

        select case (Me%SinksSourcesModel)

            case(BenthicCeQualW2Model)

                call GetCEQUALW2PropIndex(Me%ObjCEQUALW2,PropertyID,IndexNumber,STAT=STAT_CALL)

            case(BenthosModel)

                call GetBenthosPropIndex(Me%ObjBenthos, PropertyID, IndexNumber, STAT = STAT_CALL)

        end select
        
        if (STAT_CALL==SUCCESS_) then
            call InputData(Concentration, IndexNumber)
            Me%AddedProperties(IndexNumber) =.TRUE.
        else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
            stop 'FillMassTempSalinity3D - ModuleInterface - ERR05' 
        end if

        if (PropertyID== Temperature_) then
            call UnfoldMatrix(Concentration, Me%Temperature)
            TemperatureAdded    =.TRUE.
        end if  

        if (TemperatureAdded) then
            Ready = .TRUE.

            do i = PropLB, PropUB
                if (.NOT. Me%AddedProperties(i)) then
                    Ready = .FALSE.
                    exit 
                end if 
            end do 
            if (Ready) Me%AddedProperties = .FALSE.

        end if


        !----------------------------------------------------------------------

    end Subroutine FillMassTempSalinity2D

    !--------------------------------------------------------------------------
    
    subroutine FillMassTempSalinity1D(PropertyID, Concentration, Ready) 

        !Arguments-------------------------------------------------------------
        integer,            intent(IN )         :: PropertyID
        real, dimension(:), pointer             :: Concentration
        logical,            intent(OUT)         :: Ready

        !External--------------------------------------------------------------
        integer                                 :: numZoo, numPhyto, numDiatoms 
        integer                                 :: numSiBio, numSiDiss 
        integer                                 :: numAmmonia, numNitrate, numNitrite
        integer                                 :: numDONRefractory
        integer                                 :: numDONNonRefractory
        integer                                 :: numPartOrganicNitrogen
        integer                                 :: numPartOrganicNitrogenRef
        integer                                 :: numOxygen, numBOD
        integer                                 :: numDOPRefractory, numDOPNonRefractory
        integer                                 :: numPartOrganicPhosphorus, numInorganicPhosphorus 
        integer                                 :: numPONitrogen1, numPONitrogen2, numPONitrogen3
        integer                                 :: numPONitrogen4, numPONitrogen5
        integer                                 :: numPOPhosphorus1, numPOPhosphorus2, numPOPhosphorus3
        integer                                 :: numPOPhosphorus4, numPOPhosphorus5
        integer                                 :: numBacteria, numCiliate
        integer                                 :: numLarvae
        integer                                 :: numHeterotrophicN, numHeterotrophicC
        integer                                 :: numAutotrophicN, numAutotrophicC
        integer                                 :: numAutotrophicP, numHeterotrophicP
        integer                                 :: numAnaerobicN, numAnaerobicC, numAnaerobicP         
        integer                                 :: numLabil_OM_N, numLabil_OM_C, numLabil_OM_P
        integer                                 :: numRefractOM_N, numRefractOM_C, numRefractOM_P
        integer                                 :: numNgas, numInorganicP_soluble, numInorganicP_fix
        integer                                 :: numSol_C, numSol_N, numSol_P, numCO2, numUrea
!        integer                                 :: numAmmoniaGas, numMethane 
        integer                                 :: STAT_CALL
        logical                                 :: lZoo, lPhyto, lDiatoms 
        logical                                 :: lNitrogen, lPhosphorus, lSilica
        logical                                 :: lPompools 
        logical                                 :: lOxygen, lSalinity, lBOD
        logical                                 :: lBacteria, lCiliate, lSol_Bacteria
        logical                                 :: lLarvae
        logical                                 :: lCarbon

        !Local-----------------------------------------------------------------
        integer                                 :: i, PropLB, PropUB, IndexNumber 
        logical                                 :: TemperatureAdded     = .false.
        logical                                 :: SalinityAdded        = .false.
        logical                                 :: FishFoodAdded        = .false.
        logical                                 :: AlkalinityAdded      = .false.
        !----------------------------------------------------------------------

        PropLB  = Me%Prop%ILB
        PropUB  = Me%Prop%IUB

        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)

                call GetWQPropIndex(Me%ObjWaterQuality,                                         &
                                     Zoo                             = numZoo,                  &
                                     Larvae                          = numLarvae,               &
                                     Phyto                           = numPhyto,                &
                                     Diatoms                         = numDiatoms,              & 
                                     Ammonia                         = numAmmonia,              &
                                     Nitrate                         = numNitrate,              &
                                     Nitrite                         = numNitrite,              &
                                     BiogenicSilica                  = numSiBio,                & 
                                     DissolvedSilica                 = numSiDiss,               & 
                                     DissOrganicNitrogenRefractory   = numDONRefractory,        &
                                     DONNonRefractory                = numDONNonRefractory,     &
                                     PartOrganicNitrogen             = numPartOrganicNitrogen,  &
                                     PartOrganicNitrogenRefractory   = numPartOrganicNitrogenRef,  &
                                     PONitrogen1                     = numPONitrogen1 ,         &
                                     PONitrogen2                     = numPONitrogen2 ,         &
                                     PONitrogen3                     = numPONitrogen3 ,         &
                                     PONitrogen4                     = numPONitrogen4 ,         &
                                     PONitrogen5                     = numPONitrogen5 ,         & 
                                     Oxygen                          = numOxygen,               &
                                     BOD                             = numBOD,                  &
                                     DissOrganicPhosphorusRefractory = numDOPRefractory,        &
                                     DOPNonRefractory                = numDOPNonRefractory,     &
                                     PartOrganicPhosphorus           = numPartOrganicPhosphorus,&
                                     POPhosphorus1                   = numPOPhosphorus1 ,       &
                                     POPhosphorus2                   = numPOPhosphorus2 ,       &
                                     POPhosphorus3                   = numPOPhosphorus3 ,       &
                                     POPhosphorus4                   = numPOPhosphorus4 ,       &
                                     POPhosphorus5                   = numPOPhosphorus5 ,       &
                                     InorganicPhosphorus             = numInorganicPhosphorus,  &
                                     Bacteria                        = numBacteria,             &
                                     Ciliate                         = numCiliate,              &
                                     STAT                            = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity1D - ModuleInterface - ERR01' 

                call GetWQOptions(Me%ObjWaterQuality,                                           &
                                   Zoo            = lZoo,                                       &
                                   Larvae         = lLarvae,                                    &
                                   Phyto          = lPhyto,                                     &
                                   Diatoms        = lDiatoms,                                   & 
                                   Silica         = lSilica,                                    & 
                                   Nitrogen       = lNitrogen,                                  &
                                   Phosphorus     = lPhosphorus,                                &
                                   Oxygen         = lOxygen,                                    &
                                   Salinity       = lSalinity,                                  &
                                   BOD            = lBOD,                                       & 
                                   Bacteria       = lBacteria,                                  &
                                   Ciliate        = lCiliate,                                   &
                                   Pompools       = lPompools,                                  &
                                   STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity1D - ModuleInterface - ERR02'


cd16 :          if (lPhyto) then 
cd17 :          if (PropertyID == Phytoplankton_) then
                    call InputData(Concentration,numPhyto)
                    Me%AddedProperties(numPhyto)    = .TRUE.
                end if cd17
                end if cd16

cd18 :          if (lZoo) then
cd19 :          if (PropertyID == Zooplankton_) then
                    call InputData(Concentration,numZoo)
                    Me%AddedProperties(numZoo)      = .TRUE.
                end if cd19
                end if cd18

cd118 :         if (lLarvae) then
cd119 :         if (PropertyID == Larvae_) then
                    call InputData(Concentration,numLarvae)
                    Me%AddedProperties(numLarvae)   = .TRUE.
                end if cd119
                end if cd118

cd20 :          if (lPhosphorus) then

                    if (PropertyID == POP_) then
                        call InputData(Concentration, numPartOrganicPhosphorus)
                        Me%AddedProperties(numPartOrganicPhosphorus)    = .TRUE.
                    end if

                    if (PropertyID== DOPRefractory_) then
                        call InputData(Concentration, numDOPRefractory)
                        Me%AddedProperties(numDOPRefractory)            = .TRUE.
                    end if

                    if (PropertyID== DOPNon_Refractory_) then
                        call InputData(Concentration, numDOPNonRefractory)
                        Me%AddedProperties(numDOPNonRefractory)         = .TRUE.
                    end if

cd22 :              if (PropertyID== Inorganic_Phosphorus_) then
                        call InputData(Concentration, numInorganicPhosphorus)
                        Me%AddedProperties(numInorganicPhosphorus)      = .TRUE.
                    end if cd22
                    
                end if cd20


cd9 :           if (lNitrogen) then
cd10 :              if (PropertyID== PON_) then
                        call InputData(Concentration, numPartOrganicNitrogen)
                        Me%AddedProperties(numPartOrganicNitrogen)      = .TRUE.
                    end if cd10

cd11 :              if (PropertyID== PONRefractory_) then
                        call InputData(Concentration, numPartOrganicNitrogenRef)
                        Me%AddedProperties(numPartOrganicNitrogenRef)   = .TRUE.
                    end if cd11

cd12 :              if (PropertyID== DONRefractory_) then
                        call InputData(Concentration, numDONRefractory)
                        Me%AddedProperties(numDONRefractory)            = .TRUE.
                    end if cd12

cd13 :              if (PropertyID== DONNon_Refractory_) then
                        call InputData(Concentration,numDONNonRefractory)
                        Me%AddedProperties(numDONNonRefractory)         = .TRUE.
                    end if cd13

cd14 :              if (PropertyID== Ammonia_) then
                        call InputData(Concentration,numAmmonia)
                        Me%AddedProperties(numAmmonia)                  = .TRUE.
                    end if cd14

cd15 :              if (PropertyID== Nitrate_) then
                        call InputData(Concentration,numNitrate)
                        Me%AddedProperties(numNitrate)                  = .TRUE.
                    end if cd15                    

cd1113 :            if (PropertyID== Nitrite_) then
                        call InputData(Concentration,numNitrite)
                        Me%AddedProperties(numNitrite)                  = .TRUE.
                    end if cd1113
                end if cd9


cd1110 :       if (lPompools) then

cd1111 :         if (lNitrogen) then
                    
cd1116 :            if (PropertyID== PON1_) then
                        call InputData(Concentration, numPONitrogen1)
                        Me%AddedProperties(numPONitrogen1)              = .TRUE.
                    end if cd1116                    

cd1117 :             if (PropertyID== PON2_) then
                        call InputData(Concentration, numPONitrogen2)
                        Me%AddedProperties(numPONitrogen2)              = .TRUE.
                    end if cd1117 
                    
cd1118 :            if (PropertyID== PON3_) then
                        call InputData(Concentration, numPONitrogen3)
                        Me%AddedProperties(numPONitrogen3)              = .TRUE.
                    end if cd1118 
                    
cd1119 :             if (PropertyID== PON4_) then
                        call InputData(Concentration, numPONitrogen4)
                        Me%AddedProperties(numPONitrogen4)              = .TRUE.
                    end if cd1119                                                             

cd1120 :             if (PropertyID== PON5_) then
                        call InputData(Concentration, numPONitrogen5)
                        Me%AddedProperties(numPONitrogen5)              = .TRUE.
                    end if cd1120   
                 
                 end if cd1111   

cd1112 :         if (lPhosphorus) then
                 
cd1121 :             if (PropertyID== POP1_) then
                        call InputData(Concentration, numPOPhosphorus1)
                        Me%AddedProperties(numPOPhosphorus1)            = .TRUE.
                    end if cd1121                    

cd1122 :             if (PropertyID== POP2_) then
                        call InputData(Concentration, numPOPhosphorus2)
                        Me%AddedProperties(numPOPhosphorus2)            = .TRUE.
                    end if cd1122 
                    
cd1123 :             if (PropertyID== POP3_) then
                        call InputData(Concentration, numPOPhosphorus3)
                        Me%AddedProperties(numPOPhosphorus3)            = .TRUE.
                    end if cd1123 
                    
cd1124 :             if (PropertyID== POP4_) then
                        call InputData(Concentration, numPOPhosphorus4)
                        Me%AddedProperties(numPOPhosphorus4)            = .TRUE.
                    end if cd1124                                                             

cd1125 :             if (PropertyID== POP5_) then
                        call InputData(Concentration, numPOPhosphorus5)
                        Me%AddedProperties(numPOPhosphorus5)            = .TRUE.
                    end if cd1125                 
                 
                 end if cd1112
                                                           
                end if cd1110




cd7 :           if (lBOD) then
cd8 :           if (PropertyID== BOD_) then
                    call InputData(Concentration,numBOD)
                    Me%AddedProperties(numBOD)      = .TRUE.
                end if cd8
                end if cd7



cd223 :          if (lBacteria) then
cd224 :          if (PropertyID== Bacteria_) then
                    call InputData(Concentration,numBacteria)
                    Me%AddedProperties(numBacteria) = .TRUE.
                end if cd224
                end if cd223



cd225 :          if (lCiliate) then
cd226 :          if (PropertyID== Ciliate_) then
                    call InputData(Concentration,numCiliate)
                    Me%AddedProperties(numCiliate)  = .TRUE.
                end if cd226
                end if cd225


                !The oxygen is always compute. If it's not defined in the waterproperties 
                ! module then the saturation value is compute
cd6 :           if (PropertyID== Oxygen_) then
                    call InputData(Concentration,numOxygen)
                    Me%AddedProperties(numOxygen)   = .TRUE.
                end if cd6

                if (lSalinity) then
                
cd4 :               if (PropertyID== Salinity_) then
                        call UnfoldMatrix(Concentration, Me%Salinity)
                        SalinityAdded       =.TRUE.
                    end if cd4
                end if

cd3 :           if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if cd3

cd99:           if (PropertyID== FishFood_) then
                    call UnfoldMatrix(Concentration, Me%FishFood)
                    FishFoodAdded       =.TRUE.
                end if cd99


cd216 :        if (lDiatoms) then 
cd217 :        if (PropertyID == Diatoms_) then
                    call InputData(Concentration,numDiatoms)
                    Me%AddedProperties(numDiatoms)    = .TRUE.
                end if cd217
                end if cd216


cd150 :          if (lSilica) then

                    if (PropertyID == DSilica_) then
                        call InputData(Concentration, numSiDiss)
                        Me%AddedProperties(numSiBio)    = .TRUE.
                    end if

                    if (PropertyID== BioSilica_) then
                        call InputData(Concentration, numSiBio)
                        Me%AddedProperties(numSiBio)            = .TRUE.
                    end if

                end if cd150


cd2 :           if ((SalinityAdded .or. (.not. lSalinity)) .AND. TemperatureAdded) then
                    Ready = .TRUE.

do1 :               do i = PropLB, PropUB
cd1 :               if (.NOT. Me%AddedProperties(i)) then
                        Ready = .FALSE.
                        exit do1
                    end if cd1
                    end do do1

                    if (Ready) Me%AddedProperties = .FALSE.
                end if cd2        
        
        case (SedimentQualityModel)

            call GetPropIndex(  Me%ObjSedimentQuality,                                      & 
                                HeterotrophicN                  = numHeterotrophicN,        &
                                HeterotrophicC                  = numHeterotrophicC,        &
                                AutotrophicN                    = numAutotrophicN,          &
                                AutotrophicC                    = numAutotrophicC,          &
                                AnaerobicN                      = numAnaerobicN,            &
                                AnaerobicC                      = numAnaerobicC,            &
                                Labil_OM_C                      = numLabil_OM_C,            &
                                Labil_OM_N                      = numLabil_OM_N,            &
                                RefractOM_C                     = numRefractOM_C,           &
                                RefractOM_N                     = numRefractOM_N,           &
                                Ammonia                         = numAmmonia,               &
                                Nitrate                         = numNitrate,               &
                                Ngas                            = numNgas,                  &
                                Oxygen                          = numOxygen,                &
                                HeterotrophicP                  = numHeterotrophicP,        &  
                                AutotrophicP                    = numAutotrophicP,          &
                                AnaerobicP                      = numAnaerobicP,            &
                                Labil_OM_P                      = numLabil_OM_P,            &
                                RefractOM_P                     = numRefractOM_P,           &
                                Inorganic_P_soluble             = numInorganicP_soluble,    &
                                Inorganic_P_fix                 = numInorganicP_fix,        &
                                SolC                            = numSol_C,                 &
                                SolN                            = numSol_N,                 &
                                SolP                            = numSol_P,                 &
                                CO2                             = numCO2,                   &
                                Urea                            = numUrea,                  &
!                                AmmoniaGas                      = numAmmoniaGas,             &
!                                Methane                         = numMethane,               &
                                STAT                            = STAT_CALL )          
             if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity1D - ModuleInterface - ERR03' 


             call GetSQOptions( Me%ObjSedimentQuality,                                      &
                                Nitrogen                        = lNitrogen,                &
                                Carbon                          = lCarbon,                  &
                                Phosphorus                      = lPhosphorus,              &
                                Sol_Bacteria                    = lSol_Bacteria,            &
                                Oxygen                          = lOxygen,                  &
                                STAT                            = STAT_CALL  )

             if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity1D - ModuleInterface - ERR04' 



cd27 :          if (lNitrogen) then

cd28 :              if (PropertyID== PON_                  ) then
                        call InputData(Concentration, numLabil_OM_N)
                        Me%AddedProperties(numLabil_OM_N)         = .TRUE.
                    end if cd28

cd29 :              if (PropertyID== RefreactaryOrganicN_  ) then
                        call InputData(Concentration, numRefractOM_N)
                        Me%AddedProperties(numRefractOM_N)        = .TRUE.
                    end if cd29

cd30 :              if (PropertyID== Ammonia_              ) then
                        call InputData(Concentration,numAmmonia)
                        Me%AddedProperties(numAmmonia)            = .TRUE.
                    end if cd30

cd31 :              if (PropertyID== Nitrate_              ) then
                        call InputData(Concentration,numNitrate)
                        Me%AddedProperties(numNitrate)            = .TRUE.
                    end if cd31

cd32 :              if (PropertyID== Ngas_                 ) then
                        call InputData(Concentration,numNgas)
                        Me%AddedProperties(numNgas)               = .TRUE.
                    end if cd32

cd33 :              if (PropertyID== HeterotrophicN_       ) then
                        call InputData(Concentration,numHeterotrophicN)
                        Me%AddedProperties(numHeterotrophicN)      = .TRUE.
                    end if cd33

cd34 :              if (PropertyID== AutotrophicN_         ) then
                        call InputData(Concentration,numAutotrophicN)
                        Me%AddedProperties(numAutotrophicN)    = .TRUE.
                    end if cd34

cd35 :              if (PropertyID== AnaerobicN_           ) then
                        call InputData(Concentration,numAnaerobicN)
                        Me%AddedProperties(numAnaerobicN)         = .TRUE.
                    end if cd35

cd350 :             if (PropertyID== Urea_                 ) then
                        call InputData(Concentration,numUrea)
                        Me%AddedProperties(numUrea)               = .TRUE.
                    end if cd350

!cd351 :             if (PropertyID== AmmoniaGas_           ) then
!                        call InputData(Concentration,numAmmoniaGas)
!                        Me%AddedProperties(numAmmoniaGas)         = .TRUE.
!                    end if cd351

cd352 :             if (lSol_Bacteria) then

cd353 :                 if (PropertyID== SolubilizingN_   ) then
                            call InputData(Concentration,numSol_N)
                            Me%AddedProperties(numSol_N)          = .TRUE.
                        end if cd353                    
                    
                    end if cd352
                
                end if cd27


cd36 :          if (lCarbon) then

cd37 :              if (PropertyID== LabileOrganicC_       ) then
                        call InputData(Concentration, numLabil_OM_C)
                        Me%AddedProperties(numLabil_OM_C)         = .TRUE.
                    end if cd37

cd38 :              if (PropertyID== RefreactaryOrganicC_  ) then
                        call InputData(Concentration, numRefractOM_C)
                        Me%AddedProperties(numRefractOM_C)        = .TRUE.
                    end if cd38

cd39 :              if (PropertyID== HeterotrophicC_       ) then
                        call InputData(Concentration,numHeterotrophicC)
                        Me%AddedProperties(numHeterotrophicC)      = .TRUE.
                    end if cd39

cd40 :              if (PropertyID== AnaerobicC_           ) then
                        call InputData(Concentration,numAnaerobicC)
                        Me%AddedProperties(numAnaerobicC)         = .TRUE.
                    end if cd40

cd41 :              if (PropertyID== AutotrophicC_         ) then
                        call InputData(Concentration,numAutotrophicC)
                        Me%AddedProperties(numAutotrophicC)    = .TRUE.
                    end if cd41

cd410 :             if (PropertyID== CarbonDioxide_        ) then
                        call InputData(Concentration,numCO2)
                        Me%AddedProperties(numCO2)             = .TRUE.
                    end if cd410

!cd411 :             if (PropertyID== Methane_              ) then
!                        call InputData(Concentration,numMethane)
!                        Me%AddedProperties(numMethane)         = .TRUE.
!                    end if cd411

cd412 :             if (lSol_Bacteria) then

cd413 :                 if (PropertyID== SolubilizingC_   ) then
                            call InputData(Concentration,numSol_C)
                            Me%AddedProperties(numSol_C)          = .TRUE.
                        end if cd413                    
                    
                    end if cd412

                end if cd36


cd512 :         if (lPhosphorus) then

cd513:              if (PropertyID== POP_                  ) then
                        call InputData(Concentration, numLabil_OM_P)
                        Me%AddedProperties(numLabil_OM_P)         = .TRUE.
                    end if cd513

cd414 :             if (PropertyID== RefreactaryOrganicP_  ) then
                        call InputData(Concentration, numRefractOM_P)
                        Me%AddedProperties(numRefractOM_P)        = .TRUE.
                    end if cd414

cd415 :              if (PropertyID== Inorganic_Phosphorus_) then
                        call InputData(Concentration,numInorganicP_soluble)
                        Me%AddedProperties(numInorganicP_soluble) = .TRUE.
                    end if cd415

cd416 :              if (PropertyID== AdsorbedInorganicP_  ) then
                        call InputData(Concentration,numInorganicP_fix)
                        Me%AddedProperties(numInorganicP_fix)     = .TRUE.
                    end if cd416

cd417 :             if (PropertyID== AutotrophicP_         ) then
                        call InputData(Concentration,numAutotrophicP)
                        Me%AddedProperties(numAutotrophicP)       = .TRUE.
                    end if cd417                

cd418 :             if (PropertyID== HeterotrophicP_       ) then
                        call InputData(Concentration,numHeterotrophicP)
                        Me%AddedProperties(numHeterotrophicP)     = .TRUE.
                    end if cd418                

cd419 :             if (PropertyID== AnaerobicP_           ) then
                        call InputData(Concentration,numAnaerobicP)
                        Me%AddedProperties(numAnaerobicP)         = .TRUE.
                    end if cd419                 

cd420 :             if (lSol_Bacteria) then

cd421 :                 if (PropertyID== SolubilizingP_    ) then
                            call InputData(Concentration,numSol_P)
                            Me%AddedProperties(numSol_P)          = .TRUE.
                        end if cd421                    
                    
                    end if cd420
                
                end if cd512

                !The oxygen is always compute. If it's not defined in the waterproperties 
                ! module then the saturation value is compute
cd42 :          if (PropertyID== Oxygen_)then
                    call InputData(Concentration,numOxygen)
                    Me%AddedProperties(numOxygen)                 = .TRUE.
                end if cd42


cd43 :          if (PropertyID== Temperature_)then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded =.TRUE.
                end if cd43


cd44 :          if (TemperatureAdded) then
                    Ready = .TRUE.

do2  :              do i = PropLB, PropUB
cd45 :                  if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit do2
                    end if cd45
                    end do do2
                    if (Ready) Me%AddedProperties = .FALSE.
                
                end if cd44

            case(LifeModel)

                call GetLifePropIndex(Me%ObjLife,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity1D - ModuleInterface - ERR05' 
                end if
                     
                               
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

       
                if (SalinityAdded .AND. TemperatureAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if
#ifdef _BFM_  
            case(BFMModel)

                call GetBFMPropIndex(Me%ObjBFM,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity1D - ModuleInterface - ERR05a' 
                end if
                     
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

       
                if (SalinityAdded .AND. TemperatureAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if
#endif 

            case(CEQUALW2Model) 
                
                call GetCEQUALW2PropIndex(Me%ObjCEQUALW2,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity1D - ModuleInterface - ERR05' 
                end if
                     
                               
                if (PropertyID== Salinity_) then
                    call UnfoldMatrix(Concentration, Me%Salinity)
                    SalinityAdded       =.TRUE.
                end if 

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if 

                if (PropertyID== Alkalinity_) then
                    call UnfoldMatrix(Concentration, Me%Alkalinity)
                    AlkalinityAdded       =.TRUE.
                end if 

                if (SalinityAdded .AND. TemperatureAdded .AND. AlkalinityAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                    if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                    end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.
                end if
            
            case(BenthosModel)

                call GetBenthosPropIndex(Me%ObjBenthos, PropertyID, IndexNumber, STAT = STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity1D - ModuleInterface - ERR05' 
                end if

                if (PropertyID== Temperature_) then
                    call UnfoldMatrix(Concentration, Me%Temperature)
                    TemperatureAdded    =.TRUE.
                end if  

                if (TemperatureAdded) then
                    Ready = .TRUE.

                    do i = PropLB, PropUB
                        if (.NOT. Me%AddedProperties(i)) then
                            Ready = .FALSE.
                            exit 
                        end if 
                    end do 
                    if (Ready) Me%AddedProperties = .FALSE.

                end if


        end select

        !----------------------------------------------------------------------

    end subroutine FillMassTempSalinity1D
    
    !----------------------------------------------------------------------

   
    subroutine InputData3D (Concentration, nProperty)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer         :: Concentration
        integer, intent (in)                    :: nProperty 

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j, k
        integer                                 :: ILB, IUB      
        integer                                 :: JLB, JUB     
        integer                                 :: KLB, KUB    

        !----------------------------------------------------------------------

        ILB = Me%Size3D%ILB     
        IUB = Me%Size3D%IUB     

        JLB = Me%Size3D%JLB    
        JUB = Me%Size3D%JUB    

        KLB = Me%Size3D%KLB   
        KUB = Me%Size3D%KUB   
        
        !Number indexed to 3D cell in the vector 
        Index = 0

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
                Index = Index + 1
                Me%Mass(nProperty,Index) = Concentration(i,j,k)
            endif    
        enddo
        enddo
        enddo

        if ((Index)> Me%Array%IUB) stop 'InputData3D - ModuleInterface - ERR01'

        !----------------------------------------------------------------------

    end subroutine InputData3D

    !--------------------------------------------------------------------------

    subroutine InputData2D (Concentration, nProperty)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer           :: Concentration
        integer, intent (in)                    :: nProperty 

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j
        integer                                 :: ILB, IUB      
        integer                                 :: JLB, JUB     

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB     
        IUB = Me%Size2D%IUB     

        JLB = Me%Size2D%JLB    
        JUB = Me%Size2D%JUB    

        
        !Number indexed to 3D cell in the vector 
        Index = 0

        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
                Index = Index + 1
                Me%Mass(nProperty,Index) = Concentration(i,j)
            endif    
        enddo
        enddo

        if ((Index)> Me%Array%IUB) stop 'InputData2D - ModuleInterface - ERR01'

        !----------------------------------------------------------------------

    end subroutine InputData2D

    !--------------------------------------------------------------------------

    subroutine InputData1D (Concentration, nProperty)

        !Arguments-------------------------------------------------------------
        real, dimension(:), pointer             :: Concentration
        integer, intent (in)                    :: nProperty 

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i
        integer                                 :: ILB, IUB      

        !----------------------------------------------------------------------

        ILB = Me%Size1D%ILB     
        IUB = Me%Size1D%IUB     

        !Number indexed to 1D cell in the vector 
        Index = 0

        do i = ILB, IUB
            if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then
                Index = Index + 1
                Me%Mass(nProperty,Index) = Concentration(i)
            endif    
        enddo

        if ((Index)> Me%Array%IUB) stop 'InputData1D - ModuleInterface - ERR01'

        !----------------------------------------------------------------------

    end subroutine InputData1D

    !--------------------------------------------------------------------------

    subroutine UnfoldMatrix3D_R (Matrix3D, Vector)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer     :: Matrix3D
        real, dimension(:    ), pointer     :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i, j, k
        integer                             :: ILB, IUB      
        integer                             :: JLB, JUB     
        integer                             :: KLB, KUB    

        !----------------------------------------------------------------------

        ILB = Me%Size3D%ILB     
        IUB = Me%Size3D%IUB     

        JLB = Me%Size3D%JLB    
        JUB = Me%Size3D%JUB    

        KLB = Me%Size3D%KLB   
        KUB = Me%Size3D%KUB   
        
        !Number indexed to 3D cell in the vector 
        Index = 0

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
                Index         = Index + 1
                Vector(Index) = Matrix3D(i,j,k)
            endif    
        enddo
        enddo
        enddo

        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix3D_R - ModuleInterface - ERR01'
            
        !----------------------------------------------------------------------

    end subroutine UnfoldMatrix3D_R

    !--------------------------------------------------------------------------

    subroutine UnfoldMatrix3D_I (Matrix3D, Vector, Vertical1D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:,:), pointer      :: Matrix3D
        integer, dimension(:    ), pointer      :: Vector
        logical, optional, intent(IN)           :: Vertical1D

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j, k
        integer                                 :: ILB, IUB      
        integer                                 :: JLB, JUB     
        integer                                 :: KLB, KUB
        logical                                 :: Vertical1D_Aux

        !----------------------------------------------------------------------

        ILB = Me%Size3D%ILB     
        IUB = Me%Size3D%IUB     

        JLB = Me%Size3D%JLB    
        JUB = Me%Size3D%JUB    

        KLB = Me%Size3D%KLB   
        KUB = Me%Size3D%KUB   

        if (present(Vertical1D)) then
            Vertical1D_Aux = Vertical1D
        else
            Vertical1D_Aux = .false.
        endif
        
        !Number indexed to 3D cell in the vector 
        Index = 0

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
                Index         = Index + 1
                Vector(Index) = Matrix3D(i,j,k)
                if (Vertical1D_Aux .and. .not.(i==2 .and. j==2)) Vector(Index) = 0
            endif    
        enddo
        enddo
        enddo

        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix3D_I - ModuleInterface - ERR01'
            

    end subroutine UnfoldMatrix3D_I

   !----------------------------------------------------------------------

    subroutine UnfoldMatrix2D_R (Matrix2D, Vector)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer       :: Matrix2D
        real, dimension(:  ), pointer       :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i, j
        integer                             :: ILB, IUB      
        integer                             :: JLB, JUB     
          
           
        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB     
        IUB = Me%Size2D%IUB     

        JLB = Me%Size2D%JLB    
        JUB = Me%Size2D%JUB    

          
        
        !Number indexed to 3D cell in the vector 
        Index = 0

        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then

                    Index         = Index + 1
                    
             if (Me%ExternalVar%OpenPoints2D(i,j)==1) then
                
                Vector(Index) = Matrix2D(i,j)
             
             endif
                    
                    
            endif    
        enddo
        enddo

        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix2D_R - ModuleInterface - ERR01'
            

    end subroutine UnfoldMatrix2D_R
        
    !----------------------------------------------------------------------

    subroutine UnfoldMatrix2D_I (Matrix2D, Vector)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:), pointer        :: Matrix2D
        integer, dimension(:  ), pointer        :: Vector

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j
        integer                                 :: ILB, IUB      
        integer                                 :: JLB, JUB     

        !----------------------------------------------------------------------

        ILB = Me%Size2D%ILB     
        IUB = Me%Size2D%IUB     

        JLB = Me%Size2D%JLB    
        JUB = Me%Size2D%JUB    
        
        !Number indexed to 3D cell in the vector 
        Index = 0

        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
             Index         = Index + 1
             if (Me%ExternalVar%OpenPoints2D(i,j)==1) then
                
                Vector(Index) = Matrix2D(i,j)
             
             endif
            endif    
        enddo
        enddo

        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix2D_I - ModuleInterface - ERR01'
            

    end subroutine UnfoldMatrix2D_I

   !----------------------------------------------------------------------

    subroutine UnfoldMatrix1D_R (Matrix1D, Vector)

        !Arguments-------------------------------------------------------------
        real, dimension(:  ), pointer       :: Matrix1D
        real, dimension(:  ), pointer       :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i
        integer                             :: ILB, IUB      
           
        !----------------------------------------------------------------------

        ILB = Me%Size1D%ILB     
        IUB = Me%Size1D%IUB     

        !Number indexed to 1D cell in the vector 
        Index = 0

        do i = ILB, IUB
            if (Me%ExternalVar%RiverPoints1D(i)==1) then

                    Index         = Index + 1
                    
             if (Me%ExternalVar%OpenPoints1D(i)==1) then
                
                Vector(Index) = Matrix1D(i)
             
             endif
                    
                    
            endif    
        enddo

        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix1D_R - ModuleInterface - ERR01'
            

    end subroutine UnfoldMatrix1D_R
        
    !----------------------------------------------------------------------

    subroutine UnfoldMatrix1D_I (Matrix1D, Vector)

        !Arguments-------------------------------------------------------------
        integer, dimension(:  ), pointer        :: Matrix1D
        integer, dimension(:  ), pointer        :: Vector

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i
        integer                                 :: ILB, IUB      
        !----------------------------------------------------------------------

        ILB = Me%Size1D%ILB     
        IUB = Me%Size1D%IUB     

        !Number indexed to 3D cell in the vector 
        Index = 0

        do i = ILB, IUB
            if (Me%ExternalVar%RiverPoints1D(i)==1) then
             Index         = Index + 1
             if (Me%ExternalVar%OpenPoints1D(i)==1) then
                
                Vector(Index) = Matrix1D(i)
             
             endif
            endif    
        enddo

        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix1D_I - ModuleInterface - ERR01'
            

    end subroutine UnfoldMatrix1D_I

    !--------------------------------------------------------------------------
        
    integer function PropertyIndexNumber(PropertyID)

        !Arguments-------------------------------------------------------------
        integer                         :: PropertyID

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL
        integer                         :: numZoo, numPhyto, numDiatoms 
        integer                         :: numLarvae
        integer                         :: numAmmonia, numNitrate, numNitrite
        integer                         :: numSiBio, numSiDiss 
        integer                         :: numDONRefractory, numDONNonRefractory
        integer                         :: numPartOrganicNitrogen
        integer                         :: numPONitrogen1, numPONitrogen2, numPONitrogen3
        integer                         :: numPONitrogen4, numPONitrogen5
        integer                         :: numPOPhosphorus1, numPOPhosphorus2, numPOPhosphorus3
        integer                         :: numPOPhosphorus4, numPOPhosphorus5      
        integer                         :: numPartOrganicNitrogenRef
        integer                         :: numOxygen, numBOD
        integer                         :: numDOPRefractory, numDOPNonRefractory
        integer                         :: numPartOrganicPhosphorus, numInorganicPhosphorus 
        integer                         :: numBacteria, numCiliate
        integer                         :: numHeterotrophicN, numHeterotrophicC
        integer                         :: numAutotrophicN, numAutotrophicC
        integer                         :: numAutotrophicP, numHeterotrophicP
        integer                         :: numAnaerobicN, numAnaerobicC, numAnaerobicP        
        integer                         :: numLabil_OM_N, numLabil_OM_C, numLabil_OM_P
        integer                         :: numRefractOM_N, numRefractOM_C, numRefractOM_P
        integer                         :: numNgas, numInorganicP_soluble, numInorganicP_fix
        integer                         :: numSol_C, numSol_N, numSol_P, numCO2, numUrea
!        integer                         :: numAmmoniaGas, numMethane 
        
        !Local-----------------------------------------------------------------
        integer                         :: nProperty
        logical                         :: CheckName

        !Begin----------------------------------------------------------------
        

        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)
        
                call GetWQPropIndex(Me%ObjWaterQuality,                                  &
                             Zoo                              = numZoo,                  &
                             Larvae                           = numLarvae,               &
                             Phyto                            = numPhyto,                &
                             Diatoms                          = numDiatoms,              & 
                             Ammonia                          = numAmmonia,              &
                             Nitrate                          = numNitrate,              &
                             Nitrite                          = numNitrite,              &
                             BiogenicSilica                   = numSiBio,                & 
                             DissolvedSilica                  = numSiDiss,               & 
                             DissOrganicNitrogenRefractory    = numDONRefractory,        &
                             DONNonRefractory                 = numDONNonRefractory,     &
                             PartOrganicNitrogen              = numPartOrganicNitrogen,  &
                             PartOrganicNitrogenRefractory    = numPartOrganicNitrogenRef,  &
                             PONitrogen1                      = numPONitrogen1 ,         &
                             PONitrogen2                      = numPONitrogen2 ,         &
                             PONitrogen3                      = numPONitrogen3 ,         &
                             PONitrogen4                      = numPONitrogen4 ,         &
                             PONitrogen5                      = numPONitrogen5 ,         & 
                             Oxygen                           = numOxygen,               &
                             BOD                              = numBOD,                  &
                             DissOrganicPhosphorusRefractory  = numDOPRefractory,        &
                             DOPNonRefractory                 = numDOPNonRefractory,     &
                             PartOrganicPhosphorus            = numPartOrganicPhosphorus,&
                             POPhosphorus1                    = numPOPhosphorus1 ,       &
                             POPhosphorus2                    = numPOPhosphorus2 ,       &
                             POPhosphorus3                    = numPOPhosphorus3 ,       &
                             POPhosphorus4                    = numPOPhosphorus4 ,       &
                             POPhosphorus5                    = numPOPhosphorus5 ,       &
                             InorganicPhosphorus              = numInorganicPhosphorus,  &
                             Bacteria                         = numBacteria,             &
                             Ciliate                          = numCiliate,              &
                             STAT                             = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'PropertyIndexNumber - ModuleInterface - ERR01'
                   
   
cd1 :           if      (PropertyID== Phytoplankton_       ) then
                    nProperty = numPhyto

                else if (PropertyID== Zooplankton_         ) then cd1
                    nProperty = numZoo
 
                else if (PropertyID== Larvae_              ) then cd1
                    nProperty = numLarvae

                else if (PropertyID== POP_                 ) then cd1
                    nProperty = numPartOrganicPhosphorus

                else if (PropertyID== DOPRefractory_       ) then cd1
                    nProperty = numDOPRefractory

                else if (PropertyID== DOPNon_Refractory_   ) then cd1
                    nProperty = numDOPNonRefractory

                else if (PropertyID== Inorganic_Phosphorus_) then cd1
                    nProperty = numInorganicPhosphorus

                else if (PropertyID== PON_                 ) then cd1
                    nProperty = numPartOrganicNitrogen

                else if (PropertyID== PONRefractory_       ) then cd1
                    nProperty = numPartOrganicNitrogenRef

                else if (PropertyID== DONRefractory_       ) then cd1
                    nProperty = numDONRefractory

                else if (PropertyID== DONNon_Refractory_   ) then cd1
                    nProperty = numDONNonRefractory 

                else if (PropertyID== Ammonia_             ) then cd1
                    nProperty = numAmmonia

                else if (PropertyID== Nitrate_             ) then cd1
                    nProperty = numNitrate

                else if (PropertyID== Nitrite_             ) then cd1
                    nProperty = numNitrite
                    
                else if (PropertyID== PON1_                ) then cd1
                    nProperty = numPONitrogen1
                
                else if (PropertyID== PON2_                ) then cd1
                    nProperty = numPONitrogen2        
                
                else if (PropertyID== PON3_                ) then cd1
                    nProperty = numPONitrogen3    
                
                else if (PropertyID== PON4_                ) then cd1
                    nProperty = numPONitrogen4
                
                else if (PropertyID== PON5_                ) then cd1
                    nProperty = numPONitrogen5
                    
                else if (PropertyID== POP1_                ) then cd1
                    nProperty = numPOPhosphorus1
                
                else if (PropertyID== POP2_                ) then cd1
                    nProperty = numPOPhosphorus2        
                
                else if (PropertyID== POP3_                ) then cd1
                    nProperty = numPOPhosphorus3    
                
                else if (PropertyID== POP4_                ) then cd1
                    nProperty = numPOPhosphorus4
                
                else if (PropertyID== POP5_                ) then cd1
                    nProperty = numPOPhosphorus5    
                    
                else if (PropertyID== Bacteria_            ) then cd1
                    nProperty = numBacteria

                else if (PropertyID== Ciliate_             ) then cd1
                    nProperty = numCiliate

                else if (PropertyID== BOD_                 ) then cd1
                    nProperty = numBOD

                else if (PropertyID== Oxygen_              ) then cd1
                    nProperty = numOxygen
                
                else if (PropertyID== Diatoms_            ) then cd1
                    nProperty = numDiatoms 

                else if (PropertyID== DSilica_            ) then cd1
                    nProperty = numSiDiss 

                else if (PropertyID== BioSilica_           ) then cd1
                    nProperty = numSiBio

                else if (PropertyID== GrossProd_           ) then cd1 
                    CheckName = CheckPropertyName('grossprod', number = nProperty)
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR03'

                else if (PropertyID== NutrientLim_         ) then cd1   
                    CheckName = CheckPropertyName('nutrientlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR04'

                else if (PropertyID== NLim_                ) then cd1   
                    CheckName = CheckPropertyName('nitrogenlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR05'

                else if (PropertyID== PLim_                ) then cd1   
                    CheckName = CheckPropertyName('phosphoruslim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR06'

                else if (PropertyID== LightLim_            ) then cd1   
                    CheckName = CheckPropertyName('lightlim', number = nProperty) 
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR07'

                else if (PropertyID== TemperatureLim_      ) then cd1 
                    CheckName = CheckPropertyName('temperaturelim', number = nProperty)    
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR08'

                else if (PropertyID== DiaGrossProd_           ) then cd1 
                    CheckName = CheckPropertyName('diagrossprod', number = nProperty)
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR09'

                else if (PropertyID== DiaNutrientLim_         ) then cd1   
                    CheckName = CheckPropertyName('dianutrientlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR10'

                else if (PropertyID== DiaNLim_         ) then cd1   
                    CheckName = CheckPropertyName('dianitrogenlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR11'

                else if (PropertyID== DiaPLim_         ) then cd1   
                    CheckName = CheckPropertyName('diaphosphoruslim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR12'

                else if (PropertyID== DiaSiLim_         ) then cd1   
                    CheckName = CheckPropertyName('diasilicalim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR13'

                else if (PropertyID== DiaLightLim_            ) then cd1   
                    CheckName = CheckPropertyName('dialightlim', number = nProperty) 
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR14'

                else if (PropertyID== DiaTemperatureLim_      ) then cd1 
                    CheckName = CheckPropertyName('diatemperaturelim', number = nProperty)    
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR15'



                else  cd1
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and WaterQuality.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR16'
                end if cd1


            case (SedimentQualityModel)
                
                call GetPropIndex(  Me%ObjSedimentQuality,                                  & 
                                    HeterotrophicN                  = numHeterotrophicN,    &
                                    HeterotrophicC                  = numHeterotrophicC,    &
                                    AutotrophicN                    = numAutotrophicN,      &
                                    AutotrophicC                    = numAutotrophicC,      &
                                    AnaerobicN                      = numAnaerobicN,        &
                                    AnaerobicC                      = numAnaerobicC,        &
                                    Labil_OM_C                      = numLabil_OM_C,        &
                                    Labil_OM_N                      = numLabil_OM_N,        &
                                    RefractOM_C                     = numRefractOM_C,       &
                                    RefractOM_N                     = numRefractOM_N,       &
                                    Ammonia                         = numAmmonia,           &
                                    Nitrate                         = numNitrate,           &
                                    Ngas                            = numNgas,              &
                                    Oxygen                          = numOxygen,            &
                                    HeterotrophicP                  = numHeterotrophicP,    &  
                                    AutotrophicP                    = numAutotrophicP,      &
                                    AnaerobicP                      = numAnaerobicP,        &
                                    Labil_OM_P                      = numLabil_OM_P,        &
                                    RefractOM_P                     = numRefractOM_P,       &
                                    Inorganic_P_soluble             = numInorganicP_soluble,&
                                    Inorganic_P_fix                 = numInorganicP_fix,    &
                                    SolC                            = numSol_C,             &
                                    SolN                            = numSol_N,             &
                                    SolP                            = numSol_P,             &
                                    CO2                             = numCO2,               &
                                    Urea                            = numUrea,              &
!                                    AmmoniaGas                      = numAmmoniaGas,        &
!                                    Methane                         = numMethane,           &
                                    STAT                            = STAT_CALL )          
                if (STAT_CALL .NE. SUCCESS_) stop 'PropertyIndexNumber - ModuleInterface - ERR17'


                if      (PropertyID== PON_                 ) then            
                    nProperty = numLabil_OM_N
            
                else if (PropertyID== RefreactaryOrganicN_ ) then
                    nProperty = numRefractOM_N

                else if (PropertyID== Ammonia_             ) then
                    nProperty = numAmmonia

                else if (PropertyID== Nitrate_             ) then
                    nProperty = numNitrate

                else if (PropertyID== Ngas_                ) then
                    nProperty = numNgas

                else if (PropertyID== HeterotrophicN_      ) then
                    nProperty = numHeterotrophicN

                else if (PropertyID== AutotrophicN_        ) then
                    nProperty = numAutotrophicN

                else if (PropertyID== AnaerobicN_          ) then
                    nProperty = numAnaerobicN

                else if (PropertyID== Urea_                ) then
                    nProperty = numUrea

!                else if (PropertyID== AmmoniaGas_          ) then
!                    nProperty = numAmmoniaGas

                else if (PropertyID== SolubilizingN_       ) then
                    nProperty = numSol_N

                else if (PropertyID== LabileOrganicC_      ) then
                    nProperty = numLabil_OM_C

                else if (PropertyID== RefreactaryOrganicC_ ) then
                    nProperty = numRefractOM_C

                else if (PropertyID== HeterotrophicC_      ) then
                    nProperty = numHeterotrophicC

                else if (PropertyID== AnaerobicC_          ) then
                    nProperty = numAnaerobicC

                else if (PropertyID== AutotrophicC_        ) then
                    nProperty = numAutotrophicC

                else if (PropertyID== CarbonDioxide_       ) then
                    nProperty = numCO2

!                else if (PropertyID== Methane_             ) then
!                    nProperty = numMethane

                else if (PropertyID== SolubilizingC_       ) then
                    nProperty = numSol_C

                else if (PropertyID== POP_                 ) then            
                    nProperty = numLabil_OM_P
            
                else if (PropertyID== RefreactaryOrganicP_ ) then
                    nProperty = numRefractOM_P

                else if (PropertyID== Inorganic_Phosphorus_) then
                    nProperty = numInorganicP_soluble

                else if (PropertyID== AdsorbedInorganicP_  ) then
                    nProperty = numInorganicP_fix

                else if (PropertyID== HeterotrophicP_      ) then
                    nProperty = numHeterotrophicP

                else if (PropertyID== AutotrophicP_        ) then
                    nProperty = numAutotrophicN

                else if (PropertyID== AnaerobicP_          ) then
                    nProperty = numAnaerobicP

                else if (PropertyID== SolubilizingP_       ) then
                    nProperty = numSol_P

                else if (PropertyID== Oxygen_              ) then
                    nProperty = numOxygen

                else
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and SedimentQuality.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR18'

                end if

            case(LifeModel)
                call GetLifePropIndex(Me%ObjLife,PropertyID,nProperty,STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
                    
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and Life.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR19'

                end if  
#ifdef _BFM_  
            case(BFMModel)

                call GetBFMPropIndex(Me%ObjBFM, PropertyID, nProperty,STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
                    
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and BFM.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR19a'

                end if  
#endif  
            case(CEQUALW2Model, BenthicCEQUALW2Model)

                call GetCEQUALW2PropIndex(Me%ObjCEQUALW2,PropertyID,nProperty,STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
                    
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and CEQUALW2.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR20'

                end if  

            case(BenthosModel)
                
                call GetBenthosPropIndex(Me%ObjBenthos, PropertyID, nProperty, STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
                    
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and Benthos Model.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR21'

                end if

            case(MacroAlgaeModel)
                
                call GetMacroAlgaePropIndex(Me%ObjMacroAlgae, PropertyID, nProperty, STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 

                    select case(PropertyID)
                        
                        !if it's not a property it can be a rate ID number 
                        case(GrossProd_, NutrientLim_, NLim_, PLim_, LightLim_, TemperatureLim_, SalinityLim_) 

                            nProperty = PropertyID 

                        case default 

                            write(*,*) 
                            write(*,*) 'Inconsistency between Interface and MacroAlgae Model.'
                            stop       'PropertyIndexNumber - ModuleInterface - ERR22'

                    end select  
                    
                end if
                
#ifdef _PHREEQC_
            !ToDo: Need to change
            case (PhreeqCModel)
                
                select case (PropertyID)
                case (Temperature_, pH_, pE_)
                    nProperty = 0
                case default
                    call GetPhreeqCPropIndex(Me%ObjPhreeqC, PropertyID, nProperty, STAT_CALL)
                    
                    if (STAT_CALL .NE. SUCCESS_) then 
                        
                        write(*,*) 
                        write(*,*) 'Inconsistency between Interface and PhreeqC Model.'
                        stop       'PropertyIndexNumber - ModuleInterface - ERR23'

                    end if            
                end select
#endif                        
        end select

        PropertyIndexNumber = nProperty

        !----------------------------------------------------------------------

    end Function PropertyIndexNumber

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    Subroutine KillInterface(InterfaceID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                 :: InterfaceID
        integer, optional, intent(OUT)          :: STAT     

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ready_, nUsers
                                
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_) 

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mINTERFACE_, Me%InstanceID)

            if (nUsers == 0) then

                nUsers = DeassociateInstance (mTIME_,Me%ObjTime)
                if (nUsers == 0) stop 'KillInterface - ModuleInterface - ERR01'

                select case(Me%SinksSourcesModel)
                    
                    case(WaterQualityModel)

                        call KillWaterQuality(Me%ObjWaterQuality, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR02'
                    
                    case(SedimentQualityModel)

                        call KillSedimentQuality(Me%ObjSedimentQuality, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR03'

                    case(LifeModel)

                        call KillLife(Me%ObjLife, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.1'
#ifdef _BFM_  
                    case(BFMModel)

                        call KillBFM(Me%ObjBFM, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.15'
#endif  
                    case(CEQUALW2Model, BenthicCEQUALW2Model)

                        call KillCEQUALW2(Me%ObjCEQUALW2, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.2'

                    case(BenthosModel)
                        
                        call KillBenthos(Me%ObjBenthos, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.3'
                    
                    case(MacroAlgaeModel)
                        
                        call KillMacroAlgae(Me%ObjMacroAlgae, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.4'
                        
#ifdef _PHREEQC_
                    case(PhreeqCModel)
                    
                        call KillPhreeqC (Me%ObjPhreeqC, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.5'
#endif
                end select

                if(associated(Me%Salinity))then
                    deallocate(Me%Salinity, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR05'
                end if
                
                if(associated(Me%Temperature))then
                    deallocate(Me%Temperature, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR06'
                end if

                if(associated(Me%Oxygen))then
                    deallocate(Me%Oxygen, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR06A'
                end if                

                if(associated(Me%ShortWaveAverage))then
                    deallocate(Me%ShortWaveAverage, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR06B'
                end if

                
                if(associated(Me%ShortWaveTop))then
                    deallocate(Me%ShortWaveTop, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR07'
                end if
                
                
                if(associated(Me%LightExtCoefField))then
                    deallocate(Me%LightExtCoefField, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR08'
                end if
                
                if(associated(Me%Thickness))then
                    deallocate(Me%Thickness, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR09'
                end if
                
                if(associated(Me%ShearStress))then
                    deallocate(Me%ShearStress, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR09.1'
                end if

                if(associated(Me%SPMFlux))then
                    deallocate(Me%SPMFlux, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR09.2'
                end if
                
                deallocate(Me%ConcentrationIncrement, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR10'

                deallocate(Me%OpenPoints, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR11'
                
                
                if(associated(Me%FishFood))then
                    deallocate(Me%FishFood, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR12'
                end if

                if(associated(Me%Alkalinity))then
                    deallocate(Me%Alkalinity, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR13'
                end if
                
#ifdef _PHREEQC_

                if(associated(Me%SolutionVolume))then
                    deallocate(Me%SolutionVolume, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR16'
                end if
             
                if(associated(Me%pH))then
                    deallocate(Me%pH, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR17'
                end if
                
                if(associated(Me%pE))then
                    deallocate(Me%pE, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR18'
                end if
                
                if(associated(Me%SolidMass))then
                    deallocate(Me%SolidMass, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR18A'
                end if                
                
#endif                
                
                deallocate(Me%Mass, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR14'


                deallocate(Me%AddedProperties, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR15'

                !Deallocates Instance
                call DeallocateInstance

                InterfaceID = 0

                STAT_ = SUCCESS_

            end if

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_


        !----------------------------------------------------------------------

    end subroutine KillInterface


    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Interface), pointer           :: AuxInterface
        type (T_Interface), pointer           :: PreviousInterface

        !Updates pointers
        if (Me%InstanceID == FirstObjInterface%InstanceID) then
            FirstObjInterface => FirstObjInterface%Next
        else
            PreviousInterface => FirstObjInterface
            AuxInterface      => FirstObjInterface%Next
            do while (AuxInterface%InstanceID /= Me%InstanceID)
                PreviousInterface => AuxInterface
                AuxInterface      => AuxInterface%Next
            enddo

            !Now update linked list
            PreviousInterface%Next => AuxInterface%Next

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


    subroutine Ready (InterfaceID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                            :: InterfaceID
        integer                            :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (InterfaceID > 0) then
            call LocateObjInterface(InterfaceID)
            ready_ = VerifyReadLock (mINTERFACE_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjInterface (InterfaceID)

        !Arguments-------------------------------------------------------------
        integer                            :: InterfaceID

        !Local-----------------------------------------------------------------

        Me => FirstObjInterface
        do while (associated (Me))
            if (Me%InstanceID == InterfaceID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleInterface - LocateObjInterface - ERR01'

    end subroutine LocateObjInterface

    !--------------------------------------------------------------------------

end Module ModuleInterface

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

