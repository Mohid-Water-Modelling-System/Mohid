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
    !griflet
    use ModuleFunctions, only: CHUNK_I, SetMatrixValue
    use ModuleStopWatch, only: StartWatch, StopWatch
    use ModuleWaterQuality
    use ModuleSedimentQuality
    use ModuleCEQUALW2
    use ModuleLife
    use ModuleBenthos
    use ModuleMacroAlgae
    use ModuleEnterData, only: ReadFileName
    use ModuleBenthicEcology
    use ModuleWWTPQ
    use ModuleSeagrassSedimInteraction
    use ModuleSeagrassWaterInteraction
    use ModuleBivalve

#ifdef _PHREEQC_ 
    use ModulePhreeqC
#endif   

#ifdef _BFM_    
    use ModuleBFM
#endif
    !griflet
    !$ use omp_lib
    
    implicit none

    private

    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  :: ConstructInterface
    private ::      AllocateInstance
    private ::      ReadInterfaceFilesName
    private ::      StartSinksSourcesModel
    private ::      AllocateVariables  
    private ::      ConstructMapping                                              
    private ::      Check_Options
    private ::          FindProperty

    !Modifier 
    public  :: Modify_Interface
    private ::      FillMassTempSalinity
    private ::      InputData
    private ::      PropertyIndexNumber
    private ::      UnfoldMatrix
    public  :: SetSOD
    public  :: SetSettlementOnInterface
    public  :: UpdateMassDimensions

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
    
    private :: ConstructMapping_3D
    private :: ConstructMapping_2D
    private :: ConstructMapping_1D
    interface ConstructMapping
        module procedure ConstructMapping_1D
        module procedure ConstructMapping_2D
        module procedure ConstructMapping_3D
    end interface ConstructMapping

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
    private :: FillMassFromWater2D
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
    
    private :: GetRateFlux1D
    private :: GetRateFlux2D
    private :: GetRateFlux3D
    interface  GetRateFlux
        module procedure GetRateFlux1D
        module procedure GetRateFlux2D
        module procedure GetRateFlux3D
    end interface  GetRateFlux
    
    private :: UnfoldMatrix1D_I
    private :: UnfoldMatrix1D_R
    private :: UnfoldMatrix2D_I
    private :: UnfoldMatrix2D_R
    private :: UnfoldMatrix2D_R8
    private :: UnfoldMatrix3D_I
    private :: UnfoldMatrix3D_R
    private :: UnfoldMatrix3D_R8
    interface  UnfoldMatrix
        module procedure UnfoldMatrix1D_I
        module procedure UnfoldMatrix1D_R
        module procedure UnfoldMatrix2D_I
        module procedure UnfoldMatrix2D_R
        module procedure UnfoldMatrix2D_R8
        module procedure UnfoldMatrix3D_I
        module procedure UnfoldMatrix3D_R
        module procedure UnfoldMatrix3D_R8
    end interface  UnfoldMatrix


    !Type----------------------------------------------------------------------
    type       T_External
        type(T_Time)                            :: Now, BeginTime, EndTime
        integer, dimension(:      ), pointer    :: RiverPoints1D    => null()
        integer, dimension(:, :   ), pointer    :: WaterPoints2D    => null()
        integer, dimension(:, :, :), pointer    :: WaterPoints3D    => null()
        integer, dimension(:      ), pointer    :: OpenPoints1D     => null()
        integer, dimension(:, :   ), pointer    :: OpenPoints2D     => null()
        integer, dimension(:, :, :), pointer    :: OpenPoints3D     => null()
        real,    dimension(:      ), pointer    :: DWZ1D            => null()
        real,    dimension(:, :, :), pointer    :: DWZ              => null(), &
                                                   ShearStress      => null(), &
                                                   SPMFlux          => null()
        real(8), dimension(:, :, :), pointer    :: SedimCellVol     => null()
        logical                                 :: Vertical1D = .false.
    end type T_External

    type      T_Interface
        private
        integer                                 :: InstanceID           = null_int !initialization: Jauch
        character(PathLength)                   :: FileName             = null_str !initialization: Jauch        
        character(StringLength)                 :: SinksSourcesModel    = null_str !initialization: Jauch        
        type(T_Size3D)                          :: Size3D
        type(T_Size2D)                          :: Size2D
        type(T_Size1D)                          :: Size1D
        type(T_Size1D)                          :: Array
        type(T_Size1D)                          :: Prop
        type(T_External)                        :: ExternalVar
        
        !griflet: Optimizing interface
        !griflet: start
        integer, pointer, dimension(:,:,:)      :: IJK2Index                => null()
        integer, pointer, dimension(:,:  )      :: IJ2Index                 => null()
        integer, pointer, dimension(:    )      :: I2Index                  => null()
        
        integer, pointer, dimension(:    )      :: Index2I                  => null()
        integer, pointer, dimension(:    )      :: Index2J                  => null()
        integer, pointer, dimension(:    )      :: Index2K                  => null()
        !griflet: end
        
        real,    pointer, dimension(:,:  )      :: Mass                     => null()
        real,    pointer, dimension(:,:  )      :: ConcentrationIncrement   => null()
        real,    pointer, dimension(:    )      :: Salinity                 => null()
        real,    pointer, dimension(:    )      :: Alkalinity               => null()
        real,    pointer, dimension(:    )      :: Temperature              => null()
        real,    pointer, dimension(:    )      :: Oxygen                   => null()
        real,    pointer, dimension(:    )      :: SOD                      => null()
        logical                                 :: UseSOD = .false.         
        real,    pointer, dimension(:    )      :: ShortWaveTop             => null()
        real,    pointer, dimension(:    )      :: ShortWaveAverage         => null()
        real,    pointer, dimension(:    )      :: LightExtCoefField        => null()
        real,    pointer, dimension(:    )      :: Thickness                => null()
        real,    pointer, dimension(:    )      :: SPMFlux                  => null()
        real,    pointer, dimension(:    )      :: ShearStress              => null()
        real,    pointer, dimension(:    )      :: FishFood                 => null()
        real,    pointer, dimension(:    )      :: WaterPercentage          => null()
        real,    pointer, dimension(:    )      :: DissolvedToParticulate   => null()
        real,    pointer, dimension(:    )      :: SoilDryDensity           => null()
        real,    pointer, dimension(:    )      :: pH                       => null()
        real(8), pointer, dimension(:    )      :: WaterVolume1D            => null() !Array with volumes from unfold matrix 
        real,    pointer, dimension(:    )      :: CellArea1D               => null()!Array with volumes from unfold matrix
        real,    pointer, dimension(:    )      :: VelocityModulus1D
        

#ifdef _PHREEQC_        
        real,    pointer, dimension(:    )      :: pE                       => null()
        real,    pointer, dimension(:    )      :: WaterMass                => null()
        real,    pointer, dimension(:    )      :: SolidMass                => null()
        integer, pointer, dimension(:    )      :: PhreeqCID                => null()
#endif        
        real,    pointer, dimension(:    )      :: IonicStrength            => null()
        real,    pointer, dimension(:    )      :: PhosphorusAdsortionIndex => null()
        real,    pointer, dimension(:    )      :: WindVelocity             => null()
        integer, pointer, dimension(:    )      :: OpenPoints               => null()
        logical, pointer, dimension(:    )      :: AddedProperties          => null()
        !real,    pointer, dimension(:    )      :: WaterVol
        real,    pointer, dimension(:,:  )      :: MassInKgFromWater        => null()
        real,    pointer, dimension(:,:  )      :: WaterMassInKgIncrement   => null()
        real,    pointer, dimension(:    )      :: Sediment                 => null()
        real,    pointer, dimension(:    )      :: BottomSWRadiationAverage => null() !Isabella
        real,    pointer, dimension(:    )      :: ShearStress2D            => null()
        real,    pointer, dimension(:    )      :: UptakeNH4NO3w            => null()
        real,    pointer, dimension(:    )      :: UptakePO4w               => null()
        real,    pointer, dimension(:    )      :: UptakeNH4s               => null()
        real,    pointer, dimension(:    )      :: UptakePO4s               => null()
        real,    pointer, dimension(:    )      :: LightFactor              => null()
        real,    pointer, dimension(:    )      :: NintFactor               => null() !Isabella
        real,    pointer, dimension(:    )      :: NintFactorR              => null() !Isabella
        real,    pointer, dimension(:    )      :: PintFactorR              => null() !Isabella
        real,    pointer, dimension(:    )      :: PintFactor               => null() !Isabella 
        real,    pointer, dimension(:    )      :: RootsMort                => null() !Isabella
        real,    pointer, dimension(:    )      :: SeagOccupation           => null()
        real(8), pointer, dimension(:    )      :: SedimCellVol             => null() ! 3d   !Isabella
        real,    pointer, dimension(:    )      :: MacrOccupation           => null()
        
        real,    pointer, dimension(:    )      :: SettlementProbability
      
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
        
        !Instance of ModuleBenthicEcology
        integer                                 :: ObjBenthicEcology    = 0
        
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
        
        !Instance of ModuleWWTPQ
        integer                                 :: ObjWWTPQ             = 0
        
                
        !Instance of ModuleSeagrassesRoots
        integer                                 :: ObjSeagrassSedimInteraction   = 0
        
       !Instance of ModuleSeagrassWaterInteraction
        integer                                 :: ObjSeagrassWaterInteraction   = 0
        
        !Instance of ModuleBivalve
        integer                                 :: ObjBivalve           = 0

        !Collection of instances                
        type(T_Interface),          pointer     :: Next => null()

    end type T_Interface

        
    !Global Module Variables
    type (T_Interface), pointer                 :: FirstObjInterface    => null()
    type (T_Interface), pointer                 :: Me                   => null()

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
#ifdef _PHREEQC_                                    
                                    PhreeqCDatabase,                       &
                                    PhreeqCDatabaseAux,                    &
                                    PhreeqCModelID,                        &                                 
#endif
                                    Size3D,                                &
                                    Vertical1D,                            &
                                    BivalveID,                             &
                                    STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: InterfaceID
        integer                                                 :: TimeID
        character(len=StringLength)                             :: SinksSourcesModel
        integer, dimension(:), pointer, optional                :: PropertiesList
        real, intent (OUT)                                      :: DT
        integer, dimension(:,:,:), pointer                      :: WaterPoints3D
        
#ifdef _PHREEQC_
        character(LEN=*), optional                              :: PhreeqCDatabase
        character(LEN=*), optional                              :: PhreeqCDatabaseAux
        integer, intent(OUT), optional                          :: PhreeqCModelID
#endif

        type(T_Size3D)                                          :: Size3D
        logical,intent (IN),  optional                          :: Vertical1D
        integer,intent (OUT), optional                          :: BivalveID
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
            Me%Array%ILB                    =  1
            Me%Array%IUB                    =  sum(Me%ExternalVar%WaterPoints3D)

            if (present(Vertical1D)) then
                Me%ExternalVar%Vertical1D = Vertical1D
            else
                Me%ExternalVar%Vertical1D = .false.
            endif
            

            !Path to data file
            call ReadInterfaceFilesName

            
            call StartSinksSourcesModel(DT &
#ifdef _PHREEQC_
                , PhreeqCDatabase, PhreeqCDatabaseAux, PhreeqCModelID &
#endif                
            )

!            !Start sinks and sources model
!            call StartSinksSourcesModel(DT)

            if (present(BivalveID)) then
                BivalveID = Me%ObjBivalve
            endif

            !Verify model DT's
            call CheckDT

            if (present(PropertiesList)) then
                !Verify compute options
                call Check_Options(PropertiesList)
            endif

            !Allocate one-index variables global to whole module
            call AllocateVariables

            call ConstructMapping(Me%Size3D)

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
 
             !Verify model DT's
            call CheckDT

            !Verify compute options
            call Check_Options(PropertiesList)

            !Allocate variables global to whole module
            call AllocateVariables

            call ConstructMapping(Me%Size2D)

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
                                    BivalveID,                             &
                                    STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: InterfaceID
        integer                                                 :: TimeID
        character(len=StringLength)                             :: SinksSourcesModel
        integer, dimension(:    ), pointer                      :: PropertiesList
        real,   intent (OUT)                                    :: DT
        integer, dimension(:), pointer                          :: RiverPoints1D
        type(T_Size1D)                                          :: Size1D
        integer,intent (OUT), optional                          :: BivalveID
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
            
            if (present(BivalveID)) then
                BivalveID = Me%ObjBivalve
            endif


            !Verify model DT's
            call CheckDT

            !Verify compute options
            call Check_Options(PropertiesList)

            !Allocate variables global to whole module
            call AllocateVariables
            
            call ConstructMapping(Me%Size1D)

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
        
        !griflet: optimize performance
        !griflet: start
        allocate(Me%Index2I(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR01'

        allocate(Me%Index2J(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR02'

        allocate(Me%Index2K(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR03'
        !griflet: end
        
        allocate(Me%Mass(PropLB:PropUB, ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR04'

        allocate(Me%ConcentrationIncrement(PropLB:PropUB, ArrayLB:ArrayUB), &
                                                     STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR05'

        allocate(Me%Temperature(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR06'

        allocate(Me%OpenPoints(ArrayLB:ArrayUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR07'
        
        allocate(Me%AddedProperties(PropLB:PropUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR08'
        
        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR10'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR11'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR12'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR13'

                allocate(Me%FishFood(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR14'
                       
                Me%Salinity           = FillValueReal
                Me%FishFood           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal

            case (SedimentQualityModel)
                
                allocate (Me%WaterPercentage(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR20'
                
                Me%WaterPercentage = FillValueReal

                allocate (Me%DissolvedToParticulate(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR21'
                
                Me%DissolvedToParticulate = FillValueReal
         
                allocate (Me%SoilDryDensity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR22'
                
                Me%SoilDryDensity = FillValueReal

                allocate (Me%Salinity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR23'
                
                Me%Salinity = FillValueReal

                allocate (Me%pH(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR24'
                
                Me%pH = FillValueReal

                allocate (Me%IonicStrength(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR25'
                
                Me%IonicStrength = FillValueReal

                allocate (Me%PhosphorusAdsortionIndex(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR26'
                
                Me%PhosphorusAdsortionIndex = FillValueReal

                allocate (Me%WindVelocity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR27'
                
                Me%WindVelocity = FillValueReal

                allocate (Me%Oxygen(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR28'
                
                Me%WindVelocity = FillValueReal

            case (CEQUALW2Model)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR30'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR31'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR32'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR33'

                allocate(Me%Alkalinity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR34'
                       
                Me%Salinity           = FillValueReal
                Me%Alkalinity         = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal

            case (LifeModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR40'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR41'

                allocate(Me%ShortWaveAverage(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR42'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR43'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR44'

                allocate(Me%Alkalinity(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR45'
                       
                Me%Salinity           = FillValueReal
                Me%Alkalinity         = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveAverage   = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
                
#ifdef _BFM_    
            case (BFMModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR50'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR51'
        
                allocate(Me%ShortWaveRadiation(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR52'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR53'

                Me%Salinity           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
#endif

            case (BenthicCEQUALW2Model)                
                allocate(Me%Oxygen(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR60'
                
                allocate(Me%SOD (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR61'
            
            case (BenthosModel )
            
                allocate(Me%Oxygen(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR70'
                
            case (BenthicEcologyModel )
            
                allocate(Me%WaterMassInKgIncrement(PropLB:PropUB, ArrayLB:ArrayUB), &
                                                     STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR80'
                 
                allocate (Me%Sediment(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR81'
                             
                allocate(Me%WaterVolume1D(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR82'
                
                allocate(Me%CellArea1D(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR83'
                
                allocate(Me%MassinKgFromWater(PropLB:PropUB, ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR84'
                
                allocate(Me%BottomSWRadiationAverage(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR85'  
                
                allocate(Me%ShearStress2D(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR86'
                
                 allocate(Me%UptakeNH4NO3w(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR87'    
        
               allocate(Me%UptakePO4w(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR88'
                
               allocate(Me%UptakeNH4s(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR89'   
 
               allocate(Me%LightFactor(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR90'
                
                allocate(Me%UptakePO4s(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR91'
                    
                
                Me%WaterMassInKgIncrement         = FillValueReal
                Me%MassinKgFromWater              = FillValueReal
                Me%Sediment                       = FillValueReal
                Me%WaterVolume1D                  = FillValueReal
                Me%CellArea1D                     = FillValueReal
                Me%BottomSWRadiationAverage       = FillValueReal
                Me%ShearStress2D                  = FillValueReal
                Me%UptakeNH4NO3w                  = FillValueReal
                Me%UptakePO4w                     = FillValueReal
                Me%UptakeNH4s                     = FillValueReal
                Me%LightFactor                    = FillValueReal
                Me%UptakePO4s                     = FillValueReal

                
             case (SeagrassSedimInteractionModel)
                    
              !Me%Temperature               = FillValueReal
              
              allocate(Me%SedimCellVol (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR100'
              
            
              allocate(Me%NintFactorR (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR101'
              
              allocate(Me%RootsMort (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR102'
              
              allocate(Me%PintFactorR (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR103'
              
              Me%SedimCellVol      = FillValueReal
              Me%NintFactorR       = FillValueReal
              Me%PintFactorR       = FillValueReal
              Me%RootsMort         = FillValueReal

   case (SeagrassWaterInteractionModel)
                    
              !Me%Temperature               = FillValueReal
              
              allocate(Me%WaterVolume1D (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR104'
              
              allocate(Me%PintFactor (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR105'
            
              allocate(Me%NintFactor (ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR106'
              
              allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR107'
        
              allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR108'

              allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR109'
              
              
              allocate(Me%SeagOccupation(ArrayLB:ArrayUB), STAT = STAT_CALL)
              if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR110' 
            
            
              Me%WaterVolume1D       = FillValueReal
              Me%PintFactor          = FillValueReal
              Me%NintFactor          = FillValueReal
              Me%LightExtCoefField   = FillValueReal
              Me%ShortWaveTop        = FillValueReal
              Me%Thickness           = FillValueReal
              Me%SeagOccupation            = FillValueReal

            case (MacroAlgaeModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR120'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR121'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR122'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR123'
              
                allocate(Me%ShearStress(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR124'
                
                allocate(Me%MacrOccupation(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR75' 
                
                allocate(Me%SPMFlux(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR125'

                Me%Salinity           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
                Me%ShearStress        = FillValueReal
                Me%SPMFlux            = FillValueReal
                Me%MacrOccupation     = FillValueReal


#ifdef _PHREEQC_
            case (PhreeqCModel)
                                               
                allocate (Me%WaterVolume1D(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR130'
                
                allocate (Me%WaterMass(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR131'

                allocate (Me%pH(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR132'
                
                allocate (Me%pE(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR133'

                allocate (Me%Temperature(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR134'
                
                allocate (Me%SolidMass(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR135'
                
                allocate (Me%PhreeqCID(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR136'

                Me%WaterVolume1D = FillValueReal
                Me%WaterMass     = FillValueReal
                Me%pH            = FillValueReal
                Me%pE            = FillValueReal
                Me%Temperature   = FillValueReal
                Me%SolidMass     = FillValueReal
                Me%PhreeqCID     = 0
#endif

            case (WWTPQModel)
                
                !GRiflet: alocar aqui os arrays
                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR140'
        
                allocate(Me%LightExtCoefField(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR141'
        
                allocate(Me%ShortWaveTop(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR142'

                allocate(Me%Thickness(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR143'

                allocate(Me%FishFood(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR144'
                       
                Me%Salinity           = FillValueReal
                Me%FishFood           = FillValueReal
                Me%LightExtCoefField  = FillValueReal
                Me%ShortWaveTop       = FillValueReal
                Me%Thickness          = FillValueReal
            
            case (BivalveModel)

                allocate(Me%Salinity (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR150'

                allocate(Me%WaterVolume1D (ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR151'

                allocate(Me%CellArea1D(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR152'
                
                allocate(Me%VelocityModulus1D(ArrayLB:ArrayUB), STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'AllocateVariables - ModuleInterface - ERR153'
     
                Me%Salinity             = FillValueReal
                Me%WaterVolume1D        = FillValueReal
                Me%CellArea1D           = FillValueReal
                Me%VelocityModulus1D    = FillValueReal
                
            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'AllocateVariables - ModuleInterface - ERR09'
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
                
           case(BenthicEcologyModel)

                Message  = trim('Benthic Ecology Data File')
                call ReadFileName('BENTHICECOLOGY_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR08.1'
                

#ifdef _PHREEQC_
                
            case (PhreeqCModel)
            
                Message = trim('PhreeqC Data File')
                call ReadFileName('PHREEQC_DATA', Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR09' 
            
#endif
            case (WWTPQModel)

                Message  =trim('WWTPQ Data File')
                call ReadFileName('WWTPQ_DATA', Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - ModuleInterface - ERR10'
                                
            case(SeagrassSedimInteractionModel) ! Isabella

                Message  = trim('SeagrassesRoots Data File')
                call ReadFileName('SEAGRASSESROOTS_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR11'
                
            case(SeagrassWaterInteractionModel) ! Isabella

                Message  = trim('SeagrassWaterInteraction Data File')
                call ReadFileName('SEAGRASSESLEAVES_DATA',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR12'     
    
            case(BivalveModel) 
                
                Message  = trim('Bivalve Data File')
                call ReadFileName('BIV_DAT',Me%FileName, Message = Message, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'ReadInterfaceFilesName - Module Interface - ERR13' 
                
            case default
            
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'ReadInterfaceFilesName - Module Interface - ERR14' 
                
        end select

        Me%FileName = trim(Me%FileName)

        !----------------------------------------------------------------------

    end subroutine ReadInterfaceFilesName 

    !--------------------------------------------------------------------------
    
    subroutine StartSinksSourcesModel (DT &
#ifdef _PHREEQC_
                    , PhreeqCDatabase, PhreeqCDatabaseAux, PhreeqCModelID &
#endif                    
                    )

        !Arguments-------------------------------------------------------------
        real, intent(OUT)                       :: DT
#ifdef _PHREEQC_
        character(LEN=*), optional              :: PhreeqCDatabase
        character(LEN=*), optional              :: PhreeqCDatabaseAux
        integer, intent(INOUT), optional        :: PhreeqCModelID
#endif
        
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
                
                
                 case(BenthicEcologyModel)


                !Construct Benthos Model
                call ConstructBenthicEcology (Me%ObjBenthicEcology,                               &
                                   FileName = Me%FileName,                      &
                                   ILB      = Me%Array%ILB,                     &
                                   IUB      = Me%Array%IUB,                     &
                                   STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR023'
 
                !Number of properties involved
                call GetBenthicEcologySize (Me%ObjBenthicEcology,                             &
                                     PropLB = PropLB,                           &
                                     PropUB = PropUB,                           &
                                     STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR024'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from BenthicEcology model to exit as argument to WaterProperties
                call GetDTBenthicEcology(Me%ObjBenthicEcology, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR025'
                
                

#ifdef _PHREEQC_

            case (PhreeqCModel)
            
                !Construct PhreeqC Model
                call StartPhreeqC (Me%ObjPhreeqC,                    &
                                   FileName    = Me%FileName,        &
                                   Database    = PhreeqCDatabase,    &
                                   DatabaseAux = PhreeqCDatabaseAux, &
                                   STAT = STAT_CALL)
                    
                call GetPhreeqCDT(Me%ObjPhreeqC, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR22'
                 
                if (present(PhreeqCModelID)) then
                    PhreeqCModelID = Me%ObjPhreeqC
                endif
                                
#endif
            case (WWTPQModel)

                !Construct WaterQuality Model
                call StartWWTPQ(Me%ObjWWTPQ,                     &
                                       Me%FileName,                            &
                                       STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR23'
 
                !Construct mass fluxes between properties
                call Construct_WWTPQRateFlux(Me%ObjWWTPQ,                  &
                                          Me%Array%ILB,                        &
                                          Me%Array%IUB,                        &
                                          STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR24'
               
                !Number of properties involved
                call GetWWTPQSize(Me%ObjWWTPQ,                   &
                                         PropLB = PropLB,                      &
                                         PropUB = PropUB,                      &
                                         STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR25'

                !Number of properties involved
                Me%Prop%ILB     = PropLB
                Me%Prop%IUB     = PropUB

                !Get DT from WaterQuality model to exit as argument to WaterProperties
                call GetDTWWTPQM(Me%ObjWWTPQ, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR26'
                
                
            case(SeagrassSedimInteractionModel)

                !Construct SeagrassSedimInteractionModel Model
                
                call StartSeagrassSedimInteraction (Me%ObjSeagrassSedimInteraction,                                 &
                                      FileName = Me%FileName,                           &
                                      ILB      = Me%Array%ILB,                          &
                                      IUB      = Me%Array%IUB,                          &
                                      STAT     = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR030'
 
                !Number of properties involved
                call GetSeagrassSedimInteractionSize (Me%ObjSeagrassSedimInteraction,                               &
                                        PropLB = PropLB,                                &
                                        PropUB = PropUB,                                &
                                        STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR031'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from SeagrassesRoots model to exit as argument to WaterProperties
                call GetDTSeagrassSedimInteraction(Me%ObjSeagrassSedimInteraction, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR032'                
                
                
                
            case(SeagrassWaterInteractionModel)

                !Construct SeagrassWaterInteraction Model
                
                call StartSeagrassWaterInteraction (Me%ObjSeagrassWaterInteraction,                                 &
                                      FileName = Me%FileName,                           &
                                      ILB      = Me%Array%ILB,                          &
                                      IUB      = Me%Array%IUB,                          &
                                      STAT     = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR033'
 
                !Number of properties involved
                call GetSeagrassWaterInteractionSize (Me%ObjSeagrassWaterInteraction,                               &
                                        PropLB = PropLB,                                &
                                        PropUB = PropUB,                                &
                                        STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR034'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from SeagrassWaterInteraction model to exit as argument to WaterProperties
                call GetDTSeagrassWaterInteraction(Me%ObjSeagrassWaterInteraction, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR035'           
                

            case(BivalveModel) 


                call GetComputeCurrentTime(Me%ObjTime,                  &
                                           Me%ExternalVar%Now,          &
                                           STAT = STAT_CALL)         
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR118'

                call GetComputeTimeLimits(Me%ObjTime, EndTime = Me%ExternalVar%EndTime,     &      
                                                      BeginTime = Me%ExternalVar%BeginTime, &
                                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR119'
                                                      
                !Construct Bivalve Model
                call ConstructBivalve(Me%ObjBivalve,                            &
                                      FileName = Me%FileName,                   &
                                      BeginTime= Me%ExternalVar%BeginTime,      &
                                      EndTime  = Me%ExternalVar%EndTime,        &
                                      ArraySize= Me%Array,                      &
                                      STAT     = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR120'
 
                !Number of properties involved
                call GetBivalveSize(Me%ObjBivalve,                              &
                                     PropLB = PropLB,                           &
                                     PropUB = PropUB,                           &
                                     STAT   = STAT_CALL)                            
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR121'

                !Number of properties involved
                Me%Prop%ILB = PropLB
                Me%Prop%IUB = PropUB

                !Get DT from Bivalve model to exit as argument to WaterProperties
                call GetDTBivalve(Me%ObjBivalve, DTSecond = DT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSinksSourcesModel - ModuleInterface - ERR122'

            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'StartSinksSourcesModel - ModuleInterface - ERR00' 
        end select

    end subroutine StartSinksSourcesModel    
    
    !--------------------------------------------------------------------------

    subroutine CheckDT
    
        !Local-----------------------------------------------------------------
        real                                                 :: error_aux, aux_factor 
        real                                                 :: run_period, dt_lag
        character(256)                                       :: model_name
        real                                                 :: dt
        type(T_Time)                                         :: end_time
        integer                                              :: status
        
        !----------------------------------------------------------------------
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = status)                    
        if (status /= SUCCESS_) &
            stop 'CheckDT - ModuleInterface - ERR010' 

        call GetComputeTimeLimits(Me%ObjTime, EndTime = end_time, STAT = status)
        if (status /= SUCCESS_) &
            stop 'CheckDT - ModuleInterface - ERR020' 
               
        dt = InterfaceDT()
        
        if (status /= SUCCESS_) then
            write(*,*)
            write(*,*) 'Error when trying to get DT from '//trim(model_name)
            stop 'CheckDT - ModuleInterface - ERR030'
        endif
        
        !Run period in seconds
        run_period = end_time - Me%ExternalVar%Now

        !The run period must be a multiple of the SedimentQuality DT
        aux_factor = run_period / dt
        error_aux  = aux_factor - int(aux_factor)
        
        if (error_aux /= 0) then
            dt_lag = int(error_aux * dt)
            write(*,*) 
            write(*,*) 'DTSECONDS is not multiple of the run period.'
            write(*,*) trim(model_name)//' wont be computed in the last', dt_lag, ' seconds.'
            write(*,*) 'CheckDT - ModuleInterface - WRN010.'
        endif   
        
        call null_time (Me%ExternalVar%Now)     
        !----------------------------------------------------------------------
    
    end subroutine CheckDT

    !--------------------------------------------------------------------------


    subroutine Check_Options(PropertiesList)

        !Arguments-------------------------------------------------------------
        integer, dimension(:), pointer                      :: PropertiesList   

        !External--------------------------------------------------------------
!        type(T_Time)                                         :: EndTime
        logical                                              :: Zoo, Phyto
        logical                                              :: Diatoms 
        logical                                              :: Nitrogen, Phosphorus
        logical                                              :: Silica  
        logical                                              :: Oxygen, BOD
        logical                                              :: Carbon, Sol_Bacteria
        logical                                              :: Bacteria, Ciliate 
        logical                                              :: Larvae
        logical                                              :: Pompools
!        real                                                 :: DT
        integer                                              :: STAT_CALL

        !Local-----------------------------------------------------------------
!        real                                                 :: ErrorAux, auxFactor 
!        real                                                 :: RunPeriod, Dtlag
!        type(T_Size1D)                                       :: Size 
        integer, dimension(:), pointer                       :: CEQUALW2List
        integer, dimension(:), pointer                       :: MacroAlgaeList
        integer                                              :: i,PropLB, PropUB
        integer, dimension(:), pointer                       :: BenthosList, LifeList, BivalveList
        integer, dimension(:), pointer                       :: BenthicEcologyList
        integer, dimension(:), pointer                       :: SeagrassSedimInteractionList
        integer, dimension(:), pointer                       :: SeagrassWaterInteractionList
#ifdef _BFM_  
        integer, dimension(:), pointer                       :: BFMList
#endif
        !----------------------------------------------------------------------
        
!        call GetComputeCurrentTime(Me%ObjTime,                  &
!                                   Me%ExternalVar%Now,          &
!                                   STAT = STAT_CALL)                    
!        if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR00' 
!
!        !Get end time
!        call GetComputeTimeLimits(Me%ObjTime, EndTime = EndTime, STAT = STAT_CALL)
!        if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR01' 
            
        select case (Me%SinksSourcesModel)

            case (WaterQualityModel)
!                !Get water quality model time step
!                call GetDTWQM(Me%ObjWaterQuality, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL .NE. SUCCESS_)stop 'Check_Options - ModuleInterface - ERR02' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the WQ DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'Water Quality wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN01.'
!                endif 

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

!                Larvae is not a waterproperty but lagrangian
!cd3 :           if (Larvae) then
!                    if (.not.FindProperty(PropertiesList, Larvae_))                 &
!                        stop 'WQM needs property Larvae - Check_Options'
!                end if cd3

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

!                !Get SedimentQuality model time step
!                call GetDTSedimentQuality(Me%ObjSedimentQuality, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR04' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'SedimentQuality wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN02.'
!                endif
                
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
            
                if (.not.FindProperty(PropertiesList, AutotrophicPop_))                     &
                    stop 'SQM needs property "autotrophic microorganism population" - Check_Options'         

                if (.not.FindProperty(PropertiesList, HeterotrophicPop_))                     &
                    stop 'SQM needs property "heterotrophic microorganism population" - Check_Options'    

                if (.not.FindProperty(PropertiesList, AnaerobicPop_))                     &
                    stop 'SQM needs property "anaerobic microorganism population" - Check_Options'    
                
                if (Sol_Bacteria) then
                    if (.not.FindProperty(PropertiesList, SolPop_))                     &
                        stop 'SQM needs property "anaerobic microorganism population" - Check_Options'    
                endif
                
            case(CEQUALW2Model, BenthicCEQUALW2Model)
                
!                !Get CEQUALW2 model time step
!                call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR05' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'CEQUALW2 wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN03.'
!                endif
                
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

!                !Get Life model time step
!                call GetDTLife(Me%ObjLife, DT = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR09' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'Life wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN03.'
!                endif
                
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

!                !Get BFM model time step
!                call GetDTBFM(Me%ObjBFM, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR09a' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'Life wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN03a.'
!                endif
                
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

!                !Get Benthos model time step
!                call GetDTBenthos(Me%ObjBenthos, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR13' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'Benthos wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN04.'
!                endif
                
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

!                !Get MacroAlgae model time step
!                call GetDTMacroAlgae(Me%ObjMacroAlgae, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR17' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'MacroAlgae wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN05.'
!                endif
                
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

                case(BenthicEcologyModel)

!                !Get Benthos model time step
!                call GetDTBenthos(Me%ObjBenthos, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR13' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'Benthos wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN04.'
!                endif
                
                call GetBenthicEcologyPropertyList (Me%ObjBenthicEcology, BenthicEcologyList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR21'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, BenthicEcologyList(i))) then
                        write(*,*) 'Property ',GetPropertyName(BenthicEcologyList(i)),' not found'
                        stop 'Properties lists benthic ecology inconsistent  - Check_Options- ModuleInterface- ERR22'    
                    end if
                end do

                call UngetBenthicEcology (Me%ObjBenthicEcology, BenthicEcologyList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR23'


                case(SeagrassSedimInteractionModel)

                !Get SeagrassSedimInteraction model time step
                !call GetDTSeagrassSedimInteraction(Me%ObjSeagrassSedimInteraction, DTSecond = DT, STAT = STAT_CALL)
                !if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR13' 
        
                !Run period in seconds
                !RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the ? DT
                !auxFactor = RunPeriod / DT

                !ErrorAux  = auxFactor - int(auxFactor)
                
                !if (ErrorAux /= 0) then
                !    Dtlag = int(ErrorAux * DT)
                !    write(*,*) 
                !    write(*,*) 'DTSECONDS is not multiple of the run period.'
                !    write(*,*) 'SeagrassSedimInteraction wont be computed in the last', Dtlag, ' seconds.'
                !    write(*,*) 'Check_Options - ModuleInterface - WRN08.'
                !endif
                
                call GetSeagrassSedimInteractionPropertyList (Me%ObjSeagrassSedimInteraction, & 
                     SeagrassSedimInteractionList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR214'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, SeagrassSedimInteractionList(i))) then
                        write(*,*) 'Property ',GetPropertyName(SeagrassSedimInteractionList(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR215'    
                    end if
                end do

                call UngetSeagrassSedimInteraction (Me%ObjSeagrassSedimInteraction, SeagrassSedimInteractionList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR216'
                

              case(SeagrassWaterInteractionModel)

                !Get SeagrassWaterInteraction model time step
                !call GetDTSeagrassWaterInteraction(Me%ObjSeagrassWaterInteraction, DTSecond = DT, STAT = STAT_CALL)
                !if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR217' 
        
                !Run period in seconds
                !RunPeriod = EndTime - Me%ExternalVar%Now

                !The run period must be a multiple of the ? DT
                !auxFactor = RunPeriod / DT

                !ErrorAux  = auxFactor - int(auxFactor)
                
                !if (ErrorAux /= 0) then
                !    Dtlag = int(ErrorAux * DT)
                !    write(*,*) 
                !    write(*,*) 'DTSECONDS is not multiple of the run period.'
                !    write(*,*) 'SeagrassWaterInteraction wont be computed in the last', Dtlag, ' seconds.'
                !    write(*,*) 'Check_Options - ModuleInterface - WRN09.'
                !endif
                
                call GetSeagrassWaterInteractionPropertyList (Me%ObjSeagrassWaterInteraction, &
                     SeagrassWaterInteractionList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR218'
                
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, SeagrassWaterInteractionList(i))) then
                        write(*,*) 'Property ',GetPropertyName(SeagrassWaterInteractionList(i)),' not found'
                        stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR219'    
                    end if
                end do

                call UnGetSeagrassWaterInteraction (Me%ObjSeagrassWaterInteraction, SeagrassWaterInteractionList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR220' 


#ifdef _PHREEQC_

            case (PhreeqCModel)

!                !Get PhreeqC model time step
!                call GetPhreeqCDT(Me%ObjPhreeqC, DTSecond = DT, STAT = STAT_CALL)
!                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR21' 
!        
!                !Run period in seconds
!                RunPeriod = EndTime - Me%ExternalVar%Now
!
!                !The run period must be a multiple of the SedimentQuality DT
!                auxFactor = RunPeriod / DT
!
!                ErrorAux  = auxFactor - int(auxFactor)
!                
!                if (ErrorAux /= 0) then
!                    Dtlag = int(ErrorAux * DT)
!                    write(*,*) 
!                    write(*,*) 'DTSECONDS is not multiple of the run period.'
!                    write(*,*) 'PhreeqC wont be computed in the last', Dtlag, ' seconds.'
!                    write(*,*) 'Check_Options - ModuleInterface - WRN06.'
!                endif
                
                if (.NOT. FindProperty(PropertiesList, Temperature_)) &
                    stop 'PhreeqC needs property "temperature" - Check_Options'

                if (.NOT. FindProperty(PropertiesList, pH_)) &
                    stop 'PhreeqC needs property "ph" - Check_Options'

                if (.NOT. FindProperty(PropertiesList, pE_)) &
                    stop 'PhreeqC needs property "pe" - Check_Options'
                
                
                !Store number of properties involved
                Me%Prop%ILB = 1
                
                if (FindProperty(PropertiesList, SoilDryDensity_)) then               
                    Me%Prop%IUB = SIZE(PropertiesList) - 4
                else
                    Me%Prop%IUB = SIZE(PropertiesList) - 3
                endif
                
#endif

             case (WWTPQModel)

                !Sinks and sources compute options 
                call GetWWTPQOptions(Me%ObjWWTPQ,       Phyto            = Phyto,           &
                                                        Nitrogen         = Nitrogen,        &
                                                        Diatoms          = Diatoms,         &
                                                        Silica           = Silica,          & 
                                                        Phosphorus       = Phosphorus,      &
                                                        Zoo              = Zoo,             &
                                                        BOD              = BOD,             &
                                                        Oxygen           = Oxygen,          &
                                                        Bacteria         = Bacteria,        &
                                                        Ciliate          = Ciliate,         &
                                                        Pompools         = Pompools,        &
                                                        STAT             = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR30' 


                if (Phyto) then 
                    if (.not.FindProperty(PropertiesList, Phytoplankton_))          &
                        stop 'WWTPQM needs property phytoplankton - Check_Options'
                end if

               if (Zoo) then
                    if (.not.FindProperty(PropertiesList, Zooplankton_))            &
                        stop 'WWTPQM needs property zooplankton - Check_Options'
                end if

               if (Phosphorus) then

                    if (.not.FindProperty(PropertiesList, POP_))                    &
                        stop 'WWTPQM needs property particulate organic phosphorus - Check_Options'

                    if (.not.FindProperty(PropertiesList, DOPRefractory_))          &
                        stop 'WWTPQM needs property dissolved refractory organic phosphorus - Check_Options'

                    if (.not. FindProperty(PropertiesList, DOPNon_Refractory_))     &
                        stop 'WWTPQM needs property dissolved non-refractory organic phosphorus - Check_Options'

                    if (.not.FindProperty(PropertiesList, Inorganic_Phosphorus_))   &
                        stop 'WWTPQM needs property inorganic phosphorus - Check_Options'
                        
                end if

                if (Nitrogen) then

                    if (.not.FindProperty(PropertiesList, PON_))                    &
                        stop 'WWTPQM needs property particulate organic nitrogen - Check_Options'
                        
                    if (Bacteria .and. .not.FindProperty(PropertiesList, PONRefractory_))          &
                        stop 'WWTPQM needs property particulate refractory organic nitrogen - Check_Options'

                    if (.not.FindProperty(PropertiesList, DONRefractory_))          &
                        stop 'WWTPQM needs property dissolved refractory organic nitrogen - Check_Options'

                    if (.not.FindProperty(PropertiesList, DONNon_Refractory_))     &
                        stop 'WWTPQM needs property dissolved non-refractory organic nitrogen - Check_Options'

                    if (.not.FindProperty(PropertiesList, Ammonia_))                &
                        stop 'WWTPQM needs property ammonia - Check_Options'

                    if (.not.FindProperty(PropertiesList, Nitrate_))                &
                        stop 'WWTPQM needs property nitrate - Check_Options'

                    if (.not.FindProperty(PropertiesList, Nitrite_))                &
                        stop 'WWTPQM needs property nitrite - Check_Options'
                        
                end if

                if (BOD) then
                    if (.not.FindProperty(PropertiesList, BOD_))                    &
                        stop 'WWTPQM needs property biochemical oxygen demand - Check_Options'
                end if


                if (Bacteria) then
                    if (.not.FindProperty(PropertiesList, Bacteria_))               &
                        stop 'WWTPQM needs property bacteria - Check_Options'
                end if

                if (Ciliate) then
                    if (.not.FindProperty(PropertiesList, Ciliate_))                &
                        stop 'WWTPQM needs property ciliate - Check_Options'
                end if

                if (Diatoms) then
                    if (.not.FindProperty(PropertiesList, Diatoms_))                &
                        stop 'WWTPQM needs property diatoms - Check_Options'
                end if

                if (Silica) then
                    if (.not.FindProperty(PropertiesList, DSilica_))                &
                        stop 'WWTPQM needs property diatoms - Check_Options'

                    if (.not.FindProperty(PropertiesList, BioSilica_))               &
                        stop 'WWTPQM needs property diatoms - Check_Options'
                end if

                !If one wants to run just the age property oxygen is not needed - Frank 09-2003
                if (Phosphorus .or. BOD .or. Ciliate .or. Nitrogen .or. Bacteria) then
                    if (.not.FindProperty(PropertiesList, Oxygen_))                     &
                        stop 'WWTPQM needs property oxygen - Check_Options'
                endif

               if(.NOT.(Pompools)) then

                if (FindProperty(PropertiesList, PON1_) .OR. FindProperty(PropertiesList, PON2_) .OR. &
                    FindProperty(PropertiesList, PON3_) .OR. FindProperty(PropertiesList, PON4_) .OR. &
                    FindProperty(PropertiesList, PON5_) .OR. FindProperty(PropertiesList, POP1_) .OR. &
                    FindProperty(PropertiesList, POP2_) .OR. FindProperty(PropertiesList, POP3_) .OR. &
                    FindProperty(PropertiesList, POP4_) .OR. FindProperty(PropertiesList, POP5_))     &
                        stop 'WWTPQM needs the POM pools option activated - Check_Options'

               end if   

                if (Pompools) then
                 
                 if (Nitrogen) then

                    if (.not.FindProperty(PropertiesList, PON1_))                    &
                        stop 'WWTPQM with POM pools needs property PON1 - Check_Options'

                    if (.not.FindProperty(PropertiesList, PON2_))                    &
                        stop 'WWTPQM with POM pools needs property PON2 - Check_Options'    

                    if (.not.FindProperty(PropertiesList, PON3_))                    &
                        stop 'WWTPQM with POM pools needs property PON3 - Check_Options'

                    if (.not.FindProperty(PropertiesList, PON4_))                    &
                        stop 'WWTPQM with POM pools needs property PON4 - Check_Options'

                    if (.not.FindProperty(PropertiesList, PON5_))                    &
                        stop 'WWTPQM with POM pools needs property PON5 - Check_Options'
                            
                  end if

                  if (Phosphorus) then

                     if (.not.FindProperty(PropertiesList, POP1_))                   &
                        stop 'WWTPQM with POM pools needs property POP1 - Check_Options'

                    if (.not.FindProperty(PropertiesList, POP2_))                    &
                        stop 'WWTPQM with POM pools needs property POP2 - Check_Options'    

                    if (.not.FindProperty(PropertiesList, POP3_))                    &
                        stop 'WWTPQM with POM pools needs property POP3 - Check_Options'

                    if (.not.FindProperty(PropertiesList, POP4_))                    &
                        stop 'WWTPQM with POM pools needs property POP4 - Check_Options'

                    if (.not.FindProperty(PropertiesList, POP5_))                    &
                        stop 'WWTPQM with POM pools needs property POP5 - Check_Options'

                  end if

                  if (.NOT. Zoo) then                  
                    stop 'WWTPQM with POM pools needs property Zoo - Check_Options'                    
                  end if 

                end if
                
            case(BivalveModel) 
                
                !Size%ILB = 1; Size%IUB = 1
                
                !Get number of simulated properties 
                call GetBivalveSize(Me%ObjBivalve, PropLB, PropUB, STAT = STAT_CALL)
                if(STAT_CALL .NE. SUCCESS_) stop 'Check_Options - ModuleInterface - ERR40'

                call GetBivalvePropertyList(Me%ObjBivalve, BivalveList, STAT_CALL)
                if(STAT_CALL .NE. SUCCESS_) stop 'Check_Options - ModuleInterface - ERR50'
                              
                !Number of properties involved
                PropLB = Me%Prop%ILB
                PropUB = Me%Prop%IUB

                do i = PropLB, PropUB
                    if (.not.FindProperty(PropertiesList, BivalveList(i))) then
                        write(*,*) 'Property ',GetPropertyName(BivalveList(i)),' not found in the Bivalve list'
                        write(*,*) 'Please check (Water Properties file) if keyword BIVALVE should be on.'
                              stop 'Properties lists inconsistent  - Check_Options- ModuleInterface- ERR60'    
                    end if
                end do

                call UngetBivalve (Me%ObjBivalve, BivalveList, STAT=STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR70'

           case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                if (STAT_CALL /= SUCCESS_) stop 'Check_Options - ModuleInterface - ERR80'

        end select

!        call null_time   (Me%ExternalVar%Now)

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
!        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                         :: NLB, NUB
        integer                                         :: Index
        integer                                         :: i, j, k
        real,    dimension(:), pointer                  :: RateFlux
        !$ integer                                      :: CHUNK
        real,               dimension(:),       pointer :: LocalRateFlux

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)   
         
!        ILB     = Me%Size3D%ILB     
!        IUB     = Me%Size3D%IUB     
!        JLB     = Me%Size3D%JLB    
!        JUB     = Me%Size3D%JUB    
!        KLB     = Me%Size3D%KLB   
!        KUB     = Me%Size3D%KUB   

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

                    case (WWTPQModel)

                        call GetWWTPQPropRateFlux( Me%ObjWWTPQ,                     &
                                                nFirstProp, nSecondProp,                &
                                                RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR04'
                        
                    case (SeagrassWaterInteractionModel)
                          
                        call GetSeagrassWaterInteractionRateFlux(Me%ObjSeagrassWaterInteraction,      &
                                                   nFirstProp, nSecondProp,                   &
                                                   RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR06' 
                    
                    
                    case (SeagrassSedimInteractionModel)
           
                        call GetSeagrassSedimInteractionRateFlux(Me%ObjSeagrassSedimInteraction,  &
                                                   nFirstProp, nSecondProp,                   &
                                                   RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR07'

                end select

            elseif (present(RateIndex)) then


                call GetCEQUALW2RateFlux( Me%ObjCeQualW2,                               &
                                          RateIndex, RateFlux,                          &
                                          STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR09'


            endif

            
!            !Number indexed to 3D cell in the vector 
!            Index = 0 
!
!            do k = KLB, KUB
!            do j = JLB, JUB
!            do i = ILB, IUB
!
!                if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
!
!                    Index = Index + 1
!                    RateFlux3D(i, j, k) = RateFlux (Index)
!
!                end if
!            end do
!            end do
!            end do
            
            !griflet: start
            NLB = Me%Array%ILB
            NUB = Me%Array%IUB
            !$ CHUNK = CHUNK_I(NLB, NUB)
            !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalRateFlux)
            LocalRateFlux => RateFlux
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do Index = NLB, NUB
                i = Me%Index2I(Index)
                j = Me%Index2J(Index)
                k = Me%Index2K(Index)
                RateFlux3D(i, j, k) = LocalRateFlux (Index)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            !griflet: stop

            select case (Me%SinksSourcesModel)

                case (WaterQualityModel)

                    call UnGetWQPropRateFlux(Me%ObjWaterQuality, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR06' 
                
                case (SedimentQualityModel)

                    call UngetPropRateFlux(Me%ObjSedimentQuality, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR07' 

                case(CEQUALW2Model, BenthicCEQUALW2Model)
               
                    call UnGetCEQUALW2RateFlux(Me%ObjCeQualW2, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR08'
                
                case (MacroAlgaeModel)

                    call UnGetMacroAlgaeRateFlux(Me%ObjMacroAlgae, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR09' 

                case (WWTPQModel)

                    call UnGetWWTPQPropRateFlux(Me%ObjWWTPQ, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR10' 
                    
                case (SeagrassSedimInteractionModel)
                     
                    call UnGetSeagrassSedimInteractionRateFlux(Me%ObjSeagrassSedimInteraction, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR11' 
                
                case (SeagrassWaterInteractionModel)
                  
                    call UnGetSeagrassWaterInteractionRateFlux(Me%ObjSeagrassWaterInteraction, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux3D - ModuleInterface - ERR12' 
                
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
!        integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: NLB, NUB
        integer                                         :: Index
        integer                                         :: i, j, STAT_CALL
        real,    dimension(:), pointer                  :: RateFlux
        !$ integer                                      :: CHUNK
        real,    dimension(:), pointer                  :: LocalRateFlux

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)   

cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

!            ILB     = Me%Size2D%ILB     
!            IUB     = Me%Size2D%IUB     
!            JLB     = Me%Size2D%JLB    
!            JUB     = Me%Size2D%JUB    

            nullify (RateFlux  )
           ! nullify (RateFlux2D)
            RateFlux2D=0.

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
                    
                    
                    case (BenthicEcologyModel)

                    call GetBenthicEcologyRateFlux(Me%ObjBenthicEcology,                              &
                                            nFirstProp, nSecondProp,                    &
                                            RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux2D - ModuleInterface - ERR01.1'
                
            end select

!           !Number indexed to 2D cell in the vector 
!           Index = 0 
!
!           do j = JLB, JUB
!           do i = ILB, IUB
!
!                if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then
!
!                    Index = Index + 1
!                    RateFlux2D(i, j) = RateFlux (Index)
!
!                end if
!            end do
!            end do

            !griflet: start
            NLB = Me%Array%ILB
            NUB = Me%Array%IUB
            !$ CHUNK = CHUNK_I(NLB, NUB)
            !$OMP PARALLEL PRIVATE(Index,i,j,LocalRateFlux)
            LocalRateFlux => RateFlux
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do Index = NLB, NUB
                i = Me%Index2I(Index)
                j = Me%Index2J(Index)
                RateFlux2D(i, j) = LocalRateFlux (Index)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
            !griflet: stop

            select case (Me%SinksSourcesModel)

                case (BenthosModel)

                    call UnGetBenthosRateFlux(Me%ObjBenthos, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux2D - ModuleInterface - ERR02'
                    
                    
              case (BenthicEcologyModel)

                    call UnGetBenthicEcologyRateFlux(Me%ObjBenthicEcology, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux2D - ModuleInterface - ERR02.2'
                
            end select

            STAT_ = SUCCESS_
        else

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetRateFlux2D

    !--------------------------------------------------------------------------

    subroutine GetRateFlux1D(InterfaceID,                            & 
                             FirstProp,                              &
                             SecondProp,                             &
                             RateFlux1D,                             &
                             RiverPoints1D,                          &
                             RateIndex,                              &
                             STAT)

        !Arguments--------------------------------------------------------------
        integer                                         :: InterfaceID
        real,               dimension(:),   pointer     :: RateFlux1D
        integer,            dimension(:),   pointer     :: RiverPoints1D
        integer, optional,  intent(IN)                  :: FirstProp
        integer, optional,  intent(IN)                  :: SecondProp
        integer, optional,  intent(IN)                  :: RateIndex        
        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        integer                                         :: ready_         

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nFirstProp, nSecondProp          
        integer                                         :: NLB, NUB
        integer                                         :: Index
        integer                                         :: i, STAT_CALL
        real,    dimension(:), pointer                  :: RateFlux
        !$ integer                                      :: CHUNK
        real,    dimension(:), pointer                  :: LocalRateFlux

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)   

cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            nullify (RateFlux  )
            RateFlux1D=0.
            
            if (present(FirstProp)) then
                
                !Number indexed to each property 
                nFirstProp  = PropertyIndexNumber(FirstProp )
                nSecondProp = PropertyIndexNumber(SecondProp)

                !WaterPoints2D
                Me%ExternalVar%RiverPoints1D => RiverPoints1D

                select case (Me%SinksSourcesModel)

                    case (WaterQualityModel)

                        call GetWQPropRateFlux( Me%ObjWaterQuality,                     &
                                                nFirstProp, nSecondProp,                &
                                                RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux1D - ModuleInterface - ERR01'

                    case (BenthosModel)

                        call GetBenthosRateFlux(Me%ObjBenthos,                              &
                                                nFirstProp, nSecondProp,                    &
                                                RateFlux, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux1D - ModuleInterface - ERR01'

                    case(LifeModel)

                end select

            elseif (present(RateIndex)) then


                call GetCEQUALW2RateFlux( Me%ObjCeQualW2,                               &
                                          RateIndex, RateFlux,                          &
                                          STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux1D - ModuleInterface - ERR09'


            endif

            NLB = Me%Array%ILB
            NUB = Me%Array%IUB
            !$ CHUNK = CHUNK_I(NLB, NUB)
            !$OMP PARALLEL PRIVATE(Index,i,LocalRateFlux)
            LocalRateFlux => RateFlux
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do Index = NLB, NUB
                i = Me%Index2I(Index)
                RateFlux1D(i) = LocalRateFlux (Index)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL

            select case (Me%SinksSourcesModel)

                case (WaterQualityModel)

                    call UnGetWQPropRateFlux(Me%ObjWaterQuality, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux1D - ModuleInterface - ERR06' 
                
                case(CEQUALW2Model)
               
                    call UnGetCEQUALW2RateFlux(Me%ObjCeQualW2, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux1D - ModuleInterface - ERR08' 
                
                case (BenthosModel)

                    call UnGetBenthosRateFlux(Me%ObjBenthos, RateFlux, STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'GetRateFlux1D - ModuleInterface - ERR02'
                    
                
            end select

            STAT_ = SUCCESS_
        else

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetRateFlux1D

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
                                  
                                  NintFac3D, NintFac3DR, PintFac3D,                     &
                                  RootsMort, PintFac3DR,                                &
                                  SedimCellVol3D,                                       &
                                  CellArea, WaterVolume,  VelocityModulus,              &
                                  SeagOccupation,MacrOccupation,                                    &
#ifdef _PHREEQC_
                                  WaterMass, SolidMass, pE, Temperature,                & 
                                  PhreeqCID,                                            &
#endif
                                  WindVelocity,  Oxygen, DTProp,STAT)
                                 
        !Arguments-------------------------------------------------------------
        integer                                         :: InterfaceID
        integer, intent(IN)                             :: PropertyID
        real,              dimension(:,:,:), pointer    :: Concentration
        real,    optional, dimension(:,:,:), pointer    :: ShortWaveAverage
        real,    optional, dimension(:,:,:), pointer    :: ShortWaveTop
        real,    optional, dimension(:,:,:), pointer    :: LightExtCoefField
        real,    optional, dimension(:,:,:), pointer    :: WaterPercentage
        real,    optional, dimension(:,:,:), pointer    :: DissolvedToParticulate3D
        real,    optional, dimension(:,:,:), pointer    :: SoilDryDensity
        real,    optional, dimension(:,:,:), pointer    :: Salinity
        real,    optional, dimension(:,:,:), pointer    :: pH
        
        real,    optional, dimension(:,:,:), pointer    :: NintFac3D  !Isabella
        real,    optional, dimension(:,:,:), pointer    :: NintFac3DR  !Isabella
        real,    optional, dimension(:,:,:), pointer    :: PintFac3D  !Isabella
        real,    optional, dimension(:,:,:), pointer    :: PintFac3DR  !Isabella
        real(8),    optional, dimension(:,:,:), pointer    :: SedimCellVol3D  !Isabella
        !real(8),    optional, dimension(:,:,:), pointer    :: WaterCellVol3D  !Isabella
        real,    optional, dimension(:,:,:), pointer    :: RootsMort
        real(8), optional, dimension(:,:,:), pointer    :: WaterVolume
        real,    optional, dimension(:,:,:), pointer    :: VelocityModulus
        real,    optional, dimension(:,:  ), pointer    :: CellArea
        real,    optional, dimension(:,:,:), pointer    :: MacrOccupation
        real,    optional, dimension(:,:,:), pointer    :: SeagOccupation
               
#ifdef _PHREEQC_
        real,    optional, dimension(:,:,:), pointer    :: WaterMass
        real,    optional, dimension(:,:,:), pointer    :: SolidMass
        real,    optional, dimension(:,:,:), pointer    :: pE
        real,    optional, dimension(:,:,:), pointer    :: Temperature
        integer, optional, dimension(:,:,:), pointer    :: PhreeqCID
#endif
        real,    optional, dimension(:,:,:), pointer    :: IonicStrength
        real,    optional, dimension(:,:,:), pointer    :: PhosphorusAdsortionIndex
        real,    optional, dimension(:,:,:), pointer    :: WindVelocity
        real,    optional, dimension(:,:,:), pointer    :: Oxygen
        integer,           dimension(:,:,:), pointer    :: WaterPoints3D
        integer, optional, dimension(:,:,:), pointer    :: OpenPoints3D
        real,    optional, dimension(:,:,:), pointer    :: DWZ, ShearStress, SPMFlux
        real,    optional, intent(IN)                   :: DTProp
        integer, optional, intent(OUT)                  :: STAT

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
        integer                                         :: NLB, NUB
!        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB     
        integer                                         :: PropLB, PropUB, ArrayLB, ArrayUB 
        real                                            :: DTProp_
        logical                                         :: Increment
        type(T_Size2D)                                  :: PropArraySize
        !$ integer                                      :: CHUNK
        integer                                         :: LocalnProperty
        real,              dimension(:,:,:), pointer    :: LocalConcentration
        real,              dimension(:,:), pointer      :: LocalConcInc        
        real,              dimension(:,:), pointer      :: LocalMass
        
!        !DEBUG purposes--------------------------------------------------------
!        real :: old_value, new_value
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterface", "Modify_Interface3D")

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

!            ILB     = Me%Size3D%ILB
!            IUB     = Me%Size3D%IUB
!            JLB     = Me%Size3D%JLB
!            JUB     = Me%Size3D%JUB
!            KLB     = Me%Size3D%KLB
!            KUB     = Me%Size3D%KUB

!            !$ CHUNK = CHUNK_J(JLB,JUB)
            
            PropLB  = Me%Prop%ILB
            PropUB  = Me%Prop%IUB

            ArrayLB = Me%Array%ILB
            ArrayUB = Me%Array%IUB
            
            PropArraySize%ILB = PropLB
            PropArraySize%IUB = PropUB
            PropArraySize%JLB = ArrayLB
            PropArraySize%JUB = ArrayUB

            if(present(OpenPoints3D ))Me%ExternalVar%OpenPoints3D     => OpenPoints3D
            if(present(DWZ          ))Me%ExternalVar%DWZ              => DWZ
            if(present(ShearStress  ))Me%ExternalVar%ShearStress      => ShearStress
            if(present(SPMFlux      ))Me%ExternalVar%SPMFlux          => SPMFlux
            
             if(present(SeagOccupation))then
                call UnfoldMatrix(SeagOccupation, Me%SeagOccupation)
            end if
            
            if(present(MacrOccupation))then
                call UnfoldMatrix(MacrOccupation, Me%MacrOccupation)
            end if
            
           
            if(present(NintFac3D))then
                call UnfoldMatrix(NintFac3D, Me%NintFactor)
            end if
            
            if(present(NintFac3DR))then
                call UnfoldMatrix(NintFac3DR, Me%NintFactorR)
            end if
            
            if(present(PintFac3DR))then
                call UnfoldMatrix(PintFac3DR, Me%PintFactorR)
            end if
            
            if(present(RootsMort))then
                call UnfoldMatrix(RootsMort, Me%RootsMort)
            end if
            
            if(present(PintFac3D))then
                call UnfoldMatrix(PintFac3D, Me%PintFactor)
            end if
                 
            if(present(WaterVolume))then
                call UnfoldMatrix(WaterVolume, Me%WaterVolume1D) 
            end if 
            
            if(present(VelocityModulus))then
                call UnfoldMatrix(VelocityModulus, Me%VelocityModulus1D) 
            end if 
            
            if(present(SedimCellVol3D))then
                call UnfoldMatrix(SedimCellVol3D, Me%SedimCellVol) 
            end if  
            
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

                    call UnfoldMatrix(Me%ExternalVar%OpenPoints3D, Me%OpenPoints)

                    !Stores the concentration before changing them
                    call SetMatrixValue(Me%ConcentrationIncrement, PropArraySize, Me%Mass)

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
                            
                            !if forcing with oxygen send it to sediment quality
                            if (present(Oxygen)) then
                                call UnfoldMatrix(Oxygen, Me%Oxygen)
                                
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
                                                     Me%Oxygen,                        &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR055'
                                
                            else
                            
                                call SedimentQuality(Me%ObjSedimentQuality,            &
                                                     Me%Temperature,                   &
                                                     Me%Mass,                          &
                                                     Me%WaterPercentage,               &
                                                     Me%DissolvedToParticulate,        &
                                                     Me%Array%ILB,                     &
                                                     Me%Array%IUB,                     &
                                                     Me%Salinity,                      &
                                                     Me%SoilDryDensity,                &
                                                     Me%pH,                            &
                                                     Me%IonicStrength,                 &
                                                     Me%PhosphorusAdsortionIndex,      &
                                                     Me%WindVelocity,                  &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR06'
                            
                            endif
                            
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
                            call UnfoldMatrix(MacrOccupation,        Me%MacrOccupation  )

                            call ModifyMacroAlgae(ObjMacroAlgaeID       = Me%ObjMacroAlgae,       &
                                                  Temperature           = Me%Temperature,         &
                                                  Salinity              = Me%Salinity,            &
                                                  OpenPoints            = Me%OpenPoints,          &
                                                  ShearStress           = Me%ShearStress,         &
                                                  SPMDepositionFlux     = Me%SPMFlux,             &
                                                  SWRadiation           = Me%ShortWaveTop,        &
                                                  SWLightExctintionCoef = Me%LightExtCoefField,   &
                                                  Thickness             = Me%Thickness,           &
                                                  Occupation            = Me%MacrOccupation,      &
                                                  Mass                  = Me%Mass,                &
                                                  STAT                  = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR08'

#ifdef _PHREEQC_
                        case(PhreeqcModel)
                         
                            call UnfoldMatrix(WaterVolume, Me%WaterVolume1D)
                            call UnfoldMatrix(WaterMass, Me%WaterMass)
                            call UnfoldMatrix(pH, Me%pH)
                            call UnfoldMatrix(pE, Me%pE)
                            call UnfoldMatrix(Temperature, Me%Temperature)
                            call UnfoldMatrix(PhreeqCID, Me%PhreeqCID)
                            
                            if (present(SolidMass)) then
                            
                                call UnfoldMatrix(SolidMass, Me%SolidMass)
                                
                                call ModifyPhreeqC(PhreeqCID = Me%ObjPhreeqC,      &
                                                   PropertiesValues = Me%Mass,     & 
                                                   WaterVolume = Me%WaterVolume1D, &
                                                   WaterMass = Me%WaterMass,       &  
                                                   Temperature = Me%Temperature,   &
                                                   pH = Me%pH,                     &
                                                   pE = Me%pE,                     &
                                                   SolidMass = Me%SolidMass,       &
                                                   CellsArrayLB = Me%Array%ILB,    &
                                                   CellsArrayUB = Me%Array%IUB,    &
                                                   OpenPoints = Me%OpenPoints,     &
                                                   ModelID = Me%PhreeqCID,         & 
                                                   STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) &
                                    stop 'Modify_Interface3D - ModuleInterface - ERR14'

                            else

                                call ModifyPhreeqC(PhreeqCID = Me%ObjPhreeqC,      &
                                                   PropertiesValues = Me%Mass,     & 
                                                   WaterVolume = Me%WaterVolume1D, &
                                                   WaterMass = Me%WaterMass,       &  
                                                   Temperature = Me%Temperature,   &
                                                   pH = Me%pH,                     &
                                                   pE = Me%pE,                     &
                                                   CellsArrayLB = Me%Array%ILB,    &
                                                   CellsArrayUB = Me%Array%IUB,    &
                                                   OpenPoints = Me%OpenPoints,     & 
                                                   ModelID = Me%PhreeqCID,         & 
                                                   STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) &
                                    stop 'Modify_Interface3D - ModuleInterface - ERR14'
                            endif
#endif
                        case(WWTPQModel)

                            call UnfoldMatrix(ShortWaveTop,Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField, Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ,Me%Thickness)


                            call WWTPQ(Me%ObjWWTPQ,                   &
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
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR15'
                            
                         
                      case(SeagrassSedimInteractionModel)
                         !call UnfoldMatrix(Me%ExternalVar%SedimCellVol3D, Me%SedimCellVol)
                         call ModifySeagrassSedimInteraction(ObjSeagrassSedimInteractionID  = Me%ObjSeagrassSedimInteraction,   &
                                                  OpenPoints            = Me%OpenPoints,          &
                                                  Mass                  = Me%Mass,                &
                                                  NintFactorR           = Me%NintFactorR,          &
                                                  SedimCellVol          = Me%SedimCellVol,         &
                                                  RootsMort             = Me%RootsMort,            &
                                                  PintFactorR           = Me%PintFactorR,          &
                                                  STAT                  = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR09'

                      case(SeagrassWaterInteractionModel)
                      
                      
                      call UnfoldMatrix(ShortWaveTop,               Me%ShortWaveTop  )
                      call UnfoldMatrix(LightExtCoefField,          Me%LightExtCoefField   )
                      call UnfoldMatrix(SeagOccupation,             Me%SeagOccupation  )
                      call UnfoldMatrix(Me%ExternalVar%DWZ,         Me%Thickness)
                     
                      
                            
                            
                      call ModifySeagrassWaterInteraction(ObjSeagrassWaterInteractionID= Me%ObjSeagrassWaterInteraction,  &
                                                OpenPoints             = Me%OpenPoints,           &
                                                Mass                   = Me%Mass,                 &
                                                Temperature            = Me%Temperature,          &
                                                NintFactor             = Me%NintFactor,           &
                                                PintFactor             = Me%PintFactor,           &
                                                SWRadiation            = Me%ShortWaveTop,         &
                                                SWLightExctintionCoef  = Me%LightExtCoefField,    &
                                                Thickness              = Me%Thickness,            &
                                                Occupation             = Me%SeagOccupation,       &
                                                WaterCellVol           = Me%WaterVolume1D,        &
                                                STAT                   = STAT_CALL)
                            
                        case(BivalveModel) 
                        
                            call UnfoldMatrix(CellArea, Me%CellArea1D) 


                            call GetComputeCurrentTime(Me%ObjTime,                  &
                                                       Me%ExternalVar%Now,          &
                                                       STAT = STAT_CALL)                    
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR17' 

                            
                            call ModifyBivalve(Me%ObjBivalve,                        &
                                               Me%Temperature,                       &
                                               Me%Salinity,                          &
                                               Me%Mass,                              &
                                               Me%OpenPoints,                        &
                                               WaterVolumeIN = Me%WaterVolume1D,     &
                                               CellAreaIN    = Me%CellArea1D,        &
                                               VelocityIN    = Me%VelocityModulus1D, &
                                               TimeIN        = Me%ExternalVar%Now,   &
                                    STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface3D - ModuleInterface - ERR17a'
                            
                    end select

                    !griflet
                    !$ CHUNK = CHUNK_I(PropLB, PropUB)
                    !$OMP PARALLEL PRIVATE(prop,index,LocalMass,LocalConcInc)
                    LocalMass => Me%Mass
                    LocalConcInc => Me%ConcentrationIncrement
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do7 :               do prop  = PropLB,  PropUB
do6 :               do index = ArrayLB, ArrayUB
                        LocalConcInc(prop, index) = LocalMass(prop, index) - &
                                                                 LocalConcInc(prop, index)
                    end do do6
                    end do do7
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL

                end if cd4

            elseif (Increment) then cd5

                Index = 0

                nProperty = PropertyIndexNumber(PropertyID)

                if (nProperty < 0) &
                    stop 'Modify_Interface3D - ModuleInterface - ERR80'

                DT = InterfaceDT()

!               !griflet
!               !$OMP PARALLEL PRIVATE(k,j,i,Index)
!do1 :          do k = KLB, KUB
!               !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!do2 :          do j = JLB, JUB
!do3 :          do i = ILB, IUB
!cd2 :                       if (Me%ExternalVar%WaterPoints3D(i, j, k) == 1) then
!                                Index = Index + 1
!                                !Concentrations are only actualized in OpenPoints because of instability
!                                !in waterpoints that are not openpoints
!                                if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then
!                                    Concentration(i, j, k) = Concentration( i, j, k)      + &
!                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
!                                end if
!                            end if cd2
!               end do do3
!               end do do2
!               !$OMP END DO NOWAIT
!               end do do1
!               !$OMP END PARALLEL

                select case (Me%SinksSourcesModel)
                
                    case (BivalveModel)
                    
                        !griflet: start
                        NLB = Me%Array%ILB
                        NUB = Me%Array%IUB
                        !$ CHUNK = CHUNK_I(NLB, NUB)
                        !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalnProperty,LocalConcInc,LocalConcentration)
                        LocalnProperty = nProperty
                        LocalConcInc => Me%ConcentrationIncrement
                        LocalConcentration => Concentration
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

                        do Index = NLB, NUB
                            i = Me%Index2I(Index)
                            j = Me%Index2J(Index)
                            k = Me%Index2K(Index)

                                LocalConcentration(i, j, k) = LocalConcentration( i, j, k)      + &
                                                        LocalConcInc(LocalnProperty, Index) * DTProp / DT 
                                                        
                                if(abs(LocalConcentration(i, j, k)) .le. AlmostZero)then
                                    LocalConcentration(i,j,k) = 0.0
                                end if                       

                        enddo
                        !$OMP END DO NOWAIT
                        !$OMP END PARALLEL
                        !griflet: stop
                    
                    case default
                    
                        !griflet: start
                        NLB = Me%Array%ILB
                        NUB = Me%Array%IUB
                        !$ CHUNK = CHUNK_I(NLB, NUB)
                        !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalnProperty,LocalConcInc,LocalConcentration)
                        LocalnProperty = nProperty
                        LocalConcInc => Me%ConcentrationIncrement
                        LocalConcentration => Concentration
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)

                        do Index = NLB, NUB
                            i = Me%Index2I(Index)
                            j = Me%Index2J(Index)
                            k = Me%Index2K(Index)
                            if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then
                                LocalConcentration(i, j, k) = LocalConcentration( i, j, k)      + &
                                        LocalConcInc(LocalnProperty, Index) * DTProp / DT 
                                        
                                if(abs(LocalConcentration(i, j, k)) .le. AlmostZero)then
                                    LocalConcentration(i,j,k) = 0.0
                                end if

                            end if
                        enddo
                        !$OMP END DO NOWAIT
                        !$OMP END PARALLEL
                        !griflet: stop

                end select
                       
            end if cd5

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_
            
        if (MonitorPerformance) call StopWatch ("ModuleInterface", "Modify_Interface3D")

        !----------------------------------------------------------------------

    end subroutine Modify_Interface3D

    !--------------------------------------------------------------------------

    subroutine Modify_Interface2D(InterfaceID, PropertyID, Concentration,               &
                                  WaterPoints2D, OpenPoints2D, Oxygen2D,                &
                                  MassInKgFromWater, Sediment,                          &
                                  WaterVolume2D,CellArea2D,ShortWave2D,                 &
                                  ShearStress2D,UptakeNH4s2D,UptakeNH4NO3w2D,           &
                                  UptakePO4w2D,LightFactor2D, UptakePO4s2D,             &
                                  DTProp, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: InterfaceID
        integer,            intent(IN)                  :: PropertyID
        real, optional,     intent(IN)                  :: DTProp
        real,              dimension(:,:  ), pointer    :: Concentration         
        integer,           dimension(:,:  ), pointer    :: WaterPoints2D
        integer, optional, dimension(:,:  ), pointer    :: OpenPoints2D
        real   , optional, dimension(:,:  ), pointer    :: Oxygen2D
        real   , optional, dimension(:,:  ), pointer    :: ShortWave2D
        real   , optional, dimension(:,:  ), pointer    :: Sediment
        real(8), optional, dimension(:,:  ), pointer    :: WaterVolume2D
        real   , optional, dimension(:,:  ), pointer    :: CellArea2D
        real   , optional, dimension(:,:  ), pointer    :: ShearStress2D
        real   , optional, dimension(:,:  ), pointer    :: UptakeNH4s2D
        real   , optional, dimension(:,:  ), pointer    :: UptakePO4s2D
        real   , optional, dimension(:,:  ), pointer    :: UptakeNH4NO3w2D
        real   , optional, dimension(:,:  ), pointer    :: UptakePO4w2D
        real   , optional, dimension(:,:  ), pointer    :: LightFactor2D
        real   , optional, dimension(:,:  ), pointer    :: MassInKgFromWater
        integer, optional,  intent(OUT)                 :: STAT

        !External--------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL
        logical                                         :: ReadyToCompute

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: nProperty
        integer                                         :: Index 
        integer                                         :: prop, JulDay
        !integer                                         :: ILB, IUB, JLB, JUB
        integer                                         :: i, j  
        integer                                         :: NLB, NUB
        integer                                         :: PropLB, PropUB, ArrayLB, ArrayUB 
        real                                            :: DTProp_, DT
        logical                                         :: Increment
        
        !$ integer                                      :: CHUNK
        integer                                         :: LocalnProperty
        real,              dimension(:,:), pointer      :: LocalConcentration
        real,              dimension(:,:), pointer      :: LocalConcInc        
        real,              dimension(:,:), pointer      :: LocalMass

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterface", "Modify_Interface2D")

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

!            ILB     = Me%Size2D%ILB
!            IUB     = Me%Size2D%IUB
!            JLB     = Me%Size2D%JLB
!            JUB     = Me%Size2D%JUB
!
!            !$ CHUNK = CHUNK_J(JLB,JUB)
            
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
            
            if(present(WaterVolume2D))then
                call UnfoldMatrix(WaterVolume2D, Me%WaterVolume1D)
            end if 
            
             if(present(CellArea2D))then
                call UnfoldMatrix(CellArea2D, Me%CellArea1D)
            end if  
            
            if(present(ShortWave2D))then
                call UnfoldMatrix(ShortWave2D, Me%BottomSWRadiationAverage)
            end if  
            
            if(present(ShearStress2D))then
                call UnfoldMatrix(ShearStress2D, Me%ShearStress2D)
            end if  
            
            if(present(Sediment))then
                call UnfoldMatrix(Sediment, Me%Sediment)
            end if 
            
            if(present(UptakeNH4s2D))then
                call UnfoldMatrix(UptakeNH4s2D, Me%UptakeNH4s)
            end if   
            
            if(present(UptakeNH4NO3w2D))then
                call UnfoldMatrix(UptakeNH4NO3w2D, Me%UptakeNH4NO3w)
            end if   
            
           if(present(UptakePO4w2D))then
                call UnfoldMatrix(UptakePO4w2D, Me%UptakePO4w)
            end if   
           
            if(present(UptakePO4s2D))then
                call UnfoldMatrix(UptakePO4s2D, Me%UptakePO4s)
            end if 
            
           if(present(LightFactor2D))then
                call UnfoldMatrix(LightFactor2D, Me%LightFactor)
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
                
                
                if(present(MassInKgFromWater) )then
                
                call FillMassFromWater2D(PropertyID, MassInKgFromWater)
                
                endif

cd4 :           if (ReadyToCompute) then

                    !$ CHUNK = CHUNK_I(PropLB, PropUB)

                    call UnfoldMatrix(Me%ExternalVar%OpenPoints2D, Me%OpenPoints)

                    !Stores the concentration before changing them
                    Me%ConcentrationIncrement = Me%Mass
                    
                    
                    if(Me%SinksSourcesModel==BenthicEcologyModel)then
                    !Stores the concentration before changing them
                    Me%WaterMassInKgIncrement = Me%MassInKgFromWater
                    endif

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
                            
                            
                            
                        case(BenthicEcologyModel)
                        
                        call GetComputeCurrentTime(Me%ObjTime,                  &
                                                       Me%ExternalVar%Now,          &
                                                       STAT = STAT_CALL)                    
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERRXX'  
                        
                        call JulianDay(Me%ExternalVar%Now, JulDay)
                 
                        
                            call ModifyBenthicEcology  (Me%ObjBenthicEcology,             &
                                                 Me%Temperature,                          &
                                                 Me%WaterVolume1D,                        &
                                                 Me%CellArea1D,                           &
                                                 Me%MassInKgFromWater,                    &
                                                 Me%Sediment,                             &
                                                 Me%BottomSWRadiationAverage,             &
                                                 Me%ShearStress2D,                        &
                                                 Me%UptakeNH4s,                           &
                                                 Me%UptakeNH4NO3w,                        &
                                                 Me%UptakePO4w,                           &
                                                 Me%UptakePO4s,                           &
                                                 Me%LightFactor,                          &
                                                 JulDay,                                  &
                                                 Me%OpenPoints,                           &
                                                 Me%Mass,                                 &
                                                 STAT = STAT_CALL)
                                                 
                                                 
                                   
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface2D - ModuleInterface - ERR03'


                            !$OMP PARALLEL PRIVATE(prop, index)
                            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
                            do prop  = PropLB,  PropUB
                            do index = ArrayLB, ArrayUB
                                Me%WaterMassInKgIncrement(prop, index) = Me%MassInKgFromWater(prop, index) - &
                                                                         Me%WaterMassInKgIncrement(prop, index)
                            end do 
                            end do 
                            !$OMP END DO NOWAIT
                            !$OMP END PARALLEL

                    end select

                    !$OMP PARALLEL PRIVATE(prop,index,LocalMass,LocalConcInc)
                    LocalMass => Me%Mass
                    LocalConcInc => Me%ConcentrationIncrement
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do7 :               do prop  = PropLB,  PropUB
do6 :               do index = ArrayLB, ArrayUB
                        LocalConcInc(prop, index) = LocalMass(prop, index) - &
                                                                 LocalConcInc(prop, index)
                    end do do6
                    end do do7
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL
                
                end if cd4

            elseif (Increment) then cd5

                Index = 0

                nProperty = PropertyIndexNumber(PropertyID)

                if (nProperty < 0) &
                    stop 'Modify_Interface2D - ModuleInterface - ERR80'

                DT = InterfaceDT()

!                !griflet
!                !$OMP PARALLEL PRIVATE(j,i,Index)            
!                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
!                do j = JLB, JUB
!                do i = ILB, IUB
!                    if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then
!                        Index = Index + 1
!                        
!                        if (Me%ExternalVar%OpenPoints2D(i, j) == 1) then
!                            Concentration(i, j) = Concentration(i, j)      + &
!                            Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
!                        end if
!
!                    end if
!                end do
!                end do
!                !$OMP END DO NOWAIT
!                !$OMP END PARALLEL
                
               if(Me%SinksSourcesModel==BenthicEcologyModel)then
                 
                        NLB = Me%Array%ILB
                        NUB = Me%Array%IUB
                        do Index = NLB, NUB
                            i = Me%Index2I(Index)
                            j = Me%Index2J(Index)
                            if (Me%ExternalVar%OpenPoints2D(i, j) == 1) then
                                MassInKgFromWater(i, j) = MassInKgFromWater(i, j)      + &
                                Me%WaterMassInKgIncrement(nProperty, Index) * DTProp / DT 
                            end if
                        enddo 
                
                endif
             
             
             
                !griflet: start
                NLB = Me%Array%ILB
                NUB = Me%Array%IUB
                !$ CHUNK = CHUNK_I(NLB, NUB)
                !$OMP PARALLEL PRIVATE(Index,i,j,LocalnProperty,LocalConcInc,LocalConcentration)
                LocalnProperty = nProperty
                LocalConcInc => Me%ConcentrationIncrement
                LocalConcentration => Concentration
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do Index = NLB, NUB
                    i = Me%Index2I(Index)
                    j = Me%Index2J(Index)
                    if (Me%ExternalVar%OpenPoints2D(i, j) == 1) then
                        LocalConcentration(i, j) = LocalConcentration( i, j)      + &
                                LocalConcInc(LocalnProperty, Index) * DTProp / DT 
                    end if
                enddo
                !$OMP END DO NOWAIT
                !$OMP END PARALLEL
                !griflet: stop
                        
            end if cd5

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_
            
        if (MonitorPerformance) call StopWatch ("ModuleInterface", "Modify_Interface2D")

        !----------------------------------------------------------------------

    end subroutine Modify_Interface2D
    
    !--------------------------------------------------------------------------

    subroutine Modify_Interface1D(InterfaceID, PropertyID, Concentration,            &
                                  RiverPoints1D, OpenPoints1D, DWZ,                  &
                                  ShortWaveAverage, ShortWaveTop, LightExtCoefField, &
                                  WaterPercentage, DissolvedToParticulate1D,         &
                                  SoilDryDensity, Salinity, pH, IonicStrength,       &
                                  PhosphorusAdsortionIndex, WindVelocity,            &
                                  Oxygen1D, WaterVolume,                             &
                                  DTProp, STAT)

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
        !integer                                         :: ILB, IUB
        integer                                         :: NLB, NUB
        integer                                         :: PropLB, PropUB, ArrayLB, ArrayUB 
        real                                            :: DTProp_
        logical                                         :: Increment
        !$ integer                                      :: CHUNK
        integer                                         :: LocalnProperty
        real,              dimension(:), pointer        :: LocalConcentration
        real,              dimension(:,:), pointer      :: LocalConcInc

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleInterface", "Modify_Interface1D")

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

!            ILB     = Me%Size1D%ILB
!            IUB     = Me%Size1D%IUB
!
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

                            if (present(Oxygen1D)) then
                                call UnfoldMatrix(Oxygen1D, Me%Oxygen)

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
                                                     Me%Oxygen,                         &
                                                     STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR055'
                            
                            else
                            
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
                            
                            endif
                            
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

                        case(WWTPQModel)

                            call UnfoldMatrix(ShortWaveTop,         Me%ShortWaveTop)
                            call UnfoldMatrix(LightExtCoefField,    Me%LightExtCoefField)
                            call UnfoldMatrix(Me%ExternalVar%DWZ1D, Me%Thickness)


                            call WWTPQ(Me%ObjWWTPQ,                   &
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
                            if (STAT_CALL /= SUCCESS_) stop 'Modify_Interface1D - ModuleInterface - ERR10'
                            
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

                DT = InterfaceDT()

!do3 :          do i = ILB, IUB
!cd2 :              if (Me%ExternalVar%RiverPoints1D(i) == 1) then
!                                Index = Index + 1
!                                !Concentrations are only actualized in OpenPoints because of instability
!                                !in waterpoints that are not openpoints
!                                if (Me%ExternalVar%OpenPoints1D(i) == 1) then
!                                    Concentration(i) = Concentration( i)      + &
!                                    Me%ConcentrationIncrement(nProperty, Index) * DTProp / DT 
!                                end if
!                   end if cd2
!               end do do3
                        !griflet: start
                        NLB = Me%Array%ILB
                        NUB = Me%Array%IUB
                        !$ CHUNK = CHUNK_I(NLB,NUB)
                        !$OMP PARALLEL PRIVATE(Index,i,LocalConcentration,LocalConcInc,LocalnProperty)
                        LocalConcentration => Concentration
                        LocalConcInc => Me%ConcentrationIncrement
                        LocalnProperty = nProperty
                        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                        do Index = NLB, NUB
                            i = Me%Index2I(Index)
                            !Concentrations are only actualized in OpenPoints because of instability
                            !in waterpoints that are not openpoints
                            if (Me%ExternalVar%OpenPoints1D(i) == 1) then
                                LocalConcentration(i) = LocalConcentration( i)      + &
                                LocalConcInc(LocalnProperty, Index) * DTProp / DT 
                            end if
                        enddo
                        !$OMP END DO NOWAIT
                        !$OMP END PARALLEL
                        !griflet: stop

            end if cd5

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))STAT = STAT_
            
        if (MonitorPerformance) call StopWatch ("ModuleInterface", "Modify_Interface1D")

        !----------------------------------------------------------------------

    end subroutine Modify_Interface1D

    !--------------------------------------------------------------------------
    
    subroutine UpdateMassDimensions(InterfaceID, PropertiesList, STAT)
    
        !Arguments-------------------------------------------------------------        
        integer                                 :: InterfaceID
        integer, dimension(:), pointer          :: PropertiesList
        integer, optional                       :: STAT
        
        !External--------------------------------------------------------------        
        integer                                 :: STAT_, ready_,STAT_CALL

        !Local-----------------------------------------------------------------
        integer                                 :: i
        integer                                 :: PropLB, PropUB
        integer                                 :: ArrayLB, ArrayUB
        integer, dimension(:), pointer          :: BivalveList

        !----------------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModuleInterface", "UpdateMass")

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            ArrayLB = Me%Array%ILB
            ArrayUB = Me%Array%IUB

            !Get number of simulated properties 
            call GetBivalveSize(Me%ObjBivalve, PropLB, PropUB, STAT = STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_) stop 'UpdateMassDimensions - ModuleInterface - ERR01'

            call GetBivalvePropertyList(Me%ObjBivalve, BivalveList, STAT_CALL)
            if(STAT_CALL .NE. SUCCESS_) stop 'UpdateMassDimensions - ModuleInterface - ERR02'
                          
            !Number of properties involved
            Me%Prop%ILB = PropLB 
            Me%Prop%IUB = PropUB

            do i = PropLB, PropUB
                if (.not.FindProperty(PropertiesList, BivalveList(i))) then
                    write(*,*) 'Property ',GetPropertyName(BivalveList(i)),' not found in the Bivalve list'
                    write(*,*) 'Please check subroutine Construct_CohortPropertiesFromCohort if BIVALVE is on.'
                          stop 'Properties lists inconsistent  - UpdateMassDimensions - ModuleInterface - ERR03'    
                end if
            end do

            call UngetBivalve (Me%ObjBivalve, BivalveList, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'UpdateMassDimensions - ModuleInterface - ERR04'
            
            deallocate(Me%Mass)
            allocate(Me%Mass(PropLB:PropUB, ArrayLB:ArrayUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'UpdateMassDimensions - ModuleInterface - ERR05'

            deallocate(Me%ConcentrationIncrement)
            allocate(Me%ConcentrationIncrement(PropLB:PropUB, ArrayLB:ArrayUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'UpdateMassDimensions - ModuleInterface - ERR06'
            
            deallocate(Me%AddedProperties)
            allocate(Me%AddedProperties(PropLB:PropUB), STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)stop 'UpdateMassDimensions - ModuleInterface - ERR07'
            
            Me%Mass                       = FillValueReal
            Me%ConcentrationIncrement     = FillValueReal
            Me%AddedProperties            = .false.
              
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_
            
        if (MonitorPerformance) call StopWatch ("ModuleInterface", "UpdateMass")

    end subroutine UpdateMassDimensions
    
    
    !--------------------------------------------------------------------------
    
    subroutine SetSettlementOnInterface (InterfaceID, SpeciesIDNumber, SettlementProbability, STAT)
    
        !Arguments-------------------------------------------------------------        
        integer                                 :: InterfaceID
        integer, intent(in)                     :: SpeciesIDNumber
        real   , dimension(:,:,:), pointer      :: SettlementProbability        
        integer, optional                       :: STAT
        
        !External--------------------------------------------------------------        
        integer                                 :: STAT_, ready_,STAT_CALL
    
        !----------------------------------------------------------------------
       
        if (MonitorPerformance) call StartWatch ("ModuleInterface", "SetSettlementOnInterface")

        STAT_ = UNKNOWN_

        call Ready(InterfaceID, ready_)    

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            allocate(Me%SettlementProbability(Me%Array%ILB:Me%Array%IUB))

            call UnfoldMatrix(SettlementProbability, Me%SettlementProbability)
            
            call SetSettlementOnBivalve(Me%ObjBivalve, SpeciesIDNumber, Me%SettlementProbability, STAT = STAT_CALL)
            
            deallocate(Me%SettlementProbability)
            nullify   (Me%SettlementProbability)
            
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_
            
        if (MonitorPerformance) call StopWatch ("ModuleInterface", "SetSettlementOnInterface")

    end subroutine SetSettlementOnInterface
    
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
        integer                                 :: numAutotrophicPop, numHeterotrophicPop, numAnaerobicPop
        integer                                 :: numSolPop
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
        
!#ifdef _PHREEQC_
!        logical                                 :: pHAdded = .false.
!        logical                                 :: pEAdded = .false.
!!        logical                                 :: SoilDryDensityAdded = .false.
!#endif        
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


cd1110 :        if (lPompools) then

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
                                    AutotrophicPop                  = numAutotrophicPop,        &
                                    HeterotrophicPop                = numHeterotrophicPop,      &
                                    AnaerobicPop                    = numAnaerobicPop,          &
                                    SolPop                          = numSolPop,                & 
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



    cd27 :      if (lNitrogen) then

    cd28 :          if (PropertyID== PON_                  ) then
                        call InputData(Concentration, numLabil_OM_N)
                        Me%AddedProperties(numLabil_OM_N)         = .TRUE.
                    end if cd28

    cd29 :          if (PropertyID== RefreactaryOrganicN_  ) then
                        call InputData(Concentration, numRefractOM_N)
                        Me%AddedProperties(numRefractOM_N)        = .TRUE.
                    end if cd29

    cd30 :          if (PropertyID== Ammonia_              ) then
                        call InputData(Concentration,numAmmonia)
                        Me%AddedProperties(numAmmonia)            = .TRUE.
                    end if cd30

    cd31 :          if (PropertyID== Nitrate_              ) then
                        call InputData(Concentration,numNitrate)
                        Me%AddedProperties(numNitrate)            = .TRUE.
                    end if cd31

    cd32 :          if (PropertyID== Ngas_                 ) then
                        call InputData(Concentration,numNgas)
                        Me%AddedProperties(numNgas)               = .TRUE.
                    end if cd32

    cd33 :          if (PropertyID== HeterotrophicN_       ) then
                        call InputData(Concentration,numHeterotrophicN)
                        Me%AddedProperties(numHeterotrophicN)      = .TRUE.
                    end if cd33

    cd34 :          if (PropertyID== AutotrophicN_         ) then
                        call InputData(Concentration,numAutotrophicN)
                        Me%AddedProperties(numAutotrophicN)    = .TRUE.
                    end if cd34

    cd35 :          if (PropertyID== AnaerobicN_           ) then
                        call InputData(Concentration,numAnaerobicN)
                        Me%AddedProperties(numAnaerobicN)         = .TRUE.
                    end if cd35

    cd350 :         if (PropertyID== Urea_                 ) then
                        call InputData(Concentration,numUrea)
                        Me%AddedProperties(numUrea)               = .TRUE.
                    end if cd350

    !cd351 :         if (PropertyID== AmmoniaGas_           ) then
    !                    call InputData(Concentration,numAmmoniaGas)
    !                    Me%AddedProperties(numAmmoniaGas)         = .TRUE.
    !                end if cd351

    cd352 :         if (lSol_Bacteria) then

    cd353 :             if (PropertyID== SolubilizingN_   ) then
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

cd4210 :        if (PropertyID== AutotrophicPop_)then
                    call InputData(Concentration,numAutotrophicPop)
                    Me%AddedProperties(numAutotrophicPop)         = .TRUE.
                end if cd4210

cd422 :          if (PropertyID== HeterotrophicPop_)then
                    call InputData(Concentration,numHeterotrophicPop)
                    Me%AddedProperties(numHeterotrophicPop)       = .TRUE.
                end if cd422

cd423 :          if (PropertyID== AnaerobicPop_)then
                    call InputData(Concentration,numAnaerobicPop)
                    Me%AddedProperties(numAnaerobicPop)           = .TRUE.
                end if cd423
                
                if (lSol_Bacteria) then
    cd424 :          if (PropertyID== SolPop_)then
                        call InputData(Concentration,numSolPop)
                        Me%AddedProperties(numSolPop)           = .TRUE.
                    end if cd424
                endif
                
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

                call GetPhreeqCPropIndex(Me%ObjPhreeqC, PropertyID, IndexNumber, STAT=STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) = .true.
                    
                else if ((STAT_CALL .NE. SUCCESS_) .AND. (STAT_CALL .NE. NOT_FOUND_ERR_)) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR07' 
                end if

                Ready = .true.

                do i = PropLB, PropUB                       
                    if (.NOT. Me%AddedProperties(i)) then                        
                        Ready = .false.
                        exit                            
                    end if                         
                end do 
            
                if (Ready) Me%AddedProperties = .false.

!                select case (PropertyID)
!                    case (Temperature_)
!                        call UnfoldMatrix(Concentration, Me%Temperature)
!                        TemperatureAdded = .TRUE.                    
!                    case (pH_)
!                        call UnfoldMatrix(Concentration, Me%pH)
!                        pHAdded = .TRUE.                                        
!                    case (pE_)
!                        call UnfoldMatrix(Concentration, Me%pE)
!                        pEAdded = .TRUE. 
!!                    case (SoilDryDensity_)
!!                        !SoilDryDensity will be used only in these cases:
!!                        if ((PhreeqCSimOptions%Exchanger .EQ. 1) .OR. (PhreeqCSimOptions%SolidPhase)) then
!!                            call UnfoldMatrix(Concentration, Me%SoilDryDensity)
!!                        endif
!!                        SoilDryDensityAdded = .TRUE.                                                            
!                    case default
!                        call GetPhreeqCPropIndex(Me%ObjPhreeqC, PropertyID, IndexNumber, STAT=STAT_CALL)
!                        if (STAT_CALL == SUCCESS_) then
!                        
!                            call InputData(Concentration, IndexNumber)
!                            Me%AddedProperties(IndexNumber) = .true.
!                            
!                        else if ((STAT_CALL .NE. SUCCESS_) .AND. (STAT_CALL .NE. NOT_FOUND_ERR_)) then
!                            stop 'FillMassTempSalinity3D - ModuleInterface - ERR07' 
!                        end if
!                end select            
!                
!                if (TemperatureAdded .AND. pHAdded .AND. pEAdded) then ! .AND. SoilDryDensityAdded) then
!                    Ready = .TRUE.
!
!                    do i = PropLB, PropUB                       
!                        if (.NOT. Me%AddedProperties(i)) then                        
!                            Ready = .false.
!                            exit                             
!                        end if                         
!                    end do 
!                
!                    if (Ready) Me%AddedProperties = .false.
!                    
!                end if
#endif
            case (WWTPQModel)

                call GetWWTPQPropIndex(Me%ObjWWTPQ,                                         &
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
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity3D - ModuleInterface - ERR08' 

                call GetWWTPQOptions(Me%ObjWWTPQ,                                           &
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
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity3D - ModuleInterface - ERR09'

                 if (lPhyto) then 
                 if (PropertyID == Phytoplankton_) then
                     call InputData(Concentration,numPhyto)
                     Me%AddedProperties(numPhyto)    = .TRUE.
                 end if 
                 end if 

                 if (lZoo) then
                 if (PropertyID == Zooplankton_) then
                     call InputData(Concentration,numZoo)
                     Me%AddedProperties(numZoo)      = .TRUE.
                 end if 
                 end if 

                 if (lLarvae) then
                 if (PropertyID == Larvae_) then
                     call InputData(Concentration,numLarvae)
                     Me%AddedProperties(numLarvae)   = .TRUE.
                 end if 
                 end if 

                 if (lPhosphorus) then

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

                     if (PropertyID== Inorganic_Phosphorus_) then
                         call InputData(Concentration, numInorganicPhosphorus)
                         Me%AddedProperties(numInorganicPhosphorus)      = .TRUE.
                     end if
                     
                 end if


                 if (lNitrogen) then
                     if (PropertyID== PON_) then
                         call InputData(Concentration, numPartOrganicNitrogen)
                         Me%AddedProperties(numPartOrganicNitrogen)      = .TRUE.
                     end if

                     if (PropertyID== PONRefractory_) then
                         call InputData(Concentration, numPartOrganicNitrogenRef)
                         Me%AddedProperties(numPartOrganicNitrogenRef)   = .TRUE.
                     end if

                     if (PropertyID== DONRefractory_) then
                         call InputData(Concentration, numDONRefractory)
                         Me%AddedProperties(numDONRefractory)            = .TRUE.
                     end if

                     if (PropertyID== DONNon_Refractory_) then
                         call InputData(Concentration,numDONNonRefractory)
                         Me%AddedProperties(numDONNonRefractory)         = .TRUE.
                     end if

                     if (PropertyID== Ammonia_) then
                         call InputData(Concentration,numAmmonia)
                         Me%AddedProperties(numAmmonia)                  = .TRUE.
                     end if

                     if (PropertyID== Nitrate_) then
                         call InputData(Concentration,numNitrate)
                         Me%AddedProperties(numNitrate)                  = .TRUE.
                     end if                                                              

                  if (PropertyID== Nitrite_) then
                      call InputData(Concentration,numNitrite)
                      Me%AddedProperties(numNitrite)                  = .TRUE.
                  end if
              end if


             if (lPompools) then

               if (lNitrogen) then
                  
                  if (PropertyID== PON1_) then
                      call InputData(Concentration, numPONitrogen1)
                      Me%AddedProperties(numPONitrogen1)              = .TRUE.
                  end if                   

                   if (PropertyID== PON2_) then
                      call InputData(Concentration, numPONitrogen2)
                      Me%AddedProperties(numPONitrogen2)              = .TRUE.
                  end if 
                  
                  if (PropertyID== PON3_) then
                      call InputData(Concentration, numPONitrogen3)
                      Me%AddedProperties(numPONitrogen3)              = .TRUE.
                  end if
                  
                   if (PropertyID== PON4_) then
                      call InputData(Concentration, numPONitrogen4)
                      Me%AddedProperties(numPONitrogen4)              = .TRUE.
                  end if                                                            

                   if (PropertyID== PON5_) then
                      call InputData(Concentration, numPONitrogen5)
                      Me%AddedProperties(numPONitrogen5)              = .TRUE.
                  end if    
               
               end if  

               if (lPhosphorus) then
               
                   if (PropertyID== POP1_) then
                      call InputData(Concentration, numPOPhosphorus1)
                      Me%AddedProperties(numPOPhosphorus1)            = .TRUE.
                  end if                    

                   if (PropertyID== POP2_) then
                      call InputData(Concentration, numPOPhosphorus2)
                      Me%AddedProperties(numPOPhosphorus2)            = .TRUE.
                  end if 
                  
                   if (PropertyID== POP3_) then
                      call InputData(Concentration, numPOPhosphorus3)
                      Me%AddedProperties(numPOPhosphorus3)            = .TRUE.
                  end if 
                  
                   if (PropertyID== POP4_) then
                      call InputData(Concentration, numPOPhosphorus4)
                      Me%AddedProperties(numPOPhosphorus4)            = .TRUE.
                  end if                                                             

                   if (PropertyID== POP5_) then
                      call InputData(Concentration, numPOPhosphorus5)
                      Me%AddedProperties(numPOPhosphorus5)            = .TRUE.
                  end if                 
               
               end if
                                                         
              end if


              if (lBOD) then
              if (PropertyID== BOD_) then
                  call InputData(Concentration,numBOD)
                  Me%AddedProperties(numBOD)      = .TRUE.
              end if
              end if



              if (lBacteria) then
              if (PropertyID== Bacteria_) then
                  call InputData(Concentration,numBacteria)
                  Me%AddedProperties(numBacteria) = .TRUE.
              end if
              end if



              if (lCiliate) then
              if (PropertyID== Ciliate_) then
                  call InputData(Concentration,numCiliate)
                  Me%AddedProperties(numCiliate)  = .TRUE.
              end if
              end if

              !The oxygen is always compute. If it's not defined in the waterproperties 
              ! module then the saturation value is compute
              if (PropertyID== Oxygen_) then
                  call InputData(Concentration,numOxygen)
                  Me%AddedProperties(numOxygen)   = .TRUE.
              end if

              if (lSalinity) then
              
              if (PropertyID== Salinity_) then
                      call UnfoldMatrix(Concentration, Me%Salinity)
                      SalinityAdded       =.TRUE.
                  end if
              end if

              if (PropertyID== Temperature_) then
                  call UnfoldMatrix(Concentration, Me%Temperature)
                  TemperatureAdded    =.TRUE.
              end if

              if (PropertyID== FishFood_) then
                  call UnfoldMatrix(Concentration, Me%FishFood)
                  FishFoodAdded       =.TRUE.
              end if

              if (lDiatoms) then 
              if (PropertyID == Diatoms_) then
                  call InputData(Concentration,numDiatoms)
                  Me%AddedProperties(numDiatoms)    = .TRUE.
              end if
              end if


               if (lSilica) then

                if (PropertyID == DSilica_) then
                    call InputData(Concentration, numSiDiss)
                    Me%AddedProperties(numSiDiss)    = .TRUE.
                end if

                if (PropertyID== BioSilica_) then
                    call InputData(Concentration, numSiBio)
                    Me%AddedProperties(numSiBio)            = .TRUE.
                end if

            end if

            if ((SalinityAdded .or. (.not. lSalinity)) .AND. TemperatureAdded) then
                Ready = .TRUE.

                do i = PropLB, PropUB
                if (.NOT. Me%AddedProperties(i)) then
                         Ready = .FALSE.
                         exit
                     end if
                     end do

                     if (Ready) Me%AddedProperties = .FALSE.
                 end if       
     
     
            case(SeagrassSedimInteractionModel)
                
                call GetSeagrassSedimInteractionPropIndex(Me%ObjSeagrassSedimInteraction,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR010' 
                end if
              

                Ready = .TRUE.
                
                
                do i = PropLB, PropUB
                if (.NOT. Me%AddedProperties(i)) then
                        Ready = .FALSE.
                        exit 
                end if 
                end do 
                if (Ready) Me%AddedProperties = .FALSE.
                
            case(SeagrassWaterInteractionModel)
                
                call GetSeagrassWaterInteractionPropIndex(Me%ObjSeagrassWaterInteraction,PropertyID,IndexNumber,STAT=STAT_CALL)
                if (STAT_CALL==SUCCESS_) then
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR011' 
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
                
            case(BivalveModel)
                            
                call GetBivalvePropIndex(Me%ObjBivalve,PropertyID,IndexNumber,STAT=STAT_CALL)
                
                if (STAT_CALL==SUCCESS_) then
                
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                    
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity3D - ModuleInterface - ERR12' 
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

            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'FillMassTempSalinity3D - ModuleInterface - ERR13'

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
                
             case(BenthicEcologyModel)

                call GetBenthicEcologyPropIndex(Me%ObjBenthicEcology, PropertyID, IndexNumber, STAT = STAT_CALL)

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
    
    !----------------------------------------------------------------------
    
    subroutine FillMassFromWater2D(PropertyID, MassFromWater) 

        !Arguments-------------------------------------------------------------
        integer,                intent(IN )     :: PropertyID
        real, dimension(:,:),   pointer         :: MassFromWater

        !External--------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: IndexNumber, STAT_CALL
        integer                                 :: i,j
        integer                                 :: Index
        integer                                 :: NLB, NUB
        !----------------------------------------------------------------------

       select case (Me%SinksSourcesModel)

           case(BenthicEcologyModel)

                call GetBenthicEcologyPropIndex(Me%ObjBenthicEcology, PropertyID, IndexNumber, STAT = STAT_CALL)


        end select
        
        if (STAT_CALL==SUCCESS_) then
 

                    !griflet: new way
                    !griflet: start
                    NLB = Me%Array%ILB
                    NUB = Me%Array%IUB
                    !$OMP PARALLEL PRIVATE(Index,i,j)
                    !$OMP DO
                    do Index = NLB, NUB
                        i = Me%Index2I(Index)
                        j = Me%Index2J(Index)
                        Me%MassinKgFromWater(IndexNumber,Index) = MassFromWater(i,j)
                    enddo
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL
                    !griflet: stop
        
        
        else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
            stop 'FillMassTempSalinity2D - ModuleInterface - ERR05' 
        end if  

        !----------------------------------------------------------------------

    end Subroutine FillMassFromWater2D

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
        integer                                 :: numAutotrophicPop, numHeterotrophicPop, numAnaerobicPop
        integer                                 :: numSolPop
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
                                AutotrophicPop                  = numAutotrophicPop,        &
                                HeterotrophicPop                = numHeterotrophicPop,      &
                                AnaerobicPop                    = numAnaerobicPop,          & 
                                SolPop                          = numSolPop,                &
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


cd4210 :        if (PropertyID== AutotrophicPop_)then
                    call InputData(Concentration,numAutotrophicPop)
                    Me%AddedProperties(numAutotrophicPop)         = .TRUE.
                end if cd4210

cd422 :          if (PropertyID== HeterotrophicPop_)then
                    call InputData(Concentration,numHeterotrophicPop)
                    Me%AddedProperties(numHeterotrophicPop)       = .TRUE.
                end if cd422

cd423 :          if (PropertyID== AnaerobicPop_)then
                    call InputData(Concentration,numAnaerobicPop)
                    Me%AddedProperties(numAnaerobicPop)           = .TRUE.
                end if cd423
                
                if (lSol_Bacteria) then
    cd424 :          if (PropertyID== SolPop_)then
                        call InputData(Concentration,numSolPop)
                        Me%AddedProperties(numSolPop)            = .TRUE.
                    end if cd424
                endif
                
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

            case (WWTPQModel)

                call GetWWTPQPropIndex(Me%ObjWWTPQ,                                         &
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
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity1D - ModuleInterface - ERR08' 

                call GetWWTPQOptions(Me%ObjWWTPQ,                                           &
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
                if (STAT_CALL .NE. SUCCESS_)stop 'FillMassTempSalinity1D - ModuleInterface - ERR09'


          if (lPhyto) then 
          if (PropertyID == Phytoplankton_) then
              call InputData(Concentration,numPhyto)
              Me%AddedProperties(numPhyto)    = .TRUE.
          end if
          end if

          if (lZoo) then
          if (PropertyID == Zooplankton_) then
              call InputData(Concentration,numZoo)
              Me%AddedProperties(numZoo)      = .TRUE.
          end if
          end if

         if (lLarvae) then
         if (PropertyID == Larvae_) then
              call InputData(Concentration,numLarvae)
              Me%AddedProperties(numLarvae)   = .TRUE.
          end if
          end if

          if (lPhosphorus) then

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

              if (PropertyID== Inorganic_Phosphorus_) then
                  call InputData(Concentration, numInorganicPhosphorus)
                  Me%AddedProperties(numInorganicPhosphorus)      = .TRUE.
              end if
              
          end if


          if (lNitrogen) then
              if (PropertyID== PON_) then
                  call InputData(Concentration, numPartOrganicNitrogen)
                  Me%AddedProperties(numPartOrganicNitrogen)      = .TRUE.
              end if

              if (PropertyID== PONRefractory_) then
                  call InputData(Concentration, numPartOrganicNitrogenRef)
                  Me%AddedProperties(numPartOrganicNitrogenRef)   = .TRUE.
              end if

              if (PropertyID== DONRefractory_) then
                  call InputData(Concentration, numDONRefractory)
                  Me%AddedProperties(numDONRefractory)            = .TRUE.
              end if

              if (PropertyID== DONNon_Refractory_) then
                  call InputData(Concentration,numDONNonRefractory)
                  Me%AddedProperties(numDONNonRefractory)         = .TRUE.
              end if

              if (PropertyID== Ammonia_) then
                  call InputData(Concentration,numAmmonia)
                  Me%AddedProperties(numAmmonia)                  = .TRUE.
              end if 

              if (PropertyID== Nitrate_) then
                  call InputData(Concentration,numNitrate)
                  Me%AddedProperties(numNitrate)                  = .TRUE.
              end if                    

             if (PropertyID== Nitrite_) then
                  call InputData(Concentration,numNitrite)
                  Me%AddedProperties(numNitrite)                  = .TRUE.
              end if
          end if


        if (lPompools) then

          if (lNitrogen) then
              
             if (PropertyID== PON1_) then
                  call InputData(Concentration, numPONitrogen1)
                  Me%AddedProperties(numPONitrogen1)              = .TRUE.
              end if                    

              if (PropertyID== PON2_) then
                  call InputData(Concentration, numPONitrogen2)
                  Me%AddedProperties(numPONitrogen2)              = .TRUE.
              end if
              
             if (PropertyID== PON3_) then
                  call InputData(Concentration, numPONitrogen3)
                  Me%AddedProperties(numPONitrogen3)              = .TRUE.
              end if
              
              if (PropertyID== PON4_) then
                  call InputData(Concentration, numPONitrogen4)
                  Me%AddedProperties(numPONitrogen4)              = .TRUE.
              end if                                                              

              if (PropertyID== PON5_) then
                  call InputData(Concentration, numPONitrogen5)
                  Me%AddedProperties(numPONitrogen5)              = .TRUE.
              end if   
           
           end if   

          if (lPhosphorus) then
           
              if (PropertyID== POP1_) then
                  call InputData(Concentration, numPOPhosphorus1)
                  Me%AddedProperties(numPOPhosphorus1)            = .TRUE.
              end if                   

              if (PropertyID== POP2_) then
                  call InputData(Concentration, numPOPhosphorus2)
                  Me%AddedProperties(numPOPhosphorus2)            = .TRUE.
              end if 
              
              if (PropertyID== POP3_) then
                  call InputData(Concentration, numPOPhosphorus3)
                  Me%AddedProperties(numPOPhosphorus3)            = .TRUE.
              end if 
              
             if (PropertyID== POP4_) then
                  call InputData(Concentration, numPOPhosphorus4)
                  Me%AddedProperties(numPOPhosphorus4)            = .TRUE.
              end if                                                             

              if (PropertyID== POP5_) then
                  call InputData(Concentration, numPOPhosphorus5)
                  Me%AddedProperties(numPOPhosphorus5)            = .TRUE.
              end if                 
           
           end if
                                                     
          end if




          if (lBOD) then
          if (PropertyID== BOD_) then
              call InputData(Concentration,numBOD)
              Me%AddedProperties(numBOD)      = .TRUE.
          end if
          end if



          if (lBacteria) then
          if (PropertyID== Bacteria_) then
              call InputData(Concentration,numBacteria)
              Me%AddedProperties(numBacteria) = .TRUE.
          end if
          end if



          if (lCiliate) then
          if (PropertyID== Ciliate_) then
              call InputData(Concentration,numCiliate)
              Me%AddedProperties(numCiliate)  = .TRUE.
          end if
          end if


          !The oxygen is always compute. If it's not defined in the waterproperties 
          ! module then the saturation value is compute
          if (PropertyID== Oxygen_) then
              call InputData(Concentration,numOxygen)
              Me%AddedProperties(numOxygen)   = .TRUE.
          end if

          if (lSalinity) then
          
              if (PropertyID== Salinity_) then
                  call UnfoldMatrix(Concentration, Me%Salinity)
                  SalinityAdded       =.TRUE.
              end if
          end if

          if (PropertyID== Temperature_) then
              call UnfoldMatrix(Concentration, Me%Temperature)
              TemperatureAdded    =.TRUE.
          end if

          if (PropertyID== FishFood_) then
              call UnfoldMatrix(Concentration, Me%FishFood)
              FishFoodAdded       =.TRUE.
          end if


        if (lDiatoms) then 
        if (PropertyID == Diatoms_) then
              call InputData(Concentration,numDiatoms)
              Me%AddedProperties(numDiatoms)    = .TRUE.
          end if
          end if


          if (lSilica) then

              if (PropertyID == DSilica_) then
                  call InputData(Concentration, numSiDiss)
                  Me%AddedProperties(numSiBio)    = .TRUE.
              end if

              if (PropertyID== BioSilica_) then
                  call InputData(Concentration, numSiBio)
                  Me%AddedProperties(numSiBio)            = .TRUE.
              end if

          end if

          if ((SalinityAdded .or. (.not. lSalinity)) .AND. TemperatureAdded) then
              Ready = .TRUE.

              do i = PropLB, PropUB
              if (.NOT. Me%AddedProperties(i)) then
                  Ready = .FALSE.
                  exit
              end if
              end do 

                    if (Ready) Me%AddedProperties = .FALSE.
                end if      
        
            case(BivalveModel)
                            
                call GetBivalvePropIndex(Me%ObjBivalve,PropertyID,IndexNumber,STAT=STAT_CALL)
                
                if (STAT_CALL==SUCCESS_) then
                
                    call InputData(Concentration, IndexNumber)
                    Me%AddedProperties(IndexNumber) =.TRUE.
                    
                else if (STAT_CALL.NE.SUCCESS_ .AND. STAT_CALL.NE.NOT_FOUND_ERR_) then
                    stop 'FillMassTempSalinity1D - ModuleInterface - ERR10' 
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

            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'FillMassTempSalinity1D - ModuleInterface - ERR11'
        
        end select

        !----------------------------------------------------------------------

    end subroutine FillMassTempSalinity1D
    
    !griflet: optimizating performance
    !griflet: start
    !--------------------------------------------------------------------------
    subroutine ConstructMapping_3D(Size3D)
        
        !Arguments ------------------------------------------------------------
        type(T_Size3D)                          :: Size3D

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j, k
        integer                                 :: ILB, IUB
        integer                                 :: JLB, JUB
        integer                                 :: KLB, KUB
        
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------

        ILB = Size3D%ILB
        IUB = Size3D%IUB

        JLB = Size3D%JLB
        JUB = Size3D%JUB

        KLB = Size3D%KLB
        KUB = Size3D%KUB
        
        !griflet: optimize performance
        !griflet: start
        allocate(Me%IJK2Index(ILB:IUB,JLB:JUB,KLB:KUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructMapping_3D - ModuleInterface - ERR01-B1'    
        !griflet: end

        !Number indexed to 3D cell in the vector
        Index = 0

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
                Index = Index + 1
                Me%IJK2Index(i,j,k) = Index
                Me%Index2I(Index) = i
                Me%Index2J(Index) = j
                Me%Index2K(Index) = k
            endif
        enddo
        enddo
        enddo

        if ((Index) > Me%Array%IUB) stop 'ConstructMapping_3D - ModuleInterface - ERR01'

    end subroutine ConstructMapping_3D

   !----------------------------------------------------------------------
    subroutine ConstructMapping_2D(Size2D)
        
        !Arguments ------------------------------------------------------------
        type(T_Size2D)                          :: Size2D

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j
        integer                                 :: ILB, IUB
        integer                                 :: JLB, JUB

        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------

        ILB = Size2D%ILB
        IUB = Size2D%IUB

        JLB = Size2D%JLB
        JUB = Size2D%JUB

        !griflet: optimize performance
        !griflet: start
        allocate(Me%IJ2Index(ILB:IUB,JLB:JUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'ConstructInterface3D - ModuleInterface - ERR01-B1'    
       !griflet: end

        !Number indexed to 3D cell in the vector
        Index = 0

        do j = JLB, JUB
        do i = ILB, IUB
            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
                Index = Index + 1
                Me%IJ2Index(i,j) = Index
                Me%Index2I(Index) = i
                Me%Index2J(Index) = j
            endif
        enddo
        enddo

        if ((Index) > Me%Array%IUB) stop 'ConstructMapping_3D - ModuleInterface - ERR01'

    end subroutine ConstructMapping_2D

   !----------------------------------------------------------------------

    subroutine ConstructMapping_1D(Size1D)
        
        !Arguments ------------------------------------------------------------
        type(T_Size1D)                          :: Size1D

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i
        integer                                 :: ILB, IUB

        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------

        ILB = Size1D%ILB
        IUB = Size1D%IUB

        !griflet: optimize performance
        !griflet: start
        allocate(Me%I2Index(ILB:IUB), STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructInterface3D - ModuleInterface - ERR01-B1'    
        !griflet: end

        !Number indexed to 3D cell in the vector
        Index = 0

        do i = ILB, IUB
            if (Me%ExternalVar%RiverPoints1D(i)==1) then
                Index = Index + 1
                Me%I2Index(i) = Index
                Me%Index2I(Index) = i
            endif
        enddo

        if ((Index) > Me%Array%IUB) stop 'ConstructMapping_3D - ModuleInterface - ERR01'

    end subroutine ConstructMapping_1D

   !----------------------------------------------------------------------
   
   !griflet: end optimization
   
   !----------------------------------------------------------------------
   
    subroutine InputData3D (Concentration, nProperty)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer         :: Concentration
        integer, intent (in)                    :: nProperty 

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j, k
        real, dimension(:,:,:), pointer         :: LocalConcentration
        integer                                 :: LocalnProperty 
        !$ integer                              :: CHUNK
        !integer                                 :: ILB, IUB      
        !integer                                 :: JLB, JUB     
        !integer                                 :: KLB, KUB    
        
        integer                                 :: NLB, NUB

        !----------------------------------------------------------------------

!        ILB = Me%Size3D%ILB     
!        IUB = Me%Size3D%IUB     
!
!        JLB = Me%Size3D%JLB    
!        JUB = Me%Size3D%JUB    
!
!        KLB = Me%Size3D%KLB   
!        KUB = Me%Size3D%KUB   
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do k = KLB, KUB
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
!                Index = Index + 1
!                Me%Mass(nProperty,Index) = Concentration(i,j,k)
!            endif    
!        enddo
!        enddo
!        enddo
!
!        if ((Index)> Me%Array%IUB) stop 'InputData3D - ModuleInterface - ERR01'

        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalnProperty,LocalConcentration)
        LocalnProperty = nProperty
        LocalConcentration => Concentration
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            k = Me%Index2K(Index)
            Me%Mass(LocalnProperty,Index) = LocalConcentration(i,j,k)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop
        
        !----------------------------------------------------------------------

    end subroutine InputData3D

    !--------------------------------------------------------------------------

    subroutine InputData2D (Concentration, nProperty)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer           :: Concentration
        integer, intent (in)                    :: nProperty 

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: LocalnProperty
        real, dimension(:,:), pointer           :: LocalConcentration
        integer                                 :: i, j
        !integer                                 :: ILB, IUB      
        !integer                                 :: JLB, JUB  
        integer                                 :: NLB, NUB   
        !$ integer                              :: CHUNK

        !----------------------------------------------------------------------

        !ILB = Me%Size2D%ILB     
        !IUB = Me%Size2D%IUB     

        !JLB = Me%Size2D%JLB    
        !JUB = Me%Size2D%JUB    
        
        !!Number indexed to 3D cell in the vector 
        !Index = 0

        !do j = JLB, JUB
        !do i = ILB, IUB
        !    if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
        !        Index = Index + 1
        !        Me%Mass(nProperty,Index) = Concentration(i,j)
        !    endif
        !enddo
        !enddo

        !if ((Index)> Me%Array%IUB) stop 'InputData2D - ModuleInterface - ERR01'

        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,LocalnProperty,LocalConcentration)
        LocalnProperty = nProperty
        LocalConcentration => Concentration
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            Me%Mass(LocalnProperty,Index) = LocalConcentration(i,j)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop

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
        !integer                                 :: ILB, IUB   
        
        integer                                 :: NLB, NUB      
        !$ integer                              :: CHUNK
        real, dimension(:), pointer             :: LocalConcentration        
        !----------------------------------------------------------------------

        !ILB = Me%Size1D%ILB     
        !IUB = Me%Size1D%IUB
        
        !!Number indexed to 1D cell in the vector 
        !Index = 0

        !do i = ILB, IUB
        !    if (Me%ExternalVar%RiverPoints1D(i) == WaterPoint) then
        !        Index = Index + 1
        !        Me%Mass(nProperty,Index) = Concentration(i)
        !    endif
        !enddo

        !if ((Index)> Me%Array%IUB) stop 'InputData1D - ModuleInterface - ERR01'

        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB,NUB)
        !$OMP PARALLEL PRIVATE(Index,i,LocalConcentration)
        LocalConcentration => Concentration
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            Me%Mass(nProperty,Index) = LocalConcentration(i)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop

        !----------------------------------------------------------------------

    end subroutine InputData1D

    !--------------------------------------------------------------------------

    subroutine UnfoldMatrix3D_R (Matrix3D, Vector)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:,:,:), pointer     :: Matrix3D
        real(4), dimension(:    ), pointer     :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i, j, k
        integer                             :: NLB, NUB
        !$ integer                          :: CHUNK
        real(4), dimension(:,:,:), pointer     :: LocalMatrix3D
!        integer                             :: ILB, IUB      
!        integer                             :: JLB, JUB     
!        integer                             :: KLB, KUB    
!
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size3D%ILB     
!        IUB = Me%Size3D%IUB     
!
!        JLB = Me%Size3D%JLB    
!        JUB = Me%Size3D%JUB    
!
!        KLB = Me%Size3D%KLB   
!        KUB = Me%Size3D%KUB   
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do k = KLB, KUB
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
!                Index         = Index + 1
!                Vector(Index) = Matrix3D(i,j,k)
!            endif    
!        enddo
!        enddo
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix3D_R - ModuleInterface - ERR01'
        
        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalMatrix3D)
        LocalMatrix3D => Matrix3D
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            k = Me%Index2K(Index)
            Vector(Index) = LocalMatrix3D(i,j,k)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop
            
        !----------------------------------------------------------------------

    end subroutine UnfoldMatrix3D_R

    !--------------------------------------------------------------------------
    
        !--------------------------------------------------------------------------

    subroutine UnfoldMatrix3D_R8 (Matrix3D, Vector)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:,:,:), pointer  :: Matrix3D
        real(8), dimension(:    ), pointer  :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i, j, k
        integer                             :: NLB, NUB
        !$ integer                          :: CHUNK
        real(8), dimension(:,:,:), pointer  :: LocalMatrix3D
!        integer                             :: ILB, IUB      
!        integer                             :: JLB, JUB     
!        integer                             :: KLB, KUB    
!
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size3D%ILB     
!        IUB = Me%Size3D%IUB     
!
!        JLB = Me%Size3D%JLB    
!        JUB = Me%Size3D%JUB    
!
!        KLB = Me%Size3D%KLB   
!        KUB = Me%Size3D%KUB   
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do k = KLB, KUB
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
!                Index         = Index + 1
!                Vector(Index) = Matrix3D(i,j,k)
!            endif    
!        enddo
!        enddo
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix3D_R - ModuleInterface - ERR01'
        
        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalMatrix3D)
        LocalMatrix3D => Matrix3D
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            k = Me%Index2K(Index)
            Vector(Index) = LocalMatrix3D(i,j,k)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop
            
        !----------------------------------------------------------------------

    end subroutine UnfoldMatrix3D_R8

    !--------------------------------------------------------------------------

    subroutine UnfoldMatrix3D_I (Matrix3D, Vector)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:,:), pointer      :: Matrix3D
        integer, dimension(:    ), pointer      :: Vector

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j, k
        integer                                 :: NLB, NUB
        !$ integer                              :: CHUNK
        integer, dimension(:,:,:), pointer      :: LocalMatrix3D
        
!        integer                                 :: ILB, IUB      
!        integer                                 :: JLB, JUB     
!        integer                                 :: KLB, KUB
!        logical                                 :: Vertical1D_Aux
!
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size3D%ILB     
!        IUB = Me%Size3D%IUB     
!
!        JLB = Me%Size3D%JLB    
!        JUB = Me%Size3D%JUB    
!
!        KLB = Me%Size3D%KLB   
!        KUB = Me%Size3D%KUB   
!
!        if (present(Vertical1D)) then
!            Vertical1D_Aux = Vertical1D
!        else
!            Vertical1D_Aux = .false.
!        endif
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do k = KLB, KUB
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints3D(i,j,k)==1) then
!                Index         = Index + 1
!                Vector(Index) = Matrix3D(i,j,k)
!                if (Vertical1D_Aux .and. .not.(i==2 .and. j==2)) Vector(Index) = 0
!            endif    
!        enddo
!        enddo
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix3D_I - ModuleInterface - ERR01'
            
        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,k,LocalMatrix3D)
        LocalMatrix3D => Matrix3D
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            k = Me%Index2K(Index)
            Vector(Index) = LocalMatrix3D(i,j,k)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop
        
    end subroutine UnfoldMatrix3D_I

   !----------------------------------------------------------------------

    subroutine UnfoldMatrix2D_R (Matrix2D, Vector)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:,:), pointer       :: Matrix2D
        real(4), dimension(:  ), pointer       :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i, j
        integer                             :: NLB, NUB
        !$ integer                          :: CHUNK
        real(4), dimension(:,:), pointer       :: LocalMatrix2D
!        integer                             :: ILB, IUB      
!        integer                             :: JLB, JUB     
!          
!           
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size2D%ILB     
!        IUB = Me%Size2D%IUB     
!
!        JLB = Me%Size2D%JLB    
!        JUB = Me%Size2D%JUB    
!
!          
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
!
!                    Index         = Index + 1
!                    
!             if (Me%ExternalVar%OpenPoints2D(i,j)==1) then
!                
!                Vector(Index) = Matrix2D(i,j)
!             
!             endif
!                    
!                    
!            endif    
!        enddo
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix2D_R - ModuleInterface - ERR01'
            

        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,LocalMatrix2D)
        LocalMatrix2D => Matrix2D
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            Vector(Index) = LocalMatrix2D(i,j)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop

        
    end subroutine UnfoldMatrix2D_R
   !----------------------------------------------------------------------

    subroutine UnfoldMatrix2D_R8 (Matrix2D, Vector)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:,:), pointer    :: Matrix2D
        real(8), dimension(:  ), pointer    :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i, j
        integer                             :: NLB, NUB
        !$ integer                          :: CHUNK
        real(8), dimension(:,:), pointer    :: LocalMatrix2D
!        integer                             :: ILB, IUB      
!        integer                             :: JLB, JUB     
!          
!           
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size2D%ILB     
!        IUB = Me%Size2D%IUB     
!
!        JLB = Me%Size2D%JLB    
!        JUB = Me%Size2D%JUB    
!
!          
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
!
!                    Index         = Index + 1
!                    
!             if (Me%ExternalVar%OpenPoints2D(i,j)==1) then
!                
!                Vector(Index) = Matrix2D(i,j)
!             
!             endif
!                    
!                    
!            endif    
!        enddo
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix2D_R - ModuleInterface - ERR01'
            

        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,LocalMatrix2D)
        LocalMatrix2D => Matrix2D
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            Vector(Index) = LocalMatrix2D(i,j)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop
       
    end subroutine UnfoldMatrix2D_R8
            
    !----------------------------------------------------------------------

    subroutine UnfoldMatrix2D_I (Matrix2D, Vector)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:), pointer        :: Matrix2D
        integer, dimension(:  ), pointer        :: Vector

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i, j
        integer                                 :: NLB, NUB
        !$ integer                              :: CHUNK
        integer, dimension(:,:), pointer        :: LocalMatrix2D
!        integer                                 :: ILB, IUB      
!        integer                                 :: JLB, JUB     
!
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size2D%ILB     
!        IUB = Me%Size2D%IUB     
!
!        JLB = Me%Size2D%JLB    
!        JUB = Me%Size2D%JUB    
!        
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do j = JLB, JUB
!        do i = ILB, IUB
!            if (Me%ExternalVar%WaterPoints2D(i,j)==1) then
!             Index         = Index + 1
!             if (Me%ExternalVar%OpenPoints2D(i,j)==1) then
!                
!                Vector(Index) = Matrix2D(i,j)
!             
!             endif
!            endif    
!        enddo
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix2D_I - ModuleInterface - ERR01'
!            

        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$ CHUNK = CHUNK_I(NLB, NUB)
        !$OMP PARALLEL PRIVATE(Index,i,j,LocalMatrix2D)
        LocalMatrix2D => Matrix2D
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            j = Me%Index2J(Index)
            Vector(Index) = LocalMatrix2D(i,j)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop
        
   end subroutine UnfoldMatrix2D_I

   !----------------------------------------------------------------------

    subroutine UnfoldMatrix1D_R (Matrix1D, Vector)

        !Arguments-------------------------------------------------------------
        real, dimension(:  ), pointer       :: Matrix1D
        real, dimension(:  ), pointer       :: Vector

        !Local-----------------------------------------------------------------
        integer                             :: Index
        integer                             :: i
        integer                             :: NLB, NUB
        !$ integer                          :: CHUNK
        real, dimension(:  ), pointer       :: LocalMatrix1D        
!        integer                             :: ILB, IUB      
!           
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size1D%ILB     
!        IUB = Me%Size1D%IUB     
!
!        !Number indexed to 1D cell in the vector 
!        Index = 0
!
!        do i = ILB, IUB
!            if (Me%ExternalVar%RiverPoints1D(i)==1) then
!
!                    Index         = Index + 1
!                    
!             if (Me%ExternalVar%OpenPoints1D(i)==1) then
!                
!                Vector(Index) = Matrix1D(i)
!             
!             endif
!                    
!                    
!            endif    
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix1D_R - ModuleInterface - ERR01'
            
        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$OMP PARALLEL PRIVATE(Index,i,LocalMatrix1D)
        LocalMatrix1D => Matrix1D
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            Vector(Index) = LocalMatrix1D(i)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop

    end subroutine UnfoldMatrix1D_R
        
    !----------------------------------------------------------------------

    subroutine UnfoldMatrix1D_I (Matrix1D, Vector)

        !Arguments-------------------------------------------------------------
        integer, dimension(:  ), pointer        :: Matrix1D
        integer, dimension(:  ), pointer        :: Vector

        !Local-----------------------------------------------------------------
        integer                                 :: Index
        integer                                 :: i
        integer                                 :: NLB, NUB
        !$ integer                              :: CHUNK
        integer, dimension(:  ), pointer           :: LocalMatrix1D        
!        integer                                 :: ILB, IUB      
!        !----------------------------------------------------------------------
!
!        ILB = Me%Size1D%ILB     
!        IUB = Me%Size1D%IUB     
!
!        !Number indexed to 3D cell in the vector 
!        Index = 0
!
!        do i = ILB, IUB
!            if (Me%ExternalVar%RiverPoints1D(i)==1) then
!             Index         = Index + 1
!             if (Me%ExternalVar%OpenPoints1D(i)==1) then
!                
!                Vector(Index) = Matrix1D(i)
!             
!             endif
!            endif    
!        enddo
!
!        if ((Index) > Me%Array%IUB) stop 'UnfoldMatrix1D_I - ModuleInterface - ERR01'
!            
        !griflet: new way
        !griflet: start
        NLB = Me%Array%ILB
        NUB = Me%Array%IUB
        !$OMP PARALLEL PRIVATE(Index,i,LocalMatrix1D)
        LocalMatrix1D => Matrix1D
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do Index = NLB, NUB
            i = Me%Index2I(Index)
            Vector(Index) = LocalMatrix1D(i)
        enddo
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        !griflet: stop

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
        integer                         :: numAutotrophicPop, numHeterotrophicPop, numAnaerobicPop
        integer                         :: numSolPop

        
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
                                    AutotrophicPop                  = numAutotrophicPop,    &
                                    HeterotrophicPop                = numHeterotrophicPop,  &
                                    AnaerobicPop                    = numAnaerobicPop,      &
                                    SolPop                          = numSolPop,            &
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

                else if (PropertyID== HeterotrophicPop_    ) then
                    nProperty = numHeterotrophicPop

                else if (PropertyID== AutotrophicPop_      ) then
                    nProperty = numAutotrophicPop

                else if (PropertyID== AnaerobicPop_        ) then
                    nProperty = numAnaerobicPop

                else if (PropertyID== SolPop_              ) then
                    nProperty = numSolPop

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
                
                
                
               case(BenthicEcologyModel)
                
                call GetBenthicEcologyPropIndex(Me%ObjBenthicEcology, PropertyID, nProperty, STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
                
                select case(PropertyID)
                        
                        !if it's not a property it can be a rate ID number 
                        case(NintFactor_, PintFactor_, NintFactorR_, RootsMort_, PintFactorR_ ) 

                            nProperty = PropertyID 

                        case default 
                    
                        write(*,*) 
                        write(*,*) 'Inconsistency between Interface and Benthic ecology Model.'
                        stop       'PropertyIndexNumber - ModuleInterface - ERR22.1'
                    
               end select

                end if
                
                
                
                 case(SeagrassSedimInteractionModel)
                
                call GetSeagrassSedimInteractionPropIndex(Me%ObjSeagrassSedimInteraction, PropertyID, nProperty, STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 

                    select case(PropertyID)
                        
                        !if it's not a property it can be a rate ID number 
                        case(RootsUptakeN_) 

                            nProperty = PropertyID 
                            
                        case(RootsUptakeP_) 

                            nProperty = PropertyID 

                        case default 

                            write(*,*) 
                            write(*,*) 'Inconsistency between Interface and SeagrassSedimInteraction Model.'
                            stop       'PropertyIndexNumber - ModuleInterface - ERR27'

                    end select      
                    
                    
                    
                end if
                
                
         case(SeagrassWaterInteractionModel)
                   
                call GetSeagrassWaterInteractionPropIndex(Me%ObjSeagrassWaterInteraction, PropertyID, nProperty, STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
               
                        select case(PropertyID)
                       !if it's not a property it can be a rate ID number 
                        case(LeavesUptakeN_, LeavesUptakeP_, LeavesLightFactor_) 
                        
                        nProperty = PropertyID
                        
                        case default
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and SeagrassWaterInteraction Model.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR26'
                        end select
                    


                end if
                
#ifdef _PHREEQC_
            !ToDo: Need to change
            case (PhreeqCModel)
                
                select case (PropertyID)
                case (Temperature_, pH_, pE_, SoilDryDensity_)
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
            case (WWTPQModel)
        
                call GetWWTPQPropIndex(Me%ObjWWTPQ,                                  &
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
                if (STAT_CALL .NE. SUCCESS_) stop 'PropertyIndexNumber - ModuleInterface - ERR24'
                      
               if      (PropertyID== Phytoplankton_       ) then
                    nProperty = numPhyto

                else if (PropertyID== Zooplankton_         ) then
                    nProperty = numZoo
 
                else if (PropertyID== Larvae_              ) then 
                    nProperty = numLarvae

                else if (PropertyID== POP_                 ) then 
                    nProperty = numPartOrganicPhosphorus

                else if (PropertyID== DOPRefractory_       ) then 
                    nProperty = numDOPRefractory

                else if (PropertyID== DOPNon_Refractory_   ) then 
                    nProperty = numDOPNonRefractory

                else if (PropertyID== Inorganic_Phosphorus_) then 
                    nProperty = numInorganicPhosphorus

                else if (PropertyID== PON_                 ) then 
                    nProperty = numPartOrganicNitrogen

                else if (PropertyID== PONRefractory_       ) then 
                    nProperty = numPartOrganicNitrogenRef

                else if (PropertyID== DONRefractory_       ) then 
                    nProperty = numDONRefractory

                else if (PropertyID== DONNon_Refractory_   ) then 
                    nProperty = numDONNonRefractory 

                else if (PropertyID== Ammonia_             ) then 
                    nProperty = numAmmonia

                else if (PropertyID== Nitrate_             ) then 
                    nProperty = numNitrate

                else if (PropertyID== Nitrite_             ) then 
                    nProperty = numNitrite
                    
                else if (PropertyID== PON1_                ) then 
                    nProperty = numPONitrogen1
                
                else if (PropertyID== PON2_                ) then 
                    nProperty = numPONitrogen2        
                
                else if (PropertyID== PON3_                ) then 
                    nProperty = numPONitrogen3    
                
                else if (PropertyID== PON4_                ) then 
                    nProperty = numPONitrogen4
                
                else if (PropertyID== PON5_                ) then 
                    nProperty = numPONitrogen5
                    
                else if (PropertyID== POP1_                ) then 
                    nProperty = numPOPhosphorus1
                
                else if (PropertyID== POP2_                ) then 
                    nProperty = numPOPhosphorus2        
                
                else if (PropertyID== POP3_                ) then 
                    nProperty = numPOPhosphorus3    
                
                else if (PropertyID== POP4_                ) then 
                    nProperty = numPOPhosphorus4
                
                else if (PropertyID== POP5_                ) then 
                    nProperty = numPOPhosphorus5    
                    
                else if (PropertyID== Bacteria_            ) then 
                    nProperty = numBacteria

                else if (PropertyID== Ciliate_             ) then 
                    nProperty = numCiliate

                else if (PropertyID== BOD_                 ) then 
                    nProperty = numBOD

                else if (PropertyID== Oxygen_              ) then 
                    nProperty = numOxygen
                
                else if (PropertyID== Diatoms_            ) then 
                    nProperty = numDiatoms 

                else if (PropertyID== DSilica_            ) then 
                    nProperty = numSiDiss 

                else if (PropertyID== BioSilica_           ) then 
                    nProperty = numSiBio

                else if (PropertyID== GrossProd_           ) then  
                    CheckName = CheckPropertyName('grossprod', number = nProperty)
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR33'

                else if (PropertyID== NutrientLim_         ) then 
                    CheckName = CheckPropertyName('nutrientlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR34'

                else if (PropertyID== NLim_                ) then  
                    CheckName = CheckPropertyName('nitrogenlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR35'

                else if (PropertyID== PLim_                ) then  
                    CheckName = CheckPropertyName('phosphoruslim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR36'

                else if (PropertyID== LightLim_            ) then   
                    CheckName = CheckPropertyName('lightlim', number = nProperty) 
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR37'

                else if (PropertyID== TemperatureLim_      ) then
                    CheckName = CheckPropertyName('temperaturelim', number = nProperty)    
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR38'

                else if (PropertyID== DiaGrossProd_           ) then
                    CheckName = CheckPropertyName('diagrossprod', number = nProperty)
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR39'

                else if (PropertyID== DiaNutrientLim_         ) then 
                    CheckName = CheckPropertyName('dianutrientlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR40'

                else if (PropertyID== DiaNLim_         ) then   
                    CheckName = CheckPropertyName('dianitrogenlim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR41'

                else if (PropertyID== DiaPLim_         ) then 
                    CheckName = CheckPropertyName('diaphosphoruslim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR42'

                else if (PropertyID== DiaSiLim_         ) then  
                    CheckName = CheckPropertyName('diasilicalim', number = nProperty)     
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR43'

                else if (PropertyID== DiaLightLim_            ) then 
                    CheckName = CheckPropertyName('dialightlim', number = nProperty) 
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR44'

                else if (PropertyID== DiaTemperatureLim_      ) then
                    CheckName = CheckPropertyName('diatemperaturelim', number = nProperty)    
                    if (.NOT.CheckName) stop 'PropertyIndexNumber - ModuleInterface - ERR45'

                else
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and WWTPQ.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR46'
                end if
                
        case(BivalveModel)
         
                call GetBivalvePropIndex(Me%ObjBivalve,PropertyID,nProperty,STAT_CALL)
                
                if (STAT_CALL.NE.SUCCESS_) then 
                    
                    write(*,*) 
                    write(*,*) 'Inconsistency between Interface and Bivalve.'
                    stop       'PropertyIndexNumber - ModuleInterface - ERR47'

                end if  
                
        case default

            write(*,*) 
            write(*,*) 'Defined sinks and sources model was not recognised.'
            stop 'ReadInterfaceFilesName - Module Interface - ERR48' 

        end select

        PropertyIndexNumber = nProperty

        !----------------------------------------------------------------------

    end Function PropertyIndexNumber

    !--------------------------------------------------------------------------

    real function InterfaceDT()
    
        integer         :: STAT_CALL
        
        select case (Me%SinksSourcesModel)

            case(WaterQualityModel)                        
                call GetDTWQM(Me%ObjWaterQuality, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR09'
                
            case(SedimentQualityModel)                        
                call GetDTSedimentQuality(Me%ObjSedimentQuality, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR10'                

            case(CEQUALW2Model, BenthicCEQUALW2Model)
                call GetDTCEQUALW2(Me%ObjCEQUALW2, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR11'
                   
            case(LifeModel)
                call GetDTLife(Me%ObjLife, DT = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR12'
#ifdef _BFM_  
            case(BFMModel)
                call GetDTBFM(Me%ObjBFM, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR12a'
#endif  
            case(BenthosModel)
                call GetDTBenthos(Me%ObjBenthos, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR13'

            case(MacroAlgaeModel)
                call GetDTMacroAlgae(Me%ObjMacroAlgae, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR14'    
                
            case(BenthicEcologyModel)
                call GetDTBenthicEcology(Me%ObjBenthicEcology, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR15' 
                 
                                  
#ifdef _PHREEQC_
            case(PhreeqcModel)                         
                call GetPhreeqCDT(Me%ObjPhreeqC, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR16'                
#endif                   
            case(WWTPQModel)                        
                call GetDTWWTPQM(Me%ObjWWTPQ, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR17'
           
           
            case(SeagrassSedimInteractionModel)
                call GetDTSeagrassSedimInteraction(Me%ObjSeagrassSedimInteraction, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR18' 
           
            case(SeagrassWaterInteractionModel)
                call GetDTSeagrassWaterInteraction(Me%ObjSeagrassWaterInteraction, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR19' 

            case(BivalveModel)
                call GetDTBivalve(Me%ObjBivalve, DTSecond = InterfaceDT, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterfaceDT - ModuleInterface - ERR20' 

            case default
                write(*,*) 
                write(*,*) 'Defined sinks and sources model was not recognised.'
                stop 'InterfaceDT - ModuleInterface - ERR21'

        end select
                
    end function InterfaceDT

    !--------------------------------------------------------------------------

    character(256) function model_name()
    
        select case (Me%SinksSourcesModel)
        
            case (WaterQualityModel)            
                model_name = 'Water Quality'
                
            case (SedimentQualityModel)            
                model_name = 'Sediment Quality'
                
            case (CEQUALW2Model, BenthicCEQUALW2Model)            
                model_name = 'CEQUALW2'
            
            case (LifeModel)            
                model_name = 'Life'                
#ifdef _BFM_  
            case (BFMModel)            
                model_name = 'BFM'
#endif
            case (BenthosModel)            
                model_name = 'Benthos'
            
            case (MacroAlgaeModel)            
                model_name = 'Macro Algae'
                
            case (BenthicEcologyModel)
                model_name = 'BenthicEcology'
                                
#ifdef _PHREEQC_
            case (PhreeqCModel)             
                model_name = 'PhreeqC'
#endif
            case (WWTPQModel)            
                model_name = 'WWTPQ'
            
            case(SeagrassSedimInteractionModel)
                 model_name = 'SeagrassSedimInteraction'
            
            case(SeagrassWaterInteractionModel)
                 model_name = 'SeagrassWaterInteraction'

            case(BivalveModel)
                 model_name = 'Bivalve'
                 
        end select
        
    end function     
    
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
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR10'
                    
                    case(SedimentQualityModel)

                        call KillSedimentQuality(Me%ObjSedimentQuality, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR20'

                    case(LifeModel)

                        call KillLife(Me%ObjLife, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR30'
#ifdef _BFM_  
                    case(BFMModel)

                        call KillBFM(Me%ObjBFM, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR40'
#endif  
                    case(CEQUALW2Model, BenthicCEQUALW2Model)

                        call KillCEQUALW2(Me%ObjCEQUALW2, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR50'

                    case(BenthosModel)
                        
                        call KillBenthos(Me%ObjBenthos, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR60'
                    
                    case(MacroAlgaeModel)
                        
                        call KillMacroAlgae(Me%ObjMacroAlgae, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR70'
                        
                        deallocate(Me%MacrOccupation, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR70.1' 
                        
                     case(BenthicEcologyModel)
                        
                       call KillBenthicEcology(Me%ObjBenthicEcology, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.7'
                        
                       deallocate(Me%WaterMassInKgIncrement, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.8' 
                        
                       deallocate(Me%Sediment, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.10'
                       
                       deallocate(Me%WaterVolume1D, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.11'
                       
                       deallocate(Me%CellArea1D, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.12'
                       
                       deallocate(Me%MassinKgFromWater, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.13'
                       
                       deallocate(Me%BottomSWRadiationAverage, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.14'
                       
                       deallocate(Me%ShearStress2D, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.15'
                       
                       deallocate(Me%UptakeNH4NO3w, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.16'
                       
                       deallocate(Me%UptakePO4w, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.17'      

                       deallocate(Me%UptakeNH4s, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.18'  
                       
                       deallocate(Me%UptakePO4s, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.19'  
 
                       deallocate(Me%LightFactor, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.20'      
                        

#ifdef _PHREEQC_
                    case(PhreeqCModel)
                    
                        if(associated(Me%WaterVolume1D))then
                            deallocate(Me%WaterVolume1D, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR80'
                        end if
                    
                        call KillPhreeqC (Me%ObjPhreeqC, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR100'
#endif
                    case(WWTPQModel)

                        call KillWWTPQ(Me%ObjWWTPQ, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR110'
                   
                   case(SeagrassSedimInteractionModel)  ! Isabella
                        
                   call KillSeagrassSedimInteraction(Me%ObjSeagrassSedimInteraction, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.21'
                        
                    deallocate(Me%SedimCellVol, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.22'
                              
                    deallocate(Me%NintFactorR, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.23'
                       
                   deallocate(Me%RootsMort, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.24'
                   
                   deallocate(Me%PintFactorR, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.25'             
                   
                   case(SeagrassWaterInteractionModel)  ! Isabella
                        
                    call KillSeagrassWaterInteraction(Me%ObjSeagrassWaterInteraction, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.26'
                        
                    deallocate(Me%WaterVolume1D, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.27'
                              
                    deallocate(Me%PintFactor, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.28'                        
                        
                    deallocate(Me%NintFactor, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.29'
                                           
                   deallocate(Me%SeagOccupation, STAT = STAT_CALL)
                       if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR04.30'             

                   
                   case(BivalveModel)
                   
                        if(associated(Me%WaterVolume1D))then
                            deallocate(Me%WaterVolume1D, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR80'
                        end if

                        if(associated(Me%CellArea1D))then
                            deallocate(Me%CellArea1D, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR81'
                        end if
                        
                        if(associated(Me%VelocityModulus1D))then
                            deallocate(Me%VelocityModulus1D, STAT = STAT_CALL)
                            if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR82a'
                        end if

                            
                        call KillBivalve(Me%ObjBivalve, STAT = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'KillInterface - ModuleInterface - ERR82'

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

                if(associated(Me%WaterMass))then
                    deallocate(Me%WaterMass, STAT = STAT_CALL)
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
                
                !griflet: optimization performance
                !griflet: start
                if (associated(Me%I2Index))     deallocate(Me%I2Index)
                if (associated(Me%IJ2Index))    deallocate(Me%IJ2Index)
                if (associated(Me%IJK2Index))   deallocate(Me%IJK2Index)
                deallocate(Me%Index2I)
                deallocate(Me%Index2J)
                deallocate(Me%Index2K)
                !griflet: end

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

