!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Sand
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jan 2004
! REVISION      : Paulo Leitão - v4.0; Miguel Carmo (04/2005)
! DESCRIPTION   : Module to compute the non-cohesive sediment (sand) transport
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

Module ModuleSand

    use ModuleGlobalData
    use ModuleStopWatch,        only : StartWatch, StopWatch
    use ModuleFunctions         
    use ModuleTime              
    use ModuleHDF5,             only : ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,         &
                                       HDF5WriteData, HDF5FlushMemory, HDF5ReadData, KillHDF5
    use ModuleEnterData           
    use ModuleFillMatrix,       only : ConstructFillMatrix, GetDefaultValue, KillFillMatrix
    use ModuleGridData,         only : ConstructGridData, GetGridData, ModifyGridData,          &
                                       GetGridData2DReference, UngetGridData, KillGridData          
    use ModuleDischarges,       only : GetDischargesNumber, GetDischargesGridLocalization,      &
                                       GetDischargeWaterFlow, GetDischargeConcentration
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,           &
                                       GetTimeSerieLocation, CorrectsCellsTimeSerie,            &
                                       GetNumberOfTimeSeries, TryIgnoreTimeSerie       
    use ModuleHorizontalMap,    only : GetWaterPoints2D, GetBoundaries, GetOpenPoints2D,        &
                                       GetComputeFaces2D, UnGetHorizontalMap
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, WriteHorizontalGrid,                  &
                                       GetHorizontalGridSize, UnGetHorizontalGrid, GetXYCellZ
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes, BoxDif,         &
                                       UngetBoxDif, KillBoxDif
#ifndef _WAVES_
    use ModuleWaves,            only : GetWaves, UnGetWaves
#endif

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartSand
    private ::      AllocateInstance
    private ::      ConstructEvolution
    private ::      ConstructGlobalParameters
    private ::      ConstructClasses
    private ::      StartOutputBoxFluxes
    private ::      Open_HDF5_OutPut_File
    private ::      ConstructTimeSerie
    private ::      Read_Sand_Files_Name
    private ::      ReadInitialField


    !Selector
    public  :: GetSandDensity
    public  :: GetSandDiameters
    public  :: UnGetSand
    private ::      ReadLockExternalVar
    private ::      ReadUnLockExternalVar                     
    
    !Modifier
    public  :: ModifySand
    private ::      ComputeFluxes
    private ::          MeyerPeterTransport
    private ::          AckersTransport
    private ::          VanRijn1Transport !suspended load compute by algebric formula
    private ::          VanRijn2Transport !suspended load compute by numerical integration
    private ::          BailardTransport
    private ::          DibajniaTransport
    private ::          BijkerTransport
    private ::      ComputeEvolution
    private ::          ComputeDischarges   
    private ::      ComputeSmoothSlope
    private ::      BoundaryCondition
    private ::      ComputeResidualEvolution
    private ::      OutPutSandHDF
    private ::      OutPut_TimeSeries
    private ::      ComputeTauCritic
    private ::      OutputBoxFluxes
                
    !Destructor
    public  :: KillSand                                                     
    private ::      DeAllocateInstance
    private ::      WriteFinalState

    !Management
    private ::      Ready
    private ::          LocateObjSand 

    !Functions-----------------------------------------------------------------
    private :: FallVel

    !Interfaces----------------------------------------------------------------
    private :: UnGetSand2D_I
    private :: UnGetSand2D_R8
    private :: UnGetSand2D_R4
    interface  UnGetSand
        module procedure UnGetSand2D_I
        module procedure UnGetSand2D_R4
        module procedure UnGetSand2D_R8
    end interface  UnGetSand

    !Parameters
    character(LEN = StringLength), parameter    :: class_block_begin     = '<beginclass>'
    character(LEN = StringLength), parameter    :: class_block_end       = '<endclass>'
    character(LEN = StringLength), parameter    :: diam_block_begin      = '<<begindiam>>'
    character(LEN = StringLength), parameter    :: diam_block_end        = '<<enddiam>>'
    character(LEN = StringLength), parameter    :: percent_block_begin   = '<<beginpercent>>'
    character(LEN = StringLength), parameter    :: percent_block_end     = '<<endpercent>>'


    character(LEN = StringLength), parameter    :: D90_block_begin       = '<beginD90>'
    character(LEN = StringLength), parameter    :: D90_block_end         = '<endD90>'
    character(LEN = StringLength), parameter    :: D50_block_begin       = '<beginD50>'
    character(LEN = StringLength), parameter    :: D50_block_end         = '<endD50>'
    character(LEN = StringLength), parameter    :: D35_block_begin       = '<beginD35>'
    character(LEN = StringLength), parameter    :: D35_block_end         = '<endD35>'
    character(LEN = StringLength), parameter    :: Rock_block_begin      = '<beginrock>'
    character(LEN = StringLength), parameter    :: Rock_block_end        = '<endrock>'

    integer, parameter :: NoTransport = 0, Ackers = 1, MeyerPeter = 2, VanRijn1 = 3, & 
                          VanRijn2 = 4, Bailard = 5, Dibajnia = 6, Bijker = 7
    integer, parameter :: NullGradient = 1, Cyclic = 2

    !Selma
    integer, parameter :: Time_ = 1

    !Types---------------------------------------------------------------------

    private :: T_External
    type       T_External
        type(T_Time)                            :: Now
        real,    pointer, dimension(:,:)        :: DUX              => null()
        real,    pointer, dimension(:,:)        :: DVY              => null()
        real,    pointer, dimension(:,:)        :: DXX              => null()
        real,    pointer, dimension(:,:)        :: DYY              => null()
        real,    pointer, dimension(:,:)        :: DZX              => null()
        real,    pointer, dimension(:,:)        :: DZY              => null()
        integer, pointer, dimension(:,:)        :: ComputeFacesU2D  => null()
        integer, pointer, dimension(:,:)        :: ComputeFacesV2D  => null()
        integer, pointer, dimension(:,:)        :: OpenPoints2D     => null()
        integer, pointer, dimension(:,:)        :: WaterPoints2D    => null()
        integer, pointer, dimension(:,:)        :: BoundaryPoints2D => null()
        real                                    :: WaterDensity     = FillValueReal 
        logical                                 :: WaveTensionON    = .false. 
        real,    pointer, dimension(:,:)        :: Bathymetry       => null()
        real,    pointer, dimension(:,:)        :: InitialBathym    => null()
        real,    pointer, dimension(:,:)        :: WaveDirection    => null()
        real,    pointer, dimension(:,:)        :: Abw              => null()
        real,    pointer, dimension(:,:)        :: Ubw              => null()
        real,    pointer, dimension(:,:)        :: TauTotal         => null()
        real,    pointer, dimension(:,:)        :: CurrentRugosity  => null()
        real,    pointer, dimension(:,:)        :: WaveRugosity     => null()
        real,    pointer, dimension(:,:)        :: WaterColumn      => null()
        real,    pointer, dimension(:,:)        :: VelU             => null()
        real,    pointer, dimension(:,:)        :: VelV             => null()
        real,    pointer, dimension(:,:)        :: VelMod           => null()
        real,    pointer, dimension(:,:)        :: TauWave          => null()
        real,    pointer, dimension(:,:)        :: TauCurrent       => null()
        real,    pointer, dimension(:,:)        :: ShearVelocity    => null()    
        real,    pointer, dimension(:,:)        :: WaveHeight       => null()
        real,    pointer, dimension(:,:)        :: WavePeriod       => null()
        real                                    :: MinWaterColumn   = FillValueReal 
    end type T_External

    
    private :: T_Aceleration
    type       T_Aceleration
        logical         :: Yes      = .false.
        real            :: Coef     = FillValueReal
    end type   T_Aceleration

    private :: T_OutPut
    type       T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime      => null ()
         integer                                :: NextOutPut   = FillValueInt
         logical                                :: Yes          = .false.
    end type T_OutPut

!    private :: T_SubModel
!    type       T_SubModel
!        logical                                 :: ON
!        logical                                 :: Set
!        logical                                 :: InterPolTime = .false.
!        logical                                 :: Initial
!        real,    dimension(:,:), pointer        :: NextField, PreviousField
!        type(T_Time)                            :: NextTime, PreviousTime
!    endtype   

    private :: T_Filter
    type       T_Filter
        logical                                 :: ON       = .false. 
        integer                                 :: Scheme   = FillValueInt
        integer                                 :: Radius   = FillValueInt
        real, pointer, dimension(:,:  )         :: Field2D  => null()
    end type T_Filter

    private :: T_Evolution
    type       T_Evolution
        logical                                 :: Old          = .false. 
        real                                    :: SandDT       = FillValueReal
        real                                    :: BathymDT     = FillValueReal
        type (T_Time)                           :: NextSand, NextBatim
        logical                                 :: Bathym       = .false. 
        !Selma
        integer                                 :: BathymType   = Time_
    end type T_Evolution

    private :: T_Property
    type       T_Property
        type (T_PropertyID)                     :: ID
!        type (T_SubModel  )                     :: SubModel
        real                                    :: Scalar       = FillValueReal
        real, pointer, dimension(:,:  )         :: Field2D      => null ()
    end type T_Property

    private :: T_Classes
    type       T_Classes
        integer                                            :: Number        = FillValueInt
        character(Len=StringLength), dimension(:), pointer :: Name          => null ()
        type (T_Property),           dimension(:), pointer :: Diameter      => null ()
        type (T_Property),           dimension(:), pointer :: Percentage    => null ()
    end type T_Classes 

    private :: T_Files
    type       T_Files
        character(Len = StringLength)           :: ConstructData  = null_str
        character(Len = StringLength)           :: InitialSand    = null_str  
        character(Len = StringLength)           :: OutPutFields   = null_str
        character(Len = StringLength)           :: FinalSand      = null_str
    end type  T_Files

    private :: T_SmoothSlope    
    type       T_SmoothSlope    
        real        :: Critic = FillValueReal
        real        :: Factor = FillValueReal
        logical     :: ON     = .false.   
    end type  T_SmoothSlope    
        

    type       T_Discharges
        type(T_Time)                            :: NextCompute
        real                                    :: DT_Compute = FillValueReal
        logical                                 :: Yes        = .false.
    end type T_Discharges 

    type     T_Boxes
        logical                                 :: Yes      = .false.
        character(Len = StringLength)           :: File     =  null_str
        real(8), dimension(:,:), pointer        :: Mass     => null()
        real(8), dimension(:,:), pointer        :: FluxesX  => null()
        real(8), dimension(:,:), pointer        :: FluxesY  => null()
    end type T_Boxes




    private :: T_Sand
    type       T_Sand
        integer                                    :: InstanceID            = FillValueInt
        type (T_Size2D)                            :: Size, WorkSize
        type (T_Time)                              :: BeginTime, EndTime
        type (T_Evolution )                        :: Evolution
        type (T_Filter    )                        :: Filter
        type (T_SmoothSlope)                       :: SmoothSlope
        type (T_Property)                          :: BedRock, DZ, BatimIncrement, DZ_Residual
        type (T_Aceleration)                       :: Aceleration
        type (T_Property)                          :: D35, D50, D90
        real                                       :: SandMin               = FillValueReal
        real                                       :: Porosity              = FillValueReal
        real                                       :: Density               = FillValueReal
        real                                       :: RelativeDensity       = FillValueReal
        real                                       :: RhoSl                 = FillValueReal
        integer                                    :: Boundary              = FillValueInt
        real                                       :: TransportFactor       = FillValueReal
        real                                       :: TauMax                = FillValueReal
        logical                                    :: TimeSerie             = .false. 
        type (T_Classes)                           :: Classes
        type (T_Files  )                           :: Files
        type (T_OutPut )                           :: OutPut
        type (T_External)                          :: ExternalVar
        type (T_Discharges)                        :: Discharges
        type (T_Boxes     )                        :: Boxes
        real, dimension(:,:), pointer              :: FluxX                 => null ()
        real, dimension(:,:), pointer              :: FluxY                 => null ()
        real, dimension(:,:), pointer              :: TransportCapacity     => null ()
        real, dimension(:,:), pointer              :: TauCritic             => null ()
        real, dimension(:,:), pointer              :: Dast                  => null ()
        integer                                    :: TransportMethod       = FillValueInt
        !Instance of ModuleHDF5        
        integer                                    :: ObjHDF5               = 0
        !Instance of ModuleTimeSerie            
        integer                                    :: ObjTimeSerie          = 0
        !Instance of Module_EnterData           
        integer                                    :: ObjEnterData          = 0
        !Instance of ModuleGridData where the bathymetry is define             
        integer                                    :: ObjBathym             = 0
        !Instance of ModuleHorizontalGrid       
        integer                                    :: ObjHorizontalGrid     = 0
        !Instance of ModuleHorizontalMap        
        integer                                    :: ObjHorizontalMap      = 0 
        !Instance of ModuleTime                 
        integer                                    :: ObjTime               = 0
        !Instance of ModuleDischarges           
        integer                                    :: ObjDischarges         = 0
        !Instance of ModuleBoxDif               
        integer                                    :: ObjBoxDif             = 0             
        !Instance of ModuleWaves
        integer                                    :: ObjWaves              = 0
        !List of Sand Instances
        type(T_Sand), pointer                      :: Next                  => null ()

    end type  T_Sand

    !Global Module Variables
    type (T_Sand), pointer                         :: FirstObjSand          => null ()
    type (T_Sand), pointer                         :: Me                    => null () 


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartSand(ObjSandID,                             &
                         ObjGridDataID,                         &
                         ObjHorizontalGridID,                   &
                         ObjHorizontalMapID,                    &
                         ObjTimeID,                             &
                         ObjWavesID,                            &
                         ObjDischargesID,                       &
                         WaterDensity,                          &
                         WaveTensionON,                         &
                         STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSandID 
        integer                                         :: ObjGridDataID
        integer                                         :: ObjHorizontalGridID
        integer                                         :: ObjHorizontalMapID
        integer                                         :: ObjTimeID
        integer                                         :: ObjWavesID
        integer                                         :: ObjDischargesID
        real                                            :: WaterDensity
        logical                                         :: WaveTensionON
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSand_)) then
            nullify (FirstObjSand)
            call RegisterModule (mSand_) 
        endif

        call Ready(ObjSandID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ObjTime           = AssociateInstance (mTIME_,           ObjTimeID          )
            Me%ObjBathym         = AssociateInstance (mGRIDDATA_,       ObjGridDataID      )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  ObjHorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, ObjHorizontalGridID)

            if(ObjWavesID /= 0)then
                Me%ObjWaves      = AssociateInstance (mWAVES_,          ObjWavesID         )
            end if

            Me%ExternalVar%WaterDensity  = WaterDensity
            Me%ExternalVar%WaveTensionON = WaveTensionON

            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                             &
                                       Size        = Me%Size,                            &
                                       WorkSize    = Me%WorkSize,                        &
                                       STAT        = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'StartSand - ModuleSand - ERR10'

            call ReadLockExternalVar

            call Read_Sand_Files_Name

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSand - ModuleSand - ERR20'

            call ConstructEvolution

            call ConstructOutputTime

            call ConstructClasses

            call ConstructGlobalParameters

            call StartOutputBoxFluxes

            call ComputeTauCritic

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSand - ModuleSand - ERR30'


            if (Me%OutPut%Yes) call Open_HDF5_OutPut_File(Me%Files%OutPutFields)

            if (Me%Discharges%Yes) then
            
                if (ObjDischargesID == 0)  then                                                
                    write(*,*)'You need to define a water discharges in the hydrodynamic input' 
                    stop      'StartSand - Sand - ERR01'
                else
                    Me%ObjDischarges = AssociateInstance (mDISCHARGES_, ObjDischargesID)

                    Me%Discharges%NextCompute = Me%ExternalVar%Now
                endif

            endif


            call ReadUnLockExternalVar

            !Returns ID
            ObjSandID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleSand - StartSand - ERR40' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartSand
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Sand), pointer                         :: NewObjSand
        type (T_Sand), pointer                         :: PreviousObjSand


        !Allocates new instance
        allocate (NewObjSand)
        nullify  (NewObjSand%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSand)) then
            FirstObjSand         => NewObjSand
            Me                    => NewObjSand
        else
            PreviousObjSand      => FirstObjSand
            Me                    => FirstObjSand%Next
            do while (associated(Me))
                PreviousObjSand  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjSand
            PreviousObjSand%Next => NewObjSand
        endif

        Me%InstanceID = RegisterNewInstance (mSand_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine StartOutputBoxFluxes

        !External--------------------------------------------------------------
        integer                                             :: iflag, STAT_CALL
        integer                                             :: ILB, IUB, JLB, JUB
        logical                                             :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: ScalarOutputList
        character(len=StringLength), dimension(:), pointer  :: FluxesOutputList
        !----------------------------------------------------------------------


        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB


        !<BeginKeyword>
            !Keyword          : BOXFLUXES
            !<BeginDescription>       
            ! if exist fluxes between boxes are compute 
            !
            !<EndDescription>
            !Type             : Logical 
            !Default          : Do not have
            !File keyword     : *******
            !Multiple Options : Do not have
            !Search Type      : From File
        !<EndKeyword>
        
        call GetData(Me%Boxes%Yes,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     keyword      = 'BOXFLUXES',                                        &
                     ClientModule = 'ModuleSand',                                       &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'StartOutputBoxFluxes - ModuleSand - ERR10'

i1:     if (Me%Boxes%Yes) then

            !<BeginKeyword>
                !Keyword          : BOX_FILENAME
                !<BeginDescription>       
                ! Name file where the boxes are defined
                !
                !<EndDescription>
                !Type             : Character 
                !Default          : Do not have
                !File keyword     : *******
                !Multiple Options : Do not have
                !Search Type      : From File
            !<EndKeyword>
        
            call GetData(Me%Boxes%File,                                                 &
                         Me%ObjEnterData, iflag,                                        &
                         keyword      = 'BOX_FILENAME',                                 &
                         ClientModule = 'ModuleSand',                                   &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'StartOutputBoxFluxes - ModuleSand - ERR20'

            if (iflag .EQ. 0) then
                stop 'StartOutputBoxFluxes - ModuleSand - ERR30'    
            endif
        
            inquire(File = Me%Boxes%File, Exist = exist)
            if (exist) then
                inquire(File = Me%Boxes%File, Opened  = Opened)
                if (opened) then
                    write(*,*    ) 
                    write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Boxes%File))
                    write(*,*    ) 'Already opened.'
                    stop           'StartOutputBoxFluxes - ModuleSand - ERR40'    
                end if
            else
                write(*,*) 
                write(*,*)     'Could not find the boxes file.'
                write(*,'(A)') 'BoxFileName = ', Me%Boxes%File
                stop           'StartOutputBoxFluxes - ModuleSand - ERR50'    
            end if

            allocate(ScalarOutputList(1))
            allocate(FluxesOutputList(1))

            ScalarOutputList(1) = 'sand'
            FluxesOutputList(1) = 'sand'

            call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                         &
                             TimeID             = Me%ObjTime,                           &
                             HorizontalGridID   = Me%ObjHorizontalGrid,                 &
                             BoxesFilePath      = Me%Boxes%File,                        &
                             FluxesOutputList   = FluxesOutputList,                     &
                             ScalarOutputList   = ScalarOutputList,                     &
                             WaterPoints2D      = Me%ExternalVar%WaterPoints2D,         &
                             STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'StartOutputBoxFluxes - ModuleSand - ERR60'

            deallocate(FluxesOutputList)
            deallocate(ScalarOutputList)

            allocate(Me%Boxes%FluxesX(ILB:IUB, JLB:JUB))
            Me%Boxes%FluxesX(:,:) = 0.

            allocate(Me%Boxes%FluxesY(ILB:IUB, JLB:JUB))
            Me%Boxes%FluxesY(:,:) = 0.

            allocate(Me%Boxes%Mass   (ILB:IUB, JLB:JUB))
            Me%Boxes%Mass   (:,:) = 0.

        endif i1


    end subroutine StartOutputBoxFluxes


    !--------------------------------------------------------------------------

    subroutine ConstructEvolution

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        real    :: ModelDT
        integer :: STAT_CALL, iflag
        real(8) :: ErrorAux, auxFactor, DTaux



        !Begin-----------------------------------------------------------------

        call GetComputeTimeStep(Me%ObjTime, ModelDT, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                 &
            stop 'ConstructEvolution - ModuleSand - ERR10'



        !<BeginKeyword>
            !Keyword          : SAND_DT
            !<BeginDescription>       
               ! The time step of the SAND evolution
               !
            !<EndDescription>
            !Type             : real 
            !Default          : ModelDT
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Evolution%SandDT,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SAND_DT',                                           &
                     default      = ModelDT,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR20' 

        call GetComputeTimeLimits(Me%ObjTime,                                            &
                                  EndTime   = Me%EndTime,                                &
                                  BeginTime = Me%BeginTime,                              &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ConstructEvolution - ModuleSand - ERR30'

        if (Me%Evolution%SandDT < (ModelDT)) then

            !Sand DT  must be a submultiple of the ModelDT
            auxFactor = ModelDT / Me%Evolution%SandDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop 'ConstructEvolution - ModuleSand - ERR40'
            endif



        elseif (Me%Evolution%SandDT > (ModelDT)) then

            !Sand DT  must be a multiple of the ModelDT
            auxFactor = Me%Evolution%SandDT  / ModelDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop 'ConstructEvolution - ModuleSand - ERR50'
            endif

        endif

        ! Run period in seconds
        DTaux = Me%EndTime - Me%BeginTime

        !The run period   must be a multiple of the SAND DT
        auxFactor = DTaux / Me%Evolution%SandDT

        ErrorAux = auxFactor - int(auxFactor)
        if (ErrorAux /= 0) then
            write(*,*) 
            write(*,*) ' Time step error.'
            stop 'ConstructEvolution - ModuleSand - ERR60'
        endif


        Me%Evolution%NextSand = Me%BeginTime + Me%Evolution%SandDT


        !<BeginKeyword>
            !Keyword          : BATIM_DT
            !<BeginDescription>       
               ! The time step of the BATIM evolution
               !
            !<EndDescription>
            !Type             : real 
            !Default          : ModelDT
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Evolution%BathymDT,                                              &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BATIM_DT',                                          &
                     default      = Me%Evolution%SandDT,                                 &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR70' 

        call GetComputeTimeLimits(Me%ObjTime,                                            &
                                  EndTime   = Me%EndTime,                                &
                                  BeginTime = Me%BeginTime,                              &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ConstructEvolution - ModuleSand - ERR80'



        if (Me%Evolution%BathymDT < Me%Evolution%SandDT) then

            stop 'ConstructEvolution - ModuleSand - ERR90'

        elseif (Me%Evolution%BathymDT > Me%Evolution%SandDT) then

            !Batim DT  must be a multiple of the Sand DT
            auxFactor = Me%Evolution%BathymDT / Me%Evolution%SandDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop 'ConstructEvolution - ModuleSand - ERR100'
            endif

        endif

        ! Run period in seconds
        DTaux = Me%EndTime - Me%BeginTime

        !The run period   must be a multiple of the BATIM DT
        auxFactor = DTaux / Me%Evolution%BathymDT

        ErrorAux = auxFactor - int(auxFactor)
        if (ErrorAux /= 0) then
            write(*,*) 
            write(*,*) ' Time step error.'
            stop 'ConstructEvolution - ModuleSand - ERR120'
        endif

        Me%Evolution%NextBatim = Me%BeginTime + Me%Evolution%BathymDT

        !<BeginKeyword>
            !Keyword          : OLD
            !<BeginDescription>       
               ! Check if the user wants to start from the final condition of a previous run 
               !
            !<EndDescription>
            !Type             : logical
            !Default          : 0.1
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Evolution%OLD,                                                   &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'OLD',                                               &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR130' 



        call GetData(Me%Aceleration%Coef,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'ACELERATION',                                       &
                     default      = 1.,                                                  &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR140' 

        if (iflag==1) then
            Me%Aceleration%Yes = .true.
        else
            Me%Aceleration%Yes = .false.
        endif

    end subroutine ConstructEvolution

    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File(FileName)

        !Arguments-------------------------------------------------------------
        character(Len=*)                            :: FileName

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: HDF5_CREATE, i

        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5      (Me%ObjHDF5,                                &
                                 trim(FileName)//"5",                       &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR01'


        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,         &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR03'


        call GetGridData2Dreference(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR04' 

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = Me%ExternalVar%InitialBathym,                    &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR05'

        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR06' 


        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                 &
                              Array2D = Me%ExternalVar%WaterPoints2D,                    &
                              STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR07'

        call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%D35%ID%Name),   &
                               trim(Me%D35%ID%Units), Array2D = Me%D35%Field2D,          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR08'


        call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%D50%ID%Name),   &
                               trim(Me%D50%ID%Units), Array2D = Me%D50%Field2D,          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR09'

        call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%D90%ID%Name),   &
                               trim(Me%D90%ID%Units), Array2D = Me%D90%Field2D,          &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR10'

        call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%BedRock%ID%Name),&
                               trim(Me%BedRock%ID%Units), Array2D = Me%BedRock%Field2D,   &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR11'

        call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", "Tau Critic",            &
                              "m2/s2", Array2D = Me%TauCritic,                            &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR12'


        do i=1, Me%Classes%Number

            call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics/Classes/"//trim(Me%Classes%Name(i)),&
                                  "Diameter",  'm', Array2D = Me%Classes%Diameter(i)%Field2D,    &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR13'

            call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics/Classes/"//trim(Me%Classes%Name(i)),&
                                  "Percentage",  '%', Array2D = Me%Classes%Percentage(i)%Field2D,    &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR14'


        enddo

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR15'

        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File
   
    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !Read the name of the files need to construct and modify
    !the sand properties 

    subroutine Read_Sand_Files_Name

        !External--------------------------------------------------------------
        integer                       :: STAT_CALL 
        character(len = StringLength) :: Message

        !----------------------------------------------------------------------

        ! ---> ASCII file used to construct new properties
        Message   ='ASCII file used to construct sand instance.'
        Message   = trim(Message)

        call ReadFileName('SAND_DATA', Me%Files%ConstructData,                           &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sand_Files_Name - ModuleSand - ERR01' 


        ! ---> File in HDF format where is written instant fields of Sand properties
        Message   ='Instant fields of sand properties in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SAND_OUT', Me%Files%OutPutFields,                             &
                           Message = Message, TIME_END = Me%EndTime,                     &
                           Extension = 'sandlt', STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sand_Files_Name - ModuleSand - ERR02' 

        ! ---> Sand properties final values in HDF format
        Message   ='Sand properties final values in HDF format.'
        Message   = trim(Message)
        call ReadFileName('SAND_END', Me%Files%FinalSand,                                &
                           Message = Message, TIME_END = Me%EndTime,                     &
                           Extension = 'sandlf', STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sand_Files_Name - ModuleSand - ERR03' 


        ! ---> Sand properties initial values in HDF format
        Message   ='Sand properties initial values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SAND_INI', Me%Files%InitialSand,                              &
                           Message = Message, TIME_END = Me%ExternalVar%Now,             &
                           Extension = 'sanlf',STAT = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_   ) then
            write(*,*)  
            write(*,*) 'Inicial file not found.'
            stop 'Read_Sand_Files_Name - ModuleSand - ERR04' 

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then
            write(*,*)  
            write(*,*) 'Keyword for the inicial file not found in nomfich.dat.'
            write(*,*) 'Read_Sand_Files_Name - ModuleSand - WRN01'
            write(*,*)  

        else if (STAT_CALL .EQ. SUCCESS_              ) then
            continue
        else
            write(*,*) 
            write(*,*) 'Error calling ReadFileName.'
            stop 'Read_Sand_Files_Name - ModuleSand - ERR05' 
        end if cd1  

        !----------------------------------------------------------------------

    end subroutine Read_Sand_Files_Name

   !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !External--------------------------------------------------------------
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer                                             :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber  
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !----------------------------------------------------------------------

        !First checks out how many properties will have time series

        !Allocates PropertyList
        allocate(PropertyList(2), STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSand - ERR10'

        PropertyList(1) = 'DZ'
        PropertyList(1) = 'DZ_Residual'


        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleSand',                                       &
                     Default      = Me%Files%ConstructData,                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructTimeSerie - ModuleSand - ERR20' 


        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            trim(TimeSerieLocationFile),                                &
                            PropertyList, "srsand",                                     &
                            WaterPoints2D = Me%ExternalVar%WaterPoints2D,               & 
                            STAT          = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSand - ERR30'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSand - ERR40'

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR50'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= OUT_OF_BOUNDS_ERR_) then
                    stop 'ConstructTimeSerie - ModuleSand - ERR60'
                endif                            

                if (STAT_CALL == OUT_OF_BOUNDS_ERR_ .or. Id < 0 .or. Jd < 0) then

                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR70'

                    if (IgnoreOK) then
                        cycle
                    else
                        stop 'ConstructTimeSerie - ModuleSand - ERR80'
                    endif

                endif

                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR90'
            endif

        enddo

        
    end subroutine ConstructTimeSerie

   !-------------------------------------------------------------------------
    
    subroutine ConstructOutputTime()

        !Arguments------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: iflag

        !Begin----------------------------------------------------------------

        call GetOutPutTime(Me%ObjEnterData,                                              &
                           CurrentTime = Me%ExternalVar%Now,                             &
                           EndTime     = Me%EndTime,                                     &
                           keyword     = 'OUTPUT_TIME',                                  &
                           SearchType  = FromFile,                                       &
                           OutPutsTime = Me%OutPut%OutTime,                              &
                           OutPutsOn   = Me%OutPut%Yes,                                  &
                           STAT        = STAT_CALL)                                 
                                                                                    
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ConstructOutputTime - ModuleSand - ERR01' 

        Me%OutPut%NextOutPut = 1

        !<BeginKeyword>
            !Keyword          : TIME_SERIE
            !<BeginDescription>       
               ! 
               ! Checks out if the user pretends to write a time serie for this property
               ! 
            !<EndDescription>
            !Type             : Boolean 
            !Default          : .false.
            !File keyword     : SAND_DATA
            !Multiple Options : Do not have
            !Search Type      : FromFile

        !<EndKeyword>
        call GetData(Me%TimeSerie,                                                       &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword        = 'TIME_SERIE',                                      &
                     Default        = .false.,                                           &
                     SearchType     = FromFile,                                          &
                     ClientModule   = 'ModuleInterfaceSedimentWater',                    &
                     STAT           = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ConstructOutputTime - ModuleSand - ERR02' 

        if (Me%TimeSerie) call ConstructTimeSerie

    end subroutine ConstructOutputTime

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a        
    !sand property values in the domain    

    subroutine ConstructSandProperty(NewProperty, ExtractType)

        !Arguments-------------------------------------------------------------
        type(T_Property)                :: NewProperty
        integer                         :: ExtractType

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------

        call ConstructPropertyID(NewProperty%ID, Me%ObjEnterData, ExtractType)

        
        allocate(NewProperty%Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))


        call ConstructFillMatrix  (PropertyID           = NewProperty%ID,              &
                                   EnterDataID          = Me%ObjEnterData,             &
                                   TimeID               = Me%ObjTime,                  &
                                   HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                   ExtractType          = ExtractType,                 &
                                   PointsToFill2D       = Me%ExternalVar%WaterPoints2D,&
                                   Matrix2D             = NewProperty%Field2D,         &
                                   TypeZUV              = TypeZ_,                      &
                                   STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                     &
            stop 'Construct_SandProperty - ModuleSand - ERR02'

        call GetDefaultValue(NewProperty%ID%ObjFillMatrix, NewProperty%Scalar, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                     &
            stop 'Construct_SandProperty - ModuleSand - ERR03'

        if(NewProperty%ID%SolutionFromFile)then

            stop 'Construct_SandProperty - ModuleSand - ERR04'
            
        else

            call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)&
                stop 'Construct_SandProperty - ModuleSand - ERR05'

        end if


    end subroutine ConstructSandProperty

    !--------------------------------------------------------------------------
    subroutine ConstructClasses

        !External----------------------------------------------------------------
        integer                             :: ClientNumber1, ClientNumber2
        integer                             :: STAT_CALL
        logical                             :: BlockFound, BlockInBlockFound

        !Local-------------------------------------------------------------------
        logical, allocatable, dimension(:)  :: ClassOK
        integer                             :: ClassID

        !------------------------------------------------------------------------

        integer                         :: ILB, IUB, JLB, JUB, iflag, i
        integer                         :: WILB, WIUB, WJLB, WJUB

        !----------------------------------------------------------------------
 
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        WILB = Me%WorkSize%ILB
        WIUB = Me%WorkSize%IUB
        WJLB = Me%WorkSize%JLB
        WJUB = Me%WorkSize%JUB


        !<BeginKeyword>
            !Keyword          : CLASSES_NUMBER
            !<BeginDescription>       
               ! The number of sand classes the user wants to define
               !
            !<EndDescription>
            !Type             : integer 
            !Default          : 1
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Classes%Number,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CLASSES_NUMBER',                                    &
                     default      = 0,                                                   &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR01' 


if1:    if (Me%Classes%Number > 0) then

            allocate(Me%Classes%Name      (1:Me%Classes%Number))

            allocate(Me%Classes%Diameter  (1:Me%Classes%Number))

            allocate(Me%Classes%Percentage(1:Me%Classes%Number)) 

            allocate(ClassOK              (1:Me%Classes%Number))

            ClassOK = .false.

do1 :       do i=1, Me%Classes%Number
                call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                                            ClientNumber    = ClientNumber1,             &
                                            block_begin     = class_block_begin,         &
                                            block_end       = class_block_end,           &
                                            BlockFound      = BlockFound,                &
                                            STAT            = STAT_CALL)
cd1 :           if (STAT_CALL .EQ. SUCCESS_) then    
cd2 :               if (BlockFound) then    

                        call GetData(Me%Classes%Name(i),                                 &
                                     Me%ObjEnterData,iflag,                              &
                                     SearchType   = FromBlock,                           &
                                     keyword      = 'CLASS_NAME',                        &
                                     ClientModule = 'ModuleSand',                        &
                                     STAT         = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR01' 

                        if (iflag == 0) stop 'ConstructClasses - ModuleSand - ERR02'


                        call GetData(ClassID,                                            &
                                     Me%ObjEnterData,iflag,                              &
                                     SearchType   = FromBlock,                           &
                                     keyword      = 'CLASS_ID',                          &
                                     ClientModule = 'ModuleSand',                        &
                                     STAT         = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR01' 

                        if (iflag == 0) stop 'ConstructClasses - ModuleSand - ERR02'

                        if (ClassID < 1 .or. ClassID > Me%Classes%Number) stop 'ConstructClasses - ModuleSand - ERR03' 

                        if (ClassOK(ClassID)) stop 'ConstructClasses - ModuleSand - ERR04'

                        ClassOK(ClassID) = .true.

                        call ExtractBlockFromBlock(Me%ObjEnterData,                        &
                                                    ClientNumber      = ClientNumber2,     &
                                                    block_begin       = Diam_block_begin,  &
                                                    block_end         = Diam_block_end,    &
                                                    BlockInBlockFound = BlockInBlockFound, &
                                                    STAT              = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR05' 

                        if (.not. BlockInBlockFound) stop 'ConstructClasses - ModuleSand - ERR06' 
 
                        call ConstructSandProperty(Me%Classes%Diameter(ClassID), FromBlockInBlock)

                        call RewindBlock(Me%ObjEnterData, ClientNumber1, STAT = STAT_CALL)
                        
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR05' 

                        call ExtractBlockFromBlock(Me%ObjEnterData,                         &
                                                    ClientNumber      = ClientNumber2,      &
                                                    block_begin       = Percent_block_begin,&
                                                    block_end         = Percent_block_end,  &
                                                    BlockInBlockFound = BlockInBlockFound,  &
                                                    STAT              = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR07' 


                        if (.not. BlockInBlockFound) stop 'ConstructClasses - ModuleSand - ERR08' 

                        call ConstructSandProperty(Me%Classes%Percentage(ClassID), FromBlockInBlock)

                        call RewindBlock(Me%ObjEnterData, ClientNumber1, STAT = STAT_CALL)
                        
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR05' 

                    else cd2

                        write(*,*)  
                        write(*,*) 'Error calling ExtractBlockFromBlock. '
                        stop       'ConstructClasses - ModuleSand - ERR09'

                    end if cd2

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                    write(*,*)  
                    write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructClasses - ModuleSand - ERR10'
                else cd1
                    stop       'ConstructClasses - ModuleSand - ERR11'
                end if cd1
            end do do1

            deallocate(ClassOK)

            call Block_UnLock(Me%ObjEnterData, ClientNumber1, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR07'


            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSand - ERR07'

        endif if1

    end subroutine ConstructClasses

    !----------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine ConstructGlobalParameters

        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------
        integer                             :: iflag
        character(Len = StringLength)       :: Auxchar

        !------------------------------------------------------------------------

        if (Me%Classes%Number == 0) then

            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber    = ClientNumber,                  &
                                        block_begin     = D35_block_begin,               &
                                        block_end       = D35_block_end,                 &
                                        BlockFound      = BlockFound,                    &
                                        STAT            = STAT_CALL)
  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR01'

            if (BlockFound) then
                call ConstructSandProperty(Me%D35, FromBlock)

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR02'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR07'

            else
                stop 'ConstructGlobalParameters - ModuleSand - ERR02'
            endif

            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber    = ClientNumber,                  &
                                        block_begin     = D50_block_begin,               &
                                        block_end       = D50_block_end,                 &
                                        BlockFound      = BlockFound,                    &
                                        STAT            = STAT_CALL)
  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR03'

            if (BlockFound) then
                call ConstructSandProperty(Me%D50, FromBlock)

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL .NE. SUCCESS_) stop 'CConstructGlobalParameters - ModuleSand - ERR04'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR07'

            else
                stop 'ConstructGlobalParameters - ModuleSand - ERR04'
            endif


            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber    = ClientNumber,                  &
                                        block_begin     = D90_block_begin,               &
                                        block_end       = D90_block_end,                 &
                                        BlockFound      = BlockFound,                    &
                                        STAT            = STAT_CALL)
  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR05'

            if (BlockFound) then
                call ConstructSandProperty(Me%D90, FromBlock)

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR07'

                if (STAT_CALL .NE. SUCCESS_) stop 'CConstructGlobalParameters - ModuleSand - ERR02'
            else
                stop 'ConstructGlobalParameters - ModuleSand - ERR06'
            endif

        
        else

            !Compute D30, D50 e D90 based in the classes definition
            !call ComputeCharacteristicDiam            

        endif

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                     &
                                    ClientNumber    = ClientNumber,                      &
                                    block_begin     = Rock_block_begin,                  &
                                    block_end       = Rock_block_end,                    &
                                    BlockFound      = BlockFound,                        &
                                    STAT            = STAT_CALL)
  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR07'

        if (BlockFound) then
            call ConstructSandProperty(Me%BedRock, FromBlock)

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'CConstructGlobalParameters - ModuleSand - ERR02'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR07'
        else
            stop 'ConstructGlobalParameters - ModuleSand - ERR08'
        endif


        allocate(Me%DZ%Field2D         (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        Me%DZ%ID%Name  = 'DZ'
        Me%DZ%ID%Units = 'm'

        
        Me%DZ%Field2D(:,:) = 0.

        if (Me%Evolution%BathymDT > Me%Evolution%SandDT) then

            allocate(Me%BatimIncrement%Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            Me%BatimIncrement%ID%Name  = 'Batim Increment'
            Me%BatimIncrement%ID%Units = 'm'
        
            Me%BatimIncrement%Field2D(:,:) = 0.
        else

            Me%BatimIncrement%Field2D => Me%DZ%Field2D

        endif


        allocate(Me%DZ_Residual%Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        Me%DZ_Residual%ID%Name  = 'DZ_Residual'
        Me%DZ_Residual%ID%Units = 'm'

        if (Me%Evolution%OLD) then

            call ReadInitialField(FieldName = Me%DZ_Residual%ID%Name, Field2D = Me%DZ_Residual%Field2D)

        else

            Me%DZ_Residual%Field2D(:,:) = 0.

        endif

        !<BeginKeyword>
            !Keyword          : SAND_MIN
            !<BeginDescription>       
               ! The minimum sand thickness
               !
            !<EndDescription>
            !Type             : real
            !Default          : 0.01
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%SandMin,                                                         &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SAND_MIN',                                          &
                     default      = 0.01,                                                &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR11' 


        !<BeginKeyword>
            !Keyword          : POROSITY
            !<BeginDescription>       
               ! The minimum sand thickness
               !
            !<EndDescription>
            !Type             : real
            !Default          : 0.1
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Porosity,                                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'POROSITY',                                          &
                     default      = 0.1,                                                 &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR12' 


        !<BeginKeyword>
            !Keyword          : DENS_SAND
            !<BeginDescription>       
               ! Sand density in kg/m^3
               !
            !<EndDescription>
            !Type             : real
            !Default          : 2650. 
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Density,                                                         &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DENS_SAND',                                         &
                     default      = 2650.,                                               &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR13' 
        
        Me%RhoSl           = Me%Density - Me%ExternalVar%WaterDensity

        Me%RelativeDensity = Me%RhoSl / Me%ExternalVar%WaterDensity

        !Allocate fluxes
        allocate(Me%FluxX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%FluxX(:,:) = 0.

        allocate(Me%FluxY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%FluxY(:,:) = 0.

        allocate(Me%TransportCapacity(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%TransportCapacity(:,:) = 0.

        allocate(Me%TauCritic(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%TauCritic(:,:) = 0.

        !<BeginKeyword>
            !Keyword          : TRANSPORT_METHOD
            !<BeginDescription>       
               ! Methodology use to compute the sand transport
               !
            !<EndDescription>
            !Type             : character
            !Default          : MeyerPeter
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Auxchar,                                                            &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TRANSPORT_METHOD',                                  &
                     default      = 'MeyerPeter',                                        &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR14' 

        SELECT CASE (trim(Auxchar))

        Case ("no transport")
            Me%TransportMethod = NoTransport
        Case ("MeyerPeter")
            Me%TransportMethod = MeyerPeter
        Case ("Ackers")
            Me%TransportMethod = Ackers
        Case ("VanRijn1")
            Me%TransportMethod = VanRijn1
        Case ("VanRijn2")
            Me%TransportMethod = VanRijn2
        Case ("Bailard")
            Me%TransportMethod = Bailard
        Case ("Dibajnia")
            Me%TransportMethod = Dibajnia
        Case ("Bijker")
            Me%TransportMethod = Bijker

        Case default
            stop 'ConstructGlobalParameters - ModuleSand - ERR15' 
                
        END SELECT 

        if     ( Me%TransportMethod == VanRijn1  &
            .OR. Me%TransportMethod == VanRijn2  &
            .OR. Me%TransportMethod == Bijker )  then
    !.OR. Me%TransportMethod == VanRijn2 .OR. Me%TransportMethod == Bijker
     
            if (Me%ObjWaves == 0) stop 'ConstructGlobalParameters - ModuleSand - ERR16' 

            if (.not. Me%ExternalVar%WaveTensionON) stop 'ConstructGlobalParameters - ModuleSand - ERR20' 

            allocate (Me%Dast(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            Me%Dast(:,:) = FillValueReal

        endif

        !<BeginKeyword>
            !Keyword          : BATHYM_EVOLUTION
            !<BeginDescription>       
               ! Check if the user wants to let the bathymetry evolve due to sand transport
               !
            !<EndDescription>
            !Type             : logical
            !Default          : .true.
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Evolution%Bathym,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BATHYM_EVOLUTION',                                  &
                     default      = .true.,                                              &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR26' 


        call GetData(AuxChar,                                                            &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'FILTER_SCHEME',                                     &
                     ClientModule = 'ModuleSand',                                        &
                     Default      = 'No Filter',                                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR27' 

        Me%Filter%ON = .true.

        Select Case (trim(AuxChar))

        Case ("NO FILTER", "no fielter", "No filter", "No Filter")
            Me%Filter%Scheme = NoFilter
            Me%Filter%ON     = .false.
        Case ("MODIFY LAX", "modify lax", "Modify Lax", "Modify lax")
            Me%Filter%Scheme = ModifyLax
        Case default
            stop 'ConstructGlobalParameters - ModuleSand - ERR28' 
        end select 


        if (Me%Filter%ON) then

            !Radius in number of cells
            call GetData(Me%Filter%Radius,                                              &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'FILTER_RADIUS',                                &
                         ClientModule = 'ModuleSand',                                   &
                         Default      = 4,                                              &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR29' 

            allocate (Me%Filter%Field2D(Me%Size%ILB: Me%Size%IUB, Me%Size%JLB : Me%Size%JUB))

            Me%Filter%Field2D(:,:) = FillValueReal


        endif

        call GetData(Me%SmoothSlope%ON,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SMOOTH_SLOPE',                                      &
                     ClientModule = 'ModuleSand',                                        &
                     Default      = .false.,                                             &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR37' 


        if (Me%SmoothSlope%ON) then

            !Slope be on each exist flux perpendicular to the Slope
            call GetData(Me%SmoothSlope%Critic,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'CRITICAL_SLOPE',                               &
                         ClientModule = 'ModuleSand',                                   &
                         Default      = 0.1,                                            &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR39' 

            !The flux perpendicular to the flux is a percentage of the paralel flux
            call GetData(Me%SmoothSlope%Factor,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'FLUX_SLOPE',                                   &
                         ClientModule = 'ModuleSand',                                   &
                         Default      = 0.1,                                            &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR49' 

        endif

        call GetData(Me%Boundary,                                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BOUNDARY',                                          &
                     ClientModule = 'ModuleSand',                                        &
                     Default      = NullGradient,                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR59' 

        call GetData(Me%TransportFactor,                                                 &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TRANSPORT_FACTOR',                                  &
                     default      = 1.,                                                  &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR69' 


        call GetData(Me%TauMax,                                                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TAU_MAX',                                           &
                     default      = 10.,                                                 &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR79' 


        call GetData(Me%Discharges%Yes,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DISCHARGES',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR89' 




    end subroutine ConstructGlobalParameters

    !----------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !If the user want's to use the values of a previous   
    ! run the read the sand properties values form the final      
    ! results file of a previous run. By default this      
    ! file is in HDF format                                

    subroutine ReadInitialField(FieldName, Field2D)

        !Arguments-------------------------------------------------------------
        character (Len = StringLength)              :: FieldName
        real, dimension(:,:), pointer               :: Field2D


        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: STAT_CALL
        integer                                     :: ObjHDF5
        integer(4)                                  :: HDF5_READ
        integer                                     :: ILB, IUB, JLB, JUB

        !----------------------------------------------------------------------
        
        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 

        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 


        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 


        inquire (FILE=trim(Me%Files%InitialSand)//"5", EXIST = EXIST)

cd0:    if (EXIST) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%InitialSand)//"5", HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSand - ERR03.'

            Field2D(:,:) = FillValueReal

            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB,            &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSand - ERR03.'

            call HDF5ReadData   (ObjHDF5, "/Results",trim(FieldName),                    &
                                 Array2D = Field2D,                                      &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSand - ERR03.'

            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSand - ERR06.'
        
        else if(.not. EXIST) then cd0

                stop 'ReadInitialField; ModuleSand - ERR07.'

        endif cd0


        !----------------------------------------------------------------------

    end subroutine ReadInitialField

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetSandDiameters (ObjSandID,  D35, D50, D90, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSandID
        real, dimension(:, :),  pointer, optional       :: D35, D50, D90
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(D35)) then
                call Read_Lock(mSand_, Me%InstanceID)
                D35 => Me%D35%Field2D
            endif

            if (present(D50)) then
                call Read_Lock(mSand_, Me%InstanceID)
                D50 => Me%D50%Field2D
            endif

            if (present(D90)) then
                call Read_Lock(mSand_, Me%InstanceID)
                D90 => Me%D90%Field2D
            endif


            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetSandDiameters
    
    !--------------------------------------------------------------------------
    
    subroutine GetSandDensity (ObjSandID, SandDensity, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSandID
        real                                            :: SandDensity
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            SandDensity = Me%Density

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetSandDensity

    !--------------------------------------------------------------------------

    subroutine UnGetSand2D_I(ObjSandID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSandID
        integer, dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSand_, Me%InstanceID, "UnGetSand2D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSand2D_I

    !--------------------------------------------------------------------------

    subroutine UnGetSand2D_R8(ObjSandID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSandID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSand_, Me%InstanceID,  "UnGetSand2D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSand2D_R8

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine UnGetSand2D_R4(ObjSandID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSandID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSand_, Me%InstanceID,  "UnGetSand2D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSand2D_R4

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifySand(ObjSandID, TauTotal, CurrentRugosity, WaveRugosity,   &
                          WaterColumn, VelU, VelV, VelMod, TauWave, TauCurrent, &
                          ShearVelocity, MinWaterColumn, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSandID
        real, dimension(:,:), pointer               :: TauTotal, CurrentRugosity, WaveRugosity, &
                                                       WaterColumn, VelU, VelV, VelMod,         &
                                                       TauWave, TauCurrent, ShearVelocity
        real,    intent(IN )                        :: MinWaterColumn      
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        logical                                     :: ChangeBathym
                                                       

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call ReadLockExternalVar

            if (Me%TransportMethod /= NoTransport) then

                do while (Me%ExternalVar%Now >= Me%Evolution%NextSand) 

                    Me%ExternalVar%TauTotal        => TauTotal

                    Me%ExternalVar%CurrentRugosity => CurrentRugosity
                    Me%ExternalVar%WaveRugosity    => WaveRugosity

                    Me%ExternalVar%WaterColumn     => WaterColumn

                    Me%ExternalVar%VelU            => VelU
                    Me%ExternalVar%VelV            => VelV
                    Me%ExternalVar%VelMod          => VelMod

                    Me%ExternalVar%TauWave         => TauWave
                    Me%ExternalVar%TauCurrent      => TauCurrent
                    Me%ExternalVar%ShearVelocity   => ShearVelocity

                    Me%ExternalVar%MinWaterColumn  =  MinWaterColumn

                    call ComputeFluxes

                    call ComputeEvolution            

                    if (Me%SmoothSlope%ON) then

                        call ComputeSmoothSlope

                    endif

                    !Boundary Condition
                    call BoundaryCondition(Me%DZ%Field2D)

                    call ComputeResidualEvolution

                    if (Me%Boxes%Yes ) call OutputBoxFluxes

                    if (Me%OutPut%Yes) call OutPutSandHDF

                    if (Me%TimeSerie)  call OutPut_TimeSeries

                    if (Me%Evolution%Bathym) then

                        ChangeBathym = .false.
                        !Selma
                        if    (Me%Evolution%BathymType == Time_) then
                    
                            if (Me%ExternalVar%Now >= Me%Evolution%NextBatim) ChangeBathym = .true.

                        else ! (Me%Evolution%BathymType == SelmaBathym_) then

                        !call ComputeSelmaChange (ChangeBathym)

                        endif

                        if (ChangeBathym) then

                            if (Me%Filter%Scheme == ModifyLax) then
                                call FilterModifyLax 
                            endif

                            !Bathymetry 
                            call UnGetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR10'

                            call ModifyGridData(Me%ObjBathym, Me%BatimIncrement%Field2D, Add = .false.,  &
                                                STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR20.'

                            !Bathymetry
                            call GetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR30'

                            Me%BatimIncrement%Field2D(:,:) = 0.

                            if    (Me%Evolution%BathymType == Time_) then
                    
                                Me%Evolution%NextBatim = Me%Evolution%NextBatim + Me%Evolution%BathymDT

                            endif

                        endif

                    endif

                    !Selma

                    Me%Evolution%NextSand = Me%Evolution%NextSand + Me%Evolution%SandDT

                enddo

            endif

            call ReadUnLockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySand

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine FilterModifyLax 


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                               :: AuxSum, Beta
        integer                            :: i, j, iw, jw, Counter
        !----------------------------------------------------------------------


        Beta = .5
      
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                Counter = 0
                AuxSum  = 0.
                do jw = j - Me%Filter%Radius, j + Me%Filter%Radius
                do iw = i - Me%Filter%Radius, i + Me%Filter%Radius

                    if (jw >= Me%WorkSize%JLB .and. jw <= Me%WorkSize%JUB .and.         &
                        iw >= Me%WorkSize%ILB .and. iw <= Me%WorkSize%IUB) then
                        if(Me%ExternalVar%WaterPoints2D(iw, jw) == WaterPoint) then
                            Counter = Counter + 1
                            AuxSum  = AuxSum + Me%BatimIncrement%Field2D(iw, jw)
                        end if
                    endif

                enddo
                enddo

                Me%Filter%Field2D(i, j) = Beta * Me%BatimIncrement%Field2D(i, j) + (1. - Beta) * AuxSum / real(Counter)

            endif
                                   
        enddo
        enddo

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                Me%BatimIncrement%Field2D(i,j) = Me%Filter%Field2D(i,j)

            endif
                                   
        enddo
        enddo
       

    end subroutine FilterModifyLax
    
    !--------------------------------------------------------------------------


    subroutine ComputeFluxes
        !Local-----------------------------------------------------------------
        integer             :: i, j
        real                :: SandThickness, FluxX, FluxY 
        !----------------------------------------------------------------------

        SELECT CASE (Me%TransportMethod)

        Case (MeyerPeter)
            call MeyerPeterTransport
        Case (Ackers)
            call AckersTransport
        Case (VanRijn1)
            call VanRijn1Transport
        Case (VanRijn2)
            call VanRijn2Transport
        Case (Bailard)
            call BailardTransport
        Case (Dibajnia)
            call DibajniaTransport
        Case (Bijker)
            call BijkerTransport
        
        
        END SELECT 

        Me%FluxX(:, :) =  0.
        Me%FluxY(:, :) =  0.

        !Computes the sand fluxes (m3/s) in the middle of the cells

        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

                SandThickness = Me%BedRock%Field2D(i, j) + Me%DZ_Residual%Field2D(i, j)
                    
                if (SandThickness > Me%SandMin .and. Me%ExternalVar%VelMod(i,j) > 0.) then
                
                    FluxX = Me%ExternalVar%VelU(i,j) / Me%ExternalVar%VelMod(i,j) * Me%ExternalVar%DVY(i, j)
                    FluxY = Me%ExternalVar%VelV(i,j) / Me%ExternalVar%VelMod(i,j) * Me%ExternalVar%DUX(i, j)

                    Me%FluxX(i, j) = Me%TransportCapacity(i, j) * FluxX * Me%TransportFactor
                    Me%FluxY(i, j) = Me%TransportCapacity(i, j) * FluxY * Me%TransportFactor

                endif

            endif

        enddo
        enddo

        call BoundaryCondition(Me%FluxX)
        call BoundaryCondition(Me%FluxY)

    end subroutine ComputeFluxes
    !--------------------------------------------------------------------------

    subroutine MeyerPeterTransport
        !Local-----------------------------------------------------------------
        real :: Tr1,Tr2,Tr3, DeltaTau, Miu, CChezy, Clinha, Depth
        integer :: i, j

        !----------------------------------------------------------------------
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

                Depth = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                CChezy= 18.*ALOG10(12.*Depth/Me%ExternalVar%CurrentRugosity(I,J))
                Clinha= 18.*ALOG10(12.*Depth/Me%D90%Field2D(I,J))
                Miu   = (CChezy/Clinha)**1.5

                DeltaTau = min(Miu*Me%ExternalVar%TauTotal(I,J),Me%TauMax)-Me%TauCritic(I,J)        

! ---> If Tau < Tau critic there is no transport

                if (DeltaTau.LE.0.) Then
                    Me%TransportCapacity(i, j) = 0.
                    Cycle
                endif

                Tr1    = (Me%RelativeDensity*Gravity)**0.5
                Tr2                        = Me%D50%Field2D(I,J)**1.5
                Tr3                        = DeltaTau**1.5
                Me%TransportCapacity(i, j) = 8.*Tr1*Tr2*Tr3
            
            else
                Me%TransportCapacity(i, j) = 0.  !instrução adicionada
            endif

        enddo
        enddo


    end subroutine MeyerPeterTransport
    !--------------------------------------------------------------------------


    subroutine AckersTransport

        !Local-----------------------------------------------------------------
        real    :: Alfa, EP, VisCin2, R32, RhoSW, Dgr, Depth
        real    :: Vn, Va, Vm, Vc1, Vc, Fgr1, Fgr2, Fgr, Ggr1, Ggr, S1
        integer :: i, j
 
        !----------------------------------------------------------------------
 
! ---> Global Parameters

        Alfa   = 10.
        EP     = 1./3.
        VisCin2= WaterCinematicVisc**2.
        R32    = SQRT(32.)
        RhoSW  = Me%Density / Me%ExternalVar%WaterDensity
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

! ---> Dimensionless Grain Diameter
                
                Dgr  = Me%D35%Field2D(I,J)*(Gravity*Me%RelativeDensity/VisCin2)**EP

                if (Dgr.LT.1.) then
                  Write (*,*) i,j,Me%D35%Field2D(I,J)
                  Write (*,500)
                  stop 'AckersTransport - ModuleSand - ERR01'
                endif

                if (Dgr.GT.60.) then
                  Vn = 0.00
                  Va = 0.17
                  Vm = 1.78    ! 1.5 in 1973 version
                  Vc = 0.025
                else
                  Vn = 1.-0.56*Log10(Dgr)
                  Va = 0.14+0.23/(SQRT(Dgr))
                  Vm = 1.67+6.83/Dgr                              ! 1.34+9.66/Dgr in 1973 version
                  Vc1= 2.79*Log10(Dgr)-0.98*(Log10(Dgr))**2.-3.46 ! 2.86*Log10(Dgr)-(Log10(Dgr))**2.-3.53 in 1973 version
                  Vc = 10.**(Vc1)
                endif

                ! ---> Sediment Mobility (Fgr)


                Depth = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                Fgr1 = Me%ExternalVar%ShearVelocity(I,J)**Vn/(SQRT(Gravity*Me%D35%Field2D(I,J)*Me%RelativeDensity))

                Fgr2 = R32*Log10(Alfa*Depth/Me%D35%Field2D(I,J))

                Fgr  = Fgr1*(Me%ExternalVar%VelMod(I,J)/Fgr2)**(1-Vn)

                ! ---> Sediment Transport (Ggr e m3/s de sedimento por m3/s de fluxo de agua)

                Ggr1 = Fgr/Va
                if (Ggr1.GE.1.) then
                  Ggr  = Vc*(Ggr1-1)**Vm
                  S1   = (Me%ExternalVar%ShearVelocity(I,J)/Me%ExternalVar%VelMod(I,J))**Vn

                  Me%TransportCapacity(i, j) = Ggr*RhoSW*Me%D35%Field2D(I,J)/S1               !(m)
                  
                  ! [m3/s/m] = [m/s]*[m]
                  Me%TransportCapacity(i, j) = Me%ExternalVar%VelMod(I,J) * Me%TransportCapacity(i, j)

                else
                  Me%TransportCapacity(i, j) = 0.
                endif
            
            else
                  Me%TransportCapacity(i, j) = 0.
            endif

        enddo
        enddo

500 Format (/////T10,'A T E N C A O !!!',                        &
              //,T10,'O Valor do Diametro Adimensional e < 1',   &
              //T10,'Verifique os dados do problema ...',///)

    end subroutine AckersTransport

    !--------------------------------------------------------------------------
    
    subroutine VanRijn1Transport


        
        !Local-----------------------------------------------------------------
        real    :: MUw,Kapa,TetaCrit,TauC, TauW, D50, Abw
        real    :: MUc,Fc,Fc1, Clim, Ksc, Ksw
        real    :: CC,DelS,UwUcRel,Gamma,ApparentRoughness,Alfaw
        real    :: TauC_ef,TauW_ef,TauTot_ef,Thet1,Ttran1,Ttrans,Depth
        real    :: AltRock, VelQ, SuspLoad, BedLoad, Ca, Ca1, &
                   Ca2, Uast, Uastc, AA, MixLay, Zw, FAKTW, Res3, Res4, Beta, &
                   ZC, Res1, Res2, FAKTC, phi, phiWave, CrossDot
        integer :: i, j
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSand", "VanRijn1Transport")

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%ExternalVar%VelMod(I,J)>0.) then
                
                Depth   = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                Clim    = 0.01*Depth 
                Ksc     = Me%ExternalVar%CurrentRugosity(I,J)
                Ksw     = Me%ExternalVar%WaveRugosity   (I,J)
                TauW    = Me%ExternalVar%TauWave        (I,J)
                TauC    = Me%ExternalVar%TauCurrent     (I,J)
                D50     = Me%D50%Field2D                (I,J)
                Abw     = Me%ExternalVar%Abw            (i,j)
                AltRock = Me%BedRock%Field2D(I,J)
                VelQ    = FallVel(D50)

                !IF (Depth.LT.Hmin.OR.AltRock.LE.0.01) Then
                IF (AltRock.LE.0.01) Then
                  Cycle
                Endif

                !IF (Depth.LT.Clim) Depth = Clim

                !Compute reference concentration
                CC    = 18.*Alog10(12.*Depth/Ksc)
                Uast  = Sqrt(Gravity)*Me%ExternalVar%VelMod(I,J)/CC
                DelS  = 0.
                If (AbW.Gt.0.) Then
                   DelS = 0.216*Abw*(Abw/Ksw)**(-0.25)
                Endif
                
                UwUcRel = Me%ExternalVar%Ubw(I,J)/Me%ExternalVar%VelMod(I,J)
                If (UwUcRel.Ge.2.5) UwUcRel=2.5

                !Wave direction
                phiWave    = Me%ExternalVar%WaveDirection(i, j) * Pi/180.

                !Compute the angle between waves and currents.
                CrossDot = Me%ExternalVar%VelU(i,j) * cos(phiWave) + & 
                           Me%ExternalVar%VelV(i,j) * sin(phiWave)  

                if( abs(CrossDot / Me%ExternalVar%VelMod(i,j)) .LT. 1.0) then
                    phi = acos(CrossDot / Me%ExternalVar%VelMod(i,j)) * 180. / Pi                   
                endif     

                !Compute Apparent Roughness
                if (phi < 90) then
                    Gamma = 0.75
                else
                    Gamma = 0.75 + (0.35) * (phi - 90.) / 90.
                endif

                ApparentRoughness= Exp(Gamma*UwUcRel)*Ksc

                Fc  = 0.24*Alog10(12.*Depth/Ksc)**(-2.)
                Fc1 = 0.24*Alog10(12.*Depth/(3.*Me%D90%Field2D(I,J)))**(-2.)

                MUc = Fc1/Fc                       !efficiency factor currents
                MUw = 0.6/Me%Dast(I,J)             !efficiency factor waves

                If (DelS.Le.Ksc/10.) Then
                   Alfaw = 1.0
                Else
                   Alfaw = (Alog(30.*DelS/ApparentRoughness)/Alog(30.*DelS/Ksc))**2
                Endif

                !Compute Transport term

                TauW_ef   = TauW*MUw
                TauC_ef   = TauC*MUc*Alfaw
                TauTot_ef = TauC_ef+TauW_ef

                Thet1     = TauTot_ef / (Me%RhoSl*Gravity*Me%D50%Field2D(I,J))

                !Compute bed-shear stress parameter Ttrans

                TetaCrit = Me%TauCritic(i, j) / (Me%RhoSl*Gravity*Me%D50%Field2D(I,J))

                Ttran1   = (Thet1-TetaCrit)/TetaCrit
                Ttrans   = Max(0.0001,Ttran1)

                !Compute bed-load Transport

                Uastc   = Sqrt(TauC_ef/Me%ExternalVar%WaterDensity)
                BedLoad = 0.25*Uastc*Me%D50%Field2D(I,J)*Ttrans**(1.5)/Me%Dast(I,J)**(0.3)    !m3/m.s

                If (Abs(BedLoad).LT.1.E-30) Cycle

                !Compute suspended-load Transport

                AA     = 7.
                Mixlay = 2.* Ksc
                If (Depth/MixLay.LE.100) AA = 0.7*Sqrt(Depth/MixLay)
                ZW   = 0.
                FAKTW= 0.
                
                If (Me%ExternalVar%WavePeriod(i,j).GT.0..AND.Me%ExternalVar%WaveHeight(I,J).GT.0.)Then
                    ZW    = AA*(VelQ/Me%ExternalVar%VelMod(I,J))**.9*                                   &
                            (Me%ExternalVar%VelMod(I,J)*Me%ExternalVar%WavePeriod(i,j)/Me%ExternalVar%WaveHeight(I,J))**1.05

                    ZW    = Min(ZW,25.)
                    Res3  = (Clim/Depth)**ZW-(Clim/Depth)**1.2
                    Res4  = (1.-Clim/Depth)**ZW*(1.2-ZW)
                    If (abs(Res3).Lt.1.E-5.Or.abs(Res4).Lt.1.E-5) Then
                        FAKTW = 0.
                    Else
                        FAKTW = Res3/Res4
                    Endif
                Endif
                Kapa     = 0.4
                Beta     = 1.+2.*(VelQ/Uastc)**2.
                IF (Beta.GT.2.) Beta=2.
                ZC       = VelQ/(Kapa*Beta*Uast)
                ZC       = Min(ZC,25.)
                Res1     = (Clim/Depth)**ZC-(Clim/Depth)**1.2
                Res2     = (1.-Clim/Depth)**ZC*(1.2-ZC)
                If (abs(Res2).Lt.1.E-5.Or.abs(Res1).Lt.1.E-5) Then
                    FAKTC    = 0.
                Else
                    FAKTC    = Res1/Res2
                Endif
                Ca1      = Me%D50%Field2D(I,J)/Clim
                Ca2      = Ttrans**1.5*Me%Dast(I,J)**(-0.3)
                Ca       = 0.015*Ca1*Ca2
                SuspLoad = (FAKTC+FAKTW)*Me%ExternalVar%VelMod(I,J)*Depth*Ca    !m3/m.s


                Me%TransportCapacity(i, j) = (SuspLoad + BedLoad) 

            else

                Me%TransportCapacity(i, j) = 0.

            endif

        enddo
        enddo
    
    if (MonitorPerformance) call StopWatch ("ModuleSand", "VanRijn1Transport")

    end subroutine VanRijn1Transport

    !--------------------------------------------------------------------------

    subroutine VanRijn2Transport



        !Local-----------------------------------------------------------------
        real    :: DS, Depth
        real    :: Ksc,Ksw,VelMod,AltRock,WavePeriod,WaveHeight,D50,D90,RhoWater,RhoSandDry
        real    :: VisCin,TauW,Ubw,Abw,Dast,SuspLoad, BedLoad,Kapa,Fw,phiWave,CrossDot
        real    :: phi
        real    :: MUw,MUc,Fc,Fc1
        real    :: CC,DelS,UwUcRel,Gamma,ApparentRoughness,Alfaw
        real    :: TauC_ap,TauC_ef,TauW_ef,TauTot_ef,Thet1,Ttran1,Ttrans
        real    :: Clim,TetaCrit
        real    :: VelQ
        real    :: Uast,Uastc,Fca
        real    :: Ca,Ca1,Ca2,Beta
        real    :: Alt_Conf,dym,dxm,dyx,NN,HULP30,JTAL,YOLD
        real    :: Abr,Ebw,EmaxW,EmaxC,Esw,Esc,Es,fcc,FF,SS,NTEL,Y
        real    :: UDEL,Z,YPRIME,C,UC,TERM1,TERM2,XEND,XOLD,IT
        
        logical :: BO
        integer :: i, j

        !----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleSand", "VanRijn2Transport")

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%ExternalVar%VelMod(I,J)>0.) then
                
                Depth   = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                Ksc       =Me%ExternalVar%CurrentRugosity(i,j)
                Ksw       =Me%ExternalVar%WaveRugosity   (i,j)
                VelMod    =Me%ExternalVar%VelMod         (i,j)
                AltRock   =Me%BedRock%Field2D            (i,j)
                WavePeriod=Me%ExternalVar%WavePeriod     (i,j)
                WaveHeight=Me%ExternalVar%WaveHeight     (i,j)
                D50       =Me%D50%Field2D                (i,j)
                D90       =Me%D90%Field2D                (i,j)
                RhoWater  =Me%ExternalVar%WaterDensity
                RhoSandDry=Me%RhoSl
                VisCin    =WaterCinematicVisc
                TauW      =Me%ExternalVar%TauWave        (I,J)
                !TauC      =Me%ExternalVar%TauCurrent     (I,J)               
                Ubw       =Me%ExternalVar%Ubw            (i,j)
                Abw       =Me%ExternalVar%Abw            (i,j)
                Dast      =Me%Dast                       (i,j)
                VelQ      =FallVel(D50)

                !IF (Depth.LT.Hmin.OR.AltRock.LE.0.01) Then
                if (AltRock.LE.0.01) Then
                  cycle
                endif

                BedLoad     = 0.
                SuspLoad    = 0.
                Kapa        = 0.4
                Ds          = 5.E-2

                !Compute friction factor due to waves
                Fw=exp(-6.+5.2*(Abw/ksw)**(-0.19))
                if(Fw.GT.0.3) Fw=0.3 
                
                !Compute chezy and friction velocity
                CC    = 18.*alog10(12.*Depth/Ksc)   
                Uast  = sqrt(Gravity)*VelMod/CC  

                !Compute DelS, boundary layer width
                DelS  = 0.
                if (Abw.Gt.0.) then
                    DelS = 0.216*Abw*(Abw/Ksw)**(-0.25)
                endif            
               
                !Compute Apparent Roughness
                UwUcRel = Ubw/VelMod
                if(UwUcRel.Ge.2.5) UwUcRel=2.5
                
                    !Wave direction
                phiWave    = Me%ExternalVar%WaveDirection(i, j) * Pi/180

                    !Compute the angle between waves and currents.
                CrossDot = Me%ExternalVar%VelU(i,j) * cos(phiWave) + & 
                           Me%ExternalVar%VelV(i,j) * sin(phiWave)                 
                
                if( abs(CrossDot / Me%ExternalVar%VelMod(i,j)) .LT. 1.0) then
                    phi = acos(CrossDot / Me%ExternalVar%VelMod(i,j)) * 180. / Pi                   
                endif             
                           
                    !gamma
                if (phi < 90) then
                    Gamma = 0.75
                else
                    Gamma = 0.75 + (0.35) * (phi - 90.) / 90.
                endif
                
                ApparentRoughness=exp(Gamma*UwUcRel)*Ksc

                !Compute efficiency factor for currents and waves, MUc e MUw
                Fc   = 0.24*alog10(12.*Depth/Ksc)**(-2.)
                Fca  = 0.24*alog10(12.*Depth/ApparentRoughness)**(-2.)
                Fc1  = 0.24*alog10(12.*Depth/(3.*D90))**(-2.)

                MUc  = Fc1/Fc          !efficiency factor for currents
                MUw  = 0.6/Dast        !efficiency factor for waves

                !Compute wave parameter

                if(DelS.LE.Ksc/10.) then
                    Alfaw = 1.0
                else
                    Alfaw = (alog(30.*DelS/ApparentRoughness)/alog(30.*DelS/Ksc))**2.
                endif

                !Compute effective shear stress for waves and currents
                TauC_ap  = 0.125*RhoWater*Fca*VelMod**2.
                TauC_ef  = TauC_ap*MUc*Alfaw   
                TauW_ef  = TauW*MUw
                !TauW_ap  = 0.25*RhoWater*Fw*Ubw**2.
                !TauW_ef  = TauW_ap*MUw
                TauTot_ef= TauC_ef+TauW_ef
     
                !Compute bed-shear stress parameter Ttrans
                Thet1     = TauTot_ef / (RhoSandDry*Gravity*D50)
                TetaCrit  = Me%TauCritic(i,j) / (RhoSandDry*Gravity*D50)
                Ttran1    = (Thet1-TetaCrit)/TetaCrit
                Ttrans    = max(0.0001,Ttran1)

                !Compute reference level
                Alt_Conf = depth*0.11*(D50/depth)**0.3*(1-exp(-0.5*Ttrans))*(25.-Ttrans)
                if(Ttrans.GE.0.0.AND.Ttrans.LE.25.0) then 
                    Clim     = max(0.01*depth, 0.5*Alt_Conf)
                else
                    Clim     = max(0.01*depth, Ksc)
                endif

                !Compute BED-LOAD TRANSPORT
                Uastc    = Sqrt(TauC_ef/RhoWater)
                BedLoad  = 0.25*Uastc*D50*Ttrans**(1.5)/Dast**(0.3)    !m3/m.s

                !Compute reference concentration in the bed
                Ca1      = D50/Clim
                Ca2      = Ttrans**1.5*Dast**(-0.3)
                Ca       = 0.015*Ca1*Ca2

                !Compute Beta
                Beta     = 1.+2.*(VelQ/Uastc)**2. 
                if (Beta.GT.2.) Beta=2. 

                !Compute SUSPENDED-LOAD TRANSPORT by numerical integration of V*C over vertical
                NN  =12 
                JTAL=8
                NN  =JTAL*NN
                DYM =CA/NN
                DXM =depth/NN
                DYX =DYM/DXM
                BO  =.FALSE.
                HULP30=-1.+ alog(30.*depth/ApparentRoughness)
                if(DELS.GT.0.) then 
                UDEL=VelMod*alog(30.*DELS/ApparentRoughness)/HULP30
                endif

                !Starting the integral cicle
                ABR=max(3.*(WaveHeight/depth)-.8,1.)
                EBW=.004*Dast*ABR*DS*Ubw
                if(WavePeriod.GT.1.E-4) then
                    EMAXW=0.035*ABR*depth*WaveHeight/WavePeriod
                else
                    EMAXW=0.
                endif
                if(EMAXW.LE.EBW) EMAXW=EBW
                EMAXC=.25*kapa*Uast*depth*BETA
    
                C=CA
                Z=Clim
                if(Z.LE.DS)ESW=EBW
                if(Z.GT.DS.AND.Z.LE.0.5*depth)ESW=EBW+(EMAXW-EBW)*((Z-DS)/(0.5*depth-DS))
                if(Z.GE.0.5*depth)ESW=EMAXW
                if(Z.GE.0.5*depth)ESC=EMAXC
                if(Z.LT.0.5*depth)ESC=EMAXC-EMAXC*(1.-2.*Z/depth)**2
                ES=(ESW**2.+ESC**2.)**0.5
                fcc=0.
                if(C.GT.1.E-8) then
                    if(Z.GE.Clim) fcc=-VelQ/ES*C*(1.-C)**5.
                endif
                YPRIME=fcc
                FF=1./CA*YPRIME
                if(DELS.GT.0.)then
                    UC=UDEL*alog(30.*Clim/ksc)/alog(30.*DELS/ksc)
                endif
                if(Clim.GE.DELS) UC=VelMod*alog(30.*Clim/ApparentRoughness)/HULP30

 
                !Further Integration
                Y    =CA
                TERM1=UC*Y
                XEND =Clim
                SS   =0.
                NTEL =0
                IT   =2

                do while (.NOT.BO)
                    NTEL= NTEL+1
                    XOLD=XEND
                    YOLD=Y
                    if(-YPRIME.GT.DYX) then
                        Y=YOLD-DYM
                        if(Y.LT.2./3.*YOLD) Y=2./3.*YOLD
                        XEND=XOLD+alog(Y/YOLD)/FF
                    else 
                        XEND=XOLD+DXM
                        if(XEND.GE.depth)then
                            XEND=depth
                            BO=.TRUE.
                        endif
                        Y=exp(alog(YOLD)+(XEND-XOLD)*FF)
                    endif
                    C=Y
                    Z=XEND
                    if(Z.LE.DS)ESW=EBW
                    if(Z.GT.DS.AND.Z.LE.0.5*depth)ESW=EBW+(EMAXW-EBW)*((Z-DS)/(0.5*depth-DS))
                    if(Z.GE.0.5*depth)ESW=EMAXW
                    if(Z.GE.0.5*depth)ESC=EMAXC
                    if(Z.LT.0.5*depth)ESC=EMAXC-EMAXC*(1.-2.*Z/depth)**2
                    ES=(ESW**2.+ESC**2.)**0.5
                    fcc=0.
                    if(C.GT.1.E-8) then
                        if(Z.GE.Clim) fcc=-VelQ/ES*C*(1.-C)**5.
                    endif
                    YPRIME=fcc
                    FF=1./Y*YPRIME
                    if(DELS.GT.0.)then
                        UC=UDEL*alog(30.*XEND/ksc)/alog(30.*DELS/ksc)
                    endif
                    if(XEND.GE.DELS) UC=VelMod*alog(30.*XEND/ApparentRoughness)/HULP30
                    TERM2=UC*Y
                    SuspLoad=SuspLoad+(XEND-XOLD)*(TERM1+TERM2)/2.
                    TERM1=TERM2
    
                enddo

                Me%TransportCapacity(i, j) = SuspLoad + BedLoad

            else

                Me%TransportCapacity(i, j) = 0.

            endif

        enddo
        enddo
    
    if (MonitorPerformance) call StopWatch ("ModuleSand", "VanRijn2Transport")

    end subroutine VanRijn2Transport

    !--------------------------------------------------------------------------

    subroutine BailardTransport


        
        !Local-----------------------------------------------------------------
        real    :: BedEff,SusEff,Depth,WavePeriod,ksw,Ksc,D50,VisCin,VelMod,Ubw
        real    :: Abw,Fw,WaveHeight
        real    :: BedLoad,SuspLoad,TotLoad
        real    :: DynAng,DynAngr,TanDyn,TanBet
        real    :: Fc,CChezy
        real    :: cte1,cte2,cte3,cte4,TW,Uc,varT,delT,VelInst
        real    :: VelQ
        real    :: phi,phiWave,CrossDot
        integer :: i,j,g
        real,dimension(1:50)   ::varI

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSand", "BailardTransport")

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

        BedEff   = 0.1  !typically used
        SusEff   = 0.02 !typically used
        BedLoad  = 0.
        SuspLoad = 0.
            
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%ExternalVar%VelMod(I,J)>0.) then
                
                Depth   = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                WavePeriod=Me%ExternalVar%WavePeriod(i,j)
                WaveHeight=Me%ExternalVar%WaveHeight(i,j)
                ksw       =Me%ExternalVar%WaveRugosity(i,j)
                Ksc       =Me%ExternalVar%CurrentRugosity(i,j)
                D50       =Me%D50%Field2D(i,j)
                VisCin    =WaterCinematicVisc
                VelMod    =Me%ExternalVar%VelMod(i,j)
                Ubw       =Me%ExternalVar%Ubw(i,j)
                Abw       =Me%ExternalVar%Abw(i,j)
                VelQ      =FallVel(D50)

                !Compute the angle between waves and currents.
                phiWave    = Me%ExternalVar%WaveDirection(i, j) * Pi/180
                CrossDot   = Me%ExternalVar%VelU(i,j) * cos(phiWave) + & 
                             Me%ExternalVar%VelV(i,j) * sin(phiWave)                 
                
                if( abs(CrossDot / Me%ExternalVar%VelMod(i,j)) .LT. 1.0) then
                    phi = acos(CrossDot / Me%ExternalVar%VelMod(i,j)) * 180. / Pi                   
                endif    

                !Compute friction factor due to waves
                Fw=exp(-6.+5.2*(Abw/ksw)**(-0.19))
                if(Fw.GT.0.3) Fw=0.3 
     
                !dynamic friction factor and bottom Slope
                DynAng = 30.               ! degrees
                DynAngr= DynAng*pi/180.    ! radians
                TanDyn = Tan(DynAngr)      ! TanDyn=0.63, typically used value.
                TanBet = 0.                ! Slope

                !Compute total Transport
                cte1=0.5*Me%ExternalVar%WaterDensity*Fw*BedEff/(Me%RhoSl*Gravity*TanDyn)
                cte2=0.5*Me%ExternalVar%WaterDensity*Fw*SusEff/(Me%RhoSl*Gravity*VelQ)


                if(WavePeriod.GT.1.E-4.AND.WaveHeight.GT.1.E-4)then !Os Ubw's passam a Cos(teta)*Ubw
                    TW=WavePeriod
                    Uc=VelMod
                    varT=0                                                
                    varI(1)=(sqrt(Uc**2.0+(cos(phi)*Ubw*sin(2.*pi*varT/TW))**2.+2.*Uc*cos(phi)*Ubw*sin(2.*Pi*varT/TW)))**3. 
                    delT=TW/12
                    varT=delT
                    do g=2,50
                        if(varT.LE.TW) then                                                                       
                          varI(g)=(sqrt(abs(Uc**2.0+(cos(phi)*Ubw*sin(2.*Pi*varT/TW))**2.0+2.0*Uc*cos(phi)*&
                                  Ubw*sin(2.*Pi*varT/TW))))**3.
                          VelInst=VelInst+(varI(g)+varI(g-1))/2*delT
                          varT=varT+delT
                        else
                        EXIT
                        endif
                    enddo
                    VelInst=VelInst/TW
                    TotLoad=cte1*(2.*Uc*Ubw**2.+Uc**3)+cte2*(Uc*Ubw**2*VelInst)

                else         !situação sem ondas
                    CChezy = 18.*Alog10(12.*Depth/Ksc)
                    Fc=2.*Gravity/CChezy**2            !diferente do Fc em Van Rijn (1/4)
                    cte3=0.5*Me%ExternalVar%WaterDensity*Fc*BedEff/(Me%RhoSl*Gravity*TanDyn)
                    cte4=0.5*Me%ExternalVar%WaterDensity*Fc*SusEff/(Me%RhoSl*Gravity*VelQ)
                    TotLoad=cte3*VelMod**3.+cte4*VelMod**4.   
                endif
                
                Me%TransportCapacity(i, j) = TotLoad  

            else
                
                Me%TransportCapacity(i, j) = 0.
  
            endif
    enddo
    enddo
    
    if (MonitorPerformance) call StopWatch ("ModuleSand", "BailardTransport")

    end subroutine BailardTransport
        
    !--------------------------------------------------------------------------

    subroutine DibajniaTransport        


        
        !Local-----------------------------------------------------------------
        real    :: TAU,AbsTAU,Adw,Bdw
        real    :: Depth,WavePeriod,WaveHeight,ksw,Ksc,D50,VisCin,VelMod,Ubw
        real    :: Abw,Fw,TotLoad
        real    :: VelQ
        real    :: delT,varT,NN,INT1,INT2,Twc,Twt,Zerot,argu,TW,Uwc,Uwt,Fc
        real    :: transC,transT,OmegaC,OmegaT,OmegaC_susp,OmegaT_susp
        real    :: phiWC,phiW,phiC,denom,cte,phase_lagT,phase_lagC,critical
        real    :: phi,phiWave,CrossDot

        integer :: i,j,g
        real,dimension(1:50)   ::VelInst1,VelInst2

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSand", "DibajniaTransport")

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            Adw           =.0015
            Bdw           =.55
            TotLoad       = 0.
                
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%ExternalVar%VelMod(I,J)>0.) then
                
                Depth   = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                WavePeriod=Me%ExternalVar%WavePeriod(i,j)
                WaveHeight=Me%ExternalVar%WaveHeight(i,j)
                ksw       =Me%ExternalVar%WaveRugosity(i,j)
                D50       =Me%D50%Field2D(i,j)
                VisCin    =WaterCinematicVisc
                VelMod    =Me%ExternalVar%VelMod(i,j)
                Ksc       =Me%ExternalVar%CurrentRugosity(i,j)
                Ubw       =Me%ExternalVar%Ubw(i,j)
                Abw       =Me%ExternalVar%Abw(i,j)
                VelQ      =FallVel(D50)

                !Compute the angle between waves and currents.
                
                phiWave    = Me%ExternalVar%WaveDirection(i, j) * Pi/180

                CrossDot = Me%ExternalVar%VelU(i,j) * cos(phiWave) + & 
                           Me%ExternalVar%VelV(i,j) * sin(phiWave)                 
                
                if( abs(CrossDot / Me%ExternalVar%VelMod(i,j)) .LT. 1.0) then
                    phi = acos(CrossDot / Me%ExternalVar%VelMod(i,j)) * 180. / Pi                   
                endif    
                    
                !Compute friction factors due to waves and currents
                Fw=exp(-6.+5.2*(Abw/ksw)**(-0.19))
                if(Fw.GT.0.3) Fw=0.3
                Fc   = 0.24*Alog10(12.*Depth/Ksc)**(-2.)
               
                !Compute the half-periods
                argu=abs(-VelMod/Ubw)
                TW=WavePeriod
                if(TW.LT.1.E-6) TW=1.E-6
                if(argu.LE.1) then
                    Zerot=TW/(2.*Pi)*asin(-VelMod/Ubw) !o input de asin é -1 a 1 
                    Twc =TW/2.+2*abs(Zerot)
                    Twt =TW/2.-2*abs(Zerot)
                else
                    zeroT=0.
                    Twc=TW
                    Twt=1.E-7
                endif

                !Compute the velocity over each half-period
                INT1=0.
                NN=20
                varT=Zerot
                VelInst1(1)=(VelMod+cos(phi)*Ubw*sin(2.*pi*varT/TW))**2.
                delT=Twc/NN
                varT=delT
                do g=2,30
                    if(varT.LE.(Twc+Zerot)) then
                        VelInst1(g)=(VelMod+cos(phi)*Ubw*sin(2.*Pi*varT/TW))**2.  !VelMod*cos(teta) ou Ubw*cos(teta)
                        INT1=INT1+(VelInst1(g)+VelInst1(g-1))/2.*delT   
                        varT=varT+delT
                    else
                    EXIT
                    endif
                enddo

                INT2=0.
                VelInst2(1)=(VelMod+cos(phi)*Ubw*sin(2.*Pi*varT/TW))**2.
                delT=Twt/NN
                varT=varT+delT
                do g=2,30
                    if(varT.LE.(TW+Zerot)) then
                        VelInst2(g)=(VelMod+cos(phi)*Ubw*sin(2.*Pi*varT/TW))**2.  !VelMod*cos(teta) ou Ubw*cos(teta)
                        INT2=INT2+(VelInst2(g)+VelInst2(g-1))/2.*delT
                        varT=varT+delT
                    else
                        EXIT
                    endif
                enddo

                        Uwc=sqrt(2/Twc*INT1) !Uwc=sqrt(2/Twc*INT1+2.*(VelMod*sin(teta))**2.)
                        Uwt=sqrt(2/Twt*INT2) !Uwt=sqrt(2/Twt*INT2+2.*(VelMod*sin(teta))**2.)

                !Compute the critical parameter

                denom=Me%RelativeDensity*Gravity*D50
                phiW =0.5*Fw*Ubw**2/denom                        !substituí 0.5 por 1.4
                phiC =0.5*Fw*VelMod**2/denom                   
                phiWC=max(phiC,phiW)                             !maximum Shields parameter
                !phiWC= .3/(1.+1.2*Dast)+.055*(1-exp(-.02*Dast)) !SOULSBY, 1997
 
                    if(phiWC.LE.0.2) then
                        critical=0.03
                    elseif(phiWC.GT.0.2.AND.phiWC.LT.0.6) then
                        critical=1.-0.97*(1-((phiWC-.2)/.4)**2)**.5
                    else
                        critical=1.
                    endif

                ! ---> Compute the phase-lag parameter 

                phase_lagC=Uwc**2/(2.*Me%RelativeDensity*Gravity*VelQ*Twc)
                phase_lagT=Uwt**2/(2.*Me%RelativeDensity*Gravity*VelQ*Twt)

                ! ---> Compute the Omegas: amount of sand in movement during and after the half-periods

                cte=2.*VelQ/D50
                if(phase_lagC.LE.critical) then
                    OmegaC     = phase_lagC*cte*Twc
                    OmegaC_susp= 0.
                else
                    OmegaC     = cte*Twc
                    OmegaC_susp= (phase_lagC-1)*cte*Twc
                endif

                if(phase_lagT.LE.critical) then
                    OmegaT     = phase_lagT*cte*Twt
                    OmegaT_susp= 0.
                else
                    OmegaT     = cte*Twt
                    OmegaT_susp= (phase_lagT-1)*cte*Twt
                endif

                TransC=OmegaC**3.+OmegaT_susp**3.
                TransT=OmegaT**3.+OmegaC_susp**3.
                
                TAU=0.
                if(WavePeriod.GT.1.E-4 .AND. WaveHeight.GT.1.E-4)then
                    TAU=(Twc*Uwc*TransC-Twt*Uwt*TransT)/(Uwc+Uwt)/TW !SINAL MENOS por que a vel. inst. está ao quadrado
                else 
                    TAU=(Twc*Uwc*TransC+Twt*Uwt*TransT)/(Uwc+Uwt)/TW !mas sem ondas o sinal deve ser (+)
                endif
                    AbsTAU=abs(TAU)                                  
                    TotLoad=Adw*VelQ*D50*AbsTAU**Bdw
                                   
                    Me%TransportCapacity(i, j) = TotLoad  

            else
                
                    Me%TransportCapacity(i, j) = 0.
  
            endif
    enddo
    enddo
    
    if (MonitorPerformance) call StopWatch ("ModuleSand", "DibajniaTransport")

    end subroutine DibajniaTransport


    !--------------------------------------------------------------------------

    subroutine BijkerTransport
      
        !Local-----------------------------------------------------------------
        real    :: Depth,WavePeriod,WaveHeight,ksw,Ksc,D50,D90,VisCin,VelMod
        real    :: RhoWater,RhoSandDry,Fw,Ubw,Abw
        real    :: BedLoad,SuspLoad,TotLoad
        real    :: Tinic,TauTot,Ttrans
        real    :: Uast,CChezy,ChezySed,Clinha,Fc
        real    :: Miu,Coef
        real    :: VelQ
        real    :: TauC,qsi,Dast
        real    :: Kappa,KS,KB,AbsTautot,ZX,AAA,B1,B2
        integer :: i,j

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleSand", "BijkerTransport")

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
       
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%ExternalVar%VelMod(I,J)>0.) then
                
                Depth   = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) Depth = Me%ExternalVar%MinWaterColumn

                WavePeriod=Me%ExternalVar%WavePeriod(i,j)
                WaveHeight=Me%ExternalVar%WaveHeight(i,j)
                ksw       =Me%ExternalVar%WaveRugosity(i,j)
                Ksc       =Me%ExternalVar%CurrentRugosity(i,j)              
                D50       =Me%D50%Field2D(i,j)
                D90       =Me%D90%Field2D(i,j)
                VisCin    =WaterCinematicVisc
                VelMod    =Me%ExternalVar%VelMod(i,j)
                Ubw       =Me%ExternalVar%Ubw(i,j)
                Abw       =Me%ExternalVar%Abw(i,j)
                Dast      =Me%Dast(i,j)
                RhoWater  =Me%ExternalVar%WaterDensity
                RhoSandDry=Me%RhoSl
                VelQ      =FallVel(D50)

                Coef      = 5.
                Kappa     = 0.384
                BedLoad   = 0.
                SuspLoad  = 0.

                !Compute friction factors due to waves and currents
                Fw=exp(-6.+5.2*(Abw/ksw)**(-0.19))
                if(Fw.GT.0.3) Fw=0.3
                
                !Chezy coefficient and CLinha
                CChezy   = 18.*Alog10(12.*Depth/Ksc)
                CLinha   = 18.*Alog10(12.*Depth/D90)
                ChezySed = Gravity/CChezy**2

                !Compute friction factor for currents
                Fc   = 0.24*Alog10(12.*Depth/Ksc)**(-2.)

                !Compute the shear stress due to waves and currents
                TauC    = 0.5*RhoWater*Fc*VelMod**2.
                qsi     = sqrt(Fw/Fc)
                TauTot  =(1+0.5*(qsi*Ubw/VelMod)**2.)*TauC 
                
                !Initial mouvement
                Miu  = (CChezy/CLinha)**1.5
                Tinic= -0.27*RhoSandDry*D50*Gravity/(Miu*TauTot)

                !Transport
                Uast   = Sqrt(ChezySed)*VelMod         !m/s
                Ttrans = Coef*D50*Uast             !m^3/(m s)

                !Bed Load
                BedLoad= Ttrans*EXP(Tinic)             !m^3/(m s)

                !Suspended Load (Bhattacharya-TOW)
                AbsTautot = ABS(Tautot)
                if (AbsTautot.GT.1.E-8) then

                    ZX   = VelQ/Kappa/sqrt(Tautot/RhoWater)
                    AAA  = Depth/Ksc
                    B1   = 1.05*(ZX**0.96)/(AAA**(0.013*ZX))
                    B2   = 1-B1
                    if (abs(B2).LT.1.E-5) B2 = 1.E-5
                    KB   = (1.-B1*(0.1667**B2))/B2
                    KS   = 0.415*(B2*((AAA**B2)*alog(30.2*AAA)-3.4078419)+(1.-AAA**B2))/(B2**2.)
                    
                    SuspLoad=BedLoad*KS/KB
                
                else
                
                    SuspLoad=0.
                
                endif

                !Total Load
                TotLoad = BedLoad+SuspLoad
                Me%TransportCapacity(i, j) = TotLoad

            else
                
                Me%TransportCapacity(i, j) = 0.
  
            endif
    enddo
    enddo
    
    if (MonitorPerformance) call StopWatch ("ModuleSand", "BijkerTransport")

    end subroutine BijkerTransport

    !--------------------------------------------------------------------------
    
    subroutine ComputeTauCritic
        !Local-----------------------------------------------------------------
        integer             :: i, j
        real                :: TetaCrit, Exp
        !----------------------------------------------------------------------

        Exp            = 1./3.

        SELECT CASE (Me%TransportMethod)

        Case (MeyerPeter)
!        Case (MeyerPeter .or. Ackers)
            
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                    Me%TauCritic(i, j) = 0.047*Gravity*Me%RhoSl*Me%D50%Field2D(I,J)      ! (N/m2)

                else

                    Me%TauCritic(i, j) = FillValueReal

                endif

            enddo
            enddo

        Case (VanRijn1,VanRijn2,Bijker)
!     Case (VanRijn, Bijker)
            
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                    
                    Me%Dast(I,J) = Me%D50%Field2D(I,J)*(Me%RelativeDensity*Gravity/WaterCinematicVisc**2)**Exp
                    If (Me%Dast(I,J).GT.1.AND.Me%Dast(I,J).LE.4) Then
                       TetaCrit = 0.24/Me%Dast(I,J)
                    ElseIf (Me%Dast(I,J).GT.4.AND.Me%Dast(I,J).LE.10) Then
                       TetaCrit = 0.14*Me%Dast(I,J)**(-0.64)
                    ElseIf (Me%Dast(I,J).GT.10.AND.Me%Dast(I,J).LE.20) Then
                       TetaCrit = 0.04*Me%Dast(I,J)**(-0.1)
                    ElseIf (Me%Dast(I,J).GT.20.AND.Me%Dast(I,J).LE.150) Then
                       TetaCrit = 0.013*Me%Dast(I,J)**0.29
                    ElseIf (Me%Dast(I,J).GT.150) Then
                       TetaCrit = 0.055
                    EndIf

                    Me%TauCritic(i, j)= Me%RhoSl*Gravity*Me%D50%Field2D(I,J)*TetaCrit   ! N/m2

                else

                    Me%TauCritic(i, j)= FillValueReal

                endif

            enddo
            enddo
       
        END SELECT 



    end subroutine ComputeTauCritic
    !--------------------------------------------------------------------------

    subroutine ComputeEvolution            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area1, DT_Area2
        integer                            :: i, j
        !----------------------------------------------------------------------

        Me%DZ%Field2D(:,:) = 0.

        if (Me%Boxes%Yes) then

            Me%Boxes%FluxesX(:,:) = 0.
            Me%Boxes%FluxesY(:,:) = 0.

        endif

        if (Me%Aceleration%Yes) then

            Me%FluxX(:,:) = Me%Aceleration%Coef * Me%FluxX(:,:)
            Me%FluxY(:,:) = Me%Aceleration%Coef * Me%FluxY(:,:)

        endif
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB


            if      (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered .and. Me%FluxX(i, j) < 0.) then

                DT_Area1 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j-1) / Me%ExternalVar%DVY(i, j-1)
                DT_Area2 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                Me%DZ%Field2D(i, j-1) = Me%DZ%Field2D(i, j-1) - DT_Area1 * Me%FluxX(i, j) 
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  ) + DT_Area2 * Me%FluxX(i, j) 

                if (Me%Boxes%Yes) then
                    Me%Boxes%FluxesX(i,j) = Me%FluxX(i, j) * Me%Evolution%SandDT
                endif

            else if (Me%ExternalVar%ComputeFacesU2D(i,j+1) == Covered .and. Me%FluxX(i, j) > 0.) then
                
                DT_Area1 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j+1) / Me%ExternalVar%DVY(i, j+1)
                DT_Area2 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)
          
                Me%DZ%Field2D(i, j+1) = Me%DZ%Field2D(i, j+1) + DT_Area1 * Me%FluxX(i, j)  
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  ) - DT_Area2 * Me%FluxX(i, j) 

                if (Me%Boxes%Yes) then
                    Me%Boxes%FluxesX(i,j+1) = Me%FluxX(i, j) * Me%Evolution%SandDT
                endif

            endif

            if      (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered .and. Me%FluxY(i, j) < 0.) then

                DT_Area1 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i-1, j) / Me%ExternalVar%DVY(i-1, j)
                DT_Area2 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                Me%DZ%Field2D(i-1, j) = Me%DZ%Field2D(i-1, j) - DT_Area1 * Me%FluxY(i, j) 
                Me%DZ%Field2D(i  , j) = Me%DZ%Field2D(i  , j) + DT_Area2 * Me%FluxY(i, j) 

                if (Me%Boxes%Yes) then
                    Me%Boxes%FluxesY(i,j) = Me%FluxY(i, j) * Me%Evolution%SandDT
                endif

            else if (Me%ExternalVar%ComputeFacesV2D(i+1, j) == Covered .and. Me%FluxY(i, j) > 0.) then

                DT_Area1 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i+1, j) / Me%ExternalVar%DVY(i+1, j)
                DT_Area2 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                Me%DZ%Field2D(i+1, j) = Me%DZ%Field2D(i+1, j  ) + DT_Area1 * Me%FluxY(i, j)
                Me%DZ%Field2D(i  , j) = Me%DZ%Field2D(i  , j  ) - DT_Area2 * Me%FluxY(i, j) 

                if (Me%Boxes%Yes) then
                    Me%Boxes%FluxesY(i+1,j) = Me%FluxY(i, j) * Me%Evolution%SandDT
                endif

            endif
                                   
        enddo
        enddo

        call ComputeDischarges
        
    end subroutine ComputeEvolution
    
    !--------------------------------------------------------------------------

    subroutine ComputeDischarges            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area, BottomDensity, TicknessAdd,      &
                                              DischargeFlow, DischargeConc
        integer                            :: i, j, dis, DischargesNumber, STAT_CALL
        !----------------------------------------------------------------------

        
cd1:    if (Me%Discharges%Yes) then

            !Sinks and Sources
            call GetDischargesNumber(Me%ObjDischarges, DischargesNumber, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ComputeDischarges - ModuleSand - ERR10'


            !For all Discharges
d1:         do dis = 1, DischargesNumber
            

                call GetDischargesGridLocalization(Me%ObjDischarges,                    &
                                                   dis, Igrid = I, JGrid = J,           &
                                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ComputeDischarges - ModuleSand - ERR20'

                call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                           Me%ExternalVar%Now,                          &
                                           dis, -99., DischargeFlow,                    &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ComputeDischarges - ModuleSand - ERR30'



                call GetDischargeConcentration (Me%ObjDischarges,                       &
                                                Me%ExternalVar%Now,                     &
                                                dis, DischargeConc,                     &
                                                PropertyIDNumber = Sand_,               &        
                                                STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_)                                                &
                    stop 'ComputeDischarges - ModuleSand - ERR40'
                

                BottomDensity       = Me%Porosity * Me%ExternalVar%WaterDensity + (1. - Me%Porosity) * Me%Density

                DT_Area             = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                                     ![m3/s] * [g/m3/1000] / [kg/m3] * [s/m2] * []
                TicknessAdd         = DischargeFlow * (DischargeConc/1000.) / BottomDensity * DT_Area

                Me%DZ%Field2D(i, j) = Me%DZ%Field2D(i, j  ) + TicknessAdd

            enddo d1
    
    
        endif cd1

    end subroutine ComputeDischarges            
    

!--------------------------------------------------------------------------

    subroutine ComputeSmoothSlope            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area1, DT_Area2, dhdx, dhdy, FluxX, FluxY, SandThickness
        integer                            :: i, j
        !----------------------------------------------------------------------

        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            dhdx = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i, j-1)) /  &
                    Me%ExternalVar%DZX(i,j-1)

            if    ( Me%ExternalVar%WaterPoints2D(i,  j  ) == WaterPoint  .and. &
                    Me%ExternalVar%WaterPoints2D(i,  j-1) == WaterPoint  .and. &
                    abs(dhdx) > Me%SmoothSlope%Critic) then
                
                DT_Area1 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j-1) / Me%ExternalVar%DVY(i, j-1)
                DT_Area2 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                FluxX = Me%SmoothSlope%Factor * max(abs(Me%FluxY(i, j)), abs(Me%FluxY(i, j-1)))

                if (dhdx <  0.) then
                    FluxX = - 1. * FluxX
                endif

                if (dhdx >  0.) then
                    SandThickness = Me%BedRock%Field2D(i, j-1) + Me%DZ_Residual%Field2D(i, j-1)
                else
                    SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
                endif
                    
                if (SandThickness < Me%SandMin) FluxX = 0.

                Me%DZ%Field2D(i, j-1) = Me%DZ%Field2D(i, j-1) - DT_Area1 * FluxX
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  ) + DT_Area2 * FluxX
            endif

            dhdy = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i-1, j)) /  &
                    Me%ExternalVar%DZY(i-1,j)

            if    ( Me%ExternalVar%WaterPoints2D(i,   j) == WaterPoint  .and. &
                    Me%ExternalVar%WaterPoints2D(i-1, j) == WaterPoint  .and. &
                    abs(dhdy) > Me%SmoothSlope%Critic) then
                
                DT_Area1 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i-1, j) / Me%ExternalVar%DVY(i-1, j)
                DT_Area2 = Me%Evolution%SandDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                FluxY = Me%SmoothSlope%Factor * max(abs(Me%FluxX(i, j)), abs(Me%FluxX(i-1, j)))

                if (dhdy <  0.) then
                    FluxY = - 1. * FluxY
                endif

                if (dhdy >  0.) then
                    SandThickness = Me%BedRock%Field2D(i-1, j) + Me%DZ_Residual%Field2D(i-1, j)
                else
                    SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
                endif
                    
                if (SandThickness < Me%SandMin) FluxY = 0.

                Me%DZ%Field2D(i-1, j) = Me%DZ%Field2D(i-1, j) - DT_Area1 * FluxY
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  ) + DT_Area2 * FluxY

            endif

                                   
        enddo
        enddo
        
    end subroutine ComputeSmoothSlope
    
    !--------------------------------------------------------------------------


    subroutine BoundaryCondition(Field2D)
        !Arguments-------------------------------------------------------------
        real,   dimension(:,:), pointer    :: Field2D

        !Local-----------------------------------------------------------------
        real                               :: a1, a2, a3, a4, atotal
        integer                            :: i, j, ILB, IUB, JLB, JUB
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 

        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 


        if          (Me%Boundary == NullGradient) then

            do j=JLB, JUB
            do i=ILB, IUB

                if (Me%ExternalVar%BoundaryPoints2D(i, j) == Boundary) then
                
                    a1 = 0; a2 = 0; a3 = 0; a4 = 0

                    if (Me%ExternalVar%BoundaryPoints2D(i-1, j  ) == Not_Boundary .and.     &
                        Me%ExternalVar%WaterPoints2D   (i-1, j  ) == WaterPoint) a1 = 1

                    if (Me%ExternalVar%BoundaryPoints2D(i+1, j  ) == Not_Boundary .and.     &
                        Me%ExternalVar%WaterPoints2D   (i+1, j  ) == WaterPoint) a2 = 1

                    if (Me%ExternalVar%BoundaryPoints2D(i  , j-1) == Not_Boundary .and.     &
                        Me%ExternalVar%WaterPoints2D   (i  , j-1) == WaterPoint) a3 = 1

                    if (Me%ExternalVar%BoundaryPoints2D(i  , j+1) == Not_Boundary .and.     &
                        Me%ExternalVar%WaterPoints2D   (i  , j+1) == WaterPoint) a4 = 1

                    atotal = (a1 + a2 + a3 + a4)

                    if (atotal > 0) then
                
                        Field2D(i, j) = (a1 * Field2D(i-1, j) + a2 * Field2D(i+1, j)  +     &
                                         a3 * Field2D(i, j-1) + a4 * Field2D(i, j+1)) / atotal
                    endif
                                      

                endif

            enddo
            enddo    

        else if (Me%Boundary == Cyclic) then

            where (Me%ExternalVar%BoundaryPoints2D(ILB, :) == Boundary) 
                Field2D(ILB, :) = Field2D(IUB-1, :)
            end where

            where (Me%ExternalVar%BoundaryPoints2D(IUB, :) == Boundary) 
                Field2D(IUB, :) = Field2D(ILB+1, :)
            end where

            where (Me%ExternalVar%BoundaryPoints2D(: ,JLB) == Boundary) 
                Field2D(:, JLB) = Field2D(:, JUB-1)
            end where

            where (Me%ExternalVar%BoundaryPoints2D(:, JUB) == Boundary) 
                Field2D(:, JUB) = Field2D(:, JLB+1)
            end where

        endif


    end subroutine BoundaryCondition
    
    !--------------------------------------------------------------------------

    subroutine ComputeResidualEvolution            

        !Local-----------------------------------------------------------------
        integer                            :: i, j
        !----------------------------------------------------------------------
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                
                Me%DZ_Residual%Field2D(i, j)    = Me%DZ_Residual%Field2D(i, j) + Me%DZ%Field2D(i, j)

            endif

        enddo
        enddo    

        if (Me%Evolution%BathymDT > Me%Evolution%SandDT) then

            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                
                    Me%BatimIncrement%Field2D(i,j)  = Me%BatimIncrement%Field2D(i,j) + Me%DZ%Field2D(i, j)

                endif

            enddo
            enddo    

        endif


    end subroutine ComputeResidualEvolution
    
    !--------------------------------------------------------------------------



    subroutine OutPutSandHDF
        
        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
         
        !Local-----------------------------------------------------------------
        logical                            :: FirstTime
        integer                            :: OutPutNumber
        type (T_Time)                      :: Actual
        integer                            :: ILB, IUB, JLB, JUB
        real,    dimension(6    ), target  :: AuxTime
        real,    dimension(:    ), pointer :: TimePtr
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB

        !----------------------------------------------------------------------


        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 

        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        !Saida das diferentes propriedades
        Actual = Me%ExternalVar%Now

        FirstTime = .true.        

        OutPutNumber = Me%OutPut%NextOutPut

T1:     if (size(Me%OutPut%OutTime) >= OutPutNumber) then

TOut:       if (Actual >= Me%OutPut%OutTime(OutPutNumber)) then
            
                !Writes current time
                call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),              &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR01'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",     &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR02'

                !Writes OpenPoints
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                  &
                                     WorkJUB, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR03'

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",         &
                                     "-", Array2D = Me%ExternalVar%OpenPoints2D,             &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR04'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Bathymetry", "Bathymetry",        &
                                     "-", Array2D = Me%ExternalVar%Bathymetry,               &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR05'
       
                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%DZ%ID%Name), trim(Me%DZ%ID%Name),  &
                                     trim(Me%DZ%ID%Units), Array2D = Me%DZ%Field2D,                           &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR06'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%DZ_Residual%ID%Name), trim(Me%DZ_Residual%ID%Name),  &
                                     trim(Me%DZ_Residual%ID%Units), Array2D = Me%DZ_Residual%Field2D,   &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR07'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Capacity", "Transport  Capacity",  &
                                     "m3/s/m", Array2D = Me%TransportCapacity,               &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR07'


                call HDF5WriteData  (Me%ObjHDF5, "/Results/Tau Critic", "Tau Critic",        &
                                     "m3/s/m", Array2D = Me%TauCritic,                       &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR07'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux X", "Transport Flux X",        &
                                     "m3/s/m", Array2D = Me%FluxX,                       &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR80'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux Y", "Transport Flux Y",        &
                                     "m3/s/m", Array2D = Me%FluxY,                       &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR90'

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR08'

                Me%OutPut%NextOutPut = OutPutNumber + 1

            endif  TOut    

        endif T1


!        if (MonitorPerformance) call StopWatch ("ModuleSand", "OutPutSandHDF")


    end subroutine OutPutSandHDF

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------

        !if (MonitorPerformance) call StartWatch ("ModuleSand", "OutPut_TimeSeries")

        call WriteTimeSerie(Me%ObjTimeSerie,                     &
                            Data2D = Me%DZ%Field2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'OutPut_TimeSeries - ModuleSand - ERR01'

        call WriteTimeSerie(Me%ObjTimeSerie,                     &
                            Data2D = Me%DZ_Residual%Field2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'OutPut_TimeSeries - ModuleSand - ERR01'


!        if (MonitorPerformance) call StopWatch ("ModuleSand", "OutPut_TimeSeries")

    
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------


    subroutine OutputBoxFluxes


        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_CALL

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        Me%Boxes%Mass(:,:) = 0.

        do J = JLB, JUB
        do I = ILB, IUB
            Me%Boxes%Mass   (i,j) = Me%DZ%Field2D(i,j)      * (1.- Me%Porosity) * Me%Density * &
                                    Me%ExternalVar%DUX(i,j) * Me%ExternalVar%DVY(i,j) 

            !This fluxes are initialised and partial computed in the subroutine ComputeEvolution
            Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j)   * (1.- Me%Porosity) * Me%Density
            Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j)   * (1.- Me%Porosity) * Me%Density

        end do
        end do
                
        !Integration of the bottom changes
        call BoxDif(Me%ObjBoxDif, Me%Boxes%Mass,                                        &
                    trim(GetPropertyName (Sand_)),                                      &
                    Me%ExternalVar%WaterPoints2D,                                       &
                    STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'OutputBoxFluxes - ModuleSand - ERR01'

        !Integration of fluxes
        call BoxDif(Me%ObjBoxDif,                                                       &
                    Me%Boxes%FluxesX,                                                   &
                    Me%Boxes%FluxesY,                                                   &
                    trim(GetPropertyName (Sand_)),                                      &
                    Me%ExternalVar%WaterPoints2D,                                       &
                    STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'OutputBoxFluxes - ModuleSand - ERR20'

    end subroutine OutputBoxFluxes

    !--------------------------------------------------------------------------
    
    real function FallVel(D50)  !Compute Particle Fall Velocity 
        
        !Arguments-------------------------------------------------------------

        real, intent(IN) :: D50     
        
        !Local-----------------------------------------------------------------

        real :: Dsig,VQ1,VQ2

        !----------------------------------------------------------------------

                Dsig = 0.9*D50
                If(Dsig.GT.1.E-6.AND.Dsig.LT.100.E-6)    then
                    FallVel= (Me%RelativeDensity*Gravity*Dsig**2)/(18*WaterCinematicVisc)
                elseif(Dsig.GE.100.E-6.AND.Dsig.LT.1000.E-6) then
                    VQ1     = 10.*WaterCinematicVisc/D50
                    VQ2     = 1.+0.01*Me%RelativeDensity*Gravity*Dsig**3./WaterCinematicVisc**2.
                    FallVel= VQ1*(Sqrt(VQ2)-1)         
                elseif(Dsig.GE.1000E-6)                     then
                    FallVel= 1.1*sqrt(Me%RelativeDensity*Gravity*Dsig)
                endif

                 
        !----------------------------------------------------------------------

    end function FallVel

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillSand(ObjSandID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjSandID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, STAT_CALL    

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, i

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSandID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSand_,  Me%InstanceID)

            if (nUsers == 0) then
            
                 call ReadLockExternalVar

                if (Me%OutPut%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSand - ModuleSand - ERR10'
                endif

                !Kills the TimeSerie
                if (Me%TimeSerie) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSand - ModuleSand - ERR20'
                endif


                 call WriteFinalState

do1 :           do i=1, Me%Classes%Number
                    
                    deallocate(Me%Classes%Diameter(i)%Field2D)
                    deallocate(Me%Classes%Percentage(i)%Field2D)

                enddo do1


if1:            if (Me%Classes%Number > 0) then

                    deallocate(Me%Classes%Diameter)
                    deallocate(Me%Classes%Percentage)
                    deallocate(Me%Classes%Name)

                endif if1

                deallocate(Me%D35%Field2D)
                deallocate(Me%D50%Field2D)
                deallocate(Me%D90%Field2D)

                deallocate(Me%BedRock%Field2D)
                deallocate(Me%TauCritic)

                deallocate(Me%DZ%Field2D)
                
                if (Me%Evolution%BathymDT > Me%Evolution%SandDT) then
                
                    deallocate(Me%BatimIncrement%Field2D)

                endif

                deallocate(Me%DZ_Residual%Field2D)

                if (Me%Filter%ON) deallocate (Me%Filter%Field2D)

                !deallocate fluxes
                deallocate(Me%FluxX)
                deallocate(Me%FluxY)
                deallocate(Me%TransportCapacity)

                !if (Me%TransportMethod == VanRijn1) deallocate (Me%Dast)
                
                if   (Me%TransportMethod == VanRijn1  &
                 .OR. Me%TransportMethod == VanRijn2  &
                 .OR. Me%TransportMethod == Bijker ) deallocate (Me%Dast)

                if (Me%Boxes%Yes) then

                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSand - ModuleSand - ERR30'

                    deallocate(Me%Boxes%FluxesX)
                    nullify   (Me%Boxes%FluxesX)

                    deallocate(Me%Boxes%FluxesY)
                    nullify   (Me%Boxes%FluxesY)


                    deallocate(Me%Boxes%Mass   )
                    nullify   (Me%Boxes%Mass   )

                endif


                call ReadUnLockExternalVar

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillSand - ModuleSand - ERR40'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjBathym)
                if (nUsers == 0) stop 'KillSand - ModuleSand - ERR50'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillSand - ModuleSand - ERR60'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillSand - ModuleSand - ERR70'

#ifndef _WAVES_
                if(Me%ObjWaves /= 0)then
                    nUsers = DeassociateInstance (mWAVES_,Me%ObjWaves)
                    if (nUsers == 0) stop 'KillSand - ModuleSand - ERR80'
                end if
#endif

                if (Me%Discharges%Yes)    then
                    nUsers = DeassociateInstance(mDISCHARGES_,    Me%ObjDischarges)
                    if (nUsers == 0) stop 'KillSand - ModuleSand - ERR90'
                end if


                !Deallocates Instance
                call DeallocateInstance ()

                ObjSandID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillSand
        

    !------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !   Write the final water properties results in HDF format  !
  
    subroutine WriteFinalState

        !Local--------------------------------------------------------------
        integer                                 :: WorkILB, WorkIUB
        integer                                 :: WorkJLB, WorkJUB
        integer                                 :: STAT_CALL
        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        
        call Open_HDF5_OutPut_File(Me%Files%FinalSand)

       
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                               &
                             WorkJLB, WorkJUB,                                           &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSand - ERR01'

        !Final concentration
        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                     &
                             trim(Me%DZ_Residual%ID%Name),                               &
                             trim(Me%DZ_Residual%ID%Units),                              &
                             Array2D = Me%DZ_Residual%Field2D,                           &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSand - ERR02'

        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                     &
                             "Bathymetry",                                               &
                             "m",                                                        &
                             Array2D = Me%ExternalVar%Bathymetry,                        &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSand - ERR03'


   
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSand - ERR03'

        call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSand - ERR04'

    end subroutine WriteFinalState

    !--------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Sand), pointer          :: AuxObjSand
        type (T_Sand), pointer          :: PreviousObjSand

        !Updates pointers
        if (Me%InstanceID == FirstObjSand%InstanceID) then
            FirstObjSand => FirstObjSand%Next
        else
            PreviousObjSand => FirstObjSand
            AuxObjSand      => FirstObjSand%Next
            do while (AuxObjSand%InstanceID /= Me%InstanceID)
                PreviousObjSand => AuxObjSand
                AuxObjSand      => AuxObjSand%Next
            enddo

            !Now update linked list
            PreviousObjSand%Next => AuxObjSand%Next

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

    subroutine Ready (ObjSand_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSand_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjSand_ID > 0) then
            call LocateObjSand (ObjSand_ID)
            ready_ = VerifyReadLock (mSand_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjSand (ObjSandID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSandID

        !Local-----------------------------------------------------------------

        Me => FirstObjSand
        do while (associated (Me))
            if (Me%InstanceID == ObjSandID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleSand - LocateObjSand - ERR01'

    end subroutine LocateObjSand

    !--------------------------------------------------------------------------


    subroutine ReadLockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR01'

        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR02'

        !OpenPoints2D
        call GetOpenPoints2D (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR03'

        !BoundaryPoints2D
        call GetBoundaries(Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR05'

        !Compute faces 2D
        call GetComputeFaces2D(Me%ObjHorizontalMap,                                      &
                               ComputeFaces2DU = Me%ExternalVar%ComputeFacesU2D,         &
                               ComputeFaces2DV = Me%ExternalVar%ComputeFacesV2D,         &
                               ActualTime      = Me%ExternalVar%Now,                     &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR06'


        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                     &
                               DUX  = Me%ExternalVar%DUX,                                &
                               DVY  = Me%ExternalVar%DVY,                                &
                               DZX  = Me%ExternalVar%DZX,                                &
                               DZY  = Me%ExternalVar%DZY,                                &
                               DXX  = Me%ExternalVar%DXX,                                &
                               DYY  = Me%ExternalVar%DYY,                                &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR07'

        
        call GetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR08'

#ifndef _WAVES_
        if (Me%ObjWaves /=0) then        

            call GetWaves (WavesID       = Me%ObjWaves,                                  &
                           WavePeriod    = Me%ExternalVar%WavePeriod,                    &
                           WaveHeight    = Me%ExternalVar%WaveHeight,                    &
                           Abw           = Me%ExternalVar%Abw,                           &
                           Ubw           = Me%ExternalVar%Ubw,                           &
                           WaveDirection = Me%ExternalVar%WaveDirection,                 &
                           STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSand - ERR17'

        endif
#endif

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------


    subroutine ReadUnLockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR01'

        !OpenPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR02'

        !DXX
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid,                                  &
                                  Me%ExternalVar%DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR03'

        !DYY
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR04'


        !BoundaryPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR05'

        !Compute faces 2D V
        call UnGetHorizontalMap(Me%ObjHorizontalMap,                                     &
                               Me%ExternalVar%ComputeFacesV2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR06'

        !Compute faces 2D U
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%ComputeFacesU2D,     &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR07'


        !DUX
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DUX,               &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR08'

        !DVY
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                   &
                               Me%ExternalVar%DVY,                                       &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR09'


        !DZX
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DZX,               &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR88'

        !DZY
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                   &
                               Me%ExternalVar%DZY,                                       &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR89'


        !Bathymetry
        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR10'

#ifndef _WAVES_

        if (Me%ObjWaves /=0) then

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WavePeriod, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR18'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WaveHeight, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSand - ERR19'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%Abw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleSand - ERR20'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%Ubw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleSand - ERR21'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WaveDirection, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnlockExternalWater - ModuleSand - ERR22'


        endif
#endif

    end subroutine ReadUnLockExternalVar
    !--------------------------------------------------------------------------

end module ModuleSand

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------











