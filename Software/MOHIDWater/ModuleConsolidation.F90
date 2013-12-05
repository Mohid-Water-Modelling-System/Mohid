!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Consolidation
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : November 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute sediment consolidation and interstitial water fluxes
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
!
!   CONTINUOUS              : 0/1               [0]             !Speficies if initialization is based in previous run 
!   MIN_THICKNESS           : real            [1e-3]            !Minimum layer thickness 
!   MAX_THICKNESS           : real            [1e-2]            !Maximum layer thickness
!   CONSOLIDATION           : 0/1               [0]             !Specifies if consolidation is to be computed
!   CONSOLIDATION_DT        : sec.           [ModelDT]          !Time step for consolidation
!   DECAYMENT               : 0/1               [0]             !Layers consolidation
!   DECAYTIME               : sec.           [2160000]          !Decay factor for consolidation
!   COMPUTE_SHEAR_STRESS    : 0/1               [0]             !Compute shear stress or read from file
!   INFINITE_CSE            : real             [5 Pa]           !Maximum critical shear stress for erosion
!   SURFACE_CSE             : real           [0.4 Pa]           !Critical shear stress for erosion for the top layer
!   CSE_COEF                : real            [0.05m]           !Coeficient to compute exponential increase of 
!                                                               !critical shear stress for erosion with depth
!
!   <begin_stationary_porosity>
!   See module FillMatrix   : -                 -               !Initialization of porosity values
!   <end_stationary_porosity>  
!
!   <begin_critical_shear_stress>
!   See module FillMatrix   : -                 -               !Initialization of porosity values
!   <end_critical_shear_stress>  
!
!   OUTPUT_TIME             : sec. sec. sec.    []              !Output Time
!   TIME_SERIE_LOCATION     : char              []              !Path to time serie location file
!       TIME_SERIE          : 0/1               [0]             !Ouputs results in time series  
!   BOXFLUXES               : char              []              !If specified computes box integration
!                                                               !based on boxes file defined by this keyword
!<begin_porosity>
!   See module FillMatrix   : -                 -               !Initialization of porosity values
!<end_porosity>
!

Module ModuleConsolidation

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData,        only: ConstructEnterData, KillEnterData, GetData, GetOutPutTime,&
                                      RewindBlock, RewindBuffer, ExtractBlockFromBlock,         &
                                      ExtractBlockFromBuffer, ReadFileName, Block_Unlock
    use ModuleGridData,         only: GetGridData, UngetGridData            
    use ModuleHorizontalMap
    use ModuleHorizontalGrid,   only: GetHorizontalGrid, GetGridCellArea, UnGetHorizontalGrid,  &
                                      WriteHorizontalGrid, GetXYCellZ, GetDomainDecompositionMPI_ID,&
                                      GetDomainDecompositionON
    use ModuleGeometry,         only: GetGeometrySize, UnGetGeometry, ComputeInitialGeometry,   &
                                      GetGeometryVolumes, ReadGeometry, ComputeVerticalGeometry,& 
                                      WriteGeometry, GetGeometryAreas, GetGeometryDistances,    &
                                      GetGeometryKTop       
    use ModuleMap,              only: GetWaterPoints3D, GetOpenPoints3D, UngetMap,              & 
                                      GetComputeFaces3D, UpdateComputeFaces3D          
    use ModuleBoxDif,           only: StartBoxDif, BoxDif, KillBoxDif      
    use ModuleHDF5,             only: ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,          &
                                      HDF5WriteData, HDF5FlushMemory, HDF5ReadData, KillHDF5
    use ModuleTimeSerie,        only: StartTimeSerie, WriteTimeSerie, KillTimeSerie,    &
                                      GetTimeSerieLocation, CorrectsCellsTimeSerie,     &
                                      GetNumberOfTimeSeries, TryIgnoreTimeSerie
    use ModuleFillMatrix,       only: ConstructFillMatrix, GetDefaultValue, KillFillMatrix
    use ModuleStopWatch,        only: StartWatch, StopWatch

    implicit none

    private 

    !Subroutines-----------------------------------------------------------------

    !Constructor
    public  :: ConstructConsolidation
    private ::      AllocateInstance
    private ::      SetTimeLimits
    private ::          ReadConsolidationFilesName
    private ::          DataFileConnection
    private ::      Construct_GlobalVariables
    private ::          AllocateVariables
    private ::          Construct_Parameters
    private ::              ConstructScalar3D
    private ::              Read_Old_Parameters
    private ::      Construct_Initial_Geometry
    private ::      Construct_Outputs
    private ::          Construct_Output_Options
    private ::          Construct_Time_Serie
    private ::          Construct_Box_Time_Serie
    private ::          Open_HDF_OutPut_File

    !Selector
    public  :: UngetConsolidation
    public  :: GetConsolidationWaterVolume
    public  :: GetConsolidationDrySedVolume
    public  :: GetConsolidationWaterFluxes
    public  :: GetConsolidationVelocityZ
    public  :: GetConsolidationVelocityXY
    public  :: GetConsolidationTortuosity
    public  :: GetConsolidationTimeStep
    public  :: GetConsolidationMinThickness
    public  :: GetConsolidationWaterPercentage
    public  :: GetConsolidationDepth
    public  :: GetConsolidationOptions
    public  :: GetConsolidationKTopState
    public  :: GetConsolidationPorosity
    public  :: GetSedimentColumnFull
    public  :: GetConsolidationCriticalShear

    public  :: SetConsolidationFlux
    public  :: SetSedimentDryDensity

    !Modifier
    public  :: ModifyConsolidation
    private ::      ComputeConsolidation
    private ::          ComputeVerticalCoordinate
    private ::          ComputeVelocity
    private ::          ComputeWaterVelocity 
    private ::          New_Geometry    
    private ::          ComputePorosity
    private ::          ComputeTortuosity
    private ::          ComputeVolumes
    private ::          ComputeWaterFluxes
    private ::          ComputeShearStrength
    private ::          ComputeCellCenterDepth
    private ::          ConsolidationOutPut   
    private ::              Write_HDF_Format
    private ::              OutPut_TimeSeries
    private ::              Output_BoxTimeSeries

    !Destructor
    public  :: KillConsolidation
    private ::     DeallocateInstance                                                    
    private ::     Write_Final_Consolidation_HDF

    !Management
    private ::      Ready
    private ::          LocateObjConsolidation

    private ::      ReadLockExternalModules
    private ::      ReadUnLockExternalModules

    !Interfaces------------------------------------------------------------------
    private :: UngetConsolidation2D_I
    private :: UngetConsolidation2D
    private :: UngetConsolidation3D
    private :: UngetConsolidation3Dreal8
    private :: UngetConsolidation2Dreal8
    interface  UngetConsolidation
        module procedure UngetConsolidation2D_I
        module procedure UngetConsolidation2D
        module procedure UngetConsolidation3D
        module procedure UngetConsolidation3Dreal8
        module procedure UngetConsolidation2Dreal8
    end interface UngetConsolidation

    !Parameter-------------------------------------------------------------------
    character(len=20), parameter                :: Char_SedimentColumnFull    = 'Sediment Column Full'
    character(len=8 ), parameter                :: Char_Porosity              = 'Porosity'
    real, parameter                             :: WaterDensity               = 1000 ![kg]/[m3]

    !Types-----------------------------------------------------------------------
                                                
    type  T_Property_3D
         type(T_PropertyID)                     :: ID
         real,    pointer, dimension (:,:,:)    :: Field          => null()
         real                                   :: Scalar         = null_real
         logical                                :: Old            = .false.
    end type   T_Property_3D                      
                                                
    type       T_Files
         character(len=StringLength)            :: Initial        = null_str
         character(len=StringLength)            :: Final          = null_str
         character(len=StringLength)            :: OutPutFields   = null_str
         character(len=StringLength)            :: ConstructData  = null_str
         character(len=StringLength)            :: BoxesFileName  = null_str
         character(len=StringLength)            :: DomainFile     = null_str
    end type T_Files

    type       T_OutPut
         logical                                :: True         = .false.
         type (T_Time), dimension(:), pointer   :: OutTime
         integer                                :: NextOutPut   = null_int
         integer                                :: Number       = null_int
    end type   T_OutPut

    type       T_External
        !ObjTime
        type(T_Time)                            :: Now

        !ObjMap
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D   => null()
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D    => null()
        integer, pointer, dimension(:,:,:)      :: ComputeFacesU3D => null()
        integer, pointer, dimension(:,:,:)      :: ComputeFacesV3D => null()
        integer, pointer, dimension(:,:,:)      :: ComputeFacesW3D => null()

        !ObjGeometry
        real,    pointer, dimension(:,:,:)      :: DWZ             => null()
        real,    pointer, dimension(:,:  )      :: GridCellArea    => null()
        real,    pointer, dimension(:,:,:)      :: SZZ             => null()
        real(8), pointer, dimension(:,:,:)      :: Volume          => null()
        real(8), pointer, dimension(:,:,:)      :: VolumeOld       => null()
        integer, pointer, dimension(:,:  )      :: KTop            => null()
        real,    pointer, dimension(:,:,:)      :: AreaU           => null()
        real,    pointer, dimension(:,:,:)      :: AreaV           => null()

        !ObjSedimentProperties
        real,    pointer, dimension(:,:,:)      :: SedimentDryDensity  => null()

        !ObjInterfaceSedimentWater
        real,    pointer, dimension(:,:  )      :: ConsolidationFlux   => null()

    end type T_External

    type       T_WaterFluxes
        real(8), dimension(:,:,:), pointer      :: X  => null()
        real(8), dimension(:,:,:), pointer      :: Y  => null() 
        real(8), dimension(:,:,:), pointer      :: Z  => null()
    end type   T_WaterFluxes


    type      T_Consolidation
        private
        integer                                 :: InstanceID           = null_int
        type(T_Size3D        )                  :: Size   
        type(T_Size3D        )                  :: WorkSize   
        type(T_Files         )                  :: Files
        type(T_WaterFluxes   )                  :: WaterFluxes  
        type(T_External      )                  :: ExternalVar
        type(T_OutPut        )                  :: OutPut
        type(T_Property_3D)                     :: Porosity
        type(T_Property_3D)                     :: StationaryPorosity
        type(T_Property_3D)                     :: CriticalShearStress
        real(8), pointer, dimension(:,:,:)      :: WaterVolume               => null()
        real(8), pointer, dimension(:,:,:)      :: WaterVolumeOld            => null()
        real(8), pointer, dimension(:,:,:)      :: DrySedimentVolume         => null()
        real(8), pointer, dimension(:,:,:)      :: DrySedimentVolumeOld      => null()
        real,    pointer, dimension(:,:,:)      :: VelocityU                 => null()
        real,    pointer, dimension(:,:,:)      :: WaterVelocityU            => null()
        real,    pointer, dimension(:,:,:)      :: VelocityV                 => null()
        real,    pointer, dimension(:,:,:)      :: WaterVelocityV            => null()
        real,    pointer, dimension(:,:,:)      :: VelocityW                 => null()
        real,    pointer, dimension(:,:,:)      :: WaterVelocityW            => null()
        real,    pointer, dimension(:,:,:)      :: DrySedimentHeight         => null()
        real,    pointer, dimension(:,:,:)      :: WaterPercentage           => null()
        real,    pointer, dimension(:,:,:)      :: CellCenterDepth           => null()
        real,    pointer, dimension(:,:,:)      :: Tortuosity                => null()
        real,    pointer, dimension(:,:,:)      :: VerticalCoordinate        => null()
        real,    pointer, dimension(:,:  )      :: Elevation                 => null()
        integer, pointer, dimension(:,:  )      :: KTop                      => null()
        integer, pointer, dimension(:,:  )      :: KTopState                 => null()
        real,    pointer, dimension(:,:  )      :: TopCriticalShear          => null()
        integer, pointer, dimension(:,:  )      :: SedimentColumnFull        => null()

        real                                    :: MinLayerThickness    = null_real
        real                                    :: MaxLayerThickness    = null_real
        real                                    :: DecayTime            = null_real
        real                                    :: DT                   = null_real

        real                                    :: Surface_CSE          = null_real
        real                                    :: Infinite_CSE         = null_real
        real                                    :: CSE_Coef             = null_real

        type(T_Time)                            :: NextCompute     
        type(T_Time)                            :: BeginTime
        type(T_Time)                            :: EndTime
                                     
        logical                                 :: TimeSerie            = .false.
        logical                                 :: BoxTimeSerie         = .false.

        logical                                 :: ComputeShearStress   = .false. 
        logical                                 :: ContinuesCompute     = .false.
        logical                                 :: Consolidation        = ON
        logical                                 :: Decayment            = ON
                                                
        !Instance of ModuleTime                 
        integer                                 :: ObjTime              = 0

        !Instances of Module_EnterData          
        integer                                 :: ObjEnterData         = 0

        !Instance of ModuleBathymetry           
        integer                                 :: ObjGridData          = 0

        !Instance of ModuleHorizontalMap        
        integer                                 :: ObjHorizontalMap     = 0
                                                
        !Instance of ModuleHorizontalGrid       
        integer                                 :: ObjHorizontalGrid    = 0
                                                
        !Instance of ModuleGeometry             
        integer                                 :: ObjGeometry          = 0
                                                
        !Instance of ModuleMap                  
        integer                                 :: ObjMap               = 0

        !Instance of ModuleBoxDif
        integer                                 :: ObjBoxDif            = 0              
                                                
        !Instance of Module HDF5                
        integer                                 :: ObjHDF5              = 0
                                                
        !Instance of ModuleTimeSerie            
        integer                                 :: ObjTimeSerie         = 0
                                                
        type(T_Consolidation), pointer          :: Next
    end type T_Consolidation                    
                                                
                                                
    !Global Module Variables                    
    type (T_Consolidation), pointer             :: FirstObjConsolidation
    type (T_Consolidation), pointer             :: Me


    !----------------------------------------------------------------------------
    
    contains



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructConsolidation(ConsolidationID,                       &
                                      TimeID,                                &
                                      GridDataID,                            & 
                                      HorizontalMapID,                       &
                                      HorizontalGridID,                      &
                                      GeometryID,                            &
                                      MapID,                                 &
                                      STAT)

        !Arguments-------------------------------------------------------------

        integer                                      :: ConsolidationID
        integer                                      :: GridDataID
        integer                                      :: HorizontalGridID
        integer                                      :: TimeID
        integer                                      :: HorizontalMapID
        integer                                      :: GeometryID
        integer                                      :: MapID
        integer, optional, intent(out)               :: STAT
        
        !External--------------------------------------------------------------
        integer                                      :: ready_ 
                       
        !Local-----------------------------------------------------------------
        integer                                      :: STAT_

        !----------------------------------------------------------------------
        
        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mConsolidation_)) then
            nullify (FirstObjConsolidation)
            call RegisterModule (mConsolidation_) 
        endif

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID      )
            
            !REVIEW THIS
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            !REVIEW THIS

            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )
            Me%ObjMap            = AssociateInstance (mMAP_,            MapID           )

            call SetTimeLimits

            call ReadConsolidationFilesName

            call DataFileConnection(ON)

            call Construct_GlobalVariables

            call Construct_Initial_Geometry

            call Construct_Outputs

            call DataFileConnection(OFF)

            call null_time(Me%ExternalVar%Now)

            ConsolidationID = Me%InstanceID

            STAT_ = SUCCESS_
        
        else 
            
            stop 'ConstructConsolidation - ModuleConsolidation - ERR99' 

        end if 

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructConsolidation

    !--------------------------------------------------------------------------
   
    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_Consolidation), pointer           :: NewObjConsolidation
        type (T_Consolidation), pointer           :: PreviousObjConsolidation

        !Allocates new instance
        allocate (NewObjConsolidation)
        nullify  (NewObjConsolidation%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjConsolidation)) then
            FirstObjConsolidation           => NewObjConsolidation
            Me                              => NewObjConsolidation
        else
            PreviousObjConsolidation        => FirstObjConsolidation
            Me                              => FirstObjConsolidation%Next
            do while (associated(Me))
                PreviousObjConsolidation    => Me
                Me                          => Me%Next
            enddo
            Me                              => NewObjConsolidation
            PreviousObjConsolidation%Next   => NewObjConsolidation
        endif

        Me%InstanceID = RegisterNewInstance (mCONSOLIDATION_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    
    subroutine ReadConsolidationFilesName

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        character(len = StringLength)       :: Message

        !----------------------------------------------------------------------
        
        Message = trim('ASCII file used to construct new consolidation parameters.')

        call ReadFileName('CON_DAT',                                                     &
                           Me%Files%ConstructData,                                       & 
                           Message = Message,                                            &
                           STAT    = STAT_CALL)
        if (STAT_CALL/=SUCCESS_)stop 'ReadConsolidationFilesName - ModuleConsolidation - ERR01'

        ! ---> File in HDF format where is written instant fields of sediment properties
        Message   = trim('Instant fields of consolidation parameters in HDF format.')

        call ReadFileName('CON_HDF',                                                     &
                           Me%Files%OutPutFields,                                        &
                           Message   = Message,                                          &
                           Time_End  = Me%EndTime,                                       &
                           MPI_ID = GetDomainDecompositionMPI_ID (Me%ObjHorizontalGrid), &
                           DD_ON  = GetDomainDecompositionON     (Me%ObjHorizontalGrid), &
                           STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_)stop 'ReadConsolidationFilesName - ModuleConsolidation - ERR02'

        Message   = trim('Consolidation final values in HDF5 format.')

        call ReadFileName('CON_FIN',                                                     &
                           Me%Files%Final,                                               &
                           Message = Message,                                            &
                           Time_End  = Me%EndTime,                                       &
                           MPI_ID = GetDomainDecompositionMPI_ID (Me%ObjHorizontalGrid), &
                           DD_ON  = GetDomainDecompositionON     (Me%ObjHorizontalGrid), &
                           STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_)stop 'ReadConsolidationFilesName - ModuleConsolidation - ERR03'


        ! ---> sediment properties initial values in HDF format
        Message   = trim('Consolidation initial values in HDF format.')

        call ReadFileName('CON_INI',                                                     &
                           Me%Files%Initial,                                             &
                           Message   = Message,                                          &
                           Time_End  = Me%EndTime,                                       &
                           MPI_ID = GetDomainDecompositionMPI_ID (Me%ObjHorizontalGrid), &
                           DD_ON  = GetDomainDecompositionON     (Me%ObjHorizontalGrid), &
                           STAT   = STAT_CALL)
                                                                                               
cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_   ) then 

            call SetError(FATAL_, INTERNAL_,                                             & 
                 'Initial file not found - ReadConsolidationFilesName - ModuleConsolidation - ERR04')

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then cd1

            call SetError(WARNING_, KEYWORD_,                                            &
                          'Keyword CON_INI not found - ModuleConsolidation. ERR05 ',     &
                          Screen = .false.)

        else if (STAT_CALL .EQ. SUCCESS_             ) then cd1

            continue
        
        else cd1

            stop 'ReadConsolidationFilesName - ModuleConsolidation - ERR06'

        end if cd1  
                                                                             
        !----------------------------------------------------------------------

    end subroutine ReadConsolidationFilesName

    !--------------------------------------------------------------------------

    subroutine DataFileConnection(Connect)

        !Arguments-------------------------------------------------------------
        logical, intent(in)              :: Connect

        !External--------------------------------------------------------------
        integer                          :: STAT_CALL
        integer                          :: FATAL_
        
        !Begin-----------------------------------------------------------------

        if(Connect)then

            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                call SetError(FATAL_, INTERNAL_, "DataFileConnection - ModuleConsolidation. ERR01")
        else

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                         &
                call SetError(FATAL_, INTERNAL_, "DataFileConnection - ModuleConsolidation. ERR02")
                 
        end if


    end subroutine DataFileConnection

    
    !--------------------------------------------------------------------------


    subroutine SetTimeLimits

        !External--------------------------------------------------------------
        integer                          :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetComputeTimeLimits(Me%ObjTime,                               &
                                  EndTime   = Me%EndTime,                   &
                                  BeginTime = Me%BeginTime,                 &
                                  STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            call SetError(FATAL_, INTERNAL_, "SetTimeLimits - ModuleConsolidation. ERR01") 

        Me%NextCompute = Me%BeginTime

        !Actualizes the time
        call GetComputeCurrentTime(Me%ObjTime,Me%ExternalVar%Now, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            call SetError(FATAL_, INTERNAL_, "SetTimeLimits - ModuleConsolidation. ERR02") 

    end subroutine SetTimeLimits

    
    !----------------------------------------------------------------------

    
    subroutine Construct_GlobalVariables

        !External--------------------------------------------------------------
        integer                          :: STAT_CALL, iflag
        
        !Begin--------------------------------------------------------------

        call ReadLockExternalModules

        call GetData(Me%ContinuesCompute,                                                &
                     Me%ObjEnterData, iflag,                                             & 
                     SearchType   = FromFile,                                            &
                     keyword      ='CONTINUOUS',                                         &
                     Default      = .false.,                                             &
                     ClientModule ='ModuleConsolidation',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_GlobalVariables - ModuleConsolidation - ERR01'

        call GetData(Me%Consolidation,                                                   &
                     Me%ObjEnterData, iflag,                                             & 
                     SearchType   = FromFile,                                            &
                     keyword      ='CONSOLIDATION',                                      &
                     Default      = .true.,                                              &
                     ClientModule ='ModuleConsolidation',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_GlobalVariables - ModuleConsolidation - ERR02'


        call GetData(Me%Decayment,                                                       &
                     Me%ObjEnterData, iflag,                                             & 
                     SearchType   = FromFile,                                            &
                     keyword      ='DECAYMENT',                                          &
                     Default      = .true.,                                              &
                     ClientModule ='ModuleConsolidation',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_GlobalVariables - ModuleConsolidation - ERR02'


        call GetGeometrySize(Me%ObjGeometry, Size = Me%Size, WorkSize = Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Construct_GlobalVariables - ModuleConsolidation - ERR03'

        call AllocateVariables

        call Construct_Parameters

        call ReadUnLockExternalModules


    end subroutine Construct_GlobalVariables

    
    !---------------------------------------------------------------------------

    
    subroutine Construct_Parameters

        !Local--------------------------------------------------------------
        integer                          :: STAT_CALL, iflag
        real                             :: DT
        type(T_Property_3D), pointer     :: Scalar3D
        integer                          :: ILB, IUB, JLB, JUB, KLB, KUB

        !Begin--------------------------------------------------------------


        !----------------------------------------------------------------------
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        Scalar3D => Me%Porosity
        

        if(Me%ContinuesCompute)then
            call Read_Old_Parameters(Scalar3D, Me%SedimentColumnFull,           &
                                     Char_Porosity, Char_SedimentColumnFull) 
        else
            call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,           &
                                   block_begin = "<begin_porosity>",            &
                                   block_end   = "<end_porosity>")
        end if

        call GetComputeTimeStep(Me%ObjTime, DT, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR20'

        call GetData(Me%DT,                                                     &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      = 'CONSOLIDATION_DT',                         &
                     Default      = DT,                                         &
                     ClientModule = 'ModuleConsolidation',                      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR30'

        Me%NextCompute = Me%NextCompute + Me%DT


        call GetData(Me%MinLayerThickness,                                          &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'MIN_THICKNESS',                                &
                     Default      = 0.001,                                          &
                     ClientModule = 'ModuleConsolidation',                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR31'


        call GetData(Me%MaxLayerThickness,                                          &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'MAX_THICKNESS',                                &
                     Default      = 0.01,                                           &
                     ClientModule = 'ModuleConsolidation',                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR31'
        
        if(Me%Decayment)then
            
            call GetData(Me%DecayTime,                                              &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'DECAYTIME',                                &
                         Default      = 2160000.0,                                  &
                         ClientModule = 'ModuleConsolidation',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR40'


            Scalar3D => Me%StationaryPorosity
            call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                                   block_begin = "<begin_stationary_porosity>",     &
                                   block_end   = "<end_stationary_porosity>")

        end if

        call GetData(Me%ComputeShearStress,                                         &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'COMPUTE_SHEAR_STRESS',                         &
                     Default      = .false.,                                        &
                     ClientModule = 'ModuleConsolidation',                          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR50'

        if(Me%ComputeShearStress)then


            call GetData(Me%Surface_CSE,                                            &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'SURFACE_CSE',                              &
                         Default      = 0.4,                                        &
                         ClientModule = 'ModuleConsolidation',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR60'

            call GetData(Me%Infinite_CSE,                                           &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'INFINITE_CSE',                             &
                         Default      = 5.,                                         &
                         ClientModule = 'ModuleConsolidation',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR60'

            call GetData(Me%CSE_Coef,                                               &
                         Me%ObjEnterData, iflag,                                    &
                         SearchType   = FromFile,                                   &
                         keyword      = 'CSE_COEF',                                 &
                         Default      = 0.1,                                        &
                         ClientModule = 'ModuleConsolidation',                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR60'



            allocate(Me%CriticalShearStress%Field(ILB:IUB, JLB:JUB, KLB:KUB))
            Me%CriticalShearStress%Field(:,:,:) = 0.


        else

            Scalar3D => Me%CriticalShearStress
            call ConstructScalar3D(Scalar3D, ExtractType = FromBlock,               &
                                   block_begin = "<begin_critical_shear_stress>",   &
                                   block_end   = "<end_critical_shear_stress>")

        end if

    end subroutine Construct_Parameters

    
    !--------------------------------------------------------------------------


    subroutine ConstructScalar3D(Scalar3D, ExtractType, ClientNumber, block_begin, block_end)

        !Arguments-------------------------------------------------------------
        type(T_Property_3D), pointer        :: Scalar3D
        integer, intent(in)                 :: ExtractType
        integer, intent(in), optional       :: ClientNumber
        character(len=*)                    :: block_begin, block_end

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        logical                             :: BlockFound
        integer                             :: BlockClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        
        !----------------------------------------------------------------------
        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        select case(ExtractType)

            case(FromBlock)

                call ExtractBlockFromBuffer(Me%ObjEnterData, BlockClientNumber, &
                                            block_begin, block_end,             &
                                            BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR01'

            case(FromBlockInBlock)

                if(.not. present(ClientNumber))then
                    stop 'ConstructScalar3D - ModuleConsolidation - ERR02'
                end if
                
                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,       &
                                           block_begin, block_end,              &
                                           BlockFound, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR02'

        end select

        if(BlockFound)then

            allocate(Scalar3D%Field(ILB:IUB, JLB:JUB, KLB:KUB))

            call ConstructFillMatrix  (PropertyID           = Scalar3D%ID,                      &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       GeometryID           = Me%ObjGeometry,                   &
                                       ExtractType          = ExtractType,                      &
                                       PointsToFill3D       = Me%ExternalVar%WaterPoints3D,     &
                                       Matrix3D             = Scalar3D%Field,                   &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR03'


            call GetDefaultValue(Scalar3D%ID%ObjFillMatrix, Scalar3D%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR04'

            call KillFillMatrix(Scalar3D%ID%ObjFillMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR05'
            
            if(ExtractType == FromBlock)then
                call Block_Unlock(Me%ObjEnterData, BlockClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR06'
            end if

        else
            stop 'ConstructScalar3D - ModuleConsolidation - ERR07'
        end if


        select case(ExtractType)

            case(FromBlock)

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR08'

            case(FromBlockInBlock)

                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ConstructScalar3D - ModuleConsolidation - ERR09'

        end select

   
    end subroutine ConstructScalar3D

    !--------------------------------------------------------------------------

    subroutine Read_Old_Parameters(Scalar3D, SedimentColumnFull, Name, NameFull)

        !Arguments--------------------------------------------------------------
        type(T_Property_3D), pointer              :: Scalar3D
        integer,   dimension(:,:), pointer        :: SedimentColumnFull
        character (len=*), intent(in)             :: Name, NameFull

        !Local-----------------------------------------------------------------
        integer                                   :: IUB, JUB, ILB, JLB, KLB, KUB
        integer                                   :: WIUB, WJUB, WILB, WJLB, WKLB, WKUB
        integer                                   :: STAT_CALL
        integer                                   :: HDF5_READ
        integer                                   :: ObjHDF5 = 0
        logical                                   :: Exists = .false.

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB;   WILB = Me%WorkSize%ILB
        JLB = Me%Size%JLB;   WJLB = Me%WorkSize%JLB
        IUB = Me%Size%IUB;   WIUB = Me%WorkSize%IUB
        JUB = Me%Size%JUB;   WJUB = Me%WorkSize%JUB
        KLB = Me%Size%KLB;   WKLB = Me%WorkSize%KLB
        KUB = Me%Size%KUB;   WKUB = Me%WorkSize%KUB
        
        nullify (Scalar3D%Field)
        allocate(Scalar3D%Field(ILB:IUB, JLB:JUB, KLB:KUB))

        ObjHDF5 = 0

        inquire(FILE = trim(Me%Files%Initial)//"5", EXIST = Exists)

        if (.not. Exists) then
            write(*,*)'Could not find initial file to continue computation'
            stop 'Read_Old_Parameters - ModuleConsolidation - ERR10'
        end if

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


        call ConstructHDF5 (ObjHDF5, trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_Old_Parameters - ModuleConsolidation - ERR20'


        call HDF5SetLimits  (ObjHDF5, WILB, WIUB, WJLB, WJUB, WKLB, WKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_Old_Parameters - ModuleConsolidation - ERR30'

            
        call HDF5ReadData(ObjHDF5, "/"//trim(Name)//"/", &
                          Name, Array3D = Scalar3D%Field, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_Old_Parameters - ModuleConsolidation - ERR40'            


        call HDF5ReadData(ObjHDF5, "/"//trim(NameFull)//"/", &
                          NameFull, Array2D = SedimentColumnFull, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_Old_Parameters - ModuleConsolidation - ERR50'            


        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_Old_Parameters - ModuleConsolidation - ERR60'            

    end subroutine Read_Old_Parameters

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
     
    subroutine Construct_Outputs 

        call Construct_Output_Options 

        if(Me%TimeSerie)then

            call Construct_Time_Serie     

        end if

        call Construct_Box_Time_Serie 

        call Open_HDF_OutPut_File     

    end subroutine Construct_Outputs

    
    !----------------------------------------------------------------------------

    
    subroutine AllocateVariables

        !External----------------------------------------------------------------
        integer                         :: STAT_CALL

        !Local-------------------------------------------------------------------
        integer                         :: ILB, IUB
        integer                         :: JLB, JUB
        integer                         :: KLB, KUB
        !------------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB
        JUB = Me%Size%JUB
        KLB = Me%Size%KLB
        KUB = Me%Size%KUB
        
        !VerticalCoordinate
        allocate(Me%VerticalCoordinate(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR00") 
        Me%VerticalCoordinate(:,:,:) = null_real

        allocate(Me%KTop(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR01") 
        Me%KTop(:,:  ) = 0.

        allocate(Me%KTopState(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR01a") 
        Me%KTopState(:,:  ) = 0.


        !dummy elevation
        allocate(Me%Elevation(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR02") 
        Me%Elevation(:,:  ) = 0.


        allocate(Me%TopCriticalShear(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR03") 
        Me%TopCriticalShear(:,:  ) = 0.


        allocate(Me%SedimentColumnFull(ILB:IUB, JLB:JUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR04") 
        Me%SedimentColumnFull(:,:  ) = 0


       
        !dry sediment height
        allocate(Me%DrySedimentHeight(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR04") 
        Me%DrySedimentHeight(:,:,:) = null_real

        !horizontal consolidation velocity V (dummy)
        allocate(Me%VelocityV(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR07") 
        Me%velocityV(:,:,:) = 0.

        !horizontal consolidation velocity U (dummy)
        allocate(Me%VelocityU(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR08")
        Me%velocityU(:,:,:) = 0.
        
        !vertical consolidation velocity W
        allocate(Me%VelocityW(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR09") 
        Me%velocityW(:,:,:) = 0.

        !horizontal water expulsion velocity U (dummy)
        allocate(Me%WaterVelocityU(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR10") 
        Me%WatervelocityU(:,:,:) = 0.
        
        !horizontal water expulsion velocity V (dummy)
        allocate(Me%WaterVelocityV(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR11") 
        Me%WatervelocityV(:,:,:) = 0.
        
        !vertical water expulsion velocity W 
        allocate(Me%WaterVelocityW(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR12") 
        Me%WatervelocityW(:,:,:) = 0.

        !water flux X direction (dummy)
        allocate(Me%WaterFluxes%X(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR13") 
        Me%WaterFluxes%X(:,:,:) = 0.

        !water flux Y direction (dummy)
        allocate(Me%WaterFluxes%Y(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR14") 
        Me%WaterFluxes%Y(:,:,:) = 0.

        !vertical water flux Z direction 
        allocate(Me%WaterFluxes%Z(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR15") 
        Me%WaterFluxes%Z(:,:,:) = 0.
      
        !water volume within a cell
        allocate(Me%WaterVolume    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR17") 
        Me%WaterVolume (:,:,:)      = null_real
        
        !store water volume value 
        allocate(Me%WaterVolumeOld (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR18") 
        Me%WaterVolumeOld (:,:,:)   = null_real
        
        !dry sediment volume within a cell
        allocate(Me%DrySedimentVolume    (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR19") 
        Me%DrySedimentVolume (:,:,:) = null_real
        
        !store dry sediment volume 
        allocate(Me%DrySedimentVolumeOld (ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR20") 
        Me%DrySedimentVolumeOld (:,:,:)   = null_real

        !water percentage 
        allocate(Me%WaterPercentage(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR21") 
        Me%WaterPercentage (:,:,:) = 1.

        !CellCenterDepth
        allocate(Me%CellCenterDepth(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR22") 
        Me%CellCenterDepth (:,:,:) = null_real

        !Tortuosity 
        allocate(Me%Tortuosity(ILB:IUB, JLB:JUB, KLB:KUB), STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                              &
            call SetError(FATAL_, INTERNAL_, "AllocateVariables; ModuleConsolidation. ERR24") 
        Me%Tortuosity (:,:,:) = null_real
        
        

    end subroutine AllocateVariables  

    !--------------------------------------------------------------------------
    
    subroutine Construct_Initial_Geometry 
        
        !Local----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: i, j, k

        !Begin----------------------------------------------------------------

        !Initial sediment thickness equals the one specified in bathymetry 
        !so elevation equals zero
        Me%Elevation(:,:) = 0.

        if(Me%ContinuesCompute)then
            
            call ReadGeometry(Me%ObjGeometry, trim(Me%Files%Initial)//"5", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleConsolidation - ERR01'

        end if

        !Get WaterPoints3D
        call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Initial_Geometry - ModuleConsolidation - ERR02")

        
        !Computes sediment initial geometry
        call ComputeInitialGeometry(Me%ObjGeometry,                                      &
                                    WaterPoints3D    = Me%ExternalVar%WaterPoints3D,     &
                                    SurfaceElevation = Me%Elevation,                     &
                                    ContinuesCompute = Me%ContinuesCompute,              &
                                    ActualTime       = Me%ExternalVar%Now,               &
                                    STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleConsolidation - ERR03'
        
        !Unget WaterPoints3D
        call UnGetMap(Me%ObjMap,Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Initial_Geometry - ModuleConsolidation - ERR04")


        call ReadLockExternalModules

        Me%VerticalCoordinate(:,:,:) = Me%ExternalVar%SZZ(:,:,:)
        Me%KTop(:,:)                 = Me%ExternalVar%KTop(:,:)

        call ComputeTortuosity

        call ComputeCellCenterDepth

        if(Me%ComputeShearStress)then

            call ComputeShearStrength

        endif

        call ComputeTopShear

        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%WaterPoints3D (i,j,k) == WaterPoint) then

                Me%DrySedimentHeight   (i,j,k) = Me%ExternalVar%DWZ(i,j,k) * (1 - Me%Porosity%Field(i,j,k))

                Me%DrySedimentVolumeOld(i,j,k) = (1.0  - Me%Porosity%Field(i,j,k)) * Me%ExternalVar%Volume(i,j,k)

                Me%DrySedimentVolume   (i,j,k) = Me%DrySedimentVolumeOld(i,j,k)

                Me%WaterVolumeOld      (i,j,k) = Me%Porosity%Field(i,j,k) * Me%ExternalVar%Volume(i,j,k) 

                Me%WaterVolume         (i,j,k) = Me%WaterVolumeOld(i,j,k)


            endif
        enddo
        enddo
        enddo

        call ReadUnlockExternalModules

        !Computes OpenPoints3D
        call UpdateComputeFaces3D(Me%ObjMap, STAT = STAT_CALL)      
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleConsolidation - ERR04'


            
    end subroutine Construct_Initial_Geometry
    
    
    !--------------------------------------------------------------------------
    
    
    subroutine Construct_Output_Options 

        !External----------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-------------------------------------------------------------------
        integer                                 :: iflag

        !------------------------------------------------------------------------

        call GetData(Me%TimeSerie,                                          &
                     Me%ObjEnterData,                                       & 
                     flag         = iflag,                                  &
                     SearchType   = FromFile,                               &
                     keyword      ='TIME_SERIE',                            &
                     Default      = ON,                                     &
                     ClientModule ='ModuleConsolidation',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Output_Options - ModuleConsolidation - ERR01'


        call GetOutPutTime(Me%ObjEnterData,                                 &
                           CurrentTime   = Me%BeginTime,                    &
                           EndTime       = Me%EndTime,                      &
                           keyword       = 'OUTPUT_TIME',                   &
                           SearchType    = FromFile,                        &
                           OutPutsTime   = Me%OutPut%OutTime,               &
                           OutPutsOn     = Me%OutPut%True,                  &
                           OutPutsNumber = Me%OutPut%Number,                &
                           STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Output_Options - ModuleConsolidation - ERR02'

        if (Me%OutPut%True) Me%OutPut%NextOutPut = 1
            

    end subroutine Construct_Output_Options


    !--------------------------------------------------------------------------

    subroutine Open_HDF_OutPut_File

        !Local-----------------------------------------------------------------
        real,    pointer, dimension(:, :   )        :: Bathymetry
        real,    pointer, dimension(:, :   )        :: XX_IE, YY_IE
        integer, pointer, dimension(:, :, :)        :: WaterPoints3D
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        !Gets a pointer to Bathymetry
        call GetGridData      (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR01'

        !Gets XX_IE and YY_IE
        call GetHorizontalGrid  (Me%ObjHorizontalGrid,                      &
                                 XX_IE = XX_IE, YY_IE = YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR02'

        !Gets WaterPoints3D
        call GetWaterPoints3D   (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR03'


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5,                                &
                                 trim(Me%Files%OutPutFields)//"5",          &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR05'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR17'


        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",       &
                              Array2D = Bathymetry,                         &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR18'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = WaterPoints3D,                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR19'

        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB+1, WorkJLB,      &
                              WorkJUB+1, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR20'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionX", "m",      &
                              Array2D = XX_IE,                              &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR21'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "ConnectionY", "m",      &
                              Array2D = YY_IE,                              &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR22'


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR23'


        !Ungets the Bathymetry
        call UngetGridData (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR24'

        !Ungets the WaterPoints
        call UnGetMap        (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR25'

        !Ungets XX_IE, YY_IE
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR26'

        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF_OutPut_File - Consolidation - ERR27'

        !----------------------------------------------------------------------

    end subroutine Open_HDF_OutPut_File


    !----------------------------------------------------------------------------


    subroutine Construct_Time_Serie

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber     
        integer                                             :: nProperties
        integer                                             :: STAT_CALL
        integer, dimension(:,:,:), pointer                  :: WaterPoints3D
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile

        !Begin-----------------------------------------------------------------

       

        !Get WaterPoints3D
        call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Time_Serie - ModuleConsolidation - ERR10")
        
        
        WaterPoints3D => Me%ExternalVar%WaterPoints3D 
        
        nProperties = 8

        !Allocates PropertyList
        allocate(PropertyList(nProperties), STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Time_Serie - ModuleConsolidation - ERR20")

        !Fills up PropertyList
        PropertyList(1) =  Char_Porosity
        PropertyList(2) = 'Thickness'
        PropertyList(3) = 'WaterVelocityZ'
        PropertyList(4) = 'WaterVolume'
        PropertyList(5) = 'DrySedimentHeight'
        PropertyList(6) = 'DrySedVol'
        PropertyList(7) = 'DrySedVolOld'
        PropertyList(8) = 'Depth'
        
        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleWaterProperties',                            &
                     Default      = Me%Files%ConstructData,                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'Construct_Time_Serie - ModuleConsolidation - ERR30' 


        !Constructs TimeSerie (file extension *.src)
        call StartTimeSerie(Me%ObjTimeSerie,                                            &
                            Me%ObjTime,                                                 &
                            TimeSerieLocationFile,                                      &
                            PropertyList, "src",                                        &
                            WaterPoints3D,                                              &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, INTERNAL_, "Construct_Time_Serie - ModuleConsolidation - ERR40")


        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Time_Serie - ModuleConsolidation - ERR50")


        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleConsolidation - ERR60'

        do dn = 1, TimeSerieNumber

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleConsolidation - ERR70'

                if (Id < 0 .or. Jd < 0) then
                
                    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleConsolidation - ERR80'

                    if (IgnoreOK) then
                        cycle
                    else
                        stop 'Construct_Time_Serie - ModuleConsolidation - ERR90'
                    endif

                endif

                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_Time_Serie - ModuleConsolidation - ERR100'
            endif


        enddo


        !Unget WaterPoints3D
        call UnGetMap(Me%ObjMap,Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Time_Serie - ModuleConsolidation - ERR110")
    
        
    end subroutine Construct_Time_Serie

    !----------------------------------------------------------------------------

    subroutine Construct_Box_Time_Serie

        !External--------------------------------------------------------------
        integer                                                     :: STAT_CALL
        integer                                                     :: iflag
        logical                                                     :: Exist, Opened
 
        !Local-----------------------------------------------------------------
        character(len = StringLength), dimension(:    ), pointer    :: PropertyList
        integer,                       dimension(:,:,:), pointer    :: WaterPoints3D

        !----------------------------------------------------------------------

        !Get Waterpoints3D
        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Box_Time_Serie - Consolidation -ERR01'


        call GetData(Me%Files%BoxesFileName,                                        &
                     Me%ObjEnterData, iflag,                                        &
                     keyword      = 'BOXFLUXES',                                    &
                     ClientModule ='ModuleConsolidation',                           &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Box_Time_Serie - Consolidation -ERR02'

        if (iflag .EQ. 1) then
            inquire(FILE = Me%Files%BoxesFileName, EXIST = Exist)
            if (exist) then
                
                Me%BoxTimeSerie = ACTIVE

                inquire(FILE = Me%Files%BoxesFileName, OPENED  = opened)
                if (opened)then
                    write(*,*    ) 
                    write(*,'(A)') 'BoxesFileName = ', trim(adjustl(Me%Files%BoxesFileName))
                    write(*,*    ) 'Already opened.'
                    stop 'Construct_Box_Time_Serie - Consolidation -ERR02'
                end if 

                allocate(PropertyList(1:2), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_Box_Time_Serie - Consolidation -ERR03'

                PropertyList(1) = 'Sediment_Mass'
                PropertyList(2) = 'WaterVolume'

                call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                 &                
                                 TimeID             = Me%ObjTime,                   &
                                 HorizontalGridID   = Me%ObjHorizontalGrid,         &              
                                 BoxesFilePath      = Me%Files%BoxesFileName,       &
                                 ScalarOutputList   = PropertyList,                 & 
                                 WaterPoints3D      = WaterPoints3D,                &
                                 Size3D             = Me%Size,                      &
                                 WorkSize3D         = Me%WorkSize,                  &
                                 STAT               = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)stop 'Construct_Box_Time_Serie - Consolidation -ERR04'

                deallocate (PropertyList)
                nullify    (PropertyList)

            end if
        end if

        call UnGetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)stop 'Construct_Box_Time_Serie - Consolidation -ERR05'

        !----------------------------------------------------------------------

    end subroutine Construct_Box_Time_Serie



    !----------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine GetConsolidationPorosity(ConsolidationID, Porosity, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ConsolidationID
        real, dimension(:,:,:),  pointer            :: Porosity
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

            Porosity => Me%Porosity%Field


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine GetConsolidationPorosity


    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine GetSedimentColumnFull(ConsolidationID, SedimentColumnFull, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ConsolidationID
        integer, dimension(:,:),  pointer           :: SedimentColumnFull
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

            SedimentColumnFull => Me%SedimentColumnFull


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_


    end subroutine GetSedimentColumnFull


    !--------------------------------------------------------------------------
    
    subroutine GetConsolidationCriticalShear(ConsolidationID, TopCriticalShear, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ConsolidationID
        real, dimension(:,:  ),  pointer            :: TopCriticalShear
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

            TopCriticalShear => Me%TopCriticalShear


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetConsolidationCriticalShear
    
    
    !--------------------------------------------------------------------------

    
    subroutine GetConsolidationDepth(ConsolidationID, CellCenterDepth, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ConsolidationID
        real, dimension(:,:,:),  pointer            :: CellCenterDepth
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

            CellCenterDepth => Me%CellCenterDepth


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

    end subroutine GetConsolidationDepth
    
    
    !----------------------------------------------------------------------
    
    
    subroutine GetConsolidationWaterVolume(ConsolidationID, WaterVolume, WaterVolumeOld, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ConsolidationID
        real(8), dimension(:,:,:), optional, pointer:: WaterVolume, WaterVolumeOld
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(WaterVolume     )) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                WaterVolume      => Me%WaterVolume
            endif cd2

cd3 :       if (present(WaterVolumeOld  )) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                WaterVolumeOld   => Me%WaterVolumeOld
            endif cd3

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationWaterVolume


    !--------------------------------------------------------------------------

    
    subroutine GetConsolidationDrySedVolume(ConsolidationID, DrySedimentVolume, &
                                            DrySedimentVolumeOld, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                      :: ConsolidationID
        real(8), dimension(:,:,:), optional, pointer :: DrySedimentVolume
        real(8), dimension(:,:,:), optional, pointer :: DrySedimentVolumeOld
        integer, optional, intent(OUT)               :: STAT

        !External--------------------------------------------------------------
        integer                                      :: ready_        

        !Local-----------------------------------------------------------------
        integer                                      :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(DrySedimentVolume)) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                DrySedimentVolume      => Me%DrySedimentVolume
            endif cd2

cd3 :       if (present(DrySedimentVolumeOld)) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                DrySedimentVolumeOld   => Me%DrySedimentVolumeOld
            endif cd3

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationDrySedVolume


    !--------------------------------------------------------------------------
    
    
    subroutine GetConsolidationWaterFluxes (ConsolidationID,                  &
                                            WaterFluxX,                       &
                                            WaterFluxY,                       &
                                            WaterFluxZ,                       &
                                            STAT)

        !Arguments-------------------------------------------------------------
        integer                                      :: ConsolidationID
        real(8), dimension(:,:,:), pointer, optional :: WaterFluxX 
        real(8), dimension(:,:,:), pointer, optional :: WaterFluxY
        real(8), dimension(:,:,:), pointer, optional :: WaterFluxZ
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer :: ready_        

        !Local-----------------------------------------------------------------
        integer                                      :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                  &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(WaterFluxX)) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                WaterFluxX => Me%WaterFluxes%X
            end if cd2


cd3 :       if (present(WaterFluxY)) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                WaterFluxY => Me%WaterFluxes%Y
            end if cd3


cd4 :       if (present(WaterFluxZ)) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                WaterFluxZ => Me%WaterFluxes%Z
            end if cd4


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationWaterFluxes


    !--------------------------------------------------------------------------

   
    subroutine GetConsolidationVelocityZ (ConsolidationID, Velocity_W,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ConsolidationID
        real, dimension(:,:,:), optional, pointer   :: Velocity_W
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(Velocity_W     )) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                Velocity_W      => Me%VelocityW
            endif cd2

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationVelocityZ    

    !--------------------------------------------------------------------------

    subroutine GetConsolidationVelocityXY (ConsolidationID, Velocity_U, Velocity_V,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                   :: ConsolidationID
        real, dimension(:,:,:), optional, pointer :: Velocity_U, Velocity_V
        integer, optional, intent(OUT)            :: STAT

        !External--------------------------------------------------------------
        integer                                   :: ready_        

        !Local-----------------------------------------------------------------
        integer                                   :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(Velocity_U     )) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                Velocity_U      => Me%VelocityU
            endif cd2

cd3 :       if (present(Velocity_V     )) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)
                Velocity_V      => Me%VelocityV
            endif cd3


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationVelocityXY

    !--------------------------------------------------------------------------
    
    
    subroutine GetConsolidationTimeStep(ConsolidationID, TimeStep,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                    :: ConsolidationID
        real, optional,intent(OUT)                 :: TimeStep
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer                                    :: ready_        

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(TimeStep   )) then

            TimeStep = Me%DT

            endif cd2

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationTimeStep


    !--------------------------------------------------------------------------
    
    subroutine GetConsolidationMinThickness(ConsolidationID, MinThickness,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                    :: ConsolidationID
        real, optional,intent(OUT)                 :: MinThickness
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer                                    :: ready_        

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(MinThickness)) then

                MinThickness = Me%MinLayerThickness

            endif cd2

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationMinThickness


    !--------------------------------------------------------------------------
    
    subroutine GetConsolidationKTopState(ConsolidationID, KTopState,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                    :: ConsolidationID
        integer, dimension(:,:), pointer           :: KTopState
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer                                    :: ready_        

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

            KTopState => Me%KTopState


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_


    end subroutine GetConsolidationKTopState


    !--------------------------------------------------------------------------

    subroutine GetConsolidationOptions(ConsolidationID, ComputeWaterFlux,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                    :: ConsolidationID
        logical, optional, intent(OUT)             :: ComputeWaterFlux
        integer, optional, intent(OUT)             :: STAT

        !External--------------------------------------------------------------
        integer                                    :: ready_        

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR. &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            if (present(ComputeWaterFlux)) then

                ComputeWaterFlux = Me%Consolidation

            endif

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationOptions
    
    
    !--------------------------------------------------------------------------


    subroutine GetConsolidationWaterPercentage(ConsolidationID, WaterPercentage,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                    :: ConsolidationID
        real, dimension(:,:,:),          pointer   :: WaterPercentage
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer                                    :: ready_        

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

            WaterPercentage => Me%WaterPercentage


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationWaterPercentage
    
    
    !--------------------------------------------------------------------------
    
    
    subroutine GetConsolidationTortuosity(ConsolidationID, Tortuosity,  STAT) 

        !Arguments-------------------------------------------------------------
        integer                                    :: ConsolidationID
        real, dimension(:,:,:),optional, pointer   :: Tortuosity
        integer, optional, intent(OUT)             :: STAT


        !External--------------------------------------------------------------
        integer                                    :: ready_        

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd2 :       if (present(Tortuosity   )) then
                call Read_Lock(mCONSOLIDATION_, Me%InstanceID)

                Tortuosity => Me%Tortuosity

            endif cd2

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetConsolidationTortuosity

    !--------------------------------------------------------------------------
    
    
    subroutine SetSedimentDryDensity(ConsolidationID, SedimentDryDensity,  STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                 :: ConsolidationID
        real, dimension(:,:,:),optional, pointer:: SedimentDryDensity
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_            
        integer                                 :: ready_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            Me%ExternalVar%SedimentDryDensity => SedimentDryDensity
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    
    
    end subroutine SetSedimentDryDensity

    
    !--------------------------------------------------------------------------
    
    subroutine SetConsolidationFlux(ConsolidationID, ConsolidationFlux,  STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                 :: ConsolidationID
        real, dimension(:,:),optional, pointer  :: ConsolidationFlux
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_            
        integer                                 :: ready_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            Me%ExternalVar%ConsolidationFlux => ConsolidationFlux
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    
    
    end subroutine SetConsolidationFlux

    !--------------------------------------------------------------------------


    subroutine UngetConsolidation2D_I(ConsolidationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                          :: ConsolidationID
        integer, pointer, dimension(:,:) :: Array
        integer, optional, intent (OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                          :: ready_   

        !Local-----------------------------------------------------------------
        integer                          :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mCONSOLIDATION_, Me%InstanceID, "UngetConsolidation2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UngetConsolidation2D_I

    
    !--------------------------------------------------------------------------
    
    subroutine UngetConsolidation2D(ConsolidationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                          :: ConsolidationID
        real(4), pointer, dimension(:,:) :: Array
        integer, optional, intent (OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                          :: ready_   

        !Local-----------------------------------------------------------------
        integer                          :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mCONSOLIDATION_, Me%InstanceID, "UngetConsolidation2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UngetConsolidation2D


    !--------------------------------------------------------------------------

    
    subroutine UngetConsolidation3D(ConsolidationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: ConsolidationID
        real(4), pointer, dimension(:,:,:)  :: Array
        integer, optional, intent (OUT)     :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_   

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mCONSOLIDATION_, Me%InstanceID, "UngetConsolidation3D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UngetConsolidation3D


    !--------------------------------------------------------------------------

    
    subroutine UngetConsolidation3Dreal8(ConsolidationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: ConsolidationID
        real(8), pointer, dimension(:,:,:) :: Array
        integer, optional, intent (OUT)    :: STAT

        !External--------------------------------------------------------------
        integer                            :: ready_   

        !Local-----------------------------------------------------------------
        integer                            :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mCONSOLIDATION_, Me%InstanceID, "UngetConsolidation3Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UngetConsolidation3Dreal8


    !--------------------------------------------------------------------------

   
    subroutine UngetConsolidation2Dreal8(ConsolidationID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                          :: ConsolidationID
        real(8), pointer, dimension(:,:) :: Array
        integer, optional, intent (OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                          :: ready_   

        !Local-----------------------------------------------------------------
        integer                          :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mCONSOLIDATION_, Me%InstanceID, "UngetConsolidation2Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UngetConsolidation2Dreal8



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine ModifyConsolidation(ConsolidationID, STAT) 
    
        !Parameter-------------------------------------------------------------
        integer                                     :: ConsolidationID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, STAT_, Ready_ 

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleConsolidation", "ModifyConsolidation")

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Actualized the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyConsolidation - ModuleConsolidation - ERR01'


            if (Me%ExternalVar%Now >= Me%NextCompute) then     

                call ComputeConsolidation

                Me%NextCompute =  Me%NextCompute + Me%DT

            endif

            call null_time   (Me%ExternalVar%Now)

            STAT_ = SUCCESS_
        else cd1            
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        if (MonitorPerformance) call StopWatch ("ModuleConsolidation", "ModifyConsolidation")


    end subroutine ModifyConsolidation

    
    !--------------------------------------------------------------------------


    subroutine ComputeConsolidation
    
        !--------------------------------------------------------------------------------------------------

        call ReadLockExternalModules  


        if(Me%Consolidation)then

            call ComputeVerticalCoordinate

            !call ComputeVelocity

            !call ComputeWaterVelocity

            call New_Geometry

            !call ComputePorosity

            call ComputeVolumes

            !call ComputeCellCenterDepth

            if(Me%Decayment)then 

                !only re-computes if there is consolidation decayment
                !otherwise shear strength remains constant in time
                call ComputeShearStrength

            endif

            call ComputeTopShear

            !call ComputeWaterFluxes

            call ComputeTortuosity


        end if


        call ConsolidationOutPut      

        call ReadUnLockExternalModules


    end subroutine ComputeConsolidation
    
    !--------------------------------------------------------------------------
    subroutine ComputeVerticalCoordinate
    
        !Local-----------------------------------------------------------------
        integer                             :: i, j, WKUB, TotalWKUB
        real                                :: DZ, TopLayerThickness
        real                                :: ExcessDZ, ExcessMin

        !----------------------------------------------------------------------
        
        TotalWKUB = Me%WorkSize%KUB

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (Me%ExternalVar%WaterPoints3D (i ,j, TotalWKUB) == WaterPoint) then

                WKUB = Me%KTop(i, j)

             
                DZ = (Me%ExternalVar%ConsolidationFlux(i,j) * Me%DT) / &
                     (Me%ExternalVar%SedimentDryDensity(i,j,WKUB)    * &
                     (1.-Me%Porosity%Field(i,j,WKUB)))

                TopLayerThickness = Me%ExternalVar%DWZ(i,j,WKUB)

                if    (DZ > 0.)then !consolidation

                    ExcessDZ = TopLayerThickness + DZ - Me%MaxLayerThickness

                    if(ExcessDZ > 0.)then

                        Me%KTop(i,j) = WKUB+1

                        Me%KTopState(i,j) = 1

                        if(Me%KTop(i,j) > TotalWKUB)then
                            Me%SedimentColumnFull(i,j) = 1
                            Me%KTop(i,j) = TotalWKUB
                            write(*,*) 'Exceeded maximum number of layers in cell i=',i,'j=',j 
                            write(*,*) 'ComputeVerticalCoordinate - Consolidation - WRN10' 
!                            stop 'Exceeded maximum number of layers' 
                        end if

                        !NewLayerThickness = Me%MinLayerThickness

                        ExcessMin = Me%MinLayerThickness - ExcessDZ


                        Me%VerticalCoordinate(i,j,WKUB)             = Me%VerticalCoordinate(i,j,WKUB-1)     - &
                                                                      Me%MaxLayerThickness + ExcessMin

  
                        Me%VerticalCoordinate(i,j,WKUB+1:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB)       - &
                                                                      Me%MinLayerThickness

                        Me%DrySedimentHeight (i,j,WKUB)             = (Me%VerticalCoordinate(i,j,WKUB-1)    - &
                                                                       Me%VerticalCoordinate(i,j,WKUB  ))   * &
                                                                     (1. - Me%Porosity%Field(i,j,WKUB))
                        
                        Me%DrySedimentHeight (i,j,WKUB+1)           = (Me%VerticalCoordinate(i,j,WKUB  )    - &
                                                                       Me%VerticalCoordinate(i,j,WKUB+1))   * &
                                                                     (1. - Me%Porosity%Field(i,j,WKUB+1))
                    else

                        Me%VerticalCoordinate(i,j,WKUB:TotalWKUB)   = Me%VerticalCoordinate(i,j,WKUB) - DZ


                        Me%DrySedimentHeight (i,j,WKUB)             = Me%DrySedimentHeight (i,j,WKUB)         + &
                                                                      Me%ExternalVar%ConsolidationFlux(i,j)   * &
                                                                      Me%DT / Me%ExternalVar%SedimentDryDensity(i,j,WKUB)
                        Me%KTopState(i,j) = 0


                    end if


                elseif(DZ <  0.)then !erosion

                    ExcessDZ = Me%ExternalVar%DWZ(i,j,WKUB) - Me%MinLayerThickness

                    if(abs(DZ+0.001*DZ) .ge. ExcessDZ)then

                        Me%VerticalCoordinate(i,j,WKUB-1) = Me%VerticalCoordinate(i,j,WKUB-1)       &
                                                            - Me%MinLayerThickness
                                                            
                        Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB-1)

                        Me%KTop(i,j)      = WKUB-1

                        Me%KTopState(i,j) = -1


                        if(Me%KTop(i,j) < Me%WorkSize%KLB)then
                            stop 'Eroded all sediment layers'
                        endif

                        Me%DrySedimentHeight (i,j,WKUB)   = 0.

                        Me%DrySedimentHeight (i,j,WKUB-1) = Me%DrySedimentHeight (i,j,WKUB-1)                   + &
                                                            Me%MinLayerThickness                                * &
                                                            Me%ExternalVar%SedimentDryDensity(i,j,WKUB-1)


                    else

                        Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB) - DZ

                        Me%DrySedimentHeight (i,j,WKUB)           = Me%DrySedimentHeight (i,j,WKUB)             + &
                                                                    Me%ExternalVar%ConsolidationFlux(i,j)       * &
                                                                    Me%DT / Me%ExternalVar%SedimentDryDensity(i,j,WKUB)

                        Me%KTopState(i,j) = 0

                    endif

                end if

            end if

        enddo
        enddo
    
    
    end subroutine ComputeVerticalCoordinate

    !--------------------------------------------------------------------------
    
    !Computes consolidation velocity
    subroutine ComputeVelocity 
    
        !Local-----------------------------------------------------------------
        integer                             :: i, j, k, WKLB, WKUB
        real                                :: DZ, MaxDZ


        !----------------------------------------------------------------------
        
        WKUB = Me%WorkSize%KUB
        WKLB = Me%WorkSize%KLB

        if(Me%Decayment)then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                Me%VelocityW(i,j,WKLB)       = 0.0

                Me%WaterVelocityW(i,j,WKLB)  = 0.0

                do k = WKLB + 1 , WKUB + 1

                    if (Me%ExternalVar%WaterPoints3D (i ,j, k) == WaterPoint) then

                        Me%VelocityU(i,j,k) = 0.0

                        Me%VelocityV(i,j,k) = 0.0 

                        Me%VelocityW(i,j,k) = ((Me%ExternalVar%DWZ  (i,j,k-1))**2)  / &
                                               (Me%DrySedimentHeight(i,j,k-1))      * &
                                               (Me%StationaryPorosity%Field(i,j,k)  - &
                                                Me%Porosity%Field (i,j,k-1))        / &
                                                Me%DecayTime

                        if(k == WKUB + 1)then
                        
                            !If BoundaryConsolidationFlux is positive, sediment is ENTERING the domain(CONSOLIDATION)
                            !If BoundaryConsolidationFlux is negative, sediment is EXITING  the domain(EROSION)

                            !                        [kgsed m-2 s-1]
                            ![m/s] = [m/s] +  ----------------------------------
                            !                  [m3sed m-3bulk] * [kgsed m-3sed]

                            Me%VelocityW(i,j,k) = Me%VelocityW(i,j,k)                       + &
                                                  (Me%ExternalVar%ConsolidationFlux(i,j))   / &
                                                  ((1-Me%Porosity%Field(i,j,k-1))           * &
                                                  Me%ExternalVar%SedimentDryDensity(i,j,k-1))
                            
                            DZ    = Me%VelocityW(i,j,k) * Me%DT

                            MaxDZ = Me%MinLayerThickness - Me%ExternalVar%DWZ(i,j,k-1)
                            
                            if(dz < MaxDZ)then

                                Me%VelocityW(i,j,k) = MaxDZ / Me%DT

                            endif



                        end if

                    endif

                enddo
            enddo
            enddo


        else

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                Me%VelocityW     (i,j,WKLB)  = 0.0
                Me%WaterVelocityW(i,j,WKLB)  = 0.0

                do k = WKLB + 1, WKUB + 1

                    if (Me%ExternalVar%WaterPoints3D (i ,j, WKUB) == WaterPoint) then

                        Me%VelocityU(i,j,k) = 0.0

                        Me%VelocityV(i,j,k) = 0.0 

                        Me%VelocityW(i,j,k) = 0.0

                        if(k == WKUB + 1)then
                        
                            !If BoundaryConsolidationFlux is positive, sediment is ENTERING the domain(CONSOLIDATION)
                            !If BoundaryConsolidationFlux is negative, sediment is EXITING  the domain(EROSION)

                            !                        [kgsed m-2 s-1]
                            ![m/s] = [m/s] +  ----------------------------------
                            !                  [m3sed m-3bulk] * [kgsed m-3sed]

                            Me%VelocityW(i,j,k) = Me%VelocityW(i,j,k)                       + &
                                                  (Me%ExternalVar%ConsolidationFlux(i,j))   / &
                                                  ((1-Me%Porosity%Field(i,j,k-1))           * &
                                                  Me%ExternalVar%SedimentDryDensity(i,j,k-1))

                            DZ    = Me%VelocityW(i,j,k) * Me%DT

                            MaxDZ = Me%MinLayerThickness - Me%ExternalVar%DWZ(i,j,k-1)
                            
                            if(DZ < MaxDZ)then

                                !ExcessDZ = MaxDZ - DZ

                                Me%VelocityW(i,j,k)   = MaxDZ / Me%DT

                                !Me%VelocityW(i,j,k-1) = ExcessDZ / Me%DT

                            endif

                        end if

                    endif

                enddo
            enddo
            enddo

        endif



    end subroutine ComputeVelocity

    !--------------------------------------------------------------------------

    !Computes water expulsion velocity
    subroutine ComputeWaterVelocity 
    
        !Local-----------------------------------------------------------------
        integer                            :: i, j, k

        !----------------------------------------------------------------------

        
        do j = Me%WorkSize%JLB,     Me%WorkSize%JUB
        do i = Me%WorkSize%ILB,     Me%WorkSize%IUB
        do k = Me%WorkSize%KLB + 1, Me%WorkSize%KUB + 1

            if (Me%ExternalVar%WaterPoints3D (i, j, k) == WaterPoint) then

                Me%WaterVelocityU(i,j,k) = 0.0

                Me%WaterVelocityV(i,j,k) = 0.0 

                Me%WaterVelocityW(i,j,k) = Me%WaterVelocityW(i,j,k-1) + (-1. * Me%VelocityW(i,j,k))

            endif

        enddo
        enddo
        enddo

    end subroutine ComputeWaterVelocity
    
    !--------------------------------------------------------------------------

    subroutine New_Geometry

        !External----------------------------------------------------------------
        integer, pointer, dimension(:,:,:)  :: WaterPoints3D
        integer                             :: STAT_CALL
        
        !Local-------------------------------------------------------------------
        type(T_Time)                        :: ActualTime

        !------------------------------------------------------------------------

        call ReadUnLockExternalModules  
        
        call GetWaterPoints3D(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError(FATAL_, INTERNAL_, "New_Geometry; ModuleConsolidation. ERR01")
        
        ActualTime = Me%ExternalVar%Now

        !Compute new volume 
        call ComputeVerticalGeometry(Me%ObjGeometry,                                     &
                                     WaterPoints3D      = WaterPoints3D,                 &
                                     ActualTime         = ActualTime,                    &
                                     SZZ                = Me%VerticalCoordinate,         &
                                     KTOP               = Me%KTop,                       &
                                     STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError(FATAL_, INTERNAL_, "New_Geometry; ModuleConsolidation. ERR02")
        
        call UnGetMap(Me%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            call SetError(FATAL_, INTERNAL_, "New_Geometry; ModuleConsolidation. ERR03") 

        
        call ReadLockExternalModules   
                   

    end subroutine New_Geometry

    !----------------------------------------------------------------------------

    subroutine ComputeVolumes

        !Local-----------------------------------------------------------------
        integer                            :: i, j, k

        !----------------------------------------------------------------------

        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            if (Me%ExternalVar%WaterPoints3D(i,j,k) == WaterPoint) Then

                Me%WaterVolumeOld (i,j,k)      = Me%WaterVolume (i,j,k)

                Me%DrySedimentVolumeOld(i,j,k) = Me%DrySedimentVolume (i,j,k)

                Me%WaterVolume (i,j,k)         = Me%Porosity%Field (i,j,k) * Me%ExternalVar%Volume(i,j,k)
            
                Me%DrySedimentVolume(i,j,k)    = (1.0  - Me%Porosity%Field (i,j,k)) * Me%ExternalVar%Volume(i,j,k)

             end if
        
        end do
        end do
        end do


    end subroutine ComputeVolumes
   
    !----------------------------------------------------------------------------

    subroutine ComputeShearStrength

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, TotalKUB
        real                                        :: Depth
        
        !Begin-----------------------------------------------------------------

        TotalKUB = Me%WorkSize%KUB


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB

            if(Me%ExternalVar%WaterPoints3D (i, j, k) == WaterPoint) then

                if(Me%CellCenterDepth(i,j,k) < 0.)then
                    Depth = 0.
                else
                    Depth = Me%CellCenterDepth(i,j,k)
                endif

                Me%CriticalShearStress%Field(i,j,k) = Me%Infinite_CSE                   + &
                                                     (Me%Surface_CSE - Me%Infinite_CSE) * &
                                                      exp(-Depth/Me%CSE_Coef)

            end if

        enddo
        enddo
        enddo


    end subroutine ComputeShearStrength
    
    !-------------------------------------------------------------------------

    subroutine ComputeTopShear
    
        !Local-----------------------------------------------------------------        
        
        integer                            :: i, j

        !---------------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if(Me%ExternalVar%WaterPoints3D (i, j, Me%WorkSize%KUB) == WaterPoint) then

                Me%TopCriticalShear (i,j) = Me%CriticalShearStress%Field(i,j,Me%Ktop(i,j))

            end if
        
        enddo
        enddo


    end subroutine ComputeTopShear


    !--------------------------------------------------------------------------
    
    !Computes porosity as an imposition to avoid two pararel evolution which would lead to error 
    subroutine ComputePorosity 
    
        !Local-----------------------------------------------------------------        
        integer                            :: i, j, k

        !---------------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
     
            if(Me%ExternalVar%WaterPoints3D (i, j, Me%WorkSize%KUB) == WaterPoint) then

                do k = Me%WorkSize%KLB, Me%Ktop(i,j)

                    Me%Porosity%Field (i,j,k)    = (Me%ExternalVar%DWZ  (i,j,k)  -   &
                                                    Me%DrySedimentHeight(i,j,k)) / Me%ExternalVar%DWZ(i,j,k)

                enddo

            end if
        
        enddo
        enddo


    end subroutine ComputePorosity


    !--------------------------------------------------------------------------
    
    !Computes porosity as an imposition to avoid two pararel evolution which would lead to error 
    subroutine ComputeTortuosity
    
        !Local-----------------------------------------------------------------        
        integer                            :: i, j, k

        !---------------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if(Me%ExternalVar%WaterPoints3D (i, j, Me%WorkSize%KUB) == WaterPoint) then

                do k = Me%WorkSize%KLB, Me%KTop(i,j)

                    !Boudreau, 1996
                    Me%Tortuosity(i,j,k) = 1. - log(Me%Porosity%Field(i,j,k)**2)

                enddo

            end if

        enddo
        enddo


    end subroutine ComputeTortuosity
        
    !--------------------------------------------------------------------------
    
    subroutine ComputeWaterFluxes
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k
       
        !Begin-----------------------------------------------------------------
        


        do k = Me%WorkSize%KLB, Me%WorkSize%KUB   
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%ComputeFacesU3D(I,J,K) .EQ. Compute) then

                Me%WaterFluxes%X(i, j, k) = Me%WaterVelocityU(i, j, k) * Me%ExternalVar%AreaU(i, j, k) 

            else 

                Me%WaterFluxes%X(i, j, k) = 0.0

            end if 
 
        enddo
        enddo
        enddo

        
        do k = Me%WorkSize%KLB, Me%WorkSize%KUB   
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%ComputeFacesV3D(I,J,K) .EQ. Compute) then

                Me%WaterFluxes%Y(i, j, k) = Me%WaterVelocityV(i, j, k) * Me%ExternalVar%AreaV(i, j, k) 

            else

                Me%WaterFluxes%Y(i, j, k) = 0.0
                 
            end if

        enddo
        enddo
        enddo

        do k = Me%WorkSize%KLB, Me%WorkSize%KUB   
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                
            if (Me%ExternalVar%ComputeFacesW3D(I,J,K) .EQ. Compute) then

                Me%WaterFluxes%Z(i, j, k) = Me%WaterVelocityW(i, j, k) * Me%ExternalVar%GridCellArea(i,j) 

                if(k == Me%WorkSize%KUB) then

                    !REVIEW THIS!!!REVIEW THIS!!!REVIEW THIS!!!

                    !SedSurfaceFlux(i, j) = WaterVelocityW(i, j, k + 1) * Me%ExternalVar%GridCellArea(i,j) 

                    !REVIEW THIS!!!REVIEW THIS!!!REVIEW THIS!!!

                end if

            else 

                Me%WaterFluxes%Z(i, j, k) = null_real

            end if

        enddo
        enddo
        enddo

    end subroutine ComputeWaterFluxes

    
    !--------------------------------------------------------------------------


    subroutine ComputeCellCenterDepth
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k, TotalKUB, KTop
        
        !Begin-----------------------------------------------------------------

        TotalKUB = Me%WorkSize%KUB


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if(Me%ExternalVar%WaterPoints3D (i, j, Me%WorkSize%KUB) == WaterPoint) then

                KTop = Me%KTop(i,j)

                do k = KTop, Me%WorkSize%KLB, -1

                    !Depth at the center of each cell relatively to initial zero
                    Me%CellCenterDepth(i,j,k) =  Me%ExternalVar%SZZ(i,j,k  ) + &
                                                (Me%ExternalVar%SZZ(i,j,k-1) - &
                                                 Me%ExternalVar%SZZ(i,j,k  ))/2.
                enddo

                Me%CellCenterDepth(i,j,KTop:TotalKUB) = Me%CellCenterDepth(i,j,KTop)

            end if

        enddo
        enddo


    end subroutine ComputeCellCenterDepth

    
    !--------------------------------------------------------------------------

    
    subroutine ConsolidationOutPut

        !Local-----------------------------------------------------------------
        integer :: NextOutPut

        !----------------------------------------------------------------------

        if (Me%OutPut%True) then 

            NextOutPut = Me%OutPut%NextOutPut

            if (Me%ExternalVar%Now .GE. Me%OutPut%OutTime(NextOutPut)) then

                call Write_HDF_Format 

                Me%OutPut%NextOutPut = NextOutPut + 1

            end if 

        end if  

        if (Me%TimeSerie   )  call OutPut_TimeSeries

        if (Me%BoxTimeSerie)  call OutPut_BoxTimeSeries
                

        !----------------------------------------------------------------------

    end subroutine ConsolidationOutPut

    
    !--------------------------------------------------------------------------

    
    subroutine Write_HDF_Format

        !Local-----------------------------------------------------------------
        integer                              :: STAT_CALL
        real, dimension(6),          target  :: AuxTime
        real, dimension(:),          pointer :: TimePtr
        integer                              :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                              :: WorkKLB, WorkKUB, Index
        real                                 :: TotalSeconds

        !----------------------------------------------------------------------

        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 


        !Current output index
        Index        = Me%OutPut%NextOutPut
        TotalSeconds = Me%ExternalVar%Now - Me%OutPut%OutTime(1)

        !Writes current time
        call ExtractDate   (Me%ExternalVar%Now, AuxTime(1), AuxTime(2), AuxTime(3),     &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime
        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR00'

        call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",        &
                             Array1D = TimePtr, OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR01'

        !Writes SZZ
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR02'

        call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",                 &
                             "m", Array3D = Me%ExternalVar%SZZ, OutputNumber = Index,   &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR03'


        !Writes OpenPoints
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                              &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,                        &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR04'

        call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",              &
                             "-", Array3D = Me%ExternalVar%OpenPoints3D,                &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR05'

        call HDF5WriteData  (Me%ObjHDF5, "/Results/"//Char_Porosity,                    &
                             Char_Porosity, "-", Array3D = Me%Porosity%Field,           &
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR06'


        call HDF5WriteData  (Me%ObjHDF5, "/Results/Top Critical shear",                 &
                             "Top Critical shear", "N/m2", Array2D = Me%TopCriticalShear,&
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR12'
        
        
        call HDF5WriteData  (Me%ObjHDF5, "/Results/"//"Critical shear",                &
                             "Critical shear", "-", Array3D = Me%CriticalShearStress%Field,&
                             OutputNumber = Index, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR06'

    

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF_Format - Consolidation - ERR23'


    end subroutine Write_HDF_Format

    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D    = Me%Porosity%Field,                              &
                            STAT      = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR01')
        

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D   = Me%ExternalVar%DWZ,                              &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR02')

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D   = Me%WaterVelocityW,                               &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR03')

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D_8 = Me%WaterVolume,                                  &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR04')

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D   = Me%DrySedimentHeight,                                 &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR05')

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D_8 = Me%DrySedimentVolume,                            &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR06')

        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D_8 = Me%DrySedimentVolumeOld,                         &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR07')


        call WriteTimeSerie(Me%ObjTimeSerie,                                            &
                            Data3D   = Me%CellCenterDepth,                              &
                            STAT     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError(FATAL_, INTERNAL_,'OutPut_TimeSeries - ModuleConsolidation - ERR08')
   
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------

    subroutine OutPut_BoxTimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        
        !Integration of fluxes and time serie output
        call BoxDif(Me%ObjBoxDif,                        &
                    Me%WaterFluxes%X,                    &
                    Me%WaterFluxes%Y,                    &
                    Me%WaterFluxes%Z,                    &
                    'WaterVolume',                       &
                    Me%ExternalVar%WaterPoints3D,        &
                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'OutPut_BoxTimeSeries - Consolidation - ERR01'

        call BoxDif(Me%ObjBoxDif,                        &
                    Me%WaterVolume,                      &
                    'WaterVolume',                       &
                    Me%ExternalVar%WaterPoints3D,        &
                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'OutPut_BoxTimeSeries - Consolidation - ERR02'

        call BoxDif(Me%ObjBoxDif,                        &
                    Me%WaterFluxes%X,                    &
                    Me%WaterFluxes%Y,                    &
                    Me%WaterFluxes%Z,                    &
                    'Sediment_Mass',                     &
                    Me%ExternalVar%WaterPoints3D,        &
                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'OutPut_BoxTimeSeries - Consolidation - ERR03'

        call BoxDif(Me%ObjBoxDif,                        &
                    Me%DrySedimentVolume,                &
                    'Sediment_Mass',                     &
                    Me%ExternalVar%WaterPoints3D,        &
                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'OutPut_BoxTimeSeries - Consolidation - ERR04'

    end subroutine OutPut_BoxTimeSeries


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillConsolidation(ConsolidationID, STAT)

        !Arguments---------------------------------------------------------------
        integer                        :: ConsolidationID
        integer, optional, intent(OUT) :: STAT

        !External----------------------------------------------------------------
        integer                        :: ready_              
        integer                        :: STAT_CALL, nUsers

        !Local-------------------------------------------------------------------
        integer                        :: STAT_            

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ConsolidationID, ready_)


cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mCONSOLIDATION_,  Me%InstanceID)

            if (nUsers == 0) then

                if (Me%TimeSerie)then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillConsolidation - ModuleConsolidation - ERR01'
                end if
                
                if (Me%BoxTimeSerie) then
                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillConsolidation - ModuleConsolidation - ERR02'
                end if
                
                if (Me%Output%True) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'KillConsolidation - ModuleConsolidation - ERR03'
                end if
                
                call Write_Final_Consolidation_HDF 

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillConsolidation - ModuleConsolidation - ERR04'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjGridData)
                if (nUsers == 0) stop 'KillConsolidation - ModuleConsolidation - ERR05'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillConsolidation - ModuleConsolidation - ERR06'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillConsolidation - ModuleConsolidation - ERR07'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'KillConsolidation - ModuleConsolidation - ERR08'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjMap)
                if (nUsers == 0) stop 'KillConsolidation - ModuleConsolidation - ERR09'

                deallocate(Me%Elevation, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR10") 
                nullify(Me%Elevation)
                
                deallocate(Me%KTopState, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR10c") 
                nullify(Me%Elevation)
                
                deallocate(Me%TopCriticalShear, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR10a") 
                nullify(Me%TopCriticalShear)

                deallocate(Me%SedimentColumnFull, STAT = STAT_CALL) 
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR10b") 
                nullify(Me%SedimentColumnFull)

                deallocate(Me%VerticalCoordinate, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR11") 
                nullify(Me%VerticalCoordinate)

                deallocate(Me%KTop, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR11a") 
                nullify(Me%KTop)

                deallocate(Me%DrySedimentHeight, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR14") 
                nullify(Me%DrySedimentHeight)

                deallocate(Me%WaterVolume, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR15") 
                nullify(Me%WaterVolume)

                deallocate(Me%WaterVolumeOld, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR16") 
                nullify(Me%WaterVolumeOld)

                deallocate(Me%DrySedimentVolume, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR17") 
                nullify(Me%DrySedimentVolume)

                deallocate(Me%DrySedimentVolumeOld, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR18") 
                nullify(Me%DrySedimentVolumeOld)

                deallocate(Me%Porosity%Field, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR20") 
                nullify(Me%Porosity%Field)

                if(Me%Decayment)then
                    deallocate(Me%StationaryPorosity%Field, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR20") 
                    nullify(Me%StationaryPorosity%Field)
                end if

                deallocate(Me%Tortuosity, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR21") 
                nullify(Me%Tortuosity)

                deallocate(Me%velocityU, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR22") 
                nullify(Me%velocityU)

                deallocate(Me%velocityV, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR23") 
                nullify(Me%velocityV)

                deallocate(Me%velocityW, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR24") 
                nullify(Me%velocityW)

                deallocate(Me%WatervelocityU, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR25") 
                nullify(Me%WatervelocityU)


                deallocate(Me%WatervelocityV, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR26") 
                nullify(Me%WatervelocityV)

                deallocate(Me%WatervelocityW, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR27") 
                nullify(Me%WatervelocityW)

                deallocate(Me%WaterFluxes%X, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR28") 
                nullify(Me%WaterFluxes%X)


                deallocate(Me%WaterFluxes%Y, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR29") 
                nullify(Me%WaterFluxes%Y)


                deallocate(Me%WaterFluxes%Z, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR30") 
                nullify(Me%WaterFluxes%Z)

                deallocate(Me%WaterPercentage, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR32") 
                nullify(Me%WaterPercentage)

                deallocate(Me%CellCenterDepth, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_)                                    &
                    call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR33") 
                nullify(Me%CellCenterDepth)


                if(associated(Me%CriticalShearStress%Field))then
                    deallocate(Me%CriticalShearStress%Field, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        call SetError(FATAL_, INTERNAL_, "KillConsolidation; ModuleConsolidation. ERR34") 
                    nullify(Me%CriticalShearStress%Field)
                end if
                
                !Deallocates Instance
                call DeallocateInstance ()

                ConsolidationID  = 0
                STAT_            = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillConsolidation

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Local-----------------------------------------------------------------
        type (T_Consolidation), pointer     :: AuxObjConsolidation
        type (T_Consolidation), pointer     :: PreviousObjConsolidation

        !Updates pointers
        if (Me%InstanceID == FirstObjConsolidation%InstanceID) then
            FirstObjConsolidation => FirstObjConsolidation%Next
        else
            PreviousObjConsolidation        => FirstObjConsolidation
            AuxObjConsolidation             => FirstObjConsolidation%Next
            do while (AuxObjConsolidation%InstanceID /= Me%InstanceID)
                PreviousObjConsolidation    => AuxObjConsolidation
                AuxObjConsolidation         => AuxObjConsolidation%Next
            enddo

            !Now update linked list
            PreviousObjConsolidation%Next   => AuxObjConsolidation%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

    end subroutine DeallocateInstance

   
    !--------------------------------------------------------------------------


    subroutine Write_Final_Consolidation_HDF

        !Local-----------------------------------------------------------------
        integer                             :: ObjHDF5
        real,    pointer, dimension(:,:  )  :: Bathymetry
        real,    pointer, dimension(:,:,:)  :: SZZ
        integer, pointer, dimension(:,:,:)  :: OpenPoints3D
        integer, pointer, dimension(:,:,:)  :: WaterPoints3D
        integer                             :: WorkILB, WorkIUB
        integer                             :: WorkJLB, WorkJUB
        integer                             :: WorkKLB, WorkKUB
        integer                             :: HDF5_CREATE, STAT_CALL

        !----------------------------------------------------------------------
        
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB
        
        ObjHDF5 = 0
        
        !Gets a pointer to Bathymetry
        call GetGridData       (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR01'

        !Gets WaterPoints3D
        call GetWaterPoints3D   (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR02'
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (ObjHDF5,                                                &
                                 trim(Me%Files%Final)//"5",                              &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR03'

        !Writes geometry
        call WriteGeometry(Me%ObjGeometry,                                               &
                           ObjHDF5, ON,                                                  &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR04'
        

        call ReadLockExternalModules   

        
        SZZ          => Me%ExternalVar%SZZ
        OpenPoints3D => Me%ExternalVar%OpenPoints3D

        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                        &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR05'


        !Writes the Grid
        call HDF5WriteData   (ObjHDF5, "/Grid", "Bathymetry", "m",                       &
                              Array2D = Bathymetry,                                      &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR06'

        call HDF5WriteData   (ObjHDF5, "/Grid", "WaterPoints3D", "-",                    &
                              Array3D = WaterPoints3D,                                   &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR07'

        call WriteHorizontalGrid (HorizontalGridID = Me%ObjHorizontalGrid,               &
                                ObjHDF5 = ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR08'

        !Writes SZZ
        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                         &
                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR09'

        call HDF5WriteData  (ObjHDF5, "/Grid", "VerticalZ",                              &
                             "m", Array3D = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR10'


        !Writes OpenPoints
        call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB,                                  &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,                         &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR11'

        call HDF5WriteData  (ObjHDF5, "/Grid", "OpenPoints",                             &
                             "-", Array3D = OpenPoints3D,                                &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR12'

        !Sets limits for next write operations
        call HDF5SetLimits   (ObjHDF5, WorkILB, WorkIUB, WorkJLB,                        &
                              WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR13'

        !Writes Porosity
        call HDF5WriteData  (ObjHDF5,                                                    &
                             "/"//Char_Porosity,                                         &
                             Char_Porosity, "-",                                         &
                             Array3D = Me%Porosity%Field,                                &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR14'


        !Writes Porosity
        call HDF5WriteData  (ObjHDF5,                                                    &
                             "/"//Char_SedimentColumnFull,                               &
                             Char_SedimentColumnFull, "-",                               &
                             Array2D = Me%SedimentColumnFull,                            &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR14a'


        !Writes everything to disk
        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR15'

        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR16'

        !Ungets the Bathymetry
        call UngetGridData (Me%ObjGridData, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR17'

        !Ungets the WaterPoints
        call UnGetMap        (Me%ObjMap, WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'Write_Final_Consolidation_HDF - Consolidation - ERR18'



        call ReadUnLockExternalModules   


    end subroutine Write_Final_Consolidation_HDF

    
    !--------------------------------------------------------------------------


    subroutine ReadUnLockExternalModules

        !Local-------------------------------------------------------------------

        integer :: STAT_CALL

        !------------------------------------------------------------------------


        call UnGetMap(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR01") 


        call UnGetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesU3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR02") 


        call UnGetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesV3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR03") 


        call UnGetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesW3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR04") 


        call UnGetMap(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR05") 


        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%AreaU, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR06") 


        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%AreaV, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR07") 


        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%Volume, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR08") 

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%VolumeOld, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR09") 


        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%KTop, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR10") 


        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR11") 

        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR12") 


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            call SetError(FATAL_, INTERNAL_, "ReadUnLockExternalModules; ModuleConsolidation. ERR13") 
        
    end subroutine ReadUnLockExternalModules
    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine ReadLockExternalModules

        !Local-------------------------------------------------------------------
        integer :: STAT_CALL

        !------------------------------------------------------------------------


        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR02")

        call GetGeometryAreas(Me%ObjGeometry,                                       &
                              Me%ExternalVar%AreaU,                                 &
                              Me%ExternalVar%AreaV,                                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR04")


        call GetGeometryDistances(Me%ObjGeometry,                                   &
                                  SZZ         = Me%ExternalVar%SZZ,                 &
                                  DWZ         = Me%ExternalVar%DWZ,                 &
                                  STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR05")


        call GetGeometryVolumes(Me%ObjGeometry,                                     &
                                VolumeZ = Me%ExternalVar%Volume,                    &
                                VolumeZOld = Me%ExternalVar%VolumeOld,              &
                                STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR06")


        call GetGeometryKTop(Me%ObjGeometry,                                        &
                             KTopZ  = Me%ExternalVar%KTop,                          &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR07")

        call GetOpenPoints3D(Me%ObjMap,                                             &
                             Me%ExternalVar%OpenPoints3D,                           &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR08")


        call GetWaterPoints3D(Me%ObjMap,                                            &
                              Me%ExternalVar%WaterPoints3D,                         &
                              STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR09")


        call GetComputeFaces3D(Me%ObjMap,                                           &
                               ComputeFacesU3D = Me%ExternalVar%ComputeFacesU3D,    &
                               ComputeFacesV3D = Me%ExternalVar%ComputeFacesV3D,    &
                               ComputeFacesW3D = Me%ExternalVar%ComputeFacesW3D,    &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalModules; ModuleConsolidation. ERR10")


    end subroutine ReadLockExternalModules

  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjConsolidationID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjConsolidationID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjConsolidationID > 0) then
            call LocateObjConsolidation (ObjConsolidationID)
            ready_ = VerifyReadLock (mCONSOLIDATION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjConsolidation (ObjConsolidationID)

        !Arguments-------------------------------------------------------------
        integer                    :: ObjConsolidationID

        !Local-----------------------------------------------------------------

        Me => FirstObjConsolidation
        do while (associated (Me))
            if (Me%InstanceID == ObjConsolidationID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleConsolidation - LocateObjConsolidation - ERR01'

    end subroutine LocateObjConsolidation

    !--------------------------------------------------------------------------

end module ModuleConsolidation

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

