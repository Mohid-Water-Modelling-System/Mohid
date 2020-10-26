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
    use ModuleDrawing 
    use ModuleHDF5,             only : ConstructHDF5, GetHDF5FileAccess, HDF5SetLimits,         &
                                       HDF5WriteData, HDF5FlushMemory, HDF5ReadData,            &
                                       HDF5ReadWindow, KillHDF5
    use ModuleEnterData           
    use ModuleFillMatrix,       only : ConstructFillMatrix, GetDefaultValue, KillFillMatrix
    use ModuleGridData,         only : ConstructGridData, GetGridData, ModifyGridData,          &
                                       GetGridData2DReference, UngetGridData, KillGridData          
    use ModuleDischarges,       only : GetDischargesNumber, GetDischargesGridLocalization,      &
                                       GetDischargeWaterFlow, GetDischargeConcentration,        &
                                       GetDischargeON, CorrectsCellsDischarges
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,           &
                                       GetTimeSerieLocation, CorrectsCellsTimeSerie,            &
                                       GetNumberOfTimeSeries, TryIgnoreTimeSerie, GetTimeSerieName       
    use ModuleHorizontalMap,    only : GetWaterPoints2D, GetBoundaries, GetOpenPoints2D,        &
                                       GetComputeFaces2D, UnGetHorizontalMap
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, WriteHorizontalGrid,                  &
                                       GetHorizontalGridSize, UnGetHorizontalGrid, GetXYCellZ,  &
                                       GetDDecompMPI_ID, GetDDecompON, GetGridOutBorderPolygon, &
                                       GetDDecompParameters, GetDDecompWorkSize2D,              &
                                       RotateVectorFieldToGrid, RotateVectorGridToField,        &
                                       LocateCell1D
#ifdef _USE_MPI                                                  
    use ModuleHorizontalGrid,   only : ReceiveSendProperitiesMPI
#endif
    
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
    private ::      ConstructHybridMorph

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
    private ::          VanRijn2007Transport
    private ::          VanRijn1Transport !suspended load compute by algebric formula
    private ::          VanRijn2Transport !suspended load compute by numerical integration
    private ::          BailardTransport
    private ::          DibajniaTransport
    private ::          BijkerTransport
    private ::      ComputeEvolution
    private ::      ComputeEvolution_CentralDifferences
    private ::          ComputeDischarges   
    private ::      ComputeSmoothSlope
    private ::      BoundaryCondition
    private ::      ComputeHybridMorphEvolution
    private ::          ComputeAlongShoreFlow1D
    private ::          ComputeProfileCrossShoreMovement
    private ::          ComputeAlongShoreDZ1D
    private ::          ComputeProfileCrossShoreMovementDZ
    private ::          ComputeNewBathymetryFromNewProfiles
    private ::              InterpolateNewBathymProfile
    !Update the bathym increment using the HybridMoprh methodology 
    private ::          HybridMorphNewBathymIncrement 
    private ::      ComputeResidualEvolution
    private ::      OutPutSandHDF
    private ::      OutPut_TimeSeries
    private ::      ComputeTauCritic
    private ::      OutputBoxFluxes
                
    !Destructor
    public  :: KillSand                                                     
    private ::      DeAllocateInstance
    private ::      WriteFinalState
    private ::      KillHybridMorph

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
    character(LEN = StringLength), parameter    :: U_Average_block_begin = '<beginUaverage>'
    character(LEN = StringLength), parameter    :: U_Average_block_end   = '<endUaverage>'
    character(LEN = StringLength), parameter    :: V_Average_block_begin = '<beginVaverage>'
    character(LEN = StringLength), parameter    :: V_Average_block_end   = '<endVaverage>'
    character(LEN = StringLength), parameter    :: Mapp_DZ_block_begin   = '<beginMappDZ>'
    character(LEN = StringLength), parameter    :: Mapp_DZ_block_end     = '<endMappDZ>'    

    integer, parameter :: NoTransport = 0, Ackers = 1, MeyerPeter = 2, VanRijn1 = 3, & 
                          VanRijn2 = 4, Bailard = 5, Dibajnia = 6, Bijker = 7, VanRijn2007 = 8
    integer, parameter :: NullGradient = 1, Cyclic = 2, NullValue = 3

    !Selma
    integer, parameter :: Time_ = 1
    
    !srt(2.)
    real,    parameter :: SquareRoot2 = 1.414213562
    
    !domain side 
    integer, parameter :: ILB_ = 1
    integer, parameter :: IUB_ = 2
    integer, parameter :: JLB_ = 3
    integer, parameter :: JUB_ = 4

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
        real,    pointer, dimension(:,:)        :: VelResidualU     => null()
        real,    pointer, dimension(:,:)        :: VelResidualV     => null()

        real,    pointer, dimension(:,:)        :: VelU_Face        => null()
        real,    pointer, dimension(:,:)        :: VelV_Face        => null()
        real,    pointer, dimension(:,:)        :: VelMod           => null()
        real,    pointer, dimension(:,:)        :: TauWave          => null()
        real,    pointer, dimension(:,:)        :: TauCurrent       => null()
        real,    pointer, dimension(:,:)        :: ShearVelocity    => null()    
        real,    pointer, dimension(:,:)        :: WaveHeight       => null()
        real,    pointer, dimension(:,:)        :: WavePeriod       => null()
        real                                    :: MinWaterColumn   = FillValueReal 
    end type T_External

    private :: T_Residual
    type       T_Residual
!        logical                                 :: ON           = .false.
        type(T_Time)                            :: StartTime    
        real, dimension(:,:), pointer           :: FluxX        => null ()
        real, dimension(:,:), pointer           :: FluxY        => null ()        
        real, dimension(:,:), pointer           :: FluxZ        => null ()                
        real, dimension(:,:), pointer           :: OutFluxX     => null ()
        real, dimension(:,:), pointer           :: OutFluxY     => null ()        
    end type   T_Residual
   
    private :: T_Aceleration
    type       T_Aceleration
        logical                                 :: Yes          = .false.
        real                                    :: Coef         = FillValueReal
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
        real                                    :: DZDT         = FillValueReal        
        real                                    :: BathymDT     = FillValueReal
        type (T_Time)                           :: NextSand, NextBatim, NextDZ
        logical                                 :: Bathym       = .false. 
        !Selma
        integer                                 :: BathymType   = Time_
        real                                    :: Gama         = FillValueReal
        integer                                 :: NumericMethod= FillValueInt
        logical                                 :: AddBedRock2Bathym = .false. 
        integer                                 :: direction    = 1    
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
        character(Len = PathLength)           :: ConstructData  = null_str
        character(Len = PathLength)           :: InitialSand    = null_str  
        character(Len = PathLength)           :: OutPutFields   = null_str
        character(Len = PathLength)           :: FinalSand      = null_str
    end type  T_Files

    private :: T_SmoothSlope    
    type       T_SmoothSlope    
        real        :: Critic = FillValueReal
        real        :: Factor = FillValueReal
        logical     :: ON     = .false.   
        logical     :: DZ_Residual  = .false. 
    end type  T_SmoothSlope    
        

    type       T_Discharges
        type(T_Time)                            :: NextCompute
        real                                    :: DT_Compute = FillValueReal
        logical                                 :: Yes        = .false.
    end type T_Discharges 

    type     T_Boxes
        logical                                 :: Yes      = .false.
        character(Len = PathLength)             :: File     =  null_str
        real(8), dimension(:,:), pointer        :: Mass     => null()
        real(8), dimension(:,:), pointer        :: FluxesX  => null()
        real(8), dimension(:,:), pointer        :: FluxesY  => null()
        real(8), dimension(:,:), pointer        :: FluxesZ  => null()        
    end type T_Boxes
    
    type     T_HybridMorph
        logical                                 :: ON                       = .false. 
        !only one boundary can be the coast (1 - ILB_, 2 - IUB_, 3 - JLB_, 4 - JUB_
        integer                                 :: CoastBoundary            = FillValueInt
        !cell where the active profile starts - second water pointstarting from the coast line (CoastBoundary)
        integer, dimension(:),   pointer        :: InShoreMapping           => null() 
        !Offshore depth beyond which the depths do not change with time               
        real                                    :: ClosureDepth           = FillValueReal        
        !cell where the active profiel ends   - first cell with a depth greater than the depth of closure         
        integer, dimension(:),   pointer        :: OffShoreMapping          => null()
        real(8), dimension(:),   pointer        :: AlongShoreFlux           => null()
        real(8), dimension(:),   pointer        :: ResidualAlongShoreFlux   => null()
        real(8), dimension(:),   pointer        :: AlongShoreDZ             => null()
        real(8), dimension(:),   pointer        :: ResidualAlongShoreDZ     => null()        
        integer                                 :: DintegLongShore          = 1
        real(8), dimension(:),   pointer        :: CrossShoreVel            => null()
        real(8), dimension(:),   pointer        :: ResidualCrossShoreVel    => null()        
        real,    dimension(:,:), pointer        :: BathymetryNext           => null()
        real,    dimension(:,:), pointer        :: BathymetryPrevious       => null()
        !Difference between reference bathymetry and bathymetry at instant X (negative - erosion)
        type (T_Property)                       :: DZ_Residual
        !Distance of the cell centers to the domain side that corresponds to the coast boundary
        real(8), dimension(:,:), pointer        :: DistanceToCoastRef       => null()
        !Distance of the cell centers to the domain side that corresponds to the coast 
        ! boundary plus the run period x crossShore profile shift       
        real(8), dimension(:,:), pointer        :: DistanceToCoastInst      => null()  
        
        integer                                 :: Min1D                    = FillValueInt                 
        integer                                 :: Max1D                    = FillValueInt
        integer                                 :: NoTransportBufferCells   = FillValueInt              
    end type T_HybridMorph

    type T_MudErosion
        logical                                 :: ON
        real                                    :: InitialMudContent          = 0.4
        type (T_Property)                       :: MudContent, MudContentAux
        type (T_Property)                       :: Erosion, Residual_Erosion
        real, dimension(:,:), pointer           :: MassFlux                   => null ()        
    end type T_MudErosion
    
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
        type (T_Property)                          :: Uaverage, Vaverage
        type (T_Property)                          :: MappDZ
        logical                                    :: AverageON, MappDZON       
        logical                                    :: ResidualCurrent
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
        type (T_Residual  )                        :: Residual
        type (T_HybridMorph)                       :: HybridMorph
        logical                                    :: BiHarmonicFilter      = .false.
        real                                       :: BiHarmonicFilterCoef  = FillValueReal
        
        real, dimension(:,:), pointer              :: FluxX                 => null ()
        real, dimension(:,:), pointer              :: FluxY                 => null ()
        real, dimension(:,:), pointer              :: FluxZ                 => null ()
        real, dimension(:,:), pointer              :: OutFluxX              => null ()
        real, dimension(:,:), pointer              :: OutFluxY              => null ()
        real, dimension(:,:), pointer              :: FluxXIntegral         => null ()
        real, dimension(:,:), pointer              :: FluxYIntegral         => null ()        
        real, dimension(:,:), pointer              :: FluxZIntegral         => null ()                
        real, dimension(:,:), pointer              :: TransportCapacity     => null ()
        real, dimension(:,:), pointer              :: TauCritic             => null ()
        real, dimension(:,:), pointer              :: Dast                  => null ()
        real, dimension(:,:), pointer              :: BatGradient           => null ()        
        integer                                    :: TransportMethod       = FillValueInt
        logical                                    :: BedLoad               = .true.
        logical                                    :: SuspendedLoad         = .true.        
        logical                                    :: WaveEffect            = .true.                
        logical                                    :: BedSlopeEffects       = .false.         
                       
        type(T_MudErosion)                         :: MudErosion             
                       
        !Instance of ModuleHDF5        
        integer                                    :: ObjHDF5               = 0
        !Instance of ModuleHDF5 initialization
        integer                                    :: ObjHDF5In             = 0        
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
            
            call ConstructAverageCurrent            
            
            call ConstructMappDZ

            call StartOutputBoxFluxes

            call ComputeTauCritic
            
            
            if (Me%HybridMorph%ON) then
                call ConstructHybridMorph
            endif            

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartSand - ModuleSand - ERR30'
            
            if (Me%Evolution%OLD) then
                call KillHDF5 (Me%ObjHDF5In, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'StartSand; ModuleSand - ERR40.'
            endif

            if (Me%OutPut%Yes) call Open_HDF5_OutPut_File(Me%Files%OutPutFields)
            
            if (Me%OutPut%Yes) call OutPutSandHDF            

            if (Me%Discharges%Yes) then
            
                if (ObjDischargesID == 0)  then                                                
                    write(*,*)'You need to define a water discharges in the hydrodynamic input' 
                    stop      'StartSand - Sand - ERR50'
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
            
            stop 'ModuleSand - StartSand - ERR60' 

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
 
    !-----------------------------------------------------------------------------------    

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

            if (Me%MudErosion%ON) then            
            
                allocate(Me%Boxes%FluxesZ(ILB:IUB, JLB:JUB))
                Me%Boxes%FluxesZ(:,:) = 0.            
            
            endif    
                
            allocate(Me%Boxes%Mass   (ILB:IUB, JLB:JUB))
            Me%Boxes%Mass   (:,:) = 0.

        endif i1


    end subroutine StartOutputBoxFluxes


    !--------------------------------------------------------------------------

    subroutine ConstructEvolution

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        real            :: ModelDT
        integer         :: STAT_CALL, iflag
        real(8)         :: ErrorAux, auxFactor, DTaux
        logical         :: EXIST
        integer(4)      :: HDF5_READ        

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
            !Keyword          : DZ_DT
            !<BeginDescription>       
               ! The time step of the DZ evolution
               !
            !<EndDescription>
            !Type             : real 
            !Default          : ModelDT
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Evolution%DZDT,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DZ_DT',                                             &
                     default      = Me%Evolution%SandDT,                                 &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR130' 


        if (Me%Evolution%DZDT < Me%Evolution%SandDT) then

            stop 'ConstructEvolution - ModuleSand - ERR150'

        elseif (Me%Evolution%DZDT > Me%Evolution%SandDT) then

            !Batim DT  must be a multiple of the Sand DT
            auxFactor = Me%Evolution%DZDT / Me%Evolution%SandDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop 'ConstructEvolution - ModuleSand - ERR160'
            endif

        endif

        ! Run period in seconds
        DTaux = Me%EndTime - Me%BeginTime

        !The run period   must be a multiple of the BATIM DT
        auxFactor = DTaux / Me%Evolution%DZDT

        ErrorAux = auxFactor - int(auxFactor)
        if (ErrorAux /= 0) then
            write(*,*) 
            write(*,*) ' Time step error.'
            stop 'ConstructEvolution - ModuleSand - ERR170'
        endif

        Me%Evolution%NextDZ = Me%BeginTime + Me%Evolution%DZDT


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
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR180' 
        
        if (Me%Evolution%OLD) then
        
            inquire (FILE=trim(Me%Files%InitialSand)//"5", EXIST = EXIST)

cd0:        if (EXIST) then

                !Gets File Access Code
                call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

                Me%ObjHDF5In = 0

                !Opens HDF5 File
                call ConstructHDF5 (Me%ObjHDF5In,                                           &
                                    trim(Me%Files%InitialSand)//"5", HDF5_READ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'ReadResidualStartTime; ModuleSand - ERR10.'

            else if(.not. EXIST) then cd0

                stop 'ReadInitialField; ModuleSand - ERR70'

            endif cd0
        
        endif

        call GetData(Me%Aceleration%Coef,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'ACELERATION',                                       &
                     default      = 1.,                                                  &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR190' 

        if (iflag==1) then
            Me%Aceleration%Yes = .true.
        else
            Me%Aceleration%Yes = .false.
        endif


!        call GetData(Me%Residual%ON,                                                    &
!                     Me%ObjEnterData,iflag,                                             &
!                     SearchType   = FromFile,                                           &
!                     keyword      = 'RESIDUAL',                                         &
!                     default      = .false.,                                            &
!                     ClientModule = 'ModuleSand',                                       &
!                     STAT         = STAT_CALL)              
!        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR200' 
        
        if (.not. Me%Evolution%OLD) then
            Me%Residual%StartTime = Me%BeginTime
        endif                

        call GetData(Me%Evolution%Gama,                                                 &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GAMA',                                             &
                     default      = 0.0,                                                &
                     ClientModule = 'ModuleSand',                                       &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR210' 
        
        if (Me%Evolution%Gama < 0 .or. Me%Evolution%Gama > 1) then
            stop 'ConstructEvolution - ModuleSand - ERR220' 
        endif
        ! integer, parameter :: UpwindOrder1 = 1, UpwindOrder2 = 2, UpwindOrder3 = 3, P2_TVD = 4, CentralDif = 5, LeapFrog = 6
        call GetData(Me%Evolution%NumericMethod,                                        &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NUMERIC_METHOD',                                   &
                     default      = UpwindOrder1,                                       &
                     ClientModule = 'ModuleSand',                                       &
                     STAT         = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ConstructEvolution - ModuleSand - ERR230'
        
        
        if (Me%Evolution%NumericMethod /= UpwindOrder1      .and.                       &
            Me%Evolution%NumericMethod /= UpwindOrder2      .and.                       &        
            Me%Evolution%NumericMethod /= CentralDif   ) then
            stop 'ConstructEvolution - ModuleSand - ERR240'
        endif                            
        
        if (Me%Evolution%NumericMethod == CentralDif) then
            if (Me%Evolution%Gama==0) then
                write(*,*) 'Central Difference mass balance and no smoothing = unstable'
                write(*,*) 'ConstructEvolution - ModuleSand - WRN250'
            endif
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
        
        
        if (Me%Evolution%Bathym) then

            call GetGridData2Dreference(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR04' 

        else

            call GetGridData           (Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR04' 
        
        endif

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
        
        if (Me%AverageON) then

            call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%Uaverage%ID%Name),&
                                   trim(Me%Uaverage%ID%Units), Array2D = Me%Uaverage%Field2D,   &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR110'

            call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%Vaverage%ID%Name),&
                                   trim(Me%Vaverage%ID%Units), Array2D = Me%Vaverage%Field2D,   &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR120'
            
         endif  
         
        if (Me%MappDZON) then

            call HDF5WriteData   (Me%ObjHDF5, "/SandCaracteristics", trim(Me%MappDZ%ID%Name),&
                                   trim(Me%MappDZ%ID%Units), Array2D = Me%MappDZ%Field2D,   &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSand - ERR130'

         endif           

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
            stop 'Read_Sand_Files_Name - ModuleSand - ERR10' 


        ! ---> File in HDF format where is written instant fields of Sand properties
        Message   ='Instant fields of sand properties in HDF format.'
        Message   = trim(Message)
        
        !if (GetDDecompON    (Me%ObjHorizontalGrid)) then
        !    write(*,*) 'Module Sand not ready to run in domain decomposition mode'
        !    stop 'Read_Sand_Files_Name - ModuleSand - ERR20' 
        !endif

        call ReadFileName('SAND_OUT', Me%Files%OutPutFields,                             &
                           Message = Message, TIME_END = Me%EndTime,                     &
                           Extension = 'sandlt',                                         &
                           MPI_ID    = GetDDecompMPI_ID(Me%ObjHorizontalGrid),&
                           DD_ON     = GetDDecompON    (Me%ObjHorizontalGrid),&
                           STAT      = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sand_Files_Name - ModuleSand - ERR02' 

        ! ---> Sand properties final values in HDF format
        Message   ='Sand properties final values in HDF format.'
        Message   = trim(Message)
        call ReadFileName('SAND_END', Me%Files%FinalSand,                                &
                           Message = Message, TIME_END = Me%EndTime,                     &
                           Extension = 'sandlf',                                         &
                           MPI_ID    = GetDDecompMPI_ID(Me%ObjHorizontalGrid),           &
                           DD_ON     = GetDDecompON    (Me%ObjHorizontalGrid),           &
                           STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sand_Files_Name - ModuleSand - ERR03' 


        ! ---> Sand properties initial values in HDF format
        Message   ='Sand properties initial values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SAND_INI', Me%Files%InitialSand,                              &
                           Message   = Message, TIME_END = Me%ExternalVar%Now,           &
                           Extension = 'sanlf',                                          &
                           STAT      = STAT_CALL)


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
        character(len=PathLength)                           :: TimeSerieLocationFile
        integer                                             :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber  
        character(len=StringLength)                         :: TimeSerieName
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        type (T_Polygon), pointer                           :: ModelDomainLimit


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
            
        call GetGridOutBorderPolygon(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     Polygon          = ModelDomainLimit,               &
                                     STAT             = STAT_CALL)           
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructTimeSerie - ModuleSand - ERR30' 

        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            trim(TimeSerieLocationFile),                                &
                            PropertyList, "srsand",                                     &
                            WaterPoints2D = Me%ExternalVar%WaterPoints2D,               &
                            ModelDomain   = ModelDomainLimit,                           & 
                            STAT          = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSand - ERR40'
        
        call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,               &
                                 Polygon          = ModelDomainLimit,                   &
                                 STAT             = STAT_CALL)                          
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR50'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSand - ERR60'

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR70'

        do dn = 1, TimeSerieNumber
        
            call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR80'
            
            if (IgnoreOK) cycle        

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR90'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR100'  
                                  
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                
                if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= OUT_OF_BOUNDS_ERR_) then
                    stop 'ConstructTimeSerie - ModuleSand - ERR90'
                endif                            

                !if (STAT_CALL == OUT_OF_BOUNDS_ERR_ .or. Id < 0 .or. Jd < 0) then

                !    call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                !    if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR70'

                !    if (IgnoreOK) then
                !        cycle
                !    else
                !        stop 'ConstructTimeSerie - ModuleSand - ERR80'
                !    endif

                !endif

                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR110'
            endif
            
            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI   = Id,                             &
                                      LocalizationJ   = Jd,                             & 
                                      STAT     = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSand - ERR120'

            if (Me%ExternalVar%WaterPoints2D(Id, Jd) /= WaterPoint) then
                
                 write(*,*) 'Time Serie in a land cell - ',trim(TimeSerieName),' - ',' ModuleSand'

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

    subroutine ConstructSandProperty(NewProperty, ExtractType, ClientNumber)

        !Arguments-------------------------------------------------------------
        type(T_Property)                :: NewProperty
        integer                         :: ExtractType
        integer                         :: ClientNumber

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        
        call ConstructPropertyID(NewProperty%ID, Me%ObjEnterData, ExtractType)

        
        allocate(NewProperty%Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        NewProperty%Field2D(:,:) = 0.


        call ConstructFillMatrix  (PropertyID           = NewProperty%ID,               &
                                   EnterDataID          = Me%ObjEnterData,              &
                                   TimeID               = Me%ObjTime,                   &
                                   HorizontalGridID     = Me%ObjHorizontalGrid,         &
                                   ExtractType          = ExtractType,                  &
                                   PointsToFill2D       = Me%ExternalVar%WaterPoints2D, &
                                   Matrix2D             = NewProperty%Field2D,          &
                                   TypeZUV              = TypeZ_,                       &
                                   ClientID             = ClientNumber,                 & 
                                   STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
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
 
                        call ConstructSandProperty(Me%Classes%Diameter(ClassID), FromBlockInBlock, ClientNumber2)

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

                        call ConstructSandProperty(Me%Classes%Percentage(ClassID), FromBlockInBlock, ClientNumber2)

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
    
    subroutine ConstructAverageCurrent
    
        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL, iflag
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------    
        
        Me%AverageON = .true. 
    
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                    ClientNumber    = ClientNumber,                  &
                                    block_begin     = U_Average_block_begin,         &
                                    block_end       = U_Average_block_end,           &
                                    BlockFound      = BlockFound,                    &
                                    STAT            = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR10'

        if (BlockFound) then
            call ConstructSandProperty(Me%Uaverage, FromBlock, ClientNumber)

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR20'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR30'

        else
            Me%AverageON = .false. 
        endif


        call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                    ClientNumber    = ClientNumber,                  &
                                    block_begin     = V_Average_block_begin,         &
                                    block_end       = V_Average_block_end,           &
                                    BlockFound      = BlockFound,                    &
                                    STAT            = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR40'

        if (BlockFound) then
            call ConstructSandProperty(Me%Vaverage, FromBlock, ClientNumber)
        else
            Me%AverageON = .false. 
        endif

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR50'
        
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR60' 
    
        !<BeginKeyword>
            !Keyword          : RESIDUAL_CURRENT
            !<BeginDescription>       
               ! If in transport is used the residual current for define the transport direction
               !
            !<EndDescription>
            !Type             : logical
            !Default          : .false. 
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%ResidualCurrent,                                                 &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'RESIDUAL_CURRENT',                                  &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL /= SUCCESS_) stop 'ConstructAverageCurrent - ModuleSand - ERR80' 
        

        if (Me%AverageON .and. Me%ResidualCurrent) then
            stop 'ConstructAverageCurrent - ModuleSand - ERR90' 
        endif

    end subroutine ConstructAverageCurrent


    !----------------------------------------------------------------------------
    
    subroutine ConstructMappDZ
    
        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------    
        
        Me%MappDZON = .true. 
    
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                    ClientNumber    = ClientNumber,                  &
                                    block_begin     = Mapp_DZ_block_begin,           &
                                    block_end       = Mapp_DZ_block_end,             &
                                    BlockFound      = BlockFound,                    &
                                    STAT            = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ConstructMappDZ - ModuleSand - ERR10'

        if (BlockFound) then
            call ConstructSandProperty(Me%MappDZ, FromBlock, ClientNumber)

        else
            Me%MappDZON = .false. 
        endif

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructMappDZ - ModuleSand - ERR20'

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructMappDZ - ModuleSand - ERR30'


    end subroutine ConstructMappDZ

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
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR05'
        
        !<BeginKeyword>
            !Keyword          : ADD_BEDROCK_2_BATHYM
            !<BeginDescription>       
               ! Check if the user wants to add the bedrock to the batymetry
               !
            !<EndDescription>
            !Type             : logical
            !Default          : .true.
            !File keyword     : SAND_DATA
            !Search Type      : FromFile
        !<EndKeyword>

        call GetData(Me%Evolution%AddBedRock2Bathym,                                     &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'ADD_BEDROCK_2_BATHYM',                              &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR08'        
                

        if (Me%Classes%Number == 0) then

            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber    = ClientNumber,                  &
                                        block_begin     = D35_block_begin,               &
                                        block_end       = D35_block_end,                 &
                                        BlockFound      = BlockFound,                    &
                                        STAT            = STAT_CALL)
  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR10'

            if (BlockFound) then
                call ConstructSandProperty(Me%D35, FromBlock, ClientNumber)

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR20'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR30'

            else
                stop 'ConstructGlobalParameters - ModuleSand - ERR40'
            endif

            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber    = ClientNumber,                  &
                                        block_begin     = D50_block_begin,               &
                                        block_end       = D50_block_end,                 &
                                        BlockFound      = BlockFound,                    &
                                        STAT            = STAT_CALL)
  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR50'

            if (BlockFound) then
                call ConstructSandProperty(Me%D50, FromBlock, ClientNumber)

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL .NE. SUCCESS_) stop 'CConstructGlobalParameters - ModuleSand - ERR60'

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR70'

            else
                stop 'ConstructGlobalParameters - ModuleSand - ERR04'
            endif


            call ExtractBlockFromBuffer(Me%ObjEnterData,                                 &
                                        ClientNumber    = ClientNumber,                  &
                                        block_begin     = D90_block_begin,               &
                                        block_end       = D90_block_end,                 &
                                        BlockFound      = BlockFound,                    &
                                        STAT            = STAT_CALL)
  
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR80'

            if (BlockFound) then
                call ConstructSandProperty(Me%D90, FromBlock, ClientNumber)

                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR90'

                if (STAT_CALL .NE. SUCCESS_) stop 'CConstructGlobalParameters - ModuleSand - ERR100'
            else
                stop 'ConstructGlobalParameters - ModuleSand - ERR110'
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
  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR120'

        if (BlockFound) then
            call ConstructSandProperty(Me%BedRock, FromBlock, ClientNumber)

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'CConstructGlobalParameters - ModuleSand - ERR130'

            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR140'
            
        else
            stop 'ConstructGlobalParameters - ModuleSand - ERR180'
        endif


        allocate(Me%DZ%Field2D         (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        Me%DZ%ID%Name  = 'DZ'
        Me%DZ%ID%Units = 'm'

        
        Me%DZ%Field2D(:,:) = 0.

!        if (Me%Evolution%BathymDT > Me%Evolution%DZDT) then
!
!            allocate(Me%BatimIncrement%Field2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
!
!            Me%BatimIncrement%ID%Name  = 'Batim Increment'
!            Me%BatimIncrement%ID%Units = 'm'
!        
!            Me%BatimIncrement%Field2D(:,:) = 0.
!        else
!
!            Me%BatimIncrement%Field2D => Me%DZ%Field2D
!
!        endif


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
                     default      = 0.0001,                                              &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR190'


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
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR200'
        
        if (Me%Porosity < 0. .or. Me%Porosity > 1.) then
            stop 'ConstructGlobalParameters - ModuleSand - ERR210' 
        endif


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
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR220'
        
        Me%RhoSl           = Me%Density - Me%ExternalVar%WaterDensity

        Me%RelativeDensity = Me%RhoSl / Me%ExternalVar%WaterDensity

        !Allocate fluxes
        allocate(Me%FluxX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%FluxX(:,:) = 0.

        allocate(Me%FluxY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%FluxY(:,:) = 0.
        
        allocate(Me%OutFluxX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%OutFluxX(:,:) = 0.

        allocate(Me%OutFluxY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%OutFluxY(:,:) = 0.
        
        
        allocate(Me%FluxXIntegral(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%FluxXIntegral(:,:) = 0.

        allocate(Me%FluxYIntegral(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%FluxYIntegral(:,:) = 0.


        allocate(Me%Residual%FluxX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%Residual%FluxX(:,:) = 0.

        allocate(Me%Residual%FluxY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%Residual%FluxY(:,:) = 0.

        allocate(Me%Residual%OutFluxX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%Residual%OutFluxX(:,:) = 0.

        allocate(Me%Residual%OutFluxY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%Residual%OutFluxY(:,:) = 0.
        
        if (Me%MudErosion%ON) then
            !Allocate Mud fluxes
            allocate(Me%FluxZ(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%FluxZ(:,:) = 0.


            allocate(Me%FluxZIntegral(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%FluxZIntegral(:,:) = 0.

            allocate(Me%Residual%FluxZ(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            Me%Residual%FluxZ(:,:) = 0.

        endif
        
        if (Me%Evolution%Old) then
        
            call ReadResidualStartTime()
            
            call ReadInitialField(FieldName = "Residual Transport Flux X", Field2D = Me%Residual%OutFluxX)
            
            call ReadInitialField(FieldName = "Residual Transport Flux Y", Field2D = Me%Residual%OutFluxY)    
            
            call RotateVectorFieldToGrid(HorizontalGridID  = Me%ObjHorizontalGrid,      &
                                         VectorInX         = Me%Residual%OutFluxX,      &
                                         VectorInY         = Me%Residual%OutFluxY,      &
                                         VectorOutX        = Me%Residual%FluxX,         &
                                         VectorOutY        = Me%Residual%FluxY,         &
                                         WaterPoints2D     = Me%ExternalVar%WaterPoints2D,&
                                         RotateX           = .true.,                    &
                                         RotateY           = .true.,                    &
                                         STAT              = STAT_CALL)
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR225'
                        
            if (Me%MudErosion%ON) then

                call ReadInitialField(FieldName = "Residual Mud Erosion", Field2D = Me%Residual%FluxZ)
                
                allocate(Me%MudErosion%MassFlux(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
                Me%MudErosion%MassFlux(:,:) = 0.
                
            endif    
        
        endif
        

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
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR230'

        SELECT CASE (trim(Auxchar))

        Case ("no transport")
            Me%TransportMethod = NoTransport
        Case ("MeyerPeter")
            Me%TransportMethod = MeyerPeter
        Case ("Ackers")
            Me%TransportMethod = Ackers
        Case ("VanRijn2007")
            Me%TransportMethod = VanRijn2007
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
            stop 'ConstructGlobalParameters - ModuleSand - ERR240'
                
        END SELECT 

        if     ( Me%TransportMethod == VanRijn1                                         &
            .OR. Me%TransportMethod == VanRijn2                                         &
            .OR. Me%TransportMethod == Bijker)  then
    !.OR. Me%TransportMethod == VanRijn2 .OR. Me%TransportMethod == Bijker
     
            if (Me%ObjWaves == 0) stop 'ConstructGlobalParameters - ModuleSand - ERR250'

            if (.not. Me%ExternalVar%WaveTensionON) stop 'ConstructGlobalParameters - ModuleSand - ERR260'

        endif

        allocate (Me%Dast(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        Me%Dast(:,:) = FillValueReal


        call GetData(AuxChar,                                                            &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'FILTER_SCHEME',                                     &
                     ClientModule = 'ModuleSand',                                        &
                     Default      = 'No Filter',                                         &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR280'

        Me%Filter%ON = .true.

        Select Case (trim(AuxChar))

        Case ("NO FILTER", "no fielter", "No filter", "No Filter")
            Me%Filter%Scheme = NoFilter
            Me%Filter%ON     = .false.
        Case ("MODIFY LAX", "modify lax", "Modify Lax", "Modify lax")
            Me%Filter%Scheme = ModifyLax
        Case default
            stop 'ConstructGlobalParameters - ModuleSand - ERR290'
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
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR300'

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
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR310'


        if (Me%SmoothSlope%ON) then

            !Slope be on each exist flux perpendicular to the Slope
            call GetData(Me%SmoothSlope%Critic,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'CRITICAL_SLOPE',                               &
                         ClientModule = 'ModuleSand',                                   &
                         Default      = 0.1,                                            &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR320'

            !The flux perpendicular to the flux is a percentage of the paralel flux
            call GetData(Me%SmoothSlope%Factor,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'FLUX_SLOPE',                                   &
                         ClientModule = 'ModuleSand',                                   &
                         Default      = 0.1,                                            &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR330'

            !Smooth also dynamic componente DZ_Residual
            !If connected can be inconsistent option in areas with steep bathymetry
            call GetData(Me%SmoothSlope%DZ_Residual,                                    &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'DZ_RESIDUAL_SLOPE',                            &
                         ClientModule = 'ModuleSand',                                   &
                         Default      = .false.,                                        &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR330'            

        endif

        call GetData(Me%Boundary,                                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BOUNDARY',                                          &
                     ClientModule = 'ModuleSand',                                        &
                     Default      = NullGradient,                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR340'

        call GetData(Me%TransportFactor,                                                 &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TRANSPORT_FACTOR',                                  &
                     default      = 1.,                                                  &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR350'


        call GetData(Me%TauMax,                                                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TAU_MAX',                                           &
                     default      = 10.,                                                 &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR360'


        call GetData(Me%Discharges%Yes,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DISCHARGES',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR370'


        if (Me%TransportMethod == VanRijn2007) then
        
            call GetData(Me%BedLoad,                                                    &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'BEDLOAD',                                      &
                         default      = .true.,                                         &
                         ClientModule = 'ModuleSand',                                   &
                         STAT         = STAT_CALL)              
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR380'

            call GetData(Me%SuspendedLoad,                                              &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'SUSPENDEDLOAD',                                &
                         default      = .true.,                                         &
                         ClientModule = 'ModuleSand',                                   &
                         STAT         = STAT_CALL)              
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR390'

            call GetData(Me%WaveEffect,                                                 &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'WAVE_EFFECT',                                  &
                         default      = .true.,                                         &
                         ClientModule = 'ModuleSand',                                   &
                         STAT         = STAT_CALL)              
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR400'

        endif    
        
         call GetData(Me%BedSlopeEffects,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BEDSLOPE',                                          &
                     default      = .TRUE.,                                              &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR410'
         

         call GetData(Me%HybridMorph%ON,                                                 &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'HYBRID_MORPH',                                      &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR420'
         
         if (Me%HybridMorph%ON .and. .not. Me%Evolution%Bathym) then
            stop 'ConstructGlobalParameters - ModuleSand - ERR430'
         endif
         
!         if (Me%HybridMorph%ON) then
!            Me%Evolution%Bathym = .true.
!         endif         

         call GetData(Me%BiHarmonicFilter,                                              &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BIHARMONIC_FILTER',                                &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleSand',                                       &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR440'


         call GetData(Me%BiHarmonicFilterCoef,                                          &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BIHARMONIC_FILTER_COEF',                           &
                     default      = 0.01,                                               &
                     ClientModule = 'ModuleSand',                                       &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR450'

    end subroutine ConstructGlobalParameters

    !----------------------------------------------------------------------------
    
    subroutine ConstructHybridMorph

        !Local-------------------------------------------------------------------
        integer                             :: iflag, STAT_CALL
        integer                             :: ILB, IUB, JLB, JUB, Min1D, Max1D, i, j

        !------------------------------------------------------------------------
    
        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 

        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

         call GetData(Me%HybridMorph%CoastBoundary,                                      &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'HYBRID_MORPH_COAST_BOUND',                          &
                     default      = ILB_,                                                &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructHybridMorph - ModuleSand - ERR10'    
         
         if     (Me%HybridMorph%CoastBoundary == ILB_ .or. Me%HybridMorph%CoastBoundary == IUB_ ) then

            Min1D = JLB
            Max1D = JUB
            
            Me%HybridMorph%Min1D = Min1D+1
            Me%HybridMorph%Max1D = Max1D-1
            
         elseif (Me%HybridMorph%CoastBoundary == JLB_ .or. Me%HybridMorph%CoastBoundary == JUB_ ) then            

            Min1D = ILB
            Max1D = IUB
            
            Me%HybridMorph%Min1D = Min1D+1
            Me%HybridMorph%Max1D = Max1D-1
            
         endif
         
         call GetData(Me%HybridMorph%ClosureDepth,                                       &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CLOSURE_DEPTH',                                     &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructHybridMorph - ModuleSand - ERR20'    
         
         if (iflag == 0) then
            write(*,*) 'Need to define the CLOSURE_DEPTH'
            write(*,*) 'Offshore depth beyond which the depths do not change with time'
            stop 'ConstructHybridMorph - ModuleSand - ERR30'
         endif

         call GetData(Me%HybridMorph%DintegLongShore,                                    &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CELLS_INTEGRATION_ALONG_SHORE',                     &
                     default      = 1,                                                   &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructHybridMorph - ModuleSand - ERR40'    

         call GetData(Me%HybridMorph%NoTransportBufferCells,                             &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'NO_TRANSPORT_BUFFER_CELLS',                         &
                     default      = FillValueInt,                                        &
                     ClientModule = 'ModuleSand',                                        &
                     STAT         = STAT_CALL)              
         if (STAT_CALL /= SUCCESS_) stop 'ConstructHybridMorph - ModuleSand - ERR50'    

        if (2 * Me%HybridMorph%NoTransportBufferCells >= Me%HybridMorph%Max1D - Me%HybridMorph%Min1D) then
            stop 'ConstructHybridMorph - ModuleSand - ERR60'
        endif         

        !Allocate variables
        !1D vectors
        allocate(Me%HybridMorph%InShoreMapping        (Min1D:Max1D))
        allocate(Me%HybridMorph%OffShoreMapping       (Min1D:Max1D)) 
        allocate(Me%HybridMorph%AlongShoreFlux        (Min1D:Max1D)) 
        allocate(Me%HybridMorph%ResidualAlongShoreFlux(Min1D:Max1D)) 
        allocate(Me%HybridMorph%AlongShoreDZ          (Min1D:Max1D)) 
        allocate(Me%HybridMorph%ResidualAlongShoreDZ  (Min1D:Max1D)) 
        allocate(Me%HybridMorph%CrossShoreVel         (Min1D:Max1D)) 
        allocate(Me%HybridMorph%ResidualCrossShoreVel (Min1D:Max1D)) 

        Me%HybridMorph%InShoreMapping        (Min1D:Max1D) = FillValueInt
        Me%HybridMorph%OffShoreMapping       (Min1D:Max1D) = FillValueInt
        
        Me%HybridMorph%AlongShoreFlux        (Min1D:Max1D) = 0.
        Me%HybridMorph%ResidualAlongShoreFlux(Min1D:Max1D) = 0.
        Me%HybridMorph%AlongShoreDZ          (Min1D:Max1D) = 0.
        Me%HybridMorph%ResidualAlongShoreDZ  (Min1D:Max1D) = 0.
        Me%HybridMorph%CrossShoreVel         (Min1D:Max1D) = 0.
        Me%HybridMorph%ResidualCrossShoreVel (Min1D:Max1D) = 0.


        allocate(Me%HybridMorph%BathymetryNext        (ILB:IUB, JLB:JUB))
        allocate(Me%HybridMorph%BathymetryPrevious    (ILB:IUB, JLB:JUB))
        allocate(Me%HybridMorph%DistanceToCoastRef    (ILB:IUB, JLB:JUB))
        allocate(Me%HybridMorph%DistanceToCoastInst   (ILB:IUB, JLB:JUB))           
        
        Me%HybridMorph%BathymetryNext        (ILB:IUB, JLB:JUB) = FillValueReal
        Me%HybridMorph%BathymetryPrevious    (ILB:IUB, JLB:JUB) = FillValueReal
        Me%HybridMorph%DistanceToCoastRef    (ILB:IUB, JLB:JUB) = FillValueReal
        Me%HybridMorph%DistanceToCoastInst   (ILB:IUB, JLB:JUB) = FillValueReal
                     
        allocate(Me%HybridMorph%DZ_Residual%Field2D(ILB:IUB, JLB:JUB))
        
        Me%HybridMorph%DZ_Residual%Field2D(ILB:IUB, JLB:JUB) = 0.

        Me%HybridMorph%DZ_Residual%ID%Name  = 'DZ_Residual'
        Me%HybridMorph%DZ_Residual%ID%Units = 'm'
                        
        if      (Me%HybridMorph%CoastBoundary == ILB_) then
         
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB 

                Me%HybridMorph%DistanceToCoastRef(Me%WorkSize%ILB,j) = 0.

                do i = Me%WorkSize%ILB+1, Me%WorkSize%IUB 
             
                    Me%HybridMorph%DistanceToCoastRef(i,j) = Me%HybridMorph%DistanceToCoastRef(i-1,j) + Me%ExternalVar%DZY(i-1,j)
                    
                enddo
            enddo
                
         elseif (Me%HybridMorph%CoastBoundary == IUB_ ) then


            do j = Me%WorkSize%JLB,   Me%WorkSize%JUB 
            
                Me%HybridMorph%DistanceToCoastRef(Me%WorkSize%IUB,j) = 0.

                do i = Me%WorkSize%IUB-1, Me%WorkSize%ILB, -1
             
                    Me%HybridMorph%DistanceToCoastRef(i,j) = Me%HybridMorph%DistanceToCoastRef(i+1,j) + Me%ExternalVar%DZY(i  ,j)
                    
                enddo
            enddo
                            
         elseif (Me%HybridMorph%CoastBoundary == JLB_) then

            do i = Me%WorkSize%ILB,   Me%WorkSize%IUB 
         
                Me%HybridMorph%DistanceToCoastRef(i,Me%WorkSize%JLB) = 0.

                do j = Me%WorkSize%JLB+1, Me%WorkSize%JUB 
             
                    Me%HybridMorph%DistanceToCoastRef(i,j) = Me%HybridMorph%DistanceToCoastRef(i,j-1) + Me%ExternalVar%DZX(i,j-1)
                    
                enddo
            enddo
                         
         elseif (Me%HybridMorph%CoastBoundary == JUB_ ) then


            do i = Me%WorkSize%ILB,   Me%WorkSize%IUB 

                Me%HybridMorph%DistanceToCoastRef(i,Me%WorkSize%JUB) = 0.

                do j = Me%WorkSize%JUB-1, Me%WorkSize%JLB, -1
            
                    Me%HybridMorph%DistanceToCoastRef(i,j) = Me%HybridMorph%DistanceToCoastRef(i,j+1) + Me%ExternalVar%DZX(i  ,j)
                    
                enddo
            enddo

         endif        
         
        Me%HybridMorph%DistanceToCoastInst(:,:) = Me%HybridMorph%DistanceToCoastRef(:,:)
        
        call GetGridData2Dreference(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHybridMorph - ModuleSand - ERR70' 
        
        Me%HybridMorph%BathymetryPrevious(:,:) = Me%ExternalVar%InitialBathym(:,:)
            
        
        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHybridMorph - ModuleSand - ERR80'         

    end subroutine ConstructHybridMorph
            
    !--------------------------------------------------------------------------
    !If the user want's to use the values of a previous   
    ! run the read the sand properties values form the final      
    ! results file of a previous run. By default this      
    ! file is in HDF format                                

    subroutine ReadInitialField(FieldName, Field2D)

        !Arguments-------------------------------------------------------------
        character (Len = *)                         :: FieldName
        real, dimension(:,:), pointer               :: Field2D

        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer               :: Aux2D
        integer                                     :: STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: ILW, IUW, JLW, JUW
        type (T_Size2D)                             :: WindowLimitsJI
        logical                                     :: MasterOrSlave          

        !----------------------------------------------------------------------
        
        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 

        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 

                        
        call GetDDecompParameters(HorizontalGridID = Me%ObjHorizontalGrid,          &
                                  MasterOrSlave    = MasterOrSlave,                 &
                                  STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadInitialField; ModuleSand - ERR20'
        endif

            
ifMS:   if (MasterOrSlave) then
                
            call GetDDecompWorkSize2D(HorizontalGridID = Me%ObjHorizontalGrid,      &
                                      WorkSize         = WindowLimitsJI,            &
                                      STAT             = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ReadInitialField; ModuleSand - ERR30'
            endif

            
            ILW = WindowLimitsJI%ILB
            IUW = WindowLimitsJI%IUB

            JLW = WindowLimitsJI%JLB
            JUW = WindowLimitsJI%JUB
            
            write(*,*) 'MPI 1', ILW, IUW, JLW, JUW
                                                  
        else ifMS

            ILW = ILB 
            IUW = IUB

            JLW = JLB 
            JUW = JUB 

        endif ifMS                
            
        !Reads from HDF file the Property concentration and open boundary values
        call HDF5SetLimits  (Me%ObjHDF5In, ILW, IUW, JLW, JUW, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadInitialField; ModuleSand - ERR40'
        endif

        allocate(Aux2D(ILW:IUW,JLW:JUW))
        
        call HDF5ReadWindow (HDF5ID         = Me%ObjHDF5In,                             &
                             GroupName      = "/Results",                               &
                             Name           = trim(FieldName),                          &
                             Array2D        = Aux2D,                                    &
                             STAT           = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadInitialField; ModuleSand - ERR50'
        
        Field2D(ILB:IUB, JLB:JUB) = Aux2D(ILW:IUW, JLW:JUW)
       
        deallocate(Aux2D)                

        !----------------------------------------------------------------------

    end subroutine ReadInitialField

    !--------------------------------------------------------------------------


    subroutine ReadResidualStartTime()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                    :: STAT_CALL
        real,    dimension(:    ), pointer         :: AuxTime

        !----------------------------------------------------------------------


        call HDF5SetLimits  (Me%ObjHDF5In, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadResidualStartTime - ModuleSand - ERR30'
        
        allocate(AuxTime(6))

        call HDF5ReadData  (Me%ObjHDF5In, "/Time",                                      &
                            "Residual Start Time",                                      &
                             Array1D = AuxTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadResidualStartTime - ModuleSand - ERR40'
        
        call SetDate   (Me%Residual%StartTime, AuxTime(1), AuxTime(2), AuxTime(3),      &
                        AuxTime(4), AuxTime(5), AuxTime(6))
                  
        deallocate(AuxTime)
                
        !----------------------------------------------------------------------

    end subroutine ReadResidualStartTime

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
                          ShearVelocity, MinWaterColumn, VelU_Face, VelV_Face,  &
                          VelResidualU, VelResidualV, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSandID
        real, dimension(:,:), pointer               :: TauTotal, CurrentRugosity, WaveRugosity, &
                                                       WaterColumn, VelU, VelV, VelMod,         &
                                                       TauWave, TauCurrent, ShearVelocity,      &
                                                       VelU_Face, VelV_Face, VelResidualU, VelResidualV
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
            
            if (Me%Evolution%Bathym .and. Me%Evolution%AddBedRock2Bathym) then

                !Bathymetry 
                call UnGetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR150'

                call ModifyGridData(Me%ObjBathym, Me%BedRock%Field2D, Add = .false.,  &
                                    STAT = STAT_CALL)
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR160'

                !Bathymetry
                call GetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                if (STAT_CALL  /= SUCCESS_) stop 'ConstructGlobalParameters - ModuleSand - ERR170'
                
                Me%Evolution%AddBedRock2Bathym = .false. 
                
            endif             

            if (Me%TransportMethod /= NoTransport) then

                do while (Me%ExternalVar%Now >= Me%Evolution%NextSand) 

                    Me%ExternalVar%TauTotal        => TauTotal

                    Me%ExternalVar%CurrentRugosity => CurrentRugosity
                    Me%ExternalVar%WaveRugosity    => WaveRugosity

                    Me%ExternalVar%WaterColumn     => WaterColumn

                    Me%ExternalVar%VelU            => VelU
                    Me%ExternalVar%VelV            => VelV
                    Me%ExternalVar%VelMod          => VelMod
                    Me%ExternalVar%VelU_Face       => VelU_Face
                    Me%ExternalVar%VelV_Face       => VelV_Face

                    Me%ExternalVar%VelResidualU    => VelResidualU
                    Me%ExternalVar%VelResidualV    => VelResidualV

                    Me%ExternalVar%TauWave         => TauWave
                    Me%ExternalVar%TauCurrent      => TauCurrent
                    Me%ExternalVar%ShearVelocity   => ShearVelocity

                    Me%ExternalVar%MinWaterColumn  =  MinWaterColumn
                    
                    call ComputeFluxes
                    
                    if (Me%MudErosion%ON) then
                        call ComputeMudErosion
                    endif    
                    
                    if (Me%ExternalVar%Now >= Me%Evolution%NextDZ) then

                        if (Me%Evolution%NumericMethod     == UpwindOrder1 ) then
                            call ComputeEvolutionX            
                        elseif (Me%Evolution%NumericMethod == CentralDif   ) then
                            call ComputeEvolution_CentralDifferences
                        elseif (Me%Evolution%NumericMethod == UpwindOrder2 ) then
                            call ComputeEvolution_UpwindOrder2                            
                        endif                            
                        
                        call ComputeNumericalSmoothing                        

                        if (Me%SmoothSlope%ON) then

                            call ComputeSmoothSlope

                        endif

                        !Boundary Condition
                        !call BoundaryCondition(Me%DZ%Field2D)
                        call BoundaryCondition(Me%DZ_Residual%Field2D)
                        
                        !Mapp effect
                        if (Me%MappDZON) then
                            call ModifyMappDZ
                        endif
                        
                        
#if _USE_MPI                    
                        !MPI and Domain Decomposition is ON exchanges data along domain interfaces
                        call ReceiveSendProperitiesMPI(HorizontalGridID = Me%ObjHorizontalGrid, & 
                                                       Property2D       = Me%DZ%Field2D,    &
                                                       STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            stop 'ModifySand - ModuleSand - ERR10'
                        endif                        
#endif _USE_MPI

                        if (Me%HybridMorph%ON) then
                            call ComputeHybridMorphEvolution
                        endif
                        

                        call ComputeResidualEvolution
                        
                        Me%FluxXIntegral (:,:) = 0.
                        Me%FluxYIntegral (:,:) = 0.

                        
                        if (Me%MudErosion%ON) then
                            Me%FluxZIntegral (:,:) = 0.
                        endif    

                        Me%Evolution%NextDZ = Me%Evolution%NextDZ + Me%Evolution%DZDT

                    endif                         

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
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR20'
                            
                            if (Me%HybridMorph%ON) then                            
                                call ModifyGridData(Me%ObjBathym, Me%HybridMorph%DZ_Residual%Field2D, &
                                                    Add = .false., ResidualIncrement = .true., &
                                                    STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR25.'
                            else                                
    !                            call ModifyGridData(Me%ObjBathym, Me%BatimIncrement%Field2D, Add = .false.,  &
                                call ModifyGridData(Me%ObjBathym, Me%DZ_Residual%Field2D, Add = .false.,  &
                                                    ResidualIncrement = .true., STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR30.'
                                
                            endif                                

                            !Bathymetry
                            call GetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySand - ModuleSand - ERR40'

                            !Me%BatimIncrement%Field2D(:,:) = 0.

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

    subroutine ModifyMappDZ 


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                            :: i, j
        !----------------------------------------------------------------------
      
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint .and. Me%MappDZ%Field2D(i,j) < 0.5) then

                Me%DZ%Field2D(i, j) = 0. 

            endif
                                   
        enddo
        enddo

    end subroutine ModifyMappDZ
    
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine FilterModifyLax 


        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real                               :: AuxSum, Beta
        integer                            :: i, j, iw, jw, Counter
        !----------------------------------------------------------------------


        Beta = .1
      
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                Counter = 0
                AuxSum  = 0.
                do jw = j - Me%Filter%Radius, j + Me%Filter%Radius
                do iw = i - Me%Filter%Radius, i + Me%Filter%Radius
                
                    if (i == iw .and. j == jw) cycle

                    if (jw >= Me%WorkSize%JLB .and. jw <= Me%WorkSize%JUB .and.         &
                        iw >= Me%WorkSize%ILB .and. iw <= Me%WorkSize%IUB) then
                        if(Me%ExternalVar%WaterPoints2D(iw, jw) == WaterPoint) then
                            Counter = Counter + 1
                            !AuxSum  = AuxSum + Me%BatimIncrement%Field2D(iw, jw)
                            AuxSum  = AuxSum + Me%DZ_Residual%Field2D(iw, jw)                            
                        end if
                    endif

                enddo
                enddo

                !Me%Filter%Field2D(i, j) = Beta * (AuxSum / real(Counter) - Me%BatimIncrement%Field2D(i,j))
                Me%Filter%Field2D(i, j) = Beta * (AuxSum / real(Counter) - Me%DZ_Residual%Field2D(i,j))

            endif
                                   
        enddo
        enddo

        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                !Me%BatimIncrement%Field2D(i,j) = Me%BatimIncrement%Field2D(i,j) + Me%Filter%Field2D(i,j)
                Me%DZ_Residual%Field2D(i,j) = Me%DZ_Residual%Field2D(i,j) + Me%Filter%Field2D(i,j)
                
            endif
                                   
        enddo
        enddo
       

    end subroutine FilterModifyLax
    
    !--------------------------------------------------------------------------


    subroutine ComputeFluxes
        !Local-----------------------------------------------------------------
        integer             :: i, j
        real                :: Xaux, Yaux, XYMod, FluxX, FluxY, DT_Racio
        real                :: Xcomp, Ycomp
        !----------------------------------------------------------------------

        SELECT CASE (Me%TransportMethod)

        Case (MeyerPeter)
            call MeyerPeterTransport
        Case (Ackers)
            call AckersTransport
        Case (VanRijn2007) 
            call VanRijn2007Transport            
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

        !Computes the sand fluxes (m3/s) in the middle of the cells
        
        if (Me%ResidualCurrent) then
            if (.not.associated(Me%ExternalVar%VelResidualU)) then
                stop 'ComputeFluxes - ModuleSand - ERR10.'
            endif

            if (.not.associated(Me%ExternalVar%VelResidualV)) then
                stop 'ComputeFluxes - ModuleSand - ERR20.'
            endif
                            
        endif

        Me%FluxX  (:, :) =  0.
        Me%FluxY  (:, :) =  0.

!Compute FluxX and FluxY
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

                !SandThickness = Me%BedRock%Field2D(i, j) + Me%DZ_Residual%Field2D(i, j)
                    
                !if (SandThickness > Me%SandMin) then

                    XYMod = Me%ExternalVar%VelMod(i,j) 
                
                    if (Me%AverageON) then

                        Xaux  = Me%Uaverage%Field2D (i,j)
                        Yaux  = Me%Vaverage%Field2D (i,j)

                    else 
                    
                        if (Me%ResidualCurrent) then
                        
                            Xaux  = Me%ExternalVar%VelResidualU  (i,j)
                            Yaux  = Me%ExternalVar%VelResidualV  (i,j)
                        
                        else

                            Xaux  = Me%ExternalVar%VelU  (i,j)
                            Yaux  = Me%ExternalVar%VelV  (i,j)
                            
                        endif                            
                        
                    endif    
                    
                    if (XYMod > 0.) then
                    
                        Xcomp = Xaux / XYMod
!                        if (Xcomp < 1e-3) then
!                            Xcomp = 0.
!                            Ycomp = 1.    
!                        endif   

                        Ycomp = Yaux / XYMod
!                        if (Ycomp < 1e-3) then
!                            Ycomp = 0.
!                            Xcomp = 1.    
!                        endif                        
!                    
                        ![m]      = [m]/[m]*[m]
                        FluxX     = Xcomp * Me%ExternalVar%DVY(i, j) / (1. - Me%Porosity)
                        FluxY     = Ycomp * Me%ExternalVar%DUX(i, j) / (1. - Me%Porosity)
                        
                        ![m^3/s]        = [m^2/s] * [m]
                        Me%FluxX (i, j) = Me%TransportCapacity(i, j) * FluxX  * Me%TransportFactor
                        Me%FluxY (i, j) = Me%TransportCapacity(i, j) * FluxY  * Me%TransportFactor
                        
                !endif

                endif

            endif

        enddo
        enddo
        
        if(Me%BedSlopeEffects)then
            call ComputeBedSlopeEffects
        endif
        
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

                DT_Racio = Me%Evolution%SandDT / Me%Evolution%DZDT
                
                Me%FluxXIntegral (i, j) = Me%FluxXIntegral (i, j) + Me%FluxX (i, j) * DT_Racio
                Me%FluxYIntegral (i, j) = Me%FluxYIntegral (i, j) + Me%FluxY (i, j) * DT_Racio

            endif

        enddo
        enddo        

    end subroutine ComputeFluxes
    !--------------------------------------------------------------------------
    
   !--------------------------------------------------------------------------


    subroutine ComputeMudErosion
        !Local-----------------------------------------------------------------
        integer             :: i, j
        real                :: Area, Aux, TauCritic, DT_Racio
        !----------------------------------------------------------------------


        !Computes the mud erosion fluxes (m3/s) in the middle of the cells
        
        if (Me%ResidualCurrent) then
            stop 'ComputeFluxes - ModuleSand - ERR10.'
                            
        endif

        Me%FluxZ  (:, :) =  0.


!Compute FluxZ
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
                !if (Me%MudErosion%CohesiveBed) then
                !    BulkDensity = 
                !    TauCritic = 0.015 * ()
                !else
                    TauCritic = (1. + Me%MudErosion%MudContent%Field2D(i, j))**2.2*Me%TauCritic(i,j)
                !endif
                if (Me%ExternalVar%TauTotal(i,j) > TauCritic) then
                    ![kg/m2/s]
                    Me%MudErosion%MassFlux(i, j) = 5e-5 * (Me%ExternalVar%TauTotal(i,j)/TauCritic-1.) * &
                                                           Me%MudErosion%MudContent%Field2D(i, j)
                else
                    Me%MudErosion%MassFlux(i, j) = 0.
                endif     
                
                ![m2]   = [m]*[m]
                Area    = Me%ExternalVar%DVY(i, j)  * Me%ExternalVar%DUX(i, j) 
                ![m2]   = [m2] / [-] 
                Aux     = Area  / (1. - Me%Porosity)
                    
                ![m^3/s]        = [kg/m2/s] * [m2] / ([kg/m3] 
                Me%FluxZ (i, j) = Me%MudErosion%MassFlux(i, j) * Aux  / Me%Density 

            endif

        enddo
        enddo

        
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

                DT_Racio = Me%Evolution%SandDT / Me%Evolution%DZDT
                
                Me%FluxZIntegral (i, j) = Me%FluxZIntegral (i, j) + Me%FluxZ (i, j) * DT_Racio

            endif

        enddo
        enddo        

    end subroutine ComputeMudErosion
    !--------------------------------------------------------------------------
    
   !--------------------------------------------------------------------------

    subroutine ComputeBedSlopeEffects
      
        !Local-----------------------------------------------------------------                
        integer                 :: i, j
        real(8)                 :: Xaux, Yaux, AbsFlux, AuxFluxX, AuxFluxY
        real(8), parameter      :: PI_DBLE = 3.1415926536 !PI
        real                    :: alfa_bs, alfa_bn, phi, dhdx, dhdy
        real                    :: dzds, dzdn, alfa_s, alfa_n
        !----------------------------------------------------------------------
 
        alfa_bs = 1.0
        alfa_bn = 1.5
        
        !internal angle of friction of bed material (assumed to be 35º)
        phi = 35. * pi/180.
        
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then

                AbsFlux = (Me%FluxX(i, j)**2. + Me%FluxY(i, j)**2.)**0.5

                if (AbsFlux > 0.) then
                        
                    Xaux = Me%FluxX(i, j) / AbsFlux
                    Yaux = Me%FluxY(i, j) / AbsFlux
                        
                    dhdx = 0.
                    dhdy = 0.
                        
                    if (Me%FluxX(i, j) < 0.) then                            
                        if (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then
                                
                            dhdx = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i, j-1)) /  &
                                    Me%ExternalVar%DZX(i,j-1)
                        endif
                    else
                        if (Me%ExternalVar%ComputeFacesU2D(i,j+1) == Covered) then
                                
                            dhdx = (Me%ExternalVar%Bathymetry(i, j+1) - Me%ExternalVar%Bathymetry(i, j)) /  &
                                    Me%ExternalVar%DZX(i,j)
                        endif                    
                    endif

                    if (Me%FluxY(i, j) < 0.) then
                        if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then
                        
                                dhdy = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i-1, j)) /  &
                                        Me%ExternalVar%DZY(i-1,j)     
                        endif
                    else 
                        if (Me%ExternalVar%ComputeFacesV2D(i+1, j) == Covered) then
                                
                            dhdy = (Me%ExternalVar%Bathymetry(i+1, j) - Me%ExternalVar%Bathymetry(i, j)) /  &
                                        Me%ExternalVar%DZY(i,j)        
                        endif
                    endif
                        
                    !Longitudinal bed slope
                    dzds = dhdx * Xaux + dhdy * Yaux
                        
                    dzds = min(dzds,0.9*tan(phi))                        
                        
                    alfa_s = 1 + alfa_bs * (tan(phi) / (cos(atan(dzds)) * (tan(phi) - dzds)) - 1)
                        
                    !Transverse bed slope                            
                    dzdn = -dhdx * Yaux + dhdy * Xaux
                    
                    !alfa_n = alfa_bn * (Me%CriticalShearStress(i,j) /    &
                    !        Me%ExternalVar%ShearStress(i,j))**0.5 * dzdn
                    
                    if (Me%ExternalVar%TauTotal(i,j) > 0.) then
                        alfa_n = alfa_bn * (Me%TauCritic(i,j) / Me%ExternalVar%TauTotal(i,j))**0.5 * dzdn                     
                    else
                        alfa_n = 0.
                    endif                        
                            
                    !Adjustment of bedload transport for bed-slope effects
                    
                    AuxFluxX = alfa_s * (Me%FluxX(i, j) - alfa_n * Me%FluxY(i, j))
                    AuxFluxY = alfa_s * (Me%FluxY(i, j) + alfa_n * Me%FluxX(i, j))
                    
                    Me%FluxX(i, j) = AuxFluxX
                    Me%FluxY(i, j) = AuxFluxY
                    
                endif
            endif
        enddo
        enddo


    end subroutine ComputeBedSlopeEffects
      
   !--------------------------------------------------------------------------    
   
   !--------------------------------------------------------------------------

    subroutine ComputeBiHamonicFilter2D(Prop2D, Coef)
    
        !Argument--------------------------------------------------------------
        real,   pointer, dimension(:,:)     :: Prop2D    
        real                                :: Coef
      
        !Local-----------------------------------------------------------------                
        real,   pointer, dimension(:,:)     :: Aux2D_X, Aux2D_Y
        integer                             :: i, j

        !----------------------------------------------------------------------   

        !Biharmonic filter 

        allocate (Aux2D_X(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate (Aux2D_Y(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))        

        Aux2D_X(:,:) = 0.
        Aux2D_Y(:,:) = 0.        
         
        !Biharmonic flux X         
        
        !X direction   
        do i=Me%WorkSize%ILB+1, Me%WorkSize%IUB-1
        do j=Me%WorkSize%JLB+1, Me%WorkSize%JUB-1

            if (Me%ExternalVar%OpenPoints2D(i, j-1) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i, j  ) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i, j+1) == OpenPoint) then

               Aux2D_X(i,j) = Prop2D(i, j-1)- 2.* Prop2D(i, j) + Prop2D(i, j+1)

            endif
            
        enddo
        enddo                       
        
        call BoundaryCondition(Aux2D_X)
        
        !Y direction           
        
        do i=Me%WorkSize%ILB+1, Me%WorkSize%IUB-1
        do j=Me%WorkSize%JLB+1, Me%WorkSize%JUB-1

            if (Me%ExternalVar%OpenPoints2D(i-1, j) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i,   j) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i+1, j) == OpenPoint) then

               Aux2D_Y(i,j) = Prop2D(i-1, j)- 2.* Prop2D(i, j) + Prop2D(i+1, j)

            endif
            
        enddo
        enddo                       
        
        call BoundaryCondition(Aux2D_Y)
        
        !X direction           
        do i=Me%WorkSize%ILB+1, Me%WorkSize%IUB-1
        do j=Me%WorkSize%JLB+1, Me%WorkSize%JUB-1

            if (Me%ExternalVar%OpenPoints2D(i, j-1) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i, j  ) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i, j+1) == OpenPoint) then

                Prop2D(i,j) = Prop2D(i,j) - Coef*(Aux2D_X(i, j-1)- 2. * Aux2D_X(i, j) + Aux2D_X(i, j+1))

            endif
        enddo
        enddo           

        !Y direction           
        do i=Me%WorkSize%ILB+1, Me%WorkSize%IUB-1
        do j=Me%WorkSize%JLB+1, Me%WorkSize%JUB-1

            if (Me%ExternalVar%OpenPoints2D(i-1, j) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i, j  ) == OpenPoint .and.                  &
                Me%ExternalVar%OpenPoints2D(i+1, j) == OpenPoint) then

                Prop2D(i,j) =  Prop2D(i,j) - Coef*(Aux2D_Y(i-1, j)- 2. * Aux2D_Y(i, j) + Aux2D_Y(i+1, j))

            endif
        enddo
        enddo                       
        
        call BoundaryCondition(Prop2D)        
   
        deallocate (Aux2D_X)
        deallocate (Aux2D_Y)
   

    end subroutine ComputeBiHamonicFilter2D
      
   !--------------------------------------------------------------------------      

    subroutine MeyerPeterTransport
        !Local-----------------------------------------------------------------
        real    :: Tr1,Tr2,Tr3, DeltaTau, Miu, CChezy, Clinha, Depth
        real    :: DeltaTauNoDim
        integer :: i, j

        !----------------------------------------------------------------------
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
            
                Me%TransportCapacity(i, j) = 0.  

                Depth = Me%ExternalVar%WaterColumn(I,J)

                if (Depth < Me%ExternalVar%MinWaterColumn) then
                    Depth = Me%ExternalVar%MinWaterColumn
                endif                    

                CChezy= 18.*LOG10(12.*Depth/Me%ExternalVar%CurrentRugosity(I,J))
                Clinha= 18.*LOG10(12.*Depth/Me%D90%Field2D(I,J))
                Miu   = (CChezy/Clinha)**1.5

                !N/m2    = [kg * m/s^2 * m]  
                DeltaTau = min(Miu*Me%ExternalVar%TauTotal(I,J),Me%TauMax)-Me%TauCritic(I,J)      
                
                ![ ]     = [kg * m/s^2 * m] / [kg] / [m/s^2] / [m]  
                DeltaTauNoDim = DeltaTau / (Me%RhoSl * Gravity * Me%D50%Field2D(I,J))
                  

! ---> If Tau < Tau critic there is no transport

                if (DeltaTau>0.) Then
                    ![m/s^2]^0.5                = []^0.5*[m/s^2]^0.5  
                    Tr1                         = (Me%RelativeDensity*Gravity)**0.5
                    ![m]^1.5                    = [m]^1.5                      
                    Tr2                         = Me%D50%Field2D(I,J)**1.5
                    ![ ]                        = [ ]
                    Tr3                         = DeltaTauNoDim**1.5
                    ![ ]                        
                    Me%TransportCapacity(i, j)  = 8.*Tr1*Tr2*Tr3
                
                endif
                
            endif

        enddo
        enddo


    end subroutine MeyerPeterTransport
    !--------------------------------------------------------------------------


    subroutine AckersTransport

        !Local-----------------------------------------------------------------
        !real    :: Alfa, EP, VisCin2, R32, RhoSW, Dgr, Depth
        real    :: Alfa, EP, VisCin2, R32, Dgr, Depth
        real    :: Vn, Va, Vm, Vc1, Vc, Fgr1, Fgr2, Fgr, Ggr1, Ggr, S1
        integer :: i, j
 
        !----------------------------------------------------------------------
 
! ---> Global Parameters

        Alfa   = 10.
        EP     = 1./3.
        VisCin2= WaterCinematicVisc**2.
        R32    = SQRT(32.)
!        RhoSW  = Me%Density / Me%ExternalVar%WaterDensity
        
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

                  !Me%TransportCapacity(i, j) = Ggr*RhoSW*Me%D35%Field2D(I,J)/S1               !(m)
                  Me%TransportCapacity(i, j) = Ggr*Me%D35%Field2D(I,J)/S1               !(m)
                  
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
    
    !Leo C. van Rijn1 (2007) Unified View of Sediment Transport by Currents and Waves. 
    !Two papers in     !JOURNAL OF HYDRAULIC ENGINEERING © ASCE / JUNE 2007 
    !see references below
    
    subroutine VanRijn2007Transport

        !Local-----------------------------------------------------------------
 
        !Begin-----------------------------------------------------------------
 

        Me%TransportCapacity(:, :) = 0
        
        if (Me%BedLoad) then
            call VanRijn2007BedLoad
        endif
                    
        if (Me%SuspendedLoad) then
            call VanRijn2007SuspendedLoad
        endif
        
        

    end subroutine VanRijn2007Transport

    !--------------------------------------------------------------------------

    !Leo C. van Rijn1 (2007) Unified View of Sediment Transport by Currents and Waves. 
    !I: Initiation of Motion, Bed Roughness, and Bed-Load Transport 
    !JOURNAL OF HYDRAULIC ENGINEERING © ASCE / JUNE 2007 / 649

    subroutine VanRijn2007BedLoad

        !Local-----------------------------------------------------------------
        real    :: D50, U, H, M
        integer :: i, j
 
        !----------------------------------------------------------------------
        
        ! ---> Global Parameters
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
            
                D50 = Me%D50%Field2D(I,J)
                H   = Me%ExternalVar%WaterColumn(I,J)
                U   = Me%ExternalVar%VelMod(I,J)

                !----> Compute transport capacity                
                ! [m^2/s]
                if (H > 0.01) then
                    M                          = VanRijn2007Mobility(i, j)         
                    ! [m^2/s]                  = [m^2/s]                    + [] * [m/s]*[m]*[m/m]^1.2 * []^1.5 
                    Me%TransportCapacity(i, j) = Me%TransportCapacity(i, j) + 0.015*U*H*(D50/H)**1.2*M**1.5

                endif
            endif

        enddo
        enddo

500 Format (/////T10,'A T E N C A O !!!',                        &
              //,T10,'O Valor do Diametro Adimensional e < 1',   &
              //T10,'Verifique os dados do problema ...',///)

    end subroutine VanRijn2007BedLoad

    !--------------------------------------------------------------------------

    !Leo C. van Rijn1 (2007) Unified View of Sediment Transport by Currents and Waves. 
    !II : Suspended Transport 
    !JOURNAL OF HYDRAULIC ENGINEERING © ASCE / JUNE 2007 / 668    

    subroutine VanRijn2007SuspendedLoad

        !Local-----------------------------------------------------------------
        real    :: RhoSW, D50, Di, U, M, H
        integer :: i, j
 
        !----------------------------------------------------------------------
        
        ! ---> Global Parameters

        RhoSW  = Me%Density / Me%ExternalVar%WaterDensity
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
            
                D50 = Me%D50%Field2D(I,J)
                H   = Me%ExternalVar%WaterColumn(I,J)
                U   = Me%ExternalVar%VelMod(I,J)

                !----> Compute transport capacity                
                ! [m^2/s]
                if (H > 0.01) then
                    M                          = VanRijn2007Mobility(i, j)                
                    ![]                        = [m] * [[] * [m/s^2] / [m^2/s]^2]^(1/3)
                    Di                         = D50 * ((RhoSW - 1) * Gravity / WaterCinematicVisc**2.)**(0.33333)
                    ! [m^2/s]                  = [m^2/s]                    + [] * [m/s]*[m]*[]^2.4*[]^-0.6                     
                    Me%TransportCapacity(i, j) = Me%TransportCapacity(i, j) + 0.012*U*D50*M**2.4/Di**0.6

                endif
            endif

        enddo
        enddo

500 Format (/////T10,'A T E N C A O !!!',                        &
              //,T10,'O Valor do Diametro Adimensional e < 1',   &
              //T10,'Verifique os dados do problema ...',///)

    end subroutine VanRijn2007SuspendedLoad

    !--------------------------------------------------------------------------    
    
    !Leo C. van Rijn1 (2007) Unified View of Sediment Transport by Currents and Waves. 
    !I: Initiation of Motion, Bed Roughness, and Bed-Load Transport 
    !JOURNAL OF HYDRAULIC ENGINEERING © ASCE / JUNE 2007 / 649    
   
    real function VanRijn2007Mobility(i, j)
    
        !Arguments-------------------------------------------------------------
        integer :: i, j        
    
        !Local-----------------------------------------------------------------
        real    :: RhoSW, D50, D90, U, Ucr, H
        real    :: Ubw, Alpha, Beta, Ue, Ucrc, Ucrw, PeakPeriod
        
        !Begin-----------------------------------------------------------------

        RhoSW  = Me%Density / Me%ExternalVar%WaterDensity
    
        D50 = Me%D50%Field2D            (i,j)
        D90 = Me%D90%Field2D            (i,j)
        H   = Me%ExternalVar%WaterColumn(i,j)
        U   = Me%ExternalVar%VelMod     (i,j)
        
        if (Me%WaveEffect) then
            Ubw   = Me%ExternalVar%Ubw  (i,j)
            !Alpha = 0.4 irregular waves & 0.8 regular waves
            Alpha = 0.4
            Ue    = U + Alpha * Ubw
        else
            Ue    = U
        endif            
        
! ---> Compute critical velocity for currents based on Shields initiation of motion (see Van Rijn, 1993)          
        if     (D50 > 0.05e-3 .and. D50 <= 0.5e-3) then
            Ucrc = 0.19*(D50)**0.1*log10(4*H/D90)
        elseif (D50 > 0.5e-3 .and. D50 <= 2.0e-3) then
            Ucrc = 8.50*(D50)**0.6*log10(4*H/D90)
        elseif (D50 >= 2.0e-3) then
            write(*,*) "D50 > 2 mm - not valid "
            stop "VanRijn2007CurrentWave - ModuleSand - ERR10"
        elseif (D50 <= 0.05e-3) then
            Ucrc = 0.
        endif

! ---> Compute critical velocity for waves (see Van Rijn, 1993)          
        if (Me%WaveEffect) then

            PeakPeriod = Me%ExternalVar%WavePeriod(i,j)
            
            if     (D50 > 0.05e-3 .and. D50 <= 0.5e-3) then
                Ucrw = 0.24*((RhoSW-1)*Gravity)**0.66*D50**0.33*PeakPeriod**0.33
            elseif (D50 > 0.5e-3 .and. D50 <= 2.0e-3) then
                Ucrw = 0.95*((RhoSW-1)*Gravity)**0.57*D50**0.43*PeakPeriod**0.14
            elseif (D50 >= 2.0e-3) then
                write(*,*) "D50 > 2 mm - not valid "
                stop "VanRijn2007CurrentWave - ModuleSand - ERR10"
            elseif (D50 <= 0.05e-3) then
                Ucrw = 0.
            endif

        endif
        
        if (Me%WaveEffect) then
            ![]  = [m/s] / [m/s]
            if (U == 0) then
                Beta = 0.
            else                
                Beta = U / (U + Ubw)
            endif
            ![m/s] = []*[m/s] + []*[m/s]
            Ucr  = Beta * Ucrc + (1 - Beta) * Ucrw
        else
            Ucr  = Ucrc
        endif            
        
                              

         if (Ue > Ucr) then
            ![] = [m/s] / [m/s^2*m]^0.5
            VanRijn2007Mobility = (Ue - Ucr) / sqrt((RhoSW - 1)*Gravity*D50)
        else
            VanRijn2007Mobility = 0.
        endif      
                

    end function VanRijn2007Mobility
    

    !---------------------------------------------------------------------------    
                
    
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
            
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                    Me%TauCritic(i, j) = 0.047*Gravity*Me%RhoSl*Me%D50%Field2D(I,J)      ! (N/m2)

                else

                    Me%TauCritic(i, j) = 0. 

                endif

            enddo
            enddo
            
        case default                         

!        Case (VanRijn1,VanRijn2,Bijker, VanRijn2007)
!     Case (VanRijn, Bijker)
            
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                    
                    ! []         = [m] * ([] * m/s2 / [m2/s]^2)^(1/3)
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

                    Me%TauCritic(i, j)= 0.

                endif

            enddo
            enddo
       
        END SELECT 



    end subroutine ComputeTauCritic
    !--------------------------------------------------------------------------
        
    subroutine ComputeHybridMorphEvolution


        !Local-----------------------------------------------------------------
        integer     :: STAT_CALL
        !Begin-----------------------------------------------------------------
        
        call GetGridData2Dreference(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ComputeHybridMorphEvolution - ModuleSand - ERR10' 
        
            
        call ComputeAlongShoreFlow1D

        !call ComputeProfileCrossShoreMovement
        
        call ComputeAlongShoreDZ1D        

        call ComputeProfileCrossShoreMovementDZ
        
        call ComputeNewBathymetryFromNewProfiles
        
        !Update the bathym increment
        call HybridMorphNewBathymIncrement 
        
        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ComputeHybridMorphEvolution - ModuleSand - ERR20' 

        
    end subroutine ComputeHybridMorphEvolution
    
    
    !--------------------------------------------------------------------------    

    subroutine ComputeAlongShoreFlow1D

        !Local-----------------------------------------------------------------
        real                               :: RunPeriod
        integer                            :: i, j
        logical                            :: StartProfile, EndProfile
        !----------------------------------------------------------------------
        
        Me%HybridMorph%AlongShoreFlux(:) = 0.
        
        
        !ILB coast line 
        if (Me%HybridMorph%CoastBoundary == ILB_) then
        
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
                StartProfile = .false. 
                EndProfile   = .false.                

                do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.
                        Me%HybridMorph%InShoreMapping(j) = i
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends           
                        if (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(j) = i-1
                            EndProfile = .true.                            
                            exit
                        endif                            
                    endif                    
                    
                    if (StartProfile) then                    
                        !Integrate in the "j direction" - from "ILB" line = coast line 
                        Me%HybridMorph%AlongShoreFlux(j) = Me%HybridMorph%AlongShoreFlux(j) + Me%FluxXIntegral(i, j)
                    endif
                                        
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(j) = Me%WorkSize%IUB
                endif

                if (Me%HybridMorph%InShoreMapping(j) == FillValueInt) then
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR10'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(j) == FillValueInt) then                
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR20'
                endif
                
            enddo
            
        else if (Me%HybridMorph%CoastBoundary == IUB_) then
        
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
                StartProfile = .false. 
                EndProfile   = .false.                

                do i=Me%WorkSize%IUB, Me%WorkSize%ILB, -1
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.
                        !The first water point of the profile do not move never is fixed in time 
                        Me%HybridMorph%InShoreMapping(j) = i
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends                               
                        if  (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(j) = i+1
                            EndProfile = .true.                            
                            exit
                        endif                            
                    endif        
                    
                    if (StartProfile) then
                        !Integrate in the "j direction" - from "IUB" line = coast line 
                        Me%HybridMorph%AlongShoreFlux(j) = Me%HybridMorph%AlongShoreFlux(j) + Me%FluxXIntegral(i, j)
                    endif
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(j) = Me%WorkSize%ILB
                endif              
                
                if (Me%HybridMorph%InShoreMapping(j) == FillValueInt) then
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR30'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(j) == FillValueInt) then                
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR40'
                endif                
                
            enddo          
            
        else if (Me%HybridMorph%CoastBoundary == JLB_) then
        
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                StartProfile = .false. 
                EndProfile   = .false.

                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.

                        !The first water point of the profile do not move never is fixed in time 
                        Me%HybridMorph%InShoreMapping(i) = j                       
                        
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends                               
                        if (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(i) = j-1
                            EndProfile = .true.
                            exit
                        endif                            
                    endif                        

                    if (StartProfile) then
                        !Integrate in the "i direction" - from "JLB" column = coast line 
                        Me%HybridMorph%AlongShoreFlux(i) = Me%HybridMorph%AlongShoreFlux(i) + Me%FluxYIntegral(i, j)
                    endif
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(i) = Me%WorkSize%JUB
                endif                                  
                
                if (Me%HybridMorph%InShoreMapping(i) == FillValueInt) then
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR50'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(i) == FillValueInt) then                
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR60'
                endif
                
            enddo          
            
        else if (Me%HybridMorph%CoastBoundary == JUB_) then
        
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                StartProfile = .false.
                EndProfile   = .false.

                do j = Me%WorkSize%JUB, Me%WorkSize%JLB, -1
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.
                        !The first water point of the profile do not move never is fixed in time 
                        Me%HybridMorph%InShoreMapping(i) = j
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends                               
                        if (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(i) = j+1
                            EndProfile = .true.                            
                            exit
                        endif                            
                    endif                                           
                    
                    if (StartProfile) then
                        !Integrate in the "i direction" - from "JUB" column = coast line 
                        Me%HybridMorph%AlongShoreFlux(i) = Me%HybridMorph%AlongShoreFlux(i) + Me%FluxYIntegral(i, j)
                    endif
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(i) = Me%WorkSize%JLB
                endif              
                
                if (Me%HybridMorph%InShoreMapping(i) == FillValueInt) then
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR70'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(i) == FillValueInt) then                
                    stop 'ComputeAlongShoreFlow1D - ModuleSand - ERR80'
                endif                
                
            enddo          
            
        endif
        
        call BoundaryCondition1D(Me%HybridMorph%AlongShoreFlux, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)        
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime
        

        Me%HybridMorph%ResidualAlongShoreFlux(:) = ( Me%HybridMorph%ResidualAlongShoreFlux(:) * &
                                                     (RunPeriod -  Me%Evolution%DZDT)         + &
                                                    Me%HybridMorph%AlongShoreFlux(:) *          &
                                                    Me%Evolution%DZDT) / RunPeriod

    end subroutine ComputeAlongShoreFlow1D                                                    
    
            
    !--------------------------------------------------------------------------    

    !--------------------------------------------------------------------------    

    subroutine ComputeAlongShoreDZ1D

        !Local-----------------------------------------------------------------
        real                               :: RunPeriod
        integer                            :: i, j
        logical                            :: StartProfile, EndProfile
        !----------------------------------------------------------------------
        
        Me%HybridMorph%AlongShoreDZ(:) = 0.
        
        
        !ILB coast line 
        if (Me%HybridMorph%CoastBoundary == ILB_) then
        
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
                StartProfile = .false. 
                EndProfile   = .false.                

                do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.
                        Me%HybridMorph%InShoreMapping(j) = i
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends           
                        if (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(j) = i-1
                            EndProfile = .true.                            
                            exit
                        endif                            
                    endif                    
                    
                    if (StartProfile) then                                        
                        !Integrate in the "i direction" - from "ILB" line = coast line to closure depth   
                        Me%HybridMorph%AlongShoreDZ(j) = Me%HybridMorph%AlongShoreDZ(j) +   &
                            Me%DZ%Field2D(i, j) * Me%ExternalVar%DVY(i, j) / Me%HybridMorph%ClosureDepth
                    endif
                                        
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(j) = Me%WorkSize%IUB
                endif

                if (Me%HybridMorph%InShoreMapping(j) == FillValueInt) then
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR10'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(j) == FillValueInt) then                
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR20'
                endif
                
            enddo
            
        else if (Me%HybridMorph%CoastBoundary == IUB_) then
        
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
                StartProfile = .false. 
                EndProfile   = .false.                

                do i=Me%WorkSize%IUB, Me%WorkSize%ILB, -1
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.
                        !The first water point of the profile do not move never is fixed in time 
                        Me%HybridMorph%InShoreMapping(j) = i
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends                               
                        if  (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(j) = i+1
                            EndProfile = .true.                            
                            exit
                        endif                            
                    endif        
                    
                    if (StartProfile) then
                        !Integrate in the "i direction" - from "IUB" line = coast line to closure depth   
                        Me%HybridMorph%AlongShoreDZ(j) = Me%HybridMorph%AlongShoreDZ(j) + &
                            Me%DZ%Field2D(i, j) * Me%ExternalVar%DVY(i, j) / Me%HybridMorph%ClosureDepth
                    endif
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(j) = Me%WorkSize%ILB
                endif              
                
                if (Me%HybridMorph%InShoreMapping(j) == FillValueInt) then
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR30'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(j) == FillValueInt) then                
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR40'
                endif                
                
            enddo          
            
        else if (Me%HybridMorph%CoastBoundary == JLB_) then
        
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                StartProfile = .false. 
                EndProfile   = .false.

                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.

                        !The first water point of the profile do not move never is fixed in time 
                        Me%HybridMorph%InShoreMapping(i) = j                       
                        
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends                               
                        if (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(i) = j-1
                            EndProfile = .true.
                            exit
                        endif                            
                    endif                        

                    if (StartProfile) then
                        !Integrate in the "j direction" - from "JLB" line = coast line to closure depth   
                        Me%HybridMorph%AlongShoreDZ(i) = Me%HybridMorph%AlongShoreDZ(i) + &
                            Me%DZ%Field2D(i, j) * Me%ExternalVar%DUX(i, j) / Me%HybridMorph%ClosureDepth
                    endif
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(i) = Me%WorkSize%JUB
                endif                                  
                
                if (Me%HybridMorph%InShoreMapping(i) == FillValueInt) then
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR50'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(i) == FillValueInt) then                
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR60'
                endif
                
            enddo          
            
        else if (Me%HybridMorph%CoastBoundary == JUB_) then
        
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                StartProfile = .false.
                EndProfile   = .false.

                do j = Me%WorkSize%JUB, Me%WorkSize%JLB, -1
            
                    !Defines the first cell of the beach profile
                    if (Me%ExternalVar%WaterPoints2D(i, j)  == WaterPoint .and. .not. StartProfile) then
                        StartProfile = .true.
                        !The first water point of the profile do not move never is fixed in time 
                        Me%HybridMorph%InShoreMapping(i) = j
                    endif      
                    
                    !Found a land cell (e.g. detached breakwater) - the sand profile active area ends so the alongshore transport integration also ends
                    if (StartProfile) then
                        !If the closure depth is reached  - the sand profile active area ends so the alongshore transport integration also ends                               
                        if (Me%ExternalVar%InitialBathym(i,j) > Me%HybridMorph%ClosureDepth) then
                            Me%HybridMorph%OffShoreMapping(i) = j+1
                            EndProfile = .true.                            
                            exit
                        endif                            
                    endif                                           
                    
                    if (StartProfile) then
                        !Integrate in the "j direction" - from "JUB" line = coast line to closure depth   
                        Me%HybridMorph%AlongShoreDZ(i) = Me%HybridMorph%AlongShoreDZ(i) + &
                            Me%DZ%Field2D(i, j) * Me%ExternalVar%DUX(i, j) / Me%HybridMorph%ClosureDepth
                    endif
                enddo
                
                if (.not. EndProfile) then
                    Me%HybridMorph%OffShoreMapping(i) = Me%WorkSize%JLB
                endif              
                
                if (Me%HybridMorph%InShoreMapping(i) == FillValueInt) then
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR70'
                endif
                
                if (Me%HybridMorph%OffShoreMapping(i) == FillValueInt) then                
                    stop 'ComputeAlongShoreDZ1D - ModuleSand - ERR80'
                endif                
                
            enddo          
            
        endif
        
        call BoundaryCondition1D(Me%HybridMorph%AlongShoreDZ, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)        
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime
        

        Me%HybridMorph%ResidualAlongShoreDZ(:) = ( Me%HybridMorph%ResidualAlongShoreDZ(:) * &
                                                   (RunPeriod -  Me%Evolution%DZDT)         + &
                                                    Me%HybridMorph%AlongShoreDZ(:)        * &
                                                    Me%Evolution%DZDT) / RunPeriod

    end subroutine ComputeAlongShoreDZ1D                                                    
    
            
    !--------------------------------------------------------------------------    


    subroutine ComputeProfileCrossShoreMovement

        !Local-----------------------------------------------------------------
        real(8),    pointer, dimension(:)  :: Aux1D
        real                               :: DX1, DX2, DY1, DY2, Area1, Area2, RunPeriod, K
        integer                            :: i, j, ij, imin, imax, di, iaux
        
        !----------------------------------------------------------------------
        
        Me%HybridMorph%CrossShoreVel(:) = 0.

        !ILB or IUB coast line 
        if (Me%HybridMorph%CoastBoundary == ILB_ .or. Me%HybridMorph%CoastBoundary == IUB_) then
        
            if (Me%HybridMorph%CoastBoundary == ILB_) then
                ij =  Me%WorkSize%ILB                    
            endif                
            if (Me%HybridMorph%CoastBoundary == IUB_) then        
                ij =  Me%WorkSize%IUB
            endif
            
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB

                DX1   = Me%ExternalVar%DUX(ij, j-1)
                Area1 = Me%HybridMorph%ClosureDepth * DX1

                DX2   = Me%ExternalVar%DUX(ij, j  )
                Area2 = Me%HybridMorph%ClosureDepth * DX2
                
                ! [m/s]                         = [m/s] + [m^3/s] / [m]^2
                Me%HybridMorph%CrossShoreVel(j-1) = Me%HybridMorph%CrossShoreVel(j-1) - Me%HybridMorph%AlongShoreFlux(j) / Area1
                Me%HybridMorph%CrossShoreVel(j  ) = Me%HybridMorph%CrossShoreVel(j  ) + Me%HybridMorph%AlongShoreFlux(j) / Area2
            
            enddo          
            
        else if (Me%HybridMorph%CoastBoundary == JLB_ .or. Me%HybridMorph%CoastBoundary == JUB_) then
        
            if (Me%HybridMorph%CoastBoundary == JLB_) then
                ij = Me%WorkSize%JLB                    
            endif                
            if (Me%HybridMorph%CoastBoundary == JUB_) then        
                ij = Me%WorkSize%JUB
            endif
            
                
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                DY1   = Me%ExternalVar%DVY(i  , ij)
                Area1 = Me%HybridMorph%ClosureDepth * DY1

                DY2   = Me%ExternalVar%DVY(i+1, ij)
                Area2 = Me%HybridMorph%ClosureDepth * DY2
                
                ! [m/s]                         = [m/s] + [m^3/s] / [m]^2
                Me%HybridMorph%CrossShoreVel(i-1) = Me%HybridMorph%CrossShoreVel(i-1) - Me%HybridMorph%AlongShoreFlux(i) / Area1
                Me%HybridMorph%CrossShoreVel(i  ) = Me%HybridMorph%CrossShoreVel(i  ) + Me%HybridMorph%AlongShoreFlux(i) / Area2
                        
            enddo          
            
        endif
                
        call BoundaryCondition1D(Me%HybridMorph%CrossShoreVel, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)
        
        if (Me%HybridMorph%DintegLongShore > 1) then

            allocate (Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D))
            
            Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D) = &
                Me%HybridMorph%CrossShoreVel(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D)
            
            do j = Me%HybridMorph%Min1D+1, Me%HybridMorph%Max1D-1
            
                imin=max(Me%HybridMorph%Min1D,j- Me%HybridMorph%DintegLongShore)
                imax=min(Me%HybridMorph%Max1D,j+ Me%HybridMorph%DintegLongShore)            
                
                Me%HybridMorph%CrossShoreVel(j) = 0.
                
                di = imax-imin
                do i = imin, imax
                    Me%HybridMorph%CrossShoreVel(j) = Me%HybridMorph%CrossShoreVel(j) + Aux1D(i)/real(di)
                enddo            
                
            enddo            
            
            deallocate (Aux1D)
            
            call BoundaryCondition1D(Me%HybridMorph%CrossShoreVel, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)
                        
        endif
        
         
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime
        

        Me%HybridMorph%ResidualCrossShoreVel(:) = ( Me%HybridMorph%ResidualCrossShoreVel(:) * &
                                                     (RunPeriod -  Me%Evolution%DZDT)         + &
                                                    Me%HybridMorph%CrossShoreVel(:) *          &
                                                    Me%Evolution%DZDT) / RunPeriod        

!        !Biharmonic filter 

        allocate (Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D))
        
        Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D) = 0

        K = 1./8.

        do j = Me%HybridMorph%Min1D+1, Me%HybridMorph%Max1D-1
           Aux1D(j) = (Me%HybridMorph%ResidualCrossShoreVel(j-1)- 2.* Me%HybridMorph%ResidualCrossShoreVel(j) + &
                       Me%HybridMorph%ResidualCrossShoreVel(j+1))
        enddo            
        
        call BoundaryCondition1D(Aux1D, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)    
        
         do j = Me%HybridMorph%Min1D+1, Me%HybridMorph%Max1D-1
            Me%HybridMorph%ResidualCrossShoreVel(j) =  Me%HybridMorph%ResidualCrossShoreVel(j) - K*(Aux1D(j-1)- &
                                                       2. * Aux1D(j) + Aux1D(j+1))
        enddo                        
        
        deallocate (Aux1D)
        
        call BoundaryCondition1D(Me%HybridMorph%ResidualCrossShoreVel, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)

        !Buffer area adjacent to open boundary with no transport
        if (Me%HybridMorph%NoTransportBufferCells > 0) then
            iaux = Me%HybridMorph%NoTransportBufferCells
            Me%HybridMorph%ResidualCrossShoreVel(Me%HybridMorph%Min1D     : Me%HybridMorph%Min1D+iaux) = 0.
            Me%HybridMorph%ResidualCrossShoreVel(Me%HybridMorph%Max1D-iaux: Me%HybridMorph%Max1D    ) = 0.        
        endif

    end subroutine ComputeProfileCrossShoreMovement
    
    !--------------------------------------------------------------------------

    subroutine ComputeProfileCrossShoreMovementDZ

        !Local-----------------------------------------------------------------
        real(8),    pointer, dimension(:)  :: Aux1D
        real                               :: RunPeriod, K
        integer                            :: i, j, ij, imin, imax, di, iaux
        
        !----------------------------------------------------------------------
        
        Me%HybridMorph%CrossShoreVel(:) = 0.

        !ILB or IUB coast line 
        if (Me%HybridMorph%CoastBoundary == ILB_ .or. Me%HybridMorph%CoastBoundary == IUB_) then
        
            if (Me%HybridMorph%CoastBoundary == ILB_) then
                ij =  Me%WorkSize%ILB                    
            endif                
            if (Me%HybridMorph%CoastBoundary == IUB_) then        
                ij =  Me%WorkSize%IUB
            endif
            
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB

                ! [m/s]                         = [m/s] + [m] / [s]
                Me%HybridMorph%CrossShoreVel(j) = Me%HybridMorph%CrossShoreVel(j) + &
                    Me%HybridMorph%AlongShoreDZ(j) / Me%Evolution%DZDT
            
            enddo          
            
        else if (Me%HybridMorph%CoastBoundary == JLB_ .or. Me%HybridMorph%CoastBoundary == JUB_) then
        
            if (Me%HybridMorph%CoastBoundary == JLB_) then
                ij = Me%WorkSize%JLB                    
            endif                
            if (Me%HybridMorph%CoastBoundary == JUB_) then        
                ij = Me%WorkSize%JUB
            endif
            
                
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
                ! [m/s]                         = [m/s] + [m] / [s]
                Me%HybridMorph%CrossShoreVel(i) = Me%HybridMorph%CrossShoreVel(i) + &
                    Me%HybridMorph%AlongShoreDZ(i) / Me%Evolution%DZDT
                        
            enddo          
            
        endif
                
        call BoundaryCondition1D(Me%HybridMorph%CrossShoreVel, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)
        
        if (Me%HybridMorph%DintegLongShore > 1) then

            allocate (Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D))
            
            Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D) = &
                Me%HybridMorph%CrossShoreVel(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D)
            
            do j = Me%HybridMorph%Min1D+1, Me%HybridMorph%Max1D-1
            
                imin=max(Me%HybridMorph%Min1D,j- Me%HybridMorph%DintegLongShore)
                imax=min(Me%HybridMorph%Max1D,j+ Me%HybridMorph%DintegLongShore)            
                
                Me%HybridMorph%CrossShoreVel(j) = 0.
                
                di = imax-imin
                do i = imin, imax
                    Me%HybridMorph%CrossShoreVel(j) = Me%HybridMorph%CrossShoreVel(j) + Aux1D(i)/real(di)
                enddo            
                
            enddo            
            
            deallocate (Aux1D)
            
            call BoundaryCondition1D(Me%HybridMorph%CrossShoreVel, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)
                        
        endif
        
         
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime
        

        Me%HybridMorph%ResidualCrossShoreVel(:) = ( Me%HybridMorph%ResidualCrossShoreVel(:) * &
                                                     (RunPeriod -  Me%Evolution%DZDT)         + &
                                                    Me%HybridMorph%CrossShoreVel(:) *          &
                                                    Me%Evolution%DZDT) / RunPeriod        

!        !Biharmonic filter 

        allocate (Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D))
        
        Aux1D(Me%HybridMorph%Min1D:Me%HybridMorph%Max1D) = 0

        K = 1./8.

        do j = Me%HybridMorph%Min1D+1, Me%HybridMorph%Max1D-1
           Aux1D(j) = (Me%HybridMorph%ResidualCrossShoreVel(j-1)- 2.* Me%HybridMorph%ResidualCrossShoreVel(j) + &
                       Me%HybridMorph%ResidualCrossShoreVel(j+1))
        enddo            
        
        call BoundaryCondition1D(Aux1D, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)    
        
         do j = Me%HybridMorph%Min1D+1, Me%HybridMorph%Max1D-1
            Me%HybridMorph%ResidualCrossShoreVel(j) =  Me%HybridMorph%ResidualCrossShoreVel(j) - K*(Aux1D(j-1)- &
                                                       2. * Aux1D(j) + Aux1D(j+1))
        enddo                        
        
        deallocate (Aux1D)
        
        call BoundaryCondition1D(Me%HybridMorph%ResidualCrossShoreVel, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D)

        !Buffer area adjacent to open boundary with no transport
        if (Me%HybridMorph%NoTransportBufferCells > 0) then
            iaux = Me%HybridMorph%NoTransportBufferCells
            Me%HybridMorph%ResidualCrossShoreVel(Me%HybridMorph%Min1D     : Me%HybridMorph%Min1D+iaux) = 0.
            Me%HybridMorph%ResidualCrossShoreVel(Me%HybridMorph%Max1D-iaux: Me%HybridMorph%Max1D    ) = 0.        
        endif

    end subroutine ComputeProfileCrossShoreMovementDZ
    
    !--------------------------------------------------------------------------


    subroutine ComputeNewBathymetryFromNewProfiles

        !Local-----------------------------------------------------------------
        real(8), dimension(:)  , pointer   :: XXref, XXinst, InitProfile, NextProfile
        integer, dimension(:)  , pointer   :: WaterPoints1D
        real                               :: RunPeriod
        integer                            :: i, j
        integer                            :: ILB, IUB, JLB, JUB
        !----------------------------------------------------------------------
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime     
        
        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB
            
        !ILB coast Me%WorkSize%ILB:Me%WorkSize%IUB 
        if (Me%HybridMorph%CoastBoundary == ILB_) then

            allocate(XXref        (Me%Size%ILB:Me%Size%IUB))
            allocate(XXinst       (Me%Size%ILB:Me%Size%IUB))
            allocate(InitProfile  (Me%Size%ILB:Me%Size%IUB))
            allocate(NextProfile  (Me%Size%ILB:Me%Size%IUB))
            allocate(WaterPoints1D(Me%Size%ILB:Me%Size%IUB))            
            
            WaterPoints1D(:) = 0                
            
            do j= JLB, JUB
                !profile movement
                do i= Me%HybridMorph%InShoreMapping(j), Me%HybridMorph%OffShoreMapping(j)
                
                    Me%HybridMorph%DistanceToCoastInst(i, j) = Me%HybridMorph%DistanceToCoastRef(i, j) + &
                        Me%HybridMorph%ResidualCrossShoreVel(j) * RunPeriod
                
                enddo            
                
                do i = ILB, IUB
                    XXref        (i)  = Me%HybridMorph%DistanceToCoastRef (i, j)
                    XXinst       (i)  = Me%HybridMorph%DistanceToCoastInst(i, j)
                    InitProfile  (i)  = Me%ExternalVar%InitialBathym      (i, j)
                    NextProfile  (i)  = Me%HybridMorph%BathymetryNext     (i, j)
                    WaterPoints1D(i)  = Me%ExternalVar%WaterPoints2D      (i, j)
                enddo                    


                call InterpolateNewBathymProfile(LB             = ILB          ,        &
                                                 UB             = IUB          ,        &
                                                 XXref          = XXref        ,        & 
                                                 XXinst         = XXinst       ,        & 
                                                 InitProfile    = InitProfile  ,        & 
                                                 NextProfile    = NextProfile  ,        & 
                                                 WaterPoints1D  = WaterPoints1D)

                do i = ILB, IUB                                                                                 
                    Me%HybridMorph%BathymetryNext(i,j) = NextProfile(i)
                enddo                        

            enddo  
            
        
        !IUB coast line                         
        else if (Me%HybridMorph%CoastBoundary == IUB_) then

            allocate(XXref        (Me%Size%ILB:Me%Size%IUB))
            allocate(XXinst       (Me%Size%ILB:Me%Size%IUB))
            allocate(InitProfile  (Me%Size%ILB:Me%Size%IUB))
            allocate(NextProfile  (Me%Size%ILB:Me%Size%IUB))
            allocate(WaterPoints1D(Me%Size%ILB:Me%Size%IUB))     
            
            WaterPoints1D(:) = 0

            do j= JLB, JUB
                !profile movement
                do i= Me%HybridMorph%InShoreMapping(j), Me%HybridMorph%OffShoreMapping(j),-1
                
                    Me%HybridMorph%DistanceToCoastInst(i, j) = Me%HybridMorph%DistanceToCoastRef(i, j) + &
                        Me%HybridMorph%ResidualCrossShoreVel(j) * RunPeriod
                
                enddo            
                
                do i = ILB, IUB
                    XXref        (IUB + ILB -i)  = Me%HybridMorph%DistanceToCoastRef (i, j)
                    XXinst       (IUB + ILB -i)  = Me%HybridMorph%DistanceToCoastInst(i, j)
                    InitProfile  (IUB + ILB -i)  = Me%ExternalVar%InitialBathym      (i, j)
                    NextProfile  (IUB + ILB -i)  = Me%HybridMorph%BathymetryNext     (i, j)
                    WaterPoints1D(IUB + ILB -i)  = Me%ExternalVar%WaterPoints2D      (i, j)
                enddo                    


                call InterpolateNewBathymProfile(LB             = ILB          ,        &
                                                 UB             = IUB          ,        &
                                                 XXref          = XXref        ,        & 
                                                 XXinst         = XXinst       ,        & 
                                                 InitProfile    = InitProfile  ,        & 
                                                 NextProfile    = NextProfile  ,        & 
                                                 WaterPoints1D  = WaterPoints1D)

                do i= ILB, IUB
                    Me%HybridMorph%BathymetryNext(i,j) = NextProfile(IUB+ILB-i)
                enddo    
            enddo  
            
        !JLB coast line                         
        else if (Me%HybridMorph%CoastBoundary == JLB_) then

            allocate(XXref        (Me%Size%JLB:Me%Size%JUB))
            allocate(XXinst       (Me%Size%JLB:Me%Size%JUB))
            allocate(InitProfile  (Me%Size%JLB:Me%Size%JUB))
            allocate(NextProfile  (Me%Size%JLB:Me%Size%JUB))
            allocate(WaterPoints1D(Me%Size%JLB:Me%Size%JUB)) 
            
            WaterPoints1D(:) = 0            
        
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                !profile movement
                do j= Me%HybridMorph%InShoreMapping(i), Me%HybridMorph%OffShoreMapping(i)
                
                    Me%HybridMorph%DistanceToCoastInst(i, j) = Me%HybridMorph%DistanceToCoastRef(i, j) + &
                        Me%HybridMorph%ResidualCrossShoreVel(i) * RunPeriod
                
                enddo           
                
                do j = JLB, JUB
                    XXref        (j)  = Me%HybridMorph%DistanceToCoastRef (i, j)
                    XXinst       (j)  = Me%HybridMorph%DistanceToCoastInst(i, j)
                    InitProfile  (j)  = Me%ExternalVar%InitialBathym      (i, j)
                    NextProfile  (j)  = Me%HybridMorph%BathymetryNext     (i, j)
                    WaterPoints1D(j)  = Me%ExternalVar%WaterPoints2D      (i, j)
                enddo                        
                
                call InterpolateNewBathymProfile(LB             = JLB,          &
                                                 UB             = JUB,          & 
                                                 XXref          = XXref        ,& 
                                                 XXinst         = XXinst       ,& 
                                                 InitProfile    = InitProfile  ,& 
                                                 NextProfile    = NextProfile  ,& 
                                                 WaterPoints1D  = WaterPoints1D)

                do j = JLB, JUB
                    Me%HybridMorph%BathymetryNext(i,j) = NextProfile(j)
                enddo 
                           
            enddo  

        
        !JUB coast line                         
        else if (Me%HybridMorph%CoastBoundary == JUB_) then
        
            allocate(XXref        (Me%Size%JLB:Me%Size%JUB))
            allocate(XXinst       (Me%Size%JLB:Me%Size%JUB))
            allocate(InitProfile  (Me%Size%JLB:Me%Size%JUB))
            allocate(NextProfile  (Me%Size%JLB:Me%Size%JUB))
            allocate(WaterPoints1D(Me%Size%JLB:Me%Size%JUB)) 
            
            WaterPoints1D(:) = 0            
            
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                !profile movement
                do j= Me%HybridMorph%InShoreMapping(i), Me%HybridMorph%OffShoreMapping(i),-1
                
                    Me%HybridMorph%DistanceToCoastInst(i, j) = Me%HybridMorph%DistanceToCoastRef(i, j) + &
                        Me%HybridMorph%ResidualCrossShoreVel(i) * RunPeriod
                
                enddo     
                       
                do j = JLB, JUB
                    XXref        (JUB + JLB -j)  = Me%HybridMorph%DistanceToCoastRef (i, j)
                    XXinst       (JUB + JLB -j)  = Me%HybridMorph%DistanceToCoastInst(i, j)
                    InitProfile  (JUB + JLB -j)  = Me%ExternalVar%InitialBathym      (i, j)
                    NextProfile  (JUB + JLB -j)  = Me%HybridMorph%BathymetryNext     (i, j)
                    WaterPoints1D(JUB + JLB -j)  = Me%ExternalVar%WaterPoints2D      (i, j)
                enddo                    

                call InterpolateNewBathymProfile(LB             = JLB,          &
                                                 UB             = JUB,          &
                                                 XXref          = XXref        ,& 
                                                 XXinst         = XXinst       ,& 
                                                 InitProfile    = InitProfile  ,& 
                                                 NextProfile    = NextProfile  ,& 
                                                 WaterPoints1D  = WaterPoints1D)

                do j= JLB, JUB
                    Me%HybridMorph%BathymetryNext(i,j) = NextProfile(JUB+JLB-j)
                enddo    
            enddo  
            
        endif
        
        deallocate(XXref, XXinst, InitProfile, NextProfile, WaterPoints1D)
        
    end subroutine ComputeNewBathymetryFromNewProfiles
    
    !--------------------------------------------------------------------------
    
    subroutine InterpolateNewBathymProfile(LB, UB, XXref, XXinst, InitProfile, NextProfile, WaterPoints1D)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:)  , pointer   :: XXref, XXinst, InitProfile, NextProfile
        integer, dimension(:)  , pointer   :: WaterPoints1D
        integer                            :: LB, UB

        !Local-----------------------------------------------------------------
        real(8)                            :: dx1, dx2
        integer                            :: ij, iSouth, iW1, iW2        
        !----------------------------------------------------------------------
        
        NextProfile(:) = InitProfile(:)
        
        do ij = LB, UB
    
            call LocateCell1D (XXinst, XXref(ij), LB, UB, iSouth)
                               
                                                      
            !linear interpolation
            iW1 = WaterPoints1D(iSouth)
            iW2 = WaterPoints1D(iSouth+1)
                                    
            if (iW1+iW2 == 2) then
            
                dx1 = XXref (ij      ) - XXinst(iSouth)
                dx2 = XXinst(iSouth+1) - XXref (ij    )  
                
                if (dx1 < 0.) then
                    stop 'InterpolateNewBathymProfile - ModuleSand - ERR10'
                endif

                if (dx2 < 0.) then
                    stop 'InterpolateNewBathymProfile - ModuleSand - ERR20'
                endif
            
                NextProfile(ij)= (dx2 *  InitProfile(iSouth) + dx1 *  InitProfile(iSouth+1)) /(dx1 + dx2)
                                                      
            endif    
    
        enddo
        
    end subroutine InterpolateNewBathymProfile
    
    !--------------------------------------------------------------------------

    subroutine HybridMorphNewBathymIncrement

        !Local-----------------------------------------------------------------
        integer                            :: i, j
        !----------------------------------------------------------------------
        

    
        if (Me%ExternalVar%Now >= Me%Evolution%NextBatim) then

            do j= Me%WorkSize%JLB, Me%WorkSize%JUB
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                    Me%HybridMorph%DZ_Residual%Field2D(i, j) = Me%ExternalVar%InitialBathym     (i, j) - &
                                                               Me%HybridMorph%BathymetryNext(i, j)
                                                               
!                    Me%BatimIncrement%Field2D (i, j)         = Me%HybridMorph%BathymetryPrevious(i, j) - &
!                                                               Me%HybridMorph%BathymetryNext(i, j)
                endif
            enddo
            enddo                                                

            !Update the previous bathymetry
            Me%HybridMorph%BathymetryPrevious(:,:) = Me%HybridMorph%BathymetryNext(:,:)
            
        endif    

    end subroutine HybridMorphNewBathymIncrement

    !--------------------------------------------------------------------------

    subroutine ComputeEvolutionX            

        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer   :: FluxXY
        real(8), dimension(:,:), pointer   :: BoxFluxesXY
        integer, dimension(:,:), pointer   :: ComputeUV2D         
        real(8)                            :: Area2, RunPeriod, FluxSand, SandMin
        real(8)                            :: VolumeAvailable, SandThickness
        integer                            :: i, j, di, dj
        !----------------------------------------------------------------------
        
        Me%DZ%Field2D(:,:) = 0.
        
        if (Me%Evolution%Direction == 1) then
            Me%Evolution%Direction = 2
        else
            Me%Evolution%Direction = 1
        endif

!        if (Me%Boxes%Yes) then
!            Me%Boxes%FluxesX(:,:) = 0.
!            Me%Boxes%FluxesY(:,:) = 0.
!        endif

        if (Me%Aceleration%Yes) then

            Me%FluxXIntegral (:,:) = Me%Aceleration%Coef * Me%FluxXIntegral (:,:)
            Me%FluxYIntegral (:,:) = Me%Aceleration%Coef * Me%FluxYIntegral (:,:)
            
            if (Me%MudErosion%ON) then
                Me%FluxZIntegral (:,:) = Me%Aceleration%Coef * Me%FluxZIntegral (:,:)
            endif    

        endif
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime

        Me%Residual%FluxX(:,:) = ( Me%Residual%FluxX(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxXIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod

        Me%Residual%FluxY(:,:) = ( Me%Residual%FluxY(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxYIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod
        
        if (Me%MudErosion%ON) then
            Me%Residual%FluxZ(:,:) = ( Me%Residual%FluxZ(:,:) * (RunPeriod -  Me%Evolution%DZDT) + &
                                       Me%FluxZIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod
            
        endif    
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
          
            Area2    = Me%ExternalVar%DUX(i, j  ) * Me%ExternalVar%DVY(i, j)
            
            FluxSand = sqrt(Me%FluxXIntegral(i,j)**2. + Me%FluxYIntegral(i,j)**2.)
            
            ! FluxSand [m3/s] < Water velocity [m/s]* Sand Thickness [m] * Cell width [m]
            ! Sand grain limited by the water flow velocity
            ![m]    = [m3/s]   / [m/s] / [m]
            if (Me%ExternalVar%VelMod(i, j) > 0.) then
                SandMin = FluxSand / Me%ExternalVar%VelMod(i, j) / sqrt(Area2) 
            else
                SandMin = Me%SandMin 
            endif
            
            SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
            
            if (SandThickness < SandMin) then
                SandThickness = 0.
            endif
            
            VolumeAvailable = SandThickness * Area2            
      
            
            if (Me%Evolution%Direction == 1) then
                dj = 1
                di = 0
                FluxXY      => Me%FluxXIntegral
                BoxFluxesXY => Me%Boxes%FluxesX
                ComputeUV2D => Me%ExternalVar%ComputeFacesU2D
            else
                dj = 0
                di = 1
                FluxXY      => Me%FluxYIntegral
                BoxFluxesXY => Me%Boxes%FluxesY
                ComputeUV2D => Me%ExternalVar%ComputeFacesV2D
            endif            
            
            call ComputeOneDir(FluxXY, BoxFluxesXY, ComputeUV2D, Area2, VolumeAvailable, i, j, di, dj)

            if (Me%Evolution%Direction == 2) then
                dj = 1
                di = 0
                FluxXY      => Me%FluxXIntegral
                BoxFluxesXY => Me%Boxes%FluxesX
                ComputeUV2D => Me%ExternalVar%ComputeFacesU2D
            else
                dj = 0
                di = 1
                FluxXY      => Me%FluxYIntegral
                BoxFluxesXY => Me%Boxes%FluxesY
                ComputeUV2D => Me%ExternalVar%ComputeFacesV2D
            endif               
            
            call ComputeOneDir(FluxXY, BoxFluxesXY, ComputeUV2D, Area2, VolumeAvailable, i, j, di, dj)

        enddo
        enddo
        
        if (Me%MudErosion%ON) then
            call ComputeMudErosionFlux
        endif    
        
        call ComputeDischarges
        
        nullify(FluxXY, BoxFluxesXY, ComputeUV2D)
        
        
    end subroutine ComputeEvolutionX
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    
    subroutine ComputeMudErosionFlux

        !Arguments-------------------------------------------------------------
       

        !Local-----------------------------------------------------------------
        real(8)                            :: Area2, VolumeAvailable
        real(8)                            :: VolumeTransport
        real(8)                            :: Coef
        integer                            :: i, j         
        
        !Begin-----------------------------------------------------------------
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB  
                        
            if (Me%ExternalVar%OpenPoints2D(i, j) == 1) then
            
                VolumeTransport = Me%FluxZIntegral(i, j) * Me%Evolution%DZDT
                    
                if (VolumeTransport > VolumeAvailable) then
                    
                    if (VolumeTransport > 0) then
                        Coef = VolumeAvailable / VolumeTransport
                    else
                        Coef = 0
                    endif
            
                    Me%FluxZIntegral(i, j) = Coef * Me%FluxZIntegral(i, j) 
                        
                    VolumeTransport = VolumeAvailable

                    Area2 = Me%ExternalVar%DUX(i, j) * Me%ExternalVar%DVY(i, j)
                    
                    Me%DZ%Field2D(i, j) = Me%DZ%Field2D(i, j) - VolumeTransport / Area2
            
                    VolumeAvailable = VolumeAvailable - VolumeTransport

                    if (Me%Boxes%Yes) then
                        Me%Boxes%FluxesZ(i,j) = Me%Boxes%FluxesZ(i,j) + Me%FluxZIntegral(i, j)* Me%Evolution%DZDT
                    endif
                    
                endif
            endif

        enddo
        enddo
        
    end subroutine ComputeMudErosionFlux
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine ComputeOneDir(FluxXY, BoxFluxesXY, ComputeUV2D, Area2, VolumeAvailable, i, j, di, dj)            

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:), pointer   :: FluxXY
        real(8), dimension(:,:), pointer   :: BoxFluxesXY 
        integer, dimension(:,:), pointer   :: ComputeUV2D         
        real(8)                            :: Area2, VolumeAvailable
        integer                            :: i, j, di, dj
        !Local-----------------------------------------------------------------
        real(8)                            :: Area1
        real(8)                            :: VolumeTransport
        real(8)                            :: Coef, DZ
        
        !Begin-----------------------------------------------------------------
            
        VolumeTransport = abs(FluxXY(i, j)) * Me%Evolution%DZDT
                    
        if (VolumeTransport > VolumeAvailable) then
            if (VolumeTransport > 0) then
                Coef = VolumeAvailable / VolumeTransport
            else
                Coef = 0
            endif
            FluxXY(i, j) = Coef * FluxXY(i, j)
                        
            VolumeTransport = VolumeAvailable

        endif                       
            
            
        if      (ComputeUV2D(i,  j) == Covered  .and. FluxXY(i, j) < 0.) then
                    
            Area1    = Me%ExternalVar%DUX(i-di, j-dj) * Me%ExternalVar%DVY(i-di, j-dj)

            DZ                    = Me%DZ%Field2D(i-di, j-dj) + VolumeTransport / Area1 
            Me%DZ%Field2D(i-di, j-dj) = DZ
            DZ                    = Me%DZ%Field2D(i, j  ) - VolumeTransport / Area2
            Me%DZ%Field2D(i, j  ) = DZ
                    
            VolumeAvailable = VolumeAvailable - VolumeTransport

            if (Me%Boxes%Yes) then
                BoxFluxesXY(i,j) = BoxFluxesXY(i,j) + FluxXY(i, j) * Me%Evolution%DZDT
            endif
            
        elseif (ComputeUV2D(i+di,j+dj) == Covered  .and. FluxXY(i, j) > 0.) then

            Area1 = Me%ExternalVar%DUX(i+di, j+dj) * Me%ExternalVar%DVY(i+di, j+dj)

            DZ                    = Me%DZ%Field2D(i+di, j+dj) + VolumeTransport / Area1   
            Me%DZ%Field2D(i+di, j+dj) = DZ
            DZ                    = Me%DZ%Field2D(i, j  ) - VolumeTransport / Area2 
            Me%DZ%Field2D(i, j  ) = DZ
                    
            VolumeAvailable = VolumeAvailable - VolumeTransport                    

            if (Me%Boxes%Yes) then
                BoxFluxesXY(i+di,j+dj) = BoxFluxesXY(i+di,j+dj) + FluxXY(i, j) * Me%Evolution%DZDT
            endif
                   
        endif

        
    end subroutine ComputeOneDir
    
    !--------------------------------------------------------------------------
    

    subroutine ComputeEvolution            

        !Local-----------------------------------------------------------------
        real(8)                            :: Area1, Area2, RunPeriod, SandThickness
        real(8)                            :: VolumeTransport, VolumeAvailable, Coef, DZ
        integer                            :: i, j
        !----------------------------------------------------------------------
        
        Me%DZ%Field2D(:,:) = 0.

!        if (Me%Boxes%Yes) then
!            Me%Boxes%FluxesX(:,:) = 0.
!            Me%Boxes%FluxesY(:,:) = 0.
!        endif

        if (Me%Aceleration%Yes) then

            Me%FluxXIntegral (:,:) = Me%Aceleration%Coef * Me%FluxXIntegral (:,:)
            Me%FluxYIntegral (:,:) = Me%Aceleration%Coef * Me%FluxYIntegral (:,:)

        endif
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime

        Me%Residual%FluxX(:,:) = ( Me%Residual%FluxX(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxXIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod

        Me%Residual%FluxY(:,:) = ( Me%Residual%FluxY(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxYIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB


            if (Me%FluxXIntegral(i, j) < 0.) then 
                if      (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then

                    Area1    = Me%ExternalVar%DUX(i, j-1) * Me%ExternalVar%DVY(i, j-1)
                    Area2    = Me%ExternalVar%DUX(i, j  ) * Me%ExternalVar%DVY(i, j)
                    
                    SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
                    if (SandThickness > Me%SandMin) then
                        VolumeAvailable = SandThickness * Area2
                    else
                        VolumeAvailable = 0.
                    endif
                    
                    VolumeTransport = abs(Me%FluxXIntegral(i, j)) * Me%Evolution%DZDT
                    
                    if (VolumeTransport > VolumeAvailable) then
                        if (VolumeTransport > 0) then
                            Coef = VolumeAvailable / VolumeTransport
                        else
                            Coef = 0
                        endif
                        Me%FluxXIntegral(i, j) = Coef * Me%FluxXIntegral(i, j)
                        
                        VolumeTransport = VolumeAvailable

                    endif    

                    DZ                    = Me%DZ%Field2D(i, j-1) + VolumeTransport / Area1 
                    Me%DZ%Field2D(i, j-1) = DZ
                    DZ                    = Me%DZ%Field2D(i, j  ) - VolumeTransport / Area2
                    Me%DZ%Field2D(i, j  ) = DZ

                    if (Me%Boxes%Yes) then
                        Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j) + Me%FluxXIntegral(i, j) * Me%Evolution%DZDT
                    endif
                endif
            else 

                if (Me%ExternalVar%ComputeFacesU2D(i,j+1) == Covered) then

                    Area1 = Me%ExternalVar%DUX(i, j+1) * Me%ExternalVar%DVY(i, j+1)
                    Area2 = Me%ExternalVar%DUX(i, j  ) * Me%ExternalVar%DVY(i, j)
                    
                    SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
                    if (SandThickness > Me%SandMin) then
                        VolumeAvailable = SandThickness * Area2
                    else
                        VolumeAvailable = 0.
                    endif

                    
                    VolumeTransport = abs(Me%FluxXIntegral(i, j)) * Me%Evolution%DZDT
                    
                    if (VolumeTransport > VolumeAvailable) then
                        if (VolumeTransport > 0) then
                            Coef = VolumeAvailable / VolumeTransport
                        else
                            Coef = 0
                        endif
                        Me%FluxXIntegral(i, j) = Coef * Me%FluxXIntegral(i, j)
                        
                        VolumeTransport = VolumeAvailable

                    endif   
              
                    DZ                    = Me%DZ%Field2D(i, j+1) + VolumeTransport / Area1   
                    Me%DZ%Field2D(i, j+1) = DZ
                    DZ                    = Me%DZ%Field2D(i, j  ) - VolumeTransport / Area2 
                    Me%DZ%Field2D(i, j  ) = DZ

                    if (Me%Boxes%Yes) then
                        Me%Boxes%FluxesX(i,j+1) = Me%Boxes%FluxesX(i,j+1) + Me%FluxXIntegral(i, j) * Me%Evolution%DZDT
                    endif
                endif                    
            endif

            if (Me%FluxYIntegral(i, j) < 0.) then
                if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then

                    Area1 = Me%ExternalVar%DUX(i-1, j) * Me%ExternalVar%DVY(i-1, j)
                    Area2 = Me%ExternalVar%DUX(i, j  ) * Me%ExternalVar%DVY(i, j)
                    
                    SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
                    if (SandThickness > Me%SandMin) then
                        VolumeAvailable = SandThickness * Area2
                    else
                        VolumeAvailable = 0.
                    endif

                    
                    VolumeTransport = abs(Me%FluxYIntegral(i, j)) * Me%Evolution%DZDT
                    
                    if (VolumeTransport > VolumeAvailable) then
                        if (VolumeTransport > 0) then
                            Coef = VolumeAvailable / VolumeTransport
                        else
                            Coef = 0
                        endif
                        Me%FluxXIntegral(i, j) = Coef * Me%FluxXIntegral(i, j)
                        
                        VolumeTransport = VolumeAvailable

                    endif   

                    DZ                    = Me%DZ%Field2D(i-1, j) + VolumeTransport / Area1 
                    Me%DZ%Field2D(i-1, j) = DZ
                    DZ                    = Me%DZ%Field2D(i  , j) - VolumeTransport / Area2 
                    Me%DZ%Field2D(i, j  ) = DZ

                    if (Me%Boxes%Yes) then
                        Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j) + Me%FluxYIntegral(i, j) * Me%Evolution%DZDT
                    endif
                endif
            else 
                if (Me%ExternalVar%ComputeFacesV2D(i+1, j) == Covered) then

                    Area1 = Me%ExternalVar%DUX(i+1, j) * Me%ExternalVar%DVY(i+1, j)
                    Area2 = Me%ExternalVar%DUX(i, j  ) * Me%ExternalVar%DVY(i, j)
                    
                    SandThickness = Me%BedRock%Field2D(i, j  ) + Me%DZ_Residual%Field2D(i, j  )
                    if (SandThickness > Me%SandMin) then
                        VolumeAvailable = SandThickness * Area2
                    else
                        VolumeAvailable = 0.
                    endif

                    
                    VolumeTransport = abs(Me%FluxYIntegral(i, j)) * Me%Evolution%DZDT
                    
                    if (VolumeTransport > VolumeAvailable) then
                        if (VolumeTransport > 0) then
                            Coef = VolumeAvailable / VolumeTransport
                        else
                            Coef = 0
                        endif
                        Me%FluxXIntegral(i, j) = Coef * Me%FluxXIntegral(i, j)
                        
                        VolumeTransport = VolumeAvailable

                    endif   

                    DZ                    = Me%DZ%Field2D(i+1, j  ) + VolumeTransport / Area1
                    Me%DZ%Field2D(i+1, j) = DZ
                    DZ                    = Me%DZ%Field2D(i  , j  ) - VolumeTransport / Area2 
                    Me%DZ%Field2D(i  , j) = DZ

                    if (Me%Boxes%Yes) then
                        Me%Boxes%FluxesY(i+1,j) = Me%Boxes%FluxesY(i+1,j) + Me%FluxYIntegral(i, j) * Me%Evolution%DZDT
                    endif
                endif
            endif
            
            
        enddo
        enddo
        

        call ComputeDischarges
        
        
    end subroutine ComputeEvolution
    
    !--------------------------------------------------------------------------


    subroutine ComputeEvolution_CentralDifferences            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area1, DT_Area2, RunPeriod, Flux
        integer                            :: i, j
        !----------------------------------------------------------------------
        
        Me%DZ%Field2D(:,:) = 0.

!        if (Me%Boxes%Yes) then
!            Me%Boxes%FluxesX(:,:) = 0.
!            Me%Boxes%FluxesY(:,:) = 0.
!        endif

        if (Me%Aceleration%Yes) then

            Me%FluxXIntegral (:,:) = Me%Aceleration%Coef * Me%FluxXIntegral (:,:)
            Me%FluxYIntegral (:,:) = Me%Aceleration%Coef * Me%FluxYIntegral (:,:)

        endif
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime

        Me%Residual%FluxX(:,:) = ( Me%Residual%FluxX(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxXIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod

        Me%Residual%FluxY(:,:) = ( Me%Residual%FluxY(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxYIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB


            if (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then

                ![s/m^2]= [s] / [m] / [m] 
                DT_Area1 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j-1) / Me%ExternalVar%DVY(i, j-1)
                DT_Area2 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)
                
                ![m^3/s]             = ([m^3/s]*[m]+[m^3/s]*[m]) / ([m]+[m])
                Flux                 = (Me%FluxXIntegral  (i, j)   *  Me%ExternalVar%DUX(i, j-1) + &
                                        Me%FluxXIntegral  (i, j-1) *  Me%ExternalVar%DUX(i, j  ))/ &
                                       (Me%ExternalVar%DUX(i, j-1) +  Me%ExternalVar%DUX(i, j)) 
                ![m]
                Me%DZ%Field2D(i, j-1) = Me%DZ%Field2D(i, j-1) - DT_Area1 * Flux 
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  ) + DT_Area2 * Flux 

!                if (Me%Boxes%Yes) then
!                    Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j) + Flux * Me%Evolution%DZDT
!                endif
            endif

            if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then

                DT_Area1 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i-1, j) / Me%ExternalVar%DVY(i-1, j)
                DT_Area2 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                Flux                 = (Me%FluxXIntegral  (i, j)   *  Me%ExternalVar%DVY(i-1, j) + &
                                        Me%FluxXIntegral  (i-1, j) *  Me%ExternalVar%DVY(i, j  ))/ &
                                       (Me%ExternalVar%DVY(i-1, j) +  Me%ExternalVar%DVY(i, j)) 
                
                Me%DZ%Field2D(i-1, j) = Me%DZ%Field2D(i-1, j) - DT_Area1 * Flux 
                Me%DZ%Field2D(i  , j) = Me%DZ%Field2D(i  , j) + DT_Area2 * Flux 
                
                if (Me%Boxes%Yes) then
                    Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j) + Me%FluxYIntegral(i, j) * Me%Evolution%DZDT
                endif
            endif
            
        enddo
        enddo

        call ComputeDischarges
        
       
    end subroutine ComputeEvolution_CentralDifferences

    !--------------------------------------------------------------------------


    subroutine ComputeEvolution_UpwindOrder2            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area1, DT_Area2, RunPeriod, Flux, FluxWest, FluxEast
        integer                            :: i, j
        !----------------------------------------------------------------------
        
        Me%DZ%Field2D(:,:) = 0.

!        if (Me%Boxes%Yes) then
!            Me%Boxes%FluxesX(:,:) = 0.
!            Me%Boxes%FluxesY(:,:) = 0.
!        endif

        if (Me%Aceleration%Yes) then

            Me%FluxXIntegral (:,:) = Me%Aceleration%Coef * Me%FluxXIntegral (:,:)
            Me%FluxYIntegral (:,:) = Me%Aceleration%Coef * Me%FluxYIntegral (:,:)

        endif
        
        RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime

        Me%Residual%FluxX(:,:) = ( Me%Residual%FluxX(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxXIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod

        Me%Residual%FluxY(:,:) = ( Me%Residual%FluxY(:,:) * (RunPeriod -  Me%Evolution%DZDT)          + &
                                   Me%FluxYIntegral(:,:) * Me%Evolution%DZDT) / RunPeriod
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB


            if (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then

                ![s/m^2]= [s] / [m] / [m] 
                DT_Area1 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j-1) / Me%ExternalVar%DVY(i, j-1)
                DT_Area2 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)
                
                !Upwind second order                    
                if (Me%ExternalVar%OpenPoints2D(i,  j-2) == 1 .and.     &
                    Me%ExternalVar%OpenPoints2D(i,  j-1) == 1 .and.     &                    
                    Me%ExternalVar%OpenPoints2D(i,  j  ) == 1 .and.     &                                        
                    Me%ExternalVar%OpenPoints2D(i,  j+1) == 1) then

                    FluxWest = -1./8. * Me%FluxXIntegral(i, j-2) + 6./8. * Me%FluxXIntegral(i, j-1) &
                               +3./8. * Me%FluxXIntegral(i, j  )
                    FluxEast = -1./8. * Me%FluxXIntegral(i, j+1) + 6./8. * Me%FluxXIntegral(i, j  ) &
                               +3./8. * Me%FluxXIntegral(i, j-1)
                    
                    if (0.5*(FluxWest+FluxEast) > 0) then
                        Flux = FluxWest
                    else
                        Flux = FluxEast
                    endif

                !Upwind first order
                else
                    Flux = 0
                    if (Me%FluxXIntegral(i, j-1) > 0) then
                        Flux = Flux + Me%FluxXIntegral(i, j-1)
                    endif
                    if (Me%FluxXIntegral(i, j  ) < 0) then
                        Flux = Flux + Me%FluxXIntegral(i, j  )
                    endif
                    
                endif
                                       
                ![m]
                Me%DZ%Field2D(i, j-1) = Me%DZ%Field2D(i, j-1) - DT_Area1 * Flux 
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  ) + DT_Area2 * Flux 

!                if (Me%Boxes%Yes) then
!                    Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j) + Flux * Me%Evolution%DZDT
!                endif
            endif

            if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then

                DT_Area1 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i-1, j) / Me%ExternalVar%DVY(i-1, j)
                DT_Area2 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                Flux                 = (Me%FluxXIntegral  (i, j)   *  Me%ExternalVar%DVY(i-1, j) + &
                                        Me%FluxXIntegral  (i-1, j) *  Me%ExternalVar%DVY(i, j  ))/ &
                                       (Me%ExternalVar%DVY(i-1, j) +  Me%ExternalVar%DVY(i, j)) 
                
                Me%DZ%Field2D(i-1, j) = Me%DZ%Field2D(i-1, j) - DT_Area1 * Flux 
                Me%DZ%Field2D(i  , j) = Me%DZ%Field2D(i  , j) + DT_Area2 * Flux 
                
                if (Me%Boxes%Yes) then
                    Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j) + Me%FluxYIntegral(i, j) * Me%Evolution%DZDT
                endif
            endif
            
        enddo
        enddo

        call ComputeDischarges
        
       
    end subroutine ComputeEvolution_UpwindOrder2
    

    
    !--------------------------------------------------------------------------



    subroutine ComputeNumericalSmoothing            

        !Local-----------------------------------------------------------------
        real                               :: NumericalSmoothing
        integer                            :: i, j
        !----------------------------------------------------------------------
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB


            if (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then

                ![m]                       
                NumericalSmoothing    =  -(Me%ExternalVar%Bathymetry(i, j) -            &
                    Me%ExternalVar%Bathymetry(i, j-1)) * 0.5 * Me%Evolution%Gama
                                       
                ![m]
                Me%DZ%Field2D(i, j-1) = Me%DZ%Field2D(i, j-1)  + NumericalSmoothing
                Me%DZ%Field2D(i, j  ) = Me%DZ%Field2D(i, j  )  - NumericalSmoothing

            endif

            if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then

                ![m]                       
                NumericalSmoothing    =  -(Me%ExternalVar%Bathymetry(i, j) -            &
                    Me%ExternalVar%Bathymetry(i-1, j)) * 0.5 * Me%Evolution%Gama

                Me%DZ%Field2D(i-1, j) = Me%DZ%Field2D(i-1, j) + NumericalSmoothing
                Me%DZ%Field2D(i  , j) = Me%DZ%Field2D(i  , j) - NumericalSmoothing
                
            endif
        enddo
        enddo

        
       
    end subroutine ComputeNumericalSmoothing
    
    !--------------------------------------------------------------------------

    subroutine ComputeDischarges            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area, BottomDryDensity, TicknessAdd,   &
                                              DischargeFlow, DischargeConc
        integer                            :: i, j, DischargesNumber, STAT_CALL
        integer                            :: DischargeID 
        logical                            :: PropFromIntake, IgnoreOK, CoordinatesON
        real                               :: CoordinateX, CoordinateY
        
        
        !----------------------------------------------------------------------

        
cd1:    if (Me%Discharges%Yes) then

            !Sinks and Sources
            call GetDischargesNumber(Me%ObjDischarges, DischargesNumber, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ComputeDischarges - ModuleSand - ERR10'


            !For all Discharges
d1:         do DischargeID = 1, DischargesNumber
            
                call GetDischargeON(Me%ObjDischarges,DischargeID, IgnoreOK, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeDischarges - ModuleSand - ERR20'

                if (IgnoreOK) cycle

                call GetDischargesGridLocalization(Me%ObjDischarges,                    &
                                                   DischargeID,                             &
                                                   Igrid         = I,                       &
                                                   JGrid         = J,                       &
                                                   CoordinateX   = CoordinateX,             &
                                                   CoordinateY   = CoordinateY,             &
                                                   CoordinatesON = CoordinatesON,           &
                                                   TimeX         = Me%ExternalVar%Now,      &
                                                   STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'ComputeDischarges - ModuleSand - ERR30'

                if (CoordinatesON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordinateX, CoordinateY, I, J, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeDischarges - ModuleSand - ERR40'

                    call CorrectsCellsDischarges(Me%ObjDischarges, DischargeID, I, J, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ComputeDischarges - ModuleSand - ERR50'

                endif


                call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                           Me%ExternalVar%Now,                          &
                                           DischargeID, -99., DischargeFlow,            &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ComputeDischarges - ModuleSand - ERR60'



                call GetDischargeConcentration (Me%ObjDischarges,                       &
                                                Me%ExternalVar%Now,                     &
                                                DischargeID, DischargeConc,                     &
                                                PropertyIDNumber    = Sand_,            &
                                                PropertyFromIntake  = PropFromIntake,   &
                                                STAT                = STAT_CALL)                                                
                                                
                if (STAT_CALL/=SUCCESS_)                                                &
                    stop 'ComputeDischarges - ModuleSand - ERR70'
                

                BottomDryDensity    = (1. - Me%Porosity) * Me%Density

                DT_Area             = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                ![m]                = [m3/s] * [g/m3/1000] / [kg/m3] * [s/m2] * []
                TicknessAdd         = DischargeFlow * (DischargeConc/1000.) / BottomDryDensity * DT_Area

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

            if (Me%SmoothSlope%DZ_Residual) then
                    
                dhdx = (Me%DZ_Residual%Field2D(i, j-1) - Me%DZ_Residual%Field2D(i, j)) /  &
                        Me%ExternalVar%DZX(i,j-1)                    
                    
                dhdy = (Me%DZ_Residual%Field2D(i-1, j) - Me%DZ_Residual%Field2D(i, j)) /  &
                        Me%ExternalVar%DZY(i-1,j)                    
            else    

                if (Me%Evolution%Bathym) then

            dhdx = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i, j-1)) /  &
                    Me%ExternalVar%DZX(i,j-1)

                    dhdy = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i-1, j)) /  &
                            Me%ExternalVar%DZY(i-1,j)                
                    
                else
                    !bathymetry gradients can not feed transport when the bathymetry can not evolve
                    dhdx = 0.
                    dhdy = 0.                     
                endif
                
            endif
                

            if    ( Me%ExternalVar%WaterPoints2D(i,  j  ) == WaterPoint  .and. &
                    Me%ExternalVar%WaterPoints2D(i,  j-1) == WaterPoint  .and. &
                    abs(dhdx) > Me%SmoothSlope%Critic) then
                
                DT_Area1 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j-1) / Me%ExternalVar%DVY(i, j-1)
                DT_Area2 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                FluxX = Me%SmoothSlope%Factor * max(abs(Me%FluxY(i, j)), abs(Me%FluxY(i, j-1)), &
                                                    abs(Me%FluxX(i, j)), abs(Me%FluxX(i, j-1)))

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


            if    ( Me%ExternalVar%WaterPoints2D(i,   j) == WaterPoint  .and. &
                    Me%ExternalVar%WaterPoints2D(i-1, j) == WaterPoint  .and. &
                    abs(dhdy) > Me%SmoothSlope%Critic) then
                
                DT_Area1 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i-1, j) / Me%ExternalVar%DVY(i-1, j)
                DT_Area2 = Me%Evolution%DZDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                FluxY = Me%SmoothSlope%Factor * max(abs(Me%FluxX(i, j)), abs(Me%FluxX(i-1, j)), &
                                                    abs(Me%FluxY(i, j)), abs(Me%FluxY(i-1, j)))

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

!--------------------------------------------------------------------------

    subroutine ComputeBathymetryGradient           

        !Local-----------------------------------------------------------------
        real                               :: dhdx1, dhdy1, dhdx2, dhdy2 
        integer                            :: i, j
        !----------------------------------------------------------------------
        
        
        
        if (.not.associated(Me%BatGradient)) then
            allocate(Me%BatGradient(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))    
        endif

        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB
            
            if (Me%ExternalVar%WaterPoints2D(i,  j  ) == WaterPoint) then
        
                if (Me%ExternalVar%WaterPoints2D(i,j-1) == WaterPoint) then
        
                    dhdx1 = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i, j-1)) /  &
                            Me%ExternalVar%DZX(i,j-1)

                else
                
                    dhdx1 = 0.
                            
                endif                            
                

                if (Me%ExternalVar%WaterPoints2D(i,j+1) == WaterPoint) then
        
                    dhdx2 = (Me%ExternalVar%Bathymetry(i, j+1) - Me%ExternalVar%Bathymetry(i, j)) /  &
                            Me%ExternalVar%DZX(i,j)

                else
                
                    dhdx2 = 0.
                            
                endif                            
                
                if (Me%ExternalVar%WaterPoints2D(i-1,j) == WaterPoint) then
        
                    dhdy1 = (Me%ExternalVar%Bathymetry(i, j) - Me%ExternalVar%Bathymetry(i-1, j)) /  &
                            Me%ExternalVar%DZY(i-1,j)

                else
                
                    dhdy1 = 0.
                            
                endif                            
                

                if (Me%ExternalVar%WaterPoints2D(i+1,j) == WaterPoint) then
        
                    dhdy2 = (Me%ExternalVar%Bathymetry(i+1, j) - Me%ExternalVar%Bathymetry(i, j)) /  &
                            Me%ExternalVar%DZY(i,j)

                else
                
                    dhdy2 = 0.
                            
                endif   
                
                Me%BatGradient(i, j) = max(abs(dhdx1),abs(dhdx1),abs(dhdy1),abs(dhdy2))                 
                
            else
                            
                Me%BatGradient(i, j) = FillValueReal
                
            endif

                                   
        enddo
        enddo
        
    end subroutine ComputeBathymetryGradient
    
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
            
        else if (Me%Boundary == NullValue) then

            where (Me%ExternalVar%BoundaryPoints2D(:, :) == Boundary) 
                Field2D(:, :) = 0.
            end where            

        endif


    end subroutine BoundaryCondition
    
    !--------------------------------------------------------------------------


    subroutine BoundaryCondition1D(Field1D, ILB, IUB)
        !Arguments-------------------------------------------------------------
        real(8),   dimension(:  ), pointer    :: Field1D
        integer                               :: ILB, IUB

        !Begin-----------------------------------------------------------------

        if      (Me%Boundary == NullGradient) then

            Field1D(ILB) = Field1D(ILB+1)
            Field1D(IUB) = Field1D(IUB-1)

        else if (Me%Boundary == Cyclic      ) then

            Field1D(ILB) = Field1D(IUB-1)
            Field1D(IUB) = Field1D(ILB+1)
            
        else if (Me%Boundary == NullValue   ) then
        
            Field1D(ILB) = 0.
            Field1D(IUB) = 0.            
        
        endif


    end subroutine BoundaryCondition1D
    
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
        
        if (Me%BiHarmonicFilter) then
            call ComputeBiHamonicFilter2D(Prop2D = Me%DZ_Residual%Field2D, Coef = Me%BiHarmonicFilterCoef)                
        endif

!        if (.not. Me%HybridMorph%ON) then
!
!            if (Me%Evolution%BathymDT > Me%Evolution%DZDT) then
!
!                do j=Me%WorkSize%JLB, Me%WorkSize%JUB
!                do i=Me%WorkSize%ILB, Me%WorkSize%IUB
!
!                    if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
!                    
!                        !Me%BatimIncrement%Field2D(i,j)  = Me%BatimIncrement%Field2D(i,j) + Me%DZ%Field2D(i, j)
!
!                    endif
!
!                enddo
!                enddo    
!
!            endif
!
!        endif
        
        
    end subroutine ComputeResidualEvolution
    
    !--------------------------------------------------------------------------



    subroutine OutPutSandHDF
        
        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
         
        !Local-----------------------------------------------------------------
        logical                            :: FirstTime
        integer                            :: OutPutNumber
        type (T_Time)                      :: Actual
        integer                            :: ILB, IUB, JLB, JUB, i, j
        real,    dimension(6    ), target  :: AuxTime
        real,    dimension(:    ), pointer :: TimePtr
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        real,    dimension(:,:  ), pointer :: Array2D

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

                allocate(Array2D(ILB:IUB,JLB:JUB))
            
                !Writes current time
                call ExtractDate   (Actual, AuxTime(1), AuxTime(2), AuxTime(3),         &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR10'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",&
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR20'

                !Writes OpenPoints
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,             &
                                     WorkJUB, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR30'

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",      &
                                     "-", Array2D = Me%ExternalVar%OpenPoints2D,        &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR40'
                
                call ComputeBathymetryGradient
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Bathymetry Gradient", "Bathymetry Gradient",   &
                                     "-", Array2D = Me%BatGradient,          &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR50'
                

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Bathymetry", "Bathymetry",   &
                                     "-", Array2D = Me%ExternalVar%Bathymetry,          &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR50'
       
                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%DZ%ID%Name), trim(Me%DZ%ID%Name),  &
                                     trim(Me%DZ%ID%Units), Array2D = Me%DZ%Field2D,     &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR60'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%DZ_Residual%ID%Name), trim(Me%DZ_Residual%ID%Name),  &
                                     trim(Me%DZ_Residual%ID%Units), Array2D = Me%DZ_Residual%Field2D,   &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR70'

                do j = WorkJLB, Me%WorkSize%JUB
                do i = WorkILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then
                        Array2D(i, j) = Me%DZ_Residual%Field2D(i, j) + Me%BedRock%Field2D(i, j)
                    else
                        Array2D(i, j) = 0.
                    endif                        
                enddo
                enddo

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//"Sand Thickness",         &
                                    "Sand Thickness", "m", Array2D = Array2D,    &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR75'                                

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Capacity", "Transport  Capacity",  &
                                     "m3/s/m", Array2D = Me%TransportCapacity,               &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR80'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Tau Critic", "Tau Critic",        &
                                     "N/m2", Array2D = Me%TauCritic,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR90'
                

                call RotateVectorGridToField(HorizontalGridID  = Me%ObjHorizontalGrid,  &
                                             VectorInX         = Me%FluxX,              &
                                             VectorInY         = Me%FluxY,              &
                                             VectorOutX        = Me%OutFluxX,           &
                                             VectorOutY        = Me%OutFluxY,           &
                                             WaterPoints2D     = Me%ExternalVar%WaterPoints2D,&
                                             RotateX           = .true.,                &
                                             RotateY           = .true.,                &
                                             STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR100'


                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux X", "Transport Flux X",&
                                     "m3/s", Array2D = Me%OutFluxX,                     &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR110'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux Y", "Transport Flux Y",&
                                     "m3/s", Array2D = Me%OutFluxY,                     &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR120'

                do j = WorkJLB, Me%WorkSize%JUB
                do i = WorkILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then
                        Array2D(i, j) = sqrt(Me%OutFluxX(i,j)**2. + Me%OutFluxY(i,j)**2.)
                    else
                        Array2D(i, j) = 0.
                    endif                        
                enddo
                enddo
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux", "Transport Flux",&
                                     "m3/s", Array2D = Array2D,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR125'

                call RotateVectorGridToField(HorizontalGridID  = Me%ObjHorizontalGrid,  &
                                             VectorInX         = Me%Residual%FluxX,     &
                                             VectorInY         = Me%Residual%FluxY,     &
                                             VectorOutX        = Me%Residual%OutFluxX,  &
                                             VectorOutY        = Me%Residual%OutFluxY,  &
                                             WaterPoints2D     = Me%ExternalVar%WaterPoints2D,&
                                             RotateX           = .true.,                &
                                             RotateY           = .true.,                &
                                             STAT              = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR130'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Residual Transport Flux X",  &
                                     "Residual Transport Flux X",                       &
                                     "m3/s", Array2D = Me%Residual%OutFluxX,            &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR140'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Residual Transport Flux Y",  &
                                     "Residual Transport Flux Y",                       &
                                     "m3/s", Array2D = Me%Residual%OutFluxY,            &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR150'
                
                do j = WorkJLB, Me%WorkSize%JUB
                do i = WorkILB, Me%WorkSize%IUB
                    if (Me%ExternalVar%WaterPoints2D(i, j) == 1) then
                        Array2D(i, j) = sqrt(Me%Residual%OutFluxX(i,j)**2. + Me%Residual%OutFluxY(i,j)**2.)
                    else
                        Array2D(i, j) = 0.
                    endif                        
                enddo
                enddo
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Residual Transport Flux",    &
                                     "Residual Transport Flux",                         &
                                     "m3/s", Array2D = Array2D,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR155'                
                
                if (Me%HybridMorph%ON) then
                    
                    call HDF5WriteData  (Me%ObjHDF5,                                    &
                                         "/Results/HybridMorph/"//trim(Me%HybridMorph%DZ_Residual%ID%Name), &
                                         trim(Me%HybridMorph%DZ_Residual%ID%Name),         &
                                         trim(Me%HybridMorph%DZ_Residual%ID%Units),        &
                                         Array2D = Me%HybridMorph%DZ_Residual%Field2D,  &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR160'
                    

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/HybridMorph/Distance to Coast", &
                                         "Distance to Coast",                                &
                                         "m", Array2D = Me%HybridMorph%DistanceToCoastInst,           &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR220'                         
                        
                    call HDF5SetLimits  (Me%ObjHDF5, Me%HybridMorph%Min1D, Me%HybridMorph%Max1D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR170'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/HybridMorph/Along Shore Flux", &
                                         "Along Shore Flux",                                  &
                                         "m3/s", Array1D = Me%HybridMorph%AlongShoreFlux,     &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR180'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/HybridMorph/Residual Along Shore Flux", &
                                         "Residual Along Shore Flux",                         &
                                         "m3/s", Array1D = Me%HybridMorph%ResidualAlongShoreFlux,           &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR190'
                    
                    call HDF5WriteData  (Me%ObjHDF5, "/Results/HybridMorph/Cross Shore Velocity", &
                                         "Cross Shore Velocity",                                   &
                                         "m/s", Array1D = Me%HybridMorph%CrossShoreVel,     &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR200'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/HybridMorph/Residual Cross Shore Velocity", &
                                         "Residual Cross Shore Velocity",                                &
                                         "m/s", Array1D = Me%HybridMorph%ResidualCrossShoreVel,           &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR210'

              

                endif
                
                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSandHDF - ModuleSand - ERR230'

                Me%OutPut%NextOutPut = OutPutNumber + 1
                
                deallocate(Array2D)

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
            Me%Boxes%Mass   (i,j) = Me%DZ_Residual%Field2D(i,j) *                       &
                                    Me%ExternalVar%DUX(i,j)     * Me%ExternalVar%DVY(i,j) 

            !This fluxes are compute from the residual fluxes in m3/s
            Me%Boxes%FluxesX(i,j) = Me%Residual%FluxX(i,j)   
            Me%Boxes%FluxesY(i,j) = Me%Residual%FluxY(i,j)   

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
                
                if (Me%AverageON) then
                    deallocate(Me%Uaverage%Field2D)
                    deallocate(Me%Vaverage%Field2D)                
                endif
                
                if (Me%MappDZON) then
                    deallocate(Me%MappDZ%Field2D)
                endif                

                deallocate(Me%BedRock%Field2D)
                deallocate(Me%TauCritic)

                deallocate(Me%DZ%Field2D)
                
!                if (Me%Evolution%BathymDT > Me%Evolution%DZDT) then
!                
!                    deallocate(Me%BatimIncrement%Field2D)
!
!                endif

                deallocate(Me%DZ_Residual%Field2D)

                if (Me%Filter%ON) deallocate (Me%Filter%Field2D)

                !deallocate fluxes
                deallocate(Me%FluxX  )
                deallocate(Me%FluxY  )
                deallocate(Me%FluxXIntegral)
                deallocate(Me%FluxYIntegral)
                deallocate(Me%TransportCapacity)
                
                deallocate(Me%Residual%FluxX  )
                deallocate(Me%Residual%FluxY  )
                
                deallocate(Me%OutFluxX  )
                deallocate(Me%OutFluxY  )

                deallocate(Me%Residual%OutFluxX  )
                deallocate(Me%Residual%OutFluxY  )
                
                !if (Me%TransportMethod == VanRijn1) deallocate (Me%Dast)
                
                if   (associated(Me%Dast)) deallocate (Me%Dast)

                if (Me%Boxes%Yes) then

                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSand - ModuleSand - ERR30'

                    deallocate(Me%Boxes%FluxesX)
                    nullify   (Me%Boxes%FluxesX)

                    deallocate(Me%Boxes%FluxesY)
                    nullify   (Me%Boxes%FluxesY)


                    deallocate(Me%Boxes%Mass   )
                    nullify   (Me%Boxes%Mass   )
                    
                    if (Me%MudErosion%ON) then
                        deallocate(Me%Boxes%FluxesZ)
                        nullify   (Me%Boxes%FluxesZ)
                    endif                    

                endif
                
                if (Me%HybridMorph%ON) then
                    call KillHybridMorph
                endif
                
                if (associated(Me%BatGradient)) then
                    deallocate(Me%BatGradient)
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
    
    
    !----------------------------------------------------------------------------
    
    subroutine KillHybridMorph

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------
        !Deallocate variables
        !1D vectors
        deallocate(Me%HybridMorph%InShoreMapping        )
        deallocate(Me%HybridMorph%OffShoreMapping       ) 
        deallocate(Me%HybridMorph%AlongShoreFlux        ) 
        deallocate(Me%HybridMorph%ResidualAlongShoreFlux) 
        deallocate(Me%HybridMorph%CrossShoreVel         ) 
        deallocate(Me%HybridMorph%ResidualCrossShoreVel ) 


        deallocate(Me%HybridMorph%BathymetryNext        )
        deallocate(Me%HybridMorph%BathymetryPrevious    )
        deallocate(Me%HybridMorph%DistanceToCoastRef    )
        deallocate(Me%HybridMorph%DistanceToCoastInst   )           
        
        deallocate(Me%HybridMorph%DZ_Residual%Field2D)
        
    end subroutine KillHybridMorph
            
    !--------------------------------------------------------------------------    
    !--------------------------------------------------------------------------
    !   Write the final water properties results in HDF format  !
  
    subroutine WriteFinalState

        !Local--------------------------------------------------------------
        real,    dimension(6    ), target      :: AuxTime
        real,    dimension(:    ), pointer     :: TimePtr
        integer                                :: WorkILB, WorkIUB
        integer                                :: WorkJLB, WorkJUB
        integer                                :: STAT_CALL
        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        
        call Open_HDF5_OutPut_File(Me%Files%FinalSand)

       
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                              &
                             WorkJLB, WorkJUB,                                          &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR10'

        !Final concentration
        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                    &
                             trim(Me%DZ_Residual%ID%Name),                              &
                             trim(Me%DZ_Residual%ID%Units),                             &
                             Array2D = Me%DZ_Residual%Field2D,                          &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR20'

        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                    &
                             "Bathymetry",                                              &
                             "m",                                                       &
                             Array2D = Me%ExternalVar%Bathymetry,                       &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR30'
            
        call RotateVectorGridToField(HorizontalGridID  = Me%ObjHorizontalGrid,          &
                                     VectorInX         = Me%Residual%FluxX,             &
                                     VectorInY         = Me%Residual%FluxY,             &
                                     VectorOutX        = Me%Residual%OutFluxX,          &
                                     VectorOutY        = Me%Residual%OutFluxY,          &
                                     WaterPoints2D     = Me%ExternalVar%WaterPoints2D,  &
                                     RotateX           = .true.,                        &
                                     RotateY           = .true.,                        &
                                     STAT              = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR40'

        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                    &
                             "Residual Transport Flux X",                               &
                             "m3/s/m",                                                  &
                             Array2D = Me%Residual%OutFluxX,                            &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR50'

        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                    &
                             "Residual Transport Flux Y",                               &
                             "m3/s/m",                                                  &
                             Array2D = Me%Residual%OutFluxY,                            &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR60'

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR70'

        !Writes current time
        call ExtractDate   (Me%Residual%StartTime, AuxTime(1), AuxTime(2), AuxTime(3),  &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime
        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalState - ModuleSand - ERR80'

        call HDF5WriteData  (Me%ObjHDF5, "/Time",                                       &
                             "Residual Start Time", "YYYY/MM/DD HH:MM:SS",              &
                             Array1D = TimePtr, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFinalState - ModuleSand - ERR90'
                             
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR100'

        call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'WriteFinalState - ModuleSand - ERR110'

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











