!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Sediment
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2015
! REVISION      : Guilherme Franz
! DESCRIPTION   : Module to compute sediment evolution
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

Module ModuleSediment

    use ModuleGlobalData
    use ModuleStopWatch,        only : StartWatch, StopWatch
    use ModuleFunctions         
    use ModuleTime              
    use ModuleDrawing 
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
                                       GetNumberOfTimeSeries, TryIgnoreTimeSerie, GetTimeSerieName       
    use ModuleHorizontalMap,    only : GetWaterPoints2D, GetBoundaries, GetOpenPoints2D,        &
                                       GetComputeFaces2D, UnGetHorizontalMap
    use ModuleMap,              only:  GetWaterPoints3D, GetOpenPoints3D,                &
                                       GetLandPoints3D, UngetMap, UpdateComputeFaces3D
    use ModuleHorizontalGrid,   only : GetHorizontalGrid, WriteHorizontalGrid,                  &
                                       GetHorizontalGridSize, UnGetHorizontalGrid, GetXYCellZ,  &
                                       GetGridCellArea, GetDDecompMPI_ID, GetDDecompON,         &
                                       GetGridOutBorderPolygon, RotateVectorGridToField
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes, BoxDif,         &
                                       UngetBoxDif, KillBoxDif
    use ModuleGeometry,         only:  GetGeometrySize, UnGetGeometry,                      &
                                       GetGeometryDistances, GetGeometryKtop,               &
                                       ComputeInitialGeometry, ComputeVerticalGeometry,     &
                                       GetGeometryVolumes, ReadGeometryHDF, WriteGeometryHDF
    use ModuleFreeVerticalMovement,  only: SetSandParameters
    
#ifndef _WAVES_
    use ModuleWaves,            only : GetWaves, UnGetWaves
#endif

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructSediment
    private ::      AllocateInstance
    private ::      Read_Sediment_Files_Name
    private ::      ConstructGlobalParameters
    private ::      AllocateVariables
    private ::      ConstructEvolution
    private ::      Construct_Initial_Geometry
    private ::      ConstructClasses
    private ::          ConstructPropertyValue
    private ::          ConstructCohesiveDryDensity
    private ::          CheckPercentageConsistence
    private ::      ConstructD50Cell
    private ::      ConstructPorosity
    private ::      ConstructClassMass
   
    private ::      StartOutputBoxFluxes
    private ::      Open_HDF5_OutPut_File
    private ::      ConstructTimeSerie
    
    private ::      ReadInitialField
    private ::      ReadInitialField3D


    !Selector
    public  :: SetCohesiveFlux
    public  :: SetNonCohesiveFlux
    public  :: GetTopCriticalShear
    public  :: GetCohesiveMass
    public  :: GetCohesiveContent
    public  :: GetSandMass
    public  :: GetSandContent
    public  :: GetCriticalShearStress
    public  :: GetConcRef
    public  :: GetReferenceLevel
    !public  :: GetSandD50
    public  :: SetWaveTensionON
    public  :: UnGetSediment
    
    !Modifier
    public  :: ModifySediment
    private ::      ComputeOpenSediment
    private ::      ComputeD50Cell
    private ::      ComputeNDShearStress
    private ::      ComputeCriticalShearStress
    private ::      ComputeReferenceConcentration
    private ::      ComputeBedload
    private ::          Cohesiveness_Adjust
    private ::          CurrentOnly
    private ::          CurrentPlusWaves
    private ::          WavesOnly
    private ::      ComputeFluxes
    private ::          ComputeBedSlopeEffects
    private ::          FluxesCorrection
    private ::      ComputeEvolution
    private ::      ComputeDischarges   
    private ::      BoundaryCondition
    private ::      ComputeMass
    private ::      ComputePercentage
    private ::          ComputePorosity
    private ::      ComputeTotalDZ
    private ::      ComputeVerticalCoordinate
    private ::          New_Layer
    private ::              MaxLayersExceeded
    private ::          Activate_Layer
    private ::          Layer_Eroded
    private ::      New_Geometry
    private ::      ComputeResidualEvolution

    private ::      OutPutSedimentHDF
    private ::      OutPut_TimeSeries
    private ::      OutputBoxFluxes
                
    !Destructor
    public  :: KillSediment                                                     
    private ::      DeAllocateInstance
    private ::      WriteFinalState

    !Management
    private ::      Ready
    private ::      LocateObjSediment
    private ::      ReadLockExternalVar
    private ::      ReadUnLockExternalVar

    !Interfaces----------------------------------------------------------------
    private :: UnGetSediment2D_I
    private :: UnGetSediment2D_R8
    private :: UnGetSediment2D_R4
    private :: UnGetSediment3Dreal8
    private :: UnGetSediment3Dreal4
    interface  UnGetSediment
        module procedure UnGetSediment2D_I
        module procedure UnGetSediment2D_R4
        module procedure UnGetSediment2D_R8
        module procedure UnGetSediment3Dreal8
        module procedure UnGetSediment3Dreal4
    end interface  UnGetSediment

    !Parameters
    character(LEN = StringLength), parameter    :: class_block_begin     = '<beginsandclass>'
    character(LEN = StringLength), parameter    :: class_block_end       = '<endsandclass>'
    
    character(LEN = StringLength), parameter    :: cohesive_block_begin     = '<begincohesive>'
    character(LEN = StringLength), parameter    :: cohesive_block_end       = '<endcohesive>'
    
    character(LEN = StringLength), parameter    :: cohesivedensity_block_begin     = '<begindrydensity>'
    character(LEN = StringLength), parameter    :: cohesivedensity_block_end       = '<enddrydensity>'

    integer, parameter :: NoTransport = 0, CurrentOnly = 1, CurrentPlusWaves = 2, WavesOnly = 3
    
    integer, parameter :: NullGradient = 1, Cyclic = 2, NullValue = 3

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
        integer, pointer, dimension(:,:,:)      :: OpenPoints3D     => null()
        integer, pointer, dimension(:,:)        :: WaterPoints2D    => null()
        integer, pointer, dimension(:,:,:)      :: WaterPoints3D    => null()
        integer, pointer, dimension(:,:)        :: BoundaryPoints2D => null()
        real                                    :: WaterDensity     = FillValueReal 
        logical                                 :: WaveTensionON    = .false. 
        real,    pointer, dimension(:,:)        :: Bathymetry       => null()
        real,    pointer, dimension(:,:)        :: InitialBathym    => null()
        real,    pointer, dimension(:,:)        :: GridCellArea     => null()
        real,    pointer, dimension(:,:)        :: WaveDirection    => null()
        real,    pointer, dimension(:,:)        :: Abw              => null()
        real,    pointer, dimension(:,:)        :: Ubw              => null()
        real,    pointer, dimension(:,:)        :: ShearStress      => null()
        real,    pointer, dimension(:,:)        :: ShearStressMean  => null()
        real,    pointer, dimension(:,:)        :: VelU             => null()
        real,    pointer, dimension(:,:)        :: VelV             => null()
        real,    pointer, dimension(:,:)        :: VelMod           => null()
        real,    pointer, dimension(:,:)        :: TauWave          => null()
        real,    pointer, dimension(:,:)        :: Cphi             => null()
        real,    pointer, dimension(:,:)        :: CWphi            => null()
        real,    pointer, dimension(:,:)        :: Wphi             => null()
        real,    pointer, dimension(:,:)        :: WaveHeight       => null()
        real,    pointer, dimension(:,:)        :: WavePeriod       => null()
        !Sediment
        real,    pointer, dimension(:,:,:)      :: DWZ              => null()
        real,    pointer, dimension(:,:,:)      :: SZZ              => null()
        integer, pointer, dimension(:,:  )      :: KTop             => null()
        real(8),    pointer, dimension(:,:,:)   :: VolumeZ          => null()
        
        real,    pointer, dimension(:,:  )      :: ConsolidationFlux    => null()

    end type T_External

    private :: T_Residual
    type       T_Residual
        logical                                 :: ON           = .false.
        type(T_Time)                            :: StartTime    
        real, dimension(:,:), pointer           :: FluxX        => null ()
        real, dimension(:,:), pointer           :: FluxY        => null ()        
    end type   T_Residual

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

    private :: T_Evolution
    type       T_Evolution
        logical                                 :: Old          = .false. 
        real                                    :: SedimentDT       = FillValueReal
        real                                    :: BathymDT     = FillValueReal
        type (T_Time)                           :: NextSediment, NextBatim
        logical                                 :: Bathym       = .false. 
    end type T_Evolution
     
    private :: T_Property
    type ::       T_Property
        type (T_PropertyID)                     :: ID
!        type (T_SubModel  )                     :: SubModel
        real                                    :: Scalar       = FillValueReal
        real, pointer, dimension(:,:,:)         :: Field3D      => null ()
    end type T_Property

    public :: T_Sand
    type, extends (T_Property) :: T_Sand
        real(8)                                 :: D50                  = FillValueReal
        real                                    :: parameter_A          = FillValueReal, &
                                                   parameter_n          = FillValueReal, &
                                                   parameter_p          = FillValueReal
                                                               
        real(8), dimension(:,:), pointer        :: FluxX                 => null ()
        real(8), dimension(:,:), pointer        :: FluxY                 => null ()
        real,    dimension(:,:), pointer        :: FluxToSediment        => null ()
        real(8), dimension(:,:), pointer        :: Bedload               => null ()
        real(8), dimension(:,:), pointer        :: BedloadX              => null ()
        real(8), dimension(:,:), pointer        :: BedloadY              => null ()
        real,    dimension(:,:), pointer        :: CriticalShearStress   => null ()
        real, dimension(:,:), pointer           :: NDCriticalShearStress => null ()
        real(8), dimension(:,:), pointer        :: DM                    => null ()
        real(8), dimension(:,:,:), pointer      :: Mass                  => null ()
        real(8), dimension(:,:), pointer        :: TopPercentage         => null ()
        real, dimension(:,:), pointer           :: HidingFactor          => null ()
        real, dimension(:,:), pointer           :: ConcRef               => null ()
    !Field3D in the base class T_Property is the percentage (content by dry weight) 
    !of this class on each domain cell
    end type T_Sand
    
    public :: T_Cohesive
    type, extends (T_Property) :: T_Cohesive
        logical                                 :: Run                  = .false.
        real(8)                                 :: D50                  = FillValueReal
        real                                    :: WeakConsolidated_CSS = FillValueReal
        real                                    :: PM1_MAX              = FillValueReal
        real                                    :: PM2                  = FillValueReal
        
        real,    dimension(:,:), pointer        :: CriticalShearStress   => null ()
        real, dimension(:,:), pointer           :: PM1                   => null ()
        real(8), dimension(:,:), pointer        :: DM                    => null ()
        real(8), dimension(:,:,:), pointer      :: Mass                  => null ()
        real(8), dimension(:,:), pointer        :: TopPercentage         => null ()
        real, dimension(:,:,:), pointer         :: Porosity              => null ()
    !Field3D in the base class T_Property is the percentage (content by dry weight) 
    !of this class on each domain cell
    end type T_Cohesive
    
    public :: T_CohesiveDryDensity
    type, extends (T_Property) :: T_CohesiveDryDensity
        real                                     :: Min                  = FillValueReal
        real                                     :: Max                  = FillValueReal
    end type T_CohesiveDryDensity
    
    private :: T_Files
    type       T_Files
        character(Len = PathLength)           :: ConstructData  = null_str
        character(Len = PathLength)           :: Initial    = null_str  
        character(Len = PathLength)           :: OutPutFields   = null_str
        character(Len = PathLength)           :: Final      = null_str
    end type  T_Files
        
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
    end type T_Boxes


    private :: T_Sediment
    type       T_Sediment
        integer                                    :: InstanceID            = FillValueInt
        type (T_Size2D)                            :: Size, WorkSize
        type (T_Time)                              :: BeginTime, EndTime
        type (T_Evolution)                         :: Evolution
        real, dimension(:, :), pointer             :: BatimIncrement, DZ, DZ_Residual => null () 
        real                                       :: MorphologicalFactor
        real                                       :: Density               = FillValueReal
        real                                       :: RelativeDensity       = FillValueReal
        real, dimension(:, :), pointer             :: Dast                  => null ()
        real                                       :: DeltaDensity          = FillValueReal
        integer                                    :: Boundary              = FillValueInt
        real                                       :: TauMax                = FillValueReal
        logical                                    :: TimeSerie             = .false. 
        class (T_Sand), dimension(:), pointer      :: SandClass
        class (T_Cohesive), pointer                :: CohesiveClass
        class (T_CohesiveDryDensity), pointer      :: CohesiveDryDensity
        integer                                    :: NumberOfClasses
        real, dimension(:, :, :), pointer          :: TotalPercentage       => null ()
        real, dimension(:, :, :), pointer          :: Porosity              => null ()
        real, dimension(:, :), pointer             :: D50                   => null () 
        real, dimension(:, :), pointer             :: SandD50               => null () 
        real, dimension(:, :), pointer             :: NDShearStress         => null ()
        real, dimension(:, :), pointer             :: NDShearStressWaves    => null ()
        real, dimension(:, :), pointer             :: NDShearStressMean     => null ()
        real(8), dimension(:,:), pointer           :: FluxX                 => null ()
        real(8), dimension(:,:), pointer           :: FluxY                 => null ()
        real(8), dimension(:,:), pointer           :: Bedload               => null ()
        real(8), dimension(:, :), pointer          :: Mass                  => null ()
        real(8), dimension(:, :), pointer          :: DM                    => null ()
        real(8), dimension(:, :), pointer          :: FluxToSediment        => null ()
        real(8), dimension(:, :), pointer          :: BedloadMass           => null ()
        real(8), dimension(:, :), pointer          :: TotalFluxToSediment   => null ()
        
        real                                       :: PorositySand          = null_real
        real                                       :: ConcRefFactor         = null_real
        real                                       :: ReferenceLevel        = null_real
        
        type(T_Size3D     )                        :: Size3D
        type(T_Size3D     )                        :: WorkSize3D
        type(T_Size3D     )                        :: SedimentSize3D
        type(T_Size3D     )                        :: SedimentWorkSize3D
        
        type (T_Files  )                           :: Files
        type (T_OutPut )                           :: OutPut
        type (T_External)                          :: ExternalVar
        type (T_Discharges)                        :: Discharges
        type (T_Boxes     )                        :: Boxes
        type (T_Residual  )                        :: Residual
        
        real, pointer, dimension(:,:)              :: Elevation             => null()
        real, pointer, dimension(:,:,:)            :: VerticalCoordinate    => null()
        integer, pointer, dimension(:,:)           :: KTop                  => null()
        integer, pointer, dimension (:,:)          :: OpenSediment          => null()
        real                                       :: MinLayerThickness    = null_real
        real                                       :: MaxLayerThickness    = null_real

        integer                                    :: BedloadMethod        = FillValueInt
        integer                                    :: RefConcMethod        = FillValueInt
        logical                                    :: BedSlopeEffects      = .false.
        logical                                    :: WavesOn              = .false.
        
        !Instance of ModuleHDF5        
        integer                                    :: ObjHDF5               = 0
        !Instance of ModuleTimeSerie            
        integer                                    :: ObjTimeSerie          = 0
        !Instance of Module_EnterData           
        integer                                    :: ObjEnterData          = 0
        !Instance of ModuleGridData where the bathymetry is define             
        integer                                    :: ObjBathym             = 0
        !Instance of ModuleGeometry                                             
        !integer                                     :: ObjGeometry          = 0
        !Instance of ModuleHorizontalGrid       
        integer                                    :: ObjHorizontalGrid     = 0
        !Instance of ModuleHorizontalMap        
        integer                                    :: ObjHorizontalMap      = 0 
        !Instance of ModuleMap                  
        integer                                    :: ObjMap                = 0
        !Instance of ModuleTime                 
        integer                                    :: ObjTime               = 0
        !Instance of ModuleDischarges           
        integer                                    :: ObjDischarges         = 0
        !Instance of ModuleBoxDif               
        integer                                    :: ObjBoxDif             = 0             
        !Instance of ModuleWaves
        integer                                    :: ObjWaves              = 0
        !Instance of ModuleGridData                                             
        integer                                    :: ObjSedimentGridData   = 0
        !Instance of ModuleGeometry                                             
        integer                                    :: ObjSedimentGeometry   = 0
        !Instance of ModuleHorizontalMap                                    
        integer                                    :: ObjSedimentHorizontalMap = 0
        !Instance of ModuleMap                                                  
        integer                                    :: ObjSedimentMap        = 0
        
        integer                                    :: ObjFreeVerticalMovement = 0
        
        !List of Sediment Instances
        type(T_Sediment), pointer                      :: Next                  => null ()

    end type  T_Sediment

    !Global Module Variables
    type (T_Sediment), pointer                         :: FirstObjSediment      => null ()
    type (T_Sediment), pointer                         :: Me                    => null () 


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructSediment(ObjSedimentID,                             &
                         ObjGridDataID,                         &
                         ObjHorizontalGridID,                   &
                         ObjHorizontalMapID,                    &
                         ObjTimeID,                             &
                         ObjWavesID,                            &
                         ObjDischargesID,                       &
                         SedimentGridDataID,                    &
                         SedimentHorizontalMapID,               &
                         SedimentMapID,                         &
                         SedimentGeometryID,                    &
                         FreeVerticalMovementID,                &
                         STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSedimentID 
        integer                                         :: ObjGridDataID
        integer                                         :: ObjHorizontalGridID
        integer                                         :: ObjHorizontalMapID
        integer                                         :: ObjTimeID
        integer                                         :: ObjWavesID
        integer                                         :: ObjDischargesID
        integer                                         :: SedimentGridDataID
        integer                                         :: SedimentHorizontalMapID
        integer                                         :: SedimentMapID
        integer                                         :: SedimentGeometryID
        integer                                         :: FreeVerticalMovementID
        integer, optional, intent(OUT)                  :: STAT
           
        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL
        !real                                            :: WaterDensity
        !logical                                         :: WaveTensionON

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: n
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        !Jauch: Defined here because it will be passed as argument in the near future... :P
        !WaveTensionON = .false.

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mSediment_)) then
            nullify (FirstObjSediment)
            call RegisterModule (mSediment_) 
        endif

        call Ready(ObjSedimentID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ObjTime           = AssociateInstance (mTIME_,           ObjTimeID          )
            Me%ObjBathym         = AssociateInstance (mGRIDDATA_,       ObjGridDataID      )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  ObjHorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, ObjHorizontalGridID)
            
            !Sediment Column
            Me%ObjSedimentGridData     = AssociateInstance(mGRIDDATA_,          SedimentGridDataID      )
            Me%ObjSedimentHorizontalMap= AssociateInstance(mHORIZONTALMAP_,     SedimentHorizontalMapID )    
            Me%ObjSedimentGeometry     = AssociateInstance(mGEOMETRY_,          SedimentGeometryID      )
            Me%ObjSedimentMap          = AssociateInstance(mMAP_,               SedimentMapID           )

            if(ObjWavesID /= 0)then
                Me%ObjWaves      = AssociateInstance (mWAVES_,          ObjWavesID         )
            end if

            Me%ExternalVar%WaterDensity  = SigmaDensityReference
            !Me%ExternalVar%WaveTensionON = WaveTensionON

            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                             &
                                       Size        = Me%Size,                            &
                                       WorkSize    = Me%WorkSize,                        &
                                       STAT        = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR10'
            
            
            call GetGeometrySize(Me%ObjSedimentGeometry,                        &
                                 Size       = Me%SedimentSize3D,                &
                                 WorkSize   = Me%SedimentWorkSize3D,            &
                                 STAT       = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructSediment - ModuleSediment - ERR20'

            call ReadLockExternalVar

            call Read_Sediment_Files_Name

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR30'
            
            call ConstructGlobalParameters
            
            call AllocateVariables
            
            call ConstructEvolution
            
            call Construct_Initial_Geometry
            
            call ConstructClasses

            call ConstructOutputTime

            call ConstructD50Cell
          
            call ConstructPorosity
            
            call ConstructClassMass
            
            call ConstructCriticalShearStress
            
            call StartOutputBoxFluxes
            
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR40'


            if (Me%OutPut%Yes) call Open_HDF5_OutPut_File(Me%Files%OutPutFields)

            if (Me%Discharges%Yes) then
            
                if (ObjDischargesID == 0)  then                                                
                    write(*,*)'You need to define a water discharges in the hydrodynamic input' 
                    stop      'ConstructSediment - Sediment - ERR50'
                else
                    Me%ObjDischarges = AssociateInstance (mDISCHARGES_, ObjDischargesID)

                    Me%Discharges%NextCompute = Me%ExternalVar%Now
                endif

            endif
            
            if(FreeVerticalMovementID .ne. 0) then
               do n=1, Me%NumberOfClasses
                    call SetSandParameters(FreeVerticalMovementID, Me%SandClass(n)%ID%IDNumber, Me%SandClass(n)%D50, & 
                                            Me%RelativeDensity, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR60'
                enddo
            endif


            call ReadUnLockExternalVar

            !Returns ID
            ObjSedimentID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleSediment - ConstructSediment - ERR40' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructSediment
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Sediment), pointer                         :: NewObjSediment
        type (T_Sediment), pointer                         :: PreviousObjSediment


        !Allocates new instance
        allocate (NewObjSediment)
        nullify  (NewObjSediment%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjSediment)) then
            FirstObjSediment         => NewObjSediment
            Me                    => NewObjSediment
        else
            PreviousObjSediment      => FirstObjSediment
            Me                    => FirstObjSediment%Next
            do while (associated(Me))
                PreviousObjSediment  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjSediment
            PreviousObjSediment%Next => NewObjSediment
        endif

        Me%InstanceID = RegisterNewInstance (mSediment_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ConstructGlobalParameters

        !External----------------------------------------------------------------
        integer                             :: STAT_CALL

        !Local-------------------------------------------------------------------
        integer                             :: iflag
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB    
        
        !----------------------------------------------------------------------
        
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB

        call GetData(Me%MinLayerThickness,                                               &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MIN_THICKNESS',                                     &
                     default      = 0.01,                                                &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR11' 

        call GetData(Me%MaxLayerThickness,                                          &
                     Me%ObjEnterData, iflag,                                        &
                     SearchType   = FromFile,                                       &
                     keyword      = 'MAX_THICKNESS',                                &
                     Default      = 0.5,                                            &
                     ClientModule = 'ModuleSediment',                               &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Construct_Parameters - ModuleConsolidation - ERR31'


        call GetData(Me%PorositySand,                                                    &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'POROSITY_SAND',                                     &
                     default      = 0.3,                                                 &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR12' 

        call GetData(Me%Density,                                                         &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DENSITY_SEDIMENT',                                  &
                     default      = 2650.,                                               &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)               
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR13' 
        
        Me%RelativeDensity = Me%Density / Me%ExternalVar%WaterDensity

        Me%DeltaDensity    = Me%Density - Me%ExternalVar%WaterDensity
        
        call GetData(Me%BedloadMethod,                                                   &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BEDLOAD_METHOD',                                    &
                     default      = 1,                                                   &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR14' 

        !NoTransport = 0, CurrentOnly = 1, CurrentPlusWaves = 2, WavesOnly = 3
       
        if     ( Me%BedloadMethod == 2 &
            .OR. Me%BedloadMethod == 3 )  then
        
            if (Me%ObjWaves == 0) stop 'ConstructGlobalParameters - ModuleSediment - ERR16' 

            !if (.not. Me%ExternalVar%WaveTensionON) stop 'ConstructGlobalParameters - ModuleSediment - ERR20' 
            
            Me%WavesOn = .true. 

        endif
            
         call GetData(Me%BedSlopeEffects,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BEDSLOPE',                                          &
                     default      = .TRUE.,                                              &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR25' 
            

        call GetData(Me%Evolution%Bathym,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BATHYMETRY_EVOLUTION',                              &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR26' 


        call GetData(Me%Boundary,                                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'BOUNDARY',                                          &
                     ClientModule = 'ModuleSediment',                                        &
                     Default      = NullGradient,                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR59'

        call GetData(Me%TauMax,                                                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TAU_MAX',                                           &
                     default      = 10.,                                                 &
                     ClientModule = 'ModuleSediment',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR79' 


        call GetData(Me%Discharges%Yes,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DISCHARGES',                                        &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSediment',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR89' 
        
        call GetData(Me%Evolution%OLD,                                                   &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'OLD',                                               &
                     default      = .false.,                                             &
                     ClientModule = 'ModuleSediment',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSediment - ERR130' 

        call GetData(Me%MorphologicalFactor,                                             &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MORPHOLOGICAL_FACTOR',                              &
                     default      = 1.,                                                  &
                     ClientModule = 'ModuleSediment',                                        &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSediment - ERR140' 


        call GetData(Me%Residual%ON,                                                    &
                     Me%ObjEnterData,iflag,                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'RESIDUAL',                                         &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleSediment',                                       &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSediment - ERR150' 
        
        call GetData(Me%ConcRefFactor,                                                   &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'CONCREF_FACTOR',                                    &
                     default      = 1.,                                                  &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR160' 
        
        call GetData(Me%ReferenceLevel,                                                  &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'REFERENCE_LEVEL',                                   &
                     default      = 0.01,                                                &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR170'
        
        call GetData(Me%RefConcMethod,                                                   &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'REFCONC_METHOD',                                    &
                     default      = 2,                                                   &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR180'

        !Zyserman and Fredsoe (1994)  = 1, van Rijn (2007) = 2

    end subroutine ConstructGlobalParameters

    !----------------------------------------------------------------------------
    
    subroutine AllocateVariables
    
        !Local-----------------------------------------------------------------
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------
 
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB
         
        allocate(Me%D50(ILB:IUB, JLB:JUB))
        Me%D50(:,:) = FillValueReal
        
        allocate(Me%SandD50(ILB:IUB, JLB:JUB))
        Me%SandD50(:,:) = FillValueReal
            
        allocate(Me%NDShearStress (ILB:IUB, JLB:JUB))
        Me%NDShearStress(:,:) = FillValueReal
        
        if (Me%WavesOn) then
            allocate(Me%NDShearStressWaves (ILB:IUB, JLB:JUB))
            Me%NDShearStressWaves(:,:) = FillValueReal
            
            allocate(Me%NDShearStressMean (ILB:IUB, JLB:JUB))
            Me%NDShearStressMean(:,:) = FillValueReal
        endif
        
        allocate(Me%Dast (ILB:IUB, JLB:JUB))
        Me%Dast(:,:) = FillValueReal
                 
        allocate(Me%Elevation(ILB:IUB, JLB:JUB)) 
        Me%Elevation(:,:) = 0.
        
        allocate(Me%OpenSediment(ILB:IUB, JLB:JUB)) 
        Me%OpenSediment(:,:) = 0.
        
        allocate(Me%Porosity(ILB:IUB, JLB:JUB, KLB:KUB)) 
        Me%Porosity(:,:,:) = null_real
            
        allocate(Me%VerticalCoordinate(ILB:IUB, JLB:JUB, KLB:KUB)) 
        Me%VerticalCoordinate(:,:,:) = null_real
        
        allocate(Me%KTop(ILB:IUB, JLB:JUB)) 
        Me%KTop(:,:) = 0.
            
        allocate(Me%FluxX(ILB:IUB, JLB:JUB))
        Me%FluxX(:,:) = null_real

        allocate(Me%FluxY(ILB:IUB, JLB:JUB))
            Me%FluxY(:,:) = null_real

        allocate(Me%Bedload(ILB:IUB, JLB:JUB))
        Me%Bedload(:,:) = null_real
            
         allocate(Me%DZ(ILB:IUB, JLB:JUB))
        Me%DZ(:,:) = 0.
        
        allocate(Me%DM(ILB:IUB, JLB:JUB))
        Me%DM(:,:) = 0.
        
        allocate(Me%Mass(ILB:IUB, JLB:JUB))
        Me%Mass(:,:) = 0.
        
        allocate(Me%FluxToSediment(ILB:IUB, JLB:JUB))
        Me%FluxToSediment(:,:) = 0.
        
        allocate(Me%BedloadMass(ILB:IUB, JLB:JUB))
        Me%BedloadMass(:,:) = 0.
        
        allocate(Me%TotalFluxToSediment(ILB:IUB, JLB:JUB))
        Me%TotalFluxToSediment(:,:) = 0.
        
        if (Me%Evolution%Bathym) then
            allocate(Me%BatimIncrement(ILB:IUB, JLB:JUB))
            Me%BatimIncrement(:,:) = 0.
        endif

        allocate(Me%DZ_Residual(ILB:IUB, JLB:JUB))
        if (Me%Evolution%OLD) then
            call ReadInitialField(FieldName = "DZ_Residual", Field2D =Me%DZ_Residual)
        else
           Me%DZ_Residual(:,:) = 0.
        endif
        
        if (Me%Residual%ON) then
            allocate(Me%Residual%FluxX(ILB:IUB, JLB:JUB))
            Me%Residual%FluxX(:,:) = null_real

            allocate(Me%Residual%FluxY(ILB:IUB, JLB:JUB))
            Me%Residual%FluxY(:,:) = null_real
            
            if (Me%Evolution%Old) then
            
                call ReadResidualStartTime()
                
                call ReadInitialField(FieldName = "Bedload X", Field2D = Me%Residual%FluxX)
                
                call ReadInitialField(FieldName = "Bedload Y", Field2D = Me%Residual%FluxY)                
            
            endif
            
        endif
            
    end subroutine AllocateVariables
    
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
                     ClientModule = 'ModuleSediment',                                       &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'StartOutputBoxFluxes - ModuleSediment - ERR10'

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
                         ClientModule = 'ModuleSediment',                                   &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_)                                                &
                stop 'StartOutputBoxFluxes - ModuleSediment - ERR20'

            if (iflag .EQ. 0) then
                stop 'StartOutputBoxFluxes - ModuleSediment - ERR30'    
            endif
        
            inquire(File = Me%Boxes%File, Exist = exist)
            if (exist) then
                inquire(File = Me%Boxes%File, Opened  = Opened)
                if (opened) then
                    write(*,*    ) 
                    write(*,'(A)') 'BoxesFile = ',trim(adjustl(Me%Boxes%File))
                    write(*,*    ) 'Already opened.'
                    stop           'StartOutputBoxFluxes - ModuleSediment - ERR40'    
                end if
            else
                write(*,*) 
                write(*,*)     'Could not find the boxes file.'
                write(*,'(A)') 'BoxFileName = ', Me%Boxes%File
                stop           'StartOutputBoxFluxes - ModuleSediment - ERR50'    
            end if

            allocate(ScalarOutputList(1))
            allocate(FluxesOutputList(1))

            ScalarOutputList(1) = 'Sediment'
            FluxesOutputList(1) = 'Sediment'

            call StartBoxDif(BoxDifID           = Me%ObjBoxDif,                         &
                             TimeID             = Me%ObjTime,                           &
                             HorizontalGridID   = Me%ObjHorizontalGrid,                 &
                             BoxesFilePath      = Me%Boxes%File,                        &
                             FluxesOutputList   = FluxesOutputList,                     &
                             ScalarOutputList   = ScalarOutputList,                     &
                             WaterPoints2D      = Me%ExternalVar%WaterPoints2D,         &
                             STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'StartOutputBoxFluxes - ModuleSediment - ERR60'

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
            stop 'ConstructEvolution - ModuleSediment - ERR10'


        call GetData(Me%Evolution%SedimentDT,                                            &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'SEDIMENT_DT',                                       &
                     default      = ModelDT,                                             &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSediment - ERR20' 

        call GetComputeTimeLimits(Me%ObjTime,                                            &
                                  EndTime   = Me%EndTime,                                &
                                  BeginTime = Me%BeginTime,                              &
                                  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ConstructEvolution - ModuleSediment - ERR30'

        if (Me%Evolution%SedimentDT < (ModelDT)) then

            !Sediment DT  must be a submultiple of the ModelDT
            auxFactor = ModelDT / Me%Evolution%SedimentDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop 'ConstructEvolution - ModuleSediment - ERR40'
            endif

        elseif (Me%Evolution%SedimentDT > (ModelDT)) then

            !Sediment DT  must be a multiple of the ModelDT
            auxFactor = Me%Evolution%SedimentDT  / ModelDT

            Erroraux = auxFactor - int(auxFactor)
            if (Erroraux /= 0) then
                write(*,*) 
                write(*,*) ' Time step error.'
                stop 'ConstructEvolution - ModuleSediment - ERR50'
            endif

        endif

        ! Run period in seconds
        DTaux = Me%EndTime - Me%BeginTime

        !The run period   must be a multiple of the Sediment DT
        auxFactor = DTaux / Me%Evolution%SedimentDT

        ErrorAux = auxFactor - int(auxFactor)
        if (ErrorAux /= 0) then
            write(*,*) 
            write(*,*) ' Time step error.'
            stop 'ConstructEvolution - ModuleSediment - ERR60'
        endif


        Me%Evolution%NextSediment = Me%BeginTime + Me%Evolution%SedimentDT

        if (Me%Evolution%Bathym) then
            
            call GetData(Me%Evolution%BathymDT,                                              &
                         Me%ObjEnterData,iflag,                                              &
                         SearchType   = FromFile,                                            &
                         keyword      = 'EVOLUTION_DT',                                          &
                         default      = Me%Evolution%SedimentDT,                                 &
                         ClientModule = 'ModuleSediment',                                        &
                         STAT         = STAT_CALL)              
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructEvolution - ModuleSediment - ERR70'
        
        
            if (Me%Evolution%BathymDT < Me%Evolution%SedimentDT) then
        
                stop 'ConstructEvolution - ModuleSediment - ERR90'
        
            elseif (Me%Evolution%BathymDT > Me%Evolution%SedimentDT) then
        
                !Batim DT  must be a multiple of the Sediment DT
                auxFactor = Me%Evolution%BathymDT / Me%Evolution%SedimentDT
        
                Erroraux = auxFactor - int(auxFactor)
                if (Erroraux /= 0) then
                    write(*,*) 
                    write(*,*) ' Time step error.'
                    stop 'ConstructEvolution - ModuleSediment - ERR100'
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
                stop 'ConstructEvolution - ModuleSediment - ERR120'
            endif
        
            Me%Evolution%NextBatim = Me%BeginTime + Me%Evolution%BathymDT
            
        endif
        
        if (Me%Residual%ON .and. .not. Me%Evolution%OLD) then
            Me%Residual%StartTime = Me%BeginTime
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
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%SedimentWorkSize3D%ILB 
        WorkIUB = Me%SedimentWorkSize3D%IUB 
        WorkJLB = Me%SedimentWorkSize3D%JLB 
        WorkJUB = Me%SedimentWorkSize3D%JUB 
        WorkKLB = Me%SedimentWorkSize3D%KLB 
        WorkKUB = Me%SedimentWorkSize3D%KUB 


        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5      (Me%ObjHDF5,                                &
                                 trim(FileName)//"5",                       &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR010'
        
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,         &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR40'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR50'
        
        
        if (Me%Evolution%Bathym) then

            call GetGridData2Dreference(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR60' 

        else

            call GetGridData           (Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR70' 
        
        endif

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = Me%ExternalVar%InitialBathym,                    &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR80'

        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR90' 


        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                 &
                              Array2D = Me%ExternalVar%WaterPoints2D,                    &
                              STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR100'
        
        !Write WaterPoints3D
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                  &
                            WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR110'
                
                
        call HDF5WriteData  (Me%ObjHDF5, "/Grid", "WaterPoints3D",         &
                            "-", Array3D = Me%ExternalVar%WaterPoints3D,             &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR120'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR130'

        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File
   
    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !Read the name of the files need to construct and modify
    !the Sediment properties 

    subroutine Read_Sediment_Files_Name

        !External--------------------------------------------------------------
        integer                       :: STAT_CALL 
        character(len = StringLength) :: Message

        !----------------------------------------------------------------------

        ! ---> ASCII file used to construct new properties
        Message   ='ASCII file used to construct Sediment instance.'
        Message   = trim(Message)

        call ReadFileName('SEDIMENT_DAT', Me%Files%ConstructData,                           &
                           Message = Message, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sediment_Files_Name - ModuleSediment - ERR10' 


        ! ---> File in HDF format where is written instant fields of Sediment properties
        Message   ='Instant fields of Sediment properties in HDF format.'
        Message   = trim(Message)
        
        if (GetDDecompON    (Me%ObjHorizontalGrid)) then
            write(*,*) 'Module Sediment not ready to run in domain decomposition mode'
            stop 'Read_Sediment_Files_Name - ModuleSediment - ERR20' 
        endif

        call ReadFileName('SEDIMENT_HDF', Me%Files%OutPutFields,                             &
                           Message = Message, TIME_END = Me%EndTime,                     &
                           Extension = 'Sedt',                                         &
                           MPI_ID    = GetDDecompMPI_ID(Me%ObjHorizontalGrid),&
                           DD_ON     = GetDDecompON    (Me%ObjHorizontalGrid),&
                           STAT      = STAT_CALL)

        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sediment_Files_Name - ModuleSediment - ERR02' 

        ! ---> Sediment properties final values in HDF format
        Message   ='Sediment properties final values in HDF format.'
        Message   = trim(Message)
        call ReadFileName('SEDIMENT_FIN', Me%Files%Final,                                &
                           Message = Message, TIME_END = Me%EndTime,                     &
                           Extension = 'Sedf',                                         &
                           MPI_ID    = GetDDecompMPI_ID(Me%ObjHorizontalGrid),&
                           DD_ON     = GetDDecompON    (Me%ObjHorizontalGrid),&
                           STAT      = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            stop 'Read_Sediment_Files_Name - ModuleSediment - ERR03' 


        ! ---> Sediment properties initial values in HDF format
        Message   ='Sediment properties initial values in HDF format.'
        Message   = trim(Message)

        call ReadFileName('SEDIMENT_INI', Me%Files%Initial,                              &
                           Message   = Message, TIME_END = Me%ExternalVar%Now,           &
                           Extension = 'sedi',                                          &
                           MPI_ID    = GetDDecompMPI_ID(Me%ObjHorizontalGrid),&
                           DD_ON     = GetDDecompON    (Me%ObjHorizontalGrid),&
                           STAT      = STAT_CALL)


cd1 :   if      (STAT_CALL .EQ. FILE_NOT_FOUND_ERR_   ) then
            write(*,*)  
            write(*,*) 'Inicial file not found.'
            stop 'Read_Sediment_Files_Name - ModuleSediment - ERR04' 

        else if (STAT_CALL .EQ. KEYWORD_NOT_FOUND_ERR_) then
            write(*,*)  
            write(*,*) 'Keyword for the inicial file not found in nomfich.dat.'
            write(*,*) 'Read_Sediment_Files_Name - ModuleSediment - WRN01'
            write(*,*)  

        else if (STAT_CALL .EQ. SUCCESS_              ) then
            continue
        else
            write(*,*) 
            write(*,*) 'Error calling ReadFileName.'
            stop 'Read_Sediment_Files_Name - ModuleSediment - ERR05' 
        end if cd1  

        !----------------------------------------------------------------------

    end subroutine Read_Sediment_Files_Name

   !--------------------------------------------------------------------------

    subroutine ConstructTimeSerie

        !External--------------------------------------------------------------
        character(len=StringLength)                         :: TimeSerieLocationFile
        integer                                             :: STAT_CALL, iflag

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber 
        character(len=StringLength)                         :: TimeSerieName
        character(len=StringLength), dimension(:), pointer  :: PropertyList
        type (T_Polygon), pointer                           :: ModelDomainLimit
        integer                                             :: FixedOutputs, NumberOfOutputs, n


        !----------------------------------------------------------------------

        !First checks out how many properties will have time series

        !Allocates PropertyList
        if (Me%CohesiveClass%Run) then
            FixedOutputs = 14
        else
            FixedOutputs = 11
        endif
        
        NumberOfOutputs = FixedOutputs + Me%NumberOfClasses
        
        allocate(PropertyList(NumberOfOutputs), STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSediment - ERR10'

        PropertyList(1) = 'DZ(m)'
        PropertyList(2) = 'DZ_Residual(m)'
        PropertyList(3) = 'D50(m)'
        PropertyList(4) = 'Bedload(kg/s/m)'
        PropertyList(5) = 'Bedload_X(kg/s)'
        PropertyList(6) = 'Bedload_Y(kg/s)'
        PropertyList(7) = 'Flux_Sand(kg/s/m2)'
        PropertyList(8) = 'DM(kg)'
        PropertyList(9) ='TotalBedload_Mass(kg)'
        PropertyList(10) ='TotalFlux_Sand(kg)'
        PropertyList(11) ='ShearStress(N/m2)'
        
        if (Me%CohesiveClass%Run) then   
            PropertyList(12) = 'Flux_Cohesive(kg/s/m2)'
            PropertyList(13) ='TotalFlux_Cohesive(kg)'
            PropertyList(14) ='CriticalShear_Cohesive(N/m2)'
        endif
        
do1:    do n=1,Me%NumberOfClasses    
            PropertyList(FixedOutputs + n) = 'Ref_Concentration(kg/m3)_'//trim(adjustl(Me%SandClass(n)%ID%name))
        enddo do1
 
        call GetData(TimeSerieLocationFile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TIME_SERIE_LOCATION',                              &
                     ClientModule = 'ModuleSediment',                                       &
                     Default      = Me%Files%ConstructData,                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                    &
            stop 'ConstructTimeSerie - ModuleSediment - ERR20' 

        call GetGridOutBorderPolygon(HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     Polygon          = ModelDomainLimit,               &
                                     STAT             = STAT_CALL)           
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'Construct_Time_Serie - ModuleSediment - ERR25' 
        
        !Constructs TimeSerie
        call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                                &
                            trim(TimeSerieLocationFile),                                &
                            PropertyList, "srsed",                                     &
                            WaterPoints2D = Me%ExternalVar%WaterPoints2D,               &
                            ModelDomain   = ModelDomainLimit,                           & 
                            STAT          = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSediment - ERR30'
        
        call UngetHorizontalGrid(HorizontalGridID = Me%ObjHorizontalGrid,               &
                                 Polygon          = ModelDomainLimit,                   &
                                 STAT             = STAT_CALL)                          
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR35'

        !Deallocates PropertyList
        deallocate(PropertyList, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSediment - ERR40'

        !Corrects if necessary the cell of the time serie based in the time serie coordinates
        call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR50'

        do dn = 1, TimeSerieNumber
        
            call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR60'
            
            if (IgnoreOK) cycle        

            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      CoordX   = CoordX,                                &
                                      CoordY   = CoordY,                                & 
                                      CoordON  = CoordON,                               &
                                      STAT     = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR70'
            
            call GetTimeSerieName(Me%ObjTimeSerie, dn, TimeSerieName, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR80'  
                                  
            if (CoordON) then
                call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                
                if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= OUT_OF_BOUNDS_ERR_) then
                    stop 'ConstructTimeSerie - ModuleSediment - ERR90'
                endif                            

                call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR100'
            endif
            
            call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                      LocalizationI   = Id,                             &
                                      LocalizationJ   = Jd,                             & 
                                      STAT     = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ConstructTimeSerie - ModuleSediment - ERR120'

            if (Me%ExternalVar%WaterPoints2D(Id, Jd) /= WaterPoint) then
                
                 write(*,*) 'Time Serie in a land cell - ',trim(TimeSerieName),' - ',' ModuleSediment'

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
            stop 'ConstructOutputTime - ModuleSediment - ERR01' 

        Me%OutPut%NextOutPut = 1

        call GetData(Me%TimeSerie,                                                       &
                     Me%ObjEnterData, iflag,                                             &
                     Keyword        = 'TIME_SERIE',                                      &
                     Default        = .false.,                                           &
                     SearchType     = FromFile,                                          &
                     ClientModule   = 'ModuleInterfaceSedimentWater',                    &
                     STAT           = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'ConstructOutputTime - ModuleSediment - ERR02' 

        if (Me%TimeSerie) call ConstructTimeSerie

    end subroutine ConstructOutputTime

    !--------------------------------------------------------------------------

    subroutine ConstructClasses

        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB, iflag
        integer                         :: n
        !----------------------------------------------------------------------
 
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB


        call GetNumberOfBlocks (Me%ObjEnterData, class_block_begin, class_block_end,    &
                                FromFile,                                               &
                                Me%NumberOfClasses,                                     &
                                STAT = STAT_CALL)
                                
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR01'  
        
        allocate(Me%SandClass(Me%NumberOfClasses))
            
do1 :   do n=1, Me%NumberOfClasses
            call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                                    ClientNumber    = ClientNumber,             &
                                    block_begin     = class_block_begin,         &
                                    block_end       = class_block_end,           &
                                    BlockFound      = BlockFound,                &
                                    STAT            = STAT_CALL)
cd1 :       if (STAT_CALL .EQ. SUCCESS_) then    
cd2 :           if (BlockFound) then
    
                    call ConstructPropertyID(Me%SandClass(n)%ID, Me%ObjEnterData, FromBlock, CheckProperty = .false.)
                        Me%SandClass(n)%ID%IDNumber = RegisterDynamicProperty (Me%SandClass(n)%ID%Name)
          
                    call GetData(Me%SandClass(n)%D50,                                     &
                                    Me%ObjEnterData,iflag,                                 &
                                    SearchType   = FromBlock,                              &
                                    keyword      = 'D50',                                  &
                                    ClientModule = 'ModuleSediment',                           &
                                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR10' 
                        
                    if (iflag .NE. 1) stop 'ConstructClasses - ModuleSediment - ERR20' 
                    
                    if(Me%SandClass(n)%D50 < 65.0E-6) then
                        write(*,*) 
                        write(*,*) 'Sand classes must have diameters greater than 65.0E-6' 
                        write(*,*) 'Check property', Me%SandClass(n)%ID%Name
                        stop       'ConstructClasses - ModuleSediment - ERR30'
                    endif
                        
                    call GetData(Me%SandClass(n)%parameter_A,                             &
                                    Me%ObjEnterData,iflag,                                 &
                                    SearchType   = FromBlock,                              &
                                    keyword      = 'A',                                    &
                                    default      = 12.,                                     &
                                    ClientModule = 'ModuleSediment',                       &
                                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR040' 
                         
                    call GetData(Me%SandClass(n)%parameter_n,                             &
                                    Me%ObjEnterData,iflag,                                 &
                                    SearchType   = FromBlock,                              &
                                    keyword      = 'n',                                    &
                                    default      = 0.5,                                    &
                                    ClientModule = 'ModuleSediment',                       &
                                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR50'        
                        
                    call GetData(Me%SandClass(n)%parameter_p,                             &
                                    Me%ObjEnterData,iflag,                                 &
                                    SearchType   = FromBlock,                              &
                                    keyword      = 'p',                                    &
                                    default      = 1.,                                     &
                                    ClientModule = 'ModuleSediment',                       &
                                    STAT         = STAT_CALL)

                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR60'
                                      
                    call ConstructPropertyValue (Me%SandClass(n), FromBlock)
                        
                    !Allocate variables of classes 
                    allocate(Me%SandClass(n)%FluxX(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%FluxX(:,:) = null_real

                    allocate(Me%SandClass(n)%FluxY(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%FluxY(:,:) = null_real

                    allocate(Me%SandClass(n)%Bedload(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%Bedload(:,:) = null_real

                    allocate(Me%SandClass(n)%CriticalShearStress(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%CriticalShearStress(:,:) = FillValueReal
                        
                    allocate(Me%SandClass(n)%NDCriticalShearStress(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%CriticalShearStress(:,:) = FillValueReal
                        
                    allocate(Me%SandClass(n)%TopPercentage(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%TopPercentage(:,:) = null_real
                        
                    allocate(Me%SandClass(n)%HidingFactor(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%HidingFactor(:,:) = null_real
                        
                    allocate(Me%SandClass(n)%DM(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%DM(:,:) = 0. 
                    
                    allocate(Me%SandClass(n)%Mass(ILB:IUB, JLB:JUB, KLB:KUB))
                    Me%SandClass(n)%Mass(:,:,:) = 0.  
                    
                    allocate(Me%SandClass(n)%ConcRef(ILB:IUB, JLB:JUB))
                    Me%SandClass(n)%ConcRef(:,:) = 0.  
                    
                    if (Me%WavesOn) then
                        allocate(Me%SandClass(n)%BedloadX(ILB:IUB, JLB:JUB))
                        Me%SandClass(n)%BedloadX(:,:) = null_real
                        
                        allocate(Me%SandClass(n)%BedloadY(ILB:IUB, JLB:JUB))
                        Me%SandClass(n)%BedloadY(:,:) = null_real
                    endif
                    
                else cd2

                    write(*,*)  
                    write(*,*) 'Error calling ExtractBlockFromBlock. '
                    stop       'ConstructClasses - ModuleSediment - ERR80'

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'ConstructClasses - ModuleSediment - ERR90'
            else cd1
                stop       'ConstructClasses - ModuleSediment - ERR100'
            end if cd1
        end do do1

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR110'

        !call GetNumberOfBlocks (Me%ObjEnterData, cohesive_block_begin,          &
        !            cohesive_block_end,                                 &
        !            FromFile,                                           &
        !            NumberOfCohesive,                                   &
        !            STAT = STAT_CALL)
        !if (STAT_CALL .EQ. SUCCESS_) then            
        !
        !    if (NumberOfCohesive > 1) then
        !        write(*,*) 'Only one cohesive sediment class allowed. '
        !        stop       'ConstructClasses - ModuleSediment - ERR07'
        !    endif
        !endif
                
        call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                        ClientNumber    = ClientNumber,                          &
                        block_begin     = cohesive_block_begin,                  &
                        block_end       = cohesive_block_end,                    &
                        BlockFound      = BlockFound,                            &
                        STAT            = STAT_CALL)
        if (STAT_CALL .EQ. SUCCESS_) then
            
            allocate(Me%CohesiveClass)
            
            if (BlockFound) then

                Me%CohesiveClass%Run = .true.
                
                call ConstructPropertyID(Me%CohesiveClass%ID, Me%ObjEnterData, FromBlock)
   
                call GetData(Me%CohesiveClass%D50,                               &
                                Me%ObjEnterData,iflag,                           &
                                SearchType   = FromBlock,                        &
                                keyword      = 'D50',                            &
                                ClientModule = 'ModuleSediment',                 &
                                STAT         = STAT_CALL)
                        
                if (iflag .NE. 1) stop 'ConstructClasses - ModuleSediment - ERR120'
                
                !Critical erosion shear stress for weak consolidated cohesive sediment 
                 call GetData(Me%CohesiveClass%WeakConsolidated_CSS,             &
                                Me%ObjEnterData,iflag,                           &
                                SearchType   = FromBlock,                        &
                                keyword      = 'WEAK_CONSOLIDATED_CRITICAL',     &
                                default      = 0.5,                              &
                                ClientModule = 'ModuleSediment',                 &
                                STAT         = STAT_CALL)
                 
                 !Maximum mud percentage in which the bed starts to have a cohesive behaviour
                 call GetData(Me%CohesiveClass%PM1_MAX,                              &
                                Me%ObjEnterData,iflag,                           &
                                SearchType   = FromBlock,                        &
                                keyword      = 'PM1_MAX',                            &
                                default      = 0.3,                              &
                                ClientModule = 'ModuleSediment',                 &
                                STAT         = STAT_CALL)
                 
                !Mud percentage for a fully cohesive bed 
                call GetData(Me%CohesiveClass%PM2,                              &
                                Me%ObjEnterData,iflag,                           &
                                SearchType   = FromBlock,                        &
                                keyword      = 'PM2',                            &
                                default      = 0.6,                              &
                                ClientModule = 'ModuleSediment',                 &
                                STAT         = STAT_CALL)
                
                
                allocate(Me%CohesiveClass%CriticalShearStress(ILB:IUB, JLB:JUB))
                Me%CohesiveClass%CriticalShearStress(:,:) = FillValueReal
                
                allocate(Me%CohesiveClass%PM1(ILB:IUB, JLB:JUB))
                Me%CohesiveClass%PM1(:,:) = FillValueReal
                        
                allocate(Me%CohesiveClass%TopPercentage(ILB:IUB, JLB:JUB))
                Me%CohesiveClass%TopPercentage(:,:) = null_real
                        
                allocate(Me%CohesiveClass%DM(ILB:IUB, JLB:JUB))
                Me%CohesiveClass%DM(:,:) = 0. 
            
                allocate(Me%CohesiveClass%Mass(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%CohesiveClass%Mass(:,:,:) = 0. 
            
                allocate(Me%CohesiveClass%Porosity(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%CohesiveClass%Porosity(:,:,:) = null_real 
            
                call ConstructPropertyValue (Me%CohesiveClass, FromBlock)
            
                call ConstructCohesiveDryDensity
                
            endif
                    
        endif
            
        call Block_UnLock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR130'


        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR140'
        
        call CheckPercentageConsistence
               
    
    end subroutine ConstructClasses
    
    !-------------------------------------------------------------------------- 
    
    subroutine CheckPercentageConsistence

        !Local-------------------------------------------------------------------
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                         :: i, j, k, n
        real                            :: Percentage
        !----------------------------------------------------------------------
 
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB
    
        
        allocate(Me%TotalPercentage(ILB:IUB, JLB:JUB, KLB:KUB))
        Me%TotalPercentage(:,:,:)=0
        
        if (Me%CohesiveClass%Run) then
            
            do j=JLB, JUB
            do i=ILB, IUB
            do k=KLB, Me%KTop(i, j)
                    
                if (Me%ExternalVar%WaterPoints3D (i,j,k) == WaterPoint) then
                
                    Me%TotalPercentage (i, j, k) = Me%CohesiveClass%Field3D (i,j,k)
                                       
                endif
            enddo
            enddo
            enddo
        endif    
    
        do j=JLB, JUB
        do i=ILB, IUB
        do k=KLB, Me%KTop(i, j)
                    
            if (Me%ExternalVar%WaterPoints3D (i,j,k) == WaterPoint) then
                
                do n=1,Me%NumberOfClasses
                    
                    Percentage = Me%SandClass(n)%Field3D (i,j,k)
                        
                    Me%TotalPercentage (i, j, k) = Me%TotalPercentage (i, j, k) + Percentage  
                    
                enddo                
            endif
        enddo
        enddo
        enddo
                                    
        do j=JLB, JUB
        do i=ILB, IUB
        do k=KLB, Me%KTop(i, j)
                    
            if (Me%ExternalVar%WaterPoints3D (i,j,k) == WaterPoint) then
                if (Me%TotalPercentage (i, j, k)  > 1.001) then
                    write(*,*) 'The sum of the classes percentage is larger than 100%.'
                    write(*,*) i, j, k, Me%TotalPercentage (i, j, k)
                    stop 'ConstructClasses - ModuleSediment - ERR11'
                elseif (Me%TotalPercentage (i, j, k)  < 0.999) then
                        write(*,*) 'The sum of the classes percentage is smaller than 100%.'
                        write(*,*) i, j, k, Me%TotalPercentage (i, j, k)
                        stop 'ConstructClasses - ModuleSediment - ERR12'
                endif
            endif
        enddo
        enddo
        enddo
        
    end subroutine CheckPercentageConsistence
    
    !--------------------------------------------------------------------------

    subroutine ConstructCohesiveDryDensity

        !External----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockInBlockFound

        !Local-------------------------------------------------------------------
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB, iflag
        !----------------------------------------------------------------------
        
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB
        
        allocate(Me%CohesiveDryDensity)
        
        call GetData(Me%CohesiveDryDensity%Min,                                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MIN_DENSITY_COHESIVE',                              &
                     default      = 200.,                                                &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR12' 
        
        call GetData(Me%CohesiveDryDensity%MAX,                                          &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MAX_DENSITY_COHESIVE',                              &
                     default      = 1000.,                                               &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR12' 
        
        if (Me%Evolution%Old) then
            
            allocate(Me%CohesiveDryDensity%Field3D(ILB:IUB, JLB:JUB, KLB:KUB))
                            
            call ReadInitialField3D(FieldName = "Cohesive_DryDensity", Field3D = Me%CohesiveDryDensity%Field3D)  
            
        else
        
            call ExtractBlockFromBlock (Me%ObjEnterData,                              &
                                    ClientNumber      = ClientNumber,                 &
                                    block_begin       = cohesivedensity_block_begin,  &
                                    block_end         = cohesivedensity_block_end,    &
                                    BlockInBlockFound = BlockInBlockFound,            &
                                    STAT            = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then
    
                if (BlockInBlockFound) then 
          
                    call ConstructPropertyValue (Me%CohesiveDryDensity, FromBlockInBlock) 
                
                endif
        
            else
    
                allocate(Me%CohesiveDryDensity%Field3D(ILB:IUB, JLB:JUB, KLB:KUB))
 
                Me%CohesiveDryDensity%Field3D(:,:,:) = Me%CohesiveDryDensity%Min
                
                write(*,*) 'Block with initial dry density of cohesive sediment was not found'                
                write(*,*) 'Minimum value assumed', Me%CohesiveDryDensity%Min, 'kg/m3'  
        
            endif
            
        endif
        
        call CheckFieldConsistence (Me%CohesiveDryDensity)
    
    end subroutine ConstructCohesiveDryDensity

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    
    subroutine ConstructD50Cell
       !Local-----------------------------------------------------------------
        integer             :: i, j, n, WKUB
        integer             :: WILB, WIUB, WJLB, WJUB
        class(T_Sand), pointer :: SandClass
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%D50(:,:)= 0      
        Me%SandD50(:,:)= 0
        
        do j=WJLB, WJUB
        do i=WILB, WIUB
                    
            WKUB = Me%KTop(i, j)

            if (Me%ExternalVar%WaterPoints3D (i,j,WKUB) == WaterPoint) then
                        
                do n=1,Me%NumberOfClasses

                    SandClass => Me%SandClass(n)
                        
                    Me%SandD50(i,j) = SandClass%D50 * SandClass%Field3D(i, j, WKUB) + Me%SandD50(i,j)
                    
                enddo
            
                Me%D50(i,j) = Me%SandD50(i,j)
            
                if (Me%CohesiveClass%Run) then
                    
                    Me%D50(i,j) = Me%CohesiveClass%D50 * Me%CohesiveClass%Field3D(i, j, WKUB) + Me%D50(i,j)
                    
                endif
                    
                
            endif
        enddo
        enddo 

    end subroutine ConstructD50Cell
    !--------------------------------------------------------------------------

    subroutine ConstructPorosity
        
        !Local-----------------------------------------------------------------
        integer                     :: i, j, k
        real                        :: ymin, y, c
        integer                     :: ILB, IUB, JLB, JUB, KLB, KUB    
        !----------------------------------------------------------------------      
        
        ymin = 0.76
        
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB
            
        do j=JLB, JUB
        do i=ILB, IUB
        do k=KLB, KUB 
                
            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint .and.   &
                .not. Me%ExternalVar%VolumeZ(i, j, k) == 0.) then
                         
                if (Me%CohesiveClass%Run) then
                     
                    Me%CohesiveClass%Porosity(i,j,k) =  1 - Me%CohesiveDryDensity%Field3D(i,j,k) / Me%Density 
                     
                    if(Me%SandD50(i,j) > 0)then
        
                        !Volume fraction of the cohesive class based on total volume        
                        c = Me%CohesiveClass%Field3D(i,j,k) * (1 - Me%PorositySand)       
                               
                
                        if (c < Me%PorositySand) then
                    
                            y = c * (ymin - 1) / Me%PorositySand + 1
                
                            Me%Porosity(i,j,k) = Me%PorositySand - c * y * (1 - Me%CohesiveClass%Porosity(i,j,k)) +     & 
                                                        (1 - y) * c * Me%CohesiveClass%Porosity(i,j,k)
                
                        else
                    
                            y = (c - 1) * (1 - ymin)/(1 - Me%PorositySand) + 1
                    
                            Me%Porosity(i,j,k) = Me%PorositySand * (1 - y) + c * Me%CohesiveClass%Porosity(i,j,k)
                
                        endif                    
                    else                        
                        Me%Porosity(i,j,k) = Me%CohesiveClass%Porosity(i,j,k)                        
                    endif                        
                else                
                    Me%Porosity(i,j,k) = Me%PorositySand 
                endif
            endif
                
        enddo
        enddo
        enddo            
        
                
    end subroutine ConstructPorosity

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ConstructClassMass
    
        !Local-----------------------------------------------------------------
        integer                     :: i, j, k, n
        class(T_Sand), pointer      :: SandClass
        integer                     :: ILB, IUB, JLB, JUB, KLB, KUB
        real                        :: Area, Volume, TotalDryDensity, TotalDryMass
        !----------------------------------------------------------------------
        
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB
            
        do j=JLB, JUB
        do i=ILB, IUB
        do k=KLB, KUB
              
            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then
                    
                TotalDryDensity = (1 - Me%Porosity(i,j,k)) * Me%density
                    
                Area = Me%ExternalVar%DUX(i, j) * Me%ExternalVar%DVY(i, j)
                    
                Volume = Me%ExternalVar%DWZ(i,j,k) * Area
                    
                TotalDryMass = TotalDryDensity * Volume 
                    
                if (Me%CohesiveClass%Run) then
                 
                    Me%CohesiveClass%Mass(i,j,k) = Me%CohesiveClass%Field3D(i,j,k) * TotalDryMass
                
                endif
                
                do n=1,Me%NumberOfClasses

                    SandClass => Me%SandClass(n) 
                    
                    SandClass%Mass(i,j,k) = SandClass%Field3D(i,j,k) * TotalDryMass
                 
                enddo                
                
            endif
          
        enddo
        enddo
        enddo
       
    end subroutine ConstructClassMass
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructCriticalShearStress
    
        call New_Geometry
        
        call ComputeOpenSediment
    
        call ComputeCriticalShearStress
    
    end subroutine ConstructCriticalShearStress
    
    !--------------------------------------------------------------------------


    subroutine ConstructPropertyValue(NewProperty, ExtractType)

        !Arguments-------------------------------------------------------------
        class(T_Property)               :: NewProperty
        integer                         :: ExtractType

        !External--------------------------------------------------------------
        integer                         :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB    
        
        !----------------------------------------------------------------------
        
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB
              
        allocate(NewProperty%Field3D(ILB:IUB, JLB:JUB, KLB:KUB))
        
        if (Me%Evolution%Old) then
                            
            call ReadInitialField3D(FieldName = trim(NewProperty%ID%Name), Field3D = NewProperty%Field3D)
            
        else

            call ConstructFillMatrix  (PropertyID           = NewProperty%ID,              &
                                       EnterDataID          = Me%ObjEnterData,             &
                                       TimeID               = Me%ObjTime,                  &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                       GeometryID           = Me%ObjSedimentGeometry,      &
                                       ExtractType          = ExtractType,                 &
                                       PointsToFill3D       = Me%ExternalVar%WaterPoints3D,&
                                       Matrix3D             = NewProperty%Field3D,         &
                                       TypeZUV              = TypeZ_,                      &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                     &
                stop 'ConstructClassValue - ModuleSediment - ERR02'

            call GetDefaultValue(NewProperty%ID%ObjFillMatrix, NewProperty%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then                                                
                print *, NewProperty%ID%Name, ' was not build correctly.'
                stop 'ConstructClassValue - ModuleSediment - ERR03'
            endif

            if(NewProperty%ID%SolutionFromFile)then

                stop 'ConstructClassValue - ModuleSediment - ERR04'
            
            else

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'ConstructClassValue - ModuleSediment - ERR05'

            endif
        endif
        
        call CheckFieldConsistence (NewProperty)


    end subroutine ConstructPropertyValue

    !--------------------------------------------------------------------------
    subroutine CheckFieldConsistence(NewProperty)

        !Arguments-------------------------------------------------------------
        class(T_Property)                       :: NewProperty

        !Local-----------------------------------------------------------------
        integer                                 :: i,j,k
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------
        
        ILB = Me%SedimentSize3D%ILB
        IUB = Me%SedimentSize3D%IUB
        JLB = Me%SedimentSize3D%JLB
        JUB = Me%SedimentSize3D%JUB
        KLB = Me%SedimentSize3D%KLB
        KUB = Me%SedimentSize3D%KUB

        !Class percentage in empty layers is set as null_real
        do i = ILB, IUB
        do j = JLB, JUB
        do k = KLB, KUB
            
            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then
                               
                if (Me%ExternalVar%VolumeZ(i, j, k) == 0.) then

                    NewProperty%Field3D(i, j, k) = null_real
                    
                end if

            else

                NewProperty%Field3D(i, j, k) = null_real

            endif

        enddo
        enddo
        enddo

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------
    
    subroutine Construct_Initial_Geometry 
        
        !Local----------------------------------------------------------------
        integer                             :: STAT_CALL

        !Begin----------------------------------------------------------------
        
        call ReadUnLockExternalVar
        
        !Initial sediment thickness equals the one specified in bathymetry 
        !so elevation equals zero
        Me%Elevation(:,:) = 0.

        if(Me%Evolution%Old)then
            
            call ReadGeometryHDF(Me%ObjSedimentGeometry, trim(Me%Files%Initial)//"5", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleSediment - ERR01'

        end if

        !Get WaterPoints3D
        call GetWaterPoints3D(Me%ObjSedimentMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                             &
            call SetError (FATAL_, INTERNAL_, "Construct_Initial_Geometry - ModuleSediment - ERR02")

        
        !Computes sediment initial geometry
        call ComputeInitialGeometry(Me%ObjSedimentGeometry,                              &
                                    WaterPoints3D    = Me%ExternalVar%WaterPoints3D,     &
                                    SurfaceElevation = Me%Elevation,                     &
                                    ContinuesCompute = Me%Evolution%Old,                 &
                                    ActualTime       = Me%ExternalVar%Now,               &
                                    STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleSediment - ERR03'
        
        !Unget WaterPoints3D
        call UnGetMap(Me%ObjSedimentMap,Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleSediment - ERR04'


        !!Computes OpenPoints3D
        call UpdateComputeFaces3D(Me%ObjSedimentMap, STAT = STAT_CALL)      
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Initial_Geometry - ModuleSediment - ERR05'
        
        call ReadLockExternalVar

        Me%VerticalCoordinate(:,:,:) = Me%ExternalVar%SZZ(:,:,:)
        Me%KTop(:,:)                 = Me%ExternalVar%KTop(:,:)

            
    end subroutine Construct_Initial_Geometry
    
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    !If the user want's to use the values of a previous   
    ! run the read the Sediment properties values form the final      
    ! results file of a previous run. By default this      
    ! file is in HDF format                                

    subroutine ReadInitialField(FieldName, Field2D)

        !Arguments-------------------------------------------------------------
        character (Len = *)                         :: FieldName
        real, dimension(:,:), pointer               :: Field2D

        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: STAT_CALL
        integer                                     :: ObjHDF5
        integer(4)                                  :: HDF5_READ

        !----------------------------------------------------------------------


        WorkILB = Me%SedimentWorkSize3D%ILB 
        WorkIUB = Me%SedimentWorkSize3D%IUB 
        WorkJLB = Me%SedimentWorkSize3D%JLB 
        WorkJUB = Me%SedimentWorkSize3D%JUB 

        inquire (FILE=trim(Me%Files%Initial)//"5", EXIST = EXIST)

cd0:    if (EXIST) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSediment - ERR03.'

            Field2D(:,:) = FillValueReal

            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB,            &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSediment - ERR03.'

            call HDF5ReadData   (ObjHDF5, "/Results",trim(FieldName),                    &
                                 Array2D = Field2D,                                      &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSediment - ERR03.'
                
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField; ModuleSediment - ERR06.'
        
        else if(.not. EXIST) then cd0

                stop 'ReadInitialField; ModuleSediment - ERR07.'

        endif cd0


        !----------------------------------------------------------------------

    end subroutine ReadInitialField

    !--------------------------------------------------------------------------
    
    subroutine ReadInitialField3D(FieldName, Field3D)

        !Arguments-------------------------------------------------------------
        character (Len = *)                         :: FieldName
        real, dimension(:,:,:), pointer             :: Field3D

        !Local-----------------------------------------------------------------
        logical                                     :: EXIST
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: STAT_CALL
        integer                                     :: ObjHDF5
        integer(4)                                  :: HDF5_READ

        !----------------------------------------------------------------------


        WorkILB = Me%SedimentWorkSize3D%ILB 
        WorkIUB = Me%SedimentWorkSize3D%IUB 
        WorkJLB = Me%SedimentWorkSize3D%JLB 
        WorkJUB = Me%SedimentWorkSize3D%JUB 
        WorkKLB = Me%SedimentWorkSize3D%KLB 
        WorkKUB = Me%SedimentWorkSize3D%KUB 

        inquire (FILE=trim(Me%Files%Initial)//"5", EXIST = EXIST)

cd0:    if (EXIST) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField3D; ModuleSediment - ERR03.'

            Field3D(:,:,:) = FillValueReal

            ! Reads from HDF file the Property concentration and open boundary values
            call HDF5SetLimits  (ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, WorkKLB, WorkKUB, &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField3D; ModuleSediment - ERR03.'

            call HDF5ReadData   (ObjHDF5, "/Results3D", trim(FieldName),                 &
                                    Array3D = Field3D,                                   &
                                    STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField3D; ModuleSediment - ERR03.'
            
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadInitialField3D; ModuleSediment - ERR06.'
        
        else if(.not. EXIST) then cd0

                stop 'ReadInitialField3D; ModuleSediment - ERR07.'

        endif cd0


        !----------------------------------------------------------------------

    end subroutine ReadInitialField3D

    !--------------------------------------------------------------------------



    subroutine ReadResidualStartTime()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                    :: EXIST
        integer                                    :: STAT_CALL
        integer                                    :: ObjHDF5
        integer(4)                                 :: HDF5_READ
        real,    dimension(:    ), pointer         :: AuxTime

        !----------------------------------------------------------------------

        inquire (FILE=trim(Me%Files%Initial)//"5", EXIST = EXIST)

cd0:    if (EXIST) then

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            ObjHDF5 = 0

            !Opens HDF5 File
            call ConstructHDF5 (ObjHDF5,                                                 &
                                trim(Me%Files%Initial)//"5", HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadResidualStartTime; ModuleSediment - ERR10.'


            call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadResidualStartTime - ModuleSediment - ERR30'
            
            allocate(AuxTime(6))

            call HDF5ReadData  (ObjHDF5, "/Results",                           &
                                 "Start Time",                                 &
                                 Array1D = AuxTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadResidualStartTime - ModuleSediment - ERR40'
            
            call SetDate   (Me%Residual%StartTime, AuxTime(1), AuxTime(2), AuxTime(3), &
                            AuxTime(4), AuxTime(5), AuxTime(6))
                      
            deallocate(AuxTime)
                
            call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'ReadResidualStartTime; ModuleSediment - ERR50.'
        
        else if(.not. EXIST) then cd0

                stop 'ReadResidualStartTime; ModuleSediment - ERR60.'

        endif cd0

        !----------------------------------------------------------------------

    end subroutine ReadResidualStartTime

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    
    subroutine GetTopCriticalShear(ObjSedimentID, TopCriticalShear, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        real, dimension(:,:  ),  pointer            :: TopCriticalShear
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            call Read_Lock(mSediment_, Me%InstanceID)

            TopCriticalShear => Me%CohesiveClass%CriticalShearStress


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetTopCriticalShear
    
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    subroutine GetCohesiveMass(ObjSedimentID, CohesiveMass, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        real(8), dimension(:,:,:),  pointer         :: CohesiveMass
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)

            CohesiveMass => Me%CohesiveClass%Mass

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetCohesiveMass
    
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    subroutine GetCohesiveContent(ObjSedimentID, CohesiveClassRun, CohesiveContent, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        logical                                     :: CohesiveClassRun
        real, dimension(:,:,:),  pointer            :: CohesiveContent
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)
            
            CohesiveClassRun = Me%CohesiveClass%Run

            if(Me%CohesiveClass%Run)then
                CohesiveContent => Me%CohesiveClass%Field3D
            endif                

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetCohesiveContent
    

    !--------------------------------------------------------------------------
    
    subroutine GetSandMass(ObjSedimentID, SandClassID, SandMass, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        integer                                     :: SandClassID, n
        class(T_Sand), pointer                      :: SandClass
        real(8), dimension(:,:,:),pointer           :: SandMass
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)
            
 do1:        do n=1,Me%NumberOfClasses
               
               SandClass => Me%SandClass(n)
               
                if(SandClass%ID%IDNumber == SandClassID) then
                        
                    SandMass => SandClass%Mass
                    
                    exit do1
                    
                endif
    
            enddo do1


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetSandMass
    !--------------------------------------------------------------------------
        
    subroutine GetSandContent(ObjSedimentID, SandClassID, SandContent, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        integer                                     :: SandClassID, n
        class(T_Sand), pointer                      :: SandClass
        real, dimension(:,:,:),pointer              :: SandContent
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)
            
 do1:        do n=1,Me%NumberOfClasses
               
               SandClass => Me%SandClass(n)
               
                if(SandClass%ID%IDNumber == SandClassID) then
                        
                    SandContent => SandClass%Field3D
                    
                    exit do1
                    
                endif
    
            enddo do1


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetSandContent

    !--------------------------------------------------------------------------
    
    subroutine GetCriticalShearStress(ObjSedimentID, SandClassID, CriticalShearStress, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        integer                                     :: SandClassID, n
        class(T_Sand), pointer                      :: SandClass
        real, dimension(:,:),pointer                :: CriticalShearStress
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)
            
 do1:        do n=1,Me%NumberOfClasses
               
               SandClass => Me%SandClass(n)
               
                if(SandClass%ID%IDNumber == SandClassID) then
                        
                    CriticalShearStress => SandClass%CriticalShearStress
                    
                    exit do1
                    
                endif
    
            enddo do1

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
    
    end subroutine GetCriticalShearStress
     !--------------------------------------------------------------------------
    
    subroutine GetConcRef(ObjSedimentID, SandClassID, ConcRef, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        integer                                     :: SandClassID, n
        class(T_Sand), pointer                      :: SandClass
        real, dimension(:,:),  pointer              :: ConcRef
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)
            
do1:        do n=1,Me%NumberOfClasses
               
               SandClass => Me%SandClass(n)
               
                if(SandClass%ID%IDNumber == SandClassID) then

                    ConcRef => SandClass%ConcRef
                    
                    exit do1
                    
                endif
    
            enddo do1


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetConcRef
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    subroutine GetReferenceLevel(ObjSedimentID, ReferenceLevel, STAT) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        real                                        :: ReferenceLevel
        integer, optional, intent(OUT)              :: STAT

        !External--------------------------------------------------------------
        integer                                     :: ready_        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            !call Read_Lock(mSediment_, Me%InstanceID)

            ReferenceLevel = Me%ReferenceLevel


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
    
    end subroutine GetReferenceLevel
    
    !--------------------------------------------------------------------------
        
!    subroutine GetSandD50(ObjSedimentID, SandD50, STAT) 
!
!        !Arguments-------------------------------------------------------------
!        integer                                     :: ObjSedimentID
!        real, dimension(:,:),  pointer              :: SandD50
!        integer, optional, intent(OUT)              :: STAT
!
!        !External--------------------------------------------------------------
!        integer                                     :: ready_        
!
!        !Local-----------------------------------------------------------------
!        integer                                     :: STAT_
!
!        !----------------------------------------------------------------------
!
!        STAT_ = UNKNOWN_
!
!        call Ready(ObjSedimentID, ready_)
!        
!cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
!            (ready_ .EQ. READ_LOCK_ERR_)) then
!
!
!            call Read_Lock(mSediment_, Me%InstanceID)
!
!            SandD50 => Me%SandD50
!
!
!            STAT_ = SUCCESS_
!        else 
!            STAT_ = ready_
!        end if cd1
!
!
!        if (present(STAT)) STAT = STAT_
!    
!    end subroutine GetSandD50
    
    !--------------------------------------------------------------------------
 
    !--------------------------------------------------------------------------
    
    subroutine SetCohesiveFlux(ObjSedimentID, ConsolidationFlux,  STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                 :: ObjSedimentID
        real, dimension(:,:),optional, pointer  :: ConsolidationFlux
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_            
        integer                                 :: ready_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then

            Me%ExternalVar%ConsolidationFlux => ConsolidationFlux
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    
    
    end subroutine SetCohesiveFlux

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    subroutine SetNonCohesiveFlux(ObjSedimentID, SandClassID, FluxToSediment,  STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                 :: ObjSedimentID
        integer                                 :: SandClassID, n
        class(T_Sand), pointer                  :: SandClass
        real, dimension(:,:),optional, pointer  :: FluxToSediment
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_            
        integer                                 :: ready_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then
    
do1:        do n=1,Me%NumberOfClasses
               
               SandClass => Me%SandClass(n)
               
                if(SandClass%ID%IDNumber == SandClassID) then
                        
                    SandClass%FluxToSediment => FluxToSediment
                    
                    exit do1
                    
                endif
    
            enddo do1
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

    
    
    end subroutine SetNonCohesiveFlux

    !--------------------------------------------------------------------------
       
    subroutine SetWaveTensionON(ObjSedimentID, WaveTensionON, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                                 :: ObjSedimentID
        logical                                 :: WaveTensionON
        integer, optional, intent(OUT)          :: STAT

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_            
        integer                                 :: ready_
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)  
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Me%ExternalVar%WaveTensionON = WaveTensionON
            
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_
    
        if(Me%WavesOn) then
            if (.not. Me%ExternalVar%WaveTensionON) then                
                write(*,*)
                write(*,*) 'Define WAVETENSION: 1 in module InterfaceSedimentWater'
                stop 'SetWaveTensionON - ModuleSediment - ERR10'
            endif
        endif
        
    end subroutine SetWaveTensionON

    !--------------------------------------------------------------------------
    subroutine UnGetSediment2D_I(ObjSedimentID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSedimentID
        integer, dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSediment_, Me%InstanceID, "UnGetSediment2D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSediment2D_I

    !--------------------------------------------------------------------------

    subroutine UnGetSediment2D_R8(ObjSedimentID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSedimentID
        real(8), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSediment_, Me%InstanceID,  "UnGetSediment2D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSediment2D_R8

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine UnGetSediment2D_R4(ObjSedimentID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSedimentID
        real(4), dimension(:, :), pointer               :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mSediment_, Me%InstanceID,  "UnGetSediment2D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSediment2D_R4
    
        !--------------------------------------------------------------------------

    
    subroutine UnGetSediment3Dreal8(ObjSedimentID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: ObjSedimentID
        real(8), pointer, dimension(:,:,:) :: Array
        integer, optional, intent (OUT)    :: STAT

        !External--------------------------------------------------------------
        integer                            :: ready_   

        !Local-----------------------------------------------------------------
        integer                            :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mSediment_, Me%InstanceID, "UnGetSediment3Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UnGetSediment3Dreal8


    !--------------------------------------------------------------------------
    
    subroutine UnGetSediment3Dreal4(ObjSedimentID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: ObjSedimentID
        real(4), pointer, dimension(:,:,:) :: Array
        integer, optional, intent (OUT)    :: STAT

        !External--------------------------------------------------------------
        integer                            :: ready_   

        !Local-----------------------------------------------------------------
        integer                            :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_UnLock(mSediment_, Me%InstanceID, "UnGetSediment3Dreal8")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine UnGetSediment3Dreal4


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifySediment(ObjSedimentID, ShearStress, VelU, VelV, VelMod,    &
                              TauWave, ShearStressMean, Cphi, CWphi, Wphi, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        real, dimension(:,:), pointer               :: ShearStress,  &
                                                       VelU, VelV, VelMod,           &
                                                       TauWave, ShearStressMean, &
                                                       Cphi, CWphi, Wphi   
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        logical                                     :: ChangeBathym                                                
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call ReadLockExternalVar 

do1:            do while (Me%ExternalVar%Now >= Me%Evolution%NextSediment) 

                    Me%ExternalVar%ShearStress     => ShearStress                    

                    Me%ExternalVar%VelU            => VelU
                    Me%ExternalVar%VelV            => VelV
                    Me%ExternalVar%VelMod          => VelMod

                    Me%ExternalVar%TauWave         => TauWave
                    Me%ExternalVar%ShearStressMean => ShearStressMean
                    
                    Me%ExternalVar%Cphi            => Cphi
                    Me%ExternalVar%CWphi           => CWphi
                    Me%ExternalVar%Wphi            => Wphi
                    
                    call ComputeOpenSediment

                    !call ComputeD50Cell
                        
                    call ComputeNDShearStress
                    
                    call ComputeCriticalShearStress
                        
                    if (Me%BedloadMethod /= 0) then
                    
                        call ComputeBedload
                        
                        call ComputeFluxes
                        
                    endif
                    
                    call ComputeEvolution
                        
                    if (Me%Discharges%Yes) then 
                        call ComputeDischarges
                    endif           

                    call BoundaryCondition
                    
                    call ComputeMass
                    
                    call ComputePercentage
                    
                    !call ComputePorosity

                    call ComputeTotalDZ             
    
                    call ComputeVerticalCoordinate
                        
                    call New_Geometry
                        
                    call ComputeResidualEvolution

                    if (Me%Boxes%Yes ) call OutputBoxFluxes

                    if (Me%OutPut%Yes) call OutPutSedimentHDF

                    if (Me%TimeSerie)  call OutPut_TimeSeries

                    if (Me%Evolution%Bathym) then

                        ChangeBathym = .false.
                            
                        if (Me%ExternalVar%Now >= Me%Evolution%NextBatim) ChangeBathym = .true.
                            
                        if (ChangeBathym) then

                            !Bathymetry 
                            call UnGetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySediment - ModuleSediment - ERR10'

                            call ModifyGridData(Me%ObjBathym,Me%BatimIncrement, Add = .false.,  &
                                                STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySediment - ModuleSediment - ERR20.'

                            !Bathymetry
                            call GetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
                            if (STAT_CALL /= SUCCESS_) stop 'ModifySediment - ModuleSediment - ERR30'

                            Me%BatimIncrement(:,:) = 0.

                            Me%Evolution%NextBatim = Me%Evolution%NextBatim + Me%Evolution%BathymDT


                        endif

                    endif
                    
                    call ComputeD50Cell
                    
                    !call ComputeCriticalShearStress
                    
                    call ComputeReferenceConcentration

                    Me%Evolution%NextSediment = Me%Evolution%NextSediment + Me%Evolution%SedimentDT

                enddo do1            

            call ReadUnLockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySediment

    !--------------------------------------------------------------------------
    
    subroutine ComputeOpenSediment
       !Local-----------------------------------------------------------------
        integer             :: i, j 
        integer             :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        do j=WJLB, WJUB
        do i=WILB, WIUB
                    
            WKUB = Me%KTop(i, j)

            if (Me%ExternalVar%WaterPoints3D (i,j,WKUB) + &
                Me%ExternalVar%OpenPoints2D(i,j) > 1) then
                
                Me%OpenSediment(i,j) = 1            
            else                    
                Me%OpenSediment(i,j) = 0
            endif
                        
        enddo
        enddo 

    end subroutine ComputeOpenSediment
    
    !--------------------------------------------------------------------------
    
     !--------------------------------------------------------------------------
    
    subroutine ComputeD50Cell
       !Local-----------------------------------------------------------------
        integer             :: i, j, n, WKUB
        integer             :: WILB, WIUB, WJLB, WJUB
        class(T_Sand), pointer :: SandClass
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%D50(:,:)= 0      
        Me%SandD50(:,:)= 0
        
        do j=WJLB, WJUB
        do i=WILB, WIUB
                    
            WKUB = Me%KTop(i, j)

            if (Me%ExternalVar%WaterPoints3D (i,j,WKUB) == WaterPoint) then
                        
                do n=1,Me%NumberOfClasses

                    SandClass => Me%SandClass(n)
                        
                    Me%SandD50(i,j) = SandClass%D50 * SandClass%Field3D(i, j, WKUB) + Me%SandD50(i,j)
                    
                enddo
            
                Me%D50(i,j) = Me%SandD50(i,j)
            
                if (Me%CohesiveClass%Run) then
                    
                    Me%D50(i,j) = Me%CohesiveClass%D50 * Me%CohesiveClass%Field3D(i, j, WKUB) + Me%D50(i,j)
                    
                endif
                    
                
            endif
        enddo
        enddo 

    end subroutine ComputeD50Cell
    
    !--------------------------------------------------------------------------
    !Compute the nondimensional shear stress
    subroutine ComputeNDShearStress
       !Local-----------------------------------------------------------------
        integer             :: i, j
        integer             :: WILB, WIUB, WJLB, WJUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%NDShearStress(:,:) = 0.
        
        if (Me%WavesOn) then
            Me%NDShearStressMean(:,:) = 0.
            Me%NDShearStressWaves(:,:) = 0.
        endif
        
        do j=WJLB, WJUB
        do i=WILB, WIUB
            
            if ((Me%OpenSediment(i,j) == OpenPoint .and. Me%SandD50(i,j) > 0.)) then
                
                Me%NDShearStress(i,j) = Me%ExternalVar%ShearStress(i,j)/(Me%DeltaDensity*Gravity*Me%SandD50(i,j))
                
                if (Me%WavesOn) then
                    
                    Me%NDShearStressMean(i,j) = Me%ExternalVar%ShearStressMean(i,j)/(Me%DeltaDensity*Gravity*Me%SandD50(i,j))
                
                    Me%NDShearStressWaves(i,j) =  Me%ExternalVar%TauWave(i,j)/(Me%DeltaDensity*Gravity*Me%SandD50(i,j))
                endif                 
            endif
        enddo
        enddo
    end subroutine ComputeNDShearStress
     !--------------------------------------------------------------------------
    
    subroutine ComputeCriticalShearStress
        
        !Local-----------------------------------------------------------------
        integer             :: i, j, n
        real                :: ShieldsParameter
        class(T_Sand), pointer :: SandClass
        integer             :: WILB, WIUB, WJLB, WJUB, WKUB
        real                :: pm, pm1, pm1max, pm2, alfa, beta, tec
        real, dimension(:), pointer  :: aux
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%Dast(:,:) = FillValueReal
        
        do n=1,Me%NumberOfClasses
            
            SandClass => Me%SandClass(n)
            
            SandClass%CriticalShearStress(:,:)= FillValueReal
            SandClass%NDCriticalShearStress(:,:)= FillValueReal
        enddo
        
        if (Me%CohesiveClass%Run) then
             Me%CohesiveClass%CriticalShearStress(:,:)= FillValueReal
             
            !Critical erosion shear stress for weak consolidated cohesive sediment 
            tec = Me%CohesiveClass%WeakConsolidated_CSS
            !Maximum mud percentage in which the bed starts to have a cohesive behaviour
            pm1max = Me%CohesiveClass%PM1_MAX
            !Mud percentage for a fully cohesive bed 
            pm2 = Me%CohesiveClass%PM2
            
        endif
        
        do j=WJLB, WJUB
        do i=WILB, WIUB

            if (Me%OpenSediment(i,j) == OpenPoint) then
                
                WKUB = Me%KTop(i, j)
                
                if (Me%SandD50(i,j) > 0.) then
                    
                    !critical Shields parameter modified by van Rijn (2003, 2007)
                    Me%Dast(i,j) = Me%SandD50(i,j)*((Me%RelativeDensity-1)*Gravity/WaterCinematicVisc**2)**(1./3.)
                        
                    If (Me%Dast(i,j).LE.4) then
                        ShieldsParameter = 0.115*Me%Dast(i,j)**(-0.5)
                    elseIf (Me%Dast(i,j).GT.4.AND.Me%Dast(i,j).LE.10) then
                        ShieldsParameter = 0.14*Me%Dast(i,j)**(-0.64)
                    elseIf (Me%Dast(i,j).GT.10.AND.Me%Dast(i,j).LE.20) then
                        ShieldsParameter = 0.04*Me%Dast(i,j)**(-0.1)
                    elseIf (Me%Dast(i,j).GT.20.AND.Me%Dast(i,j).LE.150) then
                        ShieldsParameter = 0.013*Me%Dast(i,j)**0.29
                    elseIf (Me%Dast(i,j).GT.150) then
                        ShieldsParameter = 0.055
                    endif
                        
                    do n=1,Me%NumberOfClasses
                    
                        SandClass => Me%SandClass(n)
                        
                        if(SandClass%Field3D(i,j,WKUB) > 0.) then
                        
                            call ComputeHidingExposure(i,j,n,SandClass%D50)
                        
                            ! Critical Shear Stress of a pure sand bed [N/m2]
                            SandClass%CriticalShearStress(i, j)= Me%DeltaDensity*Gravity*SandClass%D50      &
                                                                 *ShieldsParameter*SandClass%HidingFactor(i,j)
                        endif
                    enddo
                    
                    if (Me%CohesiveClass%Run) then
                    
                        alfa = 2.*10**3
                        beta = 1.0
                
                        pm = Me%CohesiveClass%Field3D(i,j,WKUB)
                        
                        !Mud percentage in which the bed starts to have a cohesive behaviour
                        pm1 = alfa*Me%SandD50(i,j)  
                        pm1 = min(pm1,pm1max)
                        
                        Me%CohesiveClass%PM1(i,j) = pm1
                
                        if (pm > 0.) then
                        
                            !Non-cohesive bed
                            if (pm <= pm1) then
                                
                                allocate (aux(1:Me%NumberOfClasses))
                    
                                do n=1,Me%NumberOfClasses
                                
                                    SandClass => Me%SandClass(n)
                                    
                                    if(SandClass%Field3D(i,j,WKUB) > 0.) then
                            
                                        SandClass%CriticalShearStress(i, j)= (1 + pm)**beta*SandClass%CriticalShearStress(i, j)
                                
                                        aux(n) = SandClass%CriticalShearStress(i, j)
                                        
                                    else                                   
                                        aux(n) = 9999
                                    endif
                                enddo
                        
                                Me%CohesiveClass%CriticalShearStress(i, j) = minval(aux)
                                
                                deallocate (aux)
                        
                            !Cohesive bed    
                            elseif (pm <= pm2) then
                        
                                do n=1,Me%NumberOfClasses
                                
                                    SandClass => Me%SandClass(n)
                                    
                                    if(SandClass%Field3D(i,j,WKUB) > 0.) then
                        
                                        SandClass%CriticalShearStress(i, j) = (pm2 - pm)/(pm2 - pm1)*((1 + pm1)**beta         &
                                                                            *SandClass%CriticalShearStress(i, j) - tec) + tec
                                    endif                        
                                enddo
                            
                                Me%CohesiveClass%CriticalShearStress(i, j) = tec
                                                        
                            !Fully cohesive bed
                            elseif (pm > pm2) then
                            
                                 do n=1,Me%NumberOfClasses
                                 
                                    SandClass => Me%SandClass(n)
                                    
                                    if(SandClass%Field3D(i,j,WKUB) > 0.) then
                        
                                        SandClass%CriticalShearStress(i, j) = tec
                                    endif                        
                                 enddo
                             
                                 Me%CohesiveClass%CriticalShearStress(i, j) = tec                        
                            endif                   
                        endif
                        
                    endif !(Me%CohesiveClass%Run)
                        
                    do n=1,Me%NumberOfClasses
                    
                        SandClass => Me%SandClass(n)
                        
                        if(SandClass%Field3D(i,j,WKUB) > 0.) then
                
                            !Non dimensional critical shear stress
                            SandClass%NDCriticalShearStress(i, j)= SandClass%CriticalShearStress(i, j)/     &
                                                                    (Me%DeltaDensity*Gravity*SandClass%D50)
                            
                        endif
                    enddo    
                    
                else !(Me%SandD50(i,j) > 0.)
                
                    if (Me%CohesiveClass%Run)then
                        if(Me%CohesiveClass%Field3D(i,j,WKUB) > 0.)     &
                            Me%CohesiveClass%CriticalShearStress(i, j) = tec
                    endif
                        
                
                endif !(Me%SandD50(i,j) > 0.)                                                
                
            endif !(Me%OpenSediment(i,j) == OpenPoint)
                
        enddo
        enddo      

    end subroutine ComputeCriticalShearStress
    
        !--------------------------------------------------------------------------
    
    !Compute the hiding-exposure factor (Wu, et al. 2000)
    subroutine ComputeHidingExposure(i,j,n,D50)
        
        !Argument--------------------------------------------------------------
        integer             :: i, j, n
        real(8)             :: D50
        !----------------------------------------------------------------------
        !Local-----------------------------------------------------------------
        integer             :: aux
        class(T_Sand), pointer :: SandClass
        real                :: exposure, hiding, m
        integer             :: WKUB
        !----------------------------------------------------------------------
        
        m=-0.6          
        WKUB = Me%KTop(i, j)                        
        exposure = 0.
        hiding = 0.
                       
        do aux=1,Me%NumberOfClasses
            
            SandClass => Me%SandClass(aux)
                        
            !Compute hiding and exposure probabilities
            exposure = SandClass%Field3D(i,j,WKUB)*D50/             &
                        (D50 + SandClass%D50) + exposure
                        
            hiding = SandClass%Field3D(i,j,WKUB)*SandClass%D50/ &
                        (D50 +SandClass%D50) + hiding
                        
        enddo           
                        
        Me%SandClass(n)%HidingFactor(i,j) = (exposure/hiding)**m                            
                        

    end subroutine ComputeHidingExposure
    
    
       !--------------------------------------------------------------------------
    subroutine ComputeReferenceConcentration
       !Local-----------------------------------------------------------------
        integer                :: i, j, n
        integer                :: WILB, WIUB, WJLB, WJUB, WKUB
        class(T_Sand), pointer :: SandClass
        real                   :: T, Dast, pm
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        do n=1,Me%NumberOfClasses
            
            SandClass => Me%SandClass(n)
            SandClass%ConcRef(:,:)= 0.
        enddo
        
        if(Me%RefConcMethod == 1) then
        !The reference concentration (volumetric) is computed based on the formulation of Zyserman and Fredsoe (1994)
        !to calculate the erosion flux
            
            do j=WJLB, WJUB
            do i=WILB, WIUB
            
                if ((Me%OpenSediment(i,j) == OpenPoint .and. Me%SandD50(i,j) > 0.)) then
                    
                    WKUB = Me%KTop(i, j)
                    
                    do n=1,Me%NumberOfClasses
                                
                        SandClass => Me%SandClass(n)
                        
                         if(SandClass%Field3D(i,j,WKUB) > 0. .and. Me%NDShearStress(i,j) > 0.045) then

                            SandClass%ConcRef(i,j) = 0.331 * (Me%NDShearStress(i,j) - 0.045)**1.75 /        &
                                                    (1 + 0.331/0.46 * (Me%NDShearStress(i,j) - 0.045)**1.75)
                            
                            !The reference concentration upper limit is set to 0.2 (Amoudry, 2010)
                            SandClass%ConcRef(i,j) = min(SandClass%ConcRef(i,j), 0.2)
                            
                            ![kg m-3]
                            SandClass%ConcRef(i,j) = SandClass%ConcRef(i,j) * Me%Density * SandClass%Field3D(i,j,WKUB)
                        endif
                    enddo
                endif
                      
            enddo
            enddo
        endif
                    
                
        if(Me%RefConcMethod == 2) then
        !The reference concentration (volumetric) is computed based on the formulation of van Rijn(2007)
            
            do j=WJLB, WJUB
            do i=WILB, WIUB
            
                if ((Me%OpenSediment(i,j) == OpenPoint .and. Me%SandD50(i,j) > 0.)) then
                
                    WKUB = Me%KTop(i, j)                
                
                    do n=1,Me%NumberOfClasses
                                
                        SandClass => Me%SandClass(n)
                                    
                        if(SandClass%Field3D(i,j,WKUB) > 0. .and. &
                            Me%ExternalVar%ShearStress(i,j) > SandClass%CriticalShearStress(i, j)) then
                            
                            pm = 0.
                            
                            if (Me%CohesiveClass%Run) pm = Me%CohesiveClass%Field3D(i,j,WKUB)
                
                            T = (Me%ExternalVar%ShearStress(i,j) - SandClass%CriticalShearStress(i, j)) / &
                                SandClass%CriticalShearStress(i, j)
                                
                            Dast = SandClass%D50*((Me%RelativeDensity-1)*Gravity/WaterCinematicVisc**2)**(1./3.)
                                
                            SandClass%ConcRef(i,j) = Me%ConcRefFactor * 0.015 * (1 - pm) * SandClass%D50/Me%ReferenceLevel * &
                                                        T**1.5/Dast**0.3 
                            
                            SandClass%ConcRef(i,j) = min(SandClass%ConcRef(i,j), 0.05)

                            ![kg m-3]
                            SandClass%ConcRef(i,j) = SandClass%ConcRef(i,j) * Me%Density * SandClass%Field3D(i,j,WKUB)
                        
                        endif
                    enddo
                 
                endif
            enddo
            enddo
        endif
            
    end subroutine ComputeReferenceConcentration
     !--------------------------------------------------------------------------    

    subroutine ComputeBedload          
    
            if (Me%BedloadMethod==1) then
                call CurrentOnlyBedload
            elseif  (Me%BedloadMethod==2) then
                 call CurrentPlusWavesBedload
            elseif (Me%BedloadMethod==3) then
                 call WavesOnlyBedload
            endif

    end subroutine ComputeBedload
    !--------------------------------------------------------------------------

    subroutine CurrentOnlyBedload
    
        !Arguments-------------------------------------------------------------
        class(T_Sand), pointer :: SandClass
        
        !Local-----------------------------------------------------------------
        real    ::  A,n1,p,DeltaTau,NDBedload,Cohesiveness_Factor
        integer :: i, j, n
        integer :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%Bedload(:, :) = 0.
        
do1:    do n=1,Me%NumberOfClasses
    
            SandClass => Me%SandClass(n)
            
            A = SandClass%parameter_A
            n1 = SandClass%parameter_n
            p = SandClass%parameter_p
            
            SandClass%Bedload(:, :) = 0.
            
            do j=WJLB, WJUB
            do i=WILB, WIUB
                
                WKUB = Me%KTop(i, j)

                if (Me%OpenSediment(i,j) == OpenPoint) then
                    
                    if(SandClass%Field3D(i,j,WKUB) > 0.) then
                    
                        DeltaTau = Me%NDShearStress(i, j)-SandClass%NDCriticalShearStress(i, j)                  
                    
                        if (DeltaTau.GT.0.) then
                            
                            !Dimensionless bedload tranport rate
                            NDBedload = A*Me%NDShearStress(i, j)**n1*DeltaTau**p * SandClass%Field3D(i,j,WKUB)
                            
                            !Bedload transport rate per unit width [kg/s/m]
                            SandClass%BedLoad(i, j) = NDBedload*(gravity*(Me%RelativeDensity-1)*    &
                                Me%SandD50(i, j)**3)**(1./2.)*Me%Density
                            
                            if (Me%CohesiveClass%Run) then
                                
                                Cohesiveness_Factor =                                                         &
                                Cohesiveness_Adjust(Cohesive_Percentage = Me%CohesiveClass%Field3D(i,j,WKUB), &
                                                    PM1                 = Me%CohesiveClass%PM1(i,j),          &
                                                    PM2                 = Me%CohesiveClass%PM2)
                                
                                SandClass%BedLoad(i, j) = Cohesiveness_Factor * SandClass%BedLoad(i, j)    
                            
                            endif
                        
                            Me%Bedload(i, j) = Me%Bedload(i, j) + SandClass%Bedload(i, j)
                            
                        endif
                    endif
                endif
            enddo
            enddo
        enddo do1

    end subroutine CurrentOnlyBedload

    !--------------------------------------------------------------------------

    subroutine CurrentPlusWavesBedload
    
        !Arguments-------------------------------------------------------------
        class(T_Sand), pointer :: SandClass
        
        !Local-----------------------------------------------------------------
        real    :: A, n1, p, DeltaTau, NDBedload, Cohesiveness_Factor
        real    :: NDBedload1, NDBedload2, Cphi, CWphi, Asym_Factor
        real    :: NDBedloadNormal,NDBedloadParallel, NDBedloadX, NDBedloadY, aux
        integer :: i, j, n
        integer :: WILB, WIUB, WJLB, WJUB, WKUB
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin----------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%Bedload(:, :) = 0.
        
        Asym_Factor = 0.2
        
do1:    do n=1,Me%NumberOfClasses
    
            SandClass => Me%SandClass(n)
            
            A = SandClass%parameter_A
            n1 = SandClass%parameter_n
            p = SandClass%parameter_p
            
            SandClass%Bedload(:, :) = 0.
            
            do j=WJLB, WJUB
            do i=WILB, WIUB
                
                WKUB = Me%KTop(i, j)

                if (Me%OpenSediment(i,j) == OpenPoint) then
                    
                    if(SandClass%Field3D(i,j,WKUB) > 0.) then                        
                    
                        DeltaTau = Me%NDShearStress(i, j)-SandClass%NDCriticalShearStress(i, j)                  
                    
                        if (DeltaTau.GT.0.) then
                            
                            Cphi = Me%ExternalVar%Cphi(i,j) * pi/180. !Current angle in radians
                        
                            CWphi = Me%ExternalVar%CWphi(i,j) * pi/180. !Current-wave angle in radians
                            
                            !Dimensionless bedload tranport rate paralell to the current direction
                            NDBedload1 = A*Me%NDShearStressMean(i, j)**n1*DeltaTau**p * SandClass%Field3D(i,j,WKUB)
                            
                            NDBedload2 = A*(0.9534 + 0.1904 * cos(2*CWphi)) * &
                                         Me%NDShearStressWaves(i,j)**0.5 * Me%NDShearStressMean(i, j) + &
                                         A*0.229*Asym_Factor* Me%NDShearStressWaves(i,j)**1.5 * cos(CWphi)
                            
                            NDBedloadParallel = max(NDBedload1, abs(NDBedload2))
                            
                            !Dimensionless bedload tranport rate normal to the current direction
                            NDBedloadNormal = A*0.1907 * Me%NDShearStressWaves(i,j)**2 /        &
                                            (Me%NDShearStressWaves(i,j)**1.5 + 1.5*Me%NDShearStressMean(i,j)**1.5) *    &
                                            (Me%NDShearStressMean(i,j)*sin(2*CWphi) + 1.2*Asym_Factor*Me%NDShearStressWaves(i,j)*sin(CWphi)) 
                            
                            NDBedload = (NDBedloadParallel**2 + NDBedloadNormal**2)**0.5
                            
                            !Rotation to cartesian convention
                            NDBedloadX = NDBedloadParallel * cos(Cphi) - NDBedloadNormal * sin(Cphi)
                            NDBedloadY = NDBedloadParallel * sin(Cphi) + NDBedloadNormal * cos(Cphi)
                            
                            aux = (gravity*(Me%RelativeDensity-1)*    &
                                Me%SandD50(i, j)**3)**(1./2.)*Me%Density
                            
                            !Bedload transport rate per unit width [kg/s/m]
                            SandClass%BedLoad(i, j)  = NDBedload  * aux
                            SandClass%BedLoadX(i, j) = NDBedloadX * aux
                            SandClass%BedLoadY(i, j) = NDBedloadY * aux
                            
                            if (Me%CohesiveClass%Run) then
                                
                                Cohesiveness_Factor =                                                         &
                                Cohesiveness_Adjust(Cohesive_Percentage = Me%CohesiveClass%Field3D(i,j,WKUB), &
                                                    PM1                 = Me%CohesiveClass%PM1(i,j),          &
                                                    PM2                 = Me%CohesiveClass%PM2)
                                
                                SandClass%BedLoad(i, j)  = Cohesiveness_Factor * SandClass%BedLoad(i, j)
                                SandClass%BedLoadX(i, j) = Cohesiveness_Factor * SandClass%BedLoadX(i, j)
                                SandClass%BedLoadY(i, j) = Cohesiveness_Factor * SandClass%BedLoadY(i, j)
                            
                            endif                            
                        
                            Me%Bedload(i, j) = Me%Bedload(i, j) + SandClass%Bedload(i, j)
                            
                        endif
                    endif
                endif
            enddo
            enddo
        enddo do1
        
    end subroutine CurrentPlusWavesBedload

    !--------------------------------------------------------------------------
    
    subroutine WavesOnlyBedload
        !Arguments-------------------------------------------------------------
        class(T_Sand), pointer :: SandClass
        
        !Local-----------------------------------------------------------------
        real    :: A, DeltaTau, NDBedload,Cohesiveness_Factor
        real    :: Wphi, aux, Asym_Factor
        integer :: i, j, n
        integer :: WILB, WIUB, WJLB, WJUB, WKUB
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin----------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%Bedload(:, :) = 0.
        
        Asym_Factor = 0.2
        
do1:    do n=1,Me%NumberOfClasses
    
            SandClass => Me%SandClass(n)
            
            A = SandClass%parameter_A
            
            SandClass%Bedload(:, :) = 0.
            
            do j=WJLB, WJUB
            do i=WILB, WIUB
                
                WKUB = Me%KTop(i, j)

                if (Me%OpenSediment(i,j) == OpenPoint) then
                    
                    if(SandClass%Field3D(i,j,WKUB) > 0.) then                        
                    
                        DeltaTau = Me%NDShearStress(i, j)-SandClass%NDCriticalShearStress(i, j)                  
                    
                        if (DeltaTau.GT.0.) then
                        
                            !Wave angle referenced to the grid in radians
                            Wphi = Me%ExternalVar%Wphi(i,j) * pi/180.
                            
                            NDBedload = A*0.229*Asym_Factor* Me%NDShearStressWaves(i,j)**1.5                          

                            aux = (gravity*(Me%RelativeDensity-1)*    &
                                Me%SandD50(i, j)**3)**(1./2.)*Me%Density
                            
                            !Bedload transport rate per unit width [kg/s/m]
                            SandClass%BedLoad(i, j)  = NDBedload  * aux                            
                           
                            if (Me%CohesiveClass%Run) then
                                
                                Cohesiveness_Factor =                                                         &
                                Cohesiveness_Adjust(Cohesive_Percentage = Me%CohesiveClass%Field3D(i,j,WKUB), &
                                                    PM1                 = Me%CohesiveClass%PM1(i,j),          &
                                                    PM2                 = Me%CohesiveClass%PM2)
                                
                                SandClass%BedLoad(i, j) = Cohesiveness_Factor * SandClass%BedLoad(i, j)    
                            
                            endif
                            
                            SandClass%BedLoadX(i, j) =  SandClass%BedLoad(i, j) * cos(Wphi)
                            SandClass%BedLoadY(i, j) =  SandClass%BedLoad(i, j) * sin(Wphi) 
                                
                        endif
                        
                            Me%Bedload(i, j) = Me%Bedload(i, j) + SandClass%Bedload(i, j)                            
                    endif
                endif
            enddo
            enddo
        enddo do1
        
    end subroutine WavesOnlyBedload
    !--------------------------------------------------------------------------
    
    real function Cohesiveness_Adjust (Cohesive_Percentage, PM1, PM2) 
        
        !Arguments-------------------------------------------------------------
        real    :: Cohesive_Percentage, PM1, PM2
        
        !Local-----------------------------------------------------------------
        real    :: aux
        
        !Begin----------------------------------------------------------------
    
        if (Cohesive_Percentage > PM1) then
                                    
            if(Cohesive_Percentage > PM2) then
                !Fully cohesive bed - bedload is not possible 
                Cohesiveness_Adjust = 0.
            else
                !Cohesive bed - bedload decreases linearly with mud content                                        
                Cohesiveness_Adjust = (PM2 - Cohesive_Percentage)/(PM2 - PM1)
                
            endif
        endif
        
    end function Cohesiveness_Adjust
    !--------------------------------------------------------------------------

    subroutine ComputeFluxes
      
        !Local-----------------------------------------------------------------                
        integer                 :: i, j, n
        real(8)                 :: Xaux, Yaux
        class(T_Sand), pointer :: SandClass
        integer                :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)
            
            SandClass%FluxX  (:, :) =  0.
            SandClass%FluxY  (:, :) =  0.

            !Computes the Sediment fluxes (kg/s) in the middle of the cells

            do i=WILB, WIUB
            do j=WJLB, WJUB               

                if (Me%OpenSediment(i,j) == OpenPoint) then
                    
                    WKUB = Me%KTop(i, j)
                    
                    !if(SandClass%Field3D(i,j,WKUB) > 0.) then
                    if(SandClass%BedLoad(i, j) > 0.) then
                        
                        if (Me%BedloadMethod==1) then
                         
                            if (Me%ExternalVar%VelMod(i,j)> 0.) then                            
                           
                                Xaux = Me%ExternalVar%VelU(i,j) / Me%ExternalVar%VelMod(i,j)
                                Yaux = Me%ExternalVar%VelV(i,j) / Me%ExternalVar%VelMod(i,j)
                            
                                !kg/s
                                SandClass%FluxX(i, j) = SandClass%Bedload(i, j) * Xaux *        &
                                                        Me%ExternalVar%DVY(i, j)                        
                                                    
                                SandClass%FluxY(i, j) = SandClass%Bedload(i, j) * Yaux *        &
                                                        Me%ExternalVar%DUX(i, j)                                
                            endif
                            
                        elseif (Me%BedloadMethod==2 .or. Me%BedloadMethod==3) then
                            
                                !kg/s
                                SandClass%FluxX(i, j) = SandClass%BedloadX(i, j) * Me%ExternalVar%DVY(i, j)                        
                                                    
                                SandClass%FluxY(i, j) = SandClass%BedloadY(i, j) * Me%ExternalVar%DUX(i, j)
                        endif
                    endif
                endif               
            enddo
            enddo
        enddo
        
        if(Me%BedSlopeEffects)then
            call ComputeBedSlopeEffects
        endif
        
        call FluxesCorrection

    end subroutine ComputeFluxes
      
   !--------------------------------------------------------------------------

    subroutine ComputeBedSlopeEffects
      
        !Local-----------------------------------------------------------------                
        integer                 :: i, j, n
        real(8)                 :: Xaux, Yaux, AbsFlux
        class(T_Sand), pointer  :: SandClass
        integer                 :: WILB, WIUB, WJLB, WJUB
        real(8), parameter      :: PI_DBLE = 3.1415926536 !PI
        real                    :: alfa_bs, alfa_bn, phi, dhdx, dhdy
        real                    :: dzds, dzdn, alfa_s, alfa_n
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        alfa_bs = 1.0
        alfa_bn = 1.5
        
        !internal angle of friction of bed material (assumed to be 30º)
        phi = 30 * PI_DBLE/180.
        
        do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)

            do i=WILB, WIUB
            do j=WJLB, WJUB               

                if (Me%OpenSediment(i,j) == OpenPoint) then     
                        
                    if (SandClass%Bedload(i, j)> 0.) then
                            
                        AbsFlux = (SandClass%FluxX(i, j)**2 + SandClass%FluxY(i, j)**2)**0.5
                            
                        Xaux = SandClass%FluxX(i, j) / AbsFlux
                        Yaux = SandClass%FluxY(i, j) / AbsFlux
                            
                        dhdx = 0.
                        dhdy = 0.
                            
                        if (SandClass%FluxX(i, j) < 0.) then                            
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

                        if (SandClass%FluxY(i, j) < 0.) then
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
                            
                        alfa_n = alfa_bn * (SandClass%CriticalShearStress(i,j) /    &
                                Me%ExternalVar%ShearStress(i,j))**0.5 * dzdn
                                
                        !Adjustment of bedload transport for bed-slope effects
                        SandClass%FluxX(i, j) = alfa_s * (SandClass%FluxX(i, j) -   &
                                                alfa_n * SandClass%FluxY(i, j))
                            
                        SandClass%FluxY(i, j) = alfa_s * (SandClass%FluxY(i, j) +   &
                                                alfa_n * SandClass%FluxX(i, j))                            
                        
                    endif
                endif               
            enddo
            enddo
        enddo

    end subroutine ComputeBedSlopeEffects
      
   !--------------------------------------------------------------------------
    
    subroutine FluxesCorrection
              
        !Local-----------------------------------------------------------------                
        integer                 :: i, j, n
        class(T_Sand), pointer :: SandClass
        integer                :: WILB, WIUB, WJLB, WJUB, WKUB
        real(8)                :: correction, MassWithdrawal
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%FluxX  (:, :) =  0.
        Me%FluxY  (:, :) =  0.
        
        do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)

            do i=WILB, WIUB
            do j=WJLB, WJUB               

                if (Me%OpenSediment(i,j) == OpenPoint) then
                    
                    WKUB = Me%KTop(i, j)
                    
                    if(SandClass%Field3D(i,j,WKUB) > 0.) then                   
      
                        MassWithdrawal = (abs(SandClass%FluxX(i, j)) + abs(SandClass%FluxY(i, j))) * Me%Evolution%SedimentDT
                                
                        !Considers the flux to sediment (erosion, deposition)
                        MassWithdrawal = MassWithdrawal - SandClass%FluxToSediment(i,j) *   &
                                Me%Evolution%SedimentDT*Me%ExternalVar%GridCellArea(i,j)
                            
                        !Corrects the flux and bedload according to the mass content                             
                        If (MassWithdrawal > SandClass%Mass(i,j,WKUB)) then                            
                                
                            correction = SandClass%Mass(i,j,WKUB)/MassWithdrawal
                                
                            SandClass%FluxX(i, j) = SandClass%FluxX(i, j) * correction
                                
                            SandClass%FluxY(i, j) = SandClass%FluxY(i, j) * correction
                                
                            Me%Bedload(i, j) = Me%Bedload(i, j) - SandClass%Bedload(i, j) * (1 - correction)
                                
                            SandClass%BedLoad(i, j) = SandClass%BedLoad(i, j) * correction
                                
                        endif
                        
                        Me%FluxX(i,j) = Me%FluxX(i,j) + SandClass%FluxX(i, j)
                            
                        Me%FluxY(i,j) = Me%FluxY(i,j) + SandClass%FluxY(i, j)
                        
                    endif
                endif

            enddo
            enddo
        enddo

    end subroutine FluxesCorrection
    
    !--------------------------------------------------------------------------
      
    subroutine ComputeEvolution 
        
        !Local-----------------------------------------------------------------
        real                    :: RunPeriod
        integer                 :: i, j, n
        class(T_Sand), pointer :: SandClass      
        integer                 :: WILB, WIUB, WJLB, WJUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
          
        if (Me%Boxes%Yes) then

            Me%Boxes%FluxesX(:,:) = 0.
            Me%Boxes%FluxesY(:,:) = 0.

        endif
      
        Me%BedloadMass(:,:) = 0.
             
do1:    do n=1,Me%NumberOfClasses

           SandClass => Me%SandClass(n)
           
           SandClass%DM(:, :) = 0.
           
           if (Me%BedloadMethod /= 0) then
            
               if (Me%Residual%ON) then

                    RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime
                
                    Me%Residual%FluxX(:,:) = ( Me%Residual%FluxX(:,:) * (RunPeriod -  Me%Evolution%SedimentDT)          + &
                                               SandClass%FluxX(:,:) * Me%Evolution%SedimentDT) / RunPeriod
                                           
                    Me%Residual%FluxY(:,:) = ( Me%Residual%FluxY(:,:) * (RunPeriod -  Me%Evolution%SedimentDT)          + &
                                               SandClass%FluxY(:,:) * Me%Evolution%SedimentDT) / RunPeriod
                endif
           
                do j=WJLB, WJUB
                do i=WILB, WIUB
                    
                    if (SandClass%FluxX(i, j) < 0.) then 
                        if      (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then
                        
                            SandClass%DM(i, j-1) = SandClass%DM(i, j-1) - Me%Evolution%SedimentDT * SandClass%FluxX(i, j)
                            SandClass%DM(i, j  ) = SandClass%DM(i, j  ) + Me%Evolution%SedimentDT * SandClass%FluxX(i, j)
                         

                            if (Me%Boxes%Yes) then
                                Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j) + SandClass%FluxX(i, j) * Me%Evolution%SedimentDT
                            endif
                        endif
                    else 

                        if (Me%ExternalVar%ComputeFacesU2D(i,j+1) == Covered) then
                  
                            SandClass%DM(i, j+1) = SandClass%DM(i, j+1) + Me%Evolution%SedimentDT * SandClass%FluxX(i, j)  
                            SandClass%DM(i, j  ) = SandClass%DM(i, j  ) - Me%Evolution%SedimentDT * SandClass%FluxX(i, j) 

                            if (Me%Boxes%Yes) then
                                Me%Boxes%FluxesX(i,j+1) = Me%Boxes%FluxesX(i,j+1) + SandClass%FluxX(i, j) * Me%Evolution%SedimentDT
                            endif
                        endif                    
                    endif

                    if (SandClass%FluxY(i, j) < 0.) then
                        if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then
                            
                            SandClass%DM(i-1, j) = SandClass%DM(i-1, j) - Me%Evolution%SedimentDT * SandClass%FluxY(i, j) 
                            SandClass%DM(i  , j) = SandClass%DM(i  , j) + Me%Evolution%SedimentDT * SandClass%FluxY(i, j) 

                            if (Me%Boxes%Yes) then
                                Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j) + SandClass%FluxY(i, j) * Me%Evolution%SedimentDT
                            endif
                        endif
                    else 
                        if (Me%ExternalVar%ComputeFacesV2D(i+1, j) == Covered) then

                            SandClass%DM(i+1, j) = SandClass%DM(i+1, j) + Me%Evolution%SedimentDT * SandClass%FluxY(i, j)
                            SandClass%DM(i  , j) = SandClass%DM(i  , j) - Me%Evolution%SedimentDT * SandClass%FluxY(i, j) 

                            if (Me%Boxes%Yes) then
                                Me%Boxes%FluxesY(i+1,j) = Me%Boxes%FluxesY(i+1,j) + SandClass%FluxY(i, j) * Me%Evolution%SedimentDT
                            endif
                        endif
                    endif
                enddo
                enddo
            endif
           
            if (Me%TimeSerie)then
                !Only for output
                Me%BedloadMass(:,:) = Me%BedloadMass(:,:) + SandClass%DM(:,:)
            endif
            
        enddo do1

        

   
    end subroutine ComputeEvolution
    
    !--------------------------------------------------------------------------

    subroutine ComputeDischarges            

        !Local-----------------------------------------------------------------
        real                               :: MassAdd,      &
                                              DischargeFlow, DischargeConc
        integer                            :: i, j, dis, DischargesNumber, STAT_CALL
        integer                            :: n
        class(T_Sand), pointer            :: SandClass 
        !----------------------------------------------------------------------


            !Sinks and Sources
            call GetDischargesNumber(Me%ObjDischarges, DischargesNumber, STAT=STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ComputeDischarges - ModuleSediment - ERR10'


            !For all Discharges
d1:         do dis = 1, DischargesNumber
            

                call GetDischargesGridLocalization(Me%ObjDischarges,                    &
                                                   dis, Igrid = I, JGrid = J,           &
                                                   STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ComputeDischarges - ModuleSediment - ERR20'

                call GetDischargeWaterFlow(Me%ObjDischarges,                            &
                                           Me%ExternalVar%Now,                          &
                                           dis, -99., DischargeFlow,                    &
                                           STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'ComputeDischarges - ModuleSediment - ERR30'


do1:            do n=1,Me%NumberOfClasses

                    SandClass => Me%SandClass(n)
                    
                    call GetDischargeConcentration (Me%ObjDischarges,                       &
                                                    Me%ExternalVar%Now,                     &
                                                    dis, DischargeConc,                     &
                                                    PropertyIDNumber = SandClass%ID%IDNumber,         &        
                                                    STAT = STAT_CALL)
                    if (STAT_CALL/=SUCCESS_)                                                &
                        stop 'ComputeDischarges - ModuleSediment - ERR40'
                
                    !WKUB = Me%KTop(i, j)
                    !
                    !BottomDensity       = Me%Porosity(i,j,WKUB)*Me%ExternalVar%WaterDensity + (1. - Me%Porosity(i,j,WKUB))*Me%Density
                    !
                    !DT_Area             = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)
                    !
                    !                     ![m3/s] * [g/m3/1000] / [kg/m3] * [s/m2] * []
                    !TicknessAdd         = DischargeFlow * (DischargeConc/1000.) / BottomDensity * DT_Area
                
                                        ![m3/s] * [kg/m3] * [s]
                    MassAdd            = DischargeFlow * DischargeConc * Me%Evolution%SedimentDT

                    SandClass%DM(i, j) = SandClass%DM(i, j) + MassAdd
                    
                enddo do1

            enddo d1
        
    end subroutine ComputeDischarges            

    
    !--------------------------------------------------------------------------


    subroutine BoundaryCondition

        !Local-----------------------------------------------------------------
        real                               :: a1, a2, a3, a4, atotal
        integer                            :: i, j, ILB, IUB, JLB, JUB
        integer                            :: n
        class(T_Sand), pointer            :: SandClass 
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 

        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        
do1:    do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)


            if  (Me%Boundary == NullGradient) then

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
                    
                            SandClass%DM(i, j) = (a1 * SandClass%DM(i-1, j) + a2 * SandClass%DM(i+1, j)  +     &
                                             a3 * SandClass%DM(i, j-1) + a4 * SandClass%DM(i, j+1)) / atotal
                        endif
                                          

                    endif

                enddo
                enddo    

            else if (Me%Boundary == Cyclic) then

                where (Me%ExternalVar%BoundaryPoints2D(ILB, :) == Boundary) 
                    SandClass%DM(ILB, :) = SandClass%DM(IUB-1, :)
                end where

                where (Me%ExternalVar%BoundaryPoints2D(IUB, :) == Boundary) 
                    SandClass%DM(IUB, :) = SandClass%DM(ILB+1, :)
                end where

                where (Me%ExternalVar%BoundaryPoints2D(: ,JLB) == Boundary) 
                    SandClass%DM(:, JLB) = SandClass%DM(:, JUB-1)
                end where

                where (Me%ExternalVar%BoundaryPoints2D(:, JUB) == Boundary) 
                    SandClass%DM(:, JUB) = SandClass%DM(:, JLB+1)
                end where
                
            else if (Me%Boundary == NullValue) then

                where (Me%ExternalVar%BoundaryPoints2D(:, :) == Boundary) 
                    SandClass%DM(:, :) = 0.
                end where            

            endif
        enddo do1        
    
    end subroutine BoundaryCondition
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ComputeMass
    
        !Local-----------------------------------------------------------------
        integer                     :: i, j, n
        class(T_Sand), pointer      :: SandClass
        integer                     :: WILB, WIUB, WJLB, WJUB, WKUB
        real(8)                     :: aux
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%Mass(:,:) = 0.
        Me%DM(:,:) = 0.
        Me%FluxToSediment(:,:) = 0.
        
if1:    if (Me%CohesiveClass%Run) then
        
do1:        do j=WJLB, WJUB
do2:        do i=WILB, WIUB
              
if2:            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
                    
                    WKUB = Me%KTop(i, j)
                    
                    Me%CohesiveClass%DM(i, j) = Me%ExternalVar%ConsolidationFlux(i, j)* &
                                            Me%Evolution%SedimentDT*Me%ExternalVar%GridCellArea(i,j)
                 
                    !kg/m2
                    aux = (Me%CohesiveClass%Mass(i,j,WKUB) + Me%CohesiveClass%DM(i, j))/Me%ExternalVar%GridCellArea(i,j)
                    
if3:                if (aux < 1.e-2) then
                        
                        Me%CohesiveClass%DM(i, j) = - Me%CohesiveClass%Mass(i,j,WKUB)
                        
                        Me%CohesiveClass%Mass(i,j,WKUB) = 0.
                        
                        Me%CohesiveClass%Field3D(i,j,WKUB) = 0.
                        
                        Me%CohesiveClass%TopPercentage(i,j) = 0.
                        
                    else
                        
                        Me%CohesiveClass%Mass(i,j,WKUB) = Me%CohesiveClass%Mass(i,j,WKUB) + Me%CohesiveClass%DM(i, j)
                        
                    endif if3
 
                    Me%Mass(i,j) = Me%CohesiveClass%Mass(i,j,WKUB)
                    
                    Me%DM(i,j) = Me%CohesiveClass%DM(i, j)
                
                endif if2
            enddo do2
            enddo do1
        endif if1
                 
do3:    do j=WJLB, WJUB
do4:    do i=WILB, WIUB
             
if4:        if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then  
                
                WKUB = Me%KTop(i, j)
                
do5:            do n=1,Me%NumberOfClasses

                    SandClass => Me%SandClass(n)
                    
                    SandClass%DM(i, j) = SandClass%DM(i, j) + SandClass%FluxToSediment(i,j) *   &
                                        Me%Evolution%SedimentDT*Me%ExternalVar%GridCellArea(i,j)
                    
                    if (Me%TimeSerie)then
                        !kg/s/m2
                        Me%FluxToSediment(i,j) = Me%FluxToSediment(i,j) + SandClass%FluxToSediment(i,j)
                    
                        !kg
                        Me%TotalFluxToSediment(i,j) = Me%FluxToSediment(i,j)    *       &
                                                      Me%Evolution%SedimentDT   *       &
                                                      Me%ExternalVar%GridCellArea(i,j)
                    endif
                    
                    !kg/m2                    
                    aux = (SandClass%Mass(i,j,WKUB) + SandClass%DM(i, j))/Me%ExternalVar%GridCellArea(i,j)  
                    
if5:                if (aux < 1.e-2) then
                        
                        SandClass%DM(i, j) = - SandClass%Mass(i,j,WKUB)
                        
                        SandClass%Mass(i,j,WKUB) = 0.
                        
                        SandClass%Field3D(i, j, WKUB) = 0.
                        
                        SandClass%TopPercentage(i,j) = 0.
                        
                    else
                    
                        SandClass%Mass(i,j,WKUB) = SandClass%Mass(i,j,WKUB) + SandClass%DM(i, j)
                        
                    endif if5
                     
                    Me%Mass(i,j) = Me%Mass(i,j) + SandClass%Mass(i,j,WKUB)
                    
                    Me%DM(i,j) = Me%DM(i,j) + SandClass%DM(i, j)
                 
                enddo do5              
                
            endif if4
          
        enddo do4
        enddo do3        
       
    end subroutine ComputeMass
    
    !--------------------------------------------------------------------------
    
    subroutine ComputePercentage
        !Local-----------------------------------------------------------------
        integer                            :: i, j, n, WKUB
        integer                            :: WILB, WIUB, WJLB, WJUB
        class(T_Sand), pointer             :: SandClass
         !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
          do j=WJLB, WJUB
          do i=WILB, WIUB
               
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%Mass(i,j) > 0) then
                
                WKUB = Me%KTop(i, j)
                
                 Me%TotalPercentage (i,j,WKUB) = 0.
                
                if (Me%CohesiveClass%Run) then
                    
                        Me%CohesiveClass%Field3D(i,j,WKUB) = Me%CohesiveClass%Mass(i,j,WKUB) / Me%Mass(i, j)
                    
                        Me%CohesiveClass%TopPercentage(i,j) = Me%CohesiveClass%Field3D(i,j,WKUB)
                    
                        Me%TotalPercentage (i,j,WKUB) = Me%CohesiveClass%Field3D(i,j,WKUB)
                        
                        if(Me%CohesiveClass%Field3D(i,j,WKUB) >  0.) then
                            call ComputePorosity (i,j,WKUB)        
                        else
                            Me%Porosity(i,j,WKUB) = Me%PorositySand
                        endif
                            
                endif
                
                do n=1,Me%NumberOfClasses

                    SandClass => Me%SandClass(n)
                
                        SandClass%Field3D(i, j, WKUB) =  SandClass%Mass(i,j,WKUB) / Me%Mass(i, j)
                        
                        SandClass%TopPercentage(i,j) = SandClass%Field3D(i,j,WKUB)
               
                        Me%TotalPercentage (i,j,WKUB) = Me%TotalPercentage (i,j,WKUB) + SandClass%Field3D(i,j,WKUB)
                          
                    if  (SandClass%Field3D(i, j, WKUB)  > 1.001) then
                             write(*,*) 'The class percentage is larger than 100%.'
                             write(*,*) i, j, WKUB, SandClass%Field3D(i, j, WKUB), 'n=',n
                             stop 'ComputePercentage - ModuleSediment - ERR01'                 
                    elseif (SandClass%Field3D(i, j, WKUB)  < 0.) then
                             write(*,*) 'The class percentage is smaller than 0%.'
                             write(*,*) i, j, WKUB, SandClass%Field3D(i, j, WKUB), 'n=',n
                             stop 'ComputePercentage - ModuleSediment - ERR02'
                    endif         
                enddo
            endif
                    
          enddo
          enddo                               
   
        do j=WJLB, WJUB
        do i=WILB, WIUB
            
            WKUB = Me%KTop(i, j)
              
              if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint .and. Me%Mass(i,j) > 0) then          
                  
                    if  (Me%TotalPercentage (i, j, WKUB)  > 1.001) then
                         write(*,*) 'The sum of the classes percentage is larger than 100%.'
                         write(*,*) i, j, WKUB, Me%TotalPercentage (i, j, WKUB)
                         stop 'ComputePercentage - ModuleSediment - ERR03'                 
                    elseif (Me%TotalPercentage (i, j, WKUB)  < 0.999) then
                         write(*,*) 'The sum of the classes percentage is smaller than 100%.'
                         write(*,*) i, j, WKUB, Me%TotalPercentage (i, j, WKUB)
                         stop 'ComputePercentage - ModuleSediment - ERR04'
                    endif
              endif
              
            enddo
            enddo
        
    end subroutine ComputePercentage

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !Compute the porosity of a sand-mud mixture considering the fractional packing model 
    !from Koltermann and Gorelick (1995)
    subroutine ComputePorosity (i, j, k)
    
        !Local-----------------------------------------------------------------
        integer                     :: i, j, k
        real                        :: ymin, y, c
        !----------------------------------------------------------------------
        
        ymin = 0.76
        
        Me%CohesiveClass%Porosity(i,j,k) =  1 - Me%CohesiveDryDensity%Field3D(i,j,k) / Me%Density
        
        if(Me%SandD50(i,j) > 0)then
        
            !Volume fraction of the cohesive class based on total volume        
            c = Me%CohesiveClass%Field3D(i,j,k) * (1 - Me%Porosity (i,j,k))
        
            !WKUB = Me%KTop(i, j)
            !
            !if (k == WKUB) then
            !
            !    Me%CohesiveClass%Porosity(i,j,WKUB) = (Me%CohesiveClass%Porosity(i,j,WKUB) * 
            !                                          (Me%CohesiveClass%Mass(i, j) - Me%CohesiveClass%DM(i, j))      &
            !                                + Me%MaxPorosityCohesive * Me%CohesiveClass%DM(i, j)) / Me%CohesiveClass%Mass(i, j)
            !endif
        
                
            if (c < Me%PorositySand) then
                    
                y = c * (ymin - 1) / Me%PorositySand + 1
                
                Me%Porosity(i,j,k) = Me%PorositySand - c * y * (1 - Me%CohesiveClass%Porosity(i,j,k)) +     & 
                                            (1 - y) * c * Me%CohesiveClass%Porosity(i,j,k)
                
            else
                    
                y = (c - 1) * (1 - ymin)/(1 - Me%PorositySand) + 1
                    
                Me%Porosity(i,j,k) = Me%PorositySand * (1 - y) + c * Me%CohesiveClass%Porosity(i,j,k)
                
            endif
                
        else
            
            Me%Porosity(i,j,k) = Me%CohesiveClass%Porosity(i,j,k)
                
        endif        
       
    end subroutine ComputePorosity
    
    !--------------------------------------------------------------------------
     
    !--------------------------------------------------------------------------
    subroutine ComputeTotalDZ
    
        !Local-----------------------------------------------------------------
        integer                     :: i, j
        integer                     :: WILB, WIUB, WJLB, WJUB, WKUB
        real                        :: Area, factor
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%DZ(:,:) = 0.
        
        do j=WJLB, WJUB
        do i=WILB, WIUB
              
            if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
                
                WKUB = Me%KTop(i, j)
                
                Area = Me%ExternalVar%DUX(i, j) * Me%ExternalVar%DVY(i, j)
                
                !kg/m
                factor = Me%density*Area*(1-Me%Porosity (i,j,WKUB))
                
                !m    
                Me%DZ(i, j) = Me%DM(i, j) / factor
                
            endif
          
        enddo
        enddo          
       
    end subroutine ComputeTotalDZ
    
    !--------------------------------------------------------------------------
    
    subroutine ComputeVerticalCoordinate          

        !Local-----------------------------------------------------------------
        integer                             :: i, j, WKUB, TotalWKUB
        real(8)                             :: DZ, TopLayerThickness
        real(8)                             :: ExcessDZ
        integer                             :: WILB, WIUB, WJLB, WJUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        TotalWKUB = Me%SedimentWorkSize3D%KUB
           
        do j = WJLB, WJUB
        do i = WILB, WIUB

if1:        if (Me%ExternalVar%OpenPoints2D (i ,j) == OpenPoint) then
                    
                WKUB = Me%KTop(i, j)
                    
                DZ = Me%DZ(i, j)
                    
                TopLayerThickness = Me%ExternalVar%DWZ(i,j,WKUB)
                    
if2:            if (DZ > 0.)then !aggradation

                    ExcessDZ = TopLayerThickness + DZ - Me%MaxLayerThickness

if3:                if(ExcessDZ > 0.)then !New Layer
    
                        call New_Layer(i,j,WKUB,ExcessDZ)     
                        
                    else !Increase layer thickness
                        
                        if (Me%KTop(i,j) < 1) then !Activate layer                        
                            
                            call Activate_Layer (i,j,WKUB)
                            
                            WKUB = Me%KTop(i, j)
                                                     
                        endif
                        
                        Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB) - DZ
                            
                    endif if3                 

                elseif(DZ <  0.)then !degradation

                    ExcessDZ = Me%ExternalVar%DWZ(i,j,WKUB) - Me%MinLayerThickness

if7:                if(abs(DZ) .ge. ExcessDZ .or. Me%Mass(i,j) == 0.) then !Layer eroded
    
                        call Layer_Eroded (i,j,WKUB,DZ,ExcessDZ)
                        
                    else

                        Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB) - DZ

                    endif if7

                end if if2

            end if if1

        enddo
        enddo    
    
    end subroutine ComputeVerticalCoordinate      
    
    !--------------------------------------------------------------------------
    
    subroutine New_Layer(i,j,WKUB,ExcessDZ)
    
        !Arguments-------------------------------------------------------------
        integer                             :: i, j, WKUB
        real(8)                             :: ExcessDZ
    
        !Local-----------------------------------------------------------------
        integer                             :: n, TotalWKUB
        class(T_Sand), pointer              :: SandClass
        real(8)                             :: ExcessMin
        !----------------------------------------------------------------------
    
        TotalWKUB = Me%SedimentWorkSize3D%KUB
    
        Me%KTop(i,j) = WKUB+1
                        
        ExcessMin = Me%MinLayerThickness - ExcessDZ
 
                        
                            
        Me%VerticalCoordinate(i,j,WKUB)   = Me%VerticalCoordinate(i,j,WKUB-1)     - &
                                            Me%MaxLayerThickness + ExcessMin
            
        Me%VerticalCoordinate(i,j,WKUB+1) = Me%VerticalCoordinate(i,j,WKUB)       - &
                                            Me%MinLayerThickness
             
        !The new top layer has initially the same sediment properties than the older one 
        
        Me%Porosity(i,j,WKUB+1) = Me%Porosity(i,j,WKUB)
        
        do n=1,Me%NumberOfClasses 
                            
            SandClass => Me%SandClass(n)
            
            SandClass%Mass(i,j,WKUB+1) = SandClass%Mass(i,j,WKUB) *                     &
            (Me%VerticalCoordinate(i,j,WKUB) - Me%VerticalCoordinate(i,j,WKUB+1)) /       &
            (Me%VerticalCoordinate(i,j,WKUB-1) - Me%VerticalCoordinate(i,j,WKUB+1)) 
                                
            SandClass%Mass(i,j,WKUB) = SandClass%Mass(i,j,WKUB) *                       &
            (Me%VerticalCoordinate(i,j,WKUB-1) - Me%VerticalCoordinate(i,j,WKUB)) /     &
            (Me%VerticalCoordinate(i,j,WKUB-1) - Me%VerticalCoordinate(i,j,WKUB+1))
            
            SandClass%Field3D(i,j,WKUB+1) = SandClass%Field3D(i,j,WKUB)
                         
        enddo
                        
if4:    if(Me%CohesiveClass%Run) then
                
            Me%CohesiveClass%Mass(i,j,WKUB+1) = Me%CohesiveClass%Mass(i,j,WKUB) *                 &
            (Me%VerticalCoordinate(i,j,WKUB) - Me%VerticalCoordinate(i,j,WKUB+1)) /       &
            (Me%VerticalCoordinate(i,j,WKUB-1) - Me%VerticalCoordinate(i,j,WKUB+1)) 
            
            Me%CohesiveClass%Mass(i,j,WKUB) = Me%CohesiveClass%Mass(i,j,WKUB) *                   &
            (Me%VerticalCoordinate(i,j,WKUB-1) - Me%VerticalCoordinate(i,j,WKUB)) /     &
            (Me%VerticalCoordinate(i,j,WKUB-1) - Me%VerticalCoordinate(i,j,WKUB+1))

            Me%CohesiveDryDensity%Field3D(i,j,WKUB+1) = Me%CohesiveDryDensity%Field3D(i,j,WKUB)
            
            Me%CohesiveClass%Porosity(i,j,WKUB+1) = Me%CohesiveClass%Porosity(i,j,WKUB)
            
            Me%CohesiveClass%Field3D(i,j,WKUB+1) = Me%CohesiveClass%Field3D(i,j,WKUB)
                            
        endif if4

if5:    if(Me%KTop(i,j) > TotalWKUB)then  !Maximum number of layers exceeded
    
            call MaxLayersExceeded (i,j)
                            
            !write(*,*) 'Maximum number of layers exceeded in cell i=',i,'j=',j
            !write(*,*) 'Last sediment layer deleted'
            !write(*,*) 'ComputeVerticalCoordinate - ModuleSediment - WRN10' 
                                                                                                                                  
        else
    
            Me%VerticalCoordinate(i,j,WKUB+1:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB)       - &
                                                          Me%MinLayerThickness

        endif if5
        
    end subroutine New_Layer
                        
    !--------------------------------------------------------------------------
    
    subroutine MaxLayersExceeded(i,j)
    
        !Local-----------------------------------------------------------------
        integer                             :: i, j, k, n, TotalWKUB
        class(T_Sand), pointer              :: SandClass
        real(8)                             :: Mass
        !----------------------------------------------------------------------
        
        TotalWKUB = Me%SedimentWorkSize3D%KUB
        
        Me%KTop(i,j) = TotalWKUB
        
        k = 1 
        
        Me%TotalPercentage (i,j,k) = 0.
        
        Mass = 0.
        
        !Last sediment layer merged to maintain the number of layers constant 
                                    
        Me%VerticalCoordinate(i,j,k) = Me%VerticalCoordinate(i,j,k+1)
            
        !Sediment properties in the Last sediment layer averaged 
            
        Me%Porosity(i,j,k) = (Me%Porosity(i,j,k) * Me%ExternalVar%DWZ(i,j,k) + &
                              Me%Porosity(i,j,k+1) * Me%ExternalVar%DWZ(i,j,k+1))/  &
                              (Me%ExternalVar%DWZ(i,j,k) + Me%ExternalVar%DWZ(i,j,k+1))
            
        do n=1,Me%NumberOfClasses 
                                   
            SandClass => Me%SandClass(n)
                                            
            SandClass%Mass(i,j,k) = SandClass%Mass(i,j,k) + SandClass%Mass(i,j,k+1)
            
            Mass = Mass + SandClass%Mass(i,j,k)
        enddo
                                
if1:    if(Me%CohesiveClass%Run) then    
                
            Me%CohesiveClass%Porosity(i,j,k) = (Me%CohesiveClass%Porosity(i,j,k) * Me%CohesiveClass%Mass(i,j,k) +   &
                                                Me%CohesiveClass%Porosity(i,j,k+1) * Me%CohesiveClass%Mass(i,j,k+1))/    &
                                                (Me%CohesiveClass%Mass(i,j,k) + Me%CohesiveClass%Mass(i,j,k+1))
                
            Me%CohesiveDryDensity%Field3D(i,j,k) = (Me%CohesiveDryDensity%Field3D(i,j,k) * Me%CohesiveClass%Mass(i,j,k) +   &
                                                    Me%CohesiveDryDensity%Field3D(i,j,k+1) * Me%CohesiveClass%Mass(i,j,k+1))/   &
                                                    (Me%CohesiveClass%Mass(i,j,k) + Me%CohesiveClass%Mass(i,j,k+1))
                
            Me%CohesiveClass%Mass(i,j,k) = Me%CohesiveClass%Mass(i,j,k) + Me%CohesiveClass%Mass(i,j,k+1)
                
            Mass = Mass + Me%CohesiveClass%Mass(i,j,k)
                
            Me%CohesiveClass%Field3D(i,j,k) = Me%CohesiveClass%Mass(i,j,k) / Mass
                
            Me%TotalPercentage (i,j,k) = Me%CohesiveClass%Field3D(i,j,k)
                            
        endif if1


            do n=1,Me%NumberOfClasses 
                                   
                SandClass => Me%SandClass(n)
                
                SandClass%Field3D(i,j,k) = SandClass%Mass(i,j,k) / Mass

                Me%TotalPercentage (i,j,k) = Me%TotalPercentage (i,j,k) + SandClass%Field3D(i,j,k)
                
            enddo
                
            if (Me%TotalPercentage (i, j, k)  > 1.001) then
                write(*,*) 'The sum of the classes percentage is larger than 100%.'
                write(*,*) i, j, k, Me%TotalPercentage (i, j, k)
                stop 'MaxLayersExceeded - ModuleSediment - ERR10'
            elseif (Me%TotalPercentage (i, j, k)  < 0.999) then
                    write(*,*) 'The sum of the classes percentage is smaller than 100%.'
                    write(*,*) i, j, k, Me%TotalPercentage (i, j, k)
                    stop 'MaxLayersExceeded - ModuleSediment - ERR20'
            endif

                                       
        do k=2,TotalWKUB
                                    
            Me%VerticalCoordinate(i,j,k) = Me%VerticalCoordinate(i,j,k+1)
            
            !Sediment properties in the cells updated 
            
            Me%Porosity(i,j,k) = Me%Porosity(i,j,k+1)
            
            do n=1,Me%NumberOfClasses 
                                   
                SandClass => Me%SandClass(n)
                                    
                SandClass%Mass(i,j,k) = SandClass%Mass(i,j,k+1)
                
                SandClass%Field3D(i,j,k) = SandClass%Field3D(i,j,k+1)
                
            enddo
                                
if6:        if(Me%CohesiveClass%Run) then
                            
                Me%CohesiveClass%Mass(i,j,k) = Me%CohesiveClass%Mass(i,j,k+1)
                
                Me%CohesiveClass%Field3D(i,j,k) = Me%CohesiveClass%Field3D(i,j,k+1)
                
                Me%CohesiveClass%Porosity(i,j,k) = Me%CohesiveClass%Porosity(i,j,k+1)
                
                Me%CohesiveDryDensity%Field3D(i,j,k) = Me%CohesiveDryDensity%Field3D(i,j,k+1)
                            
            endif if6
                            
        enddo 
            
    end subroutine MaxLayersExceeded
            
    !--------------------------------------------------------------------------
            
    subroutine Activate_Layer (i,j,WKUB)
    
        !Local-----------------------------------------------------------------
        integer                             :: i, j, n, WKUB
        class(T_Sand), pointer              :: SandClass
        real(8)                             :: Mass
        !----------------------------------------------------------------------
    
        Mass = 0.
        
        Me%TotalPercentage (i,j,WKUB+1) = 0.
        
        do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)   
            
            !Compute values of the activated layer                 
            SandClass%Mass(i, j, WKUB+1) = SandClass%Mass(i, j, WKUB)
            
            Mass = Mass + SandClass%Mass(i,j,WKUB+1)
            
            !Nullify layer 0             
            SandClass%Mass(i, j, WKUB) = 0.
            SandClass%Field3D(i,j,WKUB) = 0.
            
        enddo
                            
        if(Me%CohesiveClass%Run) then
            
            !Compute values of the activated layer                               
            Me%CohesiveClass%Mass(i,j,WKUB+1) = Me%CohesiveClass%Mass(i,j,WKUB)
            
            Mass = Mass + Me%CohesiveClass%Mass(i,j,WKUB+1)
            
            Me%CohesiveClass%Field3D(i,j,WKUB+1) = Me%CohesiveClass%Mass(i,j,WKUB+1)/Mass 
            
            Me%TotalPercentage (i,j,WKUB+1) = Me%CohesiveClass%Field3D(i,j,WKUB+1)
            
            Me%CohesiveDryDensity%Field3D(i,j,WKUB+1) = Me%CohesiveDryDensity%Min
            
            Me%CohesiveClass%Porosity(i,j,WKUB+1) = 1 - Me%CohesiveDryDensity%Field3D(i,j,WKUB) / Me%Density
            
            !Nullify layer 0            
            Me%CohesiveClass%Mass(i,j,WKUB) = 0.
            Me%CohesiveClass%Field3D(i,j,WKUB) = 0.
            
        endif
        
        do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)        
            
            SandClass%Field3D(i,j,WKUB+1) = SandClass%Mass(i,j,WKUB+1)/Mass
            
            Me%TotalPercentage (i,j,WKUB+1) = Me%TotalPercentage (i,j,WKUB+1) + SandClass%Field3D(i,j,WKUB+1)
            
        enddo
        
        if (Me%TotalPercentage (i,j,WKUB+1)  > 1.001) then
            write(*,*) 'The sum of the classes percentage is larger than 100%.'
            write(*,*) i, j, WKUB+1, Me%TotalPercentage (i, j, WKUB+1)
            stop 'Activate_Layer - ModuleSediment - ERR10'
        elseif (Me%TotalPercentage (i,j,WKUB+1)  < 0.999) then
            write(*,*) 'The sum of the classes percentage is smaller than 100%.'
            write(*,*) i, j, WKUB+1, Me%TotalPercentage (i, j, WKUB+1)
            stop 'Activate_Layer - ModuleSediment - ERR20'
        endif
        
        Me%KTop(i,j) = 1  
        
        !Me%OpenSediment(i,j) = 1
        
    end subroutine Activate_Layer
        
    !--------------------------------------------------------------------------
    subroutine Layer_Eroded (i,j,WKUB,DZ,ExcessDZ)
    
        !Local-----------------------------------------------------------------
        integer                             :: i, j, n, WKUB,TotalWKUB
        class(T_Sand), pointer              :: SandClass
        real(8)                             :: DZ, ExcessDZ
        !----------------------------------------------------------------------
    
        TotalWKUB = Me%SedimentWorkSize3D%KUB
        
        
if8:    if (WKUB > 1) then
                    
            Me%VerticalCoordinate(i,j,WKUB-1) = Me%VerticalCoordinate(i,j,WKUB-1)       &
                                            - Me%ExternalVar%DWZ(i,j,WKUB) - DZ
                                                            
            Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB-1)

                        
            do n=1,Me%NumberOfClasses

                SandClass => Me%SandClass(n)        
                
                SandClass%Mass(i,j,WKUB-1) = SandClass%Mass(i,j,WKUB-1) + SandClass%Mass(i,j,WKUB)
                
                SandClass%Mass(i,j,WKUB) = 0.
                
                SandClass%Field3D(i,j,WKUB) = null_real
                
                SandClass%TopPercentage(i,j) = 0.
                
            enddo
                        
if9:        if (Me%CohesiveClass%Run) then
                
                Me%CohesiveClass%Mass(i,j,WKUB-1) = Me%CohesiveClass%Mass(i,j,WKUB-1) + Me%CohesiveClass%Mass(i,j,WKUB)
                
                Me%CohesiveClass%Mass(i,j,WKUB) = 0.
                
                Me%CohesiveClass%Field3D(i,j,WKUB) = null_real
                
                Me%CohesiveClass%TopPercentage(i,j) = 0.
                
                call ComputePorosity (i,j,WKUB-1)
                            
            endif if9
                       
        else !Eroded all sediment layers                         
                        
            Me%VerticalCoordinate(i,j,WKUB) =  Me%VerticalCoordinate(i,j,WKUB-1)       &
                                                    - Me%MinLayerThickness
                        
            Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB)
                        
            Me%DZ(i, j) = -ExcessDZ 
            
            do n=1,Me%NumberOfClasses

                SandClass => Me%SandClass(n)        
                
                SandClass%Mass(i,j,WKUB) = 0.
                
                SandClass%Field3D(i,j,WKUB) = null_real
                
                SandClass%TopPercentage(i,j) = 0.
                
            enddo
            
            if (Me%CohesiveClass%Run) then
                
                Me%CohesiveClass%Mass(i,j,WKUB) = 0.
                
                Me%CohesiveClass%Field3D(i,j,WKUB) = null_real
                
                Me%CohesiveClass%TopPercentage(i,j) = 0.
                
                Me%Porosity(i,j,WKUB) = Me%PorositySand
                
                Me%CohesiveDryDensity%Field3D(i,j,WKUB) = Me%CohesiveDryDensity%Min
                
                !To avoid errors when the layer is reactivated
                
                Me%Porosity(i,j,WKUB-1) = Me%PorositySand
                
                Me%CohesiveDryDensity%Field3D(i,j,WKUB-1) = Me%CohesiveDryDensity%Min
                            
            endif
                            
        endif if8
                    
        Me%KTop(i,j)      = WKUB-1
                        
    end subroutine Layer_Eroded
    
    !--------------------------------------------------------------------------

    subroutine New_Geometry

        !External----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer, pointer, dimension(:,:,:)  :: WaterPoints3D
        
        !Local-------------------------------------------------------------------
        type(T_Time)                        :: ActualTime

        !------------------------------------------------------------------------

        call ReadUnLockExternalVar  
        
        call GetWaterPoints3D(Me%ObjSedimentMap, WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError(FATAL_, INTERNAL_, "New_Geometry; ModuleSediment. ERR01")
        
        ActualTime = Me%ExternalVar%Now

        !Compute new volume 
        call ComputeVerticalGeometry(Me%ObjSedimentGeometry,                             &
                                     WaterPoints3D      = WaterPoints3D,                 &
                                     ActualTime         = ActualTime,                    &
                                     SZZ                = Me%VerticalCoordinate,         &
                                     KTOP               = Me%KTop,                       &
                                     STAT               = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError(FATAL_, INTERNAL_, "New_Geometry; ModuleSediment. ERR02")
        
        call UnGetMap(Me%ObjSedimentMap, WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                                     &
            call SetError(FATAL_, INTERNAL_, "New_Geometry; ModuleSediment. ERR03") 

        
        call ReadLockExternalVar  
                   

    end subroutine New_Geometry

    !----------------------------------------------------------------------------

    
    subroutine ComputeResidualEvolution            

        !Local-----------------------------------------------------------------
        integer                 :: i, j
        integer                 :: WILB, WIUB, WJLB, WJUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
           
            do j=WJLB, WJUB
            do i=WILB, WIUB

                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                    
                   Me%DZ_Residual(i, j) = Me%DZ_Residual(i, j) + Me%DZ(i, j)

                endif

            enddo
            enddo    

            if (Me%Evolution%Bathym) then

                do j=WJLB, WJUB
                do i=WILB, WIUB

                    if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                    
                       Me%BatimIncrement(i,j) = Me%BatimIncrement(i,j) + Me%DZ(i, j) * Me%MorphologicalFactor

                    endif

                enddo
                enddo    

            endif
                    
    end subroutine ComputeResidualEvolution
    
    !--------------------------------------------------------------------------


    subroutine OutPutSedimentHDF
        
        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
         
        !Local-----------------------------------------------------------------
        logical                            :: FirstTime
        integer                            :: OutPutNumber
        type (T_Time)                      :: Actual
        real,    dimension(6    ), target  :: AuxTime
        real,    dimension(:    ), pointer :: TimePtr
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                            :: WorkKLB, WorkKUB
        integer                            :: n

        !----------------------------------------------------------------------

        WorkILB = Me%SedimentWorkSize3D%ILB 
        WorkIUB = Me%SedimentWorkSize3D%IUB 
        WorkJLB = Me%SedimentWorkSize3D%JLB 
        WorkJUB = Me%SedimentWorkSize3D%JUB 
        WorkKLB = Me%SedimentWorkSize3D%KLB 
        WorkKUB = Me%SedimentWorkSize3D%KUB 

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
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR010'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",     &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR020'
                
                !Writes SZZ
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                                     WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR030'

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",                 &
                                     "m", Array3D = Me%ExternalVar%SZZ,                         &
                                     OutputNumber = OutPutNumber,STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR40'
                
                !Write OpenPoints3D
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                  &
                                     WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR50'

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",           &
                                     "-", Array3D = Me%ExternalVar%OpenPoints3D,             &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR60'
                

                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Bathymetry", "Bathymetry",      &
                                     "m", Array2D = Me%ExternalVar%Bathymetry,               &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR70'
       
                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/DZ", "DZ",                      &
                                     "m", Array2D = Me%DZ,                                   &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR80'

                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/D50", "D50",                    &
                                     "m", Array2D = Me%D50,                                  &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR90'

                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/DZ_Residual", "DZ_Residual",    &
                                     "m", Array2D =Me%DZ_Residual,                           &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR100'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Bedload", "Bedload",            &
                                    "kg/s/m", Array2D = Me%Bedload,                          &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR110'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Bedload X", "Bedload X",        &
                                    "kg/s", Array2D = Me%FluxX,                              &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR120'

                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Bedload Y", "Bedload Y",        &
                                    "kg/s", Array2D = Me%FluxY,                              &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR125'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results2D/KTop", "KTop",                  &
                                    "-", Array2D = Me%KTop,                                  &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR130'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Porosity", "Porosity",            &
                                     "%", Array3D = Me%Porosity,                             &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR140'
                 

do1:            do n=1,Me%NumberOfClasses

                    call HDF5WriteData   (Me%ObjHDF5, "/Results2D/Classes/TopPercentage/"//trim(Me%SandClass(n)%ID%Name),    &
                                          trim(Me%SandClass(n)%ID%Name),                                                     &
                                          "%", Array2D = Me%SandClass(n)%TopPercentage,                                      &
                                          OutputNumber = OutPutNumber,                                                       &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR150'
                
                    call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Classes/Bedload/"//trim(Me%SandClass(n)%ID%Name),           &
                                        trim(Me%SandClass(n)%ID%Name),                                                       &
                                        "kg/s/m", Array2D = Me%SandClass(n)%Bedload,                                         &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR160'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Classes/Critical Shear Stress/"//                           &
                                        trim(Me%SandClass(n)%ID%Name),                                                       &
                                        trim(Me%SandClass(n)%ID%Name),                                                       &
                                         "N/m2", Array2D = Me%SandClass(n)%CriticalShearStress,                              &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR170'
                    

                    call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Classes/Bedload X/"//trim(Me%SandClass(n)%ID%Name),         &
                                         trim(Me%SandClass(n)%ID%Name),                                                      &
                                         "kg/s", Array2D = Me%SandClass(n)%FluxX,                                            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR180'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results2D/Classes/Bedload Y/"//trim(Me%SandClass(n)%ID%Name),         &
                                        trim(Me%SandClass(n)%ID%Name),                                                       &
                                         "kg/s", Array2D = Me%SandClass(n)%FluxY,                                            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR190'
                    
                    call HDF5WriteData   (Me%ObjHDF5, "/Results/"//trim(Me%SandClass(n)%ID%Name),     &
                                         trim(Me%SandClass(n)%ID%Name),                                                        &
                                         "%", Array3D = Me%SandClass(n)%Field3D,                                               &
                                         OutputNumber = OutPutNumber,                                                        &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR200'
                    
               
            enddo do1

            if (Me%CohesiveClass%Run) then
                
                call HDF5WriteData   (Me%ObjHDF5, "/Results2D/Classes/TopPercentage/Cohesive Sediment",      &
                                        "Cohesive Sediment",                                               &
                                        "%", Array2D = Me%CohesiveClass%TopPercentage,                     &
                                        OutputNumber = OutPutNumber,                                       &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR210'
            
                call HDF5WriteData   (Me%ObjHDF5, "/Results/Cohesive Sediment",           &
                                        "Cohesive Sediment",                                &
                                        "%", Array3D = Me%CohesiveClass%Field3D,            &
                                        OutputNumber = OutPutNumber,                        &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'OutPutSedimentHDF - ModuleSediment - ERR220'
            
                call HDF5WriteData   (Me%ObjHDF5, "/Results/Cohesive_DryDensity",           &
                                        "Cohesive_DryDensity",                              &
                                        "kg/m3", Array3D = Me%CohesiveDryDensity%Field3D,   &
                                         OutputNumber = OutPutNumber,                       &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'OutPutSedimentHDF - ModuleSediment - ERR230'
                
                call HDF5WriteData   (Me%ObjHDF5, "/Results2D/Classes/Critical Shear Stress/Cohesive Sediment",    &
                                        "Cohesive Sediment",                                                       &
                                        "%", Array2D = Me%CohesiveClass%CriticalShearStress,                       &
                                        OutputNumber = OutPutNumber,                                               &
                                        STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR240'
            
            endif
                    
                if (Me%Residual%ON) then
                
                    call HDF5WriteData  (Me%ObjHDF5, "/Residual/Bedload X", &
                                         "Bedload X",                      &
                                         "kg/s", Array2D = Me%Residual%FluxX,            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR250'

                    call HDF5WriteData  (Me%ObjHDF5, "/Residual/Bedload Y", &
                                         "Bedload Y",                      &
                                         "kg/s", Array2D = Me%Residual%FluxY,            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR260'
                
                endif

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR270'

                Me%OutPut%NextOutPut = OutPutNumber + 1

            endif  TOut    

        endif T1


!        if (MonitorPerformance) call StopWatch ("ModuleSediment", "OutPutSedimentHDF")


    end subroutine OutPutSedimentHDF

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL, n

        !Local-----------------------------------------------------------------

        !if (MonitorPerformance) call StartWatch ("ModuleSediment", "OutPut_TimeSeries")

        !call WriteTimeSerie(Me%ObjTimeSerie,                                    &
        !                    Data2D = Me%KTop, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                              &
        !    stop 'OutPut_TimeSeries - ModuleSediment - ERR01'

        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D =Me%DZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR02'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D =Me%DZ_Residual, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR03'

        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D =Me%D50, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR04'
        
        !call WriteTimeSerie(Me%ObjTimeSerie,                                    &
        !                    Data2D =Me%Porosity, STAT = STAT_CALL)
        !if (STAT_CALL /= SUCCESS_)                                              &
        !    stop 'OutPut_TimeSeries - ModuleSediment - ERR05'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 = Me%Bedload, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR06'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%FluxX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR07'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%FluxY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR08'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%FluxToSediment, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR10'  
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%DM, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR11'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%BedloadMass, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR12'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%TotalFluxToSediment, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR14'
        
        call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D =Me%ExternalVar%ShearStress, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR15'

        if (Me%CohesiveClass%Run) then           
            call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                                Data2D =Me%ExternalVar%ConsolidationFlux, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'OutPut_TimeSeries - ModuleSediment - ERR09'  
              
            call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D_8 =Me%CohesiveClass%DM, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'OutPut_TimeSeries - ModuleSediment - ERR13'
        
            call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D =Me%CohesiveClass%CriticalShearStress, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'OutPut_TimeSeries - ModuleSediment - ERR16'
        endif
        
do1:    do n=1,Me%NumberOfClasses
    
            call WriteTimeSerie(Me%ObjTimeSerie,                                    &
                            Data2D =Me%SandClass(n)%ConcRef, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR17'
        enddo do1
        
        
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------


    subroutine OutputBoxFluxes


        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j, WKUB
        integer                                 :: STAT_CALL
        integer                                 :: n
        class(T_Sand), pointer                 :: SandClass

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        
 do1:    do n=1,Me%NumberOfClasses

            SandClass => Me%SandClass(n)

            Me%Boxes%Mass(:,:) = 0.

            do J = JLB, JUB
            do I = ILB, IUB
                
                WKUB = Me%KTop(i,j)
                
                Me%Boxes%Mass   (i,j) = SandClass%DM(i,j)      * (1.- Me%Porosity(i,j,WKUB)) * Me%Density * &
                                        Me%ExternalVar%DUX(i,j) * Me%ExternalVar%DVY(i,j) 

                !This fluxes are initialised and partial computed in the subroutine ComputeEvolution
                Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j)   * (1.- Me%Porosity(i,j,WKUB)) * Me%Density
                Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j)   * (1.- Me%Porosity(i,j,WKUB)) * Me%Density

            end do
            end do
                    
            !Integration of the bottom changes
            call BoxDif(Me%ObjBoxDif, Me%Boxes%Mass,                                        &
                        trim(GetPropertyName (mSediment_)),                                      &
                        Me%ExternalVar%WaterPoints2D,                                       &
                        STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_)                                                      &
                stop 'OutputBoxFluxes - ModuleSediment - ERR01'

            !Integration of fluxes
            call BoxDif(Me%ObjBoxDif,                                                       &
                        Me%Boxes%FluxesX,                                                   &
                        Me%Boxes%FluxesY,                                                   &
                        trim(GetPropertyName (mSediment_)),                                      &
                        Me%ExternalVar%WaterPoints2D,                                       &
                        STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                      &
                stop 'OutputBoxFluxes - ModuleSediment - ERR20'
        
        enddo do1

    end subroutine OutputBoxFluxes

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillSediment(ObjSedimentID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjSedimentID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, STAT_CALL    

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, i

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSediment_,  Me%InstanceID)

            if (nUsers == 0) then
            
                 call ReadLockExternalVar

                if (Me%OutPut%Yes) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSediment - ModuleSediment - ERR10'
                endif

                !Kills the TimeSerie
                if (Me%TimeSerie) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSediment - ModuleSediment - ERR20'
                endif
                
                call WriteFinalState

do1 :               do i=1, Me%NumberOfClasses
                        deallocate(Me%SandClass(i)%Field3D)
                        nullify(Me%SandClass(i)%Field3D)
                        deallocate(Me%SandClass(i)%CriticalShearStress)
                        nullify(Me%SandClass(i)%CriticalShearStress)
                        deallocate(Me%SandClass(i)%NDCriticalShearStress)
                        nullify(Me%SandClass(i)%NDCriticalShearStress)
                        deallocate(Me%SandClass(i)%DM)
                        nullify(Me%SandClass(i)%DM)
                        deallocate(Me%SandClass(i)%FluxX)
                        nullify(Me%SandClass(i)%FluxX)
                        deallocate(Me%SandClass(i)%FluxY) 
                        nullify(Me%SandClass(i)%FluxY) 
                        deallocate(Me%SandClass(i)%Bedload)
                        nullify(Me%SandClass(i)%Bedload)
                        
                        if(Me%WavesOn) then
                            deallocate(Me%SandClass(i)%BedloadX)
                            nullify(Me%SandClass(i)%BedloadX)
                            deallocate(Me%SandClass(i)%BedloadY)
                            nullify(Me%SandClass(i)%BedloadY)
                        endif
                        
                        deallocate(Me%SandClass(i)%TopPercentage)
                        nullify(Me%SandClass(i)%TopPercentage)
                        deallocate(Me%SandClass(i)%HidingFactor)
                        nullify(Me%SandClass(i)%HidingFactor)
                        deallocate(Me%SandClass(i)%ConcRef)
                        nullify(Me%SandClass(i)%ConcRef)                        
                    enddo do1
                    
                    deallocate(Me%SandClass)
                    nullify(Me%SandClass)
                    deallocate(Me%TotalPercentage)
                    nullify(Me%TotalPercentage)
                    
                if(Me%CohesiveClass%Run) then
                    deallocate(Me%CohesiveClass%CriticalShearStress)
                    nullify(Me%CohesiveClass%CriticalShearStress)
                    deallocate(Me%CohesiveClass%PM1)
                    nullify(Me%CohesiveClass%PM1) 
                    deallocate(Me%CohesiveClass%TopPercentage)
                    nullify(Me%CohesiveClass%TopPercentage)    
                    deallocate(Me%CohesiveClass%DM)
                    nullify(Me%CohesiveClass%DM)
                    deallocate(Me%CohesiveClass%Porosity)
                    nullify(Me%CohesiveClass%Porosity)
                endif
                
                if (Me%Evolution%Bathym) then
                    deallocate(Me%BatimIncrement)
                    nullify(Me%BatimIncrement)
                endif

                deallocate(Me%DZ)
                nullify(Me%DZ)
                deallocate(Me%DZ_Residual)
                nullify(Me%DZ_Residual)
                deallocate(Me%FluxX)
                nullify(Me%FluxX)
                deallocate(Me%FluxY)
                nullify(Me%FluxY)
                deallocate(Me%Bedload)
                nullify(Me%Bedload)
                deallocate(Me%DM)
                nullify(Me%DM)
                deallocate(Me%Mass)
                nullify(Me%Mass)
                deallocate(Me%FluxToSediment)
                nullify(Me%FluxToSediment)
                deallocate(Me%BedloadMass)
                nullify(Me%BedloadMass)
                deallocate(Me%TotalFluxToSediment)
                nullify(Me%TotalFluxToSediment)
                deallocate(Me%Dast)
                nullify(Me%Dast)
                
                deallocate (Me%D50, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &            
                    stop 'KillSediment - ModuleSediment. ERR01.' 
                nullify (Me%D50) 
                
                deallocate (Me%SandD50, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &            
                    stop 'KillSediment - ModuleSediment. ERR02.' 
                nullify (Me%SandD50) 
                        
                deallocate (Me%NDShearStress, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR03.' 
                nullify (Me%NDShearStress) 
                
                if (Me%WavesOn) then
                    deallocate (Me%NDShearStressWaves, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_)                                          &          
                        stop 'KillSediment - ModuleSediment. ERR03a.' 
                    nullify (Me%NDShearStressWaves)
                    
                    deallocate (Me%NDShearStressMean, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_)                                          &          
                        stop 'KillSediment - ModuleSediment. ERR03b.' 
                    nullify (Me%NDShearStressMean)
                endif
                
                deallocate(Me%Elevation, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR04.' 
                nullify (Me%Elevation)
                
                deallocate(Me%Porosity, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR05.' 
                nullify (Me%Porosity) 
            
                deallocate(Me%VerticalCoordinate, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR05.' 
                nullify (Me%VerticalCoordinate) 

                deallocate(Me%KTop, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR06.' 
                nullify (Me%KTop) 
                
                deallocate(Me%OpenSediment, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR07.' 
                nullify (Me%OpenSediment)
                
                if (Me%Residual%ON) then
                    deallocate(Me%Residual%FluxX)
                    nullify(Me%Residual%FluxX)
                    deallocate(Me%Residual%FluxY)
                    nullify(Me%Residual%FluxY)
                endif
                
                if (Me%Boxes%Yes) then

                    call KillBoxDif(Me%ObjBoxDif, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillSediment - ModuleSediment - ERR30'

                    deallocate(Me%Boxes%FluxesX)
                    nullify   (Me%Boxes%FluxesX)

                    deallocate(Me%Boxes%FluxesY)
                    nullify   (Me%Boxes%FluxesY)


                    deallocate(Me%Boxes%Mass   )
                    nullify   (Me%Boxes%Mass   )

                endif


                call ReadUnLockExternalVar

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR40'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjBathym)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR50'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR60'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR70'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjSedimentGridData)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR17'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjSedimentHorizontalMap)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR18'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjSedimentGeometry)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR19'

                nUsers = DeassociateInstance(mMAP_,             Me%ObjSedimentMap)
                if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR20'
#ifndef _WAVES_
                if(Me%ObjWaves /= 0)then
                    nUsers = DeassociateInstance (mWAVES_,Me%ObjWaves)
                    if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR80'
                end if
#endif

                if (Me%Discharges%Yes)    then
                    nUsers = DeassociateInstance(mDISCHARGES_,    Me%ObjDischarges)
                    if (nUsers == 0) stop 'KillSediment - ModuleSediment - ERR90'
                end if
                
                !Deallocates Instance
                call DeallocateInstance ()

                ObjSedimentID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillSediment
        
    !--------------------------------------------------------------------------
    !   Write the final sediment results in HDF format  !
  
    subroutine WriteFinalState

        !Local--------------------------------------------------------------
        real,    dimension(6    ), target      :: AuxTime
        real,    dimension(:    ), pointer     :: TimePtr
        integer                                :: WorkILB, WorkIUB
        integer                                :: WorkJLB, WorkJUB
        integer                                :: WorkKLB, WorkKUB
        integer                                :: STAT_CALL, n
        !----------------------------------------------------------------------

        !Bounds
        WorkILB = Me%SedimentWorkSize3D%ILB 
        WorkIUB = Me%SedimentWorkSize3D%IUB 
        WorkJLB = Me%SedimentWorkSize3D%JLB 
        WorkJUB = Me%SedimentWorkSize3D%JUB 
        WorkKLB = Me%SedimentWorkSize3D%KLB 
        WorkKUB = Me%SedimentWorkSize3D%KUB 
        
        call Open_HDF5_OutPut_File(Me%Files%Final)
        
        call ReadUnLockExternalVar
        
        !Writes geometry
        call WriteGeometryHDF(Me%ObjSedimentGeometry,                                    &
                              Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR01'
        
        call ReadLockExternalVar

       
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                               &
                             WorkJLB, WorkJUB,                                           &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR10'

        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                     &
                             "DZ_Residual",                                              &
                             "m",                                                        &
                             Array2D =Me%DZ_Residual,                                    &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR20'

        call HDF5WriteData  (Me%ObjHDF5, "/Results",                                     &
                             "Bathymetry",                                               &
                             "m",                                                        &
                             Array2D = Me%ExternalVar%Bathymetry,                        &
                             STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR30'

        if (Me%Residual%ON) then
            call HDF5WriteData  (Me%ObjHDF5, "/Results",                                &
                                 "Bedload X",                                           &
                                 "kg/s",                                              &
                                 Array2D = Me%Residual%FluxX,                           &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR40'

            call HDF5WriteData  (Me%ObjHDF5, "/Results",                                &
                                 "Bedload Y",                                           &
                                 "kg/s",                                              &
                                 Array2D = Me%Residual%FluxY,                           &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR50'

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR60'

            !Writes current time
            call ExtractDate   (Me%Residual%StartTime, AuxTime(1), AuxTime(2), AuxTime(3), &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalState - ModuleSediment - ERR70'

            call HDF5WriteData  (Me%ObjHDF5, "/Results",                                &
                                 "Start Time", "YYYY/MM/DD HH:MM:SS",                   &
                                 Array1D = TimePtr, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalState - ModuleSediment - ERR80'
                                 

        endif
        
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                               &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,                         &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR13'
            
        do n=1,Me%NumberOfClasses
            
            call HDF5WriteData   (Me%ObjHDF5, "/Results3D",                              &
                                    trim(Me%SandClass(n)%ID%Name),                       &
                                    "%", Array3D = Me%SandClass(n)%Field3D,              &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'WriteFinalState - ModuleSediment - ERR14'
            
        enddo
        
        if (Me%CohesiveClass%Run) then
            
            call HDF5WriteData   (Me%ObjHDF5, "/Results3D",                             &
                        trim(Me%CohesiveClass%ID%Name),                                 &
                        "%", Array3D = Me%CohesiveClass%Field3D,                        &
                        STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR15'
            
            call HDF5WriteData   (Me%ObjHDF5, "/Results3D",         &
            "Cohesive_DryDensity",                                                      &
            "kg/m3", Array3D = Me%CohesiveDryDensity%Field3D,                           &
            STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR16'
            
        endif
   
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR90'

        call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR100'

    end subroutine WriteFinalState

    !--------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Sediment), pointer          :: AuxObjSediment
        type (T_Sediment), pointer          :: PreviousObjSediment

        !Updates pointers
        if (Me%InstanceID == FirstObjSediment%InstanceID) then
            FirstObjSediment => FirstObjSediment%Next
        else
            PreviousObjSediment => FirstObjSediment
            AuxObjSediment      => FirstObjSediment%Next
            do while (AuxObjSediment%InstanceID /= Me%InstanceID)
                PreviousObjSediment => AuxObjSediment
                AuxObjSediment      => AuxObjSediment%Next
            enddo

            !Now update linked list
            PreviousObjSediment%Next => AuxObjSediment%Next

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

    subroutine Ready (ObjSediment_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSediment_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjSediment_ID > 0) then
            call LocateObjSediment (ObjSediment_ID)
            ready_ = VerifyReadLock (mSediment_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjSediment (ObjSedimentID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID

        !Local-----------------------------------------------------------------

        Me => FirstObjSediment
        do while (associated (Me))
            if (Me%InstanceID == ObjSedimentID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleSediment - LocateObjSediment - ERR01'

    end subroutine LocateObjSediment

    !--------------------------------------------------------------------------


    subroutine ReadLockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !Now
        call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR01'

        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR02'

        !OpenPoints2D
        call GetOpenPoints2D (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR03'

        !BoundaryPoints2D
        call GetBoundaries(Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR05'

        !Compute faces 2D
        call GetComputeFaces2D(Me%ObjHorizontalMap,                                      &
                               ComputeFaces2DU = Me%ExternalVar%ComputeFacesU2D,         &
                               ComputeFaces2DV = Me%ExternalVar%ComputeFacesV2D,         &
                               ActualTime      = Me%ExternalVar%Now,                     &
                               STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR06'


        call GetHorizontalGrid(Me%ObjHorizontalGrid,                                     &
                               DUX  = Me%ExternalVar%DUX,                                &
                               DVY  = Me%ExternalVar%DVY,                                &
                               DZX  = Me%ExternalVar%DZX,                                &
                               DZY  = Me%ExternalVar%DZY,                                &
                               DXX  = Me%ExternalVar%DXX,                                &
                               DYY  = Me%ExternalVar%DYY,                                &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR07'

        
        call GetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR08'

        !SZZ
        call GetGeometryDistances (Me%ObjSedimentGeometry, SZZ = Me%ExternalVar%SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR14'
        
        !DWZ
        call GetGeometryDistances (Me%ObjSedimentGeometry, DWZ = Me%ExternalVar%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR15'

        call GetGeometryKTop(Me%ObjSedimentGeometry,                                     &
                             KTopZ  = Me%ExternalVar%KTop,                               &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR16'
        
        call GetGeometryVolumes(Me%ObjSedimentGeometry,                                  &
                        VolumeZ = Me%ExternalVar%VolumeZ,                                &
                        STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR17'
        

#ifndef _WAVES_
        if (Me%ObjWaves /=0) then        

            call GetWaves (WavesID       = Me%ObjWaves,                                  &
                           WavePeriod    = Me%ExternalVar%WavePeriod,                    &
                           WaveHeight    = Me%ExternalVar%WaveHeight,                    &
                           Abw           = Me%ExternalVar%Abw,                           &
                           Ubw           = Me%ExternalVar%Ubw,                           &
                           WaveDirection = Me%ExternalVar%WaveDirection,                 &
                           STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadLockExternalVar - ModuleSediment - ERR18'

        endif
#endif

        call GetOpenPoints3D(Me%ObjSedimentMap,                                     &
                             Me%ExternalVar%OpenPoints3D,                           &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModuleSediment. ERR08")


        call GetWaterPoints3D(Me%ObjSedimentMap,                                    &
                              Me%ExternalVar%WaterPoints3D,                         &
                              STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModuleSediment. ERR19")
        
        call GetGridCellArea (Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            call SetError(FATAL_, INTERNAL_, "ReadLockExternalVar; ModuleSediment. ERR20")

    end subroutine ReadLockExternalVar

    !--------------------------------------------------------------------------

    subroutine ReadUnLockExternalVar
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        
        !WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR01'

        !OpenPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR02'

        !DXX
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid,                                  &
                                  Me%ExternalVar%DXX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR03'

        !DYY
        call UnGetHorizontalGrid (Me%ObjHorizontalGrid, Me%ExternalVar%DYY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR04'


        !BoundaryPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%BoundaryPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR05'

        !Compute faces 2D V
        call UnGetHorizontalMap(Me%ObjHorizontalMap,                                     &
                               Me%ExternalVar%ComputeFacesV2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR06'

        !Compute faces 2D U
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%ComputeFacesU2D,     &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR07'


        !DUX
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DUX,               &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR08'

        !DVY
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                   &
                               Me%ExternalVar%DVY,                                       &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR09'


        !DZX
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DZX,               &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR10'

        !DZY
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                   &
                               Me%ExternalVar%DZY,                                       &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR11'

        !Bathymetry
        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%Bathymetry, STAT_CALL)     
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR12'
        
        call UnGetGeometry(Me%ObjSedimentGeometry,Me%ExternalVar%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnLockExternalVar - ModuleSediment - ERR13'
        
        call UnGetGeometry(Me%ObjSedimentGeometry,Me%ExternalVar%SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnLockExternalVar - ModuleSediment - ERR14'


        call UnGetGeometry(Me%ObjSedimentGeometry, Me%ExternalVar%KTop, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) &
            stop 'ReadUnLockExternalVar - ModuleSediment - ERR15'  
        
        call UnGetMap(Me%ObjSedimentMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'ReadUnLockExternalVar - ModuleSediment - ERR16'  
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%GridCellArea, STAT = STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)                                    &
            stop 'ReadUnLockExternalVar - ModuleSediment - ERR17'
        
#ifndef _WAVES_

        if (Me%ObjWaves /=0) then

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WavePeriod, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR18'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WaveHeight, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR19'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%Abw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR20'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%Ubw, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR21'

            call UnGetWaves(Me%ObjWaves, Me%ExternalVar%WaveDirection, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR22'

        endif
#endif

            call UnGetGeometry(Me%ObjSedimentGeometry, Me%ExternalVar%VolumeZ, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR23'

            call UnGetMap(Me%ObjSedimentMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockExternalVar - ModuleSediment - ERR23' 

    end subroutine ReadUnLockExternalVar
    !--------------------------------------------------------------------------

end module ModuleSediment

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------