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
                                       GetGridOutBorderPolygon
    use ModuleBoxDif,           only : StartBoxDif, GetBoxes, GetNumberOfBoxes, BoxDif,         &
                                       UngetBoxDif, KillBoxDif
    use ModuleGeometry,         only:  GetGeometrySize, UnGetGeometry,                      &
                                       GetGeometryDistances, GetGeometryKtop,               &
                                       ComputeInitialGeometry, ComputeVerticalGeometry,     &
                                       GetGeometryVolumes, ReadGeometryHDF, WriteGeometryHDF
#ifndef _WAVES_
    use ModuleWaves,            only : GetWaves, UnGetWaves
#endif

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructSediment
    private ::      AllocateInstance
    private ::      AllocateVariables
    private ::      ConstructEvolution
    private ::      ConstructClasses
    private ::      ConstructGlobalParameters
    private ::      Construct_Initial_Geometry
    private ::      StartOutputBoxFluxes
    private ::      Open_HDF5_OutPut_File
    private ::      ConstructTimeSerie
    private ::      Read_Sediment_Files_Name
    private ::      ReadInitialField
    private ::      ReadInitialField3D


    !Selector
    public  :: UnGetSediment                 
    
    !Modifier
    public  :: ModifySediment
    private ::      ComputeNDShearStress
    private ::      ComputeCriticalShearStress
    private ::      ComputeTransport
    private ::          CurrentOnly
    private ::          CurrentPlusWaves
    private ::          WavesOnly
    private ::      ComputeFluxes
    private ::      ComputeEvolution
    private ::      ComputeDischarges   
    private ::      BoundaryCondition
    private ::      ComputeTotalDZ
    private ::      ComputeResidualEvolution
    private ::      ComputePercentage
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
    interface  UnGetSediment
        module procedure UnGetSediment2D_I
        module procedure UnGetSediment2D_R4
        module procedure UnGetSediment2D_R8
    end interface  UnGetSediment

    !Parameters
    character(LEN = StringLength), parameter    :: class_block_begin     = '<beginclass>'
    character(LEN = StringLength), parameter    :: class_block_end       = '<endclass>'

    integer, parameter :: NoTransport = 0, CurrentOnly = 1, CurrentPlusWaves = 2, WavesOnly = 3
    
    integer, parameter :: NullGradient = 1, Cyclic = 2, NullValue = 3

    !Selma
    integer, parameter :: Time_ = 1
    
    !srt(2.)
    real,    parameter :: SquareRoot2 = 1.414213562

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
        real,    pointer, dimension(:,:)        :: ShearStress         => null()
        real,    pointer, dimension(:,:)        :: CurrentRugosity  => null()
        real,    pointer, dimension(:,:)        :: WaveRugosity     => null()
        real,    pointer, dimension(:,:)        :: VelU             => null()
        real,    pointer, dimension(:,:)        :: VelV             => null()
        real,    pointer, dimension(:,:)        :: VelU_Face        => null()
        real,    pointer, dimension(:,:)        :: VelV_Face        => null()
        real,    pointer, dimension(:,:)        :: VelMod           => null()
        real,    pointer, dimension(:,:)        :: TauWave          => null()
        real,    pointer, dimension(:,:)        :: TauCurrent       => null()
        real,    pointer, dimension(:,:)        :: ShearVelocity    => null()    
        real,    pointer, dimension(:,:)        :: WaveHeight       => null()
        real,    pointer, dimension(:,:)        :: WavePeriod       => null()
        !Sediment
        real,    pointer, dimension(:,:,:)      :: DWZ              => null()
        real,    pointer, dimension(:,:,:)      :: SZZ              => null()
        integer, pointer, dimension(:,:  )      :: KTop             => null()
        real(8),    pointer, dimension(:,:,:)   :: VolumeZ          => null()

    end type T_External

    private :: T_Residual
    type       T_Residual
        logical                                 :: ON           = .false.
        type(T_Time)                            :: StartTime    
        real, dimension(:,:), pointer           :: FluxX        => null ()
        real, dimension(:,:), pointer           :: FluxY        => null ()        
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

    private :: T_Evolution
    type       T_Evolution
        logical                                 :: Old          = .false. 
        real                                    :: SedimentDT       = FillValueReal
        real                                    :: BathymDT     = FillValueReal
        type (T_Time)                           :: NextSediment, NextBatim
        logical                                 :: Bathym       = .false. 
        !Selma
        !integer                                 :: BathymType   = Time_
    end type T_Evolution
     
    private :: T_Property
    type ::       T_Property
        type (T_PropertyID)                     :: ID
!        type (T_SubModel  )                     :: SubModel
        real                                    :: Scalar       = FillValueReal
        real, pointer, dimension(:,:,:)         :: Field3D      => null ()
    end type T_Property

    public :: T_Class
    type, extends (T_Property) :: T_Class
        real(8)                                 :: D50                  = FillValueReal
        real                                    :: parameter_A          = FillValueReal, &
                                                   parameter_n          = FillValueReal, &
                                                   parameter_p          = FillValueReal
                                                               
        real(8), dimension(:,:), pointer        :: FluxX                 => null ()
        real(8), dimension(:,:), pointer        :: FluxY                 => null ()
        real(8), dimension(:,:), pointer        :: Bedload               => null ()
        real, dimension(:,:), pointer           :: CriticalShearStress   => null ()
        real, dimension(:,:), pointer           :: NDCriticalShearStress => null ()
        real(8), dimension(:,:), pointer        :: DZ                    => null ()
        real(8), dimension(:,:), pointer        :: TopPercentage         => null ()
        !Field3D in the base class T_Property is the percentage of this class on each domain cell
    end type T_Class 
    
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
        type (T_Evolution )                        :: Evolution
        real, dimension(:, :), pointer             :: BatimIncrement, DZ, DZ_Residual => null () 
        type (T_Aceleration)                       :: Aceleration
        real                                       :: Porosity              = FillValueReal
        real                                       :: Density               = FillValueReal
        real                                       :: RelativeDensity       = FillValueReal
        real                                       :: DeltaDensity          = FillValueReal
        integer                                    :: Boundary              = FillValueInt
        real                                       :: TauMax                = FillValueReal
        logical                                    :: TimeSerie             = .false. 
        class (T_Class), dimension(:), pointer     :: Classes
        integer                                    :: NumberOfClasses
        real, dimension(:, :, :), pointer          :: TotalPercentage       => null () 
        real, dimension(:, :), pointer             :: D50                   => null () 
        real, dimension(:, :), pointer             :: NDShearStress         => null ()
        real(8), dimension(:,:), pointer           :: FluxX                 => null ()
        real(8), dimension(:,:), pointer           :: FluxY                 => null ()
        real(8), dimension(:,:), pointer           :: Bedload               => null ()
        
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

        real, dimension(:,:), pointer              :: Dast                  => null()
        real, pointer, dimension(:,:)              :: Elevation             => null()
        real, pointer, dimension(:,:,:)            :: VerticalCoordinate    => null()
        integer, pointer, dimension(:,:)           :: KTop                  => null()
        real                                       :: MinLayerThickness    = null_real
        real                                       :: MaxLayerThickness    = null_real
        
        integer                                    :: TransportMethod       = FillValueInt
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
        integer, optional, intent(OUT)                  :: STAT
           
        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL
        real                                            :: WaterDensity
        logical                                         :: WaveTensionON

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_
        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

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
            Me%ExternalVar%WaveTensionON = WaveTensionON

            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                             &
                                       Size        = Me%Size,                            &
                                       WorkSize    = Me%WorkSize,                        &
                                       STAT        = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR01'
            
            
            call GetGeometrySize(Me%ObjSedimentGeometry,                        &
                                 Size       = Me%SedimentSize3D,                &
                                 WorkSize   = Me%SedimentWorkSize3D,            &
                                 STAT       = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)&
                stop 'ConstructGlobalVariables - ModuleSediment - ERR03'

            call ReadLockExternalVar

            call Read_Sediment_Files_Name

            !Construct enter data 
            call ConstructEnterData(Me%ObjEnterData, Me%Files%ConstructData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR20'
            
            call ConstructGlobalParameters
            
            call AllocateVariables
            
            call ConstructEvolution
            
            call Construct_Initial_Geometry

            call ConstructOutputTime

            call ConstructClasses

            call StartOutputBoxFluxes
            
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructSediment - ModuleSediment - ERR30'


            if (Me%OutPut%Yes) call Open_HDF5_OutPut_File(Me%Files%OutPutFields)

            if (Me%Discharges%Yes) then
            
                if (ObjDischargesID == 0)  then                                                
                    write(*,*)'You need to define a water discharges in the hydrodynamic input' 
                    stop      'ConstructSediment - Sediment - ERR01'
                else
                    Me%ObjDischarges = AssociateInstance (mDISCHARGES_, ObjDischargesID)

                    Me%Discharges%NextCompute = Me%ExternalVar%Now
                endif

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
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-------------------------------------------------------------------
        integer                             :: iflag
        character(Len = StringLength)       :: Auxchar
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


        call GetData(Me%Porosity,                                                        &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'POROSITY',                                          &
                     default      = 0.1,                                                 &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR12' 


        call GetData(Me%Density,                                                         &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'DENSITY',                                           &
                     default      = 2650.,                                               &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)               
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR13' 
        
        Me%RelativeDensity = Me%Density / Me%ExternalVar%WaterDensity

        Me%DeltaDensity    = Me%Density - Me%ExternalVar%WaterDensity
        
        call GetData(Me%TransportMethod,                                                 &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'TRANSPORT_METHOD',                                  &
                     default      = 1,                                                   &
                     ClientModule = 'ModuleSediment',                                    &
                     STAT         = STAT_CALL)              
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalParameters - ModuleSediment - ERR14' 

        !NoTransport = 0, CurrentOnly = 1, CurrentPlusWaves = 2, WavesOnly = 3
       
        if     ( Me%TransportMethod == 2 &
            .OR. Me%TransportMethod == 3 )  then
        
            if (Me%ObjWaves == 0) stop 'ConstructGlobalParameters - ModuleSediment - ERR16' 

            if (.not. Me%ExternalVar%WaveTensionON) stop 'ConstructGlobalParameters - ModuleSediment - ERR20' 

        endif

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

        call GetData(Me%Aceleration%Coef,                                                &
                     Me%ObjEnterData,iflag,                                              &
                     SearchType   = FromFile,                                            &
                     keyword      = 'ACELERATION',                                       &
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
            
        allocate (Me%Dast (ILB:IUB, JLB:JUB))
        Me%Dast(:,:) = FillValueReal
            
        allocate(Me%NDShearStress (ILB:IUB, JLB:JUB))
        Me%NDShearStress(:,:) = FillValueReal         
                 
        allocate(Me%Elevation(ILB:IUB, JLB:JUB)) 
        Me%Elevation(:,:) = 0.
            
        allocate(Me%VerticalCoordinate(ILB:IUB, JLB:JUB, KLB:KUB)) 
        Me%VerticalCoordinate(:,:,:) = null_real
        
        allocate(Me%KTop(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB)) 
        Me%KTop(:,:) = 0.
            
        allocate(Me%FluxX(ILB:IUB, JLB:JUB))
        Me%FluxX(:,:) = null_real

        allocate(Me%FluxY(ILB:IUB, JLB:JUB))
            Me%FluxY(:,:) = null_real

        allocate(Me%Bedload(ILB:IUB, JLB:JUB))
        Me%Bedload(:,:) = null_real
            
         allocate(Me%DZ(ILB:IUB, JLB:JUB))
        Me%DZ(:,:) = 0.

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
                
                call ReadInitialField(FieldName = "Transport Flux X", Field2D = Me%Residual%FluxX)
                
                call ReadInitialField(FieldName = "Transport Flux Y", Field2D = Me%Residual%FluxY)                
            
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
        integer                                     :: HDF5_CREATE, n

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
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR01'


        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,         &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR02'

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,        &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR03'
        
        
        if (Me%Evolution%Bathym) then

            call GetGridData2Dreference(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR04' 

        else

            call GetGridData           (Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR04' 
        
        endif

        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                    &
                              Array2D = Me%ExternalVar%InitialBathym,                    &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR05'

        call UnGetGridData(Me%ObjBathym, Me%ExternalVar%InitialBathym, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR06' 


        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                 &
                              Array2D = Me%ExternalVar%WaterPoints2D,                    &
                              STAT    = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR07'
        
        !Write WaterPoints3D
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                  &
                            WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR03'
                
                
        call HDF5WriteData  (Me%ObjHDF5, "/Grid", "WaterPoints3D",         &
                            "-", Array3D = Me%ExternalVar%WaterPoints3D,             &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR04a'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR15'

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


        !----------------------------------------------------------------------

        !First checks out how many properties will have time series

        !Allocates PropertyList
        allocate(PropertyList(2), STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'ConstructTimeSerie - ModuleSediment - ERR10'

        PropertyList(1) = 'DZ'
        PropertyList(2) = 'DZ_Residual'


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
        logical                             :: BlockFound, BlockInBlockFound

        !Local-------------------------------------------------------------------
        logical, allocatable, dimension(:)  :: ClassOK
        integer                             :: ClassID

        !------------------------------------------------------------------------

        integer                         :: ILB, IUB, JLB, JUB, KLB, KUB, iflag
        integer                         :: i, j, k, n, WKUB        
        real                            :: Percentage
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
        
            allocate(Me%Classes (Me%NumberOfClasses))
            
do1 :       do n=1, Me%NumberOfClasses
                call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                                            ClientNumber    = ClientNumber,             &
                                            block_begin     = class_block_begin,         &
                                            block_end       = class_block_end,           &
                                            BlockFound      = BlockFound,                &
                                            STAT            = STAT_CALL)
cd1 :           if (STAT_CALL .EQ. SUCCESS_) then    
cd2 :               if (BlockFound) then    
          
                        call GetData(Me%Classes(n)%D50,                                     &
                                     Me%ObjEnterData,iflag,                                 &
                                     SearchType   = FromBlock,                              &
                                     keyword      = 'D50',                                  &
                                     ClientModule = 'ModuleSediment',                           &
                                     STAT         = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR01a' 
                        
                        if (iflag .NE. 1) stop 'ConstructClasses - ModuleSediment - ERR01a' 
                        
                        call GetData(Me%Classes(n)%parameter_A,                             &
                                     Me%ObjEnterData,iflag,                                 &
                                     SearchType   = FromBlock,                              &
                                     keyword      = 'A',                                    &
                                     default      = 12.,                                     &
                                     ClientModule = 'ModuleSediment',                       &
                                     STAT         = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR01b' 
                         
                        call GetData(Me%Classes(n)%parameter_n,                             &
                                     Me%ObjEnterData,iflag,                                 &
                                     SearchType   = FromBlock,                              &
                                     keyword      = 'n',                                    &
                                     default      = 0.5,                                    &
                                     ClientModule = 'ModuleSediment',                       &
                                     STAT         = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR01c'        
                        
                        call GetData(Me%Classes(n)%parameter_p,                             &
                                     Me%ObjEnterData,iflag,                                 &
                                     SearchType   = FromBlock,                              &
                                     keyword      = 'p',                                    &
                                     default      = 1.,                                     &
                                     ClientModule = 'ModuleSediment',                       &
                                     STAT         = STAT_CALL)

                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR01d'  
                                      
                        call ConstructSedimentProperty (Me%Classes(n), FromBlock)
                        
                        
                        !Allocate fluxes
                        allocate(Me%Classes(n)%FluxX(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%FluxX(:,:) = null_real

                        allocate(Me%Classes(n)%FluxY(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%FluxY(:,:) = null_real

                        allocate(Me%Classes(n)%Bedload(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%Bedload(:,:) = null_real

                        allocate(Me%Classes(n)%CriticalShearStress(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%CriticalShearStress(:,:) = null_real
                        
                        allocate(Me%Classes(n)%NDCriticalShearStress(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%CriticalShearStress(:,:) = null_real 
                        
                        allocate(Me%Classes(n)%TopPercentage(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%TopPercentage(:,:) = null_real
                        
                        allocate(Me%Classes(n)%DZ(ILB:IUB, JLB:JUB))
                        Me%Classes(n)%DZ(:,:) = 0. 
                    else cd2

                        write(*,*)  
                        write(*,*) 'Error calling ExtractBlockFromBlock. '
                        stop       'ConstructClasses - ModuleSediment - ERR02'

                    end if cd2

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                    write(*,*)  
                    write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructClasses - ModuleSediment - ERR03'
                else cd1
                    stop       'ConstructClasses - ModuleSediment - ERR04'
                end if cd1
            end do do1

            call Block_UnLock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR05'


            call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL  /= SUCCESS_) stop 'ConstructClasses - ModuleSediment - ERR06'
            
            
            allocate(Me%TotalPercentage(ILB:IUB, JLB:JUB, KLB:KUB))
            Me%TotalPercentage(:,:,:)=0
            
do2 :       do n=1,Me%NumberOfClasses
                do j=JLB, JUB
                do i=ILB, IUB
                do k=KLB, Me%KTop(i, j)
                    
                    if (Me%ExternalVar%WaterPoints3D (i,j,k) == WaterPoint) then
                    
                        Percentage = Me%Classes(n)%Field3D (i,j,k)
                        
                        Me%TotalPercentage (i, j, k) = Me%TotalPercentage (i, j, k) + Percentage
                    endif
                enddo
                enddo
                enddo
             enddo do2                    
                        
                
                do j=JLB, JUB
                do i=ILB, IUB
                do k=KLB, Me%KTop(i, j)
                    
                    if (Me%ExternalVar%WaterPoints3D (i,j,k) == WaterPoint) then
                        if (Me%TotalPercentage (i, j, k)  > 1.001) then
                            write(*,*) 'The sum of the classes percentage is larger than 100%.'
                            write(*,*) i, j, k, Me%TotalPercentage (i, j, k)
                            stop 'ConstructClasses - ModuleSediment - ERR07'
                        elseif (Me%TotalPercentage (i, j, k)  < 0.999) then
                                write(*,*) 'The sum of the classes percentage is smaller than 100%.'
                                write(*,*) i, j, k, Me%TotalPercentage (i, j, k)
                                stop 'ConstructClasses - ModuleSediment - ERR08'
                        endif
                    endif
                enddo
                enddo
                enddo
    
    end subroutine ConstructClasses
    
     !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a        
    !Sediment property values in the domain    

    subroutine ConstructSedimentProperty(NewProperty, ExtractType)

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
        
        call ConstructPropertyID(NewProperty%ID, Me%ObjEnterData, ExtractType, CheckProperty = .false.)
        NewProperty%ID%IDNumber = RegisterDynamicProperty (NewProperty%ID%Name)
              
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
                stop 'Construct_SedimentProperty - ModuleSediment - ERR02'

            call GetDefaultValue(NewProperty%ID%ObjFillMatrix, NewProperty%Scalar, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then                                                
                print *, NewProperty%ID%Name, ' was not build correctly.'
                stop 'Construct_SedimentProperty - ModuleSediment - ERR03'
            endif

            if(NewProperty%ID%SolutionFromFile)then

                stop 'Construct_SedimentProperty - ModuleSediment - ERR04'
            
            else

                call KillFillMatrix(NewProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)&
                    stop 'Construct_SedimentProperty - ModuleSediment - ERR05'

            endif
        endif
        
        call CheckFieldConsistence (NewProperty)


    end subroutine ConstructSedimentProperty

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

        !Class percentage in empty layers is set as null
        do I = ILB, IUB
        do J = JLB, JUB
        do K = KLB, KUB
            
            if (Me%ExternalVar%WaterPoints3D(i, j, k) == WaterPoint) then
                               
                if (Me%ExternalVar%VolumeZ(i, j, k) == 0.) then

                    NewProperty%Field3D(i, j, k) = 0.
                    
                end if

            else

                NewProperty%Field3D(i, j, k) = 0.

            endif

        enddo
        enddo
        enddo

    end subroutine CheckFieldConsistence

    
    !----------------------------------------------------------------------
    
    subroutine Construct_Initial_Geometry 
        
        !Local----------------------------------------------------------------
        integer                             :: STAT_CALL
        integer                             :: i, j, k

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


        !Computes OpenPoints3D
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
    subroutine GetClasses (ObjSedimentID, Classes, STAT)
        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSedimentID
        type (T_Class), dimension(:), pointer           :: Classes
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSediment_, Me%InstanceID)
                Classes => Me%Classes
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    
    end subroutine GetClasses 
    !--------------------------------------------------------------------------

    subroutine UnGetSedimentClasses(ObjSedimentID, Classes, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjSedimentID
        type (T_Class), dimension(:), pointer           :: Classes
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Classes)
            call Read_Unlock(mSediment_, Me%InstanceID, "UnGetSedimentClasses")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetSedimentClasses

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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifySediment(ObjSedimentID, ShearStress, CurrentRugosity, WaveRugosity,   &
                          VelU, VelV, VelMod, TauWave, TauCurrent, &
                          ShearVelocity, VelU_Face, VelV_Face, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjSedimentID
        real, dimension(:,:), pointer               :: ShearStress, CurrentRugosity, WaveRugosity, &
                                                       VelU, VelV, VelMod,         &
                                                       TauWave, TauCurrent, ShearVelocity,      &
                                                       VelU_Face, VelV_Face     
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, STAT_CALL
        logical                                     :: ChangeBathym                                                
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjSedimentID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call ReadLockExternalVar 

            if (Me%TransportMethod /= 0) then

do1:            do while (Me%ExternalVar%Now >= Me%Evolution%NextSediment) 

                    Me%ExternalVar%ShearStress        => ShearStress

                    Me%ExternalVar%CurrentRugosity => CurrentRugosity
                    Me%ExternalVar%WaveRugosity    => WaveRugosity

                    Me%ExternalVar%VelU            => VelU
                    Me%ExternalVar%VelV            => VelV
                    Me%ExternalVar%VelMod          => VelMod
                    Me%ExternalVar%VelU_Face       => VelU_Face
                    Me%ExternalVar%VelV_Face       => VelV_Face

                    Me%ExternalVar%TauWave         => TauWave
                    Me%ExternalVar%TauCurrent      => TauCurrent
                    Me%ExternalVar%ShearVelocity   => ShearVelocity

                        call ComputeD50Cell
                        
                        call ComputeNDShearStress
                    
                        call ComputeCriticalShearStress
                    
                        call ComputeTransport
                        
                        call ComputeFluxes

                        call ComputeEvolution
                        
                        if (Me%Discharges%Yes) then 
                            call ComputeDischarges
                        endif           

                        call BoundaryCondition
                        
                        call ComputeTotalDZ
                        
                        call ComputePercentage
                        
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

                        Me%Evolution%NextSediment = Me%Evolution%NextSediment + Me%Evolution%SedimentDT

                enddo do1

            endif

            call ReadUnLockExternalVar

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifySediment

    !--------------------------------------------------------------------------

    
     !--------------------------------------------------------------------------
    
    subroutine ComputeD50Cell
       !Local-----------------------------------------------------------------
        integer             :: i, j, n, WKUB
        integer             :: WILB, WIUB, WJLB, WJUB
        class(T_Class), pointer :: SedimentClass
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%D50(:,:)= 0
        
do1:    do n=1,Me%NumberOfClasses

           SedimentClass => Me%Classes(n)

                do j=WJLB, WJUB
                do i=WILB, WIUB
                    
                     WKUB = Me%KTop(i, j)

                    if (Me%ExternalVar%WaterPoints3D (i,j,WKUB) == WaterPoint) then
                        
                        Me%D50(i,j) = SedimentClass%D50*SedimentClass%Field3D(i, j, WKUB) + Me%D50(i,j)
                    endif
                enddo
                enddo 
        enddo do1
    
    end subroutine ComputeD50Cell
    
    !--------------------------------------------------------------------------
    !Compute the nondimensional shear stress
    subroutine ComputeNDShearStress
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
            
            if ((Me%ExternalVar%OpenPoints3D(i, j, WKUB) == OpenPoint)) then
                Me%NDShearStress(i,j) = Me%ExternalVar%ShearStress(i,j)/(Me%DeltaDensity*Gravity*Me%D50(i,j))
            endif
        enddo
        enddo
    end subroutine ComputeNDShearStress
     !--------------------------------------------------------------------------
    
    subroutine ComputeCriticalShearStress
        
        !Local-----------------------------------------------------------------
        integer             :: i, j, n
        real                :: CriticalShearStress
        class(T_Class), pointer :: SedimentClass
        integer             :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
do1:    do n=1,Me%NumberOfClasses

           SedimentClass => Me%Classes(n)

                do j=WJLB, WJUB
                do i=WILB, WIUB
                    
                    WKUB = Me%KTop(i, j)

                    if (Me%ExternalVar%OpenPoints3D(i, j, WKUB) == OpenPoint) then
                        
                        !critical Shields parameter modified by van Rijn (2003, 2007)
                       
                        Me%Dast(I,J) = Me%D50(i,j)*((Me%RelativeDensity-1)*Gravity/WaterCinematicVisc**2)**(1./3.)
                        
                        If (Me%Dast(I,J).LE.4) then
                           CriticalShearStress = 0.115*Me%Dast(I,J)**(-0.5)
                        elseIf (Me%Dast(I,J).GT.4.AND.Me%Dast(I,J).LE.10) then
                           CriticalShearStress = 0.14*Me%Dast(I,J)**(-0.64)
                        elseIf (Me%Dast(I,J).GT.10.AND.Me%Dast(I,J).LE.20) then
                           CriticalShearStress = 0.04*Me%Dast(I,J)**(-0.1)
                        elseIf (Me%Dast(I,J).GT.20.AND.Me%Dast(I,J).LE.150) then
                           CriticalShearStress = 0.013*Me%Dast(I,J)**0.29
                        elseIf (Me%Dast(I,J).GT.150) then
                           CriticalShearStress = 0.055
                        endif
                        
                        ! [N/m2]
                        SedimentClass%CriticalShearStress(i, j)= Me%DeltaDensity*Gravity*SedimentClass%D50      &
                                                                *CriticalShearStress   

                        SedimentClass%NDCriticalShearStress(i, j)= SedimentClass%CriticalShearStress(i, j)/     &
                                                                (Me%DeltaDensity*Gravity*Me%D50(i,j))
                    else

                        SedimentClass%CriticalShearStress(i, j)= FillValueReal

                    endif

                enddo
                enddo
              
        enddo do1

    end subroutine ComputeCriticalShearStress
    
      !--------------------------------------------------------------------------

    subroutine ComputeTransport           
    
            if (Me%TransportMethod==1) then
                call CurrentOnlyTransport
            elseif  (Me%TransportMethod==2) then
                 call CurrentPlusWavesTransport
            elseif (Me%TransportMethod==3) then
                 call WavesOnlyTransport
            endif

    end subroutine ComputeTransport
    !--------------------------------------------------------------------------

    subroutine CurrentOnlyTransport
    
        !Arguments-------------------------------------------------------------
        class(T_Class), pointer :: SedimentClass
        
        !Local-----------------------------------------------------------------
        real    ::  A,n1,p,DeltaTau,NDBedload
        integer :: i, j, n
        integer :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%Bedload(:, :) = 0.
        
do1:    do n=1,Me%NumberOfClasses
    
            SedimentClass => Me%Classes(n)
            
            A = SedimentClass%parameter_A
            n1 = SedimentClass%parameter_n
            p = SedimentClass%parameter_p
            
            SedimentClass%Bedload(:, :) = 0.
            
            do j=WJLB, WJUB
            do i=WILB, WIUB
                
                WKUB = Me%KTop(i, j)

                if (Me%ExternalVar%OpenPoints3D(i, j, WKUB) == OpenPoint) then
                    
                    DeltaTau = Me%NDShearStress(i, j)-SedimentClass%NDCriticalShearStress(i, j)                  
                    
                    if (DeltaTau.GT.0.) then
                        !Dimensionless bedload tranport rate
                        NDBedload = A*Me%NDShearStress(i, j)**n1*DeltaTau**p 
                        !Volumetric bedload transport rate per unit width [m3/s/m]
                        SedimentClass%BedLoad(i, j) = NDBedload*(gravity*(Me%RelativeDensity-1)*&
                            Me%D50(i, j)**3)**(1./2.)
                        
                        Me%Bedload(i, j) = Me%Bedload(i, j) + SedimentClass%Bedload(i, j)
                    endif
                endif
            enddo
            enddo
        enddo do1

    end subroutine CurrentOnlyTransport

    !--------------------------------------------------------------------------

    subroutine CurrentPlusWavesTransport
    
        !Arguments-------------------------------------------------------------
        class(T_Class), pointer :: SedimentClass
      
        !Local-----------------------------------------------------------------
!        real    :: 
        integer :: i,j, n

        !----------------------------------------------------------------------
do1:    do n=1,Me%NumberOfClasses
    
            SedimentClass => Me%Classes(n)
            
            do j=Me%WorkSize%JLB, Me%WorkSize%JUB
            do i=Me%WorkSize%ILB, Me%WorkSize%IUB
       
                if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then 
                               
                    !WavePeriod=Me%ExternalVar%WavePeriod(i,j)
                    !WaveHeight=Me%ExternalVar%WaveHeight(i,j)
                    !ksw       =Me%ExternalVar%WaveRugosity(i,j)
                    !Ksc       =Me%ExternalVar%CurrentRugosity(i,j)              
                    !D50       =SedimentClass%D50
                    !D90       =SedimentClass%D90
                    !VisCin    =WaterCinematicVisc
                    !VelMod    =Me%ExternalVar%VelMod(i,j)
                    !Ubw       =Me%ExternalVar%Ubw(i,j)
                    !Abw       =Me%ExternalVar%Abw(i,j)
                    !Dast      =Me%Dast(i,j)
                    !RhoWater  =Me%ExternalVar%WaterDensity
                    !VelQ      =FallVel(D50)
                endif
                
            enddo
            enddo
            
        enddo do1               
        
    end subroutine CurrentPlusWavesTransport

    !--------------------------------------------------------------------------
    
    subroutine WavesOnlyTransport
        !Arguments-------------------------------------------------------------
        class(T_Class), pointer :: SedimentClass
        
    end subroutine WavesOnlyTransport
    !--------------------------------------------------------------------------

      subroutine ComputeFluxes
              
        
        !Local-----------------------------------------------------------------                
        integer                 :: i, j, n
        real                    :: Xaux, Yaux, XYMod, FluxX
        real                    :: FluxY, TopLayerThickness
        class(T_Class), pointer :: SedimentClass
        integer                :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%FluxX  (:, :) =  0.
        Me%FluxY  (:, :) =  0.
        
        do n=1,Me%NumberOfClasses

            SedimentClass => Me%Classes(n)
            
            SedimentClass%FluxX  (:, :) =  0.
            SedimentClass%FluxY  (:, :) =  0.

            !Computes the Sediment fluxes (m3/s) in the middle of the cells

            do i=WILB, WIUB
            do j=WJLB, WJUB
                
                WKUB = Me%KTop(i, j)

                if (Me%ExternalVar%OpenPoints3D(i, j, WKUB) == OpenPoint) then                    
                    
                    TopLayerThickness = Me%ExternalVar%DWZ(i,j,WKUB)
                        
                    if (TopLayerThickness > Me%MinLayerThickness) then

                            Xaux  = Me%ExternalVar%VelU  (i,j)
                            Yaux  = Me%ExternalVar%VelV  (i,j)
                            XYMod = Me%ExternalVar%VelMod(i,j)
                        
                        if (XYMod > 0.) then
                            
                            !m
                            FluxX     = Xaux / XYMod * Me%ExternalVar%DVY(i, j)
                            FluxY     = Yaux / XYMod * Me%ExternalVar%DUX(i, j)
                            
                            !m3/s
                            SedimentClass%FluxX(i, j) = SedimentClass%Bedload(i, j) * SedimentClass%Field3D(i, j, WKUB) * FluxX     &
                                                        * Me%Aceleration%Coef
                                                    
                            SedimentClass%FluxY(i, j) = SedimentClass%Bedload(i, j) * SedimentClass%Field3D(i, j, WKUB) * FluxY     &
                                                        * Me%Aceleration%Coef
                            
                            Me%FluxX(i,j) = Me%FluxX(i,j) + SedimentClass%FluxX(i, j)
                            
                            Me%FluxY(i,j) = Me%FluxY(i,j) + SedimentClass%FluxY(i, j)
                   
                        endif

                    endif

                endif

            enddo
            enddo
        enddo

      end subroutine ComputeFluxes
      
   !--------------------------------------------------------------------------   
      
    subroutine ComputeEvolution 
        
        !Local-----------------------------------------------------------------
        real                    :: DT_Area1, DT_Area2, AbsFluxXY, RunPeriod
        integer                 :: i, j, di, dj, n
        class(T_Class), pointer :: SedimentClass      
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
      
             
do1:    do n=1,Me%NumberOfClasses

           SedimentClass => Me%Classes(n)
           
           SedimentClass%DZ(:, :) = 0.
            
            if (Me%Residual%ON) then

                RunPeriod = Me%ExternalVar%Now- Me%Residual%StartTime
                
                Me%Residual%FluxX(:,:) = ( Me%Residual%FluxX(:,:) * (RunPeriod -  Me%Evolution%SedimentDT)          + &
                                           SedimentClass%FluxX(:,:) * Me%Evolution%SedimentDT) / RunPeriod
                                           
                Me%Residual%FluxY(:,:) = ( Me%Residual%FluxY(:,:) * (RunPeriod -  Me%Evolution%SedimentDT)          + &
                                           SedimentClass%FluxY(:,:) * Me%Evolution%SedimentDT) / RunPeriod
            endif
           
            do j=WJLB, WJUB
            do i=WILB, WIUB


                if (SedimentClass%FluxX(i, j) < 0.) then 
                    if      (Me%ExternalVar%ComputeFacesU2D(i,  j) == Covered ) then

                        DT_Area1 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j-1) / Me%ExternalVar%DVY(i, j-1)
                        DT_Area2 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)
                        
                        SedimentClass%DZ(i, j-1) = SedimentClass%DZ(i, j-1) - DT_Area1 * SedimentClass%FluxX(i, j)
                        SedimentClass%DZ(i, j  ) = SedimentClass%DZ(i, j  ) + DT_Area2 * SedimentClass%FluxX(i, j)
                         

                        if (Me%Boxes%Yes) then
                            Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j) + SedimentClass%FluxX(i, j) * Me%Evolution%SedimentDT
                        endif
                    endif
                else 

                    if (Me%ExternalVar%ComputeFacesU2D(i,j+1) == Covered) then

                        DT_Area1 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j+1) / Me%ExternalVar%DVY(i, j+1)
                        DT_Area2 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)
                  
                        SedimentClass%DZ(i, j+1) = SedimentClass%DZ(i, j+1) + DT_Area1 * SedimentClass%FluxX(i, j)  
                        SedimentClass%DZ(i, j  ) = SedimentClass%DZ(i, j  ) - DT_Area2 * SedimentClass%FluxX(i, j) 

                        if (Me%Boxes%Yes) then
                            Me%Boxes%FluxesX(i,j+1) = Me%Boxes%FluxesX(i,j+1) + SedimentClass%FluxX(i, j) * Me%Evolution%SedimentDT
                        endif
                    endif                    
                endif

                if (SedimentClass%FluxY(i, j) < 0.) then
                    if  (Me%ExternalVar%ComputeFacesV2D(i,   j) == Covered) then

                        DT_Area1 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i-1, j) / Me%ExternalVar%DVY(i-1, j)
                        DT_Area2 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                        SedimentClass%DZ(i-1, j) = SedimentClass%DZ(i-1, j) - DT_Area1 * SedimentClass%FluxY(i, j) 
                        SedimentClass%DZ(i  , j) = SedimentClass%DZ(i  , j) + DT_Area2 * SedimentClass%FluxY(i, j) 

                        if (Me%Boxes%Yes) then
                            Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j) + SedimentClass%FluxY(i, j) * Me%Evolution%SedimentDT
                        endif
                    endif
                else 
                    if (Me%ExternalVar%ComputeFacesV2D(i+1, j) == Covered) then

                        DT_Area1 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i+1, j) / Me%ExternalVar%DVY(i+1, j)
                        DT_Area2 = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                        SedimentClass%DZ(i+1, j) = SedimentClass%DZ(i+1, j) + DT_Area1 * SedimentClass%FluxY(i, j)
                        SedimentClass%DZ(i  , j) = SedimentClass%DZ(i  , j) - DT_Area2 * SedimentClass%FluxY(i, j) 

                        if (Me%Boxes%Yes) then
                            Me%Boxes%FluxesY(i+1,j) = Me%Boxes%FluxesY(i+1,j) + SedimentClass%FluxY(i, j) * Me%Evolution%SedimentDT
                        endif
                    endif
                endif
                
     
            enddo
            enddo
            
        enddo do1
   
    end subroutine ComputeEvolution
    
    !--------------------------------------------------------------------------

    subroutine ComputeDischarges            

        !Local-----------------------------------------------------------------
        real                               :: DT_Area, BottomDensity, TicknessAdd,      &
                                              DischargeFlow, DischargeConc
        integer                            :: i, j, dis, DischargesNumber, STAT_CALL
        integer                            :: n
        class(T_Class), pointer            :: SedimentClass 
        !----------------------------------------------------------------------

do1:    do n=1,Me%NumberOfClasses

           SedimentClass => Me%Classes(n)

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



                call GetDischargeConcentration (Me%ObjDischarges,                       &
                                                Me%ExternalVar%Now,                     &
                                                dis, DischargeConc,                     &
                                                PropertyIDNumber = SedimentClass%ID%IDNumber,         &        
                                                STAT = STAT_CALL)
                if (STAT_CALL/=SUCCESS_)                                                &
                    stop 'ComputeDischarges - ModuleSediment - ERR40'
                

                BottomDensity       = Me%Porosity * Me%ExternalVar%WaterDensity + (1. - Me%Porosity) * Me%Density

                DT_Area             = Me%Evolution%SedimentDT / Me%ExternalVar%DUX(i, j  ) / Me%ExternalVar%DVY(i, j)

                                     ![m3/s] * [g/m3/1000] / [kg/m3] * [s/m2] * []
                TicknessAdd         = DischargeFlow * (DischargeConc/1000.) / BottomDensity * DT_Area

                SedimentClass%DZ(i, j) = SedimentClass%DZ(i, j) + TicknessAdd

            enddo d1
    
        enddo do1
        
    end subroutine ComputeDischarges            

    
    !--------------------------------------------------------------------------


    subroutine BoundaryCondition

        !Local-----------------------------------------------------------------
        real                               :: a1, a2, a3, a4, atotal
        integer                            :: i, j, ILB, IUB, JLB, JUB
        integer                            :: n
        class(T_Class), pointer            :: SedimentClass 
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 

        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        
do1:    do n=1,Me%NumberOfClasses

            SedimentClass => Me%Classes(n)


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
                    
                            SedimentClass%DZ(i, j) = (a1 * SedimentClass%DZ(i-1, j) + a2 * SedimentClass%DZ(i+1, j)  +     &
                                             a3 * SedimentClass%DZ(i, j-1) + a4 * SedimentClass%DZ(i, j+1)) / atotal
                        endif
                                          

                    endif

                enddo
                enddo    

            else if (Me%Boundary == Cyclic) then

                where (Me%ExternalVar%BoundaryPoints2D(ILB, :) == Boundary) 
                    SedimentClass%DZ(ILB, :) = SedimentClass%DZ(IUB-1, :)
                end where

                where (Me%ExternalVar%BoundaryPoints2D(IUB, :) == Boundary) 
                    SedimentClass%DZ(IUB, :) = SedimentClass%DZ(ILB+1, :)
                end where

                where (Me%ExternalVar%BoundaryPoints2D(: ,JLB) == Boundary) 
                    SedimentClass%DZ(:, JLB) = SedimentClass%DZ(:, JUB-1)
                end where

                where (Me%ExternalVar%BoundaryPoints2D(:, JUB) == Boundary) 
                    SedimentClass%DZ(:, JUB) = SedimentClass%DZ(:, JLB+1)
                end where
            else if (Me%Boundary == NullValue) then

                where (Me%ExternalVar%BoundaryPoints2D(:, :) == Boundary) 
                    SedimentClass%DZ(:, :) = 0.
                end where            

            endif
        enddo do1        
    
    end subroutine BoundaryCondition
    
    !--------------------------------------------------------------------------
    subroutine ComputeTotalDZ
    
        !Local-----------------------------------------------------------------
        integer                 :: i, j, n
        class(T_Class), pointer :: SedimentClass 
        integer                 :: WILB, WIUB, WJLB, WJUB, WKUB
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
       
        Me%DZ(:, :) = 0.
        
do1:    do n=1,Me%NumberOfClasses

          SedimentClass => Me%Classes(n)
        
          do j=WJLB, WJUB
          do i=WILB, WIUB
              
              WKUB = Me%KTop(i, j)
              
             if (Me%ExternalVar%OpenPoints2D(i, j) == OpenPoint) then
                 
                Me%DZ(i, j) = Me%DZ(i, j) + SedimentClass%DZ(i, j)
                
                 if  (WKUB > 0 .and. Me%DZ(i, j)  < -Me%ExternalVar%DWZ(i,j,WKUB)) then
                    write(*,*) 'Sediment degradation is larger than layer thickness.'
                    write(*,*) i, j, WKUB, Me%DZ(i, j)
                    write(*,*) 'Increase minimum layer thickness.'
                    stop 'ComputeTotalDZ - ModuleSediment - ERR01'
                endif
                
             endif
          
          enddo
          enddo
          
        enddo do1
    end subroutine ComputeTotalDZ
    
     !--------------------------------------------------------------------------
     subroutine ComputePercentage
        !Local-----------------------------------------------------------------
        integer                            :: i, j, k, n, WKUB
        integer                            :: WILB, WIUB, WJLB, WJUB
        class(T_Class), pointer            :: SedimentClass 
        real(8)                            :: DZ
         !----------------------------------------------------------------------
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        Me%TotalPercentage (:,:,:) = 0.
        
do1:    do n=1,Me%NumberOfClasses

          SedimentClass => Me%Classes(n)
        
          do j=WJLB, WJUB
          do i=WILB, WIUB
              
              WKUB = Me%KTop(i, j)
              
            if (Me%ExternalVar%OpenPoints3D(i, j, WKUB) == OpenPoint) then                
                
                SedimentClass%Field3D(i, j, WKUB) = (Me%ExternalVar%DWZ(i,j,WKUB)*SedimentClass%Field3D(i, j, WKUB) + SedimentClass%DZ(i, j)) &
                                                    /(Me%ExternalVar%DWZ(i,j,WKUB)+Me%DZ(i,j))
                
                if  (SedimentClass%Field3D(i, j, WKUB)  > 1.001) then
                         write(*,*) 'The class percentage is larger than 100%.'
                         write(*,*) i, j, WKUB, SedimentClass%Field3D(i, j, WKUB), 'n=',n
                         !stop 'ComputePercentage - ModuleSediment - ERR01'                 
                    elseif (SedimentClass%Field3D(i, j, WKUB)  < 0.) then
                         write(*,*) 'The class percentage is smaller than 0%.'
                         write(*,*) i, j, WKUB, SedimentClass%Field3D(i, j, WKUB), 'n=',n
                         !stop 'ComputePercentage - ModuleSediment - ERR02'
                    endif
                
                SedimentClass%TopPercentage(i, j) = SedimentClass%Field3D(i, j, WKUB)
               
                Me%TotalPercentage (i, j, WKUB) = Me%TotalPercentage (i, j, WKUB) + SedimentClass%Field3D(i, j, WKUB)
                
            endif
                    
          enddo
          enddo                               
                                                 
        enddo do1

    
        do j=WJLB, WJUB
        do i=WILB, WIUB
            
            WKUB = Me%KTop(i, j)
              
              if (Me%ExternalVar%OpenPoints3D(i, j, WKUB) == OpenPoint) then          
                  
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
    
    subroutine ComputeVerticalCoordinate          

        !Local-----------------------------------------------------------------
        integer                             :: i, j, k, WKUB, TotalWKUB
        real(8)                             :: DZ, TopLayerThickness
        real(8)                             :: ExcessDZ, ExcessMin
        integer                             :: WILB, WIUB, WJLB, WJUB, n
        class(T_Class), pointer             :: SedimentClass
        !----------------------------------------------------------------------
        
        WILB = Me%SedimentWorkSize3D%ILB
        WIUB = Me%SedimentWorkSize3D%IUB
        WJLB = Me%SedimentWorkSize3D%JLB
        WJUB = Me%SedimentWorkSize3D%JUB
        
        TotalWKUB = Me%SedimentWorkSize3D%KUB
           
        do j=WJLB, WJUB
        do i=WILB, WIUB

            if (Me%ExternalVar%OpenPoints2D (i ,j) == OpenPoint) then
                    
                WKUB = Me%KTop(i, j)
                    
                DZ = Me%DZ(i, j)
                    
                TopLayerThickness = Me%ExternalVar%DWZ(i,j,WKUB)
                    
                if (DZ > 0.)then !aggradation

                    ExcessDZ = TopLayerThickness + DZ - Me%MaxLayerThickness

                    if(ExcessDZ > 0.)then !New Layer

                        Me%KTop(i,j) = WKUB+1
                        
                        ExcessMin = Me%MinLayerThickness - ExcessDZ
                        
                        !The new top layer has initially the same class percentages than the older one 
                        do n=1,Me%NumberOfClasses 
                            
                           SedimentClass => Me%Classes(n)
                                
                           SedimentClass%Field3D(i, j, WKUB+1) = SedimentClass%Field3D(i, j, WKUB)
                                
                        enddo
                        
                        if(Me%KTop(i,j) > TotalWKUB)then  !Maximum number of layers exceeded
                            
                            write(*,*) 'Maximum number of layers exceeded in cell i=',i,'j=',j
                            write(*,*) 'Last sediment layer deleted'
                            write(*,*) 'ComputeVerticalCoordinate - ModuleSediment - WRN10' 
                                                                                          
                            
                            Me%VerticalCoordinate(i,j,WKUB) = Me%VerticalCoordinate(i,j,WKUB-1)     - &
                                                              Me%MaxLayerThickness + ExcessMin
                            
                            Me%VerticalCoordinate(i,j,WKUB+1) = Me%VerticalCoordinate(i,j,WKUB)     - &
                                                              Me%MinLayerThickness
                                    
                            do k=0,TotalWKUB !Last sediment layer deleted to maintain the number of layers constant 
                                    
                                Me%VerticalCoordinate(i,j,k) = Me%VerticalCoordinate(i,j,k+1)
                                
                                do n=1,Me%NumberOfClasses !Class percentages in the cells updated 
                                   
                                   SedimentClass => Me%Classes(n)
                                    
                                   SedimentClass%Field3D(i,j,k) = SedimentClass%Field3D(i,j,k+1)
                                enddo
                                
                            enddo 
                            
                            Me%KTop(i,j) = TotalWKUB 
                                
                        else
                            
                            Me%VerticalCoordinate(i,j,WKUB)             = Me%VerticalCoordinate(i,j,WKUB-1)     - &
                                                                        Me%MaxLayerThickness + ExcessMin

  
                            Me%VerticalCoordinate(i,j,WKUB+1:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB)       - &
                                                                        Me%MinLayerThickness
                        
                        endif                                
                        
                    else !Increase layer thickness
                        
                        if (Me%KTop(i,j) < 1) then !Activate layer
                            
                            Me%KTop(i,j) = 1
                            
                            do n=1,Me%NumberOfClasses

                                SedimentClass => Me%Classes(n)        
                
                                WKUB = Me%KTop(i, j)
                
                                SedimentClass%Field3D(i, j, WKUB) = SedimentClass%DZ(i, j)/DZ
                            enddo
                            
                        endif
                        
                         Me%VerticalCoordinate(i,j,WKUB:TotalWKUB)   = Me%VerticalCoordinate(i,j,WKUB) - DZ

                    end if


            elseif(DZ <  0.)then !degradation

                ExcessDZ = Me%ExternalVar%DWZ(i,j,WKUB) - Me%MinLayerThickness

                if(abs(DZ) .ge. ExcessDZ)then !Layer eroded
                    
                    if(WKUB > 1)then
                    
                    Me%VerticalCoordinate(i,j,WKUB-1) = Me%VerticalCoordinate(i,j,WKUB-1)       &
                                                        - Me%ExternalVar%DWZ(i,j,WKUB) - DZ
                                                            
                    Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB-1)

                        
                        do n=1,Me%NumberOfClasses

                            SedimentClass => Me%Classes(n)        
                
                            SedimentClass%Field3D(i, j, WKUB-1) = (Me%ExternalVar%DWZ(i,j,WKUB-1)*SedimentClass%Field3D(i,j,WKUB-1) +     &
                                                                  (Me%ExternalVar%DWZ(i,j,WKUB) + DZ)*SedimentClass%Field3D(i,j,WKUB))    &
                                                                  /(Me%ExternalVar%DWZ(i,j,WKUB-1) + Me%ExternalVar%DWZ(i,j,WKUB) + DZ)                       
                        enddo
                       
                    else !(WKUB = 1)
                         
                        !write(*,*) 'Eroded all sediment layers in cell i=',i,'j=',j
                        !write(*,*) 'ComputeVerticalCoordinate - ModuleSediment - WRN11'
                        
                        Me%VerticalCoordinate(i,j,WKUB) =  Me%VerticalCoordinate(i,j,WKUB-1)       &
                                                             - Me%MinLayerThickness
                        
                        Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB)
                        
                        Me%DZ(i, j) = -ExcessDZ 
                            
                    endif
                    
                    Me%KTop(i,j)      = WKUB-1

                else

                    Me%VerticalCoordinate(i,j,WKUB:TotalWKUB) = Me%VerticalCoordinate(i,j,WKUB) - DZ

                endif

            end if

            end if

        enddo
        enddo    
    
    end subroutine ComputeVerticalCoordinate      
    
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
                    
                   Me%DZ_Residual(i, j)    =Me%DZ_Residual(i, j) + Me%DZ(i, j)

                endif

            enddo
            enddo    

            if (Me%Evolution%BathymDT >= Me%Evolution%SedimentDT) then

                do j=WJLB, WJUB
                do i=WILB, WIUB

                    if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
                    
                       Me%BatimIncrement(i,j)  =Me%BatimIncrement(i,j) + Me%DZ(i, j)

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
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR01'

                call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",     &
                                     Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR02'
                
                !Write OpenPoints2D
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                  &
                                     WorkJUB, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR03'

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints2D", "OpenPoints2D",       &
                                     "-", Array2D = Me%ExternalVar%OpenPoints2D,             &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR04'
                

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Bathymetry", "Bathymetry",        &
                                     "-", Array2D = Me%ExternalVar%Bathymetry,               &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR05'
       
                call HDF5WriteData  (Me%ObjHDF5, "/Results/DZ", "DZ",                        &
                                     "m", Array2D = Me%DZ,                                   &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR06'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/D50", "D50",                      &
                                     "m", Array2D = Me%D50,                                  &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR07'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/DZ_Residual", "DZ_Residual",      &
                                     "m", Array2D =Me%DZ_Residual,                           &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR08'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Bedload", "Bedload",              &
                                    "m3/s/m", Array2D = Me%Bedload,                          &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR09'
                
                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux X", "Transport Flux X",   &
                                    "m3/s/m", Array2D = Me%FluxX,                                   &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR10'

                call HDF5WriteData  (Me%ObjHDF5, "/Results/Transport Flux Y", "Transport Flux Y",   &
                                    "m3/s/m", Array2D = Me%FluxY,                                   &
                                    OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR11'
                 

do1:            do n=1,Me%NumberOfClasses

                    call HDF5WriteData   (Me%ObjHDF5, "/Results/Classes/TopPercentage/"//trim(Me%Classes(n)%ID%Name),      &
                                          trim(Me%Classes(n)%ID%Name),                                                     &
                                          "%", Array2D = Me%Classes(n)%TopPercentage,                                      &
                                          OutputNumber = OutPutNumber,                                                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR12'
                
                    call HDF5WriteData  (Me%ObjHDF5, "/Results/Classes/Bedload/"//trim(Me%Classes(n)%ID%Name),            &
                                        trim(Me%Classes(n)%ID%Name),                                                        &
                                        "m3/s/m", Array2D = Me%Classes(n)%Bedload,                                &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR013'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/Classes/Critical Shear Stress/"//trim(Me%Classes(n)%ID%Name),                    &
                                        trim(Me%Classes(n)%ID%Name),                                                        &
                                         "N/m2", Array2D = Me%Classes(n)%CriticalShearStress,                                       &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR14'
                    

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/Classes/Transport Flux X/"//trim(Me%Classes(n)%ID%Name),              &
                                         trim(Me%Classes(n)%ID%Name),                                                       &
                                         "m3/s/m", Array2D = Me%Classes(n)%FluxX,                                           &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR15'

                    call HDF5WriteData  (Me%ObjHDF5, "/Results/Classes/Transport Flux Y/"//trim(Me%Classes(n)%ID%Name),              &
                                        trim(Me%Classes(n)%ID%Name),                                                        &
                                         "m3/s/m", Array2D = Me%Classes(n)%FluxY,                                           &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR16'
                    
                    call HDF5WriteData   (Me%ObjHDF5, "/Results3D/"//trim(Me%Classes(n)%ID%Name),     &
                                         trim(Me%Classes(n)%ID%Name),                                                        &
                                         "%", Array3D = Me%Classes(n)%Field3D,                                               &
                                         OutputNumber = OutPutNumber,                                                        &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleSediment - ERR12'
                    
               
                enddo do1
                    
                if (Me%Residual%ON) then
                
                    call HDF5WriteData  (Me%ObjHDF5, "/Residual/Transport Flux X", &
                                         "Transport Flux X",                      &
                                         "m3/s/m", Array2D = Me%Residual%FluxX,            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR17'

                    call HDF5WriteData  (Me%ObjHDF5, "/Residual/Transport Flux Y", &
                                         "Transport Flux Y",                      &
                                         "m3/s/m", Array2D = Me%Residual%FluxY,            &
                                         OutputNumber = OutPutNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR18'
                
                endif

                !Writes everything to disk
                call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'OutPutSedimentHDF - ModuleSediment - ERR19'

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

        call WriteTimeSerie(Me%ObjTimeSerie,                     &
                            Data2D = Me%DZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR01'

        call WriteTimeSerie(Me%ObjTimeSerie,                     &
                            Data2D =Me%DZ_Residual, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                               &
            stop 'OutPut_TimeSeries - ModuleSediment - ERR02'

    
    end subroutine OutPut_TimeSeries

    !--------------------------------------------------------------------------


    subroutine OutputBoxFluxes


        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: i, j
        integer                                 :: STAT_CALL
        integer                                 :: n
        class(T_Class), pointer                 :: SedimentClass

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        
 do1:    do n=1,Me%NumberOfClasses

            SedimentClass => Me%Classes(n)

            Me%Boxes%Mass(:,:) = 0.

            do J = JLB, JUB
            do I = ILB, IUB
                Me%Boxes%Mass   (i,j) = SedimentClass%DZ(i,j)      * (1.- Me%Porosity) * Me%Density * &
                                        Me%ExternalVar%DUX(i,j) * Me%ExternalVar%DVY(i,j) 

                !This fluxes are initialised and partial computed in the subroutine ComputeEvolution
                Me%Boxes%FluxesX(i,j) = Me%Boxes%FluxesX(i,j)   * (1.- Me%Porosity) * Me%Density
                Me%Boxes%FluxesY(i,j) = Me%Boxes%FluxesY(i,j)   * (1.- Me%Porosity) * Me%Density

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
                        deallocate(Me%Classes(i)%Field3D)
                        nullify(Me%Classes(i)%Field3D)
                        deallocate(Me%Classes(i)%CriticalShearStress)
                        nullify(Me%Classes(i)%CriticalShearStress)
                        deallocate(Me%Classes(i)%NDCriticalShearStress)
                        nullify(Me%Classes(i)%NDCriticalShearStress)
                        deallocate(Me%Classes(i)%DZ)
                        nullify(Me%Classes(i)%DZ)
                        deallocate(Me%Classes(i)%FluxX)
                        nullify(Me%Classes(i)%FluxX)
                        deallocate(Me%Classes(i)%FluxY) 
                        nullify(Me%Classes(i)%FluxY) 
                        deallocate(Me%Classes(i)%Bedload)
                        nullify(Me%Classes(i)%Bedload)
                        deallocate(Me%Classes(i)%TopPercentage)
                        nullify(Me%Classes(i)%TopPercentage)
                    enddo do1
                    
                    deallocate(Me%Classes)
                    nullify(Me%Classes)
                    deallocate(Me%TotalPercentage)
                    nullify(Me%TotalPercentage)
                
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
                
                deallocate (Me%D50, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &            
                    stop 'KillSediment - ModuleSediment. ERR01.' 
                nullify (Me%D50) 
        
                deallocate (Me%Dast, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &           
                    stop 'KillSediment - ModuleSediment. ERR02.' 
                nullify (Me%Dast) 
                        
                deallocate (Me%NDShearStress, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR03.' 
                nullify (Me%NDShearStress) 
        
                deallocate(Me%Elevation, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR04.' 
                nullify (Me%Elevation) 
            
                deallocate(Me%VerticalCoordinate, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR05.' 
                nullify (Me%VerticalCoordinate) 

                deallocate(Me%KTop, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_)                                          &          
                    stop 'KillSediment - ModuleSediment. ERR06.' 
                nullify (Me%KTop) 
                
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
            call HDF5WriteData  (Me%ObjHDF5, "/Results",              &
                                 "Transport Flux X",                           &
                                 "m3/s/m",                                              &
                                 Array2D = Me%Residual%FluxX,                           &
                                 STAT    = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR40'

            call HDF5WriteData  (Me%ObjHDF5, "/Results",              &
                                 "Transport Flux Y",                           &
                                 "m3/s/m",                                              &
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

            call HDF5WriteData  (Me%ObjHDF5, "/Results",                               &
                                 "Start Time", "YYYY/MM/DD HH:MM:SS",          &
                                 Array1D = TimePtr, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFinalState - ModuleSediment - ERR80'
                                 

        endif
        
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB,                               &
                             WorkJLB, WorkJUB, WorkKLB, WorkKUB,                         &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            stop 'WriteFinalState - ModuleSediment - ERR13'
            
        do n=1,Me%NumberOfClasses
            
            call HDF5WriteData   (Me%ObjHDF5, "/Results3D",      &
                                    trim(Me%Classes(n)%ID%Name),                        &
                                    "%", Array3D = Me%Classes(n)%Field3D,                 &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'WriteFinalState - ModuleSediment - ERR14'
            
        enddo
   
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











