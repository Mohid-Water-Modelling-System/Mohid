!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Waves
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to compute waves processes
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

Module ModuleWaves

    use ModuleGlobalData
    use ModuleFunctions,        only : Secant
    use ModuleEnterData        
    use ModuleTime
    use ModuleHorizontalMap,    only : GetWaterPoints2D, GetOpenPoints2D, UnGetHorizontalMap
    use ModuleHorizontalGrid,   only : LocateCell, GetHorizontalGridSize, GetHorizontalGrid,    &
                                       GetGridAngle, GetCheckDistortion, GetCoordTypeList,      &
                                       GetGridCoordType,  GetGridOrigin, WriteHorizontalGrid,   &
                                       UnGetHorizontalGrid, GetXYCellZ
    use ModuleFillMatrix,       only : ConstructFillMatrix, ModifyFillMatrix,  &
                                       GetIfMatrixRemainsConstant, KillFillMatrix 
    use ModuleGeometry,         only : GetGeometryWaterColumn, UnGetGeometry
    use ModuleHDF5,             only : ConstructHDF5, HDF5SetLimits, HDF5WriteData, &
                                       HDF5FlushMemory, GetHDF5FileAccess, KillHDF5
    use ModuleGridData,         only : GetGridData, UngetGridData   
    use ModuleTimeSerie,        only : StartTimeSerie, WriteTimeSerie, KillTimeSerie,   &
                                       GetTimeSerieLocation, CorrectsCellsTimeSerie,    &
                                       GetNumberOfTimeSeries, TryIgnoreTimeSerie
    use ModuleDrawing         

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartWaves
    private ::      AllocateInstance
    private ::      ReadWavesFilesName
    private ::      ConstructGlobalVariables
    private ::      ConstructFetch
    private ::         ComputeFetchDistanceGrid
    private ::         ComputeFetchDistancePolygon
    private ::         ConstructLandArea
    private ::         ComputeFetch
    private ::      ConstructWaveParameters   
    private ::          ReadWaveParameters
    private ::          ConstructGlobalOutput
    private ::          Construct_Time_Serie      
    private ::      Open_HDF5_OutPut_File
      


    !Selector
    public  :: GetWaves
    public  :: GetWavesStress
    public  :: GetWavesOptions
    public  :: UnGetWaves
    
    public  :: SetWavesWind
    public  :: SetGeneric4DValues
                    
    
    !Modifier
    public  :: ModifyWaves
    public  ::      ComputeWaveParameters
    private ::      ComputeWaveHeightSimple  
    private ::      ComputeWaveHeightFetch   
    private ::      ComputeWavePeriodSimple  
    private ::      ComputeWavePeriodFetch   
    private ::      ComputeWaveDirection
    private ::      Output_Results_HDF
    private ::      Output_TimeSeries
    public  :: ComputeRadiationStress


    !Destructor
    public  :: KillWaves                                                     
    private ::      DeAllocateInstance
    
    !Management
    private ::      Ready
    private ::          LocateObjWaves 
    
    
    integer, parameter                                      :: Old          = 0
    integer, parameter                                      :: CEQUALW2     = 1 
    integer, parameter                                      :: Grid         = 0
    integer, parameter                                      :: Polygon      = 1

    integer, parameter                                      :: Sxx_         = 1
    integer, parameter                                      :: Syy_         = 2
    integer, parameter                                      :: Sxy_         = 3
    integer, parameter                                      :: Syx_         = 4

        
    !Interfaces----------------------------------------------------------------
                                                        
    !Types---------------------------------------------------------------------
                                                        
    type       T_External                               
        real,    dimension(:,:), pointer                    :: WindVelocityX
        real,    dimension(:,:), pointer                    :: WindVelocityY
        integer, dimension(:,:), pointer                    :: WaterPoints2D
        integer, dimension(:,:), pointer                    :: OpenPoints2D
        real,    dimension(:,:), pointer                    :: WaterColumn
        real,    dimension(:,:), pointer                    :: XX_IE, YY_IE
        real,    dimension(:,:), pointer                    :: DUX, DVY   
        real,    dimension(:,:), pointer                    :: Bathymetry
        real                                                :: GridAngle
        logical                                             :: DistortionOn         = .false.
        real                                                :: CurrentValue4D       = FillValueReal
        logical                                             :: Backtracking
    end type   T_External                               
                                                        
    type       T_OutPut                                 
         type (T_Time), pointer, dimension(:)               :: OutTime
         integer                                            :: NextOutPut
         logical                                            :: TimeSerie            = .false.
         logical                                            :: HDF                  = .false.
         integer                                            :: Number
    end type T_OutPut                                   
                                                        
    type T_WaveProperty                                 
        type(T_PropertyID)                                  :: ID
        logical                                             :: ON                   = .false.
        logical                                             :: Constant             = .false.
        integer                                             :: Source
        real, dimension(:,:),  pointer                      :: Field
        logical                                             :: OutputHDF            = .false.
        logical                                             :: TimeSerieOn          = .false.
    end type T_WaveProperty

    type       T_Waves
        integer                                             :: InstanceID
        character(PathLength)                               :: ModelName
        type(T_Time  )                                      :: BeginTime
        type(T_Time  )                                      :: EndTime
        type(T_Time  )                                      :: ActualTime
        type(T_Time  )                                      :: LastCompute
        type(T_Size2D)                                      :: Size, WorkSize
        type(T_External)                                    :: ExternalVar
        type (T_WaveProperty)                               :: RadiationStressX
        type (T_WaveProperty)                               :: RadiationStressY
        type (T_WaveProperty)                               :: WavePeriod
        type (T_WaveProperty)                               :: WaveHeight
        type (T_WaveProperty)                               :: WaveDirection
        real, dimension(:,:  ), pointer                     :: Abw, Ubw, WaveLength
        real, dimension(:,:  ), pointer                     :: En, k
        real, dimension(:,:,:), pointer                     :: Distance, Fetch
        real, dimension(:    ), pointer                     :: AngleList  
        real, dimension(:    ), pointer                     :: XX, YY 
        real, dimension(:    ), pointer                     :: MinDistanceFromPoint                        
        character(Len=PathLength)                           :: FileName
        character(Len=PathLength)                           :: OutputFile
        logical                                             :: ParametersON            = .false.
        integer                                             :: ObjTime                 = 0
        integer                                             :: ObjHorizontalGrid       = 0
        integer                                             :: ObjGeometry             = 0
        integer                                             :: ObjHorizontalMap        = 0
        integer                                             :: ObjEnterData            = 0
        integer                                             :: ObjHDF5                 = 0
        integer                                             :: ObjGridData             = 0
        integer                                             :: ObjTimeSerie            = 0
        integer                                             :: WaveGen_type            = Old
        integer                                             :: TotalDirections         = 8
        integer                                             :: DistanceType            = Grid
        integer                                             :: N, NNE, NE, ENE, E, ESE 
        integer                                             :: SE, SSE, S, SSW, SW, WSW 
        integer                                             :: W, WNW, NW, NNW
        real                                                :: WaveHeightParameter     = 1.
        real                                                :: WavePeriodParameter     = 1.
        type (T_OutPut)                                     :: OutPut
        type(T_Waves), pointer                              :: Next
        type(T_PointF), pointer                             :: Point
        type(T_Polygon), pointer                            :: LandArea
    end type   T_Waves
    
    
    !Global Module Variables
    type (T_Waves), pointer                                 :: FirstObjWaves
    type (T_Waves), pointer                                 :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
   !-----------------------------------------------------------------------
   
    subroutine StartWaves(ModelName, WavesID, TimeID, HorizontalMapID, HorizontalGridID, &
                          GridDataID, GeometryID, STAT)

        !Arguments---------------------------------------------------------------
        character(Len=*)                                :: ModelName
        integer                                         :: WavesID
        integer                                         :: TimeID
        integer                                         :: HorizontalMapID
        integer                                         :: HorizontalGridID
        integer                                         :: GridDataID
        integer                                         :: GeometryID
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         
      
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !Begin ------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mWaves_)) then
            nullify (FirstObjWaves)
            call RegisterModule (mWaves_) 
        endif

        call Ready(WavesID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ModelName = ModelName

            !External Modules
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID           )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID  )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_,       GridDataID       )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID       )
                        
                        
            call ReadWavesFilesName

            call ConstructGlobalVariables

            call ConstructWaveParameters
            
            if (Me%OutPut%HDF) call Open_HDF5_OutPut_File

            call ConstructFetch

            !Returns ID
            WavesID     = Me%InstanceID
            STAT_       = SUCCESS_
            
        else cd0
            
            stop 'ModuleWaves - StartWaves - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        
    end subroutine StartWaves
 
    !--------------------------------------------------------------------------
    

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Waves), pointer                         :: NewObjWaves
        type (T_Waves), pointer                         :: PreviousObjWaves

        !Begin-----------------------------------------------------------------

        !Allocates new instance
        allocate (NewObjWaves)
        nullify  (NewObjWaves%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjWaves)) then
            FirstObjWaves         => NewObjWaves
            Me                    => NewObjWaves
        else
            PreviousObjWaves      => FirstObjWaves
            Me                    => FirstObjWaves%Next
            do while (associated(Me))
                PreviousObjWaves  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjWaves
            PreviousObjWaves%Next => NewObjWaves
        endif

        Me%InstanceID = RegisterNewInstance (mWaves_)


    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------
    
    subroutine ReadWavesFilesName

        !External--------------------------------------------------------------
        integer                     :: STAT_CALL
      
        !Begin-----------------------------------------------------------------

        call ReadFileName('WAVES_DAT', Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL/=SUCCESS_) stop 'ReadWavesFilesName - ModuleWaves - ERR01'
                                
        call ReadFileName('WAVES_HDF', Me%OutputFile, STAT = STAT_CALL)               
        if (STAT_CALL/=SUCCESS_) stop 'ReadWavesFilesName - ModuleWaves - ERR10'
                                             
    end subroutine ReadWavesFilesName

    !--------------------------------------------------------------------------
    
    subroutine ConstructGlobalVariables

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
      
        
        call GetHorizontalGridSize(Me%ObjHorizontalGrid,                                &
                                   Size        = Me%Size,                               &
                                   WorkSize    = Me%WorkSize,                           &
                                   STAT        = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalVariables - ModuleWaves - ERR10'

        !Gets time
        call GetComputeTimeLimits(Me%ObjTime,                                           &
                                  BeginTime = Me%BeginTime,                             &
                                  EndTime   = Me%EndTime,                               &
                                  STAT      = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGlobalVariables - ModuleWaves - ERR20'

        !Actualize the time
        Me%ActualTime = Me%BeginTime  
        
        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ObjTime, Me%ExternalVar%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGlobalVariables - ModuleAtmosphere - ERR20'        
        

    end subroutine ConstructGlobalVariables
    
    !--------------------------------------------------------------------------

    subroutine ConstructWaveParameters

        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL, iflag

        !Begin-----------------------------------------------------------------
        
        
        !WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D,        &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructWaveParameters - ModuleWaves - ERR01'

        
        call ConstructEnterData     (Me%ObjEnterData, Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructWaveParameters - ModuleWaves - ERR010'
        
                
        call GetData(Me%WavePeriod%ON,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword    = 'WAVE_PERIOD',                                        &
                     Default    = .True.,                                               &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleWave',                                        &
                     STAT       = STAT_CALL)            

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructWaveParameters - ModuleWaves - ERR20'

        if (Me%WavePeriod%ON) then

            Me%WavePeriod%ID%Name = GetPropertyName(MeanWavePeriod_)
                     
            call ReadWaveParameters(WaveProperty = Me%WavePeriod,                       &
                                    BeginBlock   = "<begin_waveperiod>",                &
                                    EndBlock     = "<end_waveperiod>")    
        endif

        call GetData(Me%WaveHeight%ON,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword    = 'WAVE_HEIGHT',                                        &
                     Default    = .True.,                                               &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleWave',                                        &
                     STAT       = STAT_CALL)            

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructWaveParameters - ModuleWaves - ERR30'

        if (Me%WaveHeight%ON) then

            Me%WaveHeight%ID%Name = GetPropertyName(SignificantWaveHeight_)

            call ReadWaveParameters(WaveProperty = Me%WaveHeight,                       &
                                    BeginBlock   = "<begin_waveheight>",                &
                                    EndBlock     = "<end_waveheight>")    
        endif

        call GetData(Me%WaveDirection%ON,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword    = 'WAVE_DIRECTION',                                     &
                     Default    = .True.,                                               &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleWave',                                        &
                     STAT       = STAT_CALL)            

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructWaveParameters - ModuleWaves - ERR40'

        if (Me%WaveDirection%ON) then

            Me%WaveDirection%ID%Name = GetPropertyName(MeanWaveDirection_)

            call ReadWaveParameters(WaveProperty = Me%WaveDirection,                    &
                                    BeginBlock   = "<begin_wavedirection>",             &
                                    EndBlock     = "<end_wavedirection>")    
        endif

        call GetData(Me%RadiationStressX%ON,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword    = 'RADIATION_TENSION_X',                                &
                     Default    = .False.,                                              &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleWave',                                        &
                     STAT       = STAT_CALL)            

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructWaveParameters - ModuleWaves - ERR50'

        if (Me%RadiationStressX%ON) then

            Me%RadiationStressX%ID%Name = GetPropertyName(WaveStressX_)

            call ReadWaveParameters(WaveProperty = Me%RadiationStressX,                 &
                                    BeginBlock   = "<begin_radiationstress_x>",         &
                                    EndBlock     = "<end_radiationstress_x>")    

        endif

        call GetData(Me%RadiationStressY%ON,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword    = 'RADIATION_TENSION_Y',                                &
                     Default    = .False.,                                              &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleWave',                                        &
                     STAT       = STAT_CALL)            

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructWaveParameters - ModuleWaves - ERR60'

        if (Me%RadiationStressY%ON) then

            Me%RadiationStressY%ID%Name = GetPropertyName(WaveStressY_)

            call ReadWaveParameters(WaveProperty = Me%RadiationStressY,                 &
                                    BeginBlock   = "<begin_radiationstress_y>",         &
                                    EndBlock     = "<end_radiationstress_y>")    

        endif

        if (.not. Me%RadiationStressX%ON .EQV. Me%RadiationStressY%ON)                  &
            stop 'ConstructWaveParameters - ModuleWaves - ERR70'

        call GetData(Me%WaveGen_type,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     Keyword    = 'WAVEGEN_TYPE',                                       &
                     Default    = Old,                                                  &
                     SearchType = FromFile,                                             &
                     ClientModule ='ModuleWave',                                        &
                     STAT       = STAT_CALL)            

        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ConstructWaveParameters - ModuleWaves - ERR80'
        
        if (Me%WaveGen_type.eq.CEQUALW2) then

            call GetData(Me%WaveHeightParameter,                                        &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword    = 'WAVE_HEIGHT_PARAMETER',                          &
                         Default    = 1.,                                               &
                         SearchType = FromFile,                                         &
                         ClientModule ='ModuleWave',                                    &
                         STAT       = STAT_CALL)            

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructWaveParameters - ModuleWaves - ERR81'
        
            call GetData(Me%WavePeriodParameter,                                        &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword    = 'WAVE_PERIOD_PARAMETER',                          &
                         Default    = 1.,                                               &
                         SearchType = FromFile,                                         &
                         ClientModule ='ModuleWave',                                    &
                         STAT       = STAT_CALL)            

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructWaveParameters - ModuleWaves - ERR82'
                
            call GetData(Me%TotalDirections,                                            &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword    = 'WINDROSE_DIRECTIONS',                            &
                         Default    = 8,                                                &
                         SearchType = FromFile,                                         &
                         ClientModule ='ModuleWave',                                    &
                         STAT       = STAT_CALL)            

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructWaveParameters - ModuleWaves - ERR83'
            
            call GetData(Me%DistanceType,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         Keyword    = 'DISTANCE_TO_LAND_METHOD',                        &
                         Default    = Grid,                                             &
                         SearchType = FromFile,                                         &
                         ClientModule ='ModuleWave',                                    &
                         STAT       = STAT_CALL)            

            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructWaveParameters - ModuleWaves - ERR90'
        
        endif
        
        if (Me%WaveHeight%ON .and. Me%WavePeriod%ON) then

            Me%ParametersON = .true.

            allocate(Me%Abw(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%Abw(:,:)        =  FillValueReal

            allocate(Me%Ubw(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%Ubw(:,:)        =  FillValueReal

            allocate(Me%WaveLength(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
            Me%WaveLength(:,:) =  FillValueReal

        endif

        call ConstructGlobalOutput
        
        
        call Construct_Time_Serie
        

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructWaveParameters - ModuleWaves - ERR090'

        !Ungets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap,  Me%ExternalVar%WaterPoints2D,     &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructWaveParameters - ModuleWaves - ERR0100'


    end subroutine ConstructWaveParameters

    !--------------------------------------------------------------------------

    subroutine ReadWaveParameters(WaveProperty, BeginBlock, EndBlock)

        !Arguments-------------------------------------------------------------
        type (T_WaveProperty)               :: WaveProperty
        character(LEN = *)                  :: BeginBlock, EndBlock

        !Local-----------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL, iflag
        logical                             :: BlockFound

        !Begin-----------------------------------------------------------------

        !Searches for Wave Height
        call ExtractBlockFromBuffer (Me%ObjEnterData, ClientNumber,                             &
                                     BeginBlock, EndBlock,                                      &
                                     BlockFound, STAT = STAT_CALL)
        if (BlockFound) then

            !Allocates Variables
            allocate (WaveProperty%Field       (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            WaveProperty%Field      (:,:) = null_real

            call ConstructFillMatrix  (PropertyID           = WaveProperty%ID,                  &
                                       EnterDataID          = Me%ObjEnterData,                  &
                                       TimeID               = Me%ObjTime,                       &
                                       HorizontalGridID     = Me%ObjHorizontalGrid,             &
                                       ExtractType          = FromBlock,                        &
                                       PointsToFill2D       = Me%ExternalVar%WaterPoints2D,     &
                                       Matrix2D             = WaveProperty%Field,               &
                                       TypeZUV              = TypeZ_,                           &
                                       STAT                 = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadWaveParameters - ModuleWaves - ERR01'


            call GetIfMatrixRemainsConstant(FillMatrixID    = WaveProperty%ID%ObjFillMatrix,    &
                                            RemainsConstant = WaveProperty%Constant,            &
                                            STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadWaveParameters - ModuleWaves - ERR10'


            if(.not. WaveProperty%ID%SolutionFromFile)then
                call KillFillMatrix(WaveProperty%ID%ObjFillMatrix, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadWaveParameters - ModuleWaves - ERR20'
            end if

            call GetData(WaveProperty%OutputHDF,                                    &
                         Me%ObjEnterData, iflag,                                    &
                         Keyword        = 'OUTPUT_HDF',                             &
                         Default        = .false.,                                  &
                         SearchType     = FromBlock,                                &
                         ClientModule   = 'ModuleWave',                             &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ReadWaveParameters - ModuleWaves - ERR30'

            
            call GetData(WaveProperty%TimeSerieON,                                  &
                     Me%ObjEnterData, iflag,                                        &
                     keyword    = 'TIME_SERIE',                                     &
                     Default    = .false.,                                          &
                     SearchType = FromBlock,                                        &
                     ClientModule ='ModuleWave',                                    &
                     STAT       = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ConstructWaveParameters - ModuleWaves - ERR40'

            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadWaveParameters - ModuleWaves - ERR50'

        else

            write (*,*)'Block ',trim(BeginBlock),' ',trim(EndBlock),' not found'
            stop 'ReadWaveParameters - ModuleWaves - ERR60'

        endif

        call RewindBuffer (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadWaveParameters - ModuleWaves - ERR70'

        Me%LastCompute = Me%ActualTime

    end subroutine ReadWaveParameters

    !--------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: HDF5_CREATE
        
        !Begin----------------------------------------------------------------------

        if (Me%Output%HDF) then
        
            WorkILB = Me%WorkSize%ILB 
            WorkIUB = Me%WorkSize%IUB 
            WorkJLB = Me%WorkSize%JLB 
            WorkJUB = Me%WorkSize%JUB 

            !Gets a pointer to Bathymetry
            call GetGridData(Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT_CALL)     
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR01'
    
            !WaterPoints2D
            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D,    &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR10'
    
            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

            !Opens HDF5 File
            call ConstructHDF5(Me%ObjHDF5, trim(Me%OutputFile)//"5", HDF5_CREATE,       &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
               stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR20'

            !Write the Horizontal Grid
            call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
               stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR30'
    
            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB,       &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
               stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR40'

            !Writes the Grid
            call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",               &
                                  Array2D = Me%ExternalVar%Bathymetry,                  &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
               stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR50'

            call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",            &
                                  Array2D = Me%ExternalVar%WaterPoints2D,               &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
               stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR60'

            !Ungets the Bathymetry
            call UngetGridData (Me%ObjGridData, Me%ExternalVar%Bathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR70'

            !Ungets WaterPoints2D
            call UnGetHorizontalMap(Me%ObjHorizontalMap,  Me%ExternalVar%WaterPoints2D, &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR80'
        
            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                  &
               stop 'Open_HDF5_OutPut_File - ModuleWaves - ERR90'
        
        endif


    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine ConstructGlobalOutput 

        !Local-----------------------------------------------------------------
        logical                                     :: OutputON
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        

        nullify(Me%OutPut%OutTime)

        OutputON = OFF

        if ((Me%WavePeriod%OutputHDF).or.(Me%WaveHeight%OutputHDF).or.                  &
            (Me%WaveDirection%OutputHDF).or.(Me%RadiationStressX%OutputHDF).or.         &
            (Me%RadiationStressY%OutputHDF)) then
    
            OutputON = ON
            Me%Output%HDF = .true.

        endif

        if(OutputON)then

            call GetOutPutTime(Me%ObjEnterData,                                         &
                               CurrentTime     = Me%ActualTime,                         &
                               EndTime         = Me%EndTime,                            &
                               keyword         = 'OUTPUT_TIME',                         &
                               SearchType      = FromFile,                              &
                               OutPutsTime     = Me%OutPut%OutTime,                     &
                               OutPutsOn       = Me%OutPut%HDF,                         &
                               OutPutsNumber   = Me%OutPut%Number,                      &
                               STAT            = STAT_CALL)                           
                                                                                  
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ConstructGlobalOutput - ModuleWaves - ERR01'              
                                                                                
            if (Me%OutPut%HDF) then
                Me%OutPut%NextOutPut = 1
            else
                write(*,*)'Keyword OUTPUT_TIME must be defined if at least'
                write(*,*)'one property has HDF format outputs.'
                stop 'ConstructGlobalOutput - ModuleWaves - ERR10'
            endif 

        endif
 

    end subroutine ConstructGlobalOutput

    !------------------------------------------------------------------------

    subroutine Construct_Time_Serie

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        character(len=StringLength), dimension(:), pointer  :: PropertyList

        !Local-----------------------------------------------------------------
        real                                                :: CoordX, CoordY
        logical                                             :: CoordON, IgnoreOK
        integer                                             :: dn, Id, Jd, TimeSerieNumber
        integer                                             :: STAT_CALL
        integer                                             :: nProperties    = 0
        integer                                             :: PositionInList = 0
        integer                                             :: STATUS
        integer                                             :: iflag
        character(len=StringLength)                         :: TimeSerieLocationFile
        
        !Begin-----------------------------------------------------------------
        
        if ((Me%WaveHeight%TimeSerieOn).or.(Me%WavePeriod%TimeSerieOn).or.              &
            (Me%WaveDirection%TimeSerieOn).or.(Me%RadiationStressX%TimeSerieOn).or.     &
            (Me%RadiationStressY%TimeSerieOn)) then

            Me%OutPut%TimeSerie = .true.

            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D,    &
                                    STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Construct_TimeSerie - ModuleWaves - ERR10'
    
            if (Me%WaveHeight%TimeSerieOn) then
                nProperties = nProperties + 1
            endif
            if (Me%WavePeriod%TimeSerieOn) then
                nProperties = nProperties + 1
            endif
            if (Me%WaveDirection%TimeSerieOn) then
                nProperties = nProperties + 1
            endif
            if (Me%RadiationStressX%TimeSerieOn) then
                nProperties = nProperties + 1
            endif
            if (Me%RadiationStressY%TimeSerieOn) then
                nProperties = nProperties + 1
            endif

            !Allocates PropertyList
            allocate(PropertyList(nProperties), STAT = STATUS)

            if (STATUS /= SUCCESS_)                                                     &
                stop 'Construct_TimeSerie - ModuleWaves - ERR20'

            !Fills up PropertyList
            if (Me%WavePeriod%TimeSerieOn) then
                PositionInList = PositionInList + 1
                PropertyList (PositionInList) = 'Wave Period'
            endif
            if (Me%WaveHeight%TimeSerieOn) then
                PositionInList = PositionInList + 1
                PropertyList (PositionInList) = 'Wave Height'
            endif
            if (Me%WaveDirection%TimeSerieOn) then
                PositionInList = PositionInList + 1
                PropertyList (PositionInList) = 'Wave Direction'
            endif
            if (Me%RadiationStressX%TimeSerieOn) then
                PositionInList = PositionInList + 1
                PropertyList (PositionInList) = 'RadiationStressX'
            endif
            if (Me%RadiationStressY%TimeSerieOn) then
                PositionInList = PositionInList + 1
                PropertyList (PositionInList) = 'RadiationStressY'
            endif
            

            call GetData(TimeSerieLocationFile,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'TIME_SERIE_LOCATION',                          &
                         ClientModule = 'ModuleWaves',                                  &
                         Default      = Me%FileName,                                    &
                         STAT         = STATUS)
            if (STATUS .NE. SUCCESS_)                                                   &
                stop 'Construct_TimeSerie - ModuleWaves - ERR30'

            !Constructs TimeSerie
            call StartTimeSerie(Me%ObjTimeSerie, Me%ObjTime,                            &
                                TimeSerieLocationFile,                                  &
                                PropertyList, "srv",                                    &
                                WaterPoints2D = Me%ExternalVar%WaterPoints2D,           &
                                ModelName     = Me%ModelName,                           &
                                STAT = STATUS)

            if (STATUS /= SUCCESS_)                                                     &
                stop 'Construct_TimeSerie - ModuleWaves - ERR40'


            !Deallocates PropertyList
            deallocate(PropertyList, STAT = STATUS)

            if (STATUS /= SUCCESS_)                                                     &
                stop 'Construct_TimeSerie - ModuleWaves - ERR50'


            !Corrects if necessary the cell of the time serie based in the time serie coordinates
            call GetNumberOfTimeSeries(Me%ObjTimeSerie, TimeSerieNumber, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_TimeSerie - ModuleWaves - ERR60'

            do dn = 1, TimeSerieNumber

                call GetTimeSerieLocation(Me%ObjTimeSerie, dn,                              &  
                                          CoordX   = CoordX,                                &
                                          CoordY   = CoordY,                                & 
                                          CoordON  = CoordON,                               &
                                          STAT     = STAT_CALL)
                if (CoordON) then
                    call GetXYCellZ(Me%ObjHorizontalGrid, CoordX, CoordY, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_TimeSerie - ModuleWaves - ERR70'

                    if (Id < 0 .or. Jd < 0) then
                
                        call TryIgnoreTimeSerie(Me%ObjTimeSerie, dn, IgnoreOK, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'Construct_TimeSerie - ModuleWaves - ERR80'

                        if (IgnoreOK) then
                            cycle
                        else
                            stop 'Construct_TimeSerie - ModuleWaves - ERR90'
                        endif

                    endif

                    call CorrectsCellsTimeSerie(Me%ObjTimeSerie, dn, Id, Jd, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_TimeSerie - ModuleWaves - ERR100'
                endif


            enddo



            call UnGetHorizontalMap(Me%ObjHorizontalMap,  Me%ExternalVar%WaterPoints2D, &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Construct_TimeSerie - ModuleWaves - ERR110'
        
        endif

               
    end subroutine Construct_Time_Serie

    !------------------------------------------------------------------------

    subroutine ConstructFetch
    
        !External----------------------------------------------------------------
        integer                                         :: STAT_CALL
                
        !Local-------------------------------------------------------------------
    
        !Begin-------------------------------------------------------------------
        
        if (Me%Wavegen_type.eq.CEQUALW2) then
            
            write(*,*) 'Constructing Wind Fetch...'
            
            if (Me%TotalDirections.ne.8.and.Me%TotalDirections.ne.16) then
                write(*,*)'The model only runs with 8 or 16 Wind Rose directions '
                write(*,*)'Please select one of those on Waves data file'
                stop 'ConstructFetch - ModuleWaves - ERR22'

            endif

            if (Me%DistanceType.ne.Grid.and.Me%DistanceType.ne.Polygon) then
                write(*,*)'The model only has two methods to compute distance do land '
                write(*,*)'Please select 0 (zero) for grid method or 1 (one) for graphical'
                write(*,*)'method on DISTANCE_TO_LAND_METHOD keyword on Waves data file'
                stop 'ConstructFetch - ModuleWaves - ERR23'
            
            endif
            
            !Gets WaterPoints2D
            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D,    &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR01'
            
            !Gets information about grid: regular or unregular
            call GetCheckDistortion (Me%ObjHorizontalGrid, Me%ExternalVar%DistortionOn, &
                                        STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructFetch - ModuleWaves - ERR01a'

            !Gets Angle
            call GetGridAngle (Me%ObjHorizontalGrid, Angle = Me%ExternalVar%GridAngle,  &
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR010'
            
            !Gets XX_IE, YY_IE        
            call GetHorizontalGrid (Me%ObjHorizontalGrid, XX_IE = Me%ExternalVar%XX_IE, &
                                     YY_IE = Me%ExternalVar%YY_IE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR020'
    
            
            call GetHorizontalGrid(Me%ObjHorizontalGrid, DUX = Me%ExternalVar%DUX,      &        
                                   DVY = Me%ExternalVar%DVY, STAT = STAT_CALL)
                                   
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR021'

            if (Me%DistanceType.eq.Grid) then
            
                allocate(Me%XX(Me%Size%JLB:Me%Size%JUB))
                Me%XX(:) =  FillValueReal
                allocate(Me%YY(Me%Size%ILB:Me%Size%IUB))
                Me%YY(:) =  FillValueReal

            endif
            
            !allocate Distance, Fetch and Angles (angles with xx for all directions)
            allocate(Me%Fetch(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,1:8))
            Me%Fetch(:,:,:) =  FillValueReal
            
            allocate(Me%Distance(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,       &
                     1:Me%TotalDirections))
            Me%Distance(:,:,:) =  FillValueReal

            allocate(Me%AngleList(1:Me%TotalDirections))
            Me%AngleList(:) = FillValueReal
                
            if (Me%DistanceType.eq.Polygon) then
                               
                allocate(Me%Point)
                allocate(Me%MinDistanceFromPoint(1:Me%TotalDirections))
                Me%MinDistanceFromPoint(:) = FillValueReal
            
            endif

            call ComputeAngleList

            !Compute Distance to land and Fetch
            if(Me%DistanceType.eq.Polygon) then

                call ConstructLandArea
                call ComputeFetchDistancePolygon
            
            elseif(Me%DistanceType.eq.Grid) then
                
                call ComputeFetchDistanceGrid
            
            endif
            
            call ComputeFetch
            
            !deallocate lists not used in modify
            deallocate(Me%Distance)
            deallocate(Me%AngleList)

            if(Me%DistanceType.eq.Polygon) then
                deallocate(Me%Point)
                deallocate(Me%MinDistanceFromPoint)
            endif

            if (Me%DistanceType.eq.Grid) then
                deallocate(Me%XX)
                deallocate(Me%YY)
            endif

            !Ungets WaterPoints2D
            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D,  &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR030'
            
            !Unget XX_IE
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%XX_IE,        &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR040'

            !Unget YY_IE
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%YY_IE,        &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR050'

            !Unget DUX
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DUX,          &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR060'

            !Unget DVY
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%ExternalVar%DVY,          &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFetch - ModuleWaves - ERR070'
        
            write(*,*) 'Finished Constructing Wind Fetch...'
            
        endif
        
    end subroutine ConstructFetch
       
    !--------------------------------------------------------------------------

    subroutine ComputeAngleList

        !Define angles with xx positive axe, according to 8 or 16 wind rose directions                
        if (Me%TotalDirections.eq.8) then

            Me%W    = 1
            Me%SW   = 2
            Me%S    = 3
            Me%SE   = 4
            Me%E    = 5
            Me%NE   = 6
            Me%N    = 7
            Me%NW   = 8

            Me%AngleList(Me%W ) = 180.
            Me%AngleList(Me%SW) = 225.
            Me%AngleList(Me%S ) = 270.
            Me%AngleList(Me%SE) = 315.
            Me%AngleList(Me%E ) = 0.
            Me%AngleList(Me%NE) = 45.
            Me%AngleList(Me%N ) = 90.
            Me%AngleList(Me%NW) = 135.

        else if (Me%TotalDirections.eq.16) then

            Me%N    = 13
            Me%NNE  = 12
            Me%NE   = 11
            Me%ENE  = 10
            Me%E    = 9
            Me%ESE  = 8
            Me%SE   = 7
            Me%SSE  = 6
            Me%S    = 5
            Me%SSW  = 4
            Me%SW   = 3
            Me%WSW  = 2
            Me%W    = 1
            Me%WNW  = 16
            Me%NW   = 15
            Me%NNW  = 14
            
            Me%AngleList(Me%N  ) = 90. 
            Me%AngleList(Me%NNE) = 67.5
            Me%AngleList(Me%NE ) = 45.
            Me%AngleList(Me%ENE) = 22.5
            Me%AngleList(Me%E  ) = 0. 
            Me%AngleList(Me%ESE) = 337.5
            Me%AngleList(Me%SE ) = 315. 
            Me%AngleList(Me%SSE) = 292.5
            Me%AngleList(Me%S  ) = 270.
            Me%AngleList(Me%SSW) = 247.5
            Me%AngleList(Me%SW ) = 225.
            Me%AngleList(Me%WSW) = 202.5
            Me%AngleList(Me%W  ) = 180. 
            Me%AngleList(Me%WNW) = 157.5
            Me%AngleList(Me%NW ) = 135.
            Me%AngleList(Me%NNW) = 112.5
            
        endif

    
    end subroutine ComputeAngleList


    !--------------------------------------------------------------------------
    
    subroutine ComputeFetchDistanceGrid
    
        !Local-----------------------------------------------------------------
        integer                         :: i, j, x, y, runXX, runYY 
        integer                         :: ILB, IUB, JLB, JUB, Direction
        real                            :: PosX, PosY 
        real                            :: IterationDistance                   = 0.       
        real                            :: DistanceToLand                      = 0. 
        real(8)                         :: Angle
        real                            :: VerticalDistance, HorizontalDistance
                                                
        !Begin-----------------------------------------------------------------
        

        if (abs(Me%ExternalVar%GridAngle).gt.90.) then
            stop 'ComputeFetchDistanceGrid - ModuleWaves - ERR01'
        endif
 
        !Compute Distance to land
        !XX_IE and YY_IE are used instead of XX and YY to make the model compatible
        !with non metric coordinates

        !only regular grids
        if (Me%ExternalVar%DistortionOn) then

            write(*,*)'The WaveGen_Type selected in Wave data file does not support'
            write(*,*)'unregular grids. Choose WaveGen_Type:0 instead.'
            stop 'ComputeFetchDistanceGrid - ModuleWaves - ERR10'
        
        else

            IUB = Me%WorkSize%IUB
            ILB = Me%WorkSize%ILB
            JUB = Me%WorkSize%JUB
            JLB = Me%WorkSize%JLB
        
            do j=JLB, JUB
            do i=ILB, IUB
            
            
                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                    do Direction = Me%W, Me%TotalDirections

                        x = j
                        y = i
                        Angle = Me%AngleList(Direction) * Pi / 180.
                        DistanceToLand = 0.
                                                
                        !center of cell
                        PosX = Me%ExternalVar%XX_IE(i,j) + (Me%ExternalVar%DUX(i,j))/2.
                        PosY = Me%ExternalVar%YY_IE(i,j) + (Me%ExternalVar%DVY(i,j))/2.
                                                                           
                        do while (Me%ExternalVar%WaterPoints2D(y, x) == WaterPoint)
                        
                            IterationDistance = (min(Me%ExternalVar%DUX(y,x),           &
                                                    Me%ExternalVar%DVY(y,x)))/2.
                        
                            DistanceToLand = DistanceToLand + IterationDistance

                            VerticalDistance = sin(Angle)*IterationDistance
                            HorizontalDistance = cos(Angle)*IterationDistance
                            
                            !new position
                            PosY = PosY + VerticalDistance
                            PosX = PosX + HorizontalDistance

                            if ((PosX .lt. Me%ExternalVar%XX_IE(y,JLB  )) .or.  &
                                (PosY .lt. Me%ExternalVar%YY_IE(ILB,  x)) .or.  &
                                (PosX .gt. Me%ExternalVar%XX_IE(y,JUB+1)) .or.  &
                                (PosY .gt. Me%ExternalVar%YY_IE(IUB+1,x))) then
                                exit
                            endif
                            
                            !To drop one dimension (XX_IE and YY_IE are 2D and XX and 
                            !YY are 1D). XX and YY are used in the locatecell function. 
                            !In regular grids XX is independent of y and YY independent of x
                            do runXX = JLB, JUB+1
                                Me%XX(runXX) = Me%ExternalVar%XX_IE(y,runXX)
                            enddo

                            do runYY = ILB, IUB+1
                                Me%YY(runYY) = Me%ExternalVar%YY_IE(runYY,x)
                            enddo

                            call LocateCell(XX = Me%XX, YY = Me%YY, XPos = PosX,            &
                                            YPos = PosY, ILB = ILB, IUB = IUB+1,            &
                                            JLB = JLB, JUB = JUB+1, ILower = y, JLower = x)

                            if(y < 0 .or. x < 0) exit 

                        enddo
                
                        Me%Distance(i,j, Direction) = DistanceToLand

                    enddo
                endif
        
            enddo
            enddo
        
        endif

    end subroutine ComputeFetchDistanceGrid

    !--------------------------------------------------------------------------
    
    subroutine ComputeFetchDistancePolygon
        
        !Local-----------------------------------------------------------------
        integer                                 :: i, j, k, ILB, IUB, JLB, JUB
        
        !Begin-----------------------------------------------------------------
        
        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB

        !Compute Distance To Land (To Polygon Border)
        do j=JLB, JUB
        do i=ILB, IUB
            
            !Coordinates at the cell center
            Me%Point%X = Me%ExternalVar%XX_IE(i,j) + (Me%ExternalVar%DUX(i,j))/2.                 
            Me%Point%Y = Me%ExternalVar%YY_IE(i,j) + (Me%ExternalVar%DVY(i,j))/2.
                         
            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                Call PointDistanceToPolygon (Point = Me%Point, Polygons =        &
                                             Me%LandArea,                        &
                                             AngleList = Me%AngleList,           &
                                             MinDistanceFromPoint =              &
                                             Me%MinDistanceFromPoint)              
                
                !Update the distance to land list with the data computed for that point
                !in Module Drawing               
                do k = 1, Me%TotalDirections
                    Me%Distance (i,j,k) = Me%MinDistanceFromPoint(k)
                enddo
            
            endif

        enddo
        enddo

    end subroutine ComputeFetchDistancePolygon
    
    !--------------------------------------------------------------------------
    
    subroutine ComputeFetch
       
        !Local-----------------------------------------------------------------
        integer                                 :: i, j, ILB, IUB, JLB, JUB
        real                                    :: Weight22, Weight45, SumWeight

        !Begin-----------------------------------------------------------------

        !Fetch is computed in 8 directions. The difference is in the way of 
        !calculation, if 8 or 16 wind rose directions are selected.
        !With 8 wind rose directions, for each fetch direction 3 distances and 3 angles are used 
        !(the wind direction, one 45 degrees to the right and one 45 degrees to the left.
        !With 16 wind rose directions, for each fetch direction 5 distances and 5 angles are used 
        !(the wind direction, 22.5 and 45 degrees to the right and the same to the left.

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        
        !The weight is given according to cosine of the angles (angle with wind direction). 
        !The weight along the wind direction (0 deg.) is maximum (1) reducing when angle icreases.
        Weight22 = cos (22.5 * Pi/180.)
        Weight45 = cos(45. * Pi/180.)

        if (Me%TotalDirections.eq.8) then
        
            !The sum of weights with 8 wind rose directions
            SumWeight = (2. * Weight45 + 1.)
        
            do j=JLB, JUB
            do i=ILB, IUB
            
                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                    Me%Fetch(i,j, Me%N) = (Weight45 * Me%Distance(i,j, Me%NW)           &
                                            + Me%Distance(i,j, Me%N)                    &
                                            + Weight45 * Me%Distance(i,j, Me%NE)) / SumWeight
                    
                    Me%Fetch(i,j, Me%NE) = (Weight45 * Me%Distance(i,j, Me%N)           &
                                            + Me%Distance(i,j, Me%NE)                   &
                                            + Weight45 * Me%Distance(i,j, Me%E)) / SumWeight
                   
                    Me%Fetch(i,j, Me%E) = (Weight45 * Me%Distance(i,j, Me%NE)           &
                                            + Me%Distance(i,j, Me%E)                    &
                                            + Weight45 * Me%Distance(i,j, Me%SE)) / SumWeight
                    
                    Me%Fetch(i,j, Me%SE) = (Weight45 * Me%Distance(i,j, Me%E)           &
                                            + Me%Distance(i,j, Me%SE)                   &
                                            + Weight45 * Me%Distance(i,j, Me%S)) / SumWeight
                   
                    Me%Fetch(i,j, Me%S) = (Weight45 * Me%Distance(i,j, Me%SE)           &
                                            + Me%Distance(i,j, Me%S)                    &
                                            + Weight45 * Me%Distance(i,j, Me%SW)) / SumWeight
                    
                    Me%Fetch(i,j, Me%SW) = (Weight45 * Me%Distance(i,j, Me%S)           &
                                            + Me%Distance(i,j, Me%SW)                   &
                                            + Weight45 * Me%Distance(i,j, Me%W)) / SumWeight
                    
                    Me%Fetch(i,j, Me%W) = (Weight45 * Me%Distance(i,j, Me%SW)           &
                                            + Me%Distance(i,j, Me%W)                    &
                                            + Weight45 * Me%Distance(i,j, Me%NW)) / SumWeight
                    
                    Me%Fetch(i,j, Me%NW) = (Weight45 * Me%Distance(i,j, Me%W)           &
                                            + Me%Distance(i,j, Me%NW)                   &
                                            + Weight45 * Me%Distance(i,j, Me%N)) / SumWeight
                endif

            enddo
            enddo
        
        else if (Me%TotalDirections.eq.16) then

            !The sum of weights with 16 wind rose directions
            SumWeight = (2.*Weight22 + 2.*Weight45 + 1.)
        
            do j=JLB, JUB
            do i=ILB, IUB
            
                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                    !The 8 fetch directions (N-NW) now are classified from 1 to 16

                    !N
                    Me%Fetch(i,j, 7) = (Weight45 * Me%Distance(i,j, Me%NW )             &
                                      + Weight22 * Me%Distance(i,j, Me%NNW)             &
                                      +            Me%Distance(i,j, Me%N  )             &
                                      + Weight22 * Me%Distance(i,j, Me%NNE)             &
                                      + Weight45 * Me%Distance(i,j, Me%NE )) / SumWeight

                    !NE
                    Me%Fetch(i,j, 6) = (Weight45 * Me%Distance(i,j, Me%N)         &
                                                + Weight22 * Me%Distance(i,j, Me%NNE)   &
                                                + Me%Distance(i,j, Me%NE)               &
                                                + Weight22 * Me%Distance(i,j, Me%ENE)   &
                                                + Weight45 * Me%Distance(i,j, Me%E)) / SumWeight
                    
                    !E
                    Me%Fetch(i,j, 5) = (Weight45 * Me%Distance(i,j, Me%NE) &
                                                + Weight22 * Me%Distance(i,j, Me%ENE)   &
                                                + Me%Distance(i,j, Me%E)                &
                                                + Weight22 * Me%Distance(i,j, Me%ESE)   &
                                                + Weight45 * Me%Distance(i,j, Me%SE)) / SumWeight
                    !SE
                    Me%Fetch(i,j, 4) = (Weight45 * Me%Distance(i,j, Me%E)         &
                                                + Weight22 * Me%Distance(i,j, Me%ESE)   &
                                                + Me%Distance(i,j, Me%SE)               &
                                                + Weight22 * Me%Distance(i,j, Me%SSE)   &
                                                + Weight45 * Me%Distance(i,j, Me%S)) / SumWeight
                    !S
                    Me%Fetch(i,j, 3) = (Weight45 * Me%Distance(i,j, Me%SE)         &
                                                + Weight22 * Me%Distance(i,j, Me%SSE)   &
                                                + Me%Distance(i,j, Me%S)                &
                                                + Weight22 * Me%Distance(i,j, Me%SSW)   &
                                                + Weight45 * Me%Distance(i,j, Me%SW)) / SumWeight
                    !SW
                    Me%Fetch(i,j, 2) = (Weight45 * Me%Distance(i,j, Me%S)         &
                                                + Weight22 * Me%Distance(i,j, Me%SSW)   &
                                                + Me%Distance(i,j, Me%SW)               &
                                                + Weight22 * Me%Distance(i,j, Me%WSW)   &
                                                + Weight45 * Me%Distance(i,j, Me%W)) / SumWeight
                    !W
                    Me%Fetch(i,j, 1) = (Weight45 * Me%Distance(i,j, Me%SW)         &
                                                + Weight22 * Me%Distance(i,j, Me%WSW)   &
                                                + Me%Distance(i,j, Me%W)                &
                                                + Weight22 * Me%Distance(i,j, Me%WNW)   &
                                                + Weight45 * Me%Distance(i,j, Me%NW)) / SumWeight
                    !NW
                    Me%Fetch(i,j, 8) = (Weight45 * Me%Distance(i,j, Me%W)         &
                                                + Weight22 * Me%Distance(i,j, Me%WNW)   &
                                                + Me%Distance(i,j, Me%NW)               &
                                                + Weight22 * Me%Distance(i,j, Me%NNW)   &
                                                + Weight45 * Me%Distance(i,j, Me%N)) / SumWeight

                endif

            enddo
            enddo

        endif

    end subroutine ComputeFetch

    !--------------------------------------------------------------------------
    
    subroutine ConstructLandArea

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        character(len=line_length)          :: LandAreaFilesPath
        character(len=line_length)          :: FullBufferLine
        integer                             :: ClientNumber
        integer                             :: CurrentLineNumber, StartLine, EndLine
        logical                             :: BlockFound
        integer                             :: CoordType, i
        integer                             :: GEOG, SIMPLE_GEOG
        real                                :: Xorig, Yorig
        real                                :: radians, EarthRadius
        real                                :: Rad_Lat, CosenLat, Radius
        type(T_Polygon), pointer            :: Polygon

        !Begin-----------------------------------------------------------------
        
        call ConstructEnterData     (Me%ObjEnterData, Me%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR00'
        
        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber,                                       &
                                    '<BeginLandAreaFiles>',                             &
                                    '<EndLandAreaFiles>',                               &
                                    BlockFound,                                         &
                                    STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR10'

        if (.not.BlockFound) stop 'ConstructLandArea - ModuleWaves - ERR11'
          
        call GetBlockSize(Me%ObjEnterData,                                              &
                          ClientNumber,                                                 &
                          StartLine,                                                    &
                          EndLine,                                                      &
                          FromBlock_,                                                   &
                          STAT = STAT_CALL)
       
        do CurrentLineNumber = StartLine + 1 , EndLine - 1
            
            call GetFullBufferLine(Me%ObjEnterData,                                     &
                                   CurrentLineNumber,                                   &
                                   FullBufferLine,                                      &
                                   STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR20'

            !Get Polygon file path
            call GetData(LandAreaFilesPath,                                             &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType  = FromBlock_,                                      &
                         Buffer_Line = CurrentLineNumber,                               &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR30'

            !Construct PolygonArea collection
            call New(Me%LandArea,  trim(LandAreaFilesPath))
                
        enddo

        call GetCoordTypeList (GEOG = GEOG, SIMPLE_GEOG = SIMPLE_GEOG)

        call GetGridCoordType(Me%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR35'

        call GetGridOrigin(Me%ObjHorizontalGrid, Xorig = Xorig, Yorig = Yorig, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR36'

        if    (CoordType == GEOG)then
            
            write(*,*)"Cannot compute waves fetch distances in GEOG coordinates."
            write(*,*)"Conversion of polygon to grid metric coordinates needs to be programmed!"
            write(*,*)"Meanwhile you can choose SIMPLE_GEOG coordinates, if possible"
            stop 'ConstructLandArea - ModuleWaves - ERR37'

        elseif(CoordType == SIMPLE_GEOG) then

            radians      = 180.0 / PI
            EarthRadius  = 6378000.

            Polygon => Me%LandArea

            do while(associated(Polygon))

                do i = 1, Polygon%Count

                    Rad_Lat = (Polygon%VerticesF(i)%Y)/ radians 
                    CosenLat= cos(Rad_Lat)
                    Radius  = CosenLat * EarthRadius
                    
                    Polygon%VerticesF(i)%X = (Polygon%VerticesF(i)%X - Xorig) * Radius / radians
                    Polygon%VerticesF(i)%Y = (Polygon%VerticesF(i)%Y - Yorig) * EarthRadius / radians

                enddo

                call SetLimits(Polygon)

                Polygon => Polygon%Next
            
            enddo

        end if

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR40'

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR50'

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructLandArea - ModuleWaves - ERR060'

    end subroutine ConstructLandArea

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetWaves (WavesID, WavePeriod, WaveHeight, Abw, Ubw, WaveLength, &
                         WaveDirection, LastCompute, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: WavesID
        real, dimension(:,:),  pointer, optional        :: WaveHeight, WavePeriod
        real, dimension(:,:),  pointer, optional        :: Abw, Ubw, WaveLength, WaveDirection
        type(T_Time)        ,           optional        :: LastCompute
        integer, intent(OUT),           optional        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%WaveHeight%ON .and. Me%WavePeriod%ON) then


                if (present(WavePeriod)) then
                    call Read_Lock(mWaves_, Me%InstanceID)
                    WavePeriod => Me%WavePeriod%Field
                endif

                if (present(WaveHeight)) then
                    call Read_Lock(mWaves_, Me%InstanceID)
                    WaveHeight => Me%WaveHeight%Field
                endif

                if (present(Abw)) then
                    call Read_Lock(mWaves_, Me%InstanceID)
                    Abw => Me%Abw
                endif

                if (present(Ubw)) then
                    call Read_Lock(mWaves_, Me%InstanceID)
                    Ubw => Me%Ubw
                endif

                if (present(WaveLength)) then
                    call Read_Lock(mWaves_, Me%InstanceID)
                    WaveLength => Me%WaveLength
                endif

                if (present(WaveDirection)) then
                    call Read_Lock(mWaves_, Me%InstanceID)
                    WaveDirection => Me%WaveDirection%Field
                endif


                if (present(LastCompute)) then 
                    LastCompute = Me%LastCompute
                endif

                STAT_ = SUCCESS_

            endif

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine GetWaves

    !--------------------------------------------------------------------------

    subroutine GetWavesOptions (WavesID, WaveParametersON, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: WavesID
        logical, intent(OUT),           optional        :: WaveParametersON
        integer, intent(OUT),           optional        :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(WaveParametersON)) then

                WaveParametersON = Me%ParametersON

            endif


            STAT_ = SUCCESS_


        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetWavesOptions

    !--------------------------------------------------------------------------

    subroutine GetWavesStress (WavesID, RadiationStressX, RadiationStressY,             &
                               LastCompute, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: WavesID
        real, dimension(:,:),  pointer                  :: RadiationStressX, RadiationStressY
        type(T_Time)        , optional                  :: LastCompute
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                    &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            if (Me%RadiationStressX%ON .and. Me%RadiationStressY%ON) then

                call Read_Lock(mWaves_, Me%InstanceID)
                RadiationStressX => Me%RadiationStressX%Field

                call Read_Lock(mWaves_, Me%InstanceID)
                RadiationStressY => Me%RadiationStressY%Field

                if (present(LastCompute)) LastCompute = Me%LastCompute


                STAT_ = SUCCESS_
            endif

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetWavesStress
    
    !--------------------------------------------------------------------------

    subroutine UnGetWaves(WavesID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: WavesID
        real, dimension(:, :), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mWaves_, Me%InstanceID,  "UnGetWaves")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetWaves

    !--------------------------------------------------------------------------

    subroutine SetWavesWind(WavesID, WindX, WindY, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: WavesID
        real, pointer, dimension(:,:)               :: WindX, WindY
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            if (associated(WindX))         Me%ExternalVar%WindVelocityX => WindX
            if (associated(WindY))         Me%ExternalVar%WindVelocityY => WindY
            
            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            
    end subroutine SetWavesWind

    !--------------------------------------------------------------------------

    subroutine SetGeneric4DValues(WavesID, CurrentValue4D, STAT)

        !Arguments---------------------------------------------------------------
        integer                                     :: WavesID
        real                                        :: CurrentValue4D
        integer,            optional, intent(OUT)   :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_    

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_) 
        
        if (ready_ .EQ. IDLE_ERR_)then

            Me%ExternalVar%CurrentValue4D = CurrentValue4D
            
            STAT_ = SUCCESS_  

        else
            STAT_ = ready_
        end if


        if (present(STAT))STAT = STAT_
            

    end subroutine SetGeneric4DValues

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
   
    !--------------------------------------------------------------------------
   
    subroutine ModifyWaves(WavesID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: WavesID
        integer, intent(OUT), optional              :: STAT
      
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%ActualTime, STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR01'

            !WaterPoints2D
            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D,            &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR10'
          
            !Modifies Wave Height
            if (Me%WaveHeight%ON) then
                if (.not. Me%WaveHeight%Constant) then
                    if (Me%WaveHeight%ID%SolutionFromFile) then
                        call ModifyFillMatrix (FillMatrixID    = Me%WaveHeight%ID%ObjFillMatrix,&
                                               Matrix2D         = Me%WaveHeight%Field,          &
                                               PointsToFill2D   = Me%ExternalVar%WaterPoints2D, &
                                               Generic_4D_Value = Me%ExternalVar%CurrentValue4D,&
                                               STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR20'
                    else
                     
                        if (Me%Wavegen_type.eq.Old) then   
                        
                            call ComputeWaveHeightSimple

                        else if (me%Wavegen_type.eq.CEQUALW2) then

                            call ComputeWaveHeightFetch

                        else
                        
                            stop 'ModifyWaves - ModuleWaves - ERR30'
                    
                        endif

                   endif
                endif
            endif


            !Modifies Wave Period
            if (Me%WavePeriod%ON) then
                if (.not. Me%WavePeriod%Constant) then
                    if (Me%WavePeriod%ID%SolutionFromFile) then
                        call ModifyFillMatrix (FillMatrixID     = Me%WavePeriod%ID%ObjFillMatrix,&
                                               Matrix2D         = Me%WavePeriod%Field,           &
                                               PointsToFill2D   = Me%ExternalVar%WaterPoints2D,  &
                                               Generic_4D_Value = Me%ExternalVar%CurrentValue4D, &
                                               STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR40'
                    else
                        
                        
                        if (Me%Wavegen_type.eq.Old) then   
                                      
                            call ComputeWavePeriodSimple

                        else if (me%Wavegen_type.eq.CEQUALW2) then

                            call ComputeWavePeriodFetch

                        else
                   
                            stop 'ModifyWaves - ModuleWaves - ERR50'
                    
                        endif
                       
                   endif
                endif
            endif


            !Modifies Wave Direction
            if (Me%WaveDirection%ON) then
                if (.not. Me%WaveDirection%Constant) then  
                    if (Me%WaveDirection%ID%SolutionFromFile) then
                        call ModifyFillMatrix (FillMatrixID     = Me%WaveDirection%ID%ObjFillMatrix,&
                                               Matrix2D         = Me%WaveDirection%Field,           &
                                               PointsToFill2D   = Me%ExternalVar%WaterPoints2D,     &
                                               Generic_4D_Value = Me%ExternalVar%CurrentValue4D,    &
                                               STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR60'
                    else
                        call ComputeWaveDirection
                    endif
                endif
            endif


            if (Me%ParametersON .and. (.not. Me%WaveHeight%Constant .or. .not. Me%WavePeriod%Constant)) then

                call ComputeWaveParameters

            endif


            !Modifies Radiation Tension
            if (Me%RadiationStressX%ON) then
                if (.not. Me%RadiationStressX%Constant) then
                    if (Me%RadiationStressX%ID%SolutionFromFile) then
                        call ModifyFillMatrix (FillMatrixID     = Me%RadiationStressX%ID%ObjFillMatrix, &
                                               Matrix2D         = Me%RadiationStressX%Field,            &
                                               PointsToFill2D   = Me%ExternalVar%WaterPoints2D,         &
                                               Generic_4D_Value = Me%ExternalVar%CurrentValue4D,        &
                                               STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR70'
                    else
                        !call ComputeRadiationStress
                    endif
                endif
            endif

            if (Me%RadiationStressY%ON) then
                if (.not. Me%RadiationStressY%Constant) then
                    if (Me%RadiationStressY%ID%SolutionFromFile) then
                        call ModifyFillMatrix (FillMatrixID     = Me%RadiationStressY%ID%ObjFillMatrix, &
                                               Matrix2D         = Me%RadiationStressY%Field,            &
                                               PointsToFill2D   = Me%ExternalVar%WaterPoints2D,         &
                                               Generic_4D_Value = Me%ExternalVar%CurrentValue4D,        &
                                               STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR80'
                    else
                        !call ComputeRadiationStress
                    endif
                endif
            endif

            
            Me%LastCompute = Me%ActualTime

            !Ungets WaterPoints2D
            call UnGetHorizontalMap(Me%ObjHorizontalMap,  Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyWaves - ModuleWaves - ERR90'
            
            
            if(Me%OutPut%HDF) then                                 
                call OutPut_Results_HDF
            endif
            
            if (Me%OutPut%TimeSerie) then
                call Output_TimeSeries
            endif

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyWaves

    !--------------------------------------------------------------------------
  
    subroutine ComputeWaveParameters(WavesID, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN ), optional              :: WavesID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        real                                        :: COEFA, COEFB, WAVN, OMEG, C0, Celerity
        integer                                     :: STAT_, ready_
        integer                                     :: i, j, ILB, IUB, JLB, JUB, STAT_CALL
        logical                                     :: PrivateSub

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        ready_ = FillValueInt

        if (present(WavesID)) then
            call Ready(WavesID, ready_)
            PrivateSub = .false.
        else
            PrivateSub = .true.
        endif

        if (ready_ .EQ. IDLE_ERR_ .or. PrivateSub) then

            IUB = Me%WorkSize%IUB
            ILB = Me%WorkSize%ILB
            JUB = Me%WorkSize%JUB
            JLB = Me%WorkSize%JLB


            call GetGeometryWaterColumn(Me%ObjGeometry,                                 &
                                        WaterColumn = Me%ExternalVar%WaterColumn,       &
                                        STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveParameters - ModuleWaves - ERR01.'

            if (.not. PrivateSub) then
                call GetWaterPoints2D(Me%ObjHorizontalMap,                              &
                                      WaterPoints2D = Me%ExternalVar%WaterPoints2D,     &
                                      STAT = STAT_CALL)
            endif

            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveParameters - ModuleWaves - ERR10.'



            do j=JLB, JUB
            do i=ILB, IUB
            
                if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then


cd2:                if (Me%WaveHeight%Field       (i,j) .lt. 0.1 .or.                   &
                        Me%ExternalVar%WaterColumn(i,j) .lt. 0.1 .or.                   &
                        Me%WavePeriod%Field       (i,j) < 1e-3) then 
                        
                        Me%Ubw(i, j)              = 0.001
                        Me%Abw(i, j)              = 0.001
                        Me%WaveLength(i, j)       = 0.

                    else cd2

                        COEFA = Gravity * Me%WavePeriod%Field(i,j) / (2.*PI)
                                                  
                        COEFB = 2. * PI / Me%WavePeriod%Field(i,j)

                        C0    = sqrt(Gravity*Me%ExternalVar%WaterColumn(i,j))

                        call Secant(Celerity, C0, Me%ExternalVar%WaterColumn(i,j), COEFA, COEFB)

                        Me%WaveLength(i, j)      = Celerity   * Me%WavePeriod%Field(i,j)
                        OMEG                     = 2.  * PI   / Me%WavePeriod%Field(i,j)
                        WAVN                     = 2.  * PI   / Me%WaveLength(i, j)

                        !To avoid sinh results larger then 1e100
                        !if (WAVN * Me%ExternalVar%WaterColumn(i,j) .gt. 230) then
                        !To avoid sinh results larger then 1e10
                        if (WAVN * Me%ExternalVar%WaterColumn(i,j) .gt. 23.7) then
                            Me%Ubw(i, j)             = 0
                        else
                            Me%Ubw(i, j)             = 0.5 * OMEG * Me%WaveHeight%Field(i,j) / &
                                                  sinh(WAVN * Me%ExternalVar%WaterColumn(i,j))
                        endif

                        !To avoid overflow errors
                        if (Me%Ubw(i, j) < 1e-6)  Me%Ubw(i, j) = 0.
                        
                        Me%Abw(i, j) = Me%Ubw(i, j) / OMEG

                    endif cd2

                endif

            enddo
            enddo
            
            !Ungets geometry
            call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%WaterColumn, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveParameters; ModuleWaves. ERR20'

            if (.not. PrivateSub) then
                !Ungets WaterPoints2D
                call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveParameters - ModuleWaves - ERR30'
            endif


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_


    end subroutine ComputeWaveParameters

    !--------------------------------------------------------------------------
    
    subroutine ComputeWaveHeightSimple

        !Local-----------------------------------------------------------------
        real                                    :: Wind
        integer                                 :: i, j, ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------
       
        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB

        do j=JLB, JUB
        do i=ILB, IUB
            
            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then


                Wind = sqrt(Me%ExternalVar%WindVelocityX(i,j)**2. + &
                            Me%ExternalVar%WindVelocityY(i,j)**2.)  

                Me%WaveHeight%Field(i, j) = 0.243 * Wind * Wind / Gravity

            end if
        
        end do
        end do

            
    end subroutine ComputeWaveHeightSimple
    
    !--------------------------------------------------------------------------

    subroutine ComputeWaveHeightFetch

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real                                    :: Wind, LeftDirection, RightDirection 
        real                                    :: COEF1, COEF2, U2, WindX, WindY 
        integer                                 :: i, j, ILB, IUB, JLB, JUB, WindDirection
        real                                    :: WindAngle
        
        !Begin-----------------------------------------------------------------
        
        
        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
    
        call GetGeometryWaterColumn(Me%ObjGeometry,                                     &
                                    WaterColumn = Me%ExternalVar%WaterColumn,           &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveHeightFetch - ModuleWaves. ERR01.'

        do j=JLB, JUB
        do i=ILB, IUB

            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
        
               
                WindX = Me%ExternalVar%WindVelocityX(i,j)
                WindY = Me%ExternalVar%WindVelocityY(i,j)
                !Angle in degrees in triginometric referencial
                WindAngle = atan2(WindY,WindX) * 180. / Pi
                
                !only positive angles for next operations
                if (WindAngle .lt. 0.0) then
                    WindAngle = WindAngle + 360.
                endif
                
                !Wind directions are fused into 8 major directions (N,NE,E,SE,S,SW,W,NW)
                !Numbers are used instead of Me%[direction] because fetch as only 8 positions
                !no matter if 8 or 16 wind rose directions are used
                if ((WindAngle .gt. 337.5) .and. (WindAngle .le. 360.) .Or. &
                    (WindAngle .ge. 0.   ) .and. (WindAngle .le. 22.5)) then
                    
                    !if trigonometric angle is near zero then wind comes from West
                    WindDirection = Me%W

                else 
    
                    WindDirection = Me%W + 1
                    LeftDirection = 22.5
                    RightDirection = 67.5
            
                    do while (WindDirection.le.8)

                        if ((WindAngle .gt. LeftDirection).And.&
                            (WindAngle .le. RightDirection)) then
                            exit 
                        else
                            LeftDirection  = RightDirection
                            RightDirection = RightDirection + 45.
                            WindDirection  = WindDirection + 1
                        endif

                    enddo
                endif

                if (WindDirection.gt.8) then
                    stop 'ComputeWaveHeightFetch - ModuleWaves - ERR10'
                endif
        
                !compute
                Wind = sqrt(WindX**2. + WindY**2.)  

                U2 = Wind * Wind         

                COEF1 = 0.53 * (Gravity * Me%ExternalVar%WaterColumn(i, j) / U2) ** 0.75   

                COEF2 = 0.0125 * (Gravity * Me%Fetch(i,j, WindDirection) / U2)**0.42


                !WaveHeight
                Me%WaveHeight%Field(i, j) = Me%WaveHeightParameter * 0.283 * U2 /   &
                                            Gravity *                               &
                                            Tanh(COEF1) * Tanh(COEF2 / Tanh(COEF1))    
     
            endif

        enddo
        enddo
    
        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWaveHeightFetch - ModuleWaves. ERR20'

   
    
                    
    end subroutine ComputeWaveHeightFetch
    
    !--------------------------------------------------------------------------

    subroutine ComputeWavePeriodSimple

        !Local-----------------------------------------------------------------
        real                                    :: Wind
        integer                                 :: i, j, ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------


        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB

        do j=JLB, JUB
        do i=ILB, IUB
            
            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                Wind = sqrt(Me%ExternalVar%WindVelocityX(i,j)**2. + &
                            Me%ExternalVar%WindVelocityY(i,j)**2.)

                Me%WavePeriod%Field(i, j) = Wind * 8.13 / Gravity

            endif
        
        end do
        end do

            
    end subroutine ComputeWavePeriodSimple

    !--------------------------------------------------------------------------

    subroutine ComputeWavePeriodFetch

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real                                    :: Wind, LeftDirection, RightDirection 
        real                                    :: COEF3, COEF4, U2, WindX, WindY
        integer                                 :: i, j, ILB, IUB, JLB, JUB, WindDirection
        real                                    :: WindAngle

       !Begin-----------------------------------------------------------------
        

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
    
        call GetGeometryWaterColumn(Me%ObjGeometry,                                     &
                                    WaterColumn = Me%ExternalVar%WaterColumn,           &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWavePeriodFetch - ModuleWaves. ERR01.'

        do j=JLB, JUB
        do i=ILB, IUB
            
            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then
        
                WindX = Me%ExternalVar%WindVelocityX(i,j)
                WindY = Me%ExternalVar%WindVelocityY(i,j)
                !Angle in degrees in triginometric referencial
                WindAngle = atan2(WindY,WindX) * 180. / Pi

                !only positive angles for next operations
                if (WindAngle .lt. 0.0) then
                    WindAngle = WindAngle + 360.
                endif

                !Wind directions are fused into 8 major directions (N,NE,E,SE,S,SW,W,NW)
                !Numbers are used instead of Me%[direction] because fetch as only 8 positions
                !no matter if 8 or 16 wind rose directions are used
                if ((WindAngle .gt. 337.5) .and. (WindAngle .le. 360.) .Or. &
                    (WindAngle .ge. 0.   ) .and. (WindAngle .le. 22.5)) then
                    
                    !if trigonometric angle is near zero then wind comes from West                    
                    WindDirection = Me%W

                else 
    
                    WindDirection = Me%W + 1
                    LeftDirection = 22.5
                    RightDirection = 67.5
            
                    do while (WindDirection.le.8)

                        if ((WindAngle .gt. LeftDirection).And.&
                            (WindAngle .le. RightDirection)) then
                            exit 
                        else
                            LeftDirection  = RightDirection
                            RightDirection = RightDirection + 45.
                            WindDirection  = WindDirection + 1
                        endif

                    enddo
                endif
        
                if (WindDirection.gt.8) then
                    stop 'ComputeWavePeriodFetch - ModuleWaves - ERR10'
                endif

                Wind = sqrt(WindX**2. +  WindY**2.)  

                U2 = Wind * Wind              

                COEF3 = 0.833 * (Gravity * Me%ExternalVar%WaterColumn(i,j) / U2)**0.375

                COEF4 = 0.077 * (Gravity * Me%Fetch(i,j, WindDirection) / U2)**0.25

                !Wave Period
                Me%WavePeriod%Field(i, j) = Me%WavePeriodParameter * 2 * Pi * Wind / &
                                            Gravity * 1.2 *                          &
                                            Tanh(COEF3) * Tanh(COEF4 / Tanh(COEF3))

            endif

        enddo
        enddo
    
        call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ComputeWavePeriodFetch - ModuleWaves. ERR20'
    
           
    end subroutine ComputeWavePeriodFetch

    !--------------------------------------------------------------------------

    subroutine ComputeWaveDirection

        !Local-----------------------------------------------------------------
        integer                                 :: i, j, ILB, IUB, JLB, JUB

        !Begin-----------------------------------------------------------------


        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB

        do j=JLB, JUB
        do i=ILB, IUB
            
            if (Me%ExternalVar%WaterPoints2D(i, j) == WaterPoint) then

                Me%WaveDirection%Field(i, j) = atan2(Me%ExternalVar%WindVelocityY(i,j), &
                                                     Me%ExternalVar%WindVelocityX(i,j)) &
                                                     * 180. / Pi

            end if
        
        end do
        end do


            
    end subroutine ComputeWaveDirection

    !--------------------------------------------------------------------------
    
    !This subroutine compute the radiation stresses based in X
    !Xia et al. - Coastal Engineering 51 (2004) 309-321.
    real function ComputeRadiationStress(D, z, Om, Option, i, j)

        !Arguments-------------------------------------------------------------
        real                                    :: D, z, Om
        integer                                 :: Option, i, j 

        !Local-----------------------------------------------------------------
        real                                    :: S, E, k, O

        !Begin-----------------------------------------------------------------

        E = Me%En(i,j)
        k = Me%k(i,j)
        O = Me%WaveDirection%Field(i,j)


        if     (Option == Sxx_) then

            S = E * k /  sinh(2*k*D) * (cosh(2*k*(z+D))+1) * cos(O-Om)**2 - &
                E * k /  sinh(2*k*D) * (cosh(2*k*(z+D))-1) - E*Z/D**2  + &
                E *(k*(z+D)*sinh(k*(z+D))) / (D*cosh(k*D))             - &
                E / D * (1 - cosh(k*(z+D))/cosh(k*D))

        elseif (Option == Syy_) then

            S = E * k /  sinh(2*k*D) * (cosh(2*k*(z+D))+1) * sin(O-Om)**2 - &
                E * k /  sinh(2*k*D) * (cosh(2*k*(z+D))-1) - E*Z/D**2  + &
                E *(k*(z+D)*sinh(k*(z+D))) / (D*cosh(k*D))             - &
                E / D * (1 - cosh(k*(z+D))/cosh(k*D))

        elseif (Option == Sxy_ .or. Option == Syx_) then

            S = E * k /  sinh(2*k*D) * (cosh(2*k*(z+D))+1) * sin(O-Om)*cos(O-Om)

        else

            stop 'ComputeRadiationStress - ModuleWaves - ERR10'

        endif

        ComputeRadiationStress = S


    end function ComputeRadiationStress
            
    !-----------------------------------------------------
    
    subroutine Output_Results_HDF

        !External--------------------------------------------------------------
        integer                            :: STAT_CALL
        real                               :: Year, Month, Day, Hour, Minute, Second
         
        !Local-----------------------------------------------------------------

        integer                            :: OutPutNumber
        integer, dimension(6    )          :: TimeAux
        real,    dimension(6    ), target  :: AuxTime
        real,    dimension(:    ), pointer :: TimePtr
        integer                            :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        real(8)                            :: AuxPeriod, TotalTime
        type(T_Time)                       :: Aux
 
        !----------------------------------------------------------------------
        
        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 
        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB
        
        OutPutNumber = Me%OutPut%NextOutPut


TOut:   if (Me%ActualTime >= Me%OutPut%OutTime(OutPutNumber)) then

            !Gets OpenPoints2D
            call GetOpenPoints2D (Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D,     &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaves - ERR01'

            if (Me%ExternalVar%BackTracking) then
                OutPutNumber = Me%OutPut%Number - OutPutNumber + 1 
            endif 
            
            
            if (Me%ExternalVar%BackTracking) then  
                TotalTime = Me%EndTime      - Me%BeginTime                  
                AuxPeriod = Me%ActualTime   - Me%BeginTime
                AuxPeriod = TotalTime       - AuxPeriod
                
                Aux = Me%BeginTime + AuxPeriod
            else
                Aux = Me%ActualTime
            endif


            call ExtractDate(Aux,                                                       &
                             Year = Year, Month  = Month,  Day    = Day,                &
                             Hour = Hour, Minute = Minute, Second = Second)

            TimeAux(1) = int(Year  )
            TimeAux(2) = int(Month )
            TimeAux(3) = int(Day   )
            TimeAux(4) = int(Hour  )
            TimeAux(5) = int(Minute)
            TimeAux(6) = int(Second)

            !Writes current time
            call ExtractDate   (Me%ActualTime,                                          &
                                AuxTime(1), AuxTime(2), AuxTime(3),                     &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
        
            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleWaves - ERR10'

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS",    &
                                 Array1D = TimePtr, OutputNumber = OutPutNumber,        &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleWaves - ERR20'
        
            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5,                                           &
                                  Me%WorkSize%ILB,                                      &
                                  Me%WorkSize%IUB,                                      &
                                  Me%WorkSize%JLB,                                      &
                                  Me%WorkSize%JUB,                                      &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaves - ERR30'

            !Writes OpenPoints
            call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints2D",        &
                                 "-", Array2D = Me%ExternalVar%OpenPoints2D,            &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaves - ERR40'
        
            if (Me%WavePeriod%OutputHDF) then

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%WavePeriod%ID%Name),&
                                     trim(Me%WavePeriod%ID%Name), "s",                  &
                                     Array2D      = Me%WavePeriod%Field,                &
                                     OutputNumber = OutPutNumber,                       &
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OutPut_Results_HDF - ModuleWaves - ERR50'
            endif

            if (Me%WaveHeight%OutputHDF) then

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%WaveHeight%ID%Name),&
                                     trim(Me%WaveHeight%ID%Name), "m",                  &
                                     Array2D      = Me%WaveHeight%Field,                &
                                     OutputNumber = OutPutNumber,                       &
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OutPut_Results_HDF - ModuleWaves - ERR60'
        
            endif

            if (Me%WaveDirection%OutputHDF) then

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%WaveDirection%ID%Name),&
                                     trim(Me%WaveDirection%ID%Name), "o",               &
                                     Array2D      = Me%WaveDirection%Field,             &
                                     OutputNumber = OutPutNumber,                       &
                                     STAT         = STAT_CALL)                      
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OutPut_Results_HDF - ModuleWaves - ERR70'                 
                                                                                    
            endif                                                                   
                                                                                    
            if (Me%RadiationStressX%OutputHDF) then                                 
                                                                                    
                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%RadiationStressX%ID%Name),&
                                     trim(Me%RadiationStressX%ID%Name), "Pa",           &
                                     Array2D      = Me%RadiationStressX%Field,          &
                                     OutputNumber = OutPutNumber,                       &
                                     STAT         = STAT_CALL)                      
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OutPut_Results_HDF - ModuleWaves - ERR80'                 
        
            endif

            if (Me%RadiationStressY%OutputHDF) then

                call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%RadiationStressY%ID%Name),&
                                     trim(Me%RadiationStressY%ID%Name), "Pa",           &
                                     Array2D      = Me%RadiationStressY%Field,          &
                                     OutputNumber = OutPutNumber,                       &
                                     STAT         = STAT_CALL)                      
                if (STAT_CALL /= SUCCESS_)                                              &
                    stop 'OutPut_Results_HDF - ModuleWaves - ERR90'                 
        
            endif

            !Writes everything to disk
            call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) &
                stop 'OutPut_Results_HDF - ModuleWaves - ERR100'


            Me%OutPut%NextOutPut = Me%OutPut%NextOutPut + 1

            !UnGets OpenPoints2D
            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%ExternalVar%OpenPoints2D,   &
                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleWaves - ERR110'

        endif  TOut    


    end subroutine Output_Results_HDF
      
    !--------------------------------------------------------------------------

    subroutine OutPut_TimeSeries

        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
       
        !Begin-----------------------------------------------------------------
        
              
        if (Me%WavePeriod%TimeSerieOn) then                                        
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D  = Me%WavePeriod%Field,                          &
                                STAT    = STAT_CALL)                              
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OutPut_TimeSeries - ModuleWaves - ERR1'                    
        endif                                                                     
                                                                                
                                                                                
        if (Me%WaveHeight%TimeSerieOn) then                                     
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D  = Me%WaveHeight%Field,                          &
                                STAT    = STAT_CALL)                              
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OutPut_TimeSeries - ModuleWaves - ERR10'                    
        endif                                                                                                       
                                                                                
                                                                                  
        if (Me%WaveDirection%TimeSerieOn) then                                    
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D  = Me%WaveDirection%Field,                       &
                                STAT    = STAT_CALL)                              
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OutPut_TimeSeries - ModuleWaves - ERR20'                    
        endif                                                                   
                                                                                
                                                                                
        if (Me%RadiationStressX%TimeSerieOn) then                               
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D  = Me%RadiationStressX%Field,                    &
                                STAT    = STAT_CALL)                            
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OutPut_TimeSeries - ModuleWaves - ERR30'                  
        endif                                                                   
                                                                                
                                                                                
        if (Me%RadiationStressY%TimeSerieOn) then                               
            call WriteTimeSerie(Me%ObjTimeSerie,                                        &
                                Data2D  = Me%RadiationStressY%Field,                    &
                                STAT    = STAT_CALL)                            
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'OutPut_TimeSeries - ModuleWaves - ERR40'                  
        endif                                                                   
                                                                                
                                                                                
    end subroutine OutPut_TimeSeries                                            
                                                                                
                                                                                
    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillWaves(WavesID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: WavesID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers, STAT_CALL           

        !Begin-------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(WavesID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mWaves_,  Me%InstanceID)

            if (nUsers == 0) then

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillWaves - ModuleWaves - ERR10'

                nUsers = DeassociateInstance(mHORIZONTALMAP_,   Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillWaves - ModuleWaves - ERR20'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillWaves - ModuleWaves - ERR30'

                nUsers = DeassociateInstance(mGEOMETRY_,        Me%ObjGeometry)
                if (nUsers == 0) stop 'KillWaves - ModuleWaves - ERR40'

                nUsers = DeassociateInstance(mGRIDDATA_,        Me%ObjGridData)
                if (nUsers == 0) stop 'KillWaves - ModuleWaves - ERR50'

                
                if (Me%WavePeriod%ON      ) then
                
                    if (Me%WavePeriod%ID%SolutionFromFile) then
                        call KillFillMatrix(Me%WavePeriod%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillWaves - ModuleWaves - ERR60'
                    endif

                    deallocate(Me%WavePeriod%Field) 

                endif

                if (Me%WaveHeight%ON      ) then

                    if (Me%WaveHeight%ID%SolutionFromFile) then
                        call KillFillMatrix(Me%WaveHeight%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillWaves - ModuleWaves - ERR70'
                    endif

                    deallocate(Me%WaveHeight%Field)

                endif

                if (Me%WaveDirection%ON      ) then
                
                    if (Me%WaveDirection%ID%SolutionFromFile) then
                        call KillFillMatrix(Me%WaveDirection%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillWaves - ModuleWaves - ERR80'
                    endif

                    deallocate(Me%WaveDirection%Field) 

                endif


                if (Me%RadiationStressX%ON) then

                    deallocate(Me%RadiationStressX%Field)

                    if (Me%RadiationStressX%ID%SolutionFromFile) then
                        call KillFillMatrix(Me%RadiationStressX%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillWaves - ModuleWaves - ERR90'
                    endif
                endif

                if (Me%RadiationStressY%ON) then

                    deallocate(Me%RadiationStressY%Field)

                    if (Me%RadiationStressY%ID%SolutionFromFile) then
                        call KillFillMatrix(Me%RadiationStressY%ID%ObjFillMatrix, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'KillWaves - ModuleWaves - ERR100'
                    endif
                
                endif

                if (Me%ParametersON) then

                    deallocate(Me%Abw)
                    deallocate(Me%Ubw)
                    deallocate(Me%WaveLength)

                endif

                if (Me%Wavegen_type.eq.CEQUALW2) then

                    deallocate(Me%Fetch)

                endif
                
                !Kills the TimeSerie
                if (Me%Output%TimeSerie) then
                    call KillTimeSerie(Me%ObjTimeSerie, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'KillWaves - ModuleWaves - ERR110'
                endif
                
                if (associated(Me%OutPut%OutTime)) then
                    deallocate(Me%OutPut%OutTime, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) &
                        stop 'KillWaves - ModuleWaves - ERR120'
                    nullify   (Me%OutPut%OutTime)
                end if
                
                if (Me%Output%HDF) then
                    call KillHDF5 (Me%ObjHDF5, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'KillWaves - ModuleWaves - ERR130'
                endif

                !Deallocates Instance
                call DeallocateInstance ()

                WavesID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

    end subroutine KillWaves
        
    !------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Waves), pointer          :: AuxObjWaves
        type (T_Waves), pointer          :: PreviousObjWaves
        
        !Begin-----------------------------------------------------------------

        !Updates pointers
        if (Me%InstanceID == FirstObjWaves%InstanceID) then
            FirstObjWaves => FirstObjWaves%Next
        else
            PreviousObjWaves => FirstObjWaves
            AuxObjWaves      => FirstObjWaves%Next
            do while (AuxObjWaves%InstanceID /= Me%InstanceID)
                PreviousObjWaves => AuxObjWaves
                AuxObjWaves      => AuxObjWaves%Next
            enddo

            !Now update linked list
            PreviousObjWaves%Next => AuxObjWaves%Next

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

    subroutine Ready (ObjWaves_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjWaves_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjWaves_ID > 0) then
            call LocateObjWaves (ObjWaves_ID)
            ready_ = VerifyReadLock (mWaves_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjWaves (WavesID)

        !Arguments-------------------------------------------------------------
        integer                                     :: WavesID

        !Local-----------------------------------------------------------------

        Me => FirstObjWaves
        do while (associated (Me))
            if (Me%InstanceID == WavesID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleWaves - LocateObjWaves - ERR01'

    end subroutine LocateObjWaves

    !--------------------------------------------------------------------------

end module ModuleWaves

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------







