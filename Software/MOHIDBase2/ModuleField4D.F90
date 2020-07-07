!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Field4D
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : November 2011
! REVISION      : Paulo Leitão
! DESCRIPTION   : Module to Read property fields from self describing files like NetCDF and HDF5
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

Module ModuleField4D

#ifndef _NOT_IEEE_ARITHMETIC
    use ieee_arithmetic
#endif


    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions,        only : InterpolateMatrix2DInTime,                       &
                                       InterpolateMatrix3DInTime,                       &
                                       InterpolateAngle2DInTime,                        &
                                       InterpolateAngle3DInTime,                        &
                                       SetMatrixValue, ConstructPropertyID,             &
                                       FillMatrix2D, LinearInterpolation,               &
                                       FillMatrix3D, InterpolateProfileR8,              &
                                       CheckAlternativeTidalCompNames,                  &
                                       AmpPhase_To_Complex
    use ModuleDrawing,          only : ArrayPolygonWindow           
    use ModuleHorizontalGrid,   only : GetHorizontalGridSize, ConstructHorizontalGrid,  &
                                       WriteHorizontalGrid,                             &
                                       GetXYInsideDomain, GetXYCellZ, GetZCoordinates,  &
                                       KillHorizontalGrid, UnGetHorizontalGrid
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData, KillGridData
    use ModuleHorizontalMap 
    use ModuleGeometry,         only : ConstructGeometry, GetGeometrySize,              &
                                       ComputeInitialGeometry, GetGeometryDistances,    &
                                       ComputeVerticalGeometry, UnGetGeometry, KillGeometry
    use ModuleMap 
    use ModuleHDF5,             only : ConstructHDF5, HDF5ReadWindow,                   &
                                       GetHDF5FileAccess, GetHDF5GroupNumberOfItems,    &
                                       HDF5SetLimits, GetHDF5ArrayDimensions, KillHDF5, &
                                       HDF5WriteData, HDF5FlushMemory, HDF5WriteData,   &
                                       GetHDF5GroupExist, GetHDF5DataSetExist,          &
                                       GetHDF5ArrayDim, HDF5ReadData,                   &
                                       GetHDF5AllDataSetsOK
#ifndef _NO_NETCDF                                       
    ! Manages NetCDF files
    use ModuleNetCDF,           only : GetNCDFFileAccess, ConstructNETCDF,              &
                                       NETCDFReadGrid2D, NETCDFReadTime,                &
                                       NETCDFGetDimensions, NETCDFReadData,             &
                                       NETCDFReadVert, NETCDFWithVert
#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif
#endif
    use ModuleTimeSerie
    use ModuleTask2000,         only : Task2000Level               
                                       
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructField4D
    private ::      AllocateInstance
    private ::      ConstructFile
    private ::      ReadGridFromFile
    private ::      ReadBathymFromFile    
    private ::      ReadMap2DFromFile
    private ::      ReadGeometryFromFile    
    private ::      ReadMap3DFromFile    
    private ::      ReadOptions
    private ::          Generic4thDimension
    private ::      ConstructPropertyField
    private ::          HDF5TimeInstant 
    private ::          ReadValues2D
    private ::          ReadValues3D

    !Selector
    public  :: GetField4DTimeLimits
    public  :: GetField4DNumberOfInstants
    public  :: GetField4DGeneric4DValue    
    public  :: GetField4DInstant
    public  :: GetField4DSize2D
    public  :: GetField4DSize3D   
    public  :: GetField4DHarmonicsON
    public  :: GetField4DHarmonicsNumber
    public  :: GetField4DHarmonicsName
    public  :: GetField4DBathymID
    public  :: GetField4DGridID
    public  :: GetField4DIsAngle
    
    !Modifier
    public  :: ModifyField4D
    public  :: ModifyField4DXYZ
    private ::      ModifyInput2D
    private ::      ModifyInput3D
    public  :: GetBathymXY
    private ::      InterpolateBathym
    public  ::      Interpolater3D

 
    !Destructor
    public  :: KillField4D                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjField4D 
    
    !Interfaces----------------------------------------------------------------


    !Parameter-----------------------------------------------------------------
    integer, parameter                              :: Time_4D_         = 1
    integer, parameter                              :: Generic_4D_      = 2

    integer, parameter                              :: Dim2D            = 2
    integer, parameter                              :: Dim3D            = 3
    integer, parameter                              :: DimUnknown       = -99
    

    !Variable from file
    integer, parameter                              :: None             = 1

    !type of values
    integer, parameter                              :: InterpolatedValues = 1
    integer, parameter                              :: AccumulatedValues  = 2
    integer, parameter                              :: OriginalValues     = 3
    
    
    !Property types
    integer, parameter                              :: Scalar_          = 1
    integer, parameter                              :: VectorX_         = 2
    integer, parameter                              :: VectorY_         = 3

    character(len=9)                                :: char_amplitude_  = 'amplitude'
    character(len=5)                                :: char_phase_      = 'phase'    
    character(len=8)                                :: char_residual_   = 'residual'        

    !Parameter-----------------------------------------------------------------

    !Types---------------------------------------------------------------------

    type T_Harmonics
        logical                                                    :: ON             = .false.
        integer                                                    :: Number         = null_int
        real                                                       :: TimeReference  = null_real
        real                                                       :: ReferenceValue = null_real 
        character(Len=WaveNameLength), dimension(:),       pointer :: WaveGroupName  => null()
        character(Len=WaveNameLength), dimension(:),       pointer :: WaveName       => null()
        real,                          dimension(:,:,:,:), pointer :: Phase3D        => null()
        real,                          dimension(:,:,:,:), pointer :: Phase3DReal    => null()
        real,                          dimension(:,:,:,:), pointer :: Phase3DImag    => null()
        real,                          dimension(:,:,:,:), pointer :: Amplitude3D    => null()
        real,                          dimension(:,:,:  ), pointer :: Phase2D        => null()
        real,                          dimension(:,:,:  ), pointer :: Phase2DReal    => null()
        real,                          dimension(:,:,:  ), pointer :: Phase2DImag    => null()        
        real,                          dimension(:,:,:  ), pointer :: Amplitude2D    => null()
        real,                          dimension(:,:,:  ), pointer :: Residual3D     => null()
        real,                          dimension(:,:    ), pointer :: Residual2D     => null()
        
        logical                                                    :: TideStateON    = .false.
        real                                                       :: TideStateDT    = null_real 

        
        logical                                                    :: Extract           = .false.
        logical                                                    :: ExtractAmp        = .false.
        logical                                                    :: ExtractPhaseReal  = .false.
        character(Len=WaveNameLength)                              :: ExtractWave       = null_str
        character(Len=StringLength)                                :: FieldNameDim      = null_str
        logical                                                    :: SlowStartON       = .false.
        real                                                       :: SlowStartPeriod   = null_real 
        real                                                       :: DT                = null_real 
        !Xmin, Xmax, Ymin, Ymax
        real, dimension(4)                                         :: NoReadArea        = null_real
    end type T_Harmonics     


    type T_DefaultNames
        character(Len=StringLength) :: bat          = null_str   
        character(Len=StringLength) :: lon_stag     = null_str     
        character(Len=StringLength) :: lat_stag     = null_str 
        character(Len=StringLength) :: mask         = null_str 
        character(Len=StringLength) :: depth_stag   = null_str 
    end type T_DefaultNames


    !Generic 4D
    type T_Generic4D
        logical                            :: ON                    = .false.
        logical                            :: ReadFromTimeSerie     = .false.
        integer                            :: ObjTimeSerie          = null_int
        integer                            :: TimeSerieColumn       = null_int
        real                               :: CurrentValue          = null_real
        integer                            :: Instant               = null_int
        logical                            :: InstantON             = .false. 
    end type 
    
    type T_ExternalVar
        integer, dimension(:,:  ), pointer :: WaterPoints2D         => null()
        integer, dimension(:,:,:), pointer :: WaterPoints3D         => null()
        real,    dimension(:,:  ), pointer :: Bathymetry            => null()
    end type T_ExternalVar
    
    type       T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime          => null()
         integer                                :: TotalOutputs     = null_int
         integer                                :: NextOutPut       = null_int
         logical                                :: Yes              = .false.
         logical                                :: Run_End          = .false.
    end type T_OutPut
    

    type T_File
        character(PathLength)                           :: FileName             = null_str
        character(PathLength)                           :: FieldName            = null_str
        logical                                         :: FieldNameArgument    = .false.
        type (T_Time)                                   :: StartTime
        type (T_Time)                                   :: EndTime        
        integer                                         :: Obj                  = null_int
        
        logical                                         :: FileListON           = .false. 
        character(PathLength), dimension(:), pointer    :: FileNameList         => null()         
        integer                                         :: FilesNumber          = null_int
        integer, dimension(:), pointer                  :: ObjList              => null()
        integer, dimension(:), pointer                  :: ListNInst            => null()        
        integer, dimension(:), pointer                  :: ObjListInst          => null()
        integer, dimension(:), pointer                  :: ListInst             => null()

        
        integer                                         :: NumberOfInstants     = null_int        
        type (T_Time), dimension(:), pointer            :: InstantsDates        => null()    
        
        logical                                         :: CyclicTimeON         = .false.
        logical                                         :: TimeON               = .false.        
        integer                                         :: Form                 = HDF5_
        character(StringLength)                         :: LonStagName          = null_str
        character(StringLength)                         :: LatStagName          = null_str
        character(StringLength)                         :: DepthStagName        = null_str
        character(StringLength)                         :: MaskName             = null_str
        character(StringLength)                         :: BathymName           = null_str
        character(StringLength)                         :: DataSetVert          = "Vertical"
        type (T_DefaultNames)                           :: DefaultNames        
        
        logical                                         :: ReadSZZ              = .false. 
        integer                                         :: SZZLast              = null_int
        
    end type T_File
    
    type T_PropField
        character(PathLength)                       :: FieldName            = null_str
        character(PathLength)                       :: VGroupPath           = null_str
        real                                        :: MultiplyingFactor    = null_real
        logical                                     :: HasMultiplyingFactor = .false.
        real                                        :: AddingFactor         = null_real
        real                                        :: MinValue             = null_real
        real                                        :: MaxValue             = null_real
        logical                                     :: MinValueON           = .false.
        logical                                     :: MaxValueON           = .false.
        logical                                     :: HasAddingFactor      = .false.
        logical                                     :: From2Dto3D           = .false.
        logical                                     :: From3Dto2D           = .false.
        type (T_Time)                               :: NextTime
        type (T_Time)                               :: PreviousTime
        type (T_PropertyID)                         :: ID
        type(T_Generic4D)                           :: Generic4D
        type (T_Harmonics)                          :: Harmonics
        !2D/3D
        integer                                     :: SpaceDim             = null_int
        !Z/U/V
        integer                                     :: TypeZUV              = null_int
        
        logical                                     :: ChangeInTime         = .false.
        logical                                     :: Extrapolate          = .false.
        integer                                     :: ExtrapolateMethod    = null_int
        integer                                     :: InterpolMethod       = null_int
        logical                                     :: DiscardFillValues    = .false.        
        logical                                     :: Zdepths              = .false. 
       
        integer                                     :: ValuesType           = null_int
        real                                        :: Next4DValue          = null_real
        real                                        :: Previous4DValue      = null_real
        integer                                     :: NextInstant          = null_int
        integer                                     :: PreviousInstant      = null_int
        real, dimension(:,:  ), pointer             :: PreviousField2D      => null()
        real, dimension(:,:  ), pointer             :: NextField2D          => null()
        real, dimension(:,:,:), pointer             :: PreviousField3D      => null()
        real, dimension(:,:,:), pointer             :: NextField3D          => null()
        
        real                                        :: MinForDTDecrease     = AllmostZero
        real                                        :: DefaultValue         = null_real
        real                                        :: PredictedDT          = -null_real
        real                                        :: DTForNextEvent       = -null_real
        
        type(T_PropField), pointer                  :: Next                 => null()    
    end type T_PropField            


    private :: T_Field4D
    type       T_Field4D
        integer                                     :: InstanceID   = null_int
        type (T_Size2D)                             :: Size2D
        type (T_Size2D)                             :: WorkSize2D
        type (T_Size3D)                             :: Size3D
        type (T_Size3D)                             :: WorkSize3D
        
        real,    dimension(:, :   ), pointer        :: Matrix2D     => null()  
        real,    dimension(:, :, :), pointer        :: Matrix3D     => null()  
        real,    dimension(:, :, :), pointer        :: Depth3D      => null()  
        
        type (T_Time)                               :: StartTime
        type (T_Time)                               :: EndTime
        real                                        :: DT                   = null_real
        type (T_Time)                               :: CurrentTimeExt
        type (T_Time)                               :: CurrentTimeInt
        type (T_ExternalVar)                        :: ExternalVar
        type (T_PropField), pointer                 :: FirstPropField
        type (T_File      )                         :: File
        
        type (T_OutPut)                             :: OutPut
        
        logical                                     :: Backtracking         = .false.
        
        integer                                     :: ObjEnterData         = 0 
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjBathymetry        = 0    
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjGeometry          = 0    
        integer                                     :: ObjMap               = 0   
        
        integer                                     :: ObjHDF5Out           = 0

        logical                                     :: BuildHorizontalGrid  = .false.
        logical                                     :: BuildBathymetry      = .false.
        logical                                     :: BuildHorizontalMap   = .false.
        logical                                     :: BuildGeometry        = .false.
        logical                                     :: BuildMap             = .false.
        
        logical                                     :: Extrapolate          = .false. 
        integer                                     :: ExtrapolateMethod    = null_int        
        logical                                     :: DiscardFillValues    = .true.     
        

        integer                                     :: MaskDim              = Dim3D
        real                                        :: LatReference         = null_real
        real                                        :: LonReference         = null_real
        logical                                     :: ReadWindowXY         = .false. 
        logical                                     :: WindowWithData       = .false. 
        logical                                     :: ReadWindowJI         = .false. 
        logical                                     :: ReadWindow           = .false. 
        integer                                     :: ClientID             = null_int
        real,    dimension(2,2)                     :: WindowLimitsXY       = null_real
        type (T_Size2D)                             :: WindowLimitsJI        
        
        logical                                     :: CheckHDF5_File       = .false.         
        
        type(T_Field4D), pointer                    :: Next                 => null()
    end type  T_Field4D

    !Global Module Variables
    type (T_Field4D), pointer                       :: FirstObjField4D      => null()
    type (T_Field4D), pointer                       :: Me                   => null()

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructField4D(Field4DID, EnterDataID, ExtractType, TimeID,  FileName, &
                                MaskDim, HorizontalGridID, BathymetryID,                &
                                HorizontalMapID, GeometryID, MapID, LatReference,       &
                                LonReference, WindowLimitsXY, WindowLimitsJI,           &
                                Extrapolate, ExtrapolateMethod, PropertyID, ClientID,   &
                                FileNameList, FieldName, OnlyReadGridFromFile,          &
                                DiscardFillValues, CheckHDF5_File, STAT)

        !Arguments---------------------------------------------------------------
        integer,                                        intent(INOUT) :: Field4DID
        integer,                                        intent(IN )   :: EnterDataID
        integer,                                        intent(IN )   :: ExtractType        
        integer,                                        intent(IN )   :: TimeID
        character(*),                         optional, intent(IN )   :: FileName
        integer,                              optional, intent(IN )   :: MaskDim                
        integer,                              optional, intent(IN )   :: HorizontalGridID
        integer,                              optional, intent(IN )   :: BathymetryID
        integer,                              optional, intent(IN )   :: HorizontalMapID
        integer,                              optional, intent(IN )   :: GeometryID
        integer,                              optional, intent(IN )   :: MapID
        real,                                 optional, intent(IN )   :: LatReference
        real,                                 optional, intent(IN )   :: LonReference        
        real,    dimension(1:2,1:2),          optional, intent(IN )   :: WindowLimitsXY
        type (T_Size2D)            ,          optional, intent(IN )   :: WindowLimitsJI
        logical,                              optional, intent(IN )   :: Extrapolate
        integer,                              optional, intent(IN )   :: ExtrapolateMethod
        type (T_PropertyID),                  optional, intent(IN )   :: PropertyID
        integer,                              optional, intent(IN )   :: ClientID
        character(*), dimension(:), pointer,  optional, intent(IN )   :: FileNameList             
        character(*),                         optional, intent(IN )   :: FieldName             
        logical,                              optional, intent(IN )   :: OnlyReadGridFromFile 
        logical,                              optional, intent(IN )   :: DiscardFillValues
        logical,                              optional, intent(IN )   :: CheckHDF5_File  
        integer,                              optional, intent(OUT)   :: STAT     
        
        !Local-------------------------------------------------------------------
        type (T_PropField), pointer                           :: NewPropField        
        integer                                               :: ready_, STAT_, nUsers, STAT_CALL
        logical                                               :: OnlyReadGridFromFile_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mField4D_)) then
            nullify (FirstObjField4D)
            call RegisterModule (mField4D_) 
        endif
        
        call Ready(Field4DID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            
            if (present(ClientID)) then
                Me%ClientID = ClientID
            else
                Me%ClientID = FillValueInt
            endif
            
            Me%File%FileListON      = .false.
            
            if (present(FileNameList)) then
                if (associated(FileNameList)) then
            
                    Me%File%FileListON      = .true. 
                    Me%File%FilesNumber     = SIZE(FileNameList)  
                    Me%File%FileName        = trim(FileNameList(1))                       
                    
                    allocate(Me%File%FileNameList(Me%File%FilesNumber))
                    allocate(Me%File%ObjList     (Me%File%FilesNumber))
                    allocate(Me%File%ListNInst   (Me%File%FilesNumber))  
                    
                    Me%File%ObjList     (:) = 0
                    Me%File%FileNameList(:) = FileNameList(:)
                    Me%File%ListNInst   (:) = 0
                                    
                endif
            endif                
            
            if (.not. Me%File%FileListON) then
                
                if (present(FileName)) then
                    Me%File%FileName     = trim(FileName)

                else                 
                    stop 'ConstructField4D - ModuleField4D - ERR20' 
                endif
            endif                
            
            
            
            !Set field name
            if (present(FieldName)) then
                Me%File%FieldNameArgument = ON
                Me%File%FieldName         = trim(FieldName)
            endif
            
            if (present(MaskDim)) then
                Me%MaskDim  = MaskDim
            else
                Me%MaskDim  = DimUnknown
            endif
            
            Me%WindowWithData    = .true. 

            if (present(CheckHDF5_File)) then
                Me%CheckHDF5_File  = CheckHDF5_File
            else
                Me%CheckHDF5_File  = .false.
            endif
        
            call ConstructFile(ExtractType)
        
            if (present(Extrapolate)) then
                Me%Extrapolate = Extrapolate
            else
                Me%Extrapolate = .true. 
            endif
            
            if (present(ExtrapolateMethod)) then
                Me%ExtrapolateMethod = ExtrapolateMethod
            else
                Me%ExtrapolateMethod = ExtrapolAverage_
            endif   
            
            if (present(DiscardFillValues)) then
                Me%DiscardFillValues = DiscardFillValues
            else
                Me%DiscardFillValues = .true. 
            endif            
            
            
            
            
            
            if (present(HorizontalGridID)) then
                Me%ObjHorizontalGrid        = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
                
                Me%BuildHorizontalGrid      = .false.
            else
                Me%BuildHorizontalGrid      = .true.
                            
                if (present(LatReference)) then
                    Me%LatReference = LatReference
                else
                    stop 'ConstructField4D - ModuleField4D - ERR30' 
                endif

                if (present(LonReference)) then
                    Me%LonReference = LonReference
                else
                    stop 'ConstructField4D - ModuleField4D - ERR40' 
                endif
                
                if (present(WindowLimitsXY)) then
                    Me%WindowLimitsXY(1:2,1:2) = WindowLimitsXY(1:2,1:2)
                    Me%ReadWindowXY            = .true.
                else
                    Me%ReadWindowXY            = .false.
                endif
                
                if (present(WindowLimitsJI)) then
                    Me%WindowLimitsJI          = WindowLimitsJI
                    Me%ReadWindowJI            = .true.
                else
                    Me%ReadWindowJI            = .false.
                endif
                
                if (Me%ReadWindowXY .and. Me%ReadWindowJI) then                
                    stop 'ConstructField4D - ModuleField4D - ERR50' 
                endif
                
                if (Me%ReadWindowXY .or. Me%ReadWindowJI) then
                    Me%ReadWindow = .true.
                else
                    Me%ReadWindow = .false.
                endif
                
               
                call ReadGridFromFile()
                
            endif
            
            if (present(OnlyReadGridFromFile)) then
                OnlyReadGridFromFile_ = OnlyReadGridFromFile
            else
                OnlyReadGridFromFile_ = .false.
            endif                
                   
OG:         if (.not. OnlyReadGridFromFile_) then
           
                if (Me%File%NumberOfInstants == 1) then
                    Me%WindowWithData = .false.
                endif
                
wwd:            if (Me%WindowWithData) then            
                
                    call GetHorizontalGridSize(HorizontalGridID = Me%ObjHorizontalGrid,         &
                                               Size             = Me%Size2D,                    &
                                               WorkSize         = Me%WorkSize2D,                &
                                               STAT             = STAT_CALL) 
                    if (STAT_CALL/=SUCCESS_) then
                        stop 'ConstructField4D - ModuleField4D - ERR60' 
                    endif
                    
                    
                    if (present(BathymetryID)) then
                        Me%ObjBathymetry            = AssociateInstance (mGRIDDATA_,       BathymetryID    )
                        
                        Me%BuildBathymetry      = .false.
                    else
                        Me%BuildBathymetry      = .true.

                        call ReadBathymFromFile
                    endif

                    if (present(HorizontalMapID)) then
                        Me%ObjHorizontalMap         = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
                        
                        Me%BuildHorizontalMap      = .false.
                    else
                        Me%BuildHorizontalMap      = .true.


                        call ReadMap2DFromFile
                    endif
                    
                    call GetWaterPoints2D(HorizontalMapID   = Me%ObjHorizontalMap,              &
                                          WaterPoints2D     = Me%ExternalVar%WaterPoints2D,     &
                                          STAT              = STAT_CALL) 
                    if (STAT_CALL/=SUCCESS_) then
                        stop 'ConstructField4D - ModuleField4D - ERR70' 
                    endif 

                    if (Me%MaskDim == Dim3D) then

                        if (present(GeometryID)) then
                            Me%ObjGeometry          = AssociateInstance (mGEOMETRY_,       GeometryID      )
                        
                            Me%BuildGeometry      = .false.
                        else
                            Me%BuildGeometry      = .true.

                            call ReadGeometryFromFile
                        endif
                        
                        if (present(MapID)) then
                            Me%ObjMap               = AssociateInstance (mMAP_,            MapID           )
                            
                            call GetWaterPoints3D(Map_ID            = Me%ObjMap,                    &
                                                  WaterPoints3D     = Me%ExternalVar%WaterPoints3D, &
                                                  STAT              = STAT_CALL) 
                            if (STAT_CALL/=SUCCESS_) then
                                stop 'ConstructField4D - ModuleField4D - ERR80' 
                            endif 
                        
                            Me%BuildMap      = .false.
                        else
                            Me%BuildMap      = .true.

                            call ReadMap3DFromFile
                        endif
                   
                    endif
                    
                    call AllocatePropertyField  (NewPropField, Me%MaskDim)
                    
                    if (present(PropertyID)) then
                        NewPropField%ID = PropertyID
                    else
                        call ConstructPropertyID (NewPropField%ID, Me%ObjEnterData, ExtractType)
                    endif
                    
                    call ReadOptions            (NewPropField, ExtractType)     
                    
                    if (NewPropField%Harmonics%ON) then
                        call ConstructPropertyFieldHarmonics (NewPropField)
                    else
                        call ConstructPropertyField          (NewPropField)
                    endif
                    
                    if (Me%OutPut%Yes) then
                        call Open_HDF5_OutPut_File(NewPropField)
                    endif
                    
                    call UnGetHorizontalMap(HorizontalMapID = Me%ObjHorizontalMap,              &
                                            Array           = Me%ExternalVar%WaterPoints2D,     &
                                            STAT            = STAT_CALL) 
                    if (STAT_CALL/=SUCCESS_) then
                        stop 'ConstructField4D - ModuleField4D - ERR90' 
                    endif 
                    
                    if (Me%MaskDim == Dim3D) then
                    
                        call UnGetMap(Map_ID          = Me%ObjMap,                              &
                                      Array           = Me%ExternalVar%WaterPoints3D,           &
                                      STAT            = STAT_CALL) 
                        if (STAT_CALL/=SUCCESS_) then
                            stop 'ConstructField4D - ModuleField4D - ERR100' 
                        endif 
                    
                    endif
                
                endif wwd
                
            endif OG                
                         
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructField4D - ModuleField4D - ERR110' 
             
            
            !Returns ID
            Field4DID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructField4D - ModuleField4D - ERR120' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructField4D

    !----------------------------------------------------------------------
    subroutine Open_HDF5_OutPut_File(NewPropField)
        
        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                 :: NewPropField

        !Local-----------------------------------------------------------------
        real,    pointer, dimension(:, :   )        :: Bathymetry
        character (Len = PathLength)                :: FileName
        type(T_Size2D)                              :: WorkSize2D
        integer                                     :: STAT_CALL
        integer                                     :: WorkILB, WorkIUB
        integer                                     :: WorkJLB, WorkJUB
        integer                                     :: WorkKLB, WorkKUB
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------
        !Bounds

        FileName = trim(Me%File%Filename)//trim(Me%FirstPropField%FieldName)//".hdf5"

        WorkILB = Me%WorkSize2D%ILB 
        WorkIUB = Me%WorkSize2D%IUB 

        WorkJLB = Me%WorkSize2D%JLB 
        WorkJUB = Me%WorkSize2D%JUB 

        WorkSize2D%ILB = Me%WorkSize2D%ILB            
        WorkSize2D%IUB = Me%WorkSize2D%IUB 
                            
        WorkSize2D%JLB = Me%WorkSize2D%JLB 
        WorkSize2D%JUB = Me%WorkSize2D%JUB 

        !Gets a pointer to Bathymetry
        call GetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR10'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        Me%ObjHDF5Out = 0

        !Opens HDF5 File
        call ConstructHDF5      (Me%ObjHDF5Out, trim(FileName),                         &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR20'
        
        !Write the Horizontal Grid
        call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5Out,                   &
                                 WorkSize = WorkSize2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR30'
        

        !Sets limits for next write operations
        call HDF5SetLimits   (Me%ObjHDF5Out, WorkILB, WorkIUB, WorkJLB,                 &
                              WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR40'


        !Writes the Grid
        call HDF5WriteData   (Me%ObjHDF5Out, "/Grid", "Bathymetry", "m",                &
                              Array2D = Bathymetry,                                     &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR50'
        
i0:     if      (NewPropField%SpaceDim == Dim2D)then

            call HDF5WriteData   (Me%ObjHDF5Out, "/Grid", "WaterPoints2D", "-",         &
                                  Array2D = Me%ExternalVar%WaterPoints2D,               &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR60'


        elseif  (NewPropField%SpaceDim == Dim3D) then

            WorkKLB = Me%WorkSize3D%KLB
            WorkKUB = Me%WorkSize3D%KUB
        
            !Sets limits for next write operations
            call HDF5SetLimits   (Me%ObjHDF5Out, WorkILB, WorkIUB, WorkJLB,                       &
                                  WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR70'
            
            call HDF5WriteData   (Me%ObjHDF5Out, "/Grid", "WaterPoints3D", "-",         &
                                  Array3D = Me%ExternalVar%WaterPoints3D,               &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR80'
        
        endif i0
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR90'


        !Ungets the Bathymetry
        call UngetGridData (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_OutPut_File - ModuleField4D - ERR100'

        !----------------------------------------------------------------------

    end subroutine Open_HDF5_OutPut_File

    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine ReadGridFromFile()

        !Arguments-------------------------------------------------------------

                                                    
        !Local-----------------------------------------------------------------
#ifndef _NO_NETCDF            
        real(8),   pointer, dimension(:,:)      :: LatR8, LonR8, LatStagR8, LonStagR8
#endif        
        real,      pointer, dimension(:,:)      :: Lat, Lon, LatStag, LonStag        
        real,      pointer, dimension(:,:)      :: LatStagW, LonStagW
        real,   pointer, dimension(:  )         :: XXDummy, YYDummy
        integer, dimension(:,:), pointer        :: WindowDomain        
        integer                                 :: Imax, Jmax, STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB, i, j

        !Begin-----------------------------------------------------------------
        
        nullify(XXDummy)
        nullify(YYDummy)

       !Get grid horizontal space dimensions
        if      (Me%File%Form == HDF5_  ) then
            call GetHDF5ArrayDimensions (HDF5ID = Me%File%Obj, GroupName = "/Grid",     &
                                        ItemName = trim(Me%File%BathymName),            &
                                        Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
                                        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR10'
#ifndef _NO_NETCDF                       
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, STAT = STAT_CALL)
            
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR20'            
#endif            
        endif

       
        allocate(Lat     (0:Imax+1,0:Jmax+1))
        allocate(Lon     (0:Imax+1,0:Jmax+1))
        allocate(LatStag(0:Imax+1,0:Jmax+1))
        allocate(LonStag(0:Imax+1,0:Jmax+1))
        
            
       !Read horizontal grid
        if      (Me%File%Form == HDF5_  ) then
            
            call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = 1, IUB = Imax+1,           &
                                                       JLB = 1, JUB = Jmax+1,           &
                                 STAT   = STAT_CALL)                                     
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR30'
                                        
            call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                            &
                              GroupName     = "/Grid",                                  &
                              Name          = trim(Me%File%LatStagName),                &
                              Array2D       = LatStag,                                  &
                              STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR40'
            
            call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                            &
                              GroupName     = "/Grid",                                  &
                              Name          = trim(Me%File%LonStagName),                &
                              Array2D       = LonStag,                                  &
                              STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR50'
            
#ifndef _NO_NETCDF            
        else if (Me%File%Form == NetCDF_) then
        
            allocate(LatR8    (0:Imax+1,0:Jmax+1))
            allocate(LonR8    (0:Imax+1,0:Jmax+1))
            allocate(LatStagR8(0:Imax+1,0:Jmax+1))
            allocate(LonStagR8(0:Imax+1,0:Jmax+1))        
            
            LatR8    (:,:) = FillValueReal
            LonR8    (:,:) = FillValueReal      
            LatStagR8(:,:) = FillValueReal
            LonStagR8(:,:) = FillValueReal

            call NETCDFReadGrid2D(NCDFID    = Me%File%Obj,                              &
                                  Lat       = LatR8,                                    &
                                  Lon       = LonR8,                                    &
                                  Lat_Stag  = LatStagR8,                                &
                                  Lon_Stag  = LonStagR8,                                &
                                  STAT      = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR60'
            
            do j=0,Jmax+1
            do i=0,Imax+1            
                Lat     (i,j)= real(LatR8    (i,j)) 
                Lon     (i,j)= real(LonR8    (i,j)) 
                LatStag (i,j)= real(LatStagR8(i,j)) 
                LonStag (i,j)= real(LonStagR8(i,j)) 
            enddo
            enddo

            deallocate(LatR8    )
            deallocate(LonR8    )
            deallocate(LatStagR8)
            deallocate(LonStagR8)        
#endif
        endif
        
        
irw:    if (Me%ReadWindow) then

            !WindowLimitsXY(2,1) = South
            !WindowLimitsXY(2,2) = North
            !WindowLimitsXY(1,1) = West
            !WindowLimitsXY(1,2) = East

iXY:        if (Me%ReadWindowXY) then
                allocate (WindowDomain(2,2))
                
                WindowDomain(:,:) = FillValueInt
            
                call ArrayPolygonWindow(XX              = LonStag,                          &
                                        YY              = LatStag,                          &
                                        WIn             = Me%WindowLimitsXY,                &
                                        ILB             = 1,                                &
                                        IUB             = Imax+1,                           &
                                        JLB             = 1,                                &
                                        JUB             = Jmax+1,                           &
                                        WOut            = WindowDomain,                     &
                                        WindowWithData  = Me%WindowWithData)
                                        
                if (.not. Me%WindowWithData) then
                    
                    do j=0,Jmax+1
                    do i=0,Imax+1            
                        Lon     (i,j)= Lon     (i,j) - 360.
                        LonStag (i,j)= LonStag (i,j) - 360. 
                    enddo
                    enddo                

                    call ArrayPolygonWindow(XX              = LonStag,                          &
                                            YY              = LatStag,                          &
                                            WIn             = Me%WindowLimitsXY,                &
                                            ILB             = 1,                                &
                                            IUB             = Imax+1,                           &
                                            JLB             = 1,                                &
                                            JUB             = Jmax+1,                           &
                                            WOut            = WindowDomain,                     &
                                            WindowWithData  = Me%WindowWithData)
                
                    if (.not. Me%WindowWithData) then
                        write(*,*) 'Input file do not intersect the model domain'
                        write(*,*) 'ReadGridFromFile - ModuleField4D - WRN70'
                    endif
                
                endif
                
                ILB = max(WindowDomain(1,1) - 3,1   )
                IUB = min(WindowDomain(1,2) + 3,Imax)
                JLB = max(WindowDomain(2,1) - 3,1   )
                JUB = min(WindowDomain(2,2) + 3,Jmax)
                
                deallocate (WindowDomain)
                
            endif iXY                            
            
iJI:        if (Me%ReadWindowJI) then

                ILB = Me%WindowLimitsJI%ILB
                if (ILB < 1) then
                    stop 'ReadGridFromFile - ModuleField4D - ERR80'
                endif

                IUB = Me%WindowLimitsJI%IUB
                if (IUB > imax) then
                    stop 'ReadGridFromFile - ModuleField4D - ERR90'
                endif

                JLB = Me%WindowLimitsJI%JLB
                if (JLB < 1) then
                    stop 'ReadGridFromFile - ModuleField4D - ERR100'
                endif

                JUB = Me%WindowLimitsJI%JUB
                if (JUB > jmax) then
                    stop 'ReadGridFromFile - ModuleField4D - ERR110'
                endif
                
            endif iJI
            
            
wwd1:       if (Me%WindowWithData) then
                
                allocate(LatStagW(ILB-1:IUB+1,JLB-1:JUB+1))
                LatStagW(ILB-1:IUB+1,JLB-1:JUB+1) = LatStag(ILB-1:IUB+1,JLB-1:JUB+1)

                allocate(LonStagW(ILB-1:IUB+1,JLB-1:JUB+1))
                LonStagW(ILB-1:IUB+1,JLB-1:JUB+1) = LonStag(ILB-1:IUB+1,JLB-1:JUB+1)
            
            endif wwd1

        else irw
        
            ILB = 1
            IUB = Imax
            JLB = 1
            JUB = Jmax
            
            LatStagW => LatStag
            LonStagW => LonStag
            
        endif irw
        
        if (Me%WindowWithData) then
        
            !Builds horizontal grid object
            call ConstructHorizontalGrid(HorizontalGridID   = Me%ObjHorizontalGrid,         &
                                         LatitudeConn       = LatStagW,                     &
                                         LongitudeConn      = LonStagW,                     &
                                         XX                 = XXDummy,                      &
                                         YY                 = YYDummy,                      &
                                         Latitude           = Me%LatReference,              &
                                         Longitude          = Me%LonReference,              &
                                         ILB                = ILB,                          &
                                         IUB                = IUB,                          &
                                         JLB                = JLB,                          &
                                         JUB                = JUB,                          &
                                         STAT               = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR800'

        endif    
        
        deallocate(Lat    )
        deallocate(Lon    )
        deallocate(LatStag)
        deallocate(LonStag)
        if (Me%ReadWindowXY .and. Me%WindowWithData) then
            deallocate(LatStagW)
            deallocate(LonStagW)
        endif
        nullify(XXDummy)
        nullify(YYDummy)

        nullify(Lat    )
        nullify(Lon    )
        nullify(LatStag)
        nullify(LonStag)
        nullify(LatStagW)
        nullify(LonStagW)

    end subroutine ReadGridFromFile         


    !--------------------------------------------------------------------------
    
    subroutine ReadBathymFromFile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        real,      pointer, dimension(:,:)      :: Bathym       
        integer                                 :: ILB, IUB, JLB, JUB,STAT_CALL

        !Begin-----------------------------------------------------------------
        
        allocate(Bathym  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB        
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        
            
       !Read horizontal grid
        if      (Me%File%Form == HDF5_  ) then
            
            call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,            &
                                                       JLB = JLB, JUB = JUB,            &
                                 STAT   = STAT_CALL)                                                    
                                        
            call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                              &
                              GroupName     = "/Grid",                                  &
                              Name          = trim(Me%File%BathymName),                 &
                              Array2D       = Bathym,                                   &
                              STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadBathymFromFile - ModuleField4D - ERR10'                                    
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                Array2D         = Bathym,                               &
                                Name            = trim(Me%File%BathymName),             &
                                ILB             = ILB,                                  &
                                IUB             = IUB,                                  &
                                JLB             = JLB,                                  &
                                JUB             = JUB,                                  &
                                STAT            = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)stop 'ReadBathymFromFile - ModuleField4D - ERR20'
#endif                  
        endif
 
        
        !Builds horizontal grid object
        call ConstructGridData(GridDataID           = Me%ObjBathymetry,                 &
                               HorizontalGridID     = Me%ObjHorizontalGrid,             &
                               TimeID               = Me%ObjTime,                       &
                               InMatrix2D           = Bathym,                           &
                               STAT                 = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_)stop 'ReadBathymFromFile - ModuleField4D - ERR30'


        deallocate(Bathym)
        nullify   (Bathym)

    end subroutine ReadBathymFromFile         

    !--------------------------------------------------------------------------
    
    subroutine ReadMap2DFromFile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer,   pointer, dimension(:,:)      :: Mask2D
        integer,   pointer, dimension(:,:,:)    :: Mask3D
        integer                                 :: ILB, IUB, JLB, JUB, Kmax, STAT_CALL, KLB, KUB
        integer                                 :: ArrayHDF_Dim

        !Begin-----------------------------------------------------------------
        
        allocate(Mask2D  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB        
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        
        
        if (Me%MaskDim == Dim3D) then                                    

            ArrayHDF_Dim = FillValueInt    
        
           !Get grid vertical space dimensions
            if      (Me%File%Form == HDF5_  ) then

                call GetHDF5ArrayDimensions (HDF5ID = Me%File%Obj, GroupName = "/Grid",     &
                                            ItemName = trim(Me%File%MaskName),              &
                                            Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR10'
                
                ArrayHDF_Dim= GetHDF5ArrayDim(HDF5ID = Me%File%Obj, GroupName = "/Grid",   &
                                              ItemName = trim(Me%File%MaskName), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR15'
                
#ifndef _NO_NETCDF               
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, KUB = Kmax, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR20'
#endif                
            endif
            
            allocate(Mask3D  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB,1:Kmax))  
        
                   !Read horizontal grid
            if      (Me%File%Form == HDF5_  ) then
            
            
                if      (ArrayHDF_Dim == 3) then
                
                    call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,        &
                                                               JLB = JLB, JUB = JUB,        &
                                                               KLB = 1, KUB = Kmax,         &
                                         STAT   = STAT_CALL)                                
                                         
                    if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR30'
                    
                    
                                                
                    call HDF5ReadWindow(HDF5ID      = Me%File%Obj,                          &
                                      GroupName     = "/Grid",                              &
                                      Name          = trim(Me%File%MaskName),               &
                                      Array3D       = Mask3D,                               &
                                      STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR40'
                    
            elseif (ArrayHDF_Dim == 2) then

                    call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,        &
                                                               JLB = JLB, JUB = JUB,        &
                                         STAT   = STAT_CALL)                                
                    if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR42'
                    
                    call HDF5ReadWindow(HDF5ID      = Me%File%Obj,                          &
                                      GroupName     = "/Grid",                              &
                                      Name          = trim(Me%File%MaskName),               &
                                      Array2D       = Mask2D,                               &
                                      STAT          = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR44'
                    
                    
            
            else
            
                stop 'ReadMap2DFromFile - ModuleField4D - ERR46'
            
            endif                    
            
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                    Array3D         = Mask3D,                               &
                                    Name            = trim(Me%File%MaskName),               &
                                    ILB             = ILB,                                  &
                                    IUB             = IUB,                                  &
                                    JLB             = JLB,                                  &
                                    JUB             = JUB,                                  &
                                    KLB             = 1,                                    &
                                    KUB             = Kmax,                                 &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR50'
#endif                      
            endif
            
            if (ArrayHDF_Dim /= 2) then
                Mask2D(:,:) = Mask3D(:,:,Kmax)
            endif                
            
            deallocate(Mask3D)
            nullify   (Mask3D)
            

        else        
                    
           !Read horizontal grid
            if      (Me%File%Form == HDF5_  ) then

                KLB = 1
                KUB = 1

                if (trim(Me%File%MaskName) == 'WaterPoints3D') then
                    call GetHDF5ArrayDimensions (HDF5ID = Me%File%Obj, GroupName = "/Grid",     &
                                                ItemName = trim(Me%File%MaskName),              &
                                                Kmax = Kmax, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR60'
                    
                    KUB = Kmax
                    
                endif
                
                
                call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,        &
                                                           JLB = JLB, JUB = JUB,        &
                                                           KLB = KLB, KUB = KUB,        &
                                     STAT   = STAT_CALL)                                   
                                     
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR70'
                
                call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                            &
                                  GroupName     = "/Grid",                                  &
                                  Name          = trim(Me%File%MaskName),                   &
                                  Array2D       = Mask2D,                                   &
                                  STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR80'
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                    Array2D         = Mask2D,                               &
                                    Name            = trim(Me%File%MaskName),               &
                                    ILB             = ILB,                                  &
                                    IUB             = IUB,                                  &
                                    JLB             = JLB,                                  &
                                    JUB             = JUB,                                  &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR90'
#endif    
            endif
        
        endif 
        
        !Builds horizontal grid object
        call ConstructHorizontalMap(HorizontalMapID      = Me%ObjHorizontalMap,         &
                                    GridDataID           = Me%ObjBathymetry,            &
                                    HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                    Points               = Mask2D,                      &
                                    STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR100'


        deallocate(Mask2D)
        nullify   (Mask2D)

    end subroutine ReadMap2DFromFile         

 
     !--------------------------------------------------------------------------
    
    subroutine ReadGeometryFromFile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                 :: ILB, IUB, JLB, JUB, Kmax, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB        
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
                                  
        
       !Get grid vertical space dimensions
        if      (Me%File%Form == HDF5_  ) then
            call GetHDF5ArrayDimensions (HDF5ID = Me%File%Obj, GroupName = "/Grid",     &
                                        ItemName = trim(Me%File%MaskName),              &
                                        Kmax = Kmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadGeometryFromFile - ModuleField4D - ERR10'
#ifndef _NO_NETCDF           
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFGetDimensions (NCDFID = Me%File%Obj, KUB = Kmax, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)stop 'ReadGeometryFromFile - ModuleField4D - ERR20'
#endif            
        endif
            
           
        call ConstructGeometry(GeometryID       = Me%ObjGeometry,                       &
                               GridDataID       = Me%ObjBathymetry,                     &
                               HorizontalGridID = Me%ObjHorizontalGrid,                 &
                               HorizontalMapID  = Me%ObjHorizontalMap,                  &
                               Kmax             = Kmax,                                 &
                               STAT             = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ReadGeometryFromFile - ModuleField4D - ERR60'
        
        call GetGeometrySize(GeometryID         = Me%ObjGeometry,                       &
                             Size               = Me%Size3D,                            &
                             WorkSize           = Me%WorkSize3D,                        &
                             STAT               = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)stop 'ReadGeometryFromFile - ModuleField4D - ERR70'
        

    end subroutine ReadGeometryFromFile         

 
    
    subroutine ReadMap3DFromFile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
         real,      pointer, dimension(:,:,:)    :: SZZ
         integer,   pointer, dimension(:,:,:)    :: mask
         integer,   pointer, dimension(:,:  )    :: Array2D    
         integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL 
         logical                                 :: Exist, Exist1, Exist2
         integer                                 :: ArrayHDF_Dim, k  
         

        !Begin-----------------------------------------------------------------
        
        allocate(SZZ  (Me%Size3D%ILB:Me%Size3D%IUB,Me%Size3D%JLB:Me%Size3D%JUB,Me%Size3D%KLB:Me%Size3D%KUB))  
        allocate(mask (Me%Size3D%ILB:Me%Size3D%IUB,Me%Size3D%JLB:Me%Size3D%JUB,Me%Size3D%KLB:Me%Size3D%KUB))  
        

        ILB = Me%WorkSize3D%ILB
        IUB = Me%WorkSize3D%IUB        
        JLB = Me%WorkSize3D%JLB
        JUB = Me%WorkSize3D%JUB
        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB
        
       !Read horizontal grid
        if      (Me%File%Form == HDF5_  ) then
            
            call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,            &
                                                       JLB = JLB, JUB = JUB,            &
                                                       KLB = KLB-1, KUB = KUB,          &
                                 STAT   = STAT_CALL)                                
                                 
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR10'
            
            call GetHDF5DataSetExist (HDF5ID = Me%File%Obj, DataSetName = "/Grid/VerticalZ/Vertical_00001",&
                                      Exist  = Exist1, STAT = STAT_CALL)                                
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR20'
            

            call GetHDF5DataSetExist (HDF5ID = Me%File%Obj, DataSetName = "/Grid/VerticalZ/VerticalZ_00001",&
                                      Exist  = Exist2, STAT = STAT_CALL)                                
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR25'
            
            Me%File%DataSetVert = "Vertical"            
            
            if (Exist1) then
                
                Exist = .true. 
                
            else
            
                if (Exist2) then

                    Me%File%DataSetVert = "VerticalZ"
                    Exist = .true. 
                
                else
                    Exist = .false. 
                endif

            endif

            
            if (Exist) then
                                        
                call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                        &
                                    GroupName     = "/Grid/VerticalZ",                  &
                                    Name          = trim(Me%File%DataSetVert)//"_00001", &
                                    Array3D       = SZZ,                                &
                                    OffSet3       = 0,                                  &       
                                    STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR30'
                
                Me%File%ReadSZZ = .true.
                
            else
                if (KUB==1) then
                    SZZ(:,:,0) = 1.
                    SZZ(:,:,1) = 0.
                else
                    stop 'ReadMap3DFromFile - ModuleField4D - ERR40'
                endif
                
                Me%File%ReadSZZ = .false.
                
            endif
                            
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        

            call NETCDFReadVert(NCDFID          = Me%File%Obj,                          &
                                VerticalZ       = SZZ,                                  &
                                ILB             = ILB,                                  &
                                IUB             = IUB,                                  &
                                JLB             = JLB,                                  &
                                JUB             = JUB,                                  &
                                KLB             = KLB,                                  &
                                KUB             = KUB,                                  &
                                STAT            = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR50'
#endif                  
        endif


       if      (Me%File%Form == HDF5_  ) then
       
            ArrayHDF_Dim = GetHDF5ArrayDim(HDF5ID       = Me%File%Obj,                  &
                                           GroupName    = "/Grid",                      &
                                           ItemName     = trim(Me%File%MaskName),       &
                                           STAT         = STAT_CALL)
                                           
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR60'       
            
            if          (ArrayHDF_Dim == 3) then
            
                call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,        &
                                                           JLB = JLB, JUB = JUB,        &
                                                           KLB = KLB, KUB = KUB,        &
                                     STAT   = STAT_CALL)                                
                                     
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR70'

                call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                        &
                                  GroupName     = "/Grid",                              &
                                  Name          = trim(Me%File%MaskName),               &
                                  Array3D       = Mask,                                 &
                                  STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR80'
                
            else if     (ArrayHDF_Dim == 2) then
            
                call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,        &
                                                           JLB = JLB, JUB = JUB,        &
                                     STAT   = STAT_CALL)                                
                                     
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR90'
                
                allocate(Array2D(ILB:IUB, JLB:JUB))

                call HDF5ReadWindow(HDF5ID        = Me%File%Obj,                        &
                                  GroupName     = "/Grid",                              &
                                  Name          = trim(Me%File%MaskName),               &
                                  Array2D       = Array2D,                              &
                                  STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR100'            
            
                do k=KLB,KUB
                    Mask(ILB:IUB, JLB:JUB,k) = Array2D(ILB:IUB, JLB:JUB)
                enddo                    
                
                deallocate(Array2D)
                
            else
            
            endif
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array3D         = Mask,                             &
                                    Name            = trim(Me%File%MaskName),           &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    KLB             = KLB,                              &
                                    KUB             = KUB,                              &
                                    STAT            = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR110'
#endif                
        endif
        
        call ConstructMap   (Map_ID             = Me%ObjMap,                            &
                            GeometryID          = Me%ObjGeometry,                       &
                            HorizontalMapID     = Me%ObjHorizontalMap,                  &
                            WaterPoints3D       = Mask,                                 &
                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR120'
        
        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                    &
                              WaterPoints3D     = Me%ExternalVar%WaterPoints3D, &
                              STAT              = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR130'

        call ComputeInitialGeometry(GeometryID      = Me%ObjGeometry,                   &
                                    WaterPoints3D   = Me%ExternalVar%WaterPoints3D,     &
                                    SZZ             = SZZ,                              &
                                    STAT            = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR140'

    end subroutine ReadMap3DFromFile         
 
 
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_Field4D), pointer                         :: NewObjField4D
        type (T_Field4D), pointer                         :: PreviousObjField4D


        !Allocates new instance
        allocate (NewObjField4D)
        nullify  (NewObjField4D%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjField4D)) then
            FirstObjField4D    => NewObjField4D
            Me                    => NewObjField4D
        else
            PreviousObjField4D => FirstObjField4D
            Me                    => FirstObjField4D%Next
            do while (associated(Me))
                PreviousObjField4D  => Me
                Me                     => Me%Next
            enddo
            Me                          => NewObjField4D
            PreviousObjField4D%Next  => NewObjField4D
        endif

        Me%InstanceID = RegisterNewInstance (mField4D_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadOptions(PropField, ExtractType)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField                
        integer                                         :: ExtractType
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL
        integer                                         :: iflag
        logical                                         :: LastGroupEqualField
        integer                                         :: ExtractTypeBlock
        real                                            :: DT       
        !---------------------------------------------------------------------
        
        call GetData(PropField%DefaultValue,                                            &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'DEFAULTVALUE',                                   &
                     default        = 0.,                                               &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR05'             
        
        !Property ValuesType (Interpolated, accumulate, original value) 
        call GetData(PropField%ValuesType,                                              &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'VALUES_TYPE',                                    &
                     default        = InterpolatedValues,                               &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR10'       

        !When to shut down DT?
        call GetData(PropField%MinForDTDecrease, Me%ObjEnterData,  iflag,               &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'MIN_FOR_DT_DECREASE',                            &
                     default        = AllMostZero,                                      &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR20'

        call GetData(PropField%Generic4D%ON,                                            &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = '4D',                                               &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR30'


        if (PropField%Generic4D%ON) call Generic4thDimension(PropField, ExtractType)
        
        PropField%Generic4D%InstantON = .false. 
        

        call GetData(PropField%VGroupPath,                                              &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'VGROUP_PATH',                                      &
                     default      = "/Results",                                         &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR40'

        call GetData(PropField%MultiplyingFactor,                                       &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'MULTIPLYING_FACTOR',                               &
                     default      = 1.,                                                 &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR50'
        
        if (iflag == 1)then
            PropField%HasMultiplyingFactor = .true.
        end if

        call GetData(PropField%AddingFactor,                                            &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'ADDING_FACTOR',                                    &
                     default      = 0.,                                                 &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR60'
        
        if (iflag == 1)then
            PropField%HasAddingFactor = .true.
        end if
        
        !Field NAme was setted by argument?
        if (.not. Me%File%FieldNameArgument) then
            call GetData(PropField%FieldName,                                           &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'FIELD_NAME',                                   &
                         default      = trim(PropField%ID%Name),                        &
                         ClientModule = 'ModuleField4D',                                &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR70'
        
            if (iflag == 0) then
                call GetData(PropField%FieldName,                                       &
                             Me%ObjEnterData , iflag,                                   &
                             SearchType   = ExtractType,                                &
                             keyword      = 'HDF_FIELD_NAME',                           &
                             default      = trim(PropField%ID%Name),                    &
                             ClientModule = 'ModuleField4D',                            &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR80'
            endif
        else
            PropField%FieldName = Me%File%FieldName
        endif
        
        call GetData(LastGroupEqualField,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'LAST_GROUP_EQUAL_FIELD',                           &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR90'

        call GetData(PropField%From2Dto3D,                                              &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FROM_2D_TO_3D',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR100'
        

        call GetData(PropField%From3Dto2D,                                              &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FROM_3D_TO_2D',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR110'

        
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime      = Me%StartTime,                             &
                           EndTime          = Me%EndTime,                               &
                           keyword          = 'OUTPUT_TIME',                            &
                           SearchType       = ExtractType,                              &
                           OutPutsTime      = Me%OutPut%OutTime,                        &
                           OutPutsOn        = Me%OutPut%Yes,                            &
                           OutPutsNumber    = Me%OutPut%TotalOutputs,                   &
                           STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR120'        

        Me%OutPut%NextOutPut = 1
        
       call GetData(PropField%Harmonics%ON,                                             &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'HARMONICS',                                        &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR130'    
        
        if (PropField%Harmonics%ON) then
        
            PropField%From2Dto3D = .true.
        
            call GetData(PropField%Harmonics%Extract,                                   &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'EXTRACT_HARMONICS',                            &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleField4D',                                &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR140'    
            
            call GetData(PropField%Harmonics%FieldNameDim,                              &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'HARMONICS_FIELD_DIM',                          &
                         default      = char_residual_,                                 &
                         ClientModule = 'ModuleField4D',                                &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR150'    
        
            call GetData(PropField%Harmonics%TimeReference,                             &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'TIME_REF',                                     &
                         default      = 0.,                                             &
                         ClientModule = 'ModuleField4D',                                &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR160'
            
            if (iflag == 1) PropField%Harmonics%TimeReference = - PropField%Harmonics%TimeReference
        
            call GetData(PropField%Harmonics%ReferenceValue,                            &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'REF_VALUE',                                    &
                         default      = 0.,                                             &
                         ClientModule = 'ModuleField4D',                                &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR170'            


            call GetData(PropField%Harmonics%TideStateON,                               &
                         Me%ObjEnterData , iflag,                                       &
                         SearchType   = ExtractType,                                    &
                         keyword      = 'TIDE_STATE_ON',                                &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleField4D',                                &
                         STAT         = STAT_CALL)                                      
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR180'                
            
            
            if (PropField%Harmonics%TideStateON) then
                call GetData(PropField%Harmonics%TideStateDT,                           &
                             Me%ObjEnterData , iflag,                                   &
                             SearchType   = ExtractType,                                &
                             keyword      = 'TIDE_STATE_DT',                            &
                             default      = 1800.,                                      &
                             ClientModule = 'ModuleField4D',                            &
                             STAT         = STAT_CALL)                                      
                if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR190'
            endif 
            
            call GetData(PropField%Harmonics%SlowStartPeriod,                           &
                            Me%ObjEnterData , iflag,                                    &
                            SearchType   = ExtractType,                                 &
                            keyword      = 'SLOWSTART',                                 &
                            default      =  FillValueReal,                              &
                            ClientModule = 'ModuleField4D',                             &
                            STAT         = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR191'     
            
            if (iflag == 1) then
                PropField%Harmonics%SlowStartON = .true.
                if (PropField%Harmonics%SlowStartPeriod > Me%EndTime - Me%StartTime) then
                    write(*,*) 'SlowStartPeriod > RunPeriod'
                    write(*,*) 'Property name =', trim(PropField%ID%Name)
                endif
            else
                PropField%Harmonics%SlowStartON = .false.  
            endif
            
            call GetData(PropField%Harmonics%DT,                                        &
                            Me%ObjEnterData , iflag,                                    &
                            SearchType   = ExtractType,                                 &
                            keyword      = 'HARMONICS_DT',                              &
                            default      =  900.,                                       &
                            ClientModule = 'ModuleField4D',                             &
                            STAT         = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR192'  
            
            call GetComputeTimeStep(TimeID = Me%ObjTime, DT = DT, STAT = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR193'
                
            if (PropField%Harmonics%DT < DT) then
                PropField%Harmonics%DT = DT
            endif            

            !1 - Xmin, 2 - Xmax, 3 - Ymin, 4 - Ymax             
            call GetData(PropField%Harmonics%NoReadArea,                                &
                            Me%ObjEnterData , iflag,                                    &
                            SearchType   = ExtractType,                                 &
                            keyword      = 'HARMONICS_NO_READ_AREA',                    &
                            ClientModule = 'ModuleField4D',                             &
                            STAT         = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR194'  
            
            if (iflag==0) then
                PropField%Harmonics%NoReadArea(:) =   null_real
            elseif (iflag==4) then                    
                if (PropField%Harmonics%NoReadArea(1) > PropField%Harmonics%NoReadArea(2)) then
                    stop 'ReadOptions - ModuleField4D - ERR195'  
                endif                    
                if (PropField%Harmonics%NoReadArea(3) > PropField%Harmonics%NoReadArea(4)) then
                    stop 'ReadOptions - ModuleField4D - ERR196'  
                endif                                    
            elseif (iflag/=4) then    
                stop 'ReadOptions - ModuleField4D - ERR197'  
            endif
            
        endif
        
        if (LastGroupEqualField)                                                        &
            PropField%VGroupPath=trim(PropField%VGroupPath)//"/"//trim(PropField%FieldName)            

        call GetData(PropField%SpaceDim,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPACE_DIM',                                        &
                     default      = Me%MaskDim,                                         &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR200'
        
        if (iflag == 0 .and. Me%File%Form == HDF5_) then
        
            if (PropField%Harmonics%ON) then
                call GetHDF5ArrayDimensions (Me%File%Obj, trim(PropField%VGroupPath),   &
                                  PropField%Harmonics%FieldNameDim,                     &
                                  NDim = PropField%SpaceDim, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadOptions - ModuleField4D - ERR210'
            else                
                call GetHDF5ArrayDimensions (Me%File%Obj, trim(PropField%VGroupPath),   &
                                  trim(PropField%FieldName), OutputNumber = 1,          &
                                  NDim = PropField%SpaceDim, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadOptions - ModuleField4D - ERR220'
            endif                
        endif
        if (PropField%Harmonics%ON) then

            if      (ExtractType == FromFile_        ) then
                    
                ExtractTypeBlock = FromBlock_

            elseif  (ExtractType == FromBlock_       ) then   

                ExtractTypeBlock = FromBlockInBlock_
        
            elseif  (ExtractType == FromBlockInBlock_) then
            
                ExtractTypeBlock = FromBlockInBlockInBlock_
            
            else
            
                stop 'ReadOptions - ModuleField4D - ERR230'
            
            endif
            
            call ReadHarmonicWaves(PropField, ExtractTypeBlock)
            
        endif   

        call GetData(PropField%ChangeInTime,                                            &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'CHANGE_IN_TIME',                                   &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR240'


        call GetData(PropField%MinValue,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'MIN_VALUE',                                        &
                     default      = FillValueReal,                                      &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR250'
        
        PropField%MinValueON = .false.
        
        if (iflag==1) then
            PropField%MinValueON = .true.
        endif

        call GetData(PropField%MaxValue,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'MAX_VALUE',                                        &
                     default      = -FillValueReal,                                     &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR260'
        
        PropField%MaxValueON = .false.
        
        if (iflag==1) then
            PropField%MaxValueON = .true.
        endif
        call GetData(PropField%Extrapolate,                                             &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'EXTRAPOLATE',                                      &
                     default      = Me%Extrapolate,                                     &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR270'

        !ExtrapolAverage_ = 1, ExtrapolNearstCell_ = 2
        call GetData(PropField%ExtrapolateMethod,                                       &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'EXTRAPOLATE_METHOD',                               &
                     default      = Me%ExtrapolateMethod,                               &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR280'        
        
        !Bilinear2D_         = 1, NearestNeighbor2D_  = 2
        call GetData(PropField%InterpolMethod,                                          &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'INTERPOLATE_METHOD',                               &
                     default      = Bilinear2D_,                                        &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR290'        
        
        if (PropField%InterpolMethod /= Bilinear2D_ .and.                               &
            PropField%InterpolMethod /= NearestNeighbor2D_) then
            stop 'ReadOptions - ModuleField4D - ERR300'
        endif            
        
        
        call GetData(PropField%Zdepths,                                                 &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'Z_DEPTHS',                                         &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR310'             
        
        call GetData(PropField%DiscardFillValues,                                       &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'DISCARD_FILLVALUES',                               &
                     default      = Me%DiscardFillValues,                               &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR320'        
        

        

        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ObjTime, Me%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR340' 
        
        

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadHarmonicWaves(PropField, ExtractType)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                 :: PropField                
        integer                                     :: ExtractType

        !Local-----------------------------------------------------------------
        character(LEN = StringLength   ), parameter :: block_begin1 = '<beginharmonics>'
        character(LEN = StringLength   ), parameter :: block_end1   = '<endharmonics>'

        character(LEN = StringLength   ), parameter :: block_begin2 = '<<beginharmonics>>'
        character(LEN = StringLength   ), parameter :: block_end2   = '<<endharmonics>>'
        
        character(LEN = StringLength   ), parameter :: block_begin3 = '<<<beginharmonics>>>'
        character(LEN = StringLength   ), parameter :: block_end3   = '<<<endharmonics>>>'
        
        
        integer                                     :: STAT_CALL, ClientNumber, FirstLine, LastLine, i
        logical                                     :: BlockFound
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB, iflag, NW
        
        !Begin-----------------------------------------------------------------
        
        ClientNumber = Me%ClientID

        if      (ExtractType == FromBlock_   ) then
                
            call ExtractBlockFromBuffer(EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = ClientNumber,             &
                                        block_begin         = block_begin1,             &
                                        block_end           = block_end1,               &
                                        BlockFound          = BlockFound,               &
                                        FirstLine           = FirstLine,                &
                                        LastLine            = LastLine,                 &
                                        STAT                = STAT_CALL)
                                        
            if (STAT_CALL /= SUCCESS_) stop 'ReadHarmonicWaves - ModuleField4D - ERR10'
            if (.not. BlockFound     ) stop 'ReadHarmonicWaves - ModuleField4D - ERR20'
            

        elseif  (ExtractType == FromBlockInBlock_  ) then   

            call ExtractBlockFromBlock (EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = ClientNumber,             &
                                        block_begin         = block_begin2,             &
                                        block_end           = block_end2,               &
                                        BlockInBlockFound   = BlockFound,               &
                                        FirstLine           = FirstLine,                &
                                        LastLine            = LastLine,                 &
                                        STAT                = STAT_CALL)
                                        
            if (STAT_CALL /= SUCCESS_) stop 'ReadHarmonicWaves - ModuleField4D - ERR30'
            if (.not. BlockFound     ) stop 'ReadHarmonicWaves - ModuleField4D - ERR40'

        elseif  (ExtractType == FromBlockInBlockInBlock_  ) then   

            call ExtractBlockFromBlockFromBlock(EnterDataID              = Me%ObjEnterData, &
                                                ClientNumber             = ClientNumber,    &
                                                block_begin              = block_begin3,    &
                                                block_end                = block_end3,      &
                                                BlockInBlockInBlockFound = BlockFound,      &
                                                FirstLine                = FirstLine,       &
                                                LastLine                 = LastLine,        &
                                                STAT                     = STAT_CALL)
                                        
            if (STAT_CALL /= SUCCESS_) stop 'ReadHarmonicWaves - ModuleField4D - ERR130'
            if (.not. BlockFound     ) stop 'ReadHarmonicWaves - ModuleField4D - ERR140'

    
        else
        
            stop 'ReadHarmonicWaves - ModuleField4D - ERR50'
        
        endif
        
        PropField%Harmonics%Number = LastLine - FirstLine - 1 
        NW                         = PropField%Harmonics%Number

        allocate(PropField%Harmonics%WaveName     (PropField%Harmonics%Number))

        allocate(PropField%Harmonics%WaveGroupName(PropField%Harmonics%Number))        
        
i0:     if(PropField%SpaceDim == Dim2D) then

            ILB = Me%Size2D%ILB
            IUB = Me%Size2D%IUB
            JLB = Me%Size2D%JLB
            JUB = Me%Size2D%JUB

            allocate(PropField%Harmonics%Phase2D    (ILB:IUB, JLB:JUB, 1:NW))
            allocate(PropField%Harmonics%Phase2DReal(ILB:IUB, JLB:JUB, 1:NW))            
            allocate(PropField%Harmonics%Phase2DImag(ILB:IUB, JLB:JUB, 1:NW))                        
            allocate(PropField%Harmonics%Amplitude2D(ILB:IUB, JLB:JUB, 1:NW))
            allocate(PropField%Harmonics%Residual2D (ILB:IUB, JLB:JUB))
            
            PropField%Harmonics%Phase2D    (:,:, :) = FillValueReal
            PropField%Harmonics%Phase2DReal(:,:, :) = FillValueReal            
            PropField%Harmonics%Phase2DImag(:,:, :) = FillValueReal            
            PropField%Harmonics%Amplitude2D(:,:, :) = FillValueReal
            PropField%Harmonics%Residual2D (:,:)    = FillValueReal        

        else i0

            ILB = Me%Size3D%ILB
            IUB = Me%Size3D%IUB
            JLB = Me%Size3D%JLB
            JUB = Me%Size3D%JUB
            KLB = Me%Size3D%KLB
            KUB = Me%Size3D%KUB

            allocate(PropField%Harmonics%Phase3D    (ILB:IUB, JLB:JUB, KLB:KUB, 1:NW))
            allocate(PropField%Harmonics%Phase3DReal(ILB:IUB, JLB:JUB, KLB:KUB, 1:NW))
            allocate(PropField%Harmonics%Phase3DImag(ILB:IUB, JLB:JUB, KLB:KUB, 1:NW))            
            allocate(PropField%Harmonics%Amplitude3D(ILB:IUB, JLB:JUB, KLB:KUB, 1:NW))
            allocate(PropField%Harmonics%Residual3D (ILB:IUB, JLB:JUB, KLB:KUB))

            PropField%Harmonics%Phase3D    (:,:,:,:) = FillValueReal
            PropField%Harmonics%Phase3DReal(:,:,:,:) = FillValueReal
            PropField%Harmonics%Phase3DImag(:,:,:,:) = FillValueReal            
            PropField%Harmonics%Amplitude3D(:,:,:,:) = FillValueReal
            PropField%Harmonics%Residual3D (:,:,:  ) = FillValueReal
            
        endif i0
        
        do i = 1, PropField%Harmonics%Number
        
            call GetData(PropField%Harmonics%WaveGroupName(i),                          &
                         Me%ObjEnterData,  iflag, Buffer_Line  = FirstLine + i,         &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadHarmonicWaves - ModuleField4D - ERR60'
            
            call CheckAlternativeTidalCompNames (TidalName      = PropField%Harmonics%WaveGroupName(i), &
                                                 MohidTidalName = PropField%Harmonics%WaveName(i))
            
        enddo
        
        if      (ExtractType == FromBlock_   ) then
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadHarmonicWaves - ModuleField4D - ERR70'
        endif
        


    end subroutine ReadHarmonicWaves

    !--------------------------------------------------------------------------

                
    !--------------------------------------------------------------------------


    subroutine Generic4thDimension(PropField, ExtractType)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer        :: PropField
        integer                            :: ExtractType

        !Local-----------------------------------------------------------------
        character(len = PathLength)        :: Filename
        integer                            :: STAT_CALL, iflag

        !----------------------------------------------------------------------

        call GetData(PropField%Generic4D%ReadFromTimeSerie,                             &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = '4D_TIME_SERIE',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'Generic4thDimension - ModuleField4D - ERR10'


        if (PropField%Generic4D%ReadFromTimeSerie) then

            call GetData(Filename, Me%ObjEnterData, iflag,                              &
                         keyword        = 'GENERIC_4D_FILENAME',                        &  
                         SearchType     = ExtractType,                                  &
                         ClientModule   = 'ModuleField4D',                              &
                         default        = "******.***",                                 &
                         STAT           = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'Generic4thDimension  - ModuleField4D - ERR20'
            if (iflag == 0) stop 'Generic4thDimension  - ModuleField4D - ERR30'

            call GetData(PropField%Generic4D%TimeSerieColumn, Me%ObjEnterData, iflag,   &
                         keyword        = 'TIME_SERIE_COLUMN',                          &  
                         SearchType     = ExtractType,                                  &
                         ClientModule   = 'ModuleField4D',                              &
                         default        = 2,                                            &
                         STAT           = STAT_CALL)            
            if (STAT_CALL /= SUCCESS_) stop 'Generic4thDimension  - ModuleField4D - ERR40'

            !Starts Time Serie
            call StartTimeSerieInput(PropField%Generic4D%ObjTimeSerie,                  &
                                     FileName,                                          &
                                     CheckDates =.false.,                               &
                                     STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Generic4thDimension - ModuleField4D - ERR50'

        endif


     end subroutine Generic4thDimension
   !----------------------------------------------------------------------------

    subroutine ConstructFile(ExtractType)

        !Arguments-------------------------------------------------------------
        integer,                intent(IN )                     :: ExtractType        
        
        !Local-----------------------------------------------------------------
#ifndef _NO_NETCDF
        type (T_Time)                                           :: AuxTime
        real,    dimension(6)                                   :: InitialDate
        real(8), dimension(:), pointer                          :: Instants
        integer                                                 :: NCDF_READ
#endif        
        integer                                                 :: STAT_CALL, i, HDF5_READ, iflag
        logical                                                 :: exist, exist3D, exist2D, exist2D_2
        integer                                                 :: n, j, k, iaux

        !Begin-----------------------------------------------------------------



        inquire (file=trim(Me%File%FileName), exist = exist)
        if (.not. exist) then
            write(*,*)'Could not find file '//trim(Me%File%FileName)
            stop 'ConstructFile - ModuleField4D - ERR10'
        endif
        
        i = len_trim(Me%File%FileName)
        
        if      (Me%File%FileName(i-4:i) == ".hdf5") then
            Me%File%Form = HDF5_
#ifndef _NO_NETCDF            
        else if (Me%File%FileName(i-2:i) == ".nc"  ) then
            Me%File%Form = NetCDF_
#endif            
        else
            stop 'ConstructFile - ModuleField4D - ERR20'
        endif
        
        if (Me%File%Form == HDF5_) then
        
            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
            
            if (Me%File%FileListON) then

                do n = 1, Me%File%FilesNumber                            

                    call ConstructHDF5 (Me%File%ObjList(n), trim(Me%File%FileNameList(n)), HDF5_READ, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR30'
                    
                    if (Me%CheckHDF5_File) then
        
                        call GetHDF5AllDataSetsOK (HDF5ID   = Me%File%ObjList(n),       &
                                                   STAT     = STAT_CALL)                                      
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR35'
            
                    endif
                    
                enddo
                
                Me%File%Obj = Me%File%ObjList(1)
            else
            
                call ConstructHDF5 (Me%File%Obj, trim(Me%File%FileName), HDF5_READ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR40'
            
                if (Me%CheckHDF5_File) then
        
                    call GetHDF5AllDataSetsOK (HDF5ID   = Me%File%Obj,                  &
                                               STAT     = STAT_CALL)                                      
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR45'
            
                endif
            
            endif
            
            call GetHDF5GroupExist (HDF5ID      = Me%File%Obj,                          &
                                    GroupName   = "/Time",                              &
                                    Exist       = Me%File%TimeON,                       &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR50'
            
            
            if (Me%File%TimeON) then

flo:            if (Me%File%FileListON) then            
                
                    Me%File%NumberOfInstants = 0
                    iaux                     = 0
                
                    do n = 1, Me%File%FilesNumber
            
                        call GetHDF5GroupNumberOfItems(Me%File%ObjList(n), "/Time", iaux, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) then
                            write(*,*) 'FileName = ',trim(Me%File%FileName)
                            stop 'ConstructFile - ModuleField4D - ERR60'
                        endif
                        
                        Me%File%NumberOfInstants = Me%File%NumberOfInstants + iaux
                        Me%File%ListNInst(n) = iaux
                        
                    enddo                        

                    allocate(Me%File%InstantsDates(Me%File%NumberOfInstants))
                    allocate(Me%File%ObjListInst  (Me%File%NumberOfInstants))
                    allocate(Me%File%ListInst     (Me%File%NumberOfInstants))

                    k = 0

                    do n=1, Me%File%FilesNumber
                        do j = 1, Me%File%ListNInst(n) 
                            k = k + 1
                            Me%File%InstantsDates(k) = HDF5TimeInstant(j, HDF5ID = Me%File%ObjList(n))
                            Me%File%ObjListInst  (k) =  Me%File%ObjList(n)
                            Me%File%ListInst     (k) =  j
                            
                            
                            if (k > 1) then
                                if (Me%File%InstantsDates(k) < Me%File%InstantsDates(k-1)) then
                                    stop  'ConstructFile - ModuleField4D - ERR70'
                                endif
                            endif
                        enddo
                    enddo    
                    
                    if (k /= Me%File%NumberOfInstants) then
                        stop 'ConstructFile - ModuleField4D - ERR80'
                    endif
                    
                else flo
                
                    call GetHDF5GroupNumberOfItems(Me%File%Obj, "/Time", &
                                                   Me%File%NumberOfInstants, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'ConstructFile - ModuleField4D - ERR90'
                    endif

                    allocate(Me%File%InstantsDates(Me%File%NumberOfInstants))
                    
                    do i=1, Me%File%NumberOfInstants
                        Me%File%InstantsDates(i) = HDF5TimeInstant(i)
                    enddo

                endif flo

                Me%File%StartTime = Me%File%InstantsDates(1)
                Me%File%EndTime   = Me%File%InstantsDates(Me%File%NumberOfInstants)

            endif            
            
            Me%File%DefaultNames%bat        = 'Bathymetry'
            Me%File%DefaultNames%lon_stag   = 'Longitude'
            Me%File%DefaultNames%lat_stag   = 'Latitude'
 
            call GetHDF5DataSetExist (Me%File%Obj, '/Grid/WaterPoints2D', exist2D, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR100'
                       
            call GetHDF5DataSetExist (Me%File%Obj, '/Grid/WaterPoints', exist2D_2, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR110'           
            
            call GetHDF5DataSetExist (Me%File%Obj, '/Grid/WaterPoints3D', exist3D, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR120'
            
            if (Me%MaskDim == DimUnknown) then
            

                
                if (exist3d) then
                    Me%MaskDim = Dim3D
                endif
                if (exist2D) then
                    Me%MaskDim = Dim2D
                endif
                
                if (exist2D_2) then
                    Me%MaskDim = Dim2D
                endif
                                
                if (exist2D .and. exist3D) then
                    Me%MaskDim = Dim2D
                endif
                
                if (exist2D_2 .and. exist3D) then
                    Me%MaskDim = Dim3D
                endif
                                
            endif                
            
            Me%File%DefaultNames%mask   = 'WaterPoints'
            
            if (Me%MaskDim == Dim3D) then
                Me%File%DefaultNames%mask   = 'WaterPoints3D'
            else
                if (exist2D) then
                    Me%File%DefaultNames%mask   = 'WaterPoints2D'
                endif
                if (exist2D_2) then
                    Me%File%DefaultNames%mask   = 'WaterPoints'
                    write(*,*) Me%File%DefaultNames%mask                    
                endif                                    
            endif

            if (exist3D) then
                Me%File%DefaultNames%mask   = 'WaterPoints3D'
            endif

            if (exist2D) then
                Me%File%DefaultNames%mask   = 'WaterPoints2D'
            endif

            if (exist2D_2) then
                Me%File%DefaultNames%mask   = 'WaterPoints'
            endif
            
            if (exist2D_2 .and. exist3D) then
                Me%File%DefaultNames%mask   = 'WaterPoints'
            endif
            
            

            
            Me%File%DefaultNames%depth_stag = 'VerticalZ'
            
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then

            Me%File%DefaultNames%bat        = 'bathymetry'
            Me%File%DefaultNames%lon_stag   = 'lon_staggered'
            Me%File%DefaultNames%lat_stag   = 'lat_staggered'
            Me%File%DefaultNames%mask       = 'mask'
            Me%File%DefaultNames%depth_stag = 'depth_staggered'
            
        
            call GetNCDFFileAccess (NCDF_READ = NCDF_READ)
        
            call ConstructNETCDF(NCDFID = Me%File%Obj, FileName = trim(Me%File%FileName),&
                                 Access = NCDF_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR130'
            
            if (Me%MaskDim == DimUnknown) then
               
                if (NETCDFWithVert(Me%File%Obj)) then
                    Me%MaskDim = Dim3D
                else
                    Me%MaskDim = Dim2D
                endif
                
            endif             
            
            call NETCDFReadTime(NCDFID = Me%File%Obj, InitialDate = InitialDate,        &
                                nInstants = Me%File%NumberOfInstants,                   &
                                Instants = Instants, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR180'
            
            call SetDate(Time1 = AuxTime, Year = InitialDate(1), Month = InitialDate(2), &
                                           Day = InitialDate(3), Hour  = InitialDate(4), &
                                        Minute = InitialDate(5), Second= InitialDate(6))
            
            allocate(Me%File%InstantsDates(Me%File%NumberOfInstants))
            
            do i=1, Me%File%NumberOfInstants
                Me%File%InstantsDates(i) = AuxTime + Instants(i)
            enddo            
            
            Me%File%StartTime = Me%File%InstantsDates(1                       )
            Me%File%EndTime   = Me%File%InstantsDates(Me%File%NumberOfInstants)
#endif            
        endif

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = Me%StartTime,                 &
                                                EndTime = Me%EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR190'
        
        call GetComputeTimeStep(Me%ObjTime, DT = Me%DT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR195'
        
        call GetData(Me%File%LonStagName,                                               &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'LONG_GRID',                                      &
                     default        = Me%File%DefaultNames%lon_stag,                    &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR200'
        
        call GetData(Me%File%LatStagName,                                               &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'LAT_GRID',                                       &
                     default        = Me%File%DefaultNames%lat_stag,                    &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR210'

        call GetData(Me%File%DepthStagName,                                             &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'DEPTH_GRID',                                     &
                     default        = Me%File%DefaultNames%depth_stag,                  &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR220'

        call GetData(Me%File%BathymName,                                                &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'BATHYM_GRID',                                    &
                     default        = Me%File%DefaultNames%bat,                         &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR230'         
        

        call GetData(Me%File%MaskName,                                                  &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'MASK_GRID',                                      &
                     default        = Me%File%DefaultNames%mask,                        &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR240'        
 
               
    end subroutine ConstructFile

   !---------------------------------------------------------------------------

    subroutine AllocatePropertyField(NewPropField, MaskDim)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                     :: NewPropField, AuxPropField        
        integer                                         :: MaskDim
        !External--------------------------------------------------------------
 
        !Local-----------------------------------------------------------------
 

        !Begin-----------------------------------------------------------------

        !Add to list 
        allocate(NewPropField)
        
        nullify (NewPropField%Next)
        
        nullify(NewPropField%PreviousField2D, NewPropField%NextField2D)
        nullify(NewPropField%PreviousField3D, NewPropField%NextField3D)
        
        NewPropField%SpaceDim = MaskDim

        !Add a new property to the reading list 
        if (.not.associated(Me%FirstPropField)) then
            Me%FirstPropField => NewPropField
        else
            AuxPropField => Me%FirstPropField
            do while(associated(AuxPropField%Next)) 
                AuxPropField => AuxPropField%Next
            enddo
            AuxPropField%Next => NewPropField
            nullify(AuxPropField)
        endif        
        

    end subroutine AllocatePropertyField
    
    !---------------------------------------------------------------------------    

    subroutine ConstructPropertyField(NewPropField)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                     :: NewPropField        

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB, i, j, k
        logical                                         :: FoundSecondInstant

        !Begin-----------------------------------------------------------------


i0:     if(NewPropField%SpaceDim == Dim2D)then

            ILB = Me%Size2D%ILB
            IUB = Me%Size2D%IUB
            JLB = Me%Size2D%JLB
            JUB = Me%Size2D%JUB

            allocate(NewPropField%PreviousField2D (ILB:IUB, JLB:JUB))
            allocate(NewPropField%NextField2D     (ILB:IUB, JLB:JUB))

            NewPropField%PreviousField2D(:,:) = FillValueReal
            NewPropField%NextField2D    (:,:) = FillValueReal
            
            if (.not.Associated(Me%Matrix2D)) then
                allocate(Me%Matrix2D(ILB:IUB, JLB:JUB))
                Me%Matrix2D(:,:) = FillValueReal
            endif
            
            if (.not.Associated(Me%Matrix3D)) then
                allocate(Me%Matrix3D(ILB:IUB, JLB:JUB, 0:2))
                Me%Matrix3D(:,:,:) = FillValueReal
            endif
            

        else i0

            ILB = Me%Size3D%ILB
            IUB = Me%Size3D%IUB
            JLB = Me%Size3D%JLB
            JUB = Me%Size3D%JUB
            KLB = Me%Size3D%KLB
            KUB = Me%Size3D%KUB

            allocate(NewPropField%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate(NewPropField%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))

            NewPropField%PreviousField3D(:,:,:) = FillValueReal
            NewPropField%NextField3D    (:,:,:) = FillValueReal

            if (.not.Associated(Me%Matrix3D)) then
                allocate(Me%Matrix3D(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Matrix3D(:,:,:) = FillValueReal
            endif

            if (.not.Associated(Me%Depth3D)) then
                allocate(Me%Depth3D(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Depth3D(:,:,:) = FillValueReal
            endif

        endif i0
        
it:     if (NewPropField%ChangeInTime) then

            FoundSecondInstant = .false.

            if (Me%Backtracking) then

                NewPropField%PreviousInstant  = Me%File%NumberOfInstants
                NewPropField%NextInstant      = NewPropField%PreviousInstant
                NewPropField%PreviousTime     = Me%File%InstantsDates(NewPropField%PreviousInstant) 
                    

                if(NewPropField%PreviousTime .lt. Me%EndTime)then
                    write(*,*)
                    !write(*,*)'Could not read solution from file'
                    write(*,*) 'Filename =', trim(Me%File%FileName)                    
                    write(*,*)'First file instant greater than current time'
                    write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                    !stop      'ConstructPropertyField - ModuleField4D - ERR170'
                    
                    FoundSecondInstant        = .true. 
                    NewPropField%NextInstant  = NewPropField%PreviousInstant - 1
                    NewPropField%NextTime     = Me%File%InstantsDates(NewPropField%NextInstant) 
                end if

                            
                if(Me%File%StartTime .gt. Me%StartTime)then
                    write(*,*)
                    !write(*,*)'Could not read solution from file'
                    write(*,*) 'Filename =', trim(Me%File%FileName)                    
                    write(*,*)'Last instant in file lower than simulation ending time'
                    write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                    !stop      'ConstructPropertyField - ModuleField4D - ERR180'
                end if
            
            else
                NewPropField%PreviousInstant  = 1
                NewPropField%NextInstant      = NewPropField%PreviousInstant
                NewPropField%PreviousTime     = Me%File%InstantsDates(NewPropField%PreviousInstant) 
                    

                if(NewPropField%PreviousTime .gt. Me%StartTime)then
                    write(*,*)
                    !write(*,*)'Could not read solution from file'
                    write(*,*) 'Filename =', trim(Me%File%FileName)                    
                    write(*,*)'First file instant greater than current time'
                    write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                    !stop      'ConstructPropertyField - ModuleField4D - ERR190'
                    
                    FoundSecondInstant        = .true. 
                    NewPropField%NextInstant  = NewPropField%PreviousInstant + 1
                    NewPropField%NextTime     = Me%File%InstantsDates(NewPropField%NextInstant) 
                    
                end if

                            
                if(Me%File%EndTime .lt. Me%EndTime)then
                    write(*,*)
                    !write(*,*)'Could not read solution from file'
                    write(*,*) 'Filename =', trim(Me%File%FileName)                    
                    write(*,*)'Last instant in file lower than simulation ending time'
                    write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                    !stop      'ConstructPropertyField - ModuleField4D - ERR200'
                end if
                
            endif

        
            !if number of instants greater than 1 then 
            !find first and second instants
    d2:     do while(.not. FoundSecondInstant)

               NewPropField%PreviousInstant  = NewPropField%NextInstant
                                
               if (Me%Backtracking) then

                    NewPropField%NextInstant      = NewPropField%NextInstant - 1

                    NewPropField%NextTime         = Me%File%InstantsDates(NewPropField%NextInstant)

                    if(NewPropField%PreviousTime .ge. Me%EndTime .and. NewPropField%NextTime .le. Me%EndTime) then
                        FoundSecondInstant  = .true.
                        exit
                    end if

                    NewPropField%PreviousTime   = NewPropField%NextTime

                    if(NewPropField%NextInstant .lt. 1) then
                        write(*,*)
                        !write(*,*)'Could not read solution from file'
                        write(*,*) 'Filename =', trim(Me%File%FileName)                    
                        write(*,*)'Could not find second instant in file'
                        write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                        !stop      'ConstructPropertyField - ModuleField4D - ERR210'
                    end if

               
               else
                    NewPropField%NextInstant      = NewPropField%NextInstant + 1

                    NewPropField%NextTime         = Me%File%InstantsDates(NewPropField%NextInstant)

                    if(NewPropField%PreviousTime .le. Me%StartTime .and. NewPropField%NextTime .ge. Me%StartTime) then
                        FoundSecondInstant  = .true.
                        exit
                    end if

                    NewPropField%PreviousTime            = NewPropField%NextTime

                    if(NewPropField%NextInstant .gt. Me%File%NumberOfInstants) then
                        write(*,*)
                        !write(*,*)'Could not read solution from file'
                        write(*,*) 'Filename =', trim(Me%File%FileName)                    
                        write(*,*)'Could not find second instant in file'
                        write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                        !stop      'ConstructPropertyField - ModuleField4D - ERR220'
                    end if

                endif
            end do d2

    i4:     if(NewPropField%SpaceDim == Dim2D)then

                ILB = Me%Size2D%ILB
                IUB = Me%Size2D%IUB
                JLB = Me%Size2D%JLB
                JUB = Me%Size2D%JUB

                allocate(NewPropField%PreviousField2D (ILB:IUB, JLB:JUB))
                allocate(NewPropField%NextField2D     (ILB:IUB, JLB:JUB))

                NewPropField%PreviousField2D(:,:) = FillValueReal
                NewPropField%NextField2D    (:,:) = FillValueReal

                call ReadValues2D(NewPropField, Previous = .true. )
                call ReadValues2D(NewPropField, Previous = .false.)
                
                !limit maximum values
                do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                
                    if (abs(NewPropField%PreviousField2D(i,j)) > abs(FillValueReal))          &
                        NewPropField%PreviousField2D(i,j) = FillValueReal
                    
                    if (abs(NewPropField%NextField2D    (i,j)) > abs(FillValueReal))          &
                        NewPropField%NextField2D    (i,j) = FillValueReal
                enddo
                enddo   
                                
                if (NewPropField%PreviousInstant /= NewPropField%NextInstant) then
                
                    if (NewPropField%ID%IsAngle) then                
                        call InterpolateAngle2DInTime (ActualTime       = Me%StartTime,                &
                                                       Size             = Me%WorkSize2D,               &
                                                       Time1            = NewPropField%PreviousTime,   &
                                                       Matrix1          = NewPropField%PreviousField2D,&
                                                       Time2            = NewPropField%NextTime,       &
                                                       Matrix2          = NewPropField%NextField2D,    &
                                                       MatrixOut        = Me%Matrix2D,                 &
                                                       PointsToFill2D   = Me%ExternalVar%WaterPoints2D)

                    else
                        !Interpolates the two matrixes in time
                        call InterpolateMatrix2DInTime(ActualTime       = Me%StartTime,                &
                                                       Size             = Me%WorkSize2D,               &
                                                       Time1            = NewPropField%PreviousTime,   &
                                                       Matrix1          = NewPropField%PreviousField2D,&
                                                       Time2            = NewPropField%NextTime,       &
                                                       Matrix2          = NewPropField%NextField2D,    &
                                                       MatrixOut        = Me%Matrix2D,                 &
                                                       PointsToFill2D   = Me%ExternalVar%WaterPoints2D)
                    endif                                                       
                else

                    Me%Matrix2D(:,:)  = NewPropField%PreviousField2D(:,:)

                endif

            else i4

                ILB = Me%Size3D%ILB
                IUB = Me%Size3D%IUB
                JLB = Me%Size3D%JLB
                JUB = Me%Size3D%JUB
                KLB = Me%Size3D%KLB
                KUB = Me%Size3D%KUB

                allocate(NewPropField%PreviousField3D (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate(NewPropField%NextField3D     (ILB:IUB, JLB:JUB, KLB:KUB))

                NewPropField%PreviousField3D(:,:,:) = FillValueReal
                NewPropField%NextField3D    (:,:,:) = FillValueReal

                call ReadValues3D(NewPropField, Previous = .true. )
                call ReadValues3D(NewPropField, Previous = .false.)

                !limit maximum values
                do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                
                    if (abs(NewPropField%PreviousField3D(i,j,k)) > abs(FillValueReal))        &
                        NewPropField%PreviousField3D(i,j,k) = FillValueReal
                    
                    if (abs(NewPropField%NextField3D    (i,j,k)) > abs(FillValueReal))        &
                        NewPropField%NextField3D    (i,j,k) = FillValueReal
                enddo
                enddo
                enddo                


                if (NewPropField%PreviousInstant /= NewPropField%NextInstant) then
                
                    if (NewPropField%ID%IsAngle) then                                

                        call InterpolateAngle3DInTime (ActualTime       = Me%StartTime,                 &
                                                       Size             = Me%WorkSize3D,                &
                                                       Time1            = NewPropField%PreviousTime,    &
                                                       Matrix1          = NewPropField%PreviousField3D, &
                                                       Time2            = NewPropField%NextTime,        &
                                                       Matrix2          = NewPropField%NextField3D,     &
                                                       MatrixOut        = Me%Matrix3D,                  &
                                                       PointsToFill3D   = Me%ExternalVar%WaterPoints3D)
                    
                    else

                        call InterpolateMatrix3DInTime(ActualTime       = Me%StartTime,                 &
                                                       Size             = Me%WorkSize3D,                &
                                                       Time1            = NewPropField%PreviousTime,    &
                                                       Matrix1          = NewPropField%PreviousField3D, &
                                                       Time2            = NewPropField%NextTime,        &
                                                       Matrix2          = NewPropField%NextField3D,     &
                                                       MatrixOut        = Me%Matrix3D,                  &
                                                       PointsToFill3D   = Me%ExternalVar%WaterPoints3D)
                    endif
                else

                    !Prev and next are equal (last instant?)
                    Me%Matrix3D(:,:,:)  = NewPropField%NextField3D(:,:,:)

                endif

            end if i4

        endif it

    end subroutine ConstructPropertyField

    !--------------------------------------------------------------------------

    subroutine ConstructPropertyFieldHarmonics(NewPropField)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                     :: NewPropField        

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB, JLB, JUB, KLB, KUB
        type (T_Time)                                   :: CurrentTime

        !Begin-----------------------------------------------------------------
        
        if (Me%Backtracking) then
            CurrentTime = Me%EndTime
        else
            CurrentTime = Me%StartTime
        endif            

i0:     if(NewPropField%SpaceDim == Dim2D)then

            ILB = Me%Size2D%ILB
            IUB = Me%Size2D%IUB
            JLB = Me%Size2D%JLB
            JUB = Me%Size2D%JUB
            

            allocate(NewPropField%PreviousField2D (ILB:IUB, JLB:JUB))
            allocate(NewPropField%NextField2D     (ILB:IUB, JLB:JUB))

            NewPropField%PreviousField2D(:,:) = FillValueReal
            NewPropField%NextField2D    (:,:) = FillValueReal            

            if (.not.Associated(Me%Matrix2D)) then
                allocate(Me%Matrix2D(ILB:IUB, JLB:JUB))
                Me%Matrix2D(:,:) = FillValueReal
            endif
            
            if (.not.Associated(Me%Matrix3D)) then
                allocate(Me%Matrix3D(ILB:IUB, JLB:JUB, 0:2))
                Me%Matrix3D(:,:,:) = FillValueReal
            endif            

            call ReadValues2DHarmonics(NewPropField)
            
            if (Me%Backtracking) then            
                stop 'ModuleField4D - ConstructPropertyFieldHarmonics -ERR10'
            endif
            
            NewPropField%PreviousTime = Me%StartTime
            NewPropField%NextTime     = Me%StartTime + NewPropField%Harmonics%DT
            
            if (.not.NewPropField%Harmonics%Extract) then
            
                call FromHarmonics2Field2D(NewPropField, CurrentTime = NewPropField%PreviousTime) 
                NewPropField%PreviousField2D(:,:) = Me%Matrix2D(:,:)
                
                call FromHarmonics2Field2D(NewPropField, CurrentTime = NewPropField%NextTime)                 
                NewPropField%NextField2D    (:,:) = Me%Matrix2D(:,:)
                !call ModifyInputHarmonics2D(NewPropField)
            endif

        else i0

            ILB = Me%Size3D%ILB
            IUB = Me%Size3D%IUB
            JLB = Me%Size3D%JLB
            JUB = Me%Size3D%JUB
            KLB = Me%Size3D%KLB
            KUB = Me%Size3D%KUB

            if (.not.Associated(Me%Matrix3D)) then
                allocate(Me%Matrix3D(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Matrix3D(:,:,:) = FillValueReal
            endif

            if (.not.Associated(Me%Depth3D)) then
                allocate(Me%Depth3D(ILB:IUB, JLB:JUB, KLB:KUB))
                Me%Depth3D(:,:,:) = FillValueReal
            endif
            
            call ReadValues3DHarmonics(NewPropField)
            
            if (.not.NewPropField%Harmonics%Extract) then            
                call FromHarmonics2Field3D(NewPropField, CurrentTime) 
            endif                
            
        endif i0
        
    end subroutine ConstructPropertyFieldHarmonics

    !--------------------------------------------------------------------------


    character(len=StringLength) function NCDFName(Number)

        !Arguments-------------------------------------------------------------
        integer                            :: Number

        !Local-----------------------------------------------------------------
        character(len = StringLength)      :: Name
        integer                            :: i

        !----------------------------------------------------------------------

        Name = GetPropertyName (Number)

        Select case (Name)
    
            case("Bathymetry")
                NCDFName        = "bathymetry"
            case("WaterPoints2D", "WaterPoints3D", "MappingPoints2D", "WaterPoints")
                NCDFName        = "mask"
            case("OpenPoints2D", "OpenPoints3D")
                NCDFName        = "mask"
            case("temperature")
                NCDFName        = "temperature"
            case("salinity")
                NCDFName        = "salinity"
            case("density")
                NCDFName        = "sea_water_density"
            case("oxygen")
                NCDFName        = "dissolved_oxygen"
            case("dissolved_oxygen_percent_saturation")
                NCDFName        = "dissolved_oxygen_percent_saturation"
            case("velocity_U")
                NCDFName        = "u"
            case("velocity_V")
                NCDFName        = "v"
            case("velocity_W")
                NCDFName        = "w"
            case("velocity_modulus")
                NCDFName        = "vm"
            case("water_level")
                NCDFName        = "ssh"
            case("wind_velocity_X")
                NCDFName        = "x_wind"
            case("wind_velocity_Y")
                NCDFName        = "y_wind"
            case("air_temperature")
                NCDFName        = "air_temperature"
            case("atmospheric_pressure")
                NCDFName        = "air_pressure"
            case("short_wave_solar_radiation_extinction")
                NCDFName        = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
            case("short_wave_solar_radiation")
                NCDFName        = "short_wave_solar_radiation"
            case("phytoplankton")
                NCDFName        = "phytoplankton"
            case("zooplankton")
                NCDFName        = "zooplankton"
            case("nitrate")
                NCDFName        = "nitrate"
            case("ammonia")
                NCDFName        = "ammonia"
            case default
                NCDFName        = trim(adjustl(Name))
                do i=1,len_trim(NCDFName)
                    if (NCDFName(i:i)==" ") NCDFName(i:i)="_"
                enddo
        end select

    end function NCDFName
   !----------------------------------------------------------------------------


    type(T_Time) function HDF5TimeInstant(Instant, HDF5ID)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        integer, optional                       :: HDF5ID
        

        !Local-----------------------------------------------------------------
!        type(T_Time)                            :: TimeInstant
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL, HDF5ID_

        !Begin-----------------------------------------------------------------
        
        allocate(TimeVector(6))
        
        if (present(HDF5ID)) then
            HDF5ID_ = HDF5ID
        else
            HDF5ID_ = Me%File%Obj
        endif
        
        call HDF5SetLimits  (HDF5ID_, 1, 6, STAT = STAT_CALL)        

        call HDF5ReadWindow (HDF5ID         = HDF5ID_,                                  &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleField4D - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

                                     
        deallocate(TimeVector)

    end function HDF5TimeInstant

    
    !--------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    subroutine BacktrackingTime

        !Arguments------------------------------------------------------------

        !Local----------------------------------------------------------------
        real                                            :: TotalTime, AuxPeriod

        !Begin----------------------------------------------------------------


        TotalTime = Me%EndTime        - Me%StartTime                  
        AuxPeriod = Me%CurrentTimeExt - Me%StartTime
        AuxPeriod = TotalTime         - AuxPeriod
        
        Me%CurrentTimeInt = Me%StartTime + AuxPeriod

    end subroutine BacktrackingTime
    
    !-------------------------------------------------------------------------    
    !-------------------------------------------------------------------------
    
    logical function ReadNewField(PropField,n)

        !Arguments------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer      , intent(OUT)                      :: n
        !Local----------------------------------------------------------------
        logical                                         :: ReadNewField_

        !Begin----------------------------------------------------------------

        ReadNewField_ = .false.
        if (Me%BackTracking) then  
            if (Me%CurrentTimeInt .le. PropField%NextTime) ReadNewField_ = .true.
        else            
            if (Me%CurrentTimeInt .ge. PropField%NextTime) ReadNewField_ = .true.
        endif
        
        ReadNewField = ReadNewField_

        if (ReadNewField_)then

            n = 0
            
            do 

                if (Me%BackTracking) then  
                    if (Me%CurrentTimeInt .gt. PropField%NextTime) exit
                else            
                    if (Me%CurrentTimeInt .lt. PropField%NextTime) exit
                endif            
                
                PropField%PreviousInstant  = PropField%NextInstant
                
                if (Me%BackTracking) then
                    if(PropField%NextInstant .gt. 1)then
                        PropField%NextInstant  = PropField%NextInstant - 1
                    else
                        exit
                    endif
                else
                    if(PropField%NextInstant .lt. Me%File%NumberOfInstants)then
                        PropField%NextInstant  = PropField%NextInstant + 1
                    else
                        exit
                    endif
                endif
                
                PropField%PreviousTime     = PropField%NextTime
                PropField%NextTime         = Me%File%InstantsDates(PropField%NextInstant)

                n = n + 1
                
                
            enddo
            
            if (Me%BackTracking) then
                if(Me%CurrentTimeInt .lt. PropField%NextTime)then
                    write(*,*)
                    write(*,*)'----------Backtracking mode-----------'
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ReadNewField - ModuleField4D - ERR20'
                end if
            else
                if(Me%CurrentTimeInt .gt. PropField%NextTime)then
                    write(*,*)
                    write(*,*)'Could not read solution from HDF5 file'
                    write(*,*)'Time instants inconsistency.'
                    stop      'ReadNewField - ModuleField4D - ERR30'
                end if
            endif    
    
        endif 
        
    end function ReadNewField

    !-------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------
    subroutine ReadValues2D (NewPropField, Previous)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
        logical                                 :: Previous
       
        !Local-----------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:), pointer           :: Field
        integer                                 :: STAT_CALL, Imax, Jmax, i, j, ILB, IUB, JLB, JUB
        integer                                 :: Obj, iaux

        !Begin-----------------------------------------------------------------
        
        if (NewPropField%Generic4D%InstantON) then
            Instant = NewPropField%Generic4D%Instant
        else
            if (Previous) then
                Instant =  NewPropField%PreviousInstant        
            else
                Instant =  NewPropField%NextInstant
            endif   
        endif
        
        if (Previous) then
            Field    => NewPropField%PreviousField2D
        else
            Field    => NewPropField%NextField2D        
        endif
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        
        if      (Me%File%Form == HDF5_  ) then

            if (Me%File%FileListON) then
            
                Obj  = Me%File%ObjListInst(Instant)
                iaux = Me%File%ListInst   (Instant)
                 
            else

                Obj  = Me%File%Obj
                iaux = Instant
            
            endif            
        
            call GetHDF5ArrayDimensions(Obj, trim(NewPropField%VGroupPath),         &
                              trim(NewPropField%FieldName), OutputNumber = iaux,    &
                              Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR10'                                   
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, STAT = STAT_CALL)
#endif
        endif            

        if ((Imax < IUB - ILB + 1) .or.                                                &
            (Jmax < JUB - JLB + 1)) then
            
            write (*,*) trim(NewPropField%FieldName)
            write (*,*) 'miss match between the input file and model domain'
            stop 'ReadValues2D - ModuleField4D - ERR20'                                   

        endif
          


        if      (Me%File%Form == HDF5_  ) then

            call HDF5SetLimits  (Obj, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR30'
            
            call HDF5ReadWindow(HDF5ID        = Obj,                                    &
                                GroupName     = trim(NewPropField%VGroupPath),          &
                                Name          = trim(NewPropField%FieldName),           &
                                Array2D       = Field,                                  &
                                OutputNumber  = iaux,                                   &        
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR40'
            
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                Array2D         = Field,                                &
                                Name            = trim(NewPropField%FieldName),         &
                                nInstant        = Instant,                              &
                                ILB             = ILB,                                  &
                                IUB             = IUB,                                  &
                                JLB             = JLB,                                  &
                                JUB             = JUB,                                  &
                                STAT            = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR50'            
#endif            
        endif


        if(NewPropField%HasMultiplyingFactor)then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                Field(i,j) = Field(i,j) * NewPropField%MultiplyingFactor
            enddo
            enddo
        end if

        if(NewPropField%HasAddingFactor)then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                Field(i,j) = Field(i,j) + NewPropField%AddingFactor
            enddo
            enddo
        end if
        
        if(NewPropField%MinValueON)then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                if (Field(i,j) < NewPropField%MinValue) then
                    Field(i,j) = NewPropField%MinValue
                endif
            enddo
            enddo
        end if        

        if(NewPropField%MaxValueON)then
            do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                if (Field(i,j) > NewPropField%MaxValue) then
                    Field(i,j) = NewPropField%MaxValue
                endif
            enddo
            enddo
        end if        

        
        nullify(Field)

    end subroutine ReadValues2D
    
    
    !--------------------------------------------------------------------------

    
    subroutine ReadValues3D (NewPropField, Previous)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
        logical                                 :: Previous
       
        !Local-----------------------------------------------------------------
        integer                                 :: Instant
        real,    dimension(:,:,:), pointer      :: Field, Aux3D, FieldAux, SZZ
        integer, dimension(:,:,:), pointer      :: WaterPoints3D     
        real                                    :: HT   
        integer                                 :: Imax, Jmax, Kmax
        integer                                 :: STAT_CALL, i, j, k
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: SILB, SIUB, SJLB, SJUB, SKLB, SKUB        
        integer                                 :: Obj, iaux, kbottom
        integer                                 :: nItems, iVert

        !Begin-----------------------------------------------------------------
        
        if (NewPropField%Generic4D%InstantON) then
            Instant = NewPropField%Generic4D%Instant
        else
            if (Previous) then
                Instant =  NewPropField%PreviousInstant        
            else
                Instant =  NewPropField%NextInstant
            endif   
        endif
        
        if (Previous) then
            Field    => NewPropField%PreviousField3D
        else
            Field    => NewPropField%NextField3D        
        endif
                    
        ILB = Me%WorkSize3D%ILB
        IUB = Me%WorkSize3D%IUB
        JLB = Me%WorkSize3D%JLB
        JUB = Me%WorkSize3D%JUB
        
        SILB = Me%Size3D%ILB
        SIUB = Me%Size3D%IUB
        SJLB = Me%Size3D%JLB
        SJUB = Me%Size3D%JUB

        if (NewPropField%From2Dto3D) then
            KLB = 1
            KUB = 1
            SKLB = Me%Size3D%KLB
            SKUB = Me%Size3D%KUB
            
        else
 
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB
            SKLB = Me%Size3D%KLB
            SKUB = Me%Size3D%KUB
            
        endif             
        
         
        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                            &
                              WaterPoints3D     = WaterPoints3D,                        &
                              STAT              = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'ReadValues3D - ModuleField4D - ERR100'          
         
        if (NewPropField%From2Dto3D .or. NewPropField%From3Dto2D) then
            allocate(Aux3D(SILB:SIUB, SJLB:SJUB, SKLB:SKUB))   
            Aux3D(:,:,:) = 0.      
        endif            
            

        if (NewPropField%From2Dto3D) then
            FieldAux => Aux3D
        else
            FieldAux => Field
        endif
        

        if      (Me%File%Form == HDF5_  ) then

            if (Me%File%FileListON) then
            
                Obj  = Me%File%ObjListInst(Instant)
                iaux = Me%File%ListInst   (Instant)
                 
            else

                Obj  = Me%File%Obj
                iaux = Instant
            
            endif                         
        
            call GetHDF5ArrayDimensions(Obj, trim(NewPropField%VGroupPath),             &
                              trim(NewPropField%FieldName), OutputNumber = iaux,        &
                              Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR10'                                   
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, KUB = Kmax, STAT = STAT_CALL)
#endif
        endif   
        
        

        if ((Imax < IUB - ILB + 1) .or.                                                &
            (Jmax < JUB - JLB + 1) .or.                                                &
            (Kmax < KUB - KLB + 1)) then
            
            write (*,*) trim(NewPropField%FieldName)
            write (*,*) 'miss match between the input file and model domain'
            stop 'ReadValues3D - ModuleField4D - ERR30'                                   

        endif        
        
        if      (.not. Me%File%Form == HDF5_  .and. NewPropField%From3Dto2D) then
            stop 'ReadValues3D - ModuleField4D - ERR300'                                   
        endif
        
        
        if      (Me%File%Form == HDF5_  ) then
        
            call HDF5SetLimits  (Obj, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR40'
            
            
            call HDF5ReadWindow(HDF5ID        = Obj,                                    &
                                GroupName     = trim(NewPropField%VGroupPath),          &
                                Name          = trim(NewPropField%FieldName),           &
                                Array3D       = FieldAux,                               &
                                OutputNumber  = iaux,                                   &
                                STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR50'
            

            
            if (Me%File%ReadSZZ .and.  Me%File%SZZLast < iaux) then
            
                allocate(SZZ(ILB-1:IUB+1, JLB-1:JUB+1, KLB-1:KUB+1))
                SZZ(:,:,:) = 0.
                
                call HDF5SetLimits  (Obj, ILB, IUB, JLB, JUB, KLB-1, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR60'    
                
                call GetHDF5GroupNumberOfItems (HDF5ID        = Obj,                    &
                                                GroupName     = "/Grid/VerticalZ",      &
                                                nItems        = nItems,                 &
                                                STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ReadValues3D - ModuleField4D - ERR70'            
                endif
                
                !Special case where there is only one  "/Grid/VerticalZ/Vertical_00001"
                !Assumed stationay condition for the grid
                if (nItems == 1) then
                    iVert = 1
                else
                    iVert = iaux
                endif
                
                call HDF5ReadWindow(HDF5ID        = Obj,                                &
                                    GroupName     = "/Grid/VerticalZ",                  &
                                    Name          = trim(Me%File%DataSetVert),          &
                                    Array3D       = SZZ,                                &
                                    OutputNumber  = ivert,                              &
                                    OffSet3       = 0,                                  &                                    
                                    STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'ReadValues3D - ModuleField4D - ERR70'            
                endif
                
                Me%File%SZZLast = iaux
                
                                                        
                call ComputeVerticalGeometry(GeometryID         = Me%ObjGeometry,       &
                                             WaterPoints3D      = WaterPoints3D,        &
                                             SZZ                = SZZ,                  &
                                             STAT               = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR80'
                
                deallocate(SZZ)
                
            endif            
            

#ifndef _NO_NETCDF                                                   
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                Array3D         = FieldAux,                             &
                                Name            = trim(NewPropField%FieldName),         &
                                nInstant        = Instant,                              &
                                ILB             = ILB,                                  &
                                IUB             = IUB,                                  &
                                JLB             = JLB,                                  &
                                JUB             = JUB,                                  &
                                KLB             = KLB,                                  &
                                KUB             = KUB,                                  &
                                STAT            = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR90'            
#endif            
        endif

        
        if (NewPropField%From2Dto3D) then    
           
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j =               JLB,               JUB
            do i =               ILB,               IUB
                Field(i,j,k) = FieldAux(i,j,1)
            enddo
            enddo
            enddo
         endif        
         
        if (NewPropField%From3Dto2D) then   
        
            call GetGeometryDistances(  GeometryID      = Me%ObjGeometry,               &
                                                SZZ     = SZZ,                          &
                                        STAT            = STAT_CALL)                                     
            if (STAT_CALL /= SUCCESS_) stop 'ReadValues3D - ModuleValida4D - ERR110'
            
           
            do j = JLB, JUB
            do i = ILB, IUB

                kbottom = -99
                do k = KLB, KUB
                    if (WaterPoints3D (i,j,k)== 1) then
                        if (kbottom< 0) then
                            kbottom = k
                            HT      = (SZZ(i,j,kbottom-1) - SZZ(i,j,KUB))
                        endif                            
                        if (HT > 0) then
                            Aux3D(i,j,1) = Aux3D(i,j,1) + FieldAux(i,j,k) * (SZZ(i,j,k-1)-SZZ(i,j,k)) / HT
                        endif                            
                    endif                        
                enddo
            
            enddo
            enddo
            
            FieldAux(:,:,KUB) = Aux3D(:,:,1)
            
            call UnGetGeometry       (  GeometryID      = Me%ObjGeometry,               &
                                        Array           = SZZ,                          &
                                        STAT            = STAT_CALL)                                     
            if (STAT_CALL /= SUCCESS_) stop 'ReadValues3D - ModuleValida4D - ERR120'        
            

            
            
         endif               
         
         nullify(FieldAux)           

        if(NewPropField%HasMultiplyingFactor)then
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                Field(i,j,k) = Field(i,j,k) * NewPropField%MultiplyingFactor
            enddo
            enddo
            enddo
            
        end if

        if(NewPropField%HasAddingFactor)then
            
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                Field(i,j,k) = Field(i,j,k) + NewPropField%AddingFactor
            enddo
            enddo
            enddo

        end if
        
        if(NewPropField%MinValueON)then
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                if (Field(i,j,k) < NewPropField%MinValue) then
                    Field(i,j,k) = NewPropField%MinValue
                endif
            enddo
            enddo
            enddo
        end if        

        if(NewPropField%MaxValueON)then
            do k = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
                if (Field(i,j,k) > NewPropField%MaxValue) then
                    Field(i,j,k) = NewPropField%MaxValue
                endif
            enddo
            enddo
            enddo
        end if        

        if (NewPropField%From2Dto3D .or. NewPropField%From3Dto2D) then
            deallocate(Aux3D)   
        endif       


        call UnGetMap(Map_ID            = Me%ObjMap,                                &
                      Array             = WaterPoints3D,                            &
                      STAT              = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'ReadValues3D - ModuleField4D - ERR130'         

        nullify(Field)
    end subroutine ReadValues3D


    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------
    subroutine ReadValues2DHarmonics (NewPropField)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
       
        !Local-----------------------------------------------------------------
        real, dimension(:,:  ), pointer         :: Field
        integer                                 :: STAT_CALL, Imax, Jmax, ILB, IUB, JLB, JUB, NW, N
        integer                                 :: SILB, SIUB, SJLB, SJUB
        integer                                 :: i, j
        character(len=StringLength)             :: GroupName, FieldName

        !Begin-----------------------------------------------------------------
        
        ILB  = Me%WorkSize2D%ILB
        IUB  = Me%WorkSize2D%IUB
        JLB  = Me%WorkSize2D%JLB
        JUB  = Me%WorkSize2D%JUB
        
        SILB = Me%Size2D%ILB
        SIUB = Me%Size2D%IUB
        SJLB = Me%Size2D%JLB
        SJUB = Me%Size2D%JUB
        
        NW  = NewPropField%Harmonics%Number
        
        allocate(Field(SILB:SIUB,SJLB:SJUB))
        
d1:     do N =1, NW       

            GroupName = trim(NewPropField%VGroupPath)//"/"//trim(NewPropField%Harmonics%WaveGroupName(N))
            
            FieldName = char_amplitude_ 
        
            if      (Me%File%Form == HDF5_  ) then
                call GetHDF5ArrayDimensions(Me%File%Obj, GroupName, FieldName,          &
                                            Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR10'                                   
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, STAT = STAT_CALL)
#endif
            endif            

            if ((Imax < IUB - ILB + 1) .or.                                                &
                (Jmax < JUB - JLB + 1)) then
                
                write (*,*) "GroupName =",trim(GroupName)
                write (*,*) "FieldName =",trim(FieldName)
                write (*,*) 'miss match between the input file and model domain'
                stop 'ReadValues2DHarmonics - ModuleField4D - ERR20'                                   

            endif
          
            if      (Me%File%Form == HDF5_  ) then

                call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR30'
                
                     
                call HDF5ReadWindow(Me%File%Obj, GroupName, FieldName,                  &
                                    Array2D = Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR40'
                
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                FieldName = trim(NewPropField%Harmonics%WaveGroupName(N))//"/"//char_amplitude_
            
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array2D         = Field,                            &
                                    Name            = FieldName,                        &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR50'            
#endif            
            endif           
            
            NewPropField%Harmonics%Amplitude2D(:,:,N) = Field(:,:)
            
        enddo d1
        
d2:     do N =1, NW       

            GroupName = trim(NewPropField%VGroupPath)//"/"//trim(NewPropField%Harmonics%WaveGroupName(N))
            
            FieldName = char_phase_ 
        
            if      (Me%File%Form == HDF5_  ) then
                call GetHDF5ArrayDimensions(Me%File%Obj, GroupName, FieldName,          &
                                            Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR60'                                   
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, STAT = STAT_CALL)
#endif
            endif            

            if ((Imax < IUB - ILB + 1) .or.                                                &
                (Jmax < JUB - JLB + 1)) then
                
                write (*,*) "GroupName =",trim(GroupName)
                write (*,*) "FieldName =",trim(FieldName)
                write (*,*) 'miss match between the input file and model domain'
                stop 'ReadValues2DHarmonics - ModuleField4D - ERR70'                                   

            endif
          
            if      (Me%File%Form == HDF5_  ) then

                call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR80'
                
                     
                call HDF5ReadWindow(Me%File%Obj, GroupName, FieldName,                  &
                                    Array2D = Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR90'
                
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                FieldName = trim(NewPropField%Harmonics%WaveGroupName(N))//"/"//char_phase_
            
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array2D         = Field,                            &
                                    Name            = FieldName,                        &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR100'            
#endif            
            endif         
            
            do j = JLB, JUB
            do i = ILB, IUB
                NewPropField%Harmonics%Phase2D(i,j,N) = Field(i,j) / 360.
                call AmpPhase_To_Complex(Amplitude = 1.,                                       &  
                                         Phase     = NewPropField%Harmonics%Phase2D    (i,j,N),&
                                         Sreal     = NewPropField%Harmonics%Phase2DReal(i,j,N),& 
                                         Simag     = NewPropField%Harmonics%Phase2DImag(i,j,N))
            enddo
            enddo
            
        enddo d2        
        
        GroupName = trim(NewPropField%VGroupPath)
        
        FieldName = NewPropField%Harmonics%FieldNameDim 
    
        if      (Me%File%Form == HDF5_  ) then
            call GetHDF5ArrayDimensions(Me%File%Obj, GroupName, FieldName,          &
                                        Imax = Imax, Jmax = Jmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR60'                                   
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, STAT = STAT_CALL)
#endif
        endif            

        if ((Imax < IUB - ILB + 1) .or.                                                &
            (Jmax < JUB - JLB + 1)) then
            
            write (*,*) "GroupName =",trim(GroupName)
            write (*,*) "FieldName =",trim(FieldName)
            write (*,*) 'miss match between the input file and model domain'
            stop 'ReadValues2DHarmonics - ModuleField4D - ERR70'                                   

        endif
            
        if (trim(FieldName) == trim(char_residual_)) then
      
            if      (Me%File%Form == HDF5_  ) then

                call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR80'
            
                 
                call HDF5ReadWindow(Me%File%Obj, GroupName, FieldName,                  &
                                    Array2D = Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR90'
            
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
        
                FieldName = NewPropField%Harmonics%FieldNameDim
        
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array2D         = Field,                            &
                                    Name            = FieldName,                        &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues2DHarmonics - ModuleField4D - ERR100'            
#endif            
            endif           
        
            NewPropField%Harmonics%Residual2D(:,:) = Field(:,:) 
            
        else
            
            NewPropField%Harmonics%Residual2D(:,:) = 0.
        
        endif
        

        deallocate(Field)

    end subroutine ReadValues2DHarmonics
    
    
    !--------------------------------------------------------------------------

    subroutine FromHarmonics2Field2D(NewPropField, CurrentTime)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
        type (T_Time), optional                 :: CurrentTime
       
        !Local-----------------------------------------------------------------
        real, dimension(:  ),   pointer         :: Amplitude, Phase
        real, dimension(:,:),   pointer         :: Field
        real, dimension(:,:),   pointer         :: CoordX, CoordY        
        integer                                 :: STAT_CALL, ILB, IUB, JLB, JUB, NW, i, j, n
        real                                    :: T1, T2, T3, StateDT, DT_Run, Coef

        !Begin-----------------------------------------------------------------
       
        Field => Me%Matrix2D

        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        NW  = NewPropField%Harmonics%Number   
        
if1:    if (NewPropField%Harmonics%Extract) then

            do n = 1, NW
                if (trim(NewPropField%Harmonics%WaveName(n)) == trim(NewPropField%Harmonics%ExtractWave)) then
                    exit
                endif
            enddo
            
            if (NewPropField%Harmonics%ExtractAmp) then
                Field(ILB:IUB,JLB:JUB) = NewPropField%Harmonics%Amplitude2D(ILB:IUB,JLB:JUB,n)
            else
                if (NewPropField%Harmonics%ExtractPhaseReal) then
                    Field(ILB:IUB,JLB:JUB) = NewPropField%Harmonics%Phase2DReal(ILB:IUB,JLB:JUB,n)
                else
                    Field(ILB:IUB,JLB:JUB) = NewPropField%Harmonics%Phase2DImag(ILB:IUB,JLB:JUB,n)
                endif                     
            endif                
        
        else   if1
    
            call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR10'
    
            
            do i = ILB,IUB
            do j = JLB,JUB
                
                if (CoordX(i, j) > NewPropField%Harmonics%NoReadArea(1) .and.           &
                    CoordX(i, j) < NewPropField%Harmonics%NoReadArea(2) .and.           &
                    CoordY(i, j) > NewPropField%Harmonics%NoReadArea(3) .and.           &
                    CoordY(i, j) < NewPropField%Harmonics%NoReadArea(4)) then
                    Field(i, j) = NewPropField%DefaultValue
                    cycle
                endif
                Amplitude   => NewPropField%Harmonics%Amplitude2D(i, j, :)
                Phase       => NewPropField%Harmonics%Phase2D    (i, j, :)   

#ifndef _NOT_IEEE_ARITHMETIC
                do n=1,NW
                    if (ieee_is_nan(Amplitude(n))) Amplitude(n) = FillValueReal
                    if (ieee_is_nan(Phase    (n))) Phase    (n) = FillValueReal
                enddo                 
#endif
                if (sum(Amplitude(1:NW))>0.) then             
                
                    
                    call Task2000Level(WaterLevel       = T2,                                   &
                                        TimeReference    = NewPropField%Harmonics%TimeReference,&
                                        NWaves           = NewPropField%Harmonics%Number,       &
                                        WaveAmplitude    = Amplitude,                           &
                                        WavePhase        = Phase,                               &
                                        WaveName         = NewPropField%Harmonics%WaveName,     & 
                                        time_            = CurrentTime,                         &
                                        STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR20'

                    Field(i, j) = T2 + NewPropField%Harmonics%Residual2D(i, j)          
                    
                    if (NewPropField%Harmonics%TideStateON) then
                    
                        StateDT = NewPropField%Harmonics%TideStateDT

                        call Task2000Level(WaterLevel       = T1,                                   &
                                            TimeReference    = NewPropField%Harmonics%TimeReference,&
                                            NWaves           = NewPropField%Harmonics%Number,       &
                                            WaveAmplitude    = Amplitude,                           &
                                            WavePhase        = Phase,                               &
                                            WaveName         = NewPropField%Harmonics%WaveName,     & 
                                            time_            = CurrentTime - StateDT,               &
                                            STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR30'
                        
             
                        call Task2000Level(WaterLevel       = T3,                                   &
                                            TimeReference    = NewPropField%Harmonics%TimeReference,&
                                            NWaves           = NewPropField%Harmonics%Number,       &
                                            WaveAmplitude    = Amplitude,                           &
                                            WavePhase        = Phase,                               &
                                            WaveName         = NewPropField%Harmonics%WaveName,     & 
                                            time_            = CurrentTime + StateDT,               &
                                            STAT             = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR40'
                    
                        
                        !low tide 
                        if     (T3 >= T2 .and. T2 <= T1) then
                            Field(i,j) = 1
                        !flood                            
                        elseif (T3 >= T2 .and. T2 >= T1) then        
                            Field(i,j) = 2
                        !high tide 
                        elseif (T3 <= T2 .and. T2 >= T1) then
                            Field(i,j) = 3
                        !ebb tide 
                        elseif (T3 <= T2 .and. T2 <= T1) then
                            Field(i,j) = 4
                        endif                            
                    endif
                    
                else                
                    
                    Field(i, j) = FillValueReal
                    
                endif                
                
                nullify(Amplitude)
                nullify(Phase    )            

            enddo
            enddo
            
            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR50'

            call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR60'
            
            if(NewPropField%HasMultiplyingFactor)then
                do j = JLB, JUB
                do i = ILB, IUB
                    Field(i,j) = Field(i,j) * NewPropField%MultiplyingFactor
                enddo
                enddo
            end if

            if(NewPropField%HasAddingFactor)then
                do j = JLB, JUB
                do i = ILB, IUB
                    Field(i,j) = Field(i,j) + NewPropField%AddingFactor
                enddo
                enddo
            end if
        
            if(NewPropField%MinValueON)then
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Field(i,j) < NewPropField%MinValue) then
                        Field(i,j) = NewPropField%MinValue
                    endif
                enddo
                enddo
            end if        

            if(NewPropField%MaxValueON)then
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Field(i,j) > NewPropField%MaxValue) then
                        Field(i,j) = NewPropField%MaxValue
                    endif
                enddo
                enddo
            end if        
        
        
    iso:    if (NewPropField%Harmonics%SlowStartON) then
        
                DT_Run = CurrentTime -  Me%StartTime

                if (DT_Run < NewPropField%Harmonics%SlowStartPeriod) then

                    Coef = (DT_Run / NewPropField%Harmonics%SlowStartPeriod) ** 2. 

                    do j = JLB, JUB
                    do i = ILB, IUB
                        Field(i,j) = Field(i,j) * Coef
                    enddo
                    enddo

                endif        
                
            endif  iso        
            

            
        endif if1                        
    

        if (associated(Me%Matrix3D)) then
            Me%Matrix3D(ILB:IUB,JLB:JUB,1) = Field(ILB:IUB,JLB:JUB)
        endif            
        
        nullify(Field)
    
    end subroutine FromHarmonics2Field2D

    !--------------------------------------------------------------------------    
 
     
    !--------------------------------------------------------------------------
    subroutine ReadValues3DHarmonics (NewPropField)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
       
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:)  , pointer       :: Field
        integer                                 :: STAT_CALL, Imax, Jmax, Kmax
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB, NW, N
        integer                                 :: SILB, SIUB, SJLB, SJUB, SKLB, SKUB
        integer                                 :: i, j, k
        character(len=StringLength)             :: GroupName, FieldName

        !Begin-----------------------------------------------------------------
        
        ILB  = Me%WorkSize3D%ILB
        IUB  = Me%WorkSize3D%IUB
        JLB  = Me%WorkSize3D%JLB
        JUB  = Me%WorkSize3D%JUB

        SILB = Me%Size3D%ILB
        SIUB = Me%Size3D%IUB
        SJLB = Me%Size3D%JLB
        SJUB = Me%Size3D%JUB

        if (NewPropField%From2Dto3D) then
            KLB = 1
            KUB = 1
            SKLB = Me%Size3D%KLB
            SKUB = Me%Size3D%KUB
            
        else
 
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB
            SKLB = Me%Size3D%KLB
            SKUB = Me%Size3D%KUB
            
        endif        
        
        
        NW  = NewPropField%Harmonics%Number
        
        allocate(Field(SILB:SIUB,SJLB:SJUB,SKLB:SKUB))
        
d1:     do N =1, NW       

            GroupName = trim(NewPropField%VGroupPath)//"/"//trim(NewPropField%Harmonics%WaveGroupName(N))
            
            FieldName = char_amplitude_ 
        
            if      (Me%File%Form == HDF5_  ) then
                call GetHDF5ArrayDimensions(Me%File%Obj, GroupName, FieldName,          &
                                            Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR10'                                   
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, KUB = Kmax, STAT = STAT_CALL)
#endif
            endif            

            if ((Imax < IUB - ILB + 1) .or.                                             &
                (Jmax < JUB - JLB + 1) .or.                                             &
                (Kmax < KUB - KLB + 1)) then
                
                write (*,*) "GroupName =",trim(GroupName)
                write (*,*) "FieldName =",trim(FieldName)
                write (*,*) 'miss match between the input file and model domain'
                stop 'ReadValues3DHarmonics - ModuleField4D - ERR20'                                   

            endif
          
            if      (Me%File%Form == HDF5_  ) then

                call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR30'
                
                     
                call HDF5ReadWindow(Me%File%Obj, GroupName, FieldName,                  &
                                    Array3D = Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR40'
                
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                FieldName = trim(NewPropField%Harmonics%WaveGroupName(N))//"/"//char_amplitude_
            
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array3D         = Field,                            &
                                    Name            = FieldName,                        &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    KLB             = KLB,                              &
                                    KUB             = KUB,                              &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR50'            
#endif            
            endif           
            
            
            
            if (NewPropField%From2Dto3D) then
                
                do k=Me%Size3D%KLB, Me%Size3D%KUB
                    NewPropField%Harmonics%Amplitude3D(:,:,k,N) = Field(:,:,1)
                enddo
                
            else
                NewPropField%Harmonics%Amplitude3D(:,:,:,N) = Field(:,:,:)                
            endif        
                        
            
        enddo d1
        
d2:     do N =1, NW       

            GroupName = trim(NewPropField%VGroupPath)//"/"//trim(NewPropField%Harmonics%WaveGroupName(N))
            
            FieldName = char_phase_ 
        
            if      (Me%File%Form == HDF5_  ) then
                call GetHDF5ArrayDimensions(Me%File%Obj, GroupName, FieldName,          &
                                            Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR60'                                   
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, KUB = Kmax, STAT = STAT_CALL)
#endif
            endif            

            if ((Imax < IUB - ILB + 1) .or.                                             &
                (Jmax < JUB - JLB + 1) .or.                                             &
                (Kmax < KUB - KLB + 1)) then
                
                write (*,*) "GroupName =",trim(GroupName)
                write (*,*) "FieldName =",trim(FieldName)
                write (*,*) 'miss match between the input file and model domain'
                stop 'ReadValues3DHarmonics - ModuleField4D - ERR70'                                   

            endif
          
            if      (Me%File%Form == HDF5_  ) then

                call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR80'
                
                     
                call HDF5ReadWindow(Me%File%Obj, GroupName, FieldName,                  &
                                    Array3D = Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR90'
                
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
            
                FieldName = trim(NewPropField%Harmonics%WaveGroupName(N))//"/"//char_phase_
            
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array3D         = Field,                            &
                                    Name            = FieldName,                        &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    KLB             = KLB,                              &
                                    KUB             = KUB,                              &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR100'            
#endif            
            endif           
            
           
            if (NewPropField%From2Dto3D) then
                
                do k=Me%Size3D%KLB, Me%Size3D%KUB
                    NewPropField%Harmonics%Phase3D(:,:,k,N) = Field(:,:,1) / 360.
                enddo
                
            else
                NewPropField%Harmonics%Phase3D(:,:,:,N) = Field(:,:,:) / 360.
            endif     
            
            do k = KLB, KUB            
            do j = JLB, JUB
            do i = ILB, IUB
                call AmpPhase_To_Complex(Amplitude = 1.,                                         &  
                                         Phase     = NewPropField%Harmonics%Phase3D    (i,j,k,N),&
                                         Sreal     = NewPropField%Harmonics%Phase3DReal(i,j,k,N),& 
                                         Simag     = NewPropField%Harmonics%Phase3DImag(i,j,k,N))
            enddo
            enddo
            enddo            
            
        enddo d2        
        
        GroupName = trim(NewPropField%VGroupPath)
        
        FieldName = NewPropField%Harmonics%FieldNameDim 
        
        if (trim(FieldName) == trim(char_residual_)) then        
    
            if      (Me%File%Form == HDF5_  ) then
                call GetHDF5ArrayDimensions(Me%File%Obj, GroupName, FieldName,          &
                                            Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR60'                                   
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
        
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, JUB = Jmax, IUB = Imax, KUB = Kmax, STAT = STAT_CALL)
#endif
            endif            

            if ((Imax < IUB - ILB + 1) .or.                                             &
                (Jmax < JUB - JLB + 1) .or.                                             &
                (Kmax < KUB - KLB + 1)) then
            
                write (*,*) "GroupName =",trim(GroupName)
                write (*,*) "FieldName =",trim(FieldName)
                write (*,*) 'miss match between the input file and model domain'
                stop 'ReadValues3DHarmonics - ModuleField4D - ERR70'                                   

            endif
      
            if      (Me%File%Form == HDF5_  ) then

                call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR80'
            
                 
                call HDF5ReadWindow(Me%File%Obj, GroupName, FieldName,                  &
                                    Array3D = Field, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR90'
            
#ifndef _NO_NETCDF
            else if (Me%File%Form == NetCDF_) then
        
                FieldName = NewPropField%Harmonics%FieldNameDim
        
                call NETCDFReadData(NCDFID          = Me%File%Obj,                      &
                                    Array3D         = Field,                            &
                                    Name            = FieldName,                        &
                                    ILB             = ILB,                              &
                                    IUB             = IUB,                              &
                                    JLB             = JLB,                              &
                                    JUB             = JUB,                              &
                                    KLB             = KLB,                              &
                                    KUB             = KUB,                              &
                                    STAT            = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_)stop 'ReadValues3DHarmonics - ModuleField4D - ERR100'            
#endif            
            endif           
        
       
            if (NewPropField%From2Dto3D) then
            
                do k=Me%Size3D%KLB, Me%Size3D%KUB
                    NewPropField%Harmonics%Residual3D(:,:,k) = Field(:,:,1)
                enddo
            
            else
                NewPropField%Harmonics%Residual3D(:,:,:)   = Field(:,:,:)
            endif     
        
        else
            
            NewPropField%Harmonics%Residual3D(:,:,:) = 0.
            
        endif
            
        deallocate(Field)

    end subroutine ReadValues3DHarmonics
    
    
    !--------------------------------------------------------------------------

    subroutine FromHarmonics2Field3D(NewPropField, CurrentTime)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
        type (T_Time)                           :: CurrentTime
       
        !Local-----------------------------------------------------------------
        real, dimension(:    ),   pointer       :: Amplitude, Phase
        real, dimension(:,:  )  , pointer       :: CoordX, CoordY
        real, dimension(:,:,:)  , pointer       :: Field
        integer                                 :: STAT_CALL, ILB, IUB, JLB, JUB, KLB, KUB, NW, i, j, k, n

        !Begin-----------------------------------------------------------------
       
        Field       => Me%Matrix3D

        ILB = Me%WorkSize3D%ILB
        IUB = Me%WorkSize3D%IUB
        JLB = Me%WorkSize3D%JLB
        JUB = Me%WorkSize3D%JUB
        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB        
        NW  = NewPropField%Harmonics%Number   
        
        call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field3D - ModuleField4D - ERR10'        
        

        do j = JLB,JUB
        do i = ILB,IUB
            
            if (CoordX(i, j) > NewPropField%Harmonics%NoReadArea(1) .and.           &
                CoordX(i, j) < NewPropField%Harmonics%NoReadArea(2) .and.           &
                CoordY(i, j) > NewPropField%Harmonics%NoReadArea(3) .and.           &
                CoordY(i, j) < NewPropField%Harmonics%NoReadArea(4)) then
                Field(i, j, :) = NewPropField%DefaultValue
                cycle
            endif            
                    
            do k = KLB,KUB                            

                Amplitude   => NewPropField%Harmonics%Amplitude3D(i, j, k, :)
                Phase       => NewPropField%Harmonics%Phase3D    (i, j, k, :)    
            

#ifndef _NOT_IEEE_ARITHMETIC
                do n=1,NW
                    if (ieee_is_nan(Amplitude(n))) Amplitude(n) = FillValueReal
                    if (ieee_is_nan(Phase    (n))) Phase    (n) = FillValueReal
                enddo                 
#endif                
                if (sum(Amplitude(1:NW))>0.) then             
        
                    call Task2000Level(WaterLevel       = Field(i, j, k),                           &
                                       TimeReference    = NewPropField%Harmonics%TimeReference,     &
                                       NWaves           = NewPropField%Harmonics%Number,            &
                                       WaveAmplitude    = Amplitude,                                &
                                       WavePhase        = Phase,                                    &
                                       WaveName         = NewPropField%Harmonics%WaveName,          & 
                                       time_            = CurrentTime,                              &
                                       STAT             = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'FromHarmonics2Field3D - ModuleField4D - ERR20' 
                 
                    Field(i, j, k) = Field(i, j, k) + NewPropField%Harmonics%Residual3D(i, j, k)
                
                else                
                
                    Field(i, j, k) = FillValueReal
                
                endif   
            enddo                         
        enddo
        enddo
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR30'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, CoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FromHarmonics2Field2D - ModuleField4D - ERR40'        
    
        if(NewPropField%HasMultiplyingFactor)then
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Field(i,j,k) = Field(i,j,k) * NewPropField%MultiplyingFactor
            enddo
            enddo
            enddo
        end if

        if(NewPropField%HasAddingFactor)then
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                Field(i,j,k) = Field(i,j,k) + NewPropField%AddingFactor
            enddo
            enddo
            enddo
        end if
        
        if(NewPropField%MinValueON)then
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (Field(i,j,k) < NewPropField%MinValue) then
                    Field(i,j,k) = NewPropField%MinValue
                endif
            enddo
            enddo
            enddo
        end if        

        if(NewPropField%MaxValueON)then
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB                
                if (Field(i,j,k) > NewPropField%MaxValue) then
                    Field(i,j,k) = NewPropField%MaxValue
                endif
            enddo
            enddo
            enddo
        end if        

        
        nullify(Field)
    
    end subroutine FromHarmonics2Field3D

    !--------------------------------------------------------------------------    

 
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetField4DTimeLimits (Field4DID, StartTime, EndTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Field4DID
        type(T_Time)                                    :: StartTime, EndTime
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            StartTime = Me%File%StartTime
            EndTime   = Me%File%EndTime
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DTimeLimits
    
    !--------------------------------------------------------------------------
    subroutine GetField4DNumberOfInstants (Field4DID, NumberOfInstants, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Field4DID
        integer                                         :: NumberOfInstants
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        logical                                         :: TimeExist
        integer                                         :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call GetHDF5GroupExist (HDF5ID      = Me%File%Obj,                          &
                                    GroupName   = "/Time",                              &
                                    Exist       = TimeExist,                            &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DNumberOfInstants - ModuleField4D - ERR10'
            
            if (TimeExist) then            
            
                NumberOfInstants = Me%File%NumberOfInstants

            else
                NumberOfInstants = 1
            endif                
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DNumberOfInstants
    
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    type(T_Time) function GetField4DInstant (Field4DID, Instant, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Field4DID
        integer                                         :: Instant
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_, STAT_CALL
        logical                                         :: TimeExist

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call GetHDF5GroupExist (HDF5ID      = Me%File%Obj,                          &
                                    GroupName   = "/Time",                              &
                                    Exist       = TimeExist,                            &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DInstant - ModuleField4D - ERR10'
            
            if (TimeExist) then            
                GetField4DInstant = Me%File%InstantsDates(Instant)
            else
                GetField4DInstant = Me%StartTime
            endif            
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end function GetField4DInstant
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    real function GetField4DGeneric4DValue (Field4DID, Instant, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Field4DID
        integer                                         :: Instant
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        real,   dimension(:), pointer                   :: AuxVector
        integer                                         :: STAT_, ready_, STAT_CALL
        logical                                         :: Generic4DExist
        

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            call GetHDF5GroupExist (HDF5ID      = Me%File%Obj,                          &
                                    GroupName   = "/Generic4D",                         &
                                    Exist       = Generic4DExist,                       &
                                    STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DGeneric4DValue - ModuleField4D - ERR10'
            
            if (Generic4DExist) then  

                call HDF5SetLimits  (Me%File%Obj, 1, 1, STAT = STAT_CALL)

                allocate(AuxVector(1))
                                      
                call HDF5ReadData   (HDF5ID         = Me%File%Obj,                      &
                                     GroupName      = "/Generic4D",                     &
                                     Name           = "Generic4D",                      &
                                     Array1D        = AuxVector,                        &
                                     OutputNumber   = Instant,                          &
                                     STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    write(*,*) 'ObjHDF5=', Me%File%Obj 
                    write(*,*) 'FileName=', trim(Me%File%FileName)
                    write(*,*) 'STAT_CALL= ', STAT_CALL
                    stop 'GetField4DGeneric4DValue - ModuleField4D - ERR20'
                endif

                GetField4DGeneric4DValue = AuxVector(1)
     
                deallocate(AuxVector)                
                
            else
                write(*,*) "GroupName = /Generic4D does not exist in File=",trim(Me%File%FileName)
                stop 'GetField4DGeneric4DValue - ModuleField4D - ERR30'
            endif            
            
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end function GetField4DGeneric4DValue
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine GetField4DSize2D (Field4DID, Size2D, WorkSize2D, STAT)

        !Arguments-------------------------------------------------------------
        integer,        intent(IN)                      :: Field4DID
        type(T_Size2D), intent(OUT), optional           :: Size2D, WorkSize2D
        integer,        intent(OUT), optional           :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(Size2D    )) Size2D      = Me%Size2D
            if (present(WorkSize2D)) WorkSize2D  = Me%WorkSize2D
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DSize2D

    !--------------------------------------------------------------------------

    subroutine GetField4DSize3D (Field4DID, Size3D, WorkSize3D, STAT)

        !Arguments-------------------------------------------------------------
        integer,        intent(IN)                      :: Field4DID
        type(T_Size3D), intent(OUT), optional           :: Size3D, WorkSize3D
        integer,        intent(OUT), optional           :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            if (present(Size3D    )) Size3D      = Me%Size3D
            if (present(WorkSize3D)) WorkSize3D  = Me%WorkSize3D
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DSize3D
    

    !--------------------------------------------------------------------------
    
    subroutine GetField4DHarmonicsON(Field4DID, PropertyIDNumber, HarmonicsON, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN )                            :: Field4DID
        integer, intent(IN )                            :: PropertyIDNumber
        logical, intent(OUT)                            :: HarmonicsON
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DHarmonicsON - ModuleField4D - ERR10'
            
            HarmonicsON = PropField%Harmonics%ON
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DHarmonicsON

    !--------------------------------------------------------------------------
    
    subroutine GetField4DHarmonicsNumber(Field4DID, PropertyIDNumber, HarmonicsNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN )                            :: Field4DID
        integer, intent(IN )                            :: PropertyIDNumber
        integer, intent(OUT)                            :: HarmonicsNumber
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DHarmonicsNumber - ModuleField4D - ERR10'
            
            HarmonicsNumber = PropField%Harmonics%Number
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DHarmonicsNumber
    
    !--------------------------------------------------------------------------
    
    subroutine GetField4DHarmonicsName(Field4DID, PropertyIDNumber, HarmonicsName, STAT)

        !Arguments-------------------------------------------------------------
        integer,                                              intent(IN   ) :: Field4DID
        integer,                                              intent(IN   ) :: PropertyIDNumber
        character(Len=WaveNameLength), dimension(:), pointer, intent(INOUT) :: HarmonicsName
        integer, optional,                                    intent(OUT  ) :: STAT

        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DHarmonicsName - ModuleField4D - ERR10'
            
            HarmonicsName(:) = PropField%Harmonics%WaveName(:)
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DHarmonicsName
    
    !--------------------------------------------------------------------------    

    !--------------------------------------------------------------------------

    subroutine GetField4DBathymID (Field4DID, BathymID, STAT)

        !Arguments-------------------------------------------------------------
        integer,        intent(IN)                      :: Field4DID
        integer,        intent(OUT)                     :: BathymID
        integer,        intent(OUT), optional           :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            BathymID = Me%ObjBathymetry
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DBathymID
    

    !--------------------------------------------------------------------------

    subroutine GetField4DGridID (Field4DID, GridID, STAT)

        !Arguments-------------------------------------------------------------
        integer,        intent(IN)                      :: Field4DID
        integer,        intent(OUT)                     :: GridID
        integer,        intent(OUT), optional           :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            
            GridID = Me%ObjHorizontalGrid
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DGridID
    

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    
    subroutine GetField4DIsAngle(Field4DID, PropertyIDNumber, IsAngle, STAT)

        !Arguments-------------------------------------------------------------
        integer,                            intent(IN ) :: Field4DID
        integer,                            intent(IN ) :: PropertyIDNumber
        logical,                            intent(OUT) :: IsAngle
        integer, optional,                  intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                           &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetField4DIsAngle - ModuleField4D - ERR10'
            
            IsAngle = PropField%ID%IsAngle
            
            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetField4DIsAngle
    
    !-------------------------------------------------------------------------- 
 
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyField4D(Field4DID, PropertyIDNumber, CurrentTime, Matrix2D, Matrix3D, Instant, STAT)

        !Arguments-------------------------------------------------------------
        integer,                 intent(IN)             :: Field4DID
        integer,                 intent(IN)             :: PropertyIDNumber
        type (T_Time),           intent(IN)             :: CurrentTime        
        real,    dimension(:, :),    pointer, optional  :: Matrix2D
        real,    dimension(:, :, :), pointer, optional  :: Matrix3D
        integer,                 intent(IN),  optional  :: Instant
        integer,                 intent(OUT), optional  :: STAT

        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL
        integer                                         :: k
        integer                                         :: SizeI1, SizeJ1, SizeK1
        integer                                         :: SizeI2, SizeJ2, SizeK2        
        integer                                         :: ILB1,IUB1,JLB1,JUB1,KLB1,KUB1
        integer                                         :: ILB2,IUB2,JLB2,JUB2,KLB2,KUB2
        logical                                         :: CorrectTimeFrame        
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
    
            call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyField4D - ModuleField4D - ERR10'        

            Me%CurrentTimeExt = CurrentTime
        
            if (Me%BackTracking) then  
                call BacktrackingTime
            else   
                Me%CurrentTimeInt = Me%CurrentTimeExt
            endif

            if (PropField%Generic4D%ON) then
            
                if (present(Instant)) then
                
                    PropField%Generic4D%Instant   = Instant
                    PropField%Generic4D%InstantON = .true. 

                    if      (PropField%SpaceDim == Dim2D) then
                    
                        call ReadValues2D (PropField, Previous = .false.)
                        
                        Me%Matrix2D = PropField%NextField2D
                    
                    else if (PropField%SpaceDim == Dim3D) then
                    
                        call ReadValues3D (PropField, Previous = .false.)
                    
                        Me%Matrix3D = PropField%NextField3D
                    
                    endif            
                    
                else
                
                    stop 'ModifyField4D - ModuleField4D - ERR20'
                    
                endif                    
                
            else            
            
                CorrectTimeFrame = .true.
                
                if (Me%CurrentTimeInt < Me%File%StartTime) CorrectTimeFrame = .false.  
                if (Me%CurrentTimeInt > Me%File%EndTime  ) CorrectTimeFrame = .false.              
                
                if (PropField%Harmonics%ON) CorrectTimeFrame = .true.

                if (CorrectTimeFrame ) then

                    if      (PropField%SpaceDim == Dim2D) then
                        if (PropField%Harmonics%ON) then
                            !call FromHarmonics2Field2D(PropField, Me%CurrentTimeInt)                        
                            call ModifyInputHarmonics2D(PropField)
                        else
                            call ModifyInput2D        (PropField) 
                        endif                        
                    else if (PropField%SpaceDim == Dim3D) then
                        if (PropField%Harmonics%ON) then
                            call FromHarmonics2Field3D(PropField, Me%CurrentTimeInt)                        
                        else
                            call ModifyInput3D        (PropField) 
                        endif                        

                    endif
                    
                endif     
                
            endif                           
                
            if (Me%Output%Yes) then
                call WriteOutput(PropField, PropertyIDNumber)
            endif
            
            if (present(Matrix2D)) then
                Matrix2D(:,:) = Me%Matrix2D(:,:)
            endif

            if (present(Matrix3D)) then
                
                SizeI1  = size(Matrix3D, 1)
                SizeJ1  = size(Matrix3D, 2)
                SizeK1  = size(Matrix3D, 3)
                
                SizeI2 = size(Me%Matrix3D, 1)
                SizeJ2 = size(Me%Matrix3D, 2)                
                SizeK2 = size(Me%Matrix3D, 3)
                
                ILB1   = LBound(Matrix3D, 1)
                JLB1   = LBound(Matrix3D, 2)
                KLB1   = LBound(Matrix3D, 3) 
                
                IUB1   = UBound(Matrix3D, 1)
                JUB1   = UBound(Matrix3D, 2)                
                KUB1   = UBound(Matrix3D, 3)                 
                
                ILB2   = LBound(Me%Matrix3D, 1)
                JLB2   = LBound(Me%Matrix3D, 2)
                KLB2   = LBound(Me%Matrix3D, 3) 
                
                IUB2   = UBound(Me%Matrix3D, 1)
                JUB2   = UBound(Me%Matrix3D, 2)
                KUB2   = UBound(Me%Matrix3D, 3)                                
                
                
                if (SizeI1 /= SizeI2) then
                    stop 'ModifyField4D - ModuleField4D - ERR30'
                endif                        
                
                if (SizeJ1 /= SizeJ2) then
                    stop 'ModifyField4D - ModuleField4D - ERR40'
                endif                
                
                if      (PropField%From3Dto2D) then
                    !Need to be a matrix3D with one layer
                    if (SizeK1 /= 3) then
                        stop 'ModifyField4D - ModuleField4D - ERR50'
                    endif                        
                    
                    Matrix3D(:,:,1) = Me%Matrix3D(:,:, Me%WorkSize3D%KUB)
                    
                elseif (PropField%From2Dto3D) then    
                    
                    
                    do k = KLB1, KUB1
                        Matrix3D(ILB1:IUB1,JLB1:JUB1,k) = Me%Matrix3D(ILB2:IUB2,JLB2:JUB2,1)                    
                    enddo                        
                    
                    
                    
                else
                    if (SizeK1  /= SizeK2) then
                        stop 'ModifyField4D - ModuleField4D - ERR60'
                    endif                        
                    
                    
                    Matrix3D(ILB1:IUB1,JLB1:JUB1,KLB1:KUB1) = Me%Matrix3D(ILB2:IUB2,JLB2:JUB2,KLB2:KUB2)
                    
                endif                        
            endif
                                
            
            STAT_ = SUCCESS_
            
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyField4D

    !--------------------------------------------------------------------------
    
    subroutine WriteOutput(PropField, PropertyIDNumber)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer,                 intent(IN)             :: PropertyIDNumber

        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:),  pointer              :: SZZ
        real, dimension(6), target                      :: AuxTime
        real,    dimension(:    ), pointer              :: TimePtr   
        integer                                         :: WorkILB, WorkIUB,            &
                                                           WorkJLB, WorkJUB,            &
                                                           WorkKLB, WorkKUB
        integer                                         :: OutPutNumber, STAT_CALL 

        !Begin-----------------------------------------------------------------

        OutPutNumber = Me%OutPut%NextOutput
        
        if (OutPutNumber<=Me%OutPut%TotalOutputs) then
        if (Me%CurrentTimeExt>=Me%OutPut%OutTime(Me%OutPut%NextOutput)) then

            if (Me%BackTracking) then
                OutPutNumber = Me%OutPut%TotalOutputs - OutPutNumber + 1 
            endif 

            WorkILB = Me%WorkSize2D%ILB 
            WorkIUB = Me%WorkSize2D%IUB 

            WorkJLB = Me%WorkSize2D%JLB 
            WorkJUB = Me%WorkSize2D%JUB 

            
            !Writes current time
            call ExtractDate   (Me%CurrentTimeInt, AuxTime(1), AuxTime(2), AuxTime(3),  &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime
            call HDF5SetLimits  (Me%ObjHDF5Out, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR10'

            call HDF5WriteData  (Me%ObjHDF5Out, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                                 Array1D = TimePtr, OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR20'



            if (PropField%SpaceDim == Dim3D) then

                WorkKLB = Me%WorkSize3D%KLB
                WorkKUB = Me%WorkSize3D%KUB

                !Writes SZZ
                call HDF5SetLimits  (Me%ObjHDF5Out, WorkILB, WorkIUB, WorkJLB,          &
                                     WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR30'

                call GetGeometryDistances(  GeometryID      = Me%ObjGeometry,           &
                                            SZZ             = SZZ,                      &
!For Me%CurrentTimeInt a stationary geometry is assumed                                    
!                                    ActualTime      = Me%CurrentTimeInt,                             &
                                    STAT            = STAT_CALL)                                     
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR40'


                call HDF5WriteData  (Me%ObjHDF5Out, "/Grid/VerticalZ", trim(Me%File%DataSetVert),   &
                                     "m", Array3D =    SZZ,                             &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR50'
                
                call UnGetGeometry(  GeometryID      = Me%ObjGeometry,                  &
                                     Array           = SZZ,                             &
                                     STAT            = STAT_CALL)                                     
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR60'

                

                !Writes property 3D
                call HDF5SetLimits  (Me%ObjHDF5Out, WorkILB, WorkIUB,                   &
                                     WorkJLB, WorkJUB, WorkKLB, WorkKUB,                &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR70'
                
     
                call HDF5WriteData  (Me%ObjHDF5Out, "/Results/"//trim(GetPropertyName(PropertyIDNumber)), &
                                     trim(GetPropertyName(PropertyIDNumber)),           &
                                     "m", Array3D = Me%Matrix3D,                        &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR80'

            else

                !Writes property 2D
                call HDF5SetLimits  (Me%ObjHDF5Out, WorkILB, WorkIUB,                   &
                                     WorkJLB, WorkJUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR90'
                
     
                call HDF5WriteData  (Me%ObjHDF5Out, "/Results/"//trim(GetPropertyName(PropertyIDNumber)), &
                                     trim(GetPropertyName(PropertyIDNumber)),           &
                                     "m", Array2D = Me%Matrix2D,                        &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteOutput - ModuleField4D - ERR100'
            
            endif
                
            Me%OutPut%NextOutput = Me%OutPut%NextOutput + 1
            
        endif
        endif
            
    end subroutine WriteOutput
    
    !--------------------------------------------------------------------------
    subroutine ModifyField4DXYZ(Field4DID, PropertyIDNumber, CurrentTime, X, Y, Z,      &
                                Field, NoData, WaveName, ExtractAmplitudes,             &
                                ExtractPhaseReal, InterpolationDT, Instant, STAT)

        !Arguments-------------------------------------------------------------
        integer,                            intent(IN)             :: Field4DID
        integer,                            intent(IN)             :: PropertyIDNumber
        type (T_Time),                      intent(IN), optional   :: CurrentTime 
        real,    dimension(:),   pointer,   intent(IN)             :: X, Y
        real,    dimension(:),   pointer,   intent(IN),  optional  :: Z
        real,    dimension(:),   pointer,   intent(INOUT)          :: Field
        logical, dimension(:),   pointer,   intent(INOUT)          :: NoData
        character(len=*),                   intent(IN ), optional  :: WaveName
        logical,                            intent(IN ), optional  :: ExtractAmplitudes
        logical,                            intent(IN ), optional  :: ExtractPhaseReal
        real,                               intent(OUT), optional  :: InterpolationDT        
        integer,                            intent(IN),  optional  :: Instant
        integer,                            intent(OUT), optional  :: STAT
                                            
        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL
        logical                                         :: CorrectTimeFrame
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(InterpolationDT)) then                
                InterpolationDT = - FillValueReal
            endif                
        
            if (Me%WindowWithData) then
 
                call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DXYZ - ModuleField4D - ERR10'     
 
                if (present(CurrentTime)) then
                    Me%CurrentTimeExt = CurrentTime
                else
                    call null_time(Me%CurrentTimeExt)
                endif                    
            
                if (Me%BackTracking) then  
                    call BacktrackingTime
                else   
                    Me%CurrentTimeInt = Me%CurrentTimeExt
                endif
                

ifG4D:          if (PropField%Generic4D%ON) then

                    CorrectTimeFrame = .true.
                
                    if (present(Instant)) then
                    
                        PropField%Generic4D%Instant   = Instant
                        PropField%Generic4D%InstantON = .true.                         

                        if      (PropField%SpaceDim == Dim2D) then
                        
                            call ReadValues2D (PropField, Previous = .false.)
                            
                            Me%Matrix2D = PropField%NextField2D
                        
                        else if (PropField%SpaceDim == Dim3D) then
                        
                            call ReadValues3D (PropField, Previous = .false.)
                        
                            Me%Matrix3D = PropField%NextField3D
                        
                        endif            
                        
                    else
                    
                        stop 'ModifyField4D - ModuleField4D - ERR20'
                        
                    endif   
                
                else ifG4D
                
                    CorrectTimeFrame = .true.
                    
                    if (Me%CurrentTimeInt < Me%File%StartTime) CorrectTimeFrame = .false.  
                    if (Me%CurrentTimeInt > Me%File%EndTime  ) CorrectTimeFrame = .false.     
                    
                    if (PropField%Harmonics%ON) CorrectTimeFrame = .true.                                         
                    
CTF:                if (CorrectTimeFrame) then        
                
                        if (PropField%Harmonics%Extract) then
                        
                            if (present(WaveName)) then
                                PropField%Harmonics%ExtractWave = WaveName
                            else
                                stop 'ModifyField4DXYZ - ModuleField4D - ERR20'                                                
                            endif


                            if (present(ExtractAmplitudes)) then
                                PropField%Harmonics%ExtractAmp = ExtractAmplitudes
                            else
                                stop 'ModifyField4DXYZ - ModuleField4D - ERR30'                                                
                            endif

                            if (present(ExtractPhaseReal)) then
                                PropField%Harmonics%ExtractPhaseReal = ExtractPhaseReal
                            else
                                PropField%Harmonics%ExtractPhaseReal = .false. 
                            endif
                            
                        endif  

                        if      (PropField%SpaceDim == Dim2D) then
                        
                            if (PropField%Harmonics%ON) then
                                !call FromHarmonics2Field2D(PropField, Me%CurrentTimeInt)                        
                                call ModifyInputHarmonics2D(PropField)
                            else
                                call ModifyInput2D        (PropField) 
                            endif                        

                        else if (PropField%SpaceDim == Dim3D) then

                            if (PropField%Harmonics%ON) then
                                call FromHarmonics2Field3D(PropField, Me%CurrentTimeInt)                        
                            else
                                call ModifyInput3D        (PropField) 
                            endif                        

                        endif
                        
                    endif CTF                        
                        
                endif ifG4D    
                
                if (CorrectTimeFrame) then
                                       
                    if      (PropField%SpaceDim == Dim2D) then
                
                        call Interpolate2DCloud (PropField, X, Y, Field, NoData) 
                    
                    else if (PropField%SpaceDim == Dim3D) then
                
                        if (PropField%From2Dto3D) then
                            call Interpolate2DCloud3DMatrix (PropField, X, Y, Field, NoData) 
                        else
                
                        if (.not.present(Z)) then
                            stop 'ModifyField4DXYZ - ModuleField4D - ERR50'
                        endif

                        call Interpolate3DCloud (PropField, X, Y, Z, Field, NoData) 
                            
                        endif
                    endif
                    
                endif                    

                if (present(InterpolationDT)) then                                    
                    InterpolationDT = PropField%NextTime - PropField%PreviousTime
                endif

                                
                
                if (Me%Output%Yes) then
                    call WriteOutput(PropField, PropertyIDNumber)
                endif

            
            endif
            
            STAT_ = SUCCESS_
            
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyField4DXYZ
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    subroutine GetBathymXY(Field4DID, PropertyIDNumber, X, Y, Bathym, NoData, STAT)

        !Arguments-------------------------------------------------------------
        integer,                            intent(IN)             :: Field4DID
        integer,                            intent(IN)             :: PropertyIDNumber
        real,    dimension(:),   pointer,   intent(IN)             :: X, Y
        real,    dimension(:),   pointer,   intent(OUT)            :: Bathym
        logical, dimension(:),   pointer,   intent(INOUT)          :: NoData
        integer,                            intent(OUT), optional  :: STAT
                                            
        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (PropertyIDNumber /= bathymetry_) stop 'GetBathymXY - ModuleField4D - ERR10'
            
            call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GetBathymXY - ModuleField4D - ERR20'

            if      (PropField%SpaceDim == Dim2D) then
                call InterpolateBathym (X, Y, Bathym, NoData) 
            else if (PropField%SpaceDim == Dim3D) then
                stop 'GetBathymXY - ModuleField4D - ERR10'
            endif
            
            STAT_ = SUCCESS_
            
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBathymXY
    
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine SearchPropertyField(PropField, PropertyIDNumber, STAT)


        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer                 :: PropField
        integer         , optional, intent (IN)     :: PropertyIDNumber
        integer         , optional, intent (OUT)    :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_ 
        
        !----------------------------------------------------------------------

        STAT_  = UNKNOWN_

        PropField => Me%FirstPropField

do2 :   do while (associated(PropField)) 
if5 :       if (PropField%ID%IDNumber==PropertyIDNumber) then
                exit do2 
            else
                PropField => PropField%Next                 
            end if if5
        end do do2

       !A PropertyX was found
       if (associated(PropField)) then
            STAT_ = SUCCESS_  
        else
             write (*,*)'Property Not Found in Module Field4D ', &
                         trim(GetPropertyName(PropertyIDNumber))
            STAT_  = NOT_FOUND_ERR_  
        end if


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SearchPropertyField

    !----------------------------------------------------------------------

    subroutine Interpolate2DCloud (PropField, X, Y, Field, NoData) 
        
        !Arguments------------------------------------------------------------
        type (T_PropField),         pointer, intent(IN)     :: PropField
        real,       dimension(:),   pointer, intent(IN)     :: X, Y
        real,       dimension(:),   pointer, intent(INOUT)  :: Field
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real                                            :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        integer                                         :: MaskSW, MaskNW, MaskSE, MaskNE, MaskN, MaskS       
        real                                            :: X_W, X_E, Xv, Y_S, Y_N, Yv, PercI, PercJ  
        integer                                         :: STAT_CALL, nPoints, nP
        integer                                         :: jW, jE, iS, iN, i, j
        logical                                         :: InsideDomain

        !Begin----------------------------------------------------------------
        
        nPoints = size(X)
        

        call GetWaterPoints2D(HorizontalMapID   = Me%ObjHorizontalMap,              &
                              WaterPoints2D     = Me%ExternalVar%WaterPoints2D,     &
                              STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'Interpolate2DCloud - ModuleField4D - ERR10' 
        
        if (PropField%Extrapolate) then        
            call FillMatrix2D(Me%WorkSize2D%ILB,                             &
                              Me%WorkSize2D%IUB,                             &
                              Me%WorkSize2D%JLB,                             &
                              Me%WorkSize2D%JUB,                             &
                              Me%ExternalVar%Waterpoints2D,                  &
                              Me%Matrix2D,                                   &
                              FillGridMethod = PropField%ExtrapolateMethod)
        endif
        
dnP:    do nP = 1,nPoints      

            if (NoData(nP)) then
                InsideDomain = GetXYInsideDomain(Me%ObjHorizontalGrid, X(nP), Y(nP), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate2DCloud - ModuleValida4D - ERR10'
            
                if (.not. InsideDomain) then
                    cycle
                endif
                
                call GetXYCellZ(Me%ObjHorizontalGrid, X(nP), Y(nP), i, j, PercI, PercJ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate2DCloud - ModuleValida4D - ERR20'
                

                if (PercJ > 0.5) then
                    jW = j
                    jE = j+1
                    Xv = PercJ - 0.5
                else
                    jW = j-1
                    jE = j
                    Xv = PercJ + 0.5
                endif
                
                if (PercI > 0.5) then
                    iS = i
                    iN = i+1
                    Yv = PercI - 0.5
                else
                    iS = i-1
                    iN = i
                    Yv = PercI + 0.5
                endif            
                
                X_W = 0.
                X_E = 1
                Y_S = 0.                
                Y_N = 1.
                
                ValueSW     = Me%Matrix2D(iS, jW)
                ValueSE     = Me%Matrix2D(iS, jE)
                ValueNW     = Me%Matrix2D(iN, jW)
                ValueNE     = Me%Matrix2D(iN, jE)
                
                MaskSW      = Me%ExternalVar%Waterpoints2D(iS, jW)
                MaskSE      = Me%ExternalVar%Waterpoints2D(iS, jE)
                MaskNW      = Me%ExternalVar%Waterpoints2D(iN, jW)
                MaskNE      = Me%ExternalVar%Waterpoints2D(iN, jE)
                

                if (PropField%Extrapolate .or. .not. PropField%DiscardFillValues) then
                
                    NoData(nP) = .false.
                
                    if (ValueSW < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueSE < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueNW < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueNE < FillValueReal/1e4) NoData(nP) = .true.
                    
                    if (.not.  NoData(nP)) then
                    
                    
                        if (PropField%InterpolMethod == Bilinear2D_) then
                
                            ValueN      = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                            ValueS      = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                            
                            Field(nP)   = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)

                        elseif (PropField%InterpolMethod == NearestNeighbor2D_) then
                            if (Xv < 0.5) then
                                if (Yv < 0.5) then
                                     Field(nP) =  ValueSW                           
                                else
                                     Field(nP) =  ValueNW                                                           
                                endif
                            else
                                if (Yv < 0.5) then
                                     Field(nP) =  ValueSE                           
                                else
                                     Field(nP) =  ValueNE                                                           
                                endif                            
                            endif
                        endif
                        
                    endif
                    
                else
                
                    if (ValueSW < FillValueReal/1e4) ValueSW = 0.
                    if (ValueSE < FillValueReal/1e4) ValueSE = 0.                
                    if (ValueNW < FillValueReal/1e4) ValueNW = 0.                
                    if (ValueNE < FillValueReal/1e4) ValueNE = 0.     
                    
                    if (PropField%InterpolMethod == Bilinear2D_) then                                  
                    
                        if (MaskNW == WaterPoint .and. MaskNE == WaterPoint) then
                            ValueN = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                            MaskN  = 1
                        elseif (MaskNW == WaterPoint) then
                            ValueN = ValueNW
                            MaskN  = 1
                        elseif (MaskNE == WaterPoint) then
                            ValueN = ValueNE
                            MaskN  = 1
                        else
                            MaskN  = 0
                        endif

                        if (MaskSW == WaterPoint .and. MaskSE == WaterPoint) then
                            ValueS = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                            MaskS  = 1
                        elseif (MaskSW == WaterPoint) then
                            ValueS = ValueSW
                            MaskS  = 1
                        elseif (MaskSE == WaterPoint) then
                            ValueS = ValueSE
                            MaskS  = 1
                        else
                            MaskS  = 0
                        endif
                        
                        NoData(nP) = .false. 
                        
                        if (MaskN == WaterPoint .and. MaskS == WaterPoint) then
                            Field(nP) = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)
                        else if (MaskN == WaterPoint) then
                            Field(nP) = ValueN
                        else if (MaskS == WaterPoint) then
                            Field(nP) = ValueS
                        else
                            Field(nP)  = FillValueReal
                            NoData(nP) = .true. 
                        endif
                        
                    elseif (PropField%InterpolMethod == NearestNeighbor2D_) then

                        if (MaskNW == WaterPoint .and. MaskNE == WaterPoint) then
                            if (Xv < 0.5) then
                                ValueN = ValueNW
                             else
                                ValueN = ValueNE
                             endif
                            MaskN  = 1
                        elseif (MaskNW == WaterPoint) then
                            ValueN = ValueNW
                            MaskN  = 1
                        elseif (MaskNE == WaterPoint) then
                            ValueN = ValueNE
                            MaskN  = 1
                        else
                            MaskN  = 0
                        endif

                        if (MaskSW == WaterPoint .and. MaskSE == WaterPoint) then
                            if (Xv < 0.5) then
                                ValueS = ValueSW
                             else
                                ValueS = ValueSE
                             endif                            
                            MaskS  = 1
                        elseif (MaskSW == WaterPoint) then
                            ValueS = ValueSW
                            MaskS  = 1
                        elseif (MaskSE == WaterPoint) then
                            ValueS = ValueSE
                            MaskS  = 1
                        else
                            MaskS  = 0
                        endif
                        
                        NoData(nP) = .false. 
                        
                        if (MaskN == WaterPoint .and. MaskS == WaterPoint) then
                            if (Yv < 0.5) then
                                Field(nP) = ValueS
                            else
                                Field(nP) = ValueN
                            endif
                        else if (MaskN == WaterPoint) then
                            Field(nP) = ValueN
                        else if (MaskS == WaterPoint) then
                            Field(nP) = ValueS
                        else
                            Field(nP)  = FillValueReal
                            NoData(nP) = .true. 
                        endif

                    
                    endif                        
                    
                endif                
                
            endif
                            
        enddo dnP            

        call UnGetHorizontalMap(HorizontalMapID   = Me%ObjHorizontalMap,                &
                                Array             = Me%ExternalVar%WaterPoints2D,       &
                                STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'Interpolate2DCloud - ModuleField4D - ERR30' 



     end subroutine Interpolate2DCloud     
    !----------------------------------------------------------------------

    subroutine Interpolate2DCloud3DMatrix (PropField, X, Y, Field, NoData) 
        
        !Arguments------------------------------------------------------------
        type (T_PropField),         pointer, intent(IN)     :: PropField
        real,       dimension(:),   pointer, intent(IN)     :: X, Y
        real,       dimension(:),   pointer, intent(INOUT)  :: Field
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real                                            :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        integer                                         :: MaskSW, MaskNW, MaskSE, MaskNE, MaskN, MaskS       
        real                                            :: X_W, X_E, Xv, Y_S, Y_N, Yv, PercI, PercJ  
        integer                                         :: STAT_CALL, nPoints, nP
        integer                                         :: jW, jE, iS, iN, i, j, KUB
        logical                                         :: InsideDomain

        !Begin----------------------------------------------------------------
        
        nPoints = size(X)
        
        KUB = Me%WorkSize3D%KUB
        
        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                            &
                              WaterPoints3D     = Me%ExternalVar%WaterPoints3D,         &
                              STAT              = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'Interpolate2DCloud3DMatrix - ModuleField4D - ERR10' 

        if (PropField%Extrapolate) then        
            call FillMatrix3D(Me%WorkSize3D%ILB,                             &
                              Me%WorkSize3D%IUB,                             &
                              Me%WorkSize3D%JLB,                             &
                              Me%WorkSize3D%JUB,                             &
                              Me%WorkSize3D%KLB,                             &
                              Me%WorkSize3D%KUB,                             &
                              Me%ExternalVar%Waterpoints3D,                  &
                              Me%Matrix3D,                                   &
                              FillGridMethod = PropField%ExtrapolateMethod)
        endif         
        
dnP:    do nP = 1,nPoints      

            if (NoData(nP)) then
                InsideDomain = GetXYInsideDomain(Me%ObjHorizontalGrid, X(nP), Y(nP), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate2DCloud3DMatrix - ModuleValida4D - ERR10'
            
                if (.not. InsideDomain) then
                    cycle
                endif
                
                call GetXYCellZ(Me%ObjHorizontalGrid, X(nP), Y(nP), i, j, PercI, PercJ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate2DCloud3DMatrix - ModuleValida4D - ERR20'
                

                if (PercJ > 0.5) then
                    jW = j
                    jE = j+1
                    Xv = PercJ - 0.5
                else
                    jW = j-1
                    jE = j
                    Xv = PercJ + 0.5
                endif
                
                if (PercI > 0.5) then
                    iS = i
                    iN = i+1
                    Yv = PercI - 0.5
                else
                    iS = i-1
                    iN = i
                    Yv = PercI + 0.5
                endif            
                
                X_W = 0.
                X_E = 1
                Y_S = 0.                
                Y_N = 1.
                
                ValueSW     = Me%Matrix3D(iS, jW, KUB)
                ValueSE     = Me%Matrix3D(iS, jE, KUB)
                ValueNW     = Me%Matrix3D(iN, jW, KUB)
                ValueNE     = Me%Matrix3D(iN, jE, KUB)
                
                MaskSW      = Me%ExternalVar%Waterpoints3D(iS, jW, KUB)
                MaskSE      = Me%ExternalVar%Waterpoints3D(iS, jE, KUB)
                MaskNW      = Me%ExternalVar%Waterpoints3D(iN, jW, KUB)
                MaskNE      = Me%ExternalVar%Waterpoints3D(iN, jE, KUB)
                

                if (PropField%Extrapolate .or. .not. PropField%DiscardFillValues) then
                
                    NoData(nP) = .false.
                
                    if (ValueSW < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueSE < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueNW < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueNE < FillValueReal/1e4) NoData(nP) = .true.
                    
                    if (.not.  NoData(nP)) then
                    
                    
                        if (PropField%InterpolMethod == Bilinear2D_) then
                
                            ValueN      = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                            ValueS      = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                            
                            Field(nP)   = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)

                        elseif (PropField%InterpolMethod == NearestNeighbor2D_) then
                            if (Xv < 0.5) then
                                if (Yv < 0.5) then
                                     Field(nP) =  ValueSW                           
                                else
                                     Field(nP) =  ValueNW                                                           
                                endif
                            else
                                if (Yv < 0.5) then
                                     Field(nP) =  ValueSE                           
                                else
                                     Field(nP) =  ValueNE                                                           
                                endif                            
                            endif
                        endif
                        
                    endif
                    
                else
                
                    if (ValueSW < FillValueReal/1e4) ValueSW = 0.
                    if (ValueSE < FillValueReal/1e4) ValueSE = 0.                
                    if (ValueNW < FillValueReal/1e4) ValueNW = 0.                
                    if (ValueNE < FillValueReal/1e4) ValueNE = 0.     
                    
                    if (PropField%InterpolMethod == Bilinear2D_) then                                  
                    
                        if (MaskNW == WaterPoint .and. MaskNE == WaterPoint) then
                            ValueN = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                            MaskN  = 1
                        elseif (MaskNW == WaterPoint) then
                            ValueN = ValueNW
                            MaskN  = 1
                        elseif (MaskNE == WaterPoint) then
                            ValueN = ValueNE
                            MaskN  = 1
                        else
                            MaskN  = 0
                        endif

                        if (MaskSW == WaterPoint .and. MaskSE == WaterPoint) then
                            ValueS = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                            MaskS  = 1
                        elseif (MaskSW == WaterPoint) then
                            ValueS = ValueSW
                            MaskS  = 1
                        elseif (MaskSE == WaterPoint) then
                            ValueS = ValueSE
                            MaskS  = 1
                        else
                            MaskS  = 0
                        endif
                        
                        NoData(nP) = .false. 
                        
                        if (MaskN == WaterPoint .and. MaskS == WaterPoint) then
                            Field(nP) = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)
                        else if (MaskN == WaterPoint) then
                            Field(nP) = ValueN
                        else if (MaskS == WaterPoint) then
                            Field(nP) = ValueS
                        else
                            Field(nP)  = FillValueReal
                            NoData(nP) = .true. 
                        endif
                        
                    elseif (PropField%InterpolMethod == NearestNeighbor2D_) then

                        if (MaskNW == WaterPoint .and. MaskNE == WaterPoint) then
                            if (Xv < 0.5) then
                                ValueN = ValueNW
                             else
                                ValueN = ValueNE
                             endif
                            MaskN  = 1
                        elseif (MaskNW == WaterPoint) then
                            ValueN = ValueNW
                            MaskN  = 1
                        elseif (MaskNE == WaterPoint) then
                            ValueN = ValueNE
                            MaskN  = 1
                        else
                            MaskN  = 0
                        endif

                        if (MaskSW == WaterPoint .and. MaskSE == WaterPoint) then
                            if (Xv < 0.5) then
                                ValueS = ValueSW
                             else
                                ValueS = ValueSE
                             endif                            
                            MaskS  = 1
                        elseif (MaskSW == WaterPoint) then
                            ValueS = ValueSW
                            MaskS  = 1
                        elseif (MaskSE == WaterPoint) then
                            ValueS = ValueSE
                            MaskS  = 1
                        else
                            MaskS  = 0
                        endif
                        
                        NoData(nP) = .false. 
                        
                        if (MaskN == WaterPoint .and. MaskS == WaterPoint) then
                            if (Yv < 0.5) then
                                Field(nP) = ValueS
                            else
                                Field(nP) = ValueN
                            endif
                        else if (MaskN == WaterPoint) then
                            Field(nP) = ValueN
                        else if (MaskS == WaterPoint) then
                            Field(nP) = ValueS
                        else
                            Field(nP)  = FillValueReal
                            NoData(nP) = .true. 
                        endif

                    
                    endif                        
                    
                endif                
                
            endif
                            
        enddo dnP            

        call UnGetMap(Map_ID          = Me%ObjMap,                                      &
                      Array           = Me%ExternalVar%WaterPoints3D,                   &
                      STAT            = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'Interpolate2DCloud3DMatrix - ModuleField4D - ERR30' 



    end subroutine Interpolate2DCloud3DMatrix    

    !----------------------------------------------------------------------

    subroutine InterpolateBathym (X, Y, Bathym, NoData) 
        
        !Arguments------------------------------------------------------------
        real,       dimension(:),   pointer, intent(IN)     :: X, Y
        real,       dimension(:),   pointer, intent(INOUT)  :: Bathym
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real                                            :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        real                                            :: X_W, X_E, Xv, Y_S, Y_N, Yv, PercI, PercJ  
        integer                                         :: STAT_CALL, nPoints, nP
        integer                                         :: jW, jE, iS, iN, i, j
        integer                                         :: iT, i1, i2, i3, i4
        logical                                         :: InsideDomain

        !Begin----------------------------------------------------------------
        
        nPoints = size(X)

        call GetWaterPoints2D(HorizontalMapID   = Me%ObjHorizontalMap,                  &
                              WaterPoints2D     = Me%ExternalVar%WaterPoints2D,         &
                              STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'InterpolateBathym - ModuleField4D - ERR10' 

        call GetGridData  (GridDataID           = Me%ObjHorizontalMap,                  &
                           GridData2D           = Me%ExternalVar%Bathymetry,            &
                           STAT                 = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'InterpolateBathym - ModuleField4D - ERR20' 
        
        Me%Matrix2D(:,:) =  Me%ExternalVar%Bathymetry(:,:)         

        if (Me%Extrapolate) then        
            call FillMatrix2D(Me%WorkSize2D%ILB,                             &
                              Me%WorkSize2D%IUB,                             &
                              Me%WorkSize2D%JLB,                             &
                              Me%WorkSize2D%JUB,                             &
                              Me%ExternalVar%Waterpoints2D,                  &
                              Me%Matrix2D,                                   &
                              FillGridMethod = Me%ExtrapolateMethod)
        endif
        
dnP:    do nP = 1,nPoints      

            if (NoData(nP)) then
            
                InsideDomain = GetXYInsideDomain(Me%ObjHorizontalGrid, X(nP), Y(nP), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterpolateBathym - ModuleValida4D - ERR30'
            
                if (.not. InsideDomain) then
                    cycle
                endif
                
                call GetXYCellZ(Me%ObjHorizontalGrid, X(nP), Y(nP), i, j, PercI, PercJ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InterpolateBathym - ModuleValida4D - ERR40'
                

                if (PercJ > 0.5) then
                    jW = j
                    jE = j+1
                    Xv = PercJ - 0.5
                else
                    jW = j-1
                    jE = j
                    Xv = PercJ + 0.5
                endif
                
                if (PercI > 0.5) then
                    iS = i
                    iN = i+1
                    Yv = PercI - 0.5
                else
                    iS = i-1
                    iN = i
                    Yv = PercI + 0.5
                endif            
                
                X_W = 0.
                X_E = 1
                Y_S = 0.                
                Y_N = 1.
                
                ValueSW     = Me%Matrix2D(iS, jW)
                ValueSE     = Me%Matrix2D(iS, jE)
                ValueNW     = Me%Matrix2D(iN, jW)
                ValueNE     = Me%Matrix2D(iN, jE)
                
                i1=0;i2=0;i3=0;i4=0;
                
                if (ValueSW > -50) i1 = 1
                if (ValueSE > -50) i2 = 1 
                if (ValueNE > -50) i3 = 1
                if (ValueNW > -50) i4 = 1
                
                
                iT = i1+i2+i3+i4

                if (iT==4) then

                    ValueN      = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                    ValueS      = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                    Bathym(nP)  = LinearInterpolation (Y_S, ValueS , Y_N, ValueN,  Yv)

                    NoData(nP)  = .false. 
               
                else if (iT ==0) then
                    NoData(nP)  = .true. 
                else
                    Bathym(nP)  = (ValueSW*real(i1) + ValueSE*real(i2) + ValueNE*real(i3) + ValueNW*real(i4))/real(iT)
                    NoData(nP)  = .false. 
                endif
        
            endif
                            
        enddo dnP            

        call UnGetHorizontalMap(HorizontalMapID   = Me%ObjHorizontalMap,                &
                                Array             = Me%ExternalVar%WaterPoints2D,       &
                                STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'InterpolateBathym - ModuleField4D - ERR50' 


        call UnGetGridData(GridDataID           = Me%ObjHorizontalMap,                  &
                           Array                = Me%ExternalVar%Bathymetry,            &
                           STAT                 = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'InterpolateBathym - ModuleField4D - ERR60' 

     end subroutine InterpolateBathym     
    !----------------------------------------------------------------------


    subroutine Interpolate3DCloud (PropField, X, Y, Z, Field, NoData) 
        
        !Arguments------------------------------------------------------------
        type (T_PropField),         pointer, intent(IN)     :: PropField
        real,       dimension(:),   pointer, intent(IN)     :: X, Y, Z
        real,       dimension(:),   pointer, intent(INOUT)  :: Field
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real,   dimension(:,:,:),   pointer                 :: SZZ
        integer                                             :: STAT_CALL
        integer                                             :: i, j, k


        !Begin----------------------------------------------------------------
        
        
        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                            &
                              WaterPoints3D     = Me%ExternalVar%WaterPoints3D,         &
                              STAT              = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'Interpolate3DCloud - ModuleField4D - ERR10' 

        if (PropField%Extrapolate) then        
            call FillMatrix3D(Me%WorkSize3D%ILB,                             &
                              Me%WorkSize3D%IUB,                             &
                              Me%WorkSize3D%JLB,                             &
                              Me%WorkSize3D%JUB,                             &
                              Me%WorkSize3D%KLB,                             &
                              Me%WorkSize3D%KUB,                             &
                              Me%ExternalVar%Waterpoints3D,                  &
                              Me%Matrix3D,                                   &
                              FillGridMethod = PropField%ExtrapolateMethod)
        endif                                     
        call GetGeometryDistances(  GeometryID      = Me%ObjGeometry,                   &
                                            SZZ     = SZZ,                              &
!For Me%CurrentTimeInt a stationary geometry is assumed                                    
!                                    ActualTime      = Me%CurrentTimeInt,                      &
                                    STAT            = STAT_CALL)                                     
        if (STAT_CALL /= SUCCESS_) stop 'Interpolate3DCloud - ModuleValida4D - ERR20'
        
        
do1 :   do K = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
do2 :   do J = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
do3 :   do I = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
        
            Me%Depth3D(i,j,k) = (SZZ(i,j,k)+SZZ(i,j,k-1))/2. 

        end do do3
        end do do2
        end do do1
            
        if (PropField%Zdepths) then

    do4 :   do K = Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
    do5 :   do J = Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
    do6 :   do I = Me%WorkSize3D%ILB, Me%WorkSize3D%IUB        
            
                Me%Depth3D(i,j,k) = Me%Depth3D(i,j,k) - SZZ(i,j,Me%WorkSize3D%KUB)

            end do do6
            end do do5
            end do do4        
        
        endif            
        
        call UnGetGeometry       (  GeometryID      = Me%ObjGeometry,           &
                                    Array           = SZZ,                      &
                                    STAT            = STAT_CALL)                                     
        if (STAT_CALL /= SUCCESS_) stop 'Interpolate3DCloud - ModuleValida4D - ERR30'

        if (PropField%Extrapolate) then        
            call FillMatrix3D(Me%WorkSize3D%ILB,                             &
                              Me%WorkSize3D%IUB,                             &
                              Me%WorkSize3D%JLB,                             &
                              Me%WorkSize3D%JUB,                             &
                              Me%WorkSize3D%KLB,                             &
                              Me%WorkSize3D%KUB,                             &
                              Me%ExternalVar%Waterpoints3D,                  &
                              Me%Depth3D,                                    &
                              FillGridMethod = PropField%ExtrapolateMethod)
        endif
        
        call Interpolater3D(PropField           =  PropField,                           &
                            Matrix3D            =  Me%Matrix3D,                         &
                            Depth3D             =  Me%Depth3D,                          &
                            Mask3D              =  Me%ExternalVar%WaterPoints3D,        &
                            HorizontalGrid      =  Me%ObjHorizontalGrid,                &
                            KLB                 =  Me%WorkSize3D%KLB,                   &
                            KUB                 =  Me%WorkSize3D%KUB,                   &
                            X                   =  X,                                   &
                            Y                   =  Y,                                   &
                            Z                   =  Z,                                   &
                            Field               =  Field,                               &
                            NoData              =  NoData)
        
        call UnGetMap(Map_ID          = Me%ObjMap,                                      &
                      Array           = Me%ExternalVar%WaterPoints3D,                   &
                      STAT            = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'Interpolate3DCloud - ModuleField4D - ERR60' 

     end subroutine Interpolate3DCloud     
     
    !----------------------------------------------------------------------
    
    subroutine Interpolater3D(PropField, Matrix3D, Depth3D, Mask3D, HorizontalGrid,     &
                              KLB, KUB, X, Y, Z, Field, NoData)

        !Arguments------------------------------------------------------------
        type (T_PropField),         pointer, intent(IN)     :: PropField
        real,   dimension(:,:,:),   pointer, intent(IN)     :: Matrix3D        
        real,   dimension(:,:,:),   pointer, intent(IN)     :: Depth3D 
        integer,dimension(:,:,:),   pointer, intent(IN)     :: Mask3D
        integer                            , intent(IN)     :: HorizontalGrid
        integer                            , intent(IN)     :: KLB, KUB
        real,       dimension(:),   pointer, intent(IN)     :: X, Y, Z
        real,       dimension(:),   pointer, intent(INOUT)  :: Field
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real(8),    dimension(:),   pointer                 :: Depth1D, Matrix1D
        real                                                :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        integer                                             :: MaskSW, MaskNW, MaskSE, MaskNE, MaskN, MaskS
        real                                                :: X_W, X_E, Xv, Y_S, Y_N, Yv, PercI, PercJ  
        integer                                             :: STAT_CALL, nPoints, nP, k
        integer                                             :: jW, jE, iS, iN, i, j
        logical                                             :: InsideDomain

        !Begin----------------------------------------------------------------


        nPoints = size(X)  
        
        allocate(Depth1D (KLB-1:KUB+1))  
        allocate(Matrix1D(KLB-1:KUB+1))

dnP:    do nP = 1,nPoints      

            if (NoData(nP)) then
                InsideDomain = GetXYInsideDomain(HorizontalGrid, X(nP), Y(nP), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate3DCloud - ModuleValida4D - ERR40'
            
                if (.not. InsideDomain) then
                    cycle
                endif
                
                call GetXYCellZ(HorizontalGrid, X(nP), Y(nP), i, j, PercI, PercJ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate3DCloud - ModuleValida4D - ERR50'
                
                if (PercJ > 0.5) then
                    jW = j
                    jE = j+1
                    Xv = PercJ - 0.5
                else
                    jW = j-1
                    jE = j
                    Xv = PercJ + 0.5
                endif
                
                if (PercI > 0.5) then
                    iS = i
                    iN = i+1
                    Yv = PercI - 0.5
                else
                    iS = i-1
                    iN = i
                    Yv = PercI + 0.5
                endif            
                
                X_W = 0.
                X_E = 1
                Y_S = 0.                
                Y_N = 1.
                
                k = KUB
                
                if (k>1) then

                    Depth1D (KLB:KUB) = Depth3D   (iS, jW, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Matrix3D  (iS, jW, KLB:KUB)
                    ValueSW     = ValueAtDepthZ(Z(nP), KLB, KUB, Depth1D, Matrix1D)
                    
                    Depth1D (KLB:KUB) = Depth3D   (iS, jE, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Matrix3D  (iS, jE, KLB:KUB)
                    ValueSE     = ValueAtDepthZ(Z(nP), KLB, KUB, Depth1D, Matrix1D)

                    Depth1D (KLB:KUB) = Depth3D   (iN, jW, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Matrix3D  (iN, jW, KLB:KUB)
                    ValueNW     = ValueAtDepthZ(Z(nP), KLB, KUB, Depth1D, Matrix1D)

                    Depth1D (KLB:KUB) = Depth3D   (iN, jE, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Matrix3D  (iN, jE, KLB:KUB)
                    ValueNE     = ValueAtDepthZ(Z(nP), KLB, KUB, Depth1D, Matrix1D)
                    
                    
                    Depth1D (KLB:KUB) = Depth3D   (iS, jW, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Mask3D    (iS, jW, KLB:KUB)
                    MaskSW      = MaskAtDepthZ (Z(nP), KLB, KUB, Depth1D, Matrix1D)

                    Depth1D (KLB:KUB) = Depth3D   (iS, jE, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Mask3D    (iS, jE, KLB:KUB)
                    MaskSE      = MaskAtDepthZ (Z(nP), KLB, KUB, Depth1D, Matrix1D)

                    Depth1D (KLB:KUB) = Depth3D   (iN, jW, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Mask3D    (iN, jW, KLB:KUB)
                    MaskNW      = MaskAtDepthZ (Z(nP), KLB, KUB, Depth1D, Matrix1D)
                    
                    Depth1D (KLB:KUB) = Depth3D   (iN, jE, KLB:KUB) 
                    Matrix1D(KLB:KUB) = Mask3D    (iN, jE, KLB:KUB)
                    MaskNE      = MaskAtDepthZ (Z(nP), KLB, KUB, Depth1D, Matrix1D)
                                            
                else
                    ValueSW     = Matrix3D(iS, jW, k)
                    ValueSE     = Matrix3D(iS, jE, k)
                    ValueNW     = Matrix3D(iN, jW, k)
                    ValueNE     = Matrix3D(iN, jE, k)
                    
                    MaskSW      = Mask3D(iS, jW, k)
                    MaskSE      = Mask3D(iS, jE, k)
                    MaskNW      = Mask3D(iN, jW, k)
                    MaskNE      = Mask3D(iN, jE, k)
                    
                endif

                if (PropField%Extrapolate .or. .not. PropField%DiscardFillValues) then
                
                    NoData(nP) = .false.
                
                    if (ValueSW < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueSE < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueNW < FillValueReal/1e4) NoData(nP) = .true.
                    if (ValueNE < FillValueReal/1e4) NoData(nP) = .true.
                    
                    if (.not.  NoData(nP)) then
                
                        ValueN      = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                        ValueS      = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                        
                        Field(nP)   = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)

                    endif                    
                    
                else

                    if (ValueSW < FillValueReal/1e4) ValueSW = 0.
                    if (ValueSE < FillValueReal/1e4) ValueSE = 0.                
                    if (ValueNW < FillValueReal/1e4) ValueNW = 0.                
                    if (ValueNE < FillValueReal/1e4) ValueNE = 0.    

                    if (MaskNW == WaterPoint .and. MaskNE == WaterPoint) then
                        ValueN = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                        MaskN  = 1
                    elseif (MaskNW == WaterPoint) then
                        ValueN = ValueNW
                        MaskN  = 1
                    elseif (MaskNE == WaterPoint) then
                        ValueN = ValueNE
                        MaskN  = 1
                    else
                        MaskN  = 0
                    endif

                    if (MaskSW == WaterPoint .and. MaskSE == WaterPoint) then
                        ValueS = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                        MaskS  = 1
                    elseif (MaskSW == WaterPoint) then
                        ValueS = ValueSW
                        MaskS  = 1
                    elseif (MaskSE == WaterPoint) then
                        ValueS = ValueSE
                        MaskS  = 1
                    else
                        MaskS  = 0
                    endif
                    
                    NoData(nP) = .false. 
                    
                    if (MaskN == WaterPoint .and. MaskS == WaterPoint) then
                        Field(nP) = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)
                    else if (MaskN == WaterPoint) then
                        Field(nP) = ValueN
                    else if (MaskS == WaterPoint) then
                        Field(nP) = ValueS
                    else
                        Field(nP)  = FillValueReal
                        NoData(nP) = .true. 
                    endif
                    
                endif

            endif
                            
        enddo dnP   
        
        deallocate(Depth1D )
        deallocate(Matrix1D)

    end subroutine Interpolater3D
    
    !----------------------------------------------------------------------

    

    real(8) function ValueAtDepthZ(Z, KLB, KUB, Depth1D, Matrix1D)
    
        !Arguments-------------------------------------------------------------
        real(8), dimension(:), pointer :: Depth1D, Matrix1D
        real                           :: Z
        integer                        :: KLB, KUB

        !Local-----------------------------------------------------------------        
        integer     ::kb, Ndepths, k

        !Begin-----------------------------------------------------------------            

        do k = KUB,KLB+1,-1
            if (Depth1D (k-1)<Depth1D (k)) exit
        enddo
        
        kb = k
        
        Ndepths       = KUB - kb + 1

        ValueAtDepthZ = InterpolateProfileR8(dble(Z), Ndepths, Depth1D (kb:KUB), Matrix1D(kb:KUB))
        
    end function ValueAtDepthZ


    !----------------------------------------------------------------------

    integer function MaskAtDepthZ(Z, KLB, KUB, Depth1D, Matrix1D)
    
        !Arguments-------------------------------------------------------------
        real(8), dimension(:), pointer :: Depth1D, Matrix1D
        real                           :: Z
        integer                        :: KLB, KUB

        !Local-----------------------------------------------------------------        
        integer                        :: kb, Ndepths, k

        !Begin-----------------------------------------------------------------            

        do k = KUB,KLB+1,-1
            if (Depth1D (k-1)<Depth1D (k)) exit
        enddo
        
        kb = k
        
        Ndepths       = KUB - kb + 1

        MaskAtDepthZ = int(InterpolateProfileR8(dble(Z), Ndepths, Depth1D (kb:KUB), Matrix1D(kb:KUB)))
        
    end function MaskAtDepthZ

    !----------------------------------------------------------------------

    subroutine ModifyInput3D(PropField)
        
        !Arguments------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        
        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, n, i, j, k
        !Begin----------------------------------------------------------------

        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                    &
                              WaterPoints3D     = Me%ExternalVar%WaterPoints3D, &
                              STAT              = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) then
            stop 'ModifyInput3D - ModuleField4D - ERR10' 
        endif 


        if (ReadNewField(PropField,n))then
        
            if (n==1) then
                call SetMatrixValue(PropField%PreviousField3D, Me%WorkSize3D, PropField%NextField3D)
            else
                call ReadValues3D(PropField, Previous = .true. )
                
                !limit maximum values
                do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
                do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
                do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB

                    if (abs(PropField%PreviousField3D(i,j,k)) > abs(FillValueReal))     &
                            PropField%PreviousField3D(i,j,k) = FillValueReal

                enddo
                enddo
                enddo                
                
            endif

            call ReadValues3D(PropField, Previous = .false. )            
            
            !limit maximum values
            do k=Me%WorkSize3D%KLB, Me%WorkSize3D%KUB
            do j=Me%WorkSize3D%JLB, Me%WorkSize3D%JUB
            do i=Me%WorkSize3D%ILB, Me%WorkSize3D%IUB
            
                if (abs(PropField%NextField3D    (i,j,k)) > abs(FillValueReal))         &
                        PropField%NextField3D    (i,j,k) = FillValueReal

            enddo
            enddo
            enddo                
            
        end if

        if (PropField%PreviousInstant /= PropField%NextInstant) then
        
            if (PropField%ID%IsAngle) then                                        
                call InterpolateAngle3DInTime (ActualTime       = Me%CurrentTimeInt,        &
                                               Size             = Me%WorkSize3D,            &
                                               Time1            = PropField%PreviousTime,   &
                                               Matrix1          = PropField%PreviousField3D,&
                                               Time2            = PropField%NextTime,       &
                                               Matrix2          = PropField%NextField3D,    &
                                               MatrixOut        = Me%Matrix3D,              &
                                               PointsToFill3D   = Me%ExternalVar%WaterPoints3D)

            else
                call InterpolateMatrix3DInTime(ActualTime       = Me%CurrentTimeInt,        &
                                               Size             = Me%WorkSize3D,            &
                                               Time1            = PropField%PreviousTime,   &
                                               Matrix1          = PropField%PreviousField3D,&
                                               Time2            = PropField%NextTime,       &
                                               Matrix2          = PropField%NextField3D,    &
                                               MatrixOut        = Me%Matrix3D,              &
                                               PointsToFill3D   = Me%ExternalVar%WaterPoints3D)
            endif
        else
            
            !Prev and next are equal (last instant?)
            call SetMatrixValue(Me%Matrix3D, Me%WorkSize3D, PropField%NextField3D)

        endif

        call UnGetMap(Map_ID          = Me%ObjMap,                                      &
                      Array           = Me%ExternalVar%WaterPoints3D,                   &
                      STAT            = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) then
            stop 'ModifyInput3D - ModuleField4D - ERR20' 
        endif 



    end subroutine ModifyInput3D


    !----------------------------------------------------------------------------
    
    

    subroutine ModifyInputHarmonics2D(PropField)
        
        !Arguments------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin----------------------------------------------------------------

        
if1:    if (PropField%Harmonics%Extract) then    

            call FromHarmonics2Field2D(NewPropField = PropField)
            
            
        else if1      
        
            call GetWaterPoints2D(HorizontalMapID   = Me%ObjHorizontalMap,              &
                                  WaterPoints2D     = Me%ExternalVar%WaterPoints2D,     &
                                  STAT              = STAT_CALL) 

            if (STAT_CALL/=SUCCESS_) stop 'ModifyInputHarmonics2D - ModuleField4D - ERR10'         

            if (Me%CurrentTimeInt >= PropField%NextTime) then

                PropField%PreviousTime = PropField%NextTime
                PropField%NextTime     = PropField%NextTime + PropField%Harmonics%DT
            
                call FromHarmonics2Field2D(NewPropField = PropField, CurrentTime = PropField%NextTime) 
            
                PropField%PreviousField2D(:,:) = PropField%NextField2D(:,:)
                PropField%NextField2D    (:,:) = Me%Matrix2D          (:,:)
            

            end if
             
            call InterpolateMatrix2DInTime(ActualTime       = Me%CurrentTimeInt,            &
                                           Size             = Me%WorkSize2D,                &
                                           Time1            = PropField%PreviousTime,       &
                                           Matrix1          = PropField%PreviousField2D,    &
                                           Time2            = PropField%NextTime,           &
                                           Matrix2          = PropField%NextField2D,        &
                                           MatrixOut        = Me%Matrix2D,                  &
                                           PointsToFill2D   = Me%ExternalVar%WaterPoints2D)

            call UnGetHorizontalMap(HorizontalMapID   = Me%ObjHorizontalMap,                &
                                    Array             = Me%ExternalVar%WaterPoints2D,       &
                                    STAT              = STAT_CALL) 

            if (STAT_CALL/=SUCCESS_) stop 'ModifyInputHarmonics2D - ModuleField4D - ERR20' 

        endif if1 

    end subroutine ModifyInputHarmonics2D


    !--------------------------------------------------------------------------



    subroutine ModifyInput2D(PropField)
        
        !Arguments------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField

        !Local----------------------------------------------------------------
        integer                                         :: STAT_CALL, n, i, j

        !Begin----------------------------------------------------------------

        call GetWaterPoints2D(HorizontalMapID   = Me%ObjHorizontalMap,              &
                              WaterPoints2D     = Me%ExternalVar%WaterPoints2D,     &
                              STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'ModifyInput2D - ModuleField4D - ERR10' 

        

        if (ReadNewField(PropField,n))then

            if (n==1) then 
                call SetMatrixValue(PropField%PreviousField2D, Me%WorkSize2D, PropField%NextField2D, Me%ExternalVar%WaterPoints2D)
            else
                call ReadValues2D(PropField, Previous = .true. )
                
                !limit maximum values
                do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
                do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
                
                    if (abs(PropField%PreviousField2D(i,j)) > abs(FillValueReal))          &
                        PropField%PreviousField2D(i,j) = FillValueReal
                    
                enddo
                enddo   
                
            endif

            call ReadValues2D(PropField, Previous = .false. )
            
            !limit maximum values
            do j=Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
            do i=Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            
                if (abs(PropField%NextField2D(i,j)) > abs(FillValueReal))                  &
                    PropField%NextField2D    (i,j) = FillValueReal
            enddo
            enddo   
            

        end if

        if (PropField%PreviousInstant /= PropField%NextInstant) then

            if (PropField%ValuesType == OriginalValues) then
            
                Me%Matrix2D = PropField%NextField2D
                
                call PredictDTForHDF (PropField)
            
            else  if (PropField%ValuesType == AccumulatedValues) then       !For Rain
            
                Me%Matrix2D = PropField%NextField2D / (PropField%NextTime - PropField%PreviousTime)
    
                call PredictDTForHDF (PropField)
                
            else if (PropField%ValuesType == InterpolatedValues) then       !For interpolation

                !Interpolates the two matrixes in time
                if (PropField%ID%IsAngle) then                
                    call InterpolateAngle2DInTime (ActualTime       = Me%CurrentTimeInt,            &
                                                   Size             = Me%WorkSize2D,                &
                                                   Time1            = PropField%PreviousTime,       &
                                                   Matrix1          = PropField%PreviousField2D,    &
                                                   Time2            = PropField%NextTime,           &
                                                   Matrix2          = PropField%NextField2D,        &
                                                   MatrixOut        = Me%Matrix2D,                  &
                                                   PointsToFill2D   = Me%ExternalVar%WaterPoints2D)                
                else
                    call InterpolateMatrix2DInTime(ActualTime       = Me%CurrentTimeInt,            &
                                                   Size             = Me%WorkSize2D,                &
                                                   Time1            = PropField%PreviousTime,       &
                                                   Matrix1          = PropField%PreviousField2D,    &
                                                   Time2            = PropField%NextTime,           &
                                                   Matrix2          = PropField%NextField2D,        &
                                                   MatrixOut        = Me%Matrix2D,                  &
                                                   PointsToFill2D   = Me%ExternalVar%WaterPoints2D)
                endif
            endif
            
        else

            !Prev and next are equal (last instant?)
            if ( PropField%ValuesType == OriginalValues .or. PropField%ValuesType == InterpolatedValues) then

                call SetMatrixValue(Me%Matrix2D, Me%WorkSize2D, PropField%PreviousField2D, Me%ExternalVar%WaterPoints2D)

            else
            
                !do nothing
                
            endif

        endif

        call UnGetHorizontalMap(HorizontalMapID   = Me%ObjHorizontalMap,                &
                                Array             = Me%ExternalVar%WaterPoints2D,       &
                                STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'ModifyInput2D - ModuleField4D - ERR20' 



    end subroutine ModifyInput2D


    !--------------------------------------------------------------------------

    subroutine PredictDTForHDF(PropField)
        
        !Arguments-------------------------------------------------------------
        type(T_PropField), pointer                      :: PropField

        !Local-----------------------------------------------------------------
        type (T_Time)                                   :: Time1, Time2
        integer                                         :: i, j
        real                                            :: aux1, aux2
        logical                                         :: ValueDifferentZero
        



        Time1 = PropField%PreviousTime
        Time2 = PropField%NextTime        

        !Searches Maximum Value
        ValueDifferentZero = .false.
doj:    do j = Me%WorkSize2D%JLB, Me%WorkSize2D%JUB
        do i = Me%WorkSize2D%ILB, Me%WorkSize2D%IUB
            if(Me%ExternalVar%WaterPoints2D(i,j) == 1)then
                if (Me%Matrix2D(i, j) > PropField%MinForDTDecrease) then
                    ValueDifferentZero = .true.
                    exit doj
                endif
            endif
        enddo
        enddo doj

        if (ValueDifferentZero) then

            if (Time2 /= Me%CurrentTimeInt) then
                aux1 = Time2 - Me%CurrentTimeInt
            else
                aux1 = -FillValueReal
            endif

            !to ensure that DT does not go through two intervals
            aux2 = Time2 - Time1

            PropField%PredictedDT     = min(aux1, aux2)
            PropField%DTForNextEvent  = 0.0
        
        else

            !Can run until next Matrix will be read
            !This prediction is different from the one done by the GetTimeSerieDTForNextEvent
            !but I (Frank) assume that the datasets in the HDF always have a quite larger DT
            !then the model will use to run
            PropField%PredictedDT     = Time2 - Me%CurrentTimeInt
            PropField%DTForNextEvent  = Time2 - Me%CurrentTimeInt
        endif
                        

    end subroutine PredictDTForHDF

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillField4D(Field4DID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: Field4DID              
        integer, optional, intent(OUT)      :: STAT

        !Local-------------------------------------------------------------------
        type (T_PropField), pointer         :: PropField        
        integer                             :: ready_              
        integer                             :: STAT_, nUsers, STAT_CALL           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mField4D_,  Me%InstanceID)

            if (nUsers == 0) then
            
wwd:            if (Me%WindowWithData) then

                    if (Me%File%Obj /= 0) then
                    
                        PropField => Me%FirstPropField
                    
                        do while (associated(PropField))

                            if( associated(PropField%PreviousField2D))then
                                deallocate(PropField%PreviousField2D)
                                nullify   (PropField%PreviousField2D)
                            end if

                            if( associated(PropField%NextField2D))then
                                deallocate(PropField%NextField2D)
                                nullify   (PropField%NextField2D)
                            end if

                            if( associated(PropField%PreviousField3D))then
                                deallocate(PropField%PreviousField3D)
                                nullify   (PropField%PreviousField3D)
                            end if

                            if(associated (PropField%NextField3D))then
                                deallocate(PropField%NextField3D)
                                nullify   (PropField%NextField3D)
                            end if
                            
                            if(associated (PropField%Harmonics%WaveName))then
                                deallocate(PropField%Harmonics%WaveName)
                                nullify   (PropField%Harmonics%WaveName)
                            end if

                            if(associated (PropField%Harmonics%WaveGroupName))then
                                deallocate(PropField%Harmonics%WaveGroupName)
                                nullify   (PropField%Harmonics%WaveGroupName)
                            end if

                            
                            if(associated (PropField%Harmonics%Phase2D))then
                                deallocate(PropField%Harmonics%Phase2D)
                                nullify   (PropField%Harmonics%Phase2D)
                            end if

                            if(associated (PropField%Harmonics%Phase2DReal))then
                                deallocate(PropField%Harmonics%Phase2DReal)
                                nullify   (PropField%Harmonics%Phase2DReal)
                            end if

                            if(associated (PropField%Harmonics%Phase2DImag))then
                                deallocate(PropField%Harmonics%Phase2DImag)
                                nullify   (PropField%Harmonics%Phase2DImag)
                            end if
                            
                            if(associated (PropField%Harmonics%Amplitude2D))then
                                deallocate(PropField%Harmonics%Amplitude2D)
                                nullify   (PropField%Harmonics%Amplitude2D)
                            end if
                                                        

                            if(associated (PropField%Harmonics%Phase3D))then
                                deallocate(PropField%Harmonics%Phase3D)
                                nullify   (PropField%Harmonics%Phase3D)
                            end if

                            if(associated (PropField%Harmonics%Phase3DReal))then
                                deallocate(PropField%Harmonics%Phase3DReal)
                                nullify   (PropField%Harmonics%Phase3DReal)
                            end if

                            if(associated (PropField%Harmonics%Phase3DImag))then
                                deallocate(PropField%Harmonics%Phase3DImag)
                                nullify   (PropField%Harmonics%Phase3DImag)
                            end if

                            if(associated (PropField%Harmonics%Amplitude3D))then
                                deallocate(PropField%Harmonics%Amplitude3D)
                                nullify   (PropField%Harmonics%Amplitude3D)
                            end if

                            if(associated (PropField%Harmonics%Residual3D))then
                                deallocate(PropField%Harmonics%Residual3D)
                                nullify   (PropField%Harmonics%Residual3D)
                            end if

                            if(associated (PropField%Harmonics%Residual2D))then
                                deallocate(PropField%Harmonics%Residual2D)
                                nullify   (PropField%Harmonics%Residual2D)
                            end if
                            
                            PropField => PropField%Next

                        enddo
                        
                        if (associated(Me%Matrix2D)) then
                            deallocate(Me%Matrix2D)
                            nullify   (Me%Matrix2D)
                        endif                    

                        if (associated(Me%Matrix3D)) then
                            deallocate(Me%Matrix3D)
                            nullify   (Me%Matrix3D)
                        endif                                        

                        if (associated(Me%Depth3D)) then
                            deallocate(Me%Depth3D)
                            nullify   (Me%Depth3D)
                        endif                                        
                        
                        if (Me%File%Form == HDF5_) then
                            call KillHDF5(Me%File%Obj, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR30'
                        endif
                        
                        if (Me%OutPut%Yes) then
                            call KillHDF5 (Me%ObjHDF5Out, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR40'
                            
                            if (associated(Me%OutPut%OutTime)) then
                                deallocate(Me%OutPut%OutTime, STAT = STAT_CALL)
                                if (STAT_CALL /= SUCCESS_) &
                                    stop 'KillField4D - ModuleField4D - ERR50'
                                nullify   (Me%OutPut%OutTime)
                            end if
                            
                        endif                    
                          
                    endif
                
                    if (Me%ObjMap /= 0) then
                        if (Me%BuildMap) then
                            call KillMap(Me%ObjMap, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR60'
                        else
                            nUsers =  DeassociateInstance (mMAP_,  Me%ObjMap)
                            if (nUsers == 0) stop 'KillField4D - ModuleField4D - ERR70'
                        endif
                    endif

            
                    if (Me%ObjGeometry /= 0) then
                        if (Me%BuildGeometry) then
                            call KillGeometry(Me%ObjGeometry, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR80'
                        else                
                            nUsers = DeassociateInstance (mGEOMETRY_, Me%ObjGeometry)
                            if (nUsers == 0) stop 'KillField4D - ModuleField4D - ERR90'
                        endif
                    endif

                    if (Me%ObjHorizontalMap /= 0) then
                        if (Me%BuildHorizontalMap) then
                            call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR100'
                        else                  
                            nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                            if (nUsers == 0) stop 'KillField4D - ModuleField4D - ERR110'                          
                        endif
                    endif
                    
                    if (Me%ObjBathymetry  /= 0) then
                        if (Me%BuildBathymetry) then
                            call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR120'
                        else                   
                            nUsers = DeassociateInstance (mGRIDDATA_, Me%ObjBathymetry)
                            if (nUsers == 0) stop 'KillField4D - ModuleField4D - ERR130'
                        endif
                    endif

                    if (Me%ObjHorizontalGrid  /= 0) then
                        if (Me%BuildHorizontalGrid) then
                            call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
                            if (STAT_CALL /= SUCCESS_) stop 'KillField4D - ModuleField4D - ERR140'
                        else 
                            nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid )
                            if (nUsers == 0) stop 'KillField4D - ModuleField4D - ERR150'
                        endif
                    endif
                    
                end if wwd                    
                    
                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime           )
                if (nUsers == 0) stop 'KillField4D - ModuleField4D - ERR160'


                !Deallocates Instance
                call DeallocateInstance ()

                Field4DID = 0
                STAT_           = SUCCESS_

                
            endif           

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillField4D
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Field4D), pointer          :: AuxObjField4D
        type (T_Field4D), pointer          :: PreviousObjField4D

        !Updates pointers
        if (Me%InstanceID == FirstObjField4D%InstanceID) then
            FirstObjField4D => FirstObjField4D%Next
        else
            PreviousObjField4D => FirstObjField4D
            AuxObjField4D      => FirstObjField4D%Next
            do while (AuxObjField4D%InstanceID /= Me%InstanceID)
                PreviousObjField4D => AuxObjField4D
                AuxObjField4D      => AuxObjField4D%Next
            enddo

            !Now update linked list
            PreviousObjField4D%Next => AuxObjField4D%Next

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

    subroutine Ready (ObjField4D_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjField4D_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjField4D_ID > 0) then
            call LocateObjField4D (ObjField4D_ID)
            ready_ = VerifyReadLock (mField4D_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjField4D (ObjField4DID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjField4DID

        !Local-----------------------------------------------------------------

        Me => FirstObjField4D
        do while (associated (Me))
            if (Me%InstanceID == ObjField4DID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) then
            stop 'ModuleField4D - LocateObjField4D - ERR01'
        endif

    end subroutine LocateObjField4D

    !--------------------------------------------------------------------------

end module ModuleField4D

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------
