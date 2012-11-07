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

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions,        only : InterpolateMatrix2DInTime,                       &
                                       InterpolateMatrix3DInTime,                       &
                                       SetMatrixValue, ConstructPropertyID,             &
                                       FillMatrix2DNearestCell, LinearInterpolation,    &
                                       FillMatrix3DNearestCell, InterpolateProfileR8
    use ModuleDrawing,          only : ArrayPolygonWindow           
    use ModuleHorizontalGrid,   only : GetHorizontalGridSize, ConstructHorizontalGrid,  &
                                       WriteHorizontalGrid,                             &
                                       GetXYInsideDomain, GetXYCellZ, KillHorizontalGrid
    use ModuleGridData,         only : ConstructGridData, GetGridData, UngetGridData, KillGridData
    use ModuleHorizontalMap 
    use ModuleGeometry,         only : ConstructGeometry, GetGeometrySize,              &
                                       ComputeInitialGeometry, GetGeometryDistances,    &
                                       UnGetGeometry, KillGeometry
    use ModuleMap 
    use ModuleHDF5,             only : ConstructHDF5, HDF5ReadData,                     &
                                       GetHDF5FileAccess, GetHDF5GroupNumberOfItems,    &
                                       HDF5SetLimits, GetHDF5ArrayDimensions, KillHDF5, &
                                       HDF5WriteData, HDF5FlushMemory, HDF5WriteData,   &
                                       GetHDF5GroupExist
#ifndef _NO_NETCDF                                       
    ! Manages NetCDF files
    use ModuleNetCDF,           only : GetNCDFFileAccess, ConstructNETCDF,              &
                                       NETCDFReadGrid2D, NETCDFReadTime,                &
                                       NETCDFGetDimensions, NETCDFReadData,             &
                                       NETCDFReadVert
#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif
#endif
    use ModuleTimeSerie
     
                                       


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
    public  :: GetField4DSize2D
    public  :: GetField4DSize3D                         
    
    !Modifier
    public  :: ModifyField4D
    public  :: ModifyField4DXYZ
    private ::      ModifyInput2D
    private ::      ModifyInput3D
    public  :: GetBathymXY
    private ::      InterpolateBathym

 
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
    

    !Variable from file
    integer, parameter                              :: None             = 1

    !type of values
    integer, parameter                              :: InterpolatedValues = 1
    integer, parameter                              :: AccumulatedValues  = 2
    integer, parameter                              :: OriginalValues     = 3
    
    !File formats
    integer, parameter                              :: HDF5_            = 1
    integer, parameter                              :: NetCDF_          = 2
    
    !Property types
    integer, parameter                              :: Scalar_          = 1
    integer, parameter                              :: VectorX_         = 2
    integer, parameter                              :: VectorY_         = 3

    !Parameter-----------------------------------------------------------------

    !Types---------------------------------------------------------------------

    type T_DefaultNames
        character(Len=StringLength) :: bat       
        character(Len=StringLength) :: lon_stag  
        character(Len=StringLength) :: lat_stag  
        character(Len=StringLength) :: mask      
        character(Len=StringLength) :: depth_stag
    end type T_DefaultNames


    !Generic 4D
    type T_Generic4D
        logical                            :: ON
        logical                            :: ReadFromTimeSerie
        integer                            :: ObjTimeSerie
        integer                            :: TimeSerieColumn
        real                               :: CurrentValue

    end type 
    
    type T_ExternalVar
        integer, dimension(:,:  ), pointer :: WaterPoints2D
        integer, dimension(:,:,:), pointer :: WaterPoints3D
        real,    dimension(:,:  ), pointer :: Bathymetry
    end type T_ExternalVar
    
    type       T_OutPut
         type (T_Time), pointer, dimension(:)   :: OutTime
         integer                                :: TotalOutputs
         integer                                :: NextOutPut
         logical                                :: Yes                  = .false.
         logical                                :: Run_End              = .false.
    end type T_OutPut
    

    type T_File
        character(PathLength)                       :: FileName 
        type (T_Time)                               :: StartTime,  EndTime        
        integer                                     :: Obj = 0
        integer                                     :: NumberOfInstants
        type (T_Time), dimension(:), pointer        :: InstantsDates
        logical                                     :: CyclicTimeON = .false.
        integer                                     :: Form = HDF5_
        character(StringLength)                     :: LonStagName, LatStagName
        character(StringLength)                     :: DepthStagName, MaskName, BathymName
        type (T_DefaultNames)                       :: DefaultNames        
    end type T_File
    
    type T_PropField
        character(PathLength)                       :: FieldName, VGroupPath
        real                                        :: MultiplyingFactor
        logical                                     :: HasMultiplyingFactor = .false.
        real                                        :: AddingFactor
        logical                                     :: HasAddingFactor = .false.
        logical                                     :: From2Dto3D   = .false.        
        type (T_Time)                               :: NextTime,  PreviousTime
        type (T_PropertyID)                         :: ID
        type(T_Generic4D)                           :: Generic4D
       
        integer                                     :: SpaceDim             !2D/3D
        integer                                     :: TypeZUV              !Z/U/V
        
        logical                                     :: ChangeInTime
        
        integer                                     :: ValuesType
        real                                        :: Next4DValue     = FillValueReal
        real                                        :: Previous4DValue = FillValueReal
        integer                                     :: NextInstant, PreviousInstant 
        real, dimension(:,:  ), pointer             :: PreviousField2D, NextField2D
        real, dimension(:,:,:), pointer             :: PreviousField3D, NextField3D
        
        real                                        :: MinForDTDecrease     = AllmostZero
        real                                        :: DefaultValue
        real                                        :: PredictedDT          = -null_real
        real                                        :: DTForNextEvent       = -null_real
        type(T_PropField), pointer                  :: Next        
    end type T_PropField            


    private :: T_Field4D
    type       T_Field4D
        integer                                     :: InstanceID
        type (T_Size2D)                             :: Size2D, WorkSize2D
        type (T_Size3D)                             :: Size3D, WorkSize3D
        
        real,    dimension(:, :   ), pointer        :: Matrix2D
        real,    dimension(:, :, :), pointer        :: Matrix3D, Depth3D
        !Auxiliar matrixes in the interpolation process
        real(8), dimension(:      ), pointer        :: Matrix1D, Depth1D        
        
        type (T_Time)                               :: StartTime, EndTime
        type (T_Time)                               :: CurrentTimeExt, CurrentTimeInt        
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

        integer                                     :: MaskDim              = Dim3D
        real                                        :: LatReference
        real                                        :: LonReference
        logical                                     :: ReadWindow, WindowWithData 
        real,    dimension(2,2)                     :: WindowLimits
        type(T_Field4D), pointer                    :: Next
    end type  T_Field4D

    !Global Module Variables
    type (T_Field4D), pointer                    :: FirstObjField4D
    type (T_Field4D), pointer                    :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructField4D(Field4DID, EnterDataID, ExtractType, FileName, MaskDim,&
                                TimeID, HorizontalGridID, BathymetryID,                 &
                                HorizontalMapID, GeometryID, MapID, LatReference,       &
                                LonReference, LatMin, LatMax, LonMin, LonMax,           &
                                Extrapolate, STAT)

        !Arguments---------------------------------------------------------------
        integer,                intent(INOUT)               :: Field4DID
        integer,                intent(IN )                 :: EnterDataID
        integer,                intent(IN )                 :: ExtractType        
        character(*),           intent(IN )                 :: FileName
        integer,                intent(IN )                 :: MaskDim        
        integer,                intent(IN )                 :: TimeID
        integer,      optional, intent(IN )                 :: HorizontalGridID
        integer,      optional, intent(IN )                 :: BathymetryID
        integer,      optional, intent(IN )                 :: HorizontalMapID
        integer,      optional, intent(IN )                 :: GeometryID
        integer,      optional, intent(IN )                 :: MapID
        real,         optional, intent(IN )                 :: LatReference
        real,         optional, intent(IN )                 :: LonReference        
        real,         optional, intent(IN )                 :: LatMin
        real,         optional, intent(IN )                 :: LatMax
        real,         optional, intent(IN )                 :: LonMin
        real,         optional, intent(IN )                 :: LonMax
        logical,      optional, intent(IN )                 :: Extrapolate
        integer,      optional, intent(OUT)                 :: STAT     

        !Local-------------------------------------------------------------------
        type (T_PropField), pointer                         :: NewPropField        
        integer                                             :: ready_, STAT_, nUsers, STAT_CALL

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

            Me%File%FileName     = trim(FileName)
            Me%MaskDim           = MaskDim
            
            Me%WindowWithData    = .true. 

            call ConstructFile(ExtractType)
            
            if (present(Extrapolate)) then
                Me%Extrapolate = Extrapolate
            else
                Me%Extrapolate = .true. 
            endif
            
            if (present(HorizontalGridID)) then
                Me%ObjHorizontalGrid        = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
                
                Me%BuildHorizontalGrid      = .false.
            else
                Me%BuildHorizontalGrid      = .true.
                            
                if (present(LatReference)) then
                    Me%LatReference = LatReference
                else
                    stop 'ConstructField4D - ModuleField4D - ERR10' 
                endif

                if (present(LonReference)) then
                    Me%LonReference = LonReference
                else
                    stop 'ConstructField4D - ModuleField4D - ERR20' 
                endif
                
                Me%ReadWindow = .true.

                if (present(LatMin)) then                
                    Me%WindowLimits(2,1) = LatMin
                else
                    Me%ReadWindow = .false.
                endif

                if (present(LatMax)) then                
                    Me%WindowLimits(2,2) = LatMax
                else
                    Me%ReadWindow = .false.
                endif

                if (present(LonMin)) then                
                    Me%WindowLimits(1,1) = LonMin
                else
                    Me%ReadWindow = .false.
                endif

                if (present(LonMax)) then                
                    Me%WindowLimits(1,2) = LonMax
                else
                    Me%ReadWindow = .false.
                endif
                
                call ReadGridFromFile
                
            endif

wwd:        if (Me%WindowWithData) then            
            
                call GetHorizontalGridSize(HorizontalGridID = Me%ObjHorizontalGrid,         &
                                           Size             = Me%Size2D,                    &
                                           WorkSize         = Me%WorkSize2D,                &
                                           STAT             = STAT_CALL) 
                if (STAT_CALL/=SUCCESS_) then
                    stop 'ConstructField4D - ModuleField4D - ERR30' 
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
                    stop 'ConstructField4D - ModuleField4D - ERR40' 
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
                            stop 'ConstructField4D - ModuleField4D - ERR50' 
                        endif 
                    
                        Me%BuildMap      = .false.
                    else
                        Me%BuildMap      = .true.

                        call ReadMap3DFromFile
                    endif
               
                endif
                
                call AllocatePropertyField  (NewPropField, Me%MaskDim)
                
                call ReadOptions            (NewPropField, ExtractType)     
                
                call ConstructPropertyField (NewPropField)
                
                if (Me%OutPut%Yes) then
                    call Open_HDF5_OutPut_File(NewPropField)
                endif
                
                call UnGetHorizontalMap(HorizontalMapID = Me%ObjHorizontalMap,              &
                                        Array           = Me%ExternalVar%WaterPoints2D,     &
                                        STAT            = STAT_CALL) 
                if (STAT_CALL/=SUCCESS_) then
                    stop 'ConstructField4D - ModuleField4D - ERR60' 
                endif 
                
                if (Me%MaskDim == Dim3D) then
                
                    call UnGetMap(Map_ID          = Me%ObjMap,                              &
                                  Array           = Me%ExternalVar%WaterPoints3D,           &
                                  STAT            = STAT_CALL) 
                    if (STAT_CALL/=SUCCESS_) then
                        stop 'ConstructField4D - ModuleField4D - ERR70' 
                    endif 
                
                endif
            
            endif wwd
                         
            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'ConstructField4D - ModuleField4D - ERR60' 
             
            
            !Returns ID
            Field4DID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ConstructField4D - ModuleField4D - ERR20' 

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
    
    subroutine ReadGridFromFile

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        real(8),   pointer, dimension(:,:)      :: LatR8, LonR8, LatStagR8, LonStagR8
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
                                        
            call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
                              GroupName     = "/Grid",                                  &
                              Name          = trim(Me%File%LatStagName),                &
                              Array2D       = LatStag,                                  &
                              STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR40'
            
            call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
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
        
        if (Me%ReadWindow) then
        
            allocate (WindowDomain(2,2))
            
            WindowDomain(:,:) = FillValueInt
        
            call ArrayPolygonWindow(XX              = LonStag,                          &
                                    YY              = LatStag,                          &
                                    WIn             = Me%WindowLimits,                  &
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
                                        WIn             = Me%WindowLimits,                  &
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
            
            ILB = max(WindowDomain(1,1) - 1,1   )
            IUB = min(WindowDomain(1,2) + 1,Imax)
            JLB = max(WindowDomain(2,1) - 1,1   )
            JUB = min(WindowDomain(2,2) + 1,Jmax)
            
            deallocate (WindowDomain)
            
wwd1:        if (Me%WindowWithData) then

                allocate(LatStagW(ILB-1:IUB+1,JLB-1:JUB+1))
                LatStagW(ILB-1:IUB+1,JLB-1:JUB+1) = LatStag(ILB-1:IUB+1,JLB-1:JUB+1)

                allocate(LonStagW(ILB-1:IUB+1,JLB-1:JUB+1))
                LonStagW(ILB-1:IUB+1,JLB-1:JUB+1) = LonStag(ILB-1:IUB+1,JLB-1:JUB+1)
            
            endif wwd1
        else
        
            ILB = 1
            IUB = Imax
            JLB = 1
            JUB = Jmax
            
            LatStagW => LatStag
            LonStagW => LonStag
            
        endif
        
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
            if (STAT_CALL /= SUCCESS_) stop 'ReadGridFromFile - ModuleField4D - ERR80'

        endif    
        
        deallocate(Lat    )
        deallocate(Lon    )
        deallocate(LatStag)
        deallocate(LonStag)
        if (Me%ReadWindow .and. Me%WindowWithData) then
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
                                        
            call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
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
        integer                                 :: ILB, IUB, JLB, JUB, Kmax, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        allocate(Mask2D  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB))
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB        
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        
        
        if (Me%MaskDim == Dim3D) then                                    
        
           !Get grid vertical space dimensions
            if      (Me%File%Form == HDF5_  ) then

                call GetHDF5ArrayDimensions (HDF5ID = Me%File%Obj, GroupName = "/Grid",     &
                                            ItemName = trim(Me%File%MaskName),              &
                                            Kmax = Kmax, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR10'
#ifndef _NO_NETCDF               
            else if (Me%File%Form == NetCDF_) then
            
                call NETCDFGetDimensions (NCDFID = Me%File%Obj, KUB = Kmax, STAT = STAT_CALL)

                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR20'
#endif                
            endif
            
            allocate(Mask3D  (Me%Size2D%ILB:Me%Size2D%IUB,Me%Size2D%JLB:Me%Size2D%JUB,1:Kmax))  
            
        
           !Read horizontal grid
            if      (Me%File%Form == HDF5_  ) then
                
                call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,            &
                                                           JLB = JLB, JUB = JUB,            &
                                                           KLB = 1, KUB = Kmax,          &
                                     STAT   = STAT_CALL)                                
                                     
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR30'
                                            
                call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
                                  GroupName     = "/Grid",                                  &
                                  Name          = trim(Me%File%MaskName),                   &
                                  Array3D       = Mask3D,                                   &
                                  STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR40'
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
            
            Mask2D(:,:) = Mask3D(:,:,Kmax)
            
            deallocate(Mask3D)
            nullify   (Mask3D)
            

        else        
                    
           !Read horizontal grid
            if      (Me%File%Form == HDF5_  ) then
                
                call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,            &
                                                           JLB = JLB, JUB = JUB,            &
                                     STAT   = STAT_CALL)                                   
                                     
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR60'
                                            
                call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
                                  GroupName     = "/Grid",                                  &
                                  Name          = trim(Me%File%MaskName),                   &
                                  Array2D       = Mask2D,                                   &
                                  STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR70'
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
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR80'
#endif    
            endif
        
        endif 
        
        !Builds horizontal grid object
        call ConstructHorizontalMap(HorizontalMapID      = Me%ObjHorizontalMap,         &
                                    GridDataID           = Me%ObjBathymetry,            &
                                    HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                    Points               = Mask2D,                      &
                                    STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadMap2DFromFile - ModuleField4D - ERR90'


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
         integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB, STAT_CALL 
         logical                                 :: Exist  
         

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
                                                       KLB = KLB-1, KUB = KUB,             &
                                 STAT   = STAT_CALL)                                
                                 
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR10'
            
            call GetHDF5GroupExist (HDF5ID = Me%File%Obj, GroupName = "/Grid/VerticalZ",&
                                    Exist  = Exist, STAT = STAT_CALL)                                
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR15'
            
            if (Exist) then
                                        
                call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
                                  GroupName     = "/Grid/VerticalZ",                        &
                                  Name          = "Vertical_00001",                         &
                                  Array3D       = SZZ,                                      &
                                  STAT          = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR20'
                
            else
                if (KUB==1) then
                    SZZ(:,:,0) = 1.
                    SZZ(:,:,1) = 0.
                else
                    stop 'ReadMap3DFromFile - ModuleField4D - ERR25'
                endif
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
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR30'
#endif                  
        endif


       if      (Me%File%Form == HDF5_  ) then
            
            call HDF5SetLimits  (HDF5ID = Me%File%Obj, ILB = ILB, IUB = IUB,            &
                                                       JLB = JLB, JUB = JUB,            &
                                                       KLB = KLB, KUB = KUB,            &
                                 STAT   = STAT_CALL)                                
                                 
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR40'

            call HDF5ReadData(HDF5ID        = Me%File%Obj,                              &
                              GroupName     = "/Grid",                                  &
                              Name          = trim(Me%File%MaskName),                   &
                              Array3D       = Mask,                                     &
                              STAT          = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR50'
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
                call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                    Array3D         = Mask,                                 &
                                    Name            = trim(Me%File%MaskName),               &
                                    ILB             = ILB,                                  &
                                    IUB             = IUB,                                  &
                                    JLB             = JLB,                                  &
                                    JUB             = JUB,                                  &
                                    KLB             = KLB,                                  &
                                    KUB             = KUB,                                  &
                                    STAT            = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR60'
#endif                
        endif
        
        call ConstructMap   (Map_ID             = Me%ObjMap,                            &
                            GeometryID          = Me%ObjGeometry,                       &
                            HorizontalMapID     = Me%ObjHorizontalMap,                  &
                            WaterPoints3D       = Mask,                                 &
                            STAT                = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR70'
        
        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                    &
                              WaterPoints3D     = Me%ExternalVar%WaterPoints3D, &
                              STAT              = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR80'

        call ComputeInitialGeometry(GeometryID      = Me%ObjGeometry,                   &
                                    WaterPoints3D   = Me%ExternalVar%WaterPoints3D,     &
                                    SZZ             = SZZ,                              &
                                    STAT            = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMap3DFromFile - ModuleField4D - ERR90'

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

        !---------------------------------------------------------------------
        
        call ConstructPropertyID (PropField%ID, Me%ObjEnterData, ExtractType)        
        
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
        

        call GetData(PropField%VGroupPath,                                              &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'VGROUP_PATH',                                      &
                     default      = "/Results",                                         &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR20'

        call GetData(PropField%MultiplyingFactor,                                       &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'MULTIPLYING_FACTOR',                               &
                     default      = 1.,                                                 &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR30'
        
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
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleFiel4D - ERR40'
        
        if (iflag == 1)then
            PropField%HasAddingFactor = .true.
        end if

        call GetData(PropField%FieldName,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FIELD_NAME',                                       &
                     default      = trim(PropField%ID%Name),                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR50'

        call GetData(LastGroupEqualField,                                               &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'LAST_GROUP_EQUAL_FIELD',                           &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR60'

        if (LastGroupEqualField)                                                        &
            PropField%VGroupPath=trim(PropField%VGroupPath)//"/"//trim(PropField%FieldName)


        call GetData(PropField%From2Dto3D,                                              &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'FROM_2D_TO_3D',                                    &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleField4D - ERR70'
        
        
        call GetOutPutTime(Me%ObjEnterData,                                             &
                           CurrentTime      = Me%StartTime,                             &
                           EndTime          = Me%EndTime,                               &
                           keyword          = 'OUTPUT_TIME',                            &
                           SearchType       = ExtractType,                              &
                           OutPutsTime      = Me%OutPut%OutTime,                        &
                           OutPutsOn        = Me%OutPut%Yes,                            &
                           OutPutsNumber    = Me%OutPut%TotalOutputs,                   &
                           STAT             = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadOptions - ModuleField4D - ERR80'        

        Me%OutPut%NextOutPut = 1

        call GetData(PropField%SpaceDim,                                                &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'SPACE_DIM',                                        &
                     default      = Me%MaskDim,                                         &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleField4D - ERR80'

        call GetData(PropField%ChangeInTime,                                            &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = ExtractType,                                        &
                     keyword      = 'CHANGE_IN_TIME',                                   &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleField4D',                                    &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructHDFInput - ModuleField4D - ERR80'
                
        
        ! Check if the simulation goes backward in time or forward in time (default mode)
        call GetBackTracking(Me%ObjTime, Me%BackTracking, STAT = STAT_CALL)                    
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDFInput - ModuleField4D - ERR90' 

    end subroutine ReadOptions

    !--------------------------------------------------------------------------

   !----------------------------------------------------------------------------

    subroutine Generic4thDimension(PropField, ExtractType)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer        :: PropField
        integer                            :: ExtractType

        !Local-----------------------------------------------------------------
        character(len = StringLength)      :: Filename
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
        integer,                intent(IN )             :: ExtractType        
        
        !Local-----------------------------------------------------------------
        type (T_Time)                                   :: AuxTime
        real,    dimension(6)                           :: InitialDate
        real(8), dimension(:), pointer                  :: Instants
        integer                                         :: STAT_CALL, i, NCDF_READ, HDF5_READ, iflag
        logical                                         :: exist


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

            call ConstructHDF5 (Me%File%Obj, trim(Me%File%FileName), HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR30'

            call GetHDF5GroupNumberOfItems(Me%File%Obj, "/Time", &
                                           Me%File%NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR40'
            
            allocate(Me%File%InstantsDates(Me%File%NumberOfInstants))
            
            do i=1, Me%File%NumberOfInstants
                Me%File%InstantsDates(i) = HDF5TimeInstant(i)
            enddo

            Me%File%StartTime = Me%File%InstantsDates(1)
            Me%File%EndTime   = Me%File%InstantsDates(Me%File%NumberOfInstants)
            
            
            Me%File%DefaultNames%bat        = 'Bathymetry'
            Me%File%DefaultNames%lon_stag   = 'Longitude'
            Me%File%DefaultNames%lat_stag   = 'Latitude'
            if (Me%MaskDim == Dim3D) then
                Me%File%DefaultNames%mask   = 'WaterPoints3D'
            else
                Me%File%DefaultNames%mask   = 'WaterPoints2D'
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
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR50'
            
            call NETCDFReadTime(NCDFID = Me%File%Obj, InitialDate = InitialDate,        &
                                nInstants = Me%File%NumberOfInstants,                   &
                                Instants = Instants, STAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR60'
            
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
        if (STAT_CALL /= SUCCESS_) stop 'ConstructFile - ModuleField4D - ERR80'
        
        call GetData(Me%File%LonStagName,                                               &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'LONG_GRID',                                      &
                     default        = Me%File%DefaultNames%lon_stag,                    &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR90'
        
        call GetData(Me%File%LatStagName,                                               &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'LAT_GRID',                                       &
                     default        = Me%File%DefaultNames%lat_stag,                    &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR100'

        call GetData(Me%File%DepthStagName,                                             &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'DEPTH_GRID',                                     &
                     default        = Me%File%DefaultNames%depth_stag,                  &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR110'

        call GetData(Me%File%BathymName,                                                &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'BATHYM_GRID',                                    &
                     default        = Me%File%DefaultNames%bat,                         &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR120'         
        

        call GetData(Me%File%MaskName,                                                  &
                     Me%ObjEnterData,  iflag,                                           &
                     SearchType     = ExtractType,                                      &
                     keyword        = 'MASK_GRID',                                      &
                     default        = Me%File%DefaultNames%mask,                        &
                     ClientModule   = 'ModuleField4D',                                  &
                     STAT           = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructFile - ModuleFiel4D - ERR130'        
 
               
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

            if (.not.Associated(Me%Matrix1D)) then
                allocate(Me%Matrix1D(KLB:KUB))
                Me%Matrix1D(:) = FillValueReal
            endif

            if (.not.Associated(Me%Depth1D)) then
                allocate(Me%Depth1D(KLB:KUB))
                Me%Depth1D(:) = FillValueReal
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
                        write(*,*)'Could not read solution from file'
                        write(*,*) 'Filename =', trim(Me%File%FileName)                    
                        write(*,*)'Could not find second instant in file'
                        write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                        stop      'ConstructPropertyField - ModuleField4D - ERR210'
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
                        write(*,*)'Could not read solution from file'
                        write(*,*) 'Filename =', trim(Me%File%FileName)                    
                        write(*,*)'Could not find second instant in file'
                        write(*,*)'Matrix name: '//trim(NewPropField%FieldName)
                        stop      'ConstructPropertyField - ModuleField4D - ERR220'
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

                    !Interpolates the two matrixes in time
                    call InterpolateMatrix2DInTime(ActualTime       = Me%StartTime,                &
                                                   Size             = Me%WorkSize2D,               &
                                                   Time1            = NewPropField%PreviousTime,   &
                                                   Matrix1          = NewPropField%PreviousField2D,&
                                                   Time2            = NewPropField%NextTime,       &
                                                   Matrix2          = NewPropField%NextField2D,    &
                                                   MatrixOut        = Me%Matrix2D,                 &
                                                   PointsToFill2D   = Me%ExternalVar%WaterPoints2D)
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

                    call InterpolateMatrix3DInTime(ActualTime       = Me%StartTime,                 &
                                                   Size             = Me%WorkSize3D,                &
                                                   Time1            = NewPropField%PreviousTime,    &
                                                   Matrix1          = NewPropField%PreviousField3D, &
                                                   Time2            = NewPropField%NextTime,        &
                                                   Matrix2          = NewPropField%NextField3D,     &
                                                   MatrixOut        = Me%Matrix3D,                  &
                                                   PointsToFill3D   = Me%ExternalVar%WaterPoints3D)

                else

                    !Prev and next are equal (last instant?)
                    Me%Matrix3D(:,:,:)  = NewPropField%NextField3D(:,:,:)

                endif

            end if i4

        endif it

    end subroutine ConstructPropertyField

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


    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        

        !Local-----------------------------------------------------------------
!        type(T_Time)                            :: TimeInstant
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%File%Obj, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%File%Obj,                           &
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

        !Begin-----------------------------------------------------------------
        
        if (Previous) then
            Field   => NewPropField%PreviousField2D
            Instant =  NewPropField%PreviousInstant        
        else
            Field   => NewPropField%NextField2D        
            Instant =  NewPropField%NextInstant
        endif
        
        ILB = Me%WorkSize2D%ILB
        IUB = Me%WorkSize2D%IUB
        JLB = Me%WorkSize2D%JLB
        JUB = Me%WorkSize2D%JUB
        
        if      (Me%File%Form == HDF5_  ) then

            call GetHDF5ArrayDimensions(Me%File%Obj, trim(NewPropField%VGroupPath),         &
                              trim(NewPropField%FieldName), OutputNumber = Instant,         &
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

            call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR30'
            
                 
            call HDF5ReadData(Me%File%Obj, trim(NewPropField%VGroupPath),                   &
                              trim(NewPropField%FieldName),                                 &
                              Array2D = Field, OutputNumber = Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR40'
#ifndef _NO_NETCDF
        else if (Me%File%Form == NetCDF_) then
        
            call NETCDFReadData(NCDFID          = Me%File%Obj,                          &
                                Array2D         = Field,                          &
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
        
        nullify(Field)

    end subroutine ReadValues2D
    
    
    !--------------------------------------------------------------------------

    
    subroutine ReadValues3D (NewPropField, Previous)

        !Arguments-------------------------------------------------------------
        type (T_PropField), pointer             :: NewPropField                
        logical                                 :: Previous
       
        !Local-----------------------------------------------------------------
        integer                                 :: Instant
        real, dimension(:,:,:), pointer         :: Field, Aux3D, FieldAux
        integer                                 :: Imax, Jmax, Kmax
        integer                                 :: STAT_CALL, i, j, k, ILB, IUB, JLB, JUB, KLB, KUB

        !Begin-----------------------------------------------------------------
        
        if (Previous) then
            Field   => NewPropField%PreviousField3D
            Instant =  NewPropField%PreviousInstant        
        else
            Field   => NewPropField%NextField3D        
            Instant =  NewPropField%NextInstant
        endif

         ILB = Me%WorkSize3D%ILB
         IUB = Me%WorkSize3D%IUB
         JLB = Me%WorkSize3D%JLB
         JUB = Me%WorkSize3D%JUB

        if (NewPropField%From2Dto3D) then
            KLB = 1
            KUB = 1
            allocate(Aux3D(Me%Size3D%ILB:Me%Size3D%IUB,Me%Size3D%JLB:Me%Size3D%JUB,1:1))
            
            FieldAux => Aux3D
        else
            KLB = Me%WorkSize3D%KLB
            KUB = Me%WorkSize3D%KUB
            
            FieldAux => Field
        endif
        

        if      (Me%File%Form == HDF5_  ) then

            call GetHDF5ArrayDimensions(Me%File%Obj, trim(NewPropField%VGroupPath),         &
                              trim(NewPropField%FieldName), OutputNumber = Instant,         &
                              Imax = Imax, Jmax = Jmax, Kmax = Kmax, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues2D - ModuleField4D - ERR10'                                   
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
            stop 'ReadValues3D - ModuleField4D - ERR20'                                   

        endif
          


        if      (Me%File%Form == HDF5_  ) then

            call HDF5SetLimits  (Me%File%Obj, ILB, IUB, JLB, JUB, KLB, KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR30'
            
                 
            call HDF5ReadData(Me%File%Obj, trim(NewPropField%VGroupPath),                   &
                              trim(NewPropField%FieldName),                                 &
                              Array3D = FieldAux, OutputNumber = Instant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR40'
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
            if (STAT_CALL /= SUCCESS_)stop 'ReadValues3D - ModuleField4D - ERR50'            
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

        nullify(Field)
    end subroutine ReadValues3D


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
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyField4D(Field4DID, PropertyIDNumber, CurrentTime, Matrix2D, Matrix3D, STAT)

        !Arguments-------------------------------------------------------------
        integer,                 intent(IN)             :: Field4DID
        integer,                 intent(IN)             :: PropertyIDNumber
        type (T_Time),           intent(IN)             :: CurrentTime        
        real,    dimension(:, :),    pointer, optional  :: Matrix2D
        real,    dimension(:, :, :), pointer, optional  :: Matrix3D
        integer,                 intent(OUT), optional  :: STAT

        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL
        logical                                         :: CorrectTimeFrame
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            Me%CurrentTimeExt = CurrentTime
        
            if (Me%BackTracking) then  
                call BacktrackingTime
            else   
                Me%CurrentTimeInt = Me%CurrentTimeExt
            endif
            
            CorrectTimeFrame = .true.
            
            if (Me%CurrentTimeInt < Me%File%StartTime) CorrectTimeFrame = .false.  
            if (Me%CurrentTimeInt > Me%File%EndTime  ) CorrectTimeFrame = .false.              
            
            if (CorrectTimeFrame) then
        
                call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyField4D - ModuleField4D - ERR10'

                if (present(Matrix2D)) Me%Matrix2D => Matrix2D
                if (present(Matrix3D)) Me%Matrix3D => Matrix3D

                if      (PropField%SpaceDim == Dim2D) then
                    call ModifyInput2D (PropField) 
                else if (PropField%SpaceDim == Dim3D) then
                    call ModifyInput3D (PropField)
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


                call HDF5WriteData  (Me%ObjHDF5Out, "/Grid/VerticalZ", "Vertical",      &
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
    subroutine ModifyField4DXYZ(Field4DID, PropertyIDNumber, CurrentTime, X, Y, Z, Field, NoData, STAT)

        !Arguments-------------------------------------------------------------
        integer,                            intent(IN)             :: Field4DID
        integer,                            intent(IN)             :: PropertyIDNumber
        type (T_Time),                      intent(IN)             :: CurrentTime 
        real,    dimension(:),   pointer,   intent(IN)             :: X, Y
        real,    dimension(:),   pointer,   intent(IN),  optional  :: Z
        real,    dimension(:),   pointer,   intent(OUT)            :: Field
        logical, dimension(:),   pointer,   intent(INOUT)          :: NoData
        integer,                            intent(OUT), optional  :: STAT
                                            
        !Local-----------------------------------------------------------------
        type (T_PropField), pointer                     :: PropField
        integer                                         :: STAT_, ready_, STAT_CALL
        logical                                         :: CorrectTimeFrame
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Field4DID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (Me%WindowWithData) then
        
                Me%CurrentTimeExt = CurrentTime
            
                if (Me%BackTracking) then  
                    call BacktrackingTime
                else   
                    Me%CurrentTimeInt = Me%CurrentTimeExt
                endif
                
                CorrectTimeFrame = .true.
                
                if (Me%CurrentTimeInt < Me%File%StartTime) CorrectTimeFrame = .false.  
                if (Me%CurrentTimeInt > Me%File%EndTime  ) CorrectTimeFrame = .false.              
                
                if (CorrectTimeFrame) then        
            
                    call SearchPropertyField(PropField, PropertyIDNumber, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyField4DXYZ - ModuleField4D - ERR10'

                    if      (PropField%SpaceDim == Dim2D) then
                    
                        call Interpolate2DCloud (PropField, X, Y, Field, NoData) 
                        
                    else if (PropField%SpaceDim == Dim3D) then
                    
                        if (.not.present(Z)) then
                            stop 'ModifyField4DXYZ - ModuleField4D - ERR20'
                        endif

                        call Interpolate3DCloud (PropField, X, Y, Z, Field, NoData) 
                    endif
                    
                    if (Me%Output%Yes) then
                        call WriteOutput(PropField, PropertyIDNumber)
                    endif

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
        real,       dimension(:),   pointer, intent(OUT)    :: Field
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real                                            :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        real                                            :: X_W, X_E, Xv, Y_S, Y_N, Yv, PercI, PercJ  
        integer                                         :: STAT_CALL, nPoints, nP
        integer                                         :: jW, jE, iS, iN, i, j
        logical                                         :: InsideDomain

        !Begin----------------------------------------------------------------
        
        nPoints = size(X)
        
        call ModifyInput2D(PropField)
        

        call GetWaterPoints2D(HorizontalMapID   = Me%ObjHorizontalMap,              &
                              WaterPoints2D     = Me%ExternalVar%WaterPoints2D,     &
                              STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'Interpolate2DCloud - ModuleField4D - ERR10' 
        
        if (Me%Extrapolate) then        
            call FillMatrix2DNearestCell(Me%WorkSize2D%ILB,                             &
                                         Me%WorkSize2D%IUB,                             &
                                         Me%WorkSize2D%JLB,                             &
                                         Me%WorkSize2D%JUB,                             &
                                         Me%ExternalVar%Waterpoints2D,                  &
                                         Me%Matrix2D)                    
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
                
                ValueSW = Me%Matrix2D(iS, jW)

                ValueSE = Me%Matrix2D(iS, jE)

                ValueNW = Me%Matrix2D(iN, jW)

                ValueNE     = Me%Matrix2D(iN, jE)
                                               
                ValueN      = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                ValueS      = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                
                Field(nP)   = LinearInterpolation (Y_S,  ValueS, Y_N,  ValueN, Yv)
                
                NoData(nP)  = .false. 
                
            endif
                            
        enddo dnP            

        call UnGetHorizontalMap(HorizontalMapID   = Me%ObjHorizontalMap,                &
                                Array             = Me%ExternalVar%WaterPoints2D,       &
                                STAT              = STAT_CALL) 

        if (STAT_CALL/=SUCCESS_) stop 'Interpolate2DCloud - ModuleField4D - ERR30' 



     end subroutine Interpolate2DCloud     
    !----------------------------------------------------------------------

    !----------------------------------------------------------------------

    subroutine InterpolateBathym (X, Y, Bathym, NoData) 
        
        !Arguments------------------------------------------------------------
        real,       dimension(:),   pointer, intent(IN)     :: X, Y
        real,       dimension(:),   pointer, intent(OUT)    :: Bathym
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
            call FillMatrix2DNearestCell(Me%WorkSize2D%ILB,                             &
                                         Me%WorkSize2D%IUB,                             &
                                         Me%WorkSize2D%JLB,                             &
                                         Me%WorkSize2D%JUB,                             &
                                         Me%ExternalVar%Waterpoints2D,                  &
                                         Me%Matrix2D)                    
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
        real,       dimension(:),   pointer, intent(OUT)    :: Field
        logical,    dimension(:),   pointer, intent(INOUT)  :: NoData
        !Local----------------------------------------------------------------
        real,   dimension(:,:,:),   pointer                 :: SZZ
        real                                                :: ValueSW, ValueNW, ValueSE, ValueNE, ValueN, ValueS
        integer                                             :: MaskSW, MaskNW, MaskSE, MaskNE, MaskN, MaskS
        real                                                :: X_W, X_E, Xv, Y_S, Y_N, Yv, PercI, PercJ  
        integer                                             :: STAT_CALL, nPoints, nP, k
        integer                                             :: jW, jE, iS, iN, i, j
        logical                                             :: InsideDomain

        !Begin----------------------------------------------------------------

        nPoints = size(X)
        
        call ModifyInput3D(PropField)
        
        call GetWaterPoints3D(Map_ID            = Me%ObjMap,                            &
                              WaterPoints3D     = Me%ExternalVar%WaterPoints3D,         &
                              STAT              = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'Interpolate3DCloud - ModuleField4D - ERR10' 

        if (Me%Extrapolate) then        
            call FillMatrix3DNearestCell(Me%WorkSize3D%ILB,                             &
                                         Me%WorkSize3D%IUB,                             &
                                         Me%WorkSize3D%JLB,                             &
                                         Me%WorkSize3D%JUB,                             &
                                         Me%WorkSize3D%KLB,                             &
                                         Me%WorkSize3D%KUB,                             &
                                         Me%ExternalVar%Waterpoints3D,                  &
                                         Me%Matrix3D)                    
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
        
        call UnGetGeometry       (  GeometryID      = Me%ObjGeometry,           &
                                    Array           = SZZ,                      &
                                    STAT            = STAT_CALL)                                     
        if (STAT_CALL /= SUCCESS_) stop 'Interpolate3DCloud - ModuleValida4D - ERR30'

        if (Me%Extrapolate) then        
            call FillMatrix3DNearestCell(Me%WorkSize3D%ILB,                             &
                                         Me%WorkSize3D%IUB,                             &
                                         Me%WorkSize3D%JLB,                             &
                                         Me%WorkSize3D%JUB,                             &
                                         Me%WorkSize3D%KLB,                             &
                                         Me%WorkSize3D%KUB,                             &
                                         Me%ExternalVar%Waterpoints3D,                  &
                                         Me%Depth3D)                    
        endif
               
dnP:    do nP = 1,nPoints      

            if (NoData(nP)) then
                InsideDomain = GetXYInsideDomain(Me%ObjHorizontalGrid, X(nP), Y(nP), STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'Interpolate3DCloud - ModuleValida4D - ERR40'
            
                if (.not. InsideDomain) then
                    cycle
                endif
                
                call GetXYCellZ(Me%ObjHorizontalGrid, X(nP), Y(nP), i, j, PercI, PercJ, STAT = STAT_CALL)
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
                
                k = Me%WorkSize3D%KUB
                
                if (k>1) then
                    ValueSW     = ValueAtDepthZ(iS, jW, Z(nP))
                    ValueSE     = ValueAtDepthZ(iS, jE, Z(nP))
                    ValueNW     = ValueAtDepthZ(iN, jW, Z(nP))
                    ValueNE     = ValueAtDepthZ(iN, jE, Z(nP))
                    
                    MaskSW      = MaskAtDepthZ(iS, jW, Z(nP))        
                    MaskSE      = MaskAtDepthZ(iS, jE, Z(nP))
                    MaskNW      = MaskAtDepthZ(iN, jW, Z(nP))
                    MaskNE      = MaskAtDepthZ(iN, jE, Z(nP))
                                            
                else
                    ValueSW     = Me%Matrix3D(iS, jW, k)
                    ValueSE     = Me%Matrix3D(iS, jE, k)
                    ValueNW     = Me%Matrix3D(iN, jW, k)
                    ValueNE     = Me%Matrix3D(iN, jE, k)
                    
                    MaskSW      = Me%ExternalVar%Waterpoints3D(iS, jW, k)
                    MaskSE      = Me%ExternalVar%Waterpoints3D(iS, jE, k)
                    MaskNW      = Me%ExternalVar%Waterpoints3D(iN, jW, k)
                    MaskNE      = Me%ExternalVar%Waterpoints3D(iN, jE, k)
                    
                endif

                if (ValueSW < FillValueReal/1e4) ValueSW = 0.
                if (ValueSE < FillValueReal/1e4) ValueSE = 0.                
                if (ValueNW < FillValueReal/1e4) ValueNW = 0.                
                if (ValueNE < FillValueReal/1e4) ValueNE = 0.    
                
                if (Me%Extrapolate) then
                
                    ValueN      = LinearInterpolation (X_W, ValueNW, X_E, ValueNE, Xv)
                    ValueS      = LinearInterpolation (X_W, ValueSW, X_E, ValueSE, Xv)
                    
                    Field(nP)   = LinearInterpolation (Y_S, ValueS, Y_N, ValueN, Yv)

                    if (abs(Field(nP)) > 0.) then
                        NoData(nP) = .false. 
                    else
                        NoData(nP) = .true. 
                    endif                    
                    
                else
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

        call UnGetMap(Map_ID          = Me%ObjMap,                                      &
                      Array           = Me%ExternalVar%WaterPoints3D,                   &
                      STAT            = STAT_CALL) 
        if (STAT_CALL/=SUCCESS_) stop 'Interpolate3DCloud - ModuleField4D - ERR60' 

     end subroutine Interpolate3DCloud     
     
    !----------------------------------------------------------------------

    real(8) function ValueAtDepthZ(i, j, Z)
    
        !Arguments-------------------------------------------------------------
        real        :: Z
        integer     :: i, j

        !Local-----------------------------------------------------------------        
        integer     :: KLB, KUB, kb, Ndepths, k

        !Begin-----------------------------------------------------------------            

        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB
        
        do k = KUB,KLB+1,-1
            if (Me%Depth3D (i, j, k-1)<Me%Depth3D (i, j, k)) exit
        enddo
        
        kb = k
        
        Ndepths       = KUB - kb + 1
        
        Me%Depth1D (kb:KUB)  =   Me%Depth3D (i, j, kb:KUB) 
        Me%Matrix1D(kb:KUB)  =   Me%Matrix3D(i, j, kb:KUB)         

        ValueAtDepthZ = InterpolateProfileR8(dble(Z), Ndepths, Me%Depth1D (kb:KUB), Me%Matrix1D(kb:KUB))
        
    end function ValueAtDepthZ


    !----------------------------------------------------------------------

    integer function MaskAtDepthZ(i, j, Z)
    
        !Arguments-------------------------------------------------------------
        real        :: Z
        integer     :: i, j

        !Local-----------------------------------------------------------------        
        integer     :: KLB, KUB, kb, Ndepths, k

        !Begin-----------------------------------------------------------------            

        KLB = Me%WorkSize3D%KLB
        KUB = Me%WorkSize3D%KUB
        
        do k = KUB,KLB+1,-1
            if (Me%Depth3D (i, j, k-1)<Me%Depth3D (i, j, k)) exit
        enddo
        
        kb = k
        
        Ndepths       = KUB - kb + 1
        
        Me%Depth1D (kb:KUB)  =   Me%Depth3D                  (i, j, kb:KUB) 
        Me%Matrix1D(kb:KUB)  =   Me%ExternalVar%WaterPoints3D(i, j, kb:KUB)         

        MaskAtDepthZ = int(InterpolateProfileR8(dble(Z), Ndepths, Me%Depth1D (kb:KUB), Me%Matrix1D(kb:KUB)))
        
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

            call InterpolateMatrix3DInTime(ActualTime       = Me%CurrentTimeInt,                      &
                                           Size             = Me%WorkSize3D,            &
                                           Time1            = PropField%PreviousTime,   &
                                           Matrix1          = PropField%PreviousField3D,&
                                           Time2            = PropField%NextTime,       &
                                           Matrix2          = PropField%NextField3D,    &
                                           MatrixOut        = Me%Matrix3D,              &
                                           PointsToFill3D   = Me%ExternalVar%WaterPoints3D)

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
                call InterpolateMatrix2DInTime(ActualTime       = Me%CurrentTimeInt,                  &
                                               Size             = Me%WorkSize2D,        &
                                               Time1            = PropField%PreviousTime, &
                                               Matrix1          = PropField%PreviousField2D,&
                                               Time2            = PropField%NextTime,     &
                                               Matrix2          = PropField%NextField2D,  &
                                               MatrixOut        = Me%Matrix2D,          &
                                               PointsToFill2D   = Me%ExternalVar%WaterPoints2D)
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

                        if (associated(Me%Matrix1D)) then
                            deallocate(Me%Matrix1D)
                            nullify   (Me%Matrix1D)
                        endif                                        

                        if (associated(Me%Depth1D)) then
                            deallocate(Me%Depth1D)
                            nullify   (Me%Depth1D)
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
