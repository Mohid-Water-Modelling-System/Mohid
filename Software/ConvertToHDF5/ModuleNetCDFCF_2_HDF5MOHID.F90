!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : NetCDFCF_2_HDF5MOHID
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Pablo Carracedo - v4.0
! DESCRIPTION   : Module to convert NetCDFCF_2_HDF5MOHID files into HDF5 format.
!                 For reading into Mohid module HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleNetCDFCF_2_HDF5MOHID

    use ModuleNetCDF
    use ModuleTime
    use ModuleGlobalData
    use ModuleFunctions
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap

    ! Manages NetCDF files
#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif
    implicit none

    private 
    
    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::      ConvertNetCDFCF_2_HDF5MOHID
    private ::      ReadOptions
    private ::      KillNetCDFCF_2_HDF5MOHID

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: input_files_begin   = '<<begin_input_files>>'
    character(LEN = StringLength), parameter    :: input_files_end     = '<<end_input_files>>'

    character(LEN = StringLength), parameter    :: begin_time          = '<<begin_time>>'
    character(LEN = StringLength), parameter    :: end_time            = '<<end_time>>'
    
    character(LEN = StringLength), parameter    :: begin_grid          = '<<begin_grid>>'
    character(LEN = StringLength), parameter    :: end_grid            = '<<end_grid>>'

    character(LEN = StringLength), parameter    :: begin_field         = '<<begin_field>>'
    character(LEN = StringLength), parameter    :: end_field           = '<<end_field>>'


    integer                      , parameter    :: Real4_ = 1, Real8_ =2, Integer4_ = 3

    integer                      , parameter    :: Zonal = 1, Meridional = 2
    
    integer                      , parameter    :: sigma_ = 1, z_level = 2, hybrid = 3


    !Types---------------------------------------------------------------------
    private :: T_WindowOut
    type       T_WindowOut    
        logical     :: ON
        integer     :: ILB, IUB, JLB, JUB
    end type       T_WindowOut 
    
    private :: T_ValueIn
    type       T_ValueIn    
        integer                                       :: DataType = Real8_                
        integer                                       :: Dim
        integer,    dimension(:),       allocatable   :: CountDim
        integer                                       :: diL = 0, diU = 0, djL = 0, djU = 0
        real(4),    dimension(:),       allocatable   :: R41D
        real(4),    dimension(:,:  ),   allocatable   :: R42D        
        real(4),    dimension(:,:,:),   allocatable   :: R43D
        real(4),    dimension(:,:,:,:), allocatable   :: R44D        
        real(8),    dimension(:),       allocatable   :: R81D
        real(8),    dimension(:,:  ),   allocatable   :: R82D        
        real(8),    dimension(:,:,:),   allocatable   :: R83D
        real(8),    dimension(:,:,:,:), allocatable   :: R84D        
        integer(4), dimension(:),       allocatable   :: I41D
        integer(4), dimension(:,:  ),   allocatable   :: I42D        
        integer(4), dimension(:,:,:),   allocatable   :: I43D
        integer(4), dimension(:,:,:,:), allocatable   :: I44D        
    end type T_ValueIn    
    
    private :: T_Date
    type       T_Date
        character(len=StringLength)             :: NetCDFName, NetCDFDimName    
        type (T_ValueIn)                        :: ValueIn        
        real(8), dimension(:), pointer          :: ValueInTotal
        type(T_Time), dimension(:), pointer     :: Value1DOut
        integer                                 :: NumberInst
        integer                                 :: TotalInst        
        integer                                 :: TotalInstOut
        integer, dimension(:), pointer          :: InstOut
        type(T_Time)                            :: RefDateTimeIn, RefDateTimeOut            
        type(T_Time)                            :: FileEndTime, LasTEndTime                    
        real                                    :: UnitsFactor
        logical                                 :: RefAttribute
        character(len=StringLength)             :: RefAttributeName
        character(len=StringLength)             :: RefDateName
        real                                    :: RefDateOffSet
        logical                                 :: RefDateOffSetFromAtt
        character(len=StringLength)             :: RefDateOffSetProp
        character(len=StringLength)             :: RefDateOffSetAtt        
        real                                    :: RefDateOffSetAttFactor        
        integer                                 :: NetCDFvar, NetCDFdim
    end type  T_Date

    private :: T_LongLat
    type       T_LongLat
        character(len=StringLength)             :: NetCDFNameLat
        character(len=StringLength)             :: NetCDFNameLong        
        type (T_ValueIn)                        :: LongIn,  LatIn        
        real, dimension(:,:),     pointer       :: LongOut, LatOut, RotationX, RotationY
        integer                                 :: imax, jmax
        logical                                 :: Imposed
        real                                    :: LongOrig, dLong
        real                                    :: LatOrig,  dLat        
        logical                                 :: Starts180W
        integer                                 :: CorrectJDown, CorrectJUp, BreakJ
        integer                                 :: dij
    end type  T_LongLat


    private :: T_Depth
    type       T_Depth
        character(len=StringLength)             :: NetCDFName, NetCDFDimName
        character(len=StringLength)             :: NetCDFNameFace, NetCDFNameWL
        logical                                 :: Dim3D = .false.
        type (T_ValueIn)                        :: ValueIn        
        type (T_ValueIn)                        :: FaceValueIn
        type (T_ValueIn)                        :: WLValueIn        
        real, dimension(:,:,:),   pointer       :: Value3DOut
        logical                                 :: ON
        integer                                 :: GeoVert !1 - sigma 2 - z-level 3 - hybrid
        logical                                 :: RomsDistortion            = .true.
        real                                    :: theta_s = 0, theta_b = 0, Hc = 0
        character(len=StringLength)             :: positive = "down"
        integer                                 :: kmax
        integer                                 :: kmaxF
        logical                                 :: Interpolate
        integer                                 :: N_ZLevels
        real(8), dimension(:),       pointer    :: Zlevels
        logical                                 :: InvertLayers
        logical                                 :: NetCDFNameFaceOff = .false.
        integer                                 :: RemoveNsurfLayers = 0
        real                                    :: Offset = 0.
    end type  T_Depth

    private :: T_Bathym
    type       T_Bathym
        character(len=StringLength)             :: NetCDFName
        character(len=PathLength)               :: FileName        
        logical                                 :: FromMapping
        logical                                 :: ON = .false.
        real                                    :: Default
        type (T_ValueIn)                        :: ValueIn        
        real, dimension(:,:),     pointer       :: Value2DOut
        logical                                 :: InvertReferential = .false. 
    end type  T_Bathym

    private :: T_Mapping
    type       T_Mapping
        character(len=StringLength)             :: NetCDFName
        logical                                 :: ON           = .true.
        logical                                 :: Dim3D        = .true.
        logical                                 :: From2D_To_3D = .false. 
        real                                    :: Limit
        integer                                 :: Instant
        type (T_ValueIn)                        :: ValueIn        
        integer, dimension(:,:,:),   pointer    :: Value3DOut
        integer, dimension(:,:  ),   pointer    :: Value2DOut
    end type  T_Mapping    

    private :: T_Field
    type       T_Field
        type(T_PropertyID)                      :: ID
        integer                                 :: Dim
        logical                                 :: From2D_To_3D        
        character(len=StringLength)             :: NetCDFName
        real                                    :: Add, Multiply, MinValue
        real                                    :: NewMissingValue, OldMissingValue
        real                                    :: UnitsScale
        real                                    :: UnitsAdd        
        type (T_ValueIn)                        :: ValueIn        
        real, dimension(:,:),     pointer       :: Value2DOut
        real, dimension(:,:,:),   pointer       :: Value3DOut
        logical                                 :: Accumulated2StepGFS
        logical                                 :: Accumulated2Step
        logical                                 :: FromDir2Vector
        logical                                 :: FromMeteo2Algebric
        logical                                 :: FromCartesian2Meteo
        character(len=StringLength)             :: DirX
        character(len=StringLength)             :: DirY
        logical                                 :: ComputeIntensity, Rotation, Beaufort, WaveBeaufort
        logical                                 :: ComputeDirection  
        integer                                 :: VectorComponent
        character(len=StringLength)             :: VectorX
        character(len=StringLength)             :: VectorY        
        logical                                 :: CenterX, CenterY
        logical                                 :: ComputeRH
        character(len=StringLength)             :: TempRH, PressureRH, SpecificHumidityRH
        logical                                 :: ExtractLayer
        integer                                 :: LayerNumber
        integer                                 :: LayerDim       
        logical                                 :: AverageInDepth, Wfp
        character(len=StringLength)             :: AverageInDepthName, WfpName
        logical                                 :: Reflectivity2Precipitation
        character(len=StringLength)             :: ReflectivityName
        integer                                 :: DirectionReferential
        real                                    :: Limit
        logical                                 :: Energy2Power
        character(len=StringLength)             :: EnergyName
        integer, dimension(:), pointer          :: EnergyStartingHours        
        integer                                 :: N_EnergyStartH
        logical                                 :: AvModelStart2Inst
        character(len=StringLength)             :: AvModelStartName
        !integer, dimension(:), pointer          :: AvModelStartingHours        
        !integer                                 :: N_AvModelStartH
        logical                                 :: VerticalZ_2D
        character(len=StringLength)             :: Surface
        character(len=StringLength)             :: Bottom        
        logical                                 :: CheckMinMaxLimits
        real                                    :: MinLimit
        real                                    :: MaxLimit
        real                                    :: GridRotation
        logical                                 :: ComputeRotatedVector
    end type  T_Field

    type T_NetCDF_Out                                         
        character(len=PathLength)                           :: Name
        integer                                             :: ObjNETCDF    = 0
        character(len=StringLength)                         :: Title
        character(len=StringLength)                         :: Convention
        character(len=StringLength)                         :: Version
        character(len=StringLength)                         :: History
        character(len=StringLength)                         :: Source
        character(len=StringLength)                         :: Institution
        character(len=StringLength)                         :: References
        integer                                             :: iDate
    end type T_NetCDF_Out

    
    private :: T_NetCDFCF_2_HDF5MOHID
    type       T_NetCDFCF_2_HDF5MOHID
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit, ClientNumber
        integer                                 :: FieldInType
        character(len=PathLength)               :: FileName
        logical                                 :: OnlyOneFile          = .false.
        character(len=PathLength)               :: ReadingFileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: GeometryFileName
        type(T_NetCDF_Out)                      :: NetCDF_Out
        type(T_Size3D)                          :: Size, WorkSize
        type(T_Size2D)                          :: Size2D, WorkSize2D
        type(T_Field),  dimension(:), allocatable :: Field 
        integer                                 :: PropNumber, ReadPropNumber
        type(T_Date)                            :: Date  
        type(T_LongLat)                         :: LongLat
        type(T_Depth)                           :: Depth
        type(T_Mapping)                         :: Mapping
        type(T_Bathym)                          :: Bathym
        logical                                 :: OutHDF5, OutNetcdf, ReadInvertXY
        logical                                 :: ReadInvertLat
        logical                                 :: Nan_2_Null
        integer                                 :: OutCountProp = 0
        type(T_WindowOut)                       :: WindowOut
        logical                                 :: MeridionalSplit
        integer                                 :: MeridionalSplitColumn
    end type  T_NetCDFCF_2_HDF5MOHID

    type(T_NetCDFCF_2_HDF5MOHID), pointer             :: Me


    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertNetCDFCF_2_HDF5MOHID(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
        integer                                         :: iP, HDF5_READWRITE, STAT_CALL
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadOptions

        call ReadNetCDFCF_WriteHDF5MOHID
        
        if (associated(Me%LongLat%LongOut    )) deallocate(Me%LongLat%LongOut   )
        if (associated(Me%LongLat%LatOut     )) deallocate(Me%LongLat%LatOut    )
        
        if (associated(Me%Bathym%Value2DOut  )) deallocate(Me%Bathym%Value2DOut )
        
        do ip = 1, Me%PropNumber   
            if (associated(Me%Field(iP)%Value2DOut)) deallocate(Me%Field(iP)%Value2DOut)
            if (associated(Me%Field(iP)%Value3DOut)) deallocate(Me%Field(iP)%Value3DOut)
        enddo

        if (Me%OutHDF5) then
            !Close HDF5 File
            call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConvertNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

            !Gets File Access Code
            call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
        
            !Opens HDF5 File
            call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_READWRITE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConvertNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
        endif
        

        call WriteRotatedVector
          
        call WriteRotation

        call WriteComputeDirection

        call WriteComputeIntensity
        
        call WriteBeaufort 
        
        call WriteRelativeHumidity
        
        call WriteAverageInDepth
        
        call WriteReflectivity2Precipitation
        
        call WriteEnergy2Power
        
        call WriteAvModelStart2Inst
        
        call WriteVerticalZ_2D
       
        call WriteWaveFpToWaveTp
       
        if (Me%OutNetCDF) call WriteTimeNetCDF(DefDimTime=.false.)
    
        call KillNetCDFCF_2_HDF5MOHID

        STAT = SUCCESS_

    end subroutine ConvertNetCDFCF_2_HDF5MOHID

    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
    
    subroutine WriteWaveFpToWaveTp
        !Local-----------------------------------------------------------------
        real,   dimension(:,:  ), pointer           :: Prop2D
        real,   dimension(:,:  ), pointer           :: TpProp
        integer                                     :: i, iP, iPt
        logical                                     :: Found
        !Begin-----------------------------------------------------------------

        nullify(TpProp)
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%Wfp) then
        
                !Found property to calculated
                Found = .false.
                do iPt = 1, Me%PropNumber
                    if (trim(Me%Field(iPt)%ID%Name)==trim(Me%Field(iP)%WfpName)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteWaveFpToWaveTp - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


                do i=1, Me%Date%TotalInstOut

                    !Read 2D property

                    if      (Me%Field(iPt)%Dim/=2) then
                        stop 'WriteWaveFpToWaveTp - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                    endif
                    
                    allocate(Me%Field(iPt)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                    Prop2D => Me%Field(iPt)%Value2DOut                 
                    
                    call ReadFieldHDF5(iPt, i)
                    
                    !Allocate new 2D property
                    
                    allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                    if      (Me%Field(iP)%Dim/=2) then
                        stop 'WriteWaveFpToWaveTp - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                    endif
                    
                    allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                    TpProp => Me%Field(iP)%Value2DOut                    
                    
                    !Compute peak frequency to peak period
                    write(*,*) '  Calculating= ', trim(Me%Field(iP)%ID%Name), i
                    call ComputePeakPeriod(Prop2D,TpProp)

                    !Write TpProp 
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif

                    deallocate(Me%Field(iP)%Value2DOut)

                enddo
            endif                
        enddo    
        
    end subroutine WriteWaveFpToWaveTp

    !------------------------------------------------------------------------
    
    subroutine WriteComputeIntensity
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, iP, iPx, iPy
        logical                                     :: Found
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%ComputeIntensity) then
        
                !Found component X
                Found = .false.
                do iPx = 1, Me%PropNumber
                    if (trim(Me%Field(iPx)%ID%Name)==trim(Me%Field(iP)%VectorX)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                !Found component Y
                Found = .false.
                do iPy = 1, Me%PropNumber
                    if (trim(Me%Field(iPy)%ID%Name)==trim(Me%Field(iP)%VectorY)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


                do i=1, Me%Date%TotalInstOut
                    !Read component X

                    if      (Me%Field(iPx)%Dim==2) then
                        allocate(Me%Field(iPx)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPx)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    if      (Me%Field(iP)%Dim==2) then
                        allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB,           &
                                                         Me%Size%KLB:Me%Size%KUB))                    
                    endif

                                      
                    
                    call ReadFieldHDF5(iPx, i)

                    !Compute vector intensity
                    call ComputeVectorIntensity(iPx=iPx, Step=1, iP=iP)   
                    
                    if      (Me%Field(iPx)%Dim==2) then
                        deallocate(Me%Field(iPx)%Value2DOut)
                    else
                        deallocate(Me%Field(iPx)%Value3DOut)                    
                    endif                    
                    
                    if      (Me%Field(iPy)%Dim==2) then
                        allocate(Me%Field(iPy)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPy)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                                 
                    
                    !Read component Y                    
                    call ReadFieldHDF5(iPy, i)

                    !Compute vector intensity
                    call ComputeVectorIntensity(iPx=iPy, Step=2, iP=iP)              

                    if      (Me%Field(iPy)%Dim==2) then
                        deallocate(Me%Field(iPy)%Value2DOut)
                    else
                        deallocate(Me%Field(iPy)%Value3DOut)                    
                    endif
               
                    !Write Intensity
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                    
                    
                    if      (Me%Field(iP)%Dim==2) then
                        deallocate(Me%Field(iP)%Value2DOut)
                    else
                        deallocate(Me%Field(iP)%Value3DOut)
                    endif

                enddo
            endif                
        enddo    
        
    end subroutine WriteComputeIntensity
    
    !------------------------------------------------------------------------    

    subroutine WriteComputeDirection
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, iP, iPx, iPy
        logical                                     :: Found
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%ComputeDirection) then
        
                !Found component X
                Found = .false.
                do iPx = 1, Me%PropNumber
                    if (trim(Me%Field(iPx)%ID%Name)==trim(Me%Field(iP)%VectorX)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeDirection - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                !Found component Y
                Found = .false.
                do iPy = 1, Me%PropNumber
                    if (trim(Me%Field(iPy)%ID%Name)==trim(Me%Field(iP)%VectorY)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeDirection - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


                do i=1, Me%Date%TotalInstOut
                    !Read component X

                    if      (Me%Field(iPx)%Dim==2) then
                        allocate(Me%Field(iPx)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPx)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    if      (Me%Field(iP)%Dim==2) then
                        allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB,           &
                                                         Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    call ReadFieldHDF5(iPx, i)

                    !Compute vector direction
                    call ComputeVectorDirection(iPx=iPx, Step=1, iP=iP)   
                    
                    if      (Me%Field(iPx)%Dim==2) then
                        deallocate(Me%Field(iPx)%Value2DOut)
                    else
                        deallocate(Me%Field(iPx)%Value3DOut)                    
                    endif                    
                    
                    if      (Me%Field(iPy)%Dim==2) then
                        allocate(Me%Field(iPy)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPy)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                                 
                    
                    !Read component Y                    
                    call ReadFieldHDF5(iPy, i)

                    !Compute vector direction
                    call ComputeVectorDirection(iPx=iPy, Step=2, iP=iP)              

                    if      (Me%Field(iPy)%Dim==2) then
                        deallocate(Me%Field(iPy)%Value2DOut)
                    else
                        deallocate(Me%Field(iPy)%Value3DOut)                    
                    endif
               
                    !Write wind direction
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                    
                    
                    if      (Me%Field(iP)%Dim==2) then
                        deallocate(Me%Field(iP)%Value2DOut)
                    else
                        deallocate(Me%Field(iP)%Value3DOut)
                    endif

                enddo
            endif                
        enddo    
        
    end subroutine WriteComputeDirection
    
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------    

    subroutine WriteRotatedVector
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, iP, iPx, iPy
        logical                                     :: Found
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%ComputeRotatedVector) then
        
                !Found component X
                Found = .false.
                do iPx = 1, Me%PropNumber
                    if (trim(Me%Field(iPx)%ID%Name)==trim(Me%Field(iP)%VectorX)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeDirection - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                !Found component Y
                Found = .false.
                do iPy = 1, Me%PropNumber
                    if (trim(Me%Field(iPy)%ID%Name)==trim(Me%Field(iP)%VectorY)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeDirection - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


                do i=1, Me%Date%TotalInstOut
                    !Read component X

                    if      (Me%Field(iPx)%Dim==2) then
                        allocate(Me%Field(iPx)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPx)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    if      (Me%Field(iP)%Dim==2) then
                        allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB,           &
                                                         Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    call ReadFieldHDF5(iPx, i)

                    !Compute vector direction
                    call ModifyRotatedVector(iPx=iPx, Step=1, iP=iP)   
                    
                    if      (Me%Field(iPx)%Dim==2) then
                        deallocate(Me%Field(iPx)%Value2DOut)
                    else
                        deallocate(Me%Field(iPx)%Value3DOut)                    
                    endif                    
                    
                    if      (Me%Field(iPy)%Dim==2) then
                        allocate(Me%Field(iPy)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPy)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                                 
                    
                    !Read component Y                    
                    call ReadFieldHDF5(iPy, i)

                    !Compute vector direction
                    call ModifyRotatedVector(iPx=iPy, Step=2, iP=iP)              

                    if      (Me%Field(iPy)%Dim==2) then
                        deallocate(Me%Field(iPy)%Value2DOut)
                    else
                        deallocate(Me%Field(iPy)%Value3DOut)                    
                    endif
               
                    !Write wind direction
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                    
                    
                    if      (Me%Field(iP)%Dim==2) then
                        deallocate(Me%Field(iP)%Value2DOut)
                    else
                        deallocate(Me%Field(iP)%Value3DOut)
                    endif

                enddo
            endif                
        enddo    
        
    end subroutine WriteRotatedVector
    
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------    

    subroutine WriteVerticalZ_2D
    
        !Local-----------------------------------------------------------------
        integer                                     :: inst, i, j, iP, iPb, iPs
        integer                                     :: STAT_CALL
        logical                                     :: Found
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%VerticalZ_2D) then
        
                !Found Bottom
                Found = .false.
                do iPb = 1, Me%PropNumber
                    if (trim(Me%Field(iPb)%ID%Name)==trim(Me%Field(iP)%Bottom)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteVerticalZ_2D - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                !Found Surface
                Found = .false.
                do iPs = 1, Me%PropNumber
                    if (trim(Me%Field(iPs)%ID%Name)==trim(Me%Field(iP)%Surface)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteVerticalZ_2D - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


                do inst=1, Me%Date%TotalInstOut
                    !Read Bottom

                    allocate(Me%Field(iPb)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,0:2))
                    
                    call ReadFieldHDF5(iPb, inst)

                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        if (Me%Field(iPb)%Value2DOut(i, j) > -500) then
                            Me%Field(iP)%Value3DOut(i,j,0) = Me%Field(iPb)%Value2DOut(i, j)
                        else
                            Me%Field(iP)%Value3DOut(i,j,0) = 0
                        endif                            
                    enddo
                    enddo
                    

                    
                    allocate(Me%Field(iPs)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                                 
                    
                    !Read component Y                    
                    call ReadFieldHDF5(iPs, inst)
                    
                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                    do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                        if (Me%Field(iPs)%Value2DOut(i, j) > -500) then
                            Me%Field(iP)%Value3DOut(i,j,1) = - Me%Field(iPs)%Value2DOut(i, j)
                        else
                            Me%Field(iP)%Value3DOut(i,j,1) =   Me%Field(iPb)%Value2DOut(i, j)
                        endif                            
                    enddo
                    enddo           

                    !Write Vertical Z 2D
                    if (Me%OutHDF5) then
                        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,  &
                                                         Me%WorkSize%JLB, Me%WorkSize%JUB,  &
                                                         0, 1,                              &
                                                         STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'WriteVerticalZ_2D - ModuleNetCDFCF_2_HDF5MOHID - ERR30' 

                        call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",         &
                                             "m", Array3D = Me%Field(iP)%Value3DOut, OutputNumber = inst, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'WriteVerticalZ_2D - ModuleNetCDFCF_2_HDF5MOHID - ERR40'  
                        
                    !    if (inst == 1) then
                    !        
                    !        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,  &
                    !                                         Me%WorkSize%JLB, Me%WorkSize%JUB,  &
                    !                                         STAT = STAT_CALL)
                    !        if (STAT_CALL /= SUCCESS_) stop 'WriteVerticalZ_2D - ModuleNetCDFCF_2_HDF5MOHID - ERR30' 
                    !
                    !        call HDF5WriteData  (Me%ObjHDF5, "/Grid", "Bathymetry",         &
                    !                             "m", Array2D = Me%Field(iPb)%Value2DOut, STAT = STAT_CALL)
                    !        if (STAT_CALL /= SUCCESS_) stop 'WriteVerticalZ_2D - ModuleNetCDFCF_2_HDF5MOHID - ERR40'  
                    !    endif                                                    
                    endif
                    
                    deallocate(Me%Field(iPs)%Value2DOut)
                    deallocate(Me%Field(iPb)%Value2DOut)
                    deallocate(Me%Field(iP)%Value3DOut)

                enddo
            endif                
        enddo    
        
    end subroutine WriteVerticalZ_2D
    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    
    subroutine WriteBeaufort
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, iP, iPx
        logical                                     :: Found
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%Beaufort .or. Me%Field(iP)%WaveBeaufort) then
        
                !Found component X
                Found = .false.
                do iPx = 1, Me%PropNumber
                    if (trim(Me%Field(iPx)%ID%Name)==trim(Me%Field(iP)%VectorX)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteBeaufort - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                do i=1, Me%Date%TotalInstOut
                    !Read component X
                    if      (Me%Field(iPx)%Dim==2) then
                        allocate(Me%Field(iPx)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPx)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    if      (Me%Field(iP)%Dim==2) then
                        allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB,           &
                                                         Me%Size%KLB:Me%Size%KUB))                    
                    endif
                                                        
                    call ReadFieldHDF5(iPx, i)

                    !Compute Beaufort scale
                    call ComputeBeaufort(iPx=iPx,  iP=iP)   
                    
                    if      (Me%Field(iPx)%Dim==2) then
                        deallocate(Me%Field(iPx)%Value2DOut)
                    else
                        deallocate(Me%Field(iPx)%Value3DOut)                    
                    endif                    
                    
 
                    !Write Intensity
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                    
                    
                    if      (Me%Field(iP)%Dim==2) then
                        deallocate(Me%Field(iP)%Value2DOut)
                    else
                        deallocate(Me%Field(iP)%Value3DOut)
                    endif

                enddo
            endif                
        enddo    
        
    end subroutine WriteBeaufort
    
    !------------------------------------------------------------------------   
    !------------------------------------------------------------------------
    
    subroutine WriteRelativeHumidity
        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: T, P, SH, RH
        integer                                     :: i, iP, iPt, iPp, iPsh
        logical                                     :: Found
        !Begin-----------------------------------------------------------------

        nullify(T, P, SH, RH)
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%ComputeRH) then
        
                !Found temperature
                Found = .false.
                do iPt = 1, Me%PropNumber
                    if (trim(Me%Field(iPt)%ID%Name)==trim(Me%Field(iP)%TempRH)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                !Found pressure
                Found = .false.
                do iPp = 1, Me%PropNumber
                    if (trim(Me%Field(iPp)%ID%Name)==trim(Me%Field(iP)%PressureRH)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

                !Found specific humidity
                Found = .false.
                do iPsh = 1, Me%PropNumber
                    if (trim(Me%Field(iPsh)%ID%Name)==trim(Me%Field(iP)%SpecificHumidityRH)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

                do i=1, Me%Date%TotalInstOut
                    !Read temperature

                    if      (Me%Field(iPt)%Dim/=2) then
                        stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                    endif

                    allocate(Me%Field(iPt)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                    T => Me%Field(iPt)%Value2DOut

                    if      (Me%Field(iP)%Dim/=2) then
                        stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                    endif
                    
                    allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                    RH => Me%Field(iP)%Value2DOut                    
                    
                    call ReadFieldHDF5(iPt, i)
                    
                    !Compute relative humidity
                    call ComputeRelativeHumidity(T, P, SH, RH, step = 1) 
                    
                    deallocate(Me%Field(iPt)%Value2DOut)
                    nullify(Me%Field(iPt)%Value2DOut, T)
                    
                    if      (Me%Field(iPp)%Dim/=2) then
                        stop 'WriteComputeIntensity - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
                    endif

                    allocate(Me%Field(iPp)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                    P => Me%Field(iPp)%Value2DOut

                                 
                    
                    !Read pressure               
                    call ReadFieldHDF5(iPp, i)

                    !Compute RelativeHumidity
                    call ComputeRelativeHumidity(T, P, SH, RH, step = 2) 

                    deallocate(Me%Field(iPp)%Value2DOut)
                    nullify(Me%Field(iPp)%Value2DOut, P)
                    

                    allocate(Me%Field(iPsh)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                    SH => Me%Field(iPsh)%Value2DOut

                    !Read pressure               
                    call ReadFieldHDF5(iPsh, i)

                    !Compute RelativeHumidity
                    call ComputeRelativeHumidity(T, P, SH, RH, step = 3) 
                    
                    deallocate(Me%Field(iPsh)%Value2DOut)
                    nullify   (Me%Field(iPsh)%Value2DOut, SH)
               
                    !Write Intensity
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif

                    deallocate(Me%Field(iP)%Value2DOut)

                enddo
            endif                
        enddo    
        
    end subroutine WriteRelativeHumidity
    
    !------------------------------------------------------------------------       

    !------------------------------------------------------------------------
    
    subroutine WriteAverageInDepth
        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:), pointer           :: Prop3D
        real,   dimension(:,:  ), pointer           :: AvProp
        integer                                     :: i, iP, iPt
        logical                                     :: Found
        !Begin-----------------------------------------------------------------

        nullify(Prop3D, AvProp)
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%AverageInDepth) then
        
                !Found property to be average in depth
                Found = .false.
                do iPt = 1, Me%PropNumber
                    if (trim(Me%Field(iPt)%ID%Name)==trim(Me%Field(iP)%AverageInDepthName)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteAverageInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


                do i=1, Me%Date%TotalInstOut
                    !Read 3D property

                    if      (Me%Field(iPt)%Dim/=3) then
                        stop 'WriteAverageInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                    endif

                    allocate(Me%Field(iPt)%Value3DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))

                    Prop3D => Me%Field(iPt)%Value3DOut

                    if      (Me%Field(iP)%Dim/=2) then
                        stop 'WriteAverageInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                    endif
                    
                    allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                    AvProp => Me%Field(iP)%Value2DOut                    
                    
                    call ReadFieldHDF5(iPt, i)
                    
                    call ReadVerticalZHDF5(i)
                    
                    !Compute prop3D average in depth
                    call ComputeAverageInDepth(Prop3D, AvProp) 
                    
                    deallocate(Me%Field(iPt)%Value3DOut)
                    nullify   (Me%Field(iPt)%Value3DOut, Prop3D)
               
                    deallocate(Me%Depth%Value3DOut)
                    nullify   (Me%Depth%Value3DOut)


                    !Write prop3D average in depth
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif

                    deallocate(Me%Field(iP)%Value2DOut)

                enddo
            endif                
        enddo    
        
    end subroutine WriteAverageInDepth
    
    !------------------------------------------------------------------------       


    !------------------------------------------------------------------------
    
    subroutine WriteReflectivity2Precipitation
        !Local-----------------------------------------------------------------
        real,   dimension(:,:  ), pointer           :: Reflectivity, Precipitation
        integer                                     :: i, iP, iPt
        logical                                     :: Found
        !Begin-----------------------------------------------------------------

        nullify(Reflectivity, Precipitation)
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%Reflectivity2Precipitation) then
        
                !Found property to be average in depth
                Found = .false.
                do iPt = 1, Me%PropNumber
                    if (trim(Me%Field(iPt)%ID%Name)==trim(Me%Field(iP)%ReflectivityName)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteReflectivity2Precipitation - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


                do i=1, Me%Date%TotalInstOut
                    !Read 2D property

                    if      (Me%Field(iPt)%Dim/=2) then
                        stop 'WriteReflectivity2Precipitation - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                    endif

                    allocate(Me%Field(iPt)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                    Reflectivity => Me%Field(iPt)%Value2DOut

                    if      (Me%Field(iP)%Dim/=2) then
                        stop 'WriteAverageInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                    endif
                    
                    allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                    Precipitation => Me%Field(iP)%Value2DOut                    

                    if (i > 1) then                    

                        call ReadFieldHDF5(iPt, i-1)
                       
                    endif                       
                    
                    !Compute Precipitation from Reflectivity
                    call ComputePrecipitation(Reflectivity, Precipitation, i-1, Step1 = .true.) 
                    
                    if (i > 1) then

                        call ReadFieldHDF5(iPt, i)
                        
                        !Compute Precipitation from Reflectivity
                        call ComputePrecipitation(Reflectivity, Precipitation, i, Step1 = .false.) 
                    
                    endif
                    
                    deallocate(Me%Field(iPt)%Value2DOut)
                    nullify   (Me%Field(iPt)%Value2DOut, Reflectivity)
               

                    !Write Precipitation
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif

                    deallocate(Me%Field(iP)%Value2DOut)

                enddo
            endif                
        enddo    
        
    end subroutine WriteReflectivity2Precipitation
    
    !------------------------------------------------------------------------       

    !------------------------------------------------------------------------
    
    subroutine WriteEnergy2Power
        !Local-----------------------------------------------------------------
        integer                                     :: i, iP, iPt
        logical                                     :: Found
        !Begin-----------------------------------------------------------------

        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%Energy2Power) then
        
                !Found property to be average in depth
                Found = .false.
                do iPt = 1, Me%PropNumber
                    if (trim(Me%Field(iPt)%ID%Name)==trim(Me%Field(iP)%EnergyName)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteEnergy2Power - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


                if      (Me%Field(iPt)%Dim/=2) then
                    stop 'WriteEnergy2Power - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                endif

                allocate(Me%Field(iPt)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                if      (Me%Field(iP)%Dim/=2) then
                    stop 'WriteAverageInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                endif
                    
                allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                !Compute Power from Energy

                do i=1, Me%Date%TotalInstOut
                !Read 2D property           
                    if      (i == 1) then
                        
                        call ComputePowerFromEnergy(iPt, iP, 1, Option = "Start")    
                        
                    elseif  (i == Me%Date%TotalInstOut) then
                        
                        call ComputePowerFromEnergy(iPt, iP, Me%Date%TotalInstOut, Option = "End")                         
                        
                    elseif (ModelStartHour(Me%Field(ip)%EnergyStartingHours, Me%Field(ip)%N_EnergyStartH, i)) then
                        !Compute Power from Energy
                        call ComputePowerFromEnergy(iPt, iP, i, Option = "ReStart") 
                    else
                        !Compute Power from Energy
                        call ComputePowerFromEnergy(iPt, iP, i, Option = "Normal") 
                    endif
                    

                    !Write Power
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                        
                enddo
            
                deallocate(Me%Field(iPt)%Value2DOut)
                nullify   (Me%Field(iPt)%Value2DOut)
                deallocate(Me%Field(iP )%Value2DOut)
                nullify   (Me%Field(iP )%Value2DOut)
                
            endif                
        enddo    
        
    end subroutine WriteEnergy2Power
    
    !------------------------------------------------------------------------
    
    subroutine WriteAvModelStart2Inst
        !Local-----------------------------------------------------------------
        integer                                     :: i, iP, iPt
        logical                                     :: Found
        !Begin-----------------------------------------------------------------

        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%AvModelStart2Inst) then
        
                !Found property to be average in depth
                Found = .false.
                do iPt = 1, Me%PropNumber
                    if (trim(Me%Field(iPt)%ID%Name)==trim(Me%Field(iP)%AvModelStartName)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteAvModelStart2Inst - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


                if      (Me%Field(iPt)%Dim/=2) then
                    stop 'WriteAvModelStart2Inst - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                endif

                allocate(Me%Field(iPt)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

                if      (Me%Field(iP)%Dim/=2) then
                    stop 'WriteAverageInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                endif
                    
                allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    
                !Compute Instant value from average since model start

                do i=1, Me%Date%TotalInstOut
                !Read 2D property           
                    if      (i == 1) then
                        
                        call ComputeInstFromAvgSinceModelStart(iPt, iP, 1, Option = "Start")    
                        
                    elseif  (i == Me%Date%TotalInstOut) then
                        
                        call ComputeInstFromAvgSinceModelStart(iPt, iP,&
                                                               Me%Date%TotalInstOut, Option = "End")                         
                        
                    else
                        !Compute Power from Energy
                        call ComputeInstFromAvgSinceModelStart(iPt, iP, i, Option = "Normal") 
                    endif
                    

                    !Write Power
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                        
                enddo
            
                deallocate(Me%Field(iPt)%Value2DOut)
                nullify   (Me%Field(iPt)%Value2DOut)
                deallocate(Me%Field(iP )%Value2DOut)
                nullify   (Me%Field(iP )%Value2DOut)
                
            endif                
        enddo    
        
    end subroutine WriteAvModelStart2Inst
    
    
    !------------------------------------------------------------------------       
    logical function ModelStartHour(StartingHours, nTotal, i)

        !Arguments-------------------------------------------------------------
        integer, dimension(:), pointer              :: StartingHours
        integer                                     :: nTotal, i

        !Local-----------------------------------------------------------------
        real                                        :: hour
        integer                                     :: n 
        logical                                     :: Aux
        !Begin-----------------------------------------------------------------

        Aux = .false.
        
        do n=1, nTotal
            call ExtractDate(Time1 = HDF5TimeInstant(i), hour = hour) 
            if (int(hour) == StartingHours(n)) then
                Aux = .true.
                exit
            endif
        enddo            
                
        ModelStartHour = Aux
        
    end function ModelStartHour

    !------------------------------------------------------------------------     

    subroutine WriteRotation
    
        !Local-----------------------------------------------------------------
        real,  dimension(:,:), pointer              :: RotationX, RotationY
        integer                                     :: i, iP, iPx, iPy
        logical                                     :: Found
        integer                                     :: WorkJLB, WorkJUB, WorkILB, WorkIUB
        integer                                     :: ii, jj, STAT_CALL
        !Begin-----------------------------------------------------------------
        
        do iP = 1, Me%PropNumber
        
            if (Me%Field(iP)%Rotation) then
                    
                if (.not.associated(Me%LongLat%RotationX)) then
                
                
                    if (Me%WindowOut%ON) then
                        WorkILB = 1
                        WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
                        WorkJLB = 1
                        WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
                    else
                    
                        WorkILB = Me%WorkSize%ILB
                        WorkIUB = Me%WorkSize%IUB        
                        WorkJLB = Me%WorkSize%JLB
                        WorkJUB = Me%WorkSize%JUB       
                        
                        Me%WindowOut%ILB = Me%WorkSize%ILB
                        Me%WindowOut%IUB = Me%WorkSize%IUB
                        Me%WindowOut%JLB = Me%WorkSize%JLB
                        Me%WindowOut%JUB = Me%WorkSize%JUB

                    endif
                                 
                
                    allocate(Me%LongLat%RotationX (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%LongLat%RotationY (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))            
                    

                    call GetGridRotation(Me%ObjHorizontalGrid, RotationX, RotationY, STAT = STAT_CALL)   
                    if (STAT_CALL /= SUCCESS_)stop 'WriteRotation - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
                    
                    do jj = WorkJLB, WorkJUB
                    do ii = WorkILB, WorkIUB
                        Me%LongLat%RotationX(ii+ Me%WindowOut%ILB - 1,jj+ Me%WindowOut%JLB - 1) = RotationX(ii,jj)
                        Me%LongLat%RotationY(ii+ Me%WindowOut%ILB - 1,jj+ Me%WindowOut%JLB - 1) = RotationY(ii,jj)                
                    enddo
                    enddo
        
                    call UnGetHorizontalGrid(Me%ObjHorizontalGrid, RotationX, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'WriteRotation - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        
                    call UnGetHorizontalGrid(Me%ObjHorizontalGrid, RotationY, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'WriteRotation - ModuleNetCDFCF_2_HDF5MOHID - ERR30'                    
                endif                                        
        
                !Found component X
                Found = .false.
                do iPx = 1, Me%PropNumber
                    if (trim(Me%Field(iPx)%ID%Name)==trim(Me%Field(iP)%VectorX)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteRotation - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

                !Found component Y
                Found = .false.
                do iPy = 1, Me%PropNumber
                    if (trim(Me%Field(iPy)%ID%Name)==trim(Me%Field(iP)%VectorY)) then
                        Found = .true.
                        exit
                    endif
                enddo
                if (.not. Found) stop 'WriteRotation - ModuleNetCDFCF_2_HDF5MOHID - ERR50'


                do i=1, Me%Date%TotalInstOut
                    !Read component X

                    if      (Me%Field(iPx)%Dim==2) then
                        allocate(Me%Field(iPx)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPx)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    if      (Me%Field(iP)%Dim==2) then
                        allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                         Me%Size%JLB:Me%Size%JUB,           &
                                                         Me%Size%KLB:Me%Size%KUB))                    
                    endif

                    call ReadFieldHDF5(iPx, i)
                    
                    !Compute vector rotation
                    call ComputeVectorRotation(iPx=iPx, Step=1, iP=iP, Component = Me%Field(iP)%VectorComponent)              

  
                      if      (Me%Field(iPx)%Dim==2) then
                        deallocate(Me%Field(iPx)%Value2DOut)
                    else
                        deallocate(Me%Field(iPx)%Value3DOut)                    
                    endif

                    if      (Me%Field(iPy)%Dim==2) then
                        allocate(Me%Field(iPy)%Value2DOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
                    else
                        allocate(Me%Field(iPy)%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                                          Me%Size%JLB:Me%Size%JUB,           &
                                                          Me%Size%KLB:Me%Size%KUB))                    
                    endif
                    
                    
                    !Read component Y                    
                    call ReadFieldHDF5(iPy, i)


                    !Compute vector intensity
                    call ComputeVectorRotation(iPx=iPy, Step=2, iP=iP, Component = Me%Field(iP)%VectorComponent)              


                    if      (Me%Field(iPy)%Dim==2) then
                        deallocate(Me%Field(iPy)%Value2DOut)
                    else
                        deallocate(Me%Field(iPy)%Value3DOut)                    
                    endif
                    
                    
               
                    !Write vector rotation
                    if (Me%OutHDF5) then
                        call WriteFieldHDF5  (iP, i)    
                    endif

                    if (Me%OutNetCDF) then
                        call WriteFieldNetCDF(iP, i)                        
                    endif
                   
                    if      (Me%Field(iP)%Dim==2) then
                        deallocate(Me%Field(iP)%Value2DOut)
                    else
                        deallocate(Me%Field(iP)%Value3DOut)
                    endif
                enddo
            endif                
        enddo    
        
    end subroutine WriteRotation
    
    !------------------------------------------------------------------------    


    subroutine ComputeVectorIntensity(iPx, step, iP)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPx, step, iP
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                      :: Field_UV        
        integer                                     :: i, j, k, mask
        !Begin-----------------------------------------------------------------

        allocate(Field_UV)
        nullify(Field_UV)
        
        Field_UV => Me%Field(iPx)


        if      (Me%Field(iP)%Dim==3) then
        
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value3DOut(i,j,k) = 0.0
                endif
                
                if(Me%Mapping%Value3DOut(i,j,k) == 1)then
                    if (Step ==1) then
                        Me%Field(iP)%Value3DOut(i,j,k) = Field_UV%Value3DOut(i,j,k)**2
                    else if (Step ==2) then
                        Me%Field(iP)%Value3DOut(i,j,k) = sqrt(Me%Field(iP)%Value3DOut(i,j,k) + Field_UV%Value3DOut(i,j,k)**2)
                    endif
                endif

            enddo
            enddo
            enddo
        
        else if (Me%Field(iP)%Dim==2) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value2DOut(i,j) = 0.0
                endif
                
                if (Me%Depth%Dim3D) then
                    mask = Me%Mapping%Value3DOut(i,j,Me%WorkSize%KUB)
                else
                    mask = Me%Mapping%Value2DOut(i,j)
                endif                
                                
                if (mask == 1)then
                    if (Step ==1) then
                        Me%Field(iP)%Value2DOut(i,j) = Field_UV%Value2DOut(i,j)**2
                    else if (Step ==2) then
                        Me%Field(iP)%Value2DOut(i,j) = sqrt(Me%Field(iP)%Value2DOut(i,j) + Field_UV%Value2DOut(i,j)**2)
                    endif
                endif

            enddo
            enddo

        endif        

        nullify(Field_UV)

        

    end subroutine ComputeVectorIntensity
    
    !------------------------------------------------------------------------    

    !------------------------------------------------------------------------    


    subroutine ComputePeakPeriod(Prop2D, TpProp)
    
        !Arguments-------------------------------------------------------------
        real, dimension(:,:  ), pointer             :: Prop2D
        real, dimension(:,:  ), pointer             :: TpProp
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        !Begin-----------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if(Prop2D(i, j) > 0.) then
                    
                TpProp(i, j) =1. / Prop2D(i, j)
                    
            endif
        enddo
        enddo
            
    end subroutine ComputePeakPeriod
    
    !------------------------------------------------------------------------  
    
    !------------------------------------------------------------------------    

    !------------------------------------------------------------------------    


    subroutine ComputeAverageInDepth(Prop3D, AvProp)
    
        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer             :: Prop3D
        real, dimension(:,:  ), pointer             :: AvProp
        !Local-----------------------------------------------------------------
        real (8)                                    :: LayerTick, WaterColumn 
        integer                                     :: i, j, k
        !Begin-----------------------------------------------------------------
        
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            WaterColumn = 0.
            AvProp(i, j)= 0.
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB

                if(Me%Mapping%Value3DOut(i,j,k) == 1)then
                    
                    LayerTick    = Me%Depth%Value3DOut(i,j,k-1) - Me%Depth%Value3DOut(i,j,k)
                    WaterColumn  = WaterColumn  + LayerTick
                    AvProp(i, j) = AvProp(i, j) + LayerTick * Prop3D(i, j, k)
                    
                endif

            enddo
            if (WaterColumn > 0.) then
                AvProp(i, j)    =  AvProp(i, j) / WaterColumn
            endif                    
        enddo
        enddo
            
    end subroutine ComputeAverageInDepth
    
    !------------------------------------------------------------------------    

    subroutine ComputePowerFromEnergy(iPt, iP, inst, Option) 
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPt, iP, inst
        character(len=*)                            :: Option
        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: Power, Energy
        real                                        :: dt
        integer                                     :: iStart, iEnd, i, j
        logical                                     :: NullEnergy
        !Begin-----------------------------------------------------------------
        
    !------------------------------------------------------------------------    
        Power => Me%Field(iP)%Value2DOut    
                    
        Energy => Me%Field(iPt)%Value2DOut    
                    
                   
        Power(:,:) = 0.
        
        if      (Option == "Start")  then
            iEnd    = inst+1 
            iStart  = inst
        elseif  (Option == "Normal") then
            iEnd    = inst+1 
            iStart  = inst-1               
        elseif  (Option == "End"  .or. Option == "ReStart" ) then                
            iEnd    = inst
            iStart  = inst-1               
        endif                
            
        dt = HDF5TimeInstant(iEnd) - HDF5TimeInstant(iStart)
            
        call ReadFieldHDF5(iPt, iEnd)
                    
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%Mapping%Value2DOut(i,j) == 1) then
                Power(i, j) =  Energy(i, j) /dt
            endif
        enddo
        enddo  
        
        NullEnergy = ModelStartHour(Me%Field(ip)%EnergyStartingHours, Me%Field(ip)%N_EnergyStartH, iStart)
        
        if (.not. NullEnergy) then
            
            call ReadFieldHDF5(iPt, iStart)            
            
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if(Me%Mapping%Value2DOut(i,j) == 1) then
                    Power(i, j) =  Power(i, j) - Energy(i, j) /dt
                endif
            enddo
            enddo
            
        endif                        
            
        if  (Option == "ReStart") then

            iEnd    = inst+1 
            iStart  = inst
            
            dt = HDF5TimeInstant(iEnd) - HDF5TimeInstant(iStart)            
            
            call ReadFieldHDF5(iPt, iEnd)            
            
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if(Me%Mapping%Value2DOut(i,j) == 1) then
                   Power(i, j) =  (Power(i, j) + Energy(i, j) /dt) * 0.5
                endif
            enddo
            enddo         
            
        endif
            
    end subroutine ComputePowerFromEnergy
    
    !------------------------------------------------------------------------    
    
    !------------------------------------------------------------------------    

    subroutine ComputeInstFromAvgSinceModelStart(iPt, iP, inst, Option) 
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPt, iP, inst
        character(len=*)                            :: Option
        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: ValueInst, Average
        real                                        :: EndPeriod, StartPeriod, dt
        integer                                     :: iStart, iEnd, i, j
!        logical                                     :: NullAverage
        !Begin-----------------------------------------------------------------
        
    !------------------------------------------------------------------------    
        ValueInst => Me%Field(iP)%Value2DOut    
                    
        Average => Me%Field(iPt)%Value2DOut    
                    
                   
        ValueInst(:,:) = 0.
        
        if      (Option == "Start")  then
            iEnd    = inst+1
            iStart  = inst
        elseif  (Option == "Normal") then
            iEnd    = inst+1 
            iStart  = inst-1               
        elseif  (Option == "End"  .or. Option == "ReStart" ) then                
            iEnd    = inst
            iStart  = inst-1               
        endif                
            
        EndPeriod   = HDF5TimeInstant(iEnd  ) - HDF5TimeInstant(1)
        StartPeriod = HDF5TimeInstant(iStart) - HDF5TimeInstant(1)
        dt          = EndPeriod - StartPeriod
            
        call ReadFieldHDF5(iPt, iEnd)
                    
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if(Me%Mapping%Value2DOut(i,j) == 1) then
                ValueInst(i, j) =  Average(i, j) * EndPeriod / dt
            endif
        enddo
        enddo  
        
        
        !NullAverage = ModelStartHour(Me%Field(ip)%AvModelStartingHours, Me%Field(ip)%N_AvModelStartH, iStart)
        
        !if (.not. NullAverage) then
            
            call ReadFieldHDF5(iPt, iStart)            
            
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if(Me%Mapping%Value2DOut(i,j) == 1) then
                    ValueInst(i, j) =  ValueInst(i, j) - Average(i, j) * StartPeriod / dt
                endif
            enddo
            enddo
            
        !endif                        
            
        if  (Option == "ReStart") then

            iEnd    = inst+1 
            iStart  = inst
            
            dt = HDF5TimeInstant(iEnd) - HDF5TimeInstant(iStart)            
            
            call ReadFieldHDF5(iPt, iEnd)            
            
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB
                if(Me%Mapping%Value2DOut(i,j) == 1) then
                   ValueInst(i, j) =  (ValueInst(i, j) + Average(i, j)) * 0.5
                endif
            enddo
            enddo         
            
        endif
            
    end subroutine ComputeInstFromAvgSinceModelStart
    
    !------------------------------------------------------------------------     

    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant

        !Local-----------------------------------------------------------------
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5,                               &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
            
        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))
        deallocate(TimeVector)
            
    end function HDF5TimeInstant

    !------------------------------------------------------------------------    

    subroutine ComputePrecipitation(Reflectivity, Precipitation, it, Step1) 
    
        !Arguments-------------------------------------------------------------
        real, dimension(:,:),   pointer             :: Reflectivity, Precipitation
        integer                                     :: it
        logical                                     :: Step1
        !Local-----------------------------------------------------------------
        real (8)                                    :: aux
        real (8), save                              :: DT
        integer                                     :: i, j
        !Begin-----------------------------------------------------------------

        if (it > 0) then
            if (Step1) then
                ![s]
                DT = Me%Date%ValueInTotal(it+1)-Me%Date%ValueInTotal(it)
                ![h]
                DT = DT / 3600.
            endif
        endif            


        !it is assumed that the units of reflectivity is dBZ        
        !It assumed the units of precipitation is mm or in other words mm porecipitation between two consecutive instants
        !For the first instant precipitation is always ZERO
        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if(Me%Mapping%Value2DOut(i,j) == 1)then
                if (it > 0) then
                    ![mm/h]
                    if (Reflectivity(i, j) > 5) then
                        aux  = (10.**(Reflectivity(i, j)/10.)/200.)**(5./8.)
                    else
                        aux = 0
                    endif                        
                    !integral between two consecutive instants
                    ![mm]               = [mm/h] * [h]
                    if (Step1) then
                        Precipitation(i, j) = aux * DT * 0.5
                    else
                        Precipitation(i, j) = Precipitation(i, j) + aux * DT * 0.5
                    endif   
                else
                    Precipitation(i, j) = 0.
                endif
            endif
                
        enddo
        enddo

            
    end subroutine ComputePrecipitation
    
    !------------------------------------------------------------------------    


    subroutine ComputeVectorDirection(iPx, step, iP)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPx, step, iP
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                      :: Field_UV        
        integer                                     :: i, j, k
        real                                        :: DperR
        !Begin-----------------------------------------------------------------

        allocate(Field_UV)
        nullify(Field_UV)
        
        Field_UV => Me%Field(iPx)


        if      (Me%Field(iP)%Dim==3) then
        
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value3DOut(i,j,k) = 0.0
                endif
                
                if(Me%Mapping%Value3DOut(i,j,k) == 1)then
                    if (Step ==1) then
                        Me%Field(iP)%Value3DOut(i,j,k) = Field_UV%Value3DOut(i,j,k)
                    else if (Step ==2) then
                    
                        DperR=180/pi
                        
                        !CartesianDir_
                        !Dir = atan2(y,x)*DperR
                        Me%Field(iP)%Value3DOut(i,j,k) = atan2(Field_UV%Value3DOut(i,j,k), Me%Field(iP)%Value3DOut(i,j,k))*DperR
                        
                        if          (Me%Field(iP)%DirectionReferential == NauticalWind_   ) then
                            Me%Field(iP)%Value3DOut(i,j,k) = 270. - Me%Field(iP)%Value3DOut(i,j,k)
                        else if     (Me%Field(iP)%DirectionReferential == NauticalCurrent_) then
                            Me%Field(iP)%Value3DOut(i,j,k) = 90. -  Me%Field(iP)%Value3DOut(i,j,k)
                        endif

                        if      (Me%Field(iP)%Value3DOut(i,j,k)>360.)then
                             Me%Field(iP)%Value3DOut(i,j,k) = Me%Field(iP)%Value3DOut(i,j,k) - 360.
                        else if (Me%Field(iP)%Value3DOut(i,j,k)<0.  )then
                             Me%Field(iP)%Value3DOut(i,j,k) = Me%Field(iP)%Value3DOut(i,j,k) + 360.
                        end if
                        
                    endif
                endif

            enddo
            enddo
            enddo
        
        else if (Me%Field(iP)%Dim==2) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value2DOut(i,j) = 0.0
                endif
                                
                if(Me%Mapping%Value2DOut(i,j) == 1)then
                    if (Step ==1) then
                        Me%Field(iP)%Value2DOut(i,j) = Field_UV%Value2DOut(i,j)
                    else if (Step ==2) then
                    
                        DperR=180/pi
                        !CartesianDir_
                        !Dir = atan2(y,x)*DperR                        
                        Me%Field(iP)%Value2DOut(i,j) = atan2(Field_UV%Value2DOut(i,j),Me%Field(iP)%Value2DOut(i,j))*DperR
                        
                        if          (Me%Field(iP)%DirectionReferential == NauticalWind_   ) then
                            Me%Field(iP)%Value2DOut(i,j) = 270. - Me%Field(iP)%Value2DOut(i,j)
                        else if     (Me%Field(iP)%DirectionReferential == NauticalCurrent_) then
                            Me%Field(iP)%Value2DOut(i,j) = 90.  - Me%Field(iP)%Value2DOut(i,j)
                        endif
                        
                        
                        if      (Me%Field(iP)%Value2DOut(i,j)>360.) then
                             Me%Field(iP)%Value2DOut(i,j)=Me%Field(iP)%Value2DOut(i,j) - 360.
                        else if (Me%Field(iP)%Value2DOut(i,j)<0.  ) then
                             Me%Field(iP)%Value2DOut(i,j)=Me%Field(iP)%Value2DOut(i,j) + 360.
                        end if

                    endif
                endif

            enddo
            enddo

        endif        

        nullify(Field_UV)

        

    end subroutine ComputeVectorDirection
    !------------------------------------------------------------------------


    subroutine ModifyRotatedVector(iPx, step, iP)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPx, step, iP
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                      :: Field_UV        
        integer                                     :: i, j, k
        real                                        :: AngleX, AngleY, Xgrid, Ygrid
        
        !Begin-----------------------------------------------------------------

        allocate(Field_UV)
        nullify(Field_UV)
        
        Field_UV => Me%Field(iPx)


        if      (Me%Field(iP)%Dim==3) then
        
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value3DOut(i,j,k) = 0.0
                endif
                
                if(Me%Mapping%Value3DOut(i,j,k) == 1)then
                    if (Step ==1) then
                        Me%Field(iP)%Value3DOut(i,j,k) = Field_UV%Value3DOut(i,j,k)
                    else if (Step ==2) then
                        
                        AngleX = Me%Field(iP)%GridRotation * Pi / 180
                        AngleY = AngleX + Pi / 2.
                        
                        call FromCartesianToGrid (Xcart  = Me%Field(iP)%Value3DOut(i,j,k),&
                                                  Ycart  = Field_UV%Value3DOut(i,j,k),  &
                                                  Tetha1 = AngleX,                      &
                                                  Tetha2 = AngleY,                      &
                                                  Xgrid  = Xgrid,                       &
                                                  Ygrid  = Ygrid)
                        
                        if     (Me%Field(iP)%VectorComponent == Zonal) then

                            Me%Field(iP)%Value3DOut(i,j,k) = Xgrid
                        
                        elseif (Me%Field(iP)%VectorComponent == Meridional) then
                            
                            Me%Field(iP)%Value3DOut(i,j,k) = Ygrid
                            
                        endif
                        
                    endif
                endif

            enddo
            enddo
            enddo
        
        else if (Me%Field(iP)%Dim==2) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value2DOut(i,j) = 0.0
                endif
                                
                if(Me%Mapping%Value2DOut(i,j) == 1)then
                    if (Step ==1) then
                        Me%Field(iP)%Value2DOut(i,j) = Field_UV%Value2DOut(i,j)
                    else if (Step ==2) then
                    
                        AngleX = Me%Field(iP)%GridRotation * Pi / 180
                        AngleY = AngleX + Pi / 2.
                        
                        call FromCartesianToGrid (Xcart  = Me%Field(iP)%Value2DOut(i,j),&
                                                  Ycart  = Field_UV%Value2DOut(i,j),    & 
                                                  Tetha1 = AngleX,                      &
                                                  Tetha2 = AngleY,                      &
                                                  Xgrid  = Xgrid,&
                                                  Ygrid  = Ygrid)
                        
                        if     (Me%Field(iP)%VectorComponent == Zonal) then

                            Me%Field(iP)%Value2DOut(i,j) = Xgrid
                        
                        elseif (Me%Field(iP)%VectorComponent == Meridional) then
                            
                            Me%Field(iP)%Value2DOut(i,j) = Ygrid
                            
                        endif

                    endif
                endif

            enddo
            enddo

        endif        

        nullify(Field_UV)

        

    end subroutine ModifyRotatedVector
    !------------------------------------------------------------------------

    
    subroutine ComputeBeaufort(iPx, iP)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPx, iP
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                      :: Field_x        
        integer                                     :: i, j, k
        !Begin-----------------------------------------------------------------

        allocate(Field_x)
        nullify(Field_x)
        
        Field_x => Me%Field(iPx)


        if      (Me%Field(iP)%Dim==3) then
        
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if(Me%Mapping%Value3DOut(i,j,k) == 1)then
                    if     (Me%Field(iP)%Beaufort    ) then
                        Me%Field(iP)%Value3DOut(i,j,k) = WindBeaufortScale(Field_x%Value3DOut(i,j,k))
                    elseif (Me%Field(iP)%WaveBeaufort) then
                        Me%Field(iP)%Value3DOut(i,j,k) = WaveBeaufortScale(Field_x%Value3DOut(i,j,k))
                    endif                        
                endif

            enddo
            enddo
            enddo
        
        else if (Me%Field(iP)%Dim==2) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if(Me%Mapping%Value2DOut(i,j) == 1)then
                    if     (Me%Field(iP)%Beaufort) then
                        Me%Field(iP)%Value2DOut(i,j) = WindBeaufortScale(Field_x%Value2DOut(i,j))
                    elseif (Me%Field(iP)%WaveBeaufort) then
                        Me%Field(iP)%Value2DOut(i,j) = WaveBeaufortScale(Field_x%Value2DOut(i,j))
                    endif                                                
                endif

            enddo
            enddo

        endif        

        nullify(Field_x)

        

    end subroutine ComputeBeaufort

    !------------------------------------------------------------------------    
    
    subroutine ComputeRelativeHumidity(T, P, SH, RH, step)
    
        !Argument-----------------------------------------------------
        real, dimension(:,:), pointer               :: SH, T, P, RH
        integer                                     :: step

        !Local--------------------------------------------------------
        ! mixture(%), dewpoint mixture (%), dewpoint water pressure (Pa), specific humidity (%)
        real                                        :: w
        integer                                     :: i, j
        
        !Begin--------------------------------------------------------

        if (.not.associated(RH)) then
            stop 'ComputeRelativeHumidity - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif

        if (step == 1) then
            if (.not.associated(T)) then
                stop 'ComputeRelativeHumidity - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            endif
        else if (Step == 2) then
            if (.not.associated(P)) then
                stop 'ComputeRelativeHumidity - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
            endif
        else if (step == 3) then
            if (.not.associated(SH)) then
                stop 'ComputeRelativeHumidity - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            endif
        endif


        do j = Me%WorkSize%JLB, Me%WorkSize%JUB
        do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        
            if (step == 1) then
            
                ! from mm5tograds
                if (T(i,j) <= 0.) then
                    RH(i,j) = 6.11 * EXP (22.514 - 6150./(T(i,j) + AbsoluteZero))
                else
                    RH(i,j) = 6.112* EXP (17.67*((T(i,j) - 273.15 + AbsoluteZero)    &
                              /(T(i,j) - 29.65 + AbsoluteZero)))               
                endif
                
            else if (Step == 2) then

                RH(i,j) = 0.622 * RH(i,j) / ((P(i,j)/100.) - RH(i,j))
                
            else if (step == 3) then
                 
                w = SH(i,j) / (1 - SH(i,j))
                ! 5% < Rel. Hum. < 100%
                RH(i,j) =  min(max(100. * w /RH(i,j),5.), 100.) * 0.01
                !Me%RelativeHumidity(i,j) =  w / ws
                
            endif
            
        enddo
        enddo

    end subroutine ComputeRelativeHumidity
    !------------------------------------------------------------------------

    subroutine ComputeVectorRotation(iPx, step, iP, Component)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: iPx, step, iP, Component
        
        !Local-----------------------------------------------------------------
        type(T_Field), pointer                      :: Field_UV        
        real                                        :: AngleX, AngleY
        integer                                     :: i, j, k
        !Begin-----------------------------------------------------------------

        allocate(Field_UV)
        nullify(Field_UV)
        
        Field_UV => Me%Field(iPx)
    
        if      (Me%Field(iP)%Dim==3) then
        
            do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value3DOut(i,j,k) = 0.0
                endif
                
                if (Me%Mapping%Value3DOut(i,j,k) == 1)then
                !if (Field_UV%Value3DOut  (i,j,k) > FillValueReal/2.) then
                
                    AngleX = Me%LongLat%RotationX(i, j) 
                    AngleY = Me%LongLat%RotationY(i, j) 
                                    
                    if (Step ==1) then
                        if (Component == Zonal) then
                            Me%Field(iP)%Value3DOut(i,j,k) = Me%Field(iP)%Value3DOut(i,j,k) + &
                                                             Field_UV%Value3DOut(i,j,k)* cos(AngleX)
                        else if (Component == Meridional) then
                            Me%Field(iP)%Value3DOut(i,j,k) = Me%Field(iP)%Value3DOut(i,j,k) + &
                                                             Field_UV%Value3DOut(i,j,k)* sin(AngleX)
                        endif
                    else if (Step ==2) then
                        if (Component == Zonal) then
                            Me%Field(iP)%Value3DOut(i,j,k) = Me%Field(iP)%Value3DOut(i,j,k) + &
                                                             Field_UV%Value3DOut(i,j,k)* cos(AngleY)
                        else if (Component == Meridional) then
                            Me%Field(iP)%Value3DOut(i,j,k) = Me%Field(iP)%Value3DOut(i,j,k) + &
                                                             Field_UV%Value3DOut(i,j,k)* sin(AngleY)
                        endif
                    endif
                endif

            enddo
            enddo
            enddo
        
        else if (Me%Field(iP)%Dim==2) then

            do j = Me%WorkSize%JLB, Me%WorkSize%JUB
            do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Step ==1) then
                    Me%Field(iP)%Value2DOut(i,j) = 0.0
                endif
                                
                if(Me%Mapping%Value2DOut(i,j) == 1)then

                    AngleX = Me%LongLat%RotationX(i, j) 
                    AngleY = Me%LongLat%RotationY(i, j) 
                                    
                    if (Step ==1) then
                        if (Component == Zonal) then
                            Me%Field(iP)%Value2DOut(i,j) = Me%Field(iP)%Value2DOut(i,j) + Field_UV%Value2DOut(i,j)* cos(AngleX)
                        else if (Component == Meridional) then
                            Me%Field(iP)%Value2DOut(i,j) = Me%Field(iP)%Value2DOut(i,j) + Field_UV%Value2DOut(i,j)* sin(AngleX)
                        endif
                    else if (Step ==2) then
                        if (Component == Zonal) then
                            Me%Field(iP)%Value2DOut(i,j) = Me%Field(iP)%Value2DOut(i,j) + Field_UV%Value2DOut(i,j)* cos(AngleY)
                        else if (Component == Meridional) then
                            Me%Field(iP)%Value2DOut(i,j) = Me%Field(iP)%Value2DOut(i,j) + Field_UV%Value2DOut(i,j)* sin(AngleY)
                        endif
                    endif
                endif
            enddo
            enddo

        endif        

        nullify(Field_UV)

    end subroutine ComputeVectorRotation
    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer, dimension(4)                       :: Aux4
        logical                                     :: Exists

        !Begin-----------------------------------------------------------------

        call GetData(Me%FileName,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILENAME',                                         &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR00'
        
        if(iflag > 0)then
            !Just read one file, otherwise read a block with a list of files
            Me%OnlyOneFile = .true.
            
            inquire(FILE = Me%FileName, EXIST = Exists)
            if (.not. Exists) then
                write(*,*)'Netcdf file does not exist:'//trim(Me%FileName)
                stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR01'
            endif
        endif
        
        
        

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',                             &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


        call GetData(Me%GeometryFilename,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUT_GEOMETRY_FILENAME',                         &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
        

        call GetData(Me%PropNumber,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'PROPERTIES_NUMBER',                                &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
        if (iflag == 0)            stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
        
        call GetData(Me%OutHDF5,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'HDF5_OUT',                                         &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
        
        call GetData(Me%OutNetCDF,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NETCDF_OUT',                                       &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
        

        call GetData(Me%ReadInvertXY,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'READ_INVERT_XY',                                   &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR120'
        
        call GetData(Me%ReadInvertLat,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'READ_INVERT_LAT',                                  &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
        

        call GetData(Aux4,                                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'WINDOW_OUT',                                       &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (iflag==4) then
            Me%WindowOut%ON  = .true.
            Me%WindowOut%ILB = Aux4(1)
            Me%WindowOut%IUB = Aux4(2)
            Me%WindowOut%JLB = Aux4(3)
            Me%WindowOut%JUB = Aux4(4)
        else            
            Me%WindowOut%ON  = .false.
        endif
        
        call GetData(Me%NaN_2_Null,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NAN_2_NULL',                                       &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR140'
                

        if (Me%OutNetCDF) then
        
            call ReadOutNetCDFOptions
        
        endif

        allocate(Me%Field(1:Me%PropNumber))

        call ReadTimeOptions
        
        call ReadGridOptions
        
        call ReadFieldOptions


    end subroutine ReadOptions

    !------------------------------------------------------------------------

    subroutine ReadOutNetCDFOptions
    
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

    
            call GetData(Me%NetCDF_Out%Name,                                            &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'OUTPUT_NETCDF_FILE',                           &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

            call GetData(Me%NetCDF_Out%Title,                                           &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_TITLE',                                 &
                         Default      = 'EASYCO MODEL SOLUTIONS',                       &                         
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            call GetData(Me%NetCDF_Out%Convention,                                      &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_CONVENTION',                            &
                         Default      = 'CF-1.0',                                       &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            call GetData(Me%NetCDF_Out%Version,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_VERSION',                               &
                         Default      = '3.6.1',                                        &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

            call GetData(Me%NetCDF_Out%History,                                         &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_HISTORY',                               &
                         Default      = '-',                                            &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
            call GetData(Me%NetCDF_Out%Source,                                          &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_SOURCE',                                &
                         default      = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'

            
            call GetData(Me%NetCDF_Out%Institution,                                     &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_INSTITUTION',                           &
                         Default      = 'Instituto Superior Tecnico',                   &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            call GetData(Me%NetCDF_Out%References,                                      &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_REFERENCES',                            &
                         Default      = 'http://www.mohid.com',                         &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'

            call GetData(Me%NetCDF_Out%iDate,                                           &
                         Me%ObjEnterData,iflag,                                         &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'NETCDF_DATE',                                  &
                         ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',                   &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadOutNetCDFOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'    
            
    end subroutine ReadOutNetCDFOptions            
    
    !------------------------------------------------------------------------    
    
    subroutine ReadTimeOptions

        !Local-----------------------------------------------------------------
        real, dimension(6)                          :: AuxTime, DefaultTime
        logical                                     :: BlockFound
        integer                                     :: iflag, STAT_CALL
        !begin-----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   begin_time, end_time,                                &
                                   BlockInBlockFound = BlockFound,  STAT = STAT_CALL)


IS:     if(STAT_CALL .EQ. SUCCESS_) then
            
BF:         if (BlockFound) then


                call GetData(Me%Date%NetCDFName,                                        &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME',                              &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "time",                                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

               !The time in netcdf is assumed in days from 1950/1/1 : 0h:0m:0s (Mercator Ocean convention)
                DefaultTime(:) = 0.
                DefaultTime(1) = 1950
                DefaultTime(2) = 1
                DefaultTime(3) = 1

                call GetData(AuxTime,                                                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'REFERENCE_DATE_IN',                        &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'    
                
                if (iflag /= 6) AuxTime = DefaultTime
    

                call SetDate (Me%Date%RefDateTimeIn, Year    = AuxTime(1),              &
                                                     Month   = AuxTime(2),              &
                                                     Day     = AuxTime(3),              &
                                                     Hour    = AuxTime(4),              &
                                                     Minute  = AuxTime(5),              &
                                                     Second  = AuxTime(6))

                call StartComputeTime(Me%ObjTime, Me%Date%RefDateTimeIn, Me%Date%RefDateTimeIn, &
                                      Me%Date%RefDateTimeIn, DT = 0.0,                        &      
                                      VariableDT = .false., STAT = STAT_CALL)   
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
                                     
                call GetData(AuxTime,                                                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'REFERENCE_DATE_OUT',                       &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'  
                
                if (iflag==0) then
                   !The time in netcdf is assumed in days from 1950/1/1 : 0h:0m:0s (Mercator Ocean convention)
                    AuxTime(:) = 0.
                    AuxTime(1) = 1950
                    AuxTime(2) = 1
                    AuxTime(3) = 1                
                endif  

                call SetDate (Me%Date%RefDateTimeOut, Year    = AuxTime(1),             &
                                                      Month   = AuxTime(2),             &
                                                      Day     = AuxTime(3),             &
                                                      Hour    = AuxTime(4),             &
                                                      Minute  = AuxTime(5),             &
                                                      Second  = AuxTime(6))                                      

                call GetData(Me%Date%ValueIn%DataType,                                  &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'DATE_TYPE_IN',                             &
                             default      = Real8_,                                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'    
                
                Me%Date%ValueIn%Dim      = 1
          

                call GetData(Me%Date%UnitsFactor,                                       &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'UNITS_FACTOR',                             &
                             default      = 3600.,                                      &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'    

                call GetData(Me%Date%RefAttribute,                                      &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'REF_DATE_ATTRIBUTE',                       &
                             default      = .true.,                                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'    
                
                if (Me%Date%RefAttribute) then
        
                    call GetData(Me%Date%RefAttributeName,                              &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'REF_DATE_ATTRIBUTE_NAME',              &
                                 !CF convention
                                 default      = "units",                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'    

                    call GetData(Me%Date%RefDateName,                                   &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'REF_DATE_NAME',                        &
                                 default      = trim(null_str),                         &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR82'    
                               
                    !off-set in seconds
                    call GetData(Me%Date%RefDateOffSet,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'REF_DATE_OFF_SET',                     &
                                 !CF convention
                                 default      = 0.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR85'                        

                    !off-set from a specific property attribute 
                    call GetData(Me%Date%RefDateOffSetFromAtt,                          &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlock,                              &
                                 keyword      = 'REF_DATE_OFF_SET_FROM_ATT',            &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR86'
                    
                    if (Me%Date%RefDateOffSetFromAtt) then

                        !off-set from a specific property name 
                        call GetData(Me%Date%RefDateOffSetProp,                         &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlock,                          &
                                     keyword      = 'REF_DATE_OFF_SET_PROP',            &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR87'

                        if (iflag == 0) then
                            stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR87a'
                        endif

                        !off-set from a specific attribute name 
                        call GetData(Me%Date%RefDateOffSetAtt,                          &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlock,                          &
                                     keyword      = 'REF_DATE_OFF_SET_ATT',             &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR88'
                        
                        if (iflag == 0) then
                            stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR88a'
                        endif
                        

                        !off-set from a specific property attribute name - multiplying factor
                        call GetData(Me%Date%RefDateOffSetAttFactor,                    &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlock,                          &
                                     keyword      = 'REF_DATE_OFF_SET_ATT_FACTOR',      &
                                     default      = 3600.,                              &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR89'
                        
                    endif                        
                    
                endif
                
                call GetData(Me%Date%NetCDFDimName,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_DIM',                               &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = Me%Date%NetCDFName,                         &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'

                
                
 
        
            else BF
            
                stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'    
                
            endif BF
        else IS
        
            stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR200'    
        
        endif IS                             
                              
                              
    end subroutine ReadTimeOptions                              


    !------------------------------------------------------------------------
    
    subroutine ReadGridOptions

        !Local-----------------------------------------------------------------
        real, allocatable, dimension(:)             :: Aux1D
        logical                                     :: BlockFound, geo_sigma
        integer                                     :: iflag, STAT_CALL
        !begin-----------------------------------------------------------------

        call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                   begin_grid , end_grid ,                              &
                                   BlockInBlockFound = BlockFound,  STAT = STAT_CALL)


IS:     if(STAT_CALL .EQ. SUCCESS_) then
            
BF:         if (BlockFound) then


                call GetData(Me%LongLat%NetCDFNameLat,                                  &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_LAT',                          &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "latitude",                                 &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

                call GetData(Me%LongLat%NetCDFNameLong,                                 &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_LONG',                         &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "longitude",                                &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                
                call GetData(Me%LongLat%Starts180W,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'STARTS_180W',                              &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = .true.,                                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
                Me%LongLat%LatIn%DataType  = Real8_
                Me%LongLat%LongIn%DataType = Real8_
 
                Me%LongLat%LatIn%Dim       = 2
                Me%LongLat%LongIn%Dim      = Me%LongLat%LatIn%Dim

                !cut necessary to compute the location of cell corners by defaut is 1 cell 
                !but can be more. Needs to be more in the netcdf Delft3D files because
                !in this case the cell corners do not have position.
                call GetData(Me%LongLat%dij,                                            &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'DIJ',                                      &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = 1,                                          &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR35'                

                call GetData(Me%Depth%NetCDFName,                                       &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_DEPTH',                        &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

                if  (iflag == 0) then 
                    Me%Depth%Dim3D = .false.
                    Me%Depth%kmax = -99
                else 
                    Me%Depth%Dim3D = .true.
                endif
                
                Me%Depth%ValueIn%DataType = Real8_

                call GetData(Me%Depth%NetCDFNameWL,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_DEPTH_WL',                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR43'

                
                call GetData(Me%Bathym%NetCDFName,                                      &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_BATHYM',                       &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                
                if  (iflag == 0) then 
                    Me%Bathym%ON = .false.
                else 
                    Me%Bathym%ON = .true.
                endif  
                
                call GetData(Me%Bathym%FileName,                                        &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'BATHYM_FILENAME',                          &
                             default      = 'Batim.dat',                                &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'                

                call GetData(Me%Bathym%FromMapping,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'BATHYM_FROM_MAP',                          &
                             default      = .false.,                                    &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'     
                
                call GetData(Me%Bathym%ValueIn%DataType,                                &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'BATHYM_TYPE_IN',                           &
                             default      = Real8_,                                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
                

                call GetData(Me%Bathym%InvertReferential,                               &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'BATHYM_INVERT_REFERENTIAL',                &
                             default      = .false.,                                    &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR85'                
                
                Me%Bathym%ValueIn%Dim      = 2
                
                call GetData(Me%FieldInType,                                            &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlock,                                  &
                             keyword      = 'FIELD_TYPE_IN',                            &
                             default      = Real8_,                                     &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadTimeOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'                   
                
                
                if (.not. Me%Bathym%ON) then
                
                    call GetData(Me%Bathym%Default,                                         &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlockInBlock,                           &
                                 keyword      = 'BATHYM_DEFAULT',                           &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                                 default      = 0.,                                         &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
                
                endif
                
                call GetData(Me%Mapping%NetCDFName,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_MAPPING',                      &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'   
                
                if (iflag == 0) then
                    Me%Mapping%ON = .false.
                else
                    Me%Mapping%ON = .true.
                endif
                
                call GetData(Me%Mapping%Dim3D,                                          &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'NETCDF_NAME_MAPPING_3D',                   &
                             default      = Me%Depth%Dim3D,                             & 
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR35'   
                
                call GetData(Me%Mapping%From2D_To_3D,                                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'FROM_2D_TO_3D',                            &
                             default      = .false.,                                    & 
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR37'
                
                call GetData(Me%Mapping%Limit,                                          &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'MAPPING_LIMIT',                            &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = 0.5,                                        &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
                
                call GetData(Me%Mapping%Instant,                                        &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'MAPPING_INSTANT',                          &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = 1,                                          &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'

                Me%Mapping%ValueIn%DataType = Real4_
                
                call GetData(Me%Depth%Offset,                                           &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'DEPTH_OFFSET',                             &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = 0.,                                         &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
 
                !1 - sigma_ 2 - z-level 3 - hybrid
                call GetData(Me%Depth%GeoVert,                                          &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'GEO_VERT',                                 &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = z_level,                                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'  
                
                if (iflag == 0) then
                
                    !.true. - sigma_ false - z-level
                    call GetData(geo_sigma,                                                 &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlockInBlock,                           &
                                 keyword      = 'GEO_SIGMA',                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                                 default      = .false.,                                    &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'  
                    
                    if (geo_sigma) then
                        Me%Depth%GeoVert = sigma_
                    else
                        Me%Depth%GeoVert = z_level
                    endif

                endif     

                if (Me%Depth%GeoVert == Hybrid) then
                    Me%Depth%ValueIn%Dim      = 4
                else                
                    Me%Depth%ValueIn%Dim      = 1
                endif                    
                
                if (Me%Depth%GeoVert == sigma_ .and. .not. Me%Bathym%ON) then
                    stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR95'  
                endif             

                if (Me%Depth%GeoVert == sigma_) then
                
                    call GetData(Me%Depth%theta_s,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'THETA_S',                              &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 default      = 0.,                                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
                    
                    call GetData(Me%Depth%theta_b,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'THETA_B',                              &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 default      = 0.,                                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR110'

                    call GetData(Me%Depth%hc,                                           &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'HC',                                   &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 default      = 0.,                                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR120'
                    
                    call GetData(Me%Depth%NetCDFNameFace,                                   &
                                 Me%ObjEnterData, iflag,                                    &
                                 SearchType   = FromBlockInBlock,                           &
                                 keyword      = 'NETCDF_NAME_DEPTH_FACE',                   &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR45'
                    
                    if (iflag == 0) then
                        write(*,*) 'keyword not defined NETCDF_NAME_DEPTH_FACE'
                        Me%Depth%NetCDFNameFaceOff = .true.
                    else
                        Me%Depth%NetCDFNameFaceOff = .false.                                            
                    endif
                    
                
                endif
                
                call GetData(Me%Depth%Interpolate,                                      &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'GEO_INTERPOLATE',                          &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = .false.,                                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
                
                if (Me%Depth%Interpolate) then
                
                    allocate(Aux1D(200))

                    call GetData(Aux1D,                                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'GEO_ZLEVEL_INTERPOL',                  &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)  
                                       
                    if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_)             &
                        stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR140'   
                    if (iflag     == 0       ) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR150'   
                    
                    Me%Depth%N_ZLevels = iflag
                    
                    allocate(Me%Depth%Zlevels(Me%Depth%N_ZLevels))
                    
                    Me%Depth%Zlevels(1:Me%Depth%N_ZLevels) = Aux1D(1:Me%Depth%N_ZLevels)
                    
                    deallocate(Aux1D)

                endif

                !down or up
                call GetData(Me%Depth%Positive,                                         &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'GEO_POSITIVE',                             &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = "down",                                     &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR160'  
                
                if (trim(Me%Depth%Positive) /= "down" .and.                             &
                    trim(Me%Depth%Positive) /= "up"   .and.                             &
                    trim(Me%Depth%Positive) /= "inverse") then
                    stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR162'  
                endif                    

                !invert layer order
                call GetData(Me%Depth%InvertLayers,                                     &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'INVERT_LAYER_ORDER',                       &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = .false.,                                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR165'  



                call GetData(Me%Depth%RomsDistortion,                                   &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'ROMS_DISTORTION',                          &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = .false.,                                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR170'  


                call GetData(Me%LongLat%Imposed,                                        &
                             Me%ObjEnterData, iflag,                                    &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      = 'LONG_LAT_IMPOSED',                         &
                             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',               &
                             default      = .false.,                                    &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR180'  
                
                if (Me%LongLat%Imposed) then

                    call GetData(Me%LongLat%LongOrig,                                   &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'LONG_ORIG',                            &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR190'  
                    
                    if (iflag == 0) then
                        stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR200'  
                    endif
                
                    call GetData(Me%LongLat%dLong,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'LONG_DX',                              &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR210'  
                    
                    if (iflag == 0) then
                        stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR220'  
                    endif

                    call GetData(Me%LongLat%LatOrig,                                    &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'LAT_ORIG',                             &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR230'  
                    
                    if (iflag == 0) then
                        stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR240'  
                    endif
                
                    call GetData(Me%LongLat%dLat,                                       &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'LAT_DY',                               &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR250'  
                    
                    if (iflag == 0) then
                        stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR260'  
                    endif
                    
                
                    call GetData(Me%Depth%RemoveNsurfLayers,                            &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'REMOVE_SURFACE_N_LAYERS',              &
                                 default      = 0,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR270'  
  
                endif
  
            else BF
            
                stop 'ReadGridOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR400'    
                
            endif BF
        else IS
        
            stop 'ReadGrid2DOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR500'    
        
        endif IS                             
                              
                              
    end subroutine ReadGridOptions                              

    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------
    
    subroutine ReadFieldOptions

        !Local-----------------------------------------------------------------
        
        logical                                     :: BlockFound
        integer                                     :: iflag, STAT_CALL, ip
        logical                                     :: check
        integer, dimension(12)                      :: Aux1DInt
        real,    dimension(2)                       :: Array2
        !begin-----------------------------------------------------------------

        do ip = 1, Me%PropNumber 
                        
            call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                &
                                       begin_field, end_field,                          &
                                       BlockInBlockFound = BlockFound,  STAT = STAT_CALL)


    IS:     if(STAT_CALL .EQ. SUCCESS_) then
                    
    BF:         if (BlockFound) then

                    call GetData(Me%Field(ip)%NetCDFName,                               &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'NETCDF_NAME',                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
                    if (iflag==0)              stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

                    call GetData(Me%Field(ip)%Add,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ADD_FACTOR',                           &
                                 default      = 0.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR30'    
                    
                    call GetData(Me%Field(ip)%Multiply,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'MULTIPLY_FACTOR',                      &
                                 default      = 1.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR40'    

                    call GetData(Me%Field(ip)%UnitsScale,                               &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'UNITS_SCALE',                          &
                                 default      = 1.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR45'    
                    
                    call GetData(Me%Field(ip)%UnitsAdd,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'UNITS_ADD',                            &
                                 default      = 0.,                                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR45'    
                    

                    call GetData(Me%Field(ip)%MinValue,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'MINIMUM_VALUE',                        &
                                 default      = FillValueReal,                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

                    call GetData(Me%Field(ip)%OldMissingValue,                          &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'OLD_MISSING_VALUE',                    &
                                 default      = FillValueReal,                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR52'                    
                    
                    call GetData(Me%Field(ip)%NewMissingValue,                          &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'NEW_MISSING_VALUE',                    &
                                 default      = FillValueReal,                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR53'
                    
                    call GetData(Array2,                                                &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'MIN_MAX_LIMITS',                       &
                                 default      = FillValueReal,                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR55'
                                        
                    if (iflag == 2) then
                        Me%Field(ip)%CheckMinMaxLimits = .true.
                        Me%Field(ip)%MinLimit          = Array2(1)
                        Me%Field(ip)%MaxLimit          = Array2(2)
                    else
                        Me%Field(ip)%CheckMinMaxLimits = .false.
                    endif


                    call GetData(Me%Field(ip)%Dim,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DIM',                                  &
                                 default      = 2,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR60'    
                    
                    call GetData(Me%Field(ip)%From2D_To_3D,                             &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'FROM_2D_TO_3D',                        &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR65'    
                    
                    if (Me%Field(ip)%From2D_To_3D .and. Me%Field(ip)%Dim /= 2)  then
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR67'
                    endif
                    
                    if (Me%Field(ip)%Dim /= 2 .and. Me%Field(ip)%Dim /= 3)              &
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR70'    

                    if (Me%Field(ip)%Dim == 2) then
                        call GetData(Me%Field(ip)%ExtractLayer,                             &
                                     Me%ObjEnterData, iflag,                                &
                                     SearchType   = FromBlockInBlock,                       &
                                     keyword      = 'EXTRACT_LAYER',                        &
                                     default      = .true.,                                 &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
                        
                        call GetData(Me%Field(ip)%LayerNumber,                              &
                                     Me%ObjEnterData, iflag,                                &
                                     SearchType   = FromBlockInBlock,                       &
                                     keyword      = 'lAYER_number',                         &
                                     default      = 1,                                      &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
                        
                        call GetData(Me%Field(ip)%LayerDim,                                 &
                                     Me%ObjEnterData, iflag,                                &
                                     SearchType   = FromBlockInBlock,                       &
                                     keyword      = 'LAYER_DIM',                            &
                                     default      = 3,                                      &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR96'                        
                        
                    endif                        
                        
                    call GetData(check,                                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'CHECK_PROPERTY',                       &
                                 default      = .true.,                                 &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
                    
                    call ConstructPropertyID (Me%Field(ip)%ID,  Me%ObjEnterData, FromBlockInBlock, check)
                    
                    !Index 1 is time
                    Me%Field(ip)%ValueIn%DataType = Me%FieldInType
                    
                    if      (Me%Field(ip)%Dim == 2) then
                        Me%Field(ip)%ValueIn%Dim = 3
                    else if (Me%Field(ip)%Dim == 3) then
                        Me%Field(ip)%ValueIn%Dim = 4
                    endif

                    call GetData(Me%Field(ip)%FromDir2Vector,                           &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DIR_TO_VECTOR',                        &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
                    
                    if (Me%Field(ip)%FromDir2Vector) then

                        call GetData(Me%Field(ip)%DirX,                                 &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'DIR_X',                            &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR120'

                        call GetData(Me%Field(ip)%DirY,                                 &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'DIR_Y',                            &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
                                        
                                        
                    endif


                    call GetData(Me%Field(ip)%ComputeIntensity,                         &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'VECTOR_INTENSITY',                     &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR140'
                    
                    if (Me%Field(ip)%ComputeIntensity .and. .not. Me%OutHDF5) then
                        write(*,*) 'To compute vector intensity need to write hdf5 output file'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR150'
                    endif
                    
                    call GetData(Me%Field(ip)%ComputeDirection,                         &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'WIND_DIRECTION',                       &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR160'

                    if (iflag == 0) then
                        call GetData(Me%Field(ip)%ComputeDirection,                     &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'VECTOR_DIRECTION',                 &
                                     default      = .false.,                            &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR170'
                    endif
                    
                    if (Me%Field(ip)%ComputeDirection .and. .not. Me%OutHDF5) then
                        write(*,*) 'To compute wind direction need to write hdf5 output file'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR180'
                    endif

                    !Possible options:
                    !    NauticalWind_    = 1
                    !    NauticalCurrent_ = 2 
                    !    CartesianDir_    = 3
                    call GetData(Me%Field(ip)%DirectionReferential,                     &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DIRECTION_REFERENTIAL',                &
                                 default      = NauticalWind_,                          &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR185'
                    
                    if (Me%Field(ip)%DirectionReferential /= NauticalWind_    .and.     &
                        Me%Field(ip)%DirectionReferential /= NauticalCurrent_ .and.     &
                        Me%Field(ip)%DirectionReferential /= CartesianDir_) then
                        write(*,*) 'Not valid Direction Referential option'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR187'
                    endif                        

                    call GetData(Me%Field(ip)%Rotation,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'VECTOR_ROTATION',                      &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR190'
                    
                    if (Me%Field(ip)%Rotation .and. .not. Me%OutHDF5) then
                        write(*,*) 'To compute vector rotation need to write hdf5 output file'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR200'
                    endif

                    
                    call GetData(Me%Field(ip)%Beaufort,                                 &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'BEAUFORT_SCALE',                      &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR210'
                    
                    if (Me%Field(ip)%Beaufort .and. .not. Me%OutHDF5) then
                        write(*,*) 'To compute Beaufort scale need to write hdf5 output file'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR220'
                    endif


                    call GetData(Me%Field(ip)%WaveBeaufort,                             &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'WAVE_BEAUFORT_SCALE',                  &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR230'
                    
                    if (Me%Field(ip)%WaveBeaufort .and. .not. Me%OutHDF5) then
                        write(*,*) 'To compute wave Beaufort scale need to write hdf5 output file'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR240'
                    endif
                    
                    call GetData(Me%Field(ip)%ComputeRotatedVector,                     &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ROTATE_VECTOR',                        &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR245'
                    
                    if (Me%Field(ip)%ComputeRotatedVector .and. .not. Me%OutHDF5) then
                        write(*,*) 'To rotate a vector need to write hdf5 output file'
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR246'
                    endif                    
                    
                    if (Me%Field(ip)%ComputeIntensity .or. Me%Field(ip)%Rotation .or.   &
                        Me%Field(ip)%Beaufort .or. Me%Field(ip)%ComputeDirection .or.   &
                        Me%Field(ip)%ComputeRotatedVector) then

                        call GetData(Me%Field(ip)%VectorX,                              &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'VECTOR_X',                         &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR250'
                        endif                            

                    endif                        
                    
                    if (Me%Field(ip)%ComputeIntensity .or. Me%Field(ip)%Rotation .or.   &
                        Me%Field(ip)%Beaufort .or. Me%Field(ip)%ComputeDirection .or.   &
                        Me%Field(ip)%ComputeRotatedVector) then                  

                        call GetData(Me%Field(ip)%VectorY,                              &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'VECTOR_Y',                         &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR260'
                        endif                            
                                        
                                        
                    endif                  
                    
                    if (Me%Field(ip)%Rotation .or. Me%Field(ip)%ComputeRotatedVector) then
                    
                        call GetData(Me%Field(ip)%VectorComponent,                      &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'COMPONENT',                        &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR270'
                        endif                            
                    
                    endif                      

                    if (Me%Field(ip)%ComputeRotatedVector) then                    

                        call GetData(Me%Field(ip)%GridRotation,                         &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'GRID_ROTATION',                    &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR275'
                        endif                           
                        
                    endif

                    call GetData(Me%Field(ip)%CenterX,                              &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlockInBlock,                   &
                                 keyword      = 'CENTER_X',                         &
                                 default      = .false.,                            &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR280'

                    call GetData(Me%Field(ip)%CenterY,                              &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlockInBlock,                   &
                                 keyword      = 'CENTER_Y',                         &
                                 default      = .false.,                            &                                 
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR290'
                                        

                    call GetData(Me%Field(ip)%Limit,                                &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlockInBlock,                   &
                                 keyword      = 'LIMIT',                            &
                                 default      = -1000.,                             &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR295'



                    call GetData(Me%Field(ip)%ComputeRH,                            &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlockInBlock,                   &
                                 keyword      = 'RELATIVE_HUMIDITY',                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR300'

                    if (Me%Field(ip)%ComputeRH) then

                        call GetData(Me%Field(ip)%TempRH,                               &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'TEMPERATURE_RH',                   &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR310'

                        call GetData(Me%Field(ip)%PressureRH,                           &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'PRESSURE_RH',                      &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR320'
                        

                        call GetData(Me%Field(ip)%SpecificHumidityRH,                   &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'SPECIFIC_HUMIDITY_RH',             &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2640'

                    endif                        

                    call GetData(Me%Field(ip)%FromMeteo2Algebric,                       &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'METEO_TO_ALGEBRIC',                    &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2650'
                    
                    call GetData(Me%Field(ip)%FromCartesian2Meteo,                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'CARTESIAN_TO_METEO',                   &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2655'
                    
                    if (Me%Field(ip)%FromMeteo2Algebric .and. Me%Field(ip)%FromCartesian2Meteo) then
                        write(*,*) "Options METEO_TO_ALGEBRIC and CARTESIAN_TO_METEO can not be both true"
                        stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2657'
                    endif
                   

                    call GetData(Me%Field(ip)%ValueIn%diL,                              &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DI_LOWER',                             &
                                 default      = 0,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2660'

                    call GetData(Me%Field(ip)%ValueIn%diU,                              &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DI_UPPER',                             &
                                 default      = 0,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2670'

                    call GetData(Me%Field(ip)%ValueIn%djL,                              &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DJ_LOWER',                             &
                                 default      = 0,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2680'

                    call GetData(Me%Field(ip)%ValueIn%djU,                              &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'DJ_UPPER',                             &
                                 default      = 0,                                      &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2690'

                    call GetData(Me%Field(ip)%AverageInDepth,                           &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'AVERAGE_IN_DEPTH',                     &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2700'         
                    
                    if (Me%Field(ip)%AverageInDepth) then
                        
                        call GetData(Me%Field(ip)%AverageInDepthName,                   &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'AVERAGE_IN_DEPTH_NAME',            &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2710'
                        
                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2720'
                        endif   

                    endif           

                    call GetData(Me%Field(ip)%Accumulated2StepGFS,                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'STEP_ACCUMULATED_GFS',                 &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2730'
                    
                    call GetData(Me%Field(ip)%Accumulated2Step,                         &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'STEP_ACCUMULATED',                     &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2735'
                    
                    if (Me%Field(ip)%Accumulated2StepGFS) then
                        Me%Field(ip)%Accumulated2Step = .true.
                    endif


                    call GetData(Me%Field(ip)%Reflectivity2Precipitation,               &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'RFLECTIVITY_2_PRECIPITATION',          &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2740'
                    
                    if (Me%Field(ip)%Reflectivity2Precipitation) then
                        
                        call GetData(Me%Field(ip)%ReflectivityName,                     &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'REFLECTIVITY_NAME',                &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2750'   

                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2760'
                        endif                           

                    endif           
                    
                    call GetData(Me%Field(ip)%Energy2Power,                             &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'ENERGY_2_POWER',                       &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2780'
                    
                    if (Me%Field(ip)%Energy2Power) then
                        
                        call GetData(Me%Field(ip)%EnergyName,                           &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'ENERGY_NAME',                      &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2790'   

                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2800'
                        endif   
                        
                        call GetData(Aux1DInt,                                          &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'ENERGY_STARTING_HOURS',            &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2810'   
                        endif
                        
                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2820'
                        else
                            Me%Field(ip)%N_EnergyStartH = iflag
                            allocate(Me%Field(ip)%EnergyStartingHours(1:iflag))
                            Me%Field(ip)%EnergyStartingHours(1:iflag) = Aux1DInt(1:iflag)
                        endif                           

                    endif        
                    
                    call GetData(Me%Field(ip)%AvModelStart2Inst,                        &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'AVERAGE_MODEL_START_2_INST',           &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2830'
                    
                    if (Me%Field(ip)%AvModelStart2Inst) then
                        
                        call GetData(Me%Field(ip)%AvModelStartName,                     &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'AVERAGE_MODEL_START_NAME',         &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2840'   

                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2850'
                        endif   
                        
                        !call GetData(Aux1DInt,                                          &
                        !             Me%ObjEnterData, iflag,                            &
                        !             SearchType   = FromBlockInBlock,                   &
                        !             keyword      = 'AVERAGE_MODEL_STARTING_HOURS',     &
                        !             ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                        !             STAT         = STAT_CALL)       
                        !if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_) then
                        !    stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2860'   
                        !endif
                        !
                        !if (iflag == 0) then
                        !    stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2870'
                        !else
                        !    Me%Field(ip)%N_AvModelStartH = iflag
                        !    allocate(Me%Field(ip)%AvModelStartingHours(1:iflag))
                        !    Me%Field(ip)%AvModelStartingHours(1:iflag) = Aux1DInt(1:iflag)
                        !endif                           

                    endif                      
                    
                    call GetData(Me%Field(ip)%VerticalZ_2D,                             &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'VERTICAL_Z_2D',                        &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2880'
                    
                    if (Me%Field(ip)%VerticalZ_2D) then
                        
                        call GetData(Me%Field(ip)%Bottom,                               &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'BOTTOM_NAME',                      &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2890'   
                        endif
                        
                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2900'
                        endif   
                        
                        call GetData(Me%Field(ip)%Surface,                              &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'SURFACE_NAME',                     &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_ .and. STAT_CALL /= SIZE_ERR_) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2910'   
                        endif
                        
                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2920'
                        endif                           

                    endif        

                    call GetData(Me%Field(ip)%Wfp,                                      &
                                 Me%ObjEnterData, iflag,                                &
                                 SearchType   = FromBlockInBlock,                       &
                                 keyword      = 'WAVE_PF_To_WAVE_TP',                   &
                                 default      = .false.,                                &
                                 ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',           &
                                 STAT         = STAT_CALL)       
                    if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2930'         
                    
                    if (Me%Field(ip)%Wfp) then
                        
                        call GetData(Me%Field(ip)%WfpName,                              &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlockInBlock,                   &
                                     keyword      = 'WAVE_PF_To_WAVE_TP_NAME',          &
                                     ClientModule = 'ModuleNetCDFCF_2_HDF5MOHID',       &
                                     STAT         = STAT_CALL)       
                        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2940'
                        
                        if (iflag == 0) then
                            stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR2950'
                        endif   

                    endif 

                else BF
                
                    stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR3000'    
                    
                endif BF
            else IS
            
                stop 'ReadFieldOptions - ModuleNetCDFCF_2_HDF5MOHID - ERR4000'    
            
            endif IS            
        
        enddo                 
                              
                              
    end subroutine ReadFieldOptions                              


    !------------------------------------------------------------------------

    subroutine ReadNetCDFCF_WriteHDF5MOHID

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=PathLength)               :: InPutFile
        logical                                 :: BlockFound, LastLineON
        integer                                 :: iflag, line, FirstLine, LastLine,    &
                                                   STAT_CALL, HDF5_CREATE, iOut, status,&
                                                   ncid, iP, i
                 


        !Begin----------------------------------------------------------------
        
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        if (Me%OutHDF5) then
            !Opens HDF5 File
            call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Me%OutNetCDF) then
            
            call OpenNCDFFile

        endif
        
        Me%ReadPropNumber = Me%PropNumber
       
        do iP = 1, Me%PropNumber
                                                                                                    
            if (Me%Field(iP)%ComputeIntensity .or. Me%Field(iP)%Rotation                   .or. &
                Me%Field(iP)%Beaufort         .or. Me%Field(iP)%ComputeRH                  .or. &
                Me%Field(iP)%WaveBeaufort     .or. Me%Field(iP)%ComputeDirection           .or. &
                Me%Field(iP)%AverageInDepth   .or. Me%Field(iP)%Reflectivity2Precipitation .or. &
                Me%Field(iP)%Energy2Power     .or. Me%Field(iP)%AvModelStart2Inst          .or. &
                Me%Field(iP)%Wfp              .or. Me%Field(iP)%ComputeRotatedVector  ) then
            
                Me%ReadPropNumber = Me%ReadPropNumber - 1
            
            endif
        
        enddo
        
OF:     if(Me%OnlyOneFile)then
    
            iOut                = 0
            Me%Date%TotalInst   = 0
            Me%OutCountProp     = 0
            call null_time(Me%Date%LasTEndTime)
            
            Me%ReadingFileName  = Me%FileName
                    
            !Verifies if file exists
            status=NF90_OPEN(trim(Me%FileName),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) then
                write(*,*) 'Not able to open file=',trim(Me%FileName)
                stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            endif
                    
            write(*,*) 'Reading ', trim(Me%FileName)
        
            !Time
            call ReadTotalTime          (ncid)    
                    
            write(*,*) "Read total time"  

                        
            !Grid/Bathym and Grid/Longitude and Grid/Latitude
            call ReadWriteGrid2D  (ncid) 
            !Bathym, mapping
            call ReadWriteGrid3D  (ncid)

            write(*,*) "Grid"                                                   
                    
            allocate(Me%Date%InstOut(1:Me%Date%TotalInst))

            Me%Date%InstOut(:)   = -99
            Me%Date%TotalInstOut = 1
                        
            Me%Date%InstOut(1) = 1
            do i=2, Me%Date%TotalInst
                if (Me%Date%ValueInTotal(i) > Me%Date%ValueInTotal(i-1)) then
                    Me%Date%TotalInstOut = Me%Date%TotalInstOut + 1
                    Me%Date%InstOut(i)   = Me%Date%TotalInstOut
                endif
            enddo                            
                    
            if (Me%OutHDF5  ) call WriteTimeHDF5  
            if (Me%OutNetCDF) call WriteTimeNetCDF
              
            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            iOut              = 0
            call null_time(Me%Date%LasTEndTime)
                
            !Verifies if file exists
            status=NF90_OPEN(trim(Me%FileName),NF90_NOWRITE,ncid)
            if (status /= nf90_noerr) stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                    
            write(*,*) 'Reading ', trim(Me%FileName)

            if (Me%OutCountProp == 0) then
                
                call ReadTimeNetCDF(ncid)
                call deallocatevaluein(Me%Date%ValueIn)
                !Grid/Vertical Z
                if (Me%Depth%Dim3D) then
                    call WriteDepth    (iOut, ncid)
                    write(*,*) 'Write depth'
                endif
                    
            endif
                    
            !Results
            write(*,*) 'Read Write fields'  
                    
            call ReadWriteFields    (ncid, iOut) 
                                        
            if (Me%OutCountProp == Me%Date%NumberInst * Me%ReadPropNumber) then
                    
                iOut = iOut + Me%Date%NumberInst
                        
                Me%OutCountProp = 0

            endif
                    
            status=NF90_CLOSE(ncid)
            if (status /= nf90_noerr)  stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
        else OF
            
                           
            call ExtractBlockFromBlock(Me%ObjEnterData, Me%ClientNumber,                    &
                                       input_files_begin, input_files_end,                  &
                                       BlockInBlockFound = BlockFound,                      &
                                       FirstLine = FirstLine, LastLine = LastLine,          &
                                       STAT = STAT_CALL)

IS:         if(STAT_CALL .EQ. SUCCESS_) then

                !The block is found to exist before when reading depth
BF:             if (BlockFound) then            

                    iOut              = 0

                    Me%Date%TotalInst = 0
                    Me%OutCountProp   = 0
                    call null_time(Me%Date%LasTEndTime)
        
                    do line = FirstLine + 1, LastLine - 1

                        call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                                     Buffer_Line = line, STAT = STAT_CALL)

                        Me%ReadingFileName=InputFile
                    
                        !Verifies if file exists
                        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
                        if (status /= nf90_noerr) then
                            write(*,*) 'Not able to open file=',trim(InputFile)
                            stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                        endif
                    
                        write(*,*) 'Reading ', trim(InputFile)
                    
                        if (line == LastLine - 1) then
                            LastLineON = .true.
                        else
                            LastLineON = .false.
                        endif

                        !Time
                        call ReadTotalTime          (ncid)    
                    
                        write(*,*) "Read total time"  

                        if (line == FirstLine + 1) then
                            !Grid/Bathym and Grid/Longitude and Grid/Latitude
                            call ReadWriteGrid2D  (ncid) 
                            !Bathym, mapping
                            call ReadWriteGrid3D  (ncid)

                            write(*,*) "Grid"                                                   
                        endif 

                        if (LastLineON) then
                    
                            allocate(Me%Date%InstOut(1:Me%Date%TotalInst))

                            Me%Date%InstOut(:)   = -99
                            Me%Date%TotalInstOut = 1
                        
                            Me%Date%InstOut(1) = 1
                            do i=2, Me%Date%TotalInst
                                if (Me%Date%ValueInTotal(i) > Me%Date%ValueInTotal(i-1)) then
                                    Me%Date%TotalInstOut = Me%Date%TotalInstOut + 1
                                    Me%Date%InstOut(i)   = Me%Date%TotalInstOut
                                endif
                            enddo                            
                    
                            if (Me%OutHDF5  ) call WriteTimeHDF5  
                            if (Me%OutNetCDF) call WriteTimeNetCDF
                        
                        endif
              
                        status=NF90_CLOSE(ncid)
                        if (status /= nf90_noerr)  stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

                
                    enddo

                    iOut              = 0
              
                    call null_time(Me%Date%LasTEndTime)
                
                    do line = FirstLine + 1, LastLine - 1

                        call GetData(InputFile, EnterDataID = Me%ObjEnterData, flag = iflag,    &
                                     Buffer_Line = line, STAT = STAT_CALL)

                    
                        !Verifies if file exists
                        status=NF90_OPEN(trim(InputFile),NF90_NOWRITE,ncid)
                        if (status /= nf90_noerr) stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                    
                        write(*,*) 'Reading ', trim(InputFile)

                        !Time
                        if (Me%OutCountProp == 0) then
                            call ReadTimeNetCDF(ncid)
                            call deallocatevaluein(Me%Date%ValueIn)
                            !Grid/Vertical Z
                            if (Me%Depth%Dim3D) then
                                call WriteDepth    (iOut, ncid)
                            
                                write(*,*) 'Write depth'
                            endif
                        
                        
                    
                        endif
                    
                        !Results/XXX  
                    
                        write(*,*) 'Read Write fields'  
                    
                        call ReadWriteFields    (ncid, iOut) 
                    

                                        
                        if (Me%OutCountProp == Me%Date%NumberInst * Me%ReadPropNumber) then
                    
                            iOut = iOut + Me%Date%NumberInst
                        
                            Me%OutCountProp = 0

                        endif
                    
                        status=NF90_CLOSE(ncid)
                        if (status /= nf90_noerr)  stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

                
                    enddo

                
                else  BF
                    stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
                endif BF

                call Block_Unlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL) 

                if (STAT_CALL /= SUCCESS_)                                                  &
                    stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR70'

            else   IS

                stop 'ReadNetCDFCF_WriteHDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR80'

            end if IS

            
        endif OF


    end subroutine ReadNetCDFCF_WriteHDF5MOHID
    
    !------------------------------------------------------------------------
       
    subroutine OpenNCDFFile

        !Local-----------------------------------------------------------------
        integer                                     :: NCDF_CREATE, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call GetNCDFFileAccess(NCDF_CREATE = NCDF_CREATE)
        
        call ConstructNETCDF(NCDFID = Me%NetCDF_Out%ObjNetCDF,                          &
                             FileName  = trim(Me%NetCDF_Out%Name),                      &
                             Access    = NCDF_CREATE,                                   &
                             STAT      = STAT_CALL)            
        if (STAT_CALL /= SUCCESS_)stop 'OpenNCDFFile - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        

        call NETCDFWriteHeader(NCDFID         = Me%NetCDF_Out%ObjNetCDF     ,           &
                               Title          = Me%NetCDF_Out%Title         ,           &
                               Convention     = Me%NetCDF_Out%Convention    ,           &
                               Version        = Me%NetCDF_Out%Version       ,           &
                               History        = Me%NetCDF_Out%History       ,           &
                               iDate          = Me%NetCDF_Out%iDate         ,           &
                               Source         = Me%NetCDF_Out%Source        ,           &
                               Institution    = Me%NetCDF_Out%Institution   ,           &
                               References     = Me%NetCDF_Out%References    ,           &
                               STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OpenNCDFFile - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        
        write(*,*)
        write(*,*)'Opened ncdf file                : ', trim(Me%NetCDF_Out%Name)

    end subroutine OpenNCDFFile
   !------------------------------------------------------------------------
    subroutine ReadTotalTime(ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: ncid
        
        !Local-----------------------------------------------------------------
        real(8), dimension(:), pointer                  :: AuxT
        integer                                         :: iAux, j
        !Begin-----------------------------------------------------------------
        
        call ReadTimeNetCDF(ncid)
        
        if (Me%Date%FileEndTime > Me%Date%LasTEndTime) then
  
            if(associated(Me%Date%ValueInTotal)) then
                iAux = size(Me%Date%ValueInTotal)
            else 
                iAux = 0
            endif
            
            allocate(AuxT(iAux + Me%Date%NumberInst))
            
            if(associated(Me%Date%ValueInTotal)) then
                AuxT(1:iAux) = Me%Date%ValueInTotal(:)
            endif
            
            do j=1, Me%Date%NumberInst

                AuxT(j+iAux) =  GetNetCDFValue(Me%Date%ValueIn, Dim1 = j) 
                !To take in consideration the ROMS MeteoGalicia case where the final instants have value 0 and need to be discard.
                if (Me%Date%NumberInst >1 .and. j >1) then
                    if (AuxT(j+iAux)< AuxT(j-1+iAux)) then
                        Me%Date%NumberInst = j - 1
                        write(*,*)
                        write(*,*)'Warning : temporal inflexion point detected in file ', trim(Me%ReadingFileName)
                        write(*,*)
                        exit
                    endif
                endif

            enddo
            
            
            
            iAux = Me%Date%NumberInst+iAux
            
            allocate(Me%Date%ValueInTotal(iAux))
            
            Me%Date%ValueInTotal(1:iAux) = AuxT(1:iAux)
            
            deallocate(AuxT)

            Me%Date%TotalInst   = Me%Date%TotalInst + Me%Date%NumberInst
            Me%Date%LasTEndTime = Me%Date%FileEndTime            
        endif
        
        call DeAllocateValueIn(Me%Date%ValueIn)            
        

        
        
    end subroutine ReadTotalTime

    !------------------------------------------------------------------------
    
   !------------------------------------------------------------------------
    subroutine ImposedRegularGrid(ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: ncid
        
        !Local-----------------------------------------------------------------
        integer                                         :: status, numDims, i, j
        integer                                         :: RhVarIdLong
        integer, dimension(nf90_max_var_dims)           :: rhDimIdsLong
        
        !Begin-----------------------------------------------------------------        


        
        status=nf90_inq_varid(ncid,trim(Me%LongLat%NetCDFNameLong),RhVarIdLong)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        status = nf90_inquire_variable(ncid, RhVarIdLong, ndims = numDims)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        status = nf90_inquire_variable(ncid, RhVarIdLong, dimids = rhDimIdsLong(:numDims))
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        
        status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(1), len = Me%LongLat%jmax)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'  

        status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(2), len = Me%LongLat%imax)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'         
        
        if (Me%ReadInvertXY) then

            status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(2), len = Me%LongLat%jmax)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'        

            status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(1), len = Me%LongLat%imax)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR70' 
        
        endif
        
    
        !Build HDF5 MOHID Grid
        Me%WorkSize%ILB = Me%LongLat%dij
        Me%WorkSize%IUB = Me%LongLat%imax - 1 - Me%LongLat%dij
        Me%WorkSize%JLB = Me%LongLat%dij
        Me%WorkSize%JUB = Me%LongLat%jmax - 1 - Me%LongLat%dij

        
        !to warn the user before the model crashes
        !cant use a NetCDF with one of the dimension as 2 (or lower) because IUB or JUB would be zero (or lower).
        if ((Me%WorkSize%IUB < 1) .or. (Me%WorkSize%JUB < 1)) then
            write (*,*)
            write (*,*) 'Please use a NETCDF file with more than'
            write (*,*) '2x2 points so that the grid can be correctly extracted'
            stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
        endif

        Me%Size%ILB     = Me%WorkSize%ILB - 1
        Me%Size%IUB     = Me%WorkSize%IUB + 1
        Me%Size%JLB     = Me%WorkSize%JLB - 1
        Me%Size%JUB     = Me%WorkSize%JUB + 1    
        
        allocate(Me%LongLat%LongOut   (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%LongLat%LatOut    (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
      
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB+1
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB+1
            
            Me%LongLat%LongOut(i, j) = Me%LongLat%LongOrig + real(j-1)*Me%LongLat%dLong
            Me%LongLat%LatOut (i, j) = Me%LongLat%LatOrig  + real(i-1)*Me%LongLat%dLat
            
        enddo
        enddo     
            
    end subroutine ImposedRegularGrid                       

   !------------------------------------------------------------------------
    subroutine ReadWriteGrid2D(ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: ncid
        
        !Local-----------------------------------------------------------------
        real,  dimension(:,:), pointer              :: Lat, Long
        real,   dimension(:),   pointer             :: Dummy
        integer                                     :: STAT_CALL, i, j
        type (T_Size2D)                             :: WorkSize2D
        integer                                     :: WorkILB, WorkIUB, WorkJLB, WorkJUB

        !Begin-----------------------------------------------------------------
        
        if (Me%LongLat%Imposed) then
            call ImposedRegularGrid(ncid)
        else
            call ReadGrid2DNetCDF(ncid)
        endif            

        Me%WorkSize2D%ILB = Me%WorkSize%ILB
        Me%WorkSize2D%IUB = Me%WorkSize%IUB        
        Me%WorkSize2D%JLB = Me%WorkSize%JLB
        Me%WorkSize2D%JUB = Me%WorkSize%JUB        
        


        if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
            
            WorkSize2D%ILB = WorkILB
            WorkSize2D%IUB = WorkIUB            
            WorkSize2D%JLB = WorkJLB
            WorkSize2D%JUB = WorkJUB            
            
        else
        
            WorkILB = Me%WorkSize%ILB
            WorkIUB = Me%WorkSize%IUB        
            WorkJLB = Me%WorkSize%JLB
            WorkJUB = Me%WorkSize%JUB       
            
            Me%WindowOut%ILB = Me%WorkSize%ILB
            Me%WindowOut%IUB = Me%WorkSize%IUB
            Me%WindowOut%JLB = Me%WorkSize%JLB
            Me%WindowOut%JUB = Me%WorkSize%JUB

            WorkSize2D = Me%WorkSize2D 
        endif
     

        if (Me%WindowOut%ON) then

            allocate(Long(WorkILB-1:WorkIUB+1, WorkJLB-1:WorkJUB+1))
            allocate(Lat (WorkILB-1:WorkIUB+1, WorkJLB-1:WorkJUB+1))            

            do j = WorkJLB, WorkJUB+1
            do i = WorkILB, WorkIUB+1                
                Lat (i,j) = Me%LongLat%LatOut (i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
                Long(i,j) = Me%LongLat%LongOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)                
            enddo
            enddo
        else            
            Lat  => Me%LongLat%LatOut
            Long => Me%LongLat%LongOut
        endif

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Lat, Long, &
                                     XX  = Dummy, YY = Dummy, Latitude = 45., Longitude = -8.,    &
                                     ILB = WorkILB, IUB = WorkIUB,                &
                                     JLB = WorkJLB, JUB = WorkJUB,                &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadWriteGrid2D - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        if (Me%OutHDF5) then


            call WriteHorizontalGrid(Me%ObjHorizontalGrid, Me%ObjHDF5,                         &
                                     WorkSize = WorkSize2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadWriteGrid2D - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            
        endif
        
        if (Me%WindowOut%ON) then

            deallocate(Long)
            deallocate(Lat )          
        
        endif
        
        nullify(Long)
        nullify(Lat )          
        
    end subroutine ReadWriteGrid2D

    !------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine ReadWriteGrid3D(ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: ncid
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        call ReadBathymNetCDF(ncid)

        if (Me%OutNetCDF) then
            call WriteGrid2DNetCDF
        endif
        
        if (.not. Me%Bathym%FromMapping) then
            if (Me%OutHDF5  ) call WriteBathymHDF5        
            if (Me%OutNetCDF) call WriteBathymNetCDF       
        endif
        

        
        call ReadGrid3DNetCDF(ncid)
        
        if (Me%OutHDF5) call WriteGrid3DHDF5        
        
        if (Me%OutNetCDF) call WriteGridNetCDF  
              


            
        if (associated(Me%LongLat%LongOut    )) deallocate(Me%LongLat%LongOut   )
        if (associated(Me%LongLat%LatOut     )) deallocate(Me%LongLat%LatOut    )
        
        if (associated(Me%Depth%Value3DOut   )) deallocate(Me%Depth%Value3DOut  )

    end subroutine ReadWriteGrid3D

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine WriteGridNetCDF
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer, dimension(:,:,:), pointer  :: mask3D
        integer, dimension(:,:  ), pointer  :: mask2D        
        character(len=StringLength)         :: NCDFName, LongName, StandardName, Units
        real                                :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                             :: i, j, k, STAT_CALL
     
        !Begin-----------------------------------------------------------------



        
        if (Me%Depth%Dim3D) then
        
            allocate(mask3D(1:Me%WorkSize%JUB, 1:Me%WorkSize%IUB, 1:Me%WorkSize%KUB))
            
            do k=1,  Me%WorkSize%KUB
            do j=1,  Me%WorkSize%JUB
            do i=1,  Me%WorkSize%IUB            
                mask3D(j,i,k) = Me%Mapping%Value3DOut(i,j,k)
            enddo
            enddo
            enddo
            
            if (Me%MeridionalSplit) then
                do k=1,  Me%WorkSize%KUB
                do i=1,  Me%WorkSize%IUB            
                
                    do j=1,  Me%MeridionalSplitColumn
                        mask3D(j,i,k) = Me%Mapping%Value3DOut(i,Me%MeridionalSplitColumn+j,k)
                    enddo
                    
                    do j= Me%MeridionalSplitColumn+1,Me%WorkSize%JUB
                        mask3D(j,i,k) = Me%Mapping%Value3DOut(i,j-Me%MeridionalSplitColumn,k)
                    enddo
                
                enddo                
                enddo            
            endif
            
            call BuildAttributes("WaterPoints3D", NCDFName,                             &
                                 LongName, StandardName,                                &
                                 Units, ValidMin, ValidMax,                             &
                                 MinValue, MaxValue, MissingValue,                      &
                                 Int3D = mask3D)
    
            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,              &
                                  Name          = trim(NCDFName),                       &
                                  LongName      = trim(LongName),                       &
                                  StandardName  = trim(StandardName),                   & 
                                  Units         = trim(Units),                          &
                                  ValidMin      = ValidMin,                             &
                                  ValidMax      = ValidMax,                             &
                                  MinValue      = MinValue,                             &
                                  MaxValue      = MaxValue,                             &
                                  MissingValue  = MissingValue,                         &
                                  Array3D       = mask3D,                               &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
            
            deallocate(mask3D)

        else


            allocate(mask2D(1:Me%WorkSize%JUB, 1:Me%WorkSize%IUB))
            
            do j=1,  Me%WorkSize%JUB
            do i=1,  Me%WorkSize%IUB            
                mask2D(j,i) = Me%Mapping%Value2DOut(i,j)
            enddo
            enddo



            call BuildAttributes("WaterPoints2D", NCDFName, LongName, StandardName,     &
                                 Units, ValidMin, ValidMax,                             &
                                 MinValue, MaxValue, MissingValue,                      &
                                 Int2D = mask2d)

    
            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,              &
                                  Name          = trim(NCDFName),                       &
                                  LongName      = trim(LongName),                       &
                                  StandardName  = trim(StandardName),                   & 
                                  Units         = trim(Units),                          &
                                  ValidMin      = ValidMin,                             &
                                  ValidMax      = ValidMax,                             &
                                  MinValue      = MinValue,                             &
                                  MaxValue      = MaxValue,                             &
                                  MissingValue  = MissingValue,                         &
                                  Array2D       = mask2D,                               &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            deallocate(mask2D)

        endif
    
    
    end subroutine WriteGridNetCDF

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine WriteGrid2DNetCDF
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer    :: lon, lat, lon_stag, lat_stag
        real(8), dimension(:,:), pointer    :: SphericMercatorX_stag, SphericMercatorY_stag
        integer                             :: i, j, STAT_CALL
     
        !Begin-----------------------------------------------------------------
        
        
        allocate(lon     (Me%WorkSize%JUB,  Me%WorkSize%IUB  ), lat     (Me%WorkSize%JUB,  Me%WorkSize%IUB    ))
        allocate(lon_stag(Me%WorkSize%JUB+1,Me%WorkSize%IUB+1), lat_stag(Me%WorkSize%JUB+1,Me%WorkSize%IUB+1))

        allocate(SphericMercatorX_stag(Me%WorkSize%JUB+1,Me%WorkSize%IUB+1),            &
                 SphericMercatorY_stag(Me%WorkSize%JUB+1,Me%WorkSize%IUB+1))
        
        do j=1,  Me%WorkSize%JUB+1
        do i=1,  Me%WorkSize%IUB+1            
        
            lon_stag(j,i) = Me%LongLat%LongOut(i,j)
            lat_stag(j,i) = Me%LongLat%LatOut (i,j)
        
        enddo
        enddo
        
        do j=1,  Me%WorkSize%JUB
        do i=1,  Me%WorkSize%IUB            
        
            lon(j,i) = (lon_stag(j,i) + lon_stag(j+1,i) + lon_stag(j,i+1) + lon_stag(j+1,i+1))/4.
            lat(j,i) = (lat_stag(j,i) + lat_stag(j+1,i) + lat_stag(j,i+1) + lat_stag(j+1,i+1))/4.
        
        enddo
        enddo

        
        if (Me%Depth%Dim3D) then
            call NETCDFSetDimensions (Me%NetCDF_Out%ObjNetCDF,                          &
                                                    IUB  = Me%WorkSize%IUB,             &
                                                    JUB  = Me%WorkSize%JUB,             &
                                                    KUB  = Me%WorkSize%KUB,             &
                                                    STAT = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        else
            call NETCDFSetDimensions (Me%NetCDF_Out%ObjNetCDF,                          &
                                                    IUB  = Me%WorkSize%IUB,             &
                                                    JUB  = Me%WorkSize%JUB,             &
                                                    KUB  = 0,                           &
                                                    STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        endif

        call WGS84toGoogleMaps(lon_stag, lat_stag, Me%WorkSize2D%JLB, Me%WorkSize2D%JUB,&
                               Me%WorkSize2D%ILB, Me%WorkSize2D%IUB,                    &
                               SphericMercatorX_stag, SphericMercatorY_stag)

                
        call NETCDFWriteLatLon(Me%NetCDF_Out%ObjNetCDF, Lat, Lon, Lat_Stag, Lon_Stag,   &
                               SphericMercatorX_stag, SphericMercatorY_stag,            &
                               STAT = STAT_CALL)
                                
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

        deallocate(lon     , lat     )
        deallocate(lon_stag, lat_stag)
        deallocate(SphericMercatorX_stag, SphericMercatorY_stag)

        
    
    end subroutine WriteGrid2DNetCDF

    !------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine WriteBathymNetCDF
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer    :: bathym
        character(len=StringLength)         :: NCDFName, LongName, StandardName, Units, Positive
        real                                :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        integer                             :: i, j, STAT_CALL
     
        !Begin-----------------------------------------------------------------
        
        allocate(bathym(Me%WorkSize%JUB,Me%WorkSize%IUB))
        
        do j=1,  Me%WorkSize%JUB
        do i=1,  Me%WorkSize%IUB
            bathym(j,i) = Me%Bathym%Value2DOut(i,j)
        enddo
        enddo        
        
        call BuildAttributes("Bathymetry", NCDFName, LongName, StandardName,            &
                                           Units, ValidMin, ValidMax,                   &
                                           MinValue, MaxValue, MissingValue, Positive, bathym)
        
        call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,                  &
                              Name          = trim(NCDFName),                           &
                              LongName      = trim(LongName),                           &
                              StandardName  = trim(StandardName),                       & 
                              Units         = trim(Units),                              &
                              ValidMin      = ValidMin,                                 &
                              ValidMax      = ValidMax,                                 &
                              MinValue      = MinValue,                                 &
                              MaxValue      = MaxValue,                                 &
                              MissingValue  = MissingValue,                             &
                              Array2D       = bathym,                                   &
                              STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGridNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        deallocate(bathym)
        
    
    end subroutine WriteBathymNetCDF

    !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine ReadWriteFields(ncid, iOut)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut, ncid
        
        !Local-----------------------------------------------------------------
        integer                                         :: iP, iT, iFinal
        logical                                         :: WriteProp
        !Begin-----------------------------------------------------------------
        
        WriteProp = .false. 

        do iP = 1, Me%PropNumber
            do iT =1, Me%Date%NumberInst
        
                if (Me%Field(iP)%ComputeIntensity           .or. Me%Field(iP)%Rotation       .or.  &
                    Me%Field(iP)%Beaufort                   .or. Me%Field(iP)%WaveBeaufort   .or.  &
                    Me%Field(iP)%ComputeDirection           .or. Me%Field(iP)%AverageInDepth .or.  &
                    Me%Field(iP)%Reflectivity2Precipitation .or. Me%Field(iP)%Energy2Power   .or.  &
                    Me%Field(iP)%AvModelStart2Inst          .or. Me%Field(iP)%Wfp            .or.  &
                    Me%Field(iP)%ComputeRotatedVector) then
                
                    WriteProp       = .false.
                
                else
            
                    call ReadFieldNetCDF(ncid, iP, WriteProp, iT)
                    
                   
                endif
                
                if (WriteProp) then
                
                    iFinal = iOut + iT
                    if (Me%Date%InstOut(iFinal) >0) then 
                        call WriteFieldAllInst(Me%Date%InstOut(iFinal), iP, iT)
                    endif
                    
                    Me%OutCountProp = Me%OutCountProp + 1
                endif

                if (.not. (Me%Field(iP)%ComputeIntensity            .or. Me%Field(iP)%Rotation .or.     &
                           Me%Field(iP)%Beaufort                    .or. Me%Field(iP)%WaveBeaufort .or. &
                           Me%Field(iP)%ComputeDirection            .or. Me%Field(iP)%AverageInDepth.or.&
                           Me%Field(iP)%Reflectivity2Precipitation  .or. Me%Field(iP)%Energy2Power  .or.&
                           Me%Field(iP)%AvModelStart2Inst           .or. Me%Field(iP)%Wfp           .or.&
                           Me%Field(iP)%ComputeRotatedVector))  &
                    call DeAllocateValueIn(Me%Field(iP)%ValueIn)
            enddo
        enddo
        

   end subroutine ReadWriteFields        
                  
   !------------------------------------------------------------------------


   !------------------------------------------------------------------------
    subroutine WriteDepth(iOut, ncid)
        !Arguments-------------------------------------------------------------
        integer                                         :: iOut, ncid
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer              :: Vert3D
        real(8), dimension(:), pointer                  :: DepthAux
        real,    dimension(:), pointer                  :: DepthOut, DepthOutStag
        real(8)                                         :: Aux, DepthC, DepthC_Below
        integer                                         :: iFinal, i, j, k, iT, STAT_CALL, kin
        integer                                         :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        logical                                         :: SigmaIn
        real(8)                                         :: SumDepth
        real                                            :: Topdepth
        logical                                         :: Method2 = .false.

        !Begin-----------------------------------------------------------------

        allocate(Me%Depth%Value3DOut(Me%Size%ILB:Me%Size%IUB,           &
                                     Me%Size%JLB:Me%Size%JUB,           &
                                     Me%Size%KLB:Me%Size%KUB))

        
        
d1:     do iT =1, Me%Date%NumberInst
        
            iFinal = Me%Date%InstOut(iOut + iT)
            
            if (iFinal <0) cycle

i1:         if (Me%OutNetCDF) then

i2:             if (iFinal ==1) then
                    allocate(DepthOut    (Me%WorkSize%KLB:Me%WorkSize%KUB  ))
                    allocate(DepthOutStag(Me%WorkSize%KLB:Me%WorkSize%KUB+1))                                

i3:                 if (Me%Depth%Interpolate) then

                        DepthOut(:) =  Me%Depth%ZLevels(:)

                    else i3
d2:                     do k= Me%WorkSize%KLB, Me%WorkSize%KUB 
                            if (Me%Depth%InvertLayers) then
                                kin = Me%WorkSize%KUB - k + Me%WorkSize%KLB
                            else
                                kin = k
                            endif

                            if (Me%Depth%GeoVert == Hybrid) stop 'WriteDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR10'   
                            
                            DepthOut(k) = GetNetCDFValue(Me%Depth%ValueIn,  Dim1 = kin)
                        enddo d2
                    endif i3
                    
                    DepthOutStag(Me%WorkSize%KUB+1) = 0.
                    
d3:                 do k= Me%WorkSize%KUB, Me%WorkSize%KLB,-1  
                        Aux = 2*DepthOut(k) - DepthOutStag(k+1)
                        if (Aux < DepthOut(k)) then
                            DepthOutStag(k) = Aux 
                        else
                            if (k>1) then
                                DepthOutStag(k) = (DepthOut(k) + DepthOut(k-1))/2.
                            else
                                DepthOutStag(k) = DepthOut(k) + 100. 
                            endif
                        endif
                    enddo d3
                


                    if (Me%Depth%GeoVert == sigma_) then
                        SigmaIn = .true. 
                    else
                        SigmaIn = .false. 
                    endif
                    call NETCDFWriteVert    (Me%NetCDF_Out%ObjNetCDF, DepthOut, SigmaIn, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR20' 
                    
                    call NETCDFWriteVertStag(Me%NetCDF_Out%ObjNetCDF, DepthOutStag, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR30' 

                    deallocate(DepthOut, DepthOutStag)
                    
                endif  i2         

            endif  i1
            
            if (Me%Depth%Interpolate) allocate(DepthAux(1:Me%Depth%kmax+1))                                          
            
if23:       if (Me%Depth%GeoVert == sigma_) then            
                call ReadWaterLevelNetCDF(ncid, iT)                        
            endif if23
            
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                SumDepth = 0.
                Method2 = .false.
                do k= Me%WorkSize%KUB, Me%WorkSize%KLB, -1
            
           
if12:           if (Me%Mapping%Value3DOut(i,j,k) == 1) then

if13:               if (Me%Depth%GeoVert == sigma_) then

                        TopDepth = - GetNetCDFValue(Me%Depth%WLValueIn, Dim1 = j+1, Dim2 = i+1, Dim3 = 1) 

if20:                   if (k==Me%WorkSize%KUB) then
                            Me%Depth%Value3DOut(i, j, k) = GetCellInDepth(i, j, k+1,Me%WorkSize%KUB,iT, &
                                                            CellFace = .true., Topdepth = Topdepth)
                            SumDepth                     = 0.
                        endif if20

                        Me%Depth%Value3DOut(i, j, k-1) = GetCellInDepth(i, j, k,Me%WorkSize%KUB,iT,     &
                                                            CellFace = .true., Topdepth = Topdepth)
                        
                    else if13
                
if14:                   if (k==Me%WorkSize%KUB) then
                            Me%Depth%Value3DOut(i, j, k) = 0.
                            SumDepth                     = 0.
                        endif if14
                        
if15:                   if (.not. Me%Depth%Interpolate) then

                            if (Me%WorkSize%KUB == 1) then
                                Me%Depth%Value3DOut(i, j, k  ) = 0.
                                Me%Depth%Value3DOut(i, j, k-1) = 1.
                                cycle
                            endif
                                                
                            DepthC = GetCellInDepth(i, j, k,Me%WorkSize%KUB,iT)
                            
if16:                       if (DepthC <0 .and. k==Me%WorkSize%KUB) then
if19:                           if (Me%Mapping%Value3DOut(i,j,k-1) == 1) then
                                    DepthC_Below = GetCellInDepth(i, j, k-1,Me%WorkSize%KUB,iT)
                                else if19
                                    DepthC_Below = 0.
                                endif if19
                                SumDepth = (DepthC_Below - DepthC)/2. - DepthC
                            endif if16
                            
if17:                       if (SumDepth > 0) then
                                DepthC = DepthC + SumDepth
                                Aux    = DepthC
                            else if17
                                Aux  = 2. * DepthC - Me%Depth%Value3DOut(i, j, k)
                                
                                if (k==Me%WorkSize%KUB) then
                                    if (Me%Mapping%Value3DOut(i,j,k-1) == 1) then
                                        DepthC_Below = GetCellInDepth(i, j, k-1,Me%WorkSize%KUB,iT)
                                        if (Aux >= DepthC_Below) then
                                            Aux = (DepthC + DepthC_Below)/2.
                                        endif
                                    endif
                                endif                                    
                                                     
                            endif if17
                            
                            Me%Depth%Value3DOut(i, j, k-1) = Aux
                            
                            if (k > 1) then
                            
                                DepthC_Below = GetCellInDepth(i, j, k-1,Me%WorkSize%KUB,iT)
                                if (Aux >= DepthC_Below .or. Method2) then
                                    Me%Depth%Value3DOut(i, j, k-1) = (DepthC + DepthC_Below)/2.                                 
                                    Method2 = .true. 
                                endif 

                            endif
                            
                        else if15
                            Me%Depth%Value3DOut(i, j, k-1) = 2*Me%Depth%ZLevels(k) - Me%Depth%Value3DOut(i, j, k)   
                        endif if15
                                                                              
                    endif if13
                 
                 else if12
                    if (k==Me%WorkSize%KUB) then
                        Me%Depth%Value3DOut(i, j, k) = FillValueReal
                    endif                        
                    Me%Depth%Value3DOut(i, j, k-1) = FillValueReal   
                endif if12
                                                   
                enddo
            enddo
            enddo
            
            if (abs(Me%Depth%Offset)>0.) then
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                do k= Me%WorkSize%KUB, Me%WorkSize%KLB, -1
                    Me%Depth%Value3DOut(i, j, k) = Me%Depth%Value3DOut(i, j, k) - Me%Depth%Offset
                enddo
                enddo
                enddo
            endif
                
if32:       if (Me%Depth%GeoVert == sigma_) then            
                call DeAllocateValueIn(Me%Depth%WLValueIn)
            endif if32            
            
i4:         if (Me%Bathym%FromMapping) then
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do k= Me%WorkSize%KLB, Me%WorkSize%KUB
                    if (Me%Mapping%Value3DOut(i,j,Me%WorkSize%KUB)==0) then
                        Me%Bathym%Value2DOut(i,j) = -99.
                        exit                    
                    endif
                    if (Me%Mapping%Value3DOut(i,j,k) == 1) then
                        Me%Bathym%Value2DOut(i,j) = Me%Depth%Value3DOut(i,j,k-1)
                        exit
                    endif
                enddo
                enddo
                enddo
                
                call WriteBathymASCII
                
                if (Me%OutHDF5  ) call WriteBathymHDF5        
                if (Me%OutNetCDF) call WriteBathymNetCDF                    
                
                Me%Bathym%FromMapping = .false. 
                
            endif i4
            
            if (Me%Depth%Interpolate) deallocate(DepthAux)
            
i5:         if (Me%OutHDF5) then

                if (Me%WindowOut%ON) then
                    WorkILB = 1
                    WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
                    WorkJLB = 1
                    WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
                    
                    
                else
                    WorkILB = Me%WorkSize%ILB
                    WorkIUB = Me%WorkSize%IUB        
                    WorkJLB = Me%WorkSize%JLB
                    WorkJUB = Me%WorkSize%JUB                            
                endif
             

                if (Me%WindowOut%ON) then

                    allocate(Vert3D(WorkILB:WorkIUB, WorkJLB:WorkJUB, Me%Size%KLB:Me%Size%KUB))

                    do j = WorkJLB, WorkJUB
                    do i = WorkILB, WorkIUB                
                        Vert3D (i,j,:) = Me%Depth%Value3DOut (i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1,:)
                    enddo
                    enddo
                else            
                    Vert3D  => Me%Depth%Value3DOut
                endif
            
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB,WorkIUB,WorkJLB,WorkJUB,       &
                                                 Me%WorkSize%KLB-1, Me%WorkSize%KUB,    &
                                                 STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR70' 

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",         &
                                     "m", Array3D = Vert3D, OutputNumber = iFinal, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR80' 

                if (Me%WindowOut%ON) deallocate(Vert3D)
                
            endif i5

        enddo d1   

        deallocate(Me%Depth%Value3DOut)

        

   end subroutine WriteDepth        
                  
   !------------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine WriteFieldAllInst(iFinal, iP, iT)
        !Arguments-------------------------------------------------------------
        integer                                         :: iFinal, iP, iT
        
        !Local-----------------------------------------------------------------
        real(8), dimension(:), pointer                  :: DepthAux, ValueAux
        real(8)                                         :: Depthx
        integer                                         :: i, j, k, mask, l, kin, ic
        !Begin-----------------------------------------------------------------
        
        
        if      (Me%Field(iP)%Dim==3) then 
            allocate(Me%Field(iP)%Value3DOut(Me%Size%ILB:Me%Size%IUB,                   &
                                                Me%Size%JLB:Me%Size%JUB,                &
                                                Me%Size%KLB:Me%Size%KUB))
                                                 
        else if (Me%Field(iP)%Dim==2) then 
            allocate(Me%Field(iP)%Value2DOut(Me%Size%ILB:Me%Size%IUB,                   &
                                                Me%Size%JLB:Me%Size%JUB))

        endif             

                            
            
        if       (Me%Field(iP)%Dim==3) then
            
            if (Me%Depth%Interpolate) allocate(ValueAux(1:Me%Depth%kmax),DepthAux(1:Me%Depth%kmax))
            
            do k= Me%WorkSize%KLB, Me%WorkSize%KUB                    
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%Mapping%Value3DOut(i,j,k) == 1) then

                    if (.not. Me%Depth%Interpolate) then
                
                        if (Me%Depth%InvertLayers) then
                            kin = Me%WorkSize%KUB - k + Me%WorkSize%KLB
                        else
                            kin = k
                        endif    
                
                        Me%Field(iP)%Value3DOut(i, j, k) = GetNetCDFValue(Me%Field(iP)%ValueIn,  Dim1 = j+1, &
                                                                            Dim2 = i+1, Dim3 = kin, Dim4 = 1)
                    else
                        ic  = Me%Depth%kmax + 1
                        do l= Me%Depth%kmax,1,-1
                            
                            if (Me%Depth%InvertLayers) then
                                kin = Me%Depth%kmax - l + 1
                            else
                                kin = l
                            endif    
                                
                            Depthx = GetCellInDepth(i, j, l,Me%Depth%kmax, iT)
                                
                            if (Depthx > -100) then
                                    
                                ic = l
                                    
                                ValueAux(l)= GetNetCDFValue(Me%Field(iP)%ValueIn,  Dim1 = j+1, &
                                                            Dim2 = i+1, Dim3 = kin, Dim4 = 1)
                                                                
                                DepthAux(l)= Depthx
                            else
                                exit
                            endif
                        enddo
                            
                        if (ic > Me%Depth%kmax) then
                            stop 'WriteFieldAllInst - ModuleNetCDFCF_2_HDF5MOHID - ERR10' 
                        endif

                        Me%Field(iP)%Value3DOut(i, j, k) = InterpolateProfileR8 (Me%Depth%ZLevels(k), &
                                                            Me%Depth%kmax-ic+1, DepthAux(ic:Me%Depth%kmax), &
                                                            ValueAux(ic:Me%Depth%kmax))
                    endif
                    
                    Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%Value3DOut(i, j, k) * Me%Field(iP)%Multiply 
                    
                    Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%Value3DOut(i, j, k) + Me%Field(iP)%Add
                    
                    Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%Value3DOut(i, j, k) * Me%Field(iP)%UnitsScale
                    
                    Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%Value3DOut(i, j, k) + Me%Field(iP)%UnitsAdd 
                        
                    if (Me%Field(iP)%Value3DOut(i, j, k) < Me%Field(iP)%MinValue) then
                        Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%MinValue
                    endif                                     

                    if (Me%Field(ip)%CheckMinMaxLimits) then
                
                        if (Me%Mapping%Value3DOut(i,j,k) == 1) then
                        
                            if (Me%Field(iP)%Value3DOut(i, j, k) > Me%Field(iP)%MaxLimit) then
                            
                                write(*,*) 'Property=', trim(Me%Field(iP)%NetCDFName)
                                write(*,*) 'In cell i,j,k=',i,j,k
                                write(*,*) 'Netcdf value', Me%Field(iP)%Value3DOut(i, j, k)
                                write(*,*) 'Above Max limit=', Me%Field(iP)%MaxLimit
                            
                                stop 'WriteFieldAllInst - ModuleNetCDFCF_2_HDF5MOHID - ERR20' 
                            
                            endif
                        
                            if (Me%Field(iP)%Value3DOut(i, j, k) < Me%Field(iP)%MinLimit) then
                            
                                write(*,*) 'Property=', trim(Me%Field(iP)%NetCDFName)
                                write(*,*) 'In cell i,j,k=',i,j,k
                                write(*,*) 'Netcdf value', Me%Field(iP)%Value3DOut(i, j, k)
                                write(*,*) 'Below min limit=', Me%Field(iP)%MinLimit
                            
                                stop 'WriteFieldAllInst - ModuleNetCDFCF_2_HDF5MOHID - ERR30' 
                            
                            endif                        
                        endif

                    endif                                    
                    
                else 
                    Me%Field(iP)%Value3DOut(i, j, k) = FillValueReal
                endif                
                
            enddo
            enddo
            enddo
                
            if (Me%Field(iP)%OldMissingValue /= Me%Field(iP)%NewMissingValue) then     
                
                do k= Me%WorkSize%KLB, Me%WorkSize%KUB                    
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB        
                    
                    if (Me%Field(iP)%Value3DOut(i, j, k) == Me%Field(iP)%OldMissingValue) then
                        Me%Field(iP)%Value3DOut(i, j, k) =  Me%Field(iP)%NewMissingValue
                    endif
                    
                enddo
                enddo
                enddo
                
            endif
                
            if (Me%Depth%Interpolate) deallocate(DepthAux, ValueAux)
                
        else if  (Me%Field(iP)%Dim==2) then
            
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                
                if (Me%Depth%Dim3D) then
                    mask = Me%Mapping%Value3DOut(i,j,Me%WorkSize%KUB)
                else
                    mask = Me%Mapping%Value2DOut(i,j)
                endif
                    
                if (mask == 1) then
                    if (Me%Field(iP)%ValueIn%Dim == 3) then
                        Me%Field(iP)%Value2DOut(i, j) = GetNetCDFValue(Me%Field(iP)%ValueIn,  &
                                                            Dim1 = j+1, Dim2 = i+1, Dim3 = 1)
                    else
                        Me%Field(iP)%Value2DOut(i, j) = GetNetCDFValue(Me%Field(iP)%ValueIn,  &
                                                            Dim1 = j+1, Dim2 = i+1, Dim3 = 1, Dim4 = 1)
                    endif
                    Me%Field(iP)%Value2DOut(i, j) = Me%Field(iP)%Value2DOut(i, j) * Me%Field(iP)%Multiply + Me%Field(iP)%Add
                        
                    Me%Field(iP)%Value2DOut(i, j) = Me%Field(iP)%Value2DOut(i, j) * Me%Field(iP)%UnitsScale + Me%Field(iP)%UnitsAdd
                        
                    if (Me%Field(iP)%Value2DOut(i, j) < Me%Field(iP)%MinValue) then
                        Me%Field(iP)%Value2DOut(i, j) = Me%Field(iP)%MinValue
                    endif                              
                        
                    if (Me%Field(ip)%CheckMinMaxLimits) then                      
                    
                        if (Me%Field(iP)%Value2DOut(i, j) > Me%Field(iP)%MaxLimit) then
                            
                            write(*,*) 'Property=', trim(Me%Field(iP)%NetCDFName)
                            write(*,*) 'In cell i,j =',i,j,k
                            write(*,*) 'Netcdf value', Me%Field(iP)%Value2DOut(i, j)
                            write(*,*) 'Above Max limit=', Me%Field(iP)%MaxLimit
                            
                            stop 'WriteFieldAllInst - ModuleNetCDFCF_2_HDF5MOHID - ERR20' 
                            
                        endif
                        
                        if (Me%Field(iP)%Value2DOut(i, j) < Me%Field(iP)%MinLimit) then
                            
                            write(*,*) 'Property=', trim(Me%Field(iP)%NetCDFName)
                            write(*,*) 'In cell i,j=',i,j
                            write(*,*) 'Netcdf value', Me%Field(iP)%Value2DOut(i, j)
                            write(*,*) 'Below min limit=', Me%Field(iP)%MinLimit
                            
                            stop 'WriteFieldAllInst - ModuleNetCDFCF_2_HDF5MOHID - ERR30' 
                            
                        endif                         
                    
                    endif
                                                         
                else 
                    Me%Field(iP)%Value2DOut(i, j) = FillValueReal
                endif                
                    
            enddo
            enddo
                
           if (Me%Field(iP)%OldMissingValue /= Me%Field(iP)%NewMissingValue) then     
                
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB        
                    
                    if (Me%Field(iP)%Value2DOut(i, j) == Me%Field(iP)%OldMissingValue) then
                        Me%Field(iP)%Value2DOut(i, j)  = Me%Field(iP)%NewMissingValue
                    endif
                    
                enddo
                enddo
                
            endif
            
                
        endif

            
        if (Me%Field(iP)%CenterX) then
            
            call CenterProp(iP, CenterX = .true.)
            
        endif
            
        if (Me%Field(iP)%CenterY) then
            
            call CenterProp(iP, CenterY = .true.)            
            
        endif
            
        if (Me%OutHDF5) then
            call WriteFieldHDF5  (iP, iFinal)    
        endif

        if (Me%OutNetCDF) then
            call WriteFieldNetCDF(iP, iFinal)    
        endif
 
        if      (Me%Field(iP)%Dim==3) then
            deallocate(Me%Field(iP)%Value3DOut)
        else if (Me%Field(iP)%Dim==2) then 
            deallocate(Me%Field(iP)%Value2DOut)
        endif    
    
    end subroutine WriteFieldAllInst
!------------------------------------------------------------------------


    subroutine CenterProp(iP, CenterX, CenterY)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: iP
        logical, optional                               :: CenterX, CenterY
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D
        integer                                         :: i, j, k, mask, di, dj
        logical                                         :: CenterX_, CenterY_        
        !Begin-----------------------------------------------------------------
        
        allocate(Aux2D(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        
        if (present(CenterX)) then
            CenterX_ = CenterX
        else
            CenterX_ = .false. 
        endif
        
        if (present(CenterY)) then
            CenterY_ = CenterY
        else
            CenterY_ = .false. 
        endif
        
        if (CenterX_ .and. CenterY_) then
            stop 'CenterProp - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
    
        dj = 0
        di = 0
    
        if (CenterX_) dj = 1
        if (CenterY_) di = 1        

        if       (Me%Field(iP)%Dim==3) then

            do k= Me%WorkSize%KLB, Me%WorkSize%KUB      
            
                Aux2D(:,:) = Me%Field(iP)%Value3DOut(:, :, k)
                
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                
                    Me%Field(iP)%Value3DOut(i, j, k) = 0.
                
                    if (Me%Mapping%Value3DOut(i,j,k) == 1) then
                        if (abs(Aux2D(i,j))<abs(Me%Field(iP)%Limit)) then
                            Me%Field(iP)%Value3DOut(i, j, k) = Aux2D(i,j)/2.
                        endif                                

                        if (abs(Aux2D(i+di,j+dj))<abs(Me%Field(iP)%Limit)) then
                            Me%Field(iP)%Value3DOut(i, j, k) = Me%Field(iP)%Value3DOut(i, j, k) + Aux2D(i+di,j+dj)/2.
                        endif              

                    else 
                        Me%Field(iP)%Value3DOut(i, j, k) = FillValueReal
                    endif                
                enddo
                enddo
            enddo
            
            
        else if  (Me%Field(iP)%Dim==2) then
        
            Aux2D(:,:) = Me%Field(iP)%Value3DOut(:, :, k)

            
            do j= Me%WorkSize%JLB, Me%WorkSize%JUB
            do i= Me%WorkSize%ILB, Me%WorkSize%IUB
            
                Me%Field(iP)%Value2DOut(i, j) = 0.

                if (Me%Depth%Dim3D) then
                    mask = Me%Mapping%Value3DOut(i,j,Me%WorkSize%KUB)
                else
                    mask = Me%Mapping%Value2DOut(i,j)
                endif
                
                if (mask == 1) then

                    if (abs(Aux2D(i,j))<abs(Me%Field(iP)%Limit)) then
                        Me%Field(iP)%Value2DOut(i, j) = Aux2D(i,j)/2.
                    endif                                

                    if (abs(Aux2D(i+di,j+dj))<abs(Me%Field(iP)%Limit)) then
                        Me%Field(iP)%Value2DOut(i, j) = Me%Field(iP)%Value2DOut(i, j) + Aux2D(i+di,j+dj)/2.
                    endif              

                else 
                    Me%Field(iP)%Value2DOut(i, j) = FillValueReal
                endif 
                               
            enddo
            enddo
            
        endif
        
        deallocate(Aux2D)

    end subroutine CenterProp    

    !Computes the depth of a cell in the original grid 
    real function GetCellInDepth (i, j, l,lmax, iT, CellFace, TopDepth)

        !Arguments-------------------------------------------------------------
        integer                         :: i, j, l, lmax, iT
        logical, optional               :: CellFace
        real,    optional               :: TopDepth

        !Local-----------------------------------------------------------------
        integer                         :: lin, klb, kub, klr, kur, k
        real                            :: Aux, A, B, C
        real, dimension(:), allocatable :: Thickness
        logical                         :: CellFace_
        real                            :: TopDepth_
        !Begin-----------------------------------------------------------------

        if (Me%Depth%InvertLayers) then
            lin   = lmax - l  + 1
        else
            lin   = l
        endif              
        
        TopDepth_ = 0.
        CellFace_ = .false.        
        
        if (present(CellFace)) then
            CellFace_ = CellFace
            if (present(TopDepth)) then
                TopDepth_ = TopDepth
            endif                
        endif
        
i11:    if (CellFace_) then

            if (Me%Depth%InvertLayers) then
                lin = lin + 1
            endif

            Aux = GetNetCDFValue(Me%Depth%FaceValueIn,  Dim1 = lin)
        
        else i11
                    
i12:        if (Me%Depth%GeoVert == hybrid) then

i13:            if (Me%Depth%InvertLayers) then
                    klb = 1
                    kub = lin-1
                    klr = 1
                    kur = lin
                else i13
                    klb = lin+1
                    kub = lmax
                    klr = lin
                    kur = lmax
                endif i13          

                allocate(Thickness(klr:kur))
                do k=klr, kur
                    Thickness(k) = GetNetCDFValue(Me%Depth%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = k, Dim4 = iT)
                    if (Thickness(k) == 0) then
                        Thickness(k) = 1e-2
                    endif
                enddo
                Aux = 0
                do k=klb, kub
                    if (Thickness(k)>0) then
                        Aux = Aux + Thickness(k)
                    endif
                enddo
                if (Thickness(lin)>0) then
                    Aux = Aux+Thickness(lin)/2.
                else
                    !stop 'GetCellInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR05' 
                    Aux = FillValueReal
                endif
                deallocate(Thickness)
            else i12
                Aux = GetNetCDFValue(Me%Depth%ValueIn,  Dim1 = lin)
            endif i12
            
        endif i11            

i1:     if (Me%Depth%GeoVert == sigma_) then
            
        
i2:         if      (Me%Depth%Positive == "up"  ) then

                !Do not do anything            
                
            elseif  (Me%Depth%Positive == "down") then i2

                Aux =  - Aux 
            
            else  i2
            
                stop 'GetCellInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR10' 
                
            endif i2
            
            if ( Me%Bathym%Value2DOut(i, j)>-50) then
            
    i3:         if (Me%Depth%RomsDistortion) then
                !ROMS Stretching Function
                    !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit ocean circulation model 
                    !using a generalized topography-following coordinate system, J. Comp. Phys., 115 (1), 228-244. (PDF)
                    !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
                    A = (1.-Me%Depth%theta_b)*sinh(Me%Depth%theta_s * Aux) / sinh(Me%Depth%theta_s)
                    
                    B = tanh(Me%Depth%theta_s*(Aux+0.5))/tanh(0.5*Me%Depth%theta_s)
                    
                    C = A  + Me%Depth%theta_b * 0.5*(B-1.)
                    
                    GetCellInDepth = - (Me%Depth%Hc * Aux + ( Me%Bathym%Value2DOut(i, j) - Me%Depth%Hc) * C)+ Topdepth_ 
                
                else i3
                    
                    GetCellInDepth = - Aux * Me%Bathym%Value2DOut(i, j)
                    
                endif i3
                                                       
            else
                    GetCellInDepth = -99.
            endif                                                       
        !z level or hybrid            
        else if (Me%Depth%GeoVert == z_level .or. Me%Depth%GeoVert == hybrid) then i1

i4:         if      (Me%Depth%Positive == "up"  ) then
            
                GetCellInDepth = Me%Bathym%Value2DOut(i, j) - Aux
                
            elseif  (Me%Depth%Positive == "down") then i4
            
                GetCellInDepth = Aux
                
            elseif  (Me%Depth%Positive == "inverse") then i4                
            
                GetCellInDepth = - Aux            
            
            else i4
                stop 'GetCellInDepth - ModuleNetCDFCF_2_HDF5MOHID - ERR20' 
            endif i4       

        endif i1
            
            

    end function GetCellInDepth


    !--------------------------------------------------------------------------
            
   
    subroutine WriteTimeHDF5
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real(8)                                         :: Aux, daux        
        real,    dimension(:), pointer                  :: TimePtr
        type(T_Time)                                    :: CurrentTime
        integer                                         :: STAT_CALL, i, io

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Time HDF5 file...'

        !allocate(Me%Date%Value1DOut(1:Me%Date%NumberInst))
        

        do i=1, size(Me%Date%ValueInTotal)
        
            io = Me%Date%InstOut(i)           
            
            if (io <0) cycle

            Aux = Me%Date%ValueInTotal(i)

            CurrentTime = Me%Date%RefDateTimeOut
            
            do while (Aux > 0) 
         
                if (Aux > 1e8) then
                    daux = 1e8
                else
                    daux = Aux
                endif                                        
                    
                CurrentTime = CurrentTime + daux
                
                Aux = Aux - daux
            
            enddo
       

!           Dados para escriver uma soa vez cada date:
            call ExtractDate   (CurrentTime,                                            &
                                AuxTime(1), AuxTime(2), AuxTime(3),                     &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                                   &
                                 "Time", "YYYY/MM/DD HH:MM:SS",                         &
                                 Array1D = TimePtr,                                     &
                                 OutputNumber = io, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
        enddo
        
       

    end subroutine WriteTimeHDF5


   !---------------------------------------------------------------------------
   !------------------------------------------------------------------------
    subroutine WriteTimeNetCDF(DefDimTime)
        !Arguments-------------------------------------------------------------
        logical, optional                               :: DefDimTime
        !Local-----------------------------------------------------------------
        real(8), dimension(:), pointer                  :: Times        
        integer                                         :: STAT_CALL, i, TotalInst, TotalInstOut, io
        character(len=Stringlength)                     :: AuxChar

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Time NetCDF file...'
        
        TotalInst    = Me%Date%TotalInst
        TotalInstOut = Me%Date%TotalInstOut

        AuxChar = TimeToStringV2(Me%Date%RefDateTimeOut)
     
        
        allocate(Times(TotalInstOut))
        
        io = 1
        do i=1, TotalInst

            if ( Me%Date%InstOut(i) >0) then
                Times(io) = Me%Date%ValueInTotal(i)
                io = io + 1
            endif
            
        enddo
        
        if (present(DefDimTime)) then
            call NETCDFWriteTime(NCDFID       = Me%NetCDF_Out%ObjNetCDF,                    &
                                 InitialDate  = AuxChar,                                    &
                                 nInstants    = io,                                         &
                                 Times        = Times,                                      &
                                 DefDimTime   = DefDimTime,                                 &
                                 STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        else
        
            call NETCDFWriteTime(NCDFID       = Me%NetCDF_Out%ObjNetCDF,                    &
                                 InitialDate  = AuxChar,                                    &
                                 nInstants    = io,                                         &
                                 Times        = Times,                                      &
                                 STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        endif
        
        deallocate(Times)

        
    end subroutine WriteTimeNetCDF


   !---------------------------------------------------------------------------
   
   !------------------------------------------------------------------------
    subroutine WriteGrid3DHDF5
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer, dimension(:,:,:), pointer              :: Mask3D
        integer, dimension(:,:  ), pointer              :: Mask2D        
        integer                                         :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                         :: STAT_CALL, i, j

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Grid 3D HDF5 file...'
        

        if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
            
            
        else
            WorkILB = Me%WorkSize%ILB
            WorkIUB = Me%WorkSize%IUB        
            WorkJLB = Me%WorkSize%JLB
            WorkJUB = Me%WorkSize%JUB                            
        endif
     
        if (Me%Depth%Dim3D) then

            if (Me%WindowOut%ON) then

                allocate(Mask3D(WorkILB:WorkIUB, WorkJLB:WorkJUB, Me%Size%KLB:Me%Size%KUB))

                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Mask3D (i,j,:) = Me%Mapping%Value3DOut (i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1,:)
                enddo
                enddo
            else            
                Mask3D  => Me%Mapping%Value3DOut
            endif        

        else

            if (Me%WindowOut%ON) then

                allocate(Mask2D(WorkILB:WorkIUB, WorkJLB:WorkJUB))

                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Mask2D (i,j) = Me%Mapping%Value2DOut (i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
                enddo
                enddo
            else            
                Mask2D  => Me%Mapping%Value2DOut
            endif    
        
        endif

        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        if (Me%Depth%Dim3D) then

            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB,        &
                                             Me%WorkSize%KLB, Me%WorkSize%KUB,          &
                                             STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


            call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                   &
                                 "WaterPoints3D", "-",                                  &
                                 Array3D = Mask3D,                                      &
                                 STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            if (Me%WindowOut%ON) deallocate(Mask3D)
            nullify(Mask3D)            
       
        else
            
            if (Me%Mapping%From2D_To_3D) then
                
                allocate(Mask3D(WorkILB:WorkIUB, WorkJLB:WorkJUB, 0:2))
                
                call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB,        &
                                     1, 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR40'   

                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB
                    Mask3D(i,j,1) = Mask2D(i, j)
                enddo
                enddo
                
                call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                   &
                                     "WaterPoints3D", "-",                                  &
                                     Array3D = Mask3D,                                      &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                                
                
            else                

                call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                   &
                                     "WaterPoints2D", "-",                                  &
                                     Array2D = Mask2D,                                      &
                                     STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
                
            endif                
        
            if (Me%WindowOut%ON) deallocate(Mask2D)
            nullify(Mask2D)
           
        endif         
        

    end subroutine WriteGrid3DHDF5


   !---------------------------------------------------------------------------


    subroutine WriteBathymASCII
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(len=StringLength)                     :: Coment1, Coment2 
        real,       dimension(:,:), pointer             :: Bathym
        integer                                         :: STAT_CALL, i, j
        integer                                         :: WorkILB, WorkIUB, WorkJLB, WorkJUB

        !Begin-----------------------------------------------------------------

        write(*,*)"Writing bathymetry..."

        Coment1 = 'File generated by'
        Coment2 = 'Hidromod'

        

        if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
        else
            WorkILB = Me%WorkSize%ILB
            WorkIUB = Me%WorkSize%IUB        
            WorkJLB = Me%WorkSize%JLB
            WorkJUB = Me%WorkSize%JUB                            
        endif
     
        if (Me%WindowOut%ON) then

            allocate(Bathym(WorkILB:WorkIUB, WorkJLB:WorkJUB))

            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB                
                Bathym (i,j) = Me%Bathym%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
            enddo
            enddo
        else            
            Bathym  => Me%Bathym%Value2DOut
        endif                


        call WriteGridData(FileName         = trim(Me%Bathym%FileName),         & 
                           COMENT1          = Coment1,                          &
                           COMENT2          = Coment2,                          &
                           HorizontalGridID = Me%ObjHorizontalGrid,             &
                           FillValue        = -99.,                             &
                           Overwrite        = .true.,                           &
                           GridData2D_Real  = Bathym,                           &
                           STAT             = STAT_CALL) 
        
        if (Me%WindowOut%ON) deallocate(Bathym)
        nullify   (Bathym)
        
    end subroutine WriteBathymASCII        
   !------------------------------------------------------------------------
    subroutine WriteBathymHDF5
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real,       dimension(:,:), pointer             :: Bathym
        integer                                         :: STAT_CALL, i, j
        integer                                         :: WorkILB, WorkIUB, WorkJLB, WorkJUB

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing Grid 3D HDF5 file...'
        

        if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
        else
            WorkILB = Me%WorkSize%ILB
            WorkIUB = Me%WorkSize%IUB        
            WorkJLB = Me%WorkSize%JLB
            WorkJUB = Me%WorkSize%JUB                            
        endif
     
        if (Me%WindowOut%ON) then

            allocate(Bathym(WorkILB:WorkIUB, WorkJLB:WorkJUB))

            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB                
                Bathym (i,j) = Me%Bathym%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
            enddo
            enddo
        else            
            Bathym  => Me%Bathym%Value2DOut
        endif                

        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'


        call HDF5WriteData  (Me%ObjHDF5, "/Grid",                                       &
                             "Bathymetry", "m",                                         &
                             Array2D = Bathym,                                          &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'WriteGrid3DHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'


        if (Me%WindowOut%ON) deallocate(Bathym)
        nullify(Bathym)

    end subroutine WriteBathymHDF5


   !---------------------------------------------------------------------------
   
   
    subroutine ReadTimeNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        character(Len=StringLength)             :: ref_date
        real, dimension(6)                      :: AuxTime
        real(8)                                 :: Aux, HundredDays, Aux1
        integer                                 :: n, status, dimid, i, tmax, jmax
        integer                                 :: stat
        logical                                 :: ReadTime
        type (T_Time)                           :: CurrentTime
        real                                    :: AuxOffSet
        
        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Time NetCDF file...'

        write(*,*)'Reading ', trim(Me%Date%NetCDFDimName)

        status=NF90_INQ_DIMID(ncid,trim(Me%Date%NetCDFDimName),dimid)
        if (status /= nf90_noerr) then
            !Try to rea in String format - WRF model format
            call ReadTimeNetCDFString(ncid)            
        else
        
            status=NF90_INQUIRE_DIMENSION(ncid, dimid, len = Me%Date%NumberInst)
            if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
            call AllocateValueIn(Me%Date%ValueIn, Dim1 = Me%Date%NumberInst)

            status = nf90_inq_varid(ncid, trim(Me%Date%NetCDFName), n)
            if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            call GetNetCDFMatrix(ncid, n, Me%Date%ValueIn) 
            
            if (Me%Date%RefAttribute) then

                if (Me%Date%RefDateName ==  trim(null_str)) then
                    status=NF90_GET_ATT(ncid,n,trim(Me%Date%RefAttributeName), ref_date)
                    if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                else
                    ref_date = trim(Me%Date%RefDateName) 
                endif
                
                
                tmax = len_trim(ref_date)

                ReadTime =.false.
                
                Me%Date%UnitsFactor = 3600.

                do i=1,tmax-5
                    if (ref_date(i:i+5)== "second") then
                        ReadTime =.true.
                        Me%Date%UnitsFactor = 1.
                        exit
                    endif
                enddo
                
                do i=1,tmax-2
                    if (ref_date(i:i+2)== "day" .or. ref_date(i:i+2)== "DAY" .or.           &
                        ref_date(i:i+2)== "Day" .or. ref_date(i:i+2)== "daY" .or.           &
                        ref_date(i:i+2)== "DAy") then
                        ReadTime =.true.
                        Me%Date%UnitsFactor = 86400.
                        exit
                    endif
                enddo            

                do i=1,tmax-5
                    if (ref_date(i:i+5)== "minute") then
                        ReadTime =.true.
                        Me%Date%UnitsFactor = 60.
                        exit
                    endif
                enddo            


                do i=1,tmax-4            
                    if (ref_date(i:i+4)== "since") then
                        ref_date = ref_date(i+5:tmax)
                        exit
                    endif
                    
                enddo

                ReadTime = .false.
                do i=1,len_trim(ref_date)

                    if (ref_date(i:i) ==':') then
                        ReadTime = .true.
                    endif

                enddo  

                do i=1,len(ref_date)

    !                if (ref_date(i:i) =='_'.or.ref_date(i:i) ==':'.or. ref_date(i:i) =='-'&
    !                    .or. ref_date(i:i) =='Z'.or. ref_date(i:i) =='T') then
    !                    ref_date(i:i) = ' '
    !                endif
                    if (ichar(ref_date(i:i))>57 .or. ichar(ref_date(i:i))<48) ref_date(i:i)=' '
                    
                    !write(*,*) ichar("1"), ichar("9"), ichar("0")
                enddo  
                
                jmax = len(ref_date)-3
                
                do i=1,jmax

                    if (ref_date(i:i+3) ==' 00 ') then
                        ref_date(i:i+3) = ' 0  '
                    endif

                    if (ref_date(i:i+3) ==' 0.0') then
                        ref_date(i:i+3) = ' 0  '
                    endif          
                    
                    if (ref_date(i:i+2) ==' 01 ') then
                        ref_date(i:i+2) = '  1 '
                    endif                              

                    if (ref_date(i:i+2) ==' 02 ') then
                        ref_date(i:i+2) = '  2 '
                    endif                   
                    
                    if (ref_date(i:i+2) ==' 03 ') then
                        ref_date(i:i+2) = '  3 '
                    endif        
                    
                    if (ref_date(i:i+2) ==' 04 ') then
                        ref_date(i:i+2) = '  4 '
                    endif                   

                    if (ref_date(i:i+2) ==' 05 ') then
                        ref_date(i:i+2) = '  5 '
                    endif                   

                    if (ref_date(i:i+2) ==' 06 ') then
                        ref_date(i:i+2) = '  6 '
                    endif                   

                    if (ref_date(i:i+2) ==' 07 ') then
                        ref_date(i:i+2) = '  7 '
                    endif                   

                    if (ref_date(i:i+2) ==' 08 ') then
                        ref_date(i:i+2) = '  8 '
                    endif                   

                    if (ref_date(i:i+2) ==' 09 ') then
                        ref_date(i:i+2) = '  9 '
                    endif                   

                enddo            
                
                !ref_date(1:19) = trim(adjustl(ref_date))
                
                AuxTime(:) = 0.

                if (ReadTime) then                            
                    read(ref_date,*,iostat=stat) (AuxTime (i), i = 1, 6)
                    if (stat /= SUCCESS_) then
                        read(ref_date,*,iostat=stat) (AuxTime (i), i = 1, 5)
                        AuxTime(6) = 0
                    endif                    
                    
                else
                    read(ref_date,*) (AuxTime (i), i = 1, 3)
                endif

                            
                call SetDate (Me%Date%RefDateTimeIn, Year    = AuxTime(1),                  &
                                                     Month   = AuxTime(2),                  &
                                                     Day     = AuxTime(3),                  &
                                                     Hour    = AuxTime(4),                  &
                                                     Minute  = AuxTime(5),                  &
                                                     Second  = AuxTime(6))

            endif
            
            if (Me%Date%RefDateOffSetFromAtt) then
            
                AuxOffSet   = ReadOffSetAtt (ncid = ncid,                                   &
                                             Prop = Me%Date%RefDateOffSetProp,              &
                                             Att  = Me%Date%RefDateOffSetAtt)
            endif
                    
            
            do i=1, Me%Date%NumberInst

                Aux = GetNetCDFValue(Me%Date%ValueIn,  Dim1 = i)
                
                Aux = Aux * dble(Me%Date%UnitsFactor)
                
                HundredDays = 100*86400 
                Aux1        = Aux
                call JulianDateToGregorianDate(Me%Date%RefDateTimeIn, CurrentTime)
                
                !~3e7 anos
                if (Aux1 > 1e15) then  
                    write(*,*) 'error in the time instant =',i
                    
                    stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
                endif                 
                 

                if (Aux1 > HundredDays) then            
                    
                    do while (Aux1 > HundredDays)
                        
                        CurrentTime = CurrentTime + HundredDays
                        
                        Aux1 = Aux1 - HundredDays                 
                    enddo
                
                endif
                            
                CurrentTime = CurrentTime + Aux1
                
                !CurrentTime = CurrentTime + Me%Date%RefDateOffSet*86400
                !Date off set in seconds
                CurrentTime = CurrentTime + Me%Date%RefDateOffSet
                
                if (Me%Date%RefDateOffSetFromAtt) then
                    CurrentTime = CurrentTime + AuxOffSet 
                endif
                
                if (i==Me%Date%NumberInst) Me%Date%FileEndTime = CurrentTime
                
                Aux = CurrentTime - Me%Date%RefDateTimeOut
                
                call SetNetCDFValue(Me%Date%ValueIn,  Aux, Dim1 = i)
                
            enddo        
        endif    
 
    end subroutine ReadTimeNetCDF

   !---------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   
   
    subroutine ReadTimeNetCDFString(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        character(Len=StringLength)             :: aux_str
        real                                    :: Year, Month, Day, Hour, Minute, Second
        real(8)                                 :: Aux
        integer                                 :: status, dimid, i, numDims
        integer                                 :: StringLength
        type (T_Time)                           :: CurrentTime
        integer, dimension(nf90_max_var_dims)   :: DimidArray

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Time NetCDF file from an string array...'


        status=nf90_inq_varid(ncid,trim(Me%Date%NetCDFDimName),dimid)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDFString - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        status = nf90_inquire_variable(ncid, dimid, ndims = numDims)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDFString - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        status = nf90_inquire_variable(ncid, dimid, dimids = DimidArray(:numDims))
        if (status /= nf90_noerr) stop 'ReadTimeNetCDFString - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        
        status=NF90_INQUIRE_DIMENSION(ncid, DimidArray(2), len = Me%Date%NumberInst)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDFString - ModuleNetCDFCF_2_HDF5MOHID - ERR40'        

        status=NF90_INQUIRE_DIMENSION(ncid, DimidArray(1), len = StringLength)
        if (status /= nf90_noerr) stop 'ReadTimeNetCDFString - ModuleNetCDFCF_2_HDF5MOHID - ERR50'        

        call AllocateValueIn(Me%Date%ValueIn, Dim1 = Me%Date%NumberInst)

        !if (status /= nf90_noerr) stop 'ReadTimeNetCDFString - ModuleNetCDFCF_2_HDF5MOHID - ERR60'        
        
        do i=1, Me%Date%NumberInst
        

            status = NF90_GET_VAR(ncid,dimid, aux_str(1:StringLength), &
                                 start = (/ 1, i /))            

            !The follow format is assumed 2017-10-09_11:00:00
            read (aux_str( 1: 4),*) Year
            read (aux_str( 6: 7),*) Month
            read (aux_str( 9:10),*) Day
            read (aux_str(12:13),*) Hour
            read (aux_str(15:16),*) Minute
            read (aux_str(18:19),*) Second   
            
            call SetDate(CurrentTime, Year, Month, Day, Hour, Minute, Second)     
            
            Aux = CurrentTime - Me%Date%RefDateTimeOut
            
            call SetNetCDFValue(Me%Date%ValueIn,  Aux, Dim1 = i)
            
            if (i==Me%Date%NumberInst) Me%Date%FileEndTime = CurrentTime            
            
        enddo        
       
 
    end subroutine ReadTimeNetCDFString

   !---------------------------------------------------------------------------
   
   
   real function ReadOffSetAtt (ncid, Prop , Att)

        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        character(Len=*)                        :: Prop , Att

        
        !Local-----------------------------------------------------------------
        integer                                 :: OffSet
        integer                                 :: status, n
        
        !Begin-----------------------------------------------------------------
        
        status = nf90_inq_varid(ncid, trim(Prop), n)
        if (status /= nf90_noerr) stop 'ReadOffSetAtt - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        
        status=NF90_GET_ATT(ncid,n,trim(Att), OffSet)
        if (status /= nf90_noerr) stop 'ReadOffSetAtt - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        
        
        ReadOffSetAtt = real(OffSet)* Me%Date%RefDateOffSetAttFactor

    end function ReadOffSetAtt 
    
   
   !---------------------------------------------------------------------------
   
   
    subroutine ReadGrid2DNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        integer                                 :: status, numDims, i, j
        integer                                 :: RhVarIdLat, RhVarIdLong
        integer, dimension(nf90_max_var_dims)   :: rhDimIdsLat, rhDimIdsLong
        real(8), dimension(:), allocatable      :: Long1D, Lat1D
        real(8), dimension(:,:,:), allocatable  :: Aux3D
        real(8)                                 :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Aux1, Aux2, Aux
        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Grid NetCDF file...'


        status=nf90_inq_varid(ncid,trim(Me%LongLat%NetCDFNameLong),RhVarIdLong)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'

        status = nf90_inquire_variable(ncid, RhVarIdLong, ndims = numDims)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

        status = nf90_inquire_variable(ncid, RhVarIdLong, dimids = rhDimIdsLong(:numDims))
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        
        status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(1), len = Me%LongLat%jmax)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'        
        
        status=nf90_inq_varid(ncid,trim(Me%LongLat%NetCDFNameLat),RhVarIdLat)
        if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        
        if      (numDims == 2 .or. numDims == 3) then     
        
            if (Me%ReadInvertLat) then
                write(*,*) 'Can only invert the latitude reading if the Grid in not 2D'
                stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR45' 
            endif

            status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(2), len = Me%LongLat%imax)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50' 
            
            if (Me%ReadInvertXY) then

                status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(2), len = Me%LongLat%jmax)
                if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'        

                status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLong(1), len = Me%LongLat%imax)
                if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR70' 
            
            endif
                    
        else if (numDims == 1) then
            status = nf90_inquire_variable(ncid, RhVarIdLat, dimids = rhDimIdsLat(:numDims))
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR80'      
            
            status=NF90_INQUIRE_DIMENSION(ncid, rhDimIdsLat(1), len = Me%LongLat%imax)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR90' 
            
        else
            stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
        endif
        
        call AllocateValueIn(Me%LongLat%LongIn, Dim1 =  Me%LongLat%jmax, Dim2 = Me%LongLat%imax)
        call AllocateValueIn(Me%LongLat%LatIn,  Dim1 =  Me%LongLat%jmax, Dim2 = Me%LongLat%imax)        
        
        if      (numDims == 2) then    
            call GetNetCDFMatrix(ncid, RhVarIdLong, Me%LongLat%LongIn) 
            call GetNetCDFMatrix(ncid, RhVarIdLat,  Me%LongLat%LatIn) 
        else if (numDims == 3) then
        
            allocate(Aux3D(1: Me%LongLat%jmax,  Me%LongLat%imax,1))

            status = nf90_get_var(ncid, RhVarIdLong, Aux3D,                          &
                        start = (/ 1,       1, 1 /),                             &
                        count = (/ Me%LongLat%jmax, Me%LongLat%imax, 1 /))  
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR110' 
             
            do j=1, Me%LongLat%jmax
            do i=1, Me%LongLat%imax
                call SetNetCDFValue(Me%LongLat%LongIn, Aux3D (j, i, 1), Dim1 = j,   Dim2 = i  )
            enddo
            enddo    
                                

            status = nf90_get_var(ncid, RhVarIdLat,  Aux3D,                          &
                        start = (/ 1,       1, 1 /),                             &
                        count = (/ Me%LongLat%jmax, Me%LongLat%imax, 1 /))
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR120'          
            
            do j=1, Me%LongLat%jmax
            do i=1, Me%LongLat%imax
                call SetNetCDFValue(Me%LongLat%LatIn,  Aux3D (j, i, 1), Dim1 = j,   Dim2 = i  )
            enddo
            enddo    
                               

            deallocate(Aux3D)
            
        else if (numDims == 1) then
            allocate(Long1D(1:Me%LongLat%jmax))
            allocate(Lat1D (1:Me%LongLat%imax))
            
            status = nf90_get_var(ncid, RhVarIdLong, Long1D)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR110'      
            
            status = nf90_get_var(ncid, RhVarIdLat, Lat1D)
            if (status /= nf90_noerr) stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR120'     
            
            
             
            do j=1, Me%LongLat%jmax
            do i=1, Me%LongLat%imax
                call SetNetCDFValue(Me%LongLat%LatIn,  Lat1D (i), Dim1 = j,   Dim2 = i  )
                call SetNetCDFValue(Me%LongLat%LongIn, Long1D(j), Dim1 = j,   Dim2 = i  )
            enddo
            enddo    
                   
        endif
        
        Me%LongLat%CorrectJUp   = 0
        Me%LongLat%CorrectJDown = 0
        
        do j=1, Me%LongLat%jmax
        do i=1, Me%LongLat%imax        
            Aux = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j, Dim2 = i)
            if (Me%LongLat%Starts180W) then
                if (Aux >= 180) then
                    call SetNetCDFValue(Me%LongLat%LongIn, Aux-360., Dim1 = j, Dim2 = i)
                endif
            endif
        enddo
        enddo            
        
        do j=1, Me%LongLat%jmax-1
            Aux1 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j,   Dim2 = 1  )
            Aux2 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j+1, Dim2 = 1  )            
            if (Aux1 > Aux2) then
                write(*,*) 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - WRN40'
                write(*,*) 'Long > 180 assume Long = Long - 360.'
                write(*,*) 'Limits assumed are -180 < Long. < 180. '
                
                Aux = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j, Dim2 = Me%LongLat%imax)

                if (abs(Aux-Aux1)> 1e-5) then 
                    write(*,*) Aux,'/=', Aux1
                    write(*,*) 
                    write(*,*) 'Correction only valid for regular grids'
                    stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                endif               

                Me%LongLat%BreakJ       = j
                Me%LongLat%CorrectJUp   = Me%LongLat%jmax - Me%LongLat%BreakJ
                Me%LongLat%CorrectJDown = Me%LongLat%BreakJ

                exit
            endif
        enddo    
        
        
        !Build HDF5 MOHID Grid
        Me%WorkSize%ILB = Me%LongLat%dij
        Me%WorkSize%IUB = Me%LongLat%imax - 1 - Me%LongLat%dij
        Me%WorkSize%JLB = Me%LongLat%dij
        Me%WorkSize%JUB = Me%LongLat%jmax - 1 - Me%LongLat%dij
        
        !to warn the user before the model crashes
        !cant use a NetCDF with one of the dimension as 2 (or lower) because IUB or JUB would be zero (or lower).
        if ((Me%WorkSize%IUB < 1) .or. (Me%WorkSize%JUB < 1)) then
            write (*,*)
            write (*,*) 'Please use a NETCDF file with more than'
            write (*,*) '2x2 points so that the grid can be correctly extracted'
            stop 'ReadGrid2DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
        endif

        Me%Size%ILB     = Me%WorkSize%ILB - 1
        Me%Size%IUB     = Me%WorkSize%IUB + 1
        Me%Size%JLB     = Me%WorkSize%JLB - 1
        Me%Size%JUB     = Me%WorkSize%JUB + 1
       
        allocate(Me%LongLat%LongOut(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))
        allocate(Me%LongLat%LatOut (Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB))

        
        
        do j=Me%WorkSize%JLB, Me%WorkSize%JUB+1
        do i=Me%WorkSize%ILB, Me%WorkSize%IUB+1
        
            X1 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j,   Dim2 = i  )
            X2 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j+1, Dim2 = i  )
            X3 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j,   Dim2 = i+1)
            X4 = GetNetCDFValue(Me%LongLat%LongIn, Dim1 = j+1, Dim2 = i+1)

            Y1 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j,   Dim2 = i  )
            Y2 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j+1, Dim2 = i  )
            Y3 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j,   Dim2 = i+1)
            Y4 = GetNetCDFValue(Me%LongLat%LatIn,  Dim1 = j+1, Dim2 = i+1)
            
            Me%LongLat%LongOut(i, j) = (X1 + X2 + X3 + X4) / 4.
            if (Me%ReadInvertLat) then
                Me%LongLat%LatOut (Me%WorkSize%IUB+1+Me%WorkSize%ILB-i, j) = (Y1 + Y2 + Y3 + Y4) / 4.
            else
                Me%LongLat%LatOut (i, j) = (Y1 + Y2 + Y3 + Y4) / 4.
            endif                
            
        enddo
        enddo
        
        call DeAllocateValueIn(Me%LongLat%LongIn)
        call DeAllocateValueIn(Me%LongLat%LatIn )
        
        if (numDims == 1) then
            deallocate(Long1D)
            deallocate(Lat1D )
        endif
        
    end subroutine ReadGrid2DNetCDF

   !---------------------------------------------------------------------------
    subroutine ReadGrid3DNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        real(8)                                 :: DepthAux
        integer                                 :: status, i, j, k, mn, numDims, l, kin
        real(8)                                 :: Aux, Aux1
        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Read Grid 3D NetCDF file...'

  

        if (Me%Mapping%ON) then

            status=nf90_inq_varid(ncid,trim(Me%Mapping%NetCDFName),mn)
            if (status /= nf90_noerr) then
                write (*,*) 'unknown netcdf field = ', trim(Me%Mapping%NetCDFName)
                stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            endif

            status = nf90_inquire_variable(ncid, mn, ndims = numDims)
            if (status /= nf90_noerr) stop 'ReadGrid3DNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
                   
            Me%Mapping%ValueIn%Dim = numDims

            
            if (Me%Depth%Dim3D) then 
                if      (numDims == 3) then
                
                    if (Me%Mapping%Dim3D) then
                        call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                                 Dim2 = Me%LongLat%imax,        &
                                                                 Dim3 = Me%Depth%kmax)
                    else
                        call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                                 Dim2 = Me%LongLat%imax,        &
                                                                 Dim3 = Me%Date%NumberInst)
                    
                    endif                                                                 
                else if (numDims == 4) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax,        &
                                                             Dim3 = Me%Depth%kmax,          &
                                                             Dim4 = Me%Date%NumberInst)
                else if (numDims == 2) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax)
                endif
            else
                if      (numDims == 2) then
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax)
                else if (numDims == 3) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax,        &
                                                             Dim3 = Me%Date%NumberInst)
                                                             
                else if (numDims == 4) then      
                    call AllocateValueIn(Me%Mapping%ValueIn, Dim1 = Me%LongLat%jmax,        &
                                                             Dim2 = Me%LongLat%imax,        &
                                                             Dim3 = 1,                      &
                                                             Dim4 = Me%Date%NumberInst)
                endif        
            endif     

            call GetNetCDFMatrix(ncid, mn, Me%Mapping%ValueIn)         
        
        endif
        
        if (Me%Depth%Dim3D) then 

            allocate(Me%Mapping%Value3DOut(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
            
            Me%Mapping%Value3DOut(:,:,:) = 1

            if (Me%Mapping%ON) then
            
            
                do k= Me%WorkSize%KLB, Me%WorkSize%KUB                    
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                                    
i2:                 if (Me%Depth%Interpolate) then

                        if (Me%Mapping%Limit <= 1) then
                            Aux = Me%Mapping%Limit - 1
                        else
                            Aux = Me%Mapping%Limit + 1                    
                        endif

                        do l= Me%Depth%kmax, 1,-1
                        
                            if (Me%Depth%InvertLayers) then
                                kin =  Me%Depth%kmax - l + 1
                            else
                                kin = l
                            endif                         
                        
                            if      (Me%Mapping%ValueIn%Dim == 4) then
                                Aux1 = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = kin, &
                                                                           Dim4 = Me%Mapping%Instant)
                            else if (Me%Mapping%ValueIn%Dim == 3) then
                                if (Me%Mapping%Dim3D) then                            
                                    Aux1 = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = kin)
                                else
                                    Aux1 = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = Me%Mapping%Instant)
                                endif                                    
                            else if (Me%Mapping%ValueIn%Dim == 2) then
                                Aux1 = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1)
                            endif
                            
                            DepthAux= GetCellInDepth(i, j, l, Me%Depth%kmax, 1) 
                            
                            if (DepthAux > Me%Depth%ZLevels(k)) then
                                Aux = Aux1
                                exit
                            endif                                                            
                        enddo


                                    
                    else i2
                    
                        if (Me%Depth%InvertLayers) then
                            kin = Me%WorkSize%KLB + Me%WorkSize%KUB - k
                        else
                            kin = k
                        endif                     
                
                        if      (Me%Mapping%ValueIn%Dim == 4) then
                            Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = kin, Dim4 = Me%Mapping%Instant)
                        else if (Me%Mapping%ValueIn%Dim == 3) then
                            if (Me%Mapping%Dim3D) then                            
                                Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = kin               )
                            else                                
                                Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = Me%Mapping%Instant)
                            endif
                        else if (Me%Mapping%ValueIn%Dim == 2) then
                            Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1, Dim2 = i+1)
                        endif

                    endif i2    
                    
                    
                    if (Me%Mapping%Limit <= 1) then
                        if (Aux > Me%Mapping%Limit) then
                            Me%Mapping%Value3DOut(i,j,k) = 1
                        else
                            Me%Mapping%Value3DOut(i,j,k) = 0
                        endif
                    else
                        if (Aux < Me%Mapping%Limit) then
                            Me%Mapping%Value3DOut(i,j,k) = 1
                        else
                            Me%Mapping%Value3DOut(i,j,k) = 0
                        endif
                    endif
                enddo
                enddo
                enddo
            

                            
            endif
        else
            allocate(Me%Mapping%Value2DOut(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            
            Me%Mapping%Value2DOut(:,:) = 1
            
            if (Me%Mapping%ON) then

                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB

                    if      (Me%Mapping%ValueIn%Dim == 3) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1,   Dim2 = i+1, Dim3 = Me%Mapping%Instant)
                    else if (Me%Mapping%ValueIn%Dim == 2) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1,   Dim2 = i+1)
                    else if (Me%Mapping%ValueIn%Dim == 4) then
                        Aux = GetNetCDFValue(Me%Mapping%ValueIn,  Dim1 = j+1,   Dim2 = i+1, Dim3 = 1, Dim4 = Me%Mapping%Instant)
                    endif

                    if (Me%Mapping%Limit <= 1) then
                        if (Aux > Me%Mapping%Limit) then
                            Me%Mapping%Value2DOut(i,j) = 1
                        else
                            Me%Mapping%Value2DOut(i,j) = 0
                        endif
                    else
                        if (Aux < Me%Mapping%Limit) then
                            Me%Mapping%Value2DOut(i,j) = 1
                        else
                            Me%Mapping%Value2DOut(i,j) = 0
                        endif
                    endif
                    
                enddo
                enddo
                
            endif            
        endif        
        
   
            
        if (Me%Mapping%ON)                                                              &
            call DeAllocateValueIn(Me%Mapping%ValueIn)
        

        
        
    end subroutine ReadGrid3DNetCDF

    !------------------------------------------------------------------------
    subroutine ReadBathymNetCDF(ncid)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid
        
        !Local-----------------------------------------------------------------
        integer, dimension(nf90_max_var_dims)   :: dimIDs
        integer                                 :: status, i, j, bn, dn, k
        type (T_ValueIn)                        :: Aux
        real(8)                                 :: AuxK
        !Begin-----------------------------------------------------------------    

        !Read bathymetry
       !Read number of layers and their depth
        if (Me%Depth%Dim3D) then

            status = nf90_inq_varid(ncid, trim(Me%Depth%NetCDFName), dn)
            if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'            

            status = nf90_inquire_variable(ncid, dn, dimids = dimIDs(:Me%Depth%ValueIn%Dim))
            if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            if (Me%Depth%GeoVert == Hybrid) then

                status=NF90_INQUIRE_DIMENSION(ncid, dimIDs(3), len = Me%Depth%kmax)
                if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
                call AllocateValueIn(Me%Depth%ValueIn, Dim1 = Me%LongLat%jmax,          &
                                                       Dim2 = Me%LongLat%imax,          &
                                                       Dim3 = Me%Depth%kmax,            &
                                                       Dim4 = Me%Date%NumberInst)
           
            else

                status=NF90_INQUIRE_DIMENSION(ncid, dimIDs(1), len = Me%Depth%kmax)
                if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
                
                call AllocateValueIn(Me%Depth%ValueIn, Dim1 = Me%Depth%kmax)
                
            endif    
            
            call GetNetCDFMatrix(ncid, dn, Me%Depth%ValueIn) 
            
            if (Me%Depth%GeoVert == sigma_) then
            
                Me%Depth%FaceValueIn%Dim = 1

                call AllocateValueIn(Me%Depth%FaceValueIn, Dim1 = Me%Depth%kmax+1)
                
                if (Me%Depth%NetCDFNameFaceOff) then
                
                    Aux%Dim       = 1

                    call AllocateValueIn(Aux, Dim1 = Me%Depth%kmax)
                        
                    status = nf90_inq_varid(ncid, trim(Me%Depth%NetCDFName), bn)
                    if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'                     

                    call GetNetCDFMatrix(ncid, bn, Aux)  
                    
                    do k = 1, Me%Depth%kmax
                        AuxK = GetNetCDFValue(Aux, Dim1 = k)                   
                        call SetNetCDFValue(Me%Depth%FaceValueIn,  AuxK, Dim1 = k)                   
                    enddo
                    
                    k = Me%Depth%kmax
                    
                    AuxK = GetNetCDFValue(Aux, Dim1 = k)                   
                    call SetNetCDFValue(Me%Depth%FaceValueIn,  AuxK, Dim1 = k+1)                   
                    
                    call DeallocateValueIn(Aux)
                
                else                
                    status = nf90_inq_varid(ncid, trim(Me%Depth%NetCDFNameFace), bn)
                    if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'                     

                    call GetNetCDFMatrix(ncid, bn, Me%Depth%FaceValueIn) 
                endif                    
                
            endif                
            
        else
        
            Me%Depth%kmax = -99
        
        endif
        
        Me%WorkSize%KLB = 1
        
        if (Me%Depth%Interpolate) then
            Me%WorkSize%KUB = Me%Depth%N_ZLevels
        else
            Me%WorkSize%KUB = Me%Depth%kmax - Me%Depth%RemoveNsurfLayers
        endif

        Me%Size%KLB = Me%WorkSize%KLB - 1
        Me%Size%KUB = Me%WorkSize%KUB + 1     
                
        if (Me%Bathym%ON) then 
            
            status = nf90_inq_varid(ncid, trim(Me%Bathym%NetCDFName), bn)
            if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'      

            status = nf90_inquire_variable(ncid, bn, ndims = Me%Bathym%ValueIn%Dim)
            if (status /= nf90_noerr) stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR70'

            if      (Me%Bathym%ValueIn%Dim == 2) then
                call AllocateValueIn(Me%Bathym%ValueIn, Dim1 = Me%LongLat%jmax, Dim2 = Me%LongLat%imax)
            elseif  (Me%Bathym%ValueIn%Dim == 3) then
                call AllocateValueIn(Me%Bathym%ValueIn, Dim1 = Me%LongLat%jmax, Dim2 = Me%LongLat%imax, Dim3 = Me%Date%NumberInst)
            else
                stop 'ReadBathymNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
            endif                

            
            call GetNetCDFMatrix(ncid, bn, Me%Bathym%ValueIn) 
            
        endif
        
        
        allocate(Me%Bathym%Value2DOut(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        
        if (Me%Bathym%ON) then 
            
            if      (Me%Bathym%ValueIn%Dim == 2) then            
        
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                    Me%Bathym%Value2DOut(i,j) = GetNetCDFValue(Me%Bathym%ValueIn,  Dim1 = j+1, Dim2 = i+1)
                
                    if (Me%Bathym%Value2DOut(i,j) < -50) Me%Bathym%Value2DOut(i,j) = -99.
                enddo
                enddo
                                    
            elseif  (Me%Bathym%ValueIn%Dim == 3) then

                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                    Me%Bathym%Value2DOut(i,j) = GetNetCDFValue(Me%Bathym%ValueIn,  Dim1 = j+1, Dim2 = i+1, Dim3 = 1)
                
                    if (Me%Bathym%Value2DOut(i,j) < -98) Me%Bathym%Value2DOut(i,j) = -99.
                enddo
                enddo
            endif    
            
            if (Me%Bathym%InvertReferential) then
                
                do j= Me%WorkSize%JLB, Me%WorkSize%JUB
                do i= Me%WorkSize%ILB, Me%WorkSize%IUB
                    if (Me%Bathym%Value2DOut(i,j) > -98) then
                        Me%Bathym%Value2DOut(i,j) = - Me%Bathym%Value2DOut(i,j)
                    endif                        
                enddo
                enddo
                                    
            endif
        
        else

            Me%Bathym%Value2DOut(:,:) = Me%Bathym%Default
        
        endif
        
        if (Me%Bathym%ON)                                                               &
            call DeAllocateValueIn(Me%Bathym%ValueIn)
            
    end subroutine ReadBathymNetCDF         
                     
    
   !---------------------------------------------------------------------------
    subroutine ReadFieldNetCDF(ncid, iP, WriteProp, inst)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid, iP
        logical                                 :: WriteProp
        integer                                 :: inst
        real(8)                                 :: Aux
        
        !Local-----------------------------------------------------------------
        integer                                 :: status, pn, numDims, Layer, LayerDim
        !Begin-----------------------------------------------------------------
        pn = 0
        write(*,*) '  Reading= ', trim(Me%Field(iP)%NetCDFName), inst
        
        status = nf90_inq_varid(ncid, trim(Me%Field(iP)%NetCDFName), pn)

        if (status == nf90_noerr) then
            status = nf90_inquire_variable(ncid, pn, ndims = numDims)
            if (status /= nf90_noerr) stop 'ReadFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            Me%Field(iP)%ValueIn%Dim = numDims
        else 
            write(*,*) 'possible error causes:'
            write(*,*) '  Variable not found.=', trim(Me%Field(iP)%NetCDFName)
            write(*,*) '  The specified netCDF ID does not refer to an open netCDF dataset.'
            stop 'ReadFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
        endif
        
        !Read number of layers and their depth
        if      (Me%Field(iP)%Dim == 3) then
            call AllocateValueIn(Me%Field(iP)%ValueIn, Dim1 = Me%LongLat%jmax,          &
                                                       Dim2 = Me%LongLat%imax,          &
                                                       Dim3 = Me%Depth%kmax,            &
                                                       Dim4 = 1)                                                       
!                                                       Dim4 = Me%Date%NumberInst)
        else if (Me%Field(iP)%Dim == 2) then
            if      (Me%Field(iP)%ValueIn%Dim  == 3) then
                call AllocateValueIn(Me%Field(iP)%ValueIn, Dim1 = Me%LongLat%jmax,      &
                                                           Dim2 = Me%LongLat%imax,      &
                                                           Dim3 = 1)
            elseif (Me%Field(iP)%ValueIn%Dim  == 4) then
                call AllocateValueIn(Me%Field(iP)%ValueIn, Dim1 = Me%LongLat%jmax,      &
                                                           Dim2 = Me%LongLat%imax,      &
                                                           Dim3 = 1,                    &
                                                           Dim4 = 1)
            elseif (Me%Field(iP)%ValueIn%Dim  == 2) then
                call AllocateValueIn(Me%Field(iP)%ValueIn, Dim1 = Me%LongLat%jmax,      &
                                                           Dim2 = Me%LongLat%imax)
            else
                stop 'ReadFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
            endif
        endif     


        if (status /= nf90_noerr) then
            WriteProp = .false.
        else
            if (Me%Field(ip)%ExtractLayer) then
                Layer = Me%Field(ip)%LayerNumber
                LayerDim = Me%Field(ip)%LayerDim
                call GetNetCDFMatrix(ncid, pn, Me%Field(iP)%ValueIn, Inst, Layer, LayerDim)         
            else
                call GetNetCDFMatrix(ncid, pn, Me%Field(iP)%ValueIn, Inst)
            endif                
            WriteProp = .true.
        endif
        
            
        !Read add_ofset
        status=NF90_GET_ATT(ncid,pn,"add_offset", Aux)
        if (status == nf90_noerr) then
            Me%Field(iP)%Add = Aux
        endif            
            
        !Read scale_factor
        status=NF90_GET_ATT(ncid,pn,"scale_factor", Aux)
        if (status == nf90_noerr) then
            Me%Field(iP)%Multiply = Aux
        endif            

        
    end subroutine ReadFieldNetCDF

    !------------------------------------------------------------------------
    
   !---------------------------------------------------------------------------
    subroutine ReadWaterLevelNetCDF(ncid, inst)
        !Arguments-------------------------------------------------------------
        integer                                 :: ncid, inst
        
        !Local-----------------------------------------------------------------
        integer                                 :: status, pn, numDims
        !Begin-----------------------------------------------------------------
        pn = 0
        status = nf90_inq_varid(ncid, Me%Depth%NetCDFNameWL, pn)

        if (status == nf90_noerr) then
            status = nf90_inquire_variable(ncid, pn, ndims = numDims)
            if (status /= nf90_noerr) stop 'ReadFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            Me%Depth%WLValueIn%Dim = numDims
        else 
            write(*,*) 'possible error causes:'
            write(*,*) '  Variable not found.=', trim(Me%Depth%NetCDFNameWL)
            write(*,*) '  The specified netCDF ID does not refer to an open netCDF dataset.'
            stop 'ReadFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR100'
        endif
        

        call AllocateValueIn(Me%Depth%WLValueIn, Dim1 = Me%LongLat%jmax,      &
                                                 Dim2 = Me%LongLat%imax,      &
                                                 Dim3 = 1)

        if (status == nf90_noerr) then
            call GetNetCDFMatrix(ncid, pn, Me%Depth%WLValueIn, Inst)         
        endif
        
    end subroutine ReadWaterLevelNetCDF

    !------------------------------------------------------------------------    
    
    subroutine AllocateValueIn(ValueIn, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------        
        type(T_ValueIn)     :: ValueIn
        integer             :: Dim1
        integer, optional   :: Dim2, Dim3, Dim4
        !Local-----------------------------------------------------------------                
        integer             :: Dim1_, Dim2_, Dim3_, Dim4_
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Dim ==2) then
            if (.not.(present(Dim2))) stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        
        if (Dim ==3) then
            if (.not.(present(Dim2).and.present(Dim3))) stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        endif
        
        if (Dim ==4) then
            if (.not.(present(Dim2).and.present(Dim3).and.present(Dim4))) &
                stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        endif
        
        Dim1_ = Dim1
        if (present(Dim2)) Dim2_ = Dim2
        if (present(Dim3)) Dim3_ = Dim3
        if (present(Dim4)) Dim4_ = Dim4
        

        if (Me%ReadInvertXY .and. Dim >=2) then
            Dim2_ = Dim1
            Dim1_ = Dim2
        endif        
        
        if      (Dim==1) then
            allocate(ValueIn%CountDim(1))
            ValueIn%CountDim(1) = Dim1_
        else if (Dim==2) then
            allocate(ValueIn%CountDim(2))
            ValueIn%CountDim(1) = Dim1_
            ValueIn%CountDim(2) = Dim2_
        else if (Dim==3) then
            allocate(ValueIn%CountDim(3))
            ValueIn%CountDim(1) = Dim1_
            ValueIn%CountDim(2) = Dim2_
            ValueIn%CountDim(3) = Dim3_
        else if (Dim==4) then
            allocate(ValueIn%CountDim(4))
            ValueIn%CountDim(1) = Dim1_
            ValueIn%CountDim(2) = Dim2_
            ValueIn%CountDim(3) = Dim3_        
            ValueIn%CountDim(4) = Dim4_
        else
            write(*,*) 'Dim = ',Dim
            stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R81D(1:Dim1_))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R41D(1:Dim1_))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I41D(1:Dim1_))                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R82D(1:Dim1_,1:Dim2_))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R42D(1:Dim1_,1:Dim2_))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I42D(1:Dim1_,1:Dim2_))                        
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                allocate (ValueIn%R83D(1:Dim1_,1:Dim2_,1:Dim3_))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R43D(1:Dim1_,1:Dim2_,1:Dim3_))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I43D(1:Dim1_,1:Dim2_,1:Dim3_))                        
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                allocate (ValueIn%R84D(1:Dim1_,1:Dim2_,1:Dim3_,1:Dim4_))
                
            else if (DataTypeIn == Real4_   ) then
            
                allocate (ValueIn%R44D(1:Dim1_,1:Dim2_,1:Dim3_,1:Dim4_))            
            
            else if (DataTypeIn == Integer4_) then
            
                allocate (ValueIn%I44D(1:Dim1_,1:Dim2_,1:Dim3_,1:Dim4_))                        
            
            endif    

        else
            stop 'AllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
        endif
        
    end subroutine AllocateValueIn

    !------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    
    subroutine BuildAttributes(NameIn, NCDFName, LongName, StandardName, Units,           &
                               ValidMin, ValidMax, Min, Max, MissingValue, Positive,    &
                               Float2D, Float3D, Int2D, Int3D)

        !Arguments-------------------------------------------------------------
        character(len=*), intent(in )                       :: NameIn
        character(len=*), intent(out)                       :: NCDFName
        character(len=*), intent(out)                       :: LongName
        character(len=*), intent(out)                       :: StandardName
        character(len=*), intent(out)                       :: Units
        character(len=*), intent(out), optional             :: Positive
        real,             intent(out)                       :: ValidMin, Min
        real,             intent(out)                       :: ValidMax, Max
        real,             intent(out)                       :: MissingValue
        real,    dimension(:,:  ), pointer, optional        :: Float2D 
        real,    dimension(:,:,:), pointer, optional        :: Float3D 
        integer, dimension(:,:  ), pointer, optional        :: Int2D   
        integer, dimension(:,:,:), pointer, optional        :: Int3D   

        !Local-----------------------------------------------------------------
        character(len=StringLength)                         :: Name
        integer                                             :: i, j, k
        

        !Begin-----------------------------------------------------------------
        
        call CheckAndCorrectVarName(NameIn, Name)

        select case(trim(adjustl(Name)))

            case("Bathymetry")
                NCDFName        = "bathymetry"
                LongName        = "bathymetry"
                StandardName    = "sea_floor_depth_below_geoid"
                Units           = "m"
                ValidMin        = -50.
                ValidMax        = 11000.
                Positive        = "down"
                MissingValue    = -99.0

            case("WaterPoints2D", "WaterPoints3D", "MappingPoints2D", "WaterPoints")
                NCDFName        = "mask"
                LongName        = "mask of potential water points"
                StandardName    = "land_binary_mask"
                Units           = ""
                ValidMin        = 0
                ValidMax        = 1
            
            case("OpenPoints2D", "OpenPoints3D")
                NCDFName        = "mask"
                LongName        = "mask of effective water points at one given instant"
                StandardName    = "mask"
                Units           = ""
                ValidMin        = 0
                ValidMax        = 1

            case("temperature")
                NCDFName        = "temperature"
                LongName        = "temperature"
                StandardName    = "sea_water_temperature"
                Units           = "degC"
                ValidMin        = 0.
                ValidMax        = 50.
                MissingValue    = FillValueReal

            case("salinity")
                NCDFName        = "salinity"
                LongName        = "salinity"
                StandardName    = "sea_water_salinity"
                Units           = "1e-3"
                ValidMin        = 0.
                ValidMax        = 40.
                MissingValue    = FillValueReal

            case("density")
                NCDFName        = "sea_water_density"
                LongName        = "sea water density"
                StandardName    = "sea_water_density"
                Units           = "kg m-3"
                ValidMin        = 900.
                ValidMax        = 1200.
                MissingValue    = FillValueReal

            case("oxygen")
                NCDFName        = "dissolved_oxygen"
                LongName        = "dissolved oxygen"
                StandardName    = "oxygen"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 30.
                MissingValue    = FillValueReal

            case("dissolved_oxygen_percent_saturation")
                NCDFName        = "dissolved_oxygen_percent_saturation"
                LongName        = "dissolved oxygen percent saturation"
                StandardName    = "dissolved_oxygen_percent_saturation"
                Units           = "%"
                ValidMin        = 0.
                ValidMax        = 200.
                MissingValue    = FillValueReal

            case("velocity_U")
                NCDFName        = "u"
                LongName        = "east-west current velocity"
                StandardName    = "eastward_sea_water_velocity"
                Units           = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = FillValueReal

            case("velocity_V")
                NCDFName        = "v"
                LongName        = "north-south current velocity"
                StandardName    = "northward_sea_water_velocity"
                Units           = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = FillValueReal

            case("velocity_W")
                NCDFName        = "w"
                LongName        = "vertical current velocity"
                StandardName    = "upward_sea_water_velocity"
                Units           = "m s-1"
                ValidMin        = -2.
                ValidMax        = 2.
                MissingValue    = FillValueReal

            case("velocity_modulus")
                NCDFName        = "vm"
                LongName        = "sea water speed"
                StandardName    = "sea_water_speed"
                Units           = "m s-1"
                ValidMin        = -5.
                ValidMax        = 5.
                MissingValue    = FillValueReal

            case("water_level")
                NCDFName        = "ssh"
                LongName        = "sea water level"
                StandardName    = "sea_surface_height"
                Units           = "m"
                ValidMin        = -20.
                ValidMax        = 20.
                MissingValue    = FillValueReal

            case("wind_velocity_X")
                NCDFName        = "x_wind"
                LongName        = "east-west wind velocity"
                StandardName    = "east-west_wind_velocity"
                Units           = "m s-1"
                ValidMin        = -100.
                ValidMax        = 100.
                MissingValue    = FillValueReal

            case("wind_velocity_Y")
                NCDFName        = "y_wind"
                LongName        = "north-south wind velocity"
                StandardName    = "north-south_wind_velocity"
                Units           = "m s-1"
                ValidMin        = -100.
                ValidMax        = 100.
                MissingValue    = FillValueReal

            case("air_temperature")
                NCDFName        = "air_temperature"
                LongName        = "air temperature"
                StandardName    = "air_temperature"
                Units           = "degC"
                ValidMin        = -90.
                ValidMax        = 60.
                MissingValue    = FillValueReal

            case("atmospheric_pressure")
                NCDFName        = "air_pressure"
                LongName        = "atmospheric pressure"
                StandardName    = "atmospheric_pressure"
                Units           = "Pa"
                ValidMin        = -90.
                ValidMax        = 60.
                MissingValue    = FillValueReal

            case("short_wave_solar_radiation_extinction")
                NCDFName        = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
                LongName        = "short wave solar radiation light extinction coefficient"
                StandardName    = "volume_absorption_coefficient_of_radiative_flux_in_sea_water"
                Units           = "m-1"
                ValidMin        = 0.
                ValidMax        = 100.
                MissingValue    = FillValueReal

            case("short_wave_solar_radiation")
                NCDFName        = "short_wave_solar_radiation"
                LongName        = "short wave solar radiation"
                StandardName    = "omnidirectional_photosynthetic_spherical_irradiance_in_sea_water"
                Units           = "W m-2"
                ValidMin        = 0.
                ValidMax        = 1400.
                MissingValue    = FillValueReal

            case("phytoplankton")
                NCDFName        = "phytoplankton"
                LongName        = "phytoplankton"
                StandardName    = "phytoplankton"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("zooplankton")
                NCDFName        = "zooplankton"
                LongName        = "zooplankton"
                StandardName    = "zooplankton"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("nitrate")
                NCDFName        = "nitrate"
                LongName        = "nitrate"
                StandardName    = "nitrate"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("ammonia")
                NCDFName        = "ammonia"
                LongName        = "ammonia"
                StandardName    = "ammonia"
                Units           = "mg l-1"
                ValidMin        = 0.
                ValidMax        = 10.
                MissingValue    = FillValueReal

            case("velocity modulus")
                NCDFName        = "velocity_modulus"
                LongName        = "velocity modulus"
                StandardName    = "velocity modulus"
                Units           = "ms/s"
                ValidMin        = 0.
                ValidMax        = 20.
                MissingValue    = FillValueReal

            case("wind modulus")
                NCDFName        = "wind_modulus"
                LongName        = "wind modulus"
                StandardName    = "wind modulus"
                Units           = "ms/s"
                ValidMin        = 0.
                ValidMax        = 200.
                MissingValue    = FillValueReal
            case default

                NCDFName        = trim(adjustl(Name))
                LongName        = trim(adjustl(Name))
                StandardName    = trim(adjustl(Name))
                Units           = "unknown"
                ValidMin        =   null_real
                ValidMax        = - null_real
                MissingValue    = FillValueReal

        end select

        Min = ValidMax
        Max = ValidMin

if1:   if(present(Int2D) .or. present(Int3D))then
           
            Min = 0
            Max = 1

        elseif(present(Float2D))then if1

            do j = 1, Me%WorkSize%JUB
            do i = 1, Me%WorkSize%IUB

                if(Float2D(j,i) .gt. FillValueReal/2. .and. Float2D(j,i) .lt.  Min .and. &
                                                            Float2D(j,i) .ge. ValidMin)then

                    Min = Float2D(j,i)

                end if

                if(Float2D(j,i) > Max .and. Float2D(j,i) .le. ValidMax) then

                    Max = Float2D(j,i) 

                endif

            enddo
            enddo

        elseif(present(Float3D))then if1

            do j = 1, Me%WorkSize%JUB
            do i = 1, Me%WorkSize%IUB
            do k = 1, Me%WorkSize%KUB

                if(Float3D(j,i,k) .gt. FillValueReal/2. .and. Float3D(j,i,k) .lt.  Min .and. &
                                                              Float3D(j,i,k) .ge. ValidMin)then

                    Min = Float3D(j,i,k)

                end if

                if(Float3D(j,i,k) > Max .and. Float3D(j,i,k) .le. ValidMax) then

                    Max = Float3D(j,i,k) 

                endif

            enddo
            enddo
            enddo

        endif if1

    end subroutine BuildAttributes
    
    
    !--------------------------------------------------------------------------
    
    subroutine CheckAndCorrectVarName(obj_name, Name)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: obj_name

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: Name
        integer                                     :: i

        !Begin-----------------------------------------------------------------

        if(scan(obj_name, "_") .ne. 0 .and. scan(obj_name, "0") .ne. 0)then
            Name = obj_name(1:len_trim(obj_name)-6)
        else
            Name = trim(obj_name)
        endif

        do i = 1, len_trim(Name)
            if (Name(i:i) == ' ') then
                Name(i:i) =  '_'
            endif
        enddo

    end subroutine CheckAndCorrectVarName    
!--------------------------------------------------------------------------    
    
    subroutine DeAllocateValueIn(ValueIn)
        !Arguments-------------------------------------------------------------        
        type(T_ValueIn)     :: ValueIn
        !Local-----------------------------------------------------------------                
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        deallocate(ValueIn%CountDim)

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R81D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R41D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I41D)                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R82D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R42D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I42D)                        
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R83D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R43D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I43D)                        
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                deallocate (ValueIn%R84D)
                
            else if (DataTypeIn == Real4_   ) then
            
                deallocate (ValueIn%R44D)            
            
            else if (DataTypeIn == Integer4_) then
            
                deallocate (ValueIn%I44D)                        
            
            endif    

        else
            stop 'DeAllocateValueIn - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        
    end subroutine DeAllocateValueIn

    !------------------------------------------------------------------------
    
    subroutine GetNetCDFMatrix(ncid, n, ValueIn, inst, layer, layer_dim)
        !Arguments-------------------------------------------------------------        
        integer             :: ncid, n
        type(T_ValueIn)     :: ValueIn
        integer, optional   :: inst        
        integer, optional   :: layer
        integer, optional   :: layer_dim
        !Local-----------------------------------------------------------------                
        real(8), dimension(:,:,:,:), pointer    :: AuxR8D4
        real(8), dimension(:,:,:,:), pointer    :: AuxR4D4        
        integer, dimension(nf90_max_var_dims)   :: dimIDs
        integer                                 :: Dim, DataTypeIn, status, xdim, ydim, zdim, tdim
        integer                                 :: ILB, IUB, JLB, JUB
        integer                                 :: Dim1, Dim2, Dim3, Dim4
        integer                                 :: iSt1, iSt2, iSt3, iSt4
        
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (dim >=2) then
            JLB = 1 + ValueIn%djL
            JUB = ValueIn%CountDim(1) - ValueIn%djU            
            ILB = 1 + ValueIn%diL
            IUB = ValueIn%CountDim(2) - ValueIn%diU            
        endif
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif

        status = nf90_inquire_variable(ncid, n, dimids = dimIDs(:Dim))
        if (status /= nf90_noerr) then
            stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
 

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then

                status = NF90_GET_VAR(ncid,n,ValueIn%R81D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%R41D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%I41D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then

                status = NF90_GET_VAR(ncid,n,ValueIn%R82D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
                
            else if (DataTypeIn == Real4_   ) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%R42D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR70'
            
            else if (DataTypeIn == Integer4_) then
            
                status = NF90_GET_VAR(ncid,n,ValueIn%I42D)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR80'
            
            endif    
    
        else if (Dim==3) then

            if (present(inst)) then   
                                
                status = nf90_inquire_dimension(ncid, dimIDs(1), len = xdim)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR90'

                status = nf90_inquire_dimension(ncid, dimIDs(2), len = ydim)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR100'

            endif

            if      (DataTypeIn == Real8_   ) then

                if (present(inst)) then   
                
                    status = nf90_get_var(ncid,n,ValueIn%R83D,                          &
                            start = (/ 1,       1, inst /),                             &
                            count = (/ xdim, ydim, 1 /))                
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR110'
                    
                else
                    status = NF90_GET_VAR(ncid,n,ValueIn%R83D)
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR120'
                endif
                
            else if (DataTypeIn == Real4_   ) then

                if (present(inst)) then   
                                    
                    status = nf90_get_var(ncid,n,ValueIn%R43D,                          &
                            start = (/ 1, 1, inst /),                                   &
                            count = (/ xdim, ydim, 1 /))                
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR130'
                    
                else
                    status = NF90_GET_VAR(ncid,n,ValueIn%R43D)
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR140'
                endif
                
            else if (DataTypeIn == Integer4_) then
            
                if (present(inst)) then   
                                    
                    status = nf90_get_var(ncid,n,ValueIn%I43D,                          &
                            start = (/ 1, 1, inst /),                                   &
                            count = (/ xdim, ydim, 1 /))                
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR150'
                    
                else
                    status = NF90_GET_VAR(ncid,n,ValueIn%I43D)
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR160'
                endif
            
            endif

        else if (Dim==4) then     
        
            if (present(inst)) then  
                status = nf90_inquire_dimension(ncid, dimIDs(1), len = xdim)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR170'

                status = nf90_inquire_dimension(ncid, dimIDs(2), len = ydim)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR180'

                status = nf90_inquire_dimension(ncid, dimIDs(3), len = zdim)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR190'         
                
                status = nf90_inquire_dimension(ncid, dimIDs(4), len = tdim)
                if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR191'

                if (JUB-JLB+1 /= xdim) then
                    stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR192'
                endif

                if (IUB-ILB+1 /= ydim) then
                    stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR194'
                endif
                  
            endif
            
               

            

            if      (DataTypeIn == Real8_   ) then

                if (present(inst)) then   
                    
                    if (present(layer)) then

                        if      (layer_dim == 3) then
                            iSt1 =    1; iSt2 =    1; iSt3 =layer; iSt4 = inst;
                            Dim1 = xdim; Dim2 = ydim; Dim3 =    1; Dim4 =    1;
                            
                        elseif  (layer_dim == 4) then
                            iSt1 =    1; iSt2 =    1; iSt3 =inst;  iSt4 = layer;
                            Dim1 = xdim; Dim2 = ydim; Dim3 =    1; Dim4 =     1;
                        else
                            stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR196'
                        endif                            
                        
                        allocate(AuxR8D4(1:Dim1,1:Dim2,1:Dim3,1:Dim4))
                        
                        status = nf90_get_var(ncid,n,AuxR8D4,                           &
                                start = (/ iSt1, iSt2, iSt3,  iSt4 /),                  &
                                count = (/ Dim1, Dim2, Dim4,  Dim4 /))                
                        if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR200'                
                        
                        ValueIn%R84D(JLB:JUB,ILB:IUB,1:1,1:1) = AuxR8D4(1:Dim1,1:Dim2,1:Dim3,1:Dim4)
                    
                        deallocate(AuxR8D4)                    
                    
                    else
                    
                        allocate(AuxR8D4(1:xdim,1:ydim,1:zdim,1:1))
                        
                        status = nf90_get_var(ncid,n,AuxR8D4,                           &
                                start = (/    1,    1,    1, inst /),                   &
                                count = (/ xdim, ydim, zdim,    1 /))                
                        if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR202'                
                        
                        ValueIn%R84D(JLB:JUB,ILB:IUB,1:zdim,1:1) = AuxR8D4(1:xdim,1:ydim,1:zdim,1:1)
                    
                        deallocate(AuxR8D4)                    

                    endif                        

                else

                    status = NF90_GET_VAR(ncid,n,ValueIn%R84D)
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR210'

                endif
                
            else if (DataTypeIn == Real4_   ) then
            
                if (present(inst)) then   
                
                    if (present(layer)) then

                        allocate(AuxR4D4(1:xdim,1:ydim,1:1,1:1))
                        
                        status = nf90_get_var(ncid,n,AuxR4D4,                           &
                                start = (/    1,    1, layer, inst /),                  &
                                count = (/ xdim, ydim,     1,    1 /))                
                        if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR212'                
                        
                        ValueIn%R44D(JLB:JUB,ILB:IUB,1:1,1:1) = AuxR4D4(1:xdim,1:ydim,1:1,1:1)
                    
                        deallocate(AuxR4D4)    

                    else                                       
                                    
                        allocate(AuxR4D4(1:xdim,1:ydim,1:zdim,1:1))
                        
                        status = nf90_get_var(ncid,n,AuxR4D4,                               &
                                start = (/    1,    1,    1, inst /),                       &
                                count = (/ xdim, ydim, zdim,    1 /))                
                        if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR214'                
                        
                        ValueIn%R44D(JLB:JUB,ILB:IUB,1:zdim,1:1) = AuxR4D4(1:xdim,1:ydim,1:zdim,1:1)
                        
                        deallocate(AuxR4D4)                    
                        
                    endif                        

                else

                    status = NF90_GET_VAR(ncid,n,ValueIn%R44D)
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR216'

                endif
                            
            else if (DataTypeIn == Integer4_) then
            
                if (present(inst)) then   
                                    
                    status = nf90_get_var(ncid,n,ValueIn%I44D,                          &
                            start = (/    1,    1,    1, inst /),                       &
                            count = (/ xdim, ydim, zdim,    1 /))                
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR220'                

                else

                    status = NF90_GET_VAR(ncid,n,ValueIn%I44D)
                    if (status /= nf90_noerr) stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR230'

                endif            
            endif

        else
            stop 'GetNetCDFMatrix - ModuleNetCDFCF_2_HDF5MOHID - ERR240'
        endif
        
    end subroutine GetNetCDFMatrix

    !------------------------------------------------------------------------


    function GetNetCDFValue(ValueIn, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------  
        real(8)                         :: GetNetCDFValue      
        type(T_ValueIn)                 :: ValueIn
        integer, intent(in)             :: Dim1
        integer, intent(in), optional   :: Dim2, Dim3, Dim4
        !Local-----------------------------------------------------------------                
        integer             :: Dim1_, Dim2_, Dim3_, Dim4_
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Dim ==2) then
            if (.not.(present(Dim2))) stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        
        if (Dim ==3) then
            if (.not.(present(Dim2).and.present(Dim3))) stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        endif
        
        if (Dim ==4) then
            if (.not.(present(Dim2).and.present(Dim3).and.present(Dim4))) then
                stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            endif
        endif

        Dim1_ = Dim1
        if (present(Dim2)) Dim2_ = Dim2
        if (present(Dim3)) Dim3_ = Dim3        
        if (present(Dim4)) Dim4_ = Dim4     
        

        if (Dim >= 2) then
            if (Dim1_ >Me%LongLat%BreakJ) then
                Dim1_ = Dim1_ - Me%LongLat%CorrectJDown
            else
                Dim1_ = Dim1_ + Me%LongLat%CorrectJUp
            endif
        endif
                

        if (Me%ReadInvertXY .and. Dim >=2) then
            Dim1_ = Dim2
            Dim2_ = Dim1
        endif      
        
        if (Me%ReadInvertLat .and. Dim >=2) then
            Dim2_ = Me%LongLat%imax - Dim2_ + 1 
        endif  

        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R81D(Dim1_)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R41D(Dim1_)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I41D(Dim1_)                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R82D(Dim1_,Dim2_)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R42D(Dim1_,Dim2_)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I42D(Dim1_,Dim2_)                        
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R83D(Dim1_,Dim2_,Dim3_)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R43D(Dim1_,Dim2_,Dim3_)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I43D(Dim1_,Dim2_,Dim3_)                        
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                GetNetCDFValue = ValueIn%R84D(Dim1_,Dim2_,Dim3_,Dim4_)
                
            else if (DataTypeIn == Real4_   ) then
            
                GetNetCDFValue = ValueIn%R44D(Dim1_,Dim2_,Dim3_,Dim4_)            
            
            else if (DataTypeIn == Integer4_) then
            
                GetNetCDFValue = ValueIn%I44D(Dim1_,Dim2_,Dim3_,Dim4_)                        
            
            endif    

        else
            stop 'GetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        
        if (ISNAN(GetNetCDFValue)) then
            if (Me%NaN_2_Null) then
                GetNetCDFValue = 0.
            else
                GetNetCDFValue = FillValueReal            
            endif
        endif
        
    end function GetNetCDFValue

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------

    subroutine SetNetCDFValue(ValueIn, SetValue, Dim1, Dim2, Dim3, Dim4)
        !Arguments-------------------------------------------------------------  
        type(T_ValueIn)     :: ValueIn
        real(8)             :: SetValue                
        integer             :: Dim1
        integer, optional   :: Dim2, Dim3, Dim4
        !Local-----------------------------------------------------------------                
        integer             :: Dim1_, Dim2_, Dim3_, Dim4_
        integer             :: Dim, DataTypeIn
        !Begin-----------------------------------------------------------------        
        
        Dim         = ValueIn%Dim
        DataTypeIn  = ValueIn%DataType
        
        if (DataTypeIn /= Real8_ .and. DataTypeIn /= Real4_ .and. DataTypeIn /= Integer4_) then
            stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
        endif
        
        if (Dim ==2) then
            if (.not.(present(Dim2))) stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
        endif
        
        if (Dim ==3) then
            if (.not.(present(Dim2).and.present(Dim3))) stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        endif
        
        if (Dim ==4) then
            if (.not.(present(Dim2).and.present(Dim3).and.present(Dim4))) &
                stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        endif
        
        Dim1_ = Dim1
        if (present(Dim2)) Dim2_ = Dim2
        if (present(Dim3)) Dim3_ = Dim3        
        if (present(Dim4)) Dim4_ = Dim4                
        
                
        if (Me%ReadInvertXY .and. Dim >=2) then
            Dim1_ = Dim2
            Dim2_ = Dim1
        endif         
        
        if (Me%ReadInvertLat .and. Dim >=2) then
            Dim2_ = Me%LongLat%imax - Dim2_ + 1 
        endif          
        
        if      (Dim==1) then
        
            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R81D(Dim1_) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
               ValueIn%R41D(Dim1_)  = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I41D(Dim1_) = SetValue                        
            
            endif
        
        else if (Dim==2) then
    
            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R82D(Dim1_,Dim2_) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
                ValueIn%R42D(Dim1_,Dim2_) = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I42D(Dim1_,Dim2_) = SetValue
            
            endif    
    
        else if (Dim==3) then

            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R83D(Dim1_,Dim2_,Dim3_) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
                ValueIn%R43D(Dim1_,Dim2_,Dim3_) = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I43D(Dim1_,Dim2_,Dim3_) = SetValue
            
            endif    
        
        else if (Dim==4) then        

            if      (DataTypeIn == Real8_   ) then
                
                ValueIn%R84D(Dim1_,Dim2_,Dim3_,Dim4_) = SetValue
                
            else if (DataTypeIn == Real4_   ) then
            
                ValueIn%R44D(Dim1_,Dim2_,Dim3_,Dim4_) = SetValue
            
            else if (DataTypeIn == Integer4_) then
            
                ValueIn%I44D(Dim1_,Dim2_,Dim3_,Dim4_) = SetValue
            
            endif    

        else
            stop 'SetNetCDFValue - ModuleNetCDFCF_2_HDF5MOHID - ERR50'
        endif
        
    end subroutine SetNetCDFValue

    !------------------------------------------------------------------------
        
    subroutine WriteFieldHDF5(iP, iFinal)

        !Arguments-------------------------------------------------------------
        integer                                         :: iP, iFinal

        !Local-----------------------------------------------------------------
        real,   pointer, dimension(:,:,:)               :: Aux3D
        real,   pointer, dimension(:,:)                 :: Aux2D, Aux2DV, Aux2DAcc
        integer                                         :: STAT_CALL
        integer                                         :: WorkILB, WorkJLB, WorkKLB
        integer                                         :: WorkIUB, WorkJUB, WorkKUB
        integer                                         :: i, j, k, Mapping
        real                                            :: AuxOld
        character (len=StringLength)                    :: AccumulatedPropName
        type (T_Time)                                   :: CurrentTime 
        real(8)                                         :: Aux
        real                                            :: HoursOut
        logical                                         :: ComputeAccumlatedDif
        !Begin-----------------------------------------------------------------
        
        !Bounds



        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

        if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
        endif
        
        if (Me%Field(ip)%From2D_To_3D) then
            
            allocate(Aux3D(WorkILB:WorkIUB, WorkJLB:WorkJUB, 0:2))
            

            if (Me%WindowOut%ON) then

                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Aux3D(i,j,1) = Me%Field(iP)%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
                enddo
                enddo
            else
                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Aux3D(i,j,1) = Me%Field(iP)%Value2DOut(i,j)
                enddo
                enddo
            endif  
            
            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB  
                if (Me%Mapping%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1) == 0) then
                    Aux3D(i,j,1) = FillValueReal
                endif
            enddo
            enddo
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                 &
                                 WorkJUB, 1, 1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR01'        
            
            if (Me%Field(iP)%Accumulated2StepGFS) then
                stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR55'
            endif                

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array3D = Aux3D,                                       &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR05'
            
            deallocate(Aux3D)
            

        elseif  (Me%Field(iP)%Dim==2) then
        
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'        
            

            if (Me%WindowOut%ON) then

                allocate(Aux2D(WorkILB:WorkIUB, WorkJLB:WorkJUB))
                
                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Aux2D(i,j) = Me%Field(iP)%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
                enddo
                enddo
                
            else
                Aux2D => Me%Field(iP)%Value2DOut
            endif
            
            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB                
                if (associated(Me%Mapping%Value2DOut)) then
                    Mapping = Me%Mapping%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)
                else
                    Mapping = Me%Mapping%Value3DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1, Me%WorkSize%KUB)
                endif
                
                if (Mapping == 0) then 
                    Aux2D(i,j) = FillValueReal
                else
                    !Direction property - from meteorological convention to algebric 
                    if (Me%Field(iP)%FromMeteo2Algebric .or. Me%Field(iP)%FromCartesian2Meteo) then
                        Aux2D(i,j) = 270. - Aux2D(i,j)    
                    endif
                endif
            enddo
            enddo
            
            if (Me%Field(iP)%Accumulated2Step) then
            
                Aux         = Me%Date%ValueInTotal(iFinal)            
                CurrentTime = Me%Date%RefDateTimeOut + Aux
           

    !           Dados para escriver uma soa vez cada date:
                call ExtractDate   (CurrentTime, Hour = HoursOut)            
            
            
                AccumulatedPropName = trim(Me%Field(iP)%ID%Name)//" accumulated"
                
                ComputeAccumlatedDif = .true. 
                
                if (iFinal == 1) then
                    ComputeAccumlatedDif = .false. 
                endif
                
                if (Me%Field(iP)%Accumulated2StepGFS) then
                    ! The option Accumulated2StepGFS was developed specifically to face a problem 
                    ! identify in the GFS opendap service and grib output files. Precipitation 
                    ! in mm has a variable accumulation period. For precipitations at 
                    ! 3h, 9h, 15h and 21h the precipitation in mm is relative to an accumulation period 
                    ! of 3h the precipitations at 0h,6h,12h,18h and 24h is relative to a 6 h accumulation 
                    ! period. The approach followed to solve the problem consists in subtract the instants 
                    ! with an accumulated period of 6 h by the precipitation of the previous instant 
                    ! (accumulated period of 3 h). This way the precipitation is always the accumulated 
                    ! value since the last instant.
                
                    !Instant is 0h, 6h, 12h, 18h, 24h 
                    if (mod(HoursOut,6.) /= 0) then
                        ComputeAccumlatedDif = .false.
                    endif
                endif                    


                if (ComputeAccumlatedDif) then
                
                    allocate(Aux2DAcc(WorkILB:WorkIUB, WorkJLB:WorkJUB))
                    
                    call HDF5ReadData  (Me%ObjHDF5, "/Results/"//trim(AccumulatedPropName),&
                                         trim(AccumulatedPropName),                        &
                                         Array2D = Aux2DAcc,                               &
                                         OutputNumber = iFinal-1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                                        

                    do j = WorkJLB, WorkJUB
                    do i = WorkILB, WorkIUB  
                        AuxOld        = Aux2DAcc(i,j)
                        Aux2DAcc(i,j) = Aux2D   (i, j)
                        if (AuxOld > 0.) then                    
                            Aux2D(i,j) = Aux2D(i,j) - AuxOld
                            if (Aux2D(i,j) < 0.) Aux2D(i,j) = 0.
                        endif
                    enddo
                    enddo
                    
                    
                    call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(AccumulatedPropName),&
                                         trim(AccumulatedPropName),                         &
                                         trim(Me%Field(iP)%ID%Units),                       &
                                         Array2D = Aux2DAcc,                                &
                                         OutputNumber = iFinal, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
                    
                    deallocate(Aux2DAcc)                    
                    
                else
                
                    call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(AccumulatedPropName),&
                                         trim(AccumulatedPropName),                         &
                                         trim(Me%Field(iP)%ID%Units),                       &
                                         Array2D = Aux2D,                                   &
                                         OutputNumber = iFinal, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

                endif                    
                
            endif                

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array2D = Aux2D,                                       &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            if (Me%Field(iP)%FromDir2Vector) then
            
                allocate(Aux2DV(WorkILB:WorkIUB, WorkJLB:WorkJUB))            
                
                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    !Meteorological convention 
                    if (Me%Mapping%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1) == 1) then
                        if (Me%Field(iP)%FromMeteo2Algebric) then                    
                            Aux2DV(i,j) = cos((      Aux2D(i,j))*Pi/180.) 
                        else
                            Aux2DV(i,j) = cos((270 - Aux2D(i,j))*Pi/180.) 
                        endif
                    else
                        Aux2DV(i,j) = FillValueReal
                    endif
                enddo
                enddo

                
                call HDF5WriteData  (HDF5ID       = Me%ObjHDF5,                             &
                                     GroupName    = "/Results/"//trim(Me%Field(iP)%DirX),&
                                     Name         = trim(Me%Field(iP)%DirX),             &
                                     Units        = "-",                                    &
                                     Array2D      = Aux2DV,                                  &
                                     OutputNumber = iFinal,                                 &
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    !Meteorological convention 
                    if (Me%Mapping%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1) == 1) then
                        if (Me%Field(iP)%FromMeteo2Algebric) then                    
                            Aux2DV(i,j) = sin((      Aux2D(i,j))*Pi/180.) 
                        else                    
                            Aux2DV(i,j) = sin((270 - Aux2D(i,j))*Pi/180.) 
                        endif
                    else
                        Aux2DV(i,j) = FillValueReal
                    endif
                enddo
                enddo
                
                call HDF5WriteData  (HDF5ID       = Me%ObjHDF5,                         &
                                     GroupName    = "/Results/"//trim(Me%Field(iP)%DirY),&
                                     Name         = trim(Me%Field(iP)%DirY),            &
                                     Units        = "-",                                &
                                     Array2D      = Aux2DV,                             &
                                     OutputNumber = iFinal,                             &
                                     STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR40'

                deallocate(Aux2DV)                                    
            
                if (Me%WindowOut%ON) deallocate(Aux2D)                

            endif
            


        else if (Me%Field(iP)%Dim==3) then

            if (Me%WindowOut%ON) then

                allocate(Aux3D(WorkILB:WorkIUB, WorkJLB:WorkJUB, Me%Size%KLB:Me%Size%KUB))

                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Aux3D(i,j,:) = Me%Field(iP)%Value3DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1,:)
                enddo
                enddo
            else
                Aux3D => Me%Field(iP)%Value3DOut
            endif  
            
            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB  
            do k = WorkKLB, WorkKUB              
                if (Me%Mapping%Value3DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1,k) == 0) then
                    Aux3D(i,j,k) = FillValueReal
                endif
            enddo
            enddo
            enddo
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                 &
                                 WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR50'        
            
            if (Me%Field(iP)%Accumulated2StepGFS) then
                stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR55'
            endif                

            call HDF5WriteData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),   &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 trim(Me%Field(iP)%ID%Units),                           &
                                 Array3D = Aux3D,                                       &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR60'
            
            if (Me%WindowOut%ON) deallocate(Aux3D)
            
        else 

            stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR70'

        endif


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR90'


    end subroutine WriteFieldHDF5

    !----------------------------------------------------------------------
    
    !------------------------------------------------------------------------
        
    subroutine ReadVerticalZHDF5(iFinal)

        !Arguments-------------------------------------------------------------
        integer                                         :: iFinal

        !Local-----------------------------------------------------------------
        real,       dimension(:,:,:), pointer           :: Aux3D
        integer                                         :: STAT_CALL, i, j
        integer                                         :: WorkILB, WorkJLB, WorkKLB
        integer                                         :: WorkIUB, WorkJUB, WorkKUB

        !Begin-----------------------------------------------------------------
        
        !Bounds

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 

       if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
        else
            WorkILB = Me%WorkSize%ILB 
            WorkIUB = Me%WorkSize%IUB 
            WorkJLB = Me%WorkSize%JLB 
            WorkJUB = Me%WorkSize%JUB 
        endif
        
        call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                     &
                             WorkJUB, WorkKLB-1, WorkKUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
        
        if (.not. associated(Me%Depth%Value3DOut)) then
            allocate(Me%Depth%Value3DOut(Me%Size%ILB:Me%Size%IUB,                       &
                                         Me%Size%JLB:Me%Size%JUB,                       &
                                         Me%Size%KLB:Me%Size%KUB))
        endif
        
        if (Me%WindowOut%ON) then
            allocate(Aux3D(WorkILB:WorkIUB, WorkJLB:WorkJUB,Me%Size%KLB:Me%Size%KUB))
        else
            Aux3D => Me%Depth%Value3DOut
        endif            

        call HDF5ReadData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",                  &
                             Array3D = Aux3D, OutputNumber = iFinal, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
        
        if (Me%WindowOut%ON) then            
            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB                
                Me%Depth%Value3DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1,:)= Aux3D(i,j,:) 
            enddo
            enddo
            deallocate(Aux3D)
        endif       
        nullify(Aux3D)     

    end subroutine ReadVerticalZHDF5
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
        
    subroutine ReadFieldHDF5(iP, iFinal)

        !Arguments-------------------------------------------------------------
        integer                                         :: iP, iFinal

        !Local-----------------------------------------------------------------
        real,       dimension(:,:,:), pointer           :: Aux3D
        real,       dimension(:,:  ), pointer           :: Aux2D        
        integer                                         :: STAT_CALL, i, j
        integer                                         :: WorkILB, WorkJLB, WorkKLB
        integer                                         :: WorkIUB, WorkJUB, WorkKUB

        !Begin-----------------------------------------------------------------
        
        !Bounds



        WorkILB = Me%WorkSize%ILB 
        WorkIUB = Me%WorkSize%IUB 

        WorkJLB = Me%WorkSize%JLB 
        WorkJUB = Me%WorkSize%JUB 

        WorkKLB = Me%WorkSize%KLB 
        WorkKUB = Me%WorkSize%KUB 


       if (Me%WindowOut%ON) then
            WorkILB = 1
            WorkIUB = Me%WindowOut%IUB - Me%WindowOut%ILB + 1
            WorkJLB = 1
            WorkJUB = Me%WindowOut%JUB - Me%WindowOut%JLB + 1
        endif

        if      (Me%Field(iP)%Dim==2) then
        
            if (Me%WindowOut%ON) then

                allocate(Aux2D(WorkILB:WorkIUB, WorkJLB:WorkJUB))
            else
                Aux2D => Me%Field(iP)%Value2DOut
            endif

            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB, WorkJUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR10'        

            call HDF5ReadData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),    &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 Array2D = Aux2D,                                       &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR20'

            if (Me%WindowOut%ON) then            
                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Me%Field(iP)%Value2DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1)= Aux2D(i,j) 
                enddo
                enddo
                deallocate(Aux2D)
            endif
            nullify(Aux2D)
            

        else if (Me%Field(iP)%Dim==3) then
        
            call HDF5SetLimits  (Me%ObjHDF5, WorkILB, WorkIUB, WorkJLB,                 &
                                 WorkJUB, WorkKLB, WorkKUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
            
            if (Me%WindowOut%ON) then

                allocate(Aux3D(WorkILB:WorkIUB, WorkJLB:WorkJUB,Me%Size%KLB:Me%Size%KUB))
            else
                Aux3D => Me%Field(iP)%Value3DOut
            endif            

            call HDF5ReadData  (Me%ObjHDF5, "/Results/"//trim(Me%Field(iP)%ID%Name),    &
                                 trim(Me%Field(iP)%ID%Name),                            &
                                 Array3D = Aux3D,                                       &
                                 OutputNumber = iFinal, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            

            if (Me%WindowOut%ON) then            
                do j = WorkJLB, WorkJUB
                do i = WorkILB, WorkIUB                
                    Me%Field(iP)%Value3DOut(i+ Me%WindowOut%ILB - 1,j+ Me%WindowOut%JLB - 1,:)= Aux3D(i,j,:) 
                enddo
                enddo
                deallocate(Aux3D)
            endif       
            nullify(Aux3D)     

        else 

            stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        endif
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldHDF5 - ModuleNetCDFCF_2_HDF5MOHID - ERR90'
        

    end subroutine ReadFieldHDF5
    !------------------------------------------------------------------------
        
    subroutine WriteFieldNetCDF(iP, iFinal)

        !Arguments-------------------------------------------------------------
        integer                                         :: iP, iFinal

        !Local-----------------------------------------------------------------
        real,  dimension(:,:,:), pointer                :: Field3D
        real,  dimension(:,:  ), pointer                :: Field2D        
        integer                                         :: STAT_CALL, i, j, k 
        integer                                         :: WorkIUB, WorkJUB, WorkKUB
        character(len=StringLength)                     :: NCDFName, LongName, StandardName, Units
        real                                            :: MinValue, MaxValue, ValidMin, ValidMax, MissingValue
        !Begin-----------------------------------------------------------------
        
        !Bounds
        WorkIUB = Me%WorkSize%IUB 
        WorkJUB = Me%WorkSize%JUB 
        WorkKUB = Me%WorkSize%KUB 



        if      (Me%Field(iP)%Dim==2) then
        
            allocate(Field2D(1:WorkJUB,1:WorkIUB))

            do j = 1, WorkJUB
            do i = 1, WorkIUB
                Field2D(j,i)=Me%Field(iP)%Value2DOut(i,j)
               !Direction property - from meteorological convention to algebric 
                if (Me%Field(iP)%FromMeteo2Algebric) then
                    Field2D(j,i) = 270. - Field2D(j,i)    
                endif                
            enddo
            enddo      
            
            do j = 1, WorkJUB
            do i = 1, WorkIUB
                Field2D(j,i)=Me%Field(iP)%Value2DOut(i,j)
               !Direction property - from cartesian (or algebric) convention to meteorological
                if (Me%Field(iP)%FromCartesian2Meteo) then
                    Field2D(j,i) = 270. - Field2D(j,i)    
                endif                
            enddo
            enddo                           
             

            call BuildAttributes(trim(Me%Field(iP)%ID%Name), NCDFName, LongName, StandardName, &
                                       Units, ValidMin, ValidMax,        &
                                       MinValue, MaxValue, MissingValue, Float2D = Field2D)


            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,  &
                                  Name          = trim(NCDFName),           &
                                  LongName      = trim(LongName),           &
                                  StandardName  = trim(StandardName),       & 
                                  Units         = trim(Units),              &
                                  ValidMin      = ValidMin,                 &
                                  ValidMax      = ValidMax,                 &
                                  MinValue      = MinValue,                 &
                                  MaxValue      = MaxValue,                 &
                                  OutputNumber  = iFinal,                   &
                                  Array2D       = Field2D,                  &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR10'        
            
            if (Me%Field(iP)%FromDir2Vector) then
            
                do j = 1, WorkJUB
                do i = 1, WorkIUB
                    if (Me%Field(iP)%FromMeteo2Algebric) then
                        Field2D(j,i) = cos((      Me%Field(iP)%Value2DOut(i,j))*Pi/180.) 
                    else                        
                        !Meteorological convention 
                        Field2D(j,i) = cos((270 - Me%Field(iP)%Value2DOut(i,j))*Pi/180.) 
                    endif
                enddo
                enddo

                
                call BuildAttributes(trim(Me%Field(iP)%DirX), NCDFName, LongName, StandardName, &
                                           Units, ValidMin, ValidMax,        &
                                           MinValue, MaxValue, MissingValue, Float2D = Field2D)


                call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,  &
                                      Name          = trim(NCDFName),           &
                                      LongName      = trim(LongName),           &
                                      StandardName  = trim(StandardName),       & 
                                      Units         = trim(Units),              &
                                      ValidMin      = ValidMin,                 &
                                      ValidMax      = ValidMax,                 &
                                      MinValue      = MinValue,                 &
                                      MaxValue      = MaxValue,                 &
                                      OutputNumber  = iFinal,                   &
                                      Array2D       = Field2D,                  &
                                      STAT          = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR20'                  

                do j = 1, WorkJUB
                do i = 1, WorkIUB
                    if (Me%Field(iP)%FromMeteo2Algebric) then
                        Field2D(j,i) = sin((      Me%Field(iP)%Value2DOut(i,j))*Pi/180.) 
                    else
                        !Meteorological convention 
                        Field2D(j,i) = sin((270 - Me%Field(iP)%Value2DOut(i,j))*Pi/180.) 
                    endif
                enddo
                enddo

                call BuildAttributes(trim(Me%Field(iP)%DirY), NCDFName, LongName, StandardName, &
                                           Units, ValidMin, ValidMax,        &
                                           MinValue, MaxValue, MissingValue, Float2D = Field2D)


                call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,  &
                                      Name          = trim(NCDFName),           &
                                      LongName      = trim(LongName),           &
                                      StandardName  = trim(StandardName),       & 
                                      Units         = trim(Units),              &
                                      ValidMin      = ValidMin,                 &
                                      ValidMax      = ValidMax,                 &
                                      MinValue      = MinValue,                 &
                                      MaxValue      = MaxValue,                 &
                                      OutputNumber  = iFinal,                   &
                                      Array2D       = Field2D,                  &
                                      STAT          = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
                
            endif

            deallocate(Field2D)
            
            

        else if (Me%Field(iP)%Dim==3) then

            allocate(Field3D(1:WorkJUB,1:WorkIUB,1:WorkKUB))   

            do k = 1, WorkKUB
            do j = 1, WorkJUB
            do i = 1, WorkIUB
                Field3D(j,i,k)=Me%Field(iP)%Value3DOut(i,j,k)
            enddo
            enddo                
            enddo

        
            call BuildAttributes(trim(Me%Field(iP)%ID%Name), NCDFName,                  &
                                      LongName, StandardName,                           &
                                      Units, ValidMin, ValidMax,                        &
                                      MinValue, MaxValue, MissingValue,                 &
                                      Float3D = Field3D)


            call NETCDFWriteData (NCDFID        = Me%NetCDF_Out%ObjNetCDF,              &
                                  Name          = trim(NCDFName),                       &
                                  LongName      = trim(LongName),                       &
                                  StandardName  = trim(StandardName),                   & 
                                  Units         = trim(Units),                          &
                                  ValidMin      = ValidMin,                             &
                                  ValidMax      = ValidMax,                             &
                                  MinValue      = MinValue,                             &
                                  MaxValue      = MaxValue,                             &
                                  OutputNumber  = iFinal,                               &
                                  Array3D       = Field3D,                              &
                                  STAT          = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR40'
            
            deallocate(Field3D)        

        else 

            stop 'WriteFieldNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR50'

        endif

    end subroutine WriteFieldNetCDF


    !----------------------------------------------------------------------

    !--------------------------------------------------------------------------

    
    subroutine KillNetCDFCF_2_HDF5MOHID
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers, ip
        
        !Begin-----------------------------------------------------------------
        
            if(associated(Me%Date%ValueInTotal)) deAllocate(Me%Date%ValueInTotal)            

            if (Me%OutHDF5) then
                !Close HDF5 File
                call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR10'
            endif
            
            if (Me%OutNetCDF) then
                
                !Close NetCDF File
                call KillNETCDF(Me%NetCDF_Out%ObjNetCDF, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR20'
            
            endif
            
            call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR25'               

            nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
            if (nUsers == 0) stop 'KillNetCDFCF_2_HDF5MOHID - ModuleNetCDFCF_2_HDF5MOHID - ERR30'

            if (associated(Me%Date%Value1DOut    )) deallocate(Me%Date%Value1DOut   )
            
            if (associated(Me%LongLat%LongOut    )) deallocate(Me%LongLat%LongOut   )
            if (associated(Me%LongLat%LatOut     )) deallocate(Me%LongLat%LatOut    )
            
            if (associated(Me%Bathym%Value2DOut  )) deallocate(Me%Bathym%Value2DOut )            

            if (associated(Me%LongLat%RotationX  )) deallocate(Me%LongLat%RotationX )
            if (associated(Me%LongLat%RotationY  )) deallocate(Me%LongLat%RotationY )
           
            
            if (associated(Me%Bathym%Value2DOut  )) deallocate(Me%Bathym%Value2DOut )
            
            if (associated(Me%Mapping%Value2DOut )) deallocate(Me%Mapping%Value2DOut)
            if (associated(Me%Mapping%Value3DOut )) deallocate(Me%Mapping%Value3DOut)

            do ip = 1, Me%PropNumber   
                if (associated(Me%Field(iP)%Value2DOut          )) deallocate(Me%Field(iP)%Value2DOut          )
                if (associated(Me%Field(iP)%Value3DOut          )) deallocate(Me%Field(iP)%Value3DOut          )
                if (associated(Me%Field(ip)%EnergyStartingHours )) deallocate(Me%Field(iP)%EnergyStartingHours )
                !if (associated(Me%Field(ip)%AvModelStartingHours)) deallocate(Me%Field(iP)%AvModelStartingHours)
            enddo
            deallocate(Me%Field)

            deallocate(Me)
            nullify   (Me)


    
    end subroutine KillNetCDFCF_2_HDF5MOHID

    !--------------------------------------------------------------------------

end module ModuleNetCDFCF_2_HDF5MOHID
