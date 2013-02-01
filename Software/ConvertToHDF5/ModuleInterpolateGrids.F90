!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : InterpolateGrids
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to Interpolate grids with data written in HDF5 format files
!
!------------------------------------------------------------------------------


Module ModuleInterpolateGrids

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap
    use ModuleFunctions
    use ModuleTriangulation
    use ModuleInterpolation
    use ModuleDrawing

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartInterpolateGrids
    private ::      ReadOptions
    private ::      ConstructNewGrid
    private ::      ConstructFatherGrid
    private ::      FatherSonCommunication
    private ::      OpenAndReadFatherHDF5File
    private ::      Open_HDF5_OutPut_File
    private ::      Open_HDF5_OutPut_File3D
    private ::          BilinearInterpolation
    private ::          Spline2DInterpolation
    private ::          Triangulator
    private ::          AveragePointsInCells
    private ::              FillingCells
    private ::          VerticalInterpolation
    private ::          ConstructFillingCells
    private ::          ModifyFillingCells
    private ::          KillFillingCells


    
    !Parameters----------------------------------------------------------------
    integer, parameter                                      :: Bilinear         = 1
    integer, parameter                                      :: Spline2D         = 2
    integer, parameter                                      :: Triangulation    = 3
    integer, parameter                                      :: AverageInCells   = 4

    integer, parameter                                      :: Null_             = 0
    integer, parameter                                      :: MediumTriang     = 1
    integer, parameter                                      :: HighTriang       = 2
    integer, parameter                                      :: NearestNeighbor  = 3
    integer, parameter                                      :: NearestCell      = 4
    integer, parameter                                      :: ConstantValue_   = 5
    

    !Types---------------------------------------------------------------------

    type       T_InterpolWindow
        real                                                :: Xmin, Ymin, Xmax, Ymax
    end type   T_InterpolWindow

    type       T_FillingCells
        integer, dimension(:),  pointer                     :: Mapping
        real,   dimension(:),   pointer                     :: Xin, Yin, Zin
        real,   dimension(:),   pointer                     :: Xout,Yout,Zout
        integer                                             :: NFather, NSon, MapN
    end type   T_FillingCells

    type       T_StationaryMap
        logical                                             :: ON
        type (T_FillingCells), dimension(:), pointer        :: FillingCells
    end type   T_StationaryMap

    
    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        integer                                             :: IDNumber
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        type(T_Field),              pointer                 :: Next
    end type  T_Field

    type       T_Grid
        integer                                             :: ObjHorizontalGrid    = 0
        integer                                             :: ObjHDF5              = 0
        integer                                             :: ObjHorizontalMap     = 0
        integer                                             :: ObjBathymetry        = 0
        integer                                             :: ObjGeometry          = 0
        integer                                             :: ObjMap               = 0
        integer                                             :: Unit
        character(len=PathLength)                           :: FileName
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: GeometryFileName
        type(T_Time), dimension(:), pointer                 :: InstantsArray
        integer                                             :: NumberOfInstants
        integer                                             :: NumberOfProperties
        type(T_Size2D)                                      :: Size2D
        type(T_Size3D)                                      :: Size3D
        type(T_Size2D)                                      :: WorkSize2D
        type(T_Size3D)                                      :: WorkSize3D
        integer, dimension(:,:),    pointer                 :: WaterPoints2D
        integer, dimension(:,:, :), pointer                 :: WaterPoints3D
        real,    dimension(:,:),    pointer                 :: ConnectionX, ConnectionY
    end type  T_Grid

    type       T_InterpolateGrids
        integer                                             :: ObjEnterData         = 0
        integer                                             :: ObjTime              = 0
        integer                                             :: ObjInterpolation     = 0
        integer                                             :: ObjTriangulation     = 0
        logical                                             :: Interpolation3D, InterpolateGrid3D, NewInterpolation
        integer                                             :: TypeOfInterpolation
        type(T_Grid )                                       :: Father
        type(T_Grid )                                       :: New
        type(T_Grid )                                       :: Aux
        type(T_Time)                                        :: BeginTime, EndTime
        type(T_InterpolWindow)                              :: InterpolWindow
        logical                                             :: TimeWindow = .true.
        logical                                             :: FirstProperty3D
        real,    dimension(:,:), pointer                    :: SonCenterX, SonCenterY
        integer, dimension(:,:), pointer                    :: WaterPoints2D
        integer                                             :: Count
        real,    dimension(:  ), pointer                    :: NodeX, NodeY, NodeZ
        logical                                             :: DoNotBelieveMap, PoliIsEven, ExtrapolateProfile
        integer                                             :: Extrapolate2DFields
        real                                                :: ExtrapolateLimit
        integer                                             :: PoliDegree
        character(len=StringLength)                         :: BaseGroup
        logical                                             :: DoNotBelieveTime
        integer                                             :: NumberSubGroups      = 0
        integer                                             :: NumberSubSubGroups   = 0
        logical                                             :: ConvertAllFields     = .true.
        character(len=StringLength), dimension(:), pointer  :: FieldsToInterpolate
        integer                                             :: nFieldsToInterpolate
        
        character(len=StringLength), dimension(:), pointer  :: SubGroups
        character(len=StringLength), dimension(:), pointer  :: SubSubGroups

        real                                                :: ExtrapolateValue
        type (T_StationaryMap)                              :: StationaryMap
    end type  T_InterpolateGrids

    type(T_InterpolateGrids), pointer                       :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartInterpolateGrids(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)
        
        Me%FirstProperty3D = .true. 

        call ReadOptions(ClientNumber)

        call ConstructFatherGrid

        call ConstructNewGrid

        if(Me%Interpolation3D) then

            call ConstructAuxGrid

            call FatherSonCommunication(Me%Aux)

        else

            call FatherSonCommunication(Me%New)

        endif


        if(Me%Interpolation3D) then

            call Open_HDF5_OutPut_File3D

        else 

            call Open_HDF5_OutPut_File

        endif

        if(Me%NewInterpolation)then

            call ConstructNewInterpolation

        end if

        if (Me%StationaryMap%ON .and. Me%TypeOfInterpolation == AverageInCells) then
            call ConstructFillingCells
        endif


        call OpenAndReadFatherHDF5File

        call KillFatherGrid

        call KillInterpolateGrids

        
        if (Me%StationaryMap%ON) then
            call KillFillingCells
        endif

        deallocate(Me)
        nullify(Me)

        STAT = SUCCESS_

    end subroutine StartInterpolateGrids

    !------------------------------------------------------------------------

    subroutine ConstructNewInterpolation
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call Construct_Interpolation3D(InterpolationID          = Me%ObjInterpolation,                          &   
                                       EnterDataID              = Me%ObjEnterData,                              &
                                       FromWhere                = FromFile,                                     &
                                       TimeSonID                = Me%ObjTime,                                   &
                                       HorizontalGridSonID      = Me%New%ObjHorizontalGrid,                     &
                                       HorizontalMapSonID       = Me%New%ObjHorizontalMap,                      &    
                                       GeometrySonID            = Me%New%ObjGeometry,                           &
                                       MapSonID                 = Me%New%ObjMap,                                &
                                       BathymetrySonID          = Me%New%ObjBathymetry,                         &
                                       HorizontalGridFatherID   = Me%Father%ObjHorizontalGrid,                  &
                                       HorizontalMapFatherID    = Me%Father%ObjHorizontalMap,                   &
                                       GeometryFatherID         = Me%Father%ObjGeometry,                        &
                                       MapFatherID              = Me%Father%ObjMap,                             &
                                       STAT                     = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructNewInterpolation - ModuleInterpolateGrids - ERR10'




    end subroutine ConstructNewInterpolation
    
    !------------------------------------------------------------------------
    
    subroutine ReadOptions(ClientNumber)


        !Arguments---------------------------------------------------------------
        integer, intent(IN )                        :: ClientNumber

        !Local-----------------------------------------------------------------
        real,   dimension(4)                        :: Aux4
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1

        !Begin-----------------------------------------------------------------
        
        write(*,*)'Reading instructions...'

        call GetData(Me%Father%FileName,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FATHER_FILENAME',                  &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR10'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleInterpolateGrids - ERR20'
        end if


        call GetData(Me%Father%GridFileName,                            &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FATHER_GRID_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR30'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleInterpolateGrids - ERR40'
        end if

        call GetData(Me%New%FileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR50'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleInterpolateGrids - ERR60'
        end if

        call GetData(Me%New%GridFileName,                               &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NEW_GRID_FILENAME',                &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR70'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleInterpolateGrids - ERR80'
        end if


        call GetData(Me%TypeOfInterpolation,                            &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'TYPE_OF_INTERPOLATION',            &
                     Default      = Bilinear,                           &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR90'

        if (iflag == 0)then
            write(*,*)'Must specify type of interpolation'
            stop 'ReadOptions - ModuleInterpolateGrids - ERR100'
        end if

        if(Me%TypeOfInterpolation == Bilinear)then
            write(*,*)
            write(*,*)'Type of interpolation : Bilinear'
            write(*,*)
        elseif(Me%TypeOfInterpolation == Spline2D)then
            write(*,*)
            write(*,*)'Type of interpolation : Spline2D'
            write(*,*)
        elseif(Me%TypeOfInterpolation == Triangulation)then
            write(*,*)
            write(*,*)'Type of interpolation : Triangulation'
            write(*,*)
        elseif(Me%TypeOfInterpolation == AverageInCells)then
            write(*,*)
            write(*,*)'Type of interpolation : Average by cell'
            write(*,*)
        else
            write(*,*) 'Unknown type of interpolation'
            stop       'ReadOptions - ModuleInterpolateGrids - ERR110' 
        end if


        !Aux4(1) =  FillValueReal
        !Aux4(2) =  FillValueReal
        !Aux4(3) = -FillValueReal
        !Aux4(4) = -FillValueReal

        call GetData(Aux4,                                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'INTERPOLATION_WINDOW',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR115'

        if (iflag == 0) then
            !(required because default in GetData has to be only one number not a vector)

            Aux4(1) =  FillValueReal
            Aux4(2) =  FillValueReal
            Aux4(3) = -FillValueReal
            Aux4(4) = -FillValueReal

        endif

        Me%InterpolWindow%Xmin = Aux4(1)
        Me%InterpolWindow%Ymin = Aux4(2)
        Me%InterpolWindow%Xmax = Aux4(3)
        Me%InterpolWindow%Ymax = Aux4(4)

        call GetData(Me%BeginTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR120'

        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR130'

        if (iflag==0 .and. iflag1==0) Me%TimeWindow = .FALSE.

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%BeginTime, Me%EndTime, 60., .false., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR140'

        call GetData(Me%Interpolation3D,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'INTERPOLATION3D',                  &
                     ClientModule = 'ConvertToHDF5',                    &
                     Default      = .false.,                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR150'

        call GetData(Me%InterpolateGrid3D,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'INTERPOLATION_GRID_3D',            &
                     ClientModule = 'ConvertToHDF5',                    &
                     Default      = .false.,                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR155'
        
        
        call GetData(Me%NewInterpolation,                               &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NEW_INTERPOLATION',                &
                     ClientModule = 'ConvertToHDF5',                    &
                     Default      = .false.,                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR150'

        


!------------------------------------------------------------------------------------------------------
        if(Me%Interpolation3D) then

            call GetData(Me%Father%GeometryFileName,                        &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromBlock,                          &
                         keyword      = 'FATHER_GEOMETRY',                  &
                         ClientModule = 'ConvertToHDF5',                    &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR160'


            call GetData(Me%New%GeometryFileName,                           &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromBlock,                          &
                         keyword      = 'NEW_GEOMETRY',                     &
                         ClientModule = 'ConvertToHDF5',                    &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR170'


            call GetData(Me%Aux%GridFileName,                               &
                         Me%ObjEnterData, iflag,                            &
                         SearchType   = FromBlock,                          &
                         keyword      = 'AUX_GRID_FILENAME',                &
                         ClientModule = 'ConvertToHDF5',                    &
                         Default      = trim(Me%New%GridFileName),          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR175'

            call GetData(Me%PoliDegree,                                                 &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'POLI_DEGREE',                                  &
                         default      = 1,                                              &
                         ClientModule = 'ConvertToHDF5',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR180'

            if (Me%InterpolateGrid3D .and. Me%PoliDegree /= 1) then

                write(*,*) 'If you want to interpolate the grid is because your grid is not standard'
                write(*,*) 'Press enter'
                read (*,*)
                write(*,*) 'If the grid is not standard is because is generic and can have a high degree of distoriton'
                write(*,*) 'Press enter'
                read (*,*)
                write(*,*) 'If the grid has a high degree of distortion to interpolate '
                write(*,*) 'in the verticaly with a polinomial is dangerous' 
                write(*,*) 'Press enter'
                read (*,*)
                write(*,*) 'It is recomended to interpolate lineary in the vertical in this case POLI_DEGREE : 1'
                write(*,*) 'Press enter'
                read (*,*)

            endif

            call GetData(Me%Aux%FileName,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'AUX_OUTPUTFILENAME',                           &
                         ClientModule = 'ConvertToHDF5',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR50'


        endif


        call GetData(Me%DoNotBelieveMap,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DO_NOT_BELIEVE_MAP',                               &
                     default      = .false.,                                            &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR190'


        call GetData(Me%Extrapolate2DFields,                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'EXTRAPOLATE_2D',                                   &
                     default      =  Null_,                                              &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR190'

        call GetData(Me%ExtrapolateLimit,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'EXTRAPOLATE_LIMIT',                                &
                     default      =  FillValueReal/4.,                                  &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR190'


        call GetData(Me%ExtrapolateProfile,                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'EXTRAPOLATE_PROFILE',                              &
                     default      =  .true.,                                            &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR200'



        call GetData(Me%BaseGroup,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BASE_GROUP',                                       &
                     default      =  "/Results",                                        &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR210'

        !To enable interpolation from  HDF5Statistics

        call GetData(Me%DoNotBelieveTime,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'DO_NOT_BELIEVE_TIME',                              &
                     default      = .false.,                                            &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR190'

        
        call GetData(Me%NumberSubGroups,                                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'NUM_SUB_GROUPS',                                   &
                     default      = 0,                                                  &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR210'

        if (Me%NumberSubGroups > 0) then

            allocate(Me%SubGroups(1:Me%NumberSubGroups))
            Me%SubGroups = null_str

            call GetData(Me%SubGroups,                                                      &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      = 'SUB_GROUPS',                                       &
                         ClientModule = 'ConvertToHDF5',                                    &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR210'

            call GetData(Me%NumberSubSubGroups,                                             &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromBlock,                                          &
                         keyword      = 'NUM_SUB_SUB_GROUPS',                               &
                         default      = 0,                                                  &
                         ClientModule = 'ConvertToHDF5',                                    &
                        STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR210'

            if (Me%NumberSubSubGroups > 0) then

                allocate(Me%SubSubGroups(1:Me%NumberSubSubGroups))
                Me%SubSubGroups = null_str

                call GetData(Me%SubSubGroups,                                                   &
                             Me%ObjEnterData, iflag,                                            &
                             SearchType   = FromBlock,                                          &
                             keyword      = 'SUB_SUB_GROUPS',                                   &
                            ClientModule = 'ConvertToHDF5',                                     &
                            STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR210'

            endif
        endif


        if (Me%Extrapolate2DFields == ConstantValue_) then

            call GetData(Me%ExtrapolateValue,                                           &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'EXTRAPOLATE_VALUE',                            &
                         default      =  FillValueReal,                                 &
                         ClientModule = 'ConvertToHDF5',                                &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR220'

        endif



        call GetData(Me%StationaryMap%ON,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'STATIONARY_MAPPING',                               &
                     default      =  .false.,                                           &
                     ClientModule = 'ConvertToHDF5',                                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleInterpolateGrids - ERR230'

        
        call ReadFieldsToConvert(ClientNumber)

    !------------------------------------------------------------------------------------------------------


    end subroutine ReadOptions


    subroutine ReadFieldsToConvert(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                             :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                             :: flag, STAT_CALL
        integer                             :: StartLine, EndLine, Count
        integer                             :: CurrentLineNumber
        logical                             :: BlockFound
        character(len=StringLength)         :: PropertyName

        !Begin-----------------------------------------------------------------
        

        call ExtractBlockFromBlock(Me%ObjEnterData,             &
                                    ClientNumber,               &
                                    '<<BeginFields>>',          &
                                    '<<EndFields>>',            &
                                    BlockFound,                 &
                                    STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleInterpolateGrids - ERR01'

        if(BlockFound)then
    
            Me%ConvertAllFields = .false.
    
            call GetBlockSize(Me%ObjEnterData,                      &
                              ClientNumber,                         &
                              StartLine,                            &
                              EndLine,                              &
                              FromBlockInBlock,                     &
                              STAT = STAT_CALL)

            Me%nFieldsToInterpolate = EndLine - StartLine - 1

            allocate(Me%FieldsToInterpolate(1:Me%nFieldsToInterpolate))
        
            Count = 1

            do CurrentLineNumber = StartLine + 1 , EndLine - 1

                call GetData(PropertyName,                          &
                             Me%ObjEnterData,                       &
                             flag,                                  &
                             SearchType  = FromBlock_,              &
                             Buffer_Line = CurrentLineNumber,       &
                             STAT         = STAT_CALL)        
                if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleInterpolateGrids - ERR02'

                Me%FieldsToInterpolate(Count) = PropertyName

                Count = Count + 1

            end do

            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ReadFieldsToConvert - ModuleInterpolateGrids - ERR03'
        
        else

            Me%ConvertAllFields = .true. 

        endif

    end subroutine ReadFieldsToConvert


    !--------------------------------------------------------------------------


    subroutine OpenAndReadFatherHDF5File

        !Local-----------------------------------------------------------------
        logical                                 :: exist
        integer                                 :: STAT_CALL
        integer                                 :: HDF5_READ

        !Local-----------------------------------------------------------------
        integer                                 :: CurrentInstant
        integer                                 :: StartInstant, EndInstant
        integer                                 :: CurrentProperty 
        type(T_Field), pointer                  :: NewFatherField, NewField, AuxField
        character(len=StringLength)             :: PropertyName
        integer                                 :: Rank
        integer                                 :: Count = 1, NewCurrentInstant
        integer                                 :: n
!        real,    dimension(:,:  ), pointer      :: InValues2D, OutValues2D
        integer, dimension(:,:,:), pointer      :: WaterPoints3D        
        logical                                 :: ExistGroup
        logical                                 :: FirstTime, ConvertThisField
        character(len=StringLength)             :: SubGroupName, SubSubGroupName, RootGroup

        !Begin-----------------------------------------------------------------
        
        nullify(NewFatherField, NewField, AuxField)


        allocate(Me%Father%WaterPoints2D(Me%Father%Size2D%ILB:Me%Father%Size2D%IUB,     &
                                         Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))

        Me%Father%WaterPoints2D(:,:) = 0

        if(Me%Interpolation3D) then
            call GetWaterPoints3D(Me%New%ObjMap, Me%New%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR10'


            allocate(Me%Aux%WaterPoints2D(Me%Aux%Size3D%ILB:Me%Aux%Size3D%IUB,        &
                                          Me%Aux%Size3D%JLB:Me%Aux%Size3D%JUB))

            Me%Aux%WaterPoints2D(:,:) = 0

            allocate(Me%Father%WaterPoints3D(Me%Father%Size3D%ILB:Me%Father%Size3D%IUB, &
                                             Me%Father%Size3D%JLB:Me%Father%Size3D%JUB, &
                                             Me%Father%Size3D%KLB:Me%Father%Size3D%KUB))

            allocate(Me%Aux%WaterPoints3D   (Me%Aux%Size3D%ILB:Me%Aux%Size3D%IUB,       &
                                             Me%Aux%Size3D%JLB:Me%Aux%Size3D%JUB,       &
                                             Me%Aux%Size3D%KLB:Me%Aux%Size3D%KUB))

            call GetWaterPoints3D(Me%Father%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR20'

            Me%Father%WaterPoints3D(:,:,:) = WaterPoints3D(:,:,:)

            call UnGetMap(Me%Father%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR30'


            call GetWaterPoints3D(Me%Aux%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR40'

            Me%Aux%WaterPoints3D   (:,:,:) = WaterPoints3D(:,:,:)

            call UnGetMap(Me%Aux%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR50'

        endif
        
        !Verifies if file exists
        inquire(FILE = Me%Father%FileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'HDF5 file does not exist'
            stop 'OpenAndReadHDF5File - ModuleInterpolateGrids - ERR60'
        endif

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (Me%Father%ObjHDF5, trim(Me%Father%FileName), HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR70'


        !Get number of instants
        call GetHDF5GroupNumberOfItems(Me%Father%ObjHDF5, "/Time", &
                                       Me%Father%NumberOfInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR80'

        write(*,*)'Number of instants in available data: ', Me%Father%NumberOfInstants

        allocate(Me%Father%InstantsArray(1:Me%Father%NumberOfInstants))

        !fill array with instants
        do CurrentInstant = 1, Me%Father%NumberOfInstants

            Me%Father%InstantsArray(CurrentInstant) = HDF5TimeInstant(CurrentInstant)

        end do

        !check time window
        if (Me%TimeWindow) then
            if(Me%Father%InstantsArray(1) .gt. Me%BeginTime)then
                write(*,*)'Data available starts after speficied date.'
                stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR90'
            end if

            if(Me%Father%InstantsArray(Me%Father%NumberOfInstants) .lt. Me%EndTime)then
                write(*,*)'Data available ends before speficied date.'
                stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR100'
            end if
        else
            Me%BeginTime = Me%Father%InstantsArray(1)
            Me%EndTime   = Me%Father%InstantsArray(Me%Father%NumberOfInstants)
        endif
        


        !select time window begin
        do CurrentInstant = 1, Me%Father%NumberOfInstants

            if(Me%Father%InstantsArray(CurrentInstant) .ge. Me%BeginTime)then

                StartInstant = CurrentInstant
                
                exit

            end if

        end do

        
        !select time window end
        do CurrentInstant = StartInstant, Me%Father%NumberOfInstants

            if(Me%Father%InstantsArray(CurrentInstant) .ge. Me%EndTime)then

                EndInstant = CurrentInstant

                exit 

            end if

        end do

        Me%New%NumberOfInstants = EndInstant - StartInstant + 1

        allocate(Me%New%InstantsArray(1:Me%New%NumberOfInstants))

        Me%New%InstantsArray = Me%Father%InstantsArray(StartInstant:EndInstant)


        !check number of properties in file
        !call GetHDF5GroupNumberOfItems(Me%Father%ObjHDF5, trim(Me%BaseGroup), &
        !                               Me%Father%NumberOfProperties, STAT = STAT_CALL)

        call GetHDF5GroupExist(Me%Father%ObjHDF5, trim(Me%BaseGroup), ExistGroup, &
                              nGroup = Me%Father%NumberOfProperties, STAT = STAT_CALL)
        
        
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR110'


        if(Me%Father%NumberOfProperties == 0)then
            write(*,*)'No data available in file: '//trim(Me%Father%FileName)
            stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR120'
        end if
        

        Me%New%NumberOfProperties = Me%Father%NumberOfProperties

        FirstTime = .true. 

prop:   do CurrentProperty = 1, Me%New%NumberOfProperties
          
            !get property name
            call GetHDF5GroupID(Me%Father%ObjHDF5, trim(Me%BaseGroup),        &
                                CurrentProperty, PropertyName,        &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR130'

            if (PropertyName =='VolumeCreated') cycle prop

            if(.not. Me%ConvertAllFields)then

                ConvertThisField = .false.

                do n = 1, Me%nFieldsToInterpolate

                    if(Me%FieldsToInterpolate(n) == PropertyName)then

                        ConvertThisField = .true.

                    endif

                end do

                if(.not. ConvertThisField) cycle prop

            end if


            write(*,*)'Reading '//trim(PropertyName)//' fields'

            NewCurrentInstant = 0

instant:    do CurrentInstant = StartInstant, EndInstant

                Count = Count + 1

                NewCurrentInstant = NewCurrentInstant + 1
                

                if (FirstTime) then
                    call OutputInstants(NewCurrentInstant)
                endif
                
                SubGroupName    = ""
                SubSubGroupName = ""

                if (Me%NumberSubGroups > 0) write(SubGroupName,*) "/", trim(Me%SubGroups(CurrentInstant))

ifS:            if (Me%NumberSubSubGroups > 0) then

doSS:               do n = 1,Me%NumberSubSubGroups 

                        write(SubSubGroupName, *) "/", trim(Me%SubSubGroups(n))
                        
                        RootGroup = trim(Me%BaseGroup)//"/"//trim(PropertyName)//&
                                    trim(adjustl(SubGroupName))//trim(adjustl(SubSubGroupName))
                        
                        call InterpolateGrids (RootGroup        = RootGroup,            &
                                               PropertyName     = PropertyName,         &
                                               CurrentInstant   = CurrentInstant,       &
                                               NewCurrentInstant= NewCurrentInstant,    &
                                               Count            = Count,                &
                                               FirstProperty3D  = Me%FirstProperty3D,   &
                                               Rank = Rank)
 

                    enddo doSS


                else !ifS

                    RootGroup = trim(Me%BaseGroup)//"/"//trim(PropertyName)//&
                                trim(adjustl(SubGroupName))//trim(adjustl(SubSubGroupName))

                    call InterpolateGrids (RootGroup        = RootGroup,                &
                                           PropertyName     = PropertyName,             &
                                           CurrentInstant   = CurrentInstant,           &
                                           NewCurrentInstant= NewCurrentInstant,        &
                                           Count            = Count,                    &
                                           FirstProperty3D  = Me%FirstProperty3D,       &   
                                           Rank = Rank)
                
                endif ifS


            end do instant

            FirstTime = .false.

            if (Rank == 3 .and. Me%FirstProperty3D) Me%FirstProperty3D = .false.

        end do prop


        deallocate(Me%Father%WaterPoints2D)

        if(Me%Interpolation3D) then     

            deallocate(Me%Aux%WaterPoints2D)

            deallocate(Me%Father%WaterPoints3D)
            deallocate(Me%Aux%WaterPoints3D)


            call UnGetMap(Me%New%ObjMap, Me%Aux%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR320'

        endif        


    end subroutine OpenAndReadFatherHDF5File

    !------------------------------------------------------------------------

    !------------------------------------------------------------------------

    subroutine InterpolateGrids (RootGroup, PropertyName,                               &
                                 CurrentInstant, NewCurrentInstant,                     &
                                 Count, FirstProperty3D, Rank)

        !Arguments ------------------------------------------------------------
        character(len=*)                        :: RootGroup, PropertyName
        integer                                 :: CurrentInstant, NewCurrentInstant, Count
        logical                                 :: FirstProperty3D
        integer, intent(OUT)                    :: Rank

        !Locals ---------------------------------------------------------------
        type(T_Field), pointer                  :: NewFatherField, NewField, AuxField                   
        integer                                 :: STAT_CALL, GroupPosition
        integer, dimension(7)                   :: Dimensions
        type(T_Grid ), pointer                  :: AuxGrid        
        integer, dimension(:,:  ), pointer      :: FatherWP2D, AuxWP2D
        integer                                 :: i,j,k, WetFace, kin
        real,    dimension(:,:  ), pointer      :: SurfaceElevation, InValues2D, OutValues2D
        type(T_Field), pointer                  :: AuxSZZ, NewFatherSZZ        


        !Begin ----------------------------------------------------------------


        !Allocates new instance
        allocate (NewField)

        !Allocates new instance
        allocate (NewFatherField)


        !Allocates aux instance
        allocate (AuxField)

        !Allocate local matrices
        if (Me%Extrapolate2DFields /= Null_) then

            allocate(InValues2D (Me%Father%Size2D%ILB:Me%Father%Size2D%IUB,Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))
            allocate(OutValues2D(Me%Aux%Size2D%ILB   :Me%Aux%Size2D%IUB   ,Me%Aux%Size2D%JLB   :Me%Aux%Size2D%JUB))

        end if

        !Construct new fields for father and son
        NewFatherField%IDNumber = Count
        Count                   = Count + 1

        NewField%IDNumber       = NewFatherField%IDNumber
        !NewField%Name           = trim(PropertyName)

        AuxField%IDNumber       = NewFatherField%IDNumber
        !AuxField%Name           = trim(PropertyName)

        !Get field ID
        if (Me%DoNotBelieveTime) then
            GroupPosition = 1
        else 
            GroupPosition = CurrentInstant
        endif


        call GetHDF5GroupID(Me%Father%ObjHDF5, RootGroup,               &                    
                            GroupPosition, NewFatherField%Name,         &
                            NewFatherField%Units, Rank, Dimensions,     &
                            STAT = STAT_CALL)                                
        if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR140'
   
   
        NewField%Name  = trim(NewFatherField%Name)
        NewField%Units = trim(NewFatherField%Units)
        NewField%Date  = Me%New%InstantsArray(NewCurrentInstant)

        AuxField%Name  = trim(NewFatherField%Name)
        AuxField%Units = trim(NewFatherField%Units)
        AuxField%Date  = Me%New%InstantsArray(NewCurrentInstant)

        select case (Rank)

            case(2)

                if(Me%Interpolation3D) then     
                    AuxGrid => Me%Aux
                else
                    AuxGrid => Me%New
                endif

                !check dimensions
                if(Dimensions(1) .ne. Me%Father%WorkSize2D%IUB) then
                    write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%FileName)
                    stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR150'
                end if

                if(Dimensions(2) .ne. Me%Father%WorkSize2D%JUB) then
                    write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%FileName)
                    stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR160'
                end if 
            
                !allocate field
                nullify (NewFatherField%Values2D)
                allocate(NewFatherField%Values2D(Me%Father%Size2D%ILB:Me%Father%Size2D%IUB, &
                                                 Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))
                nullify (NewField%Values2D)
                allocate(NewField%Values2D(AuxGrid%Size2D%ILB:AuxGrid%Size2D%IUB,           &
                                           AuxGrid%Size2D%JLB:AuxGrid%Size2D%JUB))

                NewField%Values2D(:,:) = FillValueReal
            
                call HDF5SetLimits (Me%Father%ObjHDF5,                          &
                                    Me%Father%WorkSize2D%ILB,                   &
                                    Me%Father%WorkSize2D%IUB,                   &
                                    Me%Father%WorkSize2D%JLB,                   &
                                    Me%Father%WorkSize2D%JUB,                   &
                                    STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR170'

                !read field
                call HDF5ReadData(Me%Father%ObjHDF5, RootGroup,                         &
!                                  NewFatherField%Name,                                  &
                                  trim(PropertyName),                                   &
                                  Array2D      = NewFatherField%Values2D,               &
                                  OutputNumber = CurrentInstant,                        &
                                  STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR180'

                if (Me%DoNotBelieveMap) then 
                    call FieldMapping (Me%Father, Values2D = NewFatherField%Values2D)
                else
                    call GetWaterPoints2D(Me%Father%ObjHorizontalMap, FatherWP2D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR184'

                    Me%Father%WaterPoints2D(:,:) = FatherWP2D(:,:)

                    call UnGetHorizontalMap(Me%Father%ObjHorizontalMap, FatherWP2D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR188'
                endif

                if (Me%Extrapolate2DFields == NearestCell) then
                    call ExtraPol2DFieldsNearestCell(Me%Father, NewFatherField%Values2D)  
                    Me%Father%WaterPoints2D(:,:) = 1                         
                endif    
                
                if(Me%NewInterpolation)then

                    !call ModifyInterpolator(Me%ObjInterpolation, NewField%Values2D, &
                    !                        NewFatherField%Values2D, STAT = STAT_CALL)
                    !if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR1880'

                else

                    if(Me%TypeOfInterpolation == Bilinear)then

                        call BilinearInterpolation(NewFatherField, NewField, AuxGrid)

                    elseif(Me%TypeOfInterpolation == Spline2D) then

                        call Spline2DInterpolation(NewFatherField, NewField, AuxGrid)

                    elseif(Me%TypeOfInterpolation == Triangulation) then

                        call GetWaterPoints2D(AuxGrid%ObjHorizontalMap, AuxWP2D, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR183'

                        if(Me%Interpolation3D) then     
                            Me%Aux%WaterPoints2D(:,:) = AuxWP2D(:,:)
                        else
                            Me%New%WaterPoints2D => AuxWP2D
                        endif
                        
                        call Triangulator               (NewFatherField, NewField, AuxGrid)

                        call UnGetHorizontalMap(AuxGrid%ObjHorizontalMap, AuxWP2D, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR186'

                    elseif(Me%TypeOfInterpolation == AverageInCells)then

                        call AveragePointsInCells (NewFatherField, NewField, AuxGrid)

                    else

                        write(*,*) 'Unknown type of interpolation'
                        stop       'StartInterpolateGrids - ModuleInterpolateGrids - ERR190' 

                    end if

                end if


            


                if (Me%Extrapolate2DFields /= Null_) then

                    call GetWaterPoints2D(AuxGrid%ObjHorizontalMap, AuxWP2D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR193'


                    if(Me%Interpolation3D) then     
                        Me%Aux%WaterPoints2D(:,:) = AuxWP2D(:,:)
                    else
                        Me%New%WaterPoints2D => AuxWP2D
                    endif
            
                    if      (Me%Extrapolate2DFields == MediumTriang) then
            
                        call ExtraPol2DFieldsTriang(NewFatherField, AuxGrid, .false., &
                                                    NewFatherField%Values2D, NewField%Values2D)
           
                    else if (Me%Extrapolate2DFields == HighTriang) then

                        call ExtraPol2DFieldsTriang(NewFatherField, AuxGrid, .true.,&
                                                    NewFatherField%Values2D, NewField%Values2D)

                    else if (Me%Extrapolate2DFields == NearestNeighbor) then

                        call ExtraPol2DFieldsNearest(NewFatherField, AuxGrid, &
                                                     NewFatherField%Values2D, NewField%Values2D)

                    else if (Me%Extrapolate2DFields == NearestCell) then

                        call ExtraPol2DFieldsNearestCell(AuxGrid, NewField%Values2D)

                    else if (Me%Extrapolate2DFields == ConstantValue_) then

                        call ExtraPol2DFieldsConstantValue(AuxGrid, NewField%Values2D)

                    endif

                    call UnGetHorizontalMap(AuxGrid%ObjHorizontalMap, AuxWP2D, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR196'

                    if (Me%DoNotBelieveMap) then                             
                        call FieldMapping (AuxGrid, Values2D = NewField%Values2D)
                    endif

                endif

       
                call OutputFields(RootGroup, NewField, PropertyName, NewCurrentInstant)

                nullify(AuxGrid)

            case(3)

                !if(trim(PropertyName) == 'salinity' .or. trim(PropertyName) == 'temperature') then
                !check dimensions
                if(Dimensions(1) .ne. Me%Father%WorkSize3D%IUB) then
                    write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%FileName)
                    stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR200'
                end if

                if(Dimensions(2) .ne. Me%Father%WorkSize3D%JUB) then
                    write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%FileName)
                    stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR210'
                end if 

                if(Dimensions(3) .ne. Me%Father%WorkSize3D%KUB) then
                    write(*,*)'Fields size is not consistent with grid size : '//trim(Me%Father%FileName)
                    stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR220'
                end if 
            
                !allocate field
                nullify (NewFatherField%Values2D)
                allocate(NewFatherField%Values2D(Me%Father%Size3D%ILB:Me%Father%Size3D%IUB,         &
                                                 Me%Father%Size3D%JLB:Me%Father%Size3D%JUB))

                NewFatherField%Values2D(:,:) = FillValueReal

                nullify (NewFatherField%Values3D)
                allocate(NewFatherField%Values3D(Me%Father%Size3D%ILB:Me%Father%Size3D%IUB,         &
                                                 Me%Father%Size3D%JLB:Me%Father%Size3D%JUB,         &
                                                 Me%Father%Size3D%KLB:Me%Father%Size3D%KUB))

                NewFatherField%Values3D(:,:,:) = FillValueReal


                nullify (AuxField%Values2D)
                allocate(AuxField%Values2D(Me%New%Size3D%ILB:Me%New%Size3D%IUB,                     &
                                           Me%New%Size3D%JLB:Me%New%Size3D%JUB))

                nullify (AuxField%Values3D)
                allocate(AuxField%Values3D(Me%New%Size3D%ILB:Me%New%Size3D%IUB,                     &
                                           Me%New%Size3D%JLB:Me%New%Size3D%JUB,                     &
                                           Me%Father%Size3D%KLB:Me%Father%Size3D%KUB))

                AuxField%Values3D(:,:,:) = FillValueReal


                nullify (NewField%Values2D)
                allocate(NewField%Values2D(Me%New%Size3D%ILB:Me%New%Size3D%IUB,                     &
                                           Me%New%Size3D%JLB:Me%New%Size3D%JUB))

                nullify (NewField%Values3D)
                allocate(NewField%Values3D(Me%New%Size3D%ILB:Me%New%Size3D%IUB,                     &
                                           Me%New%Size3D%JLB:Me%New%Size3D%JUB,                     &
                                           Me%New%Size3D%KLB:Me%New%Size3D%KUB))

                NewField%Values3D(:,:,:) = FillValueReal

            
                call HDF5SetLimits (Me%Father%ObjHDF5,                          &
                                    Me%Father%WorkSize3D%ILB,                   &
                                    Me%Father%WorkSize3D%IUB,                   &
                                    Me%Father%WorkSize3D%JLB,                   &
                                    Me%Father%WorkSize3D%JUB,                   &
                                    Me%Father%WorkSize3D%KLB,                   &
                                    Me%Father%WorkSize3D%KUB,                   &
                                    STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR230'

                !read field
                call HDF5ReadData(Me%Father%ObjHDF5, RootGroup,                     &
                                  NewFatherField%Name,                              &
                                  Array3D      = NewFatherField%Values3D,           &
!                                  OutputNumber = CurrentInstant,                   &
                                  STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR240'

                if (Me%DoNotBelieveMap) then 
                    call FieldMapping (Me%Father, Values3D = NewFatherField%Values3D)
                endif



                if(Me%NewInterpolation)then

                    call ModifyInterpolator(Me%ObjInterpolation, NewField%Values3D, &
                                            NewFatherField%Values3D, STAT = STAT_CALl)
                    if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR1880'

                else

                    do k=Me%Father%WorkSize3D%KLB,Me%Father%WorkSize3D%KUB

                        AuxField%Values2D(:,:)       = FillValueReal
                        NewField%Values2D(:,:)       = FillValueReal

                        NewFatherField%Values2D(:,:) = NewFatherField%Values3D (:,:,k)
                        Me%Father%WaterPoints2D(:,:) = Me%Father%WaterPoints3D (:,:,k)
                        Me%Aux%WaterPoints2D   (:,:) = Me%Aux%WaterPoints3D (:,:,k)
                        
                        if (Me%Extrapolate2DFields == NearestCell) then
                            call ExtraPol2DFieldsNearestCell(Me%Father, NewFatherField%Values2D)     
                            Me%Father%WaterPoints2D(:,:) = 1           
                        endif                           

                        if(Me%TypeOfInterpolation == Bilinear)then

                            call BilinearInterpolation(NewFatherField, AuxField, Me%Aux, k)

                        elseif(Me%TypeOfInterpolation == Spline2D)then

                            call Spline2DInterpolation(NewFatherField, AuxField, Me%Aux)

                        elseif(Me%TypeOfInterpolation == Triangulation)then

                            call Triangulator (NewFatherField, AuxField, Me%Aux)

                        else if (Me%TypeOfInterpolation == AverageInCells)then

                            call AveragePointsInCells (NewFatherField, AuxField, Me%Aux, k)

                        else

                            write(*,*) 'Unknown type of interpolation'
                            stop       'StartInterpolateGrids - ModuleInterpolateGrids - ERR250' 

                        end if

                        AuxField%Values3D       (:,:,k) = AuxField%Values2D      (:,:)

                    enddo

                end if

                if (Me%DoNotBelieveMap) then 
                    call FieldMapping (Me%Aux, Values3D = AuxField%Values3D)
                endif

                if (Me%Extrapolate2DFields /= Null_) then

                do k=Me%Father%WorkSize3D%KLB,Me%Father%WorkSize3D%KUB

                    InValues2D (:,:) = NewFatherField%Values3D(:,:,k)
                    OutValues2D(:,:) = AuxField%Values3D      (:,:,k)

                    Me%Father%WaterPoints2D(:,:) = Me%Father%WaterPoints3D (:,:,k)
                    Me%Aux%WaterPoints2D   (:,:) = Me%Aux%WaterPoints3D (:,:,k)

                    if      (Me%Extrapolate2DFields == MediumTriang) then
            
                        call ExtraPol2DFieldsTriang(NewFatherField, Me%Aux, .false., &
                                                    InValues2D, OutValues2D)
           
                    else if (Me%Extrapolate2DFields == HighTriang) then

                        call ExtraPol2DFieldsTriang(NewFatherField, Me%Aux, .true., &
                                                    InValues2D, OutValues2D)

                    else if (Me%Extrapolate2DFields == NearestNeighbor) then

                        call ExtraPol2DFieldsNearest(NewFatherField, Me%Aux, &
                                                     InValues2D, OutValues2D)

                    else if (Me%Extrapolate2DFields == NearestCell) then

                        call ExtraPol2DFieldsNearestCell(Me%Aux, OutValues2D)

                    else if (Me%Extrapolate2DFields == ConstantValue_) then

                        call ExtraPol2DFieldsConstantValue(Me%Aux, OutValues2D)

                    endif

                    AuxField%Values3D      (:,:,k) = OutValues2D(:,:)

                    if (Me%DoNotBelieveMap) then 
                        call FieldMapping (Me%Aux, Values3D = AuxField%Values3D)
                    endif

                enddo

                endif


ifG3D:          if (Me%InterpolateGrid3D .and. FirstProperty3D) then

                    nullify (NewFatherSZZ)
                    allocate(NewFatherSZZ)

                    nullify (AuxSZZ)
                    allocate(AuxSZZ)
    
                    NewFatherSZZ%Name           = "VerticalZ"
                    NewFatherSZZ%IDNumber       = 0

                    AuxSZZ%Name           = "VerticalZ"
                    AuxSZZ%IDNumber       = 0


                    !allocate SZZ fields
                    nullify (NewFatherSZZ%Values3D)
                    allocate(NewFatherSZZ%Values3D(Me%Father%Size3D%ILB:Me%Father%Size3D%IUB, &
                                                   Me%Father%Size3D%JLB:Me%Father%Size3D%JUB, &
                                                   Me%Father%Size3D%KLB:Me%Father%Size3D%KUB))

                    nullify (NewFatherSZZ%Values2D)
                    allocate(NewFatherSZZ%Values2D(Me%Father%Size3D%ILB:Me%Father%Size3D%IUB, &
                                                   Me%Father%Size3D%JLB:Me%Father%Size3D%JUB))

                    nullify (AuxSZZ%Values3D)
                    allocate(AuxSZZ%Values3D(Me%New%Size3D%ILB:Me%New%Size3D%IUB,             &
                                             Me%New%Size3D%JLB:Me%New%Size3D%JUB,             &
                                             Me%Father%Size3D%KLB:Me%Father%Size3D%KUB))

                    nullify (AuxSZZ%Values2D)
                    allocate(AuxSZZ%Values2D(Me%New%Size3D%ILB:Me%New%Size3D%IUB,             &
                                             Me%New%Size3D%JLB:Me%New%Size3D%JUB))


            
                    call HDF5SetLimits  (Me%Father%ObjHDF5,                     &
                                         Me%Father%WorkSize3D%ILB,              &
                                         Me%Father%WorkSize3D%IUB,              &
                                         Me%Father%WorkSize3D%JLB,              &
                                         Me%Father%WorkSize3D%JUB,              &
                                         Me%Father%WorkSize3D%KLB-1,            &
                                         Me%Father%WorkSize3D%KUB,              &
                                         STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR270'


    
                    !read field
                    call HDF5ReadData(Me%Father%ObjHDF5, "/Grid/"//"VerticalZ", &
                                      "Vertical",                               &
                                      Array3D      = NewFatherSZZ%Values3D,     &
                                      OutputNumber = CurrentInstant,            &
                                      STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR280'

                    do k=Me%Father%WorkSize3D%KLB-1,Me%Father%WorkSize3D%KUB

                        NewFatherSZZ%Values2D  (:,:) = NewFatherSZZ%Values3D   (:,:,k)

                        do j = Me%Father%WorkSize3D%JLB, Me%Father%WorkSize3D%JUB
                        do i = Me%Father%WorkSize3D%ILB, Me%Father%WorkSize3D%IUB

                            WetFace = Me%Father%WaterPoints3D (i,j,k) + Me%Father%WaterPoints3D (i,j,k+1) 

                            if (WetFace>0) then

                                Me%Father%WaterPoints2D(i,j) = 1

                            else

                                Me%Father%WaterPoints2D(i,j) = 0

                            endif

                        enddo
                        enddo

                        do j = Me%Aux%WorkSize3D%JLB, Me%Aux%WorkSize3D%JUB
                        do i = Me%Aux%WorkSize3D%ILB, Me%Aux%WorkSize3D%IUB

                            WetFace = Me%Aux%WaterPoints3D (i,j,k) + Me%Aux%WaterPoints3D (i,j,k+1) 

                            if (WetFace>0) then

                                Me%Aux%WaterPoints2D(i,j) = 1

                            else

                                Me%Aux%WaterPoints2D(i,j) = 0

                            endif

                        enddo
                        enddo


                        if(Me%TypeOfInterpolation == Bilinear)then

                            kin = k
                            if (k==0) kin = 1

                            call BilinearInterpolation(NewFatherSZZ, AuxSZZ, Me%Aux, kin)

                        elseif(Me%TypeOfInterpolation == Spline2D)then

                            call Spline2DInterpolation(NewFatherSZZ, AuxSZZ, Me%Aux)

                        elseif(Me%TypeOfInterpolation == Triangulation)then

                            call Triangulator (NewFatherSZZ, AuxSZZ, Me%Aux)

                        else if (Me%TypeOfInterpolation == AverageInCells)then

                            call AveragePointsInCells (NewFatherSZZ, AuxSZZ, Me%Aux, kin)

                        else

                            write(*,*) 'Unknown type of interpolation'
                            stop       'StartInterpolateGrids - ModuleInterpolateGrids - ERR290' 

                        end if

                        AuxSZZ%Values3D         (:,:,k) = AuxSZZ%Values2D        (:,:)

                    enddo


                    !deAllocates new instance
                    deallocate (NewFatherSZZ%Values3D)
                    nullify    (NewFatherSZZ%Values3D)

                    deallocate (NewFatherSZZ%Values2D)
                    nullify    (NewFatherSZZ%Values2D)


                    deallocate(NewFatherSZZ)
                    nullify   (NewFatherSZZ)

                    allocate(SurfaceElevation (Me%Aux%Size3D%ILB:Me%Aux%Size3D%IUB,         &
                                               Me%Aux%Size3D%JLB:Me%Aux%Size3D%JUB))
                
                    SurfaceElevation(:,:) = - AuxSZZ%Values3D(:,:,Me%Father%Size3D%KUB)

                                                                                            
                    call ComputeInitialGeometry(GeometryID       = Me%Aux%ObjGeometry,      &
                                                WaterPoints3D    = Me%Aux%WaterPoints3D,    &
                                                SurfaceElevation = SurfaceElevation,        &
                                                ActualTime       = Me%Father%InstantsArray(CurrentInstant), &
                                                SZZ              = AuxSZZ%Values3D,         &
                                                STAT             = STAT_CALL )
                    if(STAT_CALL .ne. SUCCESS_) stop 'OpenAndReadFatherHDF5File -  ModuleInterpolateGrids - ERR300'

                    deallocate(SurfaceElevation)


                    deallocate (AuxSZZ%Values3D)
                    nullify    (AuxSZZ%Values3D)

                    deallocate (AuxSZZ%Values2D)
                    nullify    (AuxSZZ%Values2D)

                    deallocate(AuxSZZ)
                    nullify   (AuxSZZ)


                endif ifG3D

                if ((Me%DoNotBelieveMap .or. Me%Extrapolate2DFields /= Null_) .and. .not. Me%InterpolateGrid3D) then
                !In this case water points mapping was changed

                    allocate(SurfaceElevation (Me%Aux%Size3D%ILB:Me%Aux%Size3D%IUB,      &
                                               Me%Aux%Size3D%JLB:Me%Aux%Size3D%JUB))
                    SurfaceElevation(:,:) = 0.

                    call ComputeInitialGeometry(GeometryID       = Me%Aux%ObjGeometry,   &
                                                WaterPoints3D    = Me%Aux%WaterPoints3D, &
                                                SurfaceElevation = SurfaceElevation,     &
                                                ActualTime       = Me%BeginTime,         &
                                                STAT             = STAT_CALL )
                    if(STAT_CALL .ne. SUCCESS_) stop 'OpenAndReadFatherHDF5File -  ModuleInterpolateGrids - ERR260'

                    deallocate(SurfaceElevation)

                endif


                if (FirstProperty3D) then
                    call OutputGrid3D(NewCurrentInstant)
                endif

            
                call VerticalInterpolation (AuxField, NewField)

                call OutputFields3D(RootGroup, NewField, PropertyName, NewCurrentInstant)

                !end if


            case default 
        
!                        write(*,*)'Interpolation only available for 2D fields.'
                stop 'OpenAndReadFatherHDF5File - ModuleInterpolateGrids - ERR310'
        
            end select

            !DeAllocates local matrices
            if (Me%Extrapolate2DFields /= Null_) then

                deallocate(InValues2D)
                deallocate(OutValues2D)

            end if

            !deAllocates new instance
            if(associated(NewField%Values2D)) deallocate (NewField%Values2D)
            if(associated(NewField%Values3D)) deallocate (NewField%Values3D)
            deallocate (NewField)
            nullify    (NewField)

            !deAllocates new instance
            if(associated(NewFatherField%Values2D)) deallocate (NewFatherField%Values2D)
            if(associated(NewFatherField%Values3D)) deallocate (NewFatherField%Values3D)
            deallocate (NewFatherField)
            nullify    (NewFatherField)

            !deAllocates new instance
            if(associated(AuxField%Values2D)) deallocate (AuxField%Values2D)
            if(associated(AuxField%Values3D)) deallocate (AuxField%Values3D)
            deallocate (AuxField)
            nullify    (AuxField)


    end subroutine InterpolateGrids

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine FieldMapping (Grid, Values3D, Values2D)

        !Arguments-------------------------------------------------------------
        type(T_Grid)                               :: Grid
        real, dimension(:,:,:), pointer, optional  :: Values3D
        real, dimension(:,:  ), pointer, optional  :: Values2D
        !Local-----------------------------------------------------------------
        integer                         :: i, j, k
        !----------------------------------------------------------------------

        do k = Grid%WorkSize3D%KLB, Grid%WorkSize3D%KUB
        do j = Grid%WorkSize3D%JLB, Grid%WorkSize3D%JUB
        do i = Grid%WorkSize3D%ILB, Grid%WorkSize3D%IUB

            if (present(Values3D)) then

                if (Values3D(i,j,k)   < FillValueReal / 2.) then
                                
                    Grid%WaterPoints3D        (i,j,k)   = 0

                else

                    Grid%WaterPoints3D        (i,j,k)   = 1

                endif

            endif


            if (present(Values2D)) then

                if (Values2D(i,j)   < FillValueReal / 2.) then
                                
                    Grid%WaterPoints3D        (i,j,k)   = 0

                else

                    Grid%WaterPoints3D        (i,j,k)   = 1

                endif

            endif
        enddo
        enddo
        enddo

        Grid%WaterPoints2D(:,:) = Grid%WaterPoints3D(:,:,Grid%WorkSize3D%KUB)

    end subroutine FieldMapping

    !--------------------------------------------------------------------------


    subroutine VerticalInterpolation( AuxField, NewField )

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: NewField, AuxField

        !Local-------------------------------------------------------------
        !real,    dimension(:,:,:), pointer      :: Values3D
        real,    dimension(:,:,:), pointer      :: ZCellCenter
        real,    dimension(:,:,:), pointer      :: ZCellCenterAux
        real(8), dimension(:    ), pointer      :: Depth, Values
        real(8)                                 :: INDepth, dz, Error
        integer                                 :: i,j,k,NDEPTHS, Aux, PoliDegree, STAT_CALL
        logical                                 :: PoliIsEven, FoundBottom, FoundSurface

        !Begin-----------------------------------------------------------------

        allocate(Depth (Me%Father%Size3D%KLB: Me%Father%Size3D%KUB))
        allocate(Values(Me%Father%Size3D%KLB: Me%Father%Size3D%KUB))

        call GetGeometryDistances(Me%New%ObjGeometry, ZCellCenter = ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR10'

        call GetGeometryDistances(Me%Aux%ObjGeometry, ZCellCenter = ZCellCenterAux, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR20'

        do i = Me%New%WorkSize3D%ILB, Me%New%WorkSize3D%IUB
        do j = Me%New%WorkSize3D%JLB, Me%New%WorkSize3D%JUB

            if (Me%New%WaterPoints3D(i, j, Me%New%WorkSize3D%KUB   ) == 1) then
            
            if (Me%Aux%WaterPoints3D(i, j, Me%Father%WorkSize3D%KUB) == 1) then

                Aux=Me%Father%WorkSize3D%KLB

                do k=Me%Father%WorkSize3D%KUB, Me%Father%WorkSize3D%KLB, -1
                    dz = abs(ZCellCenterAux(i,j, k) - ZCellCenterAux(i,j, k-1))
                    if(Me%Aux%WaterPoints3D(i, j, k) == 1                               &
                    !In the Hycom model instead of having a bottom mapping to the bottom cells are 
                    !given the same depth
                    .and.  (dz < 0.0001 .or. dz > -FillvalueReal/10.)) then
                        Aux = k
                        exit
                    endif
                enddo

                Depth (Aux:Me%Father%WorkSize3D%KUB) = - ZCellCenterAux    (i,j, Aux:Me%Father%WorkSize3D%KUB)
                Values(Aux:Me%Father%WorkSize3D%KUB) =   AuxField%Values3D (i,j, Aux:Me%Father%WorkSize3D%KUB)

                NDEPTHS = Me%Father%WorkSize3D%KUB - Aux  + 1

                if(Me%ExtrapolateProfile)then
                    do k = Me%Father%WorkSize3D%KUB-1, Aux, -1
                        if(Values(k) < FillValueReal/10.)then
                            Values(k) = Values(k+1)
                        end if
                    enddo
                endif

                do k=Me%New%WorkSize3D%KLB,Me%New%WorkSize3D%KUB

                    if (Me%New%WaterPoints3D(i, j, k) == 1) then

                        INDepth = - ZCellCenter (i,j,k)

                        if (Me%PoliDegree == 1) then

                            NewField%Values3D (i,j,k) = InterpolateProfileR8 (INDepth, NDEPTHS, &
                                                        Depth (Aux:Me%Father%WorkSize3D%KUB),   &
                                                        Values(Aux:Me%Father%WorkSize3D%KUB),   &
                                                        FoundBottom, FoundSurface)

                            !if (FoundBottom) then
                            !    write(*,*) "Found bottom - i, j, k ", i, j, k
                            !endif

                            !if (FoundSurface) then
                            !    write(*,*) "Found surface - i, j, k ", i, j, k
                            !endif
                        else


                            if ( NDEPTHS == 1) then
                                !Uniform profile is assumed when there only one layer
                                NewField%Values3D (i,j,k) = AuxField%Values3D (i,j, k)
                            else
                                !Interpolation n degree
                                PoliDegree = min (Me%PoliDegree, NDEPTHS-1)

                                if(IsOdd(PoliDegree))then
                                    PoliIsEven = .false.
                                else
                                    PoliIsEven = .true.
                                endif


                                NewField%Values3D (i,j,k) = PolIntProfile (INDepth, NDEPTHS,                     &
                                                                           Depth (Aux:Me%Father%WorkSize3D%KUB), &
                                                                           Values(Aux:Me%Father%WorkSize3D%KUB), &
                                                                           PoliDegree, PoliIsEven, Error)

                            endif

                        endif

                        if (.not. Me%ExtrapolateProfile) then

                            if (INDepth< Depth (Me%Father%WorkSize3D%KUB)) NewField%Values3D (i,j,k) = FillValueReal
                            if (INDepth> Depth (Aux                     )) NewField%Values3D (i,j,k) = FillValueReal

                        endif

                    else

                        NewField%Values3D (i,j,k) = FillValueReal

                    endif

                enddo

            else

                stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR30'                

            endif
            endif

        enddo
        enddo

        deallocate(Depth )
        deallocate(Values)

        call UnGetGeometry(Me%New%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR40'

        call UnGetGeometry(Me%Aux%ObjGeometry, ZCellCenter, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'VerticalInterpolation - ModuleInterpolateGrids - ERR50'


    end subroutine VerticalInterpolation


    !------------------------------------------------------------------------
    

    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%Father%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = Me%Father%ObjHDF5,                        &
                             GroupName      = '/Time',                                  &
                             Name           = 'Time',                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleInterpolateGrids - ERR01'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !------------------------------------------------------------------------

    subroutine ConstructNewGrid
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: exist
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real, dimension(:, :), pointer              :: SurfaceElevation

        !Begin-----------------------------------------------------------------


        write(*,*)'Constructing new grid...'
       
        !Verifies if file exists
        inquire(FILE = Me%New%GridFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructNewGrid - ModuleInterpolateGrids - ERR01'
        endif

        call ConstructHorizontalGrid(Me%New%ObjHorizontalGrid, Me%New%GridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid - ModuleInterpolateGrids - ERR02'


        call GetHorizontalGridSize(Me%New%ObjHorizontalGrid,                            &
                                   WorkSize = Me%New%WorkSize2D,                        &
                                   Size     = Me%New%Size2D,                            &
                                   STAT     = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR03'

        call ConstructGridData      (GridDataID       = Me%New%ObjBathymetry,        &
                                     HorizontalGridID = Me%New%ObjHorizontalGrid,    &
                                     FileName         = Me%New%GridFileName,         &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR04'

        call ConstructHorizontalMap (HorizontalMapID  = Me%New%ObjHorizontalMap,     &
                                     GridDataID       = Me%New%ObjBathymetry,        &
                                     HorizontalGridID = Me%New%ObjHorizontalGrid,    &
                                     ActualTime       = Me%BeginTime,                & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR05'

        if(Me%Interpolation3D) then

            call ConstructGeometry      (GeometryID       = Me%New%ObjGeometry,         &
                                         GridDataID       = Me%New%ObjBathymetry,       &
                                         HorizontalGridID = Me%New%ObjHorizontalGrid,   &
                                         HorizontalMapID  = Me%New%ObjHorizontalMap,    &
                                         ActualTime       = Me%BeginTime,               &
                                         NewDomain        = Me%New%GeometryFileName,    &
                                         STAT             = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR06'

            call GetGeometrySize(GeometryID     = Me%New%ObjGeometry,                   &
                                 Size           = Me%New%Size3D,                        &
                                 WorkSize       = Me%New%WorkSize3D,                    &
                                 STAT           = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR07'

            call ConstructMap ( Map_ID          = Me%New%ObjMap,                        &
                                GeometryID      = Me%New%ObjGeometry,                   &
                                HorizontalMapID = Me%new%ObjHorizontalMap,              &
                                TimeID          = Me%ObjTime,                           &
                                STAT            = STAT_CALL)  
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR08'

            allocate(SurfaceElevation (Me%New%Size3D%ILB:Me%New%Size3D%IUB,             &
                                       Me%New%Size3D%JLB:Me%New%Size3D%JUB))
            SurfaceElevation(:,:) = 0

            call GetWaterPoints3D(Me%New%ObjMap, WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'

            call ComputeInitialGeometry(GeometryID      = Me%New%ObjGeometry,           &
                                        WaterPoints3D   = WaterPoints3D,                &
                                        SurfaceElevation= SurfaceElevation,             &
                                        ActualTime      = Me%BeginTime,                 &
                                        STAT            = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructNewGrid -  ModuleInterpolateGrids - ERR06'

            deallocate(SurfaceElevation)

        endif
        
    end subroutine ConstructNewGrid

    !------------------------------------------------------------------------


    subroutine ConstructFatherGrid
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: exist
        real, dimension(:, :), pointer              :: SurfaceElevation

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Constructing father grid...'

        !Verifies if file exists
        inquire(FILE = Me%Father%GridFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructFatherGrid - ModuleInterpolateGrids - ERR10'
        endif

        call ConstructHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%GridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid - ModuleInterpolateGrids - ERR20'

        call GetHorizontalGridSize(Me%Father%ObjHorizontalGrid,                         &
                                   WorkSize = Me%Father%WorkSize2D,                     &
                                   Size     = Me%Father%Size2D,                         &
                                   STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR30'

        Me%Father%Size3D%ILB = Me%Father%Size2D%ILB
        Me%Father%Size3D%JLB = Me%Father%Size2D%JLB
        Me%Father%Size3D%IUB = Me%Father%Size2D%IUB
        Me%Father%Size3D%JUB = Me%Father%Size2D%JUB

        Me%Father%Size3D%KLB = 0
        Me%Father%Size3D%KUB = 2


        Me%Father%WorkSize3D%ILB = Me%Father%WorkSize2D%ILB
        Me%Father%WorkSize3D%JLB = Me%Father%WorkSize2D%JLB
        Me%Father%WorkSize3D%IUB = Me%Father%WorkSize2D%IUB
        Me%Father%WorkSize3D%JUB = Me%Father%WorkSize2D%JUB

        Me%Father%WorkSize3D%KLB = 1
        Me%Father%WorkSize3D%KUB = 1

        call ConstructGridData      (GridDataID       = Me%Father%ObjBathymetry,        &
                                     HorizontalGridID = Me%Father%ObjHorizontalGrid,    &
                                     FileName         = Me%Father%GridFileName,         &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR40'

        call ConstructHorizontalMap (HorizontalMapID  = Me%Father%ObjHorizontalMap,     &
                                     GridDataID       = Me%Father%ObjBathymetry,        &
                                     HorizontalGridID = Me%Father%ObjHorizontalGrid,    &
                                     ActualTime       = Me%BeginTime,                   & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR50'


        if(Me%Interpolation3D) then

            call ConstructGeometry      (GeometryID       = Me%Father%ObjGeometry,          &
                                         GridDataID       = Me%Father%ObjBathymetry,        &
                                         HorizontalGridID = Me%Father%ObjHorizontalGrid,    &
                                         HorizontalMapID  = Me%Father%ObjHorizontalMap,     &
                                         ActualTime       = Me%BeginTime,                   &
                                         NewDomain        = Me%Father%GeometryFileName,     &
                                         STAT             = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR60'

            call GetGeometrySize(GeometryID     = Me%Father%ObjGeometry,                &
                                 Size           = Me%Father%Size3D,                     &
                                 WorkSize       = Me%Father%WorkSize3D,                 &
                                 STAT           = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR70'

            call ConstructMap ( Map_ID          = Me%Father%ObjMap,                 &
                                GeometryID      = Me%Father%ObjGeometry,            &
                                HorizontalMapID = Me%Father%ObjHorizontalMap,       &
                                TimeID          = Me%ObjTime,                       &
                                STAT            = STAT_CALL)  
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR80'

            allocate(SurfaceElevation (Me%Father%Size3D%ILB:Me%Father%Size3D%IUB,   &
                                       Me%Father%Size3D%JLB:Me%Father%Size3D%JUB))
            SurfaceElevation(:,:) = 0.

            call GetWaterPoints3D(Me%Father%ObjMap, Me%Father%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFatherGrid - ModuleInterpolateGrids - ERR90'

            call ComputeInitialGeometry(GeometryID      = Me%Father%ObjGeometry,    &
                                        WaterPoints3D   = Me%Father%WaterPoints3D,  &
                                        SurfaceElevation= SurfaceElevation,         &
                                        ActualTime      = Me%BeginTime,             &
                                        STAT            = STAT_CALL )
            if(STAT_CALL .ne. SUCCESS_) stop 'ConstructFatherGrid -  ModuleInterpolateGrids - ERR100'

            call UnGetMap(Me%Father%ObjMap, Me%Father%WaterPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFatherGrid - ModuleInterpolateGrids - ERR110'

            deallocate(SurfaceElevation)


        endif

    end subroutine ConstructFatherGrid

    
    !------------------------------------------------------------------------


    subroutine ConstructAuxGrid

        !Local-------------------------------------------------------------
        real,       dimension(:,:  ), pointer   :: SurfaceElevation
        type(T_Polygon),              pointer   :: Polygon
        real                                    :: MaximumValue
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------


        write(*,*)'Constructing aux grid...'

        call ConstructHorizontalGrid(Me%Aux%ObjHorizontalGrid, Me%Aux%GridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid - ModuleInterpolateGrids - ERR10'


        call GetHorizontalGridSize(Me%Aux%ObjHorizontalGrid,                            &
                                   WorkSize = Me%Aux%WorkSize2D,                        &
                                   Size     = Me%Aux%Size2D,                            &
                                   STAT     = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR20'

        call ConstructGridData      (GridDataID       = Me%Aux%ObjBathymetry,           &
                                     HorizontalGridID = Me%Aux%ObjHorizontalGrid,       &
                                     FileName         = Me%Aux%GridFileName,            &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR40'


        allocate(Polygon     )
        Polygon%Count = 5
        allocate(Polygon%VerticesF(1:Polygon%Count))
        
        Polygon%VerticesF(1)%X = Me%InterpolWindow%Xmin
        Polygon%VerticesF(1)%Y = Me%InterpolWindow%Ymin
        Polygon%VerticesF(2)%X = Me%InterpolWindow%Xmin
        Polygon%VerticesF(2)%Y = Me%InterpolWindow%Ymax
        Polygon%VerticesF(3)%X = Me%InterpolWindow%Xmax
        Polygon%VerticesF(3)%Y = Me%InterpolWindow%Ymax
        Polygon%VerticesF(4)%X = Me%InterpolWindow%Xmax
        Polygon%VerticesF(4)%Y = Me%InterpolWindow%Ymin
        Polygon%VerticesF(5)%X = Polygon%VerticesF(1)%X 
        Polygon%VerticesF(5)%Y = Polygon%VerticesF(1)%Y 
        
        call SetLimits(Polygon)
       
        
        call GetMaxValueInPolygon(Me%Father%ObjBathymetry, MaximumValue, Polygon, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR42'    
        
        deallocate(Polygon%VerticesF)        
        deallocate(Polygon     )  
        
              

        call ModifyGridData (Me%Aux%ObjBathymetry, MaximumValue, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR44'    

        call ConstructHorizontalMap (HorizontalMapID  = Me%Aux%ObjHorizontalMap,        &
                                     GridDataID       = Me%Aux%ObjBathymetry,           &
                                     HorizontalGridID = Me%Aux%ObjHorizontalGrid,       &
                                     ActualTime       = Me%BeginTime,                   & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR50'

        call ConstructGeometry      (GeometryID       = Me%Aux%ObjGeometry,             &
                                     GridDataID       = Me%Aux%ObjBathymetry,           &
                                     HorizontalGridID = Me%Aux%ObjHorizontalGrid,       &
                                     HorizontalMapID  = Me%Aux%ObjHorizontalMap,        &
                                     ActualTime       = Me%BeginTime,                   &
                                     NewDomain        = Me%Father%GeometryFileName,     &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR60'

        call GetGeometrySize(GeometryID     = Me%Aux%ObjGeometry,                       &
                             Size           = Me%Aux%Size3D,                            &
                             WorkSize       = Me%Aux%WorkSize3D,                        &
                             STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR70'

        call ConstructMap ( Map_ID          = Me%Aux%ObjMap,                            &
                            GeometryID      = Me%Aux%ObjGeometry,                       &
                            HorizontalMapID = Me%Aux%ObjHorizontalMap,                  &
                            TimeID          = Me%ObjTime,                               &
                            STAT            = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR80'

        allocate(SurfaceElevation (Me%Aux%Size3D%ILB:Me%Aux%Size3D%IUB,                 &
                                   Me%Aux%Size3D%JLB:Me%Aux%Size3D%JUB))
        SurfaceElevation(:,:) = 0.

        call GetWaterPoints3D(Me%Aux%ObjMap, Me%Aux%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructAuxGrid - ModuleInterpolateGrids - ERR90'

        call ComputeInitialGeometry(GeometryID       = Me%Aux%ObjGeometry,              &
                                    WaterPoints3D    = Me%Aux%WaterPoints3D,            &
                                    SurfaceElevation = SurfaceElevation,                &
                                    ActualTime       = Me%BeginTime,                    &
                                    STAT             = STAT_CALL )
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructAuxGrid -  ModuleInterpolateGrids - ERR100'

        call UnGetMap(Me%Aux%ObjMap, Me%Aux%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructAuxGrid - ModuleInterpolateGrids - ERR120'

        deallocate(SurfaceElevation)


    end subroutine ConstructAuxGrid


    !------------------------------------------------------------------------
    
    subroutine FatherSonCommunication(NewGrid)

        !Arguments-------------------------------------------------------------
        type (T_Grid)                               :: NewGrid

        !Local-----------------------------------------------------------------
        logical                                     :: Distorted
        integer                                     :: STAT_CALL, CoordType, CoordTypeSon
        integer                                     :: UTM, MIL_PORT, GEOG, SIMPLE_GEOG
        integer                                     :: GRID_COORD, NLRD
        !Begin-----------------------------------------------------------------

        write(*,*)'Constructing communication between grids...'

        call GetCheckDistortion(Me%Father%ObjHorizontalGrid, Distorted, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR10'

        if (Distorted .and. Me%TypeOfInterpolation == Bilinear) then
            write(*,*) 'Cannot use Bilinear interpolation in distorted grids.'
            write(*,*) 'Change to Triangulation => TYPE_OF_INTERPOLATION : 3'
            stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR11'
        endif

        if(.not. Distorted)then

            call ConstructFatherGridLocation(NewGrid%ObjHorizontalGrid, Me%Father%ObjHorizontalGrid, &
                                             OkCross = .false., OkZ = .true., OkU = .false., OkV = .false., &
                                             STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR20'

        end if

        call GetCoordTypeList (GEOG = GEOG, UTM = UTM, MIL_PORT = MIL_PORT,             &
                               SIMPLE_GEOG = SIMPLE_GEOG, GRID_COORD = GRID_COORD,      &
                               NLRD = NLRD)

        !Gets Coordinates in use
        call GetGridCoordType(Me%Father%ObjHorizontalGrid, CoordType, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR30'

        call GetGridCoordType(NewGrid%ObjHorizontalGrid, CoordTypeSon, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR40'

        if (CoordType /= CoordTypeSon) then
            Write (*,*) 'Fathergrid coordinate type is different than son grid coordinate type'
            stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR50'
        endif

        if(CoordType == UTM .or. CoordType == MIL_PORT .or.                             &
           CoordType == GRID_COORD .or. CoordType == NLRD)then

            call GetHorizontalGrid(Me%Father%ObjHorizontalGrid,                        & 
                                   XX_IE = Me%Father%ConnectionX,                       &
                                   YY_IE = Me%Father%ConnectionY,                       &
                                   STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR60'

            call GetHorizontalGrid(NewGrid%ObjHorizontalGrid,                           &
                                   XX_IE = NewGrid%ConnectionX,                         &
                                   YY_IE = NewGrid%ConnectionY,                         &
                                   STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR70'

        else

            call GetGridLatitudeLongitude(Me%Father%ObjHorizontalGrid,                  &
                                          GridLatitudeConn  = Me%Father%ConnectionY,    &
                                          GridLongitudeConn = Me%Father%ConnectionX,    &
                                          STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR80'

            call GetGridLatitudeLongitude(NewGrid%ObjHorizontalGrid,                  &
                                          GridLatitudeConn  = NewGrid%ConnectionY,    &
                                          GridLongitudeConn = NewGrid%ConnectionX,    &
                                          STAT  = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'FatherSonCommunication - ModuleInterpolateGrids - ERR90'

        end if





    end subroutine FatherSonCommunication
    
    !------------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        real,       dimension(:,:  ), pointer       :: Bathymetry
        integer,    dimension(:,:,:), pointer       :: WaterPoints3D
        integer,    dimension(:,:  ), pointer       :: WaterPoints2D

        !----------------------------------------------------------------------

        allocate(WaterPoints3D(Me%New%Size2D%ILB:Me%New%Size2D%IUB,&
                               Me%New%Size2D%JLB:Me%New%Size2D%JUB, 1), STAT=STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'


        call GetWaterPoints2D(Me%New%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'
        
        WaterPoints3D(Me%New%WorkSize2D%ILB:Me%New%WorkSize2D%IUB, Me%New%WorkSize2D%JLB:Me%New%WorkSize2D%JUB,1) = &
        WaterPoints2D(Me%New%WorkSize2D%ILB:Me%New%WorkSize2D%IUB, Me%New%WorkSize2D%JLB:Me%New%WorkSize2D%JUB)

        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        
        !Opens HDF5 File
        call ConstructHDF5(Me%New%ObjHDF5, Me%New%FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR01'
        
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB, Me%New%WorkSize2D%IUB,&
                             Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'

        call GetGridData(Me%New%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'
        
        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR03'            

        call UngetGridData(Me%New%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'

        call WriteHorizontalGrid (Me%New%ObjHorizontalGrid, Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR04'

        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB, Me%New%WorkSize2D%IUB,&
                             Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB, 1,1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR05'

        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR06'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR07'

        call UngetHorizontalMap(Me%New%ObjHorizontalMap, WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine Open_HDF5_OutPut_File


    subroutine Open_HDF5_OutPut_File3D

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        real,       dimension(:,:  ), pointer       :: Bathymetry
        !integer,    dimension(:,:,:), pointer       :: WaterPoints3D
 !       integer,    dimension(:,:  ), pointer       :: WaterPoints2D

        !----------------------------------------------------------------------


        call GetWaterPoints3D(Me%New%ObjMap, Me%New%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR10'

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        

        !Opens HDF5 File
        call ConstructHDF5(Me%New%ObjHDF5, Me%New%FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR20'
        
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB, Me%New%WorkSize2D%IUB,&
                             Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR30'

        call GetGridData(Me%New%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR40'
        
        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR50'            

        call UngetGridData(Me%New%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR60'

        call WriteHorizontalGrid (Me%New%ObjHorizontalGrid, Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR70'

        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize3D%ILB, Me%New%WorkSize3D%IUB,&
                             Me%New%WorkSize3D%JLB, Me%New%WorkSize3D%JUB, Me%New%WorkSize3D%KLB, &
                             Me%New%WorkSize3D%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR80'

        call HDF5WriteData   (Me%New%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = Me%New%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR90'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR100'

        call UngetMap(Me%New%ObjMap, Me%New%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR110'



        call GetWaterPoints3D(Me%Aux%ObjMap, Me%Aux%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR120'

        !Opens HDF5 File
!        call ConstructHDF5(Me%Aux%ObjHDF5, Me%Aux%FileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR130'
        
        
!        call HDF5SetLimits  (Me%Aux%ObjHDF5, Me%Aux%WorkSize2D%ILB, Me%Aux%WorkSize2D%IUB,&
!                             Me%Aux%WorkSize2D%JLB, Me%Aux%WorkSize2D%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR140'

        call GetGridData(Me%Aux%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR150'
        
!        call HDF5WriteData   (Me%Aux%ObjHDF5, "/Grid", "Bathymetry", "-",       &
!                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR160'            

        call UngetGridData(Me%Aux%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR170'

!        call WriteHorizontalGrid (Me%Aux%ObjHorizontalGrid, Me%Aux%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR180'

!        call HDF5SetLimits  (Me%Aux%ObjHDF5, Me%Aux%WorkSize3D%ILB, Me%Aux%WorkSize3D%IUB,&
!                             Me%Aux%WorkSize3D%JLB, Me%Aux%WorkSize3D%JUB, Me%Aux%WorkSize3D%KLB, &
!                             Me%Aux%WorkSize3D%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR190'

!        call HDF5WriteData   (Me%Aux%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
!                              Array3D = Me%Aux%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR200'
        
        !Writes everything to disk
!        call HDF5FlushMemory (Me%Aux%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR210'

        call UngetMap(Me%Aux%ObjMap, Me%Aux%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File3D - ModuleInterpolateGrids - ERR220'

    end subroutine Open_HDF5_OutPut_File3D

    !------------------------------------------------------------------------

    subroutine OutputFields(RootGroup, NewField, PropertyName, OutputNumber)

        !Arguments-------------------------------------------------------------
        character(len=*)                                :: RootGroup, PropertyName
        type(T_Field), pointer                          :: NewField
        integer                                         :: OutputNumber

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize2D%ILB, Me%New%WorkSize2D%IUB,&
                             Me%New%WorkSize2D%JLB, Me%New%WorkSize2D%JUB, &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleInterpolateGrids - ERR05'


        call HDF5WriteData(Me%New%ObjHDF5, RootGroup,                                   &
!                           NewField%Name,                                               &
                           trim(PropertyName),                                          &
                           NewField%Units,                                              &
                           Array2D      = NewField%Values2D,                            &
                           OutputNumber = OutputNumber,                                 &
                           STAT         = STAT_CALL)


        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleInterpolateGrids - ERR07'

    end subroutine OutputFields
    
    !------------------------------------------------------------------------
    
    subroutine OutputFields3D(RootGroup, NewField, PropertyName, OutputNumber)

        !Arguments-------------------------------------------------------------
        character(len=*)                                :: RootGroup, PropertyName
        type(T_Field), pointer                          :: NewField
        integer                                         :: OutputNumber
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%New%ObjHDF5, Me%New%WorkSize3D%ILB, Me%New%WorkSize3D%IUB,&
                             Me%New%WorkSize3D%JLB, Me%New%WorkSize3D%JUB, &
                             Me%New%WorkSize3D%KLB, Me%New%WorkSize3D%KUB, &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateGrids - ERR10'


        call HDF5WriteData(Me%New%ObjHDF5, RootGroup,                                   &
!                           NewField%Name,                                               &
                           trim(PropertyName),                                          &
                           NewField%Units,                                              &
                           Array3D      = NewField%Values3D,                            &
                           OutputNumber = OutputNumber,                                 &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateGrids - ERR20'


        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateGrids - ERR30'

!        call HDF5SetLimits  (Me%Aux%ObjHDF5, Me%Aux%WorkSize3D%ILB, Me%Aux%WorkSize3D%IUB,&
!                             Me%Aux%WorkSize3D%JLB, Me%Aux%WorkSize3D%JUB, &
!                             Me%Aux%WorkSize3D%KLB, Me%Aux%WorkSize3D%KUB, &
!                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateGrids - ERR40'

!        call HDF5WriteData(Me%Aux%ObjHDF5,                             &
!                           trim(Me%BaseGroup)//"/"//AuxField%Name,     &
!                           AuxField%Name,                              &
!                           AuxField%Units,                             &
!                           Array3D      = AuxField%Values3D,           &
!                           OutputNumber = OutputNumber,                &
!                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateGrids - ERR50'



        !Writes everything to disk
!        call HDF5FlushMemory (Me%Aux%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields3D - ModuleInterpolateGrids - ERR60'

    end subroutine OutputFields3D


    !------------------------------------------------------------------------
    
    subroutine OutputGrid3D(OutputNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: OutputNumber
        
        !Local-----------------------------------------------------------------
        real   , dimension(:,:,:), pointer              :: SZZ
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (Me%New%ObjHDF5,                                            &
                             Me%New%WorkSize3D%ILB,   Me%New%WorkSize3D%IUB,            &
                             Me%New%WorkSize3D%JLB,   Me%New%WorkSize3D%JUB,            &
                             Me%New%WorkSize3D%KLB-1, Me%New%WorkSize3D%KUB,            &
                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR10'

        call GetGeometryDistances(Me%New%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR20'


        call HDF5WriteData(Me%New%ObjHDF5,                                              &
                           "/Grid/VerticalZ",                                           &
                           "Vertical",                                                  &
                           "m",                                                         &
                           Array3D      = SZZ,                                          &
                           OutputNumber = OutputNumber,                                 &
                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR30'


        !Writes everything to disk
        call HDF5FlushMemory (Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR40'

        call UnGetGeometry(Me%New%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR50'


!        call HDF5SetLimits  (Me%Aux%ObjHDF5,                                            &
!                             Me%Aux%WorkSize3D%ILB,   Me%Aux%WorkSize3D%IUB,            &
!                             Me%Aux%WorkSize3D%JLB,   Me%Aux%WorkSize3D%JUB,            &
!                             Me%Aux%WorkSize3D%KLB-1, Me%Aux%WorkSize3D%KUB,            &
!                             STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR60'

        call GetGeometryDistances(Me%Aux%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR70'


!        call HDF5WriteData(Me%Aux%ObjHDF5,                                              &
!                           "/Grid/VerticalZ",                                           &
!                           "Vertical",                                                  &
!                           "m",                                                         &
!                           Array3D      = SZZ,                                          &
!                           OutputNumber = OutputNumber,                                 &
!                           STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR80'


        !Writes everything to disk
!        call HDF5FlushMemory (Me%Aux%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR90'

        call UnGetGeometry(Me%Aux%ObjGeometry, SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputGrid3D - ModuleInterpolateGrids - ERR100'


    end subroutine OutputGrid3D


    !--------------------------------------------------------------------------



    subroutine OutputInstants(CurrentInstant)
        !Arguments-------------------------------------------------------------
        integer                                         :: CurrentInstant

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL
        type(T_Time)                                    :: CurrentDate


        !Begin-----------------------------------------------------------------
        
        CurrentDate = Me%New%InstantsArray(CurrentInstant)

        call ExtractDate   (CurrentDate,                                                &
                            AuxTime(1), AuxTime(2), AuxTime(3),                         &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%New%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleInterpolateGrids - ERR10'


        call HDF5WriteData  (Me%New%ObjHDF5, "/Time",                                   &
                             "Time", "YYYY/MM/DD HH:MM:SS",                             &
                             Array1D = TimePtr,                                         &
                             OutputNumber = CurrentInstant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleInterpolateGrids - ERR20'


        if (Me%Interpolation3D) then


!            call HDF5SetLimits  (Me%Aux%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleInterpolateGrids - ERR30'


!            call HDF5WriteData  (Me%Aux%ObjHDF5, "/Time",                               &
!                                 "Time", "YYYY/MM/DD HH:MM:SS",                         &
!                                 Array1D = TimePtr,                                     &
!                                 OutputNumber = CurrentInstant, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputInstants - ModuleInterpolateGrids - ERR40'

        endif

    end subroutine OutputInstants

    
    !------------------------------------------------------------------------


    subroutine BilinearInterpolation(FatherField, NewField, NewGrid, k)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: FatherField, NewField
        type (T_Grid)                           :: NewGrid
        integer, optional                       :: k

        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ComputeZ, k_

        !----------------------------------------------------------------------

        write(*,*)'Interpolating : '//trim(FatherField%Name), FatherField%IDNumber

        call GetComputeZUV(NewGrid%ObjHorizontalGrid, ComputeZ = ComputeZ, STAT = STAT_CALL)

        if (present(k)) then

            k_ = k

        else

            k_ = 1

        endif

        if(Me%Interpolation3D)then

            call InterpolRegularGrid(NewGrid%ObjHorizontalGrid,               &
                                     Me%Father%ObjHorizontalGrid,             &
                                     FatherField%Values2D,                    &
                                     NewField%Values2D,                       &
                                     Compute = ComputeZ,                      &
                                     KUBFather = k_,                          &
                                     ComputeFather = Me%Father%WaterPoints3D, &
                                     STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)stop 'BilinearInterpolation - ModuleInterpolateGrids - ERR01'

        else

            call InterpolRegularGrid(NewGrid%ObjHorizontalGrid,               &
                                     Me%Father%ObjHorizontalGrid,             &
                                     FatherField%Values2D,                    &
                                     NewField%Values2D,                       &
                                     Compute = ComputeZ,                      &
                                     KUBFather = k_,                          &
                                     STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)stop 'BilinearInterpolation - ModuleInterpolateGrids - ERR01'

        end if


    end subroutine BilinearInterpolation

    
    !------------------------------------------------------------------------


    subroutine Spline2DInterpolation(FatherField, NewField, NewGrid)

        !Arguments-------------------------------------------------------------
        type(T_Field), pointer                  :: FatherField, NewField
        type (T_Grid)                           :: NewGrid

        !Local-----------------------------------------------------------------
        real, dimension(:  ),       pointer             :: Father_XX_Z, Father_YY_Z
        real, dimension(:  ),       pointer             :: New_XX_Z, New_YY_Z
        real                                            :: Xorig, Yorig, GridAngle
        integer                                         :: STAT_CALL, i, j
        real, dimension(:,:),       pointer             :: function_int
        real, dimension(:  ),       pointer             :: Father_X, Father_Y
        real, dimension(:  ),       pointer             :: New_X, New_Y

        !----------------------------------------------------------------------
        
        call GetHorizontalGrid(Me%Father%ObjHorizontalGrid, &
                               XX_Z = Father_XX_Z,          &
                               YY_Z = Father_YY_Z,          &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR10'

        call GetGridOrigin(Me%Father%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR20'
        
        call GetGridAngle(Me%Father%ObjHorizontalGrid, GridAngle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR30'

        allocate(Father_X(Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))
        allocate(Father_Y(Me%Father%Size2D%ILB:Me%Father%Size2D%IUB))

        do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB
        do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB

            Father_X(j) = Father_XX_Z(j)
            Father_Y(i) = Father_YY_Z(i)

            call RodaXY(Xorig, Yorig, GridAngle, Father_X(j), Father_Y(i))

        end do
        end do

        allocate(New_X   (NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))
        allocate(New_Y   (NewGrid%Size2D%ILB:NewGrid%Size2D%IUB))
        
        call GetGridOrigin(NewGrid%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR40'
        
        call GetGridAngle(NewGrid%ObjHorizontalGrid, GridAngle, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR50'

        call GetHorizontalGrid(NewGrid%ObjHorizontalGrid,                               &
                               XX_Z = New_XX_Z,                                         &
                               YY_Z = New_YY_Z,                                         &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR60'

        do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB
        do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB

            New_X(j) = New_XX_Z(j)
            New_Y(i) = New_YY_Z(i)

            call RodaXY(Xorig, Yorig, GridAngle, New_X(j), New_Y(i))

        end do
        end do

        allocate(function_int(Me%Father%Size2D%ILB:Me%Father%Size2D%IUB,                &
                              Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))


        write(*,*)'Interpolating : '//trim(FatherField%Name), FatherField%IDNumber

        call splie2(npmax   = max(Me%Father%WorkSize2D%IUB,Me%Father%WorkSize2D%JUB),   &
                    x1a     = Father_Y,                                                 &
                    x2a     = Father_X,                                                 &
                    ya      = FatherField%Values2D,                                     &
                    m       = Me%Father%WorkSize2D%IUB,                                 &
                    n       = Me%Father%WorkSize2D%JUB,                                 &
                    y2a     = function_int)
               
        do i = NewGrid%WorkSize2D%ILB,  NewGrid%WorkSize2D%IUB
        do j = NewGrid%WorkSize2D%JLB , NewGrid%WorkSize2D%JUB


            call splin2(npmax   = max(Me%Father%WorkSize2D%IUB,Me%Father%WorkSize2D%JUB),&
                        x1a     = Father_Y,                                             &
                        x2a     = Father_X,                                             &
                        ya      = FatherField%Values2D,                                 &
                        y2a     = function_int,                                         &
                        m       = Me%Father%WorkSize2D%IUB,                             &
                        n       = Me%Father%WorkSize2D%JUB,                             &
                        x1      = New_Y(i),                                             &
                        x2      = New_X(j),                                             &
                        yc      = NewField%Values2D(i,j))
           
        enddo
        enddo
        
        deallocate(Father_X, Father_Y, New_X, New_Y, function_int)
        nullify   (Father_X, Father_Y, New_X, New_Y, function_int)

        call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Father_XX_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR80'

        call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Father_YY_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR90'

        call UnGetHorizontalGrid(NewGrid%ObjHorizontalGrid, Father_XX_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR100'

        call UnGetHorizontalGrid(NewGrid%ObjHorizontalGrid, Father_YY_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Spline2DInterpolation - ModuleInterpolateGrids - ERR110'


    end subroutine Spline2DInterpolation
    
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------


    subroutine splie2(npmax,x1a,x2a,ya,m,n,y2a)
        
        !Arguments-------------------------------------------------------------
        integer                         :: npmax, m,n
        real, dimension(:,:), pointer   :: y2a
        real, dimension(:,:), pointer   :: ya
        real, dimension(:  ), pointer   :: x2a
        real, dimension(:  ), pointer   :: x1a

        !Local-----------------------------------------------------------------
        real                            :: y2tmp(npmax),ytmp(npmax)
        integer                         :: j, k
        
        !Begin-----------------------------------------------------------------
        
        !Just to avoid compiler warning
        x1a = x1a

        do j = 1, m
            
            do k = 1, n
                ytmp(k)=ya(j,k)
            end do

            call spline(npmax, x2a, ytmp, n, 1.e30, 1.e30, y2tmp)
            
            do k=1, n
                y2a(j,k) = y2tmp(k)
            enddo

        enddo
      
    end subroutine splie2

    !------------------------------------------------------------------------
    
    subroutine spline(npmax,xb,yb,n,yp1,ypn,y2)
        
        !Arguments-------------------------------------------------------------
        integer                     :: npmax, n
        real                        :: yp1,ypn,yb(npmax),y2(npmax)

        
        !Local-----------------------------------------------------------------
        real, dimension(:), pointer :: xb
        real                        :: p,qn,sig,un,uu(npmax)
        integer                     :: i, k
        !Begin-----------------------------------------------------------------

        if (yp1.gt..99e30) then
            y2(1) = 0.
            uu(1) = 0.
        else
            y2(1) = -0.5
            uu(1) = (3./(xb(2)-xb(1)))*((yb(2)-yb(1))/(xb(2)-xb(1))-yp1)
        endif
      
        do i=2,n-1
            sig     = (xb(i)-xb(i-1))/(xb(i+1)-xb(i-1))
            p       = sig*y2(i-1)+2.
            y2(i)   = (sig-1.)/p
            uu(i)   = (6.*((yb(i+1)-yb(i))/(xb(i+1)-xb(i))-(yb(i)-yb(i-1))/ &
                      (xb(i)-xb(i-1)))/(xb(i+1)-xb(i-1))-sig* uu(i-1))/p
        enddo


        if (ypn.gt..99e30) then
            qn = 0.
            un = 0.
        else
            qn = 0.5
            un = (3./(xb(n)-xb(n-1)))*(ypn-(yb(n)-yb(n-1))/(xb(n)-xb(n-1)))
        endif

        y2(n) =(un-qn*uu(n-1))/(qn*y2(n-1)+1.)
        
        do k=n-1,1,-1
            y2(k) = y2(k)*y2(k+1)+uu(k)
        enddo

    end subroutine spline

    !------------------------------------------------------------------------

    subroutine splint_local(npmax,xa,ya,y2a,n,xc,yc)
        
        !Arguments-------------------------------------------------------------
        integer                     :: npmax, n
        real                        :: xc,yc
        real, dimension(:), pointer :: xa
        real                        :: y2a(npmax),ya(npmax)
        
        !Local-----------------------------------------------------------------
        real                        :: a,b,hh
        integer                     :: klo, khi,k

        !Begin-----------------------------------------------------------------


        klo = 1
        khi = n

        do while(khi-klo.gt.1)

            k = (khi + klo)/2
            
            if(xa(k).gt.xc)then
                khi = k
            else
                klo = k
            endif

        end do

        hh = xa(khi)-xa(klo)

        if (hh.eq.0.) stop 'bad xa input in splint_local - ModuleInterpolateGrids'
        
        a   = (xa(khi)-xc)/hh
        b   = (xc-xa(klo))/hh
        
        yc  = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(hh**2)/6.
        

    end subroutine splint_local

    !------------------------------------------------------------------------

    subroutine splin2(npmax,x1a,x2a,ya,y2a,m,n,x1,x2,yc)

        !Arguments-------------------------------------------------------------
        integer                         :: npmax, m, n
        real                            :: x1,x2,yc
        real, dimension(:,:), pointer   :: ya,y2a
        real, dimension(:  ), pointer   :: x2a
        real, dimension(:  ), pointer   :: x1a
        
        !Local-----------------------------------------------------------------
        real                            :: y2tmp(npmax),ytmp(npmax),yytmp(npmax)
        integer                         :: j, k

        !Begin-----------------------------------------------------------------

        do j=1,m
            do k=1,n
                ytmp(k) = ya(j,k)
                y2tmp(k)= y2a(j,k)
            enddo

            call splint_local(npmax,x2a,ytmp,y2tmp,n,x2,yytmp(j))
        enddo

        call spline(npmax,x1a,yytmp,m,1.e30,1.e30,y2tmp)
        call splint_local(npmax,x1a,yytmp,y2tmp,m,x1,yc)
           

    end subroutine splin2
    
    
    !------------------------------------------------------------------------
    
    
    subroutine KillFatherGrid
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Killing father grid...'


        call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%ConnectionX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillFatherGrid - ModuleInterpolateGrids - ERR10'


        call UnGetHorizontalGrid(Me%Father%ObjHorizontalGrid, Me%Father%ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillFatherGrid - ModuleInterpolateGrids - ERR20'

        if(Me%Interpolation3D) then

            call KillMap( Me%Father%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR30'

            call KillGeometry( Me%Father%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR40'

        endif

        call KillHorizontalMap( Me%Father%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR50'

        call KillGridData(Me%Father%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR60'

        call KillHDF5(Me%Father%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR70'

        call KillHorizontalGrid(Me%Father%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillFatherGrid - ModuleInterpolateGrids - ERR80'


    end subroutine KillFatherGrid

    !------------------------------------------------------------------------

    subroutine KillInterpolateGrids
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------


        if (Me%Interpolation3D) then

            call UnGetHorizontalGrid(Me%Aux%ObjHorizontalGrid, Me%Aux%ConnectionX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR10'


            call UnGetHorizontalGrid(Me%Aux%ObjHorizontalGrid, Me%Aux%ConnectionY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR20'

        else

            call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid, Me%New%ConnectionX, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR30'


            call UnGetHorizontalGrid(Me%New%ObjHorizontalGrid, Me%New%ConnectionY, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR40'

        endif

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR50' 

        if(Me%Interpolation3D) then

            call KillMap(Me%New%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR60'

            call KillGeometry(Me%New%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR70'

        endif

        call KillHorizontalMap(Me%New%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR80'

        call KillGridData(Me%New%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR90'

        call KillHorizontalGrid(Me%New%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR100'
        
        call KillHDF5(Me%New%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR110'

        
        if (Me%Interpolation3D) then
            call KillMap(Me%Aux%ObjMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR120'

            call KillGeometry(Me%Aux%ObjGeometry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR130'

            call KillHorizontalMap(Me%Aux%ObjHorizontalMap, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR140'

            call KillGridData(Me%Aux%ObjBathymetry, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR150'

            call KillHorizontalGrid(Me%Aux%ObjHorizontalGrid, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR160'

!            call KillHDF5(Me%Aux%ObjHDF5, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'KillInterpolateGrids - ModuleInterpolateGrids - ERR170'

        endif

   
    end subroutine KillInterpolateGrids

 

    subroutine Triangulator (FatherField, NewField, NewGrid)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField, NewField
        type(T_Grid )                               :: NewGrid

        !Local-----------------------------------------------------------------
        real,    dimension(:), pointer              :: NodeX, NodeY, NodeZ
        real                                        :: AuxX, AuxY
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, Count, i, j
        logical                                     :: FillOutsidePoints   = .false.
        !integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        
        !Begin-----------------------------------------------------------------


        NumberOfNodes =  Sum(Me%Father%WaterPoints2D(Me%Father%WorkSize2D%ILB:Me%Father%WorkSize2D%IUB, &
                                                     Me%Father%WorkSize2D%JLB:Me%Father%WorkSize2D%JUB))

iN:     if (NumberOfNodes >= 3) then

        allocate(NodeX(NumberOfNodes))
        allocate(NodeY(NumberOfNodes))
        allocate(NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
        do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB

            AuxX = ((Me%Father%ConnectionX(i, j  ) + Me%Father%ConnectionX(i+1, j  ))/2. + &
                               (Me%Father%ConnectionX(i, j+1) + Me%Father%ConnectionX(i+1, j+1))/2.)/2.
    
            AuxY = ((Me%Father%ConnectionY(i, j  ) + Me%Father%ConnectionY(i+1, j  ))/2. + &
                               (Me%Father%ConnectionY(i, j+1) + Me%Father%ConnectionY(i+1, j+1))/2.)/2.

            if (Me%Father%WaterPoints2D(i, j) == WaterPoint .and.                   &
                AuxX > Me%InterpolWindow%Xmin .and. AuxX < Me%InterpolWindow%Xmax .and.&
                AuxY > Me%InterpolWindow%Ymin .and. AuxY < Me%InterpolWindow%Ymax) then

                Count           = Count + 1

                NodeX(Count) = AuxX
        
                NodeY(Count) = AuxY

                NodeZ(Count) = FatherField%Values2D(i, j)
                
            endif

        enddo
        enddo

        NumberOfNodes = Count

        allocate(Me%NodeX(NumberOfNodes))
        allocate(Me%NodeY(NumberOfNodes))
        allocate(Me%NodeZ(NumberOfNodes))

        Me%NodeX(1:NumberOfNodes) = NodeX(1:NumberOfNodes)
        Me%NodeY(1:NumberOfNodes) = NodeY(1:NumberOfNodes)
        Me%NodeZ(1:NumberOfNodes) = NodeZ(1:NumberOfNodes)

        deallocate(NodeX)
        deallocate(NodeY)
        deallocate(NodeZ)

in2:    if (NumberOfNodes >= 3) then

            !Constructs Triangulation
            call ConstructTriangulation (Me%ObjTriangulation,   &
                                         NumberOfNodes,         &
                                         Me%NodeX,              &
                                         Me%NodeY,              &
                                         STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR10'

            call SetHeightValues(Me%ObjTriangulation, Me%NodeZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR20'


            allocate(Me%SonCenterX(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, &
                                   NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))
    
            allocate(Me%SonCenterY(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, &
                                   NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))

            do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
            do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB
                
                !Find Son cell center
                Me%SonCenterX(i,j) = ((NewGrid%ConnectionX(i, j  ) + NewGrid%ConnectionX(i+1, j  ))/2. + &
                                      (NewGrid%ConnectionX(i, j+1) + NewGrid%ConnectionX(i+1, j+1))/2.)/2.
    
                Me%SonCenterY(i,j) = ((NewGrid%ConnectionY(i, j  ) + NewGrid%ConnectionY(i+1, j  ))/2. + &
                                      (NewGrid%ConnectionY(i, j+1) + NewGrid%ConnectionY(i+1, j+1))/2.)/2.

            enddo
            enddo


            write(*,*)'Interpolating : '//trim(FatherField%Name), FatherField%IDNumber

            do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
            do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

            if(NewGrid%WaterPoints2D(i, j) == WaterPoint) then

                    NewField%Values2D(i, j) = InterPolation(Me%ObjTriangulation,            &
                                                            Me%SonCenterX(i,j),             &
                                                            Me%SonCenterY(i,j),             &
                                                            FillOutsidePoints,              &
                                                            Default = null_real,            &
                                                            STAT    = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'Triangulator - ModuleInterpolateGrids - ERR30'

            end if

            enddo
            enddo

            call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Triangulator - ModuleInterpolateGrids - ERR40'

            deallocate(Me%SonCenterX)
            deallocate(Me%SonCenterY)

        endif in2

        deallocate(Me%NodeX)
        deallocate(Me%NodeY)
        deallocate(Me%NodeZ)

        endif iN


    end subroutine Triangulator

        !-----------------------------------------------------------------------

    subroutine AveragePointsInCells (FatherField, NewField, NewGrid, k)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField, NewField
        type(T_Grid )                               :: NewGrid
        integer, optional                           :: k

        !Local-----------------------------------------------------------------
        integer                                     :: k_
        !Begin-----------------------------------------------------------------

        if (Me%StationaryMap%ON) then
            if (present(k)) then
                k_ = k
            else
                k_ = 1
            endif
                 
            call ModifyFillingCells(FatherField, NewField, k_)

        else

            call AverageInCellsNotStationary (FatherField, NewField, NewGrid)

        endif

    end subroutine AveragePointsInCells 

        !-----------------------------------------------------------------------

    subroutine AverageInCellsNotStationary (FatherField, NewField, NewGrid)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField, NewField
        type(T_Grid )                               :: NewGrid

        !Local-----------------------------------------------------------------
        real,    dimension(:),   pointer            :: NodeX, NodeY, NodeZ
        real,    dimension(:),   pointer            :: Xaverage, Yaverage, Zaverage
        logical, dimension(:,:), pointer            :: Cellsaverage
        real                                        :: AuxX, AuxY
        integer                                     :: Naverage
        integer                                     :: NumberOfNodes, Count, i, j
        !integer,    dimension(:,:  ), pointer       :: WaterPoints2D
        
        !Begin-----------------------------------------------------------------


        NumberOfNodes =  Sum(Me%Father%WaterPoints2D(Me%Father%WorkSize2D%ILB:Me%Father%WorkSize2D%IUB, &
                                                     Me%Father%WorkSize2D%JLB:Me%Father%WorkSize2D%JUB))

        allocate(NodeX(NumberOfNodes))
        allocate(NodeY(NumberOfNodes))
        allocate(NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
        do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB

            AuxX = ((Me%Father%ConnectionX(i, j  ) + Me%Father%ConnectionX(i+1, j  ))/2. + &
                               (Me%Father%ConnectionX(i, j+1) + Me%Father%ConnectionX(i+1, j+1))/2.)/2.
    
            AuxY = ((Me%Father%ConnectionY(i, j  ) + Me%Father%ConnectionY(i+1, j  ))/2. + &
                               (Me%Father%ConnectionY(i, j+1) + Me%Father%ConnectionY(i+1, j+1))/2.)/2.

            if (Me%Father%WaterPoints2D(i, j) == WaterPoint .and.                   &
                AuxX > Me%InterpolWindow%Xmin .and. AuxX < Me%InterpolWindow%Xmax .and.&
                AuxY > Me%InterpolWindow%Ymin .and. AuxY < Me%InterpolWindow%Ymax) then

                Count           = Count + 1

                NodeX(Count) = AuxX
        
                NodeY(Count) = AuxY

                NodeZ(Count) = FatherField%Values2D(i, j)
                
            endif

        enddo
        enddo

        NumberOfNodes = Count

        allocate(Me%NodeX(NumberOfNodes))
        allocate(Me%NodeY(NumberOfNodes))
        allocate(Me%NodeZ(NumberOfNodes))

        Me%NodeX(1:NumberOfNodes) = NodeX(1:NumberOfNodes)
        Me%NodeY(1:NumberOfNodes) = NodeY(1:NumberOfNodes)
        Me%NodeZ(1:NumberOfNodes) = NodeZ(1:NumberOfNodes)

        deallocate(NodeX)
        deallocate(NodeY)
        deallocate(NodeZ)

        call FillingCells(Me%NodeX, Me%NodeY, Me%NodeZ, NumberOfNodes, NewGrid, NewField, &
                          Xaverage, Yaverage, Zaverage, Cellsaverage, Naverage)


        deallocate(Me%NodeX)
        deallocate(Me%NodeY)
        deallocate(Me%NodeZ)

        deallocate(Xaverage, Yaverage, Zaverage, Cellsaverage)


    end subroutine AverageInCellsNotStationary


    !--------------------------------------------------------------------------


    subroutine FillingCells(XPoints, YPoints, ZPoints, NPoints, NewGrid, NewField, &
                            Xaverage, Yaverage, Zaverage, Cellsaverage, Naverage)

        !Arguments-------------------------------------------------------------
        real,    dimension(:),   pointer    :: XPoints, YPoints, ZPoints
        real,    dimension(:),   pointer    :: Xaverage, Yaverage, Zaverage
        type(T_Grid )                       :: NewGrid 
        type(T_Field),           pointer    :: NewField
        integer                             :: NPoints, Naverage
        logical, dimension(:,:), pointer    :: Cellsaverage

        !Local-----------------------------------------------------------------
        real,    dimension(:),   pointer    :: AuxX, AuxY, AuxZ
        logical, dimension(:),   pointer    :: AlreadyUsed
        type (T_PointF),   pointer          :: GridPoint
        type(T_Polygon),   pointer          :: Rect
        real                                :: SumOfDepths
        integer                             :: nPointsInside, p, Pmax
        integer                             :: i, j, k
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW, XC, YC

        !Begin-----------------------------------------------------------------

        write(*,*)"Filling cells with known data..."

        allocate(AlreadyUsed(1:NPoints)) 

        allocate(Cellsaverage(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))

        allocate(GridPoint)
        allocate(Rect     )
        Rect%Count = 5
        allocate(Rect%VerticesF(1:Rect%Count))

        AlreadyUsed(:) = .false.

        Cellsaverage(:,:) = .false.

        Pmax = (NewGrid%WorkSize2D%JUB - NewGrid%WorkSize2D%JLB) * (NewGrid%WorkSize2D%IUB - NewGrid%WorkSize2D%ILB)

        allocate(AuxX(1:Pmax), AuxY(1:Pmax), AuxZ(1:Pmax))

        k = 0

        do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
        do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

iw:         if (NewGrid%WaterPoints2D(i, j) == WaterPoint) then

                XSW = NewGrid%ConnectionX(i, j)
                YSW = NewGrid%ConnectionY(i, j)
                XSE = NewGrid%ConnectionX(i, j + 1)
                YSE = NewGrid%ConnectionY(i, j + 1)
                XNE = NewGrid%ConnectionX(i + 1, j + 1)
                YNE = NewGrid%ConnectionY(i + 1, j + 1)
                XNW = NewGrid%ConnectionX(i + 1, j)
                YNW = NewGrid%ConnectionY(i + 1, j)

                XC  = (XSW + XSE + XNE + XNW) / 4.
                YC  = (YSW + YSE + YNE + YNW) / 4.                

i1:             if (XC > Me%InterpolWindow%Xmin .and. XC < Me%InterpolWindow%Xmax .and. &
                    YC > Me%InterpolWindow%Ymin .and. YC < Me%InterpolWindow%Ymax) then

                    Rect%VerticesF(1)%X = XSW
                    Rect%VerticesF(1)%Y = YSW
                    Rect%VerticesF(2)%X = XSE
                    Rect%VerticesF(2)%Y = YSE
                    Rect%VerticesF(3)%X = XNE
                    Rect%VerticesF(3)%Y = YNE
                    Rect%VerticesF(4)%X = XNW
                    Rect%VerticesF(4)%Y = YNW

                    Rect%VerticesF(5)%X = Rect%VerticesF(1)%X
                    Rect%VerticesF(5)%Y = Rect%VerticesF(1)%Y

                    call SetLimits(Rect)

                    SumOfDepths   = 0.
                    nPointsInside = 0

                    do p = 1, NPoints

                        if (.not.AlreadyUsed(p)) then

                            GridPoint%X = XPoints(p)
                            GridPoint%Y = YPoints(p)

                            if(IsPointInsidePolygon(GridPoint, Rect))then
                                SumOfDepths    = SumOfDepths + ZPoints(p)
                                nPointsInside  = nPointsInside + 1
                                AlreadyUsed(p) = .true.
                            end if

                        end if

                    end do

                   if(nPointsInside > 0)then
                        k           = k + 1
                        AuxX(k)     = XC
                        AuxY(k)     = YC
                        AuxZ(k)     = SumOfDepths / nPointsInside
                        NewField%Values2D(i, j) = AuxZ(k)
                        Cellsaverage(i, j) = .true. 
                    end if

                endif i1

            endif iw

        end do
        end do

        allocate(Xaverage(1:k), Yaverage(1:k), Zaverage(1:k))

        Naverage = k
        Xaverage(1:k) = AuxX(1:k)
        Yaverage(1:k) = AuxY(1:k)
        Zaverage(1:k) = AuxZ(1:k)

        deallocate(AuxX, AuxY, AuxZ)

        deallocate(Rect%VerticesF)

        deallocate(GridPoint)
        deallocate(Rect     )

        nullify(GridPoint)
        nullify(Rect     )


    end subroutine FillingCells

    !--------------------------------------------------------------------------


    subroutine ConstructFillingCells

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type(T_Grid ), pointer              :: NewGrid
        real,    dimension(:),   pointer    :: XPoints, YPoints
        integer, dimension(:,:,:), pointer  :: FatherWP3D, NewWP3D
        integer, dimension(:,:), pointer    :: FatherWP2D, NewWP2D
        integer, dimension(:),   pointer    :: MapAux, IPoints, JPoints, IndexAux
        real,    dimension(:),   pointer    :: AuxX, AuxY, Aux2X, Aux2Y
        logical, dimension(:),   pointer    :: AlreadyUsed
        type (T_PointF),   pointer          :: GridPoint
        type(T_Polygon),   pointer          :: Rect
        real                                :: FatherX, FatherY
        integer                             :: NPoints 
        integer                             :: nPointsInside, p, NumberOfSonNodes, NumberOfFatherNodes
        integer                             :: i, j, imap, aa, Nmax, jj, ii, l, KLB, KUB, k, STAT_CALL, ifa
        real                                :: XSW, YSW, XSE, YSE, XNE, YNE, XNW, YNW, XC, YC

        !Begin-----------------------------------------------------------------

        if(Me%Interpolation3D) then     
            NewGrid => Me%Aux
            KLB     =  Me%Father%Worksize3D%KLB
            KUB     =  Me%Father%Worksize3D%KUB

            call GetWaterPoints3D(Me%Father%ObjMap, FatherWP3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR10'

            call GetWaterPoints3D(Me%Aux%ObjMap, NewWP3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR20'

            allocate(FatherWP2D(Me%Father%Size2D%ILB:Me%Father%Size2D%IUB, &
                                Me%Father%Size2D%JLB:Me%Father%Size2D%JUB))

            allocate(NewWP2D   (Me%New%Size2D%ILB:Me%New%Size2D%IUB, &
                                Me%New%Size2D%JLB:Me%New%Size2D%JUB))

        else
            NewGrid => Me%New
            KLB     =  1
            KUB     =  1

            call GetWaterPoints2D(Me%Father%ObjHorizontalMap, FatherWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR30'

            call GetWaterPoints2D(Me%New%ObjHorizontalMap, NewWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR40'

        endif

        allocate(Me%StationaryMap%FillingCells(KLB:KUB))

dk:     do k=KLB, KUB

            if (Me%Interpolation3D) then     
                FatherWP2D(:,:) = FatherWP3D(:,:, k)
                NewWP2D   (:,:) = NewWP3D   (:,:, k)
            endif                

            NumberOfFatherNodes =  Sum(FatherWP2D(Me%Father%WorkSize2D%ILB:Me%Father%WorkSize2D%IUB, &
                                                  Me%Father%WorkSize2D%JLB:Me%Father%WorkSize2D%JUB))

            NumberOfSonNodes    =  Sum(NewWP2D(NewGrid%WorkSize2D%ILB:NewGrid%WorkSize2D%IUB, &
                                               NewGrid%WorkSize2D%JLB:NewGrid%WorkSize2D%JUB))

            allocate(MapAux(1:NumberOfSonNodes*3+NumberOfFatherNodes*2))


            allocate(XPoints(NumberOfFatherNodes))
            allocate(YPoints(NumberOfFatherNodes))
            allocate(IPoints(NumberOfFatherNodes))
            allocate(JPoints(NumberOfFatherNodes))

            allocate(GridPoint)
            allocate(Rect     )
            Rect%Count = 5
            allocate(Rect%VerticesF(1:Rect%Count))

            allocate(IndexAux(NumberOfFatherNodes))

            allocate(AlreadyUsed(NumberOfFatherNodes))

            AlreadyUsed(:) = .false.

            allocate(AuxX(1:NumberOfSonNodes))
            allocate(AuxY(1:NumberOfSonNodes))

            allocate(Aux2X(1:NumberOfFatherNodes))
            allocate(Aux2Y(1:NumberOfFatherNodes))


            NPoints = 0
            imap    = 0

            do jj = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
            do ii = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB

                FatherX = ((Me%Father%ConnectionX(ii, jj  ) + Me%Father%ConnectionX(ii+1, jj  ))/2. + &
                           (Me%Father%ConnectionX(ii, jj+1) + Me%Father%ConnectionX(ii+1, jj+1))/2.)/2.
    
                FatherY = ((Me%Father%ConnectionY(ii, jj  ) + Me%Father%ConnectionY(ii+1, jj  ))/2. + &
                           (Me%Father%ConnectionY(ii, jj+1) + Me%Father%ConnectionY(ii+1, jj+1))/2.)/2.

                if (FatherWP2D(ii, jj) == WaterPoint .and.                                      &
                    FatherX > Me%InterpolWindow%Xmin .and. FatherX < Me%InterpolWindow%Xmax .and.&
                    FatherY > Me%InterpolWindow%Ymin .and. FatherY < Me%InterpolWindow%Ymax) then

                    NPoints      = NPoints + 1

                    XPoints(NPoints) = FatherX
        
                    YPoints(NPoints) = FatherY

                    IPoints(NPoints) = ii
                    JPoints(NPoints) = jj

                endif

            enddo
            enddo

            ifa = 0

            do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
            do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

    iw:         if (NewWP2D(i, j) == WaterPoint) then

                    XSW = NewGrid%ConnectionX(i, j)
                    YSW = NewGrid%ConnectionY(i, j)
                    XSE = NewGrid%ConnectionX(i, j + 1)
                    YSE = NewGrid%ConnectionY(i, j + 1)
                    XNE = NewGrid%ConnectionX(i + 1, j + 1)
                    YNE = NewGrid%ConnectionY(i + 1, j + 1)
                    XNW = NewGrid%ConnectionX(i + 1, j)
                    YNW = NewGrid%ConnectionY(i + 1, j)

                    XC  = (XSW + XSE + XNE + XNW) / 4.
                    YC  = (YSW + YSE + YNE + YNW) / 4.                

    i1:             if (XC > Me%InterpolWindow%Xmin .and. XC < Me%InterpolWindow%Xmax .and. &
                        YC > Me%InterpolWindow%Ymin .and. YC < Me%InterpolWindow%Ymax) then

                        Rect%VerticesF(1)%X = XSW
                        Rect%VerticesF(1)%Y = YSW
                        Rect%VerticesF(2)%X = XSE
                        Rect%VerticesF(2)%Y = YSE
                        Rect%VerticesF(3)%X = XNE
                        Rect%VerticesF(3)%Y = YNE
                        Rect%VerticesF(4)%X = XNW
                        Rect%VerticesF(4)%Y = YNW

                        Rect%VerticesF(5)%X = Rect%VerticesF(1)%X
                        Rect%VerticesF(5)%Y = Rect%VerticesF(1)%Y

                        call SetLimits(Rect)

                        nPointsInside = 0
                        aa            = 0
                        do p = 1, NPoints

                            if (.not.AlreadyUsed(p)) then

                                GridPoint%X = XPoints(p)
                                GridPoint%Y = YPoints(p)

                                if(IsPointInsidePolygon(GridPoint, Rect))then
                                    nPointsInside  = nPointsInside + 1
                                    AlreadyUsed(p) = .true.
                                    aa             = aa + 1
                                    IndexAux(aa)   = p
                                end if

                            end if
                    
                        end do

                        if(nPointsInside > 0)then

                            imap         = imap + 1
                            Mapaux(imap) = i
        
                            imap         = imap + 1
                            Mapaux(imap) = j

                            imap         = imap + 1
                            Mapaux(imap) = nPointsInside

                            do aa = 1, nPointsInside
                                imap         = imap + 1
                                Mapaux(imap) = IPoints(IndexAux(aa))
                                imap         = imap + 1
                                Mapaux(imap) = JPoints(IndexAux(aa))
                                ifa          = ifa + 1
                                Aux2X (ifa)  = XPoints(IndexAux(aa))
                                Aux2Y (ifa)  = YPoints(IndexAux(aa))

                            enddo
                            l           = l + 1
                            AuxX(l)     = XC
                            AuxY(l)     = YC
                        end if

                    endif i1

                endif iw

            end do
            end do


            NumberOfFatherNodes = ifa
            NumberOfSonNodes    = l

            Me%StationaryMap%FillingCells(k)%NFather = NumberOfFatherNodes
            Me%StationaryMap%FillingCells(k)%NSon    = NumberOfSonNodes

            allocate(Me%StationaryMap%FillingCells(k)%Xin(1:NumberOfFatherNodes))
            allocate(Me%StationaryMap%FillingCells(k)%Yin(1:NumberOfFatherNodes))
            allocate(Me%StationaryMap%FillingCells(k)%Zin(1:NumberOfFatherNodes))

            Me%StationaryMap%FillingCells(k)%Xin(1:NumberOfFatherNodes) = Aux2X(1:NumberOfFatherNodes)
            Me%StationaryMap%FillingCells(k)%Yin(1:NumberOfFatherNodes) = Aux2Y(1:NumberOfFatherNodes)

            deallocate(XPoints)
            deallocate(YPoints)
            deallocate(IPoints)
            deallocate(JPoints)
            deallocate(Aux2X  )
            deallocate(Aux2Y  )

            Nmax = NumberOfSonNodes*3+NumberOfFatherNodes*2

            Me%StationaryMap%FillingCells(k)%MapN =Nmax

            if (imap /= Nmax) stop 'ConstructFillingCells - Module Interpolation Grids - ERR50'

            allocate(Me%StationaryMap%FillingCells(k)%Mapping(1:Nmax))

            Me%StationaryMap%FillingCells(k)%Mapping(1:Nmax) = Mapaux(1:Nmax)

            deallocate(Mapaux)

            allocate(Me%StationaryMap%FillingCells(k)%Xout(1:NumberOfSonNodes))
            allocate(Me%StationaryMap%FillingCells(k)%Yout(1:NumberOfSonNodes))
            allocate(Me%StationaryMap%FillingCells(k)%Zout(1:NumberOfSonNodes))

            Me%StationaryMap%FillingCells(k)%Xout(1:NumberOfSonNodes) = AuxX(1:NumberOfSonNodes)
            Me%StationaryMap%FillingCells(k)%Yout(1:NumberOfSonNodes) = AuxY(1:NumberOfSonNodes)

            deallocate(AuxX,  AuxY)

            deallocate(Rect%VerticesF)

            nullify(GridPoint)
            nullify(Rect     )

        enddo dk


        if(Me%Interpolation3D) then     
            call UnGetMap(Me%Father%ObjMap, FatherWP3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR60'

            call UnGetMap(Me%Aux%ObjMap, NewWP3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR70'

            deallocate(FatherWP2D)

            deallocate(NewWP2D)

        else
            call UnGetHorizontalMap(Me%Father%ObjHorizontalMap, FatherWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR80'

            call UnGetHorizontalMap(Me%New%ObjHorizontalMap, NewWP2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ConstructFillingCells - ModuleInterpolateGrids - ERR90'

        endif

        nullify(NewGrid)

    end subroutine ConstructFillingCells

!---------------------------------------------------------------------------------------------

    subroutine ModifyFillingCells(FatherField, NewField, k)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField, NewField
        integer                                     :: k

        !Local-----------------------------------------------------------------
        integer,  dimension(:),  pointer            :: Mapping
        integer                                     :: is, i, j, ii, jj, nPointsInside, imap, ip 

        !Begin-----------------------------------------------------------------

        Mapping => Me%StationaryMap%FillingCells(k)%Mapping

        imap  = 0

        do is = 1, Me%StationaryMap%FillingCells(k)%NSon

            imap = imap + 1
            i    = Mapping(imap)

            imap = imap + 1
            j    = Mapping(imap)

            imap = imap + 1
            nPointsInside = Mapping(imap)

            NewField%Values2D(i, j) = 0.

            do ip = 1, nPointsInside
                imap = imap + 1
                ii   = Mapping(imap)

                imap = imap + 1
                jj   = Mapping(imap)

                NewField%Values2D(i, j) = NewField%Values2D(i, j) +                     &
                                          FatherField%Values2D(ii, jj)/real(nPointsInside)
            enddo

        enddo

        nullify(Mapping)

    end subroutine ModifyFillingCells

!---------------------------------------------------------------------------------------------

    subroutine KillFillingCells
        
        !Arguments----------------------------------------------------------------------------

        !Local------------------------------------------------------------------------------
        integer                     :: KLB, KUB, k

        !Begin--------------------------------------------------------------------------------

        KLB = Me%Father%WorkSize3D%KLB
        KUB = Me%Father%WorkSize3D%KUB

        do k=KLB, KUB

            deallocate(Me%StationaryMap%FillingCells(k)%Xin)
            deallocate(Me%StationaryMap%FillingCells(k)%Yin)
            deallocate(Me%StationaryMap%FillingCells(k)%Zin)

            nullify   (Me%StationaryMap%FillingCells(k)%Xin)
            nullify   (Me%StationaryMap%FillingCells(k)%Yin)
            nullify   (Me%StationaryMap%FillingCells(k)%Zin)

            deallocate(Me%StationaryMap%FillingCells(k)%Mapping)

            nullify   (Me%StationaryMap%FillingCells(k)%Mapping)

            deallocate(Me%StationaryMap%FillingCells(k)%Xout)
            deallocate(Me%StationaryMap%FillingCells(k)%Yout)
            deallocate(Me%StationaryMap%FillingCells(k)%Zout)

            nullify   (Me%StationaryMap%FillingCells(k)%Xout)
            nullify   (Me%StationaryMap%FillingCells(k)%Yout)
            nullify   (Me%StationaryMap%FillingCells(k)%Zout)

        enddo

        deallocate(Me%StationaryMap%FillingCells)

        nullify   (Me%StationaryMap%FillingCells)


    end subroutine KillFillingCells

!-----------------------------------------------------------------------------------------

    subroutine ExtraPol2DFieldsTriang (FatherField, NewGrid, FillOutsidePoints, &
                                       InValues2D, OutValues2D)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField
        type(T_Grid )                               :: NewGrid
        logical                                     :: FillOutsidePoints
        real,         dimension(:,:), pointer       :: InValues2D, OutValues2D

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, Count, i, j
        real,       dimension(:,:),   pointer       :: LatitudeZ, LongitudeZ
        
        !Begin-----------------------------------------------------------------

        NumberOfNodes =  Sum(Me%Father%WaterPoints2D(Me%Father%WorkSize2D%ILB:Me%Father%WorkSize2D%IUB, &
                                                     Me%Father%WorkSize2D%JLB:Me%Father%WorkSize2D%JUB))

iN:     if (NumberOfNodes >= 3) then

        allocate(Me%NodeX(NumberOfNodes))
        allocate(Me%NodeY(NumberOfNodes))
        allocate(Me%NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
        do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB

            if (Me%Father%WaterPoints2D(i, j) == WaterPoint .and. InValues2D(i, j) > FillValueReal/4.) then

                Count           = Count + 1

                Me%NodeX(Count) = ((Me%Father%ConnectionX(i, j  ) + Me%Father%ConnectionX(i+1, j  ))/2. + &
                                   (Me%Father%ConnectionX(i, j+1) + Me%Father%ConnectionX(i+1, j+1))/2.)/2.
        
                Me%NodeY(Count) = ((Me%Father%ConnectionY(i, j  ) + Me%Father%ConnectionY(i+1, j  ))/2. + &
                                   (Me%Father%ConnectionY(i, j+1) + Me%Father%ConnectionY(i+1, j+1))/2.)/2.

                Me%NodeZ(Count) = InValues2D(i, j)
                

            endif

        enddo
        enddo

        !Actualize number of nodes
        NumberOfNodes = Count

            !Constructs Triangulation
        call ConstructTriangulation (Me%ObjTriangulation,   &
                                     NumberOfNodes,         &
                                     Me%NodeX,              &
                                     Me%NodeY,              &
                                     STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR10'

        call SetHeightValues(Me%ObjTriangulation, Me%NodeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR20'


        allocate(Me%SonCenterX(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, &
                               NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))
    
        allocate(Me%SonCenterY(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, &
                               NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))

        call GetGridLatitudeLongitude(NewGrid%ObjHorizontalGrid,                        &
                                      GridLongitude = LongitudeZ,                       &
                                      GridLatitude  = LatitudeZ,                        &
                                      STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR30'

        do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
        do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

            !Find Son cell center
            Me%SonCenterX(i,j)   = LongitudeZ(i,j)
            Me%SonCenterY(i,j)   = LatitudeZ(i,j)

        enddo
        enddo

        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LongitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR40'
            
        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LatitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR50'

        write(*,*)'Extrapolating : '//trim(FatherField%Name), FatherField%IDNumber

        do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
        do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

            if( NewGrid%WaterPoints2D(i, j) /= WaterPoint .or. &
               (NewGrid%WaterPoints2D(i, j) == WaterPoint .and. OutValues2D(i, j)  < FillValueReal/4)) then

                OutValues2D(i, j) = InterPolation(Me%ObjTriangulation,                  &
                                                        Me%SonCenterX(i,j),             &
                                                        Me%SonCenterY(i,j),             &
                                                        FillOutsidePoints,              &
                                                        Default = null_real,            &
                                                        STAT    = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR60'

            end if

        enddo
        enddo

        call KillTriangulation (Me%ObjTriangulation, STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ExtraPol2DFieldsTriang - ModuleInterpolateGrids - ERR70'

        deallocate(Me%SonCenterX)
        deallocate(Me%SonCenterY)

        deallocate(Me%NodeX)
        deallocate(Me%NodeY)
        deallocate(Me%NodeZ)

        endif iN


    end subroutine ExtraPol2DFieldsTriang



    subroutine ExtraPol2DFieldsNearest (FatherField, NewGrid, InValues2D, OutValues2D)

        !Arguments-------------------------------------------------------------
        type(T_Field),           pointer            :: FatherField
        type(T_Grid )                               :: NewGrid
        real,         dimension(:,:), pointer       :: InValues2D, OutValues2D

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: NumberOfNodes, Count, i, j
        real                                        :: Distance, Aux
        real,       dimension(:,:),   pointer       :: LatitudeZ, LongitudeZ
        
        !Begin-----------------------------------------------------------------


        NumberOfNodes =  Sum(Me%Father%WaterPoints2D(Me%Father%WorkSize2D%ILB:Me%Father%WorkSize2D%IUB, &
                                                     Me%Father%WorkSize2D%JLB:Me%Father%WorkSize2D%JUB))


iN:     if (NumberOfNodes > 0) then

        allocate(Me%NodeX(NumberOfNodes))
        allocate(Me%NodeY(NumberOfNodes))
        allocate(Me%NodeZ(NumberOfNodes))

        Count = 0

        do j = Me%Father%WorkSize2D%JLB, Me%Father%WorkSize2D%JUB
        do i = Me%Father%WorkSize2D%ILB, Me%Father%WorkSize2D%IUB

            if (Me%Father%WaterPoints2D(i, j) == WaterPoint .and. InValues2D(i, j) > Me%ExtrapolateLimit) then

                Count           = Count + 1

                Me%NodeX(Count) = ((Me%Father%ConnectionX(i, j  ) + Me%Father%ConnectionX(i+1, j  ))/2. + &
                                   (Me%Father%ConnectionX(i, j+1) + Me%Father%ConnectionX(i+1, j+1))/2.)/2.
        
                Me%NodeY(Count) = ((Me%Father%ConnectionY(i, j  ) + Me%Father%ConnectionY(i+1, j  ))/2. + &
                                   (Me%Father%ConnectionY(i, j+1) + Me%Father%ConnectionY(i+1, j+1))/2.)/2.

                Me%NodeZ(Count) = InValues2D(i, j)
                

            endif

        enddo
        enddo

        !Actualize number of nodes
        NumberOfNodes = Count


        allocate(Me%SonCenterX(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, &
                               NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))
    
        allocate(Me%SonCenterY(NewGrid%Size2D%ILB:NewGrid%Size2D%IUB, &
                               NewGrid%Size2D%JLB:NewGrid%Size2D%JUB))

        call GetGridLatitudeLongitude(NewGrid%ObjHorizontalGrid,                        &
                                      GridLongitude = LongitudeZ,                       &
                                      GridLatitude  = LatitudeZ,                        &
                                      STAT  = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsNearest - ModuleInterpolateGrids - ERR30'

        do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
        do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

                !Find Son cell center
                Me%SonCenterX(i,j)   = LongitudeZ(i,j)
                Me%SonCenterY(i,j)   = LatitudeZ(i,j)

        enddo
        enddo

        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LongitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsNearest - ModuleInterpolateGrids - ERR40'
            
        call UngetHorizontalGrid(NewGrid%ObjHorizontalGrid, LatitudeZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExtraPol2DFieldsNearest - ModuleInterpolateGrids - ERR50'


        write(*,*)'Extrapolating : '//trim(FatherField%Name), FatherField%IDNumber


        do j = NewGrid%WorkSize2D%JLB, NewGrid%WorkSize2D%JUB
        do i = NewGrid%WorkSize2D%ILB, NewGrid%WorkSize2D%IUB

        if ( NewGrid%WaterPoints2D(i, j) /= WaterPoint  .or. &
            (NewGrid%WaterPoints2D(i, j) == WaterPoint .and. OutValues2D(i, j) < Me%ExtrapolateLimit)) then

            Aux = - FillValueReal

            do Count = 1, NumberOfNodes

                Distance = sqrt((Me%NodeX(Count) - Me%SonCenterX(i,j))**2 +                &
                                (Me%NodeY(Count) - Me%SonCenterY(i,j))**2)

                !Distance = abs(Me%NodeX(Count) - Me%SonCenterX(i,j))       +            &
                !           abs(Me%NodeY(Count) - Me%SonCenterY(i,j))


                if (Distance < Aux .and. Me%NodeZ(Count) /= FillValueReal) then

                    OutValues2D(i, j) = Me%NodeZ(Count)

                    Aux = Distance

                endif

            enddo

        end if

        enddo
        enddo

        deallocate(Me%SonCenterX)
        deallocate(Me%SonCenterY)

        deallocate(Me%NodeX)
        deallocate(Me%NodeY)
        deallocate(Me%NodeZ)

        endif iN


    end subroutine ExtraPol2DFieldsNearest



    subroutine ExtraPol2DFieldsNearestCell (NewGrid, OutValues2D)

        !Arguments-------------------------------------------------------------
        type(T_Grid )                               :: NewGrid
        real,         dimension(:,:), pointer       :: OutValues2D

        !Local-----------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB
        integer                                     :: dij, Count, i, j, NumberOfCells
        integer                                     :: jj, ii, dijmax, dimax, djmax
        real                                        :: SumValues
        
        !Begin-----------------------------------------------------------------

        JLB = NewGrid%WorkSize2D%JLB; JUB = NewGrid%WorkSize2D%JUB;
        ILB = NewGrid%WorkSize2D%ILB; IUB = NewGrid%WorkSize2D%IUB;


        NumberOfCells =  Sum(NewGrid%WaterPoints2D(ILB:IUB, JLB:JUB))

        if (NumberOfCells > 0) then

            do j = JLB, JUB
            do i = ILB, IUB

            if (OutValues2D(i, j) < Me%ExtrapolateLimit) then

                    dimax = IUB-ILB + 1
                    djmax = JUB-JLB + 1

                    dijmax = max(dimax, djmax)
                
                    SumValues   = 0
                    Count = 0

                    do dij=1,dijmax

                        do jj=j-dij,j+dij
                        do ii=i-dij,i+dij

                            if (jj < JLB) cycle
                            if (jj > JUB) cycle
                            if (ii < ILB) cycle
                            if (ii > IUB) cycle

                            if (OutValues2D(ii, jj) > Me%ExtrapolateLimit) then
                                SumValues   = SumValues   + OutValues2D(ii, jj) 
                                Count = Count + 1
                            endif

                        enddo
                        enddo

                        if (Count > 0) exit

                    enddo

                    if (Count > 0) then

                        OutValues2D(i, j) = SumValues / real(Count)

                    else
                        stop 'ExtraPol2DFieldsNearestCell - ModuleInterpolateGrids - ERR10'
                    endif

                endif

            enddo
            enddo

        endif

    end subroutine ExtraPol2DFieldsNearestCell


    subroutine ExtraPol2DFieldsConstantValue (NewGrid, OutValues2D)

        !Arguments-------------------------------------------------------------
        type(T_Grid )                               :: NewGrid
        real,         dimension(:,:), pointer       :: OutValues2D

        !Local-----------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB
        integer                                     :: i, j, NumberOfCells
        
        !Begin-----------------------------------------------------------------

        JLB = NewGrid%WorkSize2D%JLB; JUB = NewGrid%WorkSize2D%JUB;
        ILB = NewGrid%WorkSize2D%ILB; IUB = NewGrid%WorkSize2D%IUB;


        NumberOfCells =  Sum(NewGrid%WaterPoints2D(ILB:IUB, JLB:JUB))

        if (NumberOfCells > 0) then

            do j = JLB, JUB
            do i = ILB, IUB

                if (OutValues2D(i, j) < FillValueReal/4.) OutValues2D(i, j) = Me%ExtrapolateValue

            enddo
            enddo

        endif

    end subroutine ExtraPol2DFieldsConstantValue



end module ModuleInterpolateGrids
