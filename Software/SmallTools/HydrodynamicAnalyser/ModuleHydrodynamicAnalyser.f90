!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : HydrodynamicAnalyser
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as HydrodynamicAnalyser to create new modules
!
!------------------------------------------------------------------------------


Module ModuleHydrodynamicAnalyser

    use ModuleGlobalData
    use ModuleFunctions
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalMap
    use ModuleGeometry
    use ModuleMap
    use ModuleStopWatch,            only : StartWatch, StopWatch         
    use ModuleFillMatrix,           only: ConstructFillMatrix, GetDefaultValue, KillFillMatrix, &
                                          ModifyFillMatrix


    implicit none

    private 

    !Parameter
    real,    parameter                              :: StefanBoltzmann          = 5.669e-08     ![W/m2/K4]
    integer, parameter                              :: SimpleHeight_            = 1
    integer, parameter                              :: ComplexHeight_           = 2

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartHydrodynamicAnalyser
    private ::      ReadOptions
    private ::          ReadInitialDensityField
    private ::      ConstructGrid
    private ::      AllocateMatrixes
    private ::      Open_HDF5_Input_File
    private ::          Open_HDF5_Hydro_File
    private ::          Open_HDF5_Water_File
    private ::      Open_HDF5_OutPut_File

                     
    
    !Modifier
    public  :: ModifyHydrodynamicAnalyser
    private ::      FlowProperties
    private ::      ResidualFlowProperties
    private ::          ComputeVorticity
    private ::          ComputeKineticEnergy
    private ::          ComputePerturbationPE
    private ::          ComputeBaroclinicForce
    private ::          ComputeSystemEnergy
    private ::          ComputeFlowAlongSections
    private ::      Read_HDF5_Residual_Hydro
    private ::      Read_HDF5_Hydro_File
    private ::      Read_HDF5_Water_File
    private ::      Read_HDF5_CommonData
    private ::      Write_HDF5_OutPut_File
    private ::      Write_HDF5_Residual_Hydro
    private ::      WriteEnergyDataFile


    !Destructor
    public  :: KillHydrodynamicAnalyser                                                     
    private ::      DeAllocateMatrixes


    !Paramenters---------------------------------------------------------------------
    !Input / Output
    integer, parameter :: FileOpen = 1, FileClose = 0
    integer, parameter :: EnergyBufferSize        = 1000

    real,    parameter :: Abbott          = 1.
    real,    parameter :: Leendertse      = 2.

    integer, parameter :: CurrentVelocity    = 0
    integer, parameter :: ResidualVelocity   = 1
    integer, parameter :: ResidualFlux       = 2
    integer, parameter :: ResidualFluxVel    = 3
    integer, parameter :: BaroclinicVelocity = 4

    character(len = StringLength), parameter    :: block_begin = '<beginsection>'
    character(len = StringLength), parameter    :: block_end   = '<endsection>'

    character(len = StringLength), parameter    :: beginvert   = '<<beginvertical>>'
    character(len = StringLength), parameter    :: endvert     = '<<endvertical>>'

    character(len = StringLength), parameter    :: begincells  = '<<begincells>>'
    character(len = StringLength), parameter    :: endcells    = '<<endcells>>'

    !Types---------------------------------------------------------------------

    private :: T_Coef_Baroc
    type T_Coef_Baroc

        Integer, dimension( : ), pointer   :: Kleft,Kright

        Real(8), dimension( : ), pointer   :: Depth_integ, Hcenter, Hleft,&
                                              Hright, HroLeft, HroRight, DensRight, DensLeft
    end type T_Coef_Baroc

    private :: T_InputVar
    type       T_InputVar
        real,    dimension(:, :, :),  pointer       :: Density, Salinity, Temperature, SigmaDens
        real,    dimension(:, :, :),  pointer       :: VelocityU, VelocityV, VelocityW
        real,    dimension(:, :, :),  pointer       :: ResidualVelU, ResidualVelV, ResidualFluxU, ResidualFluxV

        real,    dimension(:, :   ),  pointer       :: WaterLevel, Bathymetry, ResidualWaterLevel
        real,    dimension(:, :   ),  pointer       :: Latitude, Longitude
        real,    dimension(:, :   ),  pointer       :: Connection_X, Connection_Y
        real,    dimension(:, :   ),  pointer       :: DZX, DZY, DUX, DVY
        real,    dimension(:, :, :),  pointer       :: SZZ, DWZ, DUZ, DVZ, ZCellCenter
        real(8), dimension(:, :, :),  pointer       :: Volume_Z, Volume_U, Volume_V, Volume_W
        real,    dimension(:, :   ),  pointer       :: WaterColumn, WaterColumnU, WaterColumnV
        integer, dimension(:, :, :),  pointer       :: WaterPoints3D
        integer, dimension(:, :, :),  pointer       :: OpenPoints3D
        integer, dimension(:, :   ),  pointer       :: WaterPoints2D
        integer, dimension(:, :   ),  pointer       :: KFloor_U, KFloor_V, KFloor_Z

        logical                                     :: IniDensON, IniTempON, IniSalON, IniWaterLevelON
        real,    dimension(:, :, :),  pointer       :: IniDensity, IniSalinity, IniTemperature, IniSZZ, IniSigmaDens
        real,    dimension(:, :   ),  pointer       :: IniWaterLevel

    end type  T_InputVar

    private :: T_OutputVar
    type       T_OutputVar
        real,    dimension(:, :, :),  pointer       :: BaroclinicForceX, BaroclinicForceY, NN
        real,    dimension(:, :, :),  pointer       :: Vorticity
        real,    dimension(:, :, :),  pointer       :: BaroclinicU, BaroclinicV
        real,    dimension(:, :   ),  pointer       :: BaroclinicKE, KineticEnergy, PerturbationPE, &
                                                       PerturbationPE_V2, BarotropicU, BarotropicV
        real,    dimension(:, :, :),  pointer       :: ResVorticity
        real,    dimension(:, :, :),  pointer       :: ResBaroclinicU, ResBaroclinicV
        real,    dimension(:, :   ),  pointer       :: ResBaroclinicKE, ResKineticEnergy,&
                                                       ResBarotropicU,  ResBarotropicV
    end type  T_OutputVar

    private :: T_Transport
    type       T_Transport
        real,    dimension(:),        pointer       :: DepthMin, DepthMax, ZonalFlow, MeridionalFlow
        integer, dimension(:),        pointer       :: Icell, Jcell 
        real,    dimension(:,:),      pointer       :: CellZonalFlow, CellMeridionalFlow
        character(len = StringLength)               :: Name
        integer                                     :: CellsNumber, VertNumber
        type (T_Transport), pointer                 :: Next, Prev
                                                   
    end type   T_Transport


    type T_Energy
        integer                         :: FileID
        integer                         :: BufferCount
        real(8)                         :: PotentialEnergyReference
        logical                         :: FirstTime = .true.
        real, dimension(:), pointer     :: YearBuffer
        real, dimension(:), pointer     :: MonthBuffer
        real, dimension(:), pointer     :: DayBuffer
        real, dimension(:), pointer     :: HourBuffer
        real, dimension(:), pointer     :: MinuteBuffer
        real, dimension(:), pointer     :: SecondBuffer
        real, dimension(:), pointer     :: RelativeKEBuffer
        real, dimension(:), pointer     :: RelativePEBuffer
        real, dimension(:), pointer     :: KineticBuffer
        real, dimension(:), pointer     :: PotentialBuffer
        real, dimension(:), pointer     :: VorticityBuffer
        real, dimension(:), pointer     :: MassBuffer
        real, dimension(:), pointer     :: VolumeBuffer
        real, dimension(:), pointer     :: OpenVolumeBuffer
        real, dimension(:), pointer     :: WaterLevelBuffer
        real, dimension(:), pointer     :: BarotropicKEBuffer
        real, dimension(:), pointer     :: BaroclinicKEBuffer
        real, dimension(:), pointer     :: VelMaxBaroclinicBuffer
        real, dimension(:), pointer     :: VelMaxBuffer

        type (T_Size3D)                 :: Window
        type (T_Time)                   :: NextOutPut
        real                            :: DtOut
    end type T_Energy

    private :: T_Options
    type       T_Options
        logical                                     :: Vorticity, BaroclinicForce, FlowSection
        logical                                     :: KineticEnergy, PerturbationPE, PerturbationPE_V2
        logical                                     :: SquaredBruntVaisallaFrequency
        logical                                     :: HydroON, WaterON, InitialON, ResidualHydro
        logical                                     :: DensPressureCorrect, SameLayerPressure
        logical                                     :: Energy
        integer                                     :: DensityMethod, BaroclinicMethod, &
                                                       BaroclinicPoliDegree
        real                                        :: DZ
    end type  T_Options

    private :: T_FileName
    type       T_FileName
        character(LEN=StringLength)                 :: Geometry, Grid, InputHydroHDF5, InputWaterHDF5, OutPutHDF5, Energy, FlowSection
    end type  T_FileName

    
    private :: T_HydrodynamicAnalyser
    type       T_HydrodynamicAnalyser
        integer                                     :: ObjEnterData         = 0
        integer                                     :: ObjHorizontalGrid    = 0
        integer                                     :: ObjTime              = 0
        integer                                     :: ObjHorizontalMap     = 0
        integer                                     :: ObjBathymetry        = 0
        integer                                     :: ObjGeometry          = 0
        integer                                     :: ObjMap               = 0
        integer                                     :: ObjHDF5_InputHydro   = 0
        integer                                     :: ObjHDF5_InputWater   = 0
        integer                                     :: ObjHDF5_Output       = 0
        integer                                     :: WaterInstants        = FillValueInt
        integer                                     :: HydroInstants        = FillValueInt
        integer                                     :: OutInstants          = FillValueInt
        logical                                     :: DensityExist
        real                                        :: DT
        type (T_Time)                               :: BeginTime, EndTime, CurrentTime
        type (T_Coef_Baroc)                         :: Coef
        type (T_Size3D)                             :: Size, WorkSize
        type (T_InputVar )                          :: InputVar
        type (T_OutputVar)                          :: OutputVar
        type (T_Energy)                             :: Energy
        type (T_Options  )                          :: ComputeOptions
        type (T_FileName )                          :: FileName
        type (T_Transport),           pointer       :: FirstSection, LastSection
        type(T_HydrodynamicAnalyser), pointer       :: Next
    end type  T_HydrodynamicAnalyser

    !Global Module Variables
    type (T_HydrodynamicAnalyser), pointer          :: Me

    logical :: FirstDensityField = .true. 

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartHydrodynamicAnalyser

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: STAT_CALL

        !------------------------------------------------------------------------

        nullify (Me)
        allocate(Me)

        call ConstructEnterData (Me%ObjEnterData, FileName ="HydrodynamicAnalyser.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR10'

        call ReadOptions

        call ConstructGrid

        call AllocateMatrixes


        if (Me%ComputeOptions%PerturbationPE .or. Me%ComputeOptions%InitialON) then

            call ReadInitialDensityField

        endif

        if (Me%ComputeOptions%InitialON) then
        
            Me%OutInstants = 1

        else

            call Open_HDF5_Input_File

        endif

        call Open_HDF5_OutPut_File

        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'StartHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR40'


        call UnGetMap          (Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'StartHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR50'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'StartHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR60'


        !----------------------------------------------------------------------

    end subroutine StartHydrodynamicAnalyser
 
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        logical                                     :: VariableDT

        !Begin-----------------------------------------------------------------
        

        write(*,*)'Reading instructions...'


        call GetData(Me%FileName%Grid,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GRID_FILENAME',                                    &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR20'

        if (iflag == 0)then
            write(*,*)'Must specify name of the grid file'
            stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR30'
        end if


        call GetData(Me%FileName%Geometry,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GEOMETRY_FILENAME',                                &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR40'

        if (iflag == 0)then
            write(*,*)'Must specify name of the geometry file'
            stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR50'
        end if


        Me%ComputeOptions%HydroON = .false.
        Me%ComputeOptions%WaterON = .false.

        call GetData(Me%ComputeOptions%Vorticity,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'VORTICITY',                                        &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR60'

        if (Me%ComputeOptions%Vorticity) Me%ComputeOptions%HydroON = .true.

        call GetData(Me%ComputeOptions%BaroclinicForce,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BAROCLINIC_FORCE',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR70'

        if (Me%ComputeOptions%BaroclinicForce) Me%ComputeOptions%WaterON = .true.

        call GetData(Me%ComputeOptions%KineticEnergy,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'KINETIC_ENERGY',                                   &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR80'

        if (Me%ComputeOptions%KineticEnergy) Me%ComputeOptions%HydroON = .true.

        call GetData(Me%ComputeOptions%PerturbationPE,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POTENTIAL_ENERGY',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR100'

        if (Me%ComputeOptions%PerturbationPE) Me%ComputeOptions%WaterON = .true.

        call GetData(Me%ComputeOptions%PerturbationPE_V2,                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POTENTIAL_ENERGY_V2',                              &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR102'

if2:    if (Me%ComputeOptions%PerturbationPE_V2) then
        
            Me%ComputeOptions%WaterON = .true.

            call GetData(Me%ComputeOptions%DZ,                                          &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'POTENTIAL_ENERGY_DZ',                          &
                         default      = 10.,                                            &
                         ClientModule = 'ModuleHydrodynamicAnalyser',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR103'

        endif if2

        call GetData(Me%ComputeOptions%FlowSection,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'FLOW_SECTION',                                     &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR104'


        if (Me%ComputeOptions%FlowSection) then

            call GetData(Me%FileName%FlowSection,                                       &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'FLOW_SECTION_FILE_OUT',                        &
                         ClientModule = 'ModuleHydrodynamicAnalyser',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR105'
            if (iflag /= 1)            stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR105a'

            call GetSectionsFromFile

        endif
        
        
        call GetData(Me%ComputeOptions%InitialON,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'INITIAL_ANALYSIS',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR106'

        if (Me%ComputeOptions%InitialON .and. Me%ComputeOptions%HydroON)                &
            stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR107'
       


        if (Me%ComputeOptions%HydroON) then

            call GetData(Me%FileName%InputHydroHDF5,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'HYDRODYNAMIC_FILENAME',                        &
                         ClientModule = 'ModuleHydrodynamicAnalyser',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR110'

            if (iflag == 0)then
                write(*,*)'You do not specify name of the hydrodynamic input HDF5  file'
                stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR120'
            end if

        endif

        if (Me%ComputeOptions%WaterON .and. .not. Me%ComputeOptions%InitialON) then

            call GetData(Me%FileName%InputWaterHDF5,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'WATER_FILENAME',                               &
                         ClientModule = 'ModuleHydrodynamicAnalyser',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR130'

            if (iflag == 0)then
                write(*,*)'You do not specify name of the water properties input HDF5  file'
                stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR140'
            end if

        endif


        if (Me%ComputeOptions%WaterON .or. Me%ComputeOptions%HydroON) then

            call GetData(Me%FileName%OutPutHDF5,                                        &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'OUTPUT_FILENAME',                              &
                         ClientModule = 'ModuleHydrodynamicAnalyser',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR150'

            if (iflag == 0)then
                write(*,*)'You do not specify name of the output HDF5  file'
                stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR160'
            end if

        else

            write(*,*)'You do not specify any action to the hydrodynamic analyser tool'
            stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR170'

        endif


        call GetData(Me%ComputeOptions%DensPressureCorrect,                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PRESSURE_CORRECTION',                              &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR180'

        call GetData(Me%ComputeOptions%SameLayerPressure,                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SAME_LAYER_PRESSURE',                              &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR185'


        call GetData(Me%ComputeOptions%DensityMethod,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DENSITY_METHOD',                                   &
                     default      = UNESCOState_,                                       &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR190'

        call GetData(Me%ComputeOptions%BaroclinicMethod,                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BAROCLINIC_METHOD',                                &
                     default      = Leibniz2,                                           &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR200'

        call GetData(Me%ComputeOptions%BaroclinicPoliDegree,                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BAROCLINIC_POLIDEGREE',                            &
                     default      = 3,                                                  &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR210'

        call GetData(Me%ComputeOptions%Energy,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'ENERGY',                                        &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleHydrodynamicAnalyser',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR660'

        call null_time(Me%BeginTime  )
        call null_time(Me%EndTime    )
        call null_time(Me%CurrentTime)

        call ReadTimeKeyWords(Me%ObjEnterData, FromFile, Me%BeginTime, Me%EndTime,      &
                              Me%DT, VariableDT , "ModuleHydrodynamicAnalyser")

        Me%CurrentTime = Me%BeginTime

        if (Me%ComputeOptions%HydroON) then
            
            call GetData(Me%ComputeOptions%ResidualHydro,                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'RESIDUAL_HYDRO',                               &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleHydrodynamicAnalyser',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR230'
        endif

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%EndTime, Me%DT, .false., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleHydrodynamicAnalyser - ERR240'


    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    subroutine ReadInitialDensityField

        !Local-----------------------------------------------------------------
        real,   dimension(:,:,:), pointer           :: Field3D, SZZ
        type (T_PropertyID)                         :: IDProperty
        character(LEN = StringLength), parameter    :: prop_block_begin    = '<beginproperty>'
        character(LEN = StringLength), parameter    :: prop_block_end      = '<endproperty>'
        integer                                     :: ClientNumber, STAT_CALL
        logical                                     :: BlockFound


        !Begin-----------------------------------------------------------------

        allocate(Field3D                                                                &
                 (Me%Size%ILB:Me%Size%IUB,                                              &
                  Me%Size%JLB:Me%Size%JUB,                                              &
                  Me%Size%KLB:Me%Size%KUB))

        Field3D(:,:,:) = FillValueReal
       
        Me%InputVar%IniWaterLevelON  = .false.        
        Me%InputVar%IniDensON        = .false.
        Me%InputVar%IniTempON        = .false.
        Me%InputVar%IniSalON         = .false.


        do

            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = prop_block_begin,             &
                                        block_end       = prop_block_end,               &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
                
EB  :       if      (STAT_CALL .EQ. SUCCESS_ ) then    

BF  :           if (BlockFound) then                                                  
                    

                    call ConstructPropertyID (IDProperty, Me%ObjEnterData, FromBlock)

                    if (IDProperty%IDNumber == WaterLevel_) then

                        call ReadProperty2D(IDProperty, Me%InputVar%IniWaterLevel)

                        Me%InputVar%IniWaterLevelON = .true.

                    endif
    
                else BF


                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR20'

                    exit

                endif BF

            else  EB

                stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR30'

            endif EB

        enddo
        
        if (Me%InputVar%IniWaterLevelON) then

            call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,              &
                                        WaterPoints3D    = Me%InputVar%WaterPoints3D,   &
                                        SurfaceElevation = Me%InputVar%IniWaterLevel,   &
                                        ActualTime       = Me%CurrentTime,              &   
                                        STAT             = STAT_CALL )

            call UnGetMap(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR40'


            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR45'

            !Update the open points
            call UpdateComputeFaces3D(Me%ObjMap, Me%InputVar%IniWaterLevel,             &
                                      Me%CurrentTime, STAT = STAT_CALL)      
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR50'
            
            call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR60'

            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR65'

            call GetGeometryDistances (Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR67'

            Me%InputVar%IniSZZ(:,:,:) = SZZ(:,:,:)

            call UnGetGeometry(Me%ObjGeometry, SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR69'

        else
            
            stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR70'

        endif

        do

            call ExtractBlockFromBuffer(Me%ObjEnterData,                                &
                                        ClientNumber    = ClientNumber,                 &
                                        block_begin     = prop_block_begin,             &
                                        block_end       = prop_block_end,               &
                                        BlockFound      = BlockFound,                   &
                                        STAT            = STAT_CALL)
                
i1  :       if      (STAT_CALL .EQ. SUCCESS_ ) then    
i2  :           if (BlockFound) then                                                  
                    
                    call ConstructPropertyID (IDProperty, Me%ObjEnterData, FromBlock)

                    call ReadProperty3D(IDProperty, Field3D)

                    if (IDProperty%IDNumber == Density_) then

                        Me%InputVar%IniDensity  (:,:,:) = Field3D(:,:,:)

                        Me%InputVar%IniSigmaDens(:,:,:) = Me%InputVar%IniDensity(:,:,:) - SigmaDensityReference

                        Me%InputVar%IniDensON = .true.

                    endif

                    if (IDProperty%IDNumber == Salinity_) then

                        Me%InputVar%IniSalinity (:,:,:) = Field3D(:,:,:)

                        Me%InputVar%IniSalON = .true.

                    endif

                    if (IDProperty%IDNumber == Temperature_) then

                        Me%InputVar%IniTemperature (:,:,:) = Field3D(:,:,:)

                        Me%InputVar%IniTempON = .true.

                    endif



                else i2


                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR80'

                    exit

                endif i2

            else  i1

                stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR90'

            endif i1

        enddo
        
        deallocate(Field3D)

        if (.not. Me%InputVar%IniDensON)  then

            if (.not. Me%InputVar%IniSalON)                                             &
                stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR100'

            if (.not. Me%InputVar%IniTempON)                                            &
                stop 'ReadInitialDensityField - ModuleHydrodynamicAnalyser - ERR110'

            call ModifyDensity(Me%InputVar%IniSalinity,Me%InputVar%IniTemperature,      &
                               Me%InputVar%IniDensity, Me%InputVar%IniSigmaDens)

        endif


    end subroutine ReadInitialDensityField

    !--------------------------------------------------------------------------

    subroutine ModifyDensity(S,T,D, Sigma)

        !Arguments-------------------------------------------------------------
        real,    pointer, dimension(:,:,:)      :: S,T,D, Sigma
        !Local----------------------------------------------------------------- 
        real,    pointer, dimension(:,:,:)      :: SZZ, ZCellCenter
        real,    pointer, dimension(:    )      :: AverageLayerDepth
        integer                                 :: STAT_CALL
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k, Isum
        real                                    :: RoRef, AlphaS, S0, Depth
        
        !Begin----------------------------------------------------------------- 


        ILB = Me%WorkSize%ILB 
        IUB = Me%WorkSize%IUB 
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB 
        KLB = Me%WorkSize%KLB 
        KUB = Me%WorkSize%KUB 

        call GetGeometryDistances (Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyDensity - ModuleWaterProperties - ERR10'

        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  ZCellCenter   = ZCellCenter,                          &
                                  STAT          = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyDensity - ModuleWaterProperties - ERR20'

        if (Me%ComputeOptions%SameLayerPressure .and. Me%ComputeOptions%DensPressureCorrect) then

            allocate(AverageLayerDepth(KLB:KUB))
            do k = KLB, KUB
            AverageLayerDepth(k) = 0.
            Isum                 = 0
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                    AverageLayerDepth(k) = AverageLayerDepth(k) - ZCellCenter(i, j, k)
                    Isum                 = Isum + 1

                endif

            enddo
            enddo
                AverageLayerDepth(k) = AverageLayerDepth(k) / real(Isum)
            enddo

        endif

        select case(Me%ComputeOptions%DensityMethod)

            case (LeendertseState_) 

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

                   if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then
    
                        Sigma(i, j, k) = SigmaLeendertse (T(i, j, k), S(i, j, k))
            
                    end if
                    
                enddo
                enddo
                enddo

            case (UNESCOState_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

                    if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then
    
                        Sigma(i, j, k) = SigmaUNESCO     (T(i, j, k), S(i, j, k)) 

                    end if
                                          
                enddo
                enddo
                enddo

            case (Mel96State_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

                    if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                        Sigma(i, j, k) = SigmaUNESCO     (T(i, j, k), S(i, j, k)) 

                    end if
                                          
                enddo
                enddo
                enddo

            case (JMD95State_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

                    if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                        Sigma(i, j, k) = SigmaUNESCO     (T(i, j, k), S(i, j, k)) 

                    end if
                                          
                enddo
                enddo
                enddo

            case (Linear_)

                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB

                   if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then
    
                        !kg/m^3
                        RoRef = 1025. -  dble(SigmaDensityReference)
                        !psu
                        S0    = 33.75
                        !kg/m^3/psu
                        AlphaS= 0.78

                        Sigma(i, j, k) = RoRef + AlphaS * (S(i, j, k) - S0)

                    end if
                                          
                enddo
                enddo
                enddo

        end select
                       
        if (Me%ComputeOptions%DensPressureCorrect) then

            select case(Me%ComputeOptions%DensityMethod)

                case (UNESCOState_)

                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB

                        if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                            if (Me%ComputeOptions%SameLayerPressure) then
                                Depth = AverageLayerDepth (k) 
                            else
                                Depth =  - ZCellCenter(i,j,k)
                            endif
                
                            Sigma(i, j, k) = SigmaUNESCOPressureCorrection  (T(i, j, k), S(i, j, k),Depth, Sigma(i, j, k))

                        end if
                                              
                    enddo
                    enddo
                    enddo

                case (Mel96State_)

                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB

                        if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                            if (Me%ComputeOptions%SameLayerPressure) then
                                Depth = AverageLayerDepth (k) 
                            else
                                Depth =  - ZCellCenter(i,j,k)
                            endif
                
                            Sigma(i, j, k) = SigmaMel96PressureCorrection  (T(i, j, k), S(i, j, k),Depth, Sigma(i, j, k))

                       end if
                                              
                  enddo
                  enddo
                  enddo

                case (JMD95State_)           

                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB

                        if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then
    
                            if (Me%ComputeOptions%SameLayerPressure) then
                                Depth = AverageLayerDepth (k) 
                            else
                                Depth =  - ZCellCenter(i,j,k)
                            endif
                
                            Sigma(i, j, k) = SigmaJMD95PressureCorrection  (T(i, j, k), S(i, j, k),Depth, Sigma(i, j, k))

                        end if
                                              
                    enddo
                    enddo
                    enddo

            end select

        end if

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                D(i, j, k) = Sigma(i, j, k) + SigmaDensityReference

            else

                Sigma(i, j, k) = FillValueReal
                D(i, j, k) = FillValueReal

            end if
                                          
        enddo
        enddo
        enddo

        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                !One test to verify if the density value is not stupid
                if (D(i, j, k) > 1e5 .or.                 &
                    D(i, j, k) < 0. )  then
                    write(*,*) 'i,j,k'
                    write(*,*) i,j,k
                    write(*,*) 'density,temperature,salinity'
                    write(*,*) D(i, j, k),T(i, j, k),S(i, j, k)                     
                    stop 'ModifyDensity - ModuleWaterProperties - ERR04'   
                end if       

            end if
                                        
        enddo
        enddo
        enddo

        call UnGetGeometry(Me%ObjGeometry,SZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ModifyDensity - ModuleWaterProperties - ERR30'

        call UnGetGeometry(Me%ObjGeometry, ZCellCenter, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyDensity - ModuleWaterProperties - ERR40'


        if (Me%ComputeOptions%SameLayerPressure .and. Me%ComputeOptions%DensPressureCorrect) &
            deallocate(AverageLayerDepth)
        
    end subroutine ModifyDensity 

    !--------------------------------------------------------------------------

    subroutine ReadProperty3D(IDProperty, Field3D)

        !Arguments-------------------------------------------------------------
        real,   dimension(:,:,:), pointer           :: Field3D
        type (T_PropertyID)                         :: IDProperty

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL


        !Begin-----------------------------------------------------------------
        

           

        call ConstructFillMatrix   (PropertyID           = IDProperty,                  &
                                    EnterDataID          = Me%ObjEnterData,             &
                                    TimeID               = Me%ObjTime,                  &
                                    HorizontalGridID     = Me%ObjHorizontalGrid,        &
                                    GeometryID           = Me%ObjGeometry,              &
                                    ExtractType          = FromBlock,                   &
                                    PointsToFill3D       = Me%InputVar%WaterPoints3D,   &
                                    Matrix3D             = Field3D,                     &
                                    TypeZUV              = TypeZ_,                      &
                                    STAT                 = STAT_CALL)                   
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'ReadProperty3D - ModuleHydrodynamicAnalyser - ERR10' 




        if (IDProperty%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = IDProperty%ObjFillMatrix,           &
                                   Matrix3D       = Field3D,                            &
                                   PointsToFill3D = Me%InputVar%WaterPoints3D,          &
                                   STAT           = STAT_CALL)                          
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ReadProperty3D - ModuleHydrodynamicAnalyser - ERR20'

        end if


        call KillFillMatrix(IDProperty%ObjFillMatrix, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperty3D - ModuleHydrodynamicAnalyser - ERR30' 

    end subroutine ReadProperty3D

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine ReadProperty2D(IDProperty, Field2D)

        !Arguments-------------------------------------------------------------
        real,   dimension(:,:  ), pointer           :: Field2D
        type (T_PropertyID)                         :: IDProperty

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
        

           

        call ConstructFillMatrix   (PropertyID           = IDProperty,              &
                                    EnterDataID          = Me%ObjEnterData,         &
                                    TimeID               = Me%ObjTime,              &
                                    HorizontalGridID     = Me%ObjHorizontalGrid,    &
                                    ExtractType          = FromBlock,               &
                                    PointsToFill2D       = Me%InputVar%WaterPoints2D,&
                                    Matrix2D             = Field2D,                 &
                                    TypeZUV              = TypeZ_,                  &
                                    STAT                 = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
            stop 'ReadProperty2D - ModuleHydrodynamicAnalyser - ERR10' 




        if (IDProperty%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = IDProperty%ObjFillMatrix,       &
                                   Matrix2D       = Field2D,                        &
                                   PointsToFill2D = Me%InputVar%WaterPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ReadProperty2D - ModuleHydrodynamicAnalyser - ERR20'

        end if



        call KillFillMatrix(IDProperty%ObjFillMatrix, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperty2D - ModuleHydrodynamicAnalyser - ERR30' 

    end subroutine ReadProperty2D



    !--------------------------------------------------------------------------

    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        logical                                     :: exist
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------
       
        write(*,*)'Constructing grid...'

        !Verifies if file exists
        inquire(FILE = Me%FileName%Grid, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructGrid - ModuleHydrodynamicAnalyser - ERR10'
        endif

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%FileName%Grid, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid - ModuleHydrodynamicAnalyser - ERR20'


        call ConstructGridData      (GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     FileName         = Me%FileName%Grid,               &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleHydrodynamicAnalyser - ERR40'

        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     ActualTime       = Me%BeginTime,                   & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleHydrodynamicAnalyser - ERR50'


        call ConstructGeometry      (GeometryID       = Me%ObjGeometry,                 &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     ActualTime       = Me%BeginTime,                   &
                                     NewDomain        = Me%Filename%Geometry,           &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleHydrodynamicAnalyser - ERR60'

        call GetGeometrySize(GeometryID     = Me%ObjGeometry,                           &
                             Size           = Me%Size,                                  &
                             WorkSize       = Me%WorkSize,                              &
                             STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleHydrodynamicAnalyser - ERR70'

        call ConstructMap ( Map_ID          = Me%ObjMap,                                &
                            GeometryID      = Me%ObjGeometry,                           &
                            HorizontalMapID = Me%ObjHorizontalMap,                      &
                            TimeID          = Me%ObjTime,                               &
                            STAT            = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleHydrodynamicAnalyser - ERR80'


        call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleHydrodynamicAnalyser - ERR90'

        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleHydrodynamicAnalyser - ERR100'



    end subroutine ConstructGrid
    
    !------------------------------------------------------------------------
    subroutine Open_HDF5_Input_File
        
        !Local-----------------------------------------------------------------
        character(LEN=StringLength)                 :: PropertyName, GroupName 
        logical                                     :: Exist
        integer                                     :: m,n
        integer                                     :: STAT_CALL
        real,    dimension(:    ), pointer          :: TimePtr
        real,    dimension(6    ), target           :: AuxTime1, AuxTime2
        !----------------------------------------------------------------------
        
        !Begin-----------------------------------------------------------------

        if (Me%ComputeOptions%WaterON) then
            call Open_HDF5_Water_File
        endif

        if (Me%ComputeOptions%HydroON) then
            call Open_HDF5_Hydro_File
        endif

        if (Me%ComputeOptions%WaterON .and. Me%ComputeOptions%HydroON) then

            call HDF5SetLimits  (Me%ObjHDF5_InputWater, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_Input_File - ModuleHydrodynamicAnalyser - ERR10'

            call HDF5SetLimits  (Me%ObjHDF5_InputHydro, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_Input_File - ModuleHydrodynamicAnalyser - ERR10'

            PropertyName = "Time"

            do m=1,Me%WaterInstants

                TimePtr => AuxTime1

                call HDF5ReadData(Me%ObjHDF5_InputWater, "/Time",                                             &
                                  trim(PropertyName),                                           &
                                  Array1D      = TimePtr,                                       &
                                  OutputNumber = m, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstrucGrid - ModuleHydrodynamicAnalyser - ERR09a'

                TimePtr => AuxTime2

                call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Time",                                             &
                                  trim(PropertyName),                                           &
                                  Array1D      = TimePtr,                                       &
                                  OutputNumber = m, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File - ModuleHydrodynamicAnalyser - ERR09b'

                do n=1,6
                    
                    if (AuxTime1(n) /= AuxTime2(n)) then

                        stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR10'

                    endif

                enddo

    !            if (Me%WaterInstants /= Me%HydroInstants) then

    !                stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR10'

    !            endif

            enddo

        endif

        if (Me%ComputeOptions%HydroON) then
            Me%OutInstants = Me%HydroInstants
        endif

        if (Me%ComputeOptions%WaterON) then
            Me%OutInstants = Me%WaterInstants
        endif

        if (Me%ComputeOptions%HydroON) then

            PropertyName = GetPropertyName(VelocityW_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR20'

            if (.not. Exist) then
                write(*,*) 'There is not vertical velocity fields a null value will be assumed'
                write(*,*) 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - WARN10'
            endif

            PropertyName = GetPropertyName(VelocityU_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR30'

            if (.not. Exist) then
                write(*,*) 'There is not zonal velocity fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR40'
            endif

            PropertyName = GetPropertyName(VelocityV_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR50'

            if (.not. Exist) then
                write(*,*) 'There is not meridional velocity fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR60'
            endif

        endif


        if (Me%ComputeOptions%WaterON) then

            PropertyName = GetPropertyName(Temperature_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputWater, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR70'

            if (.not. Exist) then
                write(*,*) 'There is no temperature fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR80'
            endif

            PropertyName = GetPropertyName(Salinity_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputWater, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR90'

            if (.not. Exist) then
                write(*,*) 'There is not salinity fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleHydrodynamicAnalyser - ERR100'
            endif

        endif



    end subroutine Open_HDF5_Input_File

    !------------------------------------------------------------------------

    subroutine AllocateMatrixes
        
        !Local-----------------------------------------------------------------

        
        !Begin-----------------------------------------------------------------

        if (.not. Me%ComputeOptions%InitialON) then

            allocate(Me%InputVar%OpenPoints3D                                           &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))                                         


            Me%InputVar%OpenPoints3D(:,:,:) = Me%InputVar%WaterPoints3D(:,:,:)
                                                                                        
            allocate(Me%InputVar%SZZ                                                    &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%ZCellCenter                                                    &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%WaterLevel                                             &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB))



        endif


        if (Me%ComputeOptions%WaterON .and. .not. Me%ComputeOptions%InitialON) then

            allocate(Me%InputVar%Salinity                                               &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%Temperature                                            &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%Density                                                &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))
            
            allocate(Me%InputVar%SigmaDens                                              &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

        endif

        if (Me%ComputeOptions%HydroON) then

            allocate(Me%InputVar%VelocityU                                              &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%VelocityV                                              &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%VelocityW                                              &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))


            if (Me%ComputeOptions%ResidualHydro) then

                allocate(Me%InputVar%ResidualVelU                                       &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                allocate(Me%InputVar%ResidualVelV                                       &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                allocate(Me%InputVar%ResidualFluxU                                      &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                allocate(Me%InputVar%ResidualFluxV                                      &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))


                allocate(Me%InputVar%ResidualWaterLevel                                 &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))


            endif

        endif


        if (Me%ComputeOptions%Vorticity) then

            allocate(Me%OutputVar%Vorticity                                             &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            if (Me%ComputeOptions%ResidualHydro) then

                allocate(Me%OutputVar%ResVorticity                                      &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))
            endif

        endif

        if (Me%ComputeOptions%BaroclinicForce) then

            allocate(Me%OutputVar%BaroclinicForceX                                      &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%OutputVar%BaroclinicForceY                                      &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))



            Me%OutputVar%BaroclinicForceX(:,:,:) = FillValueReal
            Me%OutputVar%BaroclinicForceY(:,:,:) = FillValueReal
            Me%OutputVar%NN              (:,:,:) = FillValueReal

            allocate (Me%Coef%Kleft       (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%Kright      (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%Depth_integ (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%Hcenter     (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%Hleft       (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%Hright      (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%HroLeft     (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%HroRight    (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%DensLeft    (Me%Size%KLB:Me%Size%KUB))
            allocate (Me%Coef%DensRight   (Me%Size%KLB:Me%Size%KUB))



            Me%Coef%Kleft         = FillValueInt
            Me%Coef%Kright        = FillValueInt
            Me%Coef%Depth_integ   = FillValueReal
            Me%Coef%Hcenter       = FillValueReal 
            Me%Coef%Hleft         = FillValueReal
            Me%Coef%Hright        = FillValueReal
            Me%Coef%HroLeft       = FillValueReal
            Me%Coef%HroRight      = FillValueReal
            Me%Coef%DensLeft      = FillValueReal
            Me%Coef%DensRight     = FillValueReal


        endif

        if (Me%ComputeOptions%BaroclinicForce .or. Me%ComputeOptions%PerturbationPE .or.&
                                                   Me%ComputeOptions%PerturbationPE_V2) then

            allocate(Me%OutputVar%NN                                                    &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            Me%OutputVar%NN(:,:,:) = FillValueReal

        endif

        if (Me%ComputeOptions%KineticEnergy) then

            allocate(Me%OutputVar%KineticEnergy                                         &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB))


            if (Me%ComputeOptions%ResidualHydro) then

                allocate(Me%OutputVar%ResKineticEnergy                                  &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))
            endif


            if (Me%WorkSize%KUB > 1) then

                allocate(Me%OutputVar%BaroclinicKE                                      &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))

                allocate(Me%OutputVar%BarotropicU                                       &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))

                allocate(Me%OutputVar%BarotropicV                                       &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))

                allocate(Me%OutputVar%BaroclinicU                                       &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                allocate(Me%OutputVar%BaroclinicV                                       &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                if (Me%ComputeOptions%ResidualHydro) then

                    allocate(Me%OutputVar%ResBaroclinicKE                               &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%OutputVar%ResBarotropicU                                &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%OutputVar%ResBarotropicV                                &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB))
                    allocate(Me%OutputVar%ResBaroclinicU                                &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                    allocate(Me%OutputVar%ResBaroclinicV                                &
                         (Me%Size%ILB:Me%Size%IUB,                                      &
                          Me%Size%JLB:Me%Size%JUB,                                      &
                          Me%Size%KLB:Me%Size%KUB))

                endif
            endif


        endif

        if (Me%ComputeOptions%PerturbationPE .or. Me%ComputeOptions%InitialON .or.      &
            Me%ComputeOptions%PerturbationPE_V2) then

            allocate(Me%OutputVar%PerturbationPE                                       &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB))

            allocate(Me%OutputVar%PerturbationPE_V2                                     &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB))

            allocate(Me%InputVar%IniWaterLevel                                          &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB))

            allocate(Me%InputVar%IniTemperature                                         &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%IniSalinity                                            &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%IniDensity                                             &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))

            allocate(Me%InputVar%IniSigmaDens                                           &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))


            allocate(Me%InputVar%IniSZZ                                                 &
                     (Me%Size%ILB:Me%Size%IUB,                                          &
                      Me%Size%JLB:Me%Size%JUB,                                          &
                      Me%Size%KLB:Me%Size%KUB))


        endif

 

    end subroutine AllocateMatrixes

    !------------------------------------------------------------------------


    subroutine Open_HDF5_Water_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
 
        !----------------------------------------------------------------------



        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5_InputWater, Me%FileName%InputWaterHDF5, HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR10'
        

        call GetHDF5GroupNumberOfItems(Me%ObjHDF5_InputWater, "/Time",                  &
                                       Me%WaterInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR20'

    end subroutine Open_HDF5_Water_File

    !--------------------------------------------------------------------------


    !------------------------------------------------------------------------

    subroutine Open_HDF5_Hydro_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ
 
        !----------------------------------------------------------------------



        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)


        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5_InputHydro, Me%FileName%InputHydroHDF5, HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR10'
        

        call GetHDF5GroupNumberOfItems(Me%ObjHDF5_InputHydro, "/Time",                  &
                                       Me%HydroInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR20'
        

    end subroutine Open_HDF5_Hydro_File

    !--------------------------------------------------------------------------



    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        real,       dimension(:,:  ), pointer       :: Bathymetry

 
        !----------------------------------------------------------------------



        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)


        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5_OutPut, Me%FileName%OutPutHDF5, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR20'
        
        
        call HDF5SetLimits  (Me%ObjHDF5_OutPut, Me%WorkSize%ILB, Me%WorkSize%IUB,       &
                                                Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR30'

        call GetGridData(Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR40'
        
        call HDF5WriteData   (Me%ObjHDF5_OutPut, "/Grid", "Bathymetry", "-",            &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR70'

        call HDF5SetLimits  (Me%ObjHDF5_OutPut, Me%WorkSize%ILB, Me%WorkSize%IUB,       &
                                                Me%WorkSize%JLB, Me%WorkSize%JUB,       &
                                                Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR80'

        call HDF5WriteData   (Me%ObjHDF5_OutPut, "/Grid", "WaterPoints3D", "-",         &
                              Array3D = Me%InputVar%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR90'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR100'


        call UnGetGridData (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR120'


    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine ConstructEnergy

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                        :: status

        integer, dimension(:), pointer :: aux
        integer                        :: FromFile, iflag


        call GetExtractType  (FromFile = FromFile)

        allocate (aux(6))
        call GetData(aux,                                                                &
                     Me%ObjEnterData, iflag,                                &
                     SearchType = FromFile,                                              &
                     keyword    = 'ENERGY_WINDOW',                                       &
                     Default    = FillValueInt,                                          &                                           
                     ClientModule ='ModuleHydrodynamic',                                 &
                     STAT       = status)            
        
        if (status /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_, 'ConstructEnergy - ModuleHydrodynamicAnalyser - ERR02') 

        if   (iflag == 0) then 

            Me%Energy%Window%ILB = Me%WorkSize%ILB
            Me%Energy%Window%IUB = Me%WorkSize%IUB

            Me%Energy%Window%JLB = Me%WorkSize%JLB
            Me%Energy%Window%JUB = Me%WorkSize%JUB

            Me%Energy%Window%KLB = Me%WorkSize%KLB
            Me%Energy%Window%KUB = Me%WorkSize%KUB


        else if (iflag == 6) then

            Me%Energy%Window%ILB = aux(1) 
            Me%Energy%Window%IUB = aux(2) 

            Me%Energy%Window%JLB = aux(3) 
            Me%Energy%Window%JUB = aux(4) 

            Me%Energy%Window%KLB = aux(5) 
            Me%Energy%Window%KUB = aux(6) 

        else

            call SetError(FATAL_, KEYWORD_, 'ConstructEnergy - ModuleHydrodynamicAnalyser - ERR03') 

        endif

        deallocate (aux)

       
        open(unit = Me%Energy%FileID, file = Me%Filename%Energy, &
             form = 'FORMATTED', status = 'UNKNOWN')

        !Inits the buffer count
        Me%Energy%BufferCount = 0

        write (Me%Energy%FileID,*)  'Year Month Day Hour Minute Second '  // &
                                    'KEtotal PEtotal KE PE M V '          // & 
                                    'OpenV Level BarotropicKE '           // &
                                    'BaroclinicKE VelMax VelMaxBaroclinic'// &
                                    'Enstrophy'   
                        
        !Allocates the buffer
        allocate(Me%Energy%YearBuffer            (1:EnergyBufferSize))
        allocate(Me%Energy%MonthBuffer           (1:EnergyBufferSize))
        allocate(Me%Energy%DayBuffer             (1:EnergyBufferSize))
        allocate(Me%Energy%HourBuffer            (1:EnergyBufferSize))
        allocate(Me%Energy%MinuteBuffer          (1:EnergyBufferSize))
        allocate(Me%Energy%SecondBuffer          (1:EnergyBufferSize))
        allocate(Me%Energy%KineticBuffer         (1:EnergyBufferSize))
        allocate(Me%Energy%PotentialBuffer       (1:EnergyBufferSize))
        allocate(Me%Energy%VorticityBuffer       (1:EnergyBufferSize))
        allocate(Me%Energy%MassBuffer            (1:EnergyBufferSize))
        allocate(Me%Energy%VolumeBuffer          (1:EnergyBufferSize))
        allocate(Me%Energy%OpenVolumeBuffer      (1:EnergyBufferSize))
        allocate(Me%Energy%WaterLevelBuffer      (1:EnergyBufferSize))
        allocate(Me%Energy%BarotropicKEBuffer    (1:EnergyBufferSize))
        allocate(Me%Energy%BaroclinicKEBuffer    (1:EnergyBufferSize))
        allocate(Me%Energy%RelativeKEBuffer      (1:EnergyBufferSize))
        allocate(Me%Energy%RelativePEBuffer      (1:EnergyBufferSize))
        allocate(Me%Energy%VelMaxBuffer          (1:EnergyBufferSize))
        allocate(Me%Energy%VelMaxBaroclinicBuffer(1:EnergyBufferSize))

        call GetData(Me%Energy%DtOut,                                                    &
                     Me%ObjEnterData, iflag,                                             &
                     SearchType = FromFile,                                              &
                     keyword    = 'ENERGY_DT',                                           &
                     Default    = 600.,                                                  &
                     ClientModule ='ModuleHydrodynamic',                                 &
                     STAT       = status)            
        
        if (status /= SUCCESS_)                                                          &
            call SetError(FATAL_, KEYWORD_, 'ConstructEnergy - ModuleHydrodynamicAnalyser - ERR21') 

        Me%Energy%NextOutPut = Me%BeginTime


    end subroutine ConstructEnergy

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetSectionsFromFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STATUS
        logical                                     :: BlockFound
        logical                                     :: BlockVertFound, BlockCellsFound
        integer                                     :: iflag
        integer                                     :: ClientNumber

        real, dimension(:), allocatable             :: Aux2
        type (T_Transport), pointer                 :: NewSection
        integer                                     :: i, LastLine, FirstLine, Line, n, vmax

        !Begin--------------------------------------------------------------------------

        !Rewinds Buffer
        call RewindBuffer(Me%ObjEnterData, STAT = STATUS)
        if (STATUS /= SUCCESS_) stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR10"

        !Read the defined Sections
        BlockFound = .true.
        n = 0
        vmax = FillValueInt

db:     do while (BlockFound)

            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                  &
                                        block_begin, block_end, BlockFound,             &
                                        STAT = STATUS)
            if (STATUS /= SUCCESS_) stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR20"
            
Block:      if (BlockFound) then

                n = n + 1

                allocate (NewSection)
                nullify  (NewSection%Next)

                !Searches for the type of the Section        
                call GetData(NewSection%Name, Me%ObjEnterData, iflag,                   &
                             keyword       = 'ID',                                      &
                             SearchType    = FromBlock,                                 &
                             ClientModule  = 'ModuleHydrodynamicAnalyser',              &
                             STAT          = STATUS)                                    
                if (STATUS /= SUCCESS_ .or. iflag == 0)                                 &
                    stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR30"


                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,               &
                                           beginvert, endvert,                          &
                                           BlockVertFound,                              &
                                           FirstLine = FirstLine,                       &
                                           LastLine  = LastLine,                        &
                                           STAT      = STATUS)
                            
cd1 :           if (STATUS .EQ. SUCCESS_      ) then    
cd2 :               if (BlockVertFound) then 

                        NewSection%VertNumber = LastLine - FirstLine - 1

                        if (NewSection%VertNumber > vmax) vmax  = NewSection%VertNumber

                        allocate(Aux2(2))
                       
                        allocate (NewSection%DepthMin      (NewSection%VertNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR40"

                        allocate (NewSection%DepthMax      (NewSection%VertNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR50"

                        allocate (NewSection%MeridionalFlow(NewSection%VertNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR53"

                        allocate (NewSection%ZonalFlow     (NewSection%VertNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR55"

                        Line = FirstLine + 1

                        do  i = 1, NewSection%VertNumber

                            call GetData(Aux2, Me%ObjEnterData, iflag, Buffer_Line  = Line, STAT = STATUS)

                            if (STATUS /= SUCCESS_ .or. iflag /= 2)                     &
                                stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR60"

                            NewSection%DepthMin(i) = Aux2(1)
                            NewSection%DepthMax(i) = Aux2(2)

                            Line = Line + 1

                        enddo

                        deallocate(Aux2)

                    else

                        allocate (NewSection%DepthMin      (1), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR70"

                        allocate (NewSection%DepthMax      (1), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR80"

                        allocate (NewSection%MeridionalFlow(1), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR82"

                        allocate (NewSection%ZonalFlow     (1), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR84"

                        NewSection%DepthMin(1)    = 0
                        NewSection%DepthMax(1)    = 20000.
                        NewSection%VertNumber     = 1

                    endif cd2

                else if (STATUS .EQ. BLOCK_END_ERR_) then cd1

                    stop "GetSectionsFromFile - Geometry - ERR90"

                endif cd1

                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT= STATUS)
                if (STATUS /= SUCCESS_)                                                 &
                    stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR95"

                call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,               &
                                           begincells, endcells,                        &
                                           BlockCellsFound,                             &
                                           FirstLine = FirstLine,                       &
                                           LastLine  = LastLine,                        &
                                           STAT      = STATUS)
                            
cd3 :           if (STATUS .EQ. SUCCESS_      ) then    
cd4 :               if (BlockCellsFound) then 

                        NewSection%CellsNumber = LastLine - FirstLine - 1

                        allocate(Aux2(2))
                       
                        allocate (NewSection%Icell(NewSection%CellsNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR100"

                        allocate (NewSection%Jcell(NewSection%CellsNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR110"

                        allocate (NewSection%CellMeridionalFlow(NewSection%CellsNumber,NewSection%VertNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR102"

                        allocate (NewSection%CellZonalFlow     (NewSection%CellsNumber, NewSection%VertNumber), STAT = STATUS)
                        if (STATUS /= SUCCESS_)                                         &
                            stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR104"

                        Line = FirstLine + 1

                        do  i = 1, NewSection%CellsNumber

                            call GetData(Aux2, Me%ObjEnterData, iflag, Buffer_Line  = Line, STAT = STATUS)

                            if (STATUS /= SUCCESS_ .or. iflag /= 2)                     &
                                stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR120"

                            NewSection%Icell(i) = Aux2(1)
                            NewSection%Jcell(i) = Aux2(2)

                            Line = Line + 1

                        enddo

                        deallocate(Aux2)

                    else
                        write(*,*) 'The horizontal section was not defined'
                        stop "GetSectionsFromFile - ModuleHydrodynamicAnalyser - ERR130"

                    endif cd4

                else if (STATUS .EQ. BLOCK_END_ERR_) then cd3

                    stop "GetSectionsFromFile - Geometry - ERR150"

                endif cd3

                call Add_Section(NewSection)

            endif Block


        enddo db


        call Block_Unlock(Me%ObjEnterData, ClientNumber)


    end subroutine GetSectionsFromFile

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Add_Section(NewSection)

        !Arguments-------------------------------------------------------------
        type(T_Transport),         pointer     :: NewSection

        !----------------------------------------------------------------------

        ! Add to the WaterSection List a new property
        if (.not.associated(Me%FirstSection)) then
            Me%FirstSection       => NewSection
            Me%LastSection        => NewSection
        else
            NewSection%Prev       => Me%LastSection
            Me%LastSection%Next   => NewSection
            Me%LastSection        => NewSection
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Section 

    !--------------------------------------------------------------------------



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyHydrodynamicAnalyser

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, STAT_CALL
        type (T_Time)                               :: WaterInstant, HydroInstant

        !----------------------------------------------------------------------


        call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR10'

        do i = 1, Me%OutInstants

            
NI:         if (Me%ComputeOptions%InitialON) then

                Me%InputVar%Density     => Me%InputVar%IniDensity
                Me%InputVar%SigmaDens   => Me%InputVar%IniSigmaDens
                Me%InputVar%Temperature => Me%InputVar%IniTemperature
                Me%InputVar%Salinity    => Me%InputVar%IniSalinity
                
                call GetGeometryDistances (Me%ObjGeometry, SZZ = Me%InputVar%SZZ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleWaterProperties - ERR30'

                call GetGeometryDistances (Me%ObjGeometry, ZCellCenter = Me%InputVar%ZCellCenter, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleWaterProperties - ERR31'

                call GetOpenPoints3D(Me%ObjMap, Me%InputVar%OpenPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR20'

            else NI

                if (Me%ComputeOptions%WaterON) then
                    call Read_HDF5_CommonData(Me%ObjHDF5_InputWater, i, WaterInstant)
                    call Read_HDF5_Water_File(i)
                    Me%CurrentTime = WaterInstant
                endif

                if (Me%ComputeOptions%HydroON) then
                    call Read_HDF5_CommonData(Me%ObjHDF5_InputHydro, i, HydroInstant)
                    call Read_HDF5_Hydro_File(i)
                    Me%CurrentTime = HydroInstant
                endif

                if (Me%ComputeOptions%WaterON .and. Me%ComputeOptions%HydroON) then

                    if (WaterInstant /= HydroInstant) then

                        stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR40'

                    endif


                endif

                Me%InputVar%WaterLevel(:,:) = - Me%InputVar%SZZ(:,:,Me%WorkSize%KUB)

                !SZZ is imposed and not computed
                call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,              &
                                            WaterPoints3D    = Me%InputVar%WaterPoints3D,   &
                                            SurfaceElevation = Me%InputVar%WaterLevel,      &
                                            SZZ              = Me%InputVar%SZZ,             &
                                            ActualTime       = Me%CurrentTime,              &   
                                            STAT             = STAT_CALL )

                call GetGeometryDistances(Me%ObjGeometry, ZCellCenter = Me%InputVar%ZCellCenter, STAT = STAT_CALL)

                call UnGetMap(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR50'


                !Update the open points
                call UpdateComputeFaces3D(Me%ObjMap, Me%InputVar%WaterLevel,                &
                                          Me%CurrentTime, STAT = STAT_CALL)      
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR60'

                call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR70'
          
            endif NI

            call FlowProperties

            call Write_HDF5_OutPut_File(i)

            if (Me%ComputeOptions%FlowSection) then

                if (i == 1) call OutPut_ASCII_FlowSections(OpenF = .true.)

                call OutPut_ASCII_FlowSections(Instant = .true.)

            endif

        enddo

        if (Me%ComputeOptions%HydroON .and. Me%ComputeOptions%ResidualHydro) then

            call Read_HDF5_Residual_Hydro

            call ResidualFlowProperties
    
            call Write_HDF5_Residual_Hydro
            
            if (Me%ComputeOptions%FlowSection)                                          &
            call OutPut_ASCII_FlowSections(Residual = .true.)

        endif

        call OutPut_ASCII_FlowSections(CloseF = .true.)

        if (Me%ComputeOptions%InitialON) then

            call UnGetMap(Me%ObjMap, Me%InputVar%OpenPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR80'

            call UnGetGeometry (Me%ObjGeometry, Me%InputVar%SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleWaterProperties - ERR90'

            call UnGetGeometry (Me%ObjGeometry, Me%InputVar%ZCellCenter, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleWaterProperties - ERR91'

        endif


        call UnGetMap(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR100'


    end subroutine ModifyHydrodynamicAnalyser

    !--------------------------------------------------------------------------

    subroutine OutPut_ASCII_FlowSections(OpenF, Residual, CloseF, Instant)

        !Arguments-------------------------------------------------------------
        logical, optional                  :: OpenF, Residual, CloseF, Instant
         

        !Local-----------------------------------------------------------------
        real, allocatable, dimension(:)    :: AuxVector
        type(T_Transport), pointer         :: AuxSection
        character(len=1000)                :: FlowName
        character(len=StringLength)        :: Aux
        character(len=4)                   :: DepthMin, DepthMax
        logical                            :: OpenF_, Residual_, CloseF_, Instant_

        real                               :: Year_File, Month_File, Day_File,           &
                                              Hour_File, Minute_File, Second_File
        integer, save                      :: OutFile
        integer                            :: STAT_CALL, i, nFlow, n, v, OutFile2

        !----------------------------------------------------------------------

        
        if (present(OpenF)) then
            OpenF_ = OpenF
        else
            OpenF_ = .false.
        endif

        if (present(Instant)) then
            Instant_ = Instant
        else
            Instant_ = .false.
        endif

        if (present(Residual)) then
            Residual_ = Residual
        else
            Residual_ = .false.
        endif
        
        if (present(CloseF)) then
            CloseF_ = CloseF
        else
            CloseF_ = .false.
        endif

        !Open the output file
        if (OpenF_) then

            call UnitsManager(OutFile, FileOpen, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                  &
                 stop 'OutPut_ASCII_FlowSections; ModuleHydrodynamic. ERR10.'

            open(Unit = OutFile, File = Me%FileName%FlowSection,                        &
                 Form = 'FORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                  &
                 stop 'OutPut_ASCII_FlowSections; ModuleHydrodynamic. ERR20.'

            AuxSection => Me%FirstSection

            FlowName = ' '

            nFlow = 0    

            do while (associated(AuxSection))
                do i = 1, AuxSection%VertNumber
                    nFlow    = nFlow + 2
                    DepthMin = ' '
                    DepthMax = ' '
                    Aux      = ' '
                    write(DepthMin,'(I4)') int(AuxSection%DepthMin(i))
                    write(DepthMax,'(I4)') int(AuxSection%DepthMax(i))
                    Aux = trim(AuxSection%Name)//trim(adjustl(DepthMin))//'_'//trim(adjustl(DepthMax))
                    FlowName = trim(FlowName)//'"'//trim(Aux)//'_Zonal'//'"'//trim(Aux)//'_Meridional'
                 enddo
                
                AuxSection => AuxSection%Next
            enddo

            nullify(AuxSection)

            FlowName = 'Year Month Day Hour Minute Second'//trim(FlowName)

            write(OutFile,'(A)') trim(FlowName)

        endif

        if (Residual_) then
            write(OutFile,'(A)') 'Residual Flow'
            Instant_ = .true.
        endif

        !instant outputs 
        if (Instant_ ) then

            call ExtractDate(Me%CurrentTime, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

            allocate(AuxVector(nFlow))     
            
            n = 0       

            AuxSection => Me%FirstSection
            do while (associated(AuxSection))
                do i = 1, AuxSection%VertNumber
                    n = n + 1
                    AuxVector(n) = AuxSection%ZonalFlow     (i)
                    n = n + 1
                    AuxVector(n) = AuxSection%MeridionalFlow(i)
                enddo              
                AuxSection => AuxSection%Next
            enddo

            write(OutFile,'(6I6,1000f10.3)') int(Year_File), int(Month_File), int(Day_File), int(Hour_File), int(Minute_File), int(Second_File), AuxVector(1:nFlow)

            deallocate(AuxVector) 
            nullify   (AuxSection)

            if (Residual_) then

                AuxSection => Me%FirstSection

                do while (associated(AuxSection))

                    call UnitsManager(OutFile2, FileOpen, STAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                                  &
                         stop 'OutPut_ASCII_FlowSections; ModuleHydrodynamic. ERR10.'

                    open(Unit = OutFile2, File = trim(AuxSection%Name)//'.tsv',                 &
                         Form = 'FORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)

                    if (STAT_CALL /= SUCCESS_)                                                  &
                         stop 'OutPut_ASCII_FlowSections; ModuleHydrodynamic. ERR20.'

                    FlowName = 'I J '

                    do i = 1, AuxSection%VertNumber
                        DepthMin = ' '
                        DepthMax = ' '
                        Aux      = ' '
                        write(DepthMin,'(I4)') int(AuxSection%DepthMin(i))
                        write(DepthMax,'(I4)') int(AuxSection%DepthMax(i))
                        Aux = ' Zonal_'//trim(adjustl(DepthMin))//'_'//trim(adjustl(DepthMax))
                        FlowName = trim(FlowName)//trim(Aux)
                        Aux = ' Meridional_'//trim(adjustl(DepthMin))//'_'//trim(adjustl(DepthMax))
                        FlowName = trim(FlowName)//trim(Aux)
                     enddo                    

                     write(OutFile2,'(A)') trim(FlowName)

                    do i=1, AuxSection%CellsNumber
                        write(OutFile2,'(2I6,1000f10.6)') AuxSection%ICell(i), AuxSection%JCell(i), &
                            (AuxSection%CellZonalFlow(i,v), AuxSection%CellMeridionalFlow(i,v), &
                            v=1, AuxSection%VertNumber)
                    enddo

                    AuxSection => AuxSection%Next

                    call UnitsManager(OutFile2, FileClose, STAT = STAT_CALL) 

                    if (STAT_CALL /= SUCCESS_)                                                  &
                         stop 'OutPut_ASCII_FlowSections; ModuleHydrodynamic. ERR200.'

                enddo

                nullify   (AuxSection)

            endif

        endif

        if (CloseF_) then
            call UnitsManager(OutFile, FileClose, STAT = STAT_CALL) 

            if (STAT_CALL /= SUCCESS_)                                                  &
                 stop 'OutPut_ASCII_FlowSections; ModuleHydrodynamic. ERR200.'
        endif



    end subroutine OutPut_ASCII_FlowSections

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine Read_HDF5_Water_File(i)

        !Arguments-------------------------------------------------------------
        integer                                     :: i
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(LEN=StringLength)                 :: PropertyName      
        logical                                     :: SalinityExist, TemperatureExist

        !----------------------------------------------------------------------

        call HDF5SetLimits (Me%ObjHDF5_InputWater,                                      &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR10'

        PropertyName = GetPropertyName(Density_)

        call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName), &
                                Me%DensityExist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR20'


        if (Me%DensityExist)  then

            call HDF5ReadData(Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName),       &
                              trim(PropertyName),                                           &
                              Array3D      = Me%InputVar%Density,                           &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR30'

            Me%InputVar%SigmaDens(:,:,:) = Me%InputVar%Density - SigmaDensityReference

        endif



        PropertyName = GetPropertyName(Salinity_)

        call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName), &
                                SalinityExist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR40'

        if (SalinityExist) then

            !read field
            call HDF5ReadData(Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName),       &
                              trim(PropertyName),                                           &
                              Array3D      = Me%InputVar%Salinity,                          &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR60'

        endif

        PropertyName = GetPropertyName(Temperature_)

        call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName), &
                                TemperatureExist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ )stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR70'

        if (TemperatureExist) then

            call HDF5ReadData(Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName),       &
                              trim(PropertyName),                                           &
                              Array3D      = Me%InputVar%Temperature,                       &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleHydrodynamicAnalyser - ERR90'

        endif



    end subroutine Read_HDF5_Water_File

    !--------------------------------------------------------------------------

    subroutine Read_HDF5_Hydro_File(i)

        !Arguments-------------------------------------------------------------
        integer                                     :: i
        !Local-----------------------------------------------------------------
        logical                                     :: Exist
        character(LEN=StringLength)                 :: PropertyName, GroupName      
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        call HDF5SetLimits (Me%ObjHDF5_InputHydro,                                      &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR10'

        PropertyName = GetPropertyName(VelocityU_)

        !read field
        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Results/"//trim(PropertyName),       &
                          trim(PropertyName),                                           &
                          Array3D      = Me%InputVar%VelocityU,                         &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR20'

        PropertyName = GetPropertyName(VelocityV_)

        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Results/"//trim(PropertyName),       &
                          trim(PropertyName),                                           &
                          Array3D      = Me%InputVar%VelocityV,                         &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR30'

        PropertyName = GetPropertyName(VelocityW_)

        GroupName = "/Results/"//trim(PropertyName)

        call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR40'

        if (Exist) then

            call HDF5ReadData(Me%ObjHDF5_InputHydro, GroupName,                         &
                              trim(PropertyName),                                       &
                              Array3D      = Me%InputVar%VelocityW,                     &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleHydrodynamicAnalyser - ERR50'
        else

            Me%InputVar%VelocityW(:,:,:) = 0.

        endif



    end subroutine Read_HDF5_Hydro_File

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Read_HDF5_Residual_Hydro

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        call HDF5SetLimits (Me%ObjHDF5_InputHydro,                                      &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR10'

       !read field
        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Velocity/X",                &
                          "Vel. X",                                                     &
                          Array3D      = Me%InputVar%ResidualVelU,                      &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR20'

        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Velocity/Y",                &
                          "Vel. Y",                                                     &
                          Array3D      = Me%InputVar%ResidualVelV,                      &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR30'


        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Flux/X",                    &
                          "Flux X",                                                     &
                          Array3D      = Me%InputVar%ResidualFluxU,                     &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR40'

        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Flux/Y",                    &
                          "Flux Y",                                                     &
                          Array3D      = Me%InputVar%ResidualFluxV,                     &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR50'


        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Waterlevel",                &
                          "Vel. Z",                                                     &
                          Array2D      = Me%InputVar%ResidualWaterLevel,                &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR70'


    end subroutine Read_HDF5_Residual_Hydro

    !--------------------------------------------------------------------------


    subroutine Read_HDF5_CommonData(ObjHDF5, i, CurrentInstant)

        !Arguments-------------------------------------------------------------
        integer                                     :: i, ObjHDF5
        type (T_Time)                               :: CurrentInstant
        !Local-----------------------------------------------------------------
        real,    dimension(:    ), pointer          :: TimePtr
         real,    dimension(6    ), target           :: AuxTime
       integer                                     :: STAT_CALL
        character(LEN=StringLength)                 :: PropertyName
        logical                                     :: GroupExistWater, GroupExistHydro     

        !----------------------------------------------------------------------

        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR10'


        PropertyName = "Time"

        TimePtr => AuxTime

        call HDF5ReadData(ObjHDF5, "/Time",                                             &
                          trim(PropertyName),                                           &
                          Array1D      = TimePtr,                                       &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR20'


        call SetDate  (CurrentInstant,                                                  &
                       AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))


        call HDF5SetLimits (ObjHDF5,                                                    &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR30'

        PropertyName = "OpenPoints"
 
        if (Me%ComputeOptions%WaterON) then
        
            call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Grid/"//trim(PropertyName),&
                                    GroupExistWater, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR35'

        endif

        if (Me%ComputeOptions%HydroON) then
        
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, "/Grid/"//trim(PropertyName),&
                                    GroupExistHydro, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR37'

        endif


        if (GroupExistWater .or. GroupExistHydro) then

        !read field
        call HDF5ReadData(ObjHDF5, "/Grid/"//trim(PropertyName),                        &
                          trim(PropertyName),                                           &
                          Array3D      = Me%InputVar%OpenPoints3D,                      &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR40'

        endif

        call HDF5SetLimits (ObjHDF5,                                                    &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB-1,                                          &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR50'

        PropertyName = "VerticalZ"

        call HDF5ReadData(ObjHDF5, "/Grid/"//trim(PropertyName),                        &
                          "Vertical",                                                   &
                          Array3D      = Me%InputVar%SZZ,                               &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleHydrodynamicAnalyser - ERR60'


    end subroutine Read_HDF5_CommonData

    !--------------------------------------------------------------------------

    subroutine Write_HDF5_OutPut_File(i)

        !Arguments-------------------------------------------------------------
        integer                                     :: i
        !Local-----------------------------------------------------------------
        real,    dimension(6    ), target           :: AuxTime
        real,    dimension(:    ), pointer          :: TimePtr
        integer                                     :: STAT_CALL
        character(LEN=StringLength)                 :: PropertyName  

        !Begin----------------------------------------------------------------------

        call HDF5SetLimits  (Me%ObjHDF5_Output, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR10'


        PropertyName = "Time"

        TimePtr => AuxTime

        call ExtractDate (Me%CurrentTime,                                               &
                         AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))

        call HDF5WriteData(Me%ObjHDF5_Output, "/Time",                                  &
                          trim(PropertyName),                                           &
                          Units        = 'YYYY/MM/DD',                                  &
                          Array1D      = TimePtr,                                       &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR20'


        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR30'

        PropertyName = "OpenPoints"

        !write field
        call HDF5WriteData(Me%ObjHDF5_Output, "/Grid/"//trim(PropertyName),             &
                          trim(PropertyName),                                           &
                          Units        = '-',                                           &
                          Array3D      = Me%InputVar%OpenPoints3D,                      &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR40'

        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB-1,                                          &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR50'

        PropertyName = "VerticalZ"

        call HDF5WriteData(Me%ObjHDF5_Output, "/Grid/"//trim(PropertyName),             &
                          "Vertical",                                                   &
                          Units        = 'm',                                           &
                          Array3D      = Me%InputVar%SZZ,                               &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR60'


        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR70'

iV:     if (Me%ComputeOptions%Vorticity) then

            PropertyName = GetPropertyName(Vorticity_)

            !write vorticity field
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 's-1',                                     &
                              Array3D      = Me%OutputVar%Vorticity,                    &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR80'

        endif iV

iP:     if (Me%ComputeOptions%PerturbationPE) then

            PropertyName = GetPropertyName(PerturbationPE_)

            !write potential energy field 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%PerturbationPE,               &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR90'

        endif iP

iP2:    if (Me%ComputeOptions%PerturbationPE_V2) then

            PropertyName = GetPropertyName(PerturbationPE_)
            PropertyName = trim(PropertyName)//'V2'


            !write potential energy field 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%PerturbationPE_V2,            &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR90'

        endif iP2

iK:     if (Me%ComputeOptions%KineticEnergy) then

            PropertyName = GetPropertyName(KineticEnergy_)

            !write barotropic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%KineticEnergy,                &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR100'

            PropertyName = GetPropertyName(BaroclinicKE_)            
            
            !write baroclinic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%BaroclinicKE,                 &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR110'

            PropertyName = GetPropertyName(BarotropicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%BarotropicU,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR120'

            PropertyName = GetPropertyName(BarotropicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%BarotropicV,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR130'

            PropertyName = GetPropertyName(BaroclinicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%BaroclinicU,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR140'

            PropertyName = GetPropertyName(BaroclinicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%BaroclinicV,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR150'


        endif iK

iB:     if (Me%ComputeOptions%BaroclinicForce) then

            PropertyName = GetPropertyName(BaroclinicForceX_)

            !write baroclinic force in the X direction
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'kg/m3',                                     &
                              Array3D      = Me%OutputVar%BaroclinicForceX,             &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR160'

            PropertyName = GetPropertyName(BaroclinicForceY_)

            !write baroclinic force in the Y direction
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'kg/m3',                                   &
                              Array3D      = Me%OutputVar%BaroclinicForceY,             &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR170'

        endif iB

iB2:    if (Me%ComputeOptions%BaroclinicForce .or. Me%ComputeOptions%PerturbationPE .or. Me%ComputeOptions%PerturbationPE_V2) then


            PropertyName='NN'

            !write the square of Brunt-Vaisala frequency 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 's-2',                                     &
                              Array3D      = Me%OutputVar%NN,                           &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR180'


            PropertyName = GetPropertyName(Density_)

            !write the density 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'kg/m3',                                   &
                              Array3D      = Me%InputVar%Density,                       &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR190'
            
        endif iB2



        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_HDF5_OutPut_File - ModuleHydrodynamicAnalyser - ERR200'


    end subroutine Write_HDF5_OutPut_File

    !--------------------------------------------------------------------------

    subroutine Write_HDF5_Residual_Hydro

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        character(LEN=StringLength)                 :: PropertyName  

        !Begin----------------------------------------------------------------------

        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR10'

iV:     if (Me%ComputeOptions%Vorticity) then

            PropertyName = GetPropertyName(Vorticity_)

            !write vorticity field
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 's-1',                                     &
                              Array3D      = Me%OutputVar%ResVorticity,                 &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR20'

        endif iV


iK:     if (Me%ComputeOptions%KineticEnergy) then

            PropertyName = GetPropertyName(KineticEnergy_)

            !write barotropic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%ResKineticEnergy,             &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR100'

            PropertyName = GetPropertyName(BaroclinicKE_)            
            
            !write baroclinic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%ResBaroclinicKE,              &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR110'

            PropertyName = GetPropertyName(BarotropicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%ResBarotropicU,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR120'

            PropertyName = GetPropertyName(BarotropicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%ResBarotropicV,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR130'

            PropertyName = GetPropertyName(BaroclinicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%ResBaroclinicU,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR140'

            PropertyName = GetPropertyName(BaroclinicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%ResBaroclinicV,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR150'


        endif iK



        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_HDF5_Residual_Hydro - ModuleHydrodynamicAnalyser - ERR180'


    end subroutine Write_HDF5_Residual_Hydro

    !--------------------------------------------------------------------------

    subroutine FlowProperties

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        call ReadLockInputVar


        if (Me%ComputeOptions%Vorticity      ) call ComputeVorticity(Me%OutPutVar%Vorticity,&
                                                                     Me%InPutVar%VelocityU, &  
                                                                     Me%InPutVar%VelocityV)

        if (Me%ComputeOptions%KineticEnergy  ) call ComputeKineticEnergy(Me%OutPutVar%KineticEnergy, &
                                                                         Me%OutPutVar%BaroclinicKE, &
                                                                         Me%OutPutVar%BarotropicU,  &
                                                                         Me%OutPutVar%BarotropicV,  &
                                                                         Me%OutPutVar%BaroclinicU,  &
                                                                         Me%OutPutVar%BaroclinicV,  &
                                                                         Me%InPutVar%VelocityU,     &  
                                                                         Me%InPutVar%VelocityV) 

        if (Me%ComputeOptions%BaroclinicForce .or. Me%ComputeOptions%PerturbationPE .or. Me%ComputeOptions%PerturbationPE_V2) then
            call ComputeBruntVaisalaFrequency
            if (.not. Me%DensityExist) then
                call ModifyDensity(Me%InputVar%Salinity,Me%InputVar%Temperature,        &
                                   Me%InputVar%Density, Me%InputVar%SigmaDens)
            endif

        endif


        if (Me%ComputeOptions%PerturbationPE) call ComputePerturbationPE

        if (Me%ComputeOptions%PerturbationPE_V2) call ComputePerturbationPE_V2

        if (Me%ComputeOptions%BaroclinicForce) then
            call ComputeBaroclinicForce(Me%OutPutVar%BaroclinicForceX, 0 , 1,           &
                                        Me%InPutVar%KFloor_U, Me%InPutVar%DUZ, Me%InPutVar%DZX)
            call ComputeBaroclinicForce(Me%OutPutVar%BaroclinicForceY, 1 , 0,           &
                                        Me%InPutVar%KFloor_V, Me%InPutVar%DVZ, Me%InPutVar%DZY)

        endif

        if (Me%ComputeOptions%Energy) then
            call ComputeSystemEnergy
        endif

        if (Me%ComputeOptions%FlowSection) then
            call ComputeFlowAlongSections(Me%InPutVar%VelocityU,                        &  
                                          Me%InPutVar%VelocityV,                        &
                                          Me%InPutVar%SZZ)
        endif

        call ReadUnLockInputVar
        

    end subroutine FlowProperties

    !--------------------------------------------------------------------------

    subroutine ResidualFlowProperties

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        call ReadLockInputVar


        if (Me%ComputeOptions%Vorticity      ) call ComputeVorticity(Me%OutPutVar%ResVorticity,        &
                                                                         Me%InPutVar%ResidualVelU,     &  
                                                                         Me%InPutVar%ResidualVelV)

        if (Me%ComputeOptions%KineticEnergy  ) call ComputeKineticEnergy(Me%OutPutVar%ResKineticEnergy, &
                                                                         Me%OutPutVar%ResBaroclinicKE, &
                                                                         Me%OutPutVar%ResBarotropicU,  &
                                                                         Me%OutPutVar%ResBarotropicV,  &
                                                                         Me%OutPutVar%ResBaroclinicU,  &
                                                                         Me%OutPutVar%ResBaroclinicV,  &
                                                                         Me%InPutVar%ResidualVelU,     &  
                                                                         Me%InPutVar%ResidualVelV) 

        if (Me%ComputeOptions%FlowSection) then
            call ComputeFlowAlongSections(Me%InPutVar%ResidualVelU,                     &  
                                          Me%InPutVar%ResidualVelV,                     &
                                          Me%InPutVar%SZZ)
        endif

        call ReadUnLockInputVar
        

    end subroutine ResidualFlowProperties

    !--------------------------------------------------------------------------

    subroutine ReadLockInputVar

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Module - ModuleHorizontalGrid
        !Horizontal Grid properties
         call GetHorizontalGrid(Me%ObjHorizontalGrid,                                   &
                                DZX  = Me%InputVar%DZX,                                 &
                                DZY  = Me%InputVar%DZY,                                 &
                                DUX  = Me%InputVar%DUX,                                 &
                                DVY  = Me%InputVar%DVY,                                 &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleHydrodynamicAnalyser ERR10'


        !Gets WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleHydrodynamicAnalyser ERR20'        


        !Module - ModuleHydrodynamicAnalyser
        !3D Geometry properties

        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z = Me%InputVar%KFloor_Z,                                &
                               U = Me%InputVar%KFloor_U,                                &
                               V = Me%InputVar%KFloor_V,                                &
                               STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleHydrodynamicAnalyser ERR30'        


        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ = Me%InputVar%Volume_Z,                         &
                                VolumeU = Me%InputVar%Volume_U,                         &
                                VolumeV = Me%InputVar%Volume_V,                         &
                                VolumeW = Me%InputVar%Volume_W,                         &
                                STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleHydrodynamicAnalyser ERR40'        


        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  DWZ = Me%InputVar%DWZ,                                &
                                  DUZ = Me%InputVar%DUZ,                                &
                                  DVZ = Me%InputVar%DVZ,                                &
                                  STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleHydrodynamicAnalyser ERR50' 


        call GetGeometryWaterColumn(Me%ObjGeometry,                                     &
                                    WaterColumn = Me%InputVar%WaterColumn,              &
                                    WaterColumnU= Me%InputVar%WaterColumnU,             &
                                    WaterColumnV= Me%InputVar%WaterColumnV,             &
                                    STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleHydrodynamicAnalyser ERR60' 



    end subroutine ReadLockInputVar

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadUnLockInputVar

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        !Module - ModuleHorizontalGrid
        !Horizontal Grid properties
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%InputVar%DZX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR10'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%InputVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR20'


        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%InputVar%DUX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR24'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%InputVar%DVY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR26'

        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR30'        


        !Module - ModuleHydrodynamicAnalyser
        !3D Geometry properties
        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR40'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%KFloor_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR50'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%KFloor_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR60'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR70'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR80'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR90'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_W, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR100'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR110'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%DUZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR120'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%DVZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR130'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR140'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%WaterColumnU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR150'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%WaterColumnV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleHydrodynamicAnalyser ERR160'        


    end subroutine ReadUnLockInputVar

    !--------------------------------------------------------------------------


    subroutine ComputeVorticity(Vorticity, VelocityU, VelocityV)

        !Arguments-------------------------------------------------------------
        real,   dimension(:,:,:), pointer   :: Vorticity, VelocityU, VelocityV

        !Local-----------------------------------------------------------------
        real                                :: Aux
        integer                             :: IUB, ILB, JUB, JLB, KUB, KLB, i, j, k
        !Begin-----------------------------------------------------------------

        !Begin - Shorten variables name 

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB


        Vorticity (:, :, :) = 0.

      

    !------------Main cicle--------
        
dok:    do k = KLB, KUB
doj:    do j = JLB, JUB
doi:    do i = ILB, IUB
                
cd1:        if (Me%InputVar%WaterPoints3D(i-1, j  ,k)== WaterPoint .and.                &
                Me%InputVar%WaterPoints3D(i  , j  ,k)== WaterPoint) then 

                Aux = - (VelocityU(i, j, k) - VelocityU(i-1, j, k)) / Me%InPutVar%DZY(i,j)

                Vorticity (i  , j  , k) = Vorticity (i  , j  , k) + Aux / 2.
                Vorticity (i-1, j  , k) = Vorticity (i-1, j  , k) + Aux / 2.
            endif cd1

cd2:        if (Me%InputVar%WaterPoints3D(i  , j  ,k)== WaterPoint .and.                &
                Me%InputVar%WaterPoints3D(i  , j-1,k)== WaterPoint) then 

                Aux =   (VelocityV(i, j, k) - VelocityV(i , j-1, k)) / Me%InPutVar%DZX(i,j-1)

                Vorticity (i  , j  , k) = Vorticity (i  , j  , k) + Aux / 2.
                Vorticity (i  , j-1, k) = Vorticity (i  , j-1, k) + Aux / 2.
            endif cd2

        enddo doi
        enddo doj
        enddo dok

        

    end subroutine ComputeVorticity

    !--------------------------------------------------------------------------

    subroutine ComputeKineticEnergy(KineticEnergy, BaroclinicKE,                        &
                                    BarotropicU,  BarotropicV,                          &
                                    BaroclinicU,  BaroclinicV, VelocityU, VelocityV) 

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:  ), pointer :: KineticEnergy, BaroclinicKE, BarotropicU, BarotropicV
        real,    dimension(:,:,:), pointer :: BaroclinicU,  BaroclinicV, VelocityU, VelocityV

        !Local-----------------------------------------------------------------
        real                                :: Aux
        integer                             :: IUB, ILB, JUB, JLB, KUB, KLB, i, j, k
        !Begin-----------------------------------------------------------------

        !Begin - Shorten variables name 

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB


        KineticEnergy(:, :   ) = 0.

        if (Me%WorkSize%KUB > 1) then
            BaroclinicKE(:, :   ) = 0.
            BarotropicU (:, :   ) = 0.
            BarotropicV (:, :   ) = 0.
            BaroclinicU (:, :, :) = 0.
            BaroclinicV (:, :, :) = 0.
        else
            BarotropicU (:, :   ) = VelocityU(:,:,1)
            BarotropicV (:, :   ) = VelocityV(:,:,1)
        endif

i2:     if (Me%WorkSize%KUB > 1) then

d1:         do j = JLB, JUB
d2:         do i = ILB, IUB
                do k = KLB, KUB
                    if (Me%InputVar%WaterPoints3D(i, j  ,k)== WaterPoint) then 
                        
                        !Careful not to divide by zero
                        if (Me%InputVar%WaterColumn(i,j) .eq. 0) then 
                            Aux = Me%InputVar%DWZ(i, j, k)
                        else
                            Aux = Me%InputVar%DWZ(i, j, k)/Me%InputVar%WaterColumn(i,j)
                        endif
                            
                        BarotropicU(i, j) = BarotropicU(i, j)+ Aux * VelocityU(i, j, k)

                        BarotropicV(i, j) = BarotropicV(i, j)+ Aux * VelocityV(i, j, k)

                    endif
                enddo

                do k = KLB, KUB
                    if (Me%InputVar%WaterPoints3D(i, j  ,k)== WaterPoint) then 
                        BaroclinicU(i,j,k) = VelocityU(i,j,k) - BarotropicU(i, j)
                        BaroclinicV(i,j,k) = VelocityV(i,j,k) - BarotropicV(i, j)
                    endif
                enddo

            enddo d2
            enddo d1

        endif i2
            
            
            
    !------------Main cicle--------
        
doj:    do j = JLB, JUB
doi:    do i = ILB, IUB
                
cd1:        if (Me%InputVar%WaterPoints3D(i, j  ,KUB)== WaterPoint) then 

                KineticEnergy (i, j) = 0.
                BaroclinicKE  (i, j) = 0.

                do k = KLB, KUB

                    if (Me%InputVar%WaterPoints3D(i, j  ,k)== WaterPoint) then 

                        KineticEnergy (i, j) = KineticEnergy (i, j) + 0.5 * (VelocityU(i,j,k)**2+VelocityV(i,j, k)**2)

                        if (Me%WorkSize%KUB > 1) then
                          
                        !Careful not to divide by zero
                        if (Me%InputVar%WaterColumn(i,j) .eq. 0) then 
                            Aux = Me%InputVar%DWZ(i, j, k)
                        else
                            Aux = Me%InputVar%DWZ(i, j, k)/Me%InputVar%WaterColumn(i,j)
                        endif

                            BaroclinicKE(i, j) = BaroclinicKE(i, j) +0.5 * Aux * (BaroclinicU(i,j,k)**2+&
                                                                                  BaroclinicV(i,j,k)**2)
                        endif
                    endif

                enddo
            endif cd1
        enddo doi
        enddo doj


    end subroutine ComputeKineticEnergy

    !--------------------------------------------------------------------------

    subroutine ComputePerturbationPE_V2

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer  :: Temperature, Salinity, SZZ
        integer, dimension(:,:,:), pointer  :: Counter
        real,    dimension(:,:  ), pointer  :: Bathymetry
        real(8),    dimension(:    ), pointer  :: RefProf,TprofileRef, SprofileRef
        real(8),    dimension(:    ), pointer  :: ProfInst,TprofileInst, SprofileInst
        real(8)                             :: dw, dr, PPE, TotalDepth, D2, D1
        real                                :: Taux4, Saux4
        real(8)                             :: Taux, Saux, Prof, AuxP
        integer                             :: IUB, ILB, JUB, JLB, KUB, KLB, i, j, k, &
                                               nlayers, naux, Kbottom, STAT_CALL, nmin, nmax
        logical                             :: FoundBottom, FoundSurface
        !Begin-----------------------------------------------------------------
    !------------Main cicle--------


        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB


        call GetGridData(Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ComputePerturbationPE_V2 - ModuleHydrodynamicAnalyser - ERR40'

        Temperature       => Me%InputVar%Temperature
        Salinity          => Me%InputVar%Salinity
        SZZ               => Me%InputVar%SZZ

        
        Prof =  FillValueReal

dob:    do j = JLB, JUB
doa:    do i = ILB, IUB

aa:         if (Me%InputVar%WaterPoints3D(i, j  ,KUB)== WaterPoint) then 

                Kbottom = Me%InputVar%KFloor_Z(i,j)

                AuxP = (Me%InputVar%SZZ(i, j, Kbottom-1) + Me%InputVar%SZZ(i, j, Kbottom)) / 2. 

                if (Prof < AuxP) Prof = AuxP

            endif aa

        enddo doa
        enddo dob

        naux = int(Prof/Me%ComputeOptions%DZ) + 1

        allocate(RefProf(1:naux),TprofileRef(1:naux), SprofileRef(1:naux))

        allocate(Counter(ILB:IUB,JLB:JUB,1:naux))

        Counter(:,:,:) = 0

        RefProf(naux) = 0.

        do i =naux-1,1,-1
            RefProf(i) = RefProf(i+1) + Me%ComputeOptions%DZ
        enddo

        TprofileRef(1:naux) = FillValueReal
        SprofileRef(1:naux) = FillValueReal


doj:    do j = JLB, JUB
doi:    do i = ILB, IUB

                
cd1:        if (Me%InputVar%WaterPoints3D(i, j  ,KUB)== WaterPoint) then 

                if (Me%WorkSize%KUB > 1) then

                    Kbottom = Me%InputVar%KFloor_Z(i,j)

                    nlayers = KUB - Kbottom + 1

                    allocate(ProfInst(1:nlayers),TprofileInst(1:nlayers), SprofileInst(1:nlayers))

                    ProfInst    (1:nlayers) = (SZZ (i, j, Kbottom-1:KUB-1) + SZZ(i, j, Kbottom:KUB)) / 2.  
                    TprofileInst(1:nlayers) =  Temperature(i, j, Kbottom:KUB)
                    SprofileInst(1:nlayers) =  Salinity   (i, j, Kbottom:KUB)

                    do k = 1, naux

                        Taux = InterpolateProfileR8 (RefProf(k), nlayers, ProfInst, TprofileInst, FoundBottom, FoundSurface)

                        if (FoundBottom .or. FoundSurface) then 
                            cycle    
                        endif

                        TprofileRef(k) = (Taux  + TprofileRef(k) * real(Counter(i,j,k))) / (real(Counter(i,j,k)) + 1)

                        Saux = InterpolateProfileR8 (RefProf(k), nlayers, ProfInst, SprofileInst, FoundBottom, FoundSurface)

                        if (FoundBottom .or. FoundSurface) then 
                            cycle    
                        endif

                        SprofileRef(k) = (Saux  + SprofileRef(k) * real(Counter(i,j,k))) / (real(Counter(i,j,k)) + 1)

                        Counter(i,j,k) = Counter(i,j,k) + 1

                       
                    enddo

                    deallocate(ProfInst,TprofileInst, SprofileInst)

                endif
            endif cd1
        enddo doi
        enddo doj


        do i =naux,1,-1
            if (TprofileRef(i) /= FillValueReal) then
                nmax = i
                exit
            endif
        enddo


        do i =1, naux
            if (TprofileRef(i) /= FillValueReal) then
                nmin = i
                exit
            endif
        enddo

        naux = nmax-nmin+1


dol:    do j = JLB, JUB
dom:    do i = ILB, IUB

            Me%OutputVar%PerturbationPE_V2(i, j) = 0.
                
cd2:        if (Me%InputVar%WaterPoints3D(i, j  ,KUB)== WaterPoint) then 

                if (Me%WorkSize%KUB > 1) then

                    PPE        = 0.
                    TotalDepth = 0.

                    do k = KLB, KUB

                        if (Me%InputVar%WaterPoints3D(i, j  ,k)== WaterPoint) then 

                            Prof = (Me%InputVar%SZZ(i, j, k-1) + Me%InputVar%SZZ(i, j, k)) / 2. 

                            Taux = InterpolateProfileR8 (Prof, naux, RefProf(nmin:nmax), TprofileRef(nmin:nmax), FoundBottom, FoundSurface)
                            Saux = InterpolateProfileR8 (Prof, naux, RefProf(nmin:nmax), SprofileRef(nmin:nmax), FoundBottom, FoundSurface)

                            if (FoundBottom .or. FoundSurface) cycle

                            Taux4 = Taux;
                            Saux4 = Saux;

                            D1 = SigmaUNESCO (Temperature(i, j, k  ), Salinity(i, j, k))
                            D2 = SigmaUNESCO (Taux4                  , Saux4)

                            dw  =  Me%InputVar%SZZ(i, j, k-1) - Me%InputVar%SZZ(i, j, k) 
                            dr  =  D2 - D1

                            !Gill (1982)
                            if (Me%OutputVar%NN(i, j, k) > 0) then
                                PPE        = PPE + 0.5 * Gravity**2 * dr**2 * dw / SigmaDensityReference / Me%OutputVar%NN(i, j, k)
                            endif

                            TotalDepth = TotalDepth + dw
                            
                        endif
                    enddo

                    if (TotalDepth > 0) then
                        Me%OutputVar%PerturbationPE_V2  (i, j) = PPE / TotalDepth
                    endif

                endif
            endif cd2
        enddo dom
        enddo dol

        nullify(Temperature, Salinity)

        deallocate(RefProf,TprofileRef, SprofileRef)

        deallocate(Counter)

        call UnGetGridData (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ComputePerturbationPE_V2 - ModuleHydrodynamicAnalyser - ERR120'


    end subroutine ComputePerturbationPE_V2

   !--------------------------------------------------------------------------

    subroutine ComputePerturbationPE

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer  :: Temperature, Salinity, RefTemp, RefSal, RefSZZ
        real,    dimension(:    ), pointer  :: ProfIni,TprofileIni, SprofileIni
        real(8)                             :: dw, dr, PPE, TotalDepth, D2, D1
        real                                :: Tini, Sini, Prof
        integer                             :: IUB, ILB, JUB, JLB, KUB, KLB, i, j, k, nlayers
        logical                             :: FirstTime
        !Begin-----------------------------------------------------------------
    !------------Main cicle--------


        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB


        Temperature       => Me%InputVar%Temperature
        Salinity          => Me%InputVar%Salinity


        if (FirstDensityField) then

            RefTemp => Me%InputVar%IniTemperature
            RefSal  => Me%InputVar%IniSalinity
            RefSZZ  => Me%InputVar%IniSZZ 

            !RefTemp =  Me%InputVar%Temperature
            !RefSal  =  Me%InputVar%Salinity
            !RefSZZ  =  Me%InputVar%SZZ 


            FirstDensityField = .false.

        endif

        
        
doj:    do j = JLB, JUB
doi:    do i = ILB, IUB

            Me%OutputVar%PerturbationPE(i, j) = 0.

                
cd1:        if (Me%InputVar%WaterPoints3D(i, j  ,KUB)== WaterPoint) then 

                if (Me%WorkSize%KUB > 1) then

                    PPE        = 0.
                    TotalDepth = 0.

                    FirstTime = .true.

                    do k = KLB, KUB

                        if (Me%InputVar%WaterPoints3D(i, j  ,k)== WaterPoint) then 
                            
                            if (FirstTime) then
                                nlayers = kub-k+1
                                allocate(ProfIni(1:nlayers),TprofileIni(1:nlayers), SprofileIni(1:nlayers))

                                ProfIni    (1:nlayers) = (RefSZZ (i, j, k-1:KUB-1) + RefSZZ(i, j, k:KUB)) / 2.  
                                TprofileIni(1:nlayers) =  RefTemp(i, j, k:KUB)
                                SprofileIni(1:nlayers) =  RefSal (i, j, k:KUB)

                                FirstTime = .false.
                            endif


                            Prof = (Me%InputVar%SZZ(i, j, k-1) + Me%InputVar%SZZ(i, j, k)) / 2. 

                            Tini = InterpolateProfile (Prof, nlayers, ProfIni, TprofileIni)
                            Sini = InterpolateProfile (Prof, nlayers, ProfIni, SprofileIni)

                            D1 = SigmaUNESCO (Temperature(i, j, k  ), Salinity(i, j, k))
                            D2 = SigmaUNESCO (Tini,                   Sini)



                            dw  =  Me%InputVar%SZZ(i, j, k-1) - Me%InputVar%SZZ(i, j, k) 
                            dr  =  D2 - D1

                            !Gill (1982)
                            if (Me%OutputVar%NN(i, j, k) > 0) then
                                PPE        = PPE + 0.5 * Gravity**2 * dr**2 * dw / SigmaDensityReference / Me%OutputVar%NN(i, j, k)
                            endif

                            TotalDepth = TotalDepth + dw
                            
                        endif
                    enddo

                    deallocate(ProfIni,TprofileIni, SprofileIni)

                    if (TotalDepth > 0) then
                        Me%OutputVar%PerturbationPE  (i, j) = PPE / TotalDepth
                    endif

                endif
            endif cd1
        enddo doi
        enddo doj

        nullify(Temperature, Salinity)

    end subroutine ComputePerturbationPE

   !--------------------------------------------------------------------------


    subroutine ComputeBruntVaisalaFrequency

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer :: DWZ, Temperature, Salinity,ZCellCenter
        real                               :: DownON, UpON, NNup, NNdown, D1, D2
        real                               :: Depth
        integer                            :: i, j, k 
        integer                            :: IUB, ILB, JUB, JLB, KUB, KLB

    !------------initialization----


        !Begin - Shorten variables name 

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB

        Temperature       => Me%InputVar%Temperature
        Salinity          => Me%InputVar%Salinity
        ZCellCenter       => Me%InputVar%ZCellCenter
        DWZ               => Me%InputVar%DWZ
        

        !End - Shorten variables name 

        do K = KLB, KUB
        do J = JLB, JUB
        do I = ILB, IUB
            Me%OutPutVar%NN (I, J, K) = 0.
        enddo
        enddo
        enddo

    !------------Main cicle--------

do1:     do K = KLB, KUB
do2:     do J = JLB, JUB
do3:     do I = ILB, IUB
                    
cd1:        if (Me%InputVar%WaterPoints3D(i, j, k)== WaterPoint) then

                UpON   = 0.
                DownON = 0.

                if (Me%InputVar%WaterPoints3D(i   , j   , k+1)== WaterPoint) then 

                    D1 = SigmaUNESCO (Temperature(i, j, k  ), Salinity(i, j, k))
                    D2 = SigmaUNESCO (Temperature(i, j, k+1), Salinity(i, j, k+1))

                    !This depth better not include the water level!
                    !Depth = SZZ(i, j, k+1) + 0.5*DWZ(i ,j ,k+1)
                    Depth = -1.*ZCellCenter(i, j, k+1)
                    D1    = SigmaUNESCOPressureCorrection  (Temperature(i, j, k), &
                                  Salinity(i, j, k), Depth, D1 )

                    D2    = SigmaUNESCOPressureCorrection  (Temperature(i, j, k+1), &
                                  Salinity(i, j, k+1), Depth, D2 )

                    ![s-2] = [m/s2] * [M/m^3] / [M/m^3] / [m]
                    NNup = - Gravity * (D2 - D1) / SigmaDensityReference * &
                             2. / (DWZ(i, j, k+1) + DWZ(i, j, k))

                    UpON   =  1

                endif

                if (Me%InputVar%WaterPoints3D(i   , j   , k-1)== WaterPoint) then 

                    D1 = SigmaUNESCO (Temperature(i, j, k-1), Salinity(i, j, k-1))
                    D2 = SigmaUNESCO (Temperature(i, j, k  ), Salinity(i, j, k  ))

                    !This depth better not include the water level!
                    !Depth = SZZ(i, j, k-1) + 0.5*DWZ(i ,j ,k-1)
                    Depth = -1.*ZCellCenter(i, j, k-1)
                    D2    = SigmaUNESCOPressureCorrection  (Temperature(i, j, k), &
                                  Salinity(i, j, k), Depth, D2 )

                    D1    = SigmaUNESCOPressureCorrection  (Temperature(i, j, k-1), &
                                  Salinity(i, j, k-1), Depth, D1 )

                    ![s-2] = [m/s2] * [M/m^3] / [M/m^3] / [m]
                    NNdown = - Gravity * (D2 - D1) / SigmaDensityReference * &
                               2. / (DWZ(i, j, k-1) + DWZ(i, j, k))

                
                    DownON =  1

                endif

                if (DownON > 0 .or. UpON > 0) then

                    Me%OutPutVar%NN (i, j, k) = (NNup + NNdown) / 2.

                endif

            endif cd1

        enddo do3
        enddo do2
        enddo do1


    end subroutine ComputeBruntVaisalaFrequency

    !--------------------------------------------------------------------------

    subroutine ComputeBaroclinicForce(Rox3XY, di , dj, KFloor_UV, DUZ_VZ, DZX_ZY)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:,:), pointer :: Rox3XY
        integer, dimension(:,:),   pointer :: KFloor_UV 
        real,    dimension(:,:,:), pointer :: DUZ_VZ
        real,    dimension(:,:  ), pointer :: DZX_ZY 
        integer                            :: di, dj

        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer :: DWZ, SigmaDens, SZZ, AuxField

        Integer                            :: i, j, k, kbottom, ileft, jleft

        Integer, dimension( : ), pointer   :: Kleft,Kright

        Real(8), dimension( : ), pointer   :: Depth_integ, Hcenter, Hleft,&
                                              Hright, HroLeft, HroRight,  &
                                              DensRight, DensLeft

        Real(8)                            :: DAux,DAuxRight,DAuxLeft, AuxRight, AuxLeft, &
                                              ZRight, ZLeft, DRight, DLeft, DensZRight, DensZLeft

        integer                            :: IUB, ILB, JUB, JLB, KUB, KLB
        integer                            :: NRight, NLeft, kbright, kbleft, PoliDegree

        logical                            :: FoundBottomRight, FoundBottomLeft, FoundSurfaceRight, FoundSurfaceLeft
        logical                            :: PoliIsEven


    !------------initialization----


        !Begin - Shorten variables name 

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB

        SigmaDens         => Me%InputVar%SigmaDens
        DWZ               => Me%InputVar%DWZ
        SZZ               => Me%InputVar%SZZ

        Kleft             => Me%Coef%Kleft
        Kright            => Me%Coef%Kright
        Depth_integ       => Me%Coef%Depth_integ
        Hcenter           => Me%Coef%Hcenter
        Hleft             => Me%Coef%Hleft
        Hright            => Me%Coef%Hright
        HroLeft           => Me%Coef%HroLeft
        HroRight          => Me%Coef%HroRight
        DensRight         => Me%Coef%DensRight
        DensLeft          => Me%Coef%DensLeft
        

        !End - Shorten variables name 

        Rox3XY (:, :, :) = 0.

       

    !------------Main cicle--------

do4:     do J = JLB, JUB
do5:     do I = ILB, IUB
                    
cd1:        if (Me%InputVar%WaterPoints3D(i-di, j-dj, KUB)== WaterPoint .and.           &
                Me%InputVar%WaterPoints3D(i   , j   , KUB)== WaterPoint) then 

do7:            do K = KLB, KUB
                    Kleft      (K)       = FillValueInt    
                    Kright     (K)       = FillValueInt

                    Depth_integ(K)       = FillValueReal
                    Hcenter    (K)       = FillValueReal
                    Hleft      (K)       = FillValueReal
                    Hright     (K)       = FillValueReal
                    HroLeft    (K)       = FillValueReal
                    HroRight   (K)       = FillValueReal

                    DensRight  (K)       = FillValueReal
                    DensLeft   (K)       = FillValueReal 
                end do do7


                kbottom = KFloor_UV(i, j)

                ileft = i - di
                jleft = j - dj

                call calc_depth_and_Hro (Hcenter, Hleft, Hright, HroLeft, HroRight,     &
                                         DensRight, DensLeft, DWZ, SZZ, DUZ_VZ, SigmaDens,  &
                                         i, j, ileft, jleft, KUB, kbottom)

                if      (Me%ComputeOptions%BaroclinicMethod == DensityUniform  .or.     &
                         Me%ComputeOptions%BaroclinicMethod == DensityLinear) then

                    call Calc_Depth_integration(Hcenter, Hleft, Hright,                 &
                                                Kleft, Kright, Depth_integ, KUB, kbottom) 



do6:                do  k=KUB, kbottom,-1

                        Zright   = Depth_integ    (k)  - Hright (kright(k) + 1)    
                        Dright   = Hright (kright (k)) - Hright (kright(k) + 1)

                        Zleft    = Depth_integ    (k)  - Hleft  (kleft (k) + 1)    
                        Dleft    = Hleft  (kleft  (k)) - Hleft  (kleft (k) + 1)

                              

                        if      (Me%ComputeOptions%BaroclinicMethod == DensityUniform) then

                            !Constant density in each layer

                            AuxRight = Zright / Dright

                            AuxLeft  = Zleft  / Dleft


                            ![M/m^3 * m]  = [M/m^3 * m]               
                            DAuxRight=     HroRight(kright(k) + 1)                        + &
                                          (HroRight(kright(k)) - HroRight(kright(k) + 1)) * &
                                           AuxRight

                            DAuxLeft =     Hroleft (kleft (k) + 1)                        + &
                                          (Hroleft(kleft  (k)) - Hroleft(kleft(k)  + 1))  * &
                                           AuxLeft

                        else if (Me%ComputeOptions%BaroclinicMethod == DensityLinear) then

                            !Linear density evolution in each layer
                            DensZRight = (DensRight(kright(k)    ) *  Zright/2.          +  &
                                          DensRight(kright(k) + 1) * (Dright-Zright/2.)) /  &
                                          Dright

                            DensZLeft  = (DensLeft(kLeft  (k)    ) *  ZLeft/2.           +  &
                                          DensLeft(kLeft  (k) + 1) * (DLeft-ZLeft/2.))   /  &
                                          DLeft
                        

                            ![M/m^3 * m]  = [M/m^3 * m]               
                            DAuxRight  =  HroRight (kright(k) + 1) + DensZRight * Zright

                            DAuxLeft   =  Hroleft  (kleft (k) + 1) + DensZLeft  * ZLeft

                        endif

                        ! DAuxLeft (i or i -1, j-1 or j) - DAuxRight (i,j)
                        DAux     =     DAuxLeft - DAuxRight

                        ![M/m^3] =     [M/m^3 * m] / [m]              
                        DAux     =     DAux/dble(DZX_ZY(ileft, jleft))

                        !if Rox3 is positive then the baroclinic force is positive 
                        !in the pass was the opposite but to maintain coherence with the 
                        !other pressure forces now (Rox3 = (DAuxLeft - DAuxRight)/dxy) 
                        !in the pass (Rox3 = (- DAuxLeft + DAuxRight)/dxy)
                        Rox3XY(i,j,k) = real(DAux)

                enddo do6                    

            endif
                    
                    
            if (Me%ComputeOptions%BaroclinicMethod == Leibniz .or.                     &
                Me%ComputeOptions%BaroclinicMethod == Leibniz2) then

                Rox3XY(i,j,KUB +1) = 0.

                kbright = Me%InputVar%KFloor_Z(i, j)

                NRight = KUB - kbright + 1

                kbleft = Me%InputVar%KFloor_Z(ileft, jleft)


                Hright  (KUB+1) = null_real
                Hleft   (KUB+1) = null_real

                DensRight(KUB+1)   = null_real
                DensLeft (KUB+1)   = null_real                
                
                do k = KUB , kbright, -1
                   Hright  (k) = (dble(SZZ (  i,       j, k)) + dble(SZZ (    i,     j, k-1))) / 2.
                   DensRight(k) = dble(SigmaDens(i    ,j    ,k  ))
                enddo

                do k = KUB , kbleft , -1
                   Hleft   (k) = (dble(SZZ (ileft, jleft, k)) + dble(SZZ (ileft, jleft, k-1))) / 2.
                   DensLeft (k) = dble(SigmaDens(ileft,jleft,k  ))
                enddo   

                NLeft  = KUB - kbLeft + 1

do27:           do  k=KUB, kbottom,-1

                    if      (Me%ComputeOptions%BaroclinicMethod == Leibniz) then

                        !Linear interpolation
                        DensZRight = InterpolateProfileR8 (Hcenter(k), NRight, Hright(kbright:KUB),  &
                                                           DensRight(kbright:KUB), FoundBottomRight, &
                                                           FoundSurfaceRight)
                    
                        DensZLeft  = InterpolateProfileR8 (Hcenter(k), NLeft , HLeft (kbleft :KUB),  &
                                                           DensLeft (kbleft :KUB), FoundBottomLeft , &
                                                           FoundSurfaceLeft )

                        if (.not. FoundBottomRight .and. .not. FoundBottomLeft                       &
                            .and. .not.FoundSurfaceRight .and. .not. FoundSurfaceLeft) then

                            DAux       =  (DensZLeft - DensZRight) * DUZ_VZ(i, j, k)  / dble(DZX_ZY(ileft, jleft))

                        else
                            DAux = 0.
                        endif

                    else if (Me%ComputeOptions%BaroclinicMethod == Leibniz2) then

                        if ( (KUB - kbright) == 0) then
                            !Uniform profile is assumed when there only one layer
                            DensZRight = DensRight(KUB)
                        else
                            !Interpolation n degree
                            PoliDegree = min (Me%ComputeOptions%BaroclinicPoliDegree, KUB - kbright)

                            if(IsOdd(PoliDegree))then
                                PoliIsEven = .false.
                            else
                                PoliIsEven = .true.
                            endif

                            DensZRight = PolIntProfile  (Hcenter(k), NRight, Hright(kbright:KUB), &
                                                         DensRight(kbright:KUB), PoliDegree,      &
                                                         PoliIsEven)
                        endif

                        if ( (KUB - kbleft ) == 0) then
                            !Uniform profile is assumed when there only one layer
                            DensZLeft  = DensLeft (KUB)
                        else
                            !Interpolation n degree
                            PoliDegree = min (Me%ComputeOptions%BaroclinicPoliDegree, KUB - kbleft )

                            if(IsOdd(PoliDegree))then
                                PoliIsEven = .false.
                            else
                                PoliIsEven = .true.
                            endif

                            DensZLeft  = PolIntProfile  (Hcenter(k), NLeft , HLeft (kbleft :KUB), &
                                                         DensLeft (kbleft :KUB), PoliDegree,      &
                                                         PoliIsEven)
                        endif

                        DAux       =  (DensZLeft - DensZRight) * DUZ_VZ(i, j, k)  / dble(DZX_ZY(ileft, jleft))
                   
                    endif


                    Rox3XY(i,j,k) = Rox3XY(i,j,k+1) + DAux

                enddo do27

            endif


            end if cd1

        enddo do5
        enddo do4

        !Put the baroclinic force in the center of the cell

        allocate(AuxField(Me%Size%ILB:Me%Size%IUB,Me%Size%JLB:Me%Size%JUB,Me%Size%KLB:Me%Size%KUB))

        AuxField(:,:,:) = Rox3XY(:,:,:)


do12:   do  k = KLB, KUB
do11:   do  j = JLB, JUB
do10:   do  i = ILB, IUB

            if (Me%InputVar%WaterPoints3D(i, j, k)== WaterPoint) then 
                Rox3XY (i, j, k) = (AuxField(i-di, j-dj, k)+ AuxField(i, j, k)) / 2.
            endif


        enddo do10
        enddo do11
        enddo do12

        deallocate(AuxField)
        

        
        !Nullify auxiliar pointers
        nullify (Kleft)
        nullify (Kright)
        nullify (Depth_integ)
        nullify (Hcenter)
        nullify (Hleft)
        nullify (Hright)
        nullify (HroLeft)
        nullify (HroRight)
        nullify (DensRight)
        nullify (DensLeft )


        nullify (SigmaDens)          
        nullify (DWZ, SZZ)

        

    end subroutine ComputeBaroclinicForce

    !--------------------------------------------------------------------------

    subroutine calc_depth_and_Hro (Hcenter, Hleft, Hright, HroLeft, HroRight,           &
                                   DensRight, DensLeft, DWZ, SZZ, DUZ_VZ, SigmaDens,      &
                                   i, j, ileft, jleft, KUB, kbottom)

        !Arguments-------------------------------------------------------------
        integer,  intent(in )              :: i, j, ileft, jleft, KUB, kbottom
        
        real(8), pointer,  dimension ( : ) :: Hleft, Hright, Hcenter, HroLeft, HroRight, DensRight, DensLeft
        real, pointer, dimension(:, :, :)  :: DWZ, DUZ_VZ, SigmaDens, SZZ


        !Local-----------------------------------------------------------------
        integer     :: k

        !Begin-----------------------------------------------------------------



        if (Me%ComputeOptions%BaroclinicMethod == Leibniz .or. Me%ComputeOptions%BaroclinicMethod == Leibniz2) then

            Hcenter (KUB+1) = null_real           
            Hcenter (KUB  ) = (dble(SZZ (i, j, KUB)) + dble(SZZ (ileft, jleft, KUB))) / 2. + dble(DUZ_VZ( i, j, KUB) / 2.) 

            do k = KUB-1 , kbottom, -1
                Hcenter(k)  = Hcenter(k+1) + dble(DUZ_VZ( i, j, k+1) / 2.) + dble(DUZ_VZ( i, j, k) / 2.) 
            enddo
 
        else

            Hcenter (KUB+1) = 0           
            Hcenter (KUB  ) = dble(DUZ_VZ( i, j, KUB) / 2.) 

            do k = KUB-1 , kbottom, -1
                Hcenter(k)  = Hcenter(k+1) + dble(DUZ_VZ( i, j, k+1) / 2.) + dble(DUZ_VZ( i, j, k) / 2.) 
            enddo

            Hright  (KUB+1) = 0.0 !rcm 9
            Hleft   (KUB+1) = 0.0

            do k = KUB , kbottom, -1
               Hright  (k) = Hright  (k + 1) + dble(DWZ ( i, j, k))
               Hleft   (k) = Hleft   (k + 1) + dble(DWZ (ileft, jleft, k))
            enddo


            HroRight(KUB+1) = 0.0
            Hroleft (KUB+1) = 0.0


            if      (Me%ComputeOptions%BaroclinicMethod == DensityUniform) then

               
                do k = KUB , kbottom, -1
                   HroRight(k) = HroRight(k + 1) + dble(SigmaDens(i    ,j    ,k)) * dble(DWZ (    i,    j,  k))
                   Hroleft (k) = Hroleft (k + 1) + dble(SigmaDens(ileft,jleft,k)) * dble(DWZ (ileft, jleft, k))
                enddo                       

        
            else if (Me%ComputeOptions%BaroclinicMethod == DensityLinear) then

                DensRight(KUB+1)   = dble(SigmaDens(i    ,j    ,KUB    ))
                DensLeft (KUB+1)   = dble(SigmaDens(ileft,jleft,KUB    ))

                DensRight(kbottom) = dble(SigmaDens(i    ,j    ,kbottom))
                DensLeft (kbottom) = dble(SigmaDens(ileft,jleft,kbottom))

                do k = kbottom+1, KUB
                   DensRight(k) = (dble(SigmaDens(i    ,j    ,k  )) * dble(DWZ (    i,    j,  k-1)) +  &
                                   dble(SigmaDens(i    ,j    ,k-1)) * dble(DWZ (    i,    j,  k  ))) / &
                                  (dble(DWZ    (i    ,j    ,k-1)) + dble(DWZ (    i,    j,  k  )))
               
                   DensLeft (k) = (dble(SigmaDens(ileft,jleft,k  )) * dble(DWZ (ileft,jleft,  k-1)) +  &
                                   dble(SigmaDens(ileft,jleft,k-1)) * dble(DWZ (ileft,jleft,  k  ))) / &
                                  (dble(DWZ    (ileft,jleft,k-1)) + dble(DWZ (ileft,jleft,  k  )))
                enddo      

                 DensRight(kbottom : KUB+1) =  DensRight(kbottom : KUB+1)
                 DensLeft (kbottom : KUB+1) =  DensLeft (kbottom : KUB+1)

                do k = KUB , kbottom, -1
                   HroRight(k) = HroRight(k + 1) + (DensRight(k+1) + DensRight(k)) / 2. * dble(DWZ (    i,    j,  k))
                   Hroleft (k) = Hroleft (k + 1) + (DensLeft (k+1) + DensLeft (k)) / 2. * dble(DWZ (ileft, jleft, k))
                enddo   
            
            endif

        endif


    end subroutine calc_depth_and_Hro

    !-------------------------------------------------

    subroutine Calc_Depth_integration(Hcenter, Hleft, Hright,                           &
                                      Kleft, Kright, Depth_integ, KUB, kbottom)

        !Arguments-------------------------------------------------------------
        integer, intent(in )             :: KUB, kbottom
        real(8), pointer, dimension( : ) :: Hleft, Hright, Hcenter
        integer, pointer, dimension( : ) :: kleft, kright

        real(8), pointer, dimension( : ) :: Depth_integ

        !Local-----------------------------------------------------------------
        integer                          :: k    

        !Begin-----------------------------------------------------------------


!        Depth_integ = 0 !rcm 10
!        kleft       = 0
!        kright      = 0

dok:     do k=KUB,kbottom,-1

          kleft (k) = Locate_Layer(KUB, kbottom, Hcenter(k), Hleft )
          kright(k) = Locate_Layer(KUB, kbottom, Hcenter(k), Hright)

          Depth_integ(k)= min(Hcenter(k), Hleft(kleft (k)), Hright(kright(k)) )

          kleft (k) = Locate_Layer (KUB, kbottom, Depth_integ(k), Hleft )
          kright(k) = Locate_Layer (KUB, kbottom, Depth_integ(k), Hright)

        enddo dok

    end subroutine Calc_Depth_integration

    !--------------------------------------------------------------------

    function Locate_Layer (KUB, kbottom, Zpoint, Zside)
    Integer  :: Locate_Layer
        !Arguments-------------------------------------------------------------

        integer,   intent(in )                      :: KUB, kbottom
        real(8),   pointer, dimension(:)            :: Zside
        real(8), intent(in)                         :: Zpoint
        !Local-----------------------------------------------------------------
        integer                                     :: k_Locate_Layer

        !Begin-----------------------------------------------------------------


        k_Locate_Layer = KUB

        do while ( (Zside(k_Locate_Layer) < Zpoint) .and. (k_Locate_Layer > kbottom) )
            k_Locate_Layer = k_Locate_Layer -1
        enddo 

        Locate_Layer=k_Locate_Layer

    end function Locate_Layer

    !--------------------------------------------------------------------

    subroutine ComputeSystemEnergy !Frank Out99
    !Atention: thus far it only works for cartesian grids or z-level coordinates.

        !Arguments-------------------------------------------------------------       

        !Local-----------------------------------------------------------------
        real(8)                                    :: KineticEnergy,   TotalKineticEnergy        !*, ecin
        real(8)                                    :: TotalBarotropicKineticEnergy      
        real(8)                                    :: BaroclinicKineticEnergy,   TotalBaroclinicKineticEnergy      
        real(8)                                    :: PotentialEnergy, TotalPotentialEnergy      !potmat, epot
        real(8)                                    :: TotalEnstrophy
        real(8)                                    :: Mass,            TotalMass                 !rmassa, rmtot
        real(8)                                    :: Volume,          TotalVolume, OpenVolume
        real(8)                                    :: TotalLevelArea, TotalArea, Area
        real(8)                                    :: Velocity2, AuxVel
        real                                       :: Year, Month, Day, Hour, Minute, Second
        real(8)                                    :: VelMaxBaroclinic, VelMax, BaroclinicU, BaroclinicV
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                 :: WorkKLB, WorkKUB
        integer                                 :: i, j, k, iBuffer
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        integer, dimension(:, :, :), pointer    :: OpenPoints3D
        real, dimension(:, :, :),    pointer    :: Density, SZZ, DWZ
        real(8), dimension(:, :, :), pointer    :: VolumeZ
        real, dimension(:,:),        pointer    :: Bathymetry
        real, dimension(:,:),        pointer    :: DZX, DZY
        real, dimension(:, :, :),    pointer    :: CenterU, CenterV
        real, dimension(:, :, :),    pointer    :: CenterW
        real, dimension(:, :, :),    pointer    :: VelocityU, VelocityV, VelocityW
        real, dimension(:, :   ),    pointer    :: BarotropicU, BarotropicV, WaterColumnZ
        integer                                 :: VecSize
        real, dimension(:), pointer             :: NewDensity, NewVolume, NewHeight

        !Begin --------------------------------------------------------------------------

        if (Me%CurrentTime >= Me%Energy%NextOutPut) then

            Me%Energy%NextOutPut = Me%Energy%NextOutPut +      &
                                                Me%Energy%DtOut

            !Size
            IUB = Me%Size%IUB
            ILB = Me%Size%ILB

            JUB = Me%Size%JUB
            JLB = Me%Size%JLB

            KUB = Me%Size%KUB
            KLB = Me%Size%KLB

            !WorkSize
            WorkILB = Me%Energy%Window%ILB 
            WorkIUB = Me%Energy%Window%IUB 

            WorkJLB = Me%Energy%Window%JLB 
            WorkJUB = Me%Energy%Window%JUB 

            WorkKLB = Me%Energy%Window%KLB 
            WorkKUB = Me%Energy%Window%KUB 

            !Shorten variable names
            WaterPoints3D => Me%InputVar%WaterPoints3D
            OpenPoints3D  => Me%InputVar%OpenPoints3D
            Density       => Me%InputVar%Density
            VolumeZ       => Me%InputVar%Volume_Z
            SZZ           => Me%InputVar%SZZ
            DWZ           => Me%InputVar%DWZ
            DZX           => Me%InputVar%DZX
            DZY           => Me%InputVar%DZY
            WaterColumnZ  => Me%InputVar%WaterColumn
            Bathymetry    => Me%InputVar%Bathymetry
            VelocityU     => Me%InputVar%VelocityU
            VelocityV     => Me%InputVar%VelocityV
            VelocityW     => Me%InputVar%VelocityW
            
            call CenterVelocity( CenterU, CenterV, BarotropicU, BarotropicV, CenterW, Me%Size)

            !Inits variables
            TotalKineticEnergy   = 0.d0
            TotalPotentialEnergy = 0.d0
            TotalBaroclinicKineticEnergy = 0.d0
            TotalEnstrophy = 0.d0
            TotalMass            = 0.d0
            TotalVolume          = 0.d0

            OpenVolume           = 0.d0
            TotalLevelArea       = 0.d0
            TotalArea            = 0.d0
!            TotalA               = 0.d0    

            VelMax               = -99
            VelMaxBaroclinic     = -99

            !Calculates the barotropic Kinetic Energy
            TotalBarotropicKineticEnergy = CalcBarotropicKineticEnergy(Me%Energy%Window, BarotropicU, BarotropicV, CenterU, CenterV, DWZ, WaterColumnZ, WaterPoints3D, VolumeZ, Density)

            !Calculates Kinetic Energy and Potential Energy of all the sistem
            !Also calculates, TotalArea (surface), TotalMass, TotalVolume
            do i = WorkILB, WorkIUB
            do j = WorkJLB, WorkJUB
            do k = WorkKLB, WorkKUB

                if (WaterPoints3D(i, j, k) == WaterPoint) then 

                    Volume    = VolumeZ(i, j, k) 
                    Mass      = Density(i, j, k) * Volume

                    TotalVolume = TotalVolume + Volume
                    TotalMass   = TotalMass   + Mass
            
                    if (OpenPoints3D(i, j, k) == OpenPoint) then
                        Velocity2 = CenterU(i, j, k)**2. + CenterV(i, j, k)**2. + &
                                    CenterW(i, j, k)**2.
                        AuxVel = sqrt(Velocity2)

                        if (AuxVel > VelMax) VelMax = AuxVel


                        KineticEnergy        = Mass * Velocity2
                        TotalKineticEnergy   = TotalKineticEnergy + KineticEnergy


                        BaroclinicU = CenterU(i, j, k)-BarotropicU(i, j)
                        BaroclinicV = CenterV(i, j, k)-BarotropicV(i, j)

                        Velocity2 = BaroclinicU**2. + BaroclinicV**2. + &
                                    CenterW(i, j, k)**2.
                        AuxVel = sqrt(Velocity2)

                        if (AuxVel > VelMaxBaroclinic) VelMaxBaroclinic = AuxVel

                        BaroclinicKineticEnergy        = Mass * Velocity2
                        TotalBaroclinicKineticEnergy   = TotalBaroclinicKineticEnergy + BaroclinicKineticEnergy

                    endif                    

                    !PotentialEnergy      = Mass * Gravity * (Bathymetry(i, j) - (SZZ(i, j, k) + DWZ(i, j, k)/2.))
                    !Potential energy compute relatively to hydrographic zero
                    PotentialEnergy      = Mass * (- SZZ(i, j, k) - DWZ(i, j, k)/2.) 
                    TotalPotentialEnergy = TotalPotentialEnergy + PotentialEnergy
                  
                    !To calculate the graph Level/Volume
                    if (OpenPoints3D(i, j, k) == OpenPoint) then
                        OpenVolume  = OpenVolume  + Volume
                        if (k == WorkKUB) then
                            Area            = VolumeZ(i,j,k)/DWZ(i,j,k)
                            TotalArea       = TotalArea + Area
                        endif
                    endif

                endif

            enddo
            enddo
            enddo            

            TotalKineticEnergy   = 0.5 * TotalKineticEnergy
            TotalBaroclinicKineticEnergy   = 0.5 * TotalBaroclinicKineticEnergy
            TotalPotentialEnergy = Gravity * TotalPotentialEnergy

            !Computes the Enstrophy
            TotalEnstrophy = CalcEnstrophy(Me%Energy%Window, VolumeZ, CenterU, CenterV, DZY, DZX, WaterPoints3D)

            !Calculates the Reference Potential Energy of the system
            !Do I calculate it each time or only the first time?
            if (Me%Energy%FirstTime) then
                Me%Energy%FirstTime = .false.
                
                VecSize = GetNumberWaterpoints3D(Waterpoints3D, Me%Energy%Window)

!                call AdiabaticRedistributeField(APEmethod, NewDensity, NewVolume, NewHeight, &
!                                                Density, Volume, ZCellCenter, WaterpPoints3D, &
!                                                Me%Energy%Window, VecSize)

                !Stores the energy of the first iteration as reference
                Me%Energy%PotentialEnergyReference = CalcPotentialEnergy( VecSize, NewDensity, NewVolume, NewHeight)

            endif

            !Stores Data in a Buffer
            Me%Energy%BufferCount = Me%Energy%BufferCount + 1
            iBuffer = Me%Energy%BufferCount

            !Gets the current simulation time
            call ExtractDate(Me%CurrentTime,     &
                             Year   = Year,   Month = Month,  &
                             Day    = Day,    Hour  = Hour,   &
                             Minute = Minute, Second= Second)

            Me%Energy%YearBuffer(iBuffer)      = Year
            Me%Energy%MonthBuffer(iBuffer)     = Month
            Me%Energy%DayBuffer(iBuffer)       = Day
            Me%Energy%HourBuffer(iBuffer)      = Hour
            Me%Energy%MinuteBuffer(iBuffer)    = Minute
            Me%Energy%SecondBuffer(iBuffer)    = Second

            Me%Energy%RelativeKEBuffer(iBuffer)= TotalKineticEnergy / TotalMass
            Me%Energy%RelativePEBuffer(iBuffer)= (TotalPotentialEnergy -        &
                                                               Me%Energy%PotentialEnergyReference) / &
                                                               TotalMass

            Me%Energy%KineticBuffer(iBuffer)   =   TotalKineticEnergy

            Me%Energy%PotentialBuffer(iBuffer) =   TotalPotentialEnergy - Me%Energy%PotentialEnergyReference

            Me%Energy%VorticityBuffer(iBuffer) =   TotalEnstrophy

            Me%Energy%MassBuffer(iBuffer)      =   TotalMass

            Me%Energy%VolumeBuffer(iBuffer)    =   TotalVolume

            if (TotalArea > 0) then
                Me%Energy%WaterLevelBuffer(iBuffer)=   TotalLevelArea / TotalArea
            else
                Me%Energy%WaterLevelBuffer(iBuffer)=   -99.
            endif
            Me%Energy%OpenVolumeBuffer(iBuffer)    =   OpenVolume


            Me%Energy%BarotropicKEBuffer    (iBuffer)  =  TotalBarotropicKineticEnergy / TotalMass
            Me%Energy%BaroclinicKEBuffer    (iBuffer)  =  TotalBaroclinicKineticEnergy / TotalMass

            Me%Energy%VelMaxBuffer          (iBuffer)  =  VelMax
            Me%Energy%VelMaxBaroclinicBuffer(iBuffer)  =  VelMaxBaroclinic

            !If the buffer is full, writes the data file
            if (iBuffer == EnergyBufferSize) &
                call WriteEnergyDataFile

            !Deallocates local variables
            call UnCenterVelocity(CenterU, CenterV, CenterW, BarotropicU, BarotropicV)

            !Nullifies pointers
            nullify(WaterPoints3D, Density, VolumeZ, SZZ, DWZ)
            nullify(OpenPoints3D)

        endif

    end subroutine ComputeSystemEnergy

    !--------------------------------------------------------------------------
    subroutine ComputeFlowAlongSections(VelocityU, VelocityV, SZZ)

        !Arguments-------------------------------------------------------------
        real,   dimension(:,:,:), pointer   :: VelocityU, VelocityV, SZZ
        !Local-----------------------------------------------------------------
        type (T_Transport), pointer         :: AuxSection
        real                                :: DepthK, DW
        integer                             :: IUB, ILB, JUB, JLB, KUB, KLB
        integer                             :: i, j, k, n, s, v
        !Begin-----------------------------------------------------------------

        !Begin - Shorten variables name 

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB
        KUB = Me%WorkSize%KUB
        KLB = Me%WorkSize%KLB


        !------------Main cicle--------

        AuxSection => Me%FirstSection
        s = 1
        
d1:     do while (associated(AuxSection)) 

dv:         do v = 1, AuxSection%VertNumber

                AuxSection%ZonalFlow           (v) = 0.
                AuxSection%MeridionalFlow      (v) = 0.
                AuxSection%CellZonalFlow     (:,v) = 0.
                AuxSection%CellMeridionalFlow(:,v) = 0.


dn:             do n = 1, AuxSection%CellsNumber

                    i = AuxSection%ICell(n)
                    j = AuxSection%JCell(n)

dk:                 do k = KLB, KUB

                        DepthK = SZZ(i,j,k-1) - SZZ(i,j,KUB)

cd1:                    if (Me%InputVar%WaterPoints3D(i, j, k) == WaterPoint .and.              &
                            DepthK >= AuxSection%DepthMin(v) .and. DepthK <= AuxSection%DepthMax(v)) then

                            DW = SZZ(i,j,k-1) - SZZ(i,j,k)

                            !From m3/s to Sv

                            !Zonal flow
                            AuxSection%CellZonalFlow     (n,v) = VelocityU(i,j,k) * DW * Me%InputVar%DVY(i, j) / 1e6
                            AuxSection%ZonalFlow         (v)   = AuxSection%ZonalFlow(v)      + AuxSection%CellZonalFlow     (n,v)

                            !Meridional flow
                            AuxSection%CellMeridionalFlow(n,v) = VelocityV(i,j,k) * DW * Me%InputVar%DUX(i, j) / 1e6
                            AuxSection%MeridionalFlow    (v)   = AuxSection%MeridionalFlow(v) + AuxSection%CellMeridionalFlow(n,v)

                        endif cd1

                    enddo dk

                enddo dn

            enddo dv

            AuxSection => AuxSection%Next
            s = s + 1

        enddo d1

        nullify   (AuxSection)
               
    end subroutine ComputeFlowAlongSections

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    real(8) function CalcEnstrophy(WorkSize, VolumeZ, CenterU, CenterV, DZY, DZX, WaterPoints3D)
        
        !Arguments-----------------------------------------------------------------------
        type(T_Size3D)                          :: Worksize
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        real, dimension(:,:),        pointer    :: DZX, DZY
        real(8), dimension(:, :, :), pointer    :: VolumeZ
        real, dimension(:, :, :),    pointer    :: CenterU, CenterV

        !Locals----------------------------------------------------------------------------
        integer                                 :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                 :: WorkKLB, WorkKUB
        integer                                 :: i, j, k
        real(8)                                 :: Volume, Vorticity
 
        !Begin shorten variables names---------------------------------------
        WorkILB = WorkSize%ILB
        WorkIUB = WorkSize%IUB
        WorkJLB = WorkSize%JLB
        WorkJUB = WorkSize%JUB
        WorkKLB = WorkSize%KLB
        WorkKUB = WorkSize%KUB

dok:        do k = WorkKLB, WorkKUB
doj:        do j = WorkJLB, WorkJUB
doi:        do i = WorkILB, WorkIUB
    
                Vorticity = 0
                Volume    = VolumeZ(i, j, k) 
                
cd1:            if (WaterPoints3D(i-1, j  ,k)== WaterPoint .and.                &
                    WaterPoints3D(i  , j  ,k)== WaterPoint) then 

                    Vorticity = - (CenterU(i, j, k) - CenterU(i-1, j, k)) / DZY(i-1,j)

                endif cd1

cd2:            if (WaterPoints3D(i  , j  ,k)== WaterPoint .and.                &
                    WaterPoints3D(i  , j-1,k)== WaterPoint) then 

                    Vorticity =  Vorticity + (CenterV(i, j, k) - CenterV(i , j-1, k)) / DZX(i,j-1)

                endif cd2

                CalcEnstrophy = CalcEnstrophy + Volume * Vorticity**2.

            enddo doi
            enddo doj
            enddo dok

            ! CalcEnstrophy = 0.5 * CalcEnstrophy

    end function CalcEnstrophy

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    real(8) function CalcBarotropicKineticEnergy(WorkSize, BarotropicU, BarotropicV, CenterU, CenterV, DWZ, WaterColumnZ, WaterPoints3D, VolumeZ, Density)

        !Arguments-----------------------------------------------------------------------
        type(T_Size3D)                          :: Worksize
        integer, dimension(:, :, :), pointer    :: WaterPoints3D
        real, dimension(:, :, :),    pointer    :: Density, DWZ
        real(8), dimension(:, :, :), pointer    :: VolumeZ
        real, dimension(:, :, :),    pointer    :: CenterU, CenterV
        real, dimension(:, :   ),    pointer    :: BarotropicU, BarotropicV, WaterColumnZ
        
        !Locals----------------------------------------------------------------------------
        integer                                 :: WorkILB, WorkIUB, WorkJLB, WorkJUB
        integer                                 :: WorkKLB, WorkKUB
        integer                                 :: i, j, k
        real(8)                                 :: Mass, Volume, BarotropicKineticEnergy

        !Begin shorten variables names---------------------------------------
        WorkILB = WorkSize%ILB
        WorkIUB = WorkSize%IUB
        WorkJLB = WorkSize%JLB
        WorkJUB = WorkSize%JUB
        WorkKLB = WorkSize%KLB
        WorkKUB = WorkSize%KUB

            CalcBarotropicKineticEnergy = 0.d0

            !Calculates de barotropic velocity 
            do j = WorkJLB, WorkJUB
            do i = WorkILB, WorkIUB

                BarotropicU(i, j) = 0. 
                BarotropicV(i, j) = 0.
                Mass              = 0.

                if (WaterPoints3D(i, j, WorkKUB) == WaterPoint) then 

                    do k = WorkKLB, WorkKUB

                        if (WaterPoints3D(i, j, k) == WaterPoint) then

                            Volume    = VolumeZ(i, j, k) 
                            Mass      = Mass + Density(i, j, k) * Volume

                            BarotropicU(i, j) = BarotropicU(i, j) + CenterU(i, j, k) * DWZ(i, j, k) / WaterColumnZ(i, j)
                            BarotropicV(i, j) = BarotropicV(i, j) + CenterV(i, j, k) * DWZ(i, j, k) / WaterColumnZ(i, j) 

                        endif
                    enddo

                    BarotropicKineticEnergy      = Mass * (BarotropicU(i, j)**2 + BarotropicV(i, j)**2)
                    CalcBarotropicKineticEnergy = CalcBarotropicKineticEnergy + BarotropicKineticEnergy

                endif
            enddo
            enddo            

            CalcBarotropicKineticEnergy = 0.5 * CalcBarotropicKineticEnergy

    end function CalcBarotropicKineticEnergy

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine CenterVelocity( CenterU, CenterV, BarotropicU, BarotropicV, CenterW, Size)            

            !Arguments------------------------------------------------------------------------
            type(T_Size3D)                                         :: Size
            real, dimension(:,:,:),pointer                         :: CenterU, CenterV
            real, dimension(:,:,:),pointer                         :: CenterW
            real, dimension(:,:  ),pointer                         :: BarotropicU, BarotropicV

            !Locals----------------------------------------------------------------------------
            integer                                                         :: ILB, IUB, JLB, JUB, KLB, KUB
            integer                                                         :: i,j,k
            integer                                                         :: STAT_CALL
            real, dimension(:,:,:), pointer                                 :: VelocityW

            !Begin shorten variables names-----------------------------------------------------
            ILB = Size%ILB
            IUB = Size%IUB
            JLB = Size%JLB
            JUB = Size%JUB
            KLB = Size%KLB
            KUB = Size%KUB
            VelocityW => Me%InputVar%VelocityW

            !Centers horizontal velocity
            allocate(CenterU    (ILB:IUB, JLB:JUB, KLB:KUB),                            &
                     CenterV    (ILB:IUB, JLB:JUB, KLB:KUB),                            &
                     BarotropicU(ILB:IUB, JLB:JUB         ),                            &
                     BarotropicV(ILB:IUB, JLB:JUB         ),                            &
                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeSystemEnergy - ModuleHydrodynamic - ERR01'

            call CenterVelocityUV( CenterU, CenterV)

            !Centers vertical velocity
            allocate(CenterW(ILB:IUB, JLB:JUB, KLB:KUB), stat = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeSystemEnergy - ModuleHydrodynamic - ERR02'
            CenterW = 0.
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (k == KLB) then
                    CenterW(i, j, k) = VelocityW(i, j, k+1) / 2.
                elseif (k == KUB) then
                    CenterW(i, j, k) = VelocityW(i, j, k) / 2.
                else
                    CenterW(i, j, k) = (VelocityW(i, j, k+1) +  &
                                                VelocityW(i, j, k)) / 2.
                endif
            enddo
            enddo
            enddo

    end subroutine CenterVelocity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine UnCenterVelocity(CenterU, CenterV, CenterW, BarotropicU, BarotropicV)

        !Arguments----------------------------------------------------------------
        real, dimension(:,:,:),pointer                         :: CenterU, CenterV
        real, dimension(:,:,:),pointer                         :: CenterW
        real, dimension(:,:  ),pointer                         :: BarotropicU, BarotropicV

        !Locals-------------------------------------------------------------------
        integer                                                         :: STAT_CALL

   
            !Deallocates CenterVelocities
            deallocate(CenterU, CenterV, CenterW, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeSystemEnergy - ModuleHydrodynamicAnalyser - ERR04'

            deallocate(BarotropicU, BarotropicV, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ComputeSystemEnergy - ModuleHydrodynamicAnalyser - ERR04'


            nullify(CenterU, CenterV, CenterW)

            nullify(BarotropicU, BarotropicV) 

    end subroutine UnCenterVelocity

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine WriteEnergyDataFile !Frank Out99

        !Arguments-------------------------------------------------------------
        

        !Local
        integer                             :: iBuffer


        do iBuffer = 1, Me%Energy%BufferCount
            write(Me%Energy%FileID, fmt=100) Me%Energy%YearBuffer            (iBuffer), &
                                             Me%Energy%MonthBuffer           (iBuffer), &
                                             Me%Energy%DayBuffer             (iBuffer), &
                                             Me%Energy%HourBuffer            (iBuffer), &
                                             Me%Energy%MinuteBuffer          (iBuffer), &
                                             Me%Energy%SecondBuffer          (iBuffer), &
                                             Me%Energy%KineticBuffer         (iBuffer), &
                                             Me%Energy%PotentialBuffer       (iBuffer), &
                                             Me%Energy%RelativeKEBuffer      (iBuffer), &
                                             Me%Energy%RelativePEBuffer      (iBuffer), &
                                             Me%Energy%MassBuffer            (iBuffer), &
                                             Me%Energy%VolumeBuffer          (iBuffer), &
                                             Me%Energy%OpenVolumeBuffer      (iBuffer), &
                                             Me%Energy%WaterLevelBuffer      (iBuffer), &
                                             Me%Energy%BarotropicKEBuffer    (iBuffer), &
                                             Me%Energy%BaroclinicKEBuffer    (iBuffer), &                                                   
                                             Me%Energy%VelMaxBuffer          (iBuffer), &
                                             Me%Energy%VelMaxBaroclinicBuffer(iBuffer), &
                                             Me%Energy%VorticityBuffer       (iBuffer)
                                                         
            100 format(1x, f5.0, 1x, f4.0, 1x, f4.0, 1x, f4.0, 1x, f4.0, f4.0, 1x, e14.8, 1x, e14.8, 1x, e14.8, 1x, e14.8, 1x,  &
                       e14.8, 1x, e14.8, 1x, e14.8, 1x, f12.6, 1x, e14.8, 1x, e14.8, 1x, e14.8, 1x, e14.8, 1x, e14.8)
        enddo
        Me%Energy%BufferCount = 0

    end subroutine WriteEnergyDataFile

    !--------------------------------------------------------------------------

    subroutine CenterVelocityUV( CenterU, CenterV)

        !Arguments-------------------------------------------------------------
         
        real,    dimension(:, :, :), pointer    :: CenterU, CenterV

        !Local-----------------------------------------------------------------
        real,    dimension(:, :, :), pointer    :: FaceVelocityU, FaceVelocityV
        real(8), dimension(:, :, :), pointer    :: FaceUDouble, FaceVDouble
        real                                    :: VelU, VelV, AngleX, AngleY
        integer                                 :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                 :: i, j, k
        logical                                 :: Simple, Double
        
        
        if (MonitorPerformance) call StartWatch ("ModuleHydrodynamic", "CenterVelocityUV")

        
        !Bounds
        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB 

        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 

        KLB = Me%Size%KLB 
        KUB = Me%Size%KUB 
        
        AngleX = 0.
        AngleY = 0.

        Simple = .false.
        Double = .true.
        
        !Checks for the type of velocity to interpolate
        FaceVelocityU => Me%InputVar%VelocityU
        FaceVelocityV => Me%InputVar%VelocityV

        Simple = .true.


        !Interpolates CenterVelocities
        
        !if (MonitorPerformance) call StartWatch ("ModuleHydrodynamic", "DoLoop-CenterVelocityUV")

        if(Simple)then
            
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
            
                if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                    AngleX = 0.
                    AngleY = Pi/2.
                    
                    VelU = (FaceVelocityU(i, j, k)  +   FaceVelocityU(i, j+1, k)) / 2.

                    VelV = (FaceVelocityV(i, j, k)  +   FaceVelocityV(i+1, j, k)) / 2.
                
                else

                    VelU = 0.
                    VelV = 0.

                endif
                
                CenterU(i, j, k) = VelU * cos(AngleX) + VelV * cos(AngleY)
                CenterV(i, j, k) = VelU * sin(AngleX) + VelV * sin(AngleY)

            enddo
            enddo
            enddo

        elseif(Double)then

            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
            
                if (Me%InputVar%WaterPoints3D(i, j, k) == 1) then

                    AngleX = 0.
                    AngleY = Pi/2.

                    VelU = (FaceUDouble(i, j, k)  +  FaceUDouble(i, j+1, k)) / 2.

                    VelV = (FaceVDouble(i, j, k)  +  FaceVDouble(i+1, j, k)) / 2.
                              
                else

                    VelU = 0.
                    VelV = 0.

                endif

                CenterU(i, j, k) = VelU * cos(AngleX) + VelV * cos(AngleY)
                CenterV(i, j, k) = VelU * sin(AngleX) + VelV * sin(AngleY)

            enddo
            enddo
            enddo

        endif


        !if (MonitorPerformance) call StopWatch ("ModuleHydrodynamic", "DoLoop-CenterVelocityUV")

        !Nullifies FaceVelocityU, FaceVelocityV
        nullify(FaceVelocityU, FaceVelocityV)
        nullify(FaceUDouble,   FaceVDouble)

        if (MonitorPerformance) call StopWatch ("ModuleHydrodynamic", "CenterVelocityUV")


    end subroutine CenterVelocityUV

    !--------------------------------------------------------------------------

!--------------------------------------------------------------------------------------

    real(8) function CalcPotentialEnergy(Size, Density, Volume, Height)

        !Arguments-----------------------------------------------
        real, dimension(:), pointer                                 :: Density, Volume, Height
        integer                                                     :: Size
        
        !locals--------------------------------------------------
        integer                                                     :: m

        CalcPotentialEnergy = 0.
        do m = 1, Size
            CalcPotentialEnergy = CalcPotentialEnergy + Density(m)*Volume(m)*Height(m)
        enddo
        CalcPotentialEnergy = Gravity * CalcPotentialEnergy

    end function CalcPotentialEnergy

    integer function GetNumberWaterpoints3D(Waterpoints3D, Size3D)

        !Arguments---------------------------------------------
        type(T_Size3D)                                              :: Size3D
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D

        !Locals------------------------------------------------
        integer                                                     :: i,j,k,m

        m = 0;
        do i = Size3D%ILB, Size3D%IUB
        do j = Size3D%JLB, Size3D%JUB
        do k = Size3D%KLB, Size3D%KUB

            if(Waterpoints3D(i,j,k) == 1) then
                m = m + 1
            end if
            
        end do 
        end do
        end do

        GetNumberWaterpoints3D = m

    end function

    subroutine CalculateStageAreas( VolumeZ, VolumeHeight, Areas, Size3D, Waterpoints3D)
    
        !Arguments-------------------------------------------------------------------------
        type(T_Size3D)                                              :: Size3D
        real, dimension(:,:,:), pointer                             :: VolumeZ, VolumeHeight    
        real, dimension(:    ), pointer                             :: Areas    
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D        

        !Locals----------------------------------------------------------------------------
        integer                                                     :: i,j,k,m,p,pstart,pstop
        real                                                        :: Area

        m = 0;
        pstart = 1;
        pstop = 1;
        do k = Size3D%KLB, Size3D%KUB

            Area = 0.;

            do j = Size3D%JLB, Size3D%JUB
            do i = Size3D%ILB, Size3D%IUB

                if (Waterpoints3D(i,j,k) == 1) then

                    m = m + 1
                    Area = Area + VolumeZ(i,j,k)/VolumeHeight(i,j,k)

                end if

            enddo
            enddo

            pstop = m
            do p = pstart, pstop
                Areas(p) = Area
            enddo
            pstart = pstop

        enddo

    end subroutine CalculateStageAreas

    subroutine SortNumerically_3D( TargetReference, TargetAcolyte, SourceReference, SourceAcolyte, Size3D, Waterpoints3D)
        
        !Arguments-------------------------------------------------------------------------
        type(T_Size3D), pointer                                     :: Size3D
        real, dimension(:,:,:), pointer                             :: SourceReference, SourceAcolyte
        real, dimension(:    ), pointer                             :: TargetReference, TargetAcolyte
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D

        !Locals----------------------------------------------------------------------------
        integer                                                     :: i,j,k,m,p

        m = 0;
        do i = Size3D%ILB, Size3D%IUB
        do j = Size3D%JLB, Size3D%JUB
        do k = Size3D%KLB, Size3D%KUB

            if(Waterpoints3D(i,j,k) == 1) then
                m = m + 1
                TargetReference(m) = SourceReference(i,j,k)
                TargetAcolyte(m) = SourceAcolyte(i,j,k)
            end if
            
            p = m
            do while (p>1 .AND. TargetReference(p) > TargetReference(p-1)) 
                call swap(TargetReference, p, p-1)
                call swap(TargetAcolyte, p, p-1)
                p=p-1
            end do

        end do 
        end do
        end do

    end subroutine SortNumerically_3D

    subroutine DetermineCentreHeight_3D( TargetHeight, SourceHeight, Size3D, Waterpoints3D, method, VolumeZ, Areas)
        
        !Arguments-------------------------------------------------------------------------
        type(T_Size3D), pointer                                     :: Size3D
        real, dimension(:,:,:), pointer                             :: SourceHeight
        real, dimension(:    ), pointer                             :: TargetHeight
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D
        integer                                                     :: method
        real, dimension(:,:,:), pointer, optional                   :: VolumeZ
        real, dimension(:), pointer, optional                       :: Areas

        !Locals----------------------------------------------------------------------------
        integer                                                     :: m
        integer                                                     :: i,j,k
        real                                                        :: HeightEdge, HeightEdge_old
        real                                                        :: ElementHeight


        select case (method)

            case(SimpleHeight_) 
                m = 0;
                do i = Size3D%ILB, Size3D%IUB
                do j = Size3D%JLB, Size3D%JUB
                do k = Size3D%KLB, Size3D%KUB

                    if(Waterpoints3D(i,j,k) == 1) then
                        m = m + 1
                        TargetHeight(m) = SourceHeight(i,j,k)
                    end if
            
                end do 
                end do
                end do

            case(ComplexHeight_)
                HeightEdge_old  = 0.; 
                HeightEdge      = 0.; 
                m = 0;
                do i = Size3D%ILB, Size3D%IUB
                do j = Size3D%JLB, Size3D%JUB
                do k = Size3D%KLB, Size3D%KUB

                    if(Waterpoints3D(i,j,k) == 1) then
                        m = m + 1
                        ElementHeight = VolumeZ(i,j,k) / Areas(m)
                        if (m .GT. 1) then
                            HeightEdge = HeightEdge_old + ElementHeight
                        end if
                        TargetHeight(m) = HeightEdge + 0.5 * ElementHeight
                        HeightEdge_old = HeightEdge
                    end if
            
                end do 
                end do
                end do

        end select

    end subroutine DetermineCentreHeight_3D
    
    subroutine swap(Vec, m, n)

        !Arguments-------------------------------------------------------
        real, dimension(:)                              :: Vec
        integer                                         :: m, n
        
        !Locals----------------------------------------------------------
        real                                            :: a

        a = Vec(m)
        Vec(m) = Vec(n)
        Vec(n) = a

    end subroutine

!--------------------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillHydrodynamicAnalyser

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        type(T_Transport), pointer      :: AuxSection
        integer                         :: STAT_CALL

        !------------------------------------------------------------------------
           
        call DeAllocateMatrixes

        call KillMap(Me%ObjMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR10'

        call KillGeometry(Me%ObjGeometry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR20'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR30'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR40'

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR50'
        
        call KillHDF5(Me%ObjHDF5_Output, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR60'

NI:     if (.not. Me%ComputeOptions%InitialON) then

            if (Me%ComputeOptions%WaterON) then
                call KillHDF5(Me%ObjHDF5_InputWater, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR80'
            endif

            if (Me%ComputeOptions%HydroON) then
                call KillHDF5(Me%ObjHDF5_InputHydro, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillHydrodynamicAnalyser - ModuleHydrodynamicAnalyser - ERR70'
            endif

        endif NI

        if (Me%ComputeOptions%FlowSection) then

            AuxSection => Me%FirstSection

            do while(associated(AuxSection))

                deallocate(AuxSection%DepthMin      ,AuxSection%DepthMax          )
                deallocate(AuxSection%Icell         ,AuxSection%Jcell             )
                deallocate(AuxSection%ZonalFlow     ,AuxSection%MeridionalFlow    )
                deallocate(AuxSection%CellZonalFlow ,AuxSection%CellMeridionalFlow)
                AuxSection => AuxSection%Next   
            enddo

            nullify   (AuxSection)

        endif

        deallocate(Me)

        !------------------------------------------------------------------------

    end subroutine KillHydrodynamicAnalyser
        

    !------------------------------------------------------------------------

    subroutine DeAllocateMatrixes
        
        !Local-----------------------------------------------------------------

        
        !Begin-----------------------------------------------------------------

        if (.not. Me%ComputeOptions%InitialON) then

            deallocate(Me%InputVar%OpenPoints3D)                                            
                                                                                        
            deallocate(Me%InputVar%SZZ)

            deallocate(Me%InputVar%WaterLevel)

        endif


        if (Me%ComputeOptions%WaterON .and. .not. Me%ComputeOptions%InitialON) then

            deallocate(Me%InputVar%Salinity)

            deallocate(Me%InputVar%Temperature)

            deallocate(Me%InputVar%Density)

        endif

        if (Me%ComputeOptions%HydroON) then

            deallocate(Me%InputVar%VelocityU)

            deallocate(Me%InputVar%VelocityV)

            deallocate(Me%InputVar%VelocityW)

            if (Me%ComputeOptions%ResidualHydro) then

                deallocate(Me%InputVar%ResidualVelU)

                deallocate(Me%InputVar%ResidualVelV)

                deallocate(Me%InputVar%ResidualFluxU)

                deallocate(Me%InputVar%ResidualFluxV)

                deallocate(Me%InputVar%ResidualWaterLevel)

            endif

        endif


        if (Me%ComputeOptions%Vorticity) then

            deallocate(Me%OutputVar%Vorticity)

            if (Me%ComputeOptions%ResidualHydro) then

                deallocate(Me%OutputVar%ResVorticity)

            endif

        endif

        if (Me%ComputeOptions%BaroclinicForce) then

            deallocate(Me%OutputVar%BaroclinicForceX)

            deallocate(Me%OutputVar%BaroclinicForceY)

            deallocate (Me%Coef%Kleft      ) 
            deallocate (Me%Coef%Kright     ) 
            deallocate (Me%Coef%Depth_integ) 
            deallocate (Me%Coef%Hcenter    ) 
            deallocate (Me%Coef%Hleft      ) 
            deallocate (Me%Coef%Hright     ) 
            deallocate (Me%Coef%HroLeft    ) 
            deallocate (Me%Coef%HroRight   ) 
            deallocate (Me%Coef%DensLeft   ) 
            deallocate (Me%Coef%DensRight  ) 

        endif

        if (Me%ComputeOptions%BaroclinicForce .or. Me%ComputeOptions%PerturbationPE) then

            deallocate(Me%OutputVar%NN)

        endif

        if (Me%ComputeOptions%KineticEnergy) then

            deallocate(Me%OutputVar%KineticEnergy)

            if (Me%WorkSize%KUB > 1) then
                deallocate(Me%OutputVar%BaroclinicKE)
                deallocate(Me%OutputVar%BarotropicU)
                deallocate(Me%OutputVar%BarotropicV)
                deallocate(Me%OutputVar%BaroclinicU)
                deallocate(Me%OutputVar%BaroclinicV)

                if (Me%ComputeOptions%ResidualHydro) then

                    deallocate(Me%OutputVar%ResBaroclinicKE)
                    deallocate(Me%OutputVar%ResBarotropicU)
                    deallocate(Me%OutputVar%ResBarotropicV)
                    deallocate(Me%OutputVar%ResBaroclinicU)
                    deallocate(Me%OutputVar%ResBaroclinicV)

                endif
            endif

        endif

        if (Me%ComputeOptions%PerturbationPE .or. Me%ComputeOptions%InitialON) then

            deallocate(Me%OutputVar%PerturbationPE)

            deallocate(Me%InputVar%IniWaterLevel)

            deallocate(Me%InputVar%IniTemperature)

            deallocate(Me%InputVar%IniSalinity)

            deallocate(Me%InputVar%IniDensity)

            deallocate(Me%InputVar%IniSigmaDens)

            deallocate(Me%InputVar%IniSZZ) 

        endif

 

    end subroutine DeAllocateMatrixes

    !------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------
end module ModuleHydrodynamicAnalyser









