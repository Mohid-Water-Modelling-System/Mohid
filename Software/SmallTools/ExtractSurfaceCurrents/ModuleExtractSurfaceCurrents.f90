!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ExtractSurfaceCurrents
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to serve as ExtractSurfaceCurrents to create new modules
!
!------------------------------------------------------------------------------


Module ModuleExtractSurfaceCurrents

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
    use ModuleFillMatrix,           only: ConstructFillMatrix, GetDefaultValue, KillFillMatrix, &
                                          ModifyFillMatrix


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartExtractSurfaceCurrents
    private ::      ReadOptions
    private ::          ReadInitialDensityField
    private ::      ConstructGrid
    private ::      AllocateMatrixes
    private ::      Open_HDF5_Input_File
    private ::          Open_HDF5_Hydro_File
    private ::          Open_HDF5_Water_File
    private ::      Open_HDF5_OutPut_File

                     
    
    !Modifier
    public  :: ModifyExtractSurfaceCurrents
    private ::      FlowProperties
    private ::      ResidualFlowProperties
    private ::          ComputeVorticity
    private ::          ComputeKineticEnergy
    private ::          ComputePerturbationPE
    private ::          ComputeBaroclinicForce
    private ::      Read_HDF5_Residual_Hydro
    private ::      Read_HDF5_Hydro_File
    private ::      Read_HDF5_Water_File
    private ::      Read_HDF5_CommonData
    private ::      Write_HDF5_OutPut_File
    private ::      Write_HDF5_Residual_Hydro


    !Destructor
    public  :: KillExtractSurfaceCurrents                                                     
    private ::      DeAllocateMatrixes



    !Types---------------------------------------------------------------------



    private :: T_InputVar
    type       T_InputVar
        real,    dimension(:, :, :),  pointer       :: VelocityU, VelocityV
        real,    dimension(:, :, :),  pointer       :: Temperature
        real,    dimension(:, :   ),  pointer       :: WaterLevel, Bathymetry
        real,    dimension(:, :   ),  pointer       :: Latitude, Longitude
        real,    dimension(:, :   ),  pointer       :: Connection_X, Connection_Y
        integer, dimension(:, :, :),  pointer       :: WaterPoints3D
        integer, dimension(:, :, :),  pointer       :: OpenPoints3D
        integer, dimension(:, :   ),  pointer       :: WaterPoints2D

    end type  T_InputVar

    private :: T_FileName
    type       T_FileName
        character(LEN=StringLength)                 :: Geometry, Grid, InputHydroHDF5, InputWaterHDF5, OutPutHDF5
    end type  T_FileName

    
    private :: T_ExtractSurfaceCurrents
    type       T_ExtractSurfaceCurrents
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
        type (T_Options  )                          :: ComputeOptions
        type (T_FileName )                          :: FileName
        type(T_ExtractSurfaceCurrents), pointer       :: Next
    end type  T_ExtractSurfaceCurrents

    !Global Module Variables
    type (T_ExtractSurfaceCurrents), pointer          :: Me

    logical :: FirstDensityField = .true. 

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartExtractSurfaceCurrents

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                             :: STAT_CALL

        !------------------------------------------------------------------------

        nullify (Me)
        allocate(Me)

        call ConstructEnterData (Me%ObjEnterData, FileName ="ExtractSurfaceCurrents.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR10'

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
        if (STAT_CALL /= SUCCESS_) stop 'StartExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR40'


        call UnGetMap          (Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'StartExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR50'

        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'StartExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR60'


        !----------------------------------------------------------------------

    end subroutine StartExtractSurfaceCurrents
 
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
		logical										:: VariableDT

        !Begin-----------------------------------------------------------------
        

        write(*,*)'Reading instructions...'


        call GetData(Me%FileName%Grid,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GRID_FILENAME',                                    &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR20'

        if (iflag == 0)then
            write(*,*)'Must specify name of the grid file'
            stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR30'
        end if


        call GetData(Me%FileName%Geometry,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GEOMETRY_FILENAME',                                &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR40'

        if (iflag == 0)then
            write(*,*)'Must specify name of the geometry file'
            stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR50'
        end if


        Me%ComputeOptions%HydroON = .false.
        Me%ComputeOptions%WaterON = .false.

        call GetData(Me%ComputeOptions%Vorticity,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'VORTICITY',                                        &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR60'

        if (Me%ComputeOptions%Vorticity) Me%ComputeOptions%HydroON = .true.

        call GetData(Me%ComputeOptions%BaroclinicForce,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BAROCLINIC_FORCE',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR70'

        if (Me%ComputeOptions%BaroclinicForce) Me%ComputeOptions%WaterON = .true.

        call GetData(Me%ComputeOptions%KineticEnergy,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'KINETIC_ENERGY',                                   &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR80'

        if (Me%ComputeOptions%KineticEnergy) Me%ComputeOptions%HydroON = .true.

        call GetData(Me%ComputeOptions%PerturbationPE,                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POTENTIAL_ENERGY',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR100'

        if (Me%ComputeOptions%PerturbationPE) Me%ComputeOptions%WaterON = .true.

        call GetData(Me%ComputeOptions%PerturbationPE_V2,                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POTENTIAL_ENERGY_V2',                              &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR102'

if2:    if (Me%ComputeOptions%PerturbationPE_V2) then
        
            Me%ComputeOptions%WaterON = .true.

            call GetData(Me%ComputeOptions%DZ,                                          &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'POTENTIAL_ENERGY_DZ',                          &
                         default      = 10.,                                            &
                         ClientModule = 'ModuleExtractSurfaceCurrents',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR103'

        endif if2

        call GetData(Me%ComputeOptions%InitialON,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'INITIAL_ANALYSIS',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR105'

        if (Me%ComputeOptions%InitialON .and. Me%ComputeOptions%HydroON)                &
            stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR107'
       


        if (Me%ComputeOptions%HydroON) then

            call GetData(Me%FileName%InputHydroHDF5,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'HYDRODYNAMIC_FILENAME',                        &
                         ClientModule = 'ModuleExtractSurfaceCurrents',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR110'

            if (iflag == 0)then
                write(*,*)'You do not specify name of the hydrodynamic input HDF5  file'
                stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR120'
            end if

        endif

        if (Me%ComputeOptions%WaterON .and. .not. Me%ComputeOptions%InitialON) then

            call GetData(Me%FileName%InputWaterHDF5,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'WATER_FILENAME',                               &
                         ClientModule = 'ModuleExtractSurfaceCurrents',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR130'

            if (iflag == 0)then
                write(*,*)'You do not specify name of the water properties input HDF5  file'
                stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR140'
            end if

        endif


        if (Me%ComputeOptions%WaterON .or. Me%ComputeOptions%HydroON) then

            call GetData(Me%FileName%OutPutHDF5,                                        &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'OUTPUT_FILENAME',                              &
                         ClientModule = 'ModuleExtractSurfaceCurrents',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR150'

            if (iflag == 0)then
                write(*,*)'You do not specify name of the output HDF5  file'
                stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR160'
            end if

        else

            write(*,*)'You do not specify any action to the hydrodynamic analyser tool'
            stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR170'

        endif


        call GetData(Me%ComputeOptions%DensPressureCorrect,                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PRESSURE_CORRECTION',                              &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR180'

        call GetData(Me%ComputeOptions%SameLayerPressure,                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SAME_LAYER_PRESSURE',                              &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR185'


        call GetData(Me%ComputeOptions%DensityMethod,                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DENSITY_METHOD',                                   &
                     default      = UNESCOState_,                                       &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR190'

        call GetData(Me%ComputeOptions%BaroclinicMethod,                                &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BAROCLINIC_METHOD',                                &
                     default      = Leibniz2,                                           &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR200'

        call GetData(Me%ComputeOptions%BaroclinicPoliDegree,                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BAROCLINIC_POLIDEGREE',                            &
                     default      = 3,                                                  &
                     ClientModule = 'ModuleExtractSurfaceCurrents',                       &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR210'

        call null_time(Me%BeginTime  )
        call null_time(Me%EndTime    )
        call null_time(Me%CurrentTime)

		call ReadTimeKeyWords(Me%ObjEnterData, FromFile, Me%BeginTime, Me%EndTime,      &
                              Me%DT, VariableDT , "ModuleExtractSurfaceCurrents")

        Me%CurrentTime = Me%BeginTime


        if (Me%ComputeOptions%HydroON) then
            
            call GetData(Me%ComputeOptions%ResidualHydro,                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'RESIDUAL_HYDRO',                               &
                         default      = .false.,                                        &
                         ClientModule = 'ModuleExtractSurfaceCurrents',                   &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR230'
        endif

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%EndTime, Me%DT, .false., STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleExtractSurfaceCurrents - ERR240'


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
                        stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR20'

                    exit

                endif BF

            else  EB

                stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR30'

            endif EB

        enddo
        
        if (Me%InputVar%IniWaterLevelON) then

            call ComputeInitialGeometry(GeometryID       = Me%ObjGeometry,              &
                                        WaterPoints3D    = Me%InputVar%WaterPoints3D,   &
                                        SurfaceElevation = Me%InputVar%IniWaterLevel,   &
                                        ActualTime       = Me%CurrentTime,              &   
                                        STAT             = STAT_CALL )

            call UnGetMap(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR40'


            call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR45'

            !Update the open points
            call UpdateComputeFaces3D(Me%ObjMap, Me%InputVar%IniWaterLevel,             &
                                      Me%CurrentTime, STAT = STAT_CALL)      
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR50'
            
            call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR60'

            call GetWaterPoints2D(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR65'

			call GetGeometryDistances (Me%ObjGeometry, SZZ = SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR67'

			Me%InputVar%IniSZZ(:,:,:) = SZZ(:,:,:)

			call UnGetGeometry(Me%ObjGeometry, SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR69'

        else
            
            stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR70'

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
                        stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR80'

                    exit

                endif i2

            else  i1

                stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR90'

            endif i1

        enddo
        
        deallocate(Field3D)

        if (.not. Me%InputVar%IniDensON)  then

            if (.not. Me%InputVar%IniSalON)                                             &
                stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR100'

            if (.not. Me%InputVar%IniTempON)                                            &
                stop 'ReadInitialDensityField - ModuleExtractSurfaceCurrents - ERR110'

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
            stop 'ReadProperty3D - ModuleExtractSurfaceCurrents - ERR10' 




        if (IDProperty%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = IDProperty%ObjFillMatrix,           &
                                   Matrix3D       = Field3D,                            &
                                   PointsToFill3D = Me%InputVar%WaterPoints3D,          &
                                   STAT           = STAT_CALL)                          
            if (STAT_CALL /= SUCCESS_)                                                  &
                stop 'ReadProperty3D - ModuleExtractSurfaceCurrents - ERR20'

        end if


        call KillFillMatrix(IDProperty%ObjFillMatrix, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperty3D - ModuleExtractSurfaceCurrents - ERR30' 

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
            stop 'ReadProperty2D - ModuleExtractSurfaceCurrents - ERR10' 




        if (IDProperty%SolutionFromFile) then

            call ModifyFillMatrix (FillMatrixID   = IDProperty%ObjFillMatrix,       &
                                   Matrix2D       = Field2D,                        &
                                   PointsToFill2D = Me%InputVar%WaterPoints2D,      &
                                   STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                              &
                stop 'ReadProperty2D - ModuleExtractSurfaceCurrents - ERR20'

        end if



        call KillFillMatrix(IDProperty%ObjFillMatrix, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperty2D - ModuleExtractSurfaceCurrents - ERR30' 

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
            stop 'ConstructGrid - ModuleExtractSurfaceCurrents - ERR10'
        endif

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%FileName%Grid, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid - ModuleExtractSurfaceCurrents - ERR20'


        call ConstructGridData      (GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     FileName         = Me%FileName%Grid,               &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleExtractSurfaceCurrents - ERR40'

        call ConstructHorizontalMap (HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     ActualTime       = Me%BeginTime,                   & 
                                     STAT             = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleExtractSurfaceCurrents - ERR50'


        call ConstructGeometry      (GeometryID       = Me%ObjGeometry,                 &
                                     GridDataID       = Me%ObjBathymetry,               &
                                     HorizontalGridID = Me%ObjHorizontalGrid,           &
                                     HorizontalMapID  = Me%ObjHorizontalMap,            &
                                     ActualTime       = Me%BeginTime,                   &
                                     NewDomain        = Me%Filename%Geometry,           &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleExtractSurfaceCurrents - ERR60'

        call GetGeometrySize(GeometryID     = Me%ObjGeometry,                           &
                             Size           = Me%Size,                                  &
                             WorkSize       = Me%WorkSize,                              &
                             STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleExtractSurfaceCurrents - ERR70'

        call ConstructMap ( Map_ID          = Me%ObjMap,                                &
                            GeometryID      = Me%ObjGeometry,                           &
                            HorizontalMapID = Me%ObjHorizontalMap,                      &
                            TimeID          = Me%ObjTime,                               &
                            STAT            = STAT_CALL)  
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructGrid -  ModuleExtractSurfaceCurrents - ERR80'


        call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleExtractSurfaceCurrents - ERR90'

        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ConstructGrid - ModuleExtractSurfaceCurrents - ERR100'



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
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_Input_File - ModuleExtractSurfaceCurrents - ERR10'

            call HDF5SetLimits  (Me%ObjHDF5_InputHydro, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Open_HDF5_Input_File - ModuleExtractSurfaceCurrents - ERR10'

            PropertyName = "Time"

            do m=1,Me%WaterInstants

                TimePtr => AuxTime1

                call HDF5ReadData(Me%ObjHDF5_InputWater, "/Time",                                             &
                                  trim(PropertyName),                                           &
                                  Array1D      = TimePtr,                                       &
                                  OutputNumber = m, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstrucGrid - ModuleExtractSurfaceCurrents - ERR09a'

                TimePtr => AuxTime2

                call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Time",                                             &
                                  trim(PropertyName),                                           &
                                  Array1D      = TimePtr,                                       &
                                  OutputNumber = m, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File - ModuleExtractSurfaceCurrents - ERR09b'

                do n=1,6
                    
                    if (AuxTime1(n) /= AuxTime2(n)) then

                        stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR10'

                    endif

                enddo

    !            if (Me%WaterInstants /= Me%HydroInstants) then

    !                stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR10'

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
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR20'

            if (.not. Exist) then
                write(*,*) 'There is not vertical velocity fields a null value will be assumed'
                write(*,*) 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - WARN10'
            endif

            PropertyName = GetPropertyName(VelocityU_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR30'

            if (.not. Exist) then
                write(*,*) 'There is not zonal velocity fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR40'
            endif

            PropertyName = GetPropertyName(VelocityV_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR50'

            if (.not. Exist) then
                write(*,*) 'There is not meridional velocity fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR60'
            endif

        endif


        if (Me%ComputeOptions%WaterON) then

            PropertyName = GetPropertyName(Temperature_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputWater, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR70'

            if (.not. Exist) then
                write(*,*) 'There is no temperature fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR80'
            endif

            PropertyName = GetPropertyName(Salinity_)
            GroupName = "/Results/"//trim(PropertyName)
            
            call GetHDF5GroupExist (Me%ObjHDF5_InputWater, GroupName, Exist, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR90'

            if (.not. Exist) then
                write(*,*) 'There is not salinity fields to be read'                
                stop 'Open_HDF5_Input_File -  ModuleExtractSurfaceCurrents - ERR100'
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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR10'
        

        call GetHDF5GroupNumberOfItems(Me%ObjHDF5_InputWater, "/Time",                  &
                                       Me%WaterInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR20'

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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR10'
        

        call GetHDF5GroupNumberOfItems(Me%ObjHDF5_InputHydro, "/Time",                  &
                                       Me%HydroInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR20'
        

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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR20'
        
        
        call HDF5SetLimits  (Me%ObjHDF5_OutPut, Me%WorkSize%ILB, Me%WorkSize%IUB,       &
                                                Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR30'

        call GetGridData(Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR40'
        
        call HDF5WriteData   (Me%ObjHDF5_OutPut, "/Grid", "Bathymetry", "-",            &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR70'

        call HDF5SetLimits  (Me%ObjHDF5_OutPut, Me%WorkSize%ILB, Me%WorkSize%IUB,       &
                                                Me%WorkSize%JLB, Me%WorkSize%JUB,       &
                                                Me%WorkSize%KLB, Me%WorkSize%KUB, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR80'

        call HDF5WriteData   (Me%ObjHDF5_OutPut, "/Grid", "WaterPoints3D", "-",         &
                              Array3D = Me%InputVar%WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR90'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR100'


        call UnGetGridData (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR120'


    end subroutine Open_HDF5_OutPut_File

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

    subroutine ModifyExtractSurfaceCurrents

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, STAT_CALL
        type (T_Time)                               :: WaterInstant, HydroInstant

        !----------------------------------------------------------------------


        call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR10'

        do i = 1, Me%OutInstants
            
NI:         if (Me%ComputeOptions%InitialON) then

                Me%InputVar%Density     => Me%InputVar%IniDensity
                Me%InputVar%SigmaDens   => Me%InputVar%IniSigmaDens
                Me%InputVar%Temperature => Me%InputVar%IniTemperature
                Me%InputVar%Salinity    => Me%InputVar%IniSalinity
                
				call GetGeometryDistances (Me%ObjGeometry, SZZ = Me%InputVar%SZZ, STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleWaterProperties - ERR30'

				call GetGeometryDistances (Me%ObjGeometry, ZCellCenter = Me%InputVar%ZCellCenter, STAT = STAT_CALL)
				if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleWaterProperties - ERR31'

                call GetOpenPoints3D(Me%ObjMap, Me%InputVar%OpenPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR20'

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

                        stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR40'

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
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR50'


                !Update the open points
                call UpdateComputeFaces3D(Me%ObjMap, Me%InputVar%WaterLevel,                &
                                          Me%CurrentTime, STAT = STAT_CALL)      
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR60'

                call GetWaterPoints3D(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR70'
          
            endif NI

            call FlowProperties

            call Write_HDF5_OutPut_File(i)

        enddo

        if (Me%ComputeOptions%HydroON .and. Me%ComputeOptions%ResidualHydro) then

            call Read_HDF5_Residual_Hydro

            call ResidualFlowProperties
    
            call Write_HDF5_Residual_Hydro

        endif


        if (Me%ComputeOptions%InitialON) then

            call UnGetMap(Me%ObjMap, Me%InputVar%OpenPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR80'

            call UnGetGeometry (Me%ObjGeometry, Me%InputVar%SZZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleWaterProperties - ERR90'

            call UnGetGeometry (Me%ObjGeometry, Me%InputVar%ZCellCenter, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleWaterProperties - ERR91'

        endif


        call UnGetMap(Me%ObjMap, Me%InputVar%WaterPoints3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR100'


    end subroutine ModifyExtractSurfaceCurrents

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
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR10'

        PropertyName = GetPropertyName(Density_)

        call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName), &
                                Me%DensityExist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR20'


        if (Me%DensityExist)  then

            call HDF5ReadData(Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName),       &
                              trim(PropertyName),                                           &
                              Array3D      = Me%InputVar%Density,                           &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR30'

            Me%InputVar%SigmaDens(:,:,:) = Me%InputVar%Density - SigmaDensityReference

        endif



        PropertyName = GetPropertyName(Salinity_)

        call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName), &
                                SalinityExist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR40'

        if (SalinityExist) then

            !read field
            call HDF5ReadData(Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName),       &
                              trim(PropertyName),                                           &
                              Array3D      = Me%InputVar%Salinity,                          &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR60'

        endif

        PropertyName = GetPropertyName(Temperature_)

        call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName), &
                                TemperatureExist, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ )stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR70'

        if (TemperatureExist) then

            call HDF5ReadData(Me%ObjHDF5_InputWater, "/Results/"//trim(PropertyName),       &
                              trim(PropertyName),                                           &
                              Array3D      = Me%InputVar%Temperature,                       &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Water_File - ModuleExtractSurfaceCurrents - ERR90'

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
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR10'

        PropertyName = GetPropertyName(VelocityU_)

        !read field
        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Results/"//trim(PropertyName),       &
                          trim(PropertyName),                                           &
                          Array3D      = Me%InputVar%VelocityU,                         &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR20'

        PropertyName = GetPropertyName(VelocityV_)

        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Results/"//trim(PropertyName),       &
                          trim(PropertyName),                                           &
                          Array3D      = Me%InputVar%VelocityV,                         &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR30'

        PropertyName = GetPropertyName(VelocityW_)

        GroupName = "/Results/"//trim(PropertyName)

        call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, GroupName, Exist, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR40'

        if (Exist) then

            call HDF5ReadData(Me%ObjHDF5_InputHydro, GroupName,                         &
                              trim(PropertyName),                                       &
                              Array3D      = Me%InputVar%VelocityW,                     &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Hydro_File - ModuleExtractSurfaceCurrents - ERR50'
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
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR10'

       !read field
        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Velocity/X",                &
                          "Vel. X",                                                     &
                          Array3D      = Me%InputVar%ResidualVelU,                      &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR20'

        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Velocity/Y",                &
                          "Vel. Y",                                                     &
                          Array3D      = Me%InputVar%ResidualVelV,                      &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR30'


        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Flux/X",                    &
                          "Flux X",                                                     &
                          Array3D      = Me%InputVar%ResidualFluxU,                     &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR40'

        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Flux/Y",                    &
                          "Flux Y",                                                     &
                          Array3D      = Me%InputVar%ResidualFluxV,                     &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR50'


        call HDF5ReadData(Me%ObjHDF5_InputHydro, "/Residual/Waterlevel",                &
                          "Vel. Z",                                                     &
                          Array2D      = Me%InputVar%ResidualWaterLevel,                &
                          STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR70'


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
        if (STAT_CALL /= SUCCESS_) stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR10'


        PropertyName = "Time"

        TimePtr => AuxTime

        call HDF5ReadData(ObjHDF5, "/Time",                                             &
                          trim(PropertyName),                                           &
                          Array1D      = TimePtr,                                       &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR20'


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
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR30'

        PropertyName = "OpenPoints"
 
        if (Me%ComputeOptions%WaterON) then
        
            call GetHDF5GroupExist (Me%ObjHDF5_InputWater, "/Grid/"//trim(PropertyName),&
                                    GroupExistWater, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR35'

        endif

        if (Me%ComputeOptions%HydroON) then
        
            call GetHDF5GroupExist (Me%ObjHDF5_InputHydro, "/Grid/"//trim(PropertyName),&
                                    GroupExistHydro, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR37'

        endif


        if (GroupExistWater .or. GroupExistHydro) then

        !read field
        call HDF5ReadData(ObjHDF5, "/Grid/"//trim(PropertyName),                        &
                          trim(PropertyName),                                           &
                          Array3D      = Me%InputVar%OpenPoints3D,                      &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR40'

        endif

        call HDF5SetLimits (ObjHDF5,                                                    &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB-1,                                          &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR50'

        PropertyName = "VerticalZ"

        call HDF5ReadData(ObjHDF5, "/Grid/"//trim(PropertyName),                        &
                          "Vertical",                                                   &
                          Array3D      = Me%InputVar%SZZ,                               &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Read_HDF5_CommonData - ModuleExtractSurfaceCurrents - ERR60'


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
        if (STAT_CALL /= SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR10'


        PropertyName = "Time"

        TimePtr => AuxTime

        call ExtractDate (Me%CurrentTime,                                               &
                         AuxTime(1), AuxTime(2), AuxTime(3), AuxTime(4), AuxTime(5), AuxTime(6))

        call HDF5WriteData(Me%ObjHDF5_Output, "/Time",                                  &
                          trim(PropertyName),                                           &
                          Units        = 'YYYY/MM/DD',                                  &
                          Array1D      = TimePtr,                                       &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR20'


        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR30'

        PropertyName = "OpenPoints"

        !write field
        call HDF5WriteData(Me%ObjHDF5_Output, "/Grid/"//trim(PropertyName),             &
                          trim(PropertyName),                                           &
                          Units        = '-',                                           &
                          Array3D      = Me%InputVar%OpenPoints3D,                      &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR40'

        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB-1,                                          &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR50'

        PropertyName = "VerticalZ"

        call HDF5WriteData(Me%ObjHDF5_Output, "/Grid/"//trim(PropertyName),             &
                          "Vertical",                                                   &
                          Units        = 'm',                                           &
                          Array3D      = Me%InputVar%SZZ,                               &
                          OutputNumber = i, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR60'


        call HDF5SetLimits (Me%ObjHDF5_Output,                                          &
                            Me%WorkSize%ILB,                                            &
                            Me%WorkSize%IUB,                                            &
                            Me%WorkSize%JLB,                                            &
                            Me%WorkSize%JUB,                                            &
                            Me%WorkSize%KLB,                                            &
                            Me%WorkSize%KUB,                                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR70'

iV:     if (Me%ComputeOptions%Vorticity) then

            PropertyName = GetPropertyName(Vorticity_)

            !write vorticity field
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 's-1',                                     &
                              Array3D      = Me%OutputVar%Vorticity,                    &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR80'

        endif iV

iP:     if (Me%ComputeOptions%PerturbationPE) then

            PropertyName = GetPropertyName(PerturbationPE_)

            !write potential energy field 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%PerturbationPE,               &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR90'

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
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR90'

        endif iP2

iK:     if (Me%ComputeOptions%KineticEnergy) then

            PropertyName = GetPropertyName(KineticEnergy_)

            !write barotropic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%KineticEnergy,                &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR100'

            PropertyName = GetPropertyName(BaroclinicKE_)            
            
            !write baroclinic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%BaroclinicKE,                 &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR110'

            PropertyName = GetPropertyName(BarotropicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%BarotropicU,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR120'

            PropertyName = GetPropertyName(BarotropicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%BarotropicV,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR130'

            PropertyName = GetPropertyName(BaroclinicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%BaroclinicU,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR140'

            PropertyName = GetPropertyName(BaroclinicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%BaroclinicV,                  &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR150'


        endif iK

iB:     if (Me%ComputeOptions%BaroclinicForce) then

            PropertyName = GetPropertyName(BaroclinicForceX_)

            !write baroclinic force in the X direction
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'kg/m3',                                     &
                              Array3D      = Me%OutputVar%BaroclinicForceX,             &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR160'

            PropertyName = GetPropertyName(BaroclinicForceY_)

            !write baroclinic force in the Y direction
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'kg/m3',                                   &
                              Array3D      = Me%OutputVar%BaroclinicForceY,             &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR170'

        endif iB

iB2:    if (Me%ComputeOptions%BaroclinicForce .or. Me%ComputeOptions%PerturbationPE .or. Me%ComputeOptions%PerturbationPE_V2) then


            PropertyName='NN'

            !write the square of Brunt-Vaisala frequency 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 's-2',                                     &
                              Array3D      = Me%OutputVar%NN,                           &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR180'


            PropertyName = GetPropertyName(Density_)

            !write the density 
            call HDF5WriteData(Me%ObjHDF5_Output, "/Results/"//trim(PropertyName),      &
                              trim(PropertyName),                                       &
                              Units        = 'kg/m3',                                   &
                              Array3D      = Me%InputVar%Density,                       &
                              OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR190'
            
        endif iB2



        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_HDF5_OutPut_File - ModuleExtractSurfaceCurrents - ERR200'


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
        if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR10'

iV:     if (Me%ComputeOptions%Vorticity) then

            PropertyName = GetPropertyName(Vorticity_)

            !write vorticity field
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 's-1',                                     &
                              Array3D      = Me%OutputVar%ResVorticity,                 &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR20'

        endif iV


iK:     if (Me%ComputeOptions%KineticEnergy) then

            PropertyName = GetPropertyName(KineticEnergy_)

            !write barotropic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%ResKineticEnergy,             &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR100'

            PropertyName = GetPropertyName(BaroclinicKE_)            
            
            !write baroclinic kinetic energy
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm2/s2',                                   &
                              Array2D      = Me%OutputVar%ResBaroclinicKE,              &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR110'

            PropertyName = GetPropertyName(BarotropicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%ResBarotropicU,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR120'

            PropertyName = GetPropertyName(BarotropicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array2D      = Me%OutputVar%ResBarotropicV,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR130'

            PropertyName = GetPropertyName(BaroclinicVelocityU_)

            !write barotropic velocity U
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%ResBaroclinicU,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR140'

            PropertyName = GetPropertyName(BaroclinicVelocityV_)

            !write barotropic velocity V
            call HDF5WriteData(Me%ObjHDF5_Output, "/Residual/"//trim(PropertyName),     &
                              trim(PropertyName),                                       &
                              Units        = 'm/s',                                     &
                              Array3D      = Me%OutputVar%ResBaroclinicV,               &
                              STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR150'


        endif iK



        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5_OutPut, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Write_HDF5_Residual_Hydro - ModuleExtractSurfaceCurrents - ERR180'


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
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleExtractSurfaceCurrents ERR10'


        !Gets WaterPoints2D
        call GetWaterPoints2D(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleExtractSurfaceCurrents ERR20'        


        !Module - ModuleGeometry
        !3D Geometry properties

        call GetGeometryKFloor(Me%ObjGeometry,                                          &
                               Z = Me%InputVar%KFloor_Z,                                &
                               U = Me%InputVar%KFloor_U,                                &
                               V = Me%InputVar%KFloor_V,                                &
                               STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleExtractSurfaceCurrents ERR30'        


        call GetGeometryVolumes(Me%ObjGeometry,                                         &
                                VolumeZ = Me%InputVar%Volume_Z,                         &
                                VolumeU = Me%InputVar%Volume_U,                         &
                                VolumeV = Me%InputVar%Volume_V,                         &
                                VolumeW = Me%InputVar%Volume_W,                         &
                                STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleExtractSurfaceCurrents ERR40'        


        call GetGeometryDistances(Me%ObjGeometry,                                       &
                                  DWZ = Me%InputVar%DWZ,                                &
                                  DUZ = Me%InputVar%DUZ,                                &
                                  DVZ = Me%InputVar%DVZ,                                &
                                  STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleExtractSurfaceCurrents ERR50' 


        call GetGeometryWaterColumn(Me%ObjGeometry,                                     &
                                    WaterColumn = Me%InputVar%WaterColumn,              &
                                    WaterColumnU= Me%InputVar%WaterColumnU,             &
                                    WaterColumnV= Me%InputVar%WaterColumnV,             &
                                    STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_) stop 'ReadLockInputVar - ModuleExtractSurfaceCurrents ERR60' 



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
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR10'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, Me%InputVar%DZY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR20'


        !UnGets WaterPoints2D
        call UnGetHorizontalMap(Me%ObjHorizontalMap, Me%InputVar%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR30'        


        !Module - ModuleGeometry
        !3D Geometry properties
        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%KFloor_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR40'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%KFloor_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR50'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%KFloor_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR60'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_Z, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR70'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_U, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR80'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_V, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR90'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%Volume_W, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR100'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%DWZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR110'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%DUZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR120'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%DVZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR130'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%WaterColumn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR140'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%WaterColumnU, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR150'        

        call UnGetGeometry(Me%ObjGeometry, Me%InputVar%WaterColumnV, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockInputVar - ModuleExtractSurfaceCurrents ERR160'        


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
                        
                        Aux = Me%InputVar%DWZ(i, j, k)/Me%InputVar%WaterColumn(i,j)
    
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
                          
                            Aux = Me%InputVar%DWZ(i, j, k)/Me%InputVar%WaterColumn(i,j)

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
        real,    dimension(:    ), pointer  :: RefProf,TprofileRef, SprofileRef
        real,    dimension(:    ), pointer  :: ProfInst,TprofileInst, SprofileInst
        real(8)                             :: dw, dr, PPE, TotalDepth, D2, D1
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
        if (STAT_CALL /= SUCCESS_)stop 'ComputePerturbationPE_V2 - ModuleExtractSurfaceCurrents - ERR40'

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

                            D1 = SigmaUNESCO (Temperature(i, j, k  ), Salinity(i, j, k))
                            D2 = SigmaUNESCO (Taux                  , Saux)

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
        if (STAT_CALL /= SUCCESS_)stop 'ComputePerturbationPE_V2 - ModuleExtractSurfaceCurrents - ERR120'


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

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillExtractSurfaceCurrents

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                         :: STAT_CALL

        !------------------------------------------------------------------------
           
        call DeAllocateMatrixes

        call KillMap(Me%ObjMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR10'

        call KillGeometry(Me%ObjGeometry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR20'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR30'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR40'

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR50'
        
        call KillHDF5(Me%ObjHDF5_Output, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR60'

NI:     if (.not. Me%ComputeOptions%InitialON) then

            if (Me%ComputeOptions%WaterON) then
                call KillHDF5(Me%ObjHDF5_InputWater, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR80'
            endif

            if (Me%ComputeOptions%HydroON) then
                call KillHDF5(Me%ObjHDF5_InputHydro, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'KillExtractSurfaceCurrents - ModuleExtractSurfaceCurrents - ERR70'
            endif

        endif NI

        deallocate(Me)

        !------------------------------------------------------------------------

    end subroutine KillExtractSurfaceCurrents
        

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
end module ModuleExtractSurfaceCurrents









